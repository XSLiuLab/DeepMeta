import numpy as np
import pandas as pd
import torch
import torch.nn as nn  
import pickle
from tqdm import tqdm
import networkx as nx
import pickle
import os.path
import os
import argparse
import shutil
import warnings
warnings.filterwarnings('ignore', category=UserWarning, message='TypedStorage is deprecated')

from cell_net import cellNetDataset
from torch_geometric.loader import DataLoader
from torch_geometric.nn import GATv2Conv
###参数
parser = argparse.ArgumentParser()
parser.add_argument("-e", "--exp", help="file path of gene expression (tumor VS normal)")
parser.add_argument("-g", "--graph", help="dir path of pyg data")
parser.add_argument("-c", "--cell", help="file name of predicted cell (csv)")
parser.add_argument("-n", "--net", help="dir path of sample enzyme net")
parser.add_argument("-t", "--cores", help="multiprocess cores for creating pyg data")
parser.add_argument("-m", "--model", help="file path of model .pt file")
parser.add_argument("-d", "--datatype", help="test or val")
parser.add_argument("-b", "--batchsize", help="batch size")
parser.add_argument("-o", "--output", help="output file name")
parser.add_argument("-a", "--gpu", help="cpu or gpu")
args = parser.parse_args()
exp_file = args.exp
pyg_dir = args.graph
pre_file = args.cell
net_dir = args.net
cores = args.cores
model_file = args.model
out_file = args.output
pre_type = args.datatype
batch_size = int(args.batchsize)
device = args.gpu
###create dir
pyg_dir_raw = pyg_dir + "/raw/"
pyg_dir_processed = pyg_dir + "/processed/"
if not os.path.exists(pyg_dir_raw):
    os.makedirs(pyg_dir_raw)
if not os.path.exists(pyg_dir_processed):
    os.makedirs(pyg_dir_processed)

train_genes = pd.read_csv("/home/wt/meta_target/data/train_genes.csv")
###mv pre_file to raw dir
if not os.path.exists(pyg_dir_raw + os.path.basename(pre_file)):
    shutil.copy(pre_file, pyg_dir_raw)

class Net(torch.nn.Module):
    def __init__(self):
        super().__init__()
        self.conv1 = GATv2Conv(3247, 512, heads=3)
        self.conv2 = GATv2Conv(3 * 512, 512, heads=3)
        self.conv3 = GATv2Conv(3 * 512, 512, heads=3)
        self.lin1 = torch.nn.Linear(3 * 512 + 512, 1024)
        self.lin2 = torch.nn.Linear(1024, 512)
        self.lin3 = torch.nn.Linear(512, 1)
        
        self.encoder = torch.nn.Sequential(
            torch.nn.Linear(7993, 4000),
            torch.nn.ReLU(),
            torch.nn.Linear(4000, 1500),
            torch.nn.ReLU(),
            torch.nn.Linear(1500, 512)
        )

    def forward(self, x, edge_index, exp):
        x = torch.relu(self.conv1(x, edge_index))
        x = torch.relu(self.conv2(x, edge_index))
        x = torch.relu(self.conv3(x, edge_index))
        cell_end = self.encoder(exp)
        
        x = torch.cat((x, cell_end),1)
        x = torch.relu(self.lin1(x))
        x = torch.relu(self.lin2(x))
        x = self.lin3(x)
        return x
    
exp = pd.read_csv(exp_file)
cell_net = cellNetDataset(root = pyg_dir,
                          filename = os.path.basename(pre_file),
                          exp = exp,
                          data_type = pre_type, 
                          net_path = net_dir, 
                          cores = int(cores))
data_loader = DataLoader(cell_net, batch_size = batch_size, shuffle = False)

model = Net().to(device)
model.load_state_dict(torch.load(model_file, map_location=device))

model.eval()

for i, data in enumerate(tqdm(data_loader)):
    data = data.to(device)
    out_scale = model(data.x.float(), data.edge_index, data.exp.float())
    preds = np.rint(torch.sigmoid(out_scale).cpu().detach().numpy().flatten())
    preds_raw = torch.sigmoid(out_scale).cpu().detach().numpy().flatten()
    gene_name = np.array([i for batch in data.nodename for i in batch])
    if pre_type == "test":
        ys = data.y.cpu().detach().numpy()
        dt = {'preds':preds[data.y_index],'preds_raw':preds_raw[data.y_index],'gene_name':gene_name[data.y_index],"label":ys}
    else:
        dt = {'preds':preds,'preds_raw':preds_raw,'gene_name':gene_name}
    if i == 0:
        df = pd.DataFrame(dt)
        df["cell"] = data.cell[0]
    else:
        df1 = pd.DataFrame(dt)
        df1["cell"] = data.cell[0]
        df = pd.concat([df, df1])

df = df[(df.gene_name.isin(train_genes.id))]
df.to_csv(out_file)