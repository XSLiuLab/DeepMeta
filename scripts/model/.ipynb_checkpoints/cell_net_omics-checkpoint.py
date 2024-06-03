import numpy as np
import pandas as pd
import torch
import torch.nn as nn
import networkx as nx
import torch_geometric
from torch_geometric.utils.convert import from_networkx
from torch_geometric.data import Dataset, Data
import os
from multiprocessing import Pool
from functools import partial

def convert_ones_and_sample_zeros(matrix):
    # Iterate over each row
    for row in matrix:
        # Find the indices of the 1s and 0s in the row
        ones_indices = (row == 1).nonzero(as_tuple=True)[0]
        zeros_indices = (row == 0).nonzero(as_tuple=True)[0]
        # Convert all 1s to 0s
        row[ones_indices] = 0
        # Randomly sample the same number of 0s as there were 1s
        sample_size = len(ones_indices)
        sampled_indices = zeros_indices[torch.randperm(len(zeros_indices))[:sample_size]]
        # Set the sampled indices to 1
        row[sampled_indices] = 1
    return matrix

def get_cell_net(edge_list,feat,exp,mut,cnv,data_type,random=False):
    feat = feat.drop(columns=["gene","is_exp"])
    if (data_type == "val") and ("is_dep" in feat.columns):
        feat = feat.drop(columns=["is_dep"])

    cell_net_pre = nx.from_pandas_edgelist(edge_list, source="from",target="to")
    cell_net = cell_net_pre.subgraph(list(feat.id))
    cell_net = nx.Graph(cell_net)

    cell_node_feature = feat.set_index('id').to_dict('index')
    nx.set_node_attributes(cell_net, cell_node_feature)

    cell_net = nx.relabel.convert_node_labels_to_integers(cell_net, label_attribute="nodename")##保留节点基因名
    node_name = nx.get_node_attributes(cell_net, "nodename")
    #node_name = list(node_name.values())
    ##trans to pyg data
    attr_name = feat.columns
    attr_name = list(attr_name)
    attr_name.remove("id")
    
    cell_pyg = from_networkx(cell_net,group_node_attrs = attr_name)
    
    if data_type == "val":
        cell_pyg.exp = exp.repeat(cell_pyg.x.shape[0], 1)
        cell_pyg.mut = mut.repeat(cell_pyg.x.shape[0], 1)
        cell_pyg.cnv = cnv.repeat(cell_pyg.x.shape[0], 1)
    else:
        indices = torch.tensor([0])
        cell_y = torch.flatten(torch.index_select(cell_pyg.x, 1, indices))
        non_na_idx = (torch.where(torch.isnan(cell_y), 1.0, 0.0) == torch.tensor(0)).nonzero().flatten()
        cell_y_non_na = cell_y[non_na_idx]

        cell_pyg.y = cell_y_non_na
        cell_pyg.y_index = non_na_idx
        indices = torch.tensor(list(range(1,3248)))
        cell_x = torch.index_select(cell_pyg.x, 1, indices)
        if random:
            x_copy = cell_x.clone()
            x_random = convert_ones_and_sample_zeros(x_copy)
            cell_pyg.x = x_random
        else:
            cell_pyg.x = cell_x
        cell_pyg.exp = exp.repeat(cell_x.shape[0], 1)
        cell_pyg.mut = mut.repeat(cell_x.shape[0], 1)
        cell_pyg.cnv = cnv.repeat(cell_x.shape[0], 1)
    
    return cell_pyg

class cellNetDataset(Dataset):
    def __init__(self, root, filename, exp, mut, cnv, data_type, net_path, cores, random = False,transform=None, pre_transform=None):
        """
        root = Where the dataset should be stored. This folder is split
        into raw_dir (downloaded dataset) and processed_dir (processed data).
        filename file must contain two columns: cell, cell_index
        """
        self.filename = filename
        self.exp = exp
        self.mut = mut
        self.cnv = cnv
        self.data_type = data_type
        self.net_path = net_path
        self.cores = cores
        self.random = random
        super(cellNetDataset, self).__init__(root, transform, pre_transform)
        
    @property
    def raw_file_names(self):
        """ If this file exists in raw_dir, the download is not triggered.
            (The download func. is not implemented here)  
        """
        return self.filename

    @property
    def processed_file_names(self):
        """ If these files are found in processed_dir, processing is skipped"""
        self.data = pd.read_csv(self.raw_paths[0]).reset_index()
        
        return [f'data_{i}.pt' for i in list(self.data.index)]

    def download(self):
        pass##不需要下载

    def process(self):
        n = self.cores# Any number of threads
        with Pool(n) as p:
            p.map(partial(self._save_data,
                          exp = self.exp,
                          mut = self.mut,
                          cnv = self.cnv,
                          net_path = self.net_path,
                          data_type = self.data_type,
                          random = self.random), 
                  self.data.to_dict('records'))

    def _save_data(self, key_value, exp, mut, cnv, net_path, data_type, random):
        cell = key_value["cell"]
        index = key_value["cell_index"]
        dt = pd.read_csv(net_path+cell+".txt",sep="\t")
        dt_fea = pd.read_csv(net_path+cell+"_feat.txt",sep="\t")
        dt_exp = torch.tensor(exp.loc[exp["cell"]==cell].drop(columns=["cell"]).values.flatten())
        dt_mut = torch.tensor(mut.loc[mut["cell"]==cell].drop(columns=["cell"]).values.flatten())
        dt_cnv = torch.tensor(cnv.loc[cnv["cell"]==cell].drop(columns=["cell"]).values.flatten())
        pyg_net = get_cell_net(dt,dt_fea,dt_exp,dt_mut,dt_cnv,data_type,random=random)
        pyg_net.cell = cell
        torch.save(pyg_net, os.path.join(self.processed_dir, f"data_{index}.pt"))

    def len(self):
        return self.data.shape[0]

    def get(self, idx):
        """ - Equivalent to __getitem__ in pytorch
            - Is not needed for PyG's InMemoryDataset
        """
        data = torch.load(os.path.join(self.processed_dir, f'data_{idx}.pt')) 
        return data