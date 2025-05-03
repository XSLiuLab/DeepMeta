import numpy as np
import pandas as pd
import torch
import torch.nn as nn
from torch.utils.data import Dataset, DataLoader , random_split,SubsetRandomSampler
from sklearn.model_selection import train_test_split
from sklearn.metrics import mean_squared_error, r2_score
import pickle
from tqdm import tqdm
import networkx as nx
from networkx.algorithms.link_analysis import pagerank
import pickle
from sklearn.metrics import confusion_matrix, f1_score, accuracy_score, precision_score, recall_score, roc_auc_score
import warnings
warnings.filterwarnings('ignore', category=UserWarning, message='TypedStorage is deprecated')

exp = pd.read_csv("/home/wt/DeepMeta/autodl_143/cell_gene_exp_vs_normal_filter.csv")

from cell_net import cellNetDataset
cell_net = cellNetDataset(root="/home/wt/DeepMeta/autodl_143/cell_net_filter_exp/",
                          filename = "train_cell_info.csv",exp = exp,
                          data_type = "train", 
                          net_path = "/home/wt/DeepMeta/autodl_143/cell_enzyme_net/", 
                          cores = 30)

from sklearn.model_selection import KFold
import random
splits = KFold(n_splits=5,shuffle=True,random_state=2023052701)

from tqdm import tqdm

def train():
    model.train()

    total_loss = 0
    for _, data in enumerate(tqdm(train_loader)):
        data = data.to(device)
        optimizer.zero_grad()
        pred = model(data.x.float(), data.edge_index, data.exp.float())
        pred = pred[data.y_index]
        loss = loss_op(
            pred,
            data.y.float().reshape(-1, 1),
        )
        total_loss += loss.item() * data.num_graphs
        loss.backward()
        optimizer.step()
    return total_loss / len(train_loader.dataset)


@torch.no_grad()
def test():
    model.eval()

    ys, preds, preds_raw, genes = [], [], [],[]
    total_loss = 0
    for _, data in enumerate(tqdm(test_loader)):
        data = data.to(device)
        ys.append(data.y.cpu().detach().numpy())
        out = model(data.x.float(), data.edge_index, data.exp.float())
        out = out[data.y_index]
        loss = loss_op(
            out,
            data.y.float().reshape(-1, 1),
        )
        total_loss += loss.item() * data.num_graphs
        preds.append(np.rint(torch.sigmoid(out).cpu().detach().numpy()))
        preds_raw.append(torch.sigmoid(out).cpu().detach().numpy())
        gene_name = np.array([i for batch in data.nodename for i in batch])
        genes.append(gene_name[data.y_index.cpu()])

    all_preds = np.concatenate(preds).ravel()
    all_labels = np.concatenate(ys).ravel()
    all_preds_raw = np.concatenate(preds_raw).ravel()
    all_genes = np.concatenate(genes).ravel()
    
    res = pd.DataFrame({"preds":all_preds,"preds_raw":all_preds_raw,"label":all_labels,"genes":all_genes})
    calculate_metrics(all_preds, all_labels, epoch, "test")
    return total_loss / len(test_loader.dataset), res


def calculate_metrics(y_pred, y_true, epoch, type):
    print(f"\n Confusion matrix: \n {confusion_matrix(y_true, y_pred)}")
    print(f"F1 Score: {f1_score(y_true, y_pred)}")
    print(f"Accuracy: {accuracy_score(y_true, y_pred)}")
    prec = precision_score(y_true, y_pred)
    rec = recall_score(y_true, y_pred)
    print(f"Precision: {prec}")
    print(f"Recall: {rec}")
    try:
        roc = roc_auc_score(y_true, y_pred)
        print(f"ROC AUC: {roc}")
    except:
        print(f"ROC AUC: notdefined")

from torch_geometric.nn import GATv2Conv

class Net(torch.nn.Module):
    def __init__(self, h_d, h_n):
        super().__init__()
        ##h_d hiddle dim
        ##h_n heads number
        self.conv1 = GATv2Conv(3247, h_d, heads=h_n)
        self.conv2 = GATv2Conv(h_n * h_d, h_d, heads=h_n)
        self.conv3 = GATv2Conv(h_n * h_d, h_d, heads=h_n)
        self.lin1 = torch.nn.Linear(h_n * h_d + h_d, int((h_n * h_d + h_d)/2))
        self.lin2 = torch.nn.Linear(int((h_n * h_d + h_d)/2), int(int((h_n * h_d + h_d)/2)/2))
        self.lin3 = torch.nn.Linear(int(int((h_n * h_d + h_d)/2)/2), 1)
        self.relu = torch.nn.ReLU()
        self.encoder = torch.nn.Sequential(
            torch.nn.Linear(7993, 4000),
            torch.nn.ReLU(),
            torch.nn.Linear(4000, 1500),
            torch.nn.ReLU(),
            torch.nn.Linear(1500, h_d)
        )

    def forward(self, x, edge_index, exp):
        x = self.relu(self.conv1(x, edge_index))
        x = self.relu(self.conv2(x, edge_index))
        x = self.relu(self.conv3(x, edge_index))
        cell_end = self.encoder(exp)
        
        x = torch.cat((x, cell_end),1)
        x = self.relu(self.lin1(x))
        x = self.relu(self.lin2(x))
        x = self.lin3(x)
        return x

device = torch.device("cuda" if torch.cuda.is_available() else "cpu")
loss_op = torch.nn.BCEWithLogitsLoss()

hds = [256, 512]
hns = [2, 3, 4]
lrs = [0.01,0.001,0.0001]

from torch_geometric.loader import DataLoader
for i in hds:
    for j in hns:
        for lr in lrs:
            for fold, (train_idx,val_idx) in enumerate(splits.split(np.arange(len(cell_net)))):
                print('Fold {}'.format(fold + 1))
                train_sampler = SubsetRandomSampler(train_idx)
                test_sampler = SubsetRandomSampler(val_idx)
                train_loader = DataLoader(cell_net, batch_size=3, sampler=train_sampler, num_workers=10)
                test_loader = DataLoader(cell_net, batch_size=3, sampler=test_sampler, num_workers=10)
                
                model = Net(i, j).to(device)
                optimizer = torch.optim.Adam(model.parameters(), lr=lr)
                
                for epoch in range(15):
                    loss_train = train()
                    print(f"Epoch: {epoch:03d}, Train Loss: {loss_train:.4f}")
                    
                loss_test, test_pre = test()
                print(f"Epoch: {epoch:03d}, Test Loss: {loss_test:.4f}")
                test_pre.to_csv(f"/home/wt/DeepMeta/autodl_143/cv/hd_{i}_hn_{j}_lr_{lr}_fold_"+str(fold)+".csv")


###analysis in R 
#########grid search
library(dplyr)
library(yardstick)
library(ggprism)
library(ggplot2)
library(ggpubr)

get_roc <- function(dt){
  dt <- dt %>%
    mutate(truth = ifelse(label == 1, "Class1","Class2"),
           pred_label = ifelse(preds == 1, "Class1","Class2"))
  dt$truth <- factor(dt$truth)
  dt$pred_label <- factor(dt$pred_label)

  f1 <- try(f_meas(dt,truth,pred_label)[".estimate"] %>%
              unlist() %>% unname() %>% round(.,2),silent = TRUE)
  f1 <- ifelse('try-error' %in% class(f1),NA,f1)
  roc <-  roc_auc(dt, truth, preds_raw)[".estimate"] %>%
    unlist() %>% unname() %>% round(.,2)
  return(c(roc,f1))
}

all_cv <- c()
for (i in c(256, 512)){
  for (j in c(2,3,4)){
    for (z in c("0.01","0.001","0.0001")){
      all_cv <- append(all_cv,paste0("hd_",i,"_hn_",j,"_lr_",z))
    }
  }
}

cv_res <- vector("list",length(all_cv))
for (i in seq_along(cv_res)){
  sub_res <- vector("list",5)
  for (j in 0:4){
    pre <- read.csv(paste0("autodl_143/cv/",all_cv[i],"_fold_",j,".csv")) %>% select(-X)
    pre_res <- get_roc(pre)
    sub_res[[j+1]] <- data.frame(
      fold = paste0("Fold-",j),
      ROC = pre_res[1],
      F1 = pre_res[2]
    )
  }
  sub_res <- bind_rows(sub_res)
  sub_res$para <- all_cv[i]
  cv_res[[i]] <- sub_res
}

cv_res <- bind_rows(cv_res)
tt <- cv_res %>%
  filter(!is.na(F1)) %>% group_by(para) %>%
  summarise(mean_roc = mean(ROC),
            mean_f1 = mean(F1),
            range_roc = max(ROC)-min(ROC),
            range_f1 = max(F1)-min(F1)
            ) %>%
  ungroup() %>%
  mutate(roc_f1 = mean_roc + mean_f1,
         range_score = range_roc + range_f1) %>%
  arrange(-roc_f1,range_score)

xlsx::write.xlsx(cv_res,file = "data/grid_search_res.xlsx",sheetName = "all_res",append=TRUE)
xlsx::write.xlsx(tt,file = "data/grid_search_res.xlsx",sheetName = "summ_res",append=TRUE)







