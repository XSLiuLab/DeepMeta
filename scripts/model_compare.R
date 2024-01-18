library(dplyr)
library(yardstick)
library(ggprism)
library(ggplot2)
library(ggpubr)

get_roc <- function(dt,need="all"){
  dt <- dt %>% 
    mutate(truth = ifelse(label == 1, "Class1","Class2"),
           pred_label = ifelse(preds == 1, "Class1","Class2"))
  dt$truth <- factor(dt$truth)
  dt$pred_label <- factor(dt$pred_label)
  
  f1 <- try(f_meas(dt,truth,pred_label)[".estimate"] %>%
              unlist() %>% unname() %>% round(.,2),silent = TRUE)
  kappa <- try(kap(dt,truth,pred_label)[".estimate"] %>%
                 unlist() %>% unname() %>% round(.,2),silent = TRUE)
  f1 <- ifelse('try-error' %in% class(f1),NA,f1)
  kappa <- ifelse('try-error' %in% class(kappa),NA,kappa)
  if (need == "all"){
    roc <-  roc_auc(dt, truth, preds_raw)[".estimate"] %>% 
      unlist() %>% unname() %>% round(.,2)
    pr <-  pr_auc(dt, truth, preds_raw)[".estimate"] %>%
      unlist() %>% unname() %>% round(.,2)
    return(c(roc,pr,f1,kappa))
  }else{
    res <- switch(
      need,
      "roc" = roc_auc(dt, truth, preds_raw)[".estimate"] %>% 
        unlist() %>% unname() %>% round(.,2),
      "pr" = pr_auc(dt, truth, preds_raw)[".estimate"] %>%
        unlist() %>% unname() %>% round(.,2),
      "f1" = f1,
      "kappa" = kappa
    )
    return(res)
  }
  
}

get_gene_metric <- function(dt,gene_name,need){
  dt <- dt %>% 
    mutate(truth = ifelse(label == 1, "Class1","Class2"),
           pred_label = ifelse(preds == 1, "Class1","Class2")) %>% 
    filter(genes == gene_name)
  dt$truth <- factor(dt$truth)
  dt$pred_label <- factor(dt$pred_label)
  
  f1 <- try(f_meas(dt,truth,pred_label)[".estimate"] %>%
              unlist() %>% unname() %>% round(.,2),silent = TRUE)
  kappa <- try(kap(dt,truth,pred_label)[".estimate"] %>%
                 unlist() %>% unname() %>% round(.,2),silent = TRUE)
  f1 <- ifelse('try-error' %in% class(f1),NA,f1)
  kappa <- ifelse('try-error' %in% class(kappa),NA,kappa)
  
  res <- switch(
    need,
    "roc" = roc_auc(dt, truth, preds_raw)[".estimate"] %>% 
      unlist() %>% unname() %>% round(.,2),
    "pr" = pr_auc(dt, truth, preds_raw)[".estimate"] %>%
      unlist() %>% unname() %>% round(.,2),
    "f1" = f1,
    "kappa" = kappa
  )
  return(res)
}

test_dtV2 <- readRDS("data/test_dtV2.rds")
train_dtV2 <- readRDS("data/train_dtV2.rds")
all_dt <- bind_rows(test_dtV2,train_dtV2)
all_dt_summ <- all_dt %>% 
  group_by(gene) %>% 
  summarise(pos_counts = mean(is_dep == 1)) %>% 
  ungroup() %>% 
  filter(pos_counts < 0.9 & pos_counts > 0.1)
all_dt_summ <- left_join(
  all_dt_summ,
  all_dt %>% select(id,gene) %>% distinct_all()
)

all_cv <- c("own_model","random_feat","random_exp",
            "no_gnn","gnn_model_random_y",
            "ml_logistic","ml_rf","ml_svm")
cv_res <- vector("list",8)
for (i in 1:8){
  sub_res <- vector("list",10)
  for (j in 1:10){
    pre <- read.csv(paste0("data/cv/",all_cv[i],"/fold_",j-1,".csv")) %>% select(-X)
    pre_roc <- ifelse(i %in% c(1:7),get_roc(pre,need = "roc"),NA)
    pre_pr <- ifelse(i %in% c(1:7),get_roc(pre,need = "pr"),NA)
    pre_f1 <- get_roc(pre,need = "f1")
    pre_kap <- get_roc(pre,need = "kappa")
    dt <- pre %>% 
      filter(genes %in% all_dt_summ$id) %>% 
      group_by(genes) %>%
      summarise(pos_counts = sum(label == 1),neg_counts = sum(label == 0)) %>%
      ungroup() %>% filter(pos_counts>0 & neg_counts > 0)
    pre_gene_roc <- data.frame(gene = unique(dt$genes)) %>%
      rowwise() %>% 
      mutate(roc = ifelse(i == 8,NA,get_gene_metric(pre,gene_name = gene,"roc")),
             f1 = get_gene_metric(pre,gene_name = gene,"f1")) %>% 
      ungroup()
    
    sub_res[[j]] <- data.frame(
      fold = paste0("Fold-",j),
      ROC = pre_roc,
      PR = pre_pr,
      F1 = pre_f1,
      Kappa = pre_kap,
      per_roc = mean(pre_gene_roc$roc,na.rm=T),
      per_f1 = mean(pre_gene_roc$f1,na.rm=T)
    ) 
  }
  sub_res <- bind_rows(sub_res)
  sub_res$type <- all_cv[i]
  cv_res[[i]] <- sub_res
}
cv_res <- bind_rows(cv_res)
cv_res <- cv_res %>% 
  tidyr::pivot_longer(cols = c("ROC","PR","F1","Kappa","per_roc","per_f1"),
                      names_to = "Metric",values_to = "Value")
cv_res <- na.omit(cv_res)
saveRDS(cv_res,file = "data/cv_res.rds")

###plot
cv_res <- cv_res %>% 
  filter(Metric != "per_f1") %>% 
  mutate(Metric = ifelse(Metric == "per_roc","Per-ROC",Metric))
cv_1 <- cv_res %>% 
  filter(!(type %in% c("ml_logistic","ml_rf","ml_svm")))
cv_2 <- cv_res %>% 
  filter(type %in% c("own_model","ml_logistic","ml_rf","ml_svm"))

ggbarplot(cv_1, x = "type", y = "Value",fill="Metric",
          add = "mean_se", label = TRUE, 
          lab.vjust = -0.5,position = position_dodge(0.9),
          order = c("own_model","random_feat","random_exp","no_gnn",
                    "gnn_model_random_y"),
          palette = c("#2B4F7D","#3B77B0","#7BC5E3","#EEDC96","#637951"),
          lab.nb.digits = 2, xlab=F)+
  scale_x_discrete(labels=c("DeepMeta","Random Feature",
                            "Random Expression",
                            "Without GNN","Scrambled label"))+
  labs(x="Type",y='Value')
ggsave("figs/model_compare_with_net_1.pdf",width = 14,height = 5)

ggbarplot(cv_2, x = "type", y = "Value",fill="Metric",
          add = "mean_se", label = TRUE, 
          lab.vjust = -0.5,position = position_dodge(0.9),
          order = c("own_model","ml_logistic","ml_rf","ml_svm"),
          palette = c("#2B4F7D","#3B77B0","#7BC5E3","#EEDC96","#637951"),
          lab.nb.digits = 2, xlab=F)+
  scale_x_discrete(labels=c("DeepMeta","Logistic","RF","SVM"))+
  labs(x="Type",y='Value')
ggsave("figs/model_compare_with_net_2.pdf",width = 13,height = 5)

###十折中平均的 per-gene ROC
res <- vector("list",10)
for (j in 1:10){
  pre <- read.csv(paste0("data/cv/own_model","/fold_",j-1,".csv")) %>% select(-X)
  dt <- pre %>% 
    filter(genes %in% all_dt_summ$id) %>% 
    group_by(genes) %>%
    summarise(pos_counts = sum(label == 1),neg_counts = sum(label == 0)) %>%
    ungroup() %>% filter(pos_counts>0 & neg_counts > 0)
  pre_gene_roc <- data.frame(gene = unique(dt$genes)) %>%
    rowwise() %>% 
    mutate(roc = get_gene_metric(pre,gene_name = gene,"roc")) %>% 
    ungroup()
  pre_gene_roc$fold <- j
  pre_gene_roc <- left_join(
    pre_gene_roc %>% rename(id = gene),
    all_dt_summ
  )
  res[[j]] <- pre_gene_roc
}
res <- bind_rows(res)

res_gene_summ <- res %>% 
  group_by(gene) %>% 
  summarise(ave_per_roc = mean(roc)) %>% 
  ungroup()
saveRDS(res_gene_summ,file = "data/per_roc.rds")
res_gene_summ$tt <- ""
ggboxplot(res_gene_summ,x="tt",y="ave_per_roc",xlab = FALSE,
          ylab = "Average ROC")
ggsave("figs/per_roc.pdf",width = 3,height = 5)

###不同组学的模型比较
###突变
train_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/cell_net_filter_exp/raw/train_cell_info.csv")
test_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/enzyme_net/test/raw/test_cell_info.csv")

mut <- data.table::fread("/home/data/sdb/wt/model_data/OmicsSomaticMutations.csv",
                         data.table = F)
mut <- mut %>% 
  select(DepMap_ID,HugoSymbol,VariantInfo,Chrom,Pos,Ref,Alt) %>% 
  filter(VariantInfo != "SILENT") %>% 
  filter(DepMap_ID %in% c(train_dt$cell,test_dt$cell))
mut_summ <- mut %>% 
  group_by(HugoSymbol) %>% 
  summarise(cell_c=length(unique(DepMap_ID))) %>% 
  ungroup() %>% 
  filter(cell_c > 20)

mut <- mut %>% 
  filter(HugoSymbol %in% mut_summ$HugoSymbol) %>% 
  select(DepMap_ID,HugoSymbol) %>% 
  distinct_all() %>% 
  rename(cell=DepMap_ID,gene=HugoSymbol) %>% 
  mutate(value = 1) %>% 
  tidyr::pivot_wider(names_from = gene, values_from = value, values_fill = 0)

train_dt <- train_dt %>% filter(cell %in% mut$cell)
train_dt$cell_index <- (1:nrow(train_dt))-1

test_dt <- test_dt %>% filter(cell %in% mut$cell)
test_dt$cell_index <- (1:nrow(test_dt))-1

write.csv(train_dt,
          file = "data/mut_train_dt.csv",
          quote = F,row.names = F)
write.csv(test_dt,
          file = "data/mut_test_dt.csv",
          quote = F,row.names = F)
write.csv(mut,file = "/home/data/sda/wt/model_data/cell_mut.csv",
          quote = F,row.names = F)

####cnv
train_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/cell_net_filter_exp/raw/train_cell_info.csv")
test_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/enzyme_net/test/raw/test_cell_info.csv")
cnv <- data.table::fread("/home/data/sdb/wt/model_data/OmicsCNGene.csv",
                         data.table = F)
colnames(cnv)[2:ncol(cnv)] <- gsub("\\s*\\([^\\)]+\\)","",colnames(cnv)[2:ncol(cnv)])
colnames(cnv)[1] <- "cell"

gene_cnv_sd <- apply(cnv[,2:ncol(cnv)],2,sd)
gene_cnv_sd <- gene_cnv_sd[which(gene_cnv_sd > 0.2)]
cnv <- cnv %>% select(cell, names(gene_cnv_sd))

train_dt <- train_dt %>% filter(cell %in% cnv$cell)
train_dt$cell_index <- (1:nrow(train_dt))-1

test_dt <- test_dt %>% filter(cell %in% cnv$cell)
test_dt$cell_index <- (1:nrow(test_dt))-1

write.csv(train_dt,
          file = "data/cnv_train_dt.csv",
          quote = F,row.names = F)
write.csv(test_dt,
          file = "data/cnv_test_dt.csv",
          quote = F,row.names = F)
write.csv(cnv,file = "/home/data/sdc/wt/model_data/cell_cnv.csv",
          quote = F,row.names = F)

###cv
all_cv <- c("own_model","mut_model","cnv_model")
cv_res <- vector("list",3)
for (i in 1:3){
  sub_res <- vector("list",10)
  for (j in 1:10){
    pre <- read.csv(paste0("data/cv/",all_cv[i],"/fold_",j-1,".csv")) %>% select(-X)
    pre_roc <- get_roc(pre,need = "roc")
    pre_pr <- get_roc(pre,need = "pr")
    pre_f1 <- get_roc(pre,need = "f1")
    pre_kap <- get_roc(pre,need = "kappa")
    sub_res[[j]] <- data.frame(
      fold = paste0("Fold-",j),
      ROC = pre_roc,
      PR = pre_pr,
      F1 = pre_f1,
      Kappa = pre_kap
    ) 
  }
  sub_res <- bind_rows(sub_res)
  sub_res$type <- all_cv[i]
  cv_res[[i]] <- sub_res
}
cv_res <- bind_rows(cv_res)
cv_res <- cv_res %>% 
  tidyr::pivot_longer(cols = c("ROC","PR","F1","Kappa"),
                      names_to = "Metric",values_to = "Value")
cv_res <- na.omit(cv_res)
###plot
library(ggpubr)
ggbarplot(cv_res, x = "type", y = "Value",fill="Metric",
          add = "mean_se", label = TRUE, 
          lab.vjust = -0.5,position = position_dodge(0.9),
          order = c("own_model","mut_model","cnv_model"),
          palette = c("#F5A673","#11325D","#736B9D","#B783AF"),
          lab.nb.digits = 3, xlab=F)+
  scale_x_discrete(labels=c("DeepMeta","Mutation Model","CNV Model"))+
  labs(x="Type",y='Value')
ggsave("figs/model_compare_omics.pdf",width = 16,height = 6)

####RF test
pre <- read.csv("data/rf_pred_test.csv") %>% select(-X)
pre <- pre %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds == 1, "Class1","Class2"))
pre$truth <- factor(pre$truth)
pre$pred_label <- factor(pre$pred_label)

pr <-  pr_auc(pre, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
roc <-  roc_auc(pre, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)

###RF sanger
pre <- read.csv("data/rf_pred_sanger.csv") %>% select(-X)
train_dt <- data.table::fread("/home/wt/meta_target/data/train_dtV2.csv",
                              data.table = F)
seen <- pre %>% filter(cell %in% train_dt$cell)
unseen <- pre %>% filter(!(cell %in% train_dt$cell))
seen <- seen %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds == 1, "Class1","Class2"))
unseen <- unseen %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds == 1, "Class1","Class2"))
seen$truth <- factor(seen$truth)
seen$pred_label <- factor(seen$pred_label)
unseen$truth <- factor(unseen$truth)
unseen$pred_label <- factor(unseen$pred_label)

seen_pr <-  pr_auc(seen, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
seen_roc <-  roc_auc(seen, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
unseen_pr <-  pr_auc(unseen, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
unseen_roc <-  roc_auc(unseen, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)

dt <- unseen %>% 
  group_by(genes) %>%
  summarise(pos_counts = sum(label == 1),neg_counts = sum(label == 0)) %>%
  ungroup() %>% filter(pos_counts>0 & neg_counts > 0)
pre_gene_roc <- data.frame(gene = unique(dt$genes)) %>%
  rowwise() %>% 
  mutate(roc = get_gene_metric(unseen,gene_name = gene,"roc"),
         f1 = get_gene_metric(unseen,gene_name = gene,"f1")) %>% 
  ungroup()

###RF CCMA
ccma_dep <- data.table::fread("/home/data/sdc/wt/model_data/CCMA_CRISPRdependencyscore.csv",
                              data.table = F)
pre <- read.csv("data/rf_pred_ccma.csv") %>% select(-X)
pre_gene <-  data.frame(id = unique(c(pre$genes))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]])
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(genes=id)) %>% 
  tidyr::separate_longer_delim(cols = gene,delim = ",") %>% 
  left_join(.,
            ccma_dep %>% select(sample,gene,z) %>% 
              rename(cell=sample)) %>% 
  filter(!is.na(z)) %>%  
  group_by(cell,genes) %>% 
  summarise(preds = unique(preds),
            preds_raw = unique(preds_raw),
            z_score = min(z,na.rm = T)) %>% 
  ungroup() %>% 
  mutate(label = ifelse(z_score < (-1.5),1,0))

train_dt <- data.table::fread("/home/wt/meta_target/data/train_dtV2.csv",
                              data.table = F)
pre <- pre %>% filter(genes %in% train_dt$id)

ggboxplot(data=pre,x="preds",y="z_score")+
  stat_compare_means()+
  labs(x="Model Prediction",y="Dependency Z scores")+
  geom_hline(yintercept = (-1.5),color="red",linewidth=2,linetype=2)

###RF Drug
pre <- read.csv("data/rf_pred_drug.csv") %>% select(-X)
pre_gene <-  data.frame(id = unique(c(pre$genes))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]])
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(genes=id))
nc_drug <- readRDS("data/nc_drug.rds")
get_auc <- function(cell,gene,drug_dt){
  genes <- strsplit(gene,split = ",")[[1]]
  res <- vector("list",length(genes))
  for (i in 1:length(genes)){
    drug_dt_t <- drug_dt %>% filter(grepl(genes[i],target_gene))
    res[[i]] <- drug_dt_t
  }
  res <- bind_rows(res)
  if (nrow(res)==0){
    return(list(NA,NA))
  }else{
    res <- res %>% select(cell,drug_name)
    return(list(max(res[,1]),
                paste(res$drug_name[unlist(res[,1]==max(res[,1]))],collapse = ",")))
  }
}

pre_drug <- pre %>% 
  rowwise() %>% 
  mutate(max_auc = get_auc(cell,gene,nc_drug)[[1]],
         max_auc_drug = get_auc(cell,gene,nc_drug)[[2]]) %>% 
  ungroup()
pre_drug <- pre_drug %>% 
  na.omit() %>% 
  group_by(cell) %>% 
  mutate(rank_auc = rank(-max_auc)) %>% 
  ungroup()

pre_drug_top <- pre_drug %>%
  group_by(cell) %>% 
  slice_max(order_by = preds_raw, n =5) %>% 
  ungroup()
###label rank = 1
library(ggpubr)
library(ggrepel)
pre_drug_top <- pre_drug_top %>% 
  mutate(ll = ifelse(rank_auc <= 3,
                     paste0("Rank: ",rank_auc,", AUC: ",round(max_auc,2)),
                     "")) %>% 
  mutate(ll1 = ifelse(rank_auc <= 3,
                      paste0("Rank: ",rank_auc),""))
p1 <- ggscatter(pre_drug_top, x = "cell", y = "preds_raw")+
  geom_text_repel(aes(label = ll1),box.padding = 0.5,min.segment.length = 0,
                  max.overlaps = Inf) +
  geom_point(color = ifelse(pre_drug_top$ll == "", "grey50", "red"),size=4)+
  rotate_x_text(45)+
  labs(x="Cell Lines",y="Model Prediction")

own_drug <- readRDS("data/nc_pre_drug_V2.rds")
own_drug_top <- own_drug %>%
  group_by(cell) %>% 
  slice_max(order_by = preds_raw, n =5) %>% 
  ungroup()
own_drug_top <- own_drug_top %>% 
  mutate(ll = ifelse(rank_auc <= 3,
                     paste0("Rank: ",rank_auc,", AUC: ",round(max_auc,2)),
                     ""))
dt <- data.frame(
  aucs = c(pre_drug_top$max_auc[which(nchar(pre_drug_top$ll)>0)],
           own_drug_top$max_auc[which(nchar(own_drug_top$ll)>0)]),
  types = c(rep("RF",
                length(pre_drug_top$max_auc[which(nchar(pre_drug_top$ll)>0)])),
            rep("DeepMeta",
                length(own_drug_top$max_auc[which(nchar(own_drug_top$ll)>0)])))
)
p2 <- ggboxplot(dt,x="types",y="aucs",xlab = FALSE,ylab = "AUC")+
  stat_compare_means()

library(patchwork)
p1
ggsave("figs/RF_compare_drug_rank.pdf",width = 8,height = 4)
p2 
ggsave("figs/RF_compare_drug_auc.pdf",width = 4,height = 4)
