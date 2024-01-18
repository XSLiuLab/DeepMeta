###pre
library(yardstick)
library(dplyr)
library(ggplot2)
library(ggprism)
library(patchwork)
library(ggpubr)

enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
cpg_gene <- readRDS("data/cpg_gene.rds")
ensg2name <- function(ensg,mapping,exp_list=NULL,pathway){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  sym <- tt_gene %>% 
    paste(.,collapse=",")
  gene_fea <- pathway %>% select(any_of(tt_gene))
  if (nrow(gene_fea) != 0){
    if (ncol(gene_fea) == 1){
      gene_fea$final <- gene_fea[,1]
    }else{
      gene_fea$final <- rowSums(gene_fea[ , colnames(gene_fea)]) 
    }
    gene_fea <- gene_fea %>% 
      mutate(final = ifelse(final == 0,0,1))
    if (is.null(exp_list)){
      is_exp <- NA
    }else{
      is_exp <- all(tt_gene %in% exp_list)
    }
    return(list(sym,paste(gene_fea$final,collapse = ","),is_exp))
  }else{
    if (is.null(exp_list)){
      is_exp <- NA
    }else{
      is_exp <- all(tt_gene %in% exp_list)
    }
    return(list(sym,NA,is_exp))
  }
}

pre <- read.csv("data/test_preV2.csv") %>% select(-X)
pre <- pre %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds == 1, "Class1","Class2"))
pre$truth <- factor(pre$truth)
pre$pred_label <- factor(pre$pred_label)

pr <-  pr_auc(pre, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
roc <-  roc_auc(pre, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)

p1 <- roc_curve(pre, truth, preds_raw) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1) +
  theme_bw() +
  annotate(geom="text", x=0.75, y=0.7, label=paste0("ROC-AUC: ",roc),
           size=5)
ggsave("figs/test_roc.pdf",p1,width = 4,height = 4)

###sanger
sanger <- data.table::fread("/home/data/sdc/wt/model_data/sanger_gene_dependency_chronos.csv",data.table = F)
colnames(sanger)[2:ncol(sanger)] <- gsub("\\s*\\([^\\)]+\\)","",
                                         colnames(sanger)[2:ncol(sanger)])
sanger <- sanger %>% 
  tidyr::pivot_longer(cols = colnames(sanger)[2:ncol(sanger)],names_to = "gene",
                      values_to = "score")
colnames(sanger)[1] <- "cell"

sanger <- sanger %>%
  filter(!is.na(score)) %>% 
  mutate(type = case_when(
    score>0.5 ~ 1,
    score<=0.5 ~ 0
  ))

pre <- read.csv("data/sanger_preV2.csv") %>% select(-X)
pre_gene <-  data.frame(id = unique(c(pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]])
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(gene_name=id)) %>% 
  tidyr::separate_longer_delim(cols = gene,delim = ",") %>% 
  left_join(.,sanger) %>% filter(!is.na(type)) %>%  
  group_by(cell,gene_name) %>% 
  summarise(preds = unique(preds),
            preds_raw = unique(preds_raw),
            label = max(type,na.rm = T)) %>% 
  ungroup()
saveRDS(pre,file = "data/sanger_pre.rds")

pre <- readRDS("data/sanger_pre.rds")
pre <- pre %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds == 1, "Class1","Class2"))
pre$truth <- factor(pre$truth)
pre$pred_label <- factor(pre$pred_label)

train_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/cell_net_filter_exp/raw/train_cell_info.csv")
common_cell <- pre %>% filter(cell %in% train_dt$cell)
unique_cell <- pre %>% filter(!(cell %in% train_dt$cell))

roc_comm <-  roc_auc(common_cell, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
roc_unique <-  roc_auc(unique_cell, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)


p4 <- roc_curve(common_cell, truth, preds_raw) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.5, y=0.7, label=paste0("Seen Cells ROC-AUC: ",roc_comm),
           size=4)

p5 <- roc_curve(unique_cell, truth, preds_raw) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.5, y=0.7, label=paste0("Unseen Cell ROC-AUC: ",roc_unique),
           size=4)
p4+p5
ggsave("figs/sanger_roc.pdf",width = 8,height = 4)

p4
ggsave("figs/sanger_roc_seen.pdf",width = 4,height = 4)
p5
ggsave("figs/sanger_roc_unseen.pdf",width = 4,height = 4)

pdf("figs/sanger_roc_seen.pdf",height = 4,width = 4)
print(p4)
dev.off()
#####rnai
pre <- read.csv("data/rnai_preV2.csv") %>% select(-X)
train_dt <- readRDS("data/train_dtV2.rds")
pre <- pre %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"))
pre$truth <- factor(pre$truth)

seen <- pre %>% filter(cell %in% train_dt$cell)
unseen <- pre %>% filter(!(cell %in% train_dt$cell))

seen_roc <-  roc_auc(seen, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
unseen_roc <-  roc_auc(unseen, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)

p2 <- roc_curve(seen, truth, preds_raw) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.5, y=0.7, label=paste0("Seen Cells ROC-AUC: ",seen_roc),
           size=4)
p3 <- roc_curve(unseen, truth, preds_raw) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.5, y=0.7, label=paste0("Unseen Cells ROC-AUC: ",
                                                   unseen_roc),
           size=4)
p2+p3
ggsave("figs/rnai_roc.pdf",width = 8,height = 4)

###ccma
ccma_dep <- data.table::fread("/home/data/sdc/wt/model_data/CCMA_CRISPRdependencyscore.csv",
                              data.table = F)
pre <- read.csv("data/ccma_preV2.csv") %>% select(-X)
pre_gene <-  data.frame(id = unique(c(pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]])
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(gene_name=id)) %>% 
  tidyr::separate_longer_delim(cols = gene,delim = ",") %>% 
  left_join(.,
            ccma_dep %>% select(sample,gene,z) %>% 
              rename(cell=sample)) %>% 
  filter(!is.na(z)) %>%  
  group_by(cell,gene_name) %>% 
  summarise(preds = unique(preds),
            preds_raw = unique(preds_raw),
            z_score = min(z,na.rm = T)) %>% 
  ungroup() %>% 
  mutate(label = ifelse(z_score < (-1.5),1,0))
saveRDS(pre,file = "data/ccma_pre.rds")

pre <- pre %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds == 1, "Class1","Class2"))
pre$truth <- factor(pre$truth)
pre$pred_label <- factor(pre$pred_label)

roc <-  roc_auc(pre, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)

p6 <- roc_curve(pre, truth, preds_raw) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.4, y=0.85, label=paste0("ROC-AUC: ",roc),
           size=5)
p7 <- ggboxplot(data=pre,x="preds",y="z_score")+
  stat_compare_means()+
  labs(x="Model Prediction",y="Dependency Z scores")+
  geom_hline(yintercept = (-1.5),color="red",linewidth=2,linetype=2)

p6+p7
ggsave("figs/ccma_roc.pdf",width = 8,height = 4)


###附图
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

train_dt <- readRDS("data/train_dtV2.rds")

test <- read.csv("data/test_preV2.csv") %>% select(-X)
test_m <- get_roc(test)
sanger <- readRDS("data/sanger_pre.rds")
sanger_common_cell <- sanger %>% filter(cell %in% train_dt$cell)
sanger_unique_cell <- sanger %>% filter(!(cell %in% train_dt$cell))
sanger_seen_m <- get_roc(sanger_common_cell)
sanger_unseen_m <- get_roc(sanger_unique_cell)
rnai <- read.csv("data/rnai_preV2.csv") %>% select(-X)
rnai_seen <- rnai %>% filter(cell %in% train_dt$cell)
rnai_unseen <- rnai %>% filter(!(cell %in% train_dt$cell))
rnai_seen_m <- get_roc(rnai_seen)
rnai_unseen_m <- get_roc(rnai_unseen)
ccma <- readRDS("data/ccma_pre.rds")
ccma_m <- get_roc(ccma)

dt <- data.frame(
  CCMA = ccma_m,
  `RNAi-Seen` = rnai_seen_m,
  `Sanger-Seen` = sanger_seen_m,
  `RNAi-Unseen` = rnai_unseen_m,
  `Sanger-Unseen` = sanger_unseen_m,
  Test = test_m,
  Metrics = c("AUROC","AUPRC","F1","Kappa") 
  ,check.names = FALSE) %>% 
  tidyr::pivot_longer(cols = 1:6,names_to = "Type",values_to = "Value")

ggbarplot(dt, x = "Type", y = "Value",fill="Metrics",
          label = TRUE, 
          lab.vjust = -0.5,position = position_dodge(0.9),
          order = c("Test","RNAi-Seen","RNAi-Unseen","Sanger-Seen",
                    "Sanger-Unseen","CCMA"),
          palette = c("#2B4F7D","#3B77B0","#EEDC96","#637951"),
          lab.nb.digits = 2, xlab=F)
ggsave("figs/other_metrics.pdf",width = 12,height = 4.5)

####drug
pre <- read.csv("data/drug_preV2.csv") %>% select(-X)
pre_gene <-  data.frame(id = unique(c(pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]])
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(gene_name=id))
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
saveRDS(pre_drug,file = "data/nc_pre_drug_V2.rds")

pre_drug <- readRDS("data/nc_pre_drug_V2.rds")
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
ggscatter(pre_drug_top, x = "cell", y = "preds_raw")+
  geom_text_repel(aes(label = ll),box.padding = 0.5,min.segment.length = 0,
                  max.overlaps = Inf) +
  geom_point(color = ifelse(pre_drug_top$ll == "", "grey50", "red"),size=4)+
  rotate_x_text(45)+
  labs(x="Cell Lines",y="Model Prediction")
ggsave("figs/nc_drug.pdf",width = 11,height =6)

###正常样本的预测表现
pre <- read.csv("data/cell_normal_preV2.csv") %>% select(-X)
pre <- pre %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds == 1, "Class1","Class2"))
pre$truth <- factor(pre$truth)
pre$pred_label <- factor(pre$pred_label)

pr <-  pr_auc(pre, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
roc <-  roc_auc(pre, truth, preds_raw)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)

p1 <- pr_curve(pre, truth, preds_raw) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.5, y=0.85, label=paste0("PRC-AUC: ",pr),
           size=5)

p2 <- roc_curve(pre, truth, preds_raw) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.4, y=0.85, label=paste0("ROC-AUC: ",roc),
           size=5)
library(patchwork)
p1 + p2
ggsave("figs/normal_cell_pre.pdf",width = 11,height = 5)

# library(cvms)
# basic_table <- table(pre %>% select(label, preds))
# cfm <- as_tibble(basic_table)
# plot_confusion_matrix(cfm, 
#                       target_col = "label", 
#                       prediction_col = "preds",
#                       counts_col = "n")



