library(dplyr)
cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>%
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))

enz_gene_mapping <- readRDS("~/meta_target/data/enz_gene_mapping.rds")
cpg_gene <- readRDS("~/meta_target/data/cpg_gene.rds")
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
get_dep <- function(gene,dep_gene,per_cell_mode=FALSE,cel=NULL){
  tt <- strsplit(gene,",")[[1]]
  if (per_cell_mode){
    dt <- dep_gene %>% filter(gene %in% tt) %>% filter(ModelID %in% cel)
  }else{
    dt <- dep_gene %>% filter(gene %in% tt)
  }
  if (nrow(dt) == 0){
    return(NA)
  }else{
    return(max(dt$type))
  }
}

##gene exp
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

low_exp <- apply(gene_exp,1,
                 function(x){mean(x<1)}) %>% as.data.frame() %>% 
  rename(per=1)
median(low_exp$per)
##dependency data
dep_dt <- readRDS("/home/data/sdb/wt/model_data/dep_dt.rds")
dep_dt <- dep_dt %>%
  mutate(type = case_when(
    score>0.5 ~ 1,
    score<=0.5 ~ 0
  ))
dep_dt <- dep_dt %>% filter(gene %in% colnames(cpg_gene))

###筛选在网络中的基因
# ###细胞系中的基因
# enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
# ensg2name <- function(ensg,mapping){
#   tt <- strsplit(ensg," and ")[[1]]
#   tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>%
#     unique()
#   sym <- tt_gene %>%
#     paste(.,collapse=",")
#   return(sym)
# }
# cell_info <- cell_info %>% 
#   filter(ModelID %in% intersect(unique(dep_dt$ModelID),rownames(gene_exp)))
# 
# my.cluster <- parallel::makeCluster(
#   40,
#   type = "PSOCK"
# )
# #register it to be used by %dopar%
# doParallel::registerDoParallel(cl = my.cluster)
# library(foreach)
# res <- foreach(
#   i = 1:nrow(cell_info),
#   .export = c("enz_gene_mapping","cell_info","ensg2name"),
#   .packages = c("dplyr","tidyr")
# ) %dopar% {
#   cell <- cell_info$ModelID[i]
  # tissue_net <- read.table(paste0("/home/data/sdb/wt/model_data/Depmap_tINIT/enz_net/EnzGraphs/",
  #                                 paste0(gsub("-","_",cell),
  #                                        "_enzymes_based_graph.tsv")))
#   all_gene <- data.frame(id=unique(c(tissue_net$from,tissue_net$to))) %>%
#     rowwise() %>%
#     mutate(gene = ensg2name(id, enz_gene_mapping)) %>%
#     ungroup()
#   all_gene$cell <- cell
#   return(all_gene)
# }
# parallel::stopCluster(cl = my.cluster)
# 
# res <- bind_rows(res)
# res <- res %>% select(gene,cell) %>%
#   tidyr::separate_longer_delim(cols = "gene",delim = ",")
# res <- res %>% distinct_all()
# res <- res %>%
#   mutate(index = paste(cell,gene,sep = "-"))
# saveRDS(res,file = "data/all_net_genes_sample_specific.rds")

res <- readRDS("data/all_net_genes_sample_specific.rds")
dep_dt <- dep_dt %>%
  mutate(index = paste(ModelID,gene,sep = "-"))
dep_dt <- dep_dt %>% 
  filter(index %in% intersect(dep_dt$index,res$index))
##筛选阳性细胞多于3的基因
dep_summ <- dep_dt %>%
  group_by(gene) %>%
  summarise(pos_cell = sum(type == 1),
            neg_cell = sum(type == 0)) %>%
  ungroup() %>% filter(pos_cell>=3)
dep_dt <- dep_dt %>% filter(gene %in% dep_summ$gene)
cell_info <- cell_info %>% 
  filter(ModelID %in% rownames(gene_exp)) %>% 
  filter(ModelID %in% dep_dt$ModelID)
dep_dt <- dep_dt %>% filter(ModelID %in% cell_info$ModelID)

###生成训练数据
library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  80, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(cell_info),
  .export = c("cell_info","dep_dt","gene_exp","ensg2name","get_dep"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  cell <- cell_info$ModelID[i]
  tissue_net <- read.table(paste0("/home/data/sdb/wt/model_data/Depmap_tINIT/enz_net/EnzGraphs/",
                                  paste0(gsub("-","_",cell),
                                         "_enzymes_based_graph.tsv")))
  dep_cell <- dep_dt %>% filter(ModelID == cell)
  if (nrow(dep_cell) == 0){
    return(NA)
  }else{
    cell_exp <- gene_exp[cell,] %>% t() %>% as.data.frame()
    cell_exp$gene <- rownames(cell_exp)
    colnames(cell_exp)[1] <- "exp"
    cell_exp <- cell_exp$gene[which(cell_exp$exp > 1)]
    
    cell_net <- data.frame(id = unique(c(tissue_net$from,tissue_net$to))) %>% 
      rowwise() %>% 
      mutate(gene=ensg2name(id,enz_gene_mapping,cell_exp,cpg_gene)[[1]],
             fea=ensg2name(id,enz_gene_mapping,cell_exp,cpg_gene)[[2]],
             is_exp=ensg2name(id,enz_gene_mapping,cell_exp,cpg_gene)[[3]]) %>% 
      ungroup() %>% 
      filter(nchar(gene)>1 & is_exp) 
    cell_net <- cell_net %>% 
      tidyr::separate_wider_delim(cols = fea,delim = ",",names_sep="-") %>% 
      mutate(across(starts_with("fea-"),as.numeric))
    ####add dep data
    
    cell_net_dep <- cell_net %>% 
      rowwise() %>% 
      mutate(is_dep = get_dep(gene,dep_cell)) %>% 
      ungroup() %>% 
      select(1:2,is_dep,everything()) %>% 
      na.omit()
    
    write.table(tissue_net,
                file = paste0("/home/data/sdb/wt/model_data/enzyme_train_sample/",
                              cell,".txt"),sep = "\t",
                row.names = F)
    write.table(cell_net_dep,
                file = paste0("/home/data/sdb/wt/model_data/enzyme_train_sample/",
                              cell,"_feat.txt"),
                sep = "\t",row.names = F)
  }
}
parallel::stopCluster(cl = my.cluster)

write.csv(cell_info %>% select(ModelID) %>% 
            rename(cell=ModelID),
          file = "data/sample_specific_cell_info.csv",
          quote = F,row.names = F)
#####使用和 DeepMeta 相同的训练集细胞系
###结果
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

all_cv <- c("own_model","sample_specific")
cv_res <- vector("list",2)
for (i in 1:2){
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

ggbarplot(cv_res, x = "type", y = "Value",fill="Metric",
          add = "mean_se", label = TRUE, 
          lab.vjust = -0.5,position = position_dodge(0.9),
          order = c("own_model","sample_specific"),
          palette = c("#2B4F7D","#7BC5E3","#EEDC96","#637951"),
          lab.nb.digits = 2, xlab=F)+
  scale_x_discrete(labels=c("DeepMeta","Model using sample-specific enzyme network"))+
  labs(x="Type",y='Value')
ggsave("report/compare_sample_specific_model.pdf",width = 8,height = 5)
