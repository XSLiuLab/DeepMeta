library(dplyr)

###########机器学习模型先把特征拼起来
cell_exp <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv",
                              data.table = F)
rownames(cell_exp) <- cell_exp$cell
cell_exp <- cell_exp %>% select(-cell)
exp_pca <- prcomp(cell_exp,scale = TRUE)
exp_pc <- exp_pca$x %>% as.data.frame() %>% select(1:50)
colnames(exp_pc) <- paste0("exp_",colnames(exp_pc))
exp_pc$cell <- rownames(exp_pc)
###gene feature pca
train_feat <- data.table::fread("data/all_train_gene_cpg.csv",data.table = F)
rownames(train_feat) <- train_feat$id
train_feat <- train_feat %>% select(-id)
feat_summ <- apply(train_feat,2,sum) %>% as.data.frame()
colnames(feat_summ) <- "counts"
feat_summ <- feat_summ %>% filter(counts > 0)
train_feat <- train_feat[,rownames(feat_summ)]

gene_pca <- prcomp(train_feat,scale = TRUE)
gene_pc <- gene_pca$x %>% as.data.frame() %>% select(1:50)
colnames(gene_pc) <- paste0("gene_",colnames(gene_pc))
gene_pc$id <- rownames(gene_pc)
gene_pc <- gene_pc %>% select(id,everything())
data.table::fwrite(gene_pc, file = "data/gene_pca2.csv",
                   row.names = FALSE, nThread=30)
train_dtV2 <- readRDS("data/train_dtV2.rds")
train_dt <- left_join(
  train_dtV2,
  exp_pc
) %>% left_join(
  .,gene_pc
) %>% select(cell,id,is_dep,everything()) %>% select(-c(id))
data.table::fwrite(train_dt,
                   file = "/home/data/sdc/wt/model_data/ml_input_pca2.csv",
                   row.names = FALSE, nThread=30)



