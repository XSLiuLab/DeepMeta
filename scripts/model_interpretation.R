library(dplyr)
library(parallel)

lift_res <- data.table::fread("/home/data/sdc/wt/model_data/deeplift_res.csv",
                              data.table = F,nThread=30)
exp <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv",data.table = F)
colnames(lift_res) <- colnames(exp)[2:ncol(exp)]

genes <- data.table::fread("/home/data/sdc/wt/model_data/deeplift_gene_res.csv",
                           data.table = F, header = T) %>% select(-V1)
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
all_gene <- data.frame(id=unique(genes$gene))
ensg2name <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  sym <- tt_gene %>% 
    paste(.,collapse=",")
  return(sym)
}
all_gene <- all_gene %>% 
  rowwise() %>% 
  mutate(gene = ensg2name(id, enz_gene_mapping)) %>% 
  ungroup()

preds <- data.table::fread("/home/data/sdc/wt/model_data/deeplift_preds_res.csv",
                           data.table = F, header = T) %>% select(-V1)
labels <- data.table::fread("/home/data/sdc/wt/model_data/deeplift_y_res.csv",
                            data.table = F, header = T) %>% select(-V1)
cells <- data.table::fread("/home/data/sdc/wt/model_data/deeplift_cell_res.csv",
                           data.table = F, header = T) %>% select(-V1)
preds_label <- bind_cols(preds,labels,genes,cells) %>% 
  mutate(preds_lab = ifelse(preds>0.5,1,0))

all_gene_res <- lapply(all_gene$id,
                       function(x){
                         tt <- which(genes == x)
                         lift_res_tt <- lift_res[tt,]
                         
                         tt_preds_lable <- preds_label[tt,]
                         tt_y_1 <- which(tt_preds_lable$label == 1 &
                                           tt_preds_lable$preds_lab == 1)
                         if (length(tt_y_1) < 10){
                           return(NA)
                         }
                         lift_res_tt_1 <- lift_res_tt[tt_y_1,]
                         top <- apply(lift_res_tt_1,1,
                                      function(x){
                                        pos_x <- x[which(x>0)]
                                        if (length(pos_x) == 0){
                                          return(NA)
                                        }
                                        pos_name <- colnames(lift_res_tt_1)[which(x>0)]
                                        top_idx <- kit::topn(pos_x,n=50)
                                        top_gene <- pos_name[top_idx]
                                        return(top_gene)
                                      },simplify = F)
                         
                         ##基因出现的频次
                         top_gene <- unlist(top)
                         gene_freq <- table(top_gene) %>% as.data.frame() %>% 
                           mutate(genes = as.character(top_gene)) %>% 
                           select(-top_gene) %>% 
                           mutate(target_gene = all_gene$gene[which(all_gene$id == x)]) %>% 
                           mutate(pos_counts = length(tt_y_1))
                         return(gene_freq)
                       })
all_gene_res <- all_gene_res[lengths(all_gene_res)>1]
all_gene_res <- bind_rows(all_gene_res)
all_gene_res <- all_gene_res %>% 
  mutate(pre = Freq / pos_counts)
saveRDS(all_gene_res,file = "data/model_exp_global.rds")

###画个图
all_gene_res <- readRDS("data/model_exp_global.rds")
apply(exp[,2:ncol(exp)],2,median) %>% as.data.frame() -> exp_median
colnames(exp_median) <- "exp"
exp_median$gene <- rownames(exp_median)

all_gene_res_exp <- left_join(
  all_gene_res,
  exp_median %>% rename(genes = gene)
)

tt <- all_gene_res_exp %>% 
  group_by(genes) %>% 
  summarise(median_pre = median(pre),
            median_exp = unique(exp),
            median_exp2 = log10(unique(exp)+1)) %>% ungroup()

library(ggpubr)
p1 <- ggscatter(tt, x = "median_pre", y = "median_exp2",
                color = "black", shape = 21, size = 3, # Points color, shape and size
                add = "reg.line",  # Add regressin line
                add.params = list(color = "blue", fill = "lightgray"), # Customize reg. line
                conf.int = TRUE, # Add confidence interval
                cor.coef = TRUE, # Add correlation coefficient. see ?stat_cor
                cor.coeff.args = list(method = "pearson", label.x = 0.3, label.sep = "\n")
)+labs(x="Median global importance score",y="Log10 median differential expression")

tt <- tt %>% mutate(type = ifelse(median_exp > 2, "High (>2)","Low (<=2)"))
p2 <- ggboxplot(tt,x="type",y="median_pre")+stat_compare_means()+
  labs(y="Median global importance score",x="Median differential expression type")
library(patchwork)
p1 + p2 + plot_layout(widths = c(2,1))
ggsave("figs/model_exp_global.pdf",width = 12,height = 6)

###local
library(igraph)

train_dtV2 <- readRDS("data/train_dtV2.rds") 
train_dtV2 <- train_dtV2 %>% 
  mutate(index = paste(cell,id,sep = "-"))

##atten weight
files <- list.files("/home/data/sdc/wt/model_data/autodl_143/autodl-tmp/atten/",full.names = T)
res <- lapply(
  files,
  function(x){data.table::fread(x,data.table = F)}
)
res <- bind_rows(res)

res <- res %>% 
  mutate(index = paste(cell,genes_x,sep = "-")) %>% 
  filter(index %in% train_dtV2$index) %>% 
  left_join(.,
            train_dtV2 %>% select(index,gene,is_dep),
            by = "index")
res$weight <- rowMeans(res[,c("w1","w2","w3")])

enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
ensg2name <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  sym <- tt_gene %>% 
    paste(.,collapse=",")
  return(sym)
}
all_genes <- data.frame(genes_y = unique(res$genes_y))

all_genes <- all_genes %>% 
  rowwise() %>% 
  mutate(gene_y = ensg2name(genes_y, enz_gene_mapping)) %>% 
  ungroup()
res <- left_join(res,all_genes)
res_filter <- res %>% filter(preds_x == is_dep) ##提取预测正确的样本
saveRDS(res_filter,file = "/home/data/sdc/wt/model_data/atten_weightV2.rds")
rm(list = ls())

get_atten_plot <- function(pre_dt,gene_name){
  test <- pre_dt %>% dplyr::filter(gene == gene_name)
  ##选择权重大于0的top3基因
  test_summ <- test %>%
    dplyr::group_by(cell) %>%
    dplyr::slice_max(order_by = weight, n = 3) %>% dplyr::ungroup() %>% 
    dplyr::filter(weight > 0)
  all_cell_counts <- length(unique(test_summ$cell))
  gene_summ <- test_summ %>% 
    dplyr::group_by(gene_y) %>% 
    dplyr::summarise(cell_counts = length(unique(cell))) %>% 
    dplyr::ungroup() %>% 
    dplyr::mutate(weight = cell_counts/all_cell_counts) %>% 
    dplyr::mutate(source = gene_name) %>% 
    dplyr::rename(target = gene_y) %>% 
    dplyr::select(source, target, weight) %>% 
    dplyr::filter(target != gene_name)
  gene_g <- graph_from_data_frame(gene_summ,directed = F)
  V(gene_g)[gene_summ$target]$color <- "#11325D"
  V(gene_g)[gene_summ$source]$color <- "#F5A673"
  par(mar=c(0,0,0,0)+.2)
  plot(gene_g,vertex.size=8, edge.arrow.size=0.3,
       vertex.label.dist=1.5, vertex.label.cex = 0.8,
       edge.width = E(gene_g)$weight*20,layout=layout.circle)
  return(gene_summ)
}

atten_weight <- readRDS("/home/data/sdc/wt/model_data/atten_weightV2.rds")
atten_summ <- atten_weight %>% select(cell,gene,is_dep) %>% distinct_all() %>% 
  group_by(gene) %>% 
  summarise(pos_c = sum(is_dep == 1),
            neg_c = sum(is_dep == 0)) %>% ungroup()
kegg <- readRDS("data/kegg_all_pathway.rds")
kegg <- kegg %>% 
  filter(grepl("Metabolism",class) | grepl("metabolism",pathway)) %>% 
  mutate(pathway = gsub(" \\- Homo sapiens \\(human\\)","",pathway)) %>% 
  filter(class != "Human Diseases; Cancer: overview")

atten_summ <- atten_summ %>% 
  filter(pos_c > 0 & neg_c > 0) %>% 
  filter(gene %in% kegg$genes)

library(org.Hs.eg.db)
keytypes(org.Hs.eg.db)
genes <- atten_summ$gene
cols <- c("GENENAME")
gene_anno <- select(org.Hs.eg.db, keys=genes, columns=cols, keytype="SYMBOL")

view_genes <- c("CAD","ACLY","COASY","DHFR","FASN","GPX4","TYMS")
for (i in 1:7){
  tt <- get_atten_plot(atten_weight,gene_name = view_genes[i])
  dev.copy2pdf(file = paste0("figs/",view_genes[i],"_atten.pdf"),
               width = 8,height = 7,out.type = "pdf")
}
detach("package:org.Hs.eg.db", unload=TRUE)

##restart
atten_weight <- readRDS("/home/data/sdc/wt/model_data/atten_weightV2.rds")
weight <- atten_weight %>% 
  select(cell,gene,gene_y,weight) %>% distinct_all() %>% 
  group_by(gene,cell) %>%
  slice_max(order_by = weight,n = 3) %>% ungroup()
cell_counts <- atten_weight %>% 
  select(cell,gene,gene_y,weight) %>% 
  group_by(gene) %>% 
  summarise(counts = length(unique(cell))) %>% ungroup()
gene_attr <- weight %>% 
  filter(weight > 0) %>% 
  group_by(gene,gene_y) %>% 
  summarise(occour_counts = n()) %>% 
  ungroup() %>% 
  left_join(.,cell_counts) %>% 
  mutate(imp = occour_counts/counts)
saveRDS(gene_attr,file = "data/local_imp.rds")

###cor
local_imp <- readRDS("data/local_imp.rds")
all_genes <- strsplit(paste(paste(local_imp$gene,collapse = ","),
                            paste(local_imp$gene_y,collapse = ","),
                            collapse = ","), split = ",")[[1]] %>% unique()

dep_dt <- readRDS("/home/data/sdb/wt/model_data/dep_dt.rds")
dep_dt <- dep_dt %>% filter(gene %in% all_genes)
dep_dt <- dep_dt %>%
  tidyr::pivot_wider(names_from = "gene",values_from = "score") %>%
  as.data.frame()
cor_res <- WGCNA::cor(as.matrix(dep_dt), use = "pairwise.complete.obs")
cor_res <- as.data.frame(cor_res)
cor_res$gene <- rownames(cor_res)
cor_res <- cor_res %>%
  tidyr::pivot_longer(cols = 1:2438,names_to = "targt_genes",values_to = "cor")
cor_res <- cor_res %>% mutate(index = paste(gene,targt_genes,sep = "-"))

gene_attr_split <- local_imp %>%
  tidyr::separate_longer_delim(cols = c("gene"),delim = ",") %>%
  tidyr::separate_longer_delim(cols = c("gene_y"),delim = ",")
cor_res_all <- inner_join(
  cor_res,
  gene_attr_split %>% mutate(index = paste(gene,gene_y,sep = "-")) %>%
    select(index,imp)
) %>% 
  mutate(cor_type = case_when(
    cor > 0.8 ~ "High Dependency Cor \n(>0.8)",
    cor < 0.2 ~ "Low Dependency Cor \n(<0.2)",
    TRUE ~ "others"
  )) %>% filter(cor_type != "others")

library(ggpubr)
ggboxplot(data=cor_res_all,x="cor_type",y="imp",xlab = FALSE,
          ylab = "Local importance score",
          order = c("Low Dependency Cor \n(<0.2)","High Dependency Cor \n(>0.8)"))+
  stat_compare_means()
ggsave("figs/local_imp_dep_cor.pdf",height = 5,width = 5)

