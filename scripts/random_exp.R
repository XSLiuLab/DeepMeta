library(dplyr)
cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>% 
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))

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
get_dep <- function(gene,dep_gene,per_cell_mode=FALSE,cel=NULL){
  tt <- strsplit(gene,",")[[1]]
  if (per_cell_mode){
    dt <- dep_gene %>% filter(gene %in% tt) %>% filter(cell %in% cel)
  }else{
    dt <- dep_gene %>% filter(gene %in% tt)
  }
  if (nrow(dt) == 0){
    return(NA)
  }else{
    return(max(dt$type))
  }
}

net_cell_mapping <- data.frame(net=NA,cell=unique(cell_info$OncotreeLineage)) 
net_cell_mapping$net <- c("iOvarianCancer1620.xml","bone_marrow.xml",
                          "iColorectalCancer1750.xml","iSkinCancer1386.xml",
                          "iUrothelialCancer1532.xml","iLungCancer1490.xml",
                          "kidney.xml","iBreastCancer1746.xml",NA,
                          "iPancreaticCancer1613.xml","brain.xml",NA,NA,
                          "iStomachCancer1511.xml","iThyroidCancer1710.xml",NA,NA,
                          "iProstateCancer1560.xml",NA,"iHeadNeckCancer1628.xml",
                          "iEndometrialCancer1713.xml",NA,"iLiverCancer1788.xml",
                          "iCervicalCancer1611.xml",NA,NA,"adrenal_gland.xml",
                          NA,"iTestisCancer1483.xml",NA)
net_cell_mapping <- na.omit(net_cell_mapping)

gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

dep_dt <- readRDS("/home/data/sdb/wt/model_data/dep_dt.rds")
dep_dt <- dep_dt %>%
  mutate(type = case_when(
    score>0.5 ~ 1,
    score<=0.5 ~ 0
  ))
dep_dt <- dep_dt %>% filter(gene %in% colnames(cpg_gene))

cell_info <- cell_info %>% 
  select(ModelID,OncotreeLineage) %>% 
  rename(cell = OncotreeLineage) %>% 
  left_join(.,net_cell_mapping)
cell_info <- na.omit(cell_info)

res <- readRDS("data/all_net_genes.rds")
dep_dt <- dep_dt %>%
  mutate(index = paste(ModelID,gene,sep = "-")) %>%
  filter(index %in% res$index)
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

###随机化基因表达
gene_exp_random <- t(apply(gene_exp, 1, sample)) %>% as.data.frame()
colnames(gene_exp_random) <- colnames(gene_exp)
saveRDS(gene_exp_random,file = "/home/data/sdc/wt/model_data/gene_exp_random.rds")
###并行
library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(cell_info),
  .export = c("net_cell_mapping","cell_info","dep_dt","gene_exp_random",
              "ensg2name","get_dep"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  cell <- cell_info$ModelID[i]
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(gsub(".xml","",cell_info$net[i]),
                                         "_enzymes_based_graph.tsv")))
  dep_cell <- dep_dt %>% filter(ModelID == cell)
  if (nrow(dep_cell) == 0){
    return(NA)
  }else{
    cell_exp <- gene_exp_random[cell,] %>% t() %>% as.data.frame()
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
      select(1:2,is_dep,everything())
    
    write.table(tissue_net,
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_random/",cell,".txt"),sep = "\t",
                row.names = F)
    write.table(cell_net_dep,
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_random/",cell,"_feat.txt"),
                sep = "\t",row.names = F)
  }
}
parallel::stopCluster(cl = my.cluster)

all_dt <- cell_info %>% select(ModelID) %>% rename(cell=ModelID)
set.seed(2023042602)
test_cell <- sample(all_dt$cell,size = nrow(all_dt)*0.1)
test_dt <- all_dt %>% filter(cell %in% test_cell)
train_dt <- all_dt %>% filter(!(cell %in% test_cell))
train_dt$cell_index <- (1:nrow(train_dt))-1
write.csv(train_dt,
          file = "data/train_cell_info_random_exp.csv",
          quote = F,row.names = F)

###
origin_diff <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv",data.table = F)
normal_exp <- readRDS("/home/data/sdb/wt/model_data/gtex_normal_exp.rds")
normal_exp <- normal_exp[colnames(origin_diff)[2:ncol(origin_diff)],origin_diff$cell]
normal_exp_random <- apply(normal_exp, 2, sample) %>% as.data.frame()
rownames(normal_exp_random) <- rownames(normal_exp)

tumor_exp <- gene_exp_random %>% t() %>% as.data.frame()
tumor_exp <- tumor_exp[rownames(normal_exp),colnames(normal_exp)] %>% as.matrix()

diff <- tumor_exp / (log2(normal_exp_random + 1.01))
diff <- diff %>% t() %>% as.data.frame()
diff$cell <- rownames(diff)
diff <- diff %>% select(cell,everything())
write.csv(diff,file = "/home/data/sdc/wt/model_data/cell_gene_exp_vs_normal_filter_random.csv",
          quote = F,row.names = F)

###和之前的相关性
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))
gene_exp$cell <- rownames(gene_exp)
gene_exp <- gene_exp %>% 
  tidyr::pivot_longer(cols = colnames(gene_exp)[1:(ncol(gene_exp)-1)],
                      names_to = "gene",values_to = "exp")

gene_exp_random <- readRDS("/home/data/sdc/wt/model_data/gene_exp_random.rds")
gene_exp_random$cell <- rownames(gene_exp_random)
gene_exp_random <- gene_exp_random %>% 
  tidyr::pivot_longer(cols = colnames(gene_exp_random)[1:(ncol(gene_exp_random)-1)],
                      names_to = "gene",values_to = "exp")
merge_dt <- inner_join(
  gene_exp,
  gene_exp_random %>% rename(random_exp=exp),
  by = c("cell","gene")
)
merge_dt %>% group_by(gene) %>% summarise(cor=cor(exp,random_exp)) -> tt
merge_dt %>% group_by(cell) %>% summarise(cor=cor(exp,random_exp)) -> tt
library(ggpubr)
p1 <- gghistogram(tt, x = "cor", 
                  fill = "#00AFBB",
                  add = "median", rug = TRUE, 
                  title = "Correlation of gene expression")+
  labs(y="Gene Counts",x="Correlation")
p1
ggsave("figs/gene_exp_cor_real_random.pdf",width = 7, height = 5)


