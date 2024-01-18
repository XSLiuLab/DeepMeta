####预测永生化细胞
library(dplyr)
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

cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <-  cell_info %>% 
  filter((OncotreePrimaryDisease == "Non-Cancerous"))
dep_dt <- readRDS("/home/data/sdb/wt/model_data/dep_dt.rds")
dep_dt <- dep_dt %>%
  mutate(type = case_when(
    score>0.5 ~ 1,
    score<=0.5 ~ 0
  ))
cell_info <- cell_info %>% filter(ModelID %in% dep_dt$ModelID)
net_cell_mapping <- data.frame(net=NA,cell=unique(cell_info$OncotreeLineage))
net_cell_mapping$net <- c("kidney.xml","iProstateCancer1560.xml",
                          "bone_marrow.xml",NA)
net_cell_mapping <- na.omit(net_cell_mapping)
gene_exp <- data.table::fread("/home/data/sda/wt/update/data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

dep_dt <- dep_dt %>% filter(gene %in% colnames(cpg_gene))
cell_info <- cell_info %>% 
  select(ModelID,OncotreeLineage) %>% 
  rename(cell = OncotreeLineage) %>% 
  left_join(.,net_cell_mapping)
cell_info <- na.omit(cell_info)

cell_info <- cell_info %>% 
  filter(ModelID %in% rownames(gene_exp)) %>% 
  filter(ModelID %in% dep_dt$ModelID)
dep_dt <- dep_dt %>% filter(ModelID %in% cell_info$ModelID) %>% 
  rename(cell = ModelID)

###并行
library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  2, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(cell_info),
  .export = c("net_cell_mapping","cell_info","dep_dt","gene_exp",
              "ensg2name","get_dep"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  cell_cu <- cell_info$ModelID[i]
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(gsub(".xml","",cell_info$net[i]),
                                         "_enzymes_based_graph.tsv")))
  dep_cell <- dep_dt %>% filter(cell == cell_cu)
  if (nrow(dep_cell) == 0){
    return(NA)
  }else{
    cell_exp <- gene_exp[cell_cu,] %>% t() %>% as.data.frame()
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
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_normal/",cell_cu,".txt"),sep = "\t",
                row.names = F)
    write.table(cell_net_dep,
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_normal/",cell_cu,"_feat.txt"),
                sep = "\t",row.names = F)
  }
}
parallel::stopCluster(cl = my.cluster)

cell_info <- cell_info %>% select(ModelID) %>% rename(cell=ModelID)
cell_info$cell_index <- (1:nrow(cell_info))-1
write.csv(cell_info,
          file = "data/normal_cell_info.csv",
          quote = F,row.names = F)
###
cell_info <- read.csv("data/normal_cell_info.csv")
all_cells <- cell_info$cell
cell_info$normal <- c("Kidney - Cortex","Prostate")

gtex <- data.table::fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                          data.table = F,skip = 2) %>% 
  filter(!grepl("PAR",Name))
diff <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv",data.table = F)
all_genes <- colnames(diff)[2:ncol(diff)]

res <- vector("list",length(all_cells))
for (i in 1:length(all_cells)){
  which_normal <- cell_info$normal[which(cell_info$cell == all_cells[i])]
  ne <- gtex %>% select(2,which_normal) %>% 
    filter(Description %in% all_genes) %>% 
    rename(gene=Description)
  colnames(ne)[2] <- "exp"
  ne <- ne %>% 
    group_by(gene) %>% 
    summarise(exp=mean(exp)) %>% ungroup %>% as.data.frame()
  rownames(ne) <- ne$gene
  ne <- ne[all_genes,]
  colnames(ne)[2] <- all_cells[i]
  res[[i]] <- ne %>% select(-gene)
}

res <- bind_cols(res)
normal_exp <- as.matrix(res)

gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))
tumor_exp <- gene_exp[colnames(res),rownames(res)] %>% t() %>% as.matrix()

diff <- tumor_exp / (log2(normal_exp + 1.01))
diff <- diff %>% t() %>% as.data.frame()
diff$cell <- rownames(diff)
diff <- diff %>% select(cell,everything())
write.csv(diff,file = "/home/data/sdc/wt/model_data/cell_gene_exp_vs_normal_normal.csv",
          quote = F,row.names = F)



