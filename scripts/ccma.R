library(dplyr)

ccma_exp <- data.table::fread("/home/data/sdc/wt/model_data/CCMA_rnaseq_logcpm.csv",
                              data.table = F)
ccma_exp <- ccma_exp %>% 
  select(sample,gene,`log2(CPM+1)`) %>% 
  rename(exp=`log2(CPM+1)`) %>% 
  group_by(sample,gene) %>% 
  summarise(exp = mean(exp)) %>% ungroup() %>% 
  tidyr::pivot_wider(names_from = gene,values_from = exp)
saveRDS(ccma_exp,file = "/home/data/sdc/wt/model_data/ccma_all_exp.rds")

ccma_exp <- data.table::fread("/home/data/sdc/wt/model_data/CCMA_rnaseq_logcpm.csv",
                              data.table = F)
all_samples <- unique(ccma_exp$sample)
sample_info <- read.csv("data/Childhood Cancer Model Atlas_Diagnosis_表格.csv")
sample_info <- sample_info %>% 
  mutate(id = gsub("[-]","_",Cell.Line.ID)) %>% 
  mutate(id = gsub("[ ]\\(.+","",id)) %>% 
  mutate(id = gsub(" ","_",id))
sample_info <- sample_info %>% 
  filter(id %in% ccma_exp$sample)
control <- sample_info %>% filter(Diagnosis == "No malignancy, Brain")
brain_cns <- sample_info %>% filter(Category == "Brain/CNS")

###使用control 表达作为正常样本的表达得到diff 表达文件
diff_exp <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv",
                              data.table = F)
ccma_exp <- ccma_exp %>% 
  filter(gene %in% colnames(diff_exp)) %>% 
  select(sample,gene,`log2(CPM+1)`) %>% 
  rename(exp=`log2(CPM+1)`) %>% 
  group_by(sample,gene) %>% 
  summarise(exp = mean(exp)) %>% ungroup() %>% 
  tidyr::pivot_wider(names_from = gene,values_from = exp)
##没有的基因用 CCLE 的相同组织的加上
colnames(diff_exp)[which(!(colnames(diff_exp) %in% colnames(ccma_exp)))]
##NSG1
cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>% 
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))
brain_cell <- cell_info %>% 
  filter(OncotreeLineage == "CNS/Brain") %>% 
  filter(OncotreePrimaryDisease == "Embryonal Tumor")
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))
brain_cell_exp <- gene_exp[brain_cell$ModelID,] %>% select(NSG1) %>% 
  filter(!is.na(NSG1))

###中位数
ccma_exp$NSG1 <- rep(median(brain_cell_exp$NSG1,nrow(ccma_exp)))
ccma_exp <- ccma_exp %>% 
  select(sample,colnames(diff_exp)[2:ncol(diff_exp)])

###
tumor_exp <- ccma_exp %>% filter(sample %in% brain_cns$id) %>% as.data.frame()
normal_exp <- ccma_exp %>% filter(sample %in% control$id) %>% as.data.frame()
rownames(normal_exp) <- normal_exp$sample
normal_exp <- normal_exp %>% select(-sample)
normal_exp_median <- apply(normal_exp,2,median) %>% t() %>% as.data.frame()

rownames(tumor_exp) <- tumor_exp$sample
tumor_exp <- tumor_exp %>% select(-sample)
normal_exp_median <- normal_exp_median %>% select(colnames(tumor_exp))
normal_exp_median[1,] <- log2(2^(normal_exp_median[1,]) + 0.01)

diff_exp_ch <- apply(tumor_exp,1,function(x){x/normal_exp_median[1,]}) 
diff_exp_ch <- bind_rows(diff_exp_ch)
rownames(diff_exp_ch) <- rownames(tumor_exp)
diff_exp_ch$cell <- rownames(diff_exp_ch)
diff_exp_ch <- diff_exp_ch %>% select(colnames(diff_exp))
write.csv(diff_exp_ch,file = "/home/data/sdc/wt/model_data/ccma_gene_exp_vs_normal_filter.csv",
          quote = F,row.names = F)

cell_info <- data.frame(
  ModelID = diff_exp_ch$cell,
  net = "brain.xml"
)
saveRDS(cell_info,file = "data/ccma_cell_info.rds")
###########
ccma_exp <- readRDS("/home/data/sdc/wt/model_data/ccma_all_exp.rds") %>% as.data.frame()
rownames(ccma_exp) <- ccma_exp$sample
ccma_exp <- ccma_exp %>% select(-sample)
cell_info <- readRDS("data/ccma_cell_info.rds")

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

###并行
library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(cell_info),
  .export = c("cell_info","ccma_exp","ensg2name","enz_gene_mapping","cpg_gene"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  cell_cu <- cell_info$ModelID[i]
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(gsub(".xml","",cell_info$net[i]),
                                         "_enzymes_based_graph.tsv")))
  cell_exp <- ccma_exp[cell_cu,] %>% t() %>% as.data.frame()
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
  
  write.table(tissue_net,
              file = paste0("/home/data/sdc/wt/model_data/enzyme_net_ccma/",cell_cu,".txt"),sep = "\t",
              row.names = F)
  write.table(cell_net,
              file = paste0("/home/data/sdc/wt/model_data/enzyme_net_ccma/",cell_cu,"_feat.txt"),
              sep = "\t",row.names = F)
  
}
parallel::stopCluster(cl = my.cluster)

cell_info <- cell_info %>% select(ModelID) %>% rename(cell=ModelID)
cell_info$cell_index <- (1:nrow(cell_info))-1
write.csv(cell_info, file = "data/ccma_pre_info.csv", quote = F, row.names = F)

###
train_dt <- read.csv("data/ccma_pre_info.csv")
res <- vector("list",nrow(train_dt))
for (i in seq_along(res)){
  dt <- data.table::fread(paste0("/home/data/sdc/wt/model_data/enzyme_net_ccma/",
                                 train_dt$cell[i],"_feat.txt"),
                          data.table = F)
  dt <- dt %>% select(id)
  dt$cell <- train_dt$cell[i]
  res[[i]] <- dt
}
res <- bind_rows(res)
res <- na.omit(res)

ccma_dt <- res
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

ccma_feat <- data.frame(id = unique(ccma_dt$id)) %>% 
  rowwise() %>% 
  mutate(fea=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[2]]) %>% 
  ungroup() %>%
  tidyr::separate_wider_delim(cols = fea,delim = ",",names_sep="-") %>% 
  mutate(across(starts_with("fea-"),as.numeric)) %>% as.data.frame()

train_dt <- data.table::fread("/home/wt/meta_target/data/train_dtV2.csv",data.table = F)
ccma_dt <- ccma_dt %>% 
  filter(id %in% train_dt$id)
ccma_feat <- ccma_feat %>% 
  filter(id %in% ccma_dt$id)
write.csv(ccma_dt, file = "data/ccma_dt.csv",row.names = F)
write.csv(ccma_feat, file = "data/all_ccma_gene_cpg.csv",row.names = F)






