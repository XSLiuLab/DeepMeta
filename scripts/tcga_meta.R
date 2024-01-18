library(dplyr)
####相对于正常样本的表达
###exp log2(x+0.001)
tpm <- data.table::fread("/home/data/sdc/wt/TCGA/tcga_RSEM_gene_tpm.gz",data.table = F)
tpm <- tpm %>% 
  select(sample,which(as.numeric(substr(colnames(tpm),14,15)) < 10)) %>% 
  rename(id = sample)
###mapping
mapping <- data.table::fread("/home/data/sdc/wt/TCGA/probeMap_gencode.v23.annotation.gene.probemap",
                             data.table = F)
tpm <- left_join(tpm,mapping %>% select(id,gene)) %>% 
  select(-id) %>% 
  select(gene,everything())

tpm <- tpm[!duplicated(tpm$gene),]
rownames(tpm) <- tpm$gene
tpm <- tpm %>% select(-gene)
##转化成 log2(x+1)
tpm[,1:ncol(tpm)] <- apply(tpm[,1:ncol(tpm)],2,function(x){log2((2^x)-0.001 + 1)})
tpm <- t(tpm) %>% as.data.frame()
tpm$cell <- rownames(tpm)
tpm <- tpm %>% select(cell,everything())
saveRDS(tpm,file = "/home/data/sdc/wt/TCGA/tpm_log2_1.rds")

###癌症样本-癌症类型--net--normal
library(parallel)
all_samples <- data.frame(samples = unique(tpm$cell)) %>% 
  mutate(cancer = EasyBioinfo::get_cancer_type(samples,parallel = T,cores=30))
cancer_net <- data.frame(cancer = unique(all_samples$cancer))
cancer_net$net <- c("brain","brain","iCervicalCancer1611","iLungCancer1490",
                    "iColorectalCancer1750","bone_marrow","iBreastCancer1746",
                    "iStomachCancer1511",NA,"kidney","iStomachCancer1511",
                    "iProstateCancer1560","iSkinCancer1386","iEndometrialCancer1713",
                    "iHeadNeckCancer1628","kidney","iLungCancer1490",NA,
                    "iLiverCancer1788","iThyroidCancer1710",NA,"iColorectalCancer1750",
                    "iPancreaticCancer1613","iOvarianCancer1620","iTestisCancer1483",
                    NA,NA,NA,"iUrothelialCancer1532","kidney","kidney",NA,NA)
cancer_net <- na.omit(cancer_net)
cancer_net$normal <- c("Brain - Amygdala","Brain - Amygdala","Cervix - Endocervix",
                       "Lung","Colon - Sigmoid","Whole Blood",
                       "Breast - Mammary Tissue","Esophagus - Gastroesophageal Junction",
                       "Kidney - Cortex",
                       "Stomach","Prostate","Skin - Sun Exposed (Lower leg)",
                       "Uterus","Whole Blood","Kidney - Cortex","Lung","Liver",
                       "Thyroid","Colon - Sigmoid","Pancreas","Ovary","Testis",
                       "Bladder","Kidney - Cortex","Kidney - Cortex")
all_samples <- left_join(all_samples,cancer_net)
all_samples <- na.omit(all_samples)
saveRDS(all_samples,file = "data/tcga_all_sampleinfo.rds")

###
tcga_sample <- readRDS("data/tcga_all_sampleinfo.rds")
train_cell_info <- readRDS("data/all_cell_info.rds")
tcga_sample <- tcga_sample %>% 
  filter(net %in% gsub(".xml","",train_cell_info$net))
###only gene in net
id2gene <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  sym <- tt_gene %>% 
    paste(.,collapse=",")
  return(sym)
}
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")

res <- vector("list",20)
for (i in 1:20){
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",all_net[i]))
  dt <- data.frame(id = unique(c(tissue_net$from,tissue_net$to)))%>% 
    rowwise() %>% 
    mutate(gene=id2gene(id,enz_gene_mapping)) %>% 
    ungroup()
  res[[i]] <- dt
}
res <- bind_rows(res)
res <- res %>% distinct_all()
all_genes <- paste(res$gene,collapse = ",") %>% strsplit(.,",") %>% `[[`(1) %>% unique()
tpm_net <- tpm %>% select(cell,any_of(all_genes))
saveRDS(tpm_net,file = "/home/data/sdc/wt/TCGA/tpm_log2_1_net.rds")

###normal
tcga_sample <- readRDS("data/tcga_all_sampleinfo.rds")
gtex <- data.table::fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                          data.table = F,skip = 2) %>% 
  filter(!grepl("PAR",Name))
tt <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv")
all_sample <- tcga_sample$samples
all_genes <- colnames(tt)[2:ncol(tt)]

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
  i = all_sample,
  .export = c("tcga_sample","gtex","all_genes"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  which_normal <- tcga_sample$normal[which(tcga_sample$samples == i)]
  ne <- gtex %>% select(2,which_normal) %>% 
    filter(Description %in% all_genes) %>% 
    rename(gene=Description)
  colnames(ne)[2] <- "exp"
  ne <- ne %>% 
    group_by(gene) %>% 
    summarise(exp=mean(exp)) %>% ungroup %>% as.data.frame()
  rownames(ne) <- ne$gene
  ne <- ne[all_genes,]
  colnames(ne)[2] <- i
  ne <- ne %>% select(-gene)
  return(ne)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_cols(res)
normal_exp <- as.matrix(res)
saveRDS(normal_exp, file = "/home/data/sdc/wt/TCGA/normal_exp.rds")

###有一些基因在 TCGA 和 细胞系中的名字不一样，需要矫正
tcga_sample <- readRDS("data/tcga_all_sampleinfo.rds")
tt <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv")
all_sample <- tcga_sample$samples
all_genes <- colnames(tt)[2:ncol(tt)]

tpm <- readRDS("/home/data/sdc/wt/TCGA/tpm_log2_1.rds")
tumor_exp <- tpm %>% select(-cell) %>% t() %>% as.data.frame()
all_tcga_genes <- data.frame(gene=rownames(tumor_exp))
saveRDS(all_tcga_genes,file = "data/all_tcga_genes.rds")

mapping <- data.table::fread("/home/data/sdc/wt/TCGA/probeMap_gencode.v23.annotation.gene.probemap",
                             data.table = F)
update_name <- all_genes[which(!(all_genes %in% rownames(tumor_exp)))]
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
###下载基因名对应 https://www.genenames.org/download/custom/
sys_name <- data.table::fread("data/gene_name_als.txt",data.table = F)
sys_name <- sys_name %>% filter(`Approved symbol` %in% update_name) 
get_sys <- function(x){
  tt <- strsplit(x,split = ", ")[[1]]
  tt <- tt[which(tt %in% rownames(tumor_exp))]
  if (length(tt) != 0){
    return(tt)
  }else{
    return(NA)
  }
}
sys_name <- sys_name %>% 
  rowwise() %>% 
  mutate(tcga_symbol = get_sys(`Previous symbols`)) %>% ungroup()
remain <- enz_gene_mapping %>% 
  filter(symbol %in% sys_name$`Approved symbol`[which(is.na(sys_name$tcga_symbol))]) %>% 
  filter(ensembl_id %in% gsub("[.].+","",mapping$id)) %>%
  left_join(.,
            mapping %>% rename(ensembl_id = id) %>% 
              mutate(ensembl_id = gsub("[.].+","",ensembl_id)))
sys_name <- left_join(
  sys_name,
  remain %>% select(symbol,gene) %>% rename(`Approved symbol` = symbol)
) %>% mutate(tcga_symbol = ifelse(is.na(tcga_symbol),gene,tcga_symbol))

##剩下的手动查一下 genecard
sys_name$tcga_symbol[13] <- "CBS"
sys_name$tcga_symbol[75] <- "SMIM11"
sys_name$tcga_symbol[88] <- "CH507-152C13.1"
sys_name <- sys_name %>% select(`Approved symbol`,tcga_symbol) %>% 
  rename(cell_gene = `Approved symbol`)
###还少一个
sys_name <- bind_rows(sys_name,
                      data.frame(cell_gene=c("CCDC189"),
                                 tcga_symbol=c("C16orf93")))
origin_name <- data.frame(
  tcga_symbol = rownames(tumor_exp)
) %>% left_join(sys_name) %>% 
  mutate(cell_gene = ifelse(is.na(cell_gene),tcga_symbol,cell_gene))
rownames(tumor_exp) <- origin_name$cell_gene
tumor_exp <- tumor_exp[all_genes,all_sample] %>% as.matrix()
rownames(tumor_exp)[4700] <- "CBS"
rownames(tumor_exp)[7957] <- "CBSL"

normal_exp <- readRDS("/home/data/sdc/wt/TCGA/normal_exp.rds")
which(rownames(normal_exp) != rownames(tumor_exp))
diff <- tumor_exp / (log2(normal_exp + 1.01))
diff <- diff %>% t() %>% as.data.frame()
diff$cell <- rownames(diff)
diff <- diff %>% select(cell,everything())
write.csv(diff,file = "/home/data/sdc/wt/TCGA/tcga_gene_exp_vs_normal.csv",
          quote = F,row.names = F)

###pre
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
tpm_net <- readRDS("/home/data/sdc/wt/TCGA/tpm_log2_1_net.rds")
sample_info <- readRDS("data/tcga_all_sampleinfo.rds")

tpm_net <- tpm_net %>% select(-cell)
#create the cluster
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(sample_info),
  .export = c("sample_info","tpm_net","ensg2name","enz_gene_mapping",
              "cpg_gene"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(sample_info$net[i],
                                         "_enzymes_based_graph.tsv")))
  cell_exp <- tpm_net[sample_info$samples[i],] %>% t() %>% as.data.frame()
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
              file = paste0("/home/data/sdc/wt/model_data/tcga_net2/",
                            sample_info$samples[i],".txt"),
              sep = "\t",row.names = F)
  write.table(cell_net,
              file = paste0("/home/data/sdc/wt/model_data/tcga_net2/",
                            sample_info$samples[i],"_feat.txt"),
              sep = "\t",row.names = F)
}
parallel::stopCluster(cl = my.cluster)

sample_info <- sample_info %>% select(samples) %>% rename(cell=samples) 
sample_info$cell_index <- (1:nrow(sample_info))-1
write.csv(sample_info,
          file = "data/tcga_val_sample_info.csv",
          quote = F,row.names = F)
###TCGA-25-1870-01 这个样本大部分基因不表达
sample_info <- read.csv("data/tcga_val_sample_info.csv")
sample_info <- sample_info %>% filter(cell != "TCGA-25-1870-01")
sample_info$cell_index <- (1:nrow(sample_info))-1
write.csv(sample_info,
          file = "data/tcga_val_sample_info.csv",
          quote = F,row.names = F)

###res
pre <- data.table::fread("/home/data/sdc/wt/model_data/tcga_preV2.csv",
                         data.table = F)
pre_gene <-  data.frame(id = unique(c(pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]])
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(gene_name=id))
pre$cancer <- EasyBioinfo::get_cancer_type(pre$cell,parallel = T,cores = 30)
saveRDS(pre,file = "data/tcga_pre_V2.rds")

####预测 TCGA 正常样本
tpm <- data.table::fread("/home/data/sdc/wt/TCGA/tcga_RSEM_gene_tpm.gz",data.table = F)
tpm <- tpm %>% 
  select(sample,which(as.numeric(substr(colnames(tpm),14,15)) == 11)) %>% 
  rename(id = sample)
###mapping
mapping <- data.table::fread("/home/data/sdc/wt/TCGA/probeMap_gencode.v23.annotation.gene.probemap",
                             data.table = F)
tpm <- left_join(tpm,mapping %>% select(id,gene)) %>% 
  select(-id) %>% 
  select(gene,everything())

tpm <- tpm[!duplicated(tpm$gene),]
rownames(tpm) <- tpm$gene
tpm <- tpm %>% select(-gene)
##转化成 log2(x+1)
tpm[,1:ncol(tpm)] <- apply(tpm[,1:ncol(tpm)],2,function(x){log2((2^x)-0.001 + 1)})
tpm <- t(tpm) %>% as.data.frame()
tpm$cell <- rownames(tpm)
tpm <- tpm %>% select(cell,everything())
saveRDS(tpm,file = "/home/data/sdc/wt/TCGA/tpm_log2_1_normal.rds")
###癌症样本-癌症类型--net--normal
library(parallel)
all_samples <- data.frame(samples = unique(tpm$cell)) %>% 
  mutate(cancer = EasyBioinfo::get_cancer_type(samples,parallel = T,cores=30))
tumor <- readRDS("/home/data/sdc/wt/TCGA/tpm_log2_1.rds")
all_samples_tumor <- data.frame(samples = unique(tumor$cell)) %>% 
  mutate(cancer = EasyBioinfo::get_cancer_type(samples,parallel = T,cores=30))
cancer_net <- data.frame(cancer = unique(all_samples_tumor$cancer))
cancer_net <- cancer_net %>% filter(cancer %in% all_samples$cancer)

all_net <- list.files("data/meta_net/EnzGraphs/",pattern = "based_graph.tsv")
cancer_net$net <- c("brain","iCervicalCancer1611","iLungCancer1490",
                    "iColorectalCancer1750","iBreastCancer1746",
                    "iStomachCancer1511",NA,"kidney","iStomachCancer1511",
                    "iProstateCancer1560","iSkinCancer1386","iEndometrialCancer1713",
                    "iHeadNeckCancer1628","kidney","iLungCancer1490",NA,
                    "iLiverCancer1788","iThyroidCancer1710","iColorectalCancer1750",
                    "iPancreaticCancer1613",
                    NA,"iUrothelialCancer1532","kidney",NA)
cancer_net <- na.omit(cancer_net)
cancer_net$normal <- c("Brain - Amygdala","Cervix - Endocervix",
                       "Lung","Colon - Sigmoid",
                       "Breast - Mammary Tissue","Esophagus - Gastroesophageal Junction",
                       "Kidney - Cortex",
                       "Stomach","Prostate","Skin - Sun Exposed (Lower leg)",
                       "Uterus","Whole Blood","Kidney - Cortex","Lung","Liver",
                       "Thyroid","Colon - Sigmoid","Pancreas",
                       "Bladder","Kidney - Cortex")
all_samples <- left_join(
  all_samples,cancer_net
)
all_samples <- na.omit(all_samples)
saveRDS(all_samples,file = "data/tcga_all_sampleinfo_normal.rds")

tcga_sample <- readRDS("data/tcga_all_sampleinfo_normal.rds")
train_cell_info <- readRDS("data/all_cell_info.rds")
tcga_sample <- tcga_sample %>% filter(net %in% gsub(".xml","",train_cell_info$net))
###only gene in net
id2gene <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  sym <- tt_gene %>% 
    paste(.,collapse=",")
  return(sym)
}
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")

res <- vector("list",20)
for (i in 1:20){
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",all_net[i]))
  dt <- data.frame(id = unique(c(tissue_net$from,tissue_net$to)))%>% 
    rowwise() %>% 
    mutate(gene=id2gene(id,enz_gene_mapping)) %>% 
    ungroup()
  res[[i]] <- dt
}
res <- bind_rows(res)
res <- res %>% distinct_all()
all_genes <- paste(res$gene,collapse = ",") %>% strsplit(.,",") %>% `[[`(1) %>% unique()
tpm_net <- tpm %>% select(cell,any_of(all_genes))
saveRDS(tpm_net,file = "/home/data/sdc/wt/TCGA/tpm_log2_1_net_normal.rds")

#####
gtex <- data.table::fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                          data.table = F,skip = 2) %>% 
  filter(!grepl("PAR",Name))
tt <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv")
all_sample <- all_samples$samples
all_genes <- colnames(tt)[2:ncol(tt)]

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
  i = all_sample,
  .export = c("all_samples","gtex","all_genes"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  which_normal <- all_samples$normal[which(all_samples$samples == i)]
  ne <- gtex %>% select(2,which_normal) %>% 
    filter(Description %in% all_genes) %>% 
    rename(gene=Description)
  colnames(ne)[2] <- "exp"
  ne <- ne %>% 
    group_by(gene) %>% 
    summarise(exp=mean(exp)) %>% ungroup %>% as.data.frame()
  rownames(ne) <- ne$gene
  ne <- ne[all_genes,]
  colnames(ne)[2] <- i
  ne <- ne %>% select(-gene)
  return(ne)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_cols(res)
normal_exp <- as.matrix(res)

tumor_exp <- tpm %>% select(-cell) %>% t() %>% as.data.frame()
update_name <- all_genes[which(!(all_genes %in% rownames(tumor_exp)))]
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
sys_name <- data.table::fread("data/gene_name_als.txt",data.table = F)
sys_name <- sys_name %>% filter(`Approved symbol` %in% update_name) 
get_sys <- function(x){
  tt <- strsplit(x,split = ", ")[[1]]
  tt <- tt[which(tt %in% rownames(tumor_exp))]
  if (length(tt) != 0){
    return(tt)
  }else{
    return(NA)
  }
}
sys_name <- sys_name %>% 
  rowwise() %>% 
  mutate(tcga_symbol = get_sys(`Previous symbols`)) %>% ungroup()
remain <- enz_gene_mapping %>% 
  filter(symbol %in% sys_name$`Approved symbol`[which(is.na(sys_name$tcga_symbol))]) %>% 
  filter(ensembl_id %in% gsub("[.].+","",mapping$id)) %>%
  left_join(.,
            mapping %>% rename(ensembl_id = id) %>% 
              mutate(ensembl_id = gsub("[.].+","",ensembl_id)))
sys_name <- left_join(
  sys_name,
  remain %>% select(symbol,gene) %>% rename(`Approved symbol` = symbol)
) %>% mutate(tcga_symbol = ifelse(is.na(tcga_symbol),gene,tcga_symbol))

sys_name$tcga_symbol[13] <- "CBS"
sys_name$tcga_symbol[75] <- "SMIM11"
sys_name$tcga_symbol[88] <- "CH507-152C13.1"
sys_name <- sys_name %>% select(`Approved symbol`,tcga_symbol) %>% 
  rename(cell_gene = `Approved symbol`)
sys_name <- bind_rows(sys_name,
                      data.frame(cell_gene="CCDC189",
                                 tcga_symbol="C16orf93"))
origin_name <- data.frame(
  tcga_symbol = rownames(tumor_exp)
) %>% left_join(sys_name) %>% 
  mutate(cell_gene = ifelse(is.na(cell_gene),tcga_symbol,cell_gene))
rownames(tumor_exp) <- origin_name$cell_gene
tumor_exp <- tumor_exp[all_genes,all_sample] %>% as.matrix()
rownames(tumor_exp)[4700] <- "CBS"
rownames(tumor_exp)[7957] <- "CBSL"

diff <- tumor_exp / (log2(normal_exp + 1.01))
diff <- diff %>% t() %>% as.data.frame()
diff$cell <- rownames(diff)
diff <- diff %>% select(cell,everything())

write.csv(diff,file = "/home/data/sdc/wt/model_data/tcga_gene_exp_vs_normal_normal.csv",
          quote = F,row.names = F)
####
rm(list = ls())
tpm_net <- readRDS("/home/data/sdc/wt/TCGA/tpm_log2_1_net_normal.rds")
sample_info <- readRDS("data/tcga_all_sampleinfo_normal.rds")
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
cpg_gene <- readRDS("data/cpg_gene.rds")

tpm_net <- tpm_net %>% select(-cell)
#create the cluster
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(sample_info),
  .export = c("sample_info","tpm_net","ensg2name","enz_gene_mapping",
              "cpg_gene"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(sample_info$net[i],
                                         "_enzymes_based_graph.tsv")))
  cell_exp <- tpm_net[sample_info$samples[i],] %>% t() %>% as.data.frame()
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
              file = paste0("/home/data/sdc/wt/model_data/tcga_net_normal2/",
                            sample_info$samples[i],".txt"),
              sep = "\t",row.names = F)
  write.table(cell_net,
              file = paste0("/home/data/sdc/wt/model_data/tcga_net_normal2/",
                            sample_info$samples[i],"_feat.txt"),
              sep = "\t",row.names = F)
}
parallel::stopCluster(cl = my.cluster)

sample_info <- sample_info %>% select(samples) %>% rename(cell=samples) 
sample_info$cell_index <- (1:nrow(sample_info))-1
write.csv(sample_info,
          file = "data/tcga_normal_sample_info.csv",
          quote = F,row.names = F)

####
pre <- data.table::fread("data/tcga_normal_preV2.csv") %>%
  select(-V1)
pre_gene <-  data.frame(id = unique(c(pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]])
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(gene_name=id))
pre$cancer <- EasyBioinfo::get_cancer_type(pre$cell,parallel = T,cores = 30)
saveRDS(pre,file = "data/tcga_normal_pre_V2.rds")





