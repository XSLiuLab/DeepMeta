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
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

all_tissue <- list.files("data/met/origin_model/all_xml/")
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

cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>%
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))
cell_info <- cell_info %>% 
  filter(ModelID %in% rownames(gene_exp)) %>% 
  filter(ModelID %in% sanger$cell)
cell_info <- cell_info %>% 
  select(ModelID,OncotreeLineage) %>% 
  rename(cell = OncotreeLineage) %>% 
  left_join(.,net_cell_mapping)
cell_info <- na.omit(cell_info)
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
  .export = c("net_cell_mapping","cell_info","sanger","gene_exp",
              "ensg2name","get_dep"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  cell_cu <- cell_info$ModelID[i]
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(gsub(".xml","",cell_info$net[i]),
                                         "_enzymes_based_graph.tsv")))
  dep_cell <- sanger %>% filter(cell == cell_cu)
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
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_sanger/",
                              cell_cu,".txt"),sep = "\t",
                row.names = F)
    write.table(cell_net_dep,
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_sanger/",
                              cell_cu,"_feat.txt"),
                sep = "\t",row.names = F)
  }
}
parallel::stopCluster(cl = my.cluster)

write.csv(cell_info %>% select(ModelID) %>% rename(cell=ModelID),
          file = "/home/data/sdc/wt/model_data/new_model/enzyme_net_sanger/raw/cell_info.csv",
          quote = F,row.names = F)

sanger_cell <- read.csv("/home/data/sdc/wt/model_data/new_model/enzyme_net_sanger/raw/cell_info.csv")
sanger_cell$cell_index <- (1:nrow(sanger_cell))-1
write.csv(sanger_cell,
          file = "/home/data/sdc/wt/model_data/new_model/enzyme_net_sanger/raw/cell_info.csv",
          quote = F,row.names = F)

####RF
train_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/enzyme_net_sanger/raw/cell_info.csv")
res <- vector("list",nrow(train_dt))
for (i in seq_along(res)){
  dt <- data.table::fread(paste0("/home/data/sdc/wt/model_data/enzyme_net_sanger/",
                                 train_dt$cell[i],"_feat.txt"),
                          data.table = F)
  dt <- dt %>% select(id,is_dep)
  dt$cell <- train_dt$cell[i]
  res[[i]] <- dt
}
res <- bind_rows(res)
res <- na.omit(res)
saveRDS(res,file = "data/sanger_dt.rds")

sanger_dt <- readRDS("data/sanger_dt.rds")
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

sanger_feat <- data.frame(id = unique(sanger_dt$id)) %>% 
  rowwise() %>% 
  mutate(fea=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[2]]) %>% 
  ungroup() %>%
  tidyr::separate_wider_delim(cols = fea,delim = ",",names_sep="-") %>% 
  mutate(across(starts_with("fea-"),as.numeric)) %>% as.data.frame()

train_dt <- data.table::fread("/home/wt/meta_target/data/train_dtV2.csv",data.table = F)
sanger_dt <- sanger_dt %>% 
  filter(id %in% train_dt$id)
sanger_feat <- sanger_feat %>% 
  filter(id %in% sanger_dt$id)
write.csv(sanger_dt, file = "data/sanger_dt.csv",row.names = F)
write.csv(sanger_feat, file = "data/all_sanger_gene_cpg.csv",row.names = F)



