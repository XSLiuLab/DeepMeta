library(dplyr)
library(stringr)

nc_dt <- readxl::read_xlsx("data/nc_data1.xlsx",skip = 1)
nc_drug <- readxl::read_xlsx("data/nc_data2.xlsx",sheet = "AUC")
nc_cells <- c("ACH-000045","ACH-000362","ACH-000006","ACH-000301","ACH-000074",
              "ACH-000983","ACH-000168","ACH-002273","ACH-000004","ACH-000002",
              "ACH-000386","ACH-000551","ACH-000432","ACH-001647","ACH-000146")
colnames(nc_drug) <- c("drug",nc_cells)
cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>% 
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))

gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

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

cell_info <- cell_info %>% filter(ModelID %in% nc_cells)
net_cell_mapping <- net_cell_mapping %>% 
  filter(cell %in% cell_info$OncotreeLineage)
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
  i = 1:nrow(net_cell_mapping),
  .export = c("net_cell_mapping","cell_info","gene_exp",
              "ensg2name"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(gsub(".xml","",net_cell_mapping$net[i]),
                                         "_enzymes_based_graph.tsv")))
  tissue_cell <- cell_info %>% filter(OncotreeLineage == net_cell_mapping$cell[i])
  for (j in tissue_cell$ModelID){
    cell_exp <- gene_exp[j,] %>% t() %>% as.data.frame()
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
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_drug/",j,".txt"),sep = "\t",row.names = F)
    write.table(cell_net,
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_drug/",j,"_feat.txt"),
                sep = "\t",row.names = F)
  }
}
parallel::stopCluster(cl = my.cluster)

cell_info <- cell_info %>% select(ModelID) %>% rename(cell=ModelID) %>% 
  filter(cell %in% rownames(gene_exp))
cell_info$cell_index <- (1:nrow(cell_info))-1
write.csv(cell_info,
          file = "data/nc_drug_15_info.csv",
          quote = F,row.names = F)

####
###
drug_dt <- read.csv("data/nc_drug_15_info.csv")
res <- vector("list",nrow(drug_dt))
for (i in seq_along(res)){
  dt <- data.table::fread(paste0("/home/data/sdc/wt/model_data/enzyme_net_drug/",
                                 drug_dt$cell[i],"_feat.txt"),
                          data.table = F)
  dt <- dt %>% select(id)
  dt$cell <- drug_dt$cell[i]
  res[[i]] <- dt
}
res <- bind_rows(res)
res <- na.omit(res)

drug_dt <- res
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

drug_feat <- data.frame(id = unique(drug_dt$id)) %>% 
  rowwise() %>% 
  mutate(fea=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[2]]) %>% 
  ungroup() %>%
  tidyr::separate_wider_delim(cols = fea,delim = ",",names_sep="-") %>% 
  mutate(across(starts_with("fea-"),as.numeric)) %>% as.data.frame()

train_dt <- data.table::fread("/home/wt/meta_target/data/train_dtV2.csv",data.table = F)
drug_dt <- drug_dt %>% 
  filter(id %in% train_dt$id)
drug_feat <- drug_feat %>% 
  filter(id %in% drug_dt$id)
write.csv(drug_dt, file = "data/drug_dt.csv",row.names = F)
write.csv(drug_feat, file = "data/all_drug_gene_cpg.csv",row.names = F)








