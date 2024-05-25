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
cell_info <- cell_info %>% 
  select(ModelID,OncotreeLineage) %>% 
  rename(cell = OncotreeLineage) %>% 
  left_join(.,net_cell_mapping)
cell_info <- na.omit(cell_info)
saveRDS(cell_info,file = "data/all_cell_info.rds")
###筛选在网络中的基因
res <- readRDS("data/all_net_genes.rds")
dep_dt <- dep_dt %>%
  mutate(index = paste(ModelID,gene,sep = "-")) %>%
  filter(index %in% res$index)

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
  .export = c("net_cell_mapping","cell_info","dep_dt","gene_exp",
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
      select(1:2,is_dep,everything())
    
    write.table(tissue_net,
                file = paste0("/home/data/sdb/wt/model_data/enzyme_train/",
                              cell,".txt"),sep = "\t",
                row.names = F)
    write.table(cell_net_dep,
                file = paste0("/home/data/sdb/wt/model_data/enzyme_train/",
                              cell,"_feat.txt"),
                sep = "\t",row.names = F)
  }
}
parallel::stopCluster(cl = my.cluster)

write.csv(cell_info %>% select(ModelID) %>% rename(cell=ModelID),
          file = "/home/data/sdc/wt/model_data/new_model/enzyme_net/raw/cell_info.csv",quote = F,row.names = F)

all_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/enzyme_net/raw/cell_info.csv")

set.seed(2023042602)
test_cell <- sample(all_dt$cell,size = nrow(all_dt)*0.1)
test_dt <- all_dt %>% filter(cell %in% test_cell)
train_dt <- all_dt %>% filter(!(cell %in% test_cell))
train_dt$cell_index <- (1:nrow(train_dt))-1
write.csv(train_dt,
          file = "/home/data/sdc/wt/model_data/new_model/cell_net_filter_exp/raw/train_cell_info.csv",
          quote = F,row.names = F)
test_dt$cell_index <- (1:nrow(test_dt))-1
write.csv(test_dt,
          file = "/home/data/sdc/wt/model_data/new_model/enzyme_net/test/raw/test_cell_info.csv",
          quote = F,row.names = F)

####
###训练集中有哪些基因
train_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/cell_net_filter_exp/raw/train_cell_info.csv")
res <- vector("list",nrow(train_dt))
for (i in seq_along(res)){
  dt <- data.table::fread(paste0("/home/data/sdc/wt/model_data/enzyme_net/",
                                 train_dt$cell[i],"_feat.txt"),
                          data.table = F)
  dt <- dt %>% select(id,is_dep)
  dt$cell <- train_dt$cell[i]
  res[[i]] <- dt
}
res <- bind_rows(res)
res <- na.omit(res)
saveRDS(res,file = "data/train_dtV2.rds")

##
ensg2name <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  sym <- tt_gene %>% 
    paste(.,collapse=",")
  return(sym)
}
test_dt <- read.csv("/home/data/sdc/wt/model_data/new_model/enzyme_net/test/raw/test_cell_info.csv")
res <- vector("list",nrow(test_dt))
for (i in seq_along(res)){
  tt <- data.table::fread(paste0("/home/data/sdc/wt/model_data/enzyme_net/",
                                 test_dt$cell[i],"_feat.txt"),
                          data.table = F)
  tt <- tt %>% select(id,is_dep)
  tt$cell <- test_dt$cell[i]
  res[[i]] <- tt
}
res <- bind_rows(res)
res <- na.omit(res)

all_gene <- data.frame(id=unique(res$id))
all_gene <- all_gene %>% 
  rowwise() %>% 
  mutate(gene = ensg2name(id, enz_gene_mapping)) %>% 
  ungroup()
test_dt <- left_join(res,all_gene)
saveRDS(test_dt, "data/test_dtV2.rds")
###training gene
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
train_dt <- readRDS("data/train_dtV2.rds")
all_gene <- data.frame(id=unique(train_dt$id))

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
train_dt <- left_join(train_dt,all_gene)
saveRDS(train_dt, "data/train_dtV2.rds")

train_dtV2 <- readRDS("data/train_dtV2.rds")
train_gene <- train_dtV2 %>% select(id) %>% distinct_all()
write.csv(train_gene,file = "data/train_genes.csv",
          quote = F,row.names = F)

###save test
test_dtV2 <- readRDS("data/test_dtV2.rds")
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

test_feat <- data.frame(id = unique(test_dtV2$id)) %>% 
  rowwise() %>% 
  mutate(fea=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[2]]) %>% 
  ungroup() %>%
  tidyr::separate_wider_delim(cols = fea,delim = ",",names_sep="-") %>% 
  mutate(across(starts_with("fea-"),as.numeric)) %>% as.data.frame()

write.csv(test_dtV2 %>% select(-gene),
          file = "data/test_dtV2.csv",row.names = F)
write.csv(test_feat,file = "data/all_test_gene_cpg.csv",row.names = F)


