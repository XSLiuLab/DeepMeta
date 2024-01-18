library(dplyr)
library(ggpubr)
###rnai
rnai <- data.table::fread("/home/data/sdc/wt/model_data/rnai_gene_dependency.csv",
                          data.table = F)
colnames(rnai)[2:ncol(rnai)] <- gsub("\\s*\\([^\\)]+\\)","",
                                     colnames(rnai)[2:ncol(rnai)])
rnai <- rnai %>% 
  tidyr::pivot_longer(cols = colnames(rnai)[2:ncol(rnai)],names_to = "gene",
                      values_to = "score")
colnames(rnai)[1] <- "cell"

###探索下和crispr 的一致性
dep_dt <- readRDS("/home/data/sdb/wt/model_data/dep_dt.rds")
merge_dt <- inner_join(
  dep_dt %>% rename(cell = ModelID, cri_score = score),
  rnai %>% rename(rnai_score = score),
  by = c("cell","gene")
) %>% na.omit(.)
merge_dt %>% group_by(gene) %>% summarise(cor=cor(cri_score,rnai_score)) -> tt
p1 <- gghistogram(tt, x = "cor", 
                  fill = "#00AFBB",
                  add = "median", rug = TRUE, 
                  title = "Correlation of dependency score bewtten RNAi and CRISPR")+
  labs(y="Gene Counts",x="Correlation")

merge_dt <- merge_dt %>% 
  mutate(cri_type = case_when(
    cri_score > 0.5 ~ "Positive",
    cri_score <= 0.5 ~ "Negative",
    TRUE ~ "None"
  )) %>% 
  mutate(rnai_type = case_when(
    rnai_score > 0.5 ~ "Positive",
    rnai_score <= 0.5 ~ "Negative",
    TRUE ~ "None"
  ))

library(ggalluvial)
merge_summ <- merge_dt %>% 
  group_by(cri_type) %>% 
  summarise(Negative = sum(rnai_type == "Negative"),
            Positive = sum(rnai_type == "Positive"),
            None = sum(rnai_type == "None")) %>% 
  ungroup() %>% 
  tidyr::pivot_longer(cols = c("Negative","Positive","None"),names_to = "Rnai",
                      values_to = "freq")
p2 <- ggplot(data = merge_summ,
             aes(axis1 = cri_type, axis2 = Rnai, y = freq)) +
  geom_alluvium(aes(fill = Rnai),show.legend=T) +
  geom_stratum() +
  geom_text(stat = "stratum",
            aes(label = after_stat(stratum))) +
  theme_void()+
  labs(title = "Change of gene label between CRISPR and RNAi")
p1 + p2
ggsave("figs/crispr_rnai_compare.pdf",width = 15,height = 6)

rnai <- rnai %>% 
  filter(!is.na(score)) %>%
  mutate(type = case_when(
    score>0.5 ~ 1,
    score<=0.5 ~ 0
  ))

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
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

cell_info <- cell_info %>% 
  filter(ModelID %in% rownames(gene_exp)) %>% 
  filter(ModelID %in% rnai$cell)
cell_info <- cell_info %>% 
  select(ModelID,OncotreeLineage) %>% 
  rename(cell = OncotreeLineage) %>% 
  left_join(.,net_cell_mapping)
cell_info <- cell_info %>% na.omit()
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
  .export = c("net_cell_mapping","cell_info","rnai","gene_exp",
              "ensg2name","get_dep"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  cell_cu <- cell_info$ModelID[i]
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(gsub(".xml","",cell_info$net[i]),
                                         "_enzymes_based_graph.tsv")))
  dep_cell <- rnai %>% filter(cell == cell_cu)
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
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_rnai/",cell_cu,".txt"),sep = "\t",
                row.names = F)
    write.table(cell_net_dep,
                file = paste0("/home/data/sdc/wt/model_data/enzyme_net_rnai/",cell_cu,"_feat.txt"),
                sep = "\t",row.names = F)
  }
}
parallel::stopCluster(cl = my.cluster)

write.csv(cell_info %>% select(ModelID) %>% rename(cell=ModelID),
          file = "/home/data/sdc/wt/model_data/new_model/enzyme_net_rnai/raw/cell_info.csv",quote = F,row.names = F)

rnai <- read.csv("/home/data/sdc/wt/model_data/new_model/enzyme_net_rnai/raw/cell_info.csv")
rnai$cell_index <- (1:nrow(rnai))-1
write.csv(rnai,
          file = "/home/data/sdc/wt/model_data/new_model/enzyme_net_rnai/raw/cell_info.csv",
          quote = F,row.names = F)

####

