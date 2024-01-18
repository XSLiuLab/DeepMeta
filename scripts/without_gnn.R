library(dplyr)

train_dtV2 <- readRDS("data/train_dtV2.rds")
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

train_feat <- data.frame(id = unique(train_dtV2$id)) %>% 
  rowwise() %>% 
  mutate(fea=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[2]]) %>% 
  ungroup() %>%
  tidyr::separate_wider_delim(cols = fea,delim = ",",names_sep="-") %>% 
  mutate(across(starts_with("fea-"),as.numeric)) %>% as.data.frame()

write.csv(train_dtV2 %>% select(-gene),
          file = "data/train_dtV2.csv",row.names = F)
write.csv(train_feat,file = "data/all_train_gene_cpg.csv",row.names = F)


