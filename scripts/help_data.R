library(dplyr)

###
dep_dt <- data.table::fread("/home/data/sdb/wt/model_data/CRISPRGeneDependency.csv",
                            data.table = F)
colnames(dep_dt)[2:ncol(dep_dt)] <- gsub("\\s*\\([^\\)]+\\)","",
                                         colnames(dep_dt)[2:ncol(dep_dt)])
dep_dt <- dep_dt %>% 
  tidyr::pivot_longer(cols = colnames(dep_dt)[2:ncol(dep_dt)],names_to = "gene",
                      values_to = "score")
saveRDS(dep_dt,file = "/home/data/sdb/wt/model_data/dep_dt.rds")

###基因名转化
ensg2name <- function(ensg,mapping,exp){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  sym <- tt_gene %>% 
    paste(.,collapse=",")
  exp_tmp <- exp$exp[which(exp$gene %in% tt_gene)] %>% max()
  return(list(sym,exp_tmp))
}

enz_ensm <- org.Hs.eg.db::org.Hs.egENSEMBL
enz_sym <- org.Hs.eg.db::org.Hs.egSYMBOL
mapped_genes <- AnnotationDbi::mappedkeys(enz_ensm)
mapping_enz <- enz_ensm[mapped_genes] %>% as.data.frame()
mapped_genes <- AnnotationDbi::mappedkeys(enz_sym)
mapping_sym <- enz_sym[mapped_genes] %>% as.data.frame()
mapping <- inner_join(mapping_enz,mapping_sym)
saveRDS(mapping,file = "data/enz_gene_mapping.rds")

###生成酶网络
library(Met2Graph)
infiles <- list.files("data/GSM_xml/",full.names = T)
outDir <- "data/meta_net/"
for (i in infiles){
  Met2EnzGraph(i, rmMets=TRUE, outDir=outDir, outFormat="ncol")
}

####
cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>%
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))
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

# ###看看有多少酶
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
res <- vector("list",19)
for (i in seq_along(res)){
  dt <- read.table(paste0("data/meta_net/EnzGraphs/",
                          paste0(gsub(".xml","",net_cell_mapping$net[i]),
                                 "_enzymes_based_graph.tsv")))
  dt_gene <- data.frame(gene=strsplit(paste0(unique(c(dt$from,dt$to)),
                                             collapse = " and "),
                                      split = " and ")[[1]])
  res[[i]] <- dt_gene
}
res <- bind_rows(res)
all_gene <- enz_gene_mapping$symbol[which(enz_gene_mapping$ensembl_id %in% res$gene)] %>%
  unique()

# ###生成基因的特征
cpg <- fgsea::gmtPathways("data/c2.cgp.v2023.1.Hs.symbols.gmt")
res <- vector("list",length(all_gene))
for (i in seq_along(res)){
  dt <- sapply(cpg,function(x){all_gene[i] %in% x}) %>% as.data.frame()
  colnames(dt) <- all_gene[i]
  dt[,1] <- as.numeric(dt[,1])
  res[[i]] <- dt
}
res <- bind_cols(res)
dep_dt <- readRDS("/home/data/sdb/wt/model_data/dep_dt.rds")
res <- res %>% dplyr::select(any_of(unique(dep_dt$gene)))
res_summ <- apply(res,1,sum) %>% as.data.frame()
res_summ$pathway <- rownames(res_summ)
colnames(res_summ)[1] <- "gene_num"
res_summ <- res_summ %>% filter(gene_num > 0)
res <- res[res_summ$pathway,]
saveRDS(res,file = "data/cpg_gene.rds")

###
# ###细胞系中的基因
cell_info <- cell_info %>% 
  select(ModelID,OncotreeLineage) %>% 
  rename(cell = OncotreeLineage) %>% 
  left_join(.,net_cell_mapping)
cell_info <- na.omit(cell_info)

enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
ensg2name <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>%
    unique()
  sym <- tt_gene %>%
    paste(.,collapse=",")
  return(sym)
}

my.cluster <- parallel::makeCluster(
  40,
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

library(foreach)
res <- foreach(
  i = 1:nrow(cell_info),
  .export = c("enz_gene_mapping","cell_info",
              "ensg2name"),
  .packages = c("dplyr","tidyr")
) %dopar% {
  cell <- cell_info$ModelID[i]
  tissue_net <- read.table(paste0("data/meta_net/EnzGraphs/",
                                  paste0(gsub(".xml","",cell_info$net[i]),
                                         "_enzymes_based_graph.tsv")))

  all_gene <- data.frame(id=unique(c(tissue_net$from,tissue_net$to))) %>%
    rowwise() %>%
    mutate(gene = ensg2name(id, enz_gene_mapping)) %>%
    ungroup()
  all_gene$cell <- cell
  return(all_gene)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
res <- res %>% select(gene,cell) %>%
  tidyr::separate_longer_delim(cols = "gene",delim = ",")
res <- res %>% distinct_all()
res <- res %>%
  mutate(index = paste(cell,gene,sep = "-"))
saveRDS(res,file = "data/all_net_genes.rds")

###正常组织的表达
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>%
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))

gtex <- data.table::fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                          data.table = F,skip = 2) %>% 
  filter(!grepl("PAR",Name))

net_cell_mapping$normal <- c("Ovary","Whole Blood","Colon - Sigmoid",
                             "Skin - Sun Exposed (Lower leg)","Bladder",
                             "Lung","Kidney - Cortex","Breast - Mammary Tissue",
                             "Pancreas","Brain - Amygdala","Stomach","Thyroid","Prostate",
                             "Whole Blood","Uterus","Liver","Cervix - Endocervix",
                             "Adrenal Gland","Testis")
cell_info <- cell_info %>% 
  filter(OncotreeLineage %in% net_cell_mapping$cell) %>% 
  filter(ModelID %in% rownames(gene_exp))
all_cells <- cell_info$ModelID
all_genes <- intersect(colnames(gene_exp),gtex$Description)

res <- vector("list",length(all_cells))
for (i in 1:length(all_cells)){
  which_tissue <- cell_info$OncotreeLineage[which(cell_info$ModelID==all_cells[i])]
  which_normal <- net_cell_mapping$normal[which(net_cell_mapping$cell == which_tissue)]
  ne <- gtex %>% select(2,which_normal) %>% 
    filter(Description %in% all_genes) %>% 
    rename(gene=Description)
  colnames(ne)[2] <- "exp"
  ne <- ne %>% 
    group_by(gene) %>% 
    summarise(exp=mean(exp)) %>% ungroup() %>% as.data.frame()
  rownames(ne) <- ne$gene
  ne <- ne[all_genes,]
  colnames(ne)[2] <- all_cells[i]
  res[[i]] <- ne %>% select(-gene)
}
res <- bind_cols(res)

tumor_exp <- gene_exp %>% t() %>% as.data.frame()
tumor_exp <- tumor_exp[all_genes,all_cells] %>% as.matrix()
normal_exp <- as.matrix(res)
saveRDS(normal_exp,file = "/home/data/sdb/wt/model_data/gtex_normal_exp.rds")

diff <- tumor_exp / (log2(normal_exp + 1.01))
diff <- diff %>% t() %>% as.data.frame()
diff$cell <- rownames(diff)
diff <- diff %>% select(cell,everything())

sd_tumor <- apply(gene_exp,2,sd) %>% as.data.frame()
colnames(sd_tumor) <- "tumor_sd"
sd_tumor$gene <- rownames(sd_tumor)

gtex1 <- gtex[!duplicated(gtex$Description),] %>% select(-Name)
rownames(gtex1) <- gtex1$Description
gtex1 <- gtex1 %>% select(-Description)
sd_normal <- apply(gtex1,1,function(x){sd(log2(x+1))}) %>% as.data.frame()
colnames(sd_normal) <- "normal_sd"
sd_normal$gene <- rownames(sd_normal)

sd_all <- inner_join(sd_tumor,sd_normal) %>% 
  filter(tumor_sd >1 | normal_sd>1)

diff_exp <- diff %>% select(cell, any_of(sd_all$gene))
write.csv(diff_exp,file = "/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv",
          quote = F,row.names = F)

###cancer expression
gene_exp$cell <- rownames(gene_exp)
gene_exp <- gene_exp %>% select(cell,everything()) %>% 
  select(colnames(diff_exp)) %>% 
  filter(cell %in% diff_exp$cell)
write.csv(gene_exp,file = "/home/data/sdb/wt/model_data/cell_gene_exp_cancer.csv",
          quote = F,row.names = F)

####
nc_dt <- readxl::read_xlsx("data/nc_data1.xlsx",skip = 1)
nc_drug <- readxl::read_xlsx("data/nc_data2.xlsx",sheet = "AUC")
nc_cells <- c("ACH-000045","ACH-000362","ACH-000006","ACH-000301","ACH-000074",
              "ACH-000983","ACH-000168","ACH-002273","ACH-000004","ACH-000002",
              "ACH-000386","ACH-000551","ACH-000432","ACH-001647","ACH-000146")
colnames(nc_drug) <- c("drug",nc_cells)

drug_gene <- nc_dt %>%
  mutate(drug_name=`Compound name`) %>%
  mutate(target_gene = gsub("^[^\\|]*\\|","",Target)) %>%
  rowwise() %>%
  mutate(target_gene = paste(strsplit(target_gene,"\\|")[[1]] %>%
                               gsub(" ","",.),collapse = ",")) %>%
  select(drug_name,target_gene) %>%
  mutate(drug_name = stringr::str_trim(drug_name))
nc_drug <- nc_drug %>% na.omit()
nc_drug$drug[which(nc_drug$drug=="2,3',4,5'-Tetramethoxystilbene")] <- "Tetramethoxystilbene"
nc_drug$drug[which(nc_drug$drug=="AR-A-014418")] <- "AR-A 014418"
nc_drug$drug[which(nc_drug$drug=="Arachidonyl trifluoromethyl ketone")] <- "Arachidonyl Trifluoromethyl Ketone"
nc_drug$drug[which(nc_drug$drug=="Dephenhydramine")] <- "Diphenhydramine"
nc_drug$drug[which(nc_drug$drug=="LY 2183240")] <- "LY2183240"
nc_drug$drug[which(nc_drug$drug=="MKT 077")] <- "MKT-077"
nc_drug$drug[which(nc_drug$drug=="PD-151746")] <- "PD 151746"
nc_drug$drug[which(nc_drug$drug=="PD-146176")] <- "PD 146176"
nc_drug$drug[which(nc_drug$drug=="PF-543 hydrochloride")] <- "PF 543 hydrochloride"

nc_drug <- nc_drug %>%
  rowwise() %>%
  mutate(drug_name = drug_gene$drug_name[grep(drug,drug_gene$drug_name,fixed = T)[1]])

nc_drug <- left_join(
  nc_drug,drug_gene
)
saveRDS(nc_drug,file = "data/nc_drug.rds")

###kegg pathway
library(KEGGREST)
kegg_pathway <- read.table("data/rest.kegg.jp_list_pathway_hsa.txt",sep = "\t")
hsa <- lapply(kegg_pathway$V1, keggGet)  

res <- lapply(hsa,
              function(x){
                gene <- x[[1]]$GENE
                pathway_name <- x[[1]]$NAME
                pathway_class <- x[[1]]$CLASS
                if (length(gene) != 0){
                  all_genes <- stringr::str_split(gene[seq(2,length(gene),by=2)],';',
                                                  simplify = T)[,1]
                  dt <- data.frame(
                    pathway = pathway_name,
                    genes = all_genes,
                    class = pathway_class
                  )
                  return(dt)
                }else{
                  return(NA)
                }
              })
res <- res[which(lengths(res)>1)]
res <- bind_rows(res)
saveRDS(res,file = "data/kegg_all_pathway.rds")

###save mut data
mut <- data.table::fread("/home/data/sdc/wt/TCGA/mc3.v0.2.8.PUBLIC.xena.gz",
                         data.table = F)
mut$cancer <- EasyBioinfo::get_cancer_type(mut$sample,parallel = T,cores = 40)
saveRDS(mut,file = "/home/data/sdc/wt/TCGA/tcga_mut.rds")

###GSM 网络
library(Met2Graph)
infiles <- "data/GSM_xml/HMRdatabase.xml"
outDir <- "data/meta_net/EnzGraphs/"
Met2EnzGraph(infiles, rmMets=TRUE, outDir=outDir, 
             outFormat="ncol",add_rev_rxn = TRUE)

met_g <- read.table("data/meta_net/EnzGraphs/EnzGraphs/HMRdatabase_enzymes_based_graph.tsv")
enz_gene_mapping <- readRDS("data/enz_gene_mapping.rds")
ensg2name <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>%
    unique()
  sym <- tt_gene %>%
    paste(.,collapse=",")
  return(sym)
}
ids2gene <- data.frame(ensembl_id = unique(c(met_g$from,met_g$to))) %>%
  rowwise() %>%
  mutate(gene = ensg2name(ensembl_id,enz_gene_mapping)) %>% ungroup()

met_g <- left_join(
  met_g,
  ids2gene %>% rename(from = ensembl_id) %>% rename(from_gene = gene)
) %>% left_join(.,ids2gene %>% rename(to = ensembl_id) %>% 
                  rename(to_gene = gene))
saveRDS(met_g,file = "data/HGM_all_gene.rds")


