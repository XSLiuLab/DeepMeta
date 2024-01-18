library(dplyr)
###看看肿瘤增殖是否和嘧啶代谢依赖性相关
res <- readRDS("data/pancancer_meta_dep.rds")
py <- res %>% filter(pathway == "Pyrimidine metabolism")

tcga_drug <- readRDS("~/meta_target/data/tcga_drug.rds")
tcga_drug <- tcga_drug %>% 
  mutate(response_type = 
           ifelse(response %in% c("Complete Response","Partial Response"),"response",
                  "no-response"),
         drug_type = ifelse(drug.name %in% 
                              c("Capecitabine","Pemetrexed","Gemcitabine","Fluorouracil"),
                            "anti_nm","no_anti_nm")) %>% 
  filter(patient.arr %in% substr(py$sample,1,12)) %>% 
  select(patient.arr,drug_type) %>% 
  distinct_all()
anti_nm <- tcga_drug %>% 
  group_by(patient.arr) %>% 
  summarise(drug_type = ifelse("anti_nm" %in% drug_type,"anti_nm","no_anti_nm")) %>% 
  ungroup() %>% 
  rename(sample=patient.arr) %>% 
  filter(drug_type == "anti_nm")

###没有用抗嘧啶代谢药物的病人
py_no_nm <- py %>% filter(!(substr(sample,1,12) %in% anti_nm$sample))

###
CCS_genes <- readRDS("~/meta_target/data/CCS_genes.rds")
###
tumor_tpm <- data.table::fread("/home/data/sdc/wt/TCGA/tcga_RSEM_gene_tpm.gz",data.table = F)
tumor_tpm <- tumor_tpm %>% 
  select(sample,which(as.numeric(substr(colnames(tumor_tpm),14,15)) < 10)) %>% 
  rename(id = sample)
mapping <- data.table::fread("/home/data/sdc/wt/TCGA/probeMap_gencode.v23.annotation.gene.probemap",
                             data.table = F)
tumor_tpm <- left_join(tumor_tpm,mapping %>% select(id,gene)) %>% 
  select(-id) %>% 
  select(gene,everything())
tumor_tpm <- tumor_tpm[!duplicated(tumor_tpm$gene),]
tumor_tpm <- tumor_tpm %>% filter(gene %in% CCS_genes$V1)
tumor_ccs <- apply(tumor_tpm[,2:ncol(tumor_tpm)],2,sum) %>% as.data.frame()
colnames(tumor_ccs) <- "ccs"
tumor_ccs$sample <- rownames(tumor_ccs)

normal_tpm <- data.table::fread("/home/data/sdb/TCGA_data/gtex_RSEM_gene_tpm.gz",
                                data.table = F)
normal_tpm <- left_join(normal_tpm %>% rename(id=sample),
                        mapping %>% select(id,gene)) %>% 
  select(-id) %>% 
  select(gene,everything())
normal_tpm <- normal_tpm %>% filter(gene %in% CCS_genes$V1)
normal_tpm <- normal_tpm[!duplicated(normal_tpm$gene),]
normal_ccs <- apply(normal_tpm[,2:ncol(normal_tpm)],2,sum) %>% as.data.frame()
colnames(normal_ccs) <- "ccs"
normal_ccs$sample <- rownames(normal_ccs)

meta_info <- data.table::fread("/home/data/sdb/TCGA_data/GTEX_phenotype.gz",
                               data.table = F)
normal_ccs <- left_join(normal_ccs, meta_info %>% 
                          select(Sample,`_primary_site`) %>% 
                          rename(sample = Sample))

###注释癌症类型和组织类型
match_meta <- data.table::fread("/home/data/sdb/TCGA_data/TcgaTargetGTEX_phenotype.txt.gz",
                                data.table = F)
tumor_ccs$cancer <- EasyBioinfo::get_cancer_type(tumor_ccs$sample)
tumor_ccs <- left_join(
  tumor_ccs, 
  match_meta %>% 
    select(sample, `_primary_site`)
)

normal_ccs <- inner_join(
  normal_ccs,
  tumor_ccs %>% select(cancer, `_primary_site`) %>% distinct_all()
)
normal_ccs_summ <- normal_ccs %>% 
  group_by(cancer) %>% 
  summarise(median_ccs = median(ccs)) %>% ungroup()

tumor_ccs <- inner_join(tumor_ccs,normal_ccs_summ)

###add purity
purity <- data.table::fread("data/TCGA_mastercalls.abs_tables_JSedit.fixed.txt",
                            data.table = F)
tumor_ccs <- inner_join(
  tumor_ccs,
  purity %>% select(array,purity) %>% rename(sample=array)
)
tumor_ccs <- na.omit(tumor_ccs)

tumor_ccs <- tumor_ccs %>% 
  mutate(bc_ccs = (ccs - median_ccs) * purity) %>% 
  mutate(bc_ccs = scales::rescale(bc_ccs))

tumor_ccs <- inner_join(
  tumor_ccs,py
)

tumor_ccs_summ <- tumor_ccs %>% 
  group_by(cancer) %>% 
  summarise(median_ccs = median(ccs),
            median_bc = median(bc_ccs),
            median_or = median(ratio)) %>% ungroup()
###gsva
tumor_tpm <- data.table::fread("/home/data/sdc/wt/TCGA/tcga_RSEM_gene_tpm.gz",data.table = F)
tumor_tpm <- tumor_tpm %>% 
  select(sample,which(as.numeric(substr(colnames(tumor_tpm),14,15)) < 10)) %>% 
  rename(id = sample)
mapping <- data.table::fread("/home/data/sdc/wt/TCGA/probeMap_gencode.v23.annotation.gene.probemap",
                             data.table = F)
tumor_tpm <- left_join(tumor_tpm,mapping %>% select(id,gene)) %>% 
  select(-id) %>% 
  select(gene,everything())
tumor_tpm <- tumor_tpm[!duplicated(tumor_tpm$gene),]
rownames(tumor_tpm) <- tumor_tpm$gene
tumor_tpm$gene <- NULL

library(GSVA)
gs <- list(ccs = CCS_genes$V1)
gsvaPar <- gsvaParam(as.matrix(tumor_tpm), gs)
gsva.es <- gsva(gsvaPar, verbose=TRUE,
                BPPARAM = BiocParallel::SnowParam(workers = 20, type = "SOCK"))
gsva.es <- as.data.frame(t(gsva.es))
gsva.es$sample <- rownames(gsva.es)

res <- inner_join(
  py,
  gsva.es
)
saveRDS(res,file = "data/tcga_ccs.rds")

res <- readRDS("data/tcga_ccs.rds")
res <- res %>% 
  filter(sample %in% py_no_nm$sample)
res <- res %>% 
  mutate(type = ifelse(ratio > median(res$ratio),"high","low"))

library(ggplot2)
library(ggpubr)
ggboxplot(res,x="type",y="ccs",xlab = "OR_type",
          ylab = "Cell cycle signature score")+
  stat_compare_means()
ggsave("figs/ccs_compare_or.pdf",height = 4,width = 4)





