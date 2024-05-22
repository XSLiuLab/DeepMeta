library(dplyr)
library(ezcox)

kegg <- readRDS("~/meta_target/data/kegg_all_pathway.rds")
kegg <- kegg %>% 
  filter(grepl("Metabolism",class) | grepl("metabolism",pathway)) %>% 
  mutate(pathway = gsub(" \\- Homo sapiens \\(human\\)","",pathway))

get_pathway <- function(gene,pathway){
  gene <- strsplit(gene,split = ",")[[1]]
  pathway <- pathway %>% filter(genes %in% gene)
  return(c(paste(unique(pathway$pathway),collapse = ";"),
           paste(unique(pathway$class),collapse = ";")))
}

pre <- readRDS("data/tcga_pre_V2.rds")
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()

###
surv <- readRDS("data/pancancer_survial.rds")
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

tcga_drug <- tcga_drug %>% 
  group_by(patient.arr) %>% 
  summarise(drug_type = ifelse("anti_nm" %in% drug_type,"anti_nm","no_anti_nm")) %>% 
  ungroup()
  
###用了 anti-nm
surv_nm <- left_join(
  tcga_drug %>% 
    rename(sample=patient.arr) %>% 
    filter(drug_type == "anti_nm"),
  surv
) %>% left_join(.,
                py %>% mutate(sample = substr(sample,1,12))) %>% 
  distinct_all(.keep_all = TRUE)

surv_nm$cancers <- EasyBioinfo::get_cancer_type(surv_nm$sample)

cancer_summ <- surv_nm %>% 
  group_by(cancers) %>% 
  summarise(counts = n()) %>% ungroup() %>% 
  filter(counts > 10)
surv_nm_filter <- surv_nm %>% 
  filter(cancers %in% cancer_summ$cancers)
surv_nm_filter <- surv_nm_filter %>% 
  mutate(OR_type = ifelse(ratio > median(surv_nm_filter$ratio), "high","low")) 

surv_nm_filter$OR_type <- factor(surv_nm_filter$OR_type, levels = c("low","high"))

p1 <- show_forest(surv_nm_filter,covariates = "OR_type",time = "OS.time",
                  status = "OS",controls = "cancers",vars_to_show = "OR_type"
                  )
p2 <- EasyBioinfo::show_km(surv_nm_filter ,"OR_type",title="All Cancer Type")
p1
p2

dt <- surv_nm_filter %>% filter(!is.na(Tumor_stage))
surv_nm_filter <- surv_nm_filter %>% 
  rename(Dependency=OR_type)
p1 <- show_forest(surv_nm_filter,covariates = "Dependency",time = "OS.time",
                  status = "OS",controls = c("cancers","age","gender","Tumor_stage"),
                  vars_to_show = "Dependency")
##所有癌症类型
surv_list <- vector("list",7)
for (i in 1:7){
  surv_list[[i]] <- EasyBioinfo::show_km(surv_nm_filter %>% 
                                           filter(cancers == cancer_summ$cancers[i]),
                                         "OR_type",title=cancer_summ$cancers[i])
}
p7 <- survminer::arrange_ggsurvplots(surv_list,nrow = 2,ncol = 4)
pdf("figs/py_survial_all_cancer.pdf",height = 13,width = 22)
print(p7)
dev.off()

###没有用 anti-nm 药物的化疗病人
surv_no_nm <- left_join(
  tcga_drug %>% 
    rename(sample=patient.arr) %>% 
    filter(drug_type == "no_anti_nm") %>% 
    select(sample,drug_type),
  surv
) %>% left_join(.,
                py %>% mutate(sample = substr(sample,1,12))) %>% 
  distinct_all(.keep_all = TRUE)

surv_no_nm$cancers <- EasyBioinfo::get_cancer_type(surv_no_nm$sample)
cancer_summ <- surv_no_nm %>%
  group_by(cancers) %>%
  summarise(counts = n()) %>% ungroup() %>%
  filter(counts > 10)

surv_no_nm_filter <- surv_no_nm %>% 
  filter(cancers %in% cancer_summ$cancers)
surv_no_nm_filter <- surv_no_nm_filter %>% 
  mutate(OR_type = ifelse(ratio > median(surv_no_nm_filter$ratio), 
                          "high","low")) 

surv_no_nm_filter$OR_type <- factor(surv_no_nm_filter$OR_type, 
                                    levels = c("low","high"))

p3 <- show_forest(surv_no_nm_filter,covariates = "OR_type",time = "OS.time",
                  status = "OS",controls = "cancers",vars_to_show = "OR_type"
)
p4 <- EasyBioinfo::show_km(surv_no_nm_filter ,"OR_type",title="All Cancer Type")

p3
p4

p6 <- survminer::arrange_ggsurvplots(list(p2,p4))
pdf("figs/py_survial_with_control.pdf",height = 6,width = 14)
print(p6)
dev.off()

library(patchwork)
library(ggplot2)
p1 + p3
ggsave("figs/tcga_py_vs_nopy.pdf",height = 4,width = 12)

###
###所有病人
surv_all <- inner_join(
  surv,
  py %>% mutate(sample = substr(sample,1,12))
) %>% distinct_all(.keep_all = TRUE) %>% filter(!is.na(OS.time)) 

##所有病人里面去掉用了抗嘧啶药物的病人
surv_all <- surv_all %>% filter(!(sample %in% surv_nm$sample))

surv_all$cancers <- EasyBioinfo::get_cancer_type(surv_all$sample)
cancer_summ <- surv_all %>%
  group_by(cancers) %>%
  summarise(counts = n()) %>% ungroup() %>%
  filter(counts > 10)
surv_all_filter <- surv_all %>% 
  filter(cancers %in% cancer_summ$cancers)

surv_all_filter <- surv_all_filter %>% 
  mutate(OR_type = ifelse(ratio > median(surv_all_filter$ratio), "high","low"))

surv_all_filter$OR_type <- factor(surv_all_filter$OR_type, 
                                  levels = c("low","high"))

p5 <- show_forest(surv_all_filter,covariates = "OR_type",time = "OS.time",
                  status = "OS",controls = "cancers",vars_to_show = "OR_type")
p5
ggsave("figs/surv_drug_all.pdf",height = 4,width = 6)

p6 <- EasyBioinfo::show_km(surv_all_filter ,"OR_type",title="All Cancer Type")

p7 <- survminer::arrange_ggsurvplots(list(p2,p6))
pdf("figs/py_survial_with_control2.pdf",height = 6,width = 14)
print(p7)
dev.off()

dt <- surv_all_filter %>% filter(!is.na(Tumor_stage))
surv_all_filter <- surv_all_filter %>% 
  rename(Dependency=OR_type)
p2 <- show_forest(dt,covariates = "Dependency",time = "OS.time",
            status = "OS",controls = c("cancers","age","gender","Tumor_stage"),
            vars_to_show = "Dependency")

ggsave(p1,filename = "report/drug_tcga_controled_stage.pdf",width = 7,height = 4)
ggsave(p2,filename = "report/no_drug_tcga_controled_stage.pdf",width = 7,height = 4)

