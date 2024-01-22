library(dplyr)
####NC 文章药物的数据
cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_info <- cell_info %>%
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))

gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

nc_cells <- c("ACH-000045","ACH-000362","ACH-000006","ACH-000301","ACH-000074",
              "ACH-000983","ACH-000168","ACH-002273","ACH-000004","ACH-000002",
              "ACH-000386","ACH-000551","ACH-000432","ACH-001647","ACH-000146")
gene_exp <- gene_exp[intersect(rownames(gene_exp),nc_cells),] %>% 
  t() %>% as.data.frame()
gene_exp$gene_name <- rownames(gene_exp)
gene_exp <- gene_exp %>% select(gene_name,everything())

library(Prep4DeepDEP)

train_dtV2 <- readRDS("~/meta_target/data/train_dtV2.rds")
pred_genes <- paste(train_dtV2$gene,collapse = ",") %>% 
  strsplit(.,",") %>% `[[`(1) %>% unique()
Prep4DeepDEP(exp.data = gene_exp,
             mut.data = NULL,
             meth.data = NULL,
             cna.data = NULL,
             dep.data = data.frame(genes = pred_genes),
             mode = "prediction",
             filename.out = "/home/data/sdc/wt/DeepDEP/drug2")
###
deepdep_drug <- data.table::fread("/home/data/sdc/wt/DeepDEP/drug2_pred.txt",
                                  data.table = F)
deepdep_drug <- deepdep_drug %>% 
  tidyr::pivot_longer(cols = 2:15, names_to = "cell",
                      values_to = "preds")

nc_drug <- readRDS("data/nc_drug.rds")
get_auc <- function(cell,gene,drug_dt){
  genes <- strsplit(gene,split = ",")[[1]]
  res <- vector("list",length(genes))
  for (i in 1:length(genes)){
    drug_dt_t <- drug_dt %>% filter(grepl(genes[i],target_gene))
    res[[i]] <- drug_dt_t
  }
  res <- bind_rows(res)
  if (nrow(res)==0){
    return(list(NA,NA))
  }else{
    res <- res %>% select(cell,drug_name)
    return(list(max(res[,1]),
                paste(res$drug_name[unlist(res[,1]==max(res[,1]))],collapse = ",")))
  }
}

deepdep_pre_drug <- deepdep_drug %>% 
  rowwise() %>% 
  mutate(max_auc = get_auc(cell,CRISPR_GENE,nc_drug)[[1]],
         max_auc_drug = get_auc(cell,CRISPR_GENE,nc_drug)[[2]]) %>% 
  ungroup()
deepdep_pre_drug <- deepdep_pre_drug %>% na.omit()
deepdep_pre_drug <- deepdep_pre_drug %>% 
  group_by(cell) %>% 
  mutate(rank_auc = rank(-max_auc)) %>% 
  ungroup()
saveRDS(deepdep_pre_drug,file = "data/deepdep_pre_drug.rds")

deepdep_pre_drug <- readRDS("data/deepdep_pre_drug.rds")
deepdep_drug_top <- deepdep_pre_drug %>%
  group_by(cell) %>% 
  slice_min(order_by = preds, n =5) %>% 
  ungroup()
###label rank = 1
library(ggpubr)
library(ggrepel)
deepdep_drug_top <- deepdep_drug_top %>% 
  mutate(ll = ifelse(rank_auc <= 3,
                     paste0("Rank: ",rank_auc,", AUC: ",round(max_auc,2)),
                     ""))
p1 <- ggscatter(deepdep_drug_top, x = "cell", y = "preds")+
  geom_text_repel(aes(label = ll),box.padding = 0.5,min.segment.length = 0) +
  geom_point(color = ifelse(deepdep_drug_top$ll == "", "grey50", "red"),size=4)+
  rotate_x_text(45)+
  labs(x="Cell Lines",y="Model Prediction")
###
pre_drug <- readRDS("data/nc_pre_drug_V2.rds")
pre_drug_top <- pre_drug %>%
  group_by(cell) %>% 
  slice_max(order_by = preds_raw, n =5) %>% 
  ungroup()
pre_drug_top <- pre_drug_top %>% 
  mutate(ll = ifelse(rank_auc <= 3,
                     paste0("Rank: ",rank_auc,", AUC: ",round(max_auc,2)),
                     ""))
p2 <- ggscatter(pre_drug_top, x = "cell", y = "preds_raw")+
  geom_text_repel(aes(label = ll),box.padding = 0.5,min.segment.length = 0) +
  geom_point(color = ifelse(pre_drug_top$ll == "", "grey50", "red"),size=4)+
  rotate_x_text(45)+
  labs(x="Cell Lines",y="Model Prediction")

library(patchwork)
p3 <- p2 / p1
ggsave("figs/nc_drug_val_V2_compare.pdf",p3,width = 15,height = 15)

dt <- data.frame(
  aucs = c(pre_drug_top$max_auc[which(nchar(pre_drug_top$ll)>0)],
           deepdep_drug_top$max_auc[which(nchar(deepdep_drug_top$ll)>0)]),
  types = c(rep("DeepMeta",
                length(pre_drug_top$max_auc[which(nchar(pre_drug_top$ll)>0)])),
            rep("DeepDep",
                length(deepdep_drug_top$max_auc[which(nchar(deepdep_drug_top$ll)>0)])))
)
ggboxplot(dt,x="types",y="aucs",xlab = FALSE,ylab = "AUC",
          order = c("DeepDep","DeepMeta"))+
  stat_compare_means()
ggsave("figs/nc_drug_val_V2_compare_box.pdf",width = 5,height = 5)

####CCMA
ccma_exp <- readRDS("/home/data/sdc/wt/model_data/ccma_all_exp.rds") %>% 
  as.data.frame()
rownames(ccma_exp) <- ccma_exp$sample
ccma_exp <- ccma_exp %>% select(-sample)
ccma_exp <- as.data.frame(t(ccma_exp))
ccma_exp$gene_name <- rownames(ccma_exp)
ccma_exp <- ccma_exp %>% select(gene_name,everything())

library(Prep4DeepDEP)

ccma_dep <- data.table::fread("/home/data/sdc/wt/model_data/CCMA_CRISPRdependencyscore.csv",
                              data.table = F)
pre <- read.csv("data/ccma_preV2.csv") %>% 
  select(-X)
pre_gene <-  data.frame(id = unique(c(pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]]) ##load from other scripts
pred_genes <- paste(pre_gene$gene,collapse = ",") %>% 
  strsplit(.,",") %>% `[[`(1) %>% unique()

Prep4DeepDEP(exp.data = ccma_exp,
             mut.data = NULL,
             meth.data = NULL,
             cna.data = NULL,
             dep.data = data.frame(genes = pred_genes),
             mode = "prediction",
             filename.out = "/home/data/sdc/wt/DeepDEP/CCMA2")
###预测值
ccma_dep <- data.table::fread("/home/data/sdc/wt/model_data/CCMA_CRISPRdependencyscore.csv",
                              data.table = F)
ccma_pre <- data.table::fread("/home/data/sdc/wt/DeepDEP/CCMA2_pred.txt",
                              data.table = F)
ccma_pre <- ccma_pre %>% 
  tidyr::pivot_longer(cols = 2:191, names_to = "sample",
                      values_to = "preds") %>% 
  rename(gene = CRISPR_GENE)

ccma_pre_dep <- left_join(
  ccma_pre,
  ccma_dep
) %>% filter(!is.na(z)) %>% 
  group_by(sample,gene) %>% 
  summarise(preds = unique(preds),
            z_score = min(z,na.rm = T)) %>% 
  ungroup() %>% 
  mutate(label = ifelse(z_score < (-1.5),1,0)) %>% 
  mutate(preds_label = ifelse(preds < (-0.5),1,0))


ccma_pre_dep <- ccma_pre_dep %>% 
  mutate(truth = ifelse(label == 1, "Class1","Class2"),
         pred_label = ifelse(preds_label == 1, "Class1","Class2"))
ccma_pre_dep$truth <- factor(ccma_pre_dep$truth)
ccma_pre_dep$pred_label <- factor(ccma_pre_dep$pred_label)
ccma_pre_dep$preds1 <- - ccma_pre_dep$preds 

library(yardstick)
library(ggplot2)
pr <-  pr_auc(ccma_pre_dep, truth, preds1)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)
roc <-  roc_auc(ccma_pre_dep, truth, preds1)[".estimate"] %>% 
  unlist() %>% unname() %>% round(.,2)

p1 <- pr_curve(ccma_pre_dep, truth, preds1) %>%
  ggplot(aes(x = recall, y = precision)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.5, y=0.85, label=paste0("PRC-AUC: ",pr),
           size=5)

p2 <- roc_curve(ccma_pre_dep, truth, preds1) %>%
  ggplot(aes(x = 1 - specificity, y = sensitivity)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_bw()+
  annotate(geom="text", x=0.4, y=0.85, label=paste0("ROC-AUC: ",roc),
           size=5)

library(ggpubr)
p3 <- ggboxplot(data=ccma_pre_dep,x="preds_label",y="z_score")+
  stat_compare_means()+
  labs(x="Model Prediction",y="Dependency Z scores")+
  geom_hline(yintercept = (-1.5),color="red",linewidth=2,linetype=2)

###
ccma_cor <- ccma_pre_dep %>% 
  group_by(sample) %>% 
  summarise(sample_cor = cor(preds,z_score)) %>% 
  ungroup()

p4 <- gghistogram(ccma_cor, x = "sample_cor", 
                  fill = "#00AFBB",
                  add = "median", rug = TRUE, 
                  title = "Correlation of predicted and true dependency")+
  labs(y="Sample Counts",x="Correlation")

library(patchwork)
(p1 + p2)
ggsave("figs/ccma_compare_roc_prc.pdf",width = 11,height = 5)
(p3 + p4)
ggsave("figs/ccma_compare_z_pre.pdf",width = 11,height = 7)



