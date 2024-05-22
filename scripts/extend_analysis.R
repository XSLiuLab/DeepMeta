library(dplyr)
###培养基差异
cell_info <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
model_con <- data.table::fread("data/ModelCondition.csv",data.table = F)
model_con <- model_con %>% 
  filter(nchar(FormulationID)>0) %>% 
  mutate(media = gsub("[ + ].+","",FormulationID)) %>% 
  filter(media %in% c("DMEM","RPMI"))
model_con <- left_join(model_con %>% select(ModelID,media),
                       cell_info %>% select(ModelID,OncotreeLineage)) %>% 
  na.omit()

##选择样本量相对较多的 CNS/Brain 细胞系进行差异分析
cns <- model_con %>% filter(OncotreeLineage == "CNS/Brain")

gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionGenesExpectedCountProfile.csv",
                              data.table = F)
profiles <- data.table::fread("/home/data/sdb/wt/model_data/OmicsProfiles.csv",
                              data.table = F)
gene_exp <- inner_join(
  gene_exp,
  profiles %>% select(ProfileID,ModelID) %>% rename(V1=ProfileID)
) %>% select(ModelID,everything()) %>% select(-V1) %>% 
  distinct(ModelID,.keep_all = T)

rownames(gene_exp) <- gene_exp$ModelID
gene_exp <- gene_exp %>% select(-ModelID)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))
gene_exp <- gene_exp[intersect(rownames(gene_exp),cns$ModelID),] %>% 
  t()
gene_exp <- round(gene_exp) %>% as.data.frame()

model_con <- model_con %>% filter(ModelID %in% colnames(gene_exp))
DMEM <- model_con$ModelID[which(model_con$media == "DMEM")]
RPMI <- model_con$ModelID[which(model_con$media == "RPMI")]

gene_exp <- gene_exp[,c(DMEM,RPMI)]
deg_res <- EasyBioinfo::deg_deseq2(gene_exp,
                                   control_label = "DMEM",
                                   control_counts = length(DMEM),
                                   treatment_lable = "RPMI",
                                   treatment_counts = length(RPMI),
                                   parallel = TRUE,
                                   ncores = 40)
deg_res$genes <- rownames(deg_res)

EnhancedVolcano::EnhancedVolcano(deg_res,x="log2FoldChange",
                                 lab = rownames(deg_res),
                                 y="pvalue",pCutoffCol="padj",
                                 pCutoff = 0.1)

###建模去预测培养基类型
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))

diff_exp <- data.table::fread("/home/data/sdb/wt/model_data/cell_gene_exp_vs_normal_filter.csv",
                              data.table = F)
gene_exp <- gene_exp[intersect(rownames(gene_exp),model_con$ModelID),
                     intersect(colnames(gene_exp),colnames(diff_exp))]
gene_exp$cell <- rownames(gene_exp)
gene_exp <- gene_exp %>% select(cell,everything())
gene_exp <- left_join(
  gene_exp,
  model_con %>% select(ModelID,media) %>% rename(cell=ModelID)
) %>% select(cell,media,everything()) %>% 
  mutate(media = ifelse(media == "DMEM",0,1))
write.csv(gene_exp,file = "/home/data/sdb/wt/model_data/cell_media_pre.csv",
          quote = F,row.names = F)

model_con <- model_con %>% filter(ModelID %in% gene_exp$cell)
###res
library(yardstick)
library(dplyr)
library(ggplot2)
library(ggprism)
library(patchwork)
library(ggpubr)

pre <- lapply(list.files("/home/data/sdb/wt/model_data/",pattern = "cv",
                         full.names = F),function(x){
                           dt <- read.csv(paste0("/home/data/sdb/wt/model_data/",x)) %>% select(-X)
                           dt <- dt %>% 
                             mutate(truth = ifelse(label == 1, "Class1","Class2"),
                                    pred_label = ifelse(preds == 1, "Class1","Class2"))
                           dt$truth <- factor(dt$truth)
                           dt$pred_label <- factor(dt$pred_label)
                           dt$fold <- gsub("_media.csv","",x) %>% gsub("cv_","",.)
                           return(dt)
                         }) %>% bind_rows()

pre_met <- sapply(unique(pre$fold),
                  function(x){
                    dt <- pre %>% filter(fold == x)
                    pr <-  pr_auc(dt, truth, preds_raw)[".estimate"][[1]]
                    roc <-  roc_auc(dt, truth, preds_raw)[".estimate"][[1]]
                    return(data.frame(
                      pr = pr,
                      roc = roc,
                      fold = x
                    ))
                  },simplify = F) %>% bind_rows()

pre <- left_join(pre,pre_met)
pre <- pre %>% 
  mutate(fold_auroc = case_when(
    fold == "fold_0" ~ paste0("Fold 1: AUROC = ",round(roc,2)),
    fold == "fold_1" ~ paste0("Fold 2: AUROC = ",round(roc,2)),
    fold == "fold_2" ~ paste0("Fold 3: AUROC = ",round(roc,2)),
    fold == "fold_3" ~ paste0("Fold 4: AUROC = ",round(roc,2)),
    fold == "fold_4" ~ paste0("Fold 5: AUROC = ",round(roc,2)),
  )) %>% 
  mutate(fold_auprc = case_when(
    fold == "fold_0" ~ paste0("Fold 1: AUPRC = ",round(pr,2)),
    fold == "fold_1" ~ paste0("Fold 2: AUPRC = ",round(pr,2)),
    fold == "fold_2" ~ paste0("Fold 3: AUPRC = ",round(pr,2)),
    fold == "fold_3" ~ paste0("Fold 4: AUPRC = ",round(pr,2)),
    fold == "fold_4" ~ paste0("Fold 5: AUPRC = ",round(pr,2)),
  ))

roc_res <- vector("list",5)
for (i in 1:5){
  dt <- pre %>% filter(fold == paste0("fold_",i-1))
  dt_res <- roc_curve(dt, truth, preds_raw)
  dt_res$fold_auroc <- unique(dt$fold_auroc)
  roc_res[[i]] <- dt_res
}
roc_res <- bind_rows(roc_res)

ggplot(roc_res,aes(x = 1 - specificity, y = sensitivity, color = fold_auroc)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_prism()+
  theme(legend.position = c(0.7, 0.3))

pr_res <- vector("list",5)
for (i in 1:5){
  dt <- pre %>% filter(fold == paste0("fold_",i-1))
  dt_res <- pr_curve(dt, truth, preds_raw)
  dt_res$fold_auprc <- unique(dt$fold_auprc)
  pr_res[[i]] <- dt_res
}
pr_res <- bind_rows(pr_res)

ggplot(pr_res,aes(x = recall, y = precision, color = fold_auprc)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_prism()+
  theme(legend.position = c(0.7, 0.3))

################### TCGA 预测的依赖性值和细胞系预测的值的相关性分析
cell_pre <- data.table::fread("data/train_preV2.csv",data.table = F)
tcga_pre <- readRDS("data/tcga_pre_V2.rds") %>% as.data.frame()

###
###计算细胞系基因预测值之间的相关性
kegg <- readRDS("data/kegg_all_pathway.rds")
kegg <- kegg %>% 
  filter(grepl("Metabolism",class) | grepl("metabolism",pathway)) %>% 
  mutate(pathway = gsub(" \\- Homo sapiens \\(human\\)","",pathway))
get_pathway <- function(gene,pathway){
  gene <- strsplit(gene,split = ",")[[1]]
  pathway <- pathway %>% filter(genes %in% gene)
  return(c(paste(unique(pathway$pathway),collapse = ";"),
           paste(unique(pathway$class),collapse = ";")))
}
pre_gene <-  data.frame(id = unique(c(cell_pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]]) %>% ungroup()
pre <- cell_pre %>% 
  left_join(.,pre_gene %>% rename(gene_name=id))
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()
all_gene <- unique(pre_pathway$gene)

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
  i = 1:length(all_gene),
  .export = c("pre_pathway","all_gene"),
  .packages = c("dplyr")
) %dopar% {
  gene_dt <- pre_pathway %>% 
    filter(gene == all_gene[i]) %>% 
    select(cell,preds_raw) %>% 
    rename(preds2=preds_raw)
  
  gene_res <- sapply(all_gene,
                     function(x){
                       dt <- pre_pathway %>% filter(gene == x)
                       dt <- inner_join(
                         dt %>% select(cell,preds_raw),
                         gene_dt, by = "cell"
                       ) 
                       dt_res <- tryCatch({cor.test(dt$preds_raw,dt$preds2,
                                                    method = "sp")},
                                          error = function(e){
                                            return(data.frame(tt="error"))
                                          })
                       if (class(dt_res) == "data.frame"){
                         return(NA)
                       }else{
                         return(dt_res$estimate)
                       }
                     })
  gene_res <- as.data.frame(gene_res)
  gene_res$genes <- rownames(gene_res)
  colnames(gene_res)[1] <- "cor"
  gene_res$genes <- gsub(".rho","",gene_res$genes)
  gene_res$main_gene <- all_gene[i]
  return(gene_res)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
res_dt <- res %>% 
  tidyr::pivot_wider(names_from = main_gene,
                     values_from = cor) %>% as.data.frame()
rownames(res_dt) <- res_dt$genes
res_dt <- res_dt %>% select(-genes)
saveRDS(res_dt,file = "data/cell_pre_cor_all_gene_rho.rds")

#############
tcga_pre_cor <- readRDS("~/DeepMeta/data/tcga_pre_cor_all_gene_rho.rds")
cell_pre_cor <- readRDS("~/DeepMeta/data/cell_pre_cor_all_gene_rho.rds")
all_gene <- intersect(rownames(tcga_pre_cor),rownames(cell_pre_cor))
tcga_pre_cor <- tcga_pre_cor[all_gene,all_gene]
cell_pre_cor <- cell_pre_cor[all_gene,all_gene]

dt <- data.frame(genes=colnames(tcga_pre_cor),cor = NA)
for(i in 1:nrow(dt)){
  dt$cor[i] <- tryCatch({cor.test(tcga_pre_cor[,dt$genes[i]],
                                  cell_pre_cor[,dt$genes[i]],
                                  method = "sp")$estimate[[1]]},
                        error = function(e){
                          return(NA)
                        })
}
dt <- na.omit(dt)
library(ggpubr)
gghistogram(dt,x="cor",add = "median",fill = "#00AFBB",
            xlab = "Correlation between TCGA and Cell", ylab = "Counts",
            rug = TRUE,add.params = list(color = "red",size=1,linetype="dashed"))
ggsave("report/cor_bet_TCGA_cell.pdf",height = 5,width = 7)

library(ComplexHeatmap)
p1 <- Heatmap(cell_pre_cor,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,
        column_title  = "Correlation of predicted dependency in Cells",
        name = "Cells")
p2 <- Heatmap(tcga_pre_cor,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,
        column_title  = "Correlation of predicted dependency in TCGA",
        name = "TCGA")
p3 <- p1 + p2
pdf(file="report/cor_heatmap.pdf",width = 10,height = 5)
draw(p3)
dev.off()

###回归分析
dt <- readRDS("~/DeepMeta/data/TCGA_normal_dep_diff.rds")
pur <- dt %>% filter(pathway=="Pyrimidine metabolism")
pro <- data.table::fread("data/tcga_LM_propotion.csv",data.table = F)
all_cells <- colnames(pro)[2:23]
res <- vector("list",length(all_cells))
for (i in seq_along(res)){
  cell_pro <- pro %>% select(V1,all_cells[i]) %>% 
    rename(sample = V1) %>% 
    mutate(sample = gsub("[.]","-",sample))
  pur_cell <- left_join(pur,cell_pro)
  colnames(pur_cell)[8] <- "cell"
  fit <- lm(log2ratio ~ Type + cell, data = pur_cell)
  fit_res <- broom::tidy(fit,conf.int=TRUE)
  cell_res <- data.frame(
    cell_type = gsub("[.]","_",all_cells[i]),
    estimate = fit_res$estimate[2],
    low = fit_res$conf.low[2],
    high = fit_res$conf.high[2],
    stat = fit_res$statistic[2],
    p_value = fit_res$p.value[2]
  )
  res[[i]] <- cell_res
}
res <- bind_rows(res)

cell_pro <- pro %>% select(V1,2:23) %>% 
  rename(sample = V1) %>% 
  mutate(sample = gsub("[.]","-",sample))
pur_cell <- left_join(pur,cell_pro)
pur_cell <- pur_cell %>% select(5:29) %>% select(-cancer_type)
fit <- lm(log2ratio ~ . ,data = pur_cell)
fit_res <- broom::tidy(fit,conf.int=TRUE)
cell_res <- data.frame(
  cell_type = "All",
  estimate = fit_res$estimate[2],
  low = fit_res$conf.low[2],
  high = fit_res$conf.high[2],
  stat = fit_res$statistic[2],
  p_value = fit_res$p.value[2]
)
res <- bind_rows(res,cell_res)
res$` ` <- paste(rep(" ", 20), collapse = " ")
res$p_value <- format(res$p_value,digits = 3)
res <- res %>% rename(`Controled Cell Type`=cell_type,
                      `P value`=p_value)
library(forestploter)
p <- forest(data = res[,c(1,7,6)],
            est = res$estimate,
            lower = res$low,
            upper = res$high,
            ref_line = 0,
            ci_column = 2)
ggsave("report/Pyri_controled_cell.pdf",p,width = 8,height = 10)

####更多临床数据的验证
dt <- lapply(1:23,
             function(x){
               readxl::read_xlsx("data/chemotherapy_dataset.xlsx",sheet = x)
             })
dt <- bind_rows(dt)
dt <- dt %>% 
  filter(grepl("fluorouracil|capecitabine",Treatment))

all_files <- list.files("/home/data/sdb/wt/model_data/chemo/",pattern = "CTR",
                        full.names = T)
res <- vector("list",length(all_files))
for (i in seq_along(res)){
  dt <- data.table::fread(paste0(all_files[i],"/cli.inf_.csv"),data.table = F)
  dt <- dt %>% select(Resource,Source,Sample_id,Data_type,Cancer_type_level1,
                      Original_response,Response_standard,Response)
  dt$dataset <- gsub("/home/data/sdb/wt/model_data/chemo//","",all_files[i])
  if (dt$Resource[1] == "TCGA"){
    res[[i]] <- NA
  }else{
    res[[i]] <- dt
  }
}
res <- res[which(lengths(res) > 1)]
res <- bind_rows(res)

need_files <- unique(res$dataset)
dt_exp <- data.table::fread(paste0(all_files[2],"/matrix_.csv"),data.table = F)

