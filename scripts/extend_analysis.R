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
                                 pCutoff = 0.1,
                                 title = "RPMI VS DMEM",subtitle = NULL)

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

p1 <- ggplot(roc_res,aes(x = 1 - specificity, y = sensitivity, color = fold_auroc)) +
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

p2 <- ggplot(pr_res,aes(x = recall, y = precision, color = fold_auprc)) +
  geom_path() +
  coord_fixed(xlim = 0:1, ylim = 0:1)+
  theme_prism()+
  theme(legend.position = c(0.7, 0.3))

library(patchwork)
p1+p2
ggsave("report/auc_cell_media.pdf",width = 10,height = 5)
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
p1 <- gghistogram(dt,x="cor",add = "median",fill = "#00AFBB",
            xlab = "Correlation between TCGA and Cell", ylab = "Counts",
            rug = TRUE,add.params = list(color = "red",size=1,linetype="dashed"))
ggsave("report/cor_bet_TCGA_cell.pdf",height = 5,width = 7)

library(ComplexHeatmap)

p2 <- Heatmap(cell_pre_cor,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,
        column_title  = "Correlation of predicted dependency in Cells",
        name = "Cells")
p3 <- Heatmap(tcga_pre_cor,cluster_rows = F,cluster_columns = F,
        show_row_names = F,show_column_names = F,
        column_title  = "Correlation of predicted dependency in TCGA",
        name = "TCGA")
p3 <- p2 + p3
pdf(file="report/cor_heatmap.pdf",width = 10,height = 5)
draw(p3)
dev.off()

p3 <- grid.grabExpr(draw(p3)) 
p1 + p3

###回归分析考虑细胞比例
get_controled_p <- function(pathway_name){
  dt <- readRDS("~/DeepMeta/data/TCGA_normal_dep_diff.rds")
  pur <- dt %>% filter(pathway==pathway_name)
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
  res$` ` <- paste(rep(" ", 25), collapse = " ")
  res$p_value <- format(res$p_value,digits = 3)
  res <- res %>% rename(`Controled Cell Type`=cell_type,
                        `P value`=p_value)
  colnames(res)[6] <- paste0(pathway_name," P value")
  return(res)
}

res_py <- get_controled_p("Pyrimidine metabolism")
res_glu <- get_controled_p("Glutathione metabolism")
res_pu <- get_controled_p("Purine metabolism")
library(forestploter)
p1 <- forest(data = res_py[,c(1,7,6)],
            est = res_py$estimate,
            lower = res_py$low,
            upper = res_py$high,
            ref_line = 0,
            ci_column = 2)
p2 <- forest(data = res_glu[,c(1,7,6)],
             est = res_glu$estimate,
             lower = res_glu$low,
             upper = res_glu$high,
             ref_line = 0,
             ci_column = 2)
p3 <- forest(data = res_pu[,c(1,7,6)],
             est = res_pu$estimate,
             lower = res_pu$low,
             upper = res_pu$high,
             ref_line = 0,
             ci_column = 2)
ggsave("report/Pyri_controled_cell.pdf",p1,width = 8,height = 10)
ggsave("report/Pu_controled_cell.pdf",p3,width = 8,height = 10)
ggsave("report/Glu_controled_cell.pdf",p2,width = 8,height = 10)
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
saveRDS(res,file = "data/geo_chemo_meta.rds")

res <- readRDS("data/geo_chemo_meta.rds")
need_files <- unique(res$dataset)
exp_res <- vector("list",length(need_files))
for (i in seq_along(exp_res)){
  dt_exp <- data.table::fread(paste0("/home/data/sdb/wt/model_data/chemo/",
                                     need_files[i],"/matrix_.csv"),data.table = F)
  exp_res[[i]] <- dt_exp
}
all_exp <- Reduce(function(x,y){
  merge(x, y, by='V1', all.x=TRUE, all.y=TRUE)
}, exp_res)
rownames(all_exp) <- all_exp$V1
##用最低来填补缺失的基因
for(i in 2:ncol(all_exp)){
  # all_exp[is.na(all_exp[,i]), i] <- mean(all_exp[,i], na.rm = TRUE)
  all_exp[is.na(all_exp[,i]), i] <- min(all_exp[,i], na.rm = TRUE)
}
saveRDS(all_exp,file = "/home/data/sdb/wt/model_data/geo_chemo_exp.rds")

library(PreDeepMeta)
data("model_gene_order")
all_exp <- all_exp %>% filter(V1 %in% model_gene_order)
all_exp2 <- all_exp[model_gene_order,]
rownames(all_exp2) <- model_gene_order
all_exp2$V1 <- NULL
for(i in 1:ncol(all_exp2)){
  all_exp2[is.na(all_exp2[,i]), i] <- min(all_exp2[,i], na.rm = TRUE)
}
########准备输入文件
##正常样本的中位数表达
##BC ：GSE15852
gemma <- readRDS("~/Immune_state/data/gemma_samples.rds")
bc <- gemma %>% filter(ids == "GSE15852" & grepl("Normal",samples))
dt <- readRDS("/home/data/sdb/wt/GEO_GEMMA/GSE15852.rds")
dt <- dt %>% 
  select(-c(1:3)) %>% 
  filter(nchar(NCBIid) > 0) %>% 
  tidyr::separate_longer_delim(cols = NCBIid, delim = "|")
dt <- dt %>% select_if(~ !any(is.na(.))) ##去掉 NA 的列
esm_id <- genekitr::transId(unique(dt$NCBIid), transTo = "symbol")
esm_id <- esm_id[!duplicated(esm_id$input_id),]
esm_id <- esm_id[!duplicated(esm_id$symbol),]
dt <- inner_join(
  dt,
  esm_id %>% rename(NCBIid = input_id)
) %>% select(symbol,everything()) %>% select(-NCBIid) %>% as.data.frame()
###相同基因取均值
dt <- dt %>% 
  group_by(symbol) %>% 
  summarise(across(2:(ncol(dt)-1),~mean(.x,na.rm = T))) %>%
  ungroup() %>% as.data.frame()
library(PreDeepMeta)
data("model_gene_order")
rownames(dt) <- dt$symbol
dt$symbol <- NULL
dt2 <- dt[model_gene_order,]
rownames(dt2) <- model_gene_order
for(i in 1:ncol(dt2)){
  dt2[is.na(dt2[,i]), i] <- min(dt2[,i], na.rm = TRUE)
}
bc_median <- apply(dt2,1,median) %>% as.data.frame()
colnames(bc_median) <- "Breast"
bc_median$Description <- rownames(bc_median)

###coad
dt <- data.table::fread("/home/data/sdb/wt/model_data/CLX_Expression_PC1_2014Nov24.txt",
                        data.table = F)
dt <- dt %>% select(V1,grep("_N",colnames(dt)))
rownames(dt) <- dt$V1
dt$V1 <- NULL
dt2 <- dt[model_gene_order,]
rownames(dt2) <- model_gene_order
for(i in 1:ncol(dt2)){
  dt2[is.na(dt2[,i]), i] <- min(dt2[,i], na.rm = TRUE)
}
crc_median <- apply(dt2,1,median) %>% as.data.frame()
colnames(crc_median) <- "Colon"
crc_median$Description <- rownames(crc_median)

normal_exp <- left_join(
  bc_median %>% select(Description,everything()),
  crc_median
)
saveRDS(normal_exp,file = "data/array_bc_crc_normal.rds")
##
sample_info <- readRDS("data/geo_chemo_meta.rds")
normal_exp <- readRDS("data/array_bc_crc_normal.rds")
sample_info <- sample_info %>% select(Sample_id,Cancer_type_level1) %>% 
  rename(ModelID=Sample_id) %>% 
  mutate(normal=ifelse(Cancer_type_level1 == "Breast cancer",
                       "Breast",
                       "Colon"))
library(PreDeepMeta)
data("model_gene_order")
diff_exp <- PreDiffExp(tumor_exp = all_exp2, normal_exp = normal_exp,
                  tumor_normal_mapping = sample_info,
                  gene_order = model_gene_order,
                  save_file = FALSE)
write.csv(diff_exp,
          file = "/home/data/sdb/wt/model_data/geo_chemo_diff_exp.csv",
          quote = F,row.names = F)
###
data("enz_gene_mapping")
data("cpg_gene")
brca <- read.table("data/meta_net/EnzGraphs/iBreastCancer1746_enzymes_based_graph.tsv")
coad <- read.table("data/meta_net/EnzGraphs/iColorectalCancer1750_enzymes_based_graph.tsv")
sample_exp <- readRDS("/home/data/sdb/wt/model_data/geo_chemo_exp.rds")

# low_exp <- apply(sample_exp[,2:ncol(sample_exp)],2,
#                  function(x){mean(x<3)}) %>% as.data.frame() %>%
#   rename(per=1)
# median(low_exp$per)

library(foreach)
library(parallel)
my.cluster <- parallel::makeCluster(
  80,
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
res <- foreach(
  i = 1:nrow(sample_info),
  .export = c("sample_exp","brca","cpg_gene","coad",
              "enz_gene_mapping","sample_info"),
  .packages = c("dplyr","PreDeepMeta")
) %dopar% {
  if (sample_info$Cancer_type_level1[i] == "Breast cancer"){
    PreEnzymeNet(sample_exp, 
                 network = brca,
                 gene_mapping = enz_gene_mapping, 
                 gene_feature = cpg_gene,
                 cell_name = sample_info$ModelID[i], 
                 exp_cutoff = 1,
                 save_path = "/home/data/sdb/wt/model_data/geo_chemo_net/")
  }else{
    PreEnzymeNet(sample_exp, 
                 network = coad,
                 gene_mapping = enz_gene_mapping, 
                 gene_feature = cpg_gene,
                 cell_name = sample_info$ModelID[i], 
                 exp_cutoff = 1,
                 save_path = "/home/data/sdb/wt/model_data/geo_chemo_net/")
  }
}
parallel::stopCluster(cl = my.cluster)

dt <- sample_info %>% 
  select(ModelID) %>% rename(cell = ModelID) %>% 
  mutate(cell_index = row_number() - 1)
write.csv(dt,file = "data/geo_chemo_cells.csv",row.names = F,quote = F)

####预测结果
pre <- data.table::fread("data/geo_chemo_pre.csv",data.table = F)
sample_info <- readRDS("data/geo_chemo_meta.rds")
pre_gene <-  data.frame(id = unique(c(pre$gene_name))) %>% 
  rowwise() %>% 
  mutate(gene=ensg2name(id,enz_gene_mapping,NULL,cpg_gene)[[1]]) %>% ungroup()
pre <- pre %>% 
  left_join(.,pre_gene %>% rename(gene_name=id))
saveRDS(pre,file = "data/geo_chemo_pre.rds")

###计算通路富集 OR
pre <- readRDS("data/geo_chemo_pre.rds")
kegg <- readRDS("~/meta_target/data/kegg_all_pathway.rds")
kegg <- kegg %>% 
  filter(grepl("Metabolism",class) | grepl("metabolism",pathway)) %>% 
  mutate(pathway = gsub(" \\- Homo sapiens \\(human\\)","",pathway)) %>% 
  filter(!grepl("Drug",pathway))
kegg <- kegg %>% 
  mutate(pathway = case_when(
    pathway == "Citrate cycle (TCA cycle)" ~ "Citrate cycle",
    pathway == "Glycosylphosphatidylinositol (GPI)-anchor biosynthesis" ~ "Glycosylphosphatidylinositol",
    TRUE ~ pathway
  ))

get_pathway <- function(gene,pathway){
  gene <- strsplit(gene,split = ",")[[1]]
  pathway <- pathway %>% filter(genes %in% gene)
  return(c(paste(unique(pathway$pathway),collapse = ";"),
           paste(unique(pathway$class),collapse = ";")))
}
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()

get_pathway_p <- function(pathway_dt,pathway_name,sample_name){
  pdt <- pathway_dt %>% 
    filter(grepl(pathway_name,pathway) & cell == sample_name)
  npdt <- pathway_dt %>% 
    filter(!grepl(pathway_name,pathway) & cell == sample_name)
  res <- fisher.test(cbind(
    c(sum(pdt$preds == 1),sum(pdt$preds != 1)),
    c(sum(npdt$preds == 1),sum(npdt$preds != 1))
  ),alternative = "greater")
  return(c(res$p.value,res$estimate[[1]]))
}

all_samples <- unique(pre_pathway$cell)

library(doParallel)
library(foreach)

all_split <- pre_pathway %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";")

my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)
res <- foreach(
  i = all_samples,
  .export = c("pre_pathway","get_pathway_p"),
  .packages = c("dplyr")
) %dopar% {
  p_o <- get_pathway_p(pre_pathway,"Pyrimidine metabolism",i)
  dt <- data.frame(p_value = p_o[1],
                   ratio = p_o[2],
                   sample = i,
                   pathway = "Pyrimidine metabolism")
  return(dt)
}
parallel::stopCluster(cl = my.cluster)
res <- bind_rows(res)
saveRDS(res,file = "data/geo_chemo_meta_dep_pyri.rds")

####和药物反应的关联
sample_info <- readRDS("data/geo_chemo_meta.rds")
py <- readRDS("~/DeepMeta/data/geo_chemo_meta_dep_pyri.rds")
py <- left_join(py,sample_info %>% 
                  select(Sample_id,Response,Cancer_type_level1) %>%
                  rename(sample = Sample_id))
py <- py %>% 
  mutate(OR_type = case_when(
    ratio >= 4.075275 ~ "high",
    ratio <= 3.576558 ~ "low",
    TRUE ~ "others"
  )) %>% filter(OR_type != "others")
df <- py %>% 
  group_by(OR_type,Response) %>% 
  summarise(sample_counts = length(unique(sample))) %>% 
  ungroup() 

library(ggplot2)
library(ggprism)
ggplot(data=df,aes(x=OR_type,y=sample_counts,fill=Response))+
  geom_bar(stat = "identity",position="fill")+
  theme_prism()+
  labs(y="Percent of cases (%)",title = "Fisher's Exact Test, P = 0.004756")+
  scale_fill_manual(values=c("#00FDFE","#FE0000"))+
  theme(axis.title.x = element_blank())
ggsave("report/geo_chemo_fisher.pdf",width = 6,height = 6)
dt <- py %>% 
  mutate(truth = ifelse(Response == "Response", "Class1","Class2")) %>% 
  mutate(score = -log10(p_value) * ratio)
dt$truth <- factor(dt$truth)

library(yardstick)
pr <- pr_auc(dt, truth, ratio)[".estimate"][[1]]
roc <- roc_auc(dt, truth, ratio)[".estimate"][[1]]

###按照嘧啶代谢的 GSVA 来将样本分层
tpm <- data.table::fread("/home/data/sdc/wt/TCGA/tcga_RSEM_gene_tpm.gz",data.table = F)
tpm <- tpm %>% 
  select(sample,which(as.numeric(substr(colnames(tpm),14,15)) < 10)) %>% 
  rename(id = sample)
###mapping
mapping <- data.table::fread("/home/data/sdc/wt/TCGA/probeMap_gencode.v23.annotation.gene.probemap",
                             data.table = F)
tpm <- left_join(tpm,mapping %>% select(id,gene)) %>% 
  select(-id) %>% 
  select(gene,everything())

tcga_surv_nm <- readRDS("~/DeepMeta/data/tcga_surv_nm.rds")
tpm <- tpm %>% select(gene,which(substr(colnames(tpm),1,12) %in% tcga_surv_nm$sample))
tpm <- tpm[!duplicated(tpm$gene),]
rownames(tpm) <- tpm$gene
tpm <- tpm %>% select(-gene)

kegg <- readRDS("~/meta_target/data/kegg_all_pathway.rds")
kegg <- kegg %>% 
  filter(grepl("Metabolism",class) | grepl("metabolism",pathway)) %>% 
  mutate(pathway = gsub(" \\- Homo sapiens \\(human\\)","",pathway)) %>% 
  filter(!grepl("Drug",pathway))
py <- list(`Pyrimidine metabolism`=kegg$genes[which(kegg$pathway == "Pyrimidine metabolism")])
library(GSVA)
gsvaPar <- gsvaParam(as.matrix(tpm), py)
py_score <- gsva(gsvaPar)
py_score <- t(py_score) %>% as.data.frame()
colnames(py_score) <- "NES"
py_score$sample <- rownames(py_score)
saveRDS(py_score,file = "data/tcga_py_gsva.rds")

###
tcga_py_gsva <- readRDS("~/DeepMeta/data/tcga_py_gsva.rds")
surv <- readRDS("data/pancancer_survial.rds")
surv_nm <- inner_join(surv,
                      tcga_py_gsva %>% mutate(sample = substr(sample,1,12))) %>% 
  distinct_all(.keep_all = TRUE)
surv_nm <- surv_nm %>% 
  mutate(NES_type = ifelse(NES > median(surv_nm$NES), "high","low")) 
surv_nm$NES_type <- factor(surv_nm$NES_type, levels = c("low","high"))

p <- EasyBioinfo::show_km(surv_nm ,"NES_type")
pdf("report/py_survial_by_gsva.pdf",height = 6,width = 7)
print(p)
dev.off()

fit0 <- survminer::surv_fit(survival::Surv(OS.time, OS)~NES_type, data = surv_nm)
survminer::ggsurvplot(fit0, pval = T, risk.table = T,risk.table.pos = "in", 
                      data = surv_nm)

###PRISM repurposing
dt <- data.table::fread("/home/data/sdb/wt/model_data/PRISM/primary-screen-replicate-collapsed-logfold-change.csv",
                        data.table = F)
drug_info <- data.table::fread("/home/data/sdb/wt/model_data/PRISM/primary-screen-replicate-collapsed-treatment-info.csv",
                               data.table = F)

train_genes <- data.table::fread("/home/wt/meta_target/data/train_genes.csv",
                                 data.table = F,sep = ",")
enz_gene_mapping <- readRDS("~/meta_target/data/enz_gene_mapping.rds")
ensg2name <- function(ensg,mapping){
  tt <- strsplit(ensg," and ")[[1]]
  tt_gene <- mapping$symbol[which(mapping$ensembl_id %in% tt)] %>% 
    unique() 
  return(tt_gene)
}
train_genes <- train_genes %>% 
  rowwise() %>% 
  mutate(gene = paste0(ensg2name(id,enz_gene_mapping),collapse = ",")) %>% 
  ungroup()
all_genes <- paste0(train_genes$gene,collapse = ",") %>% 
  strsplit(.,",")
all_genes <- all_genes[[1]]

drug_info <- drug_info %>% 
  rowwise() %>% 
  mutate(is_train = any(strsplit(target,",")[[1]] %in% all_genes)) %>% 
  ungroup()
drug_info <- drug_info %>% filter(is_train)

test_pre <- data.table::fread("data/test_preV2.csv",data.table = F)
dt <- dt %>% 
  filter(V1 %in% test_pre$cell) %>% 
  select(V1,drug_info$column_name)
test_pre <- test_pre %>% 
  filter(cell %in% dt$V1) %>% 
  rowwise() %>% 
  mutate(gene = paste0(ensg2name(gene_name,enz_gene_mapping),collapse = ",")) %>% 
  ungroup()

get_drug_score<- function(cell_name,gene_name,drug_dt,score_dt){
  gene_drug <- drug_dt %>% 
    filter(target == gene_name)
  if (nrow(gene_drug) != 0){
    ids <- unique(gene_drug$column_name)
    gene_score <- score_dt %>% 
      filter(V1 == cell_name) %>% 
      select(ids) %>% t() %>% as.data.frame() %>% rename(score = V1)
    gene_score$ids <- rownames(gene_score)
    gene_score <- left_join(gene_score,
                            gene_drug %>% 
                              select(column_name,name) %>% rename(ids = 1),
                            by="ids")
    if (all(is.na(gene_score$score))){
      return(list(
        drug_name = NA,
        score = NA
      ))
    }else{
      return(list(
        drug_name = gene_score$name[which.min(gene_score$score)],
        score = gene_score$score[which.min(gene_score$score)]
      ))
    }
  }else{
    return(list(
      drug_name = NA,
      score = NA
    ))
  }
}
drug_info_split <- drug_info %>% 
  tidyr::separate_longer_delim(cols = target, delim = ", ")
test_pre <- test_pre %>% 
  tidyr::separate_longer_delim(cols = gene, delim = ",") %>% 
  rowwise() %>% 
  mutate(min_score = get_drug_score(cell,gene,drug_info_split,dt)$score,
         min_drug = get_drug_score(cell,gene,drug_info_split,dt)$drug_name) %>% 
  ungroup()
saveRDS(test_pre,file = "data/cell_test_PRISM.rds")

test_pre <- test_pre %>% 
  filter(!is.na(min_score))

library(ggpubr)
ggboxplot(data=test_pre,x="preds",y="min_score",
          xlab = "Prediction Label",ylab = "Drug Score")+
  stat_compare_means()
ggsave("report/PRISM_test_pre.pdf",width = 4,height = 5)

###基因分布
pre <- read.csv("data/test_preV2.csv") %>% select(-X)
pre_summ <- pre %>% 
  group_by(gene_name) %>% 
  summarise(count1 = sum(preds == 1),counts0=sum(preds==0),
            true0=sum(label == 0),true1 = sum(label == 1)) %>% 
  ungroup() %>% 
  filter(true0 != 0 & true1 != 0)

get_conf <- function(gene){
  tt <- pre %>% filter(gene_name == gene)
  return(c(sum(tt$label == 0 & tt$preds == 0),
           sum(tt$label == 0 & tt$preds == 1),
           sum(tt$label == 1 & tt$preds == 0),
           sum(tt$label == 1 & tt$preds == 1)))
}

dt <- matrix(rep(1,(556*2)^2),nrow = 556 * 2)
for (i in 1:nrow(pre_summ)){
  tmp_m <- get_conf(pre_summ$gene_name[i])
  dt[((2*i)-1):(((2*i)-1)+1),((2*i)-1)] <- tmp_m[1:2]
  dt[((2*i)-1):(((2*i)-1)+1),((2*i)-1)+1] <- tmp_m[3:4]
}
 
###多组学
train_cell_info <- read.csv("/home/data/sdc/wt/model_data/new_model/cell_net_filter_exp/raw/train_cell_info.csv")

mut <- data.table::fread("/home/data/sdb/wt/model_data/OmicsSomaticMutations.csv",
                         data.table = F)
mut <- mut %>% 
  select(DepMap_ID,HugoSymbol,VariantInfo,Chrom,Pos,Ref,Alt) %>% 
  filter(VariantInfo != "SILENT") %>% 
  filter(DepMap_ID %in% c(train_cell_info$cell))
train_cell_info <- train_cell_info %>% filter(cell %in% mut$DepMap_ID)
mut_summ <- mut %>% 
  group_by(HugoSymbol) %>% 
  summarise(cell_c=length(unique(DepMap_ID))) %>% 
  ungroup() %>% 
  filter(cell_c > 20)
mut <- mut %>% 
  filter(HugoSymbol %in% mut_summ$HugoSymbol) %>% 
  select(DepMap_ID,HugoSymbol) %>% 
  distinct_all() %>% 
  rename(cell=DepMap_ID,gene=HugoSymbol) %>% 
  mutate(value = 1) %>% 
  tidyr::pivot_wider(names_from = gene, values_from = value, values_fill = 0)
write.csv(mut,file = "/home/data/sdb/wt/model_data/mut_dt.csv",
          quote = F, row.names = F)

cnv <- data.table::fread("/home/data/sdb/wt/model_data/OmicsCNGene.csv",
                         data.table = F)
colnames(cnv)[2:ncol(cnv)] <- gsub("\\s*\\([^\\)]+\\)","",colnames(cnv)[2:ncol(cnv)])
colnames(cnv)[1] <- "cell"
gene_cnv_sd <- apply(cnv[,2:ncol(cnv)],2,sd)
gene_cnv_sd <- gene_cnv_sd[which(gene_cnv_sd > 0.2)]
cnv <- cnv %>% select(cell, names(gene_cnv_sd))
cnv <- cnv %>% filter(cell %in% train_cell_info$cell)
write.csv(cnv,file = "/home/data/sdb/wt/model_data/cnv_dt.csv",
          quote = F, row.names = F)

train_cell_info <- train_cell_info %>% mutate(cell_index = row_number()-1)
write.csv(train_cell_info,file = "/home/data/sdb/wt/model_data/train_cell_info_omics.csv",
          quote = F, row.names = F)

###res
library(dplyr)
library(yardstick)
library(ggprism)
library(ggplot2)
library(ggpubr)

get_roc <- function(dt,need="all"){
  dt <- dt %>% 
    mutate(truth = ifelse(label == 1, "Class1","Class2"),
           pred_label = ifelse(preds == 1, "Class1","Class2"))
  dt$truth <- factor(dt$truth)
  dt$pred_label <- factor(dt$pred_label)
  
  f1 <- try(f_meas(dt,truth,pred_label)[".estimate"] %>%
              unlist() %>% unname() %>% round(.,2),silent = TRUE)
  kappa <- try(kap(dt,truth,pred_label)[".estimate"] %>%
                 unlist() %>% unname() %>% round(.,2),silent = TRUE)
  f1 <- ifelse('try-error' %in% class(f1),NA,f1)
  kappa <- ifelse('try-error' %in% class(kappa),NA,kappa)
  if (need == "all"){
    roc <-  roc_auc(dt, truth, preds_raw)[".estimate"] %>% 
      unlist() %>% unname() %>% round(.,2)
    pr <-  pr_auc(dt, truth, preds_raw)[".estimate"] %>%
      unlist() %>% unname() %>% round(.,2)
    return(c(roc,pr,f1,kappa))
  }else{
    res <- switch(
      need,
      "roc" = roc_auc(dt, truth, preds_raw)[".estimate"] %>% 
        unlist() %>% unname() %>% round(.,2),
      "pr" = pr_auc(dt, truth, preds_raw)[".estimate"] %>%
        unlist() %>% unname() %>% round(.,2),
      "f1" = f1,
      "kappa" = kappa
    )
    return(res)
  }
  
}

all_cv <- c("own_model","multi_omics")
cv_res <- vector("list",2)
for (i in 1:2){
  sub_res <- vector("list",10)
  for (j in 1:10){
    pre <- read.csv(paste0("data/cv/",all_cv[i],"/fold_",j-1,".csv")) %>% select(-X)
    pre_roc <- ifelse(i %in% c(1:7),get_roc(pre,need = "roc"),NA)
    pre_pr <- ifelse(i %in% c(1:7),get_roc(pre,need = "pr"),NA)
    pre_f1 <- get_roc(pre,need = "f1")
    pre_kap <- get_roc(pre,need = "kappa")
    
    sub_res[[j]] <- data.frame(
      fold = paste0("Fold-",j),
      ROC = pre_roc,
      PR = pre_pr,
      F1 = pre_f1,
      Kappa = pre_kap
    ) 
  }
  sub_res <- bind_rows(sub_res)
  sub_res$type <- all_cv[i]
  cv_res[[i]] <- sub_res
}
cv_res <- bind_rows(cv_res)
cv_res <- cv_res %>% 
  tidyr::pivot_longer(cols = c("ROC","PR","F1","Kappa"),
                      names_to = "Metric",values_to = "Value")
ggbarplot(cv_res, x = "type", y = "Value",fill="Metric",
          add = "mean_se", label = TRUE, 
          lab.vjust = -0.5,position = position_dodge(0.9))+
  scale_x_discrete(labels=c("DeepMeta","Multi-Omics (EXP+CNV+MUT)"))+
  labs(x="",y='Value')
ggsave("report/multi_omics_compare.pdf",width = 8,height = 5)



