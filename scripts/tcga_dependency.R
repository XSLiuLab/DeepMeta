library(parallel)
library(dplyr)
library(foreach)

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

pre <- readRDS("data/tcga_pre_V2.rds") %>% as.data.frame()
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()

####看看 pathway 中 positive 基因的占比
all_pathways <- unique(strsplit(paste(pre_pathway$pathway,collapse = ";"),
                                split = ";")[[1]])
all_samples <- unique(pre_pathway$cell)

all_split <- pre_pathway %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";")

pre_pos_summ <- pre_pathway %>%
  group_by(cell) %>%
  summarise(pos_counts = sum(preds == 1)) %>% ungroup()
pre_pathways_summ <- all_split %>% 
  group_by(cell,pathway) %>% 
  summarise(pathway_pos_counts = sum(preds == 1)) %>% ungroup() %>% 
  left_join(.,pre_pos_summ) %>% 
  mutate(pos_percent = pathway_pos_counts/pos_counts) %>% 
  arrange(desc(pos_percent))

###去掉少于 2000 个样本中有的
sample_couts <- pre_pathways_summ %>%
  group_by(pathway) %>%
  summarise(counts = length(unique(cell))) %>% ungroup() %>% 
  filter(counts > 2000)

pathways_summ_summ <- pre_pathways_summ %>% 
  filter(pathway %in% sample_couts$pathway) %>% 
  group_by(pathway) %>% 
  summarise(mm = median(pos_percent)) %>% ungroup() %>% 
  arrange(desc(mm)) %>% 
  filter(mm > 0)
library(ggpubr)
pre_pathways_summ_filter <- pre_pathways_summ %>% 
  filter(pathway %in% pathways_summ_summ$pathway)
p1 <- ggboxplot(data = pre_pathways_summ_filter,x="pathway",y="pos_percent",
                order = pathways_summ_summ$pathway,xlab = FALSE,
                ylab = "Proportion of positive genes")+
  rotate_x_text()
p1
ggsave("figs/pathway_summ_tcga.pdf",width = 14,height = 8)

####
#####富集分析
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

#get_pathway_p(pre_pathway,"Pyrimidine metabolism","TCGA-HC-7210-01")

all_samples <- unique(pre_pathway$cell)

library(doParallel)
library(foreach)

all_split <- pre_pathway %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";")
all_pathways <- unique(all_split$pathway)
all_res <- vector("list",84)
for (i in 1:length(all_pathways)){
  tmp <- all_pathways[i]
  my.cluster <- parallel::makeCluster(
    60, 
    type = "PSOCK"
  )
  #register it to be used by %dopar%
  doParallel::registerDoParallel(cl = my.cluster)
  
  res <- foreach(
    i = all_samples,
    .export = c("pre_pathway","get_pathway_p","tmp"),
    .packages = c("dplyr")
  ) %dopar% {
    p_o <- get_pathway_p(pre_pathway,tmp,i)
    dt <- data.frame(p_value = p_o[1],
                     ratio = p_o[2],
                     sample = i,
                     pathway = tmp)
    return(dt)
  }
  parallel::stopCluster(cl = my.cluster)
  res <- bind_rows(res)
  all_res[[i]] <- res
  message("Complete ",i)
}
res <- bind_rows(all_res)
saveRDS(res,file = "data/pancancer_meta_dep.rds")

###
res <- readRDS("data/pancancer_meta_dep.rds")
res <- res %>% filter(pathway %in% sample_couts$pathway)
##每个样本中做 FDR 矫正
res_p <- res %>% 
  group_by(sample) %>% 
  mutate(padj = p.adjust(p_value,"fdr")) %>% 
  ungroup()

res_p_filter <- res_p %>% 
  filter(!is.infinite(ratio)) %>% 
  mutate(log2ratio=log2(ratio+1))
res_p_filter$cancer <- EasyBioinfo::get_cancer_type(res_p_filter$sample,
                                                    parallel = TRUE,cores = 30)
res_r_summ <- res_p_filter %>% 
  group_by(pathway) %>% 
  summarise(r_m = median(ratio),
            sample_counts = length(unique(sample))) %>% ungroup() %>% 
  arrange(desc(r_m)) %>% 
  filter(r_m > 1)
res_p_filter <- res_p_filter %>% filter(pathway %in% res_r_summ$pathway)

res_p_summ <- res_p %>% 
  group_by(pathway) %>% 
  summarise(per_p = round((sum(p_value<0.05)/n()) * 100,2),
            per_fdr = round((sum(padj<0.1)/n()) * 100,2)) %>% 
  ungroup() %>% 
  arrange(desc(per_p)) %>% 
  filter(pathway %in% res_r_summ$pathway)

res_all_summ <- inner_join(res_r_summ,res_p_summ)

library(ggpubr)
p2 <- ggboxplot(data = res_p_filter,x="pathway",y="log2ratio",ylab = "log2OR",
                order = res_r_summ$pathway,xlab = FALSE)+
  rotate_x_text()

p3 <- ggbarplot(data = res_all_summ,x="pathway",y="per_p",
                label = TRUE,
                lab.nb.digits=2,ylab = "Significant proportion",xlab = FALSE)+
  rotate_x_text()

library(patchwork)
p3
ggsave("figs/tcga_pathway_fisher_p.pdf",width = 13,height = 7)

p2
ggsave("figs/tcga_pathway_or.pdf",width = 13,height = 7)

###展示不同癌症类型中的 OR 值
res_p_filter <- res_p %>% 
  filter(!is.infinite(ratio)) %>% 
  mutate(log2ratio=log2(ratio+1))
res_p$cancer <- EasyBioinfo::get_cancer_type(res_p$sample,
                                             parallel = TRUE,cores = 30)
cancer_summ <- res_p %>% 
  group_by(cancer,pathway) %>% 
  summarise(median_or = median(ratio,na.rm=T)) %>% ungroup() %>% 
  tidyr::pivot_wider(names_from = cancer,values_from = median_or) %>% 
  as.data.frame()
rownames(cancer_summ) <- cancer_summ$pathway
cancer_summ <- cancer_summ %>% select(-pathway)
###把全是 0 的去掉
which0 <- which(apply(cancer_summ,1,function(x){all(x == 0)}))
cancer_summ <- cancer_summ[-which0,]
cancer_summ <- as.matrix(cancer_summ)
cancer_summ <- log2(cancer_summ+1)
cancer_summ[which(is.infinite(cancer_summ))] <- NA
library(ComplexHeatmap)
library(circlize)
col_fun = colorRamp2(c(0, max(cancer_summ)), c("white", "blue"))
p1 <- Heatmap(cancer_summ,cluster_rows = F,cluster_columns = F,
              rect_gp = gpar(col = "grey", lwd = 2),row_names_side = "left",
              show_heatmap_legend=T,row_names_gp = gpar(fontsize = 8),
              name="log2OR",na_col="black")
pdf("figs/pancancer_pathway_or_cancertype.pdf",height = 10,width = 10)
draw(p1)
dev.off()

####展示不同pathway基因预测值的分布
all_split <- pre_pathway %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";")

all_split <- all_split %>% filter(pathway %in% rownames(cancer_summ))
library(ggridges)
library(viridis)
library(hrbrthemes)

# Plot
all_split <- all_split %>% 
  mutate(pathway2 = stringr::str_replace(pathway," / ","-"))
gene_summ <- all_split %>% 
  group_by(pathway2) %>% 
  summarise(gene_counts = length(unique(gene))) %>% ungroup() %>% 
  filter(gene_counts >=10)
all_split <- all_split %>% 
  filter(pathway2 %in% gene_summ$pathway2) %>% 
  left_join(.,gene_summ) %>% 
  mutate(pathway3 = paste0(pathway2," (n=",gene_counts,")"))

ggplot(all_split, aes(x = preds_raw, y = pathway3, fill = stat(x))) +
  geom_density_ridges_gradient(scale = 1, rel_min_height = 0.01) +
  scale_fill_viridis(name = "Preds", option = "C") +
  theme_minimal() +
  theme(
    legend.position="none",
    panel.spacing = unit(0.1, "lines"),
    strip.text.x = element_text(size = 8)
  )+labs(x="Prediction",y="Pathways")
ggsave("figs/tcga_pathway_pre_dis2.pdf",width = 10,height = 8)

####用随机模拟的方法得到 P 值
##比如 Terpenoid backbone biosynthesis 在 TCGA-95-7947-01 样本中有两个基因预测值都是 1，
###那么在这个样本中随机选 2 个基因其中是 1 的比例是多少，和 100 % 进行比较，得到 P 值
###
get_random_p <- function(pre_dt,orign_value,sampling_counts){
  sampling1000 <- sapply(1:1000,
                         function(x,y,num){
                           dt <- y %>% slice_sample(n=num)
                           mean(dt$preds == 1)
                         }, pre_dt %>% select(preds),
                         sampling_counts)
  p_value <- mean(sampling1000 >= orign_value)
  return(p_value)
}
all_split <- pre_pathway %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";")

##只保留在所有样本中都有的
all_split_summ <- all_split %>% 
  group_by(pathway) %>% 
  summarise(sample_couts = length(unique(cell))) %>% ungroup()
all_split_summ <- all_split_summ %>% 
  filter(sample_couts >= 2000)
all_split <- all_split %>% 
  filter(pathway %in% all_split_summ$pathway)
pathway_samples <- all_split %>% 
  select(cell,pathway) %>% distinct_all()

###并行计算
my.cluster <- parallel::makeCluster(
  90, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:nrow(pathway_samples),
  .export = c("pre_pathway","get_random_p","pathway_samples"),
  .packages = c("dplyr")
) %dopar% {
  tt <- pre_pathway %>% filter(cell == pathway_samples$cell[i])
  tt1 <- tt %>% filter(grepl(pathway_samples$pathway[i],pathway))
  p <- get_random_p(tt,mean(tt1$preds == 1),nrow(tt1))
  dt <- data.frame(
    sample = pathway_samples$cell[i],
    pathway = pathway_samples$pathway[i],
    p_value = p
  )
  return(dt)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
saveRDS(res,file = "data/tcga_pathway_permutation.rds")

###同一个样本进行假设检验
res_p <- res %>% 
  group_by(sample) %>% 
  mutate(padj = p.adjust(p_value,"fdr")) %>% 
  ungroup()

###P value 分布
res_p_summ <- res %>% 
  group_by(pathway) %>% 
  summarise(sig_per = mean(p_value<0.05)) %>% ungroup() %>% 
  filter(sig_per > 0) %>% 
  arrange(desc(sig_per)) %>% 
  mutate(sig_per_100 = round(100 * sig_per,2))

library(ggpubr)
res_filter <- res %>% filter(pathway %in% res_p_summ$pathway)
ggboxplot(data=res_filter,x="pathway",y="p_value",order = res_p_summ$pathway)

ggbarplot(data = res_p_summ,x="pathway",y="sig_per_100",
          label = TRUE,
          lab.nb.digits=2,ylab = "Significant proportion",xlab = FALSE)+
  rotate_x_text()
ggsave("figs/tcga_pathway_permutation.pdf",width = 13,height = 7)

####对正常样本也做这样的分析
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

pre <- readRDS("data/tcga_normal_pre_V2.rds") %>% as.data.frame()
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()

all_split <- pre_pathway %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";")

##只保留 tumor 中有的
cancer_res <- readRDS("data/tcga_pathway_permutation.rds")
all_split <- all_split %>% 
  filter(pathway %in% cancer_res$pathway)
pathway_samples <- all_split %>% 
  select(cell,pathway) %>% distinct_all()

###跑前面那些一样的步骤
saveRDS(res,file = "data/tcga_pathway_permutation_normal.rds")


res_p_summ <- res %>% 
  group_by(pathway) %>% 
  summarise(sig_per = mean(p_value<0.05)) %>% ungroup() %>% 
  filter(sig_per > 0) %>% 
  arrange(desc(sig_per)) %>% 
  mutate(sig_per_100 = round(100 * sig_per,2))



