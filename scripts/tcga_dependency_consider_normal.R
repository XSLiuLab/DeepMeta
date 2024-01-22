library(dplyr)
library(ComplexHeatmap)
library(circlize)
get_pathway <- function(gene,pathway){
  gene <- strsplit(gene,split = ",")[[1]]
  pathway <- pathway %>% filter(genes %in% gene)
  return(c(paste(unique(pathway$pathway),collapse = ";"),
           paste(unique(pathway$class),collapse = ";")))
}

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

pre <- readRDS("data/tcga_pre_V2.rds") %>% as.data.frame()
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)

pre_summ <- pre %>% 
  filter(gene %in% pre_pathway$gene) %>% 
  group_by(cancer,gene) %>% 
  summarise(cancer_pos_counts = sum(preds == 1)) %>% ungroup()
cancer_summ <- pre %>% 
  filter(gene %in% pre_pathway$gene) %>% 
  group_by(cancer,gene) %>% 
  summarise(cancer_sample_counts = length(unique(cell))) %>% ungroup()
pre_summ <- left_join(pre_summ,cancer_summ) %>% 
  mutate(cancer_gene = paste(cancer,gene,sep = "_")) %>% 
  select(-cancer,-gene)

normal_pre <- readRDS("data/tcga_normal_pre_V2.rds")
pre_pathway <- data.frame(gene = unique(c(normal_pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)

normal_pre_summ <- normal_pre %>% 
  filter(gene %in% pre_pathway$gene) %>% 
  group_by(cancer,gene) %>% 
  summarise(pos_counts = sum(preds == 1)) %>% ungroup()
normal_cancer_summ <- normal_pre %>% 
  filter(gene %in% pre_pathway$gene) %>% 
  group_by(cancer, gene) %>% 
  summarise(sample_counts = length(unique(cell))) %>% ungroup()
normal_pre_summ <- left_join(normal_pre_summ,normal_cancer_summ) %>% 
  mutate(cancer_gene = paste(cancer,gene,sep = "_")) %>% 
  select(-cancer,-gene)

cancer_normal_vs <- inner_join(
  normal_pre_summ,pre_summ
)

###对正常样本也做富集分析
pre <- readRDS("data/tcga_normal_pre_V2.rds") %>% as.data.frame()
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()
all_samples <- unique(pre_pathway$cell)

library(doParallel)
library(foreach)

all_split <- pre_pathway %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";")
all_pathways <- unique(all_split$pathway)
all_res <- vector("list",82)
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
saveRDS(res,file = "data/pancancer_meta_dep_normal.rds")

###
###对单基因在每个癌症类型中进行 fisher
get_fisher <- function(normal_pos,normal_counts,cancer_pos,cancer_counts){
  dt <- fisher.test(cbind(c(cancer_pos,cancer_counts-cancer_pos),
                          c(normal_pos,normal_counts-normal_pos)),
                    alternative = "greater")
  return(c(dt$estimate,dt$p.value))
}
cancer_normal_vs <- cancer_normal_vs %>% 
  rowwise() %>% 
  mutate(or = get_fisher(pos_counts,sample_counts,
                         cancer_pos_counts,cancer_sample_counts)[1],
         p = get_fisher(pos_counts,sample_counts,
                        cancer_pos_counts,cancer_sample_counts)[2]) %>% 
  ungroup()
cancer_normal_vs$padj <- p.adjust(cancer_normal_vs$p,"fdr")

###先用 Fisher 的方法进行不同癌症类型的 P 值合并
sig <- cancer_normal_vs %>% 
  filter(padj < 0.1) %>% 
  tidyr::separate_wider_delim(cols = "cancer_gene",delim = "_",
                              names = c("cancer","gene"),
                              cols_remove = FALSE)

###
sig1 <- sig %>% 
  group_by(gene) %>% 
  summarise(combine_p = ifelse(length(unique(cancer)) == 1, padj,
                               metap::sumlog(padj)$p)) %>% ungroup()
sig1 <- sig1 %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1) %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";") %>% 
  mutate(log10P = (-log10(combine_p)))

##展示最少 2 个基因的通路
sig_summ <- sig1 %>% 
  group_by(pathway) %>% 
  summarise(gene_counts = length(unique(gene))) %>% ungroup() %>% 
  filter(gene_counts > 1) %>% 
  arrange(desc(gene_counts))
sig1 <- sig1 %>% filter(pathway %in% sig_summ$pathway)
sig1$pathway <- factor(sig1$pathway,levels = sig_summ$pathway)

ggplot(sig1,aes(x=gene,y=log10P))+
  geom_bar(stat = "identity") + 
  facet_grid(. ~  pathway, scale = "free_x", labeller = label_wrap_gen(width=15))+
  theme_prism()+
  theme(axis.text.x = element_text(angle = 45, vjust = 1, hjust=1))+
  theme(axis.title.x = element_blank())+
  ylab("-Log10P")
ggsave("figs/cancer_vs_normal_sig_2.pdf",width = 16,height = 6)

###展示所有显著的基因-癌症类型
dt <- sig %>% select(cancer, gene, or, padj) %>% 
  mutate(log10p = (-log10(padj)))
dt_p <- dt %>% 
  select(cancer,gene,log10p) %>% 
  tidyr::pivot_wider(names_from = cancer, values_from = log10p) %>%
  as.data.frame()
rownames(dt_p) <- dt_p$gene
dt_p$gene <- NULL
dt_p <- as.matrix(dt_p)

dt_or <- dt %>% 
  select(cancer,gene,or) %>% 
  mutate(or = sprintf(or, fmt = "%0.1f")) %>% 
  tidyr::pivot_wider(names_from = cancer, values_from = or) %>%
  as.data.frame()
rownames(dt_or) <- dt_or$gene
dt_or$gene <- NULL
dt_or <- as.matrix(dt_or)

library(circlize)
col_fun <- colorRamp2(c(1, 15), c("white", "red"))
p1 <- Heatmap(dt_p,cluster_rows = F,cluster_columns = F,
              na_col = "blue",col = col_fun,
              rect_gp = gpar(col = "white", lwd = 1),
              border_gp = gpar(col = "black", lty = 1),
              column_names_gp = grid::gpar(fontsize = 8),
              row_names_gp = grid::gpar(fontsize = 10),name = "-log10(Padj)",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(!is.na(dt_or[i, j])){
                  grid.text(dt_or[i, j], x, y, 
                            gp = gpar(fontsize = 10))
                }
              })
pdf(file="figs/cancer_vs_normal_sig_all.pdf",width = 8,height = 8)
draw(p1)
dev.off()

####癌症和正常样本中的 OR 值比较
cancer_dep <- readRDS("data/pancancer_meta_dep.rds") %>% 
  filter(!is.infinite(ratio)) %>% 
  mutate(log2ratio=log2(ratio+1)) %>% 
  mutate(type = "cancer")
normal_dep <- readRDS("data/pancancer_meta_dep_normal.rds") %>% 
  filter(!is.infinite(ratio)) %>% 
  mutate(log2ratio=log2(ratio+1)) %>% 
  mutate(type = "normal")

both <- bind_rows(normal_dep,cancer_dep)
both$cancer_type <- EasyBioinfo::get_cancer_type(both$sample,cores = 50, parallel = T)
##全局比较
cancer_summ <- both %>% 
  group_by(pathway,cancer_type) %>% 
  summarise(counts = length(unique(type))) %>% ungroup() %>% 
  filter(counts > 1)
cancer_summ2 <- both %>% 
  group_by(cancer_type,type) %>% 
  summarise(counts = length(unique(sample))) %>% 
  filter(type == "normal") %>% 
  filter(counts > 10)
need_cancer <- intersect(unique(cancer_summ$cancer_type),
                         unique(cancer_summ2$cancer_type))
both_summ <- both %>% 
  filter(cancer_type %in% need_cancer) %>% 
  group_by(pathway, type) %>% 
  summarise(median_or = median(log2ratio)) %>% ungroup() %>% 
  tidyr::pivot_wider(names_from = "type",values_from = "median_or") %>% 
  filter(cancer > normal)
dt <- both %>% 
  filter(cancer_type %in% need_cancer) %>% 
  filter(pathway %in% both_summ$pathway) %>% 
  rename(Type = type)
dt$Type <- factor(dt$Type,levels = c("normal","cancer"))
###按照 P 值进行排序
dt_p <- dt %>% group_by(pathway) %>% 
  do(w = wilcox.test(log2ratio ~ Type, data=., paired=FALSE, 
                     alternative = "less")) %>% 
  summarise(pathway, Wilcox = w$p.value) %>% 
  arrange(Wilcox)

ggboxplot(data = dt, x="pathway", y= "log2ratio", fill= "Type",
          ylab = "Log2OR",xlab = FALSE,
          order = dt_p$pathway,
          palette = c("#81E9B8","#FE7D81"))+
  stat_compare_means(aes(group = Type,
                         label = paste0("p = ", after_stat(p.format))),
                     method = "wilcox",
                     method.args = list(alternative = "greater") )+
  rotate_x_text(45)
ggsave("figs/pathway_or_cancer_vs_normal.pdf",width = 10, height = 8)
ggsave("figs/pathway_or_cancer_vs_normal.pdf",width = 14, height = 8)
##分癌症类型 tumor or 减去 normal or 的中位数
###分癌症类型
both_summ <- both %>% 
  group_by(pathway, cancer_type, type) %>% 
  summarise(median_or = median(log2ratio)) %>% 
  ungroup()
cancer_summ <- both %>% 
  group_by(pathway,cancer_type) %>% 
  summarise(counts = length(unique(type))) %>% ungroup() %>% 
  filter(counts > 1)

both_summ <- both_summ %>% 
  filter(cancer_type %in% cancer_summ$cancer_type) %>% 
  tidyr::pivot_wider(names_from = "type",values_from = "median_or")
both_summ <- both_summ %>% 
  mutate(diff_or = cancer - normal) 
both %>% group_by(cancer_type,type) %>% 
  summarise(counts = length(unique(sample))) %>% 
  filter(type == "normal") %>% 
  filter(counts > 10) -> tt

###只选择正常样本数量大于 10 
dt <- both_summ %>% 
  filter(cancer_type %in% tt$cancer_type) %>% 
  select(pathway, cancer_type, diff_or) %>% 
  tidyr::pivot_wider(names_from = cancer_type,values_from = diff_or) %>% 
  as.data.frame()
rownames(dt) <- dt$pathway
dt <- dt %>% select(-pathway)
###把全是 0 的去掉，也就是这个通路在所有癌症类型中都没有差异
which0 <- which(apply(dt,1,function(x){all(x == 0)}))
dt <- dt[-which0,]
dt <- as.matrix(dt)

library(circlize)
col_fun = colorRamp2(c(min(dt,na.rm = T), max(dt,na.rm = T)),
                     c("white", "red"))
p1 <- Heatmap(dt,cluster_rows = F,cluster_columns = F,
              rect_gp = gpar(col = "grey", lwd = 2),row_names_side = "left",
              show_heatmap_legend=T,row_names_gp = gpar(fontsize = 8),
              name="log2OR Cancer VS Normal",na_col="black",
              row_order = c(dt_p$pathway,
                            rownames(dt)[which(!(rownames(dt) %in% dt_p$pathway))]))
p1

pdf("figs/pancancer_pathway_or_cancertype_cancer_vs_normal.pdf",
    height = 8,width = 8)
draw(p1)
dev.off()

####癌症和正常样本中显著比例比较，基于 permutation test
cancer_res <- readRDS("data/tcga_pathway_permutation.rds")
normal_res <- readRDS("data/tcga_pathway_permutation_normal.rds")
both <- bind_rows(
  cancer_res %>% mutate(type = "cancer"),
  normal_res %>% mutate(type = "normal")
)
both$cancer_type <- EasyBioinfo::get_cancer_type(both$sample,cores = 50, 
                                                 parallel = TRUE)
cancer_summ2 <- both %>% 
  group_by(cancer_type,type) %>% 
  summarise(counts = length(unique(sample))) %>% 
  filter(type == "normal") %>% 
  filter(counts > 10)

both_summ <- both %>% 
  filter(cancer_type %in% cancer_summ2$cancer_type) %>% 
  group_by(pathway,cancer_type,type) %>% 
  summarise(sig_per = mean(p_value < 0.05)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(names_from = "type", values_from = "sig_per")
both_summ <- both_summ %>% 
  filter(cancer > normal)

both_summ <- both_summ %>% 
  tidyr::pivot_longer(cols = c("cancer","normal"),names_to = "type") %>% 
  mutate(cancer_type = paste(cancer_type, type, sep = "_"))

dt <- both_summ %>% 
  select(pathway, cancer_type, value) %>% 
  tidyr::pivot_wider(names_from = cancer_type, values_from = value) %>%
  as.data.frame()
rownames(dt) <- dt$pathway
dt$pathway <- NULL
dt <- as.matrix(dt)

dt <- t(dt)
library(circlize)
col_fun <- colorRamp2(c(0, 1), c("white", "red"))
p1 <- Heatmap(dt,cluster_rows = F,cluster_columns = F,
              na_col = "blue",col = col_fun,
              rect_gp = gpar(col = "white", lwd = 1),
              border_gp = gpar(col = "black", lty = 1),
              column_names_gp = grid::gpar(fontsize = 10),
              row_names_gp = grid::gpar(fontsize = 10),name = "Significant proportion",
              cell_fun = function(j, i, x, y, width, height, fill) {
                if(!is.na(dt[i, j])){
                  grid.text(sprintf("%.1f", dt[i, j] * 100), x, y, 
                            gp = gpar(fontsize = 9))
                }
              },
              row_split = rownames(dt) %>% gsub("_.+","",.),
              row_labels = rownames(dt) %>% gsub(".+_","",.),
              row_title_rot = 0)
p1
pdf(file="figs/permutation_tumor_vs_normal.pdf",width = 8,height = 10)
draw(p1)
dev.off()

###
cancer_res <- readRDS("data/tcga_pathway_permutation.rds")
normal_res <- readRDS("data/tcga_pathway_permutation_normal.rds")
both <- bind_rows(
  cancer_res %>% mutate(type = "cancer"),
  normal_res %>% mutate(type = "normal")
)
both$cancer_type <- EasyBioinfo::get_cancer_type(both$sample,cores = 50, 
                                                 parallel = TRUE)
cancer_summ2 <- both %>% 
  group_by(cancer_type,type) %>% 
  summarise(counts = length(unique(sample))) %>% 
  filter(type == "normal") %>% 
  filter(counts > 10)

both_summ <- both %>% 
  filter(cancer_type %in% cancer_summ2$cancer_type) %>% 
  group_by(pathway,cancer_type,type) %>% 
  summarise(sig_per = mean(p_value < 0.05)) %>% 
  ungroup() %>% 
  tidyr::pivot_wider(names_from = "type", values_from = "sig_per")

both_summ <- both_summ %>% 
  mutate(diff = cancer - normal) %>%
  filter(pathway %in% rownames(dt))

dt1 <- both_summ %>% 
  select(pathway, cancer_type, diff) %>% 
  tidyr::pivot_wider(names_from = cancer_type, values_from = diff) %>%
  as.data.frame()
rownames(dt1) <- dt1$pathway
dt1$pathway <- NULL
dt1 <- as.matrix(dt1)

col_fun = colorRamp2(c(min(dt1,na.rm = T), max(dt1,na.rm = T)),
                     c("white", "red"))
p1 <- Heatmap(dt1,cluster_rows = F,cluster_columns = F,
              rect_gp = gpar(col = "grey", lwd = 2),
              row_names_side = "left",
              show_heatmap_legend=T,row_names_gp = gpar(fontsize = 8),
              name="Significant sample proportion \nCancer VS Normal",
              na_col="black",
              row_order = c(dt_p$pathway,
                            rownames(dt1)[which(!(rownames(dt1) %in% dt_p$pathway))]))
p1

pdf("~/meta_target/figs/pancancer_pathway_proportion_cancertype_cancer_vs_normal.pdf",
    height = 8,width = 8.5)
draw(p1)
dev.off()
