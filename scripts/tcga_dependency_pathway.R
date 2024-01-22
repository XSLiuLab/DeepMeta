library(dplyr)

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

pre <- readRDS("data/tcga_pre_V2.rds")
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()

# high_var <- pre_pathway %>% 
#   group_by(gene) %>% 
#   summarise(pos_per = mean(preds == 1)) %>% 
#   ungroup() %>% 
#   filter(pos_per > 0 & pos_per < 1)
####两两基因之间的代谢依赖性预测相关性
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
#saveRDS(res_dt,file = "data/tcga_pre_cor_all_gene.rds")
saveRDS(res_dt,file = "data/tcga_pre_cor_all_gene_rho.rds")
###将代谢网络中相关性高的基因提取并展示
##plot
# library(dplyr)
# library(igraph)
# library(SteinerNet)
# library(qgraph)
# get_cor_net <- function(tcga_cor,gene_name,met_dt,met_net){
#   cor_dt <- tcga_cor %>% 
#     select(gene_name) %>% 
#     filter_at(1, all_vars(. > 0.8))
#   cor_dt$gene <- rownames(cor_dt)
#   cor_dt <- cor_dt %>% filter(gene != gene_name)
#   
#   out <- steinertree(type = "KB", terminals = c(cor_dt$gene,gene_name), 
#                      graph = met_net, color = FALSE)
#   out_dt <- as_data_frame(out[[1]])
#   out_dt <- left_join(
#     out_dt %>% rename(from_gene = from,to_gene = to),
#     met_dt
#   )
#   dt <- met_dt
#   colnames(dt)[1:2] <- c("to_gene","from_gene")
#   out_dt <- left_join(
#     out_dt,
#     dt %>% rename(meta2 = Metabolite)
#   )
#   out_dt <- out_dt %>% 
#     mutate(Metabolite = ifelse(is.na(Metabolite),meta2,Metabolite)) %>% 
#     select(-meta2)
#   ##plot
#   out_g <- graph_from_data_frame(out_dt,directed = F)
#   V(out_g)[cor_dt$gene]$color <- "#11325D"
#   V(out_g)[gene_name]$color <- "#F5A673"
#   e <- get.edgelist(out_g,names=FALSE)
#   l <- qgraph.layout.fruchtermanreingold(e,vcount=vcount(out_g))
#   par(mar=c(0,0,0,0)+.1)
#   plot(out_g,vertex.size=8, edge.arrow.size=0.3,
#        vertex.label = ifelse(V(out_g)$name %in% c(cor_dt$gene,gene_name), 
#                              V(out_g)$name, NA),
#        vertex.label.dist=1.5,layout=l, vertex.label.cex = 0.8)
#   return(out_g)
# }
# 
# tcga_pre_cor <- readRDS("data/tcga_pre_cor.rds")
# HGM_all_gene <- readRDS("data/HGM_all_gene.rds")
# HGM_all_gene <- HGM_all_gene %>% 
#   select(from_gene, to_gene, Metabolite)
# met_g <- graph_from_data_frame(HGM_all_gene)
# 
# get_cor_net(tcga_cor = tcga_pre_cor,
#             gene_name = "FH",
#             met_dt = HGM_all_gene,
#             met_net = met_g)
# dev.copy2pdf(file = "figs/cor_net_pkm.pdf",
#              width = 8,height = 6,out.type = "pdf")

#####
###通路富集分析
###
kegg_enricher <- function(kegg_dt, gene_list, pathway_name){
  gene_list <- strsplit(gene_list,split = ",") %>% unlist()
  pathway_dt <- kegg_dt %>% 
    mutate(isin = ifelse(genes %in% gene_list,"yes","no")) %>% 
    mutate(pathway_type = ifelse(pathway == pathway_name,"yes","no")) %>% 
    select(genes,pathway_type,isin) %>% distinct_all() 
  in_pathway <- pathway_dt %>% filter(pathway_type == "yes")
  not_in_pathway <- pathway_dt %>% filter(pathway_type == "no")
  
  counts <- pathway_dt$genes[pathway_dt$isin == "yes" & 
                               pathway_dt$pathway_type == "yes"] %>% 
    length()
  fisher_dt <-  fisher.test(cbind(c(sum(in_pathway$isin == "yes"),
                                    sum(in_pathway$isin == "no")),
                                  c(sum(not_in_pathway$isin == "yes"),
                                    sum(not_in_pathway$isin == "no"))),
                            alternative = "greater")
  p <- fisher_dt$p.value
  return(list(p,counts)) 
}

tcga_pre_cor <- readRDS("data/tcga_pre_cor_all_gene_rho.rds")
all_genes_base <- data.frame(genes = rownames(tcga_pre_cor),
                             split_gene = rownames(tcga_pre_cor)) %>% 
  tidyr::separate_longer_delim(cols = "split_gene",delim = ",") %>% 
  rowwise() %>% 
  mutate(pathway = paste(kegg$pathway[grep(split_gene,kegg$genes)],collapse = ";")) %>% 
  ungroup() %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";") %>% 
  select(genes,pathway) %>% distinct_all() %>% 
  filter(nchar(pathway) > 1) %>% 
  rowwise() %>% 
  mutate(cor_counts = length(which(tcga_pre_cor[,genes] > 0.8) - 1)) %>% 
  ungroup() %>% filter(cor_counts >= 3)

all_genes <- all_genes_base %>% 
  rowwise() %>% 
  mutate(counts = kegg_enricher(kegg,
                                rownames(tcga_pre_cor)[which((tcga_pre_cor[,genes]>0.8) & (rownames(tcga_pre_cor) != genes))],
                                pathway)[[2]]) %>% 
  mutate(p = kegg_enricher(kegg,
                           rownames(tcga_pre_cor)[which((tcga_pre_cor[,genes]>0.8) & (rownames(tcga_pre_cor) != genes))],
                           pathway)[[1]]) %>% ungroup()
all_genes$fdr <- p.adjust(all_genes$p,"fdr")
all_genes_summ <- all_genes %>% 
  group_by(genes) %>% 
  summarise(pathway_min = pathway[which.min(p)],
            min_p = p[which.min(p)],
            fdr_p = fdr[which.min(p)],
            cor_counts = unique(cor_counts),
            counts_p = counts[which.min(p)]) %>% ungroup()
mean(all_genes_summ$min_p < 0.05)
mean(all_genes_summ$fdr_p < 0.1)
#saveRDS(all_genes_summ,file = "data/tcga_codepency_fisher_all_genes.rds")
saveRDS(all_genes_summ,file = "data/tcga_codepency_fisher_all_genes_rho.rds")

####基于所有基因做置换检验
tcga_pre_cor <- readRDS("data/tcga_pre_cor_all_gene_rho.rds")
all_genes_base <- data.frame(genes = rownames(tcga_pre_cor),
                             split_gene = rownames(tcga_pre_cor)) %>% 
  tidyr::separate_longer_delim(cols = "split_gene",delim = ",") %>% 
  rowwise() %>% 
  mutate(pathway = paste(kegg$pathway[grep(split_gene,kegg$genes)],
                         collapse = ";")) %>% 
  ungroup() %>% 
  tidyr::separate_longer_delim(cols = "pathway",delim = ";") %>% 
  select(genes,pathway) %>% distinct_all() %>% 
  filter(nchar(pathway) > 1) %>% 
  rowwise() %>% 
  mutate(cor_counts = length(which(tcga_pre_cor[,genes]>0.8)) - 1) %>% 
  ungroup() %>% filter(cor_counts >= 3)

library(doParallel)
library(foreach)
#create the cluster
my.cluster <- parallel::makeCluster(
  60, 
  type = "PSOCK"
)
#register it to be used by %dopar%
doParallel::registerDoParallel(cl = my.cluster)

res <- foreach(
  i = 1:1000,
  .packages = c("dplyr","tidyr")
) %dopar% {
  all_genes_sample <- all_genes_base %>% 
    rowwise() %>% 
    mutate(p = kegg_enricher(kegg,
                             sample(rownames(tcga_pre_cor)[which(rownames(tcga_pre_cor) != genes)],
                                    cor_counts,
                                    replace=FALSE),
                             pathway)[[1]]) %>% ungroup()
  all_genes_sample$fdr <- p.adjust(all_genes_sample$p,"fdr")
  all_genes_sample_summ <- all_genes_sample %>% 
    group_by(genes) %>% 
    summarise(pathway_min = pathway[which.min(p)],
              min_p = p[which.min(p)],
              fdr_p = fdr[which.min(p)]) %>% ungroup()
  all_genes_sample_summ$sample_id <- i
  return(all_genes_sample_summ)
}
parallel::stopCluster(cl = my.cluster)

res <- bind_rows(res)
#saveRDS(res,file = "data/pathway_sim_all_genes.rds")
saveRDS(res,file = "data/pathway_sim_all_genes_rho.rds")

###plot
library(ggprism)
res <- readRDS("data/pathway_sim_all_genes_rho.rds")
sampling_res <- res %>% 
  group_by(sample_id) %>% 
  summarise(p_sig = mean(min_p < 0.05),
            fdr_sig = mean(fdr_p < 0.1)) %>% 
  ungroup()
all_genes_summ <- readRDS("data/tcga_codepency_fisher_all_genes_rho.rds")
p <- WVPlots::ShadedDensity(frame = sampling_res, 
                            xvar = "fdr_sig",
                            threshold = mean(all_genes_summ$fdr_p < 0.1),
                            title = "",
                            tail = "right",linecolor="red")
p$layers[[1]]$aes_params$size <- 1
p$layers[[2]]$aes_params$fill <- "blue"     #geom_ribbon
p$layers[[3]]$aes_params$colour <- "black" 
p$layers[[3]]$aes_params$size <- 1
p1 <- p + labs(x="Significant proportion in the simulation (FDR < 0.1)",
               y="Density")+
  theme_prism()

p2 <- p + labs(x="Significant proportion in the simulation (P < 0.05)",
               y="Density")+
  theme_prism()
p1
ggsave("figs/cor_density_sig_fdr.pdf",width = 7,height = 5)
p2
ggsave("figs/cor_density_sig_p.pdf",width = 7,height = 5)

# median(sampling_res$p_sig)
# [1] 0.1428571
# median(sampling_res$fdr_sig)
# [1] 0.05167173

####画图
tcga_pre_cor <- readRDS("data/tcga_pre_cor_all_gene_rho.rds")
tcga_pre_cor$main_gene <- rownames(tcga_pre_cor)
tcga_pre_cor <- tcga_pre_cor %>% 
  tidyr::pivot_longer(cols = ACSL4:PGAM4, names_to = "target_gene",
                      values_to = "cor") %>% 
  rowwise() %>% 
  mutate(main_pathway = get_pathway(main_gene,kegg)[1],
         target_pathway = get_pathway(target_gene,kegg)[1]) %>% 
  ungroup() %>% 
  filter(nchar(main_pathway)>1) %>% 
  filter(nchar(target_pathway)>1)
#saveRDS(tcga_pre_cor, file = "data/gene_cor_add_pathway.rds")
saveRDS(tcga_pre_cor, file = "data/gene_cor_add_pathway_rho.rds")

####
tcga_pre_cor <- readRDS("data/gene_cor_add_pathway_rho.rds")
tcga_pre_cor <- tcga_pre_cor %>% 
  rowwise() %>% 
  mutate(pathway_overlap = 
           ifelse(length(intersect(strsplit(main_pathway,";")[[1]],strsplit(target_pathway,";")[[1]])) > 0,"yes","no")) %>% 
  ungroup()

tcga_pre_cor_summ <- tcga_pre_cor %>% 
  na.omit() %>% 
  filter(main_gene != target_gene) %>% 
  group_by(main_gene, pathway_overlap) %>% 
  summarise(median_cor = median(cor)) %>% ungroup() 
dt <- tcga_pre_cor_summ %>% 
  tidyr::pivot_wider(names_from = pathway_overlap, values_from = median_cor) %>% 
  na.omit()
dt <- dt %>% 
  mutate(diff = log2(abs(yes)/abs(no)))

library(ggpubr)
ggdensity(dt, x = "diff", fill="#69b3a2",
          xlab = "log2(In pathway / Out pathway)", ylab = "Density")+
  geom_vline(xintercept = 0, color = "red", linewidth=1.5)
ggsave("figs/cor_diff_all_genes.pdf",width = 5,height = 3)
# mean(dt$diff > 0)
# [1] 0.5751634

###以 ACLY 为例
tcga_pre_cor <- readRDS("data/gene_cor_add_pathway_rho.rds")
tca_acly <- kegg %>% 
  filter(pathway == "Citrate cycle (TCA cycle)")

acly <- tcga_pre_cor %>% 
  tidyr::separate_longer_delim(cols = "target_pathway", delim = ";") %>% 
  filter(main_gene == "ACLY") %>% 
  filter(target_gene != "ACLY") %>% 
  mutate(type = ifelse(grepl("Citrate cycle \\(TCA cycle\\)",target_pathway),
                       "TCA cycle pathway","Other pathway")) %>% 
  select(1,2,type,cor) %>% 
  distinct_all()
ggboxplot(data=acly,x="type",y="cor",
          xlab = FALSE, ylab = "Spearman Correlation")+
  stat_compare_means()
ggsave("figs/acly_cor_diff.pdf",width = 4,height = 4)

###绘制一个pathway的图
tcga_pre_cor <- readRDS("data/gene_cor_add_pathway_rho.rds")
tca_acly <- kegg %>% 
  filter(pathway == "Citrate cycle (TCA cycle)")
HGM_all_gene <- readRDS("data/HGM_all_gene.rds")
HGM_all_gene <- HGM_all_gene %>%
  select(from_gene, to_gene, Metabolite)
if_in <- function(gene,dt){
  sum(sapply(strsplit(gene,split = ",")[[1]],
             function(x){length(grep(x,dt$genes))}))
}
HGM_all_gene <- HGM_all_gene %>% 
  filter(nchar(from_gene)>0 & nchar(to_gene)>0) %>% 
  rowwise() %>% 
  mutate(from_tca = ifelse(if_in(from_gene,tca_acly)!=0,"yes","no"),
         to_tca = ifelse(if_in(to_gene,tca_acly)!=0,"yes","no")) %>% 
  ungroup()
HGM_tca <- HGM_all_gene %>% 
  filter(from_tca == "yes"  & to_tca == "yes")

acly_cor <- tcga_pre_cor %>% 
  tidyr::separate_longer_delim(cols = "target_pathway", delim = ";") %>% 
  filter(main_gene == "ACLY") %>% 
  filter(target_gene != "ACLY") %>% 
  mutate(type = ifelse(target_pathway == "Citrate cycle (TCA cycle)",
                       "TCA cycle pathway","Other pathway")) %>% 
  filter(type == "TCA cycle pathway")

library(igraph)
library(SteinerNet)
library(qgraph)

met_g <- graph_from_data_frame(HGM_tca %>% select(1,2))
V(met_g)$color <- "white"
V(met_g)[unique(acly_cor$main_gene)]$color <- "#11325D"
V(met_g)[intersect(acly_cor$target_gene,V(met_g)$name)]$color <- "#F5A673"


plot(delete.vertices(simplify(met_g), degree(simplify(met_g))==0),
     vertex.size=8, edge.arrow.size=0.3,
     vertex.label = ifelse(V(met_g)$name %in% c(acly_cor$main_gene,
                                               acly_cor$target_gene),
                          V(met_g)$name, NA),
    vertex.label.dist=1.5, vertex.label.cex = 0.8,
    layout=layout.circle)

sim_g <- igraph::delete.vertices(igraph::simplify(met_g),
                                 igraph::degree(igraph::simplify(met_g))==0)
sim_g <- igraph::as_data_frame(sim_g)
library(ggnetwork) 
library(network)
ggnet <- network(sim_g, loops = T)
dt <- data.frame(
  nodes = ggnet %v% "vertex.names"
) %>% mutate(color = case_when(
  nodes %in% acly_cor$main_gene ~ "main",
  TRUE ~ NA
)) %>% 
  left_join(.,
           acly_cor %>% select(target_gene,cor) %>% 
             rename(nodes = target_gene)
  ) %>% 
  mutate(label = ifelse(nodes == "ACLY","ACLY",NA))
dt$cor[which(dt$nodes == "ACLY")] <- 1
ggnet %v% "color" <- dt$color
ggnet %v% "Cor" <- dt$cor
ggnet %v% "name" <- dt$label

ggplot(ggnet, aes(x = x, y = y, xend = xend, yend = yend)) +
  geom_edges(arrow = arrow(length = unit(6, "pt"), type = "closed")) +
  geom_nodes(aes(fill = Cor), size = 4, pch=21) +
  geom_nodelabel_repel(aes(label = name))+
  theme_blank()+
  scale_fill_gradient(low = "white",
                      high = "red")
ggsave("figs/acly_pathway.pdf",width = 6, height = 5)
