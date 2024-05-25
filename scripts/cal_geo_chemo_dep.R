library(dplyr)
###计算通路富集 OR
pre <- readRDS("~/DeepMeta/data/geo_chemo_pre.rds")
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
all_pathways <- unique(all_split$pathway)
all_res <- vector("list",length(all_pathways))
for (i in 1:length(all_pathways)){
  tmp <- all_pathways[i]
  my.cluster <- parallel::makeCluster(
    60, 
    type = "PSOCK"
  )
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
saveRDS(res,file = "~/DeepMeta/data/geo_chemo_meta_dep.rds")

