library(dplyr)
library(parallel)
library(caret)

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

pre <- readRDS("data/tcga_pre_V2.rds") %>% as.data.frame()
pre_pathway <- data.frame(gene = unique(c(pre$gene))) %>% 
  rowwise() %>% 
  mutate(pathway = get_pathway(gene,kegg)[1],
         class = get_pathway(gene,kegg)[2]) %>% 
  ungroup() %>% 
  filter(nchar(pathway)>1)
pre_pathway <- pre %>% 
  left_join(.,pre_pathway) %>% na.omit()

gene_screen <- function(pre_dt){
  gene_summ <- pre_dt %>% 
    group_by(gene) %>% 
    summarise(counts=mean(preds == 1)) %>%
    ungroup() %>% 
    filter(counts > 0 & counts < 1) 
  gene_summ2 <- pre_dt %>% 
    group_by(gene) %>% 
    summarise(counts=mean(type == "mut"))  %>% 
    ungroup() %>% 
    filter(counts > 0 & counts < 1) 
  all_genes <- intersect(gene_summ$gene,gene_summ2$gene)
  
  pre_dt <- pre_dt %>% filter(gene %in% all_genes)
  
  res <- vector("list",length = length(all_genes))
  for (i in 1:length(all_genes)){
    dt <- pre_dt %>% filter(gene == all_genes[i])
    dt$type <- factor(dt$type,levels = c("non-mut","mut"))
    dt_cancer <- dt %>% select(cancer)
    dummy <- caret::dummyVars(" ~ .", data=dt_cancer)
    dt_cancer <- data.frame(predict(dummy, newdata=dt_cancer))
    
    dt <- bind_cols(dt,dt_cancer)
    all_caners <- colnames(dt_cancer)
    dt_model <- glm(data = dt, 
                    as.formula(paste0("preds ~ type + ",
                                      paste(all_caners,collapse = "+"))), 
                    family =binomial(link = "logit"))
    model_summ <- jtools::summ(dt_model,exp=TRUE)$coeftable %>% as.data.frame()
    model_summ$var <- rownames(model_summ)
    model_summ <- model_summ %>% filter(var == "typemut")
    model_summ$gene <- all_genes[i]
    res[[i]] <- model_summ
  }
  res <- bind_rows(res)
  colnames(res)[1:3] <- c("OR","lower","upper")
  colnames(res)[5] <- c("P")
  colnames(res)[7] <- c("Gene")
  res <- res %>% filter(P<0.05 & OR>1)
  res$mean <- res$OR
  #res$OR <- round(res$OR,3)
  #res$P <- round(res$P,3)
  res <- res %>% 
    arrange(P,desc(OR))
  return(res)
}

##
mut <- readRDS("/home/data/sdc/wt/TCGA/tcga_mut.rds")

####ctnnb
need_gene <- "CTNNB1"
mut_pos <- c("p.S45F","p.S45P","p.S45Y","p.S45del","p.S45_P52del",
             "p.K335I","p.K335T","p.N387K","p.T41A","p.T41I",
             "p.D32A","p.D32Y","p.D32G","p.D32H","p.D32V","p.D32N",
             "p.S33A","p.S33C","p.S33F","p.S33Y","p.S33P",
             "p.G34V","p.G34E","p.G34R",
             "p.S37A","p.S37F","p.S37Y","p.S37C")
mut_gene <- mut %>% 
  filter(gene == need_gene & Amino_Acid_Change %in% mut_pos)
all_mut_gene <- mut %>% 
  filter(gene == need_gene)

ctnnb_gene_pre <- pre_pathway %>% 
  mutate(type = case_when(
    cell %in% mut_gene$sample ~ "mut",
    !(cell %in% all_mut_gene$sample) ~ "non-mut",
    TRUE ~ "other"
  )) %>% filter(type != "other") 
ctnnb_res <- gene_screen(pre_dt = ctnnb_gene_pre)

###MYC
gistic <- data.table::fread("/home/data/sdc/wt/TCGA/GISTIC.focal_data_by_genes.conf_95.txt.gz",data.table = F)
gistic_myc <- gistic %>% 
  filter(`Gene Symbol` %in% c("MYC","MYCL1","MYCN")) %>% 
  select(-c(2,3)) %>% as.data.frame()
rownames(gistic_myc) <- gistic_myc$`Gene Symbol`
gistic_myc <- gistic_myc %>% 
  select(-`Gene Symbol`) %>% 
  t() %>% 
  as.data.frame()
gistic_myc$sample <- rownames(gistic_myc)
gistic_myc <- gistic_myc %>% 
  rowwise() %>% 
  mutate(type = ifelse(any(c(MYC,MYCL1,MYCN)>0),"amp","non-amp")) %>% 
  ungroup()
gistic_myc <- gistic_myc %>% 
  mutate(sample = substr(sample,1,15))

myc_amp <- gistic_myc %>% filter(type == "amp")
all_mut_gene <- mut %>% 
  filter(gene %in% c("MYC","MYCL1","MYCN"))
all_mut_sample <- unique(c(all_mut_gene$sample,myc_amp$sample))

myc_gene_pre <- pre_pathway %>% 
  mutate(type = case_when(
    cell %in% myc_amp$sample ~ "mut",
    !(cell %in% all_mut_sample) ~ "non-mut",
    TRUE ~ "other"
  )) %>% filter(type != "other") 
myc_res <- gene_screen(pre_dt = myc_gene_pre)

###tp53
tp53_mut <- mut %>% 
  filter(gene == "TP53") %>% 
  filter(effect %in% c("Frame_Shift_Del","Frame_Shift_Ins","In_Frame_Del",
                       "In_Frame_Ins","Missense_Mutation","Nonsense_Mutation",
                       "Splice_Site") | Amino_Acid_Change == "p.T125T") %>% 
  mutate(mut_index = paste(chr,start,end,reference,alt,sep = "-"))
tp53_mut_ckb <- readRDS("data/tp53_mut_ckb.rds")
tp53_mut_ckb <- paste0("p.",tp53_mut_ckb$mut)
tp53_mut_ckb <- tp53_mut_ckb[which(tp53_mut_ckb %in% tp53_mut$Amino_Acid_Change)]

need_gene <- "TP53"
mut_pos <- tp53_mut_ckb
mut_gene <- mut %>% 
  filter(gene == need_gene & Amino_Acid_Change %in% mut_pos)
all_mut_gene <- mut %>% 
  filter(gene == need_gene)
tp53_gene_pre <- pre_pathway %>% 
  mutate(type = case_when(
    cell %in% mut_gene$sample ~ "mut",
    !(cell %in% all_mut_gene$sample) ~ "non-mut",
    TRUE ~ "other"
  )) %>% filter(type != "other") 
tp53_res <- gene_screen(pre_dt = tp53_gene_pre)

###汇总
library(forestplot)
all_res <- bind_rows(
  myc_res %>% mutate(target = "MYC"),
  tp53_res %>% mutate(target = "TP53"),
  ctnnb_res %>% mutate(target = "CTNNB1")
)
saveRDS(all_res,file = "data/merge_driver_dep2.rds")

gene_pre_res <- list(
  MYC = myc_gene_pre,
  TP53 = tp53_gene_pre,
  CTNNB = ctnnb_gene_pre
)
saveRDS(gene_pre_res,file = "/home/data/sdc/wt/TCGA/tcga_driver_pre.rds")

rm(list = ls())

all_res <- readRDS("data/merge_driver_dep.rds")
all_res_filter <- all_res %>% 
  mutate(target = ifelse(target == "CTNNB", "CTNNB1", target)) %>% 
  filter(!is.infinite(OR)) %>% 
  mutate(`Target-Gene` = paste(target,Gene,sep = "-")) %>% 
  filter(target %in% c("TP53","KRAS","MYC","CTNNB1")) %>% 
  mutate(`Target-Gene` = ifelse(
    `Target-Gene` %in% c("TP53-GPX4","CTNNB1-COASY","CTNNB1-IMPDH2",
                         "KRAS-GFPT1","KRAS-PKM","MYC-CAD","MYC-TYMS"),
    paste0(`Target-Gene`," *"), `Target-Gene`
  ))
all_res_filter <- all_res_filter %>% filter(target != "KRAS")
p <- all_res_filter |>
  forestplot(labeltext = c(`Target-Gene`, OR, P),
             xlog = TRUE,
             xlab = "OR",boxsize = 0.25,) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(`Target-Gene` = c("Target-Gene"),
                OR = c("OR"),
                P = c("P")) |>
  fp_set_zebra_style("#EFEFEF")
p
library(gridExtra)
pdf(file = "figs/merge_OR.pdf", onefile=TRUE, width = 10,height = 9)
grid.newpage()
print(p)
dev.off()

##only ctnnb1
all_res <- readRDS("data/merge_driver_dep.rds")
all_res_filter <- all_res %>% 
  filter(!is.infinite(OR)) %>% 
  filter(target %in% c("CTNNB")) %>% 
  mutate(target = "CTNNB1") %>% 
  mutate(`Target-Gene` = paste(target,Gene,sep = "-"))
p <- all_res_filter |>
  forestplot(labeltext = c(`Target-Gene`, OR, P),
             xlog = TRUE,
             xlab = "OR",boxsize = 0.25,) |>
  fp_set_style(box = "royalblue",
               line = "darkblue",
               summary = "royalblue") |> 
  fp_add_header(`Target-Gene` = c("Target-Gene"),
                OR = c("OR"),
                P = c("P")) |>
  fp_set_zebra_style("#EFEFEF")
library(gridExtra)
pdf(file = "figs/merge_OR_ctnnb1.pdf", onefile=TRUE, width = 9,height = 5)
grid.newpage()
print(p)
dev.off()
#####单个癌症类型的汇总
gene_screen_one_cancer_fisher <- function(pre_dt,need_gene){
  dt <- pre_dt %>% filter(gene == need_gene)
  dt$type <- factor(dt$type,levels = c("mut","non-mut"))
  dt$preds <- factor(dt$preds, levels = c(1,0))
  
  all_cancers <- unique(dt$cancer)
  cancer_res <- vector("list",length(all_cancers))
  names(cancer_res) <- all_cancers
  for (i in 1:length(all_cancers)){
    dt_cancer <- dt %>% filter(cancer == all_cancers[i])
    if (length(unique(as.character(dt_cancer$preds))) == 2 & 
        length(unique(as.character(dt_cancer$type))) == 2){
      dt_fisher <- fisher.test(table(dt_cancer$preds,dt_cancer$type),
                               alternative = "greater")
      
      dt_cancer_summ <- dt_cancer %>% 
        group_by(cancer,type,preds) %>% 
        summarise(counts=n()) %>% ungroup()
      dt_cancer_summ$p <- dt_fisher$p.value
      
      cancer_res[[i]] <- dt_cancer_summ
    }else{
      cancer_res[[i]] <- NA
    }
  }
  cancer_res <- bind_rows(cancer_res[which(lengths(cancer_res)>1)])
  
  return(cancer_res)
}

get_gene_split_cancer <- function(res_dt,pre_dt){
  all_genes <- unique(res_dt$Gene)
  all_genes_res <- vector("list",length(all_genes))
  for (i in 1:length(all_genes)){
    res_cancer <- gene_screen_one_cancer_fisher(pre_dt,all_genes[i])
    res_cancer$gene <- all_genes[i]
    all_genes_res[[i]] <- res_cancer
  }
  all_genes_res <- bind_rows(all_genes_res)
  return(all_genes_res)
}

all_res <- readRDS("data/merge_driver_dep.rds")
gene_pre <- readRDS("/home/data/sdc/wt/TCGA/tcga_driver_pre.rds")
cancer_res <- vector("list",4)
all_genes <- c("TP53","KRAS","MYC","CTNNB")
for (i in 1:length(all_genes)){
  dt <- all_res %>% filter(target == all_genes[i])
  dt_gene_pre <- gene_pre[[all_genes[i]]]
  dt_res <- get_gene_split_cancer(dt, dt_gene_pre)
  dt_res <- dt_res %>% mutate(target = all_genes[i])
  cancer_res[[i]] <- dt_res
}
cancer_res <- bind_rows(cancer_res)
cancer_res <- cancer_res %>%  filter(p < 0.05)
saveRDS(cancer_res,file = "data/split_cancer_res_greater.rds")


###按癌症类型画图
all_res <- readRDS("data/split_cancer_res_greater.rds")
all_res_summ <- all_res %>% 
  group_by(target,gene,cancer,type) %>% 
  summarise(percent_1 = counts[which(preds == 1)]/sum(counts),
            p_value = -log10(unique(p)),
            all_counts = sum(counts)) %>% ungroup()
##只选择 mut 里面大于 non-mut 里面的
all_res_summ2 <- all_res_summ %>% 
  group_by(target,gene,cancer) %>% 
  summarise(large=ifelse(percent_1[which(type == "mut")] > percent_1[which(type == "non-mut")],1,0)) %>% 
  ungroup() %>% filter(large == 1)

all_res_summ <- left_join(
  all_res_summ,
  all_res_summ2
)

all_res_summ <- all_res_summ %>% filter(!is.na(large))
all_res_summ <- all_res_summ %>% 
  select(-large) %>% 
  tidyr::pivot_wider(names_from = "type",
                     values_from = c("percent_1","all_counts"))
all_res_summ <- all_res_summ %>% 
  mutate(index = paste(gene,cancer,sep = "-"))

###normal
tcga_normal <- readRDS("data/tcga_normal_pre_V2.rds") %>% as.data.frame()
normal_summ <- tcga_normal %>% 
  group_by(gene,cancer) %>% 
  summarise(counts = n(),
            pos_percent = mean(preds == 1)) %>% ungroup() %>% 
  mutate(index = paste(gene,cancer,sep = "-")) %>% 
  filter(index %in% all_res_summ$index)

all_res_summ <- left_join(
  all_res_summ,
  normal_summ %>% select(-cancer,-gene)
)

# dt <- all_res_summ %>%  ##至少突变3个样本
#   filter(all_counts_mut >=3) %>% 
#   mutate(ratio = round(percent_1_mut / `percent_1_non-mut`,2),
#          `Cancer Sample Counts\n (Mut/Non-mut/Normal)` = 
#            paste(all_counts_mut,`all_counts_non-mut`,counts,sep = "/"),
#          `Positive Percentage\n (Mut/Non-mut/Normal)` = 
#            paste(round(percent_1_mut,2),round(`percent_1_non-mut`,2),
#                  round(pos_percent,2),sep = "/"),
#          p_value = round(p_value,2)) %>% 
#   select(target, gene, cancer, ratio, p_value, c(13,14))
#colnames(dt)[1:3] <- c("Target Gene","Gene","Cancer")
#saveRDS(dt,file = "data/split_cancer_webshot.rds")

dt <- all_res_summ %>%  ##至少突变3个样本
  filter(all_counts_mut >=3) %>% 
  mutate(ratio1 = round(percent_1_mut / `percent_1_non-mut`,2),
         ratio2 = round((percent_1_mut - pos_percent) / percent_1_mut,2),
         p_value = round(p_value,2)) %>% 
  select(target, gene, cancer, p_value, c(12,13))

colnames(dt)[1:3] <- c("Target Gene","Gene","Cancer")
saveRDS(dt,file = "data/split_cancer_webshot2.rds")

library(htmltools)
library(reactable)
dt <- readRDS("data/split_cancer_webshot2.rds") 
dt <- dt %>% 
  filter(`Target Gene` == "CTNNB") %>% 
  mutate(`Target Gene` = "CTNNB1")
dt$Gene <- factor(dt$Gene, levels = c("IMPDH2","TYMS","PDPK1","UAP1","ACACA","OGDH","COASY","SLC2A1"))
dt <- dt %>% arrange(Gene)
# dt <- dt %>% 
#   mutate(`Target Gene` = ifelse(`Target Gene` == "CTNNB", "CTNNB1", `Target Gene`))
# Render a bar chart with a label on the left
bar_chart <- function(label, width = "100%", height = "1rem", fill = "#00bfc4", background = NULL) {
  bar <- div(style = list(background = fill, width = width, height = height))
  chart <- div(style = list(flexGrow = 1, marginLeft = "0.5rem"), bar)
  div(style = list(display = "flex", alignItems = "center"), label, chart)
}

bar_chart_pos_neg <- function(label, width = "100%", height = "1rem", fill = "#00bfc4", background = NULL) {
  neg_chart <- div(style = list(flex = "1 1 0"))
  pos_chart <- div(style = list(flex = "1 1 0"))
  if (label < 0){
    bar <- div(style = list(background = fill, width = width, height = height, marginLeft = "0.8rem"))
    chart <- div(style = list(flexGrow = 1, display = "flex", alignItems = "center", justifyContent = "flex-end"), bar)
    neg_chart <- tagAppendChild(neg_chart, chart)
  } else {
    bar <- div(style = list(background = fill, width = width, height = height, marginRight = "0.8rem"))
    chart <- div(style = list(flexGrow = 1, display = "flex", alignItems = "center"), bar)
    pos_chart <- tagAppendChild(pos_chart, chart)
  }
  div(style = list(display = "flex"), label, neg_chart, pos_chart)
}

tt <- reactable(
  dt,
  columns = list(
    ratio1 = colDef(name = "Ratio <br> (MP / NMP)", align = "left", html = TRUE, cell = function(value) {
      width <- paste0(abs(value / max(dt$ratio1)) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart(label, width = width, fill = "#11325D", background = "#e1e1e1")
    },minWidth = 100),
    ratio2 = colDef(name = "Ratio <br> ([MP - NP] / MP)", align = "left", html = TRUE, cell = function(value) {
      width <- paste0(abs(value / max(dt$ratio2)) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart_pos_neg(label, width = width, fill = "#f5a673", background = "#e1e1e1")
    },minWidth = 120),
    p_value = colDef(name = "P Value (-log10)", align = "left", cell = function(value) {
      width <- paste0(value / max(dt$p_value) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart(label, width = width, fill = "#B783AF", background = "#e1e1e1")
    },minWidth = 100)
  ),
  defaultPageSize = 40
)
html <- "figs/dt_table_ctnnb1_add_normal.html"
library(htmlwidgets)
saveWidget(tt, html)
library(webshot2)
webshot(html, "figs/dt_table_ctnnb1_add_normal.png",vwidth=1000,zoom=2)

psycModel::html_to_pdf(file_path = html)

###
dt <- readRDS("data/split_cancer_webshot.rds") %>% 
  mutate(`Target Gene` = ifelse(`Target Gene` == "CTNNB", "CTNNB1", `Target Gene`))
dt <- dt %>% filter(`Target Gene` %in% c("TP53","MYC","CTNNB1"))

# Render a bar chart with a label on the left
bar_chart <- function(label, width = "100%", height = "1rem", fill = "#00bfc4", background = NULL) {
  bar <- div(style = list(background = fill, width = width, height = height))
  chart <- div(style = list(flexGrow = 1, marginLeft = "0.5rem", background = background), bar)
  div(style = list(display = "flex", alignItems = "center"), label, chart)
}
tt <- reactable(
  dt,
  columns = list(
    ratio = colDef(name = "Ratio", align = "left", cell = function(value) {
      width <- paste0(value / max(dt$ratio) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart(label, width = width, fill = "#11325D", background = "#e1e1e1")
    }),
    p_value = colDef(name = "P Value (-log10)", align = "left", cell = function(value) {
      width <- paste0(value / max(dt$p_value) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart(label, width = width, fill = "#B783AF", background = "#e1e1e1")
    })
  ),
  defaultPageSize = 40
)
html <- "figs/dt_table.html"
library(htmlwidgets)
saveWidget(tt, html)
library(webshot2)
webshot(html, "figs/dt_table.png",vwidth=1500,zoom=4)

###
dt <- readRDS("data/split_cancer_webshot2.rds") %>% 
  mutate(`Target Gene` = ifelse(`Target Gene` == "CTNNB", "CTNNB1", `Target Gene`))
dt <- dt %>% filter(`Target Gene` %in% c("TP53","MYC"))

dt_summ <- dt %>% 
  group_by(`Target Gene`,Gene) %>% 
  summarise(p = max(p_value),
            r1 = max(ratio1),
            r2 = max(ratio2)) %>% 
  ungroup() %>% arrange(desc(p),desc(r1),desc(r2))
dt_summ$`Target Gene` <- factor(dt_summ$`Target Gene`, levels = c("TP53","MYC"))
dt_summ <- dt_summ %>% arrange(desc(`Target Gene`))

dt$`Target Gene` <- factor(dt$`Target Gene`, levels = c("TP53","MYC"))
dt$Gene <- factor(dt$Gene,levels = dt_summ$Gene)
dt <- dt %>% arrange(desc(`Target Gene`),Gene)
tt <- reactable(
  dt,
  columns = list(
    ratio1 = colDef(name = "Ratio <br> (MP / NMP)", align = "left", html = TRUE, cell = function(value) {
      width <- paste0(abs(value / max(dt$ratio1)) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart(label, width = width, fill = "#11325D", background = "#e1e1e1")
    },minWidth = 100),
    ratio2 = colDef(name = "Ratio <br> ([MP - NP] / MP)", align = "left", html = TRUE, cell = function(value) {
      width <- paste0(abs(value / max(dt$ratio2)) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart_pos_neg(label, width = width, fill = "#f5a673", background = "#e1e1e1")
    },minWidth = 120),
    p_value = colDef(name = "P Value (-log10)", align = "left", cell = function(value) {
      width <- paste0(value / max(dt$p_value) * 100, "%")
      label <- format(value, nsmall = 2)
      bar_chart(label, width = width, fill = "#B783AF", background = "#e1e1e1")
    },minWidth = 100)
  ),
  defaultPageSize = 40
)

html <- "figs/dt_table_ratio_without_ctnnb.html"
library(htmlwidgets)
saveWidget(tt, html)
library(webshot2)
webshot(html, "figs/dt_table_ratio_without_ctnnb.png",vwidth=1400,zoom=2)
psycModel::html_to_pdf(file_path = html, scale = 0.5)

