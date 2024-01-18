library(dplyr)
library(PreDeepMeta)
##以 test 数据集为例使用 PreDeepMeta 包产生输入数据，并使用模型进行预测

###cell gene expression
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",
                              data.table = F)
rownames(gene_exp) <- gene_exp$V1
gene_exp <- gene_exp %>% select(-V1)
colnames(gene_exp) <- gsub(" [(].+","",colnames(gene_exp))
gene_exp <- as.data.frame(t(gene_exp))

###cell info 
cell_mapping <- read.csv("/home/data/sdc/wt/update/data/Model.csv")
cell_mapping <- cell_mapping %>%
  filter((OncotreeLineage != "Normal") & (OncotreePrimaryDisease != "Non-Cancerous"))
net_cell_mapping <- data.frame(origin_net=NA,
                               cell=unique(cell_mapping$OncotreeLineage))
net_cell_mapping$origin_net <- c("iOvarianCancer1620.xml","bone_marrow.xml",
                                 "iColorectalCancer1750.xml","iSkinCancer1386.xml",
                                 "iUrothelialCancer1532.xml","iLungCancer1490.xml",
                                 "kidney.xml","iBreastCancer1746.xml",NA,
                                 "iPancreaticCancer1613.xml","brain.xml",NA,NA,
                                 "iStomachCancer1511.xml","iThyroidCancer1710.xml",NA,NA,
                                 "iProstateCancer1560.xml",NA,"iHeadNeckCancer1628.xml",
                                 "iEndometrialCancer1713.xml",NA,"iLiverCancer1788.xml",
                                 "iCervicalCancer1611.xml",NA,NA,"adrenal_gland.xml",
                                 NA,"iTestisCancer1483.xml",NA)
net_cell_mapping <- na.omit(net_cell_mapping)
net_cell_mapping$normal <- c("Ovary","Whole Blood","Colon - Sigmoid",
                             "Skin - Sun Exposed (Lower leg)","Bladder",
                             "Lung","Kidney - Cortex","Breast - Mammary Tissue",
                             "Pancreas","Brain - Amygdala","Stomach","Thyroid","Prostate",
                             "Whole Blood","Uterus","Liver","Cervix - Endocervix",
                             "Adrenal Gland","Testis")
cell_mapping <- cell_mapping %>%
  select(ModelID,OncotreeLineage) %>%
  rename(cell = OncotreeLineage) %>%
  left_join(.,net_cell_mapping)
cell_mapping <- na.omit(cell_mapping)
cell_mapping <- cell_mapping %>%
  mutate(net = paste0(gsub(".xml","",origin_net),"_enzymes_based_graph.tsv"))

####
enz_gene_mapping <- readRDS("~/DeepMeta/data/enz_gene_mapping.rds")
cpg_gene <- readRDS("~/DeepMeta/data/cpg_gene.rds")

###test cell 
test_cell <- read.csv("data/test_cell_info.csv")

##
library(doParallel)
library(foreach)
my.cluster <- parallel::makeCluster(
  40, 
  type = "PSOCK"
)
doParallel::registerDoParallel(cl = my.cluster)

foreach(
  i = test_cell$cell,
  .export = c("enz_gene_mapping","cpg_gene","cell_mapping","gene_exp"),
  .packages = c("dplyr","tidyr","PreDeepMeta")
) %dopar% {
  cell_net <- cell_mapping$net[which(cell_mapping$ModelID == i)]
  cell_net <- read.table(paste0("data/meta_net/EnzGraphs/",cell_net))
  PreDeepMeta::PreEnzymeNet(gene_exp, network = cell_net,
                            gene_mapping = enz_gene_mapping, gene_feature = cpg_gene,
                            cell_name = i, save_path = "/home/data/sdb/wt/model_data/enzyme_net_test/")
}
parallel::stopCluster(cl = my.cluster)

###expression
gtex <- data.table::fread("data/GTEx_Analysis_2017-06-05_v8_RNASeQCv1.1.9_gene_median_tpm.gct",
                          data.table = F,skip = 2) %>% 
  filter(!grepl("PAR",Name)) %>% 
  select(-Name)
data("model_gene_order")

cell_mapping <- cell_mapping %>% filter(ModelID %in% test_cell$cell)
gene_exp <- gene_exp %>% select(all_of(cell_mapping$ModelID))
PreDiffExp(tumor_exp = gene_exp, normal_exp = gtex,
           tumor_normal_mapping = cell_mapping,
           gene_order = model_gene_order,
           save_file = TRUE, 
           save_path = "/home/data/sdb/wt/model_data/test_diff_exp.csv")

