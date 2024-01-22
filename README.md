### Predicting metabolic dependency using DeepMeta

We used `PreDeepMeta`​ package to prepare input data needed for DeepMeta:

```shell
###https://github.com/wt12318/PreDeepMeta
##devtools::install_github("wt12318/PreDeepMeta")
library(dplyr)
library(PreDeepMeta)
```

Here, we present the **DeepMeta** framework to predict the metabolic gene dependence based on the metabolic pathway information characterized by enzyme networks and sample status information defined by gene expression profile. Thus, DeepMeta has two inputs, the sample-specific enzyme network and the gene expression information (**Figure 1a**). We used the graph attention network (GAT) module to extract information from the sample-specific enzyme network to obtain the embedding of the metabolic enzyme gene node, and used the fully connected neuron network module to extract information from the expression profiles of cells to get the sample state embedding. Here, the features of gene node in enzyme network was implemented by binary vector which denotes the involvement of the gene in CPG (chemical and genetic perturbation) signature. These two parts of information are then merged to predict metabolic dependencies

![](https://picgo-wutao.oss-cn-shanghai.aliyuncs.com/image-20240118105406-fuurffy.png)​

Two functions `PreEnzymeNet`​ and `PreDiffExp`​ can be used for preparing DeepMeta inputs.

```R
###cell gene expression
gene_exp <- data.table::fread("/home/data/sdb/wt/model_data/OmicsExpressionProteinCodingGenesTPMLogp1.csv",data.table = F)
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

####gene mapping and CPG features, The code that generates this data is in `scripts/help_data.R`
enz_gene_mapping <- readRDS("~/DeepMeta/data/enz_gene_mapping.rds")
cpg_gene <- readRDS("~/DeepMeta/data/cpg_gene.rds")
```

We used 76 test cell lines as the example:

```R
###test cell 
test_cell <- read.csv("data/test_cell_info.csv")

##enzyme net input
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

###expression input
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

```

Then we can use python script `pred_enzyme.py`​ (which in `scripts/model`​ fold) to predict metabolic dependency:

```R
python ~/DeepMeta/scripts/model/pred_enzyme.py 
  -e /home/data/sdb/wt/model_data/test_diff_exp.csv 
  -g /home/data/sdb/wt/model_data/tmp/example_test/ 
  -c /home/wt/DeepMeta/data/test_cell_info.csv 
  -n /home/data/sdb/wt/model_data/enzyme_net_test/ 
  -t 30 
  -m /home/data/sdc/wt/model_data/new_model/enzyme_model_filterV2.pt 
  -o /home/wt/DeepMeta/data/example_test.csv 
  -d val -b 1
```

The arguments are :

* `-e`​ : file path of gene expression (tumor VS normal), output of `PreDiffExp`​
* `-g`: dir path of pyg data, the `processed`​ and `raw`​ dir will be created to save PyG processed graph data and raw sample information data
* `-c`: file name of predicted cell (csv), with columns of `cell`​ and `cell_index`​
* `-n`: dir path of sample enzyme net and gene node features, output of `PreEnzymeNet`​
* `-t`: multiprocess cores for creating pyg data
* `-m`: file path of model `.pt`​ file, users can download the DeepMeta model from [google drive](https://drive.google.com/file/d/1ZQAaSeOgmgBy-dE23i5qu5pieaAdCCUd/view?usp=drive_link)  
* `-o`: output file name
* `-d`: datatype, test or val; test mean you have the dependency data, while val mean not have the dependency data
* `-b`: batch size

The output is the csv file with predicted dependency probability (`preds_raw`​ column) and lable (using cutoff probability 0.5, `preds`​ column):

```R
dt <- read.csv("data/example_test.csv") %>% select(-X)
View(dt)
```

![](https://picgo-wutao.oss-cn-shanghai.aliyuncs.com/image-20240118111542-urhi2ix.png)​

### Reproducible analysis report

The code for training sample pre-processing and generation of auxiliary data can be found in `scripts`. The code for analysis and figures have been integrate into Rmarkdown and can be view online. Data for analysis can be found in `data`.    

### Acknowledgement

We thank ShanghaiTech University High Performance Computing Public Service Platform for computing services. We thank multi-omics facility, molecular and cell biology core facility of ShanghaiTech University for technical help. This work is supported by cross disciplinary Research Fund of Shanghai Ninth People's Hospital, Shanghai JiaoTong University School of Medicine (JYJC202227). Shanghai Science and Technology Commission (21ZR1442400), National Natural Science Foundation of China (82373149), and startup funding from ShanghaiTech University.

### Citation


