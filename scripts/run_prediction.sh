python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/model_data/cell_gene_exp_vs_normal_filter.csv -g /home/data/sda/wt/model_data/new_model/tmp/test2/ -c /home/data/sda/wt/model_data/new_model/enzyme_net/test/raw/test_cell_info.csv -n /home/data/sda/wt/model_data/enzyme_net/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/test_preV2.csv -d test -b 1

python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/model_data/cell_gene_exp_vs_normal_filter.csv -g /home/data/sda/wt/model_data/new_model/tmp/sanger2/ -c /home/data/sda/wt/model_data/new_model/enzyme_net_sanger/raw/cell_info.csv -n /home/data/sda/wt/model_data/enzyme_net_sanger/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/sanger_preV2.csv -d test -b 1

python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/model_data/cell_gene_exp_vs_normal_filter.csv -g /home/data/sda/wt/model_data/new_model/tmp/rnai2/ -c /home/data/sda/wt/model_data/new_model/enzyme_net_rnai/raw/cell_info.csv -n /home/data/sda/wt/model_data/enzyme_net_rnai/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/rnai_preV2.csv -d test -b 1

python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/model_data/ccma_gene_exp_vs_normal_filter.csv -g /home/data/sda/wt/model_data/new_model/tmp/ccma3/ -c ~/meta_target/data/ccma_pre_info.csv -n /home/data/sda/wt/model_data/enzyme_net_ccma/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/ccma_preV2.csv -d val -b 1

python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/model_data/cell_gene_exp_vs_normal_filter.csv -g /home/data/sda/wt/model_data/new_model/tmp/drug2/ -c ~/meta_target/data/nc_drug_15_info.csv -n /home/data/sda/wt/model_data/enzyme_net_drug/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/drug_preV2.csv -d val -b 1

python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/TCGA/tcga_gene_exp_vs_normal.csv -g /home/data/sda/wt/model_data/new_model/tmp/tcga2/ -c ~/meta_target/data/tcga_val_sample_info.csv -n /home/data/sda/wt/model_data/tcga_net2/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/tcga_preV2.csv -d val -b 1

python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/model_data/tcga_gene_exp_vs_normal_normal.csv -g /home/data/sda/wt/model_data/new_model/tmp/tcga_normal2/ -c ~/meta_target/data/tcga_normal_sample_info.csv -n /home/data/sda/wt/model_data/tcga_net_normal2/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/tcga_normal_preV2.csv -d val -b 1

python ~/meta_target/scripts/model/pred_enzyme.py -e /home/data/sda/wt/model_data/cell_gene_exp_vs_normal_normal.csv -g /home/data/sda/wt/model_data/new_model/tmp/cell_normal2/ -c ~/meta_target/data/normal_cell_info.csv -n /home/data/sda/wt/model_data/enzyme_net_normal/ -t 30 -m /home/data/sda/wt/model_data/new_model/enzyme_model_filterV2.pt -o /home/wt/meta_target/data/cell_normal_preV2.csv -d test -b 1

python ~/meta_target/scripts/model/PredictNewSamples.py -e /home/data/sda/wt/DeepDEP/drug2_exp_prediction.txt -f /home/data/sda/wt/DeepDEP/drug2_fingerprint_prediction.txt -o /home/data/sda/wt/DeepDEP/drug2_pred.txt

python ~/meta_target/scripts/model/PredictNewSamples.py -e /home/data/sda/wt/DeepDEP/CCMA2_exp_prediction.txt -f /home/data/sda/wt/DeepDEP/CCMA2_fingerprint_prediction.txt -o /home/data/sda/wt/DeepDEP/CCMA2_pred.txt



