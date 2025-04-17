#########grid search
library(dplyr)
library(yardstick)
library(ggprism)
library(ggplot2)
library(ggpubr)

get_roc <- function(dt){
  dt <- dt %>% 
    mutate(truth = ifelse(label == 1, "Class1","Class2"),
           pred_label = ifelse(preds == 1, "Class1","Class2"))
  dt$truth <- factor(dt$truth)
  dt$pred_label <- factor(dt$pred_label)
  
  f1 <- try(f_meas(dt,truth,pred_label)[".estimate"] %>%
              unlist() %>% unname() %>% round(.,2),silent = TRUE)
  f1 <- ifelse('try-error' %in% class(f1),NA,f1)
  roc <-  roc_auc(dt, truth, preds_raw)[".estimate"] %>% 
    unlist() %>% unname() %>% round(.,2)
  pr <-  pr_auc(dt, truth, preds_raw)[".estimate"] %>%
    unlist() %>% unname() %>% round(.,2)
  return(c(roc,pr,f1,kappa))
}

all_cv <- c()
for (i in c(256, 512)){
  for (j in c(2,3,4)){
    for (z in c("0.01","0.001","0.0001")){
      all_cv <- append(all_cv,paste0("hd_",i,"_hn_",j,"_lr_",z))
    }
  }
}

cv_res <- vector("list",length(all_cv))
for (i in seq_along(cv_res)){
  sub_res <- vector("list",5)
  for (j in 0:4){
    pre <- read.csv(paste0("data/grid_search/cv/",all_cv[i],"_fold_",j,".csv")) %>% select(-X)
    pre_res <- get_roc(pre)
    sub_res[[j+1]] <- data.frame(
      fold = paste0("Fold-",j),
      ROC = pre_res[1],
      PR = pre_res[2],
      F1 = pre_res[3],
      Kappa = pre_res[4]
    ) 
  }
  sub_res <- bind_rows(sub_res)
  sub_res$para <- all_cv[i]
  cv_res[[i]] <- sub_res
}

cv_res <- bind_rows(cv_res)
tt <- cv_res %>% 
  filter(!is.na(F1)) %>% group_by(para) %>% 
  summarise(mean_roc = mean(ROC),
            mean_pr = mean(PR),
            mean_f1 = mean(F1),
            range_roc = max(ROC)-min(ROC),
            range_pr = max(PR)-min(PR),
            range_f1 = max(F1)-min(F1)
  ) %>% 
  ungroup() %>% 
  mutate(roc_pr_f1 = mean_roc + mean_pr + mean_f1,
         range_score = range_roc + range_pr + range_f1) %>% 
  arrange(-roc_pr_f1,range_score)

xlsx::write.xlsx(cv_res,file = "data/grid_search/grid_search_res.xlsx",sheetName = "all_res",append=TRUE)
xlsx::write.xlsx(tt,file = "data/grid_search/grid_search_res.xlsx",sheetName = "summ_res",append=TRUE)