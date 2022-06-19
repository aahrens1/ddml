library("DoubleML")
library("haven")
library("mlr3")
library(data.table)
library("dplyr")

sessionInfo()

######################

dta<-read_dta("http://fmwww.bc.edu/repec/bocode/j/jtpa.dta") %>%
  zap_label() %>% zap_labels()
 

ml_g = lrn("regr.lm")
ml_m = lrn("classif.log_reg")
ml_r = lrn("classif.log_reg")
smpl <- list(list(train_ids = list(1:5602,5603:11204),
                  test_ids = list(5603:11204,1:5602)))
dta <- dta[,c("sex","age","married","black","hispanic","lnearnings","training","assignmt")]
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("sex","age","married","black","hispanic"),
                    y_col="earnings",
                    d_cols="training",
                    z_cols="assignmt"
                    )
 
### interactive
DML_obj <- NULL
DML_obj <- DoubleMLIIVM$new(dml_plr_obj,
            ml_g=ml_g,
            ml_m=ml_m,
            ml_r=ml_r
            )
DML_obj$set_sample_splitting(smpl)
fit <- DML_obj$fit()
# print(fit)
DML_obj$summary()


######################


rm(list=ls())
dta<-read_dta("http://www.stata-press.com/data/r13/cattaneo2.dta") %>%
  zap_label() %>% zap_labels()

learner <- lrn("regr.lm")
ml_g = lrn("regr.lm")
ml_m = lrn("classif.log_reg")
smpl <- list(list(train_ids = list(1:2321,2322:4642),
                  test_ids = list(2322:4642,1:2321)))
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols="mage",
                    y_col="bweight",
                    d_cols="mbsmoke"
                    )
dml_plr_obj_cl = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols="mage",
                    y_col="bweight",
                    d_cols="mbsmoke",
                    cluster_cols="medu"
                    )



### interactive
DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj,
            ml_g=ml_g,
            ml_m=ml_m
            )
DML_obj$set_sample_splitting(smpl)
fit <- DML_obj$fit()
# print(fit)
DML_obj$summary()

### interactive
DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj_cl,
            ml_g=ml_g,
            ml_m=ml_m
            )
DML_obj$set_sample_splitting(smpl)
fit <- DML_obj$fit()
# print(fit)
DML_obj$summary()


### partial 
#DML_obj <- DoubleMLIRM$new(dml_plr_obj,
DML_obj <- DoubleMLPLR$new(dml_plr_obj,
            ml_l=ml_g,
            ml_m=ml_m
            )
DML_obj$set_sample_splitting(smpl)
fit <- DML_obj$fit()
print(fit)
DML_obj$summary()


