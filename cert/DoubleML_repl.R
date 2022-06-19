library("DoubleML")
library("haven")
library("mlr3")
library("mlr3learner")
library(data.table)
library("dplyr")

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


