library("DoubleML")
library("haven")
library("mlr3")
library(data.table)
library("dplyr")

sessionInfo()

rm(list=ls())

all_results <- as.data.frame(NULL)

setwd("/Users/kahrens/MyProjects/ddml/cert")

##################################################################
###################### cattaneo2 data. 
##################################################################


dta<-read_dta("cattaneo2.dta") %>%
  zap_label() %>% zap_labels()

smpl <- list(list(train_ids = list(1:2321,2322:4642),
                  test_ids = list(2322:4642,1:2321)))

dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("prenatal1", "mmarried", "fbaby", "mage", "medu"),
                    y_col="bweight",
                    d_cols="mbsmoke"
                    )

### interactive: ATE
DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="cattaneo2",
                                model="interactive ATE",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))

### interactive: ATTE
DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg"),
            score="ATTE"
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="cattaneo2",
                                model="interactive ATTE",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))

### partial linear model
DML_obj <- NULL
DML_obj <- DoubleMLPLR$new(dml_plr_obj,
            ml_l=lrn("regr.lm"),
            ml_m=lrn("regr.lm")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="cattaneo2",
                                model="partial",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))

##################################################################
###################### jtpa. 
##################################################################

dta<-read_dta("jtpa.dta") %>%
  zap_label() %>% zap_labels()
 
smpl <- list(list(train_ids = list(1:5602,5603:11204),
                  test_ids = list(5603:11204,1:5602)))
dta <- dta[,c("sex","age","married","black","hispanic","earnings","training","assignmt")]

### interactive IV
dml_plr_obj <- NULL
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("sex","age","married","black","hispanic"),
                    y_col="earnings",
                    d_cols="training",
                    z_cols="assignmt"
                    )

DML_obj <- NULL
DML_obj <- DoubleMLIIVM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg"),
            ml_r=lrn("classif.log_reg")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="jtpa",
                                model="interactive IV",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))

### partial IV
DML_obj <- NULL
DML_obj <- DoubleMLPLIV$new(dml_plr_obj,
            ml_l=lrn("regr.lm"),
            ml_m=lrn("regr.lm"),
            ml_r=lrn("regr.lm")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="jtpa",
                                model="partial IV",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))


### interactive
dml_plr_obj <- NULL
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("sex","age","married","black","hispanic"),
                    y_col="earnings",
                    d_cols="assignmt"
                    )

DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="jtpa",
                                model="interactive ATE",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))

DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg"),
            score="ATTE"
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="jtpa",
                                model="interactive ATTE",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))



### partial
dml_plr_obj <- NULL
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("sex","age","married","black","hispanic"),
                    y_col="earnings",
                    d_cols="assignmt"
                    )

DML_obj <- NULL
DML_obj <- DoubleMLPLR$new(dml_plr_obj,
            ml_l=lrn("regr.lm"),
            ml_m=lrn("regr.lm")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="jtpa",
                                model="partial",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))



##################################################################
###################### 401k 
##################################################################

dta <- read_dta("sipp1991.dta")

id <- (1:9915 %% 3)+1
f1 <- which(id==1)
f2 <- which(id==2)
f3 <- which(id==3)
smpl <- list(list(train_ids = list(c(f1,f2),c(f2,f3),c(f1,f3)),
                  test_ids = list(f3,f1,f2)))

### partial
dml_plr_obj <- NULL
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("tw","age","inc","fsize","educ","db","marr","twoearn","pira","hown"),
                    y_col="net_tfa",
                    d_cols="e401"
                    )

DML_obj <- NULL
DML_obj <- DoubleMLPLR$new(dml_plr_obj,
            ml_l=lrn("regr.lm"),
            ml_m=lrn("regr.lm")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="401k",
                                model="partial",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))


### interactive ATE
dml_plr_obj <- NULL
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("tw","age","inc","fsize","educ","db","marr","twoearn","pira","hown"),
                    y_col="net_tfa",
                    d_cols="e401"
                    )

DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg"),
            score="ATE"
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="401k",
                                model="interactive ATE",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))

DML_obj <- NULL
DML_obj <- DoubleMLIRM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg"),
            score="ATTE"
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="401k",
                                model="interactive ATTE",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))


### interactive IV
dml_plr_obj <- NULL
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("tw","age","inc","fsize","educ","db","marr","twoearn","pira","hown"),
                    y_col="net_tfa",
                    d_cols="p401",
                    z_cols="e401"
                    )

DML_obj <- NULL
DML_obj <- DoubleMLIIVM$new(dml_plr_obj,
            ml_g=lrn("regr.lm"),
            ml_m=lrn("classif.log_reg"),
            ml_r=lrn("classif.log_reg")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="401k",
                                model="interactive IV",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))

### partial IV
dml_plr_obj <- NULL
dml_plr_obj = double_ml_data_from_data_frame(
                    df=as.data.table(dta),
                    x_cols=c("tw","age","inc","fsize","educ","db","marr","twoearn","pira","hown"),
                    y_col="net_tfa",
                    d_cols="p401",
                    z_cols="e401"
                    )

DML_obj <- NULL
DML_obj <- DoubleMLPLIV$new(dml_plr_obj,
            ml_l=lrn("regr.lm"),
            ml_m=lrn("regr.lm"),
            ml_r=lrn("regr.lm")
            )
DML_obj$set_sample_splitting(smpl)
DML_obj$fit()
DML_obj$summary()
DML_obj$coef
DML_obj$se
all_results <- rbind(all_results,
                     data.frame(application="401k",
                                model="partial IV",
                                coef = as.numeric(DML_obj$coef),
                                se = as.numeric(DML_obj$se)
                                ))
 
all_results
library(readr)
write_csv(all_results,"DoubleML_results.csv")
 