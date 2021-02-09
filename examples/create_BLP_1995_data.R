#' ddml Demo using the Berry, Levinson, Pakes (1995) data  =====================
#'
#' This R script is a demonstration of the ddml R package. We revisit the
#'     empirical example of Chernozhukov, Hansen, and Spindler (2015)
#'     (CHS2015, hereafter), which extends the instruments of Berry, Levinson,
#'     and Pakes (1995) (BLP1995, hereafter) and applies an instrument selection
#'     procedure based on the Lasso. We consider the same instrument extension
#'     and apply a selection of ensemble procedures that combines conventional
#'     linear estimators with computational alternatives including Lasso-based
#'     approaches and random forests.
#'
#' The BLP1995 data is included as an exemplary dataset with the ddml package
#'     under the name `BLP_1995'.

# Reproduction of LS and TSLS Estimates from CHS2015 ===========================

library(ddml)

# Data prepraration
n <- length(BLP_1995$model.id)
y <- log(BLP_1995$share) -  log(BLP_1995$outshr)
x1 <- as.matrix(cbind(1, BLP_1995[, c("hpwt", "air", "mpd", "space", "price")]))

#' Construct the BLP1995 instruments. Instruments are sums of product
#'     characteristics (excluding price and other potentially endogenous
#'     variables) of other products offered by the firm as well as over
#'     competing firms. The exact instrument specification here follows the
#'     approach of CHS2015.
X_ <- x1[, 1:5]; ncol_X <- ncol(X_) # exclude price
sum_other <- sum_rival <- matrix(0, n, 5)
for (i in 1:n) {
  other_ind <- (BLP_1995$firmid == BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i]) & (BLP_1995$id != BLP_1995$id[i])
  rival_ind <- (BLP_1995$firmid != BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i])
  sum_other[i, ] <- colSums(X_[other_ind, , drop = FALSE])
  sum_rival[i, ] <- colSums(X_[rival_ind, , drop = FALSE])
}#FOR
Z_ <- cbind(sum_other, sum_rival); ncol_Z <- ncol(Z_)

#' Calculate LS and TSLS estimates corresponding to price. Note that the TSLS
#'     estimates differ from those in CHS2015. This is due to a slight
#'     instrument-construction error in the code of the CHS2015.
summary(ols(y, x1), type = "HC1")$res[6, ]
summary(tsls(y, D = BLP_1995$price, Z = Z_, X = X_), type = "HC1")$res[1, ]

# Extending the BLP1995 Instruments ============================================

#' Here, we follow the extension of the instruments as in CHS2015. In
#'     particular, various interaction terms are considered.
tu = BLP_1995$trend/19;
mpdu = BLP_1995$mpd/7;
spaceu = BLP_1995$space/2;
XL_ <- as.matrix(cbind(1, BLP_1995[, c("hpwt", "air")], mpdu, spaceu, tu,
                       BLP_1995$hpwt^2, BLP_1995$hpwt^3, mpdu^2, mpdu^3,
                       spaceu^2, spaceu^3, tu^2, tu^3, BLP_1995$hpwt *
                         BLP_1995$air,  mpdu * BLP_1995$air, spaceu *
                         BLP_1995$air, tu * BLP_1995$air, BLP_1995$hpwt *
                         mpdu, BLP_1995$hpwt * spaceu, BLP_1995$hpwt * tu,
                       mpdu * spaceu,  mpdu * tu, spaceu * tu))
ncol_XL <- ncol(XL_)
sum_otherL <- sum_rivalL <- matrix(0, n, 24)
for (i in 1:n) {
  other_ind <- (BLP_1995$firmid == BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i]) & (BLP_1995$id != BLP_1995$id[i])
  rival_ind <- (BLP_1995$firmid != BLP_1995$firmid[i]) &
    (BLP_1995$cdid == BLP_1995$cdid[i])
  sum_otherL[i, ] <- colSums(XL_[other_ind, , drop = FALSE])
  sum_rivalL[i, ] <- colSums(XL_[rival_ind, , drop = FALSE])
}#FOR
ZL_ <- cbind(sum_otherL,sum_rivalL); ncol_ZL <- ncol(ZL_)

# TSLS estimates on the extended set of instruments
summary(tsls(y, D = BLP_1995$price, Z = ZL_, X = XL_), type = "HC1")$res[1, ]


##########################################################
### export data as csv                                 ###
##########################################################

XL_ <- janitor::clean_names(as.data.frame(XL_))
XL_ <- XL_[,-1] # remove constant

colnames(ZL_) <- paste("z",1:ncol(ZL_),sep="")

data <- cbind(y,BLP_1995$price,ZL_,XL_)
colnames(data) <- c("share","price",colnames(ZL_),colnames(XL_))

readr::write_csv(data,"/Users/kahrens/MyProjects/ddml/cert/BLP1995.csv")


