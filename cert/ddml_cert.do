
clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/pystacked"
	adopath + "/Users/kahrens/MyProjects/ddml"
	cd "/Users/kahrens/MyProjects/ddml/cert"
}

cap log close
log using "ddml_cert", replace text

// net install ddml, from(https://raw.githubusercontent.com/aahrens1/ddml/dev/) replace

// necessary programs for cert; script exits with error if not installed
findfile pystacked.ado
findfile svmachines.ado
findfile rforest.ado

global tol = 0.0001
which ddml
which pystacked


******************************************************************************** 
**** Partially-linear model													****
******************************************************************************** 

use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear
global Y net_tfa
global D e401
global X tw age inc fsize educ db marr twoearn pira hown
set seed 42

sample 30 

ddml init partial, kfolds(2) 

ddml E[Y|X]: reg $Y $X
ddml E[Y|X]: pystacked $Y $X, type(reg) method(rf gradboost)
ddml E[Y|X]: svmachines $Y $X, type(svr)
ddml E[D|X]: reg $D $X
ddml E[D|X]: pystacked $D $X, type(reg) method(rf gradboost)
ddml E[D|X]: svmachines $D $X, type(svr)
// add eqn if parsnip installed
cap findfile parsnip2.ado
local hasparsnip=_rc
cap findfile parsnip.ado
local hasparsnip=`hasparsnip'+_rc
if `hasparsnip'==0 {
	ddml E[D|X]: parsnip2 $D $X, model(linear_reg) engine(glmnet) penalty(.5) clearR
}
ddml desc

ddml crossfit 

ddml estimate, robust
local b1=_b[e401]
cap drop Y_
cap drop D_
gen double Y_ = $Y - Y2_pystacked_1
gen double D_ = $D - D2_pystacked_1
reg Y_ D_, robust
local b2=_b[D_]
assert reldif(`b1',`b2')<$tol

ddml estimate, robust allcombos

ddml extract, show(mse)
ddml extract, show(n)

ddml drop 

// check that everything was deleted
cap ddml desc
assert _rc==3259

// pystacked as a single learner with multiple base learners;
// shortstacking; multiple reps; ddml extract
ddml init partial, kfolds(2) reps(2)

ddml E[Y|X]: pystacked $Y $X, type(reg) method(rf gradboost)
ddml E[D|X]: pystacked $D $X, type(reg) method(rf gradboost)
ddml desc

ddml crossfit, shortstack

ddml estimate, robust

ddml extract, show(pystacked)
ddml extract, show(weights)

******************************************************************************** 
**** Partially-linear IV model										        ****
******************************************************************************** 


use https://statalasso.github.io/dta/AJR.dta, clear
global Y logpgp95
global D avexpr
global Z logem4
global X lat_abst edes1975 avelf temp* humid* steplow-oilres
set seed 42

ddml init iv, kfolds(30)

ddml E[Y|X]: reg $Y $X
ddml E[Y|X], vtype(none): rforest $Y $X, type(reg)
ddml E[D|X]: reg $D $X
ddml E[D|X], vtype(none): rforest $D $X, type(reg)
ddml E[Z|X]: reg $Z $X
ddml E[Z|X], vtype(none): rforest $Z $X, type(reg)

ddml crossfit
ddml estimate, robust
local a = _b[ave]

cap drop Y_
cap drop D_
cap drop Z_
gen double Y_ = $Y - Y2_rforest_1
gen double D_ = $D - D2_rforest_1
gen double Z_ = $Z - Z2_rforest_1
ivreg Y_ (D_ = Z_) , robust
local b = _b[D_]
assert reldif(`a',`b')<$tol


******************************************************************************** 
**** Partially-linear IV model with multiple treatments				        ****
******************************************************************************** 


use https://statalasso.github.io/dta/AJR.dta, clear
global Y logpgp95
global D1 avexpr
global D2 democ1
global Z1 logem4
global Z2 lat_abst 
global X edes1975 avelf temp* humid* steplow-oilres
set seed 42

ddml init iv, kfolds(30)

ddml E[Y|X]: reg $Y $X
ddml E[D|X]: reg $D1 $X
ddml E[D|X]: reg $D2 $X
ddml E[Z|X]: reg $Z1 $X
ddml E[Z|X]: reg $Z2 $X

ddml crossfit

ddml estimate, robust
local a1 = _b[ave]
local a2 = _b[demo]

ddml extract, show(mse)
ddml extract, show(n)

cap drop Y_
cap drop D_1 D_2
cap drop Z_1 Z_2
gen double Y_ = $Y - Y1_reg_1
gen double D_1 = $D1 - D1_reg_1
gen double D_2 = $D2 - D2_reg_1
gen double Z_1 = $Z1 - Z1_reg_1
gen double Z_2 = $Z2 - Z2_reg_1
ivreg Y_ (D_1 D_2 = Z_1 Z_2) , robust
local b1 = _b[D_1]
local b2 = _b[D_2]
assert reldif(`a1',`b1')<$tol
assert reldif(`a2',`b2')<$tol


******************************************************************************** 
**** Partially-linear model with repetitions								****
******************************************************************************** 

use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear
global Y net_tfa
global D e401
global X tw age inc fsize educ db marr twoearn pira hown
set seed 42

sample 30 

ddml init partial, kfolds(2) reps(3)

ddml E[Y|X]: reg $Y $X
ddml E[Y|X]: pystacked $Y $X, type(reg) method(rf)
ddml E[D|X]: reg $D $X
ddml E[D|X]: pystacked $D $X, type(reg) method(rf)

ddml crossfit 
ddml estimate, robust

ddml estimate, robust allcombos  
ddml estimate, spec(3) rep(1) notable replay 
cap drop bhat
gen bhat=.
replace bhat = _b[e401] if _n==1
ddml estimate, spec(3) rep(2) notable replay 
replace bhat = _b[e401] if _n==2
ddml estimate, spec(3) rep(3) notable replay 
replace bhat = _b[e401] if _n==3

ddml estimate, spec(mse) rep(md) notable replay 
local bmedian =_b[e401]
sum bhat, detail
assert reldif(`r(p50)',`bmedian')<$tol

ddml estimate, spec(mse) rep(mn) notable replay 
local bmean =_b[e401]
sum bhat, detail
assert reldif(`r(mean)',`bmean')<$tol

******************************************************************************** 
**** Partially-linear model with 2 treatments								****
******************************************************************************** 

use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear
global Y net_tfa
global D1 e401
global D2 educ
global X tw age inc fsize db marr twoearn pira hown
set seed 42

sample 30 

ddml init partial, kfolds(2)

ddml E[Y|X]: reg $Y $X
ddml E[Y|X]: pystacked $Y $X, type(reg) method(rf lassocv)
ddml E[D|X]: reg $D1 $X
ddml E[D|X]: pystacked $D1 $X, type(reg) method(rf lassocv)
ddml E[D|X]: reg $D2 $X
ddml E[D|X]: pystacked $D2 $X, type(reg) method(rf lassocv )

ddml crossfit 

ddml estimate, robust 
local a1 = _b[e401]
local a2 = _b[educ]
local a3 = _b[_cons]
cap drop Y_
cap drop D_1
cap drop D_2
gen Y_ = $Y - Y2_pystacked_1
gen D_1 = $D1 - D2_pystacked_1
gen D_2 = $D2 - D4_pystacked_1
reg Y_ D_1 D_2,robust
local b1 = _b[D_1]
local b2 = _b[D_2]
local b3 = _b[_cons]
assert reldif(`a1',`b1')<$tol
assert reldif(`a2',`b2')<$tol
assert reldif(`a3',`b3')<$tol

ddml extract, show(mse)
ddml extract, show(n)


******************************************************************************** 
**** Interactive model--ATE and ATET estimation.							****
******************************************************************************** 

webuse cattaneo2, clear
global Y bweight
global D mbsmoke
global X mage prenatal1 mmarried fbaby mage medu
set seed 42

sample 30

ddml init interactive, kfolds(5) reps(5)

ddml E[Y|X,D]: reg $Y $X
ddml E[Y|X,D]: pystacked $Y $X, type(reg) method(gradboost rf)
ddml E[D|X]: logit $D $X
ddml E[D|X]: pystacked $D $X, type(class) method(gradboost rf)
ddml E[D|X]: svmachines $D $X, type(svc) 

ddml crossfit, shortstack

ddml estimate
ddml estimate, atet trim(0)

ddml extract, show(mse)
ddml extract, show(n)

******************************************************************************** 
**** Interactive IV model--LATE estimation.     							****
******************************************************************************** 


use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta,clear
keep in 1/1000
global Y earnings
global D training
global Z assignmt
global X sex age married black hispanic
set seed 42

ddml init interactiveiv, kfolds(5) reps(2)

ddml E[Y|X,Z]: reg $Y $X
ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), type(reg) m(lassocv rf )
ddml E[D|X,Z]: logit $D $X
ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), type(class) m(lassocv rf)
ddml E[Z|X]: logit $Z $X
ddml E[Z|X]: pystacked $Z c.($X)# #c($X), type(class) m(lassocv rf)

ddml crossfit, shortstack
ddml estimate

ddml extract, show(ssweights)
ddml extract, show(mse)
ddml extract, show(n)

******************************************************************************** 
**** Flexible IV							    							****
******************************************************************************** 

use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear
global Y share
global D price
global X hpwt air mpd space
global Z sum*
set seed 42

ddml init fiv

ddml E[Y|X]: reg $Y $X
ddml E[Y|X]: pystacked $Y $X, type(reg)

ddml E[D|Z,X], learner(Dhat_reg): reg $D $X $Z
ddml E[D|Z,X], learner(Dhat_pystacked): pystacked $D $X $Z, type(reg)

ddml E[D|X], learner(Dhat_reg) vname($D): reg {D} $X
ddml E[D|X], learner(Dhat_pystacked) vname($D): pystacked {D} $X, type(reg)
 
ddml crossfit
ddml estimate
local a1 = _b[price]

ddml extract, show(mse)
ddml extract, show(n)

cap drop Ytilde
cap drop Dtilde
cap drop Zopt
gen Ytilde = $Y - Y2_pystacked_1
gen Dtilde = $D - Dhat_pystacked_h_1
gen Zopt = Dhat_pystacked_1 - Dhat_pystacked_h_1
 
ivreg Ytilde (Dtilde=Zopt) 
local b1 = _b[Dtil]
assert reldif(`a1',`b1')<$tol

log close
