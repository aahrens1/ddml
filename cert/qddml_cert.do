clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/pystacked"
	cd "/Users/kahrens/MyProjects/ddml/cert"
}

cap log close
log using "qddml_cert", replace text

global tol = 0.0001
which ddml
which pystacked
	
**** Partially linear model.

use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear
global Y net_tfa
global D e401
global X tw age inc fsize educ db marr twoearn pira hown

set seed 42
ddml init partial, kfolds(2)
ddml E[Y|X]: pystacked $Y $X, type(reg) method(rf)
ddml E[D|X]: pystacked $D $X, type(reg) method(rf)
ddml crossfit
ddml estimate  
local b1 = _b[$D]

set seed 42
qddml $Y $D ($X), kfolds(2) model(partial) cmd(pystacked) cmdopt(type(reg) method(rf)) noreg
local b2 = _b[$D]
ddml estimate 

assert reldif(`b1',`b2')<$tol

**** Partially linear IV model.

use https://statalasso.github.io/dta/AJR.dta, clear
global Y logpgp95
global D avexpr
global Z logem4
global X lat_abst edes1975 avelf temp* humid* steplow-oilres

set seed 42
ddml init iv, kfolds(30)
ddml E[Y|X], vtype(none): rforest $Y $X, type(reg)
ddml E[D|X], vtype(none): rforest $D $X, type(reg)
ddml E[Z|X], vtype(none): rforest $Z $X, type(reg)
ddml crossfit
ddml estimate, robust
local b1 = _b[$D]

set seed 42
qddml $Y ($X) ($D=$Z), kfolds(30) model(iv) cmd(rforest) cmdopt(type(reg)) vtype(none) robust noreg
local b2 = _b[$D]

assert reldif(`b1',`b2')<$tol

**** Interactive model--ATE and ATET estimation.

webuse cattaneo2, clear
global Y bweight
global D mbsmoke
global X mage prenatal1 mmarried fbaby mage medu
set seed 42

ddml init interactive, kfolds(5) reps(5)
ddml E[Y|X,D]: pystacked $Y $X, type(reg) method(gradboost)
ddml E[D|X]: pystacked $D $X, type(class) method(gradboost)
ddml crossfit
ddml estimate
local b1 = _b[$D]

ddml estimate, atet
local b1_atet=_b[$D]

set seed 42
qddml $Y $D ($X), kfolds(5) reps(5) model(interactive) cmd(pystacked) ycmdopt(type(reg) method(gradboost)) dcmdopt(type(class) method(gradboost)) noreg
local b2 = _b[$D]

ddml estimate, atet
local b2_atet=_b[$D]

assert reldif(`b1',`b2')<$tol
assert reldif(`b1_atet',`b2_atet')<$tol


**** Interactive IV model--LATE estimation.

use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta,clear
global Y earnings
global D training
global Z assignmt
global X sex age married black hispanic
set seed 42

ddml init interactiveiv, kfolds(5)
ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), type(reg) m(lassocv)
ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), type(class) m(lassocv)
ddml E[Z|X]: pystacked $Z c.($X)# #c($X), type(class) m(lassocv)
ddml crossfit
ddml estimate
local b1 = _b[$D]

set seed 42
qddml $Y (c.($X)# #c($X)) ($D=$Z), kfolds(5) model(interactiveiv) cmd(pystacked) ycmdopt(type(reg) m(lassocv)) dcmdopt(type(class) m(lassocv)) zcmdopt(type(class) m(lassocv)) noreg
local b2 = _b[$D]

assert reldif(`b1',`b2')<$tol

**** FIV

use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear
global Y share
global D price
global X hpwt air mpd space
global Z sum*

set seed 42
ddml init fiv
ddml E[Y|X]: pystacked $Y $X, type(reg)
ddml E[D|Z,X], learner(Dhat_pystacked): pystacked $D $X $Z, type(reg)
ddml E[D|X], learner(Dhat_pystacked) vname($D): pystacked {D} $X, type(reg)
ddml crossfit
ddml estimate
local b1 = _b[$D]

set seed 42
qddml $Y ($X) ($D=$Z), model(fiv) cmd(pystacked) cmdopt(type(reg)) noreg
local b2 = _b[$D]

assert reldif(`b1',`b2')<$tol

log close
