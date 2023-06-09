clear all

if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
}

cap cd "/Users/kahrens/MyProjects/ddml/cert"
cap cd "C:\LocalStore\ecomes\Documents\GitHub\ddml\cert"

cap log close
log using "ddml_cert_interactive", replace text

which ddml
mata: whichddml()

webuse cattaneo2, clear
keep in 1/1000

gen double mage_sq = mage^2

which ddml
which pystacked

// necessary programs for cert; script exits with error if not installed
findfile pystacked.ado

set seed 123

********************************************************************************
*** interactive model	 													 ***
********************************************************************************

global Y bweight
global D mbsmoke
global X c.(mmarried mage mage_sq fbaby medu)#c.(mmarried mage mage_sq fbaby medu)

*** pystacked, no SS

*** initialise ddml and select model; 
ddml init interactive, kfolds(2) reps(2)
ddml E[Y|X,D]: pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[D|X]: pystacked $D $X, type(class) method(logit gradboost)
ddml crossfit
ddml estimate
*** replay
ddml estimate, mname(m0) spec(st) rep(1) replay notable
*** append, estimate, replay
ddml sample, append(1)
ddml crossfit
ddml estimate
*** replay
ddml estimate, mname(m0) spec(st) rep(1) replay notable
ddml estimate, mname(m0) spec(st) rep(2) replay notable
ddml estimate, mname(m0) spec(st) rep(mn) replay notable
ddml estimate, mname(m0) spec(st) rep(md) replay notable

*** pystacked, SS and PS

*** initialise ddml and select model; 
ddml init interactive, kfolds(2) reps(2)
ddml E[Y|X,D]: pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[D|X]: pystacked $D $X, type(class) method(logit gradboost)
ddml crossfit, shortstack poolstack
ddml estimate
*** replay
ddml estimate, mname(m0) spec(st) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
ddml estimate, mname(m0) spec(ps) rep(1) replay notable

*** append, estimate, replay
ddml sample, append(1)
ddml crossfit, shortstack poolstack
ddml estimate
*** replay
ddml estimate, mname(m0) spec(st) rep(1) replay notable
ddml estimate, mname(m0) spec(st) rep(2) replay notable
ddml estimate, mname(m0) spec(st) rep(mn) replay notable
ddml estimate, mname(m0) spec(st) rep(md) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(2) replay notable
ddml estimate, mname(m0) spec(ss) rep(mn) replay notable
ddml estimate, mname(m0) spec(ss) rep(md) replay notable
ddml estimate, mname(m0) spec(ps) rep(1) replay notable
ddml estimate, mname(m0) spec(ps) rep(2) replay notable
ddml estimate, mname(m0) spec(ps) rep(mn) replay notable
ddml estimate, mname(m0) spec(ps) rep(md) replay notable

*** multiple learners, no SS

*** initialise ddml and select model; 
ddml init interactive, kfolds(2) reps(2)
ddml E[Y|X,D]: pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[Y|X,D]: reg $Y $X
ddml E[D|X]: pystacked $D $X, type(class) method(logit gradboost)
ddml E[D|X]: logit $D $X
ddml crossfit
ddml estimate
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
*** allcombos
ddml estimate, allcombos
forvalues i=1/8 {
    ddml estimate, mname(m0) spec(`i') rep(1) replay notable
}
*** append, estimate, replay
ddml sample, append(1)
ddml crossfit
ddml estimate
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
ddml estimate, mname(m0) spec(mse) rep(2) replay notable
ddml estimate, mname(m0) spec(mse) rep(mn) replay notable
ddml estimate, mname(m0) spec(mse) rep(md) replay notable
*** allcombos
ddml estimate, allcombos
forvalues i=1/4 {
	forvalues r=1/2 {
		ddml estimate, mname(m0) spec(`i') rep(`r') replay notable
	}
}
ddml estimate, mname(m0) spec(mse) rep(mn) replay notable
ddml estimate, mname(m0) spec(mse) rep(md) replay notable

*** multiple learners, SS

*** initialise ddml and select model; 
ddml init interactive, kfolds(2) reps(2)
ddml E[Y|X,D]: pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[Y|X,D]: reg $Y $X
ddml E[D|X]: pystacked $D $X, type(class) method(logit gradboost)
ddml E[D|X]: logit $D $X
ddml crossfit, shortstack
ddml estimate
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
*** allcombos
ddml estimate, allcombos
forvalues i=1/8 {
    ddml estimate, mname(m0) spec(`i') rep(1) replay notable
}
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
*** append, estimate, replay
ddml sample, append(1)
ddml crossfit, shortstack
ddml estimate
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
ddml estimate, mname(m0) spec(mse) rep(2) replay notable
ddml estimate, mname(m0) spec(mse) rep(mn) replay notable
ddml estimate, mname(m0) spec(mse) rep(md) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(2) replay notable
ddml estimate, mname(m0) spec(ss) rep(mn) replay notable
ddml estimate, mname(m0) spec(ss) rep(md) replay notable
*** allcombos
ddml estimate, allcombos
forvalues i=1/8 {
    forvalues r=1/2 {
		ddml estimate, mname(m0) spec(`i') rep(`r') replay notable
	}
}
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
ddml estimate, mname(m0) spec(mse) rep(2) replay notable
ddml estimate, mname(m0) spec(mse) rep(mn) replay notable
ddml estimate, mname(m0) spec(mse) rep(md) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(2) replay notable
ddml estimate, mname(m0) spec(ss) rep(mn) replay notable
ddml estimate, mname(m0) spec(ss) rep(md) replay notable

*** ddml overlap
ddml init interactive, kfolds(2) reps(2)
ddml E[Y|X,D]: pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[D|X]: pystacked $D $X, type(class) method(logit gradboost)
ddml crossfit
ddml estimate
ddml overlap
ddml overlap, replist(1)
ddml overlap, pslist(D1_pystacked_L1 D1_pystacked_L2)
ddml overlap, name(triangle, replace)							///
	title("Propensity score: triangle kernel")
ddml overlap, kernel(epanechnikov) name(epanechnikov, replace)	///
	title("Propensity score: epanechnikov kernel")

log close


