clear all

if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
}

cap cd "/Users/kahrens/MyProjects/ddml/cert"
cap cd "C:\LocalStore\ecomes\Documents\GitHub\ddml\cert"

cap log close
log using "ddml_cert_partial", replace text

which ddml
mata: whichddml()

use https://statalasso.github.io/dta/AJR.dta, clear

// necessary programs for cert; script exits with error if not installed
findfile pystacked.ado

set seed 123

******************************************************************************** 
**** Partially linear model.												 ***
******************************************************************************** 

global Y logpgp95
global X lat_abst edes1975 temp* humid* steplow-oilres
global D1 avexpr 
global D2 democ1 
global D3 avelf

*** pystacked, no SS

*** initialise ddml and select model; 
ddml init partial, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X , type(reg)
ddml E[D|X]: pystacked $D1 $X , type(reg)
ddml E[D|X]: pystacked $D2 $X , type(reg)
ddml crossfit
ddml estimate, robust
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

*** pystacked, SS

*** initialise ddml and select model; 
ddml init partial, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X , type(reg)
ddml E[D|X]: pystacked $D1 $X , type(reg)
ddml E[D|X]: pystacked $D2 $X , type(reg)
ddml crossfit, shortstack poolstack
*** estimation of parameter of interest
ddml estimate, robust
*** replay
ddml estimate, mname(m0) spec(st) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
*** append, estimate, replay
ddml sample, append(1)
ddml crossfit, shortstack
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

*** multiple learners, no SS
ddml init partial, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X , type(reg)
ddml E[Y|X]: reg $Y $X 
ddml E[D|X]: pystacked $D1 $X , type(reg)
ddml E[D|X]: reg $D1 $X 
ddml crossfit
ddml estimate, robust
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
*** allcombos
ddml estimate, allcombos
forvalues i=1/4 {
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
ddml init partial, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X , type(reg)
ddml E[Y|X]: reg $Y $X 
ddml E[D|X]: pystacked $D1 $X , type(reg)
ddml E[D|X]: reg $D1 $X 
ddml crossfit, shortstack
ddml estimate, robust
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
*** allcombos
ddml estimate, allcombos shortstack
forvalues i=1/4 {
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
ddml estimate, allcombos shortstack
forvalues i=1/4 {
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

log close
