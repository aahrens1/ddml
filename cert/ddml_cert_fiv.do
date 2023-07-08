clear all

if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
}

cap cd "/Users/kahrens/MyProjects/ddml/cert"
cap cd "C:\LocalStore\ecomes\Documents\GitHub\ddml\cert"

cap log close
log using "ddml_cert_fiv", replace text

which ddml
mata: whichddml()

use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear

// necessary programs for cert; script exits with error if not installed
findfile pystacked.ado

set seed 123

******************************************************************************** 
**** Flexible IV							    							****
******************************************************************************** 

global Y share
global D price
global X hpwt air mpd space
global Z sum*

// single learner = pystacked
ddml init fiv, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X, type(reg)
ddml E[D|Z,X], learner(Dhat_pystacked): pystacked $D $X $Z, type(reg)
ddml E[D|X], learner(Dhat_pystacked) vname($D): pystacked {D} $X, type(reg)
ddml crossfit
ddml estimate
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

// multiple learners = pystacked+gradboost, pystacked+lassocv, pystacked+rf
ddml init fiv, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X, type(reg) m(gradboost)
ddml E[Y|X]: pystacked $Y $X, type(reg) m(lassocv)
ddml E[Y|X]: pystacked $Y $X, type(reg) m(rf)
ddml E[D|Z,X], learner(Dhat_psgradboost): pystacked $D $X $Z, type(reg) m(gradboost)
ddml E[D|X], learner(Dhat_psgradboost) vname($D): pystacked {D} $X, type(reg) m(gradboost)
ddml E[D|Z,X], learner(Dhat_pslassocv): pystacked $D $X $Z, type(reg) m(lassocv)
ddml E[D|X], learner(Dhat_pslassocv) vname($D): pystacked {D} $X, type(reg) m(lassocv)
ddml E[D|Z,X], learner(Dhat_rf): pystacked $D $X $Z, type(reg) m(rf)
ddml E[D|X], learner(Dhat_rf) vname($D): pystacked {D} $X, type(reg) m(rf)
ddml crossfit
ddml estimate
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
*** append, estimate, replay
ddml sample, append(1)
ddml crossfit
ddml estimate
*** replay
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
ddml estimate, mname(m0) spec(mse) rep(2) replay notable
ddml estimate, mname(m0) spec(mse) rep(mn) replay notable
ddml estimate, mname(m0) spec(mse) rep(md) replay notable


// as above, with shortstacking
// multiple learners = pystacked+gradboost, pystacked+lassocv, pystacked+rf
ddml init fiv, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X, type(reg) m(gradboost)
ddml E[Y|X]: pystacked $Y $X, type(reg) m(lassocv)
ddml E[Y|X]: pystacked $Y $X, type(reg) m(rf)
ddml E[D|Z,X], learner(Dhat_psgradboost): pystacked $D $X $Z, type(reg) m(gradboost)
ddml E[D|X], learner(Dhat_psgradboost) vname($D): pystacked {D} $X, type(reg) m(gradboost)
ddml E[D|Z,X], learner(Dhat_pslassocv): pystacked $D $X $Z, type(reg) m(lassocv)
ddml E[D|X], learner(Dhat_pslassocv) vname($D): pystacked {D} $X, type(reg) m(lassocv)
ddml E[D|Z,X], learner(Dhat_rf): pystacked $D $X $Z, type(reg) m(rf)
ddml E[D|X], learner(Dhat_rf) vname($D): pystacked {D} $X, type(reg) m(rf)
ddml crossfit, shortstack
ddml estimate
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
*** append, estimate, replay
ddml sample, append(1)
ddml crossfit, shortstack
ddml estimate
*** replay
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
forvalues i=1/27 {
    ddml estimate, mname(m0) spec(`i') rep(1) replay notable
}
ddml estimate, mname(m0) spec(mse) rep(1) replay notable
ddml estimate, mname(m0) spec(mse) rep(2) replay notable
ddml estimate, mname(m0) spec(mse) rep(mn) replay notable
ddml estimate, mname(m0) spec(mse) rep(md) replay notable
ddml estimate, mname(m0) spec(ss) rep(1) replay notable
ddml estimate, mname(m0) spec(ss) rep(2) replay notable
ddml estimate, mname(m0) spec(ss) rep(mn) replay notable
ddml estimate, mname(m0) spec(ss) rep(md) replay notable

// poolstacking not supported with multiple calls to pystacked

log close
