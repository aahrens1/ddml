clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pylearn2"
	cd "/Users/kahrens/MyProjects/ddml/examples"
}

** Note that you need to run the file "lddml.do" to install
** required Mata structures.
** This only needs to be done once.
//run "/Users/kahrens/MyProjects/ddml/lddml.do"

cap log close
log using "log_late.txt", replace text  
 
use "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta",clear   

which ddml
which pystacked
  
********************************************************************************
*** interactive model	 													 ***
********************************************************************************

*** initialise ddml and select model; 
ddml init late

gen lnearnings = log(earnings) 
global Y lnearnings
global D training
global Z assignmt 
global X sex-age4554

*** specify supervised machine learners for E[Y|X,Z=1] & E[Y|X,Z=0] ("yeq")
*** and E[D|X,Z=0] & E[D|X,Z=1] ("deq")

* y-equation:
ddml yeq, gen(lassoy): lasso2 $Y $X, lic(aicc) postres
ddml yeq, gen(rlassoy): rlasso $Y $X 
ddml yeq, gen(pystackedy): pystacked $Y $X, type(reg)

* d-equation:
ddml deq, gen(lassod): lasso2 $D $X, lic(aicc) postres
ddml deq, gen(rlassod): rlasso $D $X
ddml deq, gen(pystackedd): pystacked $D $X, type(reg)

* z-equation:
ddml zeq, gen(lassoz): lasso2 $Z $X, lic(aicc) postres
ddml zeq, gen(rlassoz): rlasso $Z $X
ddml zeq, gen(pystackedz): pystacked $Z $X, type(reg)

*** cross-fitting and display mean-squared prediction error
ddml crossfit, kfolds(2) tabfold

ddml desc

*** estimation of parameter of interest
ddml estimate 

*** now, do the same using one-line command
qddml $Y ($X) ($D=$Z), model(late)
 
cap log close
