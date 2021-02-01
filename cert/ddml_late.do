clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	cd "/Users/kahrens/Dropbox (PP)/ddml"
}
//use Data/jtpa.dta //
use "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta",clear   

which ddml
  
********************************************************************************
*** interactive model	 													 ***
********************************************************************************

*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init late

gen lnearnings = log(earnings) 
global Y lnearnings
global D training
global Z assignmt 
global X sex-age4554

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")

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
ddml estimate, show(all) 
 
