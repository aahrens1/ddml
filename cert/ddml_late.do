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
ddml init late, mname(myest)

gen lnearnings = log(earnings) 
global Y lnearnings
global D training
global Z assignmt 
global X sex-age4554

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(lassoy) vname($Y) mname(myest): lasso2 $Y $X, lic(aicc) postres
ddml yeq, gen(rlassoy) vname($Y) mname(myest): rlasso $Y $X 
//ddml yeq pfory: (pyadaboost $Y $X, type(regress))

* d-equation:
ddml deq, gen(lassod) vname($D) mname(myest): lasso2 $D $X, lic(aicc) postres
ddml deq, gen(rlassod) vname($D) mname(myest):  rlasso $D $X

* z-equation:
ddml zeq, gen(lassoz) vname($Z) mname(myest): lasso2 $Z $X, lic(aicc) postres
ddml zeq , gen(rlassoz) vname($Z) mname(myest): rlasso $Z $X

*** cross-fitting and display mean-squared prediction error
ddml crossfit, kfolds(2) mname(myest) tabfold

//_ddml_late, yvar($Y) dvar($D) zvar($Z) y0tilde(lassoy) y1tilde(lassoy) ///
//					d0tilde(lassod) d1tilde(lassod) ztilde(lassoz)

*** estimation of parameter of interest
ddml estimate, show(all) mname(myest)
 
