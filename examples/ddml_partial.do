clear all

if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pylearn2"
	cd "/Users/kahrens/MyProjects/ddml/examples"
}

cap log close
log using "log_partial.txt", replace text  
 
use https://statalasso.github.io/dta/AJR.dta, clear

which ddml
which pystacked
 
********************************************************************************
*** partial model														     ***
********************************************************************************

global Y logpgp95
global X lat_abst edes1975 avelf temp* humid* steplow-oilres
global D avexpr 

*** initialise ddml and select model; 
ddml init partial 

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(lassoy): lasso2 $Y $X, lic(aicc) postres
ddml yeq, gen(rigy): rlasso $Y $X
ddml yeq, gen(pystackedy): pystacked $Y $X, type(regress) method(rf gradboost)

* d-equation:
ddml deq, gen(lassod1): lasso2 $D $X, lic(aicc) postres
ddml deq, gen(rigd1): rlasso $D $X
ddml deq, gen(pys1): pystacked $D $X, type(reg) method(rf gradboost)

*** cross-fitting and display mean-squared prediction error
ddml crossfit 
ddml desc 

*** estimation of parameter of interest
ddml estimate 

*** now using the one-line command
qddml $Y $D ($X), model(partial) cmdopt(method(rf gradboost))

cap log close
