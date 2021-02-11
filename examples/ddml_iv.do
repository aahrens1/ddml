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
log using "log_iv.txt", replace text  
 
use https://statalasso.github.io/dta/AJR.dta
  
which ddml
which pystacked

********************************************************************************
*** iv model																 ***
********************************************************************************

global Y logpgp95
global X edes1975 avelf temp* humid* steplow-oilres
global D avexpr 
global Z logem4

*** initialise ddml and select model
ddml init iv 

*** specify supervised machine learners for E[Y|X] ("yeq"), E[D|X] ("deq")
*** as well as E[Z|X] ("zeq")
* y-equation:
ddml yeq, gen(lassoy): lasso2 $Y $X, lic(aicc) postres
ddml yeq, gen(rigy): rlasso $Y $X
ddml yeq, gen(pystackedy): pystacked $Y $X, type(reg) method(rf gradboost)

* d-equation:
ddml deq, gen(lassod1): lasso2 $D $X, lic(aicc) postres
ddml deq, gen(rigd1): rlasso $D $X
ddml deq, gen(pystackedd): pystacked $D $X, type(reg) method(rf gradboost)

* z-equation:
ddml zeq, gen(lassoz): lasso2 $Z $X, lic(aicc) postres
ddml zeq, gen(pystackedd): pystacked $Z $X, type(reg) method(rf gradboost)
	
*** cross-fitting and display mean-squared prediction error
ddml crossfit
ddml desc

*** estimation of parameter of interest
ddml estimate

*** now, using one-line command:
qddml $Y ($X) ($D = $Z), model(iv) cmdopt(method(rf gradboost))

cap log close
