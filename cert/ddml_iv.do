 

clear all

use Data/AJR.dta // https://statalasso.github.io/dta/AJR.dta
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	cd "/Users/kahrens/Dropbox (PP)/ddml"
}
  
which ddml

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
ddml yeq, gen(pystackedy): pystacked $Y $X, type(reg)

* d-equation:
ddml deq, gen(lassod1): lasso2 $D $X, lic(aicc) postres
ddml deq, gen(rigd1): rlasso $D $X
ddml deq, gen(pystackedd): pystacked $D $X, type(reg)

* z-equation:
ddml zeq, gen(lassoz): lasso2 $Z $X, lic(aicc) postres
ddml zeq, gen(pystackedd): pystacked $Z $X, type(reg)
	
*** cross-fitting and display mean-squared prediction error
ddml crossfit

ddml desc

*** estimation of parameter of interest
ddml estimate
