clear all

if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	cd "/Users/kahrens/Dropbox (PP)/ddml"
}
// use Data/AJR.dta
use https://statalasso.github.io/dta/AJR.dta, clear


which ddml
 
********************************************************************************
*** partial model														     ***
********************************************************************************

*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init partial, mname(myest)

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(lassoy) mname(myest) vname(logpgp95) : lasso2 logpgp95 edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml yeq, gen(rigy) mname(myest) vname(logpgp95) : rlasso logpgp95 edes1975 avelf temp* humid* steplow-oilres
	
* d-equation:
ddml deq, gen(lassod1) mname(myest) vname(avexpr) : lasso2 avexpr edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml deq, gen(rigd1) mname(myest) vname(lat_abst) : rlasso avexpr edes1975 avelf temp* humid* steplow-oilres
ddml deq, gen(lassod2) mname(myest) vname(avexpr) : lasso2 lat_abst edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml deq, gen(rigd2) mname(myest) vname(lat_abst) : rlasso lat_abst edes1975 avelf temp* humid* steplow-oilres

ddml sample if democ1 < ., mname(myest) vars(logpgp95 avexpr lat_abst edes1975 avelf temp* humid* steplow-oilres)

ddml desc, mname(myest)

*** cross-fitting and display mean-squared prediction error
ddml crossfit, mname(myest)

*** estimation of parameter of interest
ddml estimate, mname(myest)
