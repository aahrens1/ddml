 

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
 

*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init iv, mname(myest)

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(lassoy) mname(myest) vname(logpgp95) : ///
		lasso2 logpgp95 edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml yeq, gen(rigy) mname(myest) vname(logpgp95) : ///
		rlasso logpgp95 edes1975 avelf temp* humid* steplow-oilres
	
* d-equation:
ddml deq, gen(lassod1) mname(myest) vname(avexpr) : lasso2 avexpr edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml deq, gen(rigd1) mname(myest) vname(avexpr) : rlasso avexpr edes1975 avelf temp* humid* steplow-oilres
ddml deq, gen(rigdemo) mname(myest) vname(democ00a) : rlasso democ00a edes1975 avelf temp* humid* steplow-oilres
ddml deq, gen(rigdemo2) mname(myest) vname(democ00a) : reg democ00a edes1975 avelf temp* humid* steplow-oilres

* z-equation:
ddml zeq, gen(lassoz) vname(logem4) mname(myest): ///
	lasso2 logem4   temp* humid* steplow-oilres, lic(aicc) postres
ddml zeq, gen(rigz) vname(logem4) mname(myest): ///
	rlasso logem4  temp* humid* steplow-oilres
 
ddml desc, mname(myest)

	
*** cross-fitting and display mean-squared prediction error
ddml crossfit, mname(myest)

*** estimation of parameter of interest
ddml estimate,  mname(myest) // show(all)
