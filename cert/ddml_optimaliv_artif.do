 

clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	cd "/Users/kahrens/Dropbox (PP)/ddml"
}
  
which ddml

********************************************************************************
*** gen data																 ***
********************************************************************************

set obs 1000
gen w = rnormal()
gen xall = 0
gen zall = 0
forvalues i= 1(1)20 {
	gen x`i'=rnormal()
	replace xall = xall + x`i'*(0.9)^(`i')
}
forvalues i= 1(1)20 {
	gen z`i'=rnormal()
	replace zall = zall + z`i'*(0.9)^(`i')
}
gen d1 = rnormal() + zall + w +xall
gen d2 = rnormal() + zall + w + xall
gen e = rnormal() 
gen y = w + d1 + d2 + e + xall

reg y d* x*

********************************************************************************
*** iv model																 ***
********************************************************************************

 
*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init optimaliv, mname(myest)

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(yt) mname(myest) vname(y) : lasso2 y x*, lic(aicc) postres
ddml yeq, gen(yt2) mname(myest) vname(y) : rlasso y x*
	
* d-equation:
ddml deq, gen(d1t1) mname(myest) vname(d1) : lasso2 d1 x*, lic(aicc) postres
ddml deq, gen(d1t2) mname(myest) vname(d1) : rlasso d1 x*
ddml deq, gen(d2t1) mname(myest) vname(d2) : rlasso d2 x*
ddml deq, gen(d2t2) mname(myest) vname(d2) : reg d2 x*

* z-equation:
ddml dheq, gen(d1H1) mname(myest) vname(d1) : lasso2 d1 z* x*, lic(aicc) postres
ddml dheq, gen(d1H2) mname(myest) vname(d1) : rlasso d1 z* x*
ddml dheq, gen(d2H1) mname(myest) vname(d2) : rlasso d2 z* x*
ddml dheq, gen(d2H2) mname(myest) vname(d2) : reg d2 z* x*
	
ddml desc, mname(myest)
	
*** cross-fitting and display mean-squared prediction error
ddml crossfit, mname(myest)

*** estimation of parameter of interest
ddml estimate,  mname(myest) // show(all)
