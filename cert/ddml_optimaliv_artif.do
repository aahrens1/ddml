 

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
gen wall = 0
gen xall = 0
gen zall = 0
gen d1 = rnormal()
gen d2 = rnormal()
local p = 20
forvalues i= 1(1)`p' {
	gen w`i'=rnormal()
	gen x`i'=rnormal()
	gen z`i'=rnormal()
	replace wall = wall + w`i'
	replace xall = xall + x`i'
	replace d1 = d1 + x`i'*(0.9)^(`i') + w`i'*(0.9)^(`i') + z`i'*(0.9)^(21-`i')
	replace d2 = d2 + x`i'*(0.9)^(`i') + w`i'*(0.9)^(`i') + z`i'*(0.9)^(`i')
	replace zall = zall + z`i'*(0.9)^(`i')
}
gen e = rnormal() 
gen y = d1 + d2 + e + xall + 2*wall


********************************************************************************
*** iv estimation 															 ***
********************************************************************************

order y d1 d2 w* z* x*

reg y d1 d2 x1-x20 w1-w20
reg y d1 d2 x1-x20

ivreg2 y x1-x20 (d1 d2 = z1-z20), first

ivreg2 y x1-x20 (d1 = z1-z20), first

********************************************************************************
*** ddml 																	 ***
********************************************************************************
 
*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init optimaliv, mname(myest)

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(yt) mname(myest) vname(y) : lasso2 y x*, lic(aicc) postres
ddml yeq, gen(yt2) mname(myest) vname(y) : rlasso y x1-x20
	
* d-equation:
ddml deq, gen(d1t1) mname(myest) vname(d1) : lasso2 d1 x*, lic(aicc) postres
ddml deq, gen(d1t2) mname(myest) vname(d1) : rlasso d1 x1-x20
ddml deq, gen(d2t1) mname(myest) vname(d2) : rlasso d2 x*
ddml deq, gen(d2t2) mname(myest) vname(d2) : reg d2 x*

* z-equation:
ddml dheq, gen(d2H1) mname(myest) vname(d2) : rlasso d2 z* x*
ddml dheq, gen(d2H2) mname(myest) vname(d2) : reg d2 z* x*
ddml dheq, gen(d1H1) mname(myest) vname(d1) : lasso2 d1 z* x*, lic(aicc) postres
ddml dheq, gen(d1H2) mname(myest) vname(d1) : rlasso d1 z1-z20 x1-x20

		
*** cross-fitting and display mean-squared prediction error
ddml crossfit, mname(myest)

ddml desc, mname(myest)


*** estimation of parameter of interest
ddml estimate,  mname(myest)   // show(all)

