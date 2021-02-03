
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
ddml init optimaliv 

*** specify supervised machine learners for E[Y|X] ("yeq"), E[D|X] ("deq")
*** as well as E[D|ZX] ("dheq")

* y-equation:
ddml yeq, gen(yt1) : lasso2 y x1-x20, lic(aicc) postres
ddml yeq, gen(yt2) : rlasso y x1-x20
ddml yeq, gen(yt3) : pystacked y x1-x20, type(reg)

* d-equation:
ddml deq, gen(d1t1): lasso2 d1 x1-x20, lic(aicc) postres
ddml deq, gen(d1t2): rlasso d1 x1-x20
ddml deq, gen(d1t3): pystacked d1 x1-x20, type(reg)
ddml deq, gen(d2t1): rlasso d2 x1-x20
ddml deq, gen(d2t2): reg d2 x1-x20
ddml deq, gen(d2t3): pystacked d2 x1-x20, type(reg)

* dh-equation:
ddml dheq, gen(d2H1): rlasso d2 z1-z20 x1-x20
ddml dheq, gen(d2H2): reg d2 z1-z20 x1-x20
ddml dheq, gen(d2H3): pystacked d2 z1-z20 x1-x20, type(reg)
ddml dheq, gen(d1H1): lasso2 d1 z1-z20 x1-x20, lic(aicc) postres
ddml dheq, gen(d1H2): rlasso d1 z1-z20 x1-x20
ddml dheq, gen(d1H3): pystacked d1 z1-z20 x1-x20, type(reg)

*** cross-fitting and display mean-squared prediction error
ddml crossfit

ddml desc

*** estimation of parameter of interest
ddml estimate 
