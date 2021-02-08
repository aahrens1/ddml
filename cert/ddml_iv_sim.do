 
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
gen z1 = rnormal()
gen z2 = rnormal()
gen w = rnormal()
gen xall = 0
forvalues i= 1(1)20 {
	gen x`i'=rnormal()
	replace xall = xall + x`i'*(0.9)^(`i')
}
gen d1 = rnormal() + z1*0.8 + z2*0.2  + w + xall
gen d2 = rnormal() + z1*0.2 + z2*0.8  + w + xall
gen e = rnormal() 
gen y = w + d1 + d2 + e + xall

reg y d* x*

********************************************************************************
*** iv model																 ***
********************************************************************************

 
*** initialise ddml and select model; 
ddml init iv

*** specify supervised machine learners for E[Y|X] ("yeq"), E[D|X] ("deq")
*** as well as E[Z|X] ("zeq")* y-equation:
ddml yeq, gen(yt1): lasso2 y x*, lic(aicc) postres
ddml yeq, gen(yt2): rlasso y x*
ddml yeq, gen(yt3): pystacked y x*, type(reg)
	
* d-equation:
ddml deq, gen(d1t1): lasso2 d1 x*, lic(aicc) postres
ddml deq, gen(d1t2): rlasso d1 x*
ddml deq, gen(d1t3): pystacked d1 x*, type(reg)
ddml deq, gen(d2t1): rlasso d2 x*
ddml deq, gen(d2t2): reg d2 x*
ddml deq, gen(d2t2): pystacked d2 x*, type(reg)

* z-equation:
ddml zeq, gen(z1t1): lasso2 z1 x*, lic(aicc) postres
ddml zeq, gen(z1t2): rlasso z1 x*
ddml zeq, gen(z1t3): pystacked z1 x*, type(reg)
ddml zeq, gen(z2t1): lasso2 z2 x*, lic(aicc) postres
ddml zeq, gen(z2t2): rlasso z2 x*
ddml zeq, gen(z2t3): pystacked z2 x*, type(reg)

*** cross-fitting and display mean-squared prediction error
ddml crossfit 
	
ddml desc 

*** estimation of parameter of interest
ddml estimate // show(all)

*** now, using one-line command:
qddml y (x*) (d1 d2 = z1 z2), model(iv)
