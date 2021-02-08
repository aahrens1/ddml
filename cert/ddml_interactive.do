clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pylearn2"
	cd "/Users/kahrens/Dropbox (PP)/ddml"
}
webuse cattaneo2
//use Data/cattaneo2.dta
  
which ddml
  
********************************************************************************
*** interactive model	 													 ***
********************************************************************************

*** initialise ddml and select model; 
ddml init interactive

gen double mage_sq = mage^2

global Y bweight
global D mbsmoke
global X c.(mmarried mage mage_sq fbaby medu)#c.(mmarried mage mage_sq fbaby medu)

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(lassoy): lasso2 $Y $X, lic(aicc) postres
ddml yeq, gen(rlassoy): rlasso $Y $X 
ddml yeq, gen(pystackedy): pystacked $Y $X, type(reg) 

* d-equation:
ddml deq, gen(lassod): lasso2 $D $X, lic(aicc) postres
ddml deq, gen(rlassod): rlasso $D $X
ddml deq, gen(pystackedd): pystacked $D $X, type(reg)

*** cross-fitting and display mean-squared prediction error
ddml crossfit, kfolds(2) 

ddml desc

*** estimation of parameter of interest
ddml estimate, show(all) 

*** now, do the same using one-line command
qddml $Y $D ($X), model(interactive)
