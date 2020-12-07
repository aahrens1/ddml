clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	cd "/Users/kahrens/Dropbox (PP)/ddml"
}
webuse cattaneo2
//use Data/cattaneo2.dta
  
which ddml
  
********************************************************************************
*** interactive model	 													 ***
********************************************************************************

*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init interactive, mname(myest)

gen double mage_sq = mage^2

global Y bweight
global D mbsmoke
global X c.(mmarried mage mage_sq fbaby medu)#c.(mmarried mage mage_sq fbaby medu)

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(lassoy) vname($Y) mname(myest): lasso2 $Y $X, lic(aicc) postres
ddml yeq, gen(rlassoy) vname($Y) mname(myest): rlasso $Y $X 
//ddml2 yeq, gen(rfy) vname($Y) mname(myest): rforest $Y $X, type(reg)
//ddml2 yeq, gen(pfory) vname($Y) mname(myest): pyforest $Y $X, type(regress)

* d-equation:
ddml deq, gen(lassod) vname($D) mname(myest): lasso2 $D $X, lic(aicc) postres
ddml deq, gen(rigd) vname($D) mname(myest): rlasso $D $X

*** cross-fitting and display mean-squared prediction error
ddml crossfit, kfolds(2) mname(myest)

//_ddml_ate, yvar($Y) dvar($D) y0tilde(lassoy) y1tilde(rlassoy) dtilde(lassod)  
			
*** estimation of parameter of interest
ddml estimate, show(all) mname(myest)
 
pdslasso $Y $D ($X)
