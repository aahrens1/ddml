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

global Y logpgp95
global X edes1975 avelf temp* humid* steplow-oilres
global D1 avexpr 
global D2 lat_abst

*** initialise ddml and select model; 
* currently only the partial linear model is supported
ddml init partial 

*** specify supervised machine learners for E[Y|X] ("yeq") and E[D|X] ("deq")
* y-equation:
ddml yeq, gen(lassoy): lasso2 $Y $X, lic(aicc) postres
ddml yeq, gen(rigy): rlasso $Y $X
ddml yeq, gen(pystackedy): pystacked $Y $X, type(regress) 

* d-equation:
ddml deq, gen(lassod1): lasso2 $D1 $X, lic(aicc) postres
ddml deq, gen(rigd1): rlasso $D1 $X
ddml deq, gen(pys1): pystacked $D1 $X, type(reg) 
ddml deq, gen(lassod2): lasso2 $D2 $X, lic(aicc) postres
ddml deq, gen(rigd2): rlasso $D2 $X
ddml deq, gen(pystackedd): pystacked $D2 $X, type(regress)  

ddml sample if democ1 < ., vars($Y $D1 $D2 $X) 

*** cross-fitting and display mean-squared prediction error
ddml crossfit 

ddml desc 

*** estimation of parameter of interest
ddml estimate 
