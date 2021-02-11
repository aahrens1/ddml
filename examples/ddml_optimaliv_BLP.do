clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pylearn2"
	cd "/Users/kahrens/MyProjects/ddml/examples"
}

** Note that you need to run the file "lddml.do" to install
** required Mata structures.
** This only needs to be done once.
//run "/Users/kahrens/MyProjects/ddml/lddml.do"

cap log close
log using "log_BLP.txt", replace text  
 
which ddml
which pystacked 

insheet using "https://statalasso.github.io/dta/BLP1995-prepared.csv", comma clear

********************************************************************************
*** ddml 																	 ***
********************************************************************************

global Y share
global D price
global X hpwt-spaceu_tu
global Z z1-z28
 
*** initialise ddml and select model
set seed 123
ddml init optimaliv 

*** specify supervised machine learners for E[Y|X] ("yeq"), E[D|X] ("deq")
*** as well as E[D|Z,X] ("zeq")

* "y"-equation:
ddml yeq, gen(yt1) : rlasso $Y $X
ddml yeq, gen(yt2) : pystacked $Y $X, type(reg) seed(99)
ddml yeq, gen(yt3) : pyvote $Y $X, type(reg)

* "d"-equation:
ddml deq, gen(d1t1): rlasso $D $X
ddml deq, gen(d1t2): pystacked $D $X, type(reg)
ddml deq, gen(d1t3): pyvote $D $X, type(reg)

* "dheq"-equation:
ddml dheq, gen(d2H1): rlasso $D $X $Z
ddml dheq, gen(d2H2): pystacked $D $X $Z, type(reg)
ddml dheq, gen(d2H3): pyvote $D $X $Z, type(reg)
 
*** cross-fitting and display mean-squared prediction error
ddml crossfit

*** get an overview of "equations" and MSE
ddml desc

*** estimation of parameter of interest
ddml estimate 

*** now, do the same using the one-line command (qddml) 
*** .. which uses only pystacked
qddml $Y ($X) ($D = $Z), model(optimaliv)  
*** use only rlasso:
qddml $Y ($X) ($D = $Z), model(optimaliv) cmd(rlasso)

cap log close
