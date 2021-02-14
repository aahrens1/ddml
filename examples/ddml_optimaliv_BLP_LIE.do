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
log using "log_BLP_LIE.txt", replace text  
 
which ddml
which pystacked 

insheet using "https://statalasso.github.io/dta/BLP1995-prepared.csv", comma clear

********************************************************************************
*** LIE-conform DDML-IV														 ***
********************************************************************************

global Y share
global D price
global X hpwt-spaceu_tu
global Z z1-z28
 
*** initialise ddml and select model
set seed 123
ddml init optimaliv 

*** specify supervised machine learners for E[Y|X] ("yeq"), E[D^|X] ("deq")
*** as well as D^=E[D|Z,X] ("dheq")

* "y"-equation:
ddml yeq, gen(yt1) : rlasso $Y $X
ddml yeq, gen(yt2) : pystacked $Y $X, type(reg) seed(99)
ddml yeq, gen(yt3) : pyvote $Y $X, type(reg)

* "dheq"-equation:
ddml dheq, gen(d2H1) noprefix: rlasso $D $X $Z
ddml dheq, gen(d2H2) noprefix: pystacked $D $X $Z, type(reg)
ddml dheq, gen(d2H3) noprefix: pyvote $D $X $Z, type(reg)
  
* "d"-equation:
/* In this step, we use estimate of E[D|X,Z] from "dheq"-step as LHS variable.
	Notes:
		- We need to specify variable name explicitly ($D) as it doesn't 
		correspond to LHS variable (d2H1 etc).
		- Previous step uses "noprefix" option as by default a prefix
		is added to all predicted values. 
		- Currently "dheq" needs to be added before "deq"; otherwise an error
		occurs. (Will be fixed.)
*/ 
ddml deq, gen(d1t1) vname($D): rlasso d2H1 $X
ddml deq, gen(d1t2) vname($D): pystacked d2H2 $X, type(reg)
ddml deq, gen(d1t3) vname($D): pyvote d2H3 $X, type(reg) 
 
*** cross-fitting and display mean-squared prediction error
ddml desc
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
