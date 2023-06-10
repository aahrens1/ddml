clear all

if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
}

cap cd "/Users/kahrens/MyProjects/ddml/cert"
cap cd "C:\LocalStore\ecomes\Documents\GitHub\ddml\cert"

cap log close
log using "ddml_cert_crossfit", replace text

which ddml
mata: whichddml()
which crossfit

use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta, clear
global X sex age married black hispanic

set seed 42
crossfit, estring(reg earnings $X) gen(yhat) kfolds(3)
sum earnings yhat_1

set seed 42
crossfit, estring(pystacked earnings $X) gen(yhat) kfolds(3)
sum earnings yhat*

set seed 42
crossfit, estring(reg earnings $X) gen(yhat) kfolds(3) reps(5)
sum earnings yhat*

// check that norandom is equivalent to provided fold identifier
count
gen fid = _n<=(r(N)/2)
set seed 42
crossfit, estring(reg earnings $X) gen(noran) kfolds(2) norandom
set seed 42
crossfit, estring(reg earnings $X) gen(foldv) foldvar(fid)
assert noran==foldv

log close
