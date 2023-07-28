clear all

if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/ddml"
	adopath + "/Users/kahrens/MyProjects/pystacked"
}

cap cd "/Users/kahrens/MyProjects/ddml/cert"
cap cd "C:\LocalStore\ecomes\Documents\GitHub\ddml\cert"

cap log close
log using "ddml_cert_misc", replace text

which ddml
mata: whichddml()

use https://statalasso.github.io/dta/AJR.dta, clear

// necessary programs for cert; script exits with error if not installed
findfile pystacked.ado

set seed 123

******************************************************************************** 
**** ddml describe															****
******************************************************************************** 

global Y logpgp95
global X lat_abst edes1975 temp* humid* steplow-oilres
global D1 avexpr 
global D2 democ1 

*** pystacked, no SS

*** initialise ddml and select model; 
ddml init partial, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X , type(reg)
ddml E[D|X]: pystacked $D1 $X , type(reg)
ddml E[D|X]: pystacked $D2 $X , type(reg)

ddml describe
ddml describe, learners
ddml describe, sample
ddml describe, all
// no error, just message that results aren't available
ddml describe, crossfit
ddml describe, estimates
// crossfit
ddml crossfit
ddml describe, crossfit
ddml describe, estimates
// estimate
ddml estimate
ddml describe, all
ddml describe, crossfit
ddml describe, estimates

*** multiple learners, SS
ddml init partial, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X , type(reg)
ddml E[Y|X]: reg $Y $X 
ddml E[D|X]: pystacked $D1 $X , type(reg)
ddml E[D|X]: reg $D1 $X 

ddml describe
ddml describe, learners
ddml describe, sample
ddml describe, all
// no error, just message that results aren't available
ddml describe, crossfit
ddml describe, estimates
// crossfit
ddml crossfit
ddml describe, crossfit
ddml describe, estimates
// estimate
ddml estimate
ddml describe, all
ddml describe, crossfit
ddml describe, estimates


******************************************************************************** 
**** ddml extract															****
******************************************************************************** 

global Y logpgp95
global X lat_abst edes1975 temp* humid* steplow-oilres
global D1 avexpr 
global D2 democ1 

*** initialise ddml and select model; 
ddml init partial, kfolds(2) reps(2)
ddml E[Y|X]: pystacked $Y $X , type(reg)
ddml E[D|X]: pystacked $D1 $X , type(reg)
ddml E[D|X]: pystacked $D2 $X , type(reg)
ddml crossfit, shortstack poolstack
*** estimation of parameter of interest
ddml estimate, robust

// ddml extract - basic usage with show(.) option
ddml extract, show(pystacked)
ddml extract, show(stweights)
ddml extract, show(ssweights)
ddml extract, show(psweights)
ddml extract, show(mse)
ddml extract, show(n)

// ddml extract - details
// mStruct has two AAs:
//     (1) eqnAA with all equations; (2) estAA with estimation results
//     estAA itself has AAs. each est results is in an AA, plus some misc AAs.
// eStruct has two AAs:
//     (1) lrnAA with all learners; (2) resAA with all results for learner
// info can be extracted using Mata's AA utilities or using ddml extract
// nb: all keys are strings but ddml extract does not expect them to be in ""

// all keys for everything (lots of output if many learners)
* ddml extract, keys

// working with the estimation results in the model AA
mata: (m0.estAA).keys()		// estimation results
// AA with estimation results for spec st, resample 1
mata: (m0.estAA).get(("st","1"))
mata: AA_e2_r1 = (m0.estAA).get(("st","1"))
mata: AA_e2_r1.keys()
// estimated beta for spec st, resample 1. same output.
// key1 = spec st, key2 = resample 1, subkey1 = b, subkey2 = post
mata: AA_e2_r1.get(("b","post"))
ddml extract, mname(m0) key1(st) key2(1) subkey1(b) subkey2(post)
// save as a Mata object or a Stata object
ddml extract bvec, mname(m0) key1(st) key2(1) subkey1(b) subkey2(post)
mata: bvec
ddml extract bvec, mname(m0) key1(st) key2(1) subkey1(b) subkey2(post) stata
mat list r(bvec)
// scalar example
ddml extract D4_mse, mname(m0) key1(st) key2(1) subkey1(D4_pystacked_mse) subkey2(scalar)
mata: D4_mse
ddml extract D4_mse, mname(m0) key1(st) key2(1) subkey1(D4_pystacked_mse) subkey2(scalar) stata
di r(D4_mse)
// string/local example
ddml extract dnames, mname(m0) key1(st) key2(1) subkey1(dnames) subkey2(local)
mata: dnames
ddml extract dnames, mname(m0) key1(st) key2(1) subkey1(dnames) subkey2(local) stata
di r(dnames)

// working with the eqn AAs
// all eStructs on the mStruct AA
mata: (m0.eqnAA).keys()
// all keys for avexpr eqn
ddml extract, keys vname(avexpr)
// extract the eStruct for avexpr
ddml extract avexpr_eqn, vname(avexpr)
mata: avexpr_eqn
mata: avexpr_eqn.vname
mata: avexpr_eqn.nlearners
ddml extract, keys ename(avexpr_eqn)
// extract info from AA with all learners. same output.
// 2 keys: key1=learner, key2=main estimation string
mata: (avexpr_eqn.lrnAA).get(("D1_pystacked","est_main"))
ddml extract, ename(avexpr_eqn) key1(D1_pystacked) key2(est_main)
ddml extract, mname(m0) vname(avexpr) key1(D1_pystacked) key2(est_main)
// extract info from AA with all estimation results. same output.
// 3 keys: key1=learner, key2=output category, key3=resample number
mata: (avexpr_eqn.resAA).get(("D1_pystacked","MSE_folds","2"))
ddml extract, ename(avexpr_eqn) key1(D1_pystacked) key2(MSE_folds) key3(2)
ddml extract, mname(m0) vname(avexpr) key1(D1_pystacked) key2(MSE_folds) key3(2)
// save as a Mata object or a Stata object
ddml extract MSE_folds, mname(m0) vname(avexpr) key1(D1_pystacked) key2(MSE_folds) key3(2)
mata: MSE_folds
ddml extract MSE_folds, mname(m0) vname(avexpr) key1(D1_pystacked) key2(MSE_folds) key3(2) stata
mat list r(MSE_folds)


******************************************************************************** 
**** ddml export															****
******************************************************************************** 

// partially linear, pystacked-only, 2 reps

use https://statalasso.github.io/dta/AJR.dta, clear
global Y logpgp95
global X lat_abst edes1975 temp* humid* steplow-oilres
global D1 avexpr 
global D2 democ1 

*** initialise ddml and select model; 
ddml init partial, kfolds(2) reps(2) mname(m_export)
ddml E[Y|X], mname(m_export): pystacked $Y $X , type(reg)
ddml E[D|X], mname(m_export): pystacked $D1 $X , type(reg)
ddml E[D|X], mname(m_export): pystacked $D2 $X , type(reg)
ddml crossfit, shortstack poolstack mname(m_export)
*** estimation of parameter of interest
ddml estimate, robust mname(m_export)

gen stata_id = _n

ddml export using cert_ddml_export.csv, mname(m_export) replace addvars(stata_id)
local numvars = r(numvars)
fvexpand m_export_sample-stata_id
assert `numvars'== `: word count `r(varlist)''

// partially linear IV, pystacked-only, no SS or PS, 1 rep, prefix

use https://statalasso.github.io/dta/AJR.dta, clear
global Y logpgp95
global X edes1975 temp* humid* steplow-oilres
global D avexpr  
global Z1 lat_abst
global Z2 logem4

*** initialise ddml and select model; 
ddml init iv, kfolds(2) reps(1) mname(m_export) prefix
ddml E[Y|X], mname(m_export): pystacked $Y $X , type(reg)
ddml E[D|X], mname(m_export): pystacked $D $X , type(reg)
ddml E[Z|X], mname(m_export): pystacked $Z1 $X , type(reg)
ddml E[Z|X], mname(m_export): pystacked $Z2 $X , type(reg)
ddml crossfit, mname(m_export)
ddml estimate, mname(m_export)

gen stata_id = _n

ddml export using cert_ddml_export.csv, mname(m_export) replace addvars(stata_id)
local vreplist = r(vreplist)
local numvars = r(numvars)
local vreplist : subinstr local vreplist "m_export_" "m_export_", all count(local nvreplist)
// number of vars = total minus sample, 2 fold varis and stata_id=4
assert `nvreplist' == `numvars'-4
fvexpand m_export_sample-stata_id
assert `numvars' == `: word count `r(varlist)''

// interactive model, pystacked-only

webuse cattaneo2, clear
keep in 1/1000
gen double mage_sq = mage^2

global Y bweight
global D mbsmoke
global X c.(mmarried mage mage_sq fbaby medu)#c.(mmarried mage mage_sq fbaby medu)

ddml init interactive, kfolds(2) reps(2) mname(m_export)
ddml E[Y|X,D], mname(m_export): pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[D|X], mname(m_export): pystacked $D $X, type(class) method(logit gradboost)
ddml crossfit, shortstack poolstack mname(m_export)
ddml estimate, mname(m_export)
gen stata_id = _n

ddml export using cert_ddml_export.csv, mname(m_export) replace addvars(stata_id)
local numvars = r(numvars)
fvexpand m_export_sample-stata_id
assert `numvars'== `: word count `r(varlist)''

// interactive model, multiple learners

webuse cattaneo2, clear
keep in 1/1000
gen double mage_sq = mage^2

global Y bweight
global D mbsmoke
global X c.(mmarried mage mage_sq fbaby medu)#c.(mmarried mage mage_sq fbaby medu)

ddml init interactive, kfolds(2) reps(2) mname(m_export)
ddml E[Y|X,D], mname(m_export): pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[Y|X,D], mname(m_export): reg $Y $X
ddml E[D|X], mname(m_export): pystacked $D $X, type(class) method(logit gradboost)
ddml E[D|X], mname(m_export): logit $D $X
ddml crossfit, shortstack mname(m_export)
ddml estimate, mname(m_export)
gen stata_id = _n

ddml export using cert_ddml_export.csv, mname(m_export) replace addvars(stata_id)
local numvars = r(numvars)
fvexpand m_export_sample-stata_id
assert `numvars'== `: word count `r(varlist)''

// interactiveiv (LATE), pystacked-only

use "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta",clear   
keep in 1/5000
gen lnearnings = log(earnings) 

global Y lnearnings
global D training
global Z assignmt 
global X sex-age4554

ddml init interactiveiv, kfolds(2) reps(2) mname(m_export)
ddml E[Y|X,Z], mname(m_export): pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[D|X,Z], mname(m_export): pystacked $D $X, type(class) method(logit gradboost)
ddml E[Z|X], mname(m_export): pystacked $Z $X, type(class) method(logit gradboost)
ddml crossfit, mname(m_export)
ddml estimate, mname(m_export)
gen stata_id = _n

ddml export using cert_ddml_export.csv, mname(m_export) replace addvars(stata_id)
local numvars = r(numvars)
fvexpand m_export_sample-stata_id
assert `numvars'== `: word count `r(varlist)''

// interactiveiv (LATE), multiple learners

use "http://fmwww.bc.edu/repec/bocode/j/jtpa.dta",clear   
keep in 1/5000
gen lnearnings = log(earnings) 

global Y lnearnings
global D training
global Z assignmt 
global X sex-age4554

ddml init interactiveiv, kfolds(2) reps(2) mname(m_export)
ddml E[Y|X,Z], mname(m_export): pystacked $Y $X, type(reg) method(ols gradboost)
ddml E[Y|X,Z], mname(m_export): reg $Y $X
ddml E[D|X,Z], mname(m_export): pystacked $D $X, type(class) method(logit gradboost)
ddml E[D|X,Z], mname(m_export): logit $D $X
ddml E[Z|X], mname(m_export): pystacked $Z $X, type(class) method(logit gradboost)
ddml E[Z|X], mname(m_export): logit $Z $X
ddml crossfit, shortstack mname(m_export)
ddml estimate, mname(m_export)
gen stata_id = _n

ddml export using cert_ddml_export.csv, mname(m_export) replace addvars(stata_id)
local numvars = r(numvars)
fvexpand m_export_sample-stata_id
assert `numvars'== `: word count `r(varlist)''

// flexible iv, pystacked-only

use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear

global Y share
global D price
global X hpwt air mpd space
global Z sum*

ddml init fiv, kfolds(2) reps(2) mname(m_export)
ddml E[Y|X], mname(m_export): pystacked $Y $X, type(reg)
ddml E[D|Z,X], learner(Dhat_pystacked) mname(m_export): pystacked $D $X $Z, type(reg)
ddml E[D|X], learner(Dhat_pystacked) vname($D) mname(m_export): pystacked {D} $X, type(reg)
ddml crossfit, mname(m_export)
ddml estimate, mname(m_export)
gen stata_id = _n

ddml export using cert_ddml_export.csv, mname(m_export) replace addvars(stata_id)
local numvars = r(numvars)
fvexpand m_export_sample-stata_id
assert `numvars'== `: word count `r(varlist)''

log close
