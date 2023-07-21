clear all
 
if ("`c(username)'"=="kahrens") {
	adopath + "/Users/kahrens/MyProjects/pystacked"
	cd "/Users/kahrens/MyProjects/ddml/cert"
}

cap log close
log using "qddml_cert", replace text

* global tol = 0.0001
which ddml
which pystacked

**** Partially linear model w/pystacked integration

use https://github.com/aahrens1/ddml/raw/master/data/sipp1991.dta, clear
global Y net_tfa
global D e401
global X tw age inc fsize educ db marr twoearn pira hown

set seed 42
ddml init partial, kfolds(2)
ddml E[Y|X]: pystacked $Y $X, type(reg)
ddml E[D|X]: pystacked $D $X, type(reg)
// ddml reports only standard stacking by default
ddml crossfit, shortstack
ddml estimate
global b1_ss = _b[$D]
ddml estimate, mname(m0) spec(st) rep(1) notable replay
global b1_std = _b[$D]

set seed 42
// qddml reports only short-stacking by default
qddml $Y $D ($X), kfolds(2) model(partial) stdstack shortstack mname(m0q)
global b2_ss = _b[$D]
ddml estimate, mname(m0q) spec(st) rep(1) notable replay
global b2_std = _b[$D]

assert reldif($b1_ss, $b2_ss)  < 10e-10
assert reldif($b1_std,$b2_std) < 10e-10

**** Partially linear IV model, rforest

use https://statalasso.github.io/dta/AJR.dta, clear
global Y logpgp95
global D avexpr
global Z logem4
global X lat_abst edes1975 avelf temp* humid* steplow-oilres

set seed 42
ddml init iv, kfolds(30)
ddml E[Y|X], vtype(none): rforest $Y $X, type(reg)
ddml E[D|X], vtype(none): rforest $D $X, type(reg)
ddml E[Z|X], vtype(none): rforest $Z $X, type(reg)
ddml crossfit
ddml estimate, robust
global b1 = _b[$D]

set seed 42
qddml $Y ($X) ($D=$Z), kfolds(30) model(iv) mname(m0q)		///
	cmd(rforest) cmdopt(type(reg)) vtype(none) robust
global b2 = _b[$D]

assert reldif($b1,$b2) < 10e-10

**** Interactive model--ATE and ATET estimation w/pystacked integration

webuse cattaneo2, clear
keep in 1/1000
global Y bweight
global D mbsmoke
global X mage prenatal1 mmarried fbaby mage medu
global pystacked_y_options type(reg) method(rf gradboost)
global pystacked_d_options type(class) method(rf gradboost)

set seed 42
ddml init interactive, kfolds(2) reps(2)
ddml E[Y|X,D]: pystacked $Y $X, $pystacked_y_options
ddml E[D|X]: pystacked $D $X, $pystacked_d_options
// ddml reports only standard stacking by default
ddml crossfit, shortstack
ddml estimate
global b1_md_ss			= _b[$D]
ddml estimate, mname(m0) spec(st) rep(md) notable replay
global b1_md_std		= _b[$D]
ddml estimate, atet
global b1_atet_md_ss	= _b[$D]
ddml estimate, mname(m0) spec(st) rep(md) notable replay
global b1_atet_md_std	= _b[$D]

set seed 42
// qddml reports only short-stacking by default
qddml $Y $D ($X), kfolds(2) reps(2) model(interactive)		///
	mname(m0q)												///
	stdstack shortstack										///
	pystacked_y($pystacked_y_options)						///
	pystacked_d($pystacked_d_options)
global b2_md_ss			= _b[$D]
ddml estimate, mname(m0q) spec(st) rep(md) notable replay
global b2_md_std		= _b[$D]
ddml estimate, atet mname(m0q)
global b2_atet_md_ss	= _b[$D]
ddml estimate, mname(m0q) spec(st) rep(md) notable replay
global b2_atet_md_std	= _b[$D]

assert reldif($b1_md_ss,       $b2_md_ss)       < 10e-10
assert reldif($b1_md_std,      $b2_md_std)      < 10e-10
assert reldif($b1_atet_md_ss,  $b2_atet_md_ss)  < 10e-10
assert reldif($b1_atet_md_std, $b2_atet_md_std) < 10e-10

**** Interactive IV model w/pystacked integration

use http://fmwww.bc.edu/repec/bocode/j/jtpa.dta,clear
global Y earnings
global D training
global Z assignmt
global X sex age married black hispanic
global pystacked_y_options type(reg) m(lassocv)
global pystacked_d_options type(class) m(lassocv)
global pystacked_z_options type(class) m(lassocv)

set seed 42
ddml init interactiveiv, kfolds(5)
ddml E[Y|X,Z]: pystacked $Y c.($X)# #c($X), $pystacked_y_options
ddml E[D|X,Z]: pystacked $D c.($X)# #c($X), $pystacked_d_options
ddml E[Z|X]: pystacked $Z c.($X)# #c($X), $pystacked_z_options
ddml crossfit
ddml estimate
global b1 = _b[$D]

set seed 42
qddml $Y (c.($X)# #c($X)) ($D=$Z), kfolds(5) model(interactiveiv)	///
	mname(m0q)														///
	pystacked_y($pystacked_y_options)								///
	pystacked_d($pystacked_d_options)								///
	pystacked_z($pystacked_z_options)
global b2 = _b[$D]

assert reldif($b1,$b2) < 10e-10

**** FIV

use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear
global Y share
global D price
global X hpwt air mpd space
global Z sum*

// standard stacking by pystacked
set seed 42
ddml init fiv
ddml E[Y|X]: pystacked $Y $X, type(reg)
ddml E[D|Z,X], learner(Dhat_pystacked): pystacked $D $X $Z, type(reg)
ddml E[D|X], learner(Dhat_pystacked) vname($D): pystacked {D} $X, type(reg)
ddml crossfit
ddml estimate
global b1 = _b[$D]

set seed 42
// enforce standard stacking using the cmd(.) option
qddml $Y ($X) ($D=$Z), model(fiv) cmd(pystacked) cmdopt(type(reg))
global b2 = _b[$D]

assert reldif($b1,$b2) < 10e-10

ddml extract, show(pystacked)
ddml extract, show(stweights)

// short-stacking only
use https://github.com/aahrens1/ddml/raw/master/data/BLP.dta, clear
global Y share
global D price
global X hpwt air mpd space
global Z sum*

set seed 42
ddml init fiv
ddml E[Y|X], learner(Y1_ols): pystacked $Y $X, type(reg) m(ols)
ddml E[Y|X], learner(Y2_lassocv): pystacked $Y $X, type(reg) m(lassocv)
ddml E[Y|X], learner(Y3_gradboost): pystacked $Y $X, type(reg) m(gradboost)
ddml E[D|Z,X], learner(D1_ols): pystacked $D $X $Z, type(reg) m(ols)
ddml E[D|X], learner(D1_ols) vname($D): pystacked {D} $X, type(reg) m(ols)
ddml E[D|Z,X], learner(D2_lassocv): pystacked $D $X $Z, type(reg) m(lassocv)
ddml E[D|X], learner(D2_lassocv) vname($D): pystacked {D} $X, type(reg) m(lassocv)
ddml E[D|Z,X], learner(D3_gradboost): pystacked $D $X $Z, type(reg) m(gradboost)
ddml E[D|X], learner(D3_gradboost) vname($D): pystacked {D} $X, type(reg) m(gradboost)
ddml crossfit, shortstack
ddml estimate
global b1 = _b[$D]
ddml extract, show(ssweights)
mat Yss1 = r(Y_share_ss)
mat Dss1 = r(D_price_ss)

set seed 42
// short-stacking only by default
qddml $Y ($X) ($D=$Z), model(fiv)
global b2 = _b[$D]
ddml extract, show(ssweights)
mat Yss2 = r(Y_share_ss)
mat Dss2 = r(D_price_ss)

assert reldif($b1,$b2) < 10e-10
assert mreldif(Yss1,Yss2) < 10e-10
assert mreldif(Dss1,Dss2) < 10e-10

log close
