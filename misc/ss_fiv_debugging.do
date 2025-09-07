cap cd "C:\LocalStore\ecomes\GitHub\ddml\misc"

cap log close
log using debug_log_`c(hostname)'.txt, text replace

clear all

which ddml, all
which crossfit, all

global reps 1
global kfolds 3

cap prog drop crossfit
cap prog drop _ddml_estimate_linear
cap prog drop _ddml_crossfit
cap prog drop ddml
cap prog drop _ddml_extract

// want to control any randomization by learners
// since non-pystacked and pystacked code order differs
global pyopt norandom noshuffle
//global L1 ols
//global opt1
global L1 rf
global opt1 random_state(1) max_features(1)
//global L2 ols
//global L2 lassoic
global L2 ridgecv
global opt2 random_state(1)


sysuse auto, clear
drop if rep78==.
gen byte one=1
global Y mpg
global D price
global Z rep78-displacement
global X gear_ratio foreign



**************** H learner is regress, no Xs, LIE enforced *******************

cap drop ps*
cap drop nps*

// Non-pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(nps) prefix
ddml E[Y|X], mname(nps)								: reg $Y
ddml E[D|Z,X], learner(D1) mname(nps)				: pystacked $D $Z, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D1) vname($D) mname(nps)		: reg {D}
ddml E[D|Z,X], learner(D2) mname(nps)				: pystacked $D $Z, method($L2) $pyopt
ddml E[D|X],   learner(D2) vname($D) mname(nps)		: reg {D}
// force branch to non-pystacked routine
ddml crossfit, mname(nps) crossfitother shortstack
ddml estimate, mname(nps)
mat b_nps=e(b)
mat V_nps=e(V)

// pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(ps) prefix
ddml E[Y|X], mname(ps)								: reg $Y
ddml E[D|Z,X], learner(D) mname(ps)					: pystacked $D $Z, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D) vname($D) mname(ps)		: reg {D}
ddml crossfit, mname(ps) shortstack 
ddml estimate, mname(ps)
mat b_ps=e(b)
mat V_ps=e(V)

// estimation results
assert mreldif(b_nps,b_ps)<10e-7
assert mreldif(V_nps,V_ps)<10e-7

// ss weights
ddml extract, show(ssweights) mname(nps)
mat ssw_nps=r(nps_D_price_ss)
ddml extract, show(ssweights) mname(ps)
mat ssw_ps=r(ps_D_price_ss)
// short-stacking h is redundant because the two h learners to be stacked
// by the non-pystacked code are identical:
assert nps_D1_h_ss_1==nps_D2_h_ss_1
// and the nnls solution with the ssw_h weights is arbitrary.
// the pystacked code anticipates this and skips the nnls.
// hence we check only first two (non-h) rows.
mat ssw_nps = ssw_nps[1..2,.]
mat ssw_ps = ssw_ps[1..2,.]
assert mreldif(ssw_nps,ssw_ps)<10e-7

// Y
sum nps_Y1_reg_1
sum ps_Y1_reg_1
assert reldif(nps_Y1_reg_1,ps_Y1_reg_1)<10e-7

// D learners for E[D|X,Z] - predicted values match
sum nps_D1_1 ps_D_L1_1
sum nps_D2_1 ps_D_L2_1
assert reldif(nps_D1_1,ps_D_L1_1)<10e-7
assert reldif(nps_D2_1,ps_D_L2_1)<10e-7

// h learner is regress for all D learners, so predicted values match w/LIE
sum nps_D1_h_1 ps_D_h_L1_1
sum nps_D2_h_1 ps_D_h_L2_1
assert reldif(nps_D1_h_1,ps_D_h_L1_1)<10e-7
assert reldif(nps_D2_h_1,ps_D_h_L2_1)<10e-7

// shortstack - all match
sum nps_Y_mpg_ss_1 ps_Y_mpg_ss_1
assert reldif(nps_Y_mpg_ss_1,ps_Y_mpg_ss_1)<10e-7

sum nps_D_price_ss_1 ps_D_price_ss_1
sum ps_D_price_h_ss_1 nps_D_price_h_ss_1
assert reldif(nps_D_price_ss_1,ps_D_price_ss_1)<10e-7
assert reldif(nps_D_price_h_ss_1,ps_D_price_h_ss_1)<10e-7


**************** H learner is regress, LIE enforced *******************

cap drop ps*
cap drop nps*

// Non-pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(nps) prefix
ddml E[Y|X], mname(nps)								: reg $Y $X
ddml E[D|Z,X], learner(D1) mname(nps)				: pystacked $D $X $Z, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D1) vname($D) mname(nps)		: reg {D} $X
ddml E[D|Z,X], learner(D2) mname(nps)				: pystacked $D $X $Z, method($L2) $pyopt
ddml E[D|X],   learner(D2) vname($D) mname(nps)		: reg {D} $X
// force branch to non-pystacked routine
ddml crossfit, mname(nps) crossfitother shortstack
ddml estimate, mname(nps)
mat b_nps=e(b)
mat V_nps=e(V)

// pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(ps) prefix
ddml E[Y|X], mname(ps)								: reg $Y $X
ddml E[D|Z,X], learner(D) mname(ps)					: pystacked $D $X $Z, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D) vname($D) mname(ps)		: reg {D} $X
ddml crossfit, mname(ps) shortstack 
ddml estimate, mname(ps)
mat b_ps=e(b)
mat V_ps=e(V)

// estimation results
assert mreldif(b_nps,b_ps)<10e-7
assert mreldif(V_nps,V_ps)<10e-7

// ss weights
ddml extract, show(ssweights) mname(nps)
mat ssw_nps=r(nps_D_price_ss)
ddml extract, show(ssweights) mname(ps)
mat ssw_ps=r(ps_D_price_ss)
// short-stacking h is redundant because the two h learners to be stacked
// by the non-pystacked code are identical:
assert nps_D1_h_ss_1==nps_D2_h_ss_1
// and the nnls solution with the ssw_h weights is arbitrary.
// the pystacked code anticipates this and skips the nnls.
// hence we check only first two (non-h) rows.
mat ssw_nps = ssw_nps[1..2,.]
mat ssw_ps = ssw_ps[1..2,.]
assert mreldif(ssw_nps,ssw_ps)<10e-7

// Y
sum nps_Y1_reg_1
sum ps_Y1_reg_1
assert reldif(nps_Y1_reg_1,ps_Y1_reg_1)<10e-7

// D learners for E[D|X,Z] - predicted values match
sum nps_D1_1 ps_D_L1_1
sum nps_D2_1 ps_D_L2_1
assert reldif(nps_D1_1,ps_D_L1_1)<10e-7
assert reldif(nps_D2_1,ps_D_L2_1)<10e-7

// h learner is regress for all D learners, so predicted values match w/LIE
sum nps_D1_h_1 ps_D_h_L1_1
sum nps_D2_h_1 ps_D_h_L2_1
assert reldif(nps_D1_h_1,ps_D_h_L1_1)<10e-7
assert reldif(nps_D2_h_1,ps_D_h_L2_1)<10e-7

// shortstack - all match
sum nps_Y_mpg_ss_1 ps_Y_mpg_ss_1
assert reldif(nps_Y_mpg_ss_1,ps_Y_mpg_ss_1)<10e-7

sum nps_D_price_ss_1 ps_D_price_ss_1
sum ps_D_price_h_ss_1 nps_D_price_h_ss_1
assert reldif(nps_D_price_ss_1,ps_D_price_ss_1)<10e-7
assert reldif(nps_D_price_h_ss_1,ps_D_price_h_ss_1)<10e-7


**************** H learner is pystacked, no LIE *******************

cap drop ps*
cap drop nps*

// Non-pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(nps) prefix
ddml E[Y|X], mname(nps)								: reg $Y $X
ddml E[D|Z,X], learner(D1) mname(nps)				: pystacked $D $X $Z, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D1) vname($D) mname(nps)		: pystacked {D} $X, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|Z,X], learner(D2) mname(nps)				: pystacked $D $X $Z, method($L2) $pyopt
ddml E[D|X],   learner(D2) vname($D) mname(nps)		: pystacked {D} $X, method($L2) $pyopt
// force branch to non-pystacked routine
ddml crossfit, mname(nps) crossfitother shortstack noenforce
ddml estimate, mname(nps)
mat b_nps=e(b)
mat V_nps=e(V)

// pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(ps) prefix
ddml E[Y|X], mname(ps)								: reg $Y $X
ddml E[D|Z,X], learner(D) mname(ps)					: pystacked $D $X $Z, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D) vname($D) mname(ps)		: pystacked {D} $X, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml crossfit, mname(ps) shortstack noenforce
ddml estimate, mname(ps)
mat b_ps=e(b)
mat V_ps=e(V)

// estimation results
assert mreldif(b_nps,b_ps)<10e-7
assert mreldif(V_nps,V_ps)<10e-7

// ss weights
ddml extract, show(ssweights) mname(nps)
mat ssw_nps=r(nps_D_price_ss)
ddml extract, show(ssweights) mname(ps)
mat ssw_ps=r(ps_D_price_ss)
assert mreldif(ssw_nps,ssw_ps)<10e-7

// Y
sum nps_Y1_reg_1
sum ps_Y1_reg_1
assert reldif(nps_Y1_reg_1,ps_Y1_reg_1)<10e-7

// all D learners w/no LIE - all match
sum nps_D1_1 ps_D_L1_1
sum nps_D2_1 ps_D_L2_1
assert reldif(nps_D1_1,ps_D_L1_1)<10e-7
assert reldif(nps_D2_1,ps_D_L2_1)<10e-7

sum nps_D1_h_1 ps_D_h_L1_1
sum nps_D2_h_1 ps_D_h_L2_1
assert reldif(nps_D1_h_1,ps_D_h_L1_1)<10e-7
assert reldif(nps_D2_h_1,ps_D_h_L2_1)<10e-7

// shortstack w/no LIE - all match
sum nps_Y_mpg_ss_1 ps_Y_mpg_ss_1
assert reldif(nps_Y_mpg_ss_1,ps_Y_mpg_ss_1)<10e-7

sum nps_D_price_ss_1 ps_D_price_ss_1
sum ps_D_price_h_ss_1 nps_D_price_h_ss_1
assert reldif(nps_D_price_ss_1,ps_D_price_ss_1)<10e-7
assert reldif(nps_D_price_h_ss_1,ps_D_price_h_ss_1)<10e-7


**************** H learner is pystacked, no LIE, nostd stack *******************

cap drop ps*
cap drop nps*

// Non-pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(nps) prefix
ddml E[Y|X], mname(nps)								: reg $Y $X
ddml E[D|Z,X], learner(D1) mname(nps)				: pystacked $D $X $Z, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D1) vname($D) mname(nps)		: pystacked {D} $X, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|Z,X], learner(D2) mname(nps)				: pystacked $D $X $Z, method($L2) $pyopt
ddml E[D|X],   learner(D2) vname($D) mname(nps)		: pystacked {D} $X, method($L2) $pyopt
// force branch to non-pystacked routine
ddml crossfit, mname(nps) crossfitother shortstack noenforce
ddml estimate, mname(nps)
mat b_nps=e(b)
mat V_nps=e(V)

// pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(ps) prefix
ddml E[Y|X], mname(ps)								: reg $Y $X
ddml E[D|Z,X], learner(D) mname(ps)					: pystacked $D $X $Z, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D) vname($D) mname(ps)		: pystacked {D} $X, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml crossfit, mname(ps) shortstack noenforce nostd
ddml estimate, mname(ps)
mat b_ps=e(b)
mat V_ps=e(V)

// estimation results
assert mreldif(b_nps,b_ps)<10e-7
assert mreldif(V_nps,V_ps)<10e-7

// ss weights
ddml extract, show(ssweights) mname(nps)
mat ssw_nps=r(nps_D_price_ss)
ddml extract, show(ssweights) mname(ps)
mat ssw_ps=r(ps_D_price_ss)
assert mreldif(ssw_nps,ssw_ps)<10e-7

// Y
sum nps_Y1_reg_1
sum ps_Y1_reg_1
assert reldif(nps_Y1_reg_1,ps_Y1_reg_1)<10e-7

// all D learners w/no LIE
sum nps_D1_1 ps_D_L1_1
sum nps_D2_1 ps_D_L2_1
assert reldif(nps_D1_1,ps_D_L1_1)<10e-7
assert reldif(nps_D2_1,ps_D_L2_1)<10e-7

sum nps_D1_h_1 ps_D_h_L1_1
sum nps_D2_h_1 ps_D_h_L2_1
assert reldif(nps_D1_h_1,ps_D_h_L1_1)<10e-7
assert reldif(nps_D2_h_1,ps_D_h_L2_1)<10e-7

// shortstack w/no LIE
sum nps_Y_mpg_ss_1 ps_Y_mpg_ss_1
assert reldif(nps_Y_mpg_ss_1,ps_Y_mpg_ss_1)<10e-7

sum nps_D_price_ss_1 ps_D_price_ss_1
sum ps_D_price_h_ss_1 nps_D_price_h_ss_1
assert reldif(nps_D_price_ss_1,ps_D_price_ss_1)<10e-7
assert reldif(nps_D_price_h_ss_1,ps_D_price_h_ss_1)<10e-7


********* H learner is pystacked, nostd stack, LIE - ss won't match ************

cap drop ps*
cap drop nps*

// Non-pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(nps) prefix
ddml E[Y|X], mname(nps)								: reg $Y $X
ddml E[D|Z,X], learner(D1) mname(nps)				: pystacked $D $X $Z, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D1) vname($D) mname(nps)		: pystacked {D} $X, method($L1) cmdopt1($opt1) $pyopt
ddml E[D|Z,X], learner(D2) mname(nps)				: pystacked $D $X $Z, method($L2) $pyopt
ddml E[D|X],   learner(D2) vname($D) mname(nps)		: pystacked {D} $X, method($L2) $pyopt
// force branch to non-pystacked routine
ddml crossfit, mname(nps) crossfitother shortstack
ddml estimate, mname(nps)
mat b_nps=e(b)
mat V_nps=e(V)

// pystacked code
cap prog drop crossfit
cap prog drop _ddml_estimate_linear
set seed 42
ddml init fiv, kfolds($kfolds) reps($reps) mname(ps) prefix
ddml E[Y|X], mname(ps)								: reg $Y $X
ddml E[D|Z,X], learner(D) mname(ps)					: pystacked $D $X $Z, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml E[D|X],   learner(D) vname($D) mname(ps)		: pystacked {D} $X, method($L1 $L2) cmdopt1($opt1) $pyopt
ddml crossfit, mname(ps) shortstack nostd
ddml estimate, mname(ps)
mat b_ps=e(b)
mat V_ps=e(V)

// estimation results - won't match
cap noi assert mreldif(b_nps,b_ps)<10e-7
assert _rc==9
cap noi assert mreldif(V_nps,V_ps)<10e-7
assert _rc==9

// ss weights
ddml extract, show(ssweights) mname(nps)
mat ssw_nps=r(nps_D_price_ss)
ddml extract, show(ssweights) mname(ps)
mat ssw_ps=r(ps_D_price_ss)
cap noi assert mreldif(ssw_nps,ssw_ps)<10e-7
assert _rc==9

// Y
sum nps_Y1_reg_1
sum ps_Y1_reg_1
assert reldif(nps_Y1_reg_1,ps_Y1_reg_1)<10e-7

// D learners for E[D|X,Z] - predicted values match
sum nps_D1_1 ps_D_L1_1
sum nps_D2_1 ps_D_L2_1
assert reldif(nps_D1_1,ps_D_L1_1)<10e-7
assert reldif(nps_D2_1,ps_D_L2_1)<10e-7

// h predicted values do not match
sum nps_D1_h_1 ps_D_h_L1_1
sum nps_D2_h_1 ps_D_h_L2_1
cap noi assert reldif(nps_D1_h_1,ps_D_h_L1_1)<10e-7
assert _rc==9
cap noi assert reldif(nps_D2_h_1,ps_D_h_L2_1)<10e-7
assert _rc==9

// Y shortstack
sum nps_Y_mpg_ss_1 ps_Y_mpg_ss_1
assert reldif(nps_Y_mpg_ss_1,ps_Y_mpg_ss_1)<10e-7

// D shortstacked learners for E[D|X,Z] - predicted values match
sum nps_D_price_ss_1 ps_D_price_ss_1
assert reldif(nps_D_price_ss_1,ps_D_price_ss_1)<10e-7
// h shortstack - won't match
sum nps_D_price_h_ss_1 ps_D_price_h_ss_1
cap noi assert reldif(nps_D_price_h_ss_1,ps_D_price_h_ss_1)<10e-7
assert _rc==9


cap log close
