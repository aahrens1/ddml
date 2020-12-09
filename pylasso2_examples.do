clear all
cap prog drop pylasso2
cap cd "C:\Users\ecomes\Documents\GitHub\lassopack\lassopack_v141"
sysuse auto, clear
drop if rep78==.
gen double price000 = price/1000
gen double price0000 = price/10000

foreach var of varlist price mpg-foreign {
	qui sum `var'
	qui gen double sd_`var' = (`var'-r(mean))/(r(sd) * sqrt((r(N)-1)/r(N)))
	di "SD(`var') = " _col(20) %13.8f r(sd) * sqrt((r(N)-1)/r(N))
}

// SDs match
// pylasso2 sd_price mpg-foreign, lambda(.1) standardize
// Std vars have SD=1
// pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1) standardize

// match - prestandardized variables
pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1) nocons
mat list e(b)
mat b_pylasso2=e(b)
qui lasso2 sd_price sd_mpg-sd_foreign, lglmnet lambda(.1) nocons
mat list e(sbetaAll)
mat b_lasso2=e(sbetaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// match - standardized coefs, programs do the standardizing
pylasso2 price mpg-foreign, lambda(100) stdcoef
mat list e(b)
mat b_pylasso2=e(b)
qui lasso2 price mpg-foreign, lglmnet lambda(100) stdcoef
mat list e(sbetaAll)
mat b_lasso2=e(sbetaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// standardizing already standardized X changes nothing
pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1) unitloadings
pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1)

// normalizing already standardized X changes nothing
pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1) unitloadings
pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1) unitloadings normalize

// with a pre-standardized y, standardized is equiv to std+std coefs
pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1) unitloadings nocons
pylasso2 sd_price mpg-foreign, lambda(.1) stdcoef
mat b_pylasso2=e(b)
// match
qui lasso2 sd_price mpg-foreign, lglmnet lambda(.1) stdcoef
mat list e(sbetaAll)
mat b_lasso2=e(sbetaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// match - standardization and then unstandardized coefs
pylasso2 price mpg-foreign, lambda(1000)
mat list e(b)
mat b_pylasso2=e(b)
lasso2 price mpg-foreign, lglmnet lambda(1000)
mat list e(betaAll)
mat b_lasso2=e(betaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// match - no standardization at all
// need to fix constant
pylasso2 price mpg-foreign, lambda(1000) unitloadings
mat list e(b)
mat b_pylasso2=e(b)
lasso2 price mpg-foreign, lglmnet lambda(1000) unitloadings
mat list e(betaAll)
mat b_lasso2=e(betaAll)
// note lower tolerance
assert mreldif(b_pylasso2,b_lasso2) < 1e-6

// training
cap drop training
gen training=(_n<=50)
// in std units
pylasso2 price mpg-foreign if _n<=50, lambda(1000) stdcoef
ereturn list
pylasso2 price mpg-foreign, lambda(1000) training(training) stdcoef
ereturn list
// in natural units, no standardization
pylasso2 price mpg-foreign if _n<=50, lambda(1000) unitloadings
ereturn list
pylasso2 price mpg-foreign, lambda(1000) training(training) unitloadings
ereturn list
// in natural units with standardization loadings
pylasso2 price mpg-foreign if _n<=50, lambda(1000)
ereturn list
pylasso2 price mpg-foreign, lambda(1000) training(training)
ereturn list

**** elastic net ****

// only prestandardization of dep var - match
pylasso2 sd_price mpg-foreign, lambda(1) alpha(0.5) unitloadings
mat b_pylasso2=e(b)
lasso2 sd_price mpg-foreign, lglmnet lambda(1) alpha(0.5) unitloadings
mat b_lasso2=e(betaAll)
// note lower tolerance
assert mreldif(b_pylasso2,b_lasso2) < 1e-7

// prestandardization of dep var, programs standardize X - match
pylasso2 sd_price mpg-foreign, lambda(1) alpha(0.5)
mat b_pylasso2=e(b)
lasso2 sd_price mpg-foreign, lglmnet lambda(1) alpha(0.5)
mat b_lasso2=e(betaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// programs do the standardizing - match
pylasso2 price mpg-foreign, lambda(100) alpha(0.5)
mat list e(b)
mat b_pylasso2=e(b)
qui lasso2 price mpg-foreign, lglmnet lambda(100) alpha(0.5)
mat list e(betaAll)
mat b_lasso2=e(betaAll)
// note lower tolerance
assert mreldif(b_pylasso2,b_lasso2) < 1e-7

// no standardization - match
pylasso2 price mpg-foreign, lambda(100) alpha(0.5) unitloadings
mat b_pylasso2=e(b)
lasso2 price mpg-foreign, lglmnet lambda(100) alpha(0.5) unitloadings
mat b_lasso2=e(betaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// how it should behave - change scale of depvar and lambda
lasso2 price mpg-foreign, lglmnet lambda(100) alpha(0.5) unitloadings
lasso2 price000 mpg-foreign, lglmnet lambda(0.100) alpha(0.5) unitloadings

// how sklearn behaves
pylasso2 price mpg-foreign, lambda(100) alpha(0.5) unitloadings
pylasso2 price000 mpg-foreign, lambda(0.100) alpha(0.5) unitloadings
pylasso2 price0000 mpg-foreign, lambda(0.0100) alpha(0.5) unitloadings

*** ridge ***

// programs do the standardizing
pylasso2 price mpg-foreign, lambda(100) alpha(0) stdcoef
mat list e(b)
mat b_pylasso2=e(b)
qui lasso2 price mpg-foreign, lglmnet lambda(100) alpha(0) stdcoef
mat list e(sbetaAll)
mat b_lasso2=e(sbetaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// only standardization of dep var - match
pylasso2 sd_price mpg-foreign, lambda(1) alpha(0) unitloadings
mat b_pylasso2=e(b)
lasso2 sd_price mpg-foreign, lglmnet lambda(1) alpha(0) unitloadings
mat b_lasso2=e(betaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// no standardization - match
pylasso2 price mpg-foreign, lambda(100) alpha(0) unitloadings
mat b_pylasso2=e(b)
lasso2 price mpg-foreign, lglmnet lambda(100) alpha(0) unitloadings
mat b_lasso2=e(betaAll)
assert mreldif(b_pylasso2,b_lasso2) < 1e-8

// how it should behave - change scale of depvar and lambda
lasso2 price mpg-foreign, lglmnet lambda(100) alpha(0) unitloadings
lasso2 price000 mpg-foreign, lglmnet lambda(0.100) alpha(0) unitloadings

// how sklearn behaves
pylasso2 price mpg-foreign, lambda(100) alpha(0) unitloadings
pylasso2 price000 mpg-foreign, lambda(0.100) alpha(0) unitloadings
pylasso2 price0000 mpg-foreign, lambda(0.0100) alpha(0) unitloadings

************* speed check *************

timer clear
timer on 1
forvalues i=1/100 {
	qui pylasso2 price mpg-foreign, lambda(100) stdcoef
}
timer off 1
timer on 2
forvalues i=1/100 {
	qui lasso2 price mpg-foreign, lglmnet lambda(100) stdcoef
}
timer off 2
timer list

timer clear
timer on 1
forvalues i=1/100 {
	qui pylasso2 sd_price sd_mpg-sd_foreign, lambda(.1) unitloadings
}
timer off 1
timer on 2
forvalues i=1/100 {
	qui lasso2 sd_price sd_mpg-sd_foreign, lglmnet lambda(.1) unitloadings
}
timer off 2
timer list
