sysuse auto, clear
cap prog drop pylasso2

// lasso - match
lasso2 mpg rep78-foreign, lglmnet lambda(1)
pylasso2 mpg rep78-foreign, verbose(1) lambda(1)
pylasso2 mpg rep78-foreign, verbose(1) lambda(1) normalize

// lasso lars - match
lasso2 mpg rep78-foreign, lglmnet lambda(1)
pylasso2 mpg rep78-foreign, verbose(1) lambda(1) lars
pylasso2 mpg rep78-foreign, verbose(1) lambda(1) lars normalize

// elastic net
lasso2 mpg rep78-foreign, lglmnet lambda(1) alpha(0.4)
pylasso2 mpg rep78-foreign, verbose(1) lambda(1) alpha(0.4)
pylasso2 mpg rep78-foreign, verbose(1) lambda(1) alpha(0.4) normalize

// ridge - match
lasso2 mpg rep78-foreign, lglmnet lambda(1) alpha(0)
pylasso2 mpg rep78-foreign, verbose(1) lambda(1) alpha(0)
pylasso2 mpg rep78-foreign, verbose(1) lambda(1) alpha(0) normalize

// CV

// lasso - match
pylasso2 mpg rep78-foreign, verbose(1) cv
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')
pylasso2 mpg rep78-foreign, verbose(1) cv normalize
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')

// elastic net
pylasso2 mpg rep78-foreign, verbose(1) cv cvalpha(0.4)
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(0.4)
pylasso2 mpg rep78-foreign, verbose(1) cv normalize cvalpha(0.4)
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(0.4)

// elastic net - cv over alpha and lambda
pylasso2 mpg rep78-foreign, verbose(1) cv cvalpha(0.1 0.5 0.9)
local lambda = el(e(lambda),1,1)
local alpha = el(e(alpha),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(`alpha')
pylasso2 mpg rep78-foreign, verbose(1) cv normalize cvalpha(0.1 0.5 0.9)
local lambda = el(e(lambda),1,1)
local alpha = el(e(alpha),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(`alpha')

// lasso lars - match
pylasso2 mpg rep78-foreign, verbose(1) lars cv
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')
pylasso2 mpg rep78-foreign, verbose(1) lars cv normalize
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')

// ridge - match
pylasso2 mpg rep78-foreign, verbose(1) cv alpha(0)
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(0)
pylasso2 mpg rep78-foreign, verbose(1) cv normalize alpha(0)
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(0)

// IC

// lasso bic - match
pylasso2 mpg rep78-foreign, verbose(1) ic(bic)
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')
pylasso2 mpg rep78-foreign, verbose(1) ic(bic) normalize
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')

// lasso aic - match
pylasso2 mpg rep78-foreign, verbose(1) ic(aic)
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')
pylasso2 mpg rep78-foreign, verbose(1) ic(aic) normalize
local lambda = el(e(lambda),1,1)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda')


**********************

// match
pylasso2 price mpg = rep78-foreign, lambda(100 1) verbose(1)
pylasso2 price mpg = rep78-foreign, lambda(100 1) verbose(1) normalize

// elastic net, 2 dep vars
pylasso2 price mpg = rep78-foreign, verbose(1) cv cvalpha(0.01 0.1 0.3 0.5 0.7 0.9 0.99)
local lambda = el(e(lambda),1,1)
local alpha = el(e(alpha),1,1)
lasso2 price rep78-foreign, lglmnet lambda(`lambda') alpha(`alpha')
pylasso2 price mpg = rep78-foreign, verbose(1) cv cvalpha(0.01 0.1 0.3 0.5 0.7 0.9 0.99)
local lambda = el(e(lambda),1,2)
local alpha = el(e(alpha),1,2)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(`alpha')

// elastic net, 2 dep vars, normalize
pylasso2 price mpg = rep78-foreign, verbose(1) cv cvalpha(0.01 0.1 0.3 0.5 0.7 0.9 0.99) normalize
local lambda = el(e(lambda),1,1)
local alpha = el(e(alpha),1,1)
lasso2 price rep78-foreign, lglmnet lambda(`lambda') alpha(`alpha')
pylasso2 price mpg = rep78-foreign, verbose(1) cv cvalpha(0.01 0.1 0.3 0.5 0.7 0.9 0.99) normalize
local lambda = el(e(lambda),1,2)
local alpha = el(e(alpha),1,2)
lasso2 mpg rep78-foreign, lglmnet lambda(`lambda') alpha(`alpha')

