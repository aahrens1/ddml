program define pylasso2, eclass
version 16.0

	if indexnot("=","`0'") {
		// single y variable
		gettoken yvars 0: 0
		tokenize "`0'", parse(",")
		local xvars "`1'"
		mac shift 1
		local 0 `"`*'"'
	}
	else {
		// multiple y variables
		tokenize "`0'", parse("=")
		local yvars "`1'"
		// skip the = sign
		tokenize "`3'", parse(",")
		local xvars "`1'"
		mac shift 1
		local 0 `"`*'"'
	}

	syntax [if] [in] [aweight fweight],				 ///
				[									 ///
					/// lambda(real 1)				 /// penalty
					lambda(string)					 /// penalties, one per y - can be numlist or matrix
					///alpha(real 1)				 /// elastic net parameter, alpha=1 is lasso, alpha=0 is ridge
					alpha(string)					 /// elastic net parameters, one per y - alpha=1 is lasso, alpha=0 is ridge
					cvalpha(string)					 /// range of elastic net alphas over which to cross-validate
					lars							 ///
					ic(string)						 /// bic, aic or empty
					n_jobs(integer -1)		 	 	 /// number of processors to use when computing stuff - default is all
					random_state(integer 0) 	 	 /// seed used by random number generator
					verbose(integer 0)	 			 	 	 /// controls verbosity
					warm_start(string asis)	 	 	 /// when set to true, reuse solution of previous call to fit
					prediction(string asis) 	     /// variable name to save predictions
					standardize                  	 /// standardize features
					NORMalize						 /// use scikit's normalize (implies no std-by-hand)
					UNITLoadings					 /// don't prestandardize
					STDCoef							 /// report standardized coefficients
					TOLOpt(real 1e-10)				 /// tol for optimization (ElasticNet default is 1e-4)
					NOCONStant						 /// standard option
					cv								 /// cross-validate
				]
	
	local varlist `yvars' `xvars'

	timer on 40

	**** n_jobs from pyforest ****
	* n_jobs: number of processors to use
	if "`n_jobs'"=="" local n_jobs -1
	if `n_jobs'<1 & `n_jobs'!=-1 {
		di as error "Syntax error: num_jobs() must be positive integer or -1."
		di as error " num_jobs() specifies number of processors to use; the default -1 means all."
		di as error " If not -1, this has to be a positive integer. But you should probably not mess around with this."
		exit 1
	}

	if "`warm_start'"=="" local warm_start False

	local normflag		= ("`normalize'"~="")
	local stdcoefflag	= ("`stdcoef'"~="")			// return coef estimates in std units
	local stdflag		= ("`unitloadings'"=="") | `stdcoefflag'
	if `normflag' {
		local stdflag	= 0							// sklearn's normalize means no std-by-hand
	}
	local consmodel		= ("`noconstant'"=="")
	local icflag		= ("`ic'"~="")
	local larsflag		= ("`lars'"~="") | `icflag'	// sklearn has IC selection only for LARS
	local cvflag		= ("`cv'"~="")
	if `stdcoefflag' | `consmodel'==0 {
		local consflag	= 0
	}
	else {
		local consflag = 1
	}

	marksample touse
	markout `touse' `yvars' `xvars'
	qui count if `touse'
	local N		= r(N)

	* Deal with factor and time-series vars

	fvexpand `yvars' if `touse'
	local yvars `r(varlist)'
	// now create a varlist with temps etc. that can be passed to Python
	fvrevar `yvars' if `touse'
	local yvars_t `r(varlist)'
	
	fvexpand `xvars' if `touse'
	local xvars `r(varlist)'
	// now create a varlist with temps etc. that can be passed to Python
	fvrevar `xvars' if `touse'
	local xvars_t `r(varlist)'

	local nyvars	: word count `yvars'
	local nxvars	: word count `xvars'
	
	// process lambda - can be a matrix or a numlist
	tempname lmat
	cap mat list `lambda'
	if "`lambda'"=="" {
		mat `lmat' = J(1,`nyvars',1)
	}
	else if _rc==0 {
		mat `lmat' = `lambda'
		}
	else {
		mata: st_matrix("`lmat'",strtoreal(tokens("`lambda'")))
		if colsof(`lmat')==1 & `nyvars'>1 {
			mat `lmat' = J(1,`nyvars',el(`lmat',1,1))
		}
	}

	// process alpha - can be a matrix or a numlist
	tempname amat
	cap mat list `alpha'
	if "`alpha'"=="" {
		// default is lasso
		mat `amat' = J(1,`nyvars',1)
	}
	else if _rc==0 {
		mat `amat' = `alpha'
		}
	else {
		mata: st_matrix("`amat'",strtoreal(tokens("`alpha'")))
		if colsof(`amat')==1 & `nyvars'>1 {
			mat `amat' = J(1,`nyvars',el(`amat',1,1))
		}
	}

	// ridgeflag=1 if all eqns are ridge
	mata: st_local("ridgeflag",strofreal(!sum(st_matrix("`amat'"):~=0)))
	
	// process cvalpha - can be a matrix or a numlist
	tempname cvamat
	cap mat list `cvalpha'
	if "`cvalpha'"=="" {
		// default is lasso
		mat `cvamat' = 1
	}
	else if _rc==0 {
		mat `cvamat' = `cvalpha'
		}
	else {
		mata: st_matrix("`cvamat'",strtoreal(tokens("`cvalpha'")))
	}
	
	timer off 40
	timer on 50
	
	python: run_elastic_net(					///
		"`amat'",								///
		`ridgeflag',							///
		"`lmat'",								///
		"`cvamat'",								///
		`larsflag',								///
		`icflag',								///
		"`ic'",									///
		"`yvars_t'",							///
		"`xvars_t'",							///
		"`touse'",								///
		`n_jobs',								///
		`random_state',							///
		`verbose',								///
		`warm_start',							///
		 `stdflag',								///
		 `stdcoefflag',							///
		 `consflag',							///
		 `tolopt',								///
		 `normflag',							///
		 `cvflag')

	timer off 50
	timer on 60

	tempname b lambdamat sk_lambdamat alphamat sk_alphamat r2mat
	// vectors are returned from Python as col vectors; transpose
	mat `b' = r(beta)'
	mat `lambdamat'		= r(lambda)'
	mat `sk_lambdamat'	= r(sk_lambda)'
	mat `alphamat'		= r(alpha)'
	mat `sk_alphamat'	= r(sk_alpha)'
	mat `r2mat'			= r(r2)'

	// colnames for b
	if `consflag' {
		local cnames `xvars' _cons
	}
	else {
		local cnames `xvars'
	}
	if `nyvars'==1 {
		mat colnames `b' = `cnames'
	}
	else {
		foreach eqname in `yvars' {
			foreach cn in `cnames' {
				local fullcolnames `fullcolnames' `eqname':`cn'
			}
		}
		mat colnames `b' = `fullcolnames'
	}
	
	// colnames for lambda and alpha
	mat colnames `lambdamat'	= `yvars'
	mat colnames `sk_lambdamat'	= `yvars'
	mat colnames `alphamat'		= `yvars'
	mat colnames `sk_alphamat'	= `yvars'
	mat colnames `r2mat'		= `yvars'

	*** fix colnames of beta vectors to include omitted "o." notation
	// trick is to use _ms_findomitted utility but give it
	// diag(beta) as vcv matrix where it looks for zeros
	// also build in fv info
	tempname tempvmat
	mat `tempvmat'	= diag(`b')
	_ms_findomitted	`b' `tempvmat'
	_ms_build_info	`b' if `touse'

	ereturn post `b', depname(`depname') esample(`touse') obs(`N') properties(b)

	// ereturn scalar mae				= `mae'
	// ereturn scalar rmse				= `rmse'

	ereturn local cmd				pylasso2
	ereturn local predict			pylasso2_p
	ereturn local depvar			`yvars'
	
	ereturn scalar cons				= `consflag'
	ereturn scalar std				= `stdflag'
	ereturn scalar stdcoef			= `stdcoefflag'
	ereturn scalar nyvars			= `nyvars'

	ereturn matrix lambda			= `lambdamat', copy
	ereturn matrix sk_lambda		= `sk_lambdamat', copy
	ereturn matrix alpha			= `alphamat', copy
	ereturn matrix sk_alpha			= `sk_alphamat', copy
	ereturn matrix r2				= `r2mat', copy

	di
	ereturn di, noomitted vsquish nofvlabel
	
	tokenize `yvars'
	di _col(13) "lambda" _col(22) "lambda (sk)" _col(35) "alpha" _col(42) "alpha (sk)" _col(55) "R-sq"
	forvalues i=1/`nyvars' {
		local lvalue		= el(`lambdamat',1,`i')
		local sk_lvalue		= el(`sk_lambdamat',1,`i')
		local avalue		= el(`alphamat',1,`i')
		local sk_avalue		= el(`sk_alphamat',1,`i')
		local r2value		= el(`r2mat',1,`i')
		di "``i''" _c
		di _col(10) %10.4f `lvalue' _c
		di _col(20) %10.4f `sk_lvalue' _c
		di _col(30) %10.4f `avalue' _c
		di _col(40) %10.4f `sk_avalue' _c
		di _col(50) %10.4f `r2value' _c
		di
	}

	timer off 60

end


*===============================================================================
* Python helper function
*===============================================================================

version 16.0
python:

#-------------------------------------------------------------------------------
# Import required packages and attempt to install w/ Pip if that fails
#-------------------------------------------------------------------------------

# Import required Python modules (pandas, scikit-learn, sfi)
from pandas import DataFrame
from sfi import Data,Matrix,Scalar,SFIToolkit
from sklearn import metrics, preprocessing
import numpy as np
from sklearn.linear_model import ElasticNet,Ridge,RidgeCV,LassoLars,ElasticNetCV,LassoLarsCV,LassoLarsIC
from sklearn.pipeline import make_pipeline

# To pass objects to Stata
import __main__

# Set random seed
import random
random.seed(50)

#-------------------------------------------------------------------------------
# Define Python function: run_elastic_net
#-------------------------------------------------------------------------------

# MS: "lambda" probably a reserved word

def run_elastic_net(amat,ridgeflag,lmat,cvamat,larsflag,icflag,ic,yvars,xvars,touse,n_jobs,random_state,verbose,warm_start,stdflag,stdcoefflag,consflag,tolopt,normflag,cvflag):

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 51")
	
	# Load data
	x = Data.get(xvars,selectvar=touse)
	# y = Data.get(yvars,selectvar=touse)
	y = np.array(Data.get(yvars,selectvar=touse))
	if y.ndim==1:
		y = np.array([y]).T
	SFIToolkit.stata("timer off 51")

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 52")
	
	# If standardizing, scale features to mean zero, std dev one
	if stdflag==1:
		SFIToolkit.stata("timer on 91")
		x_mean = np.mean(x,axis=0)
		SFIToolkit.stata("timer off 91")
		SFIToolkit.stata("timer on 92")
		x_std = np.std(x,axis=0)
		SFIToolkit.stata("timer off 92")
		# If SD is zero, replace with 1
		x_std[x_std==0] = 1
		SFIToolkit.stata("timer on 93")
		x = (x-x_mean)/x_std
		SFIToolkit.stata("timer on 93")

	SFIToolkit.stata("timer off 52")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 53")

	# Misc
	n = len(y)
	y_mean = np.mean(y,axis=0)
	y_std = np.std(y,axis=0)

	# Always standardize y
	y = (y-y_mean)/y_std
	
	# lambda supplied as a Stata matrix; convert to Python vector
	pylambda = Matrix.get(lmat)
	pylambda = pylambda[0]
	
	# alpha supplied as a Stata matrix; convert to a Python vector
	pyalpha = Matrix.get(amat)
	pyalpha = pyalpha[0]
	pycvalpha = Matrix.get(cvamat)
	pycvalpha = pycvalpha[0]
	
	# Alway need to rescale lambda penalty
	lpenalty = pylambda/y_std
	# Initialize
	l1_ratio = pyalpha
	cvl1_ratio = pycvalpha
	# More rescaling
	if ridgeflag==1 and normflag==0:
		lpenalty = lpenalty*n
	elif ridgeflag==0 and normflag==1 and cvflag==0:
		a = lpenalty*l1_ratio*1/np.sqrt(n)
		b = lpenalty*(1.0 - np.array(l1_ratio))*1/n
		lpenalty = a+b
		l1_ratio = a/(a+b)
	# Disabled - elastic net + normalize + CV over alpha (and code doesn't work anyway)
	#elif ridgeflag==0 and normflag==1 and cvflag==1:
		#a = lpenalty*cvl1_ratio*1/np.sqrt(n)
		#b = lpenalty*(1.0 - np.array(cvl1_ratio))*1/n
		#lpenalty = a+b
		#cvl1_ratio = a/(a+b)

	SFIToolkit.stata("timer off 53")

	if ridgeflag==1:
	
		# All ridge
		# Ridge handled separately because it can estimate multi-output with separate lambda penalties
		
		if cvflag==1:
			model = RidgeCV(fit_intercept=consflag, normalize=normflag, alpha_per_target=True)
		else:		
			# Ridge uses a different definition of the penalty
			model = Ridge(alpha=lpenalty, fit_intercept=consflag, normalize=normflag, tol=tolopt)
		
		SFIToolkit.stata("timer on 54")
		# Estimate model
		model.fit(x, y)
		SFIToolkit.stata("timer off 54")
		
		for col in range(y.shape[1]):

			# Start Stata timer
			SFIToolkit.stata("timer on 56")
		
			# Add constant, unstandardize, etc.
			if stdflag==0:
				if consflag==1:
					b = np.append(model.coef_[col]*y_std.T[col],(model.intercept_[col]*y_std.T[col] + y_mean.T[col]))
				else:
					b = model.coef_[col]
			else:
				# estimation used standardization
				if stdcoefflag==0:
					# Unstandardize coefficients
					b = model.coef_[col] / x_std * y_std.T[col]
					if consflag==1:
						# Get constant; use prestandardized means
						predict_nostd = np.sum(x_mean * b)
						b = np.append(b,y_mean.T[col] - predict_nostd)
				else:
					b = model.coef_[col]
	
			if col==0:
				ball = b
			else:
				ball = np.append(ball,b)

		if cvflag==1:
			cvlambda = model.alpha_
			cvalpha = [0]*y.shape[1]
			# rescale to >0
			r2 = -model.best_score_
		else:
			# Returns single R-sq = 'uniform_average' so not suitable.
			# r2 = model.score(x,y)
			r2 = metrics.r2_score(y,model.predict(x),multioutput='raw_values')
		
	else:
	
		# All estimators that don't accept multi-output with separate tuning parameters
	
		for col in range(y.shape[1]):
			# equivalent
			# print(y[:,col])
			# print(y.T[col])
	
			##############################################################
	
			# Initialize model object
			if icflag==1:
				model = LassoLarsIC(criterion=ic,fit_intercept=consflag, normalize=normflag)
			elif cvflag==1 and larsflag==1:
				model = LassoLarsCV(fit_intercept=consflag, normalize=normflag)
			elif cvflag==1:
				model = ElasticNetCV(l1_ratio=cvl1_ratio, random_state=random_state, fit_intercept=consflag, normalize=normflag, tol=tolopt)
			elif larsflag==1:
				model = LassoLars(alpha=lpenalty[col], random_state=random_state, fit_intercept=consflag, normalize=normflag)
			elif l1_ratio[col]>0:
				model = ElasticNet(alpha=lpenalty[col], l1_ratio=l1_ratio[col], random_state=random_state, fit_intercept=consflag, normalize=normflag, tol=tolopt)
			else:
				# Ridge uses a different definition of the penalty
				model = Ridge(alpha=lpenalty[col]*n, fit_intercept=consflag, normalize=normflag, tol=tolopt)

			##############################################################
			
			# Start Stata timer
			SFIToolkit.stata("timer on 54")
			
			# Estimate model
			model.fit(x, y.T[col])

			SFIToolkit.stata("timer off 54")

			##############################################################
		
			# Start Stata timer
			SFIToolkit.stata("timer on 56")
		
			# Add constant, unstandardize, etc.
			if stdflag==0:
				if consflag==1:
					b = np.append(model.coef_*y_std.T[col],(model.intercept_*y_std.T[col] + y_mean.T[col]))
				else:
					b = model.coef_
			else:
				# estimation used standardization
				if stdcoefflag==0:
					# Unstandardize coefficients
					b = model.coef_ / x_std * y_std.T[col]
					if consflag==1:
						# Get constant; use prestandardized means
						predict_nostd = np.sum(x_mean * b)
						b = np.append(b,y_mean.T[col] - predict_nostd)
				else:
					b = model.coef_
	
			if col==0:
				ball = b
				if cvflag==1 or icflag==1:
					cvlambda = model.alpha_
					if len(cvl1_ratio)>0 and larsflag==0:
						cvalpha = model.l1_ratio_
					else:
						cvalpha = cvl1_ratio
			else:
				ball = np.append(ball,b)
				if cvflag==1 or icflag==1:
					cvlambda = np.append(cvlambda,model.alpha_)
					if len(cvl1_ratio)>0 and larsflag==0:
						cvalpha = np.append(cvalpha,model.l1_ratio_)
					else:
						cvalpha = np.append(cvalpha,cvl1_ratio)

			if col==0:
				r2 = model.score(x, y.T[col])
			else:
				r2 = np.append(r2,model.score(x, y.T[col]))
				
			SFIToolkit.stata("timer off 56")
	
	##############################################################
	
	# Start Stata timer
	# SFIToolkit.stata("timer on 55")
	# Get in-sample prediction
	# model_predict = model.predict(x)
	# SFIToolkit.stata("timer off 55")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 57")

	#b = np.array([b])
	#ball = np.array([ball])
	Matrix.store("r(beta)",ball)
	
	if cvflag==1 or icflag==1:
		if ridgeflag==1 and normflag==0:
			pylambda = cvlambda*y_std/n
			pyalpha = cvalpha
		elif ridgeflag==1:
			pylambda = cvlambda*y_std
			pyalpha = cvalpha
		elif normflag==1:
			a = np.array(cvlambda)*np.array(cvalpha)
			b = cvlambda - a
			a = a*np.sqrt(n)
			b = b*n
			pylambda = (a+b)
			pylambda = pylambda*y_std
			pyalpha = a/(a+b)
		else:
			pylambda = cvlambda*y_std
			pyalpha = cvalpha
		Matrix.store("r(lambda)", pylambda)
		Matrix.store("r(sk_lambda)", cvlambda)
		Matrix.store("r(alpha)",pyalpha)
		Matrix.store("r(sk_alpha)",cvalpha)
	else:
		Matrix.store("r(lambda)", pylambda)
		Matrix.store("r(sk_lambda)", lpenalty)
		Matrix.store("r(alpha)", pyalpha)
		Matrix.store("r(sk_alpha)", l1_ratio)

	Matrix.store("r(r2)",r2)

	# Get mae, rmse
	# SFIToolkit.stata("timer on 96")
	# mae = metrics.mean_absolute_error(y.T[col], model_predict)
	# mse = metrics.mean_squared_error(y.T[col], model_predict)
	# SFIToolkit.stata("timer off 96")
	# rmse = np.sqrt(mse)
	# if stdcoefflag==0:
	# 	mae = mae * y_std.T[col]
	# 	rmse = rmse * y_std.T[col]
	# Scalar.setValue("r(mae)", mae, vtype='visible')
	# Scalar.setValue("r(rmse)", rmse, vtype='visible')

	# Unstandardize predictions if required
	# if stdflag==1 and prediction != "":
	#	SFIToolkit.stata("timer on 97")
	#	model_predict = model_predict * y_std.T[col] + y_mean.T[col]
	#	SFIToolkit.stata("timer off 97")

	SFIToolkit.stata("timer off 57")

	##############################################################

	if verbose==1:
		print(model)

	##############################################################

	# Pass objects back to __main__ namespace to interact w/ later
	__main__.model_object  = model
	__main__.model_xvars = xvars
	# If number of y vars > 1, store for latest
	__main__.model_y_std = y_std[col]
	__main__.model_y_mean = y_mean[col]
	if stdflag==1:
		__main__.model_x_std = x_std
		__main__.model_x_mean = x_mean

	
end
