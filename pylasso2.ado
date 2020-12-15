program define pylasso2, eclass
version 16.0
syntax varlist(min=2 fv) [if] [in] [aweight fweight],	 ///
				[									 ///
					alpha(real 1)					 /// elastic net parameter, alpha=1 is lasso, alpha=0 is ridge
					lambda(real 1)					 /// penalty
					lars							 ///
					n_jobs(integer -1)		 	 	 /// number of processors to use when computing stuff - default is all
					random_state(integer -1) 	 	 /// seed used by random number generator
					verbose		 			 	 	 /// controls verbosity
					warm_start(string asis)	 	 	 /// when set to true, reuse solution of previous call to fit
					prediction(string asis) 	     /// variable name to save predictions
					standardize                  	 /// standardize features
					normalize						 /// use scikit's normalize
					UNITLoadings					 /// don't prestandardize
					STDCoef							 /// report standardized coefficients
					TOLOpt(real 1e-10)				 /// tol for optimization (ElasticNet default is 1e-4)
					NOCONStant						 /// standard option
				]

	timer on 40

	**** n_jobs and random_state code from pyforest ****
	* n_jobs: number of processors to use
	if "`n_jobs'"=="" local n_jobs -1
	if `n_jobs'<1 & `n_jobs'!=-1 {
		di as error "Syntax error: num_jobs() must be positive integer or -1."
		di as error " num_jobs() specifies number of processors to use; the default -1 means all."
		di as error " If not -1, this has to be a positive integer. But you should probably not mess around with this."
		exit 1
	}
	if "`random_state'"=="-1" local random_state None
	if "`random_state'"!="" & "`random_state'"!="None" {
		if `random_state'<1 {
			di as error "Syntax error: random_state should be a positive integer."
			exit 1
		}
		set seed `random_state'
	}

	if "`verbose'"=="" local verbose 0

	if "`warm_start'"=="" local warm_start False

	local stdcoefflag	= ("`stdcoef'"~="")		// return coef estimates in std units
	local stdflag		= ("`unitloadings'"=="") | `stdcoefflag'
	local normflag		= ("`normalize'"~="")
	local consmodel		= ("`noconstant'"=="")
	local larsflag		= ("`lars'"~="")
	if `stdcoefflag' | `consmodel'==0 {
		local consflag	= 0
	}
	else {
		local consflag = 1
	}
	
	marksample touse
	qui count if `touse'
	local N		= r(N)

	* Deal with factor and time-series vars
	// first expand and unabbreviate
	fvexpand `varlist' if `touse'
	local varlist `r(varlist)'
	// now create a varlist with temps etc. that can be passed to Python
	fvrevar `varlist' if `touse'
	local varlist_t `r(varlist)'

	* Pass varlist into varlists called yvar and xvars
	gettoken yvar xvars : varlist
	gettoken yvar_t xvars_t : varlist_t
	
	// create an empty prediction variable
	if "`prediction'"~="" {
		qui gen double `prediction' = .
	}

	timer off 40
	timer on 50
	
	python: run_elastic_net(					///
		`alpha',								///
		`lambda',								///
		`larsflag',								///
		"`yvar_t'",								///
		"`xvars_t'",							///
		"`touse'",								///
		`n_jobs',								///
		`random_state',							///
		`verbose',								///
		`warm_start',							///
		"`prediction'",							///
		 `stdflag',								///
		 `stdcoefflag',							///
		 `consflag',							///
		 `tolopt',								///
		 `normflag')

	timer off 50
	timer on 60

	tempname b
	mat `b' = r(beta)
	if `consflag' {
		local cnames `xvars' _cons
	}
	else {
		local cnames `xvars'
	}
	mat colnames `b' = `cnames'
	
	local mae			= r(mae)
	local rmse			= r(rmse)

	local sk_lambda		= r(lambda)
	local sk_alpha		= r(alpha)

	*** fix colnames of beta vectors to include omitted "o." notation
	// trick is to use _ms_findomitted utility but give it
	// diag(beta) as vcv matrix where it looks for zeros
	// also build in fv info
	tempname tempvmat
	mat `tempvmat'	= diag(`b')
	_ms_findomitted	`b' `tempvmat'
	_ms_build_info	`b' if `touse'
	// rowname gets zapped for some reason after posting
	mat rownames `b' = `yvar'

	ereturn post `b', depname(`yvar') esample(`touse') obs(`N') properties(b)

	ereturn scalar mae				= `mae'
	ereturn scalar rmse				= `rmse'

	ereturn local cmd				pylasso2
	ereturn local predict			pylasso2_p
	ereturn local depvar			`yvar'
	
	ereturn scalar cons				= `consflag'
	ereturn scalar std				= `stdflag'
	ereturn scalar stdcoef			= `stdcoefflag'

	di
	ereturn di, noomitted vsquish nofvlabel
	di "sklearn lambda=" _col(20) %10.4f `sk_lambda'
	di "sklearn alpha=" _col(20) %10.4f `sk_alpha'

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
from sklearn.ensemble import RandomForestClassifier,RandomForestRegressor
from sfi import Data,Matrix,Scalar,SFIToolkit
from sklearn import metrics, preprocessing
import numpy as np
# MS
from sklearn.linear_model import ElasticNet,Ridge,LassoLars
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

def run_elastic_net(lratio,lpenalty,larsflag,yvar,xvars,touse,n_jobs,random_state,verbose,warm_start,prediction,stdflag,stdcoefflag,consflag,tolopt,normflag):

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 51")
	
	# Load data
	y = Data.get(yvar,selectvar=touse)
	x = Data.get(xvars,selectvar=touse)
	SFIToolkit.stata("timer off 51")

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 52")
	
	# If standardizing, scale features to mean zero, std dev one
	if stdflag==1:
		x_mean = np.mean(x,axis=0)
		x_std = np.std(x,axis=0)
		# If SD is zero, replace with 1
		x_std[x_std==0] = 1
		x = (x-x_mean)/x_std

	SFIToolkit.stata("timer off 52")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 53")

	# Misc
	n = len(y)
	y_mean = np.mean(y)
	y_std = np.std(y)
	
	# Always standardize y
	y = (y-y_mean)/y_std
	
	# Need to rescale lambda penalty
	lpenalty = lpenalty/y_std
	
	SFIToolkit.stata("timer off 53")

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 54")

	# Initialize model object
	if larsflag==1:
		model = LassoLars(alpha=lpenalty, random_state=0, fit_intercept=consflag, normalize=normflag)
	elif lratio>0:
		model = ElasticNet(alpha=lpenalty, l1_ratio=lratio, random_state=0, fit_intercept=consflag, normalize=normflag, tol=tolopt)
	else:
		# Ridge uses a different definition of the penalty
		model = Ridge(alpha=lpenalty*n, fit_intercept=consflag, normalize=normflag, tol=tolopt)

	# Estimate model
	model.fit(x, y)
	
	SFIToolkit.stata("timer off 54")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 55")
	
	# Get in-sample prediction
	model_predict = model.predict(x)

	SFIToolkit.stata("timer off 55")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 56")
	
	# Stata coefficient vectors are row vectors 
	if stdflag==0:
		if consflag==1:
			b = np.append(model.coef_*y_std,(model.intercept_*y_std + y_mean))
		else:
			b = model.coef_
		b = np.array([b])
	else:
		# estimation used standardization
		if stdcoefflag==0:
			# Unstandardize coefficients
			b = model.coef_ / x_std * y_std
			if consflag==1:
				# Get constant; use prestandardized means
				predict_nostd = np.sum(x_mean * b)
				b = np.append(b,y_mean - predict_nostd)
				b = np.array([b])
			else:
				b = np.array([b])
		else:
			b = np.array([model.coef_])
	Matrix.store("r(beta)",b)

	SFIToolkit.stata("timer off 56")
	
	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 57")
	
	Scalar.setValue("r(alpha)", lratio)
	Scalar.setValue("r(lambda)", lpenalty)

	# Get mae, rmse
	mae = metrics.mean_absolute_error(y, model_predict)
	mse = metrics.mean_squared_error(y, model_predict)
	rmse = np.sqrt(mse)
	if stdcoefflag==0:
		mae = mae * y_std
		rmse = rmse * y_std
	Scalar.setValue("r(mae)", mae, vtype='visible')
	Scalar.setValue("r(rmse)", rmse, vtype='visible')

	# Unstandardize predictions if required
	if stdflag==1 and prediction != "":
		model_predict = model_predict * y_std + y_mean
	
	# Store in supplied variable
	if prediction != "":
		Data.store(var=prediction,obs=None,val=model_predict,selectvar=touse)

	SFIToolkit.stata("timer off 57")

	##############################################################

	# Pass objects back to __main__ namespace to interact w/ later
	__main__.model_object  = model
	__main__.model_y_std = y_std
	__main__.model_y_mean = y_mean
	__main__.model_xvars = xvars
	if stdflag==1:
		__main__.model_x_std = x_std
		__main__.model_x_mean = x_mean

	
end
