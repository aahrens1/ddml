program define pylasso2, eclass
version 16.0
syntax varlist(min=2 fv) [if] [in] [aweight fweight],	 ///
				[									 ///
					alpha(real 1)					 /// elastic net parameter, alpha=1 is lasso, alpha=0 is ridge
					lambda(real 1)					 /// penalty
					n_jobs(integer -1)		 	 	 /// number of processors to use when computing stuff - default is all
					random_state(integer -1) 	 	 /// seed used by random number generator
					verbose		 			 	 	 /// controls verbosity
					warm_start(string asis)	 	 	 /// when set to true, reuse solution of previous call to fit
				    training(varname) 	             /// training dataset identifier
					prediction(string asis) 	     /// variable name to save predictions
					standardize                  	 /// standardize features
					normalize						 /// use scikit's normalize
					UNITLoadings					 /// don't prestandardize
					STDCoef							 /// report standardized coefficients
					NOCONStant						 /// standard option
				]

	// pylearn, check
	
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
	
	* Define a temporary variable for the training sample
	tempvar training_var
	if "`training'"=="" {
		qui gen `training_var' = 1 if `touse'
	}
	else {
		qui gen `training_var' = `training' if `touse'
	}

	qui count if `training_var'==0 & `touse'
	local nonempty_test = r(N)>0
	
	// create an empty prediction variable
	if "`prediction'"~="" {
		qui gen double `prediction' = .
	}

	timer off 40
	timer on 50
	
	python: run_elastic_net(					///
		`alpha',								///
		`lambda',								///
		"`training_var' `yvar_t' `xvars_t'",	///
		"`touse'",								///
		`n_jobs',								///
		`random_state',							///
		`verbose',								///
		`warm_start',							///
		"`prediction'",							///
		"`training_var'",						///
		 `nonempty_test',						///
		 `stdflag',								///
		 `stdcoefflag',							///
		 `consflag',							///
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
	
	local training_mae	= r(training_mae)
	local training_rmse	= r(training_rmse)
	if `nonempty_test' {
		local test_mae	= r(test_mae)
		local test_rmse	= r(test_rmse)
	}
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

	ereturn scalar training_mae		= `training_mae'
	ereturn scalar training_rmse	= `training_rmse'
	if `nonempty_test' {
		ereturn scalar test_mae		= r(test_mae)
		ereturn scalar test_rmse	= r(test_rmse)
	}

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
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Ridge
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

def run_elastic_net(lratio,lpenalty,vars,touse,n_jobs,random_state,verbose,warm_start,prediction,training,nonempty_test,stdflag,stdcoefflag,consflag,normflag):

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 50")

	# Load into Pandas data frame
	df = DataFrame(Data.get(vars,selectvar=touse))
	
	SFIToolkit.stata("timer off 50")

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 51")

	colnames = []
	for var in vars.split():
		 colnames.append(var)
	df.columns = colnames

	# Split training data and test data into separate data frames
	df_train, df_test = df[df[training]==1], df[df[training]==0]
	
	# Create list of feature names
	features = df.columns[2:]
	y        = df.columns[1]

	# Split training data frame into features (x) and outcome (y)
	x           = df[features]
	x_insample  = df_train[features]
	y_insample  = df_train[y]

	SFIToolkit.stata("timer off 51")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 52")

	# Misc
	n_insample = x_insample.shape[0]
	y_mean = np.mean(y_insample)
	y_std = np.std(y_insample)
	x_mean = np.mean(x_insample,axis=0)
	x_std = np.std(x_insample,axis=0)
	
	# If SD is zero, replace with 1
	x_std[x_std==0] = 1
	
	# Always standardize y
	y_insample = (y_insample-y_mean)/y_std
	
	# Need to rescale lambda penalty
	lpenalty = lpenalty/y_std
	
	# If standardizing, scale features to mean zero, std dev one
	if stdflag==1:
		x_insample = (x_insample-x_mean)/x_std
		x = (x-x_mean)/x_std
	if nonempty_test==1:
		x_test = df_test[features]
		y_test = df_test[y]
		if stdflag==1:
			x_test = (x_test-x_mean)/x_std
			y_test = (y_test-y_mean)/y_std

	SFIToolkit.stata("timer off 52")

	##############################################################
	
	# Start Stata timer
	SFIToolkit.stata("timer on 53")

	# Initialize model object
	if lratio>0:
		model = ElasticNet(alpha=lpenalty, l1_ratio=lratio, random_state=0, fit_intercept=consflag, normalize=normflag, tol=1e-10)
	else:
		# Ridge uses a different definition of the penalty
		model = Ridge(alpha=lpenalty*n_insample, fit_intercept=consflag, normalize=normflag, tol=1e-10)

	# Train model on training data
	model.fit(x_insample, y_insample)
	
	SFIToolkit.stata("timer off 53")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 54")
	
	# Get in-sample prediction
	insample_predict = model.predict(x_insample)

	# Get full-sample prediction
	model_predict = model.predict(x)

	SFIToolkit.stata("timer off 54")

	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 55")
	
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
				# Get constant; use prestandardized data
				insample_predict_nostd = np.nansum(df_train[features] * b)/n_insample
				b = np.append(b,y_mean - insample_predict_nostd)
				b = np.array([b])
			else:
				b = np.array([b])
		else:
			b = np.array([model.coef_])
	Matrix.store("r(beta)",b)

	SFIToolkit.stata("timer off 55")
	
	##############################################################

	# Start Stata timer
	SFIToolkit.stata("timer on 56")
	
	Scalar.setValue("r(alpha)", lratio)
	Scalar.setValue("r(lambda)", lpenalty)

	# Get in sample (training sample) mae, rmse
	insample_mae = metrics.mean_absolute_error(y_insample, insample_predict)
	insample_mse = metrics.mean_squared_error(y_insample, insample_predict)
	insample_rmse = np.sqrt(insample_mse)
	if stdcoefflag==0:
		insample_mae = insample_mae * y_std
		insample_rmse = insample_rmse * y_std
	Scalar.setValue("r(training_mae)", insample_mae, vtype='visible')
	Scalar.setValue("r(training_rmse)", insample_rmse, vtype='visible')

	# If nonempty test sample, get out of sample stats
	if nonempty_test==1:
		outsample_predict = model.predict(x_test)
		outsample_mae = metrics.mean_absolute_error(y_test, outsample_predict)
		outsample_mse = metrics.mean_squared_error(y_test, outsample_predict)
		outsample_rmse = np.sqrt(outsample_mse)
		if stdflag==1 and stdcoefflag==0:
			outsample_mae = outsample_mae * y_std
			outsample_rmse = outsample_rmse * y_std
		Scalar.setValue("r(test_mae)", outsample_mae, vtype='visible')
		Scalar.setValue("r(test_rmse)", outsample_rmse, vtype='visible')

	# Unstandardize predictions if required
	if stdflag==1:
		model_predict = model_predict * y_std + y_mean
	
	# Store in supplied variable
	if prediction != "":
		Data.store(var=prediction,obs=None,val=model_predict,selectvar=touse)

	SFIToolkit.stata("timer off 56")

	##############################################################

	
end
