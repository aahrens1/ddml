program define pylasso2, eclass
version 16.0
syntax varlist(min=2) [if] [in] [aweight fweight],	 ///
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

	**** code from pyforest ****
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
	
	* Pass varlist into varlists called yvar and xvars
	gettoken yvar xvars : varlist
	local num_features : word count `xvars'
	
	marksample touse
	qui count if `touse'
	local N		= r(N)

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

	python: run_elastic_net(				///
		`alpha',							///
		`lambda',							///
		"`training_var' `yvar' `xvars'",	///
		"`touse'",							///
		`n_jobs',							///
		`random_state',						///
		`verbose',							///
		`warm_start',						///
		"`prediction'",						///
		"`training_var'",					///
		 `nonempty_test',					///
		 `stdflag',							///
		 `stdcoefflag',						///
		 `consflag',						///
		 `normflag')

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
	ereturn di, noomitted
	di "sklearn lambda=" _col(20) %10.4f `sk_lambda'
	di "sklearn alpha=" _col(20) %10.4f `sk_alpha'

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
from sfi import Data,Matrix,Scalar
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

	#-------------------------------------------------------------
	# Load data from Stata into Python
	#-------------------------------------------------------------
	
	# Load into Pandas data frame
	df = DataFrame(Data.get(vars,selectvar=touse))
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

	# Misc
	y_mean = np.nanmean(y_insample)
	y_scale = np.nanstd(y_insample)
	x_mean = np.nanmean(x_insample,axis=0)
	
	# Always standardize y
	y_insample = (y_insample-y_mean)/y_scale
	# Need to rescale lambda penalty
	lpenalty = lpenalty/y_scale
	
	# If standardizing, scale features to mean zero, std dev one
	if stdflag==1:
		scaler = preprocessing.StandardScaler().fit(x_insample)
		x_insample = scaler.transform(x_insample)
		x_insample = DataFrame(x_insample)
		x_insample.columns = features
		x = scaler.transform(x)
		x = DataFrame(x)
		x.columns = features
		x_scale = scaler.scale_
	
	Scalar.setValue("r(alpha)", lratio)
	Scalar.setValue("r(lambda)", lpenalty)

	##############################################################

	# Initialize model object
	if lratio>0:
		model = ElasticNet(alpha=lpenalty, l1_ratio=lratio, random_state=0, fit_intercept=consflag, normalize=normflag, tol=1e-10)
	else:
		model = Ridge(alpha=lpenalty, fit_intercept=consflag, normalize=normflag, tol=1e-10)

	# Train model on training data
	model.fit(x_insample, y_insample)

	# Get in-sample prediction
	pred_insample = model.predict(x_insample)

	# Get full-sample prediction
	model_predict = model.predict(x)


	# Stata coefficient vectors are row vectors 
	if stdflag==0:
		if consflag==1:
			b = np.append(model.coef_*y_scale,(model.intercept_*y_scale + y_mean))
		else:
			b = model.coef_
		b = np.array([b])
	else:
		# estimation used standardization
		if stdcoefflag==0:
			# Unstandardize coefficients
			b = model.coef_ / x_scale * y_scale
			if consflag==1:
				# Get constant; use prestandardized data
				pred_insample_nostd = np.nansum(df_train[features] * b)/np.shape(x_insample)[0]
				b = np.append(b,y_mean - pred_insample_nostd)
				b = np.array([b])
			else:
				b = np.array([b])
		else:
			b = np.array([model.coef_])
	Matrix.store("r(beta)",b)

	# Get in sample (training sample) mae, rmse
	insample_mae = metrics.mean_absolute_error(y_insample, pred_insample)
	insample_mse = metrics.mean_squared_error(y_insample, pred_insample)
	insample_rmse = np.sqrt(insample_mse)
	if stdcoefflag==0:
		insample_mae = insample_mae * y_scale
		insample_rmse = insample_rmse * y_scale
	Scalar.setValue("r(training_mae)", insample_mae, vtype='visible')
	Scalar.setValue("r(training_rmse)", insample_rmse, vtype='visible')

	# If nonempty test sample, get out of sample stats
	if nonempty_test==1:
		x_test = df_test[features]
		y_outsample = df_test[y]
		if stdflag==1:
			y_outsample = (y_outsample-y_mean)/y_scale
			x_test = (x_test-x_mean)/x_scale
		pred_outsample = model.predict(x_test)
		outsample_mae = metrics.mean_absolute_error(y_outsample, pred_outsample)
		outsample_mse = metrics.mean_squared_error(y_outsample, pred_outsample)
		outsample_rmse = np.sqrt(outsample_mse)
		if stdflag==1 and stdcoefflag==0:
			outsample_mae = outsample_mae * y_scale
			outsample_rmse = outsample_rmse * y_scale
		Scalar.setValue("r(test_mae)", outsample_mae, vtype='visible')
		Scalar.setValue("r(test_rmse)", outsample_rmse, vtype='visible')

	# Pass objects back to __main__ namespace to interact w/ later
	# Code from pyforest; not in use
	__main__.model_object  = model
	__main__.model_predict = model_predict
	__main__.features = features

	##############################################################

	
end
