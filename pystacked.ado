program define pystacked, eclass
version 16.0
syntax varlist(min=2 fv) [if] [in] [aweight fweight],	 	///
					type(string) /// classification or regression
				[									 		///
					lasso /// 
					rf /// 
					gradboost ///
				]

	// pylearn, check
	
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

	python: run_stacked(					///
		"`type'",								///
		"`lasso' `rf' `gradboost'", ///
		"`training_var' `yvar_t' `xvars_t'",	///
		"`training_var'", ///
		"`touse'")

	ereturn local predict "pystacked_p"

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
from sklearn.ensemble import GradientBoostingClassifier,GradientBoostingRegressor
from sklearn.ensemble import StackingRegressor
from sfi import Data,Matrix,Scalar,SFIToolkit
from sklearn import metrics, preprocessing
from sklearn.preprocessing import StandardScaler
from sklearn.linear_model import LassoCV,LassoLarsIC,RidgeCV
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
# Define Python function: run_stacked
#-------------------------------------------------------------------------------

def run_stacked(type,methods,vars,training,touse):

	##############################################################
	### load data  											   ###
	##############################################################	

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

	
	##############################################################
	### fitting												   ###
	##############################################################

	methods = methods.split()
	est_list = []
	if type=="regress":
		if "lasso" in methods:
			est_list.append(('lasso',make_pipeline(StandardScaler(), LassoLarsIC())))
		if "rf" in methods:
			est_list.append(('rf',RandomForestRegressor()))
		if "gradboost" in methods:
			est_list.append(('rf',GradientBoostingRegressor()))
		fin_est = RidgeCV()
	elif type=="classify":
		if "rf" in methods:
			est_list.append(('rf',RandomForestClassifier()))
		if "gradboost" in methods:
			est_list.append(('rf',GradientBoostingClassifier()))
		fin_est = LogisticRegression()

	model = StackingRegressor(
	   estimators=est_list,
	   final_estimator=fin_est
	)
	
	#-------------------------------------------------------------
	# Fit model, get predictions, pass objects back to main
	#-------------------------------------------------------------

	# Train model on training data
	model.fit(x_insample, y_insample)

	# Get in-sample prediction
	pred_insample = model.predict(x_insample)

	# Get full-sample prediction
	model_predict = model.predict(x)

	# Pass objects back to __main__ namespace to interact w/ later
	__main__.model_object  = model
	__main__.model_predict = model_predict
	__main__.model_transform = model.transform(x)
	__main__.features = features
	
end
