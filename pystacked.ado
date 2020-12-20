program define pystacked, eclass
version 16.0
syntax varlist(min=2 fv) [if] [in] [aweight fweight],	 	///
					type(string) /// classification or regression
				[									 		///
					LASSO /// 
					RF /// 
					GRADboost ///
					SVM /// sklearn.svm.SVM
					LINSvm /// sklearn.svm.LinearSVR
					svmopt(string asis) ///
					lassoopt(string asis) ///
					rfopt(string asis) ///
					gradboostopt(string asis) ///
					linsvmopt(string asis) ///
					seed(integer 0) ///
					printopt ///
				]

	// pylearn, check

	*** lasso
	if ("`lassoopt'"!="") {
		local lasso lasso
	}
	parse_LassoIC `lassoopt'  
	local lassoopt `r(optstr)'

	*** random forest
	if ("`rfopt'"!="") {
		local rf rf
	}
	if "`rf'"!="" {
		_pyparse , type(`type') method(rf)
		local rfopt `r(optstr)'
	}

	*** SVM
	if ("`svmopt'"!="") {
		local svm svm
	}
	if "`svm'"!="" {
		_pyparse , type(`type') method(svm)
		local svmopt `r(optstr)'
	}

	*** gradboost 
	if ("`gradboostopt'"!="") {
		local gradboost gradboost
	}
	if "`gradboost'"!="" {
		_pyparse , type(`type') method(gradboost)
		local gradboostopt `r(optstr)'
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

	tempvar idvar
	gen byte `idvar'=_n

	python: run_stacked(					///
		"`type'",								///
		"`lasso' `rf' `gradboost'", ///
		"`training_var' `yvar_t' `xvars_t'",	///
		"`training_var'", ///
		"`lassoopt'", ///
		"`gradboostopt'", ///
		"`rfopt'", ///
		"`svmopt'", ///
		"`touse'", ///
		"`idvar'", ///
		`seed', ///
		1)

	if ("`lasso'"!=""&"`printopt'"!="") di as text "`lassoopt'"
	if ("`gradboost'"!=""&"`printopt'"!="") di as text "`gradboostopt'"
	if ("`rf'"!=""&"`printopt'"!="") di as text "`rfopt'"
	if ("`svm'"!=""&"`printopt'"!="") di as text "`svmopt'"

	ereturn local predict "pystacked_p"
	if ("`lasso'"!="") ereturn local lassoopt `lassoopt'
	if ("`gradboost'"!="") ereturn local gradboostopt `gradboostopt'
	if ("`rf'"!="") ereturn local rfopt `rfopt'
	if ("`svm'"!="") ereturn local svmopt `svmopt'
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
from sklearn.linear_model import ElasticNet
from sklearn.linear_model import Ridge
from sklearn.pipeline import make_pipeline
from sklearn.svm import SVC

# To pass objects to Stata
import __main__


import random

#-------------------------------------------------------------------------------
# Define Python function: run_stacked
#-------------------------------------------------------------------------------

def run_stacked(type,methods,vars,training,lassoopt,gradboostopt,rfopt,svmopt,touse,idvar,seed,save_transform):

	# Set random seed
	if seed>0:
		random.seed(seed)

	##############################################################
	### load data  											   ###
	##############################################################	

	# Load into Pandas data frame
	df = DataFrame(Data.get(vars,selectvar=touse))
	id = Data.get(idvar,selectvar=touse)

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

	lassoopt = eval(lassoopt)
	gradboostopt = eval(gradboostopt)
	rfopt = eval(rfopt)
	svmopt = eval(svmopt)

	methods = methods.split()
	est_list = []
	if type=="regress":
		if "lasso" in methods:
			est_list.append(('lasso',make_pipeline(StandardScaler(), LassoLarsIC(**lassoopt))))
		if "rf" in methods:
			est_list.append(('rf',RandomForestRegressor(**rfopt)))
		if "gradboost" in methods:
			est_list.append(('gradboost',GradientBoostingRegressor(**gradboostopt)))
		if "svm" in methods:
			est_list.append(('svm',make_pipeline(StandardScaler(), SVR(**svmopt))))
		fin_est = RidgeCV()
	elif type=="classify":
		if "lasso" in methods:
			est_list.append(('lasso',make_pipeline(StandardScaler(), LassoLarsIC(**lassoopt))))
		if "rf" in methods:
			est_list.append(('rf',RandomForestClassifier(**rfopt)))
		if "gradboost" in methods:
			est_list.append(('gradboost',GradientBoostingClassifier(**gradboostopt)))
		if "svm" in methods:
			est_list.append(('svm',make_pipeline(StandardScaler(), SVM(**svmopt))))
		fin_est = LogisticRegression()

	model = StackingRegressor(
	   estimators=est_list,
	   final_estimator=fin_est
	)
	
	for i in range(len(est_list)):
		print(est_list[i][1].get_params())

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
	__main__.model_touse = touse
	__main__.model_id = id
	__main__.features = features
	__main__.methods = methods
	if save_transform==1:
		__main__.model_transform = model.transform(x)
	
end
