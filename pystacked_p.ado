program define pystacked_p, rclass
	version 16.0
	syntax anything(id="argument name" name=arg) [if] [in], [ ///
															pr /// not implemented yet
															xb /// default
															r /// not implemented yet
															TRANSForm ///
															]
	
	* Count number of variables
	local numVars : word count `arg'
	if `numVars'!=1 {
		di as error "Error: More than 1 prediction variable specified"
		exit 1
	}
	
	* Define locals prediction, features
	local predict_var "`arg'"
	
	/* Check to see if variable exists
	cap confirm new variable `predict_var'
	if _rc>0 {
		di as error "Error: prediction variable `predict_var' could not be created - probably already exists in dataset."
		di as error "Choose another name for the prediction."
		exit 1
	}
	*/

	python: get_dim()
	local n_methods = r(n_methods)

	if "`transform'"!="" {
		local transform_tempvars
		forvalues i = 1/`n_methods' {
			tempvar pred_var`i'
			local transform_tempvars `transform_tempvars' `pred_var`i''
		}
	}
	
	if "`xb'"!="" {
		* Keep only if touse
		tempvar temp_predict
		qui gen double `temp_predict' = .
	}

	* Get predictions
	python: post_prediction("`temp_predict'","`transform_tempvars'")
	
	if "`xb'"!="" {
		qui gen double `predict_var'= `temp_predict'  `if' `in'
		label var `predict_var' "optimal predictions"
	}

	if "`transform'"!="" {
		local j = 1
		foreach var of varlist `transform_tempvars' {
			qui gen double `predict_var'`j'= `var' `if' `in'
			local m: word `j' of `methods'
			label var `predict_var'`j' "predictions `m'"
			local j =`j'+1
		}
	}
	
end

python:

# Import SFI, always with stata 16
from sfi import Data,Matrix,Scalar,Macro
from pandas import DataFrame

def post_prediction(pred_var,transf_vars):

	# Start with a working flag
	Scalar.setValue("r(import_success)", 1, vtype='visible')

	# Import model from Python namespace
	try:
		from __main__ import model_predict as pred
		from __main__ import model_object as model
		from __main__ import model_touse as touse
		from __main__ import model_id as id
		from __main__ import model_transform as transform
		from __main__ import methods as methods
	except ImportError:
		print("Error: Could not find estimation results. Run a pylearn command before loading this.")
		Scalar.setValue("r(import_success)", 0, vtype='visible')
		return

	if transf_vars!="":
		transf_vars = transf_vars.split()
		ncol = transform.shape[1]
		for j in range(ncol):
			Data.addVarDouble(transf_vars[j])
			Data.setVarLabel(transf_vars[j],"Prediction"+" "+methods[j])
			Data.store(var=transf_vars[j],val=transform[:,j],obs=None,selectvar=touse)

	if pred_var!="":
		Data.store(var=pred_var,val=pred,obs=None,selectvar=touse)

def get_dim():

	# Start with a working flag
	Scalar.setValue("r(import_success)", 1, vtype='visible')

	# Import model from Python namespace
	try:
		from __main__ import methods as methods
	except ImportError:
		print("Error: Could not find estimation results. Run a pylearn command before loading this.")
		Scalar.setValue("r(import_success)", 0, vtype='visible')
		return

	Scalar.setValue("r(n_methods)", len(methods))
	Macro.setLocal("methods"," ".join(methods))

end


