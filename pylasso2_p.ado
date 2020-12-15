program define pylasso2_p, rclass
	version 16.0
	syntax namelist(min=1 max=2) [if] [in], [		///
											xb		/// default
											Resid	/// not implemented yet
											]
	
	* get var type & name
	tokenize `namelist'
	if "`2'"=="" {					//  only new varname provided
		local predictvar `1'
		qui gen `predictvar' = .
	}
	else {							//  datatype also provided
		local vtype `1'
		local predictvar `2'
		qui gen `vtype' `predictvar' = .
	}
	*

	local command=e(cmd)
	if ("`command'"~="pylasso2") {
		di as err "error: -pylasso2_p- supports only the -pylasso2- command"
		exit 198
	}
	*
	
	marksample touse, novarlist

	* Get predictions (fitted values)
	python: post_prediction("`predictvar'",`e(std)')
	
	* Replace with residuals if requested
	if "`resid'"~="" {
		qui replace `predictvar' = `e(depvar)' - `predictvar'
	}
	
	* Set any obs not to be used to missing
	qui replace `predictvar' = . if `touse'==0

end

python:

# Import SFI, always with stata 16
from sfi import Data,Matrix,Scalar,Macro
from pandas import DataFrame
import numpy as np

def post_prediction(pred_var,stdflag):

	# Start with a working flag
	Scalar.setValue("r(import_success)", 1, vtype='visible')

	# Import model from Python namespace
	try:
		from __main__ import model_object as model
		from __main__ import model_xvars as xvars
		from __main__ import model_x_std as x_std
		from __main__ import model_x_mean as x_mean
		from __main__ import model_y_std as y_std
		from __main__ import model_y_mean as y_mean

	except ImportError:
		print("Error: Could not find pylasso2 estimation results.")
		Scalar.setValue("r(import_success)", 0, vtype='visible')
		return

	# Load data as an array
	x = Data.get(xvars,missingval=np.nan)

	# Track NaNs
	x_hasnan = np.isnan(x).any(axis=1)
	
	# Prepare data including standardization if required
	if stdflag==1:
		x = (x-x_mean)/x_std

	# Set any NaNs to zeros so that model.predict(.) won't crash
	x = np.nan_to_num(x)
	
	# Get predictions
	pred = model.predict(x)
	
	# Unstandardize and remean
	pred = pred * y_std + y_mean
	
	# Set any predictions that should be missing to missing (NaN)
	pred[x_hasnan] = np.nan

	# Put predictions into the designated Stata variable
	Data.store(var=pred_var,val=pred,obs=None)

end


