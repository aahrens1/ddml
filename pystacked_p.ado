program define pystacked_p, eclass
	version 16.0
	syntax anything(id="argument name" name=arg) [if] [in], [pr xb]
	
	* Mark sample with if/in
	marksample touse, novarlist
	
	* Count number of variables
	local numVars : word count `arg'
	if `numVars'!=1 {
		di as error "Error: More than 1 prediction variable specified"
		exit 1
	}
	
	* Define locals prediction, features
	local predict_var "`arg'"
	local features "`e(features)'"
	
	* Check to see if variable exists
	cap confirm new variable `predict_var'
	if _rc>0 {
		di as error "Error: prediction variable `predict_var' could not be created - probably already exists in dataset."
		di as error "Choose another name for the prediction."
		exit 1
	}
	
	* Get predictions
	python: post_prediction("`features'","`predict_var'")

	* Keep only if touse
	qui replace `predict_var'=. if `touse'==0
	
	
end

python:

# Import SFI, always with stata 16
from sfi import Data,Matrix,Scalar
from pandas import DataFrame

def post_prediction(vars, prediction):

	# Start with a working flag
	Scalar.setValue("import_success", 1, vtype='visible')

	# Import model from Python namespace
	try:
		from __main__ import model_predict as pred
		from __main__ import model_object as model
	except ImportError:
		print("Error: Could not find estimation results. Run a pylearn command before loading this.")
		Scalar.setValue("import_success", 0, vtype='visible')
		return

	# Load data into Pandas data frame
	df = DataFrame(Data.get(vars))
	colnames = []
	for var in vars.split():
		colnames.append(var)
	df.columns = colnames
	
	# Create list of feature names
	features = df.columns[0:]
	
	# Generate predictions (on both training and test data)
	pred    = model.predict(df[features])

	# Export predictions back to Stata
   	Data.addVarFloat(prediction)
	Data.store(prediction,None,pred)
	
end