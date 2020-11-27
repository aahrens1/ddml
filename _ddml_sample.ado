program _ddml_sample
	version 13

	syntax [if] [in] , mname(name) [ foldvar(varname) vars(varlist) kfolds(integer 2) tabfold ]
	
	marksample touse

	cap drop `mname'_esample
	if "`vars'" ~= "" {
		// set sample indicator to 0 if obs have missings
		fvunab vars : `vars'
		markout `touse' `vars'
	}

	*** gen folds
	// create foldvar
	// Stata name will be mname_fid
	cap drop `mname'_fid
	if "`foldvar'"=="" {
		tempvar uni cuni
		qui gen double `uni' = runiform() if `touse'
		qui cumul `uni' if `touse', gen(`cuni')
		qui gen int `mname'_fid = ceil(`kfolds'*`cuni') if `touse'
	}
	else {
		// check that fold var is legit
		cap count if `foldvar' < . & `touse'
		if _rc > 0 {
			di as err "error - fold variable `foldvar' does not exist or is not a valid identifier"
			exit 198
		}
		qui count if `foldvar'==. & `touse'
		if r(N)>0 {
			di as res "note - fold variable missing for some observations"
			di as res "these observations will be excluded from the estimation sample"
			qui replace `touse' = 0 if `foldvar'==.
		} 
		qui gen `mname'_fid = `foldvar'
	}

	// add sample indicator to model struct (col 1 = id, col 2 = fold id)
	qui gen byte `mname'_esample = `touse'
	mata: `mname'.strDatavars = "`vars'"
	mata: `mname'.idSample = st_data(., ("`mname'_id", "`mname'_esample"))
	// add fold id to model struct (col 1 = id, col 2 = fold id)
	mata: `mname'.idFold = st_data(., ("`mname'_id", "`mname'_fid"))

	if ("`tabfold'"!="") {
		di
		di "Overview of frequencies by fold and sample:"
		tab `mname'_fid `mname'_esample, miss
		di
	}


end
