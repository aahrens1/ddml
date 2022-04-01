program _ddml_sample
	version 13

	syntax [if] [in] , mname(name) [				///
							foldvar(varlist)		///
							reps(integer 0)			///
							NORANDOM				/// first fold ID uses obs in existing order
							vars(varlist)			///
							kfolds(integer 5)		///
							tabfold					///
							sreset					///
							]
	
	marksample touse
	
	// clear any preexisting equation results from the model struct
	mata: clear_model_results(`mname')

	if "`sreset'"~="" {
		// reset sample indicator
		cap drop `mname'_sample
		qui gen byte `mname'_sample = 1
	}
	
	if "`vars'" ~= "" {
		// set sample indicator to 0 if obs have missings
		fvunab vars : `vars'
		markout `touse' `vars'
		// add list of vars to model struct
		mata: `mname'.strDatavars = "`vars'"
	}
	
	*** gen folds
	// create foldvar
	// Stata name will be mname_fid with _m as rep extension
	cap drop `mname'_fid*
	// incompatible options
	if "`foldvar'"~="" & `reps'>0 {
		di as err "error - incompatible options, foldvar(`foldvar') and reps(`reps')"
		exit 198
	}
	if "`foldvar'"=="" {
		// if neither foldvar nor reps provided, set reps to default
		if `reps'==0 {
			local reps 1
		}
		forvalues m=1/`reps' {
			*** gen folds
			tempvar uni cuni
			if `m'==1 & "`norandom'"~="" {
				qui gen `uni' = _n
				local labtext "(original order)"
			}
			else {
				qui gen double `uni' = runiform() if `mname'_sample
				local labtext "(randomly generated)"
			}
			qui cumul `uni' if `mname'_sample, gen(`cuni')
			// epsfloat() used below so that if cuni is on boundary it goes to the lower group
			qui gen int `mname'_fid_`m' = ceil(`kfolds'*(`cuni'-epsfloat())) if `mname'_sample
			label var `mname'_fid_`m' "`labtext'"
		}
	}
	else {
		// update reps
		local reps : word count `foldvar'
		tokenize `foldvar'
		forvalues m=1/`reps' {
			local vname ``m''
			// check that fold var is legit
			cap count if `vname' < . & `touse'
			if _rc > 0 {
				di as err "error - fold variable `foldvar' does not exist or is not a valid identifier"
				exit 198
			}
			qui count if `vname'==. & `touse'
			if r(N)>0 {
				di as res "note - fold variable missing for some observations"
				di as res "these observations will be excluded from the estimation sample"
				qui replace `touse' = 0 if `vname'==.
			}
			qui egen `mname'_fid_`m' = group(`vname')
			label var `mname'_fid_`m' "(based on `vname')"
			// enforce that the number of folds is the same for all fold vars
			qui tab `mname'_fid_`m'
			if `m'==1 {
				// initialize kfolds for checking
				local kfolds = r(r)
			}
			else {
				if r(r)~=`kfolds' {
					di as err "error - fold variables must have same number of folds"
					exit 198
				}
			}
			
		}
		// update sample indicator to account for missing fold vars
		qui replace `mname'_sample = `touse'
	}

	// update model struct
	mata: `mname'.kfolds = `kfolds'
	mata: `mname'.nreps = `reps'

	forvalues m=1/`reps' {
		if ("`tabfold'"!="") {
			di
			di "Overview of frequencies by fold (sample `m'):"
			tab `mname'_fid_`m' if `mname'_sample
			di
		}
	}


end

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end
