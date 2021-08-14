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
			}
			else {
				qui gen double `uni' = runiform() if `mname'_sample
			}
			qui cumul `uni' if `mname'_sample, gen(`cuni')
			qui gen int `mname'_fid_`m' = ceil(`kfolds'*`cuni') if `mname'_sample
			// add fold id to model struct (col 1 = id, col 2 = fold id)
			// mata: `mname'.idFold = (`mname'.idFold , st_data(., ("`mname'_fid_`m'")))
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
			label var `mname'_fid_`m'
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
		// add fold id to model struct (col 1 = id, col 2 = fold id)
		//mata: `mname'.idFold = (`mname'.idFold , st_data(., ("`mname'_fid_`m'")))
	}

	// add sample indicator to model struct (col 1 = id, col 2 = fold id)
	// mata: `mname'.idSample = st_data(., ("`mname'_id", "`mname'_sample"))
	
	// update model struct
	mata: `mname'.kfolds = `kfolds'
	mata: `mname'.nreps = `reps'

	// clear previous results (incomplete)
	mata: clear_model_results(`mname',`kfolds',`reps')

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

void clear_model_results(		struct ddmlStruct m,	///
								real scalar kfolds,		///
								real scalar reps		///
								)

{
	pointer(struct eqnStruct) scalar p
	
	m.crossfitted	= 0

	for (i=1; i<=cols(m.eqnlist); i++) {
		p					= m.eqnlist[i]
		(*p).MSE			= J(0,1,0)
		(*p).MSE_h			= J(0,1,0)
		(*p).MSE0			= J(0,1,0)
		(*p).MSE1			= J(0,1,0)
		(*p).N				= J(0,1,0)
		(*p).N_h			= J(0,1,0)
		(*p).N0				= J(0,1,0)
		(*p).N0				= J(0,1,0)

		(*p).MSE_folds		= J(0,kfolds,0)
		(*p).MSE_h_folds	= J(0,kfolds,0)
		(*p).MSE0_folds		= J(0,kfolds,0)
		(*p).MSE1_folds		= J(0,kfolds,0)
		(*p).N_folds		= J(0,kfolds,0)
		(*p).N_h_folds		= J(0,kfolds,0)
		(*p).N0_folds		= J(0,kfolds,0)
		(*p).N1_folds		= J(0,kfolds,0)
	}

}

end

/*
	real colvector		MSE
	real matrix			MSE_folds		// MSE by fold; col=fold, row=resample
	real colvector      MSE_h 			// (intended for LIE)
	real matrix			MSE_h_folds		// (intended for LIE)
	real colvector 		MSE0
	real colvector 		MSE1
	real matrix			MSE0_folds		// MSE by fold; col=fold, row=resample
	real matrix			MSE1_folds		// MSE by fold; col=fold, row=resample
	real colvector		N
	real matrix			N_folds			// sample size by fold; col=fold, row=resample
	real colvector		N_h				// (intended for LIE)
	real matrix         N_h_folds		// (intended for LIE)
	real colvector		N0
	real colvector		N1
	real matrix			N0_folds		// sample size by fold; col=fold, row=resample
	real matrix			N1_folds		// sample size by fold; col=fold, row=resample
*/