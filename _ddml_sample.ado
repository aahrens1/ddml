*! ddml v1.4.2
*! last edited: 8aug2023
*! authors: aa/ms

program _ddml_sample, sortpreserve					//  sortpreserve needed for fold IDs that respect clustering
	version 16

	syntax [if] [in] , mname(name) [				///
							foldvar(varlist)		/// optional list of variables indicating folds, one per rep
							reps(integer 0)			/// default=1 below
							APPEND1					/// option abbrev is "append", allowing both "append" and "append(#)"
							append(integer 0)		///
							NORANDOM				/// first fold ID uses obs in existing order
							vars(varlist)			///
							kfolds(integer 0)		/// default=5 below
							tabfold					///
							]
	
	// incompatible options
	if "`foldvar'"~="" & `reps' {
		di as err "error - incompatible options, foldvar(`foldvar') and reps(`reps')"
		exit 198
	}
	if "`foldvar'"~="" & `kfolds' {
		di as err "error - incompatible options, foldvar(`foldvar') and kfolds(`kfolds')"
		exit 198
	}
	if `reps' & `append' {
		di as err "error - incompatible options, append(`append') and reps(`reps')"
		exit 198
	}
	if "`append1'"~="" & `append' {
		di as err "error - incompatible options, append and append(`append')"
		exit 198
	}
	if "`append1'"~="" & "`foldvar'"=="" {
		di as err "error - append option with no argument requires foldvar(varlist)"
		exit 198
	}
	
	// syntax checks and defaults
	if `kfolds'<0 {
		di as err "error - invalid kfolds(`kfolds'); must be an integer > 1"
		exit 198
	}
	else if `kfolds'==0 & "`foldvar'"=="" {
		// default number of folds unless foldvar is provided
		local kfolds=5
	}
	if `reps'<0 {
		di as err "error - invalid reps(`reps'); must be an integer > 0"
		exit 198
	}
	else if `reps'==0 & "`foldvar'"=="" {
		// default number of folds unless foldvar is provided
		local kfold	=5
		local reps	=1
	}
	else if `reps'==0 {
		local reps : word count `foldvar'
	}
	
	
	// update append macro; append1 macro not needed after this
	if "`append1'"~="" {
		// update append macro to have number of appended resamples from #foldvars
		local append : word count `foldvar'
	}
	
	// if appending, reps and kfolds = current setting for model
	if `append' {
		mata: st_local("reps", strofreal(`mname'.nreps))
		mata: st_local("kfolds", strofreal(`mname'.kfolds))
	}
	
	// clear all results (crossfits and estimation) or just estimations
	if `append'==0 {
		// clear any preexisting equation results from the model struct
		mata: clear_model_results(`mname')
	}
	else {
		// keep crossfit equations but clear estimation results
		mata: clear_model_estimation(`mname')
	}
	// set reps, firstrep, lastrep
	if `append' {
		// if appending, reps=prev reps setting for model, firstrep=reps+1, lastrep=reps+append
		local firstrep	= `reps'+1
		local lastrep	= `reps'+`append'
	}
	else {
		// will create full set of foldvars for all resamples
		local firstrep	= 1
		local lastrep	= `reps'
	}
	
	// estimation sample	
	marksample touse
	if "`vars'" ~= "" & `append'==0 {
		// set sample indicator to 0 if obs have missings, unless appending to existing model
		fvunab vars : `vars'
		markout `touse' `vars'
		// add list of vars to model struct
		mata: `mname'.strDatavars = "`vars'"
	}
	if `append' {
		// replace touse macro with pre-existing sample variable
		local touse = `mname'_sample
	}
	
	// tempvar fclustid is either defined using fclustvar or is equal to _n.
	tempvar fclustid
	mata: st_local("fclustvar", `mname'.fclustvar)
	if "`fclustvar'"=="" {
		qui gen double `fclustid'=_n
	}
	else {
		qui egen double `fclustid' = group(`fclustvar')
	}
	
	*** gen folds
	// create foldvar
	// Stata name will be mname_fid with _m as rep extension
	
	// delete existing foldvars, unless appending to existing model
	if `append'==0	cap drop `mname'_fid*
	
	if "`foldvar'"=="" {
		forvalues m=`firstrep'/`lastrep' {
			*** gen folds
			tempvar uni cuni tag
			// tag one ob per fold cluster; if no clustering, all obs are tagged
			qui egen `tag' = tag(`fclustid') if `mname'_sample
			if `m'==1 & "`norandom'"~="" {
				qui gen `uni' = _n if `mname'_sample & `tag'
				local labtext "Fold ID (original order), rep `m'"
			}
			else {
				qui gen double `uni' = runiform() if `mname'_sample & `tag'
				local labtext "Fold ID (randomly generated) rep `m'"
			}
			qui cumul `uni' if `mname'_sample, gen(`cuni')
			// create equal-sized folds (#obs or #cluster)
			qui egen long `mname'_fid_`m' = cut(`uni'), group(`kfolds')
			sort `fclustid' `tag'
			// propagate random uniforms (last ob in fcluster) within fclusters
			qui by `fclustid': replace `mname'_fid_`m'=`mname'_fid_`m'[_N] if `mname'_sample
			qui replace `mname'_fid_`m' = `mname'_fid_`m' + 1
			label var `mname'_fid_`m' "`labtext'"
		}
	}
	else {
		local pos = 1
		forvalues m=`firstrep'/`lastrep' {
			local vname : word `pos' of `foldvar'
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
			qui egen long `mname'_fid_`m' = group(`vname')
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
			local ++pos			
		}
		// update sample indicator to account for missing fold vars
		qui replace `mname'_sample = `touse'
	}

	// update model struct
	if `append' {
		local reps = `reps'+`append'
	}
	else {
		mata: `mname'.kfolds = `kfolds'
	}
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
