*! ddml v1.2
*! last edited: 8 july 2023
*! authors: aa/ms

program _ddml_estimate_ate_late, eclass sortpreserve
	version 16
	syntax [anything] [if] [in] ,					/// 
								[					///
								y0(varname)			/// for estimating by hand...
								y1(varname)			/// 
								d(varname)			/// 
								d0(varname)			/// 
								d1(varname)			/// 
								z(varname)			///
								shortstack			/// re-short-stack
								poolstack			/// re-pool-stack
								finalest(name)		/// relevant only for re-stacking
								ssfinalest(name)	///
								psfinalest(name)	///
								* ]

	if "`shortstack'"~="" {
		// restack
		if "`ssfinalest'"==""	local ssfinalest `finalest'
		_ddml_estimate_stacking `anything' `if' `in', ss finalest(`ssfinalest') `options'
	}
	if "`poolstack'"~="" {
		// restack
		if "`psfinalest'"==""	local psfinalest `finalest'
		_ddml_estimate_stacking `anything' `if' `in', ps finalest(`psfinalest') `options'
	}
	if "`y0'`y1'`d'`d0'`d1'`z'"=="" {
		// main program for estimation
		_ddml_estimate_main `anything' `if' `in', `options'
	}
	else {
		// a single user-specified estimation
		_ddml_estimate_single `anything' `if' `in', y0(`y0') y1(`y1') d(`d') d0(`d0') d1(`d1') z(`z') `options'
	}
end

// (re-)stack
program _ddml_estimate_stacking, eclass sortpreserve
	version 16
	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								ss					///
								ps					///
								finalest(name)		///
								NOIsily				///
								*					///
								]

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	if "`noisily'"==""	local qui qui
	
	// macro `ts' is either "ss" or "ts"
	// macro `typestack' is either "shortstack" or "poolstack"
	
	if "`ss'"~="" {
		// shortstacking
		local ssflag = 1
		local typestack	shortstack
		local ts ss
	}
	else if "`ps'"~="" {
		// poolstacking
		local ssflag = 0
		local typestack poolstack
		local ts ps
	}
	else {
		// error
		di as err "internal _ddml_estimate_stacking error - missing ss or ts option"
		exit 198
	}
	
	marksample touse

	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	// reps = total number of reps; crossfitted = reps done so far (=0 if none)
	mata: st_local("reps", strofreal(`mname'.nreps))
	mata: st_local("kfolds", strofreal(`mname'.kfolds))
	mata: st_local("crossfitted", strofreal(`mname'.crossfitted))

	// two sets of lists, one for variables that come in 0/1 pairs and one for the rest
	// with pystacked, will always be a single 0/1 vtilde pair and eqn struct for every variable
	// for rest, with pystacked+interactive, always a single vtilde and eqn struct
	
	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	// used for checking
	mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
	if `pystackedmulti' {
		mata: st_local("vtildeY",invtokens(`eqn'.vtlist))
		local namelist01	`nameY'
		local vtlist01		`vtildeY'
	}
	else {
		di as err "error - restacking available only if pystacked is the only learner for each conditional expectation"
		exit 198
	}

	if "`model'"=="interactive" {
		// for interactive model, D is standard, not 0/1
		local treatvar	`nameD'
		mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
		// used for checking
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		if `pystackedmulti' {
			mata: st_local("vtildeD",invtokens(`eqn'.vtlist))
			local namelist	`nameD'
			local vtlist	`vtildeD'
		}
		else {
			di as err "error - restacking available only if pystacked is the only learner for each conditional expectation"
			exit 198
		}
	}
	else {
		// for late model, D is 0/1 and Z is standard
		mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
		// used for checking
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		if `pystackedmulti' {
			mata: st_local("vtildeD",invtokens(`eqn'.vtlist))
			local namelist01	`namelist01' `nameD'
			local vtlist01		`vtlist01' `vtildeD'
		}
		else {
			di as err "error - restacking available only if pystacked is the only learner for each conditional expectation"
			exit 198
		}
		local treatvar	`nameZ'
		mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
		// used for checking
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		if `pystackedmulti' {
			mata: st_local("vtildeZ",invtokens(`eqn'.vtlist))
			local namelist	`nameZ'
			local vtlist	`vtildeZ'
		}
		else {
			di as err "error - restacking available only if pystacked is the only learner for each conditional expectation"
			exit 198
		}
	}

	tempname sweights N_folds mse_folds
	local numvts : word count `vtlist'
	// loop through standard vtildes
	forvalues i=1/`numvts' {
		local vname :	word `i' of `namelist'
		local vtilde :	word `i' of `vtlist'

		mata: `eqn' = (`mname'.eqnAA).get("`vname'")
		mata: st_local("base_est",return_learner_item(`eqn',"`vtilde'","stack_base_est"))
		mata: st_local("stype",return_learner_item(`eqn',"`vtilde'","stack_type"))
		mata: st_local("etype",`eqn'.etype)
		mata: st_local("lieflag", strofreal(`eqn'.lieflag))
		local nlearners : word count `base_est'
		// check if previously stacked; required for poolstacking
		if `ssflag' {
			mata: st_local("shortstack", `eqn'.shortstack)
			if "`shortstack'"=="" {
				// not previously shortstacked, set local and struct field to default
				local shortstack `etype'_`vname'
				mata: `eqn'.shortstack = "`shortstack'"
				local newstack 1
			}
			else {
				local newstack 0
			}
		}
		else {
			mata: st_local("poolstack", `eqn'.poolstack)
			// must previously have been poolstacked for restacking to work
			if "`poolstack'"=="" {
				di as err "error - must poolstack at crossfitting stage first"
				exit 198
			}
			else {
				local newstack 0
			}
		}

		// loop through reps
		forvalues m=1/`reps' {
			local fid `mname'_fid_`m'
			// assemble learner list
			// clear macro
			local learner_list
			forvalues j=1/`nlearners' {
				local learner_list `learner_list' `vtilde'_L`j'_`m'
			}

			// get stacking weights
			tempvar yhat
			if `ssflag' {
				// shortstacking uses crossfit predictions
				`qui' _ddml_nnls `vname' `learner_list', finalest(`finalest') stype(`stype') if `touse'
				local finalest	`e(finalest)'
				mat `sweights'	= e(b)
				qui predict double `yhat'
			}
			else {
				// poolstacking uses stacking CV predictions, stored in a mata struct
				tempname tframe y_stacking_cv
				qui frame pwf
				local cframe `r(currentframe)'
				frame create `tframe'
				frame change `tframe'
				tempname
				mata: `y_stacking_cv' = return_result_item(`eqn',"``typestack''_ps","y_stacking_cv", "`m'")
				// touse ignored at weights stage
				getmata (`fid' `vname' `learner_list')=`y_stacking_cv', force replace
				`qui' _ddml_nnls `vname' `learner_list', finalest(`finalest') stype(`stype')
				mata: `sweights' = st_matrix("e(b)")
				frame change `cframe'
				frame drop `tframe'
				mata: st_matrix("`sweights'",`sweights')
				mat colnames `sweights' =  `learner_list'
				mat score double `yhat' = `sweights' if `touse'
				cap mata: mata drop `y_stacking_cv'
				cap mata: mata drop `sweights'
			}
			
			if `newstack' {
				cap drop ``typestack''_`ts'_`m'
				qui gen double ``typestack''_`ts'_`m' = `yhat'
				label var ``typestack''_`ts'_`m' "Predicted values cond. exp. of `vname' using `typestack'ing"
			}
			else {
				if "`noisily'"~="" {
					di
					di "Existing vs new predicted values:"
					sum ``typestack''_`ts'_`m' `yhat'
				}
				qui replace ``typestack''_`ts'_`m' = `yhat'
			}
			get_stack_stats if `touse', kfolds(`kfolds') fid(`mname'_fid_`m') vname(`nameY') vhat(`yhat')
			local N				= r(N)
			local mse			= r(mse)
			mat `N_folds'		= r(N_folds)
			mat `mse_folds'		= r(mse_folds)
			
			// to store:
			mata: add_result_item(`eqn',"``typestack''_`ts'","N",            "`m'", `N')
			mata: add_result_item(`eqn',"``typestack''_`ts'","N_folds",      "`m'", st_matrix("`N_folds'"))
			mata: add_result_item(`eqn',"``typestack''_`ts'","MSE",          "`m'", `mse')
			mata: add_result_item(`eqn',"``typestack''_`ts'","MSE_folds",    "`m'", st_matrix("`mse_folds'"))
			mata: add_result_item(`eqn',"``typestack''_`ts'","`ts'_weights",   "`m'", st_matrix("`sweights'"))
			// final estimator used to stack is a learner item
			mata: add_learner_item(`eqn',"``typestack''_`ts'","`ts'_final_est", "`finalest'")
			// need base estimators as well
			mata: add_learner_item(`eqn',"``typestack''_`ts'","stack_base_est", "`base_est'")
			// replace updated eqn
			mata: (`mname'.eqnAA).put("`vname'",`eqn')
		}
	}

	tempname sweights0 sweights1 N0_folds N1_folds mse0_folds mse1_folds
	local numvts : word count `vtlist01'
	// loop through 0/1 vtildes
	forvalues i=1/`numvts' {
		local vname :	word `i' of `namelist01'
		local vtilde :	word `i' of `vtlist01'

		mata: `eqn' = (`mname'.eqnAA).get("`vname'")
		mata: st_local("base_est",return_learner_item(`eqn',"`vtilde'","stack_base_est"))
		mata: st_local("stype",return_learner_item(`eqn',"`vtilde'","stack_type"))
		mata: st_local("etype",`eqn'.etype)
		mata: st_local("lieflag", strofreal(`eqn'.lieflag))
		local nlearners : word count `base_est'
		// check if previously stacked; required for poolstacking
		if `ssflag' {
			mata: st_local("shortstack", `eqn'.shortstack)
			if "`shortstack'"=="" {
				// not previously shortstacked, set local and struct field to default
				local shortstack `etype'_`vname'
				mata: `eqn'.shortstack = "`shortstack'"
				local newstack 1
			}
			else {
				local newstack 0
			}
		}
		else {
			mata: st_local("poolstack", `eqn'.poolstack)
			// must previously have been poolstacked for restacking to work
			if "`poolstack'"=="" {
				di as err "error - must poolstack at crossfitting stage first"
				exit 198
			}
		}

		// loop through reps
		forvalues m=1/`reps' {
			local fid `mname'_fid_`m'
			// loop through treatvar=0/1
			forvalues t=0/1 {
			
				// assemble learner list
				// clear macro
				local learner_list
				forvalues j=1/`nlearners' {
					local learner_list `learner_list' `vtilde'`t'_L`j'_`m'
				}
				
				// get stacking weights
				tempvar yhat
				if `ssflag' {
					// shortstacking uses crossfit predictions
					`qui' _ddml_nnls `vname' `learner_list', finalest(`finalest') stype(`stype') if `touse' & `treatvar'==`t'
					local finalest		`e(finalest)'
					mat `sweights`t''	= e(b)
					qui predict double `yhat'
				}
				else {
					// poolstacking uses stacking CV predictions, stored in a mata struct
					tempname tframe y_stacking_cv
					qui frame pwf
					local cframe `r(currentframe)'
					frame create `tframe'
					frame change `tframe'
					tempname
					mata: `y_stacking_cv' = return_result_item(`eqn',"``typestack''_ps","y_stacking_cv`t'", "`m'")				
					// touse ignored at weights stage
					getmata (`fid' `vname' `learner_list')=`y_stacking_cv', force replace
					`qui' _ddml_nnls `vname' `learner_list', finalest(`finalest') stype(`stype')
					mata: `sweights`t'' = st_matrix("e(b)")
					frame change `cframe'
					frame drop `tframe'
					mata: st_matrix("`sweights`t''",`sweights`t'')
					mat colnames `sweights`t'' =  `learner_list'
					mat score double `yhat' = `sweights`t'' if `touse'
					cap mata: mata drop `y_stacking_cv'
					cap mata: mata drop `sweights`t''
				}
				
				if `newstack' {
					cap drop ``typestack''_`ts'`t'_`m'
					qui gen double ``typestack''_`ts'`t'_`m' = `yhat'
					label var ``typestack''_`ts'`t'_`m' "Predicted values cond. exp. of `vname' given `treatvar'=`t' using `typestack'ing"
				}
				else {
					if "`noisily'"~="" {
						di
						di "Existing vs new predicted values:"
						sum ``typestack''_`ts'`t'_`m' `yhat'
					}
					qui replace ``typestack''_`ts'`t'_`m' = `yhat'
				}
				get_stack_stats if `touse', kfolds(`kfolds') fid(`mname'_fid_`m') vname(`vname') vhat(`yhat')
				local N				= r(N)
				local mse			= r(mse)
				mat `N_folds'		= r(N_folds)
				mat `mse_folds'		= r(mse_folds)
				
				// to store:
				mata: add_result_item(`eqn',"``typestack''_`ts'","N`t'",            "`m'", `N')
				mata: add_result_item(`eqn',"``typestack''_`ts'","N`t'_folds",      "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`eqn',"``typestack''_`ts'","MSE`t'",          "`m'", `mse')
				mata: add_result_item(`eqn',"``typestack''_`ts'","MSE`t'_folds",    "`m'", st_matrix("`mse_folds'"))
				mata: add_result_item(`eqn',"``typestack''_`ts'","`ts'_weights`t'",   "`m'", st_matrix("`sweights'"))
				// final estimator used to stack is a learner item
				mata: add_learner_item(`eqn',"``typestack''_`ts'","`ts'_final_est", "`finalest'")
			}
			
			// replace updated eqn
			mata: (`mname'.eqnAA).put("`vname'",`eqn')
		}
	}

	// update flag on mstruct
	if `ssflag'		mata: `mname'.ssflag = 1
	else			mata: `mname'.psflag = 1
	// re-stacking means any previous model estimation results should be dropped
	mata: clear_model_estimation(`mname')

end

// utility for stacking results
program get_stack_stats, rclass
	version 16
	syntax [anything] [if] [in] , [ kfolds(integer 2) fid(varname) vname(varname) vhat(varname) ]
	
	marksample touse
	markout `touse' `fid' `vname' `vhat'
	
	tempname mse_folds N_folds mse_list N_list mse_folds_list N_folds_list
	
	// calculate and return mspe and sample size
	tempvar vres_sq
	// stack macros have fitted values
	qui gen double `vres_sq' = (`vname' - `vhat')^2

	// additive-type model
	qui sum `vres_sq' if `touse', meanonly
	local mse			= r(mean)
	local N				= r(N)
	forvalues k = 1(1)`kfolds' {
		qui sum `vres_sq' if `touse' & `fid'==`k', meanonly
		mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
		qui count if `touse' & `fid'==`k' & `vres_sq'<.
		mat `N_folds' = (nullmat(`N_folds'), r(N))
	}
	
	return scalar mse		= `mse'
	return scalar N			= `N'
	return mat mse_folds	= `mse_folds'
	return mat N_folds		= `N_folds'
	
end



// a single user-specified estimation
program _ddml_estimate_single, eclass sortpreserve
	version 16
	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								y0(varname)			/// for estimating by hand...
								y1(varname)			/// 
								d(varname)			/// 
								d0(varname)			/// 
								d1(varname)			/// 
								z(varname)			///
								ATET 				///
								ATEU				///
								foldid(varname)		/// needed for ATET and ATEU
								ROBust				/// has no effect
								CLUster(varname)	///
								vce(string)			///
								trim(real 0.01)		///
								* ]

	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))

	// checks
	if "`y0'"=="" | "`y1'"=="" {
		di as err "options y0(.) and y1(.) required"
		exit 198
	}
	if "`model'"=="interactive" {
		if "`z'"~="" {
			di as err "z(.) should be empty for model=interactive"
			exit 198
		}
		if "`d'"=="" {
			di as err "option d(.) missing"
			exit 198
		}
		if "`atet'`ateu'"~="" & "`foldvar'"=="" {
			di as err "option foldvar(varname) required for ATET or ATEU"
			exit 198
		}
	}
	else if "`model'"=="late" {
		if "`z'"=="" {
			di as err "option z(.) missing"
			exit 198
		}
		if "`d0'"=="" | "`d1'"=="" {
			di as err "options d0(.) and d1(.) required"
			exit 198
		}
	}
	else {
		di as err "error - unknown model `model'"
		exit 198
	}

	local ateflag=("`model'"=="interactive")
	// teffect is either ATE, ATET, ATEU or LATE
	local check_opt : word count `atet' `ateu'
	if `check_opt' > 1 {
		di as err "error - may request only one of atet or ateu"
		exit 198
	}
	else if `check_opt'==1 & ~`ateflag' {
		di as err "error - `atet'`ateu' not available with LATE (interactiveiv)"
		exit 198
	}
	if ~`ateflag' {
		// it's LATE
		local teffect LATE
	}
	else if "`atet'`ateu'"=="" {
		// default
		local teffect ATE
	}
	else if "`atet'"~="" {
		local teffect ATET
	}
	else {
		local teffect ATEU
	}
	
	if "`robust'"!=""	local vce robust
	local vce1: word 1 of `vce'
	if "`vce1'"=="cluster" {
		local clustvar : word 2 of `vce'
	}
	tempname b V
	tempvar esample
	marksample touse
	if "`model'"=="interactive" {
		qui gen byte `esample' = `y0'<. & `y1'<. & `d'<. & `touse'
		mata: ATE("`teffect'","`nameY'","`nameD'","`y0'", "`y1'", "`d'","`touse'","`b'","`V'","`clustvar'","`foldid'",`trim')
	}
	else {
		qui gen byte `esample' = `y0'<. & `y1'<. & `d0'<. & `d1'<. & `z'<. & `touse'
		mata: LATE("`nameY'","`nameD'","`nameZ'","`y0'", "`y1'", "`d0'","`d1'","`z'","`touse'","`b'","`V'","`clustvar'",`trim')
	}
	if "`clustvar'"=="" {
		// e(.) for basic robust
		local vce		robust
		local vcetype	Robust
	}
	else {
		// e(.) for cluster-robust; clustvar already defined
		local vce		cluster
		local vcetype	Robust
		local N_clust	=r(N_clust)
	}
	local lltrim	
	mat `b'				=r(b)
	mat `V'				=r(V)
	local N				=r(N)
	local lltrim		=r(lltrim)
	local ultrim		=r(ultrim)
	mat colnames `b'	=`nameD'
	mat colnames `V'	=`nameD'
	mat rownames `V'	=`nameD'

	// qui gen byte `esample'=
	ereturn post `b' `V', depname(`nameY') obs(`N') esample(`esample')
	ereturn local cmd		ddml
	ereturn local model		`model'
	ereturn local z			`z'
	ereturn local d			`d'
	ereturn local d0		`d0'
	ereturn local d1		`d1'
	ereturn local y0		`y0'
	ereturn local y1		`y1'
	ereturn local vce		`vce'
	ereturn local vcetype	`vcetype'
	ereturn scalar lltrim	=`lltrim'
	ereturn scalar ultrim	=`ultrim'
	if "`clustvar'"~="" ereturn scalar N_clust=`Nclust'
	
	di
	di as text "E[y|X,D=0]" _col(14) "= " as res "`e(y0)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	di as text "E[y|X,D=1]" _col(14) "= " as res "`e(y1)'"
	if "`e(model)'"=="interactive" {
		di as text "E[D|X]" _col(14)  "= " as res "`e(d)'"
	}
	else {
		di as text "E[D|X,Z=0]" _col(14)  "= " as res "`e(d0)'"
		di as text "E[D|X,Z=1]" _col(14)  "= " as res "`e(d1)'"
		di as text "E[Z|X]" _col(14)  "= " as res "`e(z)'"
	}
	ereturn display
	
	// report warning if clustered SEs requested but doesn't match clustered crossfitting
	mata: st_local("fclustvar",`mname'.fclustvar)
	if "`e(clustvar)'"~="" {
		if "`fclustvar'"=="" {
			di as res "Warning" as text ": crossfit folds do not respect cluster structure used for VCE."
		}
		else if "`fclustvar'"~="`e(clustvar)'" {
			di as res "Warning" as text ": cluster variable for VCE does not match cluster variable for crossfit folds."
		}
	}
	
	// warn if any values trimmed
	if e(lltrim)>0 & e(lltrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(lltrim) " propensity scores trimmed to lower limit " e(trim) "."
	}
	if e(ultrim) & e(ultrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(ultrim) " propensity scores trimmed to upper limit " 1-e(trim) "."
	}
end

// main estimation program
program _ddml_estimate_main
	version 16
	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								ATET 				///
								ATEU				///
								ROBust				/// has no effect - vcv always robust or cluster-robust
								CLUster(varname)	///
								vce(string)			///
								allcombos			/// estimate and show all combinations
								NOTable				/// suppress summary table
								clear				/// deletes all tilde-variables (to be implemented)
								spec(string)		/// specification to post/display
								REP(string)			/// resampling iteration to post/display
								replay				/// model has been estimated, just display results
								trim(real 0.01)		///
								debug				///
								noisily				///
								]
	
	if "`debug'`noisily'"==""	local qui qui
	
	** standard errors
	// local vce is the argument to the Stata option vce(.)
	// SEs are always either robust or cluster-robust
	if "`cluster'"~=""	local vce cluster `cluster'
	else				local vce robust
	
	marksample touse

	// replay existing results
	local replayflag = "`replay'"~=""
	// display summary table
	local tableflag = "`notable'"==""
	// request estimation/reporting of all combinations
	local doallcombos = "`allcombos'"~=""
	// remaining macro flags	
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("stdflag",strofreal(`mname'.stdflag))
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	mata: st_local("psflag",strofreal(`mname'.psflag))
	
	// can't estimate unless crossfitted first
	if `crossfitted'==0 {
		di as err "error - model must be crossfitted before final estimation"
		exit 198
	}
	else if `crossfitted'<`nreps' {
		di as err "error - total reps=`nreps' but only `crossfitted' reps crossfitted"
		exit 198
	}
	
	// reestimation necessary unless replay specified
	if `replayflag' {
		// estimated macro =0/1 indicating estimation results exist
		mata: st_local("estimated", strofreal(`mname'.estimated))
		// initial ncombos; will be 0 if all combos not (yet) estimated
		mata: st_local("ncombos", strofreal(`mname'.ncombos))
		// error checks
		if `estimated'==0 {
			di as err "error - replay specified but model not yet estimated"
			exit 198
		}
		if `ncombos'==0 & "`spec'"~="" & real("`spec'")<. & real("`spec'")>1 {
			di as err "error - spec(`spec') not available; add 'allcombos' to estimate all combinations"
			di as err "add 'replay' to retrieve one of the available estimations stored in memory"
			exit 198
		}
	}
	else {
		// error checks
		if `doallcombos'==0 & "`spec'"~="" & real("`spec'")<. & real("`spec'")>1 {
			di as err "error - spec(`spec') not available; add 'allcombos' to estimate all combinations"
			exit 198
		}
		mata: clear_model_estimation(`mname')
		local estimated = 0
		local ncombos = 0
	}
	
	// number of possible combos; if only one spec, will replace "mse" label
	mata: st_local("poss_combos",strofreal(return_ncombos(`mname')))
	
	if `doallcombos' & (`poss_combos'==1) {
		di as text "only one possible combination of specifications; allcombos option ignored"
		local doallcombos = 0
	}
		
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// locals used below
	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	local allpystackedmulti = 1
	foreach vname in `nameY' `nameD' `nameZ' {
		mata: `eqn' = (`mname'.eqnAA).get("`vname'")
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		local allpystackedmulti = `allpystackedmulti' & `pystackedmulti'
	}

	// reset stdflag (standard stacking) and psflag (pooled stacking)
	// if not all eqns are pystacked multi-learners
	// nb: stdflag may also be 0 if pystacked is used with voting to get short-stacked only
	if `allpystackedmulti'==0 {
		local stdflag	=0
		local psflag	=0
	}

	// numspecflag=1 unless there is no numbered (i.e. not numbers, ss or ps) spec
	// possible with all pystacked and only shortstacking
	local numspecflag	= (`allpystackedmulti'==0) | (`stdflag'==1) | (`psflag'==1)

	// default specs, in order/if available: ss, ps, mse if ncombos>1, otherwise 1 since ncombos=1
	if "`spec'"=="" {
		if `ssflag' {
			local spec "ss"
		}
		else if "`spec'"=="" & `poss_combos'>1 {
			local spec "mse"
		}
		else if `stdflag' {
			local spec st
		}
		else if "`spec'"=="" {
			local spec 1
		}
	}
	
	// allowed forms of spec and rep
	if "`spec'"=="shortstack"	local spec ss
	if "`spec'"=="poolstack"	local spec ps
	if "`spec'"=="minmse"		local spec mse
	if "`rep'"=="mean"			local rep mn
	if "`rep'"=="median"		local rep md

	// if rep not specified, default is rep=1 when nreps==1; md if nreps>1
	if "`rep'"=="" & `nreps'>1 {
		local rep md
	}
	else if "`rep'"=="" & `nreps'==1 {
		local rep 1
	}
	
	// checks
	if "`spec'"~="ss" & "`spec'"~="ps" & "`spec'"~="mse" & "`spec'"~="st" & real("`spec'")==. {
		di as err "error - invalid spec(`spec')"
		exit 198
	}
	if real("`rep'")==. {
		// rep is an integer or mn/md
		if "`rep'"~="mn" & "`rep'"~="md" {
			di as err "error - illegal rep(.) option `rep'"
			exit 198
		}
	}
	else {
		if (real("`rep'")<1) | (real("`rep'")~=int(real("`rep'"))) {
			di as err "error - illegal rep(.) option `rep'"
			exit 198
		}
	}
	// check that rep, if integer, isn't larger than nreps
	if real("`rep'")!=. {
		if `rep'>`nreps' {
			di as err "rep() cannot be larger than `nreps' in current model specification"
			exit 198
		}
	}

	local ateflag=("`model'"=="interactive")
	// teffect is either ATE, ATET, ATEU or LATE
	local check_opt : word count `atet' `ateu'
	if `check_opt' > 1 {
		di as err "error - may request only one of atet or ateu"
		exit 198
	}
	else if `check_opt'==1 & ~`ateflag' {
		di as err "error - `atet'`ateu' not available with LATE (interactiveiv)"
		exit 198
	}
	if ~`ateflag' {
		// it's LATE
		local teffect LATE
	}
	else if "`atet'`ateu'"=="" {
		// default
		local teffect ATE
	}
	else if "`atet'"~="" {
		local teffect ATET
	}
	else {
		local teffect ATEU
	}
	
	// should only ever be 1 D or 1 Z eqn
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	if `numeqnD'>1 {
		di as err "error - model `model' supports only a single treatment variable"
		exit 198
	}
	if `ateflag' & (`numeqnZ'>1) {
		di as err "error - model `model' supports only a single instrument"
		exit 198
	}

	*** shortstack names
	if `ssflag' {
		// code works for both ATE and LATE
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("shortstack", `eqn'.shortstack)
		local Y0ss	`shortstack'_ss
		local Y1ss	`shortstack'_ss
		mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
		mata: st_local("shortstack", `eqn'.shortstack)
		local Dss	`shortstack'_ss
		local D0ss	`shortstack'_ss
		local D1ss	`shortstack'_ss
		if "`model'"=="late" {
			mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
			mata: st_local("shortstack", `eqn'.shortstack)
			local Zss	`shortstack'_ss
		}
	}
	*** poolstack names
	if `psflag' {
		// code works for both ATE and LATE
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("poolstack", `eqn'.poolstack)
		local Y0ps	`poolstack'_ps
		local Y1ps	`poolstack'_ps
		mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
		mata: st_local("poolstack", `eqn'.poolstack)
		local Dps	`poolstack'_ps
		local D0ps	`poolstack'_ps
		local D1ps	`poolstack'_ps
		if "`model'"=="late" {
			mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
			mata: st_local("poolstack", `eqn'.poolstack)
			local Zps	`poolstack'_ps
		}
	}
	
	// multiple specs
	// text locals control messages and lookups
	if `replayflag' {
		local spectext	`spec'
	}
	else if `poss_combos'>1 {
		local msetext	"min-mse "
		local MSEtext	"Min MSE "
		local sttext
		local STtext
		local spectext	mse
	}
	else if `stdflag' {
		local msetext
		local MSEtext
		local sttext	"stacking "
		local STtext	"Stacking "
		local spectext	st
	}
	else {
		local msetext
		local MSEtext
		local sttext
		local STtext
		local spectext	1
	}

	************* ESTIMATE ************
	
	if `estimated'==0 {
		// enter if no estimates exist
		// Loop over resamples and estimate/save the min mse and ss model for each
		forvalues m=1/`nreps' {
		
			// text used in output below
			// multiple reps
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			
			// reset locals
			local Y0opt
			local Y1opt
			local Dopt
			local D0opt
			local D1opt
			local Zopt
			
			*** retrieve best model
			mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
			mata: st_local("Y0opt",return_learner_item(`eqn',"opt0","`m'"))
			mata: st_local("Y1opt",return_learner_item(`eqn',"opt1","`m'"))
			mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
			if `ateflag' {
				mata: st_local("Dopt",return_learner_item(`eqn',"opt","`m'"))
			}
			else {
				mata: st_local("D0opt",return_learner_item(`eqn',"opt0","`m'"))
				mata: st_local("D1opt",return_learner_item(`eqn',"opt1","`m'"))
				mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
				mata: st_local("Zopt",return_learner_item(`eqn',"opt","`m'"))
			}

			// code works for both ATE and LATE
			// estimate optimal (or only numbered) spec for this rep
			if `numspecflag' {
				local title "`STtext'`MSEtext'DDML model`stext' (`teffect')"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')				///
					y0tilde(`Y0opt') y1tilde(`Y1opt')						///
					dtilde(`Dopt') d0tilde(`D0opt') d1tilde(`D1opt')		///
					ztilde(`Zopt')											///
					spec(`spectext') rep(`m')								///
					mname(`mname')											///
					title(`title')											///
					trim(`trim')											///
					vce(`vce')												///
					te(`teffect')
			}
			// estimate shortstack for this rep
			if `ssflag' {
				// code works for both ATE and LATE
				local title "Shortstack DDML model`stext' (`teffect')"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')				///
					y0tilde(`Y0ss') y1tilde(`Y1ss')							///
					dtilde(`Dss') d0tilde(`D0ss') d1tilde(`D1ss')			///
					ztilde(`Zss')											///
					spec(ss) rep(`m'	)									///
					mname(`mname')											///
					title(`title')											///
					trim(`trim')											///
					vce(`vce')												///
					te(`teffect')
			}
			
			// estimate poolstack for this rep
			if `psflag' {
				// code works for both ATE and LATE
				local title "Poolstack DDML model`stext' (`teffect')"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')				///
					y0tilde(`Y0ps') y1tilde(`Y1ps')							///
					dtilde(`Dps') d0tilde(`D0ps') d1tilde(`D1ps')			///
					ztilde(`Zps')											///
					spec(ps) rep(`m')										///
					mname(`mname')											///
					title(`title')											///
					trim(`trim')											///
					vce(`vce')												///
					te(`teffect')
			}

		}

		// have looped over reps to get each optimal model and shortstack per rep
		// now aggregate over reps to get mean/median
		if `nreps' > 1 {
			// numbered/stacking estimates
			if `numspecflag' {
				local title "Mean over `nreps' `sttext'`msetext'resamples (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(`spectext') medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "Median over `nreps' `sttext'`msetext'resamples (`teffect')"
 				`qui' medmean_and_store, mname(`mname') spec(`spectext') medmean(md) title(`title') vce(`vce') te(`teffect')
 			}
			// shortstack
			if `ssflag' {
				local title "Shortstack DDML model (mean over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ss) medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "Shortstack DDML model (median over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ss) medmean(md) title(`title') vce(`vce') te(`teffect')
			}
			// poolstack
			if `psflag' {
				local title "Poolstack DDML model (mean over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ps) medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "Poolstack DDML model (median over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ps) medmean(md) title(`title') vce(`vce') te(`teffect')
			}
		}
		
		// Estimation of med/med/shortstack complete, mark as estimated.
		mata: `mname'.estimated = 1
		// (re-)set estimated
		mata: st_local("estimated", strofreal(`mname'.estimated))
	}
	
	******************************************
	
	if `ncombos' {
		// all combos have already been estimated, so just recover matrices
		tempname nmat bmat semat
		mata: `nmat' = (`mname'.estAA).get(("nmat","all"))
		mata: `bmat' = (`mname'.estAA).get(("bmat","all"))
		mata: `semat' = (`mname'.estAA).get(("semat","all"))
		// recover min MSE specs
		forvalues m=1/`nreps' {
			mata: check_spec(`mname',"optspec","`m'")
			mata: st_local("optspec",(`mname'.estAA).get(("optspec","`m'")))
			local optspec`m' = `optspec'
		}
	}
	else if `doallcombos' {
		// allcombos need to be estimated
		// ATE and LATE both use Y0/Y1
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
		mata: st_local("Ytilde",invtokens(`eqn'.vtlist))
		// ATE and LATE both have a D eqn
		mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameD)
		mata: st_local("Dtilde",invtokens(`eqn'.vtlist))
		if `ateflag' {
			// ATE has a single Dtilde
			_ddml_allcombos `Ytilde'- `Ytilde' - `Dtilde',		///
				addprefix("")									///
				`debug'  
			local Dlist `r(colstr3)'
		}
		else {
			// LATE has a D0/D1 and a single Z
			mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameZ)
			mata: st_local("Ztilde",invtokens(`eqn'.vtlist))
			_ddml_allcombos `Ytilde' - `Ytilde' - `Dtilde' - `Dtilde' - `Ztilde' ,	///
				`debug'																///
				addprefix("")
			local d0list `r(colstr3)'
			local d1list `r(colstr4)'
			local Zlist `r(colstr5)' 
		}
		local y0list `r(colstr1)'
		local y1list `r(colstr2)'
		local ncombos = r(ncombos)
		local tokenlen = `ncombos'*2 -1
			
		tempname nmat bmat semat
		mata: `nmat' = J(`ncombos',5,"")
		mata: `bmat' = J(`ncombos'*`nreps',`numeqnD',.)
		mata: `semat' = J(`ncombos'*`nreps',`numeqnD',.)
		
		// simplest if put into a Mata string matrix
		tokenize `y0list' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',1] = strtrim("``idx''")
		}
		tokenize `y1list' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',2] = strtrim("``idx''")
		}
		if `ateflag' {
			// ATE has a single D
			tokenize `Dlist' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',3] = strtrim("``idx''")
			}
		}
		else {
			// LATE has D0/D1 and a single Z
			tokenize `d0list' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',3] = strtrim("``idx''")
			}
			tokenize `d1list' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',4] = strtrim("``idx''")
			}
			tokenize `Zlist' , parse("-")
			forvalues i=1/`ncombos' {
				local idx = 2*`i'-1
				mata: `nmat'[`i',5] = strtrim("``idx''")
			}
		}
		
		forvalues m=1/`nreps' {
			
			// reset locals
			local Y0opt
			local Y1opt
			local Dopt
			local D0opt
			local D1opt
			local Zopt
			
			*** retrieve best model
			mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
			mata: st_local("Y0opt",return_learner_item(`eqn',"opt0","`m'"))
			mata: st_local("Y1opt",return_learner_item(`eqn',"opt1","`m'"))
			mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
			if `ateflag' {
				mata: st_local("Dopt",return_learner_item(`eqn',"opt","`m'"))
			}
			else {
				mata: st_local("D0opt",return_learner_item(`eqn',"opt0","`m'"))
				mata: st_local("D1opt",return_learner_item(`eqn',"opt1","`m'"))
				mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
				mata: st_local("Zopt",return_learner_item(`eqn',"opt","`m'"))
			}
			
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			
			forvalues i = 1/`ncombos' {
				mata: st_local("y0",`nmat'[`i',1])
				mata: st_local("y1",`nmat'[`i',2])
				if `ateflag' {
					mata: st_local("d",`nmat'[`i',3])
				}
				else {
					mata: st_local("d0",`nmat'[`i',3])
					mata: st_local("d1",`nmat'[`i',4])
					mata: st_local("z",`nmat'[`i',5])
				}
				// check if opt for this resample
				// code works for both ATE and LATE
				local isopt
				local isY0opt	: list Y0opt == y0
				local isY1opt	: list Y1opt == y1
				local isDopt	: list Dopt == d
				local isD0opt	: list D0opt == d0
				local isD1opt	: list D1opt == d1
				local isZopt	: list Zopt == z
				local title "DDML model, specification `i'`stext' (`teffect')"
				if `isY0opt' & `isY1opt' & `isDopt' & `isD0opt' & `isD1opt' & `isZopt' {
					local optspec`m' = `i'
					local isopt *
					local title `MSEtext'`title'
					// save in AA
					mata: (`mname'.estAA).put(("optspec","`m'"),"`i'")
				}
				// code works for both ATE and LATE
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')				///
					y0tilde(`y0') y1tilde(`y1')								///
					dtilde(`d') d0tilde(`d0') d1tilde(`d1')					///
					ztilde(`z')												///
					spec(`i') rep(`m')										///
					mname(`mname')											///
					title(`title')											///
					trim(`trim')											///
					vce(`vce')												///
					te(`teffect')
				
				mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(bmat)")
				mata: `semat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(semat)")
				
			}
			
			if `ssflag' {
				// code works for both ATE and LATE
				local title "Shortstack DDML model`stext' (`teffect')"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')				///
					y0tilde(`Y0ss') y1tilde(`Y1ss')							///
					dtilde(`Dss') d0tilde(`D0ss') d1tilde(`D1ss')			///
					ztilde(`Zss')											///
					spec(ss) rep(`m')										///
					mname(`mname')											///
					title(`title')											///
					trim(`trim')											///
					vce(`vce')												///
					te(`teffect')
			}
			
			if `psflag' {
				// code works for both ATE and LATE
				local title "Poolstack DDML model`stext' (`teffect')"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
					yvar(`nameY') dvar(`nameD') zvar(`nameZ')				///
					y0tilde(`Y0ps') y1tilde(`Y1ps')							///
					dtilde(`Dps') d0tilde(`D0ps') d1tilde(`D1ps')			///
					ztilde(`Zps')											///
					spec(ps) rep(`m')										///
					mname(`mname')											///
					title(`title')											///
					trim(`trim')											///
					vce(`vce')												///
					te(`teffect')
			}

		}

		// we make a copy of the MSE-optimal model for each m
		forvalues m=1/`nreps' {
			tempname Bopt
			mata: st_local("optspec",(`mname'.estAA).get(("optspec","`m'")))
			mata: `Bopt' = (`mname'.estAA).get(("`optspec'","`m'"))
			mata: (`mname'.estAA).put(("mse","`m'"),`Bopt')
			mata: mata drop `Bopt'
		}
		
		// aggregate across resamplings
		if `nreps' > 1 {
			local title "Mean over `nreps' min-mse specifications (`teffect')"
 			`qui' medmean_and_store, mname(`mname') spec(mse) medmean(mn) title(`title') vce(`vce') te(`teffect')
 			local title "Median over `nreps' min-mse specifications (`teffect')"
 			`qui' medmean_and_store, mname(`mname') spec(mse) medmean(md) title(`title') vce(`vce') te(`teffect')
			// numbered specifications
			forvalues i = 1/`ncombos' {
				local title "DDML model, specification `i' (mean over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(`i') medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "DDML model, specification `i' (median over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(`i') medmean(md) title(`title') vce(`vce') te(`teffect')
			}
			// shortstack
			if `ssflag' {
				local title "Shortstack DDML model (mean over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ss) medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "Shortstack DDML model (median over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ss) medmean(md) title(`title') vce(`vce') te(`teffect')
			}
			// poolstack
			if `psflag' {
				local title "Poolstack DDML model (mean over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ps) medmean(mn) title(`title') vce(`vce') te(`teffect')
				local title "Poolstack DDML model (median over `nreps' resamples) (`teffect')"
				`qui' medmean_and_store, mname(`mname') spec(ps) medmean(md) title(`title') vce(`vce') te(`teffect')
			}
		}
		
		// estimation of all combos complete; ncombos > 0 indicates all combos estimated
		mata: `mname'.ncombos = `ncombos'
		mata: (`mname'.estAA).put(("nmat","all"),`nmat')
		mata: (`mname'.estAA).put(("bmat","all"),`bmat')
		mata: (`mname'.estAA).put(("semat","all"),`semat')
		
	}
	
	************** REPORT RESULTS **************

	
	*** Results ***
	
	// optional table of all results
	if `tableflag' {
		di
		ddml describe, mname(`mname')
		di
		di as text "DDML estimation results (`teffect'):"
		di as text "spec  r" %14s "Y(0) learner" _c
		di as text           %14s "Y(1) learner" _c
		if `ateflag' {
			di as text           %14s "D learner" _c
		}
		else {
			di as text           %14s "D(0) learner" _c
			di as text           %14s "D(1) learner" _c
		}
		// space after SE for the (.) enclosing the SE
		di as text %10s "b" %11s "SE " _c
		if ~`ateflag' {
			di as text           %14s "Z learner" _c
		}
		di
		forvalues m=1/`nreps' {
			if `doallcombos' {
				// all combos available, so loop through
				forvalues i=1/`ncombos' {
					mata: st_local("yt0",abbrev(`nmat'[`i',1],13))
					mata: st_local("yt1",abbrev(`nmat'[`i',2],13))
					if "`optspec`m''"=="`i'" & `ncombos'>1 {
						// print * for opt (lowest MSE combo) if more than 1 combo
						di "*" _c
					}
					else {
						di " " _c
					}
					local specrep "`: di %3.0f `i' %3.0f `m''"
					local rcmd stata ddml estimate, mname(`mname') spec(`i') rep(`m') notable replay
					di %6s "{`rcmd':`specrep'}" _c
					di as res %14s abbrev("`yt0'",13) _c
					di as res %14s abbrev("`yt1'",13) _c
					if `ateflag' {
						mata: st_local("dt",`nmat'[`i',3])
						di as res %14s abbrev("`dt'",13) _c
					}
					else {
						mata: st_local("dt0",abbrev(`nmat'[`i',3],13))
						mata: st_local("dt1",abbrev(`nmat'[`i',4],13))
						di as res %14s abbrev("`dt0'",13) _c
						di as res %14s abbrev("`dt1'",13) _c
					}
					mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
					mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
					// precede with a space
					di as res " " %9.3f `b' _c
					local se = strtrim("`: di %8.3f `se''")
					local pse (`se')
					// precede with a space
					di as res " " %10s "`pse'" _c
					if ~`ateflag' {
						mata: st_local("zt",abbrev(`nmat'[`i',5],13))
						di as res %14s "`zt'" _c
					}
					di
				}
			}
			else {
				// only mse/ss specs available/reported
				if `numspecflag' {
					`qui' replay_estimate, mname(`mname') spec(`spectext') rep(`m')
					tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
					mat `btemp' = e(b)
					mat `Vtemp' = e(V)
					local specrep "`: di %4s "`spectext'" %3.0f `m''"
					local rcmd stata ddml estimate, mname(`mname') spec(`spectext') rep(`m') notable replay
					di %6s "{`rcmd':`specrep'}" _c
					mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
					mata: st_local("yt0",return_learner_item(`eqn',"opt0","`m'"))
					di as res %14s abbrev("`yt0'",13) _c
					mata: st_local("yt1",return_learner_item(`eqn',"opt1","`m'"))
					di as res %14s abbrev("`yt1'",13) _c
					if `ateflag' {
						mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
						mata: st_local("dt",return_learner_item(`eqn',"opt","`m'"))
						di as res %14s abbrev("`dt'",13) _c
					}
					else {
						mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
						mata: st_local("dt0",return_learner_item(`eqn',"opt0","`m'"))
						di as res %14s abbrev("`dt0'",13) _c
						mata: st_local("dt1",return_learner_item(`eqn',"opt1","`m'"))
						di as res %14s abbrev("`dt1'",13) _c
					
					}
					// precede with a space
					di as res " " %9.3f el(`btemp',1,1) _c
					local se = sqrt(el(`Vtemp',1,1))
					local se = strtrim("`: di %8.3f `se''")
					local pse (`se')
					// precede with a space
					di as res " " %10s "`pse'" _c
					if ~`ateflag' {
						mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
						mata: st_local("zt",return_learner_item(`eqn',"opt","`m'"))
						di as res %14s abbrev("`zt'",13) _c
					}
					di
				}
			}

			if `ssflag' {
				`qui' replay_estimate, mname(`mname') spec(ss) rep(`m')
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ss" %3.0f `m''"
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`m') notable replay
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[shortstack]" _c
				di as res %14s "[ss]" _c
				di as res %14s "[ss]" _c
				if ~`ateflag' {
					di as res %14s "[ss]" _c
				}
				// precede with a space
				di as res " " %9.3f el(`btemp',1,1) _c
				local se = sqrt(el(`Vtemp',1,1))
				local se = strtrim("`: di %8.3f `se''")
				local pse (`se')
				// precede with a space
				di as res " " %10s "`pse'" _c
				if ~`ateflag' {
					di as res %14s "[ss]" _c
				}
				di
			}

			if `psflag' {
				`qui' replay_estimate, mname(`mname') spec(ps) rep(`m')
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ps" %3.0f `m''"
				local rcmd stata ddml estimate, mname(`mname') spec(ps) rep(`m') notable replay
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[poolstack]" _c
				di as res %14s "[ps]" _c
				di as res %14s "[ps]" _c
				if ~`ateflag' {
					di as res %14s "[ps]" _c
				}
				// precede with a space
				di as res " " %9.3f el(`btemp',1,1) _c
				local se = sqrt(el(`Vtemp',1,1))
				local se = strtrim("`: di %8.3f `se''")
				local pse (`se')
				// precede with a space
				di as res " " %10s "`pse'" _c
				if ~`ateflag' {
					di as res %14s "[ps]" _c
				}
				di
			}
		}
		// footnote needed only if multiple specs possible
		if `poss_combos'>1 {
			if `doallcombos' {
				di as res "*" _c
			}
			else {
				di as text "mse" _c
			}
			di as text " = minimum MSE specification for that resample."
		}
	}
		
	if `nreps' > 1 & `tableflag' {
		di
		di as text "Mean/med Y(0) learner" _c
		di as text               %14s "Y(1) learner" _c
		if `ateflag' {
			di as text           %14s "D learner" _c
		}
		else {
			di as text           %14s "D(0) learner" _c
			di as text           %14s "D(1) learner" _c
		}
		// space after SE for the (.) enclosing the SE
		di as text               %10s "b" %11s "SE " _c
		if ~`ateflag' {
			di as text           %14s "Z learner" _c
		}
		di
		foreach medmean in mn md {
			** mean and median over mse
			if `numspecflag' {
				`qui' replay_estimate, mname(`mname') spec(`spectext') rep(`medmean')
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "mse" %3s "`medmean'"'"
				local rcmd stata ddml estimate, mname(`mname') spec(`spectext') rep(`medmean') notable replay
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[min-mse]" _c
				di as res %14s "[mse]" _c
				di as res %14s "[mse]" _c
				if ~`ateflag' {
					di as res %14s "[mse]" _c
				}
				// precede with a space
				di as res " " %9.3f el(`btemp',1,1) _c
				local se = sqrt(el(`Vtemp',1,1))
				local se = strtrim("`: di %8.3f `se''")
				local pse (`se')
				// precede with a space
				di as res " " %10s "`pse'" _c
				if ~`ateflag' {
					di as res %14s "[mse]" _c
				}
				di
			}
			if `ssflag' {
				`qui' replay_estimate, mname(`mname') spec(ss) rep(`medmean')
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ss" %3s "`medmean'"'"
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`medmean') notable replay
				di as res %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[shortstack]" _c
				di as res %14s "[ss]" _c
				di as res %14s "[ss]" _c
				if ~`ateflag' {
					di as res %14s "[ss]" _c
				}
				// precede with a space
				di as res " " %9.3f el(`btemp',1,1) _c
				local se = sqrt(el(`Vtemp',1,1))
				local se = strtrim("`: di %8.3f `se''")
				local pse (`se')
				// precede with a space
				di as res " " %10s "`pse'" _c
				if ~`ateflag' {
					di as res %14s "[ss]" _c
				}
				di
			}
			if `psflag' {
				`qui' replay_estimate, mname(`mname') spec(ps) rep(`medmean')
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ps" %3s "`medmean'"'"
				local rcmd stata ddml estimate, mname(`mname') spec(ps) rep(`medmean') notable replay
				di as res %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[poolstack]" _c
				di as res %14s "[ps]" _c
				di as res %14s "[ps]" _c
				if ~`ateflag' {
					di as res %14s "[ps]" _c
				}
				// precede with a space
				di as res " " %9.3f el(`btemp',1,1) _c
				local se = sqrt(el(`Vtemp',1,1))
				local se = strtrim("`: di %8.3f `se''")
				local pse (`se')
				// precede with a space
				di as res " " %10s "`pse'" _c
				if ~`ateflag' {
					di as res %14s "[ps]" _c
				}
				di
			}
		}
	}

	// select result to display and post
	// default is, as available: 1. ss 2. ps 3. only spec or MSE
	if `replayflag' | "`spec'"~="" {
		local specdisp `spec'
	}
	else if `ssflag' {
		local specdisp ss
	}
	else if `psflag' {
		local specdisp ps
	}
	else if `stdflag' {
		local specdisp st
	}
	else if `poss_combos'==1 {
		local specdisp 1
	}
	else {
		local specdisp mse
	}
	
	// post selected estimates; rep is the resample number (default=1)
	di
	replay_estimate, mname(`mname') spec(`specdisp') rep(`rep')
	di
	
	if `nreps' > 1 & ("`rep'"=="mn" | "`rep'"=="md") {
		tempvar bhat
		svmat e(b_resamples), names(`bhat')
		// variables in Stata will look like _000000A1, _000000A2, etc. and will disappear as temps after exit
		local dnames : colnames e(b)
		
		di as text "Summary over " `nreps' " resamples:"
		di as text %12s "D eqn" %10s "mean" %10s "min" %10s "p25" %10s "p50" %10s "p75" %10s "max"
		local i 1
		foreach vn in `dnames' {
			di as text %12s abbrev("`vn'",11) _col(15) _c
			qui sum `bhat'`i', detail
			// intersperse with spaces
			di as res " " %9.4f r(mean) _c
			di as res " " %9.4f r(min) _c
			di as res " " %9.4f r(p25) _c
			di as res " " %9.4f r(p50) _c
			di as res " " %9.4f r(p75) _c
			di as res " " %9.4f r(max)
			local ++i
		}
	}
	
	// temp Mata objects no longer needed
	foreach obj in `eqn' `nmat' `bmat' `semat' {
		cap mata: mata drop `obj'
	}

end

// adds rep number suffixes to list of varnames
program define add_suffix, sclass
	version 16
	syntax [anything] , suffix(name)

	// anything is a list of to-be-varnames that need suffix added to them
	foreach vn in `anything' {
		local vnames `vnames' `vn'`suffix'
	}
	
	sreturn local vnames `vnames'
end


********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

void ATE(   
			string scalar teffect,    // ATE, ATET or ATEU
			string scalar yvar,       // Y
			string scalar dvar,       // D
			string scalar y0tilde,    // E[Y|X,D=0]
			string scalar y1tilde,    // E[Y|X,D=1]
			string scalar dtilde,     // E[D|X]
			string scalar sample,     // sample
			string scalar outate,     // output: name of matrix to store b
			string scalar outatese,   // output: name of matrix to store V
			string scalar clustvar,   //
			string scalar foldvar, 	  //
			real scalar trim          // trim the propensity score
			)
{
	st_view(my_d0x,.,y0tilde,sample)
	st_view(my_d1x,.,y1tilde,sample)
	st_view(d,.,dvar,sample)
	st_view(y,.,yvar,sample)
	st_view(fid,.,foldvar,sample)
	// copy since we may trim it
	md_x=st_data(.,dtilde,sample)
	if (clustvar!="") {
		st_view(clustid,.,clustvar,sample)
		clustid_uni=uniqrows(clustid)
		nclust = rows(clustid_uni)
	}

	n = rows(y)
	
	// first a vector
	lltrim = md_x :< trim
	ultrim = md_x :> (1-trim)
	// trim
	md_x = md_x :* (1:-lltrim) + trim * lltrim
	md_x = md_x :* (1:-ultrim) + (1-trim) * ultrim
	// now a scalar
	lltrim = sum(lltrim)
	ultrim = sum(ultrim)

	// psi = psi_b + psi_a*theta, e.g. equation 5.62
	if (teffect=="ATE") {
		psi_b  = (d :* (y :- my_d1x) :/ md_x) :-  ((1 :- d) :* (y :- my_d0x) :/ (1 :- md_x)) :+ my_d1x :- my_d0x 
		psi_a  = J(n,1,-1) 
	}
	else {
		// ATET and ATEU
		// calculate mean of d by fold
		fid_uni = uniqrows(fid)
		folds = rows(fid_uni)
		p_hat = J(n,1,.)
		for (j=1;j<=folds;j++) {
			k=fid_uni[j,1]
			sel = selectindex(fid:==k)
			meank = mean(d[sel])
			p_hat[sel] = J(length(sel), 1, meank)
		}
		if (teffect=="ATET") {
			psi_b = (d :* (y :- my_d0x) :/ p_hat) :-  md_x :* (1 :- d) :* (y :- my_d0x) :/ (p_hat :*(1 :- md_x)) 
			psi_a = -d :/ p_hat
		}
		else if (teffect=="ATEU") {
			psi_b = ((1:-d) :* (y :- my_d1x) :/ (1:-p_hat)) :- (1:-md_x) :* d :* (y :- my_d1x) :/ ((1:-p_hat) :* md_x) 
			psi_a = (1:-d) :/ (1:-p_hat)
		}
		else {
			errprintf("internal ddml error - teffect estimator %s not defined\n",teffect)
			exit(499)
		}
	}
	theta = -mean(psi_b) / mean(psi_a)
	psi = psi_a :* theta :+ psi_b

	if (clustvar=="") {
		V =  mean(psi:^2) / (mean(psi_a):^2) / n
	}
	else {
		gamma = 0
		jhat = 0
		for (i=1;i<=nclust;i++) {
			psi_c = select(psi,clustid:==clustid_uni[i,1])
			psi_a_c = select(psi_a,clustid:==clustid_uni[i,1])
			gamma = gamma :+ 1/nclust :* sum(psi_c*psi_c')
			jhat = jhat :+  1/nclust :* sum(psi_a_c)
		}
		V = gamma / jhat:^2 / nclust
		st_numscalar("r(N_clust)",nclust)
	}

	st_numscalar("r(N)",n)
	st_matrix("r(b)",theta)
	st_matrix("r(V)",V)
	st_numscalar("r(lltrim)",lltrim)
	st_numscalar("r(ultrim)",ultrim)
}

void LATE(  string scalar yvar,      // Y
            string scalar dvar,      // D
            string scalar zvar,      // Z
            string scalar y0tilde,   // E[Y|X,Z=0]
            string scalar y1tilde,   // E[Y|X,Z=1]
            string scalar d0tilde,   // E[D|X,Z=0]
            string scalar d1tilde,   // E[D|X,Z=1]
            string scalar ztilde,    // E[Z|X]
            string scalar sample,    // sample
            string scalar outlate,   // output: name of matrix to store b
            string scalar outlatese,  // output: name of matrix to store V
            string scalar clustvar,
 			real scalar trim          // trim the propensity score
           )
{
    st_view(my_z0x,.,y0tilde,sample)
    st_view(my_z1x,.,y1tilde,sample)
    st_view(md_z0x,.,d0tilde,sample)
    st_view(md_z1x,.,d1tilde,sample)
    st_view(d,.,dvar,sample)
    st_view(y,.,yvar,sample)
    st_view(z,.,zvar,sample)
	// copy since we may trim it
    mz_x=st_data(.,ztilde,sample)
    if (clustvar!="") {
		st_view(clustid,.,clustvar,sample)
		clustid_uni=uniqrows(clustid)
		nclust = rows(clustid_uni)
	}

    n = rows(y)

	// first a vector
	lltrim = mz_x :< trim
	ultrim = mz_x :> (1-trim)
	// trim
	mz_x = mz_x :* (1:-lltrim) + trim * lltrim
	mz_x = mz_x :* (1:-ultrim) + (1-trim) * ultrim
	// now a scalar
	lltrim = sum(lltrim)
	ultrim = sum(ultrim)
	
    psi_b =  z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x 
    psi_a =  -(z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x)

	theta = -mean(psi_b) / mean(psi_a)
	psi = psi_a :* theta :+ psi_b

	if (clustvar=="") {
		V =  mean(psi:^2) :/ mean(psi_a):^2 :/ n
	}
	else {
		gamma = 0
		jhat = 0
		for (i=1;i<=nclust;i++) {
			psi_c = select(psi,clustid:==clustid_uni[i,1])
			psi_a_c = select(psi_a,clustid:==clustid_uni[i,1])
			gamma = gamma :+ 1/nclust :* sum(psi_c*psi_c')
			jhat = jhat :+  1/nclust :* sum(psi_a_c)
		}
		V = gamma / jhat:^2 / nclust
		st_numscalar("r(N_clust)",nclust)
	}

    st_numscalar("r(N)",n)
    st_matrix("r(b)",theta)
    st_matrix("r(V)",V)
	st_numscalar("r(lltrim)",lltrim)
	st_numscalar("r(ultrim)",ultrim)
}

end

* code below currently supports only a single treatment variable, but coded for multiple variables in places
program define estimate_and_store, eclass

	version 16

	syntax [anything] [if] [in] ,				///
						[						///
							yvar(varname)		///
							dvar(varname)		///
							zvar(varname)		///
							y0tilde(name)		///
							y1tilde(name)		///
							dtilde(name)		///
							d0tilde(name)		///
							d1tilde(name)		///
							ztilde(name)		///
							mname(name)			///
							spec(string)		///
							rep(string)			///
							replay				///
							title(string)		///
							medmean(string)		///
							vce(string) 		///
							TEffect(name)		/// either ATE, ATET or ATEU
							trim(real 0)	    /// value should be provided by calling program
						]
		
	mata: st_local("model",`mname'.model)
	
	marksample touse
		
	tempname A
	mata: `A' = AssociativeArray()
	mata: `A'.reinit("string",2)
	mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
	
	// 0/1 etc
	local y0		`y0tilde'0
	local y1		`y1tilde'1
	local d			`dtilde'
	local d0		`d0tilde'0
	local d1		`d1tilde'1
	local z			`ztilde'
	// add suffixes and 0/1 indicator
	local y0_m		`y0tilde'0_`rep'
	local y1_m		`y1tilde'1_`rep'
	local d_m		`dtilde'_`rep'
	local d0_m		`d0tilde'0_`rep'
	local d1_m		`d1tilde'1_`rep'
	local z_m		`ztilde'_`rep'

	// estimation samples may differ across conditional expectations
	if "`model'"=="interactive" {
		markout `touse' `y0_m' `y1_m' `d_m'
	}
	else {
		markout `touse' `y0_m' `y1_m' `d0_m' `d1_m' `z_m'
	}
	
	local vce1: word 1 of `vce'
	if "`vce1'"=="cluster" {
		local clustvar : word 2 of `vce'
	}
	if "`model'"=="interactive" {
		mata: ATE("`teffect'","`yvar'","`dvar'","`y0_m'", "`y1_m'", "`d_m'","`touse'","`b'","`V'","`clustvar'","`mname'_fid_`rep'",`trim')
	}
	else {
		mata: LATE("`yvar'","`dvar'","`zvar'","`y0_m'", "`y1_m'", "`d0_m'","`d1_m'","`z_m'","`touse'","`b'","`V'","`clustvar'",`trim')
	}
	if "`clustvar'"=="" {
		// e(.) for basic robust
		local vce		robust
		local vcetype	Robust
	}
	else {
		// e(.) for cluster-robust; clustvar already defined
		local vce		cluster
		local vcetype	Robust
		local N_clust	=r(N_clust)
	}
	
	// store post objects
	mata: `A'.put(("N","post"),`r(N)')
	mata: `A'.put(("b","post"),st_matrix("r(b)"))
	mata: `A'.put(("V","post"),st_matrix("r(V)"))
	mata: `A'.put(("depvar","post"),"`yvar'")
	
	// for calling program
	ereturn clear
	mata: st_matrix("e(bmat)",st_matrix("r(b)"))
	mata: st_matrix("e(semat)",sqrt(diagonal(st_matrix("r(V)"))'))
	
	// store locals
	local list_local title yvar dvar y0 y0_m y1 y1_m vce vcetype teffect
	if "`model'"=="interactive" {
		local list_local `list_local' d d_m
	}
	else {
		local list_local `list_local' d0 d0_m d1 d1_m z z_m
	}
	if "`clustvar'"~=""		local list_local `list_local' clustvar
	foreach obj in `list_local' {
		mata: `A'.put(("`obj'","local"),"``obj''")
	}
	// store scalars
	mata: `A'.put(("lltrim","scalar"),`r(lltrim)')
	mata: `A'.put(("ultrim","scalar"),`r(ultrim)')
	mata: `A'.put(("trim","scalar"),`trim')
	if "`clustvar'"~="" {
		mata: `A'.put(("N_clust","scalar"),`N_clust')
	}
	// additional estimation results
	tempname eqn
	mata: `eqn' = init_eStruct()
	// Y eqn results
	mata: `eqn' = (`mname'.eqnAA).get("`yvar'")
	mata: st_local("shortstack",`eqn'.shortstack)
	mata: st_local("poolstack",`eqn'.poolstack)
	// MSE
	mata: `A'.put(("`y0'_mse","scalar"),return_result_item(`eqn',"`y0tilde'","MSE0","`rep'"))
	mata: `A'.put(("`y1'_mse","scalar"),return_result_item(`eqn',"`y1tilde'","MSE1","`rep'"))
	// MSE folds
	mata: `A'.put(("`y0'_mse_folds","matrix"),return_result_item(`eqn',"`y0tilde'","MSE0_folds","`rep'"))
	mata: `A'.put(("`y1'_mse_folds","matrix"),return_result_item(`eqn',"`y1tilde'","MSE1_folds","`rep'"))
	// pystacked final est (pystacked multi only)
	mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
	if `pystackedmulti' {
	// cap because won't exist for e.g. shortstack variable
		cap mata: `A'.put(("`y0'_stack_final_est","local"), return_learner_item(`eqn',"`y0tilde'","stack_final_est"))
		cap mata: `A'.put(("`y1'_stack_final_est","local"), return_learner_item(`eqn',"`y1tilde'","stack_final_est"))
	}
	// shortstack results
	if "`spec'"=="ss" {
		mata: `A'.put(("`yvar'_ssw0","matrix"),return_result_item(`eqn',"`shortstack'_ss","ss_weights0","`rep'"))
		mata: `A'.put(("`yvar'_ssw1","matrix"),return_result_item(`eqn',"`shortstack'_ss","ss_weights1","`rep'"))
		mata: `A'.put(("`yvar'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
	}
	// poolstack results
	if "`spec'"=="ps" {
		mata: `A'.put(("`yvar'_psw0","matrix"),return_result_item(`eqn',"`poolstack'_ps","ps_weights0","`rep'"))
		mata: `A'.put(("`yvar'_psw1","matrix"),return_result_item(`eqn',"`poolstack'_ps","ps_weights1","`rep'"))
		mata: `A'.put(("`yvar'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
	}
	if "`model'"=="interactive" {
		// D eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
		mata: st_local("shortstack",`eqn'.shortstack)
		mata: st_local("poolstack",`eqn'.poolstack)
		// MSE
		mata: `A'.put(("`dtilde'_mse","scalar"),return_result_item(`eqn',"`dtilde'","MSE","`rep'"))
		// MSE folds
		mata: `A'.put(("`dtilde'_mse_folds","matrix"),return_result_item(`eqn',"`dtilde'","MSE_folds","`rep'"))
		// shortstack results
		if "`spec'"=="ss" {
			mata: `A'.put(("`dvar'_ssw","matrix"),return_result_item(`eqn',"`shortstack'_ss","ss_weights","`rep'"))
			mata: `A'.put(("`dvar'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
		}
		// poolstack results
		if "`spec'"=="ps" {
			mata: `A'.put(("`dvar'_psw","matrix"),return_result_item(`eqn',"`poolstack'_ps","ps_weights","`rep'"))
			mata: `A'.put(("`dvar'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
		}
	}
	else {
		// D
		mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
		mata: st_local("shortstack",`eqn'.shortstack)
		mata: st_local("poolstack",`eqn'.poolstack)
		// MSE, D
		mata: `A'.put(("`d0'_mse","scalar"),return_result_item(`eqn',"`d0tilde'","MSE0","`rep'"))
		mata: `A'.put(("`d1'_mse","scalar"),return_result_item(`eqn',"`d1tilde'","MSE1","`rep'"))
		// MSE folds, D
		mata: `A'.put(("`d0'_mse_folds","matrix"),return_result_item(`eqn',"`d0tilde'","MSE0_folds","`rep'"))
		mata: `A'.put(("`d1'_mse_folds","matrix"),return_result_item(`eqn',"`d1tilde'","MSE1_folds","`rep'"))
		// shortstack results, D
		if "`spec'"=="ss" {
			mata: `A'.put(("`dvar'_ssw0","matrix"),return_result_item(`eqn',"`shortstack'_ss","ss_weights0","`rep'"))
			mata: `A'.put(("`dvar'_ssw1","matrix"),return_result_item(`eqn',"`shortstack'_ss","ss_weights1","`rep'"))
			mata: `A'.put(("`dvar'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
		}
		// poolstack results, D
		if "`spec'"=="ps" {
			mata: `A'.put(("`dvar'_psw0","matrix"),return_result_item(`eqn',"`poolstack'_ps","ps_weights0","`rep'"))
			mata: `A'.put(("`dvar'_psw1","matrix"),return_result_item(`eqn',"`poolstack'_ps","ps_weights1","`rep'"))
			mata: `A'.put(("`dvar'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
		}
		// Z
		mata: `eqn' = (`mname'.eqnAA).get("`zvar'")
		mata: st_local("shortstack",`eqn'.shortstack)
		mata: st_local("poolstack",`eqn'.poolstack)
		// MSE, Z
		mata: `A'.put(("`ztilde'_mse","scalar"),return_result_item(`eqn',"`ztilde'","MSE","`rep'"))
		// MSE folds, Z
		mata: `A'.put(("`ztilde'_mse_folds","matrix"),return_result_item(`eqn',"`ztilde'","MSE_folds","`rep'"))
		// shortstack results, Z
		if "`spec'"=="ss" {
			mata: `A'.put(("`zvar'_ssw","matrix"),return_result_item(`eqn',"`shortstack'_ss","ss_weights","`rep'"))
			mata: `A'.put(("`zvar'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
		}
		// poolstack results, Z
		if "`spec'"=="ps" {
			mata: `A'.put(("`zvar'_psw","matrix"),return_result_item(`eqn',"`poolstack'_ps","ps_weights","`rep'"))
			mata: `A'.put(("`zvar'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
		}
	}

	mata: (`mname'.estAA).put(("`spec'","`rep'"),`A')
	
	// no longer needed
	foreach obj in `A' `eqn' {
		cap mata: mata drop `obj'
	}

end

// estimates and stores mean/median estimates across resamples
program medmean_and_store, eclass

	version 16

	syntax [anything] [if] [in] ,				///
						[						///
							mname(name)			///
							spec(string)		///
							title(string)		///
							medmean(string)		///
							vce(string) 		///
							TEffect(name)		/// either ATE, ATET or ATEU
						]
		
	mata: st_local("model",`mname'.model)
	
	marksample touse

	tempname b V bagg Vagg Vi
	tempname bvec sbvec bmed Vvec sVvec Vmed
	tempname nlltrim nultrim trimval
	tempvar esample
	tempname B
	
	// initialize
	mata: st_local("nameD",invtokens(`mname'.nameD))
	local K : word count `nameD'
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: `B' = AssociativeArray()
	local isodd = mod(`nreps',2)
	local medrow = ceil(`nreps'/2)
	local N = 0
	
	// bvec a misnomer - usually a vector, but can be a matrix if multiple D variables
	mata: `bvec' = J(`nreps',`K',0)
	mata: `bagg' = J(1,`K',0)
	forvalues m=1/`nreps' {
		mata: check_spec(`mname',"`spec'","`m'")
		mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
		mata: `bvec'[`m',.] = `B'.get(("b","post"))
		// row/colnames etc. - need to do this only once
		if `m'==1 {
			mata: st_local("depvar",`B'.get(("depvar","post")))
			// retrieve locals; if empty, will be ""
			local list_local y0 y0_m y1 y1_m d d_m d0 d0_m d1 d1_m z z_m yvar dvar vce vcetype clustvar
			foreach obj in `list_local' {
				mata: st_local("`obj'",`B'.get(("`obj'","local")))
			}
			// retrieve scalars (as locals)
			local list_scalar
			if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
			foreach obj in `list_scalar' {
				mata: st_local("`obj'",strofreal(`B'.get(("`obj'","scalar"))))
			}
		}
		// possible that different estimation samples have different #obs
		qui count if `mname'_sample_`m'==1
		local N = `N' + r(N)
	}
	local N = round(`N'/`nreps')
	if "`medmean'"=="mn" {
		// mean beta
		mata: `bagg' = mean(`bvec')
		mata: st_matrix("`bagg'",`bagg')
	}
	else if "`medmean'"=="md" {
		// median beta
		forvalues k=1/`K' {
			// leave order of bvec unchanged
			mata: `sbvec' = sort(`bvec',`k')
			// mata: _sort(`bvec',`k')
			if `isodd' {
				mata: `bagg'[1,`k'] = `sbvec'[`medrow',`k']
			}
			else {
				mata: `bagg'[1,`k'] = (`sbvec'[`medrow',`k'] + `sbvec'[`medrow'+1,`k'])/2
			}
		}
		mata: st_matrix("`bagg'",`bagg')
	}
	else {
		di as err "replay_estimate error - unrecognized option `medmean'"
		exit 198
	}
	
	mata: `Vagg' = J(`K',`K',0)
	mata: `Vvec' = J(`nreps',1,0)
	if "`medmean'"=="mn" {
		// harmonic mean
		// inefficient - does off-diagonals twice
		forvalues m=1/`nreps' {
			mata: check_spec(`mname',"`spec'","`m'")
			mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
			mata: `Vi' = `B'.get(("V","post"))
			forvalues j=1/`K' {
				forvalues k=1/`K' {
					// abs(.) needed?
					mata: `Vi'[`j',`k'] = `Vi'[`j',`k'] + abs((`bvec'[`m',`j'] - `bagg'[1,`j'])*(`bvec'[`m',`k'] - `bagg'[1,`k']))
				}
			}
			mata: `Vagg' = `Vagg' + 1:/`Vi'
		}
		mata: `Vagg' = `nreps' :/ `Vagg'
		mata: st_matrix("`Vagg'",`Vagg')
	}
	else if "`medmean'"=="md" {
		// median VCV
		// inefficient - does off-diagonals twice
		forvalues j=1/`K' {
			forvalues k=1/`K' {
				forvalues m=1/`nreps' {
					mata: check_spec(`mname',"`spec'","`m'")
					mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
					mata: `Vi' = `B'.get(("V","post"))
					mata: `Vvec'[`m'] = `Vi'[`j',`k']
				}
				// adjustment as per
				// https://docs.doubleml.org/stable/guide/resampling.html#repeated-cross-fitting-with-k-folds-and-m-repetition
				// (generalized to multiple D variables)
				mata: `Vvec' = `Vvec' + abs((`bvec'[.,`j'] :- `bagg'[1,`j']):*(`bvec'[.,`k'] :- `bagg'[1,`k']))
				// leave order of Vvec unchanged
				mata: `sVvec' = sort(`Vvec',1)
				// mata: _sort(`Vvec',1)
				if `isodd' {
					mata: `Vagg'[`j',`k'] = `sVvec'[`medrow',1]
				}
				else {
					mata: `Vagg'[`j',`k'] = (`sVvec'[`medrow',1] + `sVvec'[`medrow'+1,1])/2
				}
			}
		}
		mata: st_matrix("`Vagg'",`Vagg')
	}
	else {
		di as err "replay_estimate error - unrecognized option `medmean'"
		exit 198
	}
	
	// count trim instances
	mata: `nlltrim' = 0
	mata: `nultrim' = 0
	forvalues m=1/`nreps' {
		mata: check_spec(`mname',"`spec'","`m'")
		mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
		mata: `nlltrim' = `nlltrim' + (`B'.get(("lltrim","scalar"))>0)
		mata: `nultrim' = `nultrim' + (`B'.get(("ultrim","scalar"))>0)
	}
	// retrieve from last rep (will be stored in all)
	mata: `trimval' = `B'.get(("trim","scalar"))

	tempname A
	mata: `A' = AssociativeArray()
	mata: `A'.reinit("string",2)
	mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
	
	mata: `A'.put(("N","post"),`N')
	mata: `A'.put(("b","post"),`bagg')
	mata: `A'.put(("V","post"),`Vagg')
	mata: `A'.put(("depvar","post"),"`depvar'")
	mata: `A'.put(("b_resamples","matrix"),`bvec')
	
	// store locals
	local list_local title yvar dvar y0 y1 vce vcetype teffect
	if "`model'"=="interactive" {
		local list_local `list_local' d
	}
	else {
		local list_local `list_local' d0 d1 z
	}
	if "`clustvar'"~=""		local list_local `list_local' clustvar
	foreach obj in `list_local' {
		mata: `A'.put(("`obj'","local"),"``obj''")
	}
	// special case - "_m" subscript doesn't apply to mean/median over resamplings
	// so store without resample subscript
	foreach obj in title y0 y1 d d0 d1 z {
		mata: `A'.put(("`obj'_m","local"),"``obj''")
	}
	
	// store scalars
	local trim `trimval'	// hack, to fix
	local list_scalar nlltrim nultrim trim
	if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
	foreach obj in `list_scalar' {
		mata: `A'.put(("`obj'","scalar"),``obj'')
	}

/*
	// additional estimation results
	tempname eqn
	mata: `eqn' = init_eStruct()
	// Y eqn results
	mata: `eqn' = (`mname'.eqnAA).get("`yname'")
	// pystacked final est (pystacked multi only)
	mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
	if `pystackedmulti' {
		// cap because won't exist for e.g. shortstack variable
		cap mata: `A'.put(("`y'_stack_final_est","local"), return_learner_item(`eqn',"`y'","stack_final_est"))
	}
	// ss results
	if "`spec'"=="ss" {
		mata: st_local("shortstack", `eqn'.shortstack)
		mata: `A'.put(("`yname'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
	}
	// ps results
	if "`spec'"=="ps" {
		mata: st_local("poolstack", `eqn'.poolstack)
		mata: `A'.put(("`yname'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
	}

	// D eqn results - uses vtilde names in d
	forvalues i=1/`numeqnD' {
		local dname : word `i' of `dnames'
		local vtilde : word `i' of `d'
		local vtilde_h `vtilde'

		mata: `eqn' = (`mname'.eqnAA).get("`dname'")
		// pystacked final est (pystacked multi only)
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		if `pystackedmulti' {
			// cap because won't exist for e.g. shortstack variable
			cap mata: `A'.put(("`vtilde'_stack_final_est","local"), return_learner_item(`eqn',"`vtilde'","stack_final_est"))
		}
		// ss results
		if "`spec'"=="ss" {
			mata: st_local("shortstack", `eqn'.shortstack)
			mata: `A'.put(("`dname'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
		}
		// ps results
		if "`spec'"=="ps" {
			mata: st_local("poolstack", `eqn'.poolstack)
			mata: `A'.put(("`dname'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
		}
	}
	if `fivflag'==0 {
		// Z eqn results; fiv won't enter
		forvalues i=1/`numeqnZ' {
			local zname : word `i' of `znames'
			local vtilde : word `i' of `z'
			mata: `eqn' = (`mname'.eqnAA).get("`zname'")
			// pystacked final est (pystacked multi only)
			mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
			if `pystackedmulti' {
				// cap because won't exist for e.g. shortstack variable
				cap mata: `A'.put(("`vtilde'_stack_final_est","local"), return_learner_item(`eqn',"`vtilde'","stack_final_est"))
			}
			// ss results
			if "`spec'"=="ss" {
				mata: st_local("shortstack", `eqn'.shortstack)
				mata: `A'.put(("`zname'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
			}
			// ps results
			if "`spec'"=="ps" {
				mata: st_local("poolstack", `eqn'.poolstack)
				mata: `A'.put(("`zname'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
			}
		}
	}
*/


	// store AA with median/mean results
	mata: (`mname'.estAA).put(("`spec'","`medmean'"),`A')
	
	// no longer needed
	foreach obj in `A' `B' `bagg' `bvec' `sbvec' `Vagg' `Vvec' `sVvec' `Vi' `nlltrim' `nultrim' `trimval' {
		cap mata: mata drop `obj'
	}
end


* code below currently supports only a single treatment variable, but coded for multiple variables in places
program replay_estimate, eclass

	version 16

	syntax [anything] [if] [in] ,				///
						[						///
							mname(name)			///
							spec(string)		///
							rep(string)			///
						]
		
	mata: st_local("model",`mname'.model)
	
	marksample touse

	// replay
			
	tempname B keys isscalar islocal ismatrix

	mata: `B' = AssociativeArray()
	mata: check_spec(`mname',"`spec'","`rep'")
	mata: `B' = (`mname'.estAA).get(("`spec'","`rep'"))
	mata: `keys' = `B'.keys()
	mata: st_local("nentries",strofreal(rows(`keys')))
	mata: `isscalar'	= (`keys'[.,2] :== "scalar")
	mata: `islocal'		= (`keys'[.,2] :== "local")
	mata: `ismatrix'	= (`keys'[.,2] :== "matrix")
	
	tempname b V
	mata: st_matrix("`b'",`B'.get(("b","post")))
	mata: st_matrix("`V'",`B'.get(("V","post")))
	mata: st_local("N",strofreal(`B'.get(("N","post"))))
	mata: st_local("depvar",`B'.get(("depvar","post")))
	
	mata: st_local("yvar",`B'.get(("yvar","local")))
	mata: st_local("dvar",`B'.get(("dvar","local")))
	
	matrix rownames `b' = `depvar'
	matrix colnames `b' = `dvar'
 	matrix colnames `V' = `dvar'
	matrix rownames `V' = `dvar'
	
	tempvar esample
	cap gen byte `esample' = `mname'_sample_`rep'
	if _rc>0 {
		// sample variable doesn't exist; ignore
		local esample
	}
	
	ereturn clear
	ereturn post `b' `V', depname(`depvar') obs(`N') esample(`esample')
	
	ereturn local cmd ddml
	ereturn local model `model'
	ereturn local rep `rep'
	ereturn local spec `spec'
	ereturn local tmname `mname'
	ereturn local teffect `teffect'
	
	// extract and post scalars, locals, matrices
	forvalues i=1/`nentries' {
		mata: st_local("topost",strofreal(`isscalar'[`i']))
		if `topost' {
			mata: st_local("sname",substr(`keys'[`i',1],1,32))
			mata: st_numscalar("e(`sname')",`B'.get(`keys'[`i',.]))
		}
	}
	forvalues i=1/`nentries' {
		mata: st_local("topost",strofreal(`islocal'[`i']))
		if `topost' {
			mata: st_local("lname",substr(`keys'[`i',1],1,32))
			mata: st_global("e(`lname')",`B'.get(`keys'[`i',.]))
		}
	}
	forvalues i=1/`nentries' {
		mata: st_local("topost",strofreal(`ismatrix'[`i']))
		if `topost' {
			mata: st_local("tmname",substr(`keys'[`i',1],1,32))
			mata: st_matrix("e(`tmname')",`B'.get(`keys'[`i',.]))
		}
	}
	
	// no longer needed
	foreach obj in `B' `keys' `isscalar' `islocal' `ismatrix' {
		cap mata: mata drop `obj'
	}
	
	// display results
	di as text "`e(title)'"
	di as text "E[y|X,D=0]" _col(14) "= " as res "`e(y0_m)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	di as text "E[y|X,D=1]" _col(14) "= " as res "`e(y1_m)'"
	if "`e(model)'"=="interactive" {
		di as text "E[D|X]" _col(14)  "= " as res "`e(d_m)'"
	}
	else {
		di as text "E[D|X,Z=0]" _col(14)  "= " as res "`e(d0_m)'"
		di as text "E[D|X,Z=1]" _col(14)  "= " as res "`e(d1_m)'"
		di as text "E[Z|X]" _col(14)  "= " as res "`e(z_m)'"
	}
	ereturn display
	
	// report warning if clustered SEs requested but doesn't match clustered crossfitting
	mata: st_local("fclustvar",`mname'.fclustvar)
	if "`e(clustvar)'"~="" {
		if "`fclustvar'"=="" {
			di as res "Warning" as text ": crossfit folds do not respect cluster structure used for VCE."
		}
		else if "`fclustvar'"~="`e(clustvar)'" {
			di as res "Warning" as text ": cluster variable for VCE does not match cluster variable for crossfit folds."
		}
	}
	
	// warn if any values trimmed
	if e(lltrim)>0 & e(lltrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(lltrim) " propensity scores trimmed to lower limit " e(trim) "."
	}
	if e(ultrim) & e(ultrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(ultrim) " propensity scores trimmed to upper limit " 1-e(trim) "."
	}
	// for mean/median over resamples
	if e(nlltrim)>0 & e(nlltrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(nlltrim) " resamples had propensity scores trimmed to lower limit " e(trim) "."
	}
	if e(ultrim) & e(ultrim)<. {
		di as res "Warning" as text ": " _c
		di as text e(nultrim) " resamples had propensity scores trimmed to upper limit " 1-e(trim) "."
	}
	
end

