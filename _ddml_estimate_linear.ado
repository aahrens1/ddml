*! ddml v1.2
*! last edited: 27 feb 2023
*! authors: aa/ms

program _ddml_estimate_linear, eclass sortpreserve
	syntax [anything] [if] [in] ,					/// 
								[					///
								y(varname)			/// for estimating by hand...
								d(varlist)			/// 
								z(varlist)			///
								dh(varname)			///
								shortstack			///
								poolstack			///
								* ]

	if "`shortstack'"~="" {
		// restack
		_ddml_estimate_stacking `anything' `if' `in', ss `options'
	}
	if "`poolstack'"~="" {
		// restack
		_ddml_estimate_stacking `anything' `if' `in', ps `options'
	}
	if "`y'`d'`z'`dh'"=="" {
		// main program for estimation
		_ddml_estimate_main `anything' `if' `in', `options'
	}
	else {
		// a single user-specified estimation
		_ddml_estimate_single `anything' `if' `in', y(`y') d(`d') z(`z') dh(`dh') `options'
	}
end

// (re-)estimate shortstacked
program _ddml_estimate_stacking, eclass sortpreserve
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
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	// reps = total number of reps; crossfitted = reps done so far (=0 if none)
	mata: st_local("reps", strofreal(`mname'.nreps))
	mata: st_local("kfolds", strofreal(`mname'.kfolds))
	mata: st_local("crossfitted", strofreal(`mname'.crossfitted))
	// assume stacking available for all
	local stackflag = 1
	
	// with pystacked, will always be a single vtilde and eqn struct for every variable
	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	// used for checking
	mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
	if `pystackedmulti' {
		mata: st_local("vtildeY",invtokens(`eqn'.vtlist))
		local namelist	`nameY'
		local vtlist	`vtildeY'
	}
	else {
		di as err "error - restacking of `nameY' requires pystacked to be the sole learner"
		local stackflag	= 0
	}

	// can be multiple D eqns for multiple D vars, and similarly for Z
	foreach vname in `nameD' `nameZ' {
		mata: `eqn' = (`mname'.eqnAA).get("`vname'")
		// used for checking
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		if `pystackedmulti' {
			mata: st_local("vtilde",invtokens(`eqn'.vtlist))
			local namelist	`namelist' `vname'
			local vtlist	`vtlist' `vtilde'
		}
		else {
			di as err "error - restacking of `vname' requires pystacked to be the sole learner"
			local stackflag	= 0
		}
	}

	tempname sweights N_folds mse_folds
	local numvts : word count `vtlist'
	// loop through vtildes
	forvalues i=1/`numvts' {
		local vname :	word `i' of `namelist'
		local vtilde :	word `i' of `vtlist'

		mata: `eqn' = (`mname'.eqnAA).get("`vname'")
		mata: st_local("base_est",return_learner_item(`eqn',"`vtilde'","stack_base_est"))
		mata: st_local("lieflag", strofreal(`eqn'.lieflag))
		local nlearners : word count `base_est'
		// check if previously stacked; required for poolstacking
		if `ssflag' {
			mata: st_local("shortstack", `eqn'.shortstack)
			if "`shortstack'"=="" {
				// not previously shortstacked, set local and struct field to default
				local shortstack `vname'
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
				`qui' _ddml_nnls `vname' `learner_list', finalest(`finalest') if `touse'
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
				getmata (`vname' `learner_list')=`y_stacking_cv', force replace
				`qui' _ddml_nnls `vname' `learner_list', finalest(`finalest')
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
			get_stack_stats if `touse', kfolds(`kfolds') fid(`mname'_fid_`m') vname(`vname') vhat(`yhat')
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
	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								y(varname)			/// for estimating by hand...
								d(varlist)			/// 
								z(varlist)			///
								dh(varname)			///
								ROBust				///
								CLUster(varname)	///
								vce(string)			///
								NOConstant			///
								* ]

	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'

	// checks
	if "`y'"=="" {
		di as err "option y(.) missing"
		exit 198
	}
	if "`d'"=="" {
		di as err "option d(.) missing"
		exit 198
	}
	if "`model'"=="partial" {
		if "`z'"~="" {
			di as err "z(.) should be empty for model=partial"
			exit 198
		}
		local numD : word count `d'
		if `numD'~=`numeqnD' {
			di as err "error - model has `numeqnD' D variables but `numD' specified"
			exit 198
		}
	}
	else if "`model'"=="iv" {
		if "`z'"=="" {
			di as err "option z(.) missing"
			exit 198
		}
	}
	else if "`model'"=="fiv" {
		if "`z'"~="" {
			di as err "z(.) should be empty for model=fiv"
			exit 198
		}
		if "`dh'"=="" {
			di as err "option dh(.) missing"
			exit 198
		}
	}
	else {
		di as err "error - unknown model `model'"
		exit 198
	}

	// residualization
	tempvar yt
	qui gen double `yt' = `nameY' - `y'
	if "`model'"~="fiv" {
		// applies to plm and iv models
		forvalues i=1/`numeqnD' {
			tempvar d_resid
			local dname_i : word `i' of `nameD'
			local dtilde_i : word `i' of `d'
			qui gen double `d_resid' = `dname_i' - `dtilde_i'
			local d_resid_list `d_resid_list' `d_resid'
		}
		// if plm model, won't enter
		forvalues i=1/`numeqnZ' {
			tempvar z_resid
			local zname_i : word `i' of `nameZ'
			local ztilde_i : word `i' of `z'
			qui gen double `z_resid' = `zname_i' - `ztilde_i'
			local z_resid_list `z_resid_list' `z_resid'
		}
	}
	else {
		// create new D and Z
		tempvar d_resid_list z_resid_list
		qui gen double `d_resid_list'	= `nameD' - `dh'
		qui gen double `z_resid_list'	= `d'     - `dh'

	}

	if "`robust'"!=""	local vce robust
	if "`cluster'"~=""	local vce cluster `cluster'
	tempname b V
	tempvar esample
	marksample touse
	if "`model'"=="partial" {
		qui reg `yt' `d_resid_list' if `touse', vce(`vce') `noconstant'
	}
	else {
		qui reg `yt' `d_resid_list' (`z_resid_list') if `touse', vce(`vce') `noconstant'
	}

	mat `b'=e(b)
	mat `V'=e(V)
	local cnames			`nameD'
	if "`noconstant'"=="" {
		local cnames		`cnames' _cons
	}
	mat colnames `b' = `cnames'
	mat colnames `V' = `cnames'
	mat rownames `V' = `cnames'
	local N=e(N)
	local vcetype			`e(vcetype)'
	local clustvar			`e(clustvar)'
	qui gen byte `esample'=e(sample)
	ereturn post `b' `V', depname(`nameY') obs(`N') esample(`esample')
	ereturn local cmd		ddml
	ereturn local model		`model'
	ereturn local dh		`dh'
	ereturn local z			`z'
	ereturn local d			`d'
	ereturn local y			`y'
	ereturn local vce		`vce'
	ereturn local vcetype	`vcetype'
	
	di
	di as text "y-E[y|X]" _col(11) "= " as res "`e(y)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	if "`e(model)'"~="fiv" {
		di as text "D-" _c
	}
	di as text "E[D|X,Z]" _col(11)  "= " as res "`e(d)'"
	if "`e(model)'" == "iv" {
		di as text "Z-E[Z|X]" _col(11) "= " as res "`e(z)'"
	}
	else if "`e(model)'" == "fiv" {
		di as text "E[D|X]" _col(11) "= " as res "`e(dh)'"
	}
	if "`e(model)'" == "fiv" {
		di as text "Orthogonalized D = D - E[D|X]; optimal IV = E[D|X,Z] - E[D|X]."
	}
	ereturn display

end

// main estimation program
program _ddml_estimate_main
	syntax namelist(name=mname) [if] [in] ,			/// 
								[					///
								ROBust				///
								CLUster(varname)	///
								NOConstant			/// suppress constant in estimation
								SHOWConstant		/// display constant in summary table
								vce(string)			///
								allcombos			/// estimate and show all combinations (dfn changed below)
								NOTable				/// suppress summary table
								clear				/// deletes all tilde-variables (to be implemented)
								spec(string)		/// specification to post/display
								REP(string)			/// resampling iteration to post/display or mean/median
								replay				/// model has been estimated, just display results
								debug				///
								NOIsily				///
								* ]

	if "`debug'`noisily'"==""	local qui qui
	
	marksample touse
	
	// consflag
	local consflag = ("`noconstant'"=="")
	local showconsflag = ("`showconstant'"~="" & `consflag')	// ignored if nocons
	// replay existing results
	local replayflag = "`replay'"~=""
	// display summary table
	local tableflag = "`notable'"==""
	// request estimation/reporting of all combinations
	local doallcombos = "`allcombos'"~=""
	// remaining macro flags
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
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
	local ivflag	= "`model'"=="iv"
	local fivflag	= "`model'"=="fiv"

	** standard errors
	// local vce is the argument to the Stata option vce(.)
	if "`robust'"!=""	local vce robust
	if "`cluster'"~=""	local vce cluster `cluster'
	
	if ~`crossfitted' {
		di as err "ddml model not cross-fitted; call `ddml crossfit` first"
		exit 198
	}

	// default spec is ss if available, otherwise mse if ncombos>1, otherwise 1 since ncombos=1
	if "`spec'"=="" & `ssflag' {
		local spec "ss"
	}
	else if "`spec'"=="" & `ncombos'>1 {
		local spec "mse"
	}
	else if "`spec'"=="" {
		local spec 1
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
	if "`spec'"~="ss" & "`spec'"~="ps" & "`spec'"~="mse" & real("`spec'")==. {
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

	// shortstack and poolstack names
	if `ssflag' {
		local Yss `nameY'_ss
		foreach var in `nameD' {
			local Dss `Dss' `var'_ss
		}
		if `fivflag' {
			foreach var in `nameD' {
				local Zss `Zss' `var'_h_ss
			}
		}
		else {
			foreach var in `nameZ' {
				local Zss `Zss' `var'_ss
			}
		}
	}
	if `psflag' {
		local Yps `nameY'_ps
		foreach var in `nameD' {
			local Dps `Dps' `var'_ps
		}
		if `fivflag' {
			foreach var in `nameD' {
				local Zps `Zps' `var'_h_ps
			}
		}
		else {
			foreach var in `nameZ' {
				local Zps `Zps' `var'_ps
			}
		}
	}
	
	// number of possible combos; if only one spec, will replace "mse" label
	mata: st_local("poss_combos",strofreal(return_ncombos(`mname')))
	
	// multiple specs
	// text locals control messages and lookups
	if `replayflag' {
		local spectext	`spec'
	}
	else if `poss_combos'>1 {
		local msetext	"min-mse "
		local MSEtext	"Min MSE "
		local spectext	mse
		local otext		mse
	}
	else {
		local msetext
		local MSEtext
		local spectext		1
		local otext			1
		// if only one combo, enter allcombos code so estimation gets numbered
		local doallcombos	=1
	}
	
	************* ESTIMATE ************
	
	if `estimated'==0 {
		// enter if no estimates exist

		// Loop over resamples and estimate/save the min mse and ss/ps models for each
		forvalues m=1/`nreps' {
			
			// text used in output below
			// multiple reps
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			
			// reset locals
			local Yopt
			local Dopt
			local Zopt
			
			*** retrieve best model
			mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
			mata: st_local("Yopt",return_learner_item(`eqn',"opt","`m'"))
			
			foreach var in `nameD' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneDopt",return_learner_item(`eqn',"opt","`m'"))
				local Dopt `Dopt' `oneDopt'
				// DHopt is stored in list Zopt
				if "`model'"=="fiv" {
					mata: st_local("oneDHopt",return_learner_item(`eqn',"opt_h","`m'"))
					local Zopt `Zopt' `oneDHopt'
				}
			}
			
			// nameZ is empty for fiv model
			foreach var in `nameZ' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneZopt",return_learner_item(`eqn',"opt","`m'"))
				local Zopt `Zopt' `oneZopt'
			}
		
			local title "`MSEtext'DDML model`stext'"
			`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
					`noconstant' vce(`vce')								///
					y(`Yopt') yname(`nameY')							///
					d(`Dopt') dnames(`nameD')			 				///
					z(`Zopt') znames(`nameZ')							///
					mname(`mname') spec(`spectext') rep(`m')			///
					title(`title')
		
			// estimate shortstack for this rep
			if `ssflag' {
				local title "Shortstack DDML model`stext'"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
						`noconstant' vce(`vce')								///
						y(`Yss') yname(`nameY')								///
						d(`Dss') dnames(`nameD') 							///
						z(`Zss') znames(`nameZ') 							///
						mname(`mname') spec(ss) rep(`m')					///
						title(`title')
			}
			// estimate poolstack for this rep
			if `psflag' {
				local title "Poolstack DDML model`stext'"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
						`noconstant' vce(`vce')								///
						y(`Yps') yname(`nameY')								///
						d(`Dps') dnames(`nameD') 							///
						z(`Zps') znames(`nameZ')							///
						mname(`mname') spec(ps) rep(`m')					///
						title(`title')
			}
		}
		
		// have looped over reps to get each optimal model and shortstack per rep
		// now aggregate over reps to get mean/median
		if `nreps' > 1 {
			local title "Mean over `nreps' `msetext'resamples"
			`qui' medmean_and_store, mname(`mname') spec(`spectext') medmean(mn) title(`title') `noconstant'
			local title "Median over `nreps' `msetext'resamples"
			`qui' medmean_and_store, mname(`mname') spec(`spectext') medmean(md) title(`title') `noconstant'
			// shortstack
			if `ssflag' {
				local title "Shortstack DDML model (mean over `nreps' resamples)"
				`qui' medmean_and_store, mname(`mname') spec(ss) medmean(mn) title(`title') `noconstant'
				local title "Shortstack DDML model (median over `nreps' resamples)"
				`qui' medmean_and_store, mname(`mname') spec(ss) medmean(md) title(`title') `noconstant'
			}
			// poolstack
			if `psflag' {
				local title "Poolstack DDML model (mean over `nreps' resamples)"
				`qui' medmean_and_store, mname(`mname') spec(ps) medmean(mn) title(`title') `noconstant'
				local title "Poolstack DDML model (median over `nreps' resamples)"
				`qui' medmean_and_store, mname(`mname') spec(ps) medmean(md) title(`title') `noconstant'
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
		// get varlists
		_ddml_make_varlists, mname(`mname')
		local yvars `r(yvars)'
		local dvars `r(dvars)'
		local zvars `r(zvars)'
		local zpos	`r(zpos_start)'
	
		// obtain all combinations
		_ddml_allcombos `yvars' - `dvars' - `zvars' ,	///
			`debug'										///
			zpos(`zpos')		 						///
			addprefix("")
		
		local ncombos = r(ncombos)
		local tokenlen = `ncombos'*2
		local ylist `r(ystr)'
		local Dlist `r(dstr)'
		local Zlist `r(zstr)' 
		
		tempname nmat bmat semat
		mata: `nmat' = J(`ncombos',3,"")
		mata: `bmat' = J(`ncombos'*`nreps',`numeqnD'+`consflag',.)
		mata: `semat' = J(`ncombos'*`nreps',`numeqnD'+`consflag',.)
		
		// simplest if put into a Mata string matrix
		tokenize `ylist' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',1] = strtrim("``idx''")
		}
		tokenize `Dlist' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',2] = strtrim("``idx''")
		}
		tokenize `Zlist' , parse("-")
		forvalues i=1/`ncombos' {
			local idx = 2*`i'-1
			mata: `nmat'[`i',3] = strtrim("``idx''")
		}
		
		forvalues m=1/`nreps' {
			
			// reset locals
			local Yopt
			local Dopt
			local Zopt
			
			*** retrieve best model
			mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
			mata: st_local("Yopt",return_learner_item(`eqn',"opt","`m'"))
			foreach var in `nameD' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneDopt",return_learner_item(`eqn',"opt","`m'"))
				local Dopt `Dopt' `oneDopt'
				// DHopt is stored in list Zopt
				if "`model'"=="fiv" {
					mata: st_local("oneDHopt",return_learner_item(`eqn',"opt_h","`m'"))
					local Zopt `Zopt' `oneDHopt'
				}
			}
			// nameZ is empty for fiv model
			foreach var in `nameZ' {
				mata: `eqn' = (`mname'.eqnAA).get("`var'")
				mata: st_local("oneZopt",return_learner_item(`eqn',"opt","`m'"))
				local Zopt `Zopt' `oneZopt'
			}
			
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			
			forvalues i = 1/`ncombos' {
				mata: st_local("y",`nmat'[`i',1])
				mata: st_local("d",`nmat'[`i',2])
				mata: st_local("z",`nmat'[`i',3])
				// check if opt for this resample; note === for D so order doesn't matter
				local isopt
				local isYopt : list Yopt == y
				local isDopt : list Dopt === d
				local isZopt : list Zopt === z
				local title "DDML model, specification `i'`stext'"
				if `isYopt' & `isDopt' & `isZopt' {
					local optspec`m' = `i'
					local isopt *
					local title `MSEtext'`title'
					// save in AA
					mata: (`mname'.estAA).put(("optspec","`m'"),"`i'")
				}
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
						`noconstant' vce(`vce')								///
						y(`y') yname(`nameY')								///
						d(`d') dnames(`nameD')						 		///
						z(`z') znames(`nameZ')								///
						mname(`mname') spec(`i') rep(`m')					///
						title(`title')
				mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(bmat)")
				mata: `semat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(semat)")
			}
			
			if `ssflag' {
				local title "Shortstack DDML model`stext'"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
						`noconstant' vce(`vce')								///
						y(`Yss') yname(`nameY')								///
						d(`Dss') dnames(`nameD')							///
						z(`Zss') znames(`nameZ')							///
						mname(`mname') spec(ss) rep(`m')					///
						title(`title')
			}
			if `psflag' {
				local title "Poolstack DDML model`stext'"
				`qui' estimate_and_store if `mname'_sample_`m' & `touse',	///
						`noconstant' vce(`vce')								///
						y(`Yps') yname(`nameY')								///
						d(`Dps') dnames(`nameD') 							///
						z(`Zps') znames(`nameZ')							///
						mname(`mname') spec(ps) rep(`m')					///
						title(`title')
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
	
	// optional table of results
	if `tableflag' {
		di
		di as text "DDML estimation results:"
		di as text "spec  r" %14s "Y learner" _c
		forvalues j=1/`numeqnD' {
			di as text %14s "D learner" %10s "b" %10s "SE" _c
		}
		if `showconsflag' {
			di as text %10s "_cons" %10s "SE" _c
		}
		if "`model'"=="fiv" {
			forvalues j=1/`numeqnD' {
				di as text %14s "DH learner" _c
			}
		}
		forvalues j=1/`numeqnZ' {
			di as text %14s "Z learner" _c
		}
		di
		forvalues m=1/`nreps' {
			if `doallcombos' {
				// all combos available, so loop through
				forvalues i=1/`ncombos' {
					mata: st_local("yt",abbrev(`nmat'[`i',1],13))
					mata: st_local("dtlist",invtokens(abbrev(tokens(`nmat'[`i',2]),13)))
					mata: st_local("ztlist",`nmat'[`i',3])
					if "`ztlist'"~="" {
						mata: st_local("ztlist",invtokens(abbrev(tokens(`nmat'[`i',3]),13)))
					}
					if "`optspec`m''"=="`i'" & `ncombos'>1 {
						// print * for opt (lowest MSE combo) if more than 1 combo
						di "*" _c
					}
					else {
						di " " _c
					}
					local specrep "`: di %3.0f `i' %3.0f `m''"
					local rcmd stata ddml estimate, mname(`mname') spec(`i') rep(`m') notable replay `noconstant'
					di %6s "{`rcmd':`specrep'}" _c
					di as res %14s "`yt'" _c
					forvalues j=1/`numeqnD' {
						local vt : word `j' of `dtlist'
						mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
						mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
						di as res %14s "`vt'" _c
						di as res %10.3f `b' _c
						local pse (`: di %6.3f `se'')
						di as res %10s "`pse'" _c
					}
					if `showconsflag' {
						// set j by hand; cons is in last column
						local j = `numeqnD' + 1
						mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
						mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
						di as res %10.3f `b' _c
						local pse (`: di %6.3f `se'')
						di as res %10s "`pse'" _c
					}
					forvalues j=1/`numeqnZ' {
						local vt : word `j' of `ztlist'
						di as res %14s "`vt'" _c
					}
					if "`model'"=="fiv" {
						forvalues j=1/`numeqnD' {
							local vt : word `j' of `ztlist'
							di as res %14s "`vt'" _c
						}
					}
					di
				}
			}
			else {
				// only mse/ss specs available/reported
				`qui' replay_estimate, mname(`mname') spec(`spectext') rep(`m') `noconstant'
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "`otext'" %3.0f `m''"
				local rcmd stata ddml estimate, mname(`mname') spec(`spectext') rep(`m') notable replay `noconstant'
				di %6s "{`rcmd':`specrep'}" _c
				mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
				mata: st_local("yt",return_learner_item(`eqn',"opt","`m'"))
				di as res %14s abbrev("`yt'",13) _c
				forvalues j=1/`numeqnD' {
					local dd : word `j' of `nameD'
					mata: `eqn' = (`mname'.eqnAA).get("`dd'")
					mata: st_local("vt",return_learner_item(`eqn',"opt","`m'"))
					di as res %14s abbrev("`vt'",13) _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if `showconsflag' {
					// set j by hand; cons is in last column
					local j = `numeqnD' + 1
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				forvalues j=1/`numeqnZ' {
					local zz : word `j' of `nameZ'
					mata: `eqn' = (`mname'.eqnAA).get("`zz'")
					mata: st_local("vt",return_learner_item(`eqn',"opt","`m'"))
					di as res %14s abbrev("`vt'",13) _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						local dd : word `j' of `nameD'
						mata: `eqn' = (`mname'.eqnAA).get("`dd'")
						mata: st_local("vt",return_learner_item(`eqn',"opt_h","`m'"))
						di as res %14s abbrev("`vt'",13) _c
					}
				}
				di
			}
			if `ssflag' {
				`qui' replay_estimate, mname(`mname') spec(ss) rep(`m') `noconstant'
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ss" %3.0f `m''"
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`m') notable replay `noconstant'
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[shortstack]" _c
				forvalues j=1/`numeqnD' {
					di as res %14s "[ss]" _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if `showconsflag' {
					// set j by hand; cons is in last column
					local j = `numeqnD' + 1
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						di as res %14s "[ss]" _c
					}
				}
				forvalues j=1/`numeqnZ' {
					di as res %14s "[ss]" _c
				}
				di
			}
			if `psflag' {
				`qui' replay_estimate, mname(`mname') spec(ps) rep(`m') `noconstant'
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ps" %3.0f `m''"
				local rcmd stata ddml estimate, mname(`mname') spec(ps) rep(`m') notable replay `noconstant'
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[poolstack]" _c
				forvalues j=1/`numeqnD' {
					di as res %14s "[ps]" _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if `showconsflag' {
					// set j by hand; cons is in last column
					local j = `numeqnD' + 1
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						di as res %14s "[ps]" _c
					}
				}
				forvalues j=1/`numeqnZ' {
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
		di as text "Mean/med    Y learner" _c
		forvalues j=1/`numeqnD' {
			di as text %14s "D learner" %10s "b" %10s "SE" _c
		}
		if "`model'"=="fiv" {
			forvalues j=1/`numeqnD' {
				di as text %14s "DH learner" _c
			}
		}
		forvalues j=1/`numeqnZ' {
			di as text %14s "Z learner" _c
		}
		di
		foreach medmean in mn md {
			** mean and median over mse
			// force noconstant with mean/median
			`qui' replay_estimate, mname(`mname') spec(`spectext') rep(`medmean') noconstant
			tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
			mat `btemp' = e(b)
			mat `Vtemp' = e(V)
			local specrep "`: di " " %3s "`spectext'" %3s "`medmean'"'"
			// force noconstant with mean/median
			local rcmd stata ddml estimate, mname(`mname') spec(`spectext') rep(`medmean') notable replay `noconstant'
			di %6s "{`rcmd':`specrep'}" _c
			di as res %14s "[min-mse]" _c
			forvalues j=1/`numeqnD' {
				di as res %14s "[mse]" _c
				di as res %10.3f el(`btemp',1,`j') _c
				local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
				di as res %10s "`pse'" _c
			}
			if "`model'"=="fiv" {
				forvalues j=1/`numeqnD' {
					di as res %14s "[mse]" _c
				}
			}
			forvalues j=1/`numeqnZ' {
				di as res %14s "[mse]" _c
			}
			di
			if `ssflag' {
				// force noconstant with mean/median
				`qui' replay_estimate, mname(`mname') spec(ss) rep(`medmean') noconstant
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ss" %3s "`medmean'"'"
				// force noconstant with mean/median
				local rcmd stata ddml estimate, mname(`mname') spec(ss) rep(`medmean') notable replay noconstant
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[shortstack]" _c
				forvalues j=1/`numeqnD' {
					di as res %14s "[ss]" _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						di as res %14s "[ss]" _c
					}
				}
				forvalues j=1/`numeqnZ' {
					di as res %14s "[ss]" _c
				}
				di
			}
			if `psflag' {
				// force noconstant with mean/median
				`qui' replay_estimate, mname(`mname') spec(ps) rep(`medmean') noconstant
				tempname btemp Vtemp	// pre-Stata 16 doesn't allow el(e(b),1,1) etc.
				mat `btemp' = e(b)
				mat `Vtemp' = e(V)
				local specrep "`: di %4s "ps" %3s "`medmean'"'"
				// force noconstant with mean/median
				local rcmd stata ddml estimate, mname(`mname') spec(ps) rep(`medmean') notable replay noconstant
				di %6s "{`rcmd':`specrep'}" _c
				di as res %14s "[poolstack]" _c
				forvalues j=1/`numeqnD' {
					di as res %14s "[ps]" _c
					di as res %10.3f el(`btemp',1,`j') _c
					local pse (`: di %6.3f sqrt(el(`Vtemp',`j',`j'))')
					di as res %10s "`pse'" _c
				}
				if "`model'"=="fiv" {
					forvalues j=1/`numeqnD' {
						di as res %14s "[ps]" _c
					}
				}
				forvalues j=1/`numeqnZ' {
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
	// reach this point if only 1 spec and numbered, or mult specs and numbered, or 1 spec and MSE
	else if `poss_combos'==1 {
		local specdisp 1
	}
	else {
		local specdisp mse
	}

	di
	if ("`rep'"=="mn" | "`rep'"=="md") {
		// force noconstant with mean/median
		replay_estimate, mname(`mname') spec(`specdisp') rep(`rep') noconstant
	}
	else {
		replay_estimate, mname(`mname') spec(`specdisp') rep(`rep') `noconstant'
	}
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
			di %12s as text "`vn'" _col(15) _c
			qui sum `bhat'`i', detail
			di as res %10.4f r(mean) _c
			di as res %10.4f r(min) _c
			di as res %10.4f r(p25) _c
			di as res %10.4f r(p50) _c
			di as res %10.4f r(p75) _c
			di as res %10.4f r(max)
			local ++i
		}
	}
	
	// temp Mata object no longer needed
	foreach obj in `eqn' `nmat' `bmat' `semat' {
		cap mata: mata drop `obj'
	}


end

// adds model name prefixes to list of varnames
program define add_prefix, sclass
	syntax [anything] , prefix(name)

	// anything is a list of to-be-varnames that need prefix added to them
	foreach vn in `anything' {
		local vnames `vnames' `prefix'`vn' 
	}
	
	sreturn local vnames `vnames'
end

// adds rep number suffixes to list of varnames
program define add_suffix, sclass
	syntax [anything] , suffix(name)

	// anything is a list of to-be-varnames that need suffix added to them
	foreach vn in `anything' {
		local vnames `vnames' `vn'`suffix'
	}
	
	sreturn local vnames `vnames'
end

program define _ddml_make_varlists, rclass

	syntax [anything], mname(name)
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()

	// locals used below
	mata: st_local("model",`mname'.model)
	
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	local numeqnD	: word count `nameD'
	local numeqnZ	: word count `nameZ'
	
	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	local numlnrY : word count `vtlistY'
	
	if `numeqnD' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			tempname Dt_list
			local `Dt_list' `vtlistD'
			local Dorthog_lists `Dorthog_lists' `Dt_list'
		}
	}
	
	if `numeqnD' & `lieflag' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			foreach vn in `vtlistD' {
				local vtlistD_h `vtlistD_h' `vn'_h
			}
			tempname DHt_list
			local `DHt_list' `vtlistD_h'
			local DHorthog_lists `DHorthog_lists' `DHt_list'
		}
	}
	
	if `numeqnZ' {
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			tempname Zt_list
			local `Zt_list' `vtlistZ'
			local Zorthog_lists `Zorthog_lists' `Zt_list'
		}
	}

	local dash
	foreach vl in `Dorthog_lists' {
		local Dtilde `Dtilde' `dash' ``vl''
		local dash -
	}
	local dash
	foreach vl in `Zorthog_lists' {
		local Ztilde `Ztilde' `dash' ``vl''
		local dash -
	}
	local dash
	foreach vl in `DHorthog_lists' {
		local DHtilde `DHtilde' `dash' ``vl''
		local dash -
	}

	// clear from Mata
	mata: mata drop `eqn'
	
	return scalar dpos_end = `numeqnD' + 1
	return scalar dpos_start = 2
	if (`numeqnZ'>0) {
		return scalar zpos_start = `numeqnD' +2
		return scalar zpos_end = `numeqnD' + `numeqnZ' + 1
	}
	else if ("`model'"=="fiv") {		
		return scalar zpos_start = `numeqnD' +2
		// return scalar zpos_end = `numeqnD' + `numeqnDH' + 1
		// not sure this will work
		return scalar zpos_end = 2*`numeqnD' + 1
	}
	else {
		return scalar zpos_start = 0
		return scalar zpos_end = 0
	}
	return scalar numD = `numeqnD'
	return scalar numZ = `numeqnZ'
	return local yvars `vtlistY'
	return local dvars `Dtilde'
	return local zvars `Ztilde' `DHtilde'

end

// does OLS/IV and reports with substitute yname and dnames
program define estimate_and_store, eclass
	syntax [anything] [if] [in] , [				///
				y(name) yname(name)				/// predicted var, original var
				d(namelist) dnames(namelist)	/// predicted var, original var
				z(namelist) znames(namelist)	/// predicted var, original var
												/// if model=fiv, z=DH
				mname(name)						///
				spec(string)					///
				rep(string)						///
				vce(string)						///
				title(string)					///
				NOConstant						///
				*								///
				]

	local consflag = ("`noconstant'"=="")
	mata: st_local("model",`mname'.model)
	local ivflag	= "`model'"=="iv"
	local fivflag	= "`model'"=="fiv"
	local numeqnD	: word count `dnames'
	local numeqnZ	: word count `znames'
		
	marksample touse
	
	tempname A
	mata: `A' = AssociativeArray()
	mata: `A'.reinit("string",2)
	mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
	
	// add resample suffixes
	local y_m `y'_`rep'
	add_suffix `d', suffix("_`rep'")
	local d_m `s(vnames)'
	add_suffix `z', suffix("_`rep'")
	local z_m `s(vnames)'
	
	// residualize
	// y var
	tempvar y_resid
	qui gen double `y_resid' = `yname' - `y_m'
	// D vars
	forvalues i=1/`numeqnD' {
		tempvar d_resid
		local dname_i : word `i' of `dnames'
		local dtilde_i : word `i' of `d_m'
		qui gen double `d_resid' = `dname_i' - `dtilde_i'
		local d_resid_list `d_resid_list' `d_resid'
	}
	// Z vars, iv model
	if `ivflag' {
		forvalues i=1/`numeqnZ' {
			tempvar z_resid
			local zname_i : word `i' of `znames'
			local ztilde_i : word `i' of `z_m'
			qui gen double `z_resid' = `zname_i' - `ztilde_i'
			local z_resid_list `z_resid_list' `z_resid'
		}
	}
	// D and Z, fiv model (nb: awkward naming)
	if `fivflag' {
		tempvar d_resid_list z_resid_list
		qui gen double `d_resid_list' = `dnames' - `z_m'
		qui gen double `z_resid_list' = `d_m' - `z_m'
	}
	// estimate
	if ~`ivflag' & ~`fivflag' {
		qui reg `y_resid' `d_resid_list' if `touse', vce(`vce') `noconstant' `options'
	}
	else {
		// old-style regress syntax: put IVs in parentheses
		 qui reg `y_resid' `d_resid_list' (`z_resid_list') if `touse', vce(`vce') `noconstant' `options'
	}
	tempname b V
	mat `b' = e(b)
	mat `V' = e(V)
	local N = e(N)
	local vce		`e(vce)'
	local vcetype	`e(vcetype)'
	local clustvar	`e(clustvar)'
	local N_clust	=e(N_clust)
	if `fivflag' {
		local dh	`z'
		local dh_m	`z_m'
		// clear locals (IV info is stored as dh not z)
		local z
		local z_m
	}
	
	// store post objects
	mata: `A'.put(("N","post"),`N')
	mata: `A'.put(("b","post"),st_matrix("`b'"))
	mata: `A'.put(("V","post"),st_matrix("`V'"))
	mata: `A'.put(("depvar","post"),"`yname'")
	
	// for calling program
	ereturn clear
	mata: st_matrix("e(bmat)",st_matrix("`b'"))
	mata: st_matrix("e(semat)",sqrt(diagonal(st_matrix("`V'"))'))
	
	// store locals
	local list_local title y y_m d d_m yname dnames znames vce vcetype
	if "`model'"=="iv" {
		local list_local `list_local' z z_m
	}
	else if "`model'"=="fiv" {
		local list_local `list_local' z z_m dh dh_m
	}
	if "`clustvar'"~=""		local list_local `list_local' clustvar
	foreach obj in `list_local' {
		mata: `A'.put(("`obj'","local"),"``obj''")
	}
	// store scalars
	local list_scalar
	if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
	foreach obj in `list_scalar' {
		mata: `A'.put(("`obj'","scalar"),``obj'')
	}
	
	// additional estimation results
	ereturn scalar resample = `rep'
	tempname eqn
	mata: `eqn' = init_eStruct()
	// Y eqn results
	mata: `eqn' = (`mname'.eqnAA).get("`yname'")
	// MSE
	mata: `A'.put(("`y'_mse","scalar"),return_result_item(`eqn',"`y'","MSE","`rep'"))
	// MSE folds
	mata: `A'.put(("`y'_mse_folds","matrix"),return_result_item(`eqn',"`y'","MSE_folds","`rep'"))
	// pystacked final est (pystacked multi only)
	mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
	if `pystackedmulti' {
		// cap because won't exist for e.g. shortstack variable
		cap mata: `A'.put(("`y'_stack_final_est","local"), return_learner_item(`eqn',"`y'","stack_final_est"))
	}
	// ss results
	if "`spec'"=="ss" {
		mata: st_local("shortstack", `eqn'.shortstack)
		mata: `A'.put(("`yname'_ssw","matrix"), return_result_item(`eqn',"`shortstack'_ss","ss_weights","`rep'"))
		mata: `A'.put(("`yname'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
	}
	// ps results
	if "`spec'"=="ps" {
		mata: st_local("poolstack", `eqn'.poolstack)
		mata: `A'.put(("`yname'_psw","matrix"), return_result_item(`eqn',"`poolstack'_ps","ps_weights","`rep'"))
		mata: `A'.put(("`yname'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
	}

	// D eqn results - uses vtilde names in d
	forvalues i=1/`numeqnD' {
		local dname : word `i' of `dnames'
		local vtilde : word `i' of `d'
		local vtilde_h `vtilde'

		mata: `eqn' = (`mname'.eqnAA).get("`dname'")
		// MSE
		mata: `A'.put(("`vtilde'_mse","scalar"),return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
		// MSE folds
		mata: `A'.put(("`vtilde'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
		// pystacked final est (pystacked multi only)
		mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
		if `pystackedmulti' {
			// cap because won't exist for e.g. shortstack variable
			cap mata: `A'.put(("`vtilde'_stack_final_est","local"), return_learner_item(`eqn',"`vtilde'","stack_final_est"))
		}
		// ss results
		if "`spec'"=="ss" {
			mata: st_local("shortstack", `eqn'.shortstack)
			mata: `A'.put(("`dname'_ssw","matrix"), return_result_item(`eqn',"`shortstack'_ss","ss_weights","`rep'"))
			mata: `A'.put(("`dname'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
		}
		// ps results
		if "`spec'"=="ps" {
			mata: st_local("poolstack", `eqn'.poolstack)
			mata: `A'.put(("`dname'_psw","matrix"), return_result_item(`eqn',"`poolstack'_ps","ps_weights","`rep'"))
			mata: `A'.put(("`dname'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
		}
		if `fivflag' {
			// MSE
			mata: `A'.put(("`vtilde_h'_mse","scalar"),return_result_item(`eqn',"`vtilde_h'","MSE_h","`rep'"))
			// MSE folds
			mata: `A'.put(("`vtilde_h'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde_h'","MSE_h_folds","`rep'"))
			// ss weights (no h final est saved - same as main est)
			if "`spec'"=="ss" {
					mata: st_local("shortstack", `eqn'.shortstack)
					mata: `A'.put(("`dname'_h_ssw","matrix"), return_result_item(`eqn',"`shortstack'_ss","ss_weights_h","`rep'"))
			}
			// ps weights (no h final est saved - same as main est)
			if "`spec'"=="ps" {
					mata: st_local("poolstack", `eqn'.poolstack)
					mata: `A'.put(("`dname'_h_psw","matrix"), return_result_item(`eqn',"`poolstack'_ps","ps_weights_h","`rep'"))
			}
		}
	}
	if `fivflag'==0 {
		// Z eqn results; fiv won't enter
		forvalues i=1/`numeqnZ' {
			local zname : word `i' of `znames'
			local vtilde : word `i' of `z'
			mata: `eqn' = (`mname'.eqnAA).get("`zname'")
			// MSE
			mata: `A'.put(("`vtilde'_mse","scalar"),return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
			// MSE folds
			mata: `A'.put(("`vtilde'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
			// pystacked final est (pystacked multi only)
			mata: st_local("pystackedmulti", strofreal(`eqn'.pystackedmulti))
			if `pystackedmulti' {
				// cap because won't exist for e.g. shortstack variable
				cap mata: `A'.put(("`vtilde'_stack_final_est","local"), return_learner_item(`eqn',"`vtilde'","stack_final_est"))
			}
			// ss results
			if "`spec'"=="ss" {
				mata: st_local("shortstack", `eqn'.shortstack)
				mata: `A'.put(("`zname'_ssw","matrix"), return_result_item(`eqn',"`shortstack'_ss","ss_weights","`rep'"))
				mata: `A'.put(("`zname'_ss_final_est","local"), return_learner_item(`eqn',"`shortstack'_ss","ss_final_est"))
			}
			// ps results
			if "`spec'"=="ps" {
				mata: st_local("poolstack", `eqn'.poolstack)
				mata: `A'.put(("`zname'_psw","matrix"), return_result_item(`eqn',"`poolstack'_ps","ps_weights","`rep'"))
				mata: `A'.put(("`zname'_ps_final_est","local"), return_learner_item(`eqn',"`poolstack'_ps","ps_final_est"))
			}
		}
	}
	
	mata: (`mname'.estAA).put(("`spec'","`rep'"),`A')
	
	// no longer needed
	foreach obj in `A' `eqn' {
		cap mata: mata drop `obj'
	}
	
end

// estimates and stores mean/median estimates across resamples
program define medmean_and_store, eclass
	syntax [anything] [if] [in] , [								///
				mname(name)										///
				spec(string)									///
				title(string)									///
				medmean(string)									///
				NOConstant										///
				]

	local consflag = ("`noconstant'"=="")
	mata: st_local("model",`mname'.model)
	local ivflag	= "`model'"=="iv"
	local fivflag	= "`model'"=="fiv"
		
	tempname b V bagg Vagg Vi
	tempname bvec brow sbvec bmed Vvec sVvec Vmed
	tempvar esample
	tempname B
	
	// initialize
	mata: st_local("nameD",invtokens(`mname'.nameD))
	// don't aggregate the constant
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
		mata: `brow' = `B'.get(("b","post"))
		// don't aggregate the constant
		mata: `bvec'[`m',.] = `brow'[1,(1..(cols(`brow')-`consflag'))]
		// row/colnames etc. - need to do this only once
		if `m'==1 {
			mata: st_local("depvar",`B'.get(("depvar","post")))
			// retrieve locals; if empty, will be ""
			local list_local y d dh z yname dnames znames vce vcetype clustvar
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
			mata: `sbvec' = sort(`bvec',`k')
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
			// don't aggregate the constant
			mata: `Vi' = `Vi'[(1::(cols(`Vi')-`consflag')),(1..(cols(`Vi')-`consflag'))]
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
				mata: `sVvec' = sort(`Vvec',1)
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

	tempname A
	mata: `A' = AssociativeArray()
	mata: `A'.reinit("string",2)
	mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
	
	mata: `A'.put(("N","post"),`N')
	mata: `A'.put(("b","post"),`bagg')
	mata: `A'.put(("V","post"),`Vagg')
	mata: `A'.put(("depvar","post"),"`depvar'")
	
	// store locals
	local list_local title y d yname dnames vce vcetype
	if "`model'"=="iv" {
		local list_local `list_local' z
	}
	else if "`model'"=="fiv" {
		local list_local `list_local' z dh
	}
	if "`clustvar'"~=""		local list_local `list_local' clustvar
	foreach obj in `list_local' {
		mata: `A'.put(("`obj'","local"),"``obj''")
	}
	// store scalars
	local list_scalar
	if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
	foreach obj in `list_scalar' {
		mata: `A'.put(("`obj'","scalar"),``obj'')
	}
	// special case - "_m" subscript doesn't apply to mean/median over resamplings
	// so store without resample subscript
	foreach obj in y d dh z {
		mata: `A'.put(("`obj'_m","local"),"``obj''")
	}
	// special case - vector of betas available only for mean/median
	mata: `A'.put(("b_resamples","matrix"),`bvec')
	
	// additional estimation results
	local numeqnD	: word count `dnames'
	local numeqnZ	: word count `znames'
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
	
	// store AA with median/mean results
	mata: (`mname'.estAA).put(("`spec'","`medmean'"),`A')
	
	// no longer needed
	foreach obj in `A' `B' `bagg' `bvec' `brow' `sbvec' `Vagg' `Vvec' `sVvec' `Vi' {
		cap mata: mata drop `obj'
	}
	
end

// replays stored estimate
program define replay_estimate, eclass
	syntax [anything] [if] [in] , [								///
				mname(name)										///
				spec(string)									///
				rep(string)										///
				NOConstant										///
				]

	local consflag = ("`noconstant'"=="")
	mata: st_local("model",`mname'.model)
	local ivflag	= "`model'"=="iv"
	local fivflag	= "`model'"=="fiv"

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
	mata: st_local("yname",`B'.get(("yname","local")))
	mata: st_local("dnames",`B'.get(("dnames","local")))
	
	if `consflag' {
		// will be empty if no constant
		local consname "_cons"
	}
	
	matrix rownames `b' = `depvar'
	matrix colnames `b' = `dnames' `consname'
 	matrix colnames `V' = `dnames' `consname'
	matrix rownames `V' = `dnames' `consname'
	
	tempvar esample
	cap gen `esample' = `mname'_sample_`rep'
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
	ereturn local tmname `mname'	// temporary mname
	
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
	di as text "y-E[y|X]" _col(11) "= " as res "y-`e(y_m)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	if "`e(model)'"~="fiv" {
		di as text "D-" _c
	}
	di as text "E[D|X,Z]" _col(11) "= " _c
	local numeqnD : word count `e(d_m)'
	forvalues i=1/`numeqnD' {
		local Dtilde : word `i' of `e(d_m)' {
		// iv model residualizes, fiv model uses different approach
		if "`e(model)'"=="iv" di as res "D`i'-" _c
		di as res "`Dtilde' " _c
	}
	di
	if "`e(model)'" == "iv" {
		di as text "Z-E[Z|X]" _col(11) "= " _c
		local numeqnZ : word count `e(z_m)'
		forvalues i=1/`numeqnZ' {
			local Ztilde : word `i' of `e(z_m)' {
			di as res "Z`i'-`Ztilde' " _c
		}
		di
	}
	if "`e(model)'" == "fiv" {
		di as text "E[D|X]" _col(11) "= " as res "`e(dh_m)'"
		di as text "Orthogonalized D = D - E[D|X]; optimal IV = E[D|X,Z] - E[D|X]."
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

end
