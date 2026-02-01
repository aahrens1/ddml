*! ddml v1.5.0
*! last edited: 27jan2026
*! authors: aa/ms
* need to accommodate weights in parsing of estimation strings

program define crossfit, rclass sortpreserve
	version 16
	syntax [anything] [if] [in] ,						/// [anything] is renamed to vname below; currently undocumented option
							[							///
							estring(string asis)		/// estimation string
							estringh(string asis)		/// est string for E[D^|XZ]
														/// need asis option in case it includes strings
							ename(name)					/// name for Mata struct; default is "crossfit"
							pystackedmulti(integer 1)	/// if Mata eStruct provided, use pystacked-specific code if possible
							NOREPLACE					/// use already-initialized eqn in ename
							vtilde(namelist)			/// name(s) of fitted variable(s)
							Generate(namelist)			/// synonym for vtilde
							vtype(string)				/// datatype of fitted variable; default=double
							shortstack(name)			/// for interactive use
							poolstack(name)				/// for interactive use
							model(name)					/// when called by _ddml_crossfit
							crossfitother				/// force use of general (not-pystacked-specific) code
							predopt(string asis)		/// undocumented
							NOIsily						///
							*							/// other options that go to subroutines
							]

	// renaming for clarity
	local vtlist `vtilde' `generate'
	local vname `anything'
	// clear the local macros
	local vtilde
	local generate
	// used in call to subroutines
	local cfoptions `options'

	** indicator for initializing eqn struct
	local initflag	= "`noreplace'"==""

	if `initflag' {
		// for interactive use; will be only one equation (one learner)
		// start with a blank setup
		
		*** set up struct with equation/learner info
		
		if "`ename'"== "" {
			local eqn_info crossfit
		}
		else {
			local eqn_info `ename'
		}
		
		*** variable type
		if "`vtype'"==""		local vtype double
		if "`vtype'"=="none"	local vtype
	
		mata: `eqn_info' = init_eStruct()
		initialize_eqn_info,										///
							ename(`eqn_info')						///
							vname(`vname')							///
							vtlist(`vtlist')						///
							shortstack(`shortstack')				///
							poolstack(`poolstack')					///
							estring(`estring')						///
							estringh(`estringh')					///
							vtype(`vtype')							///
							predopt(`predopt')					 	///
							`noisily'
		local nlearners = r(nlearners)
		if "`vtlist'"=="" {
			// default names set by subroutine
			mata: st_local("vtlist", invtokens(`eqn_info'.vtlist))
		}
		if "`vname'"=="" {
			// depvar identified by subroutine
			mata: st_local("vname", `eqn_info'.vname)
		}
		// identify if pystacked with multiple learners
		local vtilde : word 1 of `vtlist'
		mata: st_local("est_main", return_learner_item(`eqn_info',"`vtilde'","est_main"))
		mata: st_local("est_options", return_learner_item(`eqn_info',"`vtilde'","est_options"))
		local cmd : word 1 of `est_main'
		if "`cmd'"=="pystacked" {
			// pystacked as single learner?
			tempname holdname
			_estimates hold `holdname', nullok
			`qui' di as text "calling pystacked on full sample with noestimate option..."
			`qui' `est_main' , `est_options' noestimate
			if _rc==0 {
				`qui' di as text "N=" as res e(N)
				`qui' di as text "number of learners = " as res `e(mcount)'
				local mcount = e(mcount)
				`qui' di as text "Base learners: " _c
				forvalues j=1/`mcount' {
					di as res e(method`j') " " _c
				}
				`qui' di
			}
			_estimates unhold `holdname'
		}
		else {
			local mcount 1
		}
		// pystackedmulti default = 0; update if multiple pystacked learners
		if "`cmd'"=="pystacked" & `mcount'>1 {
			mata: `eqn_info'.pystackedmulti = `mcount'
			local pystackedmulti = `mcount'
		}
		// set pystackedmulti flag (overwrite default)
		mata: st_local("pystackedmulti", strofreal(`eqn_info'.pystackedmulti))
	}
	else {
		// ename is a pre-populated eqn struct
		local eqn_info `ename'
		mata: st_local("ssflag", strofreal(`ename'.shortstack~=""))
		mata: st_local("psflag", strofreal(`ename'.poolstack~=""))
		// default is to use pystacked-specific code unless user says "no" or doesn't apply
		if `pystackedmulti' {
			// user says "yes if possible" so set to whatever eqn says
			mata: st_local("pystackedmulti", strofreal(`eqn_info'.pystackedmulti))
		}
	}

	if (`pystackedmulti' | `psflag') & ("`crossfitother'"=="") {
		// crossfitother forces branching to _crossfit_other
		// enter single learner = pystacked and number of pystacked base learners > 1
		_crossfit_pystacked	`if' `in',	///
			ename(`eqn_info')			///
			`noisily'					///
			`cfoptions'
	}
	else {
		// special handling for single pystacked multilearner with fiv
		// supported but no code longer branches here
		local stdflag=(`pystackedmulti'>1)
		// multiple ddml learners or single learner = pystacked with a single base learner
		_crossfit_other	`if' `in',		///
			ename(`eqn_info')			///
			`noisily'					///
			`cfoptions'					///
			stdflag(`stdflag')
	}
	
	return add

end


program define _crossfit_pystacked, rclass sortpreserve
	version 16
	syntax [if] [in] ,								///
							[ estring(string asis)	/// estimation string
													/// need asis option in case it includes strings
													///
							ename(name)				/// name for Mata struct; default is "crossfit"
							foldvar(varlist)		/// one fold var per resample
							kfolds(integer 5)		/// ignored if foldvars provided
							reps(integer 1)			/// ignored if foldvars provided
							firstrep(integer 1)		/// add reps starting with m=firstrep
							NORANDOM				/// first fold ID uses obs in existing order
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							NOENFORCElie			/// do not enforce LIE in fiv model
							NOIsily					///
							stdfinalest(name)		/// final estimator for standard stacking
							ssfinalest(name)		/// final estimator for short-stacking
							psfinalest(name)		/// final estimator for pooled-stacking
							finalest(name)			/// final estimator for all
							NOSTDstack				/// no standard stacking - use psytacked+voting to get learners only
							CVCbootnum(integer 500)	///
							*						/// ignored options
							]

	// syntax errors
	if "`options'"~="" {
		di as err "crossfit error: illegal option '`options''"
		exit 197
	}

	*** set sample
	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `vname'
	
	// used throughout	
	local cmd pystacked

	mata: st_local("vname", `ename'.vname)
	mata: st_local("vtlist", invtokens(`ename'.vtlist))
	mata: st_local("shortstack", `ename'.shortstack)
	mata: st_local("poolstack", `ename'.poolstack)
	mata: st_local("etype", `ename'.etype)	// "Y", "D", "DH" or "Z"

	// 3 cases: basic partially linear, interactive, fiv
	if "`etype'"=="DH"			local est_type "fiv"
	else if "`treatvar'"~=""	local est_type "interactive"
	else						local est_type "partial"
	
	** indicator for enforcing LIE in fiv model **
	local enforceflag	= "`noenforcelie'"==""
	
	** indicator for LATE model special case: perfect assignment to treatment
	// nb: D (treatment, `vname') can be continuous; code below aimed at binary treatment
	//	 Z (assignment, `treatvar') is always binary
	// (note bad naming convention)
	// perfect assignment: no one treated if not assigned
	// means that in the treatvar Z=0 subsample of the estimation sample, no one is assigned to treatment (vname D=0)
	// since no variation in D, estimation fails, so catch here as a special case
	// and set predicted values vhat of D to 0
	// applies to full sample rather than fold-by-fold
	if "`est_type'"=="interactive" & "`etype'"=="D" {
		qui count if `treatvar'==0 & `touse'
		local Z0_count = r(N)
		qui count if `treatvar'==0 & `vname'==0 & `touse'
		local D0Z0_count = r(N)
		local perfectflag = (`Z0_count'==`D0Z0_count')
	}
	else {
		local perfectflag = 0
	}

	** indicator for short-stacking
	local ssflag		= "`shortstack'"~=""
	** indicator for pooled-stacking
	local psflag		= "`poolstack'"~=""
	** indicator for standard stacking
	local stdflag		= "`nostdstack'"==""
	if `psflag' & ~`stdflag' {
		di as res "pooled stacking requires standard stacking; poolstack option ignored"
		local psflag = 0
	}
	if ~`stdflag' & ~`ssflag' {
		di as err "error - pystacked integration requires using either standard and/or short-stacking"
		exit 198
	}
	
	// unless specified, finalest sets the final estimator for all stacking methods
	// if empty, finalest will be the pystacked/ddml default
	if "`stdfinalest'"==""	local stdfinalest `finalest'
	// standard stacking final estimator goes directly to pystacked as an option
	// but only if not empty (otherwise it can conflict with the pystacked finalest(.) option)
	if "`stdfinalest'"~=""	local stdfinalestopt finalest(`stdfinalest')
	// nostd => pystacked option is voting instead of stacking
	// votetype is ignored if type=reg; relevant only for type=class
	// also reset other std final est macros to null
	if ~`stdflag' {
		local nostdstackopt voting votetype(soft)
		local stdfinalest
		local stdfinalestopt
	}

	if "`noisily'"=="" {
		local qui quietly
		local cap cap
	}
	
	*** debugging message
	`qui' di as text "entering _crossfit_pystacked..."
	
	*** set up fold/reps variables
	
	// fold id variables
	if "`foldvar'"~="" {
		// fold variables provided
		foreach fvar of varlist `foldvar' {
			tempvar fid
			qui egen `fid' = group(`fvar')
			local fidlist `fidlist' `fid'
		}
		check_foldvar, fidlist(`fidlist') touse(`touse')
		local kfolds = r(kfolds)
		local reps = r(reps)
		local lastrep = `firstrep' + `reps' - 1
	}
	else {
		// fold variables not provided, so generate
		local lastrep = `firstrep' + `reps' - 1
		forvalues m=`firstrep'/`lastrep' {
			*** gen folds
			tempvar uni cuni fid
			if `m'==1 & "`norandom'"~="" {
				qui gen `uni' = _n
			}
			else {
				qui gen double `uni' = runiform() if `touse'
			}
			qui cumul `uni' if `touse', gen(`cuni')
			qui gen int `fid' = ceil(`kfolds'*`cuni') if `touse'
			local fidlist `fidlist' `fid'
		}
	}
	
	*** initialization
	// only one learner, pystacked
	local vtilde `vtlist'
	mata: st_local("est_main", return_learner_item(`ename',"`vtilde'","est_main"))
	mata: st_local("est_options", return_learner_item(`ename',"`vtilde'","est_options"))
	mata: st_local("predopt",return_learner_item(`ename',"`vtilde'","predopt"))
	mata: st_local("vtype",return_learner_item(`ename',"`vtilde'","vtype"))
	if "`est_type'"=="fiv" {
		// FIV locals
		mata: st_local("est_main_h", return_learner_item(`ename',"`vtilde'","est_main_h"))
		mata: st_local("est_options_h", return_learner_item(`ename',"`vtilde'","est_options_h"))
		mata: st_local("predopt_h",return_learner_item(`ename',"`vtilde'","predopt_h"))
	}

	// call pystacked with noestimate option to parse and get basic model specs
	`qui' di as text "calling pystacked on full sample with noestimate option..."
	`qui' `est_main' if `touse', `est_options' noestimate
	`qui' di as text "N=" as res e(N)
	`qui' di as text "number of learners = `e(mcount)'"
	local cmd `e(cmd)'
	// holds for all reps and folds
	local nlearners	= e(mcount)
	// will be "reg" or "class"
	local stype `e(type)'
	if "`est_type'"=="fiv" {
		// check h estimation and initialize cmd_h local
		// if h estimation is also pystacked, number of learners is =1 or =nlearners
		local cmd_h : word 1 of `est_main_h'		
		if "`cmd_h'"=="pystacked" {
			// temporarily replace {D} placeholder in order to call pystacked w/noestimate
			local est_main_h_st = subinstr("`est_main_h'","{D}","`vname'",1)
			`qui' di as text "calling pystacked h cmd on full sample with noestimate option..."
			`qui' `est_main_h_st' if `touse', `est_options' noestimate
			local nlearners_h = e(mcount)
			if `nlearners_h'~=1 & `nlearners_h'~=`nlearners' {
				di as err "error - number of E[D|X] learners must be either 1 or match number of E[D|Z,X] learners"
				exit 198
			}
		}
		else {
			local nlearners_h = 1
		}
	}
	
	// stored results; initialize here
	tempname mse_list N_list mse_folds_list N_folds_list
	tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
	tempname mse0_list N0_list mse0_folds_list N0_folds_list
	tempname mse1_list N1_list mse1_folds_list N1_folds_list
	// more temps
	tempname pysw pysw0 pysw1 pysw_h
	tempname pysw_temp pysw1_temp pysw0_temp pysw_h_temp
	tempname pysm pysm0 pysm1 pysm_h
	tempname pysm_temp pysm1_temp pysm0_temp pysm_h_temp
	
	// will save base learner CV predictions in mata
	tempname y_stacking_cv y_stacking_cv0 y_stacking_cv1 y_stacking_cv_h y_stacking_cv_h_ss

	tempvar tomata fidtouse fidtouse1 fidtouse0
	tempvar vhat vhat0 vhat1 vhat_not_k hhat vres hres vres_sq vres0_sq vres1_sq hres_sq
	tempvar vhat_ss vhat_not_k_ss hhat_ss
	qui gen byte `tomata'=.
	qui gen int `fidtouse'=.
	qui gen int `fidtouse0'=.
	qui gen int `fidtouse1'=.
	qui gen `vtype' `vhat'=.
	qui gen `vtype' `vhat0'=.
	qui gen `vtype' `vhat1'=.
	qui gen `vtype' `vhat_not_k'=.
	qui gen `vtype' `hhat'=.
	qui gen double `vres'=.
	qui gen double `hres'=.
	qui gen double `vres_sq'=.
	qui gen double `vres0_sq'=.
	qui gen double `vres1_sq'=.
	qui gen double `hres_sq'=.
	qui gen `vtype' `vhat_ss'=.
	qui gen `vtype' `vhat_not_k_ss'=.
	qui gen `vtype' `hhat_ss'=.

	// needed only for fiv - temp vars for in-sample predicted values of individual learners
	if "`est_type'"=="fiv" {
		forvalues j=1/`nlearners' {
			tempvar vt_not_k_L`j'
			qui gen double `vt_not_k_L`j''=.
			local vt_not_k_list `vt_not_k_list' `vt_not_k_L`j''
		}
	}
	
	
	// notes:
	// `vtilde': name for fitted variable (predicted value)
	// `vhat_k': OOS (crossfit) stacked predicted values for fold k.
	// `vhat1_k', `vhat0_k': as above for TV.
	// `vhat': OOS (crossfit) stacked predicted values, for all folds; filled when looping over folds. tempvar.
	// `vtilde'_`m': `vhat' using user-provided varname.
	// `vtilde'_L`j'_`m': OOS (crossfit) predicted values, by learner. m is resample j is base learner number.
	// `vtilde_list': list of `vtilde'_L`j'_`m' for j=1,2,...; reset for each m loop.
	// `stacking_p_cv'[1,2,...]: in-sample CV predicted values of base learners; saved in mata for poolstacking.
	// `stacking_p'[1,2,...]: predicted values of base learners, in-sample and OOS; OOS => `vtilde'_L`j'_`m'.
	// `shortstack': name for shortstacked fitted variable (predicted value).
	// `shortstack'_ss_`m': shortstacked fitted variable of D for resample m; not a tempvar.
	// `y_stacking_cv': mata matrix with depvar in column 1 and in-sample CV base learner predicted values in rest
	//  accumulated in mata since row dimension >> rows of dataset; 0/1 variants for ate/late
	
	// loop over resamples (crossfitting, shortstacking, store results)
	forvalues m=`firstrep'/`lastrep' {
	
		// create rnames for renaming the rows of saved matrices
		local rnames `rnames' resample_`m'
		
		// initialize for each resample
		local failed
		local failed1
		local failed0
		local failed_h
		cap mat drop `pysw'
		cap mat drop `pysw_h'
		cap mat drop `pysw0'
		cap mat drop `pysw1'
		cap mat drop `pysm'
		cap mat drop `pysm_h'
		cap mat drop `pysm0'
		cap mat drop `pysm1'
		
		
		// create list of names of variables created by learners for rep m
		// need to reset for each m
		if "`est_type'"=="partial" { // case 1 - partially-linear
			local vtilde_list
			forvalues j=1/`nlearners' {
				local vtilde_list `vtilde_list' `vtilde'_L`j'_`m'
			}
		}
		else if "`est_type'"=="interactive" {	// case 2 - interactive
			local vtilde0_list
			local vtilde1_list
			forvalues j=1/`nlearners' {
				local vtilde0_list `vtilde0_list' `vtilde'0_L`j'_`m'
				local vtilde1_list `vtilde1_list' `vtilde'1_L`j'_`m'
			}
		}
		else if "`est_type'"=="fiv" { // case 3 - FIV
			local vtilde_list
			local vtilde_h_list
			forvalues j=1/`nlearners' {
				local vtilde_list   `vtilde_list'   `vtilde'_L`j'_`m'
				local vtilde_h_list `vtilde_h_list' `vtilde'_h_L`j'_`m'
			}
		}
		else {
			di as err "internal crossfit error - unknown est type `est_type'"
			exit 198
		}
			
		******************************** CROSSFITTING ************************************
		
		// cross-fitting fundamentally depends on three cases: plm, interactive, and fiv
		
		if `lastrep'>1 {
			di as text "Resample `m'..."
		}

		// list of fold IDs refers to current set of reps (not full set including appended)
		local thisrep = `m'-`firstrep'+1
		local fid : word `thisrep' of `fidlist'

		// create blank fitted variable(s), other initializations for resample m
		if "`est_type'"=="partial" { // case 1 - partially-linear
			qui replace `vhat'=.
			// base learner predicted values
			// not tempvars; name based on `vtilde' macro; used for poolstacking
			forvalues j=1/`nlearners' {
				cap drop `vtilde'_L`j'_`m'
				qui gen `vtype' `vtilde'_L`j'_`m' = .
			}
		}  
		else if "`est_type'"=="interactive" {	// case 2 - interactive
			qui replace `vhat0'=.
			qui replace `vhat1'=.
			// base learner predicted values
			// not tempvars; name based on `vtilde' macro; used for poolstacking
			forvalues j=1/`nlearners' {
				// tv=0
				cap drop `vtilde'0_L`j'_`m'
				qui gen `vtype' `vtilde'0_L`j'_`m' = .
				// tv=1
				cap drop `vtilde'1_L`j'_`m'
				qui gen `vtype' `vtilde'1_L`j'_`m' = .
			}
		}
		else if "`est_type'"=="fiv" {	// case 3 - FIV
			qui replace `vhat'=.
			qui replace `vhat_not_k'=.
			qui replace `hhat'=.
			qui replace `vhat_ss'=.			// short-stacked version
			qui replace `vhat_not_k_ss'=.	// short-stacked version
			qui replace `hhat_ss'=.			// short-stacked version
			// base learner predicted values
			// not tempvars; name based on `vtilde' macro; used for poolstacking
			forvalues j=1/`nlearners' {
				cap drop `vtilde'_L`j'_`m'
				qui gen `vtype' `vtilde'_L`j'_`m' = .
				cap drop `vtilde'_h_L`j'_`m'
				qui gen `vtype' `vtilde'_h_L`j'_`m' = .
			}
		}
		else {
			di as err "internal crossfit error"
			exit 198
		}

		
		// crossfit
		di as text "Cross-fitting fold " _c
		forvalues k = 1(1)`kfolds' {
	
			// initialize
			local notyetfailed = 1
			di as text "`k' " _c
			
			if "`est_type'"=="partial" { // case 1 - partially-linear
				
				if (`k'==1) & (`stdflag' | `nlearners'==1) {
					// initialize mata object to hold dep var and in-sample crossfit CV predictions
					// also if just one learner (so stacking=single learner)
					mata: `y_stacking_cv' = J(0,`nlearners'+3,0)
				}
			
				// estimate excluding kth fold
				_fit_pystacked,	k(`k')					///
					m(`m')								///
					nlearners(`nlearners')				///
					stdflag(`stdflag')					///
					cap(`cap')							///
					qui(`qui')							///
					est_main(`est_main')				///
					est_options(`est_options')			///
					stdfinalestopt(`stdfinalestopt')	///
					nostdstackopt(`nostdstackopt')		///
					vname(`vname')						///
					vhat(`vhat')						///
					vtilde(`vtilde_list')				///
					fid(`fid')							///
					touse(`touse')						///
					mmatname(`y_stacking_cv')			///
					etype(`etype')

				if `r(failed)'==0 {
					// no failed estimation
					local base_est			`r(base_est)'
					local stack_final_est	`r(final_est)'
					mat `pysw'				= (nullmat(`pysw'),r(pysw))
					mat `pysm'				= (nullmat(`pysm'),r(pysm))
				}
				else {
					// failed means estimation failed, no results to save
					// track so that pysw and psym can have missings inserted at end
					local failed			`failed' `k'
					if `notyetfailed'		di
					local notyetfailed = 0
					di as res "Warning: estimation failed with return code `r(rc)' when cross-fitting fold `k'."
				}
			}	// end case 1 - partially linear

			else if "`est_type'"=="interactive" {	// case 2 - interactive
	
				// pystacked learners
				if (`k'==1) & (`stdflag' | `nlearners'==1) {
					// initialize mata object to hold dep var and in-sample crossfit CV predictions
					mata: `y_stacking_cv0' = J(0,`nlearners'+3,0)
					mata: `y_stacking_cv1' = J(0,`nlearners'+3,0)
				}
	
				// outcome equation so estimate separately for treated and untreated

				// for treatvar = 1
				_fit_pystacked,	k(`k')					///
					m(`m')								///
					nlearners(`nlearners')				///
					stdflag(`stdflag')					///
					cap(`cap')							///
					qui(`qui')							///
					est_main(`est_main')				///
					est_options(`est_options')			///
					nostdstackopt(`nostdstackopt')		///
					stdfinalestopt(`stdfinalestopt')	///
					vname(`vname')						///
					vhat(`vhat1')						///
					vtilde(`vtilde1_list')				///
					treatvar(`treatvar')				///
					treatval(1)							///
					fid(`fid')							///
					touse(`touse')						///
					mmatname(`y_stacking_cv1')			///
					etype(`etype')

				if `r(failed)'==0 & `r(perfect)'==0 {
					// no failed estimation
					// no perfect assignment
					local base_est			`r(base_est)'
					local stack_final_est	`r(final_est)'
					mat `pysw1'				= (nullmat(`pysw1'),r(pysw))
					mat `pysm1'				= (nullmat(`pysm1'),r(pysm))
				}
				else {
					// failed means estimation failed, no results to save
					// perfect means perfect assignment to treatment (no variation in D)
					// estimation would have failed so we don't estimation, no results to save
					// track so that pysw and psym can have missings inserted at end
					local failed1			`failed1' `k'
					if `notyetfailed'		di
					local notyetfailed = 0
					if `r(failed)' {
						di as res "Warning: estimation failed with return code `r(rc)' when cross-fitting fold `k'."
					}
					else {
						di as res "Warning: perfect assignment to treatment when cross-fitting fold `k'."
					}
				}

				// for treatvar = 0
				_fit_pystacked,	k(`k')					///
					m(`m')								///
					nlearners(`nlearners')				///
					stdflag(`stdflag')					///
					cap(`cap')							///
					qui(`qui')							///
					est_main(`est_main')				///
					est_options(`est_options')			///
					nostdstackopt(`nostdstackopt')		///
					stdfinalestopt(`stdfinalestopt')	///
					vname(`vname')						///
					vhat(`vhat0')						///
					vtilde(`vtilde0_list')				///
					treatvar(`treatvar')				///
					treatval(0)							///
					perfect(`perfectflag')				/// special LATE case - perfect assignment
					fid(`fid')							///
					touse(`touse')						///
					mmatname(`y_stacking_cv0')			///
					etype(`etype')
// check
				if `r(failed)'==0 /* & `r(perfect)'==0 */ {
					// no failed estimation
					// no perfect assignment
					// already set in treat=1 case:
					// local base_est			`r(base_est)'
					// local stack_final_est	`r(final_est)'
					mat `pysw0'				= (nullmat(`pysw0'),r(pysw))
					mat `pysm0'				= (nullmat(`pysm0'),r(pysm))
				}
				else {
					// failed means estimation failed, no results to save
					// perfect means perfect assignment to treatment (no variation in D)
					// estimation would have failed so we don't estimate, no results to save
					// track so that pysw and psym can have missings inserted at end
					local failed0			`failed0' `k'
					if `notyetfailed'		di
					local notyetfailed = 0
					if `r(failed)' {
						di as res "Warning: estimation failed with return code `r(rc)' when cross-fitting fold `k'."
					}
					else {
						di as res "Warning: perfect assignment to non-treatment when cross-fitting fold `k'."
					}
				}
			}	// end case 2 - interactive
			
			else if "`est_type'"=="fiv" { // case 3 - FIV

				if (`k'==1) & (`stdflag' | `nlearners'==1) {
					// initialize mata object to hold dep var and in-sample crossfit CV predictions
					mata: `y_stacking_cv'		= J(0,`nlearners'+3,0)
					mata: `y_stacking_cv_h'		= J(0,`nlearners'+3,0)
					mata: `y_stacking_cv_h_ss'	= J(0,`nlearners'+3,0)		// short-stacked version
				}
				
				// initialize variable to hold in-sample predicted values
				qui replace `vhat_not_k'=.
				qui replace `vhat_not_k_ss'=.		// shortstacked version

				// Step I: estimation of E[D|Z,X]=D^
				// estimate excluding kth fold

				_fit_pystacked,	k(`k')					///
					m(`m')								///
					nlearners(`nlearners')				///
					stdflag(`stdflag')					///
					cap(`cap')							///
					qui(`qui')							///
					est_main(`est_main')				///
					est_options(`est_options')			///
					stdfinalestopt(`stdfinalestopt')	///
					nostdstackopt(`nostdstackopt')		///
					vname(`vname')						///
					vhat(`vhat')						/// where OOS fitted values are collected
					vhat_ss(`vhat_ss')					/// OOS fitted values, needed for short-stacking final nnls
					vhat_not_k(`vhat_not_k')			/// in-sample fitted values
					vhat_not_k_ss(`vhat_not_k_ss')		/// in-sample fitted values, short-stacked version
					vtilde(`vtilde_list')				/// OOS (crossfit) base learner predicted values
					vtilde_not_k(`vt_not_k_list')		/// in-sample base learner predicted values
					fid(`fid')							///
					touse(`touse')						///
					mmatname(`y_stacking_cv')			///
					etype(`etype')

				if `r(failed)'==0 {
					// no failed estimation
					local base_est			`r(base_est)'
					local stack_final_est	`r(final_est)'
					mat `pysw'				= (nullmat(`pysw'),r(pysw))
					mat `pysm'				= (nullmat(`pysm'),r(pysm))
				}
				else {
					// failed means estimation failed, no results to save
					// track so that pysw and psym can have missings inserted at end
					local failed			`failed' `k'
					if `notyetfailed'		di
					local notyetfailed = 0
					di as res "Warning: estimation failed with return code `r(rc)' when cross-fitting fold `k'."
				}

				// Step II: estimation of E[D|X]

				// E[D|X] learner may not be pystacked (e.g., if no Xs)

				if "`cmd_h'"=="pystacked" & (`nlearners'==`nlearners_h') {

					if `enforceflag' {
						// replace {D}-placeholder in estimation string with variable name
						local est_main_h_st = subinstr("`est_main_h'","{D}","`vhat_not_k'",1)
					}
					else {
						// replace {D}-placeholder with varname of original D
						local est_main_h_st = subinstr("`est_main_h'","{D}","`vname'",1)
					}

					if `stdflag' | `nlearners'==1 {

						_fit_pystacked,							///
							k(`k')								///
							m(`m')								///
							nlearners(`nlearners')				///
							stdflag(`stdflag')					///
							cap(`cap')							///
							qui(`qui')							///
							est_main(`est_main_h_st')			///
							est_options(`est_options_h')		///
							stdfinalestopt(`stdfinalestopt')	///
							nostdstackopt(`nostdstackopt')		///
							vname(`vhat_not_k')					///
							vhat(`hhat')						/// where to save (replace) stacked OOS (crossfit) predicted values
							vtilde(`vtilde_h_list')				/// where to save (replace) OOS predicted values for indiv learners
							fid(`fid')							///
							touse(`touse')						///
							mmatname(`y_stacking_cv_h')			///
							etype(`etype')
						if `r(failed)'==0 {
							// no failed estimation
							local base_est_h		`r(base_est)'
							local stack_final_est_h	`r(final_est)'
							mat `pysw_h'			= (nullmat(`pysw_h'),r(pysw))
							mat `pysm_h'			= (nullmat(`pysm_h'),r(pysm))
						}
						else {
							// failed means estimation failed, no results to save
							// track so that pysw and psym can have missings inserted at end
							local failed_h			`failed_h' `k'
							if `notyetfailed'		di
							local notyetfailed = 0
							di as res "Warning: estimation failed with return code `r(rc)' when cross-fitting fold `k'."
						}
					}

					if `ssflag' {

						// need to enter only if short-stacking and main call to pystacked is LIE compliant
						// in that case the main call to pystacked uses the stacked vhat_not_k rather than the actual vname (D)
						// here, all the individual learners use the same actual D in a single call to pystacked
						// which is why a single call works (otherwise, would have to call pystacked separately for every learner)

						// always use non-LIE compliant for short-stacking
						// reason is that LIE would require indiv learners to all use different E^[D|XZ]
						local est_main_h_ss = subinstr("`est_main_h'","{D}","`vname'",1)

						// save warning and display after SS is complete
						if `enforceflag' {
							local fiv_warn "warning - enforcement of LIE with pystacked integration not supported with short-stackeding"
						}

						_fit_pystacked,							///
							k(`k')								///
							m(`m')								///
							nlearners(`nlearners')				///
							stdflag(`stdflag')					///
							cap(`cap')							///
							qui(`qui')							///
							est_main(`est_main_h_ss')			///
							est_options(`est_options_h')		///
							stdfinalestopt(`stdfinalestopt')	///
							nostdstackopt(`nostdstackopt')		///
							vname(`vhat_not_k_ss)')				///
							vhat(`hhat_ss')						/// where to save (replace) stacked OOS (crossfit) predicted values
							vtilde(`vtilde_h_list')				/// where to save (replace) OOS predicted values for indiv learners
							fid(`fid')							///
							touse(`touse')						///
							mmatname(`y_stacking_cv_h_ss')		///
							etype(`etype')
					}

				}	// end of pystacked FIV h block
				
				else {
					// Step 1 (main call) to pystacked will usually have multiple learners within pystacked.
					// Step 2 (h call), if not pystacked, will be a single learner that applies to all Step 1 learners.
					// vhat_not_k will have the h fitted values of the stacked Step 1 learners
					// vt_not_k_list will have the h fitted values corresponding to the Step 1 learners, so need to loop through these.
					// hhat will have the 2nd step OOS h for the stacked learners
					// vtilde_h_list will have the 2nd step OOS h for the base learners

					if `enforceflag' {
						// replace {D}-placeholder in estimation string with variable name
						local est_main_h_st = subinstr("`est_main_h'","{D}","`vhat_not_k'",1)
					}
					else {
						// replace {D}-placeholder with varname of original D
						local est_main_h_st = subinstr("`est_main_h'","{D}","`vname'",1)
					}

					_fit_other,	k(`k')						///
						m(`m')								///
						cap(`cap')							///
						qui(`qui')							///
						est_main(`est_main_h_st')			///
						est_options(`est_options_h')		///
						vname(`vhat_not_k')					///
						vhat(`hhat')						/// where to save (replace) stacked OOS (crossfit) predicted values
						fid(`fid')							///
						touse(`touse')
	
					forvalues j=1/`nlearners' {
						local vn : word `j' of `vt_not_k_list'
						local vh : word `j' of `vtilde_h_list'
						if `enforceflag'	local est_main_h_st = subinstr("`est_main_h'","{D}","`vn'",1)
						else 				local est_main_h_st = subinstr("`est_main_h'","{D}","`vname'",1)

						_fit_other,	k(`k')						///
							m(`m')								///
							cap(`cap')							///
							qui(`qui')							///
							est_main(`est_main_h_st')			///
							est_options(`est_options_h')		///
							vname(`vn')							///
							vhat(`vh')							/// where to save (replace) stacked OOS (crossfit) predicted values
							fid(`fid')							///
							touse(`touse')

					}

					local base_est_h		
					local base_est_h		`r(base_est)'
					local stack_final_est_h	"n.a."
					// pystacked weights go in as missing with same dimensions as main pystacked call
					mat `pysw_h'			= `pysw' * .
					mat `pysm_h'			= (nullmat(`pysm_h'),r(pysm))

					if `ssflag' {

						if `enforceflag' {
							// for fold k: nnls of D on in-sample base learners predicted values (vt_not_k_list)
							`qui' _ddml_nnls `vname' `vt_not_k_list' if `fid'~=`k' & `touse'
							tempvar vt_ss dhat_not_k_ss
							qui predict double `vt_ss' if `touse'
							qui gen `dhat_not_k_ss'=`vt_ss' if `fid'~=`k' & `touse'
							// replace {D}-placeholder in estimation string with variable name
							local est_main_h_ss = subinstr("`est_main_h'","{D}","`dhat_not_k_ss'",1)
						}
						else {
							// replace {D}-placeholder with varname of original D
							local est_main_h_ss = subinstr("`est_main_h'","{D}","`vname'",1)
						}

						// h estimation needs to be done only once since it's the same learner each time
						`qui' `est_main_h_ss' if `fid'~=`k' & `touse', `est_options_h'
						tempvar hhat_ss_k
						qui predict double `hhat_ss_k' if `fid'==`k' & `touse'
						qui replace `hhat_ss' = `hhat_ss_k' if `fid'==`k' & `touse'

					}	// end of ss section for FIV h block

				}	// end of non-pystacked FIV h block

			}	// end of FIV block

			else {
				// not case 1, 2 (treatvar), 3 (fiv) so...
				di as err "internal crossfit error"
				exit 198
			}

		}	// end of resamples loop

		// update number of crossfit reps
		mata: `ename'.crossfitted = `lastrep'

		// last fold, report completed
		di as text "...completed cross-fitting" _c
		// if noisily, print new line
		`qui' di

		******************************** SHORTSTACKING ************************************

		// shortstacking. we need to distinguish between 3 cases again. 
		if `ssflag' {
		
			tempname ssw ssw0 ssw1 ssw_h
		
			// if short-stacking final estimator not supplied, used std stack final est
			// unless stdflag==0, which means pystacked used voting, in which case use finalest(.)
			if "`ssfinalest'"=="" & `stdflag'	local ssfinalest `stack_final_est'
			else if "`ssfinalest'"==""			local ssfinalest `finalest'

			if "`est_type'"=="partial" { // case 1 - partially-linear
				if `nlearners'>1 {
					// stack dep var against all OOS (cross-fit) learner predicted values
					`qui' di
					`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (additive model):"
					`qui' _ddml_nnls `vname' `vtilde_list' if `touse', finalest(`ssfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mat `ssw' = e(b)
					cap drop `shortstack'_ss_`m'
					mat score double `shortstack'_ss_`m' = `ssw' if `touse'
				}
				else {
					cap drop `shortstack'_ss_`m'
					qui gen double `shortstack'_ss_`m' = `vhat' if `touse'
					mat `ssw' = 1
				}
				// if no std stacking, update with default
				if "`ssfinalest'"==""			local ssfinalest `e(finalest)'

			}	// end case 1 - partially-linear
			else if "`est_type'"=="interactive" {	// case 2 - interactive
				if `nlearners'>1 {
					// treatvar == 1
					`qui' di
					`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=1):"
					`qui' _ddml_nnls `vname' `vtilde1_list' if `touse' & `treatvar'==1, finalest(`ssfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mat `ssw1' = e(b)
					cap drop `shortstack'_ss1_`m'
					mat score double `shortstack'_ss1_`m' = `ssw1' if `touse'
					// if no std stacking, update with default
					if "`ssfinalest'"==""			local ssfinalest `e(finalest)'
					
					// treatvar == 0
					// possible in LATE model that D is always =0 if Z=0 (perfect assignment to treatment)
					// check estimation sample for this and handle as a special case
					qui count if `vname'!=0 & `treatvar'==0 & `touse'
					if (r(N)==0 & "`etype'"=="D") {
						cap drop `shortstack'_ss0_`m'
						qui gen double `shortstack'_ss0_`m' = 0
						mat `ssw0'		= J(1,`nlearners',.)
					}
					else {
						`qui' di
						`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=0):"
						`qui' _ddml_nnls `vname' `vtilde0_list' if `touse' & `treatvar'==0, finalest(`ssfinalest') stype(`stype') `noisily'
						`qui' di as text "N=" as res e(N)
						mat `ssw0' = e(b)
						cap drop `shortstack'_ss0_`m'
						mat score double `shortstack'_ss0_`m' = `ssw0' if `touse'
					}
				}
				else {
					cap drop `shortstack'_ss1_`m'
					cap drop `shortstack'_ss0_`m'
					qui gen double `shortstack'_ss1_`m' = `vhat1' if `touse'
					qui gen double `shortstack'_ss0_`m' = `vhat0' if `touse'
					mat `ssw1' = 1
					mat `ssw0' = 1
				}
			}	// end case 2 - interactive
			else if "`est_type'"=="fiv" { // case 3 - FIV
			
				// Step I: apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
				if `nlearners'>1 {
					// stack D against all OOS (cross-fit) learner predicted values
					`qui' di
					`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (additive model):"
					`qui' _ddml_nnls `vname' `vtilde_list' if `touse', finalest(`ssfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mat `ssw' = e(b)
					cap drop `shortstack'_ss_`m'
					mat score double `shortstack'_ss_`m' = `ssw' if `touse'
				}
				else {
					cap drop `shortstack'_ss_`m'
					qui gen double `shortstack'_ss_`m' = `vhat' if `touse'
					mat `ssw' = 1
				}
				// if no std stacking, update with default
				if "`ssfinalest'"==""			local ssfinalest `e(finalest)'

				// Step II: apply short-stacking on fitted values of learners

				if "`cmd_h'"=="pystacked" & (`nlearners'==`nlearners_h') {

					// multiple pystacked h learners
					// di as res "Final step of SS: nnls of vhat_ss (=dhat_oosSS)"
					// di "  on fitted values of learners (vtilde_h_list = hhatSS_list = all hhatSSj)"
					// di "vhat_ss:"
					// sum `vhat_ss' `vtilde_h_list' if `touse'
					
					`qui' _ddml_nnls `vhat_ss' `vtilde_h_list' if `touse', finalest(`ssfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mat `ssw_h' = e(b)
					tempvar vtemp
					qui predict double `vtemp' if `touse'
					cap drop `shortstack'_h_ss_`m'
					qui gen double `shortstack'_h_ss_`m' = `vtemp'

				}
				else {
					// single non-pystacked non-stacking learner for h.
					// hence no need for stacking
					cap drop `shortstack'_h_ss_`m'
					qui gen double `shortstack'_h_ss_`m' = `hhat_ss'
					// ssw_h will have 1 for the first learner and zeros for the rest
					mat `ssw_h' = `ssw' * 0
					mat `ssw_h'[1,1] = 1
				}
			}	// end case 3 - fiv
			else {
				di as err "internal ddml error - unknown est_type=`est_type'"
				exit 198
			}
		}
	
		if `ssflag' & `nlearners'>1 {
			di as text "...completed short-stacking" _c
		}

		******************************** POOLSTACKING *************************************
		
		if `psflag' {
		
			tempname psw psw0 psw1 psw_h
		
			// if pool-stacking final estimator not supplied, used std stack final est
			// note that pool-stacking can take place only if std stacking was done
			// hence macro stack_final_est is the appropriate default
			if "`psfinalest'"  ==""		local psfinalest   `stack_final_est'
			if "`psfinalest_h'"==""		local psfinalest_h `stack_final_est_h'
			
			// mata object y_stacking_cv has y and predicted yhats of all learners in all crossfits
			if "`est_type'"=="partial" { // case 1 - partially-linear
				if `nlearners'>1 {
					tempname tframe
					qui frame pwf
					local cframe `r(currentframe)'
					frame create `tframe'
					frame change `tframe'
					getmata (`fid' `fidtouse' `vname' `vtilde_list')=`y_stacking_cv', force replace
					`qui' di
					`qui' di as text "Pooled-stacking, finalest=`psfinalest' (additive model):"
					`qui' _ddml_nnls `vname' `vtilde_list', finalest(`psfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mata: `psw' = st_matrix("e(b)")
					frame change `cframe'
					frame drop `tframe'
					mata: st_matrix("`psw'",`psw')
					mat colnames `psw' =  `vtilde_list'
					cap drop `poolstack'_ps_`m'
					mat score double `poolstack'_ps_`m' = `psw' if `touse'
				}
				else {
					cap drop `poolstack'_ps_`m'
					qui gen double `poolstack'_ps_`m' = `vhat' if `touse'
					mat `psw' = 1
				}
			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive

				if `nlearners'>1 {
					tempname tframe
					
					// treatvar=1
					qui frame pwf
					local cframe `r(currentframe)'
					frame create `tframe'
					frame change `tframe'
					getmata (`fid' `fidtouse' `vname' `vtilde1_list')=`y_stacking_cv1', force replace
					`qui' di
					`qui' di as text "Pooled-stacking, finalest=`psfinalest' (interactive model):"
					`qui' _ddml_nnls `vname' `vtilde1_list', finalest(`psfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mata: `psw1' = st_matrix("e(b)")
					frame change `cframe'
					frame drop `tframe'
					mata: st_matrix("`psw1'",`psw1')
					mat colnames `psw1' = `vtilde1_list'
					cap drop `poolstack'_ps1_`m'
					mat score double `poolstack'_ps1_`m' = `psw1' if `touse'
	
					// treatvar=0
					// possible in LATE model that D is always =0 if Z=0 (perfect assignment to treatment)
					// check estimation sample for this and handle as a special case (no estimation took place)
					qui count if `vname'!=0 & `treatvar'==0 & `touse'
					if (r(N)==0 & "`etype'"=="D") {
						cap drop `poolstack'_ps0_`m'
						qui gen double `poolstack'_ps0_`m' = 0
						mat `psw0'		= J(1,`nlearners',.)
					}
					else {
						qui frame pwf
						local cframe `r(currentframe)'
						frame create `tframe'
						frame change `tframe'
						getmata (`fid' `fidtouse' `vname' `vtilde0_list')=`y_stacking_cv0', force replace
						`qui' di
						`qui' di as text "Pooled-stacking, finalest=`psfinalest' (interactive model):"
						`qui' _ddml_nnls `vname' `vtilde0_list', finalest(`psfinalest') stype(`stype') `noisily'
						`qui' di as text "N=" as res e(N)
						mata: `psw0' = st_matrix("e(b)")
						frame change `cframe'
						frame drop `tframe'
						mata: st_matrix("`psw0'",`psw0')
						mat colnames `psw0' =  `vtilde0_list'
						cap drop `poolstack'_ps0_`m'
						mat score double `poolstack'_ps0_`m' = `psw0' if `touse'
					}
				}
				else {
					cap drop `poolstack'_ps1_`m'
					cap drop `poolstack'_ps0_`m'
					qui gen double `poolstack'_ps1_`m' = `vhat1' if `touse'
					qui gen double `poolstack'_ps0_`m' = `vhat0' if `touse'
					mat `psw1' = 1
					mat `psw0' = 1
				}
			}
			else if "`est_type'"=="fiv" { // case 3 - FIV

				if `nlearners'>1 {
					tempname tframe
					
					// Step I: estimation of E[D|Z,X]=D^
					qui frame pwf
					local cframe `r(currentframe)'
					frame create `tframe'
					frame change `tframe'
					getmata (`fid' `fidtouse' `vname' `vtilde_list')=`y_stacking_cv', force replace
					`qui' di
					`qui' di as text "Pooled-stacking, finalest=`psfinalest' (additive model):"
					`qui' _ddml_nnls `vname' `vtilde_list', finalest(`psfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mata: `psw' = st_matrix("e(b)")
					frame change `cframe'
					frame drop `tframe'
					mata: st_matrix("`psw'",`psw')
					mat colnames `psw' =  `vtilde_list'
					cap drop `poolstack'_ps_`m'
					mat score double `poolstack'_ps_`m' = `psw' if `touse'
	
					// Step II: estimation of E[D|X]
					qui frame pwf
					local cframe `r(currentframe)'
					frame create `tframe'
					frame change `tframe'
					getmata (`fid' `fidtouse' `vhat_not_k' `vtilde_h_list')=`y_stacking_cv_h', force replace
					`qui' di
					`qui' di as text "Pooled-stacking, finalest_h=`psfinalest_h' (additive model):"
					`qui' _ddml_nnls `vhat_not_k' `vtilde_h_list', finalest(`psfinalest_h') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					mata: `psw_h' = st_matrix("e(b)")
					frame change `cframe'
					frame drop `tframe'
					mata: st_matrix("`psw_h'",`psw_h')
					mat colnames `psw_h' =  `vtilde_h_list'
					cap drop `poolstack'_h_ps_`m'
					mat score double `poolstack'_h_ps_`m' = `psw_h' if `touse'
				}
				else {
					cap drop `poolstack'_ps_`m'
					cap drop `poolstack'_h_ps_`m'
					qui gen double `poolstack'_ps_`m' = `vtilde'_L1_`m' if `touse'
					qui gen double `poolstack'_h_ps_`m' = `vtilde'_h_L1_`m' if `touse'
					mat `psw' = 1
					mat `psw_h' = 1
				}

			}
			else {
				di as err "internal crossfit error"
				exit 198
			}
		}
		
		if `psflag' & `nlearners'>1 {
			di as text "...completed pooled-stacking" _c
		}
		
		************************************ MISC *****************************************
		
		if "`fiv_warn'"~="" {
			di
			di as res "`fiv_warn'"
			di as res "to enforce LIE with short-stacking, re-specify model without pystacked integration"
			di as res "and specify individual learners separately instead of in a single pystacked call" _c
		}

		***************************** ESTIMATION COMPLETE *********************************
		// estimation done, insert newline
		di
		******************************** STORE RESULTS ************************************

		// vtilde, mspe, etc.
		if "`est_type'"=="partial" { // case 1 - partially-linear
	
			// always label learner predicted values
			local vt_L_list
			forvalues j=1/`nlearners' {
				qui label var `vtilde'_L`j'_`m' "Pred. values E[`vname'|X] using base learner `j', rep `m'"
				local vt_L_list `vt_L_list' `vtilde'_L`j'_`m'
			}
			// always available even if no standard stacking
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est","`base_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est","`stack_final_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_type","`stype'")
			
			// rsq, mse, N by learner - save under vtilde
			rsqmse `vt_L_list' if `touse', yvar(`vname')
			mata: add_result_item(`ename',"`vtilde'","MSE_L",  "`m'", st_matrix("r(mse)"))
			mata: add_result_item(`ename',"`vtilde'","RMSE_L", "`m'", st_matrix("r(rmse)"))
			mata: add_result_item(`ename',"`vtilde'","R-sq_L", "`m'", st_matrix("r(rsq)"))
			mata: add_result_item(`ename',"`vtilde'","N_L",    "`m'", st_matrix("r(N)"))
			
			// cvc by learner - save under vname
			cvc `vt_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			mata: add_result_item(`ename',"`vname'","cvc_pval", "`m'", st_matrix("r(pmat)"))
			mata: add_result_item(`ename',"`vname'","cvc_bootnum", "`m'", `cvcbootnum')
			
			// save results relating to stacked learner if it exists
			if `stdflag' {
			
				// vtilde has fitted values
				cap drop `vtilde'_`m'
				qui gen `vtype' `vtilde'_`m' = `vhat'
				qui label var `vtilde'_`m' "Pred. values E[`vname'|X] using `cmd', rep `m'"

				// weights and MSEs will be missing values if #learners=1
				if "`failed'"~="" {
					add_m_col `pysw', flist(`failed') lastcol(`kfolds')
					mat `pysw' = r(A)
					add_m_col `pysm', flist(`failed') lastcol(`kfolds')
					mat `pysm' = r(A)
				}
				mata: add_result_item(`ename',"`vtilde'","stack_weights", "`m'", st_matrix("`pysw'"))
				mata: add_result_item(`ename',"`vtilde'","stack_MSEs",	  "`m'", st_matrix("`pysm'"))
				
				// stacked learner, full sample - mse, rsq, etc.
				rsqmse `vtilde'_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename',"`vtilde'","MSE",  "`m'", `mse')
				mata: add_result_item(`ename',"`vtilde'","RMSE", "`m'", `rmse')
				mata: add_result_item(`ename',"`vtilde'","R-sq", "`m'", `rsq')
				mata: add_result_item(`ename',"`vtilde'","N",    "`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')

				// stacked learner by fold - mse, rsq, etc.				
				rsqmse `vtilde'_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`vtilde'","MSE_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`vtilde'","N_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')   \ r(N))

			}
			
			// only one learner so it's the opt; set opt="" if no std stacking
			if `stdflag'	mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
			else			mata: add_learner_item(`ename',"opt","`m'","")
			
		}	// end case 1 - partially linear
		else if "`est_type'"=="interactive" {	// case 2 - interactive
		
			// always label learner predicted values
			local vt0_L_list
			local vt1_L_list
			forvalues j=1/`nlearners' {
				qui label var `vtilde'0_L`j'_`m' "Pred. values E[`vname'|X] given `treatvar'==0 using base learner `j', rep `m'"
				qui label var `vtilde'1_L`j'_`m' "Pred. values E[`vname'|X] given `treatvar'==1 using base learner `j', rep `m'"
				local vt0_L_list `vt0_L_list' `vtilde'0_L`j'_`m'
				local vt1_L_list `vt1_L_list' `vtilde'1_L`j'_`m'
			}
			// always available even if no standard stacking
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est","`base_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est","`stack_final_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_type","`stype'")
			
			// cvc by learner - save under vname
			cvc `vt1_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			mata: add_result_item(`ename',"`vname'","cvc_pval1", "`m'", st_matrix("r(pmat)"))
			cvc `vt0_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			mata: add_result_item(`ename',"`vname'","cvc_pval0", "`m'", st_matrix("r(pmat)"))
			// same bootnum for both
			mata: add_result_item(`ename',"`vname'","cvc_bootnum", "`m'", `cvcbootnum')
			
			// rsq, mse, N by learner - save under vtilde
			rsqmse `vt1_L_list' if `treatvar'==1 & `touse', yvar(`vname')
			mata: add_result_item(`ename',"`vtilde'","MSE_L1",  "`m'", st_matrix("r(mse)"))
			mata: add_result_item(`ename',"`vtilde'","RMSE_L1", "`m'", st_matrix("r(rmse)"))
			mata: add_result_item(`ename',"`vtilde'","R-sq_L1", "`m'", st_matrix("r(rsq)"))
			mata: add_result_item(`ename',"`vtilde'","N_L1",    "`m'", st_matrix("r(N)"))
			rsqmse `vt0_L_list' if `treatvar'==0 & `touse', yvar(`vname')
			mata: add_result_item(`ename',"`vtilde'","MSE_L0",  "`m'", st_matrix("r(mse)"))
			mata: add_result_item(`ename',"`vtilde'","RMSE_L0", "`m'", st_matrix("r(rmse)"))
			mata: add_result_item(`ename',"`vtilde'","R-sq_L0", "`m'", st_matrix("r(rsq)"))
			mata: add_result_item(`ename',"`vtilde'","N_L0",    "`m'", st_matrix("r(N)"))
			
			// save results relating to stacked learner if it exists
			if `stdflag' {
			
				cap drop `vtilde'0_`m'
				cap drop `vtilde'1_`m'
				// vtilde is predicted values
				qui gen `vtype' `vtilde'0_`m' = `vhat0'
				qui gen `vtype' `vtilde'1_`m' = `vhat1'
				qui label var `vtilde'0_`m' "Pred. values E[`vname'|X] given `treatvar'==0 using `cmd', rep `m'"
				qui label var `vtilde'1_`m' "Pred. values E[`vname'|X] given `treatvar'==1 using `cmd', rep `m'"

				// weights and MSEs will be missing values if #learners=1
				if "`failed1'"~="" {
					add_m_col `pysw1', flist(`failed1') lastcol(`kfolds')
					mat `pysw1' = r(A)
					add_m_col `pysm1', flist(`failed1') lastcol(`kfolds')
					mat `pysm1' = r(A)
				}
				if "`failed0'"~="" {
					add_m_col `pysw0', flist(`failed0') lastcol(`kfolds')
					mat `pysw0' = r(A)
					add_m_col `pysm0', flist(`failed0') lastcol(`kfolds')
					mat `pysm0' = r(A)
				}
				mata: add_result_item(`ename',"`vtilde'","stack_weights0","`m'", st_matrix("`pysw0'"))
				mata: add_result_item(`ename',"`vtilde'","stack_weights1","`m'", st_matrix("`pysw1'"))
				mata: add_result_item(`ename',"`vtilde'","stack_MSEs0","`m'",	st_matrix("`pysm0'"))
				mata: add_result_item(`ename',"`vtilde'","stack_MSEs1","`m'",	st_matrix("`pysm1'"))
				
				// calculate and return mse and sample size
				// interactive-type model, return mse etc. separately for treatvar =0 and =1
				// treatvar=0
				rsqmse `vtilde'0_`m' if `treatvar'==0 & `touse', yvar(`vname')
				local mse0			= el(r(mse),1,1)
				local rmse0			= el(r(rmse),1,1)
				local rsq0			= el(r(rsq),1,1)
				local N0			= el(r(N),1,1)
				// treatvar=1
				rsqmse `vtilde'1_`m' if `treatvar'==1 & `touse', yvar(`vname')
				local mse1			= el(r(mse),1,1)
				local rmse1			= el(r(rmse),1,1)
				local rsq1			= el(r(rsq),1,1)
				local N1			= el(r(N),1,1)
				local N				= `N0'+`N1'
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`vtilde'","MSE`t'",  "`m'", `mse`t'')
					mata: add_result_item(`ename',"`vtilde'","RMSE`t'", "`m'", `rmse`t'')
					mata: add_result_item(`ename',"`vtilde'","R-sq`t'",	"`m'", `rsq`t'')
					mata: add_result_item(`ename',"`vtilde'","N`t'",    "`m'", `N`t'')
				}
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list')   \ `N0')
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list')   \ `N1')

				// by fold and treat=0: mse, rsq, etc.				
				rsqmse `vtilde'0_`m' if `treatvar'==0 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`vtilde'","MSE0_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`vtilde'","N0_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ r(mse))
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')   \ r(N))
				
				// by fold and treat=1: mse, rsq, etc.				
				rsqmse `vtilde'1_`m' if `treatvar'==1 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`vtilde'","MSE1_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`vtilde'","N1_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ r(mse))
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')   \ r(N))

			}

			// only one learner so it's the opt; set opt="" if no std stacking
			forvalues t=0/1 {
				if `stdflag'	mata: add_learner_item(`ename',"opt`t'","`m'","`vtilde'")
				else			mata: add_learner_item(`ename',"opt`t'","`m'","")
			}
		}	// end case 2 - interactive
		
		else if "`est_type'"=="fiv" { // case 3 - FIV

			// always label learner predicted values
			local vt_L_list
			local vt_h_L_list
			forvalues j=1/`nlearners' {
				qui label var `vtilde'_L`j'_`m'   "Pred. values E[`vname'|X,Z], rep `m'"
				qui label var `vtilde'_h_L`j'_`m' "Pred. values E[`vtilde'|X], rep `m'"
				local vt_L_list `vt_L_list' `vtilde'_L`j'_`m'
				local vt_h_L_list `vt_h_L_list' `vtilde'_h_L`j'_`m'
			}
			// always available even if no standard stacking
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est_h","`base_est_h'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est_h","`stack_final_est_h'")
			mata: add_learner_item(`ename',"`vtilde'","stack_type_h","`stype_h'")
			
			// cvc by learner - save under vname
			cvc `vt_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			mata: add_result_item(`ename',"`vname'","cvc_pval", "`m'", st_matrix("r(pmat)"))
			cvc `vt_h_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			mata: add_result_item(`ename',"`vname'","cvc_pval_h", "`m'", st_matrix("r(pmat)"))
			// same bootnum for both
			mata: add_result_item(`ename',"`vname'","cvc_bootnum", "`m'", `cvcbootnum')
			
			// rsq, mse, N by learner - save under vtilde
			rsqmse `vt_L_list' if `touse', yvar(`vname')
			mata: add_result_item(`ename',"`vtilde'","MSE_L",    "`m'", st_matrix("r(mse)"))
			mata: add_result_item(`ename',"`vtilde'","RMSE_L",   "`m'", st_matrix("r(rmse)"))
			mata: add_result_item(`ename',"`vtilde'","R-sq_L",   "`m'", st_matrix("r(rsq)"))
			mata: add_result_item(`ename',"`vtilde'","N_L",      "`m'", st_matrix("r(N)"))
			rsqmse `vt_h_L_list' if `touse', yvar(`vname')
			mata: add_result_item(`ename',"`vtilde'","MSE_L_h",  "`m'", st_matrix("r(mse)"))
			mata: add_result_item(`ename',"`vtilde'","RMSE_L_h", "`m'", st_matrix("r(rmse)"))
			mata: add_result_item(`ename',"`vtilde'","R-sq_L_h", "`m'", st_matrix("r(rsq)"))
			mata: add_result_item(`ename',"`vtilde'","N_L_h",    "`m'", st_matrix("r(N)"))
			
			// save results relating to stacked learner if it exists
			if `stdflag' {
			
				cap drop `vtilde'_`m'
				cap drop `vtilde'_h_`m'
				// vtilde and vtilde_h are predicted values
				qui gen `vtype' `vtilde'_`m'   = `vhat'
				qui gen `vtype' `vtilde'_h_`m' = `hhat'
				qui label var `vtilde'_`m'   "Pred. values E[`vname'|X,Z] using `cmd', rep `m'"
				qui label var `vtilde'_h_`m' "Pred. values E[`vtilde'|X] using `cmd', rep `m'"
				
				// weights and MSEs will be missing values if #learners=1
				if "`failed'"~="" {
					add_m_col `pysw', flist(`failed') lastcol(`kfolds')
					mat `pysw' = r(A)
					add_m_col `pysm', flist(`failed') lastcol(`kfolds')
					mat `pysm' = r(A)
				}
				if "`failed_h'"~="" {
					add_m_col `pysw_h', flist(`failed_h') lastcol(`kfolds')
					mat `pysw_h' = r(A)
					add_m_col `pysm_h', flist(`failed_h') lastcol(`kfolds')
					mat `pysm_h' = r(A)
				}
				mata: add_result_item(`ename',  "`vtilde'", "stack_weights",   "`m'", st_matrix("`pysw'"))
				mata: add_result_item(`ename',  "`vtilde'", "stack_MSEs",      "`m'", st_matrix("`pysm'"))
				mata: add_learner_item(`ename', "`vtilde'", "stack_base_est",         "`base_est'")
				mata: add_learner_item(`ename', "`vtilde'", "stack_final_est",        "`stack_final_est'")
				mata: add_learner_item(`ename', "`vtilde'", "stack_type",             "`stype'")
				
				mata: add_result_item(`ename',  "`vtilde'", "stack_weights_h", "`m'", st_matrix("`pysw_h'"))
				mata: add_result_item(`ename',  "`vtilde'", "stack_MSEs_h",    "`m'", st_matrix("`pysm_h'"))
				
				// stacked learner, full sample - mse, rsq, etc.
				rsqmse `vtilde'_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename',"`vtilde'","MSE",  "`m'", `mse')
				mata: add_result_item(`ename',"`vtilde'","RMSE", "`m'", `rmse')
				mata: add_result_item(`ename',"`vtilde'","R-sq", "`m'", `rsq')
				mata: add_result_item(`ename',"`vtilde'","N",    "`m'", `N')
				rsqmse `vtilde'_h_`m' if `touse', yvar(`vname')
				local mse_h			= el(r(mse),1,1)
				local rmse_h		= el(r(rmse),1,1)
				local rsq_h			= el(r(rsq),1,1)
				local N_h			= el(r(N),1,1)
				mata: add_result_item(`ename',"`vtilde'","MSE_h",  "`m'", `mse_h')
				mata: add_result_item(`ename',"`vtilde'","RMSE_h", "`m'", `rmse_h')
				mata: add_result_item(`ename',"`vtilde'","R-sq_h", "`m'", `rsq_h')
				mata: add_result_item(`ename',"`vtilde'","N_h",    "`m'", `N_h')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list')   \ `mse')
				mat `N_list'			= (nullmat(`N_list')     \ `N')
				mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
				mat `N_h_list'			= (nullmat(`N_h_list')   \ `N_h')

				// stacked learner by fold - mse, rsq, etc.				
				rsqmse `vtilde'_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',  "`vtilde'", "MSE_folds",       "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',  "`vtilde'", "N_folds",         "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list')   \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')     \ r(N))
				rsqmse `vtilde'_h_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',  "`vtilde'", "MSE_h_folds",     "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',  "`vtilde'", "N_h_folds",       "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ r(mse))
				mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')   \ r(N))
			}
			
			// only one learner so it's the opt; set opt="" if no std stacking
			if `stdflag' {
				mata: add_learner_item(`ename',"opt",   "`m'", "`vtilde'")
				mata: add_learner_item(`ename',"opt_h", "`m'", "`vtilde'_h")
			}
			else {
				mata: add_learner_item(`ename',"opt",   "`m'", "")
				mata: add_learner_item(`ename',"opt_h", "`m'", "")
			}
			
		}	// end case 3 - fiv
		else {
			di as err "internal crossfit error"
			exit 198
		} // end vtilde, mspe, etc.

		
		// add standard stacking base learner CV predictions
		if `stdflag' {
			if "`est_type'"=="partial" { // case 1 - partially-linear
				mata: add_result_item(`ename',"`vtilde'","y_stacking_cv", "`m'", `y_stacking_cv')
			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`vtilde'","y_stacking_cv`t'", "`m'", `y_stacking_cv`t'')
				}
			}
			else if "`est_type'"=="fiv" { // case 3 - FIV
				mata: add_result_item(`ename',"`vtilde'","y_stacking_cv",   "`m'", `y_stacking_cv')
				mata: add_result_item(`ename',"`vtilde'","y_stacking_cv_h", "`m'", `y_stacking_cv_h')
			}
			else {
				di as err "internal crossfit error"
				exit 198
			}
		} // end standard stacking base learner predictions

		
		// add shortstack results
		if `ssflag' {
			if "`est_type'"=="partial" { // case 1 - partially-linear
		
				label var `shortstack'_ss_`m' "Pred. values E[`vname'|X] using shortstacking, rep `m'"
				
				// save weights as column vector
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights", "`m'", st_matrix("`ssw'")')
				// save base estimator list with rest of shortstack results
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_base_est","`base_est'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_type","`stype'")
				
				// full sample - mse, rsq, etc.
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename',"`shortstack'_ss","MSE",	"`m'", `mse')
				mata: add_result_item(`ename',"`shortstack'_ss","RMSE",	"`m'", `rmse')
				mata: add_result_item(`ename',"`shortstack'_ss","R-sq",	"`m'", `rsq')
				mata: add_result_item(`ename',"`shortstack'_ss","N",	"`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				
				// by fold - mse, rsq, etc.				
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`shortstack'_ss","N_folds",	 "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')   \ r(N))

			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive
			
				label var `shortstack'_ss0_`m'  "Pred. values E[`vname'|X] given `treatvar'==0 using shortstacking, rep `m'"
				label var `shortstack'_ss1_`m'  "Pred. values E[`vname'|X] given `treatvar'==1 using shortstacking, rep `m'"
				
				// save weights as column vectors
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights0",	"`m'", st_matrix("`ssw0'")')
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights1",	"`m'", st_matrix("`ssw1'")')
				// save base estimator list with rest of shortstack results
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_base_est","`base_est'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_type","`stype'")

				// calculate and return mspe and sample size
				// interactive-type model, return mse separately for treatvar =0 and =1
				// treatvar=0
				rsqmse `shortstack'_ss0_`m' if `treatvar'==0 & `touse', yvar(`vname')
				local mse0			= el(r(mse),1,1)
				local rmse0			= el(r(rmse),1,1)
				local rsq0			= el(r(rsq),1,1)
				local N0			= el(r(N),1,1)
				// treatvar=1
				rsqmse `shortstack'_ss1_`m' if `treatvar'==1 & `touse', yvar(`vname')
				local mse1			= el(r(mse),1,1)
				local rmse1			= el(r(rmse),1,1)
				local rsq1			= el(r(rsq),1,1)
				local N1			= el(r(N),1,1)
				local N				= `N0'+`N1'
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`shortstack'_ss","MSE`t'",  "`m'", `mse`t'')
					mata: add_result_item(`ename',"`shortstack'_ss","RMSE`t'", "`m'", `rmse`t'')
					mata: add_result_item(`ename',"`shortstack'_ss","R-sq`t'",	"`m'", `rsq`t'')
					mata: add_result_item(`ename',"`shortstack'_ss","N`t'",    "`m'", `N`t'')
				}
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list')   \ `N0')
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list')   \ `N1')

				// by fold and treat=0: mse, rsq, etc.				
				rsqmse `shortstack'_ss0_`m' if `treatvar'==0 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE0_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`shortstack'_ss","N0_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ r(mse))
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')   \ r(N))
				// by fold and treat=1: mse, rsq, etc.				
				rsqmse `shortstack'_ss1_`m' if `treatvar'==1 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE1_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`shortstack'_ss","N1_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ r(mse))
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')   \ r(N))

			}
			else if "`est_type'"=="fiv" { // case 3 - FIV
				label var `shortstack'_ss_`m'	"Pred. values E[`vname'|X,Z] using shortstacking, rep `m'"
				label var `shortstack'_h_ss_`m'  "Pred. values E[`vtilde'|X] using shortstacking, rep `m'"
				
				// save weights as column vector
				mata: add_result_item(`ename',  "`shortstack'_ss", "ss_weights",       "`m'", st_matrix("`ssw'")')
				mata: add_result_item(`ename',  "`shortstack'_ss", "ss_weights_h",     "`m'", st_matrix("`ssw_h'")')
				// save base estimator list with rest of shortstack results
				mata: add_learner_item(`ename', "`shortstack'_ss", "stack_base_est",   "`base_est'")
				mata: add_learner_item(`ename', "`shortstack'_ss", "stack_base_est_h", "`base_est_h'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename', "`shortstack'_ss", "ss_final_est",     "`ssfinalest'")
				mata: add_learner_item(`ename', "`shortstack'_ss", "ss_final_est_h",   "`ssfinalest_h'")
				mata: add_learner_item(`ename', "`shortstack'_ss", "stack_type",       "`stype'")

				// rsq, mse, N - save under vtilde
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE",    "`m'", `mse')
				mata: add_result_item(`ename', "`shortstack'_ss", "RMSE",   "`m'", `rmse')
				mata: add_result_item(`ename', "`shortstack'_ss", "R-sq",   "`m'", `rsq')
				mata: add_result_item(`ename', "`shortstack'_ss", "N",      "`m'", `N')
				rsqmse `shortstack'_h_ss_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE_h",  "`m'", `mse')
				mata: add_result_item(`ename', "`shortstack'_ss", "RMSE_h", "`m'", `rmse')
				mata: add_result_item(`ename', "`shortstack'_ss", "R-sq_h", "`m'", `rsq')
				mata: add_result_item(`ename', "`shortstack'_ss", "N_h",    "`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list')   \ `mse')
				mat `N_list'			= (nullmat(`N_list')     \ `N')
				mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
				mat `N_h_list'			= (nullmat(`N_h_list')   \ `N_h')

				// by fold - mse, rsq, etc.				
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE_folds",   "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`shortstack'_ss", "N_folds",     "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list')   \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')     \ r(N))
				rsqmse `shortstack'_h_ss_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE_h_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`shortstack'_ss", "N_h_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ r(mse))
				mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')   \ r(N))

			}
			else {
				di as err "internal crossfit error"
			}
		} // end shortstack results
	
		// add poolstack results
		if `psflag' {
			if "`est_type'"=="partial" { // case 1 - partially-linear
		
				label var `poolstack'_ps_`m' "Pred. values E[`vname'|X] using poolstacking, rep `m'"
				
				// save weights as column vector
				mata: add_result_item(`ename',"`poolstack'_ps","ps_weights",	"`m'", st_matrix("`psw'")')
				// save base estimator list with rest of poolstack results
				mata: add_learner_item(`ename',"`poolstack'_ps","stack_base_est","`base_est'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`poolstack'_ps","ps_final_est", "`psfinalest'")
				mata: add_learner_item(`ename',"`poolstack'_ps","stack_type","`stype'")
				
				// calculate and return mspe and sample size
				// poolstack macros have fitted values
				qui replace `vres_sq' = (`vname' - `poolstack'_ps_`m')^2 if `touse'
			
				// full sample - mse, rsq, etc.
				rsqmse `poolstack'_ps_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename',"`poolstack'_ps","MSE",		"`m'", `mse')
				mata: add_result_item(`ename',"`poolstack'_ps","RMSE",		"`m'", `mse')
				mata: add_result_item(`ename',"`poolstack'_ps","R-sq",		"`m'", `rsq')
				mata: add_result_item(`ename',"`poolstack'_ps","N",			"`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				
				// by fold - mse, rsq, etc.				
				rsqmse `poolstack'_ps_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`poolstack'_ps","MSE_folds",	"`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`poolstack'_ps","N_folds",	"`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')   \ r(N))
				
			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive
			
				label var `poolstack'_ps0_`m'  "Pred. values E[`vname'|X] given `treatvar'==0 using poolstacking, rep `m'"
				label var `poolstack'_ps1_`m'  "Pred. values E[`vname'|X] given `treatvar'==1 using poolstacking, rep `m'"
				
				// save weights as column vector
				mata: add_result_item(`ename',"`poolstack'_ps","ps_weights0",	"`m'", st_matrix("`psw0'")')
				mata: add_result_item(`ename',"`poolstack'_ps","ps_weights1",	"`m'", st_matrix("`psw1'")')
				// save base estimator list with rest of poolstack results
				mata: add_learner_item(`ename',"`poolstack'_ps","stack_base_est","`base_est'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`poolstack'_ps","ps_final_est", "`psfinalest'")
				mata: add_learner_item(`ename',"`poolstack'_ps","stack_type","`stype'")

				// calculate and return mspe and sample size
				// interactive-type model, return mse separately for treatvar =0 and =1
				// treatvar=0
				rsqmse `poolstack'_ps0_`m' if `treatvar'==0 & `touse', yvar(`vname')
				local mse0			= el(r(mse),1,1)
				local rmse0			= el(r(rmse),1,1)
				local rsq0			= el(r(rsq),1,1)
				local N0			= el(r(N),1,1)
				// treatvar=1
				rsqmse `poolstack'_ps1_`m' if `treatvar'==1 & `touse', yvar(`vname')
				local mse1			= el(r(mse),1,1)
				local rmse1			= el(r(rmse),1,1)
				local rsq1			= el(r(rsq),1,1)
				local N1			= el(r(N),1,1)
				local N				= `N0'+`N1'
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`poolstack'_ps","MSE`t'",  "`m'", `mse`t'')
					mata: add_result_item(`ename',"`poolstack'_ps","RMSE`t'", "`m'", `rmse`t'')
					mata: add_result_item(`ename',"`poolstack'_ps","R-sq`t'",	"`m'", `rsq`t'')
					mata: add_result_item(`ename',"`poolstack'_ps","N`t'",    "`m'", `N`t'')
				}
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list')   \ `N0')
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list')   \ `N1')

				// by fold and treat=0: mse, rsq, etc.				
				rsqmse `poolstack'_ps0_`m' if `treatvar'==0 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`poolstack'_ps","MSE0_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`poolstack'_ps","N0_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ r(mse))
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')   \ r(N))
				// by fold and treat=1: mse, rsq, etc.				
				rsqmse `poolstack'_ps1_`m' if `treatvar'==1 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename',"`poolstack'_ps","MSE1_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename',"`poolstack'_ps","N1_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ r(mse))
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')   \ r(N))
			}
			else if "`est_type'"=="fiv" { // case 3 - FIV
				label var `poolstack'_ps_`m'	"Pred. values E[`vname'|X,Z] using poolstacking, rep `m'"
				label var `poolstack'_h_ps_`m'  "Pred. values E[`vtilde'|X] using poolstacking, rep `m'"

				// calculate and return mspe and sample size
				// poolstack macros have fitted values
				qui replace `vres_sq' = (`vname' - `poolstack'_ps_`m')^2 if `touse'
				qui replace `hres_sq' = (`vname' - `poolstack'_h_ps_`m')^2 if `touse'

				// save weights as column vector
				mata: add_result_item(`ename',  "`poolstack'_ps", "ps_weights",       "`m'", st_matrix("`psw'")')
				mata: add_result_item(`ename',  "`poolstack'_ps", "ps_weights_h",     "`m'", st_matrix("`psw_h'")')
				// save base estimator list with rest of poolstack results
				mata: add_learner_item(`ename', "`poolstack'_ps", "stack_base_est",   "`base_est'")
				mata: add_learner_item(`ename', "`poolstack'_ps", "stack_base_est_h", "`base_est_h'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename', "`poolstack'_ps", "ps_final_est",     "`psfinalest'")
				mata: add_learner_item(`ename', "`poolstack'_ps", "ps_final_est_h",   "`psfinalest_h'")
				mata: add_learner_item(`ename', "`poolstack'_ps", "stack_type",       "`stype'")

				// rsq, mse, N - save under vtilde
				rsqmse `poolstack'_ps_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename', "`poolstack'_ps", "MSE",    "`m'", `mse')
				mata: add_result_item(`ename', "`poolstack'_ps", "RMSE",   "`m'", `rmse')
				mata: add_result_item(`ename', "`poolstack'_ps", "R-sq",   "`m'", `rsq')
				mata: add_result_item(`ename', "`poolstack'_ps", "N",      "`m'", `N')
				rsqmse `poolstack'_h_ps_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename', "`poolstack'_ps", "MSE_h",  "`m'", `mse')
				mata: add_result_item(`ename', "`poolstack'_ps", "RMSE_h", "`m'", `rmse')
				mata: add_result_item(`ename', "`poolstack'_ps", "R-sq_h", "`m'", `rsq')
				mata: add_result_item(`ename', "`poolstack'_ps", "N_h",    "`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list')   \ `mse')
				mat `N_list'			= (nullmat(`N_list')     \ `N')
				mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
				mat `N_h_list'			= (nullmat(`N_h_list')   \ `N_h')

				// by fold - mse, rsq, etc.				
				rsqmse `poolstack'_ps_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`poolstack'_ps", "MSE_folds",   "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`poolstack'_ps", "N_folds",     "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list')   \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')     \ r(N))
				rsqmse `poolstack'_h_ps_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`poolstack'_ps", "MSE_h_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`poolstack'_ps", "N_h_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ r(mse))
				mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')   \ r(N))

			}
			else {
				di as err "internal crossfit error"
				exit 198
			}
		} // end poolstack results
		
		// final clean up
		cap mata: mata drop `y_stacking_cv'
		cap mata: mata drop `y_stacking_cv0'
		cap mata: mata drop `y_stacking_cv1'
		cap mata: mata drop `y_stacking_cv_h'
		cap mata: mata drop `y_stacking_cv_h_ss'
		cap mata: mata drop `psw'
		cap mata: mata drop `psw1'
		cap mata: mata drop `psw0'
		cap mata: mata drop `ssw'
		cap mata: mata drop `ssw0'
		cap mata: mata drop `ssw1'

	}	// end of resampling loop

		
	******************************** RETURN RESULTS ************************************

	if "`est_type'"=="partial" { // case 1 - partially-linear
		foreach mname in mse_list N_list mse_folds_list N_folds_list {
			mat rownames ``mname''	= `rnames'
		}
		return mat mse				= `mse_list'
		return mat N				= `N_list'
		return mat mse_folds		= `mse_folds_list'
		return mat N_folds			= `N_folds_list'
	}
	else if "`est_type'"=="interactive" {	// case 2 - interactive
		foreach mname in mse0_list N0_list mse0_folds_list N0_folds_list mse1_list N1_list mse1_folds_list N1_folds_list {
			mat rownames ``mname''	= `rnames'
		}
		return mat mse0				= `mse0_list'
		return mat N0				= `N0_list'
		return mat mse1				= `mse1_list'
		return mat N1				= `N1_list'
		return mat mse0_folds		= `mse0_folds_list'
		return mat mse1_folds		= `mse1_folds_list'
		return mat N0_folds			= `N0_folds_list'
		return mat N1_folds			= `N1_folds_list'
	}
	else if "`est_type'"=="fiv" { // case 3 - FIV
		foreach mname in mse_list N_list mse_folds_list N_folds_list mse_h_list N_h_list mse_h_folds_list N_h_folds_list {
			mat rownames ``mname''	= `rnames'
		}
		return mat mse				= `mse_list'
		return mat N				= `N_list'
		return mat mse_h			= `mse_h_list'
		return mat N_h				= `N_h_list'
		return mat mse_folds		= `mse_folds_list'
		return mat mse_h_folds		= `mse_h_folds_list'
		return mat N_folds			= `N_folds_list'
		return mat N_h_folds		= `N_h_folds_list'
		// special case - pystacked + LIE not enforced with short-stacking
		return local fiv_warn		`fiv_warn'
	}
	else {
		di as err "internal crossfit error"
		exit 198
	}
	
	qui count if `touse'
	return scalar N				= r(N)
	return local cmd_list		pystacked
	return scalar perfectflag	= `perfectflag'
	
end


program define _crossfit_other, rclass sortpreserve
	version 16
	syntax [if] [in] ,								///
							[ estring(string asis)	/// estimation string
													/// need asis option in case it includes strings
													///
							ename(name)				/// name for Mata struct; default is "crossfit"
							foldvar(varlist)		/// one fold var per resample
							kfolds(integer 5)		/// ignored if foldvars provided
							reps(integer 1)			/// ignored if foldvars provided
							firstrep(integer 1)		/// add reps starting with m=firstrep
							NORANDOM				/// first fold ID uses obs in existing order
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							NOENFORCElie			/// do not enforce LIE in fiv model
							NOIsily					///
							ssfinalest(name)		/// final estimator for short-stacking
							finalest(name)			/// final estimator for both
							stdflag(integer 0)		/// flag for pystacked with multiple learners
							NOSTDstack				/// no standard stacking - pystacked only, not supported, used to catch error
													/// legal but ignored options:
							psfinalest(name)		///
							stdfinalest(name)		///
							CVCbootnum(integer 500)	///
							*						/// ignored options - will trigger syntax error
							]

	// syntax errors
	if "`options'"~="" {
		di as err "crossfit error: illegal option '`options''"
		exit 197
	}

	*** set sample
	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `vname'
	
	// final estimator choice; default is NNLS + coefs sum to 1
	if "`finalest'"==""		local finalest nnls1
	if "`ssfinalest'"==""	local ssfinalest `finalest'
	
	mata: st_local("vname", `ename'.vname)
	mata: st_local("vtlist", invtokens(`ename'.vtlist))
	mata: st_local("nlearners", strofreal(`ename'.nlearners))
	mata: st_local("shortstack", `ename'.shortstack)
	mata: st_local("poolstack", `ename'.poolstack)
	mata: st_local("etype", `ename'.etype)	// "Y", "D" or "Z"

	// 3 cases: basic partially linear, interactive, fiv
	if "`etype'"=="DH"			local est_type "fiv"
	else if "`treatvar'"~=""	local est_type "interactive"
	else						local est_type "partial"
	
	** indicator for LATE model special case: perfect assignment to treatment
	// nb: D (treatment, `vname') can be continuous; code below aimed at binary treatment
	//	 Z (assignment, `treatvar') is always binary
	// (note bad naming convention)
	// perfect assignment: no one treated if not assigned
	// means that in the treatvar Z=0 subsample of the estimation sample, no one is assigned to treatment (vname D=0)
	// since no variation in D, estimation fails, so catch here as a special case
	// and set predicted values vhat of D to 0
	// applies to full sample rather than fold-by-fold
	if "`est_type'"=="interactive" & "`etype'"=="D" {
		qui count if `treatvar'==0 & `touse'
		local Z0_count = r(N)
		qui count if `treatvar'==0 & `vname'==0 & `touse'
		local D0Z0_count = r(N)
		local perfectflag = (`Z0_count'==`D0Z0_count')
	}
	else {
		local perfectflag = 0
	}

	** indicator for enforcing LIE in fiv model **
	local enforceflag	= "`noenforcelie'"==""

	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
	** indicator for pooled-stacking
	local psflag	= "`poolstack'"~=""
	
	if "`noisily'"=="" {
		local qui quietly
		local cap cap
	}

	** error check
	if "`nostdstack'"~="" & "`est_type'"=="fiv" {
		di as err "err: option nostdstack invalid with pystacked+fiv model; standard stacking required"
		exit 198
	}
	
	*** debugging message
	`qui' di as text "entering _crossfit_other..."
	
	*** set up fold/reps variables
	
	// fold id variables
	if "`foldvar'"~="" {
		// fold variables provided
		foreach fvar of varlist `foldvar' {
			tempvar fid
			qui egen `fid' = group(`fvar')
			local fidlist `fidlist' `fid'
		}
		check_foldvar, fidlist(`fidlist') touse(`touse') `noisily'
		local kfolds = r(kfolds)
		local reps = r(reps)
		local lastrep = `firstrep' + `reps' - 1
	}
	else {
		// fold variables not provided, so generate
		local lastrep = `firstrep' + `reps' - 1
		forvalues m=`firstrep'/`lastrep' {
			*** gen folds
			tempvar uni cuni fid
			if `m'==1 & "`norandom'"~="" {
				qui gen `uni' = _n
			}
			else {
				qui gen double `uni' = runiform() if `touse'
			}
			qui cumul `uni' if `touse', gen(`cuni')
			qui gen int `fid' = ceil(`kfolds'*`cuni') if `touse'
			local fidlist `fidlist' `fid'
		}
	}
	
	*** syntax

	if `psflag' {
		di as err "error - {opt poolstack} available only with pystacked integration; also unavailable for fiv model"
		exit 198
	}
	
	// stored results; initialize here
	tempname mse_list N_list mse_folds_list N_folds_list
	tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
	tempname mse0_list N0_list mse0_folds_list N0_folds_list
	tempname mse1_list N1_list mse1_folds_list N1_folds_list
	// pystacked support for fiv/lie
	tempname pysw pysw_h pysw_temp pysm pysm_h pysm_temp
	// initialize temp vars here
	tempvar dhat_oosSS dhatSS hhatSS vres_sq vres0_sq vres1_sq hres_sq dres_sq
	qui gen double `dhat_oosSS'=.
	qui gen double `dhatSS'=.
	qui gen double `hhatSS'=.
	qui gen double `vres_sq'=.
	qui gen double `vres0_sq'=.
	qui gen double `vres1_sq'=.
	qui gen double `hres_sq'=.
	qui gen double `dres_sq'=.
	forvalues j=1/`nlearners' {
		tempvar vhat`j' vres`j' vhat0`j' vres0`j' vhat1`j' vres1`j' dhat`j' hhat`j' hhatSS`j'
		qui gen double `vhat`j''=.
		qui gen double `vres`j''=.
		qui gen double `vhat0`j''=.
		qui gen double `vres0`j''=.
		qui gen double `vhat1`j''=.
		qui gen double `vres1`j''=.
		qui gen double `dhat`j''=.
		qui gen double `hhat`j''=.
		qui gen double `hhatSS`j''=.
		forvalues k=1/`kfolds' {
			tempvar dhat_`j'_`k'
			qui gen double `dhat_`j'_`k''=.
		}
	}
	forvalues k = 1(1)`kfolds' {
		tempvar dhat_isSS_`k'
		qui gen double `dhat_isSS_`k''=.
	}
	
	// loop over resamples (crossfitting, shortstacking, store results)
	forvalues m=`firstrep'/`lastrep' {
	
		// create rnames for renaming the rows of saved matrices
		local rnames `rnames' resample_`m'
		
		// initialize for each resample
		local failed
		local failed0
		local failed1
		local failed_lie_h
		cap mat drop `pysw'
		cap mat drop `pysw_h'
		cap mat drop `pysm'
		cap mat drop `pysm_h'
	
		******************************** CROSSFITTING ************************************
		
		// cross-fitting fundamentally depends on three cases: partially-linear, interactive, FIV
	
		if `lastrep'>1 {
			di as text "Resample `m'..."
		}

		// list of fold IDs refers to current set of reps (not full set including appended)
		local thisrep = `m'-`firstrep'+1
		local fid : word `thisrep' of `fidlist'

		// create or (re-)initialize blank fitted variable(s)
		if "`est_type'"=="partial" { // case 1 - partially-linear
			if `ssflag' {
				cap drop `shortstack'_ss_`m'
				qui gen double `shortstack'_ss_`m'=.
			}
			forvalues j=1/`nlearners' {
				qui replace `vhat`j''=.
				qui replace `vres`j''=.
			}
		}
		else if "`est_type'"=="interactive" {	// case 2 - interactive
			if `ssflag' {
				cap drop `shortstack'_ss0_`m'
				cap drop `shortstack'_ss1_`m'
				qui gen double `shortstack'_ss0_`m'=.
				qui gen double `shortstack'_ss1_`m'=.
			}
			forvalues j=1/`nlearners' {
				qui replace `vhat0`j''=.
				qui replace `vres0`j''=.	
				qui replace `vhat1`j''=.
				qui replace `vres1`j''=.	
			}
		}
		else if "`est_type'"=="fiv" { // case 3 - FIV
			// out-of-sample predicted values for E[D|XZ] 
			forvalues j=1/`nlearners' {
				qui replace `dhat`j''=.  // using learner i
			}
			// out-of-sample predicted values for E[D|X] 
			forvalues j=1/`nlearners' {
				qui replace `hhat`j''=.  // using learner i in both steps
			}
			// predicted values for E[D|ZX] for each k & learner i
			forvalues k=1/`kfolds' {
				forvalues j=1/`nlearners' {
					qui replace `dhat_`j'_`k''=.
				}
			}
			if `ssflag' { // with short-stacking
				// for final results
				cap drop `shortstack'_ss_`m'
				cap drop `shortstack'_h_ss_`m'
				qui gen double `shortstack'_ss_`m'=.
				qui gen double `shortstack'_h_ss_`m'=.			
				// in-sample predicted values for E[D|ZX] for each k: short-stacked
				forvalues k=1/`kfolds' {		
					qui replace `dhat_isSS_`k''=.
				}
				// out-of-sample predicted values for E[D|ZX] from short-stacking by k
				// NB: this is different from "dhatSS", which are from applying constrained regression to full sample
				qui replace `dhat_oosSS'=. 
				// out-of-sample predicted values for E[D|X] 
				forvalues j=1/`nlearners' {
					qui replace `hhatSS`j''=. // using short-stacking for E[D|XZ] & then using learner i for E[D|X]
				}
				// short-stacking predicted values
				qui replace `dhatSS'=. // this will become `shortstack'
				qui replace `hhatSS'=. // this will become `shortstack'_h_ss_
			} 
		}
		else {
			di as err "internal crossfit error"
			exit 198
		}
		
		// crossfit
		di as text "Cross-fitting fold " _c
		forvalues k = 1(1)`kfolds' {
	
			// initialize
			local notyetfailed = 1
			di as text "`k' " _c
	
			forvalues j=1/`nlearners' {
				local vtilde : word `j' of `vtlist'
				mata: st_local("est_main", return_learner_item(`ename',"`vtilde'","est_main"))
				mata: st_local("est_options", return_learner_item(`ename',"`vtilde'","est_options"))
				mata: st_local("predopt",return_learner_item(`ename',"`vtilde'","predopt"))
				mata: st_local("vtype",return_learner_item(`ename',"`vtilde'","vtype"))
				if "`est_type'"=="fiv" { // case 3 - FIV
					// LIE locals
					local hhat `vtilde'_h
					mata: st_local("est_main_h", return_learner_item(`ename',"`vtilde'","est_main_h"))
					mata: st_local("est_options_h", return_learner_item(`ename',"`vtilde'","est_options_h"))
					mata: st_local("predopt_h",return_learner_item(`ename',"`vtilde'","predopt_h"))
				}
				if "`est_type'"=="partial" { // case 1 - partially-linear
				
					// need to declare here since needed for predict
					tempvar vhat_k
					
					// estimate excluding kth fold
					`cap' `est_main' if `fid'!=`k' & `touse', `est_options'

					local rc = _rc
					if `rc'>0 {
						// estimation failed
						if `notyetfailed'		di
						local notyetfailed = 0
						di as res "Warning: estimation failed with return code `rc' when cross-fitting fold `k'."
						qui gen `vtype' `vhat_k' = . if `touse', `predopt'
						qui replace `vhat`j'' = `vhat_k' if `fid'==`k' & `touse'
						qui replace `vres`j'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
						local failed `failed' `k'
					}
					else {
					
						local cmd `e(cmd)'
						if "`cmd'"=="" {
							// macro e(cmd) may be missing, so deduce from command line
							local cmd : word 1 of `est_main'
						}
						
						// special treatment of single call to pystacked with possible multiple learners
						if "`e(cmd)'"=="pystacked" {
							local base_est `e(base_est)'
							local stack_final_est `e(finalest)'
							local stype `e(type)'

							// save pystacked weights and MSEs
							qui pystacked, table(rmspe)		// create rmspe matrix
							mat `pysw_temp' = e(weights)
							mat `pysm_temp' = e(rmspe)
							mat `pysm_temp' = `pysm_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
						 	if `stdflag' & e(mcount)>1 {
								mat `pysw' = (nullmat(`pysw'),`pysw_temp')
								mata: st_replacematrix("`pysm_temp'",st_matrix("`pysm_temp'"):^2)
								mat `pysm' = (nullmat(`pysm'),`pysm_temp')
							}
							else {
								local stack_final_est "n.a."
								local stype "n.a."
								mat `pysw' = `pysw_temp' * .
								mat `pysm' = `pysm_temp' * .
							}
						}
						
						// get fitted values for kth fold	
						qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'
			
						// get predicted values
						qui replace `vhat`j'' = `vhat_k' if `fid'==`k' & `touse'
						qui replace `vres`j'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
					}
					
					// no longer needed
					drop `vhat_k'
				}

				else if "`est_type'"=="interactive" {	// case 2 - interactive

					// need to declare here since needed for predict
					tempvar vhat_k
		
					// outcome equation so estimate separately
		
					// for treatvar = 1
					// estimate excluding kth fold

					`cap' `est_main' if `fid'!=`k' & `treatvar' == 1 & `touse', `est_options'

					local rc = _rc
					if `rc'>0 {
						// estimation failed
						if `notyetfailed'		di
						local notyetfailed = 0
						di as res "Warning: estimation failed with return code `rc' when cross-fitting fold `k'."
						qui gen `vtype' `vhat_k' = . if `touse', `predopt'
						local failed1 `failed1' `k'
					}
					else {
						local cmd `e(cmd)'
						if "`cmd'"=="" {
							// macro e(cmd) may be missing, so deduce from command line
							local cmd : word 1 of `est_main'
						}
			
						// get fitted values for kth fold	
						qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'
					}
					qui replace `vhat1`j'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres1`j'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
					
					// no longer needed
					drop `vhat_k'

					// for treatvar = 0

					// special LATE case - perfect assignment
					if `perfectflag' {
						// fitted values for kth fold are all zeros with perfect assignment and treatvar=0
						`qui' gen double `vhat_k'=0  if `fid'==`k' & `touse' 	
					} 
					else {

						// estimate excluding kth fold
						`cap' `est_main' if `fid'!=`k' & `treatvar' == 0 & `touse', `est_options'

						local rc = _rc
						if `rc'>0 {
							// estimation failed
							if `notyetfailed'		di
							local notyetfailed = 0
							di as res "Warning: estimation failed with return code `rc' when cross-fitting fold `k'."
							qui gen `vtype' `vhat_k' = . if `touse', `predopt'
							local failed0 `failed0' `k'
						}
						else {
				
							// get fitted values for kth fold	
							qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'					
						}
					}
					qui replace `vhat0`j'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres0`j'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
					
					// no longer needed
					drop `vhat_k'
					
				}
		
				else if "`est_type'"=="fiv" { // case 3 - FIV
		
					// need to declare here since needed for predict
					tempvar vhat_k // stores predicted values for E[D|Z,X] temporarily
					tempvar vtil_k // stores predicted values for E[D|X] temporarily
		
					// Step I: estimation of E[D|Z,X]=D^
					// estimate excluding kth fold
					`cap' `est_main' if `fid'!=`k' & `touse', `est_options'
					local rc = _rc
					if `rc'>0 {
						// estimation failed
						if `notyetfailed'		di
						local notyetfailed = 0
						di as res "Warning: estimation of E[D|Z,X] failed with return code `rc' when cross-fitting fold `k'."
						qui gen `vtype' `vhat_k' = . if `touse', `predopt'
						// get *combined* out-of-sample predicted values
						qui replace `dhat`j'' = `vhat_k' if `fid'==`k' & `touse'
						// get predicted values in wide format for each i and k
						qui replace `dhat_`j'_`k'' = `vhat_k' if `touse'
						local failed `failed' `k'
					}
					else {
						if "`cmd'"=="" {
							// macro e(cmd) may be missing, so deduce from command line
							local cmd : word 1 of `est_main'
						}
						
						// special treatment of single call to pystacked with possible multiple learners
						if "`e(cmd)'"=="pystacked" {
							local base_est `e(base_est)'
							local stack_final_est `e(finalest)'
							local stype `e(type)'

							// save pystacked weights and MSEs
							qui pystacked, table(rmspe)		// create rmspe matrix
							mat `pysw_temp' = e(weights)
							mat `pysm_temp' = e(rmspe)
							mat `pysm_temp' = `pysm_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
							if  `stdflag' & e(mcount)>1 {
								mat `pysw' = (nullmat(`pysw'),`pysw_temp')
								mata: st_replacematrix("`pysm_temp'",st_matrix("`pysm_temp'"):^2)
								mat `pysm' = (nullmat(`pysm'),`pysm_temp')
							}
							else {
								local stack_final_est "n.a."
								local stype "n.a."
								mat `pysw' = `pysw_temp' * .
								mat `pysm' = `pysm_temp' * .
							}
						}
			
						// get fitted values (in and out of sample)
						qui predict `vtype' `vhat_k' if `touse', `predopt'
			
						// get *combined* out-of-sample predicted values
						qui replace `dhat`j'' = `vhat_k' if `fid'==`k' & `touse'

						// get predicted values in wide format for each j and k
						qui replace `dhat_`j'_`k'' = `vhat_k' if `touse'

					}
		
					// Step II: estimation of E[D|X]
		
					if `enforceflag' {
						// replace {D}-placeholder in estimation string with variable name
						local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_`j'_`k''",1)
					}
					else {
						// replace {D}-placeholder with varname of original D
						local est_main_h_k = subinstr("`est_main_h'","{D}","`vname'",1)
					}

					// estimation	
					`cap' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
					local rc = _rc
					if `rc'>0 {
						// estimation failed
						if `notyetfailed'		di
						local notyetfailed = 0
						di as res "Warning: estimation of E[D|X] failed with return code `rc' when cross-fitting fold `k'."
						qui gen `vtype' `vtil_k' = . if `touse', `predopt_h'
						// get *combined* out-of-sample predicted values
						qui replace `hhat`j'' = `vtil_k' if `fid'==`k' & `touse'
						local failed_lie_h `failed_lie_h' `k'
					}
					else {
						local cmd_h `e(cmd)'
			
						// special treatment of single call to pystacked with multiple learners
						// h estimation
						if "`e(cmd_h)'"=="pystacked" &`stdflag' & e(mcount)>1 {
							local base_est_h `e(base_est)'
							local stack_final_est_h `e(finalest)'
							local stype_h `e(type)'
							assert e(mcount)>1 & e(mcount)<.
							// save pystacked weights and MSEs
							qui pystacked, table(rmspe)		// create rmspe matrix
							mat `pysw_temp' = e(weights)
							mat `pysw_h' = (nullmat(`pysw_h'),`pysw_temp')
							mat `pysm_temp' = e(rmspe)
							mat `pysm_temp' = `pysm_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
							mata: st_replacematrix("`pysm_temp'",st_matrix("`pysm_temp'"):^2)
							mat `pysm_h' = (nullmat(`pysm_h'),`pysm_temp')
						}
						else if `stdflag' & "`e(cmd)'"=="pystacked" {
							// main command is pystacked but h command may not be pystacked
							local base_est_h "`e(cmd)'"
							local stack_final_est_h "n.a."
							local stype_h "n.a."
							mat `pysw_h' = `pysw' * .
							mat `pysm_h' = `pysm' * .
						}
						
						// get fitted values  
						qui predict `vtype' `vtil_k' if `touse', `predopt_h'
			
						// get *combined* out-of-sample predicted values
						qui replace `hhat`j'' = `vtil_k' if `fid'==`k' & `touse'

					}
					
					// no longer needed
					drop `vhat_k'
					drop `vtil_k'
				}
				
				if `k'==1 & `m'==1 {
					local cmd_list   `cmd_list'   `cmd'
					local cmd_h_list `cmd_h_list' `cmd_h'
				}
			}
		}
	
		// update number of crossfit reps
		mata: `ename'.crossfitted = `lastrep'
		// last fold, report completed
		di as text "...completed cross-fitting" _c

		******************************** SHORTSTACKING ************************************

		// shortstacking. we need to distinguish between 3 cases again. 

		if `ssflag' {
		
			tempname ssw ssw1 ssw0 ssw_h
			
			if "`est_type'"=="partial" { // case 1 - partially-linear
				local vhats
				forvalues j=1/`nlearners' {
					local vhats `vhats' `vhat`j''
				}
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (additive model):"
				if `nlearners'>1 {
					`qui' _ddml_nnls `vname' `vhats', finalest(`ssfinalest') `noisily'
					`qui' di as text "N=" as res e(N)
					mat `ssw' = e(b)
					// declare here since needed for predict
					tempvar vtemp
					qui predict double `vtemp'
					qui replace `shortstack'_ss_`m' = `vtemp'
					// no longer needed
					drop `vtemp'
				}
				else {
					qui replace `shortstack'_ss_`m' = `vhat1'
					mat `ssw' = 1
				}
			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive
				local vhats1 
				local vhats0
				if `nlearners'>1 {
					// declare here since needed for predict
					tempvar vtemp
					
					// treatvar == 1
					forvalues j=1/`nlearners' {
						local vhats1 `vhats1' `vhat1`j''
					}
					`qui' di
					`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=1):"
					`qui' _ddml_nnls `vname' `vhats1' if `treatvar'==1, finalest(`ssfinalest') `noisily'
					`qui' di as text "N=" as res e(N)
					mat `ssw1' = e(b)
					qui predict double `vtemp'
					qui replace `shortstack'_ss1_`m'=`vtemp'
					// no longer needed
					drop `vtemp'
						
					// treatvar == 0
					// if late and Z==0, we know Dhat should be 0 (eps = appx 0)
					if "`etype'"=="D" {
						cap drop `shortstack'_ss0_`m'
						qui gen double `shortstack'_ss0_`m' = 0
						mat `ssw0' = J(1,`nlearners',.)
					}
					else {
						forvalues j=1/`nlearners' {
							local vhats0 `vhats0' `vhat0`j''
						}
						`qui' di
						`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=0):"
						`qui' _ddml_nnls `vname' `vhats0' if `treatvar'==0, finalest(`ssfinalest') `noisily'
						`qui' di as text "N=" as res e(N)
						mat `ssw0' = e(b)
						qui predict double `vtemp'
						qui replace `shortstack'_ss0_`m'=`vtemp'
					}
					// no longer needed
					cap drop `vtemp'
				}
				else {
					qui replace `shortstack'_ss1_`m' = `vhat11' if `touse'
					qui replace `shortstack'_ss0_`m' = `vhat01' if `touse'
					mat `ssw1' = 1
					mat `ssw0' = 1
				}
			}
			else if "`est_type'"=="fiv" { // case 3 - FIV

				// apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
				local dhats
				forvalues j=1/`nlearners' {
					local dhats `dhats' `dhat`j''
				}
				`qui' di
				`qui' di as text "Short-stacking, finalest=`ssfinalest' (LIE, OOS E[D|XZ]):"

				`qui' _ddml_nnls `vname' `dhats' if `touse'
				mat `ssw'= e(b)
				// declare here since needed for predict
				tempvar vtemp
				qui predict double `vtemp' if `touse'
				qui replace `dhatSS'=`vtemp' 
				// no longer needed
				drop `vtemp'

				// apply short-stacking to in-sample predicted values of E[D|XZ] *for each k*
				forvalues k = 1(1)`kfolds' {
					local dhats_is
					forvalues j=1/`nlearners' {
						local dhats_is `dhats_is' `dhat_`j'_`k''
					}
					// declare here since needed for predict
					tempvar vtemp

					`qui' di
					`qui' di as text "Short-stacking, finalest=`ssfinalest' (LIE, in-sample E[D|XZ] fold `k':"
					`qui' _ddml_nnls `vname' `dhats_is' if `fid'!=`k' & `touse' 
					qui predict double `vtemp'
					qui replace `dhat_isSS_`k'' = `vtemp' if `fid'!=`k' & `touse'
					qui replace `dhat_oosSS' = `vtemp' if `fid'==`k' & `touse'

					// no longer needed
					drop `vtemp'
				}

				// need to cross-fit stacked in-sample predicted values against X
				forvalues k = 1(1)`kfolds' {

					forvalues j=1/`nlearners' {
						local vtilde : word `j' of `vtlist'
						mata: st_local("est_main_h", return_learner_item(`ename',"`vtilde'","est_main_h"))
						mata: st_local("est_options_h", return_learner_item(`ename',"`vtilde'","est_options_h"))
						mata: st_local("predopt_h",return_learner_item(`ename',"`vtilde'","predopt_h"))
						mata: st_local("vtype",return_learner_item(`ename',"`vtilde'","vtype"))				
	
						if `enforceflag' {
							// replace {D}-placeholder in estimation string with variable name
							local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_isSS_`k''",1)
						}
						else {
							// replace {D}-placeholder with varname of original D
							local est_main_h_k = subinstr("`est_main_h'","{D}","`vname'",1)
						}
						
						// estimation
						`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
						local cmd_h `e(cmd)'
					
						// get fitted values  
						// declare here since needed for predict
						tempvar vtemp
						qui predict double `vtemp' if `touse', `predopt_h'
						// get out-of-sample predicted values
						qui replace `hhatSS`j'' = `vtemp' if `fid'==`k' & `touse'
						// no longer needed
						drop `vtemp'
					}
				}
				// final stacking for E[D|X]
				local hhatSS_list
				forvalues j=1/`nlearners' {
					local hhatSS_list `hhatSS_list' `hhatSS`j''
				}

				`qui' di
				`qui' di as text "Short-stacking, finalest=`ssfinalest' (LIE, E[D|X]):"
				`qui' _ddml_nnls `dhat_oosSS' `hhatSS_list'
				mat `ssw_h' = e(b)
				// declare here since needed for predict
				tempvar vtemp

				qui predict double `vtemp'
				qui replace `hhatSS'=`vtemp'
				qui replace `shortstack'_ss_`m'=`dhatSS'
				qui replace `shortstack'_h_ss_`m'=`hhatSS'
				// no longer needed
				drop `vtemp'

			}
			else {
				di as err "internal crossfit error"
				exit 198
			}
		}
		else if `ssflag' {
			// single learner case, so shortstack vars are just copies of learner vars
			
			if "`est_type'"=="partial" { // case 1 - partially-linear
				qui replace `shortstack'_ss_`m' = `vhat1'
			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive
				qui replace `shortstack'_ss1_`m'=`vhat11'
				qui replace `shortstack'_ss0_`m'=`vhat01'
			}
			else if "`est_type'"=="fiv" { // case 3 - FIV
				qui replace `shortstack'_ss_`m'=`dhat1'
				label var `shortstack'_ss_`m' "Predicted values E[D|Z,X] of `vname' using shortstacking, rep `m'"
				qui replace `shortstack'_h_ss_`m'=`hhat1'
				label var `shortstack'_h_ss_`m' "Predicted values E[D^|X] of `vname' using shortstacking, rep `m'"
			}
		}
	
		if `ssflag' & `nlearners'>1 {
			// last fold, insert new line
			di as text "...completed short-stacking"
		}
		else {
			di
		}
		
		******************************** STORE RESULTS ************************************
		
		// reset locals
		local vt_L_list
		local vt0_L_list
		local vt1_L_list
		local vt_h_L_list
		
		forvalues j=1/`nlearners' {
			
			local vtilde	: word `j' of `vtlist'
			local cmd		: word `j' of `cmd_list'
			local cmd_h		: word `j' of `cmd_h_list'

			// vtilde, mspe, etc.
			if "`est_type'"=="partial" { // case 1 - partially-linear
		
				// not a temp var so can overwrite
				cap drop `vtilde'_`m'
				// vtilde is predicted values
				mata: st_local("vtype", return_learner_item(`ename',"`vtilde'","vtype"))
				qui gen `vtype' `vtilde'_`m' = `vhat`j''
				qui label var `vtilde'_`m' "Pred. values E[`vname'|X] using `cmd', rep `m'"
				local vt_L_list `vt_L_list' `vtilde'_`m'
				
				// rsq, mse, N
				rsqmse `vtilde'_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename', "`vtilde'", "MSE",  "`m'", `mse')
				mata: add_result_item(`ename', "`vtilde'", "RMSE", "`m'", `rmse')
				mata: add_result_item(`ename', "`vtilde'", "R-sq", "`m'", `rsq')
				mata: add_result_item(`ename', "`vtilde'", "N",    "`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list')   \ `N')
				
				// by fold - mse, rsq, etc.				
				rsqmse `vtilde'_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`vtilde'", "MSE_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`vtilde'", "N_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')   \ r(N))

				if `j'==1 {
					local mse_opt		= `mse'
					mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
				}
				else if `mse' < `mse_opt' {
					// overwrite with new opt
					local mse_opt		= `mse'
					mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
				}
			
			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive
			
				cap drop `vtilde'0_`m'
				cap drop `vtilde'1_`m'
				// vtilde is predicted values
				mata: st_local("vtype", return_learner_item(`ename',"`vtilde'","vtype"))
				qui gen `vtype' `vtilde'0_`m' = `vhat0`j''
				qui gen `vtype' `vtilde'1_`m' = `vhat1`j''
				qui label var `vtilde'0_`m' "Pred. values E[`vname'|X] given `treatvar'==0 using `cmd', rep `m'"
				qui label var `vtilde'1_`m' "Pred. values E[`vname'|X] given `treatvar'==1 using `cmd', rep `m'"
				local vt0_L_list `vt0_L_list' `vtilde'0_`m'
				local vt1_L_list `vt1_L_list' `vtilde'1_`m'
		
				// calculate and return mspe and sample size
				// interactive-type model, return mse separately for treatvar =0 and =1
				// treatvar=1
				rsqmse `vtilde'1_`m' if `treatvar'==1 & `touse', yvar(`vname')
				local mse1			= el(r(mse),1,1)
				local rmse1			= el(r(rmse),1,1)
				local rsq1			= el(r(rsq),1,1)
				local N1			= el(r(N),1,1)
				mata: add_result_item(`ename', "`vtilde'", "MSE1",  "`m'", `mse1')
				mata: add_result_item(`ename', "`vtilde'", "RMSE1", "`m'", `rmse1')
				mata: add_result_item(`ename', "`vtilde'", "R-sq1", "`m'", `rsq1')
				mata: add_result_item(`ename', "`vtilde'", "N1",    "`m'", `N1')
				// treatvar=0
				rsqmse `vtilde'0_`m' if `treatvar'==0 & `touse', yvar(`vname')
				local mse0			= el(r(mse),1,1)
				local rmse0			= el(r(rmse),1,1)
				local rsq0			= el(r(rsq),1,1)
				local N0			= el(r(N),1,1)
				mata: add_result_item(`ename', "`vtilde'", "MSE0",  "`m'", `mse0')
				mata: add_result_item(`ename', "`vtilde'", "RMSE0", "`m'", `rmse0')
				mata: add_result_item(`ename', "`vtilde'", "R-sq0", "`m'", `rsq0')
				mata: add_result_item(`ename', "`vtilde'", "N0",    "`m'", `N0')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list')   \ `N0')
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list')   \ `N1')
				
				// by fold - mse, rsq, etc.				
				rsqmse `vtilde'1_`m' if `treatvar'==1 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`vtilde'", "MSE1_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`vtilde'", "N1_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ r(mse))
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')   \ r(N))
				rsqmse `vtilde'0_`m' if `treatvar'==0 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`vtilde'", "MSE0_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`vtilde'", "N0_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ r(mse))
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')   \ r(N))
				
				forvalues t=0/1 {
					if `j'==1 {
						local mse`t'_opt		= `mse`t''
						mata: add_learner_item(`ename',"opt`t'","`m'","`vtilde'")
					}
					else if `mse`t'' < `mse`t'_opt' {
						// overwrite with new opt
						local mse`t'_opt		= `mse`t''
						mata: add_learner_item(`ename',"opt`t'","`m'","`vtilde'")
					}
				}
				
			}
			else if "`est_type'"=="fiv" { // case 3 - FIV
	
				cap drop `vtilde'_`m'
				cap drop `vtilde'_h_`m'
				
				mata: st_local("vtype", return_learner_item(`ename',"`vtilde'","vtype"))
				qui gen `vtype' `vtilde'_`m' = `dhat`j''
				qui label var `vtilde'_`m' "Pred. values E[`vname'|X,Z], rep `m'"
				qui gen `vtype' `vtilde'_h_`m' = `hhat`j''
				qui label var `vtilde'_h_`m' "Pred. values E[`vtilde'|X], rep `m'"
				local vt_L_list `vt_L_list' `vtilde'_`m'
				local vt_h_L_list `vt_h_L_list' `vtilde'_h_`m'

				if `ssflag' {
					// intermediate dhat predicted values by learner
					cap drop `vtilde'_h_ss_`m'
					qui gen `vtype' `vtilde'_h_ss_`m' = `hhatSS`j''
					qui label var `vtilde'_h_ss_`m' "Pred. values (ss) E[`vtilde'|X], rep `m'"
				}

				// rsq, mse, N
				rsqmse `vtilde'_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename', "`vtilde'", "MSE",  "`m'", `mse')
				mata: add_result_item(`ename', "`vtilde'", "RMSE", "`m'", `rmse')
				mata: add_result_item(`ename', "`vtilde'", "R-sq", "`m'", `rsq')
				mata: add_result_item(`ename', "`vtilde'", "N",    "`m'", `N')
				rsqmse `vtilde'_h_`m' if `touse', yvar(`vname')
				local mse_h			= el(r(mse),1,1)
				local rmse_h		= el(r(rmse),1,1)
				local rsq_h			= el(r(rsq),1,1)
				local N_h			= el(r(N),1,1)
				mata: add_result_item(`ename', "`vtilde'", "MSE_h",  "`m'", `mse_h')
				mata: add_result_item(`ename', "`vtilde'", "RMSE_h", "`m'", `rmse_h')
				mata: add_result_item(`ename', "`vtilde'", "R-sq_h", "`m'", `rsq_h')
				mata: add_result_item(`ename', "`vtilde'", "N_h",    "`m'", `N_h')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list')   \ `mse')
				mat `N_list'			= (nullmat(`N_list')     \ `N')
				mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
				mat `N_h_list'			= (nullmat(`N_h_list')   \ `N_h')
				
				// by fold - mse, rsq, etc.				
				rsqmse `vtilde'_`m' if `touse',   yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`vtilde'", "MSE_folds",   "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`vtilde'", "N_folds",     "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list')   \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')     \ r(N))
				rsqmse `vtilde'_h_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`vtilde'", "MSE_h_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`vtilde'", "N_h_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ r(mse))
				mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')   \ r(N))
				
				// optimal D and H
				if `j'==1 {
					local mse_opt		= `mse'
					mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
				}
				else if `mse' < `mse_opt' {
					// overwrite with new opt
					local mse_opt		= `mse'
					mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
				}
				if `j'==1 {
					local mse_h_opt		= `mse_h'
					mata: add_learner_item(`ename',"opt_h","`m'","`vtilde'_h")
				}
				else if `mse_h' < `mse_h_opt' {
					// overwrite with new opt
					local mse_h_opt		= `mse_h'
					mata: add_learner_item(`ename',"opt_h","`m'","`vtilde'_h")
				}
				
			}

		}
		
		// cvc
		if "`est_type'"=="partial" { // case 1 - partially-linear
			cvc `vt_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			tempname pmat
			mat `pmat' = r(pmat)
			mata: add_result_item(`ename',"`vname'","cvc_pval", "`m'", st_matrix("`pmat'"))
			mata: add_result_item(`ename',"`vname'","cvc_bootnum", "`m'", `cvcbootnum')
		}
		else if "`est_type'"=="interactive" {	// case 2 - interactive
			tempname pmat1
			cvc `vt1_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			mat `pmat1' = r(pmat)
			mata: add_result_item(`ename',"`vname'","cvc_pval1", "`m'", st_matrix("`pmat1'"))
			// LATE special case - perfect assignment (when cvc should be missing)
			tempname pmat0
			if `perfectflag'==0 {
				cvc `vt0_L_list' if `touse', yvar(`vname') foldvar(`fid')
				mat `pmat0' = r(pmat)
			}
			else {
				mat `pmat0' = `pmat1' * .
			}
			mata: add_result_item(`ename',"`vname'","cvc_pval0", "`m'", st_matrix("`pmat0'"))
			// same bootnum for both
			mata: add_result_item(`ename',"`vname'","cvc_bootnum", "`m'", `cvcbootnum')
		}
		else if "`est_type'"=="fiv" { // case 3 - FIV
			cvc `vt_L_list' if `touse', yvar(`vname') foldvar(`fid') bootnum(`cvcbootnum')
			tempname pmat
			mat `pmat' = r(pmat)
			mata: add_result_item(`ename',"`vname'","cvc_pval", "`m'", st_matrix("`pmat'")) bootnum(`cvcbootnum')
			cvc `vt_h_L_list' if `touse', yvar(`vname') foldvar(`fid')
			tempname pmat_h
			mat `pmat_h' = r(pmat)
			mata: add_result_item(`ename',"`vname'","cvc_pval_h", "`m'", st_matrix("`pmat_h'"))
			// same bootnum for both
			mata: add_result_item(`ename',"`vname'","cvc_bootnum", "`m'", `cvcbootnum')
		}

		// add shortstack results
		if `ssflag' {
			if "`est_type'"=="partial" { // case 1 - partially-linear
		
				label var `shortstack'_ss_`m' "Pred. values E[`vname'|X] using shortstacking, rep `m'"
				
				// save weights as column vector
				mata: add_result_item(`ename', "`shortstack'_ss",  "ss_weights",     "`m'", st_matrix("`ssw'")')
				// base est will be n.a.
				mata: add_learner_item(`ename', "`shortstack'_ss", "stack_base_est", "`base_est'")
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename', "`shortstack'_ss", "ss_final_est",   "`ssfinalest'")
				
				// full sample - mse, rsq, etc.
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE",  "`m'", `mse')
				mata: add_result_item(`ename', "`shortstack'_ss", "RMSE", "`m'", `rmse')
				mata: add_result_item(`ename', "`shortstack'_ss", "R-sq", "`m'", `rsq')
				mata: add_result_item(`ename', "`shortstack'_ss", "N",    "`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list')   \ `N')
				
				// by fold - mse, rsq, etc.				
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`shortstack'_ss", "N_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')   \ r(N))
				
			}
			else if "`est_type'"=="interactive" {	// case 2 - interactive
			
				label var `shortstack'_ss0_`m'  "Pred. values E[`vname'|X] given `treatvar'==0 using shortstacking, rep `m'"
				label var `shortstack'_ss1_`m'  "Pred. values E[`vname'|X] given `treatvar'==1 using shortstacking, rep `m'"
				// save weights as column vector
				mata: add_result_item(`ename',  "`shortstack'_ss", "ss_weights0",  "`m'", st_matrix("`ssw0'")')
				mata: add_result_item(`ename',  "`shortstack'_ss", "ss_weights1",  "`m'", st_matrix("`ssw1'")')
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename', "`shortstack'_ss", "ss_final_est", "`ssfinalest'")
				
				// rsq, mse, N
				rsqmse `shortstack'_ss1_`m' if `treatvar'==1 & `touse', yvar(`vname')
				local mse1			= el(r(mse),1,1)
				local rmse1			= el(r(rmse),1,1)
				local rsq1			= el(r(rsq),1,1)
				local N1			= el(r(N),1,1)
				mata: add_result_item(`ename',"`shortstack'_ss","MSE1",  "`m'", `mse1')
				mata: add_result_item(`ename',"`shortstack'_ss","RMSE1", "`m'", `rmse1')
				mata: add_result_item(`ename',"`shortstack'_ss","R-sq1", "`m'", `rsq1')
				mata: add_result_item(`ename',"`shortstack'_ss","N1",    "`m'", `N1')
				rsqmse `shortstack'_ss0_`m' if `treatvar'==0 & `touse', yvar(`vname')
				local mse0			= el(r(mse),1,1)
				local rmse0			= el(r(rmse),1,1)
				local rsq0			= el(r(rsq),1,1)
				local N0			= el(r(N),1,1)
				mata: add_result_item(`ename',"`shortstack'_ss","MSE0",  "`m'", `mse0')
				mata: add_result_item(`ename',"`shortstack'_ss","RMSE0", "`m'", `rmse0')
				mata: add_result_item(`ename',"`shortstack'_ss","R-sq0", "`m'", `rsq0')
				mata: add_result_item(`ename',"`shortstack'_ss","N0",    "`m'", `N0')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list')   \ `N0')
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list')   \ `N1')
				
				// by fold - mse, rsq, etc.				
				rsqmse `shortstack'_ss1_`m' if `treatvar'==1 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE1_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`shortstack'_ss", "N1_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ r(mse))
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')   \ r(N))
				rsqmse `shortstack'_ss0_`m' if `treatvar'==0 & `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE0_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`shortstack'_ss", "N0_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ r(mse))
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')   \ r(N))
	

			}
			else if "`est_type'"=="fiv" { // case 3 - FIV
	
				label var `shortstack'_ss_`m' "Pred. values E[`vname'|Z,X] using shortstacking, rep `m'"
				label var `shortstack'_h_ss_`m' "Pred. values E[Dhat|X] of `vname' using shortstacking, rep `m'"
				
				// save weights as column vector
				mata: add_result_item(`ename',  "`shortstack'_ss", "ss_weights",   "`m'", st_matrix("`ssw'")')
				mata: add_result_item(`ename',  "`shortstack'_ss", "ss_weights_h", "`m'", st_matrix("`ssw_h'")')
				// base est will be n.a.
				mata: add_learner_item(`ename', "`shortstack'_ss", "stack_base_est", "`base_est'")
				// final estimator used to stack is a learner item; same in "_h" estimation
				mata: add_learner_item(`ename', "`shortstack'_ss", "ss_final_est", "`ssfinalest'")

				// rsq, mse, N
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename',"`shortstack'_ss","MSE",	"`m'", `mse')
				mata: add_result_item(`ename',"`shortstack'_ss","RMSE",	"`m'", `rmse')
				mata: add_result_item(`ename',"`shortstack'_ss","R-sq",	"`m'", `rsq')
				mata: add_result_item(`ename',"`shortstack'_ss","N",	"`m'", `N')
				rsqmse `shortstack'_h_ss_`m' if `touse', yvar(`vname')
				local mse			= el(r(mse),1,1)
				local rmse			= el(r(rmse),1,1)
				local rsq			= el(r(rsq),1,1)
				local N				= el(r(N),1,1)
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_h",  "`m'", `mse')
				mata: add_result_item(`ename',"`shortstack'_ss","RMSE_h", "`m'", `rmse')
				mata: add_result_item(`ename',"`shortstack'_ss","R-sq_h", "`m'", `rsq')
				mata: add_result_item(`ename',"`shortstack'_ss","N_h",    "`m'", `N')
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_list'			= (nullmat(`mse_list')   \ `mse')
				mat `N_list'			= (nullmat(`N_list')     \ `N')
				mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
				mat `N_h_list'			= (nullmat(`N_h_list')   \ `N_h')
				
				// by fold - mse, rsq, etc.				
				rsqmse `shortstack'_ss_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE_folds",   "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`shortstack'_ss", "N_folds",     "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_folds_list'	= (nullmat(`mse_folds_list')   \ r(mse))
				mat `N_folds_list'		= (nullmat(`N_folds_list')     \ r(N))
				rsqmse `shortstack'_h_ss_`m' if `touse', yvar(`vname') foldvar(`fid') kfolds(`kfolds')
				mata: add_result_item(`ename', "`shortstack'_ss", "MSE_h_folds", "`m'", st_matrix("r(mse)"))
				mata: add_result_item(`ename', "`shortstack'_ss", "N_h_folds",   "`m'", st_matrix("r(N)"))
				// returned as r(.) macros by crossfit; not saved in mata struct
				mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ r(mse))
				mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')   \ r(N))
			}
		}
		
		// save results relating to stacked learner if it exists
		if `stdflag' {
			if "`failed'"~="" {
				add_m_col `pysw', flist(`failed') lastcol(`kfolds')
				mat `pysw' = r(A)
				add_m_col `pysm', flist(`failed') lastcol(`kfolds')
				mat `pysm' = r(A)
			}
			mata: add_result_item(`ename',"`vtilde'","stack_weights",   "`m'", st_matrix("`pysw'"))
			mata: add_result_item(`ename',"`vtilde'","stack_MSEs",	  "`m'", st_matrix("`pysm'"))
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est","`base_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est","`stack_final_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_type","`stype'")
		}
		if `stdflag' & "`est_type'"=="fiv" {
			if "`failed_lie_h'"~="" {
				add_m_col `pysw_h', flist(`failed_lie_h') lastcol(`kfolds')
				mat `pysw_h' = r(A)
				add_m_col `pysm_h', flist(`failed_lie_h') lastcol(`kfolds')
				mat `pysm_h' = r(A)
			}
			mata: add_result_item(`ename',"`vtilde'","stack_weights_h",   "`m'", st_matrix("`pysw_h'"))
			mata: add_result_item(`ename',"`vtilde'","stack_MSEs_h",	  "`m'", st_matrix("`pysm_h'"))
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est_h","`base_est_h'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est_h","`stack_final_est_h'")
			mata: add_learner_item(`ename',"`vtilde'","stack_type_h","`stype_h'")
		}
	
	}		// end of resampling loop
	
		
	******************************** RETURN RESULTS ************************************
	if "`est_type'"=="partial" { // case 1 - partially-linear
		foreach mname in mse_list N_list mse_folds_list N_folds_list {
			mat rownames ``mname''	= `rnames'
		}
		return mat mse				= `mse_list'
		return mat N				= `N_list'
		return mat mse_folds		= `mse_folds_list'
		return mat N_folds			= `N_folds_list'
	}
	else if "`est_type'"=="interactive" {
		foreach mname in mse0_list N0_list mse0_folds_list N0_folds_list mse1_list N1_list mse1_folds_list N1_folds_list {
			mat rownames ``mname''	= `rnames'
		}
		return mat mse0				= `mse0_list'
		return mat N0				= `N0_list'
		return mat mse1				= `mse1_list'
		return mat N1				= `N1_list'
		return mat mse0_folds		= `mse0_folds_list'
		return mat mse1_folds		= `mse1_folds_list'
		return mat N0_folds			= `N0_folds_list'
		return mat N1_folds			= `N1_folds_list'
	}
	else if "`est_type'"=="fiv" {
		foreach mname in mse_h_list N_h_list mse_h_folds_list N_h_folds_list {
			mat rownames ``mname''	= `rnames'
		}
		return mat mse_h			= `mse_h_list'
		return mat N_h				= `N_h_list'
		return mat mse_h_folds		= `mse_h_folds_list'
		return mat N_h_folds		= `N_h_folds_list'
	}

	qui count if `touse'
	return scalar N				= r(N)
	return local cmd_list		`cmd_list'
	return local cmd_h_list		`cmd_h_list'
	return scalar perfectflag	= `perfectflag'
	
end

program define check_foldvar, rclass
	syntax [anything], fidlist(varlist) touse(varname) [ NOIsily ]
	if "`noisily'"=="" local qui quietly
	local reps : word count `fidlist'
	tokenize `fidlist'
	forvalues m=1/`reps' {
		local fid ``m''
		// check that fold var is legit
		qui count if `fid'==. & `touse'
		if r(N)>0 {
			di as res "note - fold variable missing for some observations"
			di as res "these observations will be excluded from the estimation sample"
			qui replace `touse' = 0 if `vname'==.
		}
		// enforce that the number of folds is the same for all fold vars
		`qui' tab `fid'
		if `m'==1 {
			// initialize kfolds in resample 1 for checking vs resamplings 2,3,...
			local kfolds = r(r)
			if `kfolds'==1 {
				di as err "error - fold variable identifies only one group"
				exit 198
			}
		}
		else {
			if r(r)~=`kfolds' {
				di as err "error - fold variables must have same number of folds"
				exit 198
			}
		}
	}
	return scalar reps		= `reps'
	return scalar kfolds	= `kfolds'
end

// for interactive use
program define initialize_eqn_info, rclass

	syntax [anything] [if] [in] ,					/// 
							[						///
							ename(name)				/// name of mata struct
							vname(varname)			/// name of dep var to be orthogonalized
							vtlist(string)			/// names of corresponding tilde variables
							shortstack(name)		/// name of short-stacking variable; may be empty
							poolstack(name)			/// name of pooled-stacking variable; may be empty
							estring(string asis)	/// names of estimation strings
													/// need asis option in case it includes strings
							estringh(string asis)	/// names of LIE estimation strings
													/// need asis option in case it includes strings
							vtype(string)			///
							predopt(string asis)	///
							NOIsily					///
							]
	
	if "`noisily'"=="" {
		local qui quietly
	}
	
	// name for temp mata object
	tempname t
	
	parse_estring, vtlist(`vtlist') ename(`ename') estring(`estring') vtype(`vtype') predopt(`predopt') `noisily'
	if "`vtlist'"=="" {
		// parse_estring set the default vtilde names
		local vtlist `r(vtlist)'
	}
	if "`vname'"=="" {
		// parse_estring set the default vname
		local vname `r(vname)'
	}
	
	mata: `ename'.vname = "`vname'"
	mata: `ename'.shortstack = "`shortstack'"	// may be empty string
	mata: `ename'.poolstack = "`poolstack'"		// may be empty string
		
	if "`estringh'"~="" {
		parse_estring, vtlist(`vtlist') ename(`ename') estring(`estringh') h vtype(`vtype') predopt(`predopt') `noisily'
	}
	
	// if there were duplicate vtilde names in vtlist, info for the last learner would be stored in the AA
	// so we remove any duplicates in the vtilde key list
	local vtlist : list uniq vtlist
	// in two steps, to accommodate singleton lists (which are otherwise string scalars and not matrices
	mata: `t' = tokens("`vtlist'")
	mata: `ename'.vtlist	= `t'
	
	local nlearners : word count `vtlist'
	mata: `ename'.nlearners	= `nlearners'

	return scalar nlearners = `nlearners'
	
	// no longer needed so clear from Mata
	mata: mata drop `t'

end

program define parse_estring, rclass

	syntax [anything] [if] [in] ,					/// 
							[						///
							vtlist(namelist)		///
							ename(name)				/// name of mata struct
							estring(string asis)	/// names of estimation strings
													/// need asis option in case it includes strings
							h						/// indicates LIE eqn
							vtype(string)			///
							predopt(string asis)	///
							NOIsily					///
							]
	
	if "`noisily'"=="" {
		local qui quietly
	}

	// used for temporary Mata object
	tempname t
	
	local doparse = 1
	local vnum = 1
	local hasvtlist = "`vtlist'"~=""
	while `doparse' {
		
		tokenize `"`estring'"', parse("||")
		mata: `t' = "`1'"
		// used below with syntax command
		local 0 `"`1'"'
		
		// catch special case - a single | appears inside the estimation string
		if "`2'"=="|" & "`3'"~="|" {
			mata: `t' = "`1' `2' `3'"
			// used below with syntax command
			local 0 `"`1' `2' `3'"'
			mac shift 2
			local estring `*'
			tokenize `"`estring'"', parse("||")
		}		
		
		syntax [anything] [if] [in] , [*]
		local est_main `anything'
		local est_options `options'
		
		if `hasvtlist' {
			local vtilde : word `vnum' of `vtlist'
		}
		else {
			// assign default name
			local vtilde : word 1 of `1'
			local vtilde Y`vnum'_`vtilde'
			local vtlist `vtlist' `vtilde'
		}
		
		if "`h'"=="" {
			local nvname : word 2 of `1'
			if "`vname'"=="" {
				local vname `nvname'
			}
			// check
			if "`vname'"~="`nvname'" {
				di as err "warning - conflicting depvar names, `vname' and `nvname'"
			}
		}
		
		if "`h'"=="" {
			mata: add_learner_item(`ename',"`vtilde'","estring","`0'")
			mata: add_learner_item(`ename',"`vtilde'","est_main","`est_main'")
			mata: add_learner_item(`ename',"`vtilde'","est_options","`est_options'")
			mata: add_learner_item(`ename',"`vtilde'","predopt","`predopt'")
			mata: add_learner_item(`ename',"`vtilde'","vtype","`vtype'")
		}
		else {
			mata: add_learner_item(`ename',"`vtilde'","estring_h","`0'")
			mata: add_learner_item(`ename',"`vtilde'","est_main_h","`est_main'")
			mata: add_learner_item(`ename',"`vtilde'","est_options_h","`est_options'")
			mata: add_learner_item(`ename',"`vtilde'","predopt_h","`predopt'")
			mata: add_learner_item(`ename',"`vtilde'","vtype","`vtype'")
		}

		if "`2'"~="|" & "`3'"~="|" {
			// done parsing
			local doparse = 0
		}

		mac shift 3
		
		local estring `*'
		
		local ++vnum
	}
	
	return local vtlist `vtlist'
	return local vname `vname'
	
	// no longer needed so clear from Mata
	mata: mata drop `t'

end

prog define add_m_col, rclass
	version 16.0
	syntax name, flist(numlist) lastcol(integer)
	
	tempname A A_miss
	mat `A' = `namelist'
	
	mat `A_miss' = J(rowsof(`A'),1,.)
	foreach c in `flist' {
		if `c'==1 {
			mat `A' = `A_miss',`A'
		}
		else if `c'==`lastcol' {
			mat `A' = `A',`A_miss'
		}
		else {
			mat `A' = `A'[.,1..`c'-1], `A_miss', `A'[.,`c'...] 
		}
	}
	return mat A=`A'

end

prog define _fit_pystacked, rclass
	version 16.0
	syntax [anything], [				///
		k(integer 0)					///
		m(integer 0)					///
		nlearners(integer 1)			///
		stdflag(integer 0)				/// flag for doing standard stacking
		cap(string)						///
		qui(string)						///
		est_main(string)				///
		est_options(string)				///
		NOSTDstackopt(string)			///
		stdfinalestopt(string)			///
		vname(varname)					///
		vhat(varname)					/// where to save (replace) stacked OOS (crossfit) predicted values
		vhat_ss(varname)				/// short-stacked version
		vhat_not_k(varname)				/// where to save in-sample predicted values
		vhat_not_k_ss(varname)			/// short-stacked version
		vtilde(varlist)					/// where to save (replace) OOS (crossfit) base learner predicted values
		vtilde_not_k(varlist)			/// where to save in-sample base learner predicted values
		treatvar(varname)				///
		treatval(integer 1)				///
		perfect(integer 0)				///
		fid(varname)					///
		touse(varname)					///
		mmatname(name)					/// name of Mata matrix for storing the in-sample base learner predicted values
		etype(string)					/// will be "D" or blank
						]

	tempname pysw pysm
	tempvar tomata fidtouse stacking_p y_hat toest
	// initialize
	mat `pysw' = .
	mat `pysm' = .
	local failed = 0
	local rc = 0

	// estimation is either all obs excluding the fold ID
	// or also for the treated or untreated sample (treatvar=1 or treatvar=0)
	if "`treatvar'"=="" {
		qui gen `toest' = `fid'!=`k' & `touse'
	}
	else {
		// treatment effects / interactive model
		qui gen `toest' = `fid'!=`k' & `touse' & (`treatvar'==`treatval')
	}
	if `perfect' {
		// create predicted values = treatval (0 or 1) for kth (OOS) fold
		qui replace `vhat' = `treatval' if `fid'==`k' & `touse'
		// standard stacking case
		// OOS (crossfit) base learner predicted values
		tokenize `vtilde'
		forvalues j=1/`nlearners' {
			qui replace ``j'' = `treatval' if `fid'==`k' & `touse'
		}
		if `stdflag' | (`nlearners'==1) {
			// in-sample CV base learner predicted values - temp vars, for mata
			// initialize
			local st_cv_list
			forvalues j=1/`nlearners' {
				tempvar stacking_p_cv`j'
				qui gen double `stacking_p_cv`j'' = `treatval'
				local st_cv_list `st_cv_list' `stacking_p_cv`j''
			}
			qui gen `tomata' = `fid'!=`k' & `touse'
			qui gen `fidtouse' = (`k'*`touse')
			mata: `mmatname' = `mmatname' \ st_data(., "`fid' `fidtouse' `vname' `st_cv_list'", "`tomata'")
		}
		// empty weights and RMSPEs
		mat `pysw' = J(`nlearners',1,.)
		mat `pysm' = J(`nlearners',1,.)
	}
	else {
		// main estimation block
		`cap' `est_main' if `toest', `est_options' `nostdstackopt' `stdfinalestopt'

		local rc = _rc
		if `rc'==0 {
			// estimation successful, _rc=0
			`qui' di as text "N=" as res e(N)
			local base_est `e(base_est)'
			local finalest `e(finalest)'
			// check
			assert "`e(cmd)'"=="pystacked"
			assert e(mcount)<.
			
			// *** pystacked complains if data change after estimation and before predict ***
			// *** so make sure all predict commands are done first ***
			
			// always save crossfit base learner predicted values
			// OOS (crossfit) base learner predicted values
// confirm
			qui predict double `stacking_p' /* if `fid'==`k' & `touse' */, basexb

			if `stdflag' | (`nlearners'==1) {
				// standard stacking case; also enter if just 1 learner
				// get predicted values; used for both OOS and in-sample
				qui predict `vtype' `y_hat' if `touse'
				// in-sample CV base learner predicted values
				tempvar stacking_p_cv
				qui predict double `stacking_p_cv', basexb cv
				// save pystacked weights and MSEs
				qui pystacked, table(rmspe)		// create rmspe matrix
				mat `pysw' = e(weights)
				mat `pysm' = e(rmspe)
				mat `pysm' = `pysm'[2...,"RMSPE_cv".."RMSPE_cv"]
				mata: st_replacematrix("`pysm'",st_matrix("`pysm'"):^2)
			}
			
			// OOS (crossfit) base learner predicted values
			tokenize `vtilde'
			forvalues j=1/`nlearners' {
				qui replace ``j'' = `stacking_p'`j' if `fid'==`k' & `touse'
			}
			// in-sample base learner predicted values (if supplied)
			if "`vtilde_not_k'"~="" {
				tokenize `vtilde_not_k'
				forvalues j=1/`nlearners' {
					qui replace ``j'' = `stacking_p'`j' if `fid'~=`k' & `touse'
				}
			}

			tempvar y_hat_2
			if `nlearners' > 1 {
				qui _ddml_nnls `vname' `stacking_p'* if `fid'!=`k' & `touse'
				qui predict double `y_hat_2'
			}
			else {
				// singler learner case, no NNLS needed. "1" on end added by pystacked (learner number)
				qui gen double `y_hat_2' = `stacking_p'1
			}

			if "`vhat_not_k_ss'"~="" {
				qui replace `vhat_not_k_ss' = `y_hat_2' if `toest'
			}

			// vhat is OOS stacked predicted values for fold k based on
			// pystacked nnls of vname (original D in Step on) on base learner CV predicted values
			// vhat_ss is OOS stacked predicted values for fold k based on
			// nnls of vname (original D in Step 1) on the base learner predicted values.

			if "`vhat_ss'"~="" {
				qui replace `vhat_ss' = `y_hat_2' if `fid'==`k' & `touse'
			}

			if `stdflag' | (`nlearners'==1) {
				// insert previously-created stacked out-of-sample predicted values
				qui replace `vhat' = `y_hat' if `fid'==`k' & `touse'
				// insert in-sample predicted values if requested
				if "`vhat_not_k'"~="" {
					qui replace `vhat_not_k' = `y_hat' if `fid'~=`k' & `touse'
				}
				// in-sample CV base learner predicted values
				// accumulated in mata along with corresponding values of dep var
				qui gen `tomata' = `fid'!=`k' & `touse'
				qui gen `fidtouse' = (`k'*`touse')
				fvexpand `stacking_p_cv'*
				mata: `mmatname' = `mmatname' \ st_data(., "`fid' `fidtouse' `vname' `r(varlist)'", "`tomata'")
			}
		}
		else {
			// estimation failed with _rc>0 - all in-sample and OOS predicted values are missings
			local failed = 1
			// predicted values = missing values for kth (OOS) fold (though should be missing already)
			qui replace `vhat' = . if `fid'==`k'
			// predicted values = missing values for in-sample estimation (though should be missing already)
			if "`vhat_not_k'"~="" {
				qui replace `vhat_not_k' = .
			}
			// OOS (crossfit) base learner predicted values
			tokenize `vtilde'
			forvalues j=1/`nlearners' {
				qui replace ``j'' = . if `fid'==`k' & `touse'
			}
			if `stdflag' | (`nlearners'==1) {
				// standard stacking case
				// in-sample CV base learner predicted values - temp vars, for mata
				// initialize
				local st_cv_list
				forvalues j=1/`nlearners' {
					tempvar stacking_p_cv`j'
					qui gen double `stacking_p_cv`j'' = .
					local st_cv_list `st_cv_list' `stacking_p_cv`j''
				}
				qui gen `tomata' = `fid'!=`k' & `touse'
				qui gen `fidtouse' = (`k'*`touse')
				mata: `mmatname' = `mmatname' \ st_data(., "`fid' `fidtouse' `vname' `st_cv_list'", "`tomata'")
			}
		}
	}

	return local base_est	`base_est'
	return local final_est	`finalest'
	return local rc			= `rc'
	return local failed		= `failed'
	return local perfect	= `perfect'
	return mat pysw			= `pysw'
	return mat pysm			= `pysm'
end


prog define _fit_other, rclass
	version 16.0
	syntax [anything], [				///
		k(integer 0)					///
		m(integer 0)					///
		cap(string)						///
		qui(string)						///
		est_main(string)				///
		est_options(string)				///
		vname(varname)					///
		vhat(varname)					/// where to save (replace) stacked OOS (crossfit) predicted values
		treatvar(varname)				///
		treatval(integer 1)				///
		fid(varname)					///
		touse(varname)					///
						]

	tempvar fidtouse y_hat toest
	// initialize
	local failed = 0
	local perfect = 0
	local rc = 0

	// estimation is either all obs excluding the fold ID
	// or also for the treated or untreated sample (treatvar=1 or treatvar=0)
	if "`treatvar'"=="" {
		qui gen `toest' = `fid'!=`k' & `touse'
	}
	else {
		// treatment effects / interactive model
		qui gen `toest' = `fid'!=`k' & `touse' & (`treatvar'==`treatval')
		// LATE model special cases:
		// nb: D (treatment, `vname') can be continuous; code below aimed at binary treatment
		//	 Z (assignment, `treatvar') is always binary
		// possible that in the Z=1 subsample of the estimation sample, all are assigned to treatment (D=1)
		// possible that in the Z=0 subsample of the estimation sample, no one is assigned to treatment (D=0)
		// main case: perfect assignment to (non-)treatment: D=1 for all obs when Z=1, or D=0 for all obs when Z=0
		// if no variation in D, estimation may fail, so catch here as a special case
		qui count if `touse'
		local N = r(N)
		qui count if `treatvar'==0 & `vname'==0 & `touse'
		local Z0D0 = (r(N)==`N')
		qui count if `treatvar'==1 & `vname'==1 & `touse'
		local Z1D1 = (r(N)==`N')
		local perfect = (`Z0D0' | `Z1D1') & "`etype'"=="D"
	}

	if `perfect' {
		// create predicted values = treatval (0 or 1) for kth (OOS) fold
		qui replace `vhat' = `treatval' if `fid'==`k' & `touse'
	}
	else {

		`cap' `est_main' if `toest', `est_options'

		local rc = _rc
		if `rc'>0 {
			// estimation failed - all in-sample and OOS predicted values are missings
			local failed = 1
			// predicted values = missing values for kth (OOS) fold (though should be missing already)
			qui replace `vhat' = . 
		}
		else {
			`qui' di as text "N=" as res e(N)
			local base_est `e(cmd)'
			qui predict `vtype' `y_hat' if `fid'==`k'
			qui replace `vhat' = `y_hat' if `fid'==`k'
		}
	}
	
	// mse
	tempvar vres vres_sq
	qui gen double `vres' = `vname' - `vhat' if `fid'==`k'
	qui gen double `vres_sq' = `vres'^2
	qui sum `vres_sq'
	
	return scalar mse		=r(mean)
	return local base_est	`e(cmd)'
	return local rc			= `rc'
	return local failed		= `failed'
	return local perfect	= `perfect'

end

program define rsqmse, rclass

	version 16
	syntax varlist [if] [in],			   ///
					yvar(varname)			///
					[						///
					foldvar(varname)		///
					kfolds(integer 0)		///
					]

	marksample touse
	markout `touse' `yvar'
	
	if "`foldvar'"=="" {
		tempvar foldvar
		qui gen byte `foldvar'=1
		local kfolds 1
	}
	
	tempname rsqmat msemat rmsemat Nmat rsqmat_i msemat_i Nmat_i rmsemat_i
	tempvar resid_sq ydm_sq
	qui gen double `resid_sq' = .
	qui gen double `ydm_sq' = .
	
	forvalues k = 1(1)`kfolds' {
		cap mat drop `Nmat_i'
		cap mat drop `msemat_i'
		cap mat drop `rmsemat_i'
		cap mat drop `rsqmat_i'
		foreach v of varlist `varlist' {
			qui replace `resid_sq' = (`yvar' - `v')^2 if `touse' & `foldvar'==`k'
			// should re-calculate mean each time in case #obs of v changes
			sum `yvar' if `touse' & `v'<. & `foldvar'==`k', meanonly
			qui replace `ydm_sq' = (`yvar' - r(mean))^2 if `touse' & `v'<. & `foldvar'==`k'
			sum `resid_sq' if `touse' & `v'<. & `foldvar'==`k', meanonly
			local numerator = r(sum)
			mat `msemat_i' = nullmat(`msemat_i') \ r(mean)
			mat `rmsemat_i' = nullmat(`rmsemat_i') \ sqrt(r(mean))
			mat `Nmat_i' = nullmat(`Nmat_i') \ r(N)
			sum `ydm_sq' if `touse' & `v'<. & `foldvar'==`k', meanonly
			local denominator = r(sum)
			mat `rsqmat_i' = nullmat(`rsqmat_i') \ 1-`numerator'/`denominator'
		}
		mat `Nmat'		= nullmat(`Nmat') , `Nmat_i'
		mat `msemat'	= nullmat(`msemat'), `msemat_i'
		mat `rmsemat'	= nullmat(`rmsemat'), `rmsemat_i'
		mat `rsqmat'	= nullmat(`rsqmat'), `rsqmat_i'
	}
	
	mat rownames `rsqmat'	= `varlist'
	mat rownames `msemat'	= `varlist'
	mat rownames `msemat'	= `varlist'
	mat rownames `Nmat'		= `varlist'
	return matrix rsq		= `rsqmat'
	return matrix mse		= `msemat'
	return matrix rmse		= `rmsemat'
	return matrix N			= `Nmat'
					
end

program define cvc, rclass

	version 16
	syntax varlist [if] [in],			   ///
					yvar(varname)		   ///
					foldvar(varname)		///
					[					   ///
					bootnum(integer 500)	///
					]
	
	marksample touse
	markout `touse' `yvar' `foldvar'
	local resid = "`resid'"!=""
	local nlearners : word count `varlist'

	tempname pmat
	
	if `nlearners' == 1 {
		// trivial case, only one learner supplied, return missing value
		mat `pmat' = .
	}
	else {
	
		foreach v of varlist `varlist' {
		
			local yhat1 `v'
			local yhat2 : list varlist - v
			
			// catch perfect assignment, when yhat1 is all zeros
			qui count if `yhat1'~=0 & `touse'
			if r(N) > 0 {
				mata: cvc_calc("`yvar'","`yhat1'","`yhat2'","`foldvar'","`touse'",`bootnum')
				local pval = r(pval)
			}
			else {
				// all zeros (perfect assignment), cvc pval = missing
				local pval = .
			}
			mat `pmat' = nullmat(`pmat') \ `pval'
	
		}
	}
	
	mat rownames `pmat'   = `varlist'

	return matrix pmat	= `pmat'

end

mata: 
void cvc_calc(   
			string scalar yvar,
			string scalar yhat1,
			string scalar yhat2,
			string scalar foldvar,
			string scalar touse,
			real scalar bootnum
			)
{

	st_view(Y,.,yvar,touse)
	st_view(Yhat1,.,yhat1,touse) 
	st_view(Yhat2,.,yhat2,touse) 
	st_view(fid,.,foldvar,touse) 
	
	Nt=rows(Y)
	N2=cols(Yhat2)

	// fitted values always supplied, not residuals	
	Yhat1 = Y :- Yhat1
	Yhat2 = Y :- Yhat2
	
	loss1 = mean((Yhat1):^2)
	loss2 = mean((Yhat2):^2)
	
	if (N2>1) {
		Yhat1 = Yhat1 :* J(Nt,N2,1)
	}
	// step 1:
	zeta = ((Yhat1):^2) :- ((Yhat2):^2)

	// step 2: calculate mean of zeta by fold
	fid_uni = uniqrows(fid)
	folds = rows(fid_uni)
	zeta_m = J(Nt,N2,.)
	for (j=1;j<=folds;j++) {
		k=fid_uni[j,1]
		sel = selectindex(fid:==k)
		meank = mean(zeta[sel,.])
		zeta_m[sel,.] = J(length(sel), N2, 1) :* meank
	}

	// step 3:
	zeta_til = zeta:-zeta_m

	// step 4:
	zeta_sd = sqrt(diagonal(variance(zeta_til)))'

	// step 5:
	Tx = max(sqrt(Nt)*mean(zeta):/zeta_sd)

	// step 6:
	Txb=J(bootnum,1,.)
	for (b=1; b<=bootnum; b++) {
		bw=rnormal(Nt,1,0,1)
		Txb[b]=max((1/sqrt(Nt))*sum((zeta_til:/zeta_sd):*bw))
	}

	// step 7:
	Pval=mean(Txb:>Tx)

	st_numscalar("r(pval)",Pval)
	st_numscalar("r(Nt)",Nt)
	st_numscalar("r(loss1)",loss1)
	st_matrix("r(loss2)",loss2)
	st_numscalar("r(folds)",folds)
}
end
