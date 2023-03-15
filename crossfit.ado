*! crossfit v0.6
*! last edited: 21feb2023
*! authors: aa/ms
* need to accommodate weights in parsing of estimation strings

program define crossfit, rclass sortpreserve
	// minimum Stata is version 14, with support for associative arrays
	version 14
	syntax [anything] [if] [in] ,					/// [anything] is renamed to vname below; currently undocumented option
							[ estring(string asis)	/// estimation string
							estringh(string asis)	/// est string for E[D^|XZ]		
													/// need asis option in case it includes strings
							ename(name)				/// name for Mata struct; default is "crossfit"
							NOREPLACE				/// use already-initialized eqn in ename
							vtilde(namelist)		/// name(s) of fitted variable(s)
							Generate(namelist)		/// synonym for vtilde
							vtype(string)			/// datatype of fitted variable; default=double
							vtype_h(string)			/// datatype of fitted variable; default=double
							shortstack(name)		/// for interactive use
							poolstack(name)			/// for interactive use
							predopt(string asis)	/// undocumented
							NOIsily					///
							*						/// other options that go to subroutines
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
							vtype_h(`vtype_h')						///
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
		}
	}
	else {
		// ename is a pre-populated eqn struct
		local eqn_info `ename'
	}
	mata: st_local("pystackedmulti", strofreal(`eqn_info'.pystackedmulti))
	
	if `pystackedmulti' {
		// enter single learner = pystacked and number of pystacked base learners > 1
		_crossfit_pystacked	`if' `in',	///
			ename(`eqn_info')			///
			vtype(`vtype')				///
			`noisily'					///
			`cfoptions'			
	}
	else {
		// multiple ddml learners or single learner = pystacked with a single base learner
		_crossfit_other	`if' `in',		///
			ename(`eqn_info')			///
			vtype(`vtype')				///
			`noisily'					///
			`cfoptions'			
	}
	
	return add

end


program define _crossfit_pystacked, rclass sortpreserve
	// minimum Stata is version 16, with support for Python
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
							resid					/// return residuals (prediction errors); default is predicted values
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							NOIsily					///
							allowallzero			/// in the LATE model: allow D
							ssfinalest(name)		/// final estimator for short-stacking
							psfinalest(name)		/// final estimator for pooled-stacking
							finalest(name)			/// final estimator for both
							*						/// ignored options
							]

	// used throughout	
	local cmd pystacked

	// final estimator choice; default is NNLS + coefs sum to 1
	if "`finalest'"==""		local finalest nnls1
	if "`ssfinalest'"==""	local ssfinalest `finalest'
	if "`psfinalest'"==""	local psfinalest `finalest'
	
	mata: st_local("vname", `ename'.vname)
	mata: st_local("vtlist", invtokens(`ename'.vtlist))
	mata: st_local("shortstack", `ename'.shortstack)
	mata: st_local("poolstack", `ename'.poolstack)
	mata: st_local("lieflag", strofreal(`ename'.lieflag))

	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
	** indicator for pooled-stacking
	local psflag	= "`poolstack'"~=""
	** indicator for interactive model
	local tvflag	= "`treatvar'"~=""
	** indicator for residuals vs predicted values
	local residflag	= "`resid'"~=""

	// LIE => we want predicted values not resids
	if `lieflag' & `residflag' {
		di as res "resid option ignored"
		local residflag=0
	}
	
	if "`noisily'"=="" {
		local qui quietly
	}

	*** set sample
	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `vname'
	
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
	if `lieflag' {
		// LIE locals
		mata: st_local("est_main_h", return_learner_item(`ename',"`vtilde'","est_main_h"))
		mata: st_local("est_options_h", return_learner_item(`ename',"`vtilde'","est_options_h"))
		mata: st_local("vtype_h",return_learner_item(`ename',"`vtilde'","vtype_h"))
		mata: st_local("predopt_h",return_learner_item(`ename',"`vtilde'","predopt_h"))
		local cmd_h : word 1 of `est_main_h'
		if "`cmd_h'"=="pystacked" {
			// call pystacked with noestimate option to parse and get basic model specs
			// replace placeholder with depvar for parsing purposes
			local est_main_h_k = subinstr("`est_main_h'","{D}","`vname'",1)
			`qui' di as text "calling pystacked {D} on full sample with noestimate option..."
			`qui' `est_main_h_k' if `touse', `est_options_h' noestimate
			`qui' di as text "N=" as res e(N)
			`qui' di as text "number of learners = `e(mcount)'"
			if e(mcount) > 1 {
				local nlearners_h= e(mcount)
				// will be "reg" or "class"
				local stype_h `e(type)'
				local pystacked_h_flag = 1
			}
			else {
				local pystacked_h_flag = 0
				local nlearners_h = 1
			}
		}
		else {
			local pystacked_h_flag = 0
			local nlearners_h = 1
		}
	}
	
	// stored results; initialize here
	tempname mse_list N_list mse_folds_list N_folds_list
	tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
	tempname mse0_list N0_list mse0_folds_list N0_folds_list
	tempname mse1_list N1_list mse1_folds_list N1_folds_list
	
	// will save base learner CV predictions in mata
	tempname y_stacking_cv y_stacking_cv0 y_stacking_cv1 d_xz_stack_cv dhat_x_stack_cv d_dhat_j_insample
	
	// notes:
	// `vtilde': name for fitted variable (predicted value or residual).
	// `vhat_k': OOS (crossfit) stacked predicted values for fold k.
	// `vhat1_k', `vhat0_k': as above for TV.
	// `vhat': OOS (crossfit) stacked predicted values, for all folds; filled when looping over folds. tempvar.
	// `vres': OOS (crossfit) stacked residuals, for all folds; filled when looping over folds. tempvar.
	// `vtilde'_`m': `vhat' or `vres' using user-provided varname.
	// `vtilde'_L`j'_`m': OOS (crossfit) predicted values, by learner. m is resample j is base learner number.
	// `vtilde_list': list of `vtilde'_L`j'_`m' for j=1,2,...; reset for each m loop.
	// not (yet?) implemented - making resid option apply to `vtilde'_L`j'_`m'
	// `stacking_p_cv'[1,2,...]: in-sample CV predicted values of base learners; saved in mata for poolstacking.
	// `stacking_p'[1,2,...]: predicted values of base learners, in-sample and OOS; OOS => `vtilde'_L`j'_`m'.
	// `shortstack': name for shortstacked fitted variable (predicted value or residual).
	// `shortstack'_ss_`m': shortstacked fitted variable of D for resample m; not a tempvar.
	// `shortstack'_h_ss_`m': shortstacked fitted variable of Dhat for resample m; not a tempvar.
	// `hhat_k': OOS (crossfit) stacked predicted values for Dhat (Step II) for fold K.
	// `hhat': OOS (crossfit) stacked predicted values for Dhat (Step II), all folds. tempvar.
	// `vtilde'_h_`m' = `hhat' using user-provided varname.
	// `vtilde'_h_L`j'_`m': OOS (crossfit) predicted values, by learner. m is resample j is base learner number.
	// `vtilde_h_list': list of `vtilde'_h_L`j'_`m' for j=1,2,...; reset for each m loop.
	
	// loop over resamples (crossfitting, shortstacking, store results)
	forvalues m=`firstrep'/`lastrep' {
	
		// create rnames for renaming the rows of saved matrices
		local rnames `rnames' resample_`m'
		// create list of names of variables created by learners for rep m
		// need to reset for each m
		if ~`tvflag' {
			// PLM and LIE
			local vtilde_list
			forvalues j=1/`nlearners' {
				local vtilde_list `vtilde_list' `vtilde'_L`j'_`m'
			}
		}
		else {
			// interactive only
			local vtilde0_list
			local vtilde1_list
			forvalues j=1/`nlearners' {
				local vtilde0_list `vtilde0_list' `vtilde'0_L`j'_`m'
				local vtilde1_list `vtilde1_list' `vtilde'1_L`j'_`m'
			}
		}
		if `lieflag' {
			// LIE only
			local vtilde_h_list
			forvalues j=1/`nlearners_h' {
				local vtilde_h_list `vtilde_h_list' `vtilde'_h_L`j'_`m'
			}
		}
		******************************** CROSSFITTING ************************************
		
		// cross-fitting fundamentally depends on three cases 
		// case 1: `tvflag'==0 & `lieflag'==0
		// case 2: `tvflag'==1 & `lieflag'==0
		// case 3: `lieflag'==1
		// nb: `tvflag'==1 & `lieflag'==1 is impossible
		
		if `lastrep'>1 {
			di as text "Resample `m'..."
		}

		// list of fold IDs refers to current set of reps (not full set including appended)
		local thisrep = `m'-`firstrep'+1
		local fid : word `thisrep' of `fidlist'

		// create blank fitted variable(s), other initializations for resample m
		if ~`tvflag' & ~`lieflag' { // case 1
			tempvar vhat vres
			qui gen `vtype' `vhat'=.
			qui gen `vtype' `vres'=.
			// base learner predicted values
			// not tempvars; name based on `vtilde' macro; used for poolstacking
			forvalues j=1/`nlearners' {
				cap drop `vtilde'_L`j'_`m'
				qui gen `vtype' `vtilde'_L`j'_`m' = .
			}
		}  
		else if `tvflag' & ~`lieflag' { // case 2
			tempvar vhat0 vres0
			tempvar vhat1 vres1
			qui gen `vtype' `vhat0'=.
			qui gen `vtype' `vres0'=.
			qui gen `vtype' `vhat1'=.
			qui gen `vtype' `vres1'=.
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
		else if `lieflag' { // case 3
			// out-of-sample predicted values for E[D|XZ] 
			tempvar vhat
			qui gen `vtype' `vhat'=.
			// out-of-sample predicted values for E[D|X] 
			tempvar hhat
			qui gen `vtype_h' `hhat'=.  // using learner i in both steps
			// base learner predicted values
			// not tempvars; name based on `vtilde' macro; used for poolstacking
			forvalues j=1/`nlearners' {
				cap drop `vtilde'_L`j'_`m'
				qui gen `vtype' `vtilde'_L`j'_`m' = .
			}
			forvalues j=1/`nlearners_h' {
				cap drop `vtilde'_h_L`j'_`m'
				qui gen `vtype_h' `vtilde'_h_L`j'_`m' = .
			}
		}
		else {
			di as err "internal crossfit error"
			exit 198
		}
		
	
		// will save pystacked weights and MSEs
		tempname pysw pysw0 pysw1 pyswh
		tempname pysw_temp pysw1_temp pysw0_temp pyswh_temp
		tempname pysm pysm0 pysm1 pysmh
		tempname pysm_temp pysm1_temp pysm0_temp pysmh_temp
		tempvar tomata
		qui gen byte `tomata'=.
		
		// crossfit
		di as text "Cross-fitting fold " _c
		forvalues k = 1(1)`kfolds' {
	
			di as text "`k' " _c
			
			// needed for each fold loop
			// predicted values for fold k based on estimation for not-k
			tempvar vhat_k vhat1_k vhat0_k
			
			if ~`tvflag' & ~`lieflag' { // case 1
				
				if (`k'==1) & `psflag' {
					// initialize mata object to hold dep var and in-sample crossfit CV predictions for pooled-stacking
					mata: `y_stacking_cv' = J(0,`nlearners'+1,0)
				}
				
				// estimate excluding kth fold
				`qui' `est_main' if `fid'!=`k' & `touse', `est_options'
				`qui' di as text "N=" as res e(N)
				local base_est `e(base_est)'
				local stack_final_est `e(finalest)'
				// check
				assert "`e(cmd)'"=="pystacked"
				assert e(mcount)>1 & e(mcount)<.
				
				// save pystacked weights and MSEs
				qui pystacked, table(rmspe)		// create rmspe matrix
				if (`k'==1) {
					mat `pysw' = e(weights)
					mat `pysm_temp' = e(rmspe)
					mat `pysm' = `pysm_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
					mata: st_replacematrix("`pysm'",st_matrix("`pysm'"):^2)
				}
				else {
					mat `pysw_temp' = e(weights)
					mat `pysw' = (`pysw',`pysw_temp')
					mat `pysm_temp' = e(rmspe)
					mat `pysm_temp' = `pysm_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
					mata: st_replacematrix("`pysm_temp'",st_matrix("`pysm_temp'"):^2)
					mat `pysm' = (`pysm',`pysm_temp')
				}
	
				// get fitted values and residuals for kth fold	
				qui predict `vtype' `vhat_k' if `fid'==`k' & `touse'
				// get base learner predicted values
				tempvar stacking_p_cv stacking_p
				// in-sample CV base learner predicted values
				qui predict double `stacking_p_cv', basexb cv
				// OOS (crossfit) base learner predicted values
				qui predict double `stacking_p' if `fid'==`k' & `touse', basexb
				// OOS (crossfit) base learner predicted values
				forvalues j=1/`nlearners' {
					qui replace `vtilde'_L`j'_`m' = `stacking_p'`j' if `fid'==`k' & `touse'
				}
	
				// get stacked out-of-sample predicted values
				qui replace `vhat' = `vhat_k' if `fid'==`k' & `touse'
				qui replace `vres' = `vname' - `vhat_k' if `fid'==`k' & `touse'

				if `psflag' {				
					// in-sample CV base learner predicted values
					// accumulated in mata along with corresponding values of dep var
					qui replace `tomata' = `fid'!=`k' & `touse'
					fvexpand `stacking_p_cv'*
					mata: `y_stacking_cv' = `y_stacking_cv' \ st_data(., "`vname' `r(varlist)'", "`tomata'")
				}
			}
	
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
	
				// pystacked learners
				if (`k'==1) & `psflag' {
					// initialize mata object to hold dep var and in-sample crossfit CV predictions for pooled-stacking
					mata: `y_stacking_cv0' = J(0,`nlearners'+1,0)
					mata: `y_stacking_cv1' = J(0,`nlearners'+1,0)
				}
	
				// outcome equation so estimate separately for treated and untreated
	
				// for treatvar = 1
				// estimate excluding kth fold
				`qui' `est_main' if `fid'!=`k' & `treatvar' == 1 & `touse', `est_options'
				`qui' di as text "N=" as res e(N)
				local base_est `e(base_est)'
				local stack_final_est `e(finalest)'
				// check
				assert "`e(cmd)'"=="pystacked"
				assert e(mcount)>1 & e(mcount)<.
				
				// save pystacked weights and MSEs if #learners>1
				qui pystacked, table(rmspe)		// create rmspe matrix
				if (`k'==1) {
					mat `pysw1' = e(weights)
					mat `pysm1_temp' = e(rmspe)
					mat `pysm1' = `pysm1_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
					mata: st_replacematrix("`pysm1'",st_matrix("`pysm1'"):^2)
				}
				else {
					mat `pysw1_temp' = e(weights)
					mat `pysw1' = (`pysw1',`pysw1_temp')
					mat `pysm1_temp' = e(rmspe)
					mat `pysm1_temp' = `pysm1_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
					mata: st_replacematrix("`pysm1_temp'",st_matrix("`pysm1_temp'"):^2)
					mat `pysm1' = (`pysm1',`pysm1_temp')
				}
	
				// get fitted values for kth fold	
				qui predict `vtype' `vhat1_k' if `fid'==`k' & `touse'
				// get base learner predicted values
				tempvar stacking_p_cv stacking_p
				// in-sample CV base learner predicted values
				qui predict double `stacking_p_cv', basexb cv
				// OOS (crossfit) base learner predicted values
				qui predict double `stacking_p' if `fid'==`k' & `touse', basexb
				// OOS (crossfit) base learner predicted values
				forvalues j=1/`nlearners' {
					qui replace `vtilde'1_L`j'_`m' = `stacking_p'`j' if `fid'==`k' & `touse'
				}
	
				// get stacked out-of-sample predicted values
				qui replace `vhat1' = `vhat1_k' if `fid'==`k' & `touse'
				qui replace `vres1' = `vname' - `vhat1_k' if `fid'==`k' & `touse'
				
				if `psflag' {
					// in-sample CV base learner predicted values
					// accumulated in mata along with corresponding values of dep var
					qui replace `tomata' = `fid'!=`k' & `touse'
					fvexpand `stacking_p_cv'*
					mata: `y_stacking_cv1' = `y_stacking_cv1' \ st_data(., "`vname' `r(varlist)'", "`tomata'")
				}
				
				// for treatvar = 0

				// first, we need to account for D always = 0 if Z=0
				// check if dvar is always 0
				qui count if `treatvar'==0 & `vname'!=0 & `fid'==`k' & `touse' 

				if (`r(N)'==0 & "`allowallzero'"!="") {
					// create fitted values = eps (appx 0) for kth fold
					`qui' gen `vhat0_k'=10e-12  if `fid'==`k' & `touse' 
					// OOS (crossfit) base learner predicted values
					forvalues j=1/`nlearners' {
						qui replace `vtilde'0_L`j'_`m' = 10e-12 if `fid'==`k' & `touse'
					}
					// don't need in-sample CV base learner predicted values for poolstacking
				} 
				else {
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `treatvar' == 0 & `touse', `est_options'
					`qui' di as text "N=" as res e(N)
					local base_est `e(base_est)'
					local stack_final_est `e(finalest)'
		
					// save pystacked weights and MSEs if #learners>1
					if e(mcount)>1 & e(mcount)<. {
						qui pystacked, table(rmspe)		// create rmspe matrix
						if (`k'==1) {
							mat `pysw0' = e(weights)
							mat `pysm0_temp' = e(rmspe)
							mat `pysm0' = `pysm0_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
							mata: st_replacematrix("`pysm0'",st_matrix("`pysm0'"):^2)
						}
						else {
							mat `pysw0_temp' = e(weights)
							mat `pysw0' = (`pysw0',`pysw0_temp')
							mat `pysm0_temp' = e(rmspe)
							mat `pysm0_temp' = `pysm0_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
							mata: st_replacematrix("`pysm0_temp'",st_matrix("`pysm0_temp'"):^2)
							mat `pysm0' = (`pysm0',`pysm0_temp')
						}
					}

					// get fitted values for kth fold	
					qui predict `vtype' `vhat0_k' if `fid'==`k' & `touse'
					// get base learner predicted values
					tempvar stacking_p_cv stacking_p
					// in-sample CV base learner predicted values
					qui predict double `stacking_p_cv', basexb cv
					// OOS (crossfit) base learner predicted values
					qui predict double `stacking_p' if `fid'==`k' & `touse', basexb
					// OOS (crossfit) base learner predicted values
					forvalues j=1/`nlearners' {
						qui replace `vtilde'0_L`j'_`m' = `stacking_p'`j' if `fid'==`k' & `touse'
					}
				}
				
				// get stacked out-of-sample predicted values
				qui replace `vhat0' = `vhat0_k' if `fid'==`k' & `touse'
				qui replace `vres0' = `vname' - `vhat0_k' if `fid'==`k' & `touse'
				
				// in-sample CV base learner predicted values
				// accumulated in mata along with corresponding values of dep var
				if `psflag' {
					qui replace `tomata' = `fid'!=`k' & `touse'
					fvexpand `stacking_p_cv'*
					mata: `y_stacking_cv0' = `y_stacking_cv0' \ st_data(., "`vname' `r(varlist)'", "`tomata'")
				}	
			}
	
			else if `lieflag' { // case 3
				
				if (`k'==1) & `psflag' {
					// initialize mata object to hold D and in-sample crossfit CV predictions for pooled-stacking
					mata: `d_xz_stack_cv' = J(0,`nlearners'+1,0)
				}
				if (`k'==1) & `ssflag' {
					// initialize mata object to hold D and in-sample predictions for short-stacking
					mata: `d_dhat_j_insample' = J(0,`nlearners'+1,0)
				}
					
				// Step I: estimation of E[D|XZ]=D^
				// estimate excluding kth fold
				`qui' di
				`qui' di as text "Estimating E[D|XZ]=d(X,Z)"
				`qui' `est_main' if `fid'!=`k' & `touse', `est_options'
				`qui' di as text "N=" as res e(N)
				local base_est `e(base_est)'
				local stack_final_est `e(finalest)'
				// check
				assert "`e(cmd)'"=="pystacked"
				assert e(mcount)>1 & e(mcount)<.
		
				// save pystacked weights and MSEs
				qui pystacked, table(rmspe)		// create rmspe matrix
				if (`k'==1) {
					mat `pysw' = e(weights)
					mat `pysm_temp' = e(rmspe)
					mat `pysm' = `pysm_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
					mata: st_replacematrix("`pysm_temp'",st_matrix("`pysm_temp'"):^2)
				}
				else {
					mat `pysw_temp' = e(weights)
					mat `pysw' = (`pysw',`pysw_temp')
					mat `pysm_temp' = e(rmspe)
					mat `pysm_temp' = `pysm_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
					mata: st_replacematrix("`pysm_temp'",st_matrix("`pysm_temp'"):^2)
					mat `pysm' = (`pysm',`pysm_temp')
				}
				
				// get fitted values (in and out of sample)
				// stores predicted values for E[D|ZX] as a tempvar
				// OOS vhat_k stored below as vhat; in-sample vhat_k used in Step II
				qui predict `vtype' `vhat_k' if `touse'
				// get base learner predicted values
				tempvar stacking_p_cv stacking_p
				// in-sample CV base learner predicted values
				qui predict double `stacking_p_cv', basexb cv
				// crossfit base learner predicted values, OOS (fid=k) and in-sample (fid~=k)
				qui predict double `stacking_p' if `touse', basexb
				// OOS (crossfit) base learner predicted values
				forvalues j=1/`nlearners' {
					qui replace `vtilde'_L`j'_`m' = `stacking_p'`j' if `fid'==`k' & `touse'
				}
	
				// vhat = OOS Dhat stacked predicted values
				qui replace `vhat' = `vhat_k' if `fid'==`k' & `touse'
				
				if `psflag' {
					// in-sample CV base learner predicted values
					// accumulated in mata along with corresponding values of dep var
					qui replace `tomata' = `fid'!=`k' & `touse'
					fvexpand `stacking_p_cv'*
					mata: `d_xz_stack_cv' = `d_xz_stack_cv' \ st_data(., "`vname' `r(varlist)'", "`tomata'")
				}
				if `ssflag' {
					// in-sample base learner predicted values, used for short-stacking
					// accumulated in mata along with corresponding values of D
					qui replace `tomata' = `fid'!=`k' & `touse'
					fvexpand `stacking_p'*
					mata: `d_dhat_j_insample' = `d_dhat_j_insample' \ st_data(., "`vname' `r(varlist)'", "`tomata'")
				}

				// Step II: estimation of E[D^|X]

				// replace {D}-placeholder in estimation string with variable name
				// vhat_k has in-sample and OOS values; estimation restricted to OOS below
				local est_main_h_k = subinstr("`est_main_h'","{D}","`vhat_k'",1)
	
				// estimation	
				`qui' di
				`qui' di as text "Estimating E[D^|X]=h(X)"
				`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
				`qui' di as text "N=" as res e(N)
				local cmd_h `e(cmd)'
				if `pystacked_h_flag' {
					local base_est_h `e(base_est)'
					local stack_final_est_h `e(finalest)'
					// check
					assert "`cmd_h'"=="pystacked"
					assert e(mcount)>1 & e(mcount)<.
					if (`k'==1) & `psflag' {
						// initialize mata object to hold dep var and in-sample crossfit CV predictions for pooled-stacking
						mata: `dhat_x_stack_cv' = J(0,`nlearners_h'+1,0)
					}
					
					// save pystacked weights and MSEs
					qui pystacked, table(rmspe)		// create rmspe matrix
					if (`k'==1) {
						// initialize
						mat `pyswh' = e(weights)
						mat `pysmh_temp' = e(rmspe)
						mat `pysmh' = `pysmh_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
						mata: st_replacematrix("`pysmh'",st_matrix("`pysmh'"):^2)
					}
					else {
						mat `pyswh_temp' = e(weights)
						mat `pyswh' = (`pyswh',`pyswh_temp')
						mat `pysmh_temp' = e(rmspe)
						mat `pysmh_temp' = `pysmh_temp'[2...,"RMSPE_cv".."RMSPE_cv"]
						mata: st_replacematrix("`pysmh_temp'",st_matrix("`pysmh_temp'"):^2)
						mat `pysmh' = (`pysmh',`pysmh_temp')
					}
		
					// get fitted values (in and out of sample)
					tempvar hhat_k // stores predicted values for E[D^|X] temporarily
					qui predict `vtype_h' `hhat_k' if `touse'
					// get base learner predicted values
					tempvar stacking_p_cv_h stacking_p_h
					// in-sample CV base learner predicted values
					qui predict double `stacking_p_cv_h', basexb cv
					// crossfit base learner predicted values, OOS (fid=k) and in-sample (fid~=k)
					qui predict double `stacking_p_h' if `touse', basexb
					// OOS (crossfit) base learner predicted values
					forvalues j=1/`nlearners_h' {
						qui replace `vtilde'_h_L`j'_`m' = `stacking_p_h'`j' if `fid'==`k' & `touse'
					}
	
					// get stacked out-of-sample predicted values
					qui replace `hhat' = `hhat_k' if `fid'==`k' & `touse'
	
					if `psflag' {
						// in-sample CV base learner predicted values
						// accumulated in mata along with corresponding values of {Dhat}
						qui replace `tomata' = `fid'!=`k' & `touse'
						fvexpand `stacking_p_cv_h'*
						mata: `dhat_x_stack_cv' = `dhat_x_stack_cv' \ st_data(., "`vhat' `r(varlist)'", "`tomata'")
					}
				}
				else {

					// get fitted values  
					tempvar hhat_k // stores predicted values for E[D^|X] temporarily
					qui predict `vtype_h' `hhat_k' if `touse', `predopt_h'
		
					// get *combined* out-of-sample predicted values
					qui replace `hhat`j'' = `hhat_k' if `fid'==`k' & `touse'
			
				}
			}
			
			if `k'==1 & `m'==1 {
				local cmd_h_list `cmd_h_list' `cmd_h'
			}
		}
	
		// last fold
		di as text "...completed cross-fitting" _c
		// if noisily, print new line
		`qui' di

		******************************** SHORTSTACKING ************************************

		// shortstacking. we need to distinguish between 3 cases again. 
		if `ssflag' {
	
			if ~`tvflag' & ~`lieflag' { // case 1
				// stack dep var against all OOS (cross-fit) learner predicted values
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (additive model):"
				`qui' get_stack_weights `vname' `vtilde_list' if `touse', finalest(`ssfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw
				mat `ssw' = e(b)
				cap drop `shortstack'_ss_`m'
				mat score `vtype' `shortstack'_ss_`m' = `ssw' if `touse'
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss_`m' = `vname' - `shortstack'_ss_`m'
				}
			}
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=1):"
				`qui' get_stack_weights `vname' `vtilde1_list' if `touse' & `treatvar'==1, finalest(`ssfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw1
				mat `ssw1' = e(b)
				cap drop `shortstack'_ss1_`m'
				mat score `vtype' `shortstack'_ss1_`m' = `ssw1' if `touse'
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss1_`m' = `vname' - `shortstack'_ss1_`m'
				}
				// treatvar == 0
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=0):"
				`qui' get_stack_weights `vname' `vtilde0_list' if `touse' & `treatvar'==0, finalest(`ssfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw0
				mat `ssw0' = e(b)
				cap drop `shortstack'_ss0_`m'
				mat score `vtype' `shortstack'_ss0_`m' = `ssw0' if `touse'
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss0_`m' = `vname' - `shortstack'_ss0_`m'
				}
			}
			else if `lieflag' {
				// apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (LIE, OOS E[D|XZ]):"
				`qui' get_stack_weights `vname' `vtilde_list' if `touse', finalest(`ssfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw
				mat `ssw'= e(b)
				cap drop `shortstack'_ss_`m'
				mat score `vtype' `shortstack'_ss_`m' = `ssw' if `touse'
				if `pystacked_h_flag' {
					// apply short-stacking to in-sample predicted values of E[D|XZ]=D^
					tempname tframe
					qui frame pwf
					local cframe `r(currentframe)'
					frame create `tframe'
					frame change `tframe'
					getmata (`vname' `vtilde_h_list')=`d_dhat_j_insample', force replace
					`qui' di
					`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (LIE, OOS E[D^|X]):"
					`qui' get_stack_weights `vname' `vtilde_h_list', finalest(`ssfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					tempname ssw_h
					mata: `ssw_h' = st_matrix("e(b)")
					frame change `cframe'
					frame drop `tframe'
					mata: st_matrix("`ssw_h'",`ssw_h')
					mat colnames `ssw_h' =  `vtilde_h_list'
					cap drop `shortstack'_h_ss_`m'
					mat score `vtype_h' `shortstack'_h_ss_`m' = `ssw_h' if `touse'
				}
			}
			else {
				di as err "internal crossfit error"
				exit 198
			}
		}
	
		if `ssflag' {
			di as text "...completed short-stacking" _c
		}

		******************************** POOLSTACKING *************************************
		
		if `psflag' {
			// mata object y_stacking_cv has y and predicted yhats of all learners in all crossfits
			if ~`tvflag' & ~`lieflag' { // case 1
				tempname tframe
				qui frame pwf
				local cframe `r(currentframe)'
				frame create `tframe'
				frame change `tframe'
				getmata (`vname' `vtilde_list')=`y_stacking_cv', force replace
				`qui' di
				`qui' di as text "Pooled-stacking, finalest=`psfinalest' (additive model):"
				`qui' get_stack_weights `vname' `vtilde_list', finalest(`psfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname psw
				mata: `psw' = st_matrix("e(b)")
				frame change `cframe'
				frame drop `tframe'
				mata: st_matrix("`psw'",`psw')
				mat colnames `psw' =  `vtilde_list'
				cap drop `poolstack'_ps_`m'
				mat score `vtype' `poolstack'_ps_`m' = `psw' if `touse'
				if `residflag' {
					// vtilde is the residual
					qui replace `poolstack'_ps_`m' = `vname' - `poolstack'_ps_`m'
				}
			}
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
			
				// treatvar=1
				tempname tframe
				qui frame pwf
				local cframe `r(currentframe)'
				frame create `tframe'
				frame change `tframe'
				getmata (`vname' `vtilde1_list')=`y_stacking_cv1', force replace
				`qui' di
				`qui' di as text "Pooled-stacking, finalest=`psfinalest' (additive model):"
				`qui' get_stack_weights `vname' `vtilde1_list', finalest(`psfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname psw1
				mata: `psw1' = st_matrix("e(b)")
				frame change `cframe'
				frame drop `tframe'
				mata: st_matrix("`psw1'",`psw1')
				mat colnames `psw1' = `vtilde1_list'
				cap drop `poolstack'_ps1_`m'
				mat score `vtype' `poolstack'_ps1_`m' = `psw1' if `touse'
				if `residflag' {
					// vtilde is the residual
					qui replace `poolstack'_ps1_`m' = `vname' - `poolstack'_ps1_`m'
				}

				// treatvar=0
				tempname tframe
				qui frame pwf
				local cframe `r(currentframe)'
				frame create `tframe'
				frame change `tframe'
				getmata (`vname' `vtilde0_list')=`y_stacking_cv0', force replace
				`qui' di
				`qui' di as text "Pooled-stacking, finalest=`psfinalest' (additive model):"
				`qui' get_stack_weights `vname' `vtilde0_list', finalest(`psfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname psw0
				mata: `psw0' = st_matrix("e(b)")
				frame change `cframe'
				frame drop `tframe'
				mata: st_matrix("`psw0'",`psw0')
				mat colnames `psw0' =  `vtilde0_list'
				cap drop `poolstack'_ps0_`m'
				mat score `vtype' `poolstack'_ps0_`m' = `psw0' if `touse'
				if `residflag' {
					// vtilde is the residual
					qui replace `poolstack'_ps0_`m' = `vname' - `poolstack'_ps0_`m'
				}
			
			}
			else if `lieflag' {
				tempname tframe
				qui frame pwf
				local cframe `r(currentframe)'
				frame create `tframe'
				
				frame change `tframe'
				getmata (`vname' `vtilde_list')=`d_xz_stack_cv', force replace
				`qui' di
				`qui' di as text "Pooled-stacking, finalest=`psfinalest' (LIE, Step I, D on X,Z):"
				`qui' get_stack_weights `vname' `vtilde_list', finalest(`psfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname psw
				mata: `psw' = st_matrix("e(b)")
				frame change `cframe'
				mata: st_matrix("`psw'",`psw')
				mat colnames `psw' =  `vtilde_list'
				cap drop `poolstack'_ps_`m'
				mat score `vtype' `poolstack'_ps_`m' = `psw' if `touse'
				if `residflag' {
					// vtilde is the residual
					qui replace `poolstack'_ps_`m' = `vname' - `poolstack'_ps_`m'
				}
				if `pystacked_h_flag' {
					frame change `tframe'
					// clear Step I data
					clear
					getmata (`vname' `vtilde_h_list')=`dhat_x_stack_cv', force replace
					`qui' di
					`qui' di as text "Pooled-stacking, finalest=`psfinalest' (LIE, Step II, Dhat on X):"
					`qui' get_stack_weights `vname' `vtilde_h_list', finalest(`psfinalest') stype(`stype') `noisily'
					`qui' di as text "N=" as res e(N)
					tempname psw_h
					mata: `psw_h' = st_matrix("e(b)")
					frame change `cframe'
					frame drop `tframe'
					mata: st_matrix("`psw_h'",`psw_h')
					fvexpand `stacking_p_h'*
					mat colnames `psw' =  `r(varlist)'
					cap drop `poolstack'_h_ps_`m'
					mat score `vtype' `poolstack'_h_ps_`m' = `psw' if `touse'
					if `residflag' {
						// vtilde is the residual
						qui replace `poolstack'_h_ps_`m' = `vname' - `poolstack'_ps_`m'
					}
				}
				else {
					frame drop `tframe'
				}

			}			
		}
		if `psflag' {
			di as text "...completed pooled-stacking" _c
		}
		
		// residualize pystacked base learners (linear models)
		// use general foreach to allow for an empty list
		if `residflag' {
			foreach v in `vtilde_list' `vtilde_h_list' {
				qui replace `v' = `vname' - `v'
			}
		}
		************************************************************************************
		
		// clean up
		cap mata: mata drop `y_stacking_cv'
		cap mata: mata drop `y_stacking_cv0' `y_stacking_cv1'
		cap mata: mata drop `d_xz_stack_cv'
		cap mata: mata drop `dhat_x_stack_cv'
		cap mata: mata drop `d_dhat_j_insample'
		cap mata: mata drop `psw'
		cap mata: mata drop `psw_h'
		cap mata: mata drop `psw1' `psw0'
		cap mata: mata drop `ssw'
		cap mata: mata drop `ssw_h'
		

		***************************** ESTIMATION COMPLETE *********************************
		// estimation done, insert newline
		di
		******************************** STORE RESULTS ************************************

		// vtilde, mspe, etc.
		if ~`tvflag' & ~`lieflag' {
	
			cap drop `vtilde'_`m'
			if "`resid'"=="" {
				// vtilde is predicted values
				qui gen `vtype' `vtilde'_`m' = `vhat'
				qui label var `vtilde'_`m' "Predicted values cond. exp. of `vname' using `cmd'"
				if `ssflag' qui label var `shortstack'_ss_`m' "Predicted values cond. exp. of `vname' using shortstacking"
				if `psflag' qui label var `poolstack'_ps_`m' "Predicted values cond. exp. of `vname' using poolstacking"
			}
			else {
				// vtilde is residuals
				qui gen `vtype' `vtilde'_`m' = `vres'
				qui label var `vtilde'_`m' "Residuals cond. exp. of `vname' using `cmd'"
				if `ssflag' qui label var `shortstack'_ss_`m' "Residuals cond. exp. of `vname' using shortstacking"
				if `psflag' qui label var `poolstack'_ps_`m' "Residuals cond. exp. of `vname' using poolstacking"
			}
			forvalues j=1/`nlearners' {
				qui label var `vtilde'_L`j'_`m' "Predicted values cond. exp. of `vname' using base learner `j'"
			}
	
			// calculate and return mspe and sample size
			tempvar vres_sq
			qui gen double `vres_sq' = `vres'^2 if `touse'

			// additive-type model
			qui sum `vres_sq' if `touse', meanonly
			local mse			= r(mean)
			local N				= r(N)
			tempname mse_folds N_folds
			forvalues k = 1(1)`kfolds' {
				qui sum `vres_sq' if `touse' & `fid'==`k', meanonly
				mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
				qui count if `touse' & `fid'==`k' & `vres_sq'<.
				mat `N_folds' = (nullmat(`N_folds'), r(N))
			}
		
			mat `mse_list'			= (nullmat(`mse_list') \ `mse')
			mat `N_list'			= (nullmat(`N_list') \ `N')
			mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
			mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')
		
			mata: add_result_item(`ename',"`vtilde'","N",         "`m'", `N')
			mata: add_result_item(`ename',"`vtilde'","N_folds",   "`m'", st_matrix("`N_folds'"))
			mata: add_result_item(`ename',"`vtilde'","MSE",       "`m'", `mse')
			mata: add_result_item(`ename',"`vtilde'","MSE_folds", "`m'", st_matrix("`mse_folds'"))
			
			// weights and MSEs will be missing values if #learners=1
			mata: add_result_item(`ename',"`vtilde'","stack_weights",   "`m'", st_matrix("`pysw'"))
			mata: add_result_item(`ename',"`vtilde'","stack_MSEs",      "`m'", st_matrix("`pysm'"))
			
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est","`base_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est","`stack_final_est'")
			// only one learner so it's the opt			
			local mse_opt		= `mse'
			mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
			
		}
		else if `tvflag' & ~`lieflag' {
		
			cap drop `vtilde'0_`m'
			cap drop `vtilde'1_`m'
			if "`resid'"=="" {
				// vtilde is predicted values
				qui gen `vtype' `vtilde'0_`m' = `vhat0'
				qui gen `vtype' `vtilde'1_`m' = `vhat1'
				qui label var `vtilde'0_`m' "Predicted values cond. exp. of `vname' given `treatvar'==0 using `cmd'"
				qui label var `vtilde'1_`m' "Predicted values cond. exp. of `vname' given `treatvar'==1 using `cmd'"
				if `ssflag' qui label var `shortstack'_ss0_`m'  "Predicted values cond. exp. of `vname' given `treatvar'==0 using shortstacking"
				if `ssflag' qui label var `shortstack'_ss1_`m'  "Predicted values cond. exp. of `vname' given `treatvar'==1 using shortstacking"
				if `psflag' qui label var `poolstack'_ps0_`m'  "Predicted values cond. exp. of `vname' given `treatvar'==0 using poolstacking"
				if `psflag' qui label var `poolstack'_ps1_`m'  "Predicted values cond. exp. of `vname' given `treatvar'==1 using poolstacking"
			}
			else {
				// vtilde is residuals
				qui gen `vtype' `vtilde'0_`m' = `vres0'
				qui gen `vtype' `vtilde'1_`m' = `vres1'
				qui label var `vtilde'0_`m' "Residuals cond. exp. of `vname' given `treatvar'==0 using `cmd'"
				qui label var `vtilde'1_`m' "Residuals cond. exp. of `vname' given `treatvar'==1 using `cmd'"
				if `ssflag' qui label var `shortstack'_ss0_`m'  "Residuals cond. exp. of `vname' given `treatvar'==0 using shortstacking"
				if `ssflag' qui label var `shortstack'_ss1_`m'  "Residuals cond. exp. of `vname' given `treatvar'==1 using shortstacking"
				if `psflag' qui label var `poolstack'_ps0_`m'  "Residuals cond. exp. of `vname' given `treatvar'==0 using poolstacking"
				if `psflag' qui label var `poolstack'_ps1_`m'  "Residuals cond. exp. of `vname' given `treatvar'==1 using poolstacking"
			}
			forvalues j=1/`nlearners' {
				qui label var `vtilde'0_L`j'_`m' "Predicted values cond. exp. of `vname' given `treatvar'==0 using base learner `j'"
				qui label var `vtilde'1_L`j'_`m' "Predicted values cond. exp. of `vname' given `treatvar'==1 using base learner `j'"
			}
	
			// calculate and return mspe and sample size
			tempvar vres0_sq vres1_sq
			// vtilde has fitted values
			qui gen double `vres0_sq' = (`vname' - `vhat0')^2 if `treatvar' == 0 & `touse'
			qui gen double `vres1_sq' = (`vname' - `vhat1')^2 if `treatvar' == 1 & `touse'
	
			// interactive-type model, return mse separately for treatvar =0 and =1
			qui sum `vres0_sq' if `treatvar' == 0 & `touse', meanonly
			local mse0			= r(mean)
			local N0			= r(N)
			qui sum `vres1_sq' if `treatvar' == 1 & `touse', meanonly
			local mse1			= r(mean)
			local N1			= r(N)
			local N				= `N0'+`N1'
			tempname mse0_folds N0_folds mse1_folds N1_folds
			forvalues k = 1(1)`kfolds' {
				qui sum `vres0_sq' if `treatvar' == 0 & `touse' & `fid'==`k', meanonly
				mat `mse0_folds' = (nullmat(`mse0_folds'), r(mean))
				qui sum `vres1_sq' if `treatvar' == 1 & `touse' & `fid'==`k', meanonly
				mat `mse1_folds' = (nullmat(`mse1_folds'), r(mean))
				qui count if `treatvar' == 0 & `touse' & `fid'==`k' & `vres0_sq'<.
				mat `N0_folds' = (nullmat(`N0_folds'), r(N))
				qui count if `treatvar' == 1 & `touse' & `fid'==`k' & `vres1_sq'<.
				mat `N1_folds' = (nullmat(`N1_folds'), r(N))
			}

			mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
			mat `N0_list'			= (nullmat(`N0_list') \ `N0')
			mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ `mse0_folds')
			mat `N0_folds_list'		= (nullmat(`N0_folds_list')\ `N0_folds')
			
			mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
			mat `N1_list'			= (nullmat(`N1_list') \ `N1')
			mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ `mse1_folds')
			mat `N1_folds_list'		= (nullmat(`N1_folds_list')\ `N1_folds')

			forvalues t=0/1 {
				mata: add_result_item(`ename',"`vtilde'","N`t'",         "`m'", `N`t'')
				mata: add_result_item(`ename',"`vtilde'","N`t'_folds",   "`m'", st_matrix("`N`t'_folds'"))
				mata: add_result_item(`ename',"`vtilde'","MSE`t'",       "`m'", `mse`t'')
				mata: add_result_item(`ename',"`vtilde'","MSE`t'_folds", "`m'", st_matrix("`mse`t'_folds'"))
			}
			
			// weights and MSEs will be missing values if #learners=1
			mata: add_result_item(`ename',"`vtilde'","stack_weights0","`m'", st_matrix("`pysw0'"))
			mata: add_result_item(`ename',"`vtilde'","stack_weights1","`m'", st_matrix("`pysw1'"))
			mata: add_result_item(`ename',"`vtilde'","stack_MSEs0","`m'",    st_matrix("`pysm0'"))
			mata: add_result_item(`ename',"`vtilde'","stack_MSEs1","`m'",    st_matrix("`pysm1'"))
			
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est","`base_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est","`stack_final_est'")
			forvalues t=0/1 {
				local mse`t'_opt		= `mse`t''
				mata: add_learner_item(`ename',"opt`t'","`m'","`vtilde'")
			}
		}
		else if `lieflag' {

			cap drop `vtilde'_`m'
			cap drop `vtilde'_h_`m'
			
			qui gen `vtype' `vtilde'_`m' = `vhat'
			qui label var `vtilde'_`m' "Predicted values E[`vname'|X,Z]"
			qui gen `vtype' `vtilde'_h_`m' = `hhat'
			qui label var `vtilde'_h_`m' "Predicted values E[`vtilde'_`m'|X]"
			if `ssflag' qui label var `shortstack'_ss_`m'  "Predicted values E[`vname'|X,Z] using shortstacking"
			if `psflag' qui label var `poolstack'_ps_`m'  "Predicted values E[`vname'|X,Z] using poolstacking"
			forvalues j=1/`nlearners' {
				qui label var `vtilde'_L`j'_`m' "Predicted values E[`vname'|X,Z] using base learner `j'"
			}
			if `pystacked_h_flag' {
				if `ssflag' qui label var `shortstack'_h_ss_`m'  "Predicted values E[`vtilde'_`m'|X] using shortstacking"
				if `psflag' qui label var `poolstack'_h_ps_`m'  "Predicted values E[`vtilde'_`m'|X] using poolstacking"
				forvalues j=1/`nlearners_h' {
					qui label var `vtilde'_h_L`j'_`m' "Predicted values Predicted values E[`vtilde'_`m'|X] using base learner `j'"
				}
			}
			// calculate and return mspe and sample size
			tempvar hres dres hres_sq dres_sq
			// vtilde has fitted values
			qui gen double `dres_sq' = (`vname' - `vhat')^2 if `touse'
			qui gen double `hres_sq' = (`vname' - `hhat')^2 if `touse'

			qui sum `dres_sq' if `touse', meanonly
			local mse			= r(mean)
			local N				= r(N)
			tempname mse_folds N_folds
			forvalues k = 1(1)`kfolds' {
				qui sum `dres_sq' if `touse' & `fid'==`k', meanonly
				mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
				qui count if `touse' & `fid'==`k' & `dres_sq'<.
				mat `N_folds' = (nullmat(`N_folds'), r(N))
			}
			mat `mse_list'			= (nullmat(`mse_list') \ `mse')
			mat `N_list'			= (nullmat(`N_list') \ `N')
			mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
			mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')
			
			mata: add_result_item(`ename',"`vtilde'","N",         "`m'", `N')
			mata: add_result_item(`ename',"`vtilde'","N_folds",   "`m'", st_matrix("`N_folds'"))
			mata: add_result_item(`ename',"`vtilde'","MSE",       "`m'", `mse')
			mata: add_result_item(`ename',"`vtilde'","MSE_folds", "`m'", st_matrix("`mse_folds'"))

			qui sum `hres_sq' if `touse', meanonly
			local mse_h			= r(mean)
			local N_h			= r(N)	
			tempname mse_h_folds N_h_folds
			forvalues k = 1(1)`kfolds' {
				qui sum `hres_sq' if `touse' & `fid'==`k', meanonly
				mat `mse_h_folds' = (nullmat(`mse_h_folds'), r(mean))
				qui count if `touse' & `fid'==`k' & `hres_sq'<.
				mat `N_h_folds' = (nullmat(`N_h_folds'), r(N))
			}
			mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
			mat `N_h_list'			= (nullmat(`N_h_list') \ `N_h')
			mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ `mse_h_folds')
			mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')\ `N_h_folds')
			
			mata: add_result_item(`ename',"`vtilde'","N_h",	        "`m'", `N_h')
			mata: add_result_item(`ename',"`vtilde'","N_h_folds",   "`m'", st_matrix("`N_h_folds'"))
			mata: add_result_item(`ename',"`vtilde'","MSE_h",       "`m'", `mse_h')
			mata: add_result_item(`ename',"`vtilde'","MSE_h_folds", "`m'", st_matrix("`mse_h_folds'"))
			
			// weights and MSEs will be missing values if #learners=1
			mata: add_result_item(`ename',"`vtilde'","stack_weights",   "`m'", st_matrix("`pysw'"))
			mata: add_result_item(`ename',"`vtilde'","stack_MSEs",      "`m'", st_matrix("`pysm'"))
			mata: add_result_item(`ename',"`vtilde'","stack_weights_h", "`m'", st_matrix("`pyswh'"))
			mata: add_result_item(`ename',"`vtilde'","stack_MSEs_h",    "`m'", st_matrix("`pysmh'"))
			
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est",   "`base_est'")					
			mata: add_learner_item(`ename',"`vtilde'","stack_base_est_h", "`base_est_h'")					
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est",  "`stack_final_est'")
			mata: add_learner_item(`ename',"`vtilde'","stack_final_est_h","`stack_final_est_h'")
			// optimal D and H
			local mse_opt		= `mse'
			mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
			local mse_h_opt		= `mse_h'
			mata: add_learner_item(`ename',"opt_h","`m'","`vtilde'_h")
		}

		
		// add shortstack results
		if `ssflag' {
			if ~`tvflag' & ~`lieflag' {	// case 1
		
				// calculate and return mspe and sample size
				tempvar vres_sq
				// shortstack macro has the residuals
				qui gen double `vres_sq' = (`shortstack'_ss_`m')^2 if `touse'
			
				// additive-type model
				qui sum `vres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `vres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
			
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')

				mata: add_result_item(`ename',"`shortstack'_ss","N",            "`m'", `N')
				mata: add_result_item(`ename',"`shortstack'_ss","N_folds",      "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","MSE",          "`m'", `mse')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_folds",    "`m'", st_matrix("`mse_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights",   "`m'", st_matrix("`ssw'"))
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
			}
			else if `tvflag' & ~`lieflag' {	// case 2
			
				// calculate and return mspe and sample size
				tempvar vres0_sq vres1_sq
				// shortstack macros have fitted values
				qui gen double `vres0_sq' = (`shortstack'_ss0_`m')^2 if `treatvar' == 0 & `touse'
				qui gen double `vres1_sq' = (`shortstack'_ss1_`m')^2 if `treatvar' == 1 & `touse'
		
				// interactive-type model, return mse separately for treatvar =0 and =1
				qui sum `vres0_sq' if `treatvar' == 0 & `touse', meanonly
				local mse0			= r(mean)
				local N0			= r(N)
				qui sum `vres1_sq' if `treatvar' == 1 & `touse', meanonly
				local mse1			= r(mean)
				local N1			= r(N)
				local N				= `N0'+`N1'
				tempname mse0_folds N0_folds mse1_folds N1_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres0_sq' if `treatvar' == 0 & `touse' & `fid'==`k', meanonly
					mat `mse0_folds' = (nullmat(`mse0_folds'), r(mean))
					qui sum `vres1_sq' if `treatvar' == 1 & `touse' & `fid'==`k', meanonly
					mat `mse1_folds' = (nullmat(`mse1_folds'), r(mean))
					qui count if `treatvar' == 0 & `touse' & `fid'==`k' & `vres1_sq'<.
					mat `N0_folds' = (nullmat(`N0_folds'), r(N))
					qui count if `treatvar' == 1 & `touse' & `fid'==`k' & `vres0_sq'<.
					mat `N1_folds' = (nullmat(`N1_folds'), r(N))
				}
	
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list') \ `N0')
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ `mse0_folds')
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')\ `N0_folds')
				
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list') \ `N1')
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ `mse1_folds')
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')\ `N1_folds')
				
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`shortstack'_ss","N`t'",		        "`m'", `N`t'')
					mata: add_result_item(`ename',"`shortstack'_ss","N`t'_folds",       "`m'", st_matrix("`N`t'_folds'"))
					mata: add_result_item(`ename',"`shortstack'_ss","MSE`t'",           "`m'", `mse`t'')
					mata: add_result_item(`ename',"`shortstack'_ss","MSE`t'_folds",     "`m'", st_matrix("`mse`t'_folds'"))
					mata: add_result_item(`ename',"`shortstack'_ss","ss_weights`t'",    "`m'", st_matrix("`ssw`t''"))
				}
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
			}
			else if `lieflag' {	// case 3
	
				// calculate and return mspe and sample size
				tempvar hres dres hres_sq dres_sq
				// vtilde has fitted values
				qui gen double `dres_sq' = (`vname' - `shortstack'_ss_`m')^2 if `touse'
				if `pystacked_h_flag' qui gen double `hres_sq' = (`vname' - `shortstack'_h_ss_`m')^2 if `touse'
	
				qui sum `dres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `dres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `dres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')

				mata: add_result_item(`ename',"`shortstack'_ss","N",            "`m'", `N')
				mata: add_result_item(`ename',"`shortstack'_ss","N_folds",      "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","MSE",          "`m'", `mse')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_folds",    "`m'", st_matrix("`mse_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights",   "`m'", st_matrix("`ssw'"))
				// final estimator used to stack is a learner item; same in "_h" estimation
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
	
				if `pystacked_h_flag' {
					qui sum `hres_sq' if `touse', meanonly
					local mse_h			= r(mean)
					local N_h			= r(N)	
					tempname mse_h_folds N_h_folds
					forvalues k = 1(1)`kfolds' {
						qui sum `hres_sq' if `touse' & `fid'==`k', meanonly
						mat `mse_h_folds' = (nullmat(`mse_h_folds'), r(mean))
						qui count if `touse' & `fid'==`k' & `hres_sq'<.
						mat `N_h_folds' = (nullmat(`N_h_folds'), r(N))
					}
					mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
					mat `N_h_list'			= (nullmat(`N_h_list') \ `N_h')
					mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ `mse_h_folds')
					mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')\ `N_h_folds')
	
					mata: add_result_item(`ename',"`shortstack'_ss","N_h",          "`m'", `N_h')
					mata: add_result_item(`ename',"`shortstack'_ss","N_h_folds",    "`m'", st_matrix("`N_h_folds'"))
					mata: add_result_item(`ename',"`shortstack'_ss","MSE_h",        "`m'", `mse_h')
					mata: add_result_item(`ename',"`shortstack'_ss","MSE_h_folds",  "`m'", st_matrix("`mse_h_folds'"))
					mata: add_result_item(`ename',"`shortstack'_ss","ss_weights_h", "`m'", st_matrix("`ssw_h'"))
				}
			}
		}
	
		// add poolstack results
		if `psflag' {
			if ~`tvflag' & ~`lieflag' {	// case 1
		
				// calculate and return mspe and sample size
				tempvar vres_sq
				// poolstack macro has the residuals
				qui gen double `vres_sq' = (`poolstack'_ps_`m')^2 if `touse'
			
				// additive-type model
				qui sum `vres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `vres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
			
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')

				mata: add_result_item(`ename',"`poolstack'_ps","N",             "`m'", `N')
				mata: add_result_item(`ename',"`poolstack'_ps","N_folds",       "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`poolstack'_ps","MSE",           "`m'", `mse')
				mata: add_result_item(`ename',"`poolstack'_ps","MSE_folds",     "`m'", st_matrix("`mse_folds'"))
				mata: add_result_item(`ename',"`poolstack'_ps","ps_weights",    "`m'", st_matrix("`psw'"))
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename',"`poolstack'_ps","ps_final_est", "`psfinalest'")
				
			}
			else if `tvflag' & ~`lieflag' {	// case 2
			
				// calculate and return mspe and sample size
				tempvar vres0_sq vres1_sq
				// poolstack macros have fitted values
				qui gen double `vres0_sq' = (`poolstack'_ps0_`m')^2 if `treatvar' == 0 & `touse'
				qui gen double `vres1_sq' = (`poolstack'_ps1_`m')^2 if `treatvar' == 1 & `touse'
		
				// interactive-type model, return mse separately for treatvar =0 and =1
				qui sum `vres0_sq' if `treatvar' == 0 & `touse', meanonly
				local mse0			= r(mean)
				local N0			= r(N)
				qui sum `vres1_sq' if `treatvar' == 1 & `touse', meanonly
				local mse1			= r(mean)
				local N1			= r(N)
				local N				= `N0'+`N1'
				tempname mse0_folds N0_folds mse1_folds N1_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres0_sq' if `treatvar' == 0 & `touse' & `fid'==`k', meanonly
					mat `mse0_folds' = (nullmat(`mse0_folds'), r(mean))
					qui sum `vres1_sq' if `treatvar' == 1 & `touse' & `fid'==`k', meanonly
					mat `mse1_folds' = (nullmat(`mse1_folds'), r(mean))
					qui count if `treatvar' == 0 & `touse' & `fid'==`k' & `vres1_sq'<.
					mat `N0_folds' = (nullmat(`N0_folds'), r(N))
					qui count if `treatvar' == 1 & `touse' & `fid'==`k' & `vres0_sq'<.
					mat `N1_folds' = (nullmat(`N1_folds'), r(N))
				}
	
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list') \ `N0')
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ `mse0_folds')
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')\ `N0_folds')
				
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list') \ `N1')
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ `mse1_folds')
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')\ `N1_folds')
				
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`poolstack'_ps","N`t'",          "`m'", `N`t'')
					mata: add_result_item(`ename',"`poolstack'_ps","N`t'_folds",    "`m'", st_matrix("`N`t'_folds'"))
					mata: add_result_item(`ename',"`poolstack'_ps","MSE`t'",        "`m'", `mse`t'')
					mata: add_result_item(`ename',"`poolstack'_ps","MSE`t'_folds",  "`m'", st_matrix("`mse`t'_folds'"))
					mata: add_result_item(`ename',"`poolstack'_ps","ps_weights`t'", "`m'", st_matrix("`psw`t''"))
				}
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename',"`poolstack'_ps","ps_final_est", "`psfinalest'")
			}
			else if `lieflag' {	// case 3
	
				// calculate and return mspe and sample size
				tempvar hres dres hres_sq dres_sq
				// vtilde has fitted values
				qui gen double `dres_sq' = (`vname' - `poolstack'_ps_`m')^2 if `touse'
				if `pystacked_h_flag' qui gen double `hres_sq' = (`vname' - `poolstack'_h_ps_`m')^2 if `touse'
	
				qui sum `dres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `dres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `dres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')

				mata: add_result_item(`ename',"`poolstack'_ps","N",             "`m'", `N')
				mata: add_result_item(`ename',"`poolstack'_ps","N_folds",       "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`poolstack'_ps","MSE",           "`m'", `mse')
				mata: add_result_item(`ename',"`poolstack'_ps","MSE_folds",     "`m'", st_matrix("`mse_folds'"))
				mata: add_result_item(`ename',"`poolstack'_ps","ps_weights",    "`m'", st_matrix("`psw'"))
				// final estimator used to stack is a learner item; same in "_h" estimation
				mata: add_learner_item(`ename',"`poolstack'_ps","ps_final_est", "`psfinalest'")
				
				if `pystacked_h_flag' {
					qui sum `hres_sq' if `touse', meanonly
					local mse_h			= r(mean)
					local N_h			= r(N)	
					tempname mse_h_folds N_h_folds
					forvalues k = 1(1)`kfolds' {
						qui sum `hres_sq' if `touse' & `fid'==`k', meanonly
						mat `mse_h_folds' = (nullmat(`mse_h_folds'), r(mean))
						qui count if `touse' & `fid'==`k' & `hres_sq'<.
						mat `N_h_folds' = (nullmat(`N_h_folds'), r(N))
					}
					mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
					mat `N_h_list'			= (nullmat(`N_h_list') \ `N_h')
					mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ `mse_h_folds')
					mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')\ `N_h_folds')
	
					mata: add_result_item(`ename',"`poolstack'_ps","N_h",           "`m'", `N_h')
					mata: add_result_item(`ename',"`poolstack'_ps","N_h_folds",     "`m'", st_matrix("`N_h_folds'"))
					mata: add_result_item(`ename',"`poolstack'_ps","MSE_h",         "`m'", `mse_h')
					mata: add_result_item(`ename',"`poolstack'_ps","MSE_h_folds",   "`m'", st_matrix("`mse_h_folds'"))
					mata: add_result_item(`ename',"`poolstack'_ps","ps_weights_h",  "`m'", st_matrix("`psw_h'"))
				}
			}
		}
	}		// end of resampling loop
	
		
	******************************** RETURN RESULTS ************************************

	if ~`tvflag' {
		foreach mname in mse_list N_list mse_folds_list N_folds_list {
			mat rownames ``mname''	= `rnames'
			return mat `mname'		= ``mname''
		}
		return mat mse_folds		= `mse_folds'
		return mat N_folds			= `N_folds'
		return scalar mse			= `mse'
	}
	else {
		foreach mname in mse0_list N0_list mse0_folds_list N0_folds_list mse1_list N1_list mse1_folds_list N1_folds_list {
			mat rownames ``mname''	= `rnames'
			return mat `mname'		= ``mname''
		}
		return scalar mse0			= `mse0'
		return scalar N0			= `N0'
		return scalar mse1			= `mse1'
		return scalar N1			= `N1'
		return mat mse0_folds		= `mse0_folds'
		return mat mse1_folds		= `mse1_folds'
		return mat N0_folds			= `N0_folds'
		return mat N1_folds			= `N1_folds'
	}
	
	if `lieflag' {
		foreach mname in mse_h_list N_h_list mse_h_folds_list N_h_folds_list {
			mat rownames ``mname''	= `rnames'
			return mat `mname'		= ``mname''
		}
		return scalar N_h			= `N_h'
		return mat mse_h_folds		= `mse_h_folds'
		return mat N_h_folds		= `N_h_folds'
		return scalar mse_h			= `mse_h'
	}

	return scalar N			= `N'
	return local cmd_list	pystacked
	return local cmd_h_list	`cmd_h_list'
	
end


program define _crossfit_other, rclass sortpreserve
	// minimum Stata is version 14, with support for associative arrays
	version 14
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
							resid					/// return residuals (prediction errors); default is predicted values
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							NOIsily					///
							allowallzero			/// in the LATE model: allow D
							ssfinalest(name)		/// final estimator for short-stacking
							psfinalest(name)		/// final estimator for pooled-stacking
							finalest(name)			/// final estimator for both
							*						/// ignored options
							]

	// final estimator choice; default is NNLS + coefs sum to 1
	if "`finalest'"==""		local finalest nnls1
	if "`ssfinalest'"==""	local ssfinalest `finalest'
	if "`psfinalest'"==""	local psfinalest `finalest'
	
	mata: st_local("vname", `ename'.vname)
	mata: st_local("vtlist", invtokens(`ename'.vtlist))
	mata: st_local("nlearners", strofreal(`ename'.nlearners))
	mata: st_local("shortstack", `ename'.shortstack)
	mata: st_local("lieflag", strofreal(`ename'.lieflag))

	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
	** indicator for interactive model
	local tvflag	= "`treatvar'"~=""
	** indicator for residuals vs predicted values
	local residflag	= "`resid'"~=""

	// LIE => we want predicted values not resids
	if `lieflag' & `residflag' {
		di as res "resid option ignored"
		local resid = 0
	}
	
	if "`noisily'"=="" {
		local qui quietly
	}
	
	*** set sample
	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `vname'
	
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
	if `ssflag' & `nlearners'==1 {
		di as res "warning - shortstack option relevant only for multiple learners"
	}
	
	// stored results; initialize here
	tempname mse_list N_list mse_folds_list N_folds_list
	tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
	tempname mse0_list N0_list mse0_folds_list N0_folds_list
	tempname mse1_list N1_list mse1_folds_list N1_folds_list
		
	// loop over resamples (crossfitting, shortstacking, store results)
	forvalues m=`firstrep'/`lastrep' {
	
		// create rnames for renaming the rows of saved matrices
		local rnames `rnames' resample_`m'
	
		******************************** CROSSFITTING ************************************
		
		// cross-fitting fundamentally depends on three cases 
		// case 1: `tvflag'==0 & `lieflag'==0
		// case 2: `tvflag'==1 & `lieflag'==0
		// case 3: `lieflag'==1
		// nb: `tvflag'==1 & `lieflag'==1 is impossible
	
		if `firstrep'>1 {
			di as text "Resample `m'..."
		}

		// list of fold IDs refers to current set of reps (not full set including appended)
		local thisrep = `m'-`firstrep'+1
		local fid : word `thisrep' of `fidlist'

		// create blank fitted variable(s)
		if ~`tvflag' & ~`lieflag' { // case 1
			if `ssflag' {
				cap drop `shortstack'_ss_`m'
				qui gen `vtype' `shortstack'_ss_`m'=.
			}
			forvalues j=1/`nlearners' {
				tempvar vhat`j' vres`j'
				qui gen `vtype' `vhat`j''=.
				qui gen `vtype' `vres`j''=.
			}
		}  
		else if `tvflag' & ~`lieflag' { // case 2
			if `ssflag' {
				cap drop `shortstack'_ss0_`m'
				cap drop `shortstack'_ss1_`m'
				qui gen `vtype' `shortstack'_ss0_`m'=.
				qui gen `vtype' `shortstack'_ss1_`m'=.
			}
			forvalues j=1/`nlearners' {
				tempvar vhat0`j' vres0`j'
				tempvar vhat1`j' vres1`j'
				qui gen `vtype' `vhat0`j''=.
				qui gen `vtype' `vres0`j''=.	
				qui gen `vtype' `vhat1`j''=.
				qui gen `vtype' `vres1`j''=.	
			}
		}
		else if `lieflag' { // case 3
			// out-of-sample predicted values for E[D|XZ] 
			forvalues j=1/`nlearners' {
				tempvar dhat`j'
				qui gen `vtype' `dhat`j''=.  // using learner i
			}
			// out-of-sample predicted values for E[D|X] 
			forvalues j=1/`nlearners' {
				tempvar hhat`j'
				qui gen `vtype' `hhat`j''=.  // using learner i in both steps
			}
			// predicted values for E[D|ZX] for each k & learner i
			forvalues k=1/`kfolds' {
				forvalues j=1/`nlearners' {
					tempvar dhat_`j'_`k'
					qui gen `vtype' `dhat_`j'_`k''=.
				}
			}
			if `ssflag' { // with short-stacking
				// for final results
				cap drop `shortstack'_ss_`m'
				cap drop `shortstack'_h_ss_`m'
				qui gen `vtype' `shortstack'_ss_`m'=.
				qui gen `vtype' `shortstack'_h_ss_`m'=.			
				// in-sample predicted values for E[D|ZX] for each k: short-stacked
				forvalues k=1/`kfolds' {		
					tempvar dhat_isSS_`k' 
					qui gen `vtype' `dhat_isSS_`k''=.
				}
				// out-of-sample predicted values for E[D|ZX] from short-stacking by k
				// NB: this is different from "dhatSS", which are from applying constrained regression to full sample
				tempvar dhat_oosSS 
				qui gen `vtype' `dhat_oosSS'=. 
				// out-of-sample predicted values for E[D|X] 
				forvalues j=1/`nlearners' {
					tempvar hhatSS`j'
					qui gen `vtype' `hhatSS`j''=. // using short-stacking for E[D|XZ] & then using learner i for E[D|X]
				}
				// short-stacking predicted values
				tempvar dhatSS
				qui gen `vtype' `dhatSS'=. // this will become `shortstack'
				tempvar hhatSS
				qui gen `vtype' `hhatSS'=. // this will become `shortstack'_h_ss_
			} 
		}
		else {
			di as err "internal crossfit error"
			exit 198
		}
		
		// crossfit
		di as text "Cross-fitting fold " _c
		forvalues k = 1(1)`kfolds' {
	
			di as text "`k' " _c
	
			forvalues j=1/`nlearners' {
				local vtilde : word `j' of `vtlist'
				mata: st_local("est_main", return_learner_item(`ename',"`vtilde'","est_main"))
				mata: st_local("est_options", return_learner_item(`ename',"`vtilde'","est_options"))
				mata: st_local("predopt",return_learner_item(`ename',"`vtilde'","predopt"))
				mata: st_local("vtype",return_learner_item(`ename',"`vtilde'","vtype"))
				if `lieflag' {
					// LIE locals
					local hhat `vtilde'_h
					mata: st_local("est_main_h", return_learner_item(`ename',"`vtilde'","est_main_h"))
					mata: st_local("est_options_h", return_learner_item(`ename',"`vtilde'","est_options_h"))
					mata: st_local("predopt_h",return_learner_item(`ename',"`vtilde'","predopt_h"))
					mata: st_local("vtype_h",return_learner_item(`ename',"`vtilde'","vtype_h"))			
				}
				
				if ~`tvflag' & ~`lieflag' { // case 1
				
					tempvar vhat_k
					
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `touse', `est_options'
					local cmd `e(cmd)'
					if "`cmd'"=="" {
						// macro e(cmd) may be missing, so deduce from command line
						local cmd : word 1 of `est_main'
					}
					
					// get fitted values and residuals for kth fold	
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'
		
					// get predicted values
					qui replace `vhat`j'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres`j'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
				}
		
				else if `tvflag' & ~`lieflag' {	// case 2: interactive models
		
					// outcome equation so estimate separately
		
					// for treatvar = 1
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `treatvar' == 1 & `touse', `est_options'
					local cmd `e(cmd)'
					if "`cmd'"=="" {
						// macro e(cmd) may be missing, so deduce from command line
						local cmd : word 1 of `est_main'
					}
		
					// get fitted values for kth fold	
					tempvar vhat_k
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'
					qui replace `vhat1`j'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres1`j'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
		
					// for treatvar = 0

					// first, we need to account for D always = 0 if Z=0
					// check if dvar is always 0
					qui count if `treatvar'==0 & `vname'!=0 & `fid'==`k' & `touse' 

					if (`r(N)'==0 & "`allowallzero'"!="") {
						// get fitted values for kth fold	
						tempvar vhat_k	
						`qui' gen `vhat_k'=10e-12  if `fid'==`k' & `touse' 	
					} 
					else {
						// estimate excluding kth fold
						`qui' `est_main' if `fid'!=`k' & `treatvar' == 0 & `touse', `est_options'
			
						// get fitted values for kth fold	
						tempvar vhat_k
						qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'					

					}
					qui replace `vhat0`j'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres0`j'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
					
				}
		
				else if `lieflag' { // case 3
		
					tempvar vhat_k // stores predicted values for E[D|ZX] temporarily
					tempvar vtil_k // stores predicted values for E[D^|X] temporarily
		
					// Step I: estimation of E[D|XZ]=D^
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `touse', `est_options'
					local cmd `e(cmd)'
					if "`cmd'"=="" {
						// macro e(cmd) may be missing, so deduce from command line
						local cmd : word 1 of `est_main'
					}
		
					// get fitted values (in and out of sample)
					qui predict `vtype' `vhat_k' if `touse', `predopt'
		
					// get *combined* out-of-sample predicted values
					qui replace `dhat`j'' = `vhat_k' if `fid'==`k' & `touse'
		
					// get predicted values in wide format for each i and k
					qui replace `dhat_`j'_`k'' = `vhat_k' if `touse'
		
					// Step II: estimation of E[D^|X]
		
					// replace {D}-placeholder in estimation string with variable name
					local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_`j'_`k''",1)
		
					// estimation	
					`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
					local cmd_h `e(cmd)'
		
					// get fitted values  
					qui predict `vtype' `vtil_k' if `touse', `predopt_h'
		
					// get *combined* out-of-sample predicted values
					qui replace `hhat`j'' = `vtil_k' if `fid'==`k' & `touse'
		
				}
				
				if `k'==1 & `m'==1 {
					local cmd_list   `cmd_list'   `cmd'
					local cmd_h_list `cmd_h_list' `cmd_h'
				}
			}
		}
	
		// last fold, insert new line
		di as text "...completed cross-fitting" _c
		

		******************************** SHORTSTACKING ************************************
	

		// shortstacking. we need to distinguish between 3 cases again. 
		if `ssflag' & `nlearners'>1 {
	
			if ~`tvflag' & ~`lieflag' { // case 1
				local vhats
				forvalues j=1/`nlearners' {
					local vhats `vhats' `vhat`j''
				}
				tempvar vss
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (additive model):"
				`qui' get_stack_weights `vname' `vhats', finalest(`ssfinalest') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw
				mat `ssw' = e(b)
				tempvar vtemp
				qui predict `vtype' `vtemp'
				qui replace `shortstack'_ss_`m' = `vtemp'
					
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss_`m' = `vname' - `shortstack'_ss_`m'
				}
			}
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
				local vhats1 
				local vhats0
				// treatvar == 1
				forvalues j=1/`nlearners' {
					local vhats1 `vhats1' `vhat1`j''
				}
				tempvar vtemp
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=1):"
				`qui' get_stack_weights `vname' `vhats1' if `treatvar'==1, finalest(`ssfinalest') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw1
				mat `ssw1' = e(b)
				qui predict `vtype' `vtemp'
				qui replace `shortstack'_ss1_`m'=`vtemp'
					
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss1_`m' = `vname' - `shortstack'_ss1_`m'
				}
				// treatvar == 0
				forvalues j=1/`nlearners' {
					local vhats0 `vhats0' `vhat0`j''
				}
				tempvar vtemp
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=0):"
				`qui' get_stack_weights `vname' `vhats0' if `treatvar'==0, finalest(`ssfinalest') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw0
				mat `ssw0' = e(b)
				qui predict `vtype' `vtemp'
				qui replace `shortstack'_ss0_`m'=`vtemp'
		
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss0_`m' = `vname' - `shortstack'_ss0_`m'
				}
			}
			else if `lieflag' {
				// apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
				local dhats
				forvalues j=1/`nlearners' {
					local dhats `dhats' `dhat`j''
				}
				`qui' di
				`qui' di as text "Short-stacking, finalest=`ssfinalest' (LIE, OOS E[D|XZ]):"
				`qui' get_stack_weights `vname' `dhats' if `touse'
				tempname ssw
				mat `ssw'= e(b)
				tempvar vtemp
				qui predict `vtype' `vtemp' if `touse'
				qui replace `dhatSS'=`vtemp' 
	
				// apply short-stacking to in-sample predicted values of E[D|XZ] *for each k*
				forvalues k = 1(1)`kfolds' {
					local dhats_is
					forvalues j=1/`nlearners' {
						local dhats_is `dhats_is' `dhat_`j'_`k''
					}
					tempvar vtemp
					`qui' di
					`qui' di as text "Short-stacking, finalest=`ssfinalest' (LIE, in-sample E[D|XZ] fold `k':"
					`qui' get_stack_weights `vname' `dhats_is' if `fid'!=`k' & `touse' 
					qui predict `vtype' `vtemp'
					qui replace `dhat_isSS_`k'' = `vtemp' if `fid'!=`k' & `touse'
					qui replace `dhat_oosSS' = `vtemp' if `fid'==`k' & `touse'
				}
	
				// need to cross-fit stacked in-sample predicted values against X
				forvalues k = 1(1)`kfolds' {
					forvalues j=1/`nlearners' {
						local vtilde : word `j' of `vtlist'
						mata: st_local("est_main_h", return_learner_item(`ename',"`vtilde'","est_main_h"))
						mata: st_local("est_options_h", return_learner_item(`ename',"`vtilde'","est_options_h"))
						mata: st_local("predopt_h",return_learner_item(`ename',"`vtilde'","predopt_h"))
						mata: st_local("vtype_h",return_learner_item(`ename',"`vtilde'","vtype_h"))				
	
						// replace {D}-placeholder in estimation string with variable name
						local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_isSS_`k''",1)
			
						// estimation	
						`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
						local cmd_h `e(cmd)'
					
						// get fitted values  
						tempvar vtemp
						qui predict `vtype' `vtemp' if `touse', `predopt_h'
			
						// get out-of-sample predicted values
						qui replace `hhatSS`j'' = `vtemp' if `fid'==`k' & `touse'
					}
				}
	
				// final stacking for E[D|X]
				tempname ssw_h ssw_h_temp
				forvalues k = 1(1)`kfolds' {
					local hhatSS_list
					forvalues j=1/`nlearners' {
						local hhatSS_list `hhatSS_list' `hhatSS`j''
					}
					`qui' di
					`qui' di as text "Short-stacking, finalest=`ssfinalest' (LIE, E[D|X]):"
					`qui' get_stack_weights `dhat_oosSS' `hhatSS_list'
					if (`k'==1) {
						mat `ssw_h' = e(b)
					}
					else {
						mat `ssw_h_temp' = e(b)
						mat  `ssw_h_temp' = (`ssw_h' \ `ssw_h_temp' )
					}
					tempvar vtemp
					qui predict `vtype' `vtemp'
					qui replace `hhatSS'=`vtemp'
				}
				qui replace `shortstack'_ss_`m'=`dhatSS'
				label var `shortstack'_ss_`m' "short-stacking cross-fitted E[D|Z,X]"
				qui replace `shortstack'_h_ss_`m'=`hhatSS'
				label var `shortstack'_h_ss_`m' "short-stacking cross-fitted E[D|X]"
			}
			else {
				di as err "internal crossfit error"
				exit 198
			}
		}
		else if `ssflag' {
			// single learner case, so shortstack vars are just copies of learner vars
			
			if ~`tvflag' & ~`lieflag' { // case 1
				qui replace `shortstack'_ss_`m' = `vhat1'
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss_`m' = `vname' - `shortstack'_ss_`m'
				}
			}
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
				qui replace `shortstack'_ss1_`m'=`vhat11'
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss1_`m' = `vname' - `shortstack'_ss1_`m' if `treatvar'==1
				}
				qui replace `shortstack'_ss0_`m'=`vhat01'
				if `residflag' {
					// vtilde is the residual
					qui replace `shortstack'_ss0_`m' = `vname' - `shortstack'_ss0_`m' if `treatvar'==0
				}
			}
			else if `lieflag' {
				qui replace `shortstack'_ss_`m'=`dhat1'
				label var `shortstack'_ss_`m' "short-stacking cross-fitted E[D|Z,X]"
				qui replace `shortstack'_h_ss_`m'=`hhat1'
				label var `shortstack'_h_ss_`m' "short-stacking cross-fitted E[D|X]"
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
		
		forvalues j=1/`nlearners' {
			
			local vtilde	: word `j' of `vtlist'
			local cmd		: word `j' of `cmd_list'
			local cmd_h		: word `j' of `cmd_h_list'

			// vtilde, mspe, etc.
			if ~`tvflag' & ~`lieflag' {
		
				cap drop `vtilde'_`m'
				if "`resid'"=="" {
					// vtilde is predicted values
					qui gen `vtype' `vtilde'_`m' = `vhat`j''
					qui label var `vtilde'_`m' "Predicted values cond. exp. of `vname' using `cmd'"
				}
				else {
					// vtilde is residuals
					qui gen `vtype' `vtilde'_`m' = `vres`j''
					qui label var `vtilde'_`m' "Residuals cond. exp. of `vname' using `cmd'"
				}
		
				// calculate and return mspe and sample size
				tempvar vres_sq
				qui gen double `vres_sq' = `vres`j''^2 if `touse'

				// additive-type model
				qui sum `vres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `vres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
			
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')
			
				mata: add_result_item(`ename',"`vtilde'","N",         "`m'", `N')
				mata: add_result_item(`ename',"`vtilde'","N_folds",   "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`vtilde'","MSE",       "`m'", `mse')
				mata: add_result_item(`ename',"`vtilde'","MSE_folds", "`m'", st_matrix("`mse_folds'"))
				
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
			else if `tvflag' & ~`lieflag' {
			
				cap drop `vtilde'0_`m'
				cap drop `vtilde'1_`m'
				if "`resid'"=="" {
					// vtilde is predicted values
					qui gen `vtype' `vtilde'0_`m' = `vhat0`j''
					qui gen `vtype' `vtilde'1_`m' = `vhat1`j''
					qui label var `vtilde'0_`m' "Predicted values cond. exp. of `vname' given `treatvar'==0 using `cmd'"
					qui label var `vtilde'1_`m' "Predicted values cond. exp. of `vname' given `treatvar'==1 using `cmd'"
				}
				else {
					// vtilde is residuals
					qui gen `vtype' `vtilde'0_`m' = `vres0`j''
					qui gen `vtype' `vtilde'1_`m' = `vres1`j''
					qui label var `vtilde'0_`m' "Residuals cond. exp. of `vname' given `treatvar'==0 using `cmd'"
					qui label var `vtilde'1_`m' "Residuals cond. exp. of `vname' given `treatvar'==1 using `cmd'"
				}
		
				// calculate and return mspe and sample size
				tempvar vres0_sq vres1_sq
				// vtilde has fitted values
				qui gen double `vres0_sq' = (`vname' - `vhat0`j'')^2 if `treatvar' == 0 & `touse'
				qui gen double `vres1_sq' = (`vname' - `vhat1`j'')^2 if `treatvar' == 1 & `touse'
		
				// interactive-type model, return mse separately for treatvar =0 and =1
				qui sum `vres0_sq' if `treatvar' == 0 & `touse', meanonly
				local mse0			= r(mean)
				local N0			= r(N)
				qui sum `vres1_sq' if `treatvar' == 1 & `touse', meanonly
				local mse1			= r(mean)
				local N1			= r(N)
				local N				= `N0'+`N1'
				tempname mse0_folds N0_folds mse1_folds N1_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres0_sq' if `treatvar' == 0 & `touse' & `fid'==`k', meanonly
					mat `mse0_folds' = (nullmat(`mse0_folds'), r(mean))
					qui sum `vres1_sq' if `treatvar' == 1 & `touse' & `fid'==`k', meanonly
					mat `mse1_folds' = (nullmat(`mse1_folds'), r(mean))
					qui count if `treatvar' == 0 & `touse' & `fid'==`k' & `vres0_sq'<.
					mat `N0_folds' = (nullmat(`N0_folds'), r(N))
					qui count if `treatvar' == 1 & `touse' & `fid'==`k' & `vres1_sq'<.
					mat `N1_folds' = (nullmat(`N1_folds'), r(N))
				}
	
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list') \ `N0')
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ `mse0_folds')
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')\ `N0_folds')
				
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list') \ `N1')
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ `mse1_folds')
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')\ `N1_folds')

				forvalues t=0/1 {
					mata: add_result_item(`ename',"`vtilde'","N`t'",         "`m'", `N`t'')
					mata: add_result_item(`ename',"`vtilde'","N`t'_folds",   "`m'", st_matrix("`N`t'_folds'"))
					mata: add_result_item(`ename',"`vtilde'","MSE`t'",       "`m'", `mse`t'')
					mata: add_result_item(`ename',"`vtilde'","MSE`t'_folds", "`m'", st_matrix("`mse`t'_folds'"))
				}
				
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
			else if `lieflag' {
	
				cap drop `vtilde'_`m'
				cap drop `vtilde'_h_`m'
				
				qui gen `vtype' `vtilde'_`m' = `dhat`j''
				qui label var `vtilde'_`m' "Predicted values E[`vname'|X,Z]"
				qui gen `vtype' `vtilde'_h_`m' = `hhat`j''
				qui label var `vtilde'_h_`m' "Predicted values E[`vtilde'_`m'|X]"

				// calculate and return mspe and sample size
				tempvar hres dres hres_sq dres_sq
				// vtilde has fitted values
				qui gen double `dres_sq' = (`vname' - `dhat`j'')^2 if `touse'
				qui gen double `hres_sq' = (`vname' - `hhat`j'')^2 if `touse'
	
				qui sum `dres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `dres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `dres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')
				
				mata: add_result_item(`ename',"`vtilde'","N",         "`m'", `N')
				mata: add_result_item(`ename',"`vtilde'","N_folds",   "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`vtilde'","MSE",	      "`m'", `mse')
				mata: add_result_item(`ename',"`vtilde'","MSE_folds", "`m'", st_matrix("`mse_folds'"))
	
				qui sum `hres_sq' if `touse', meanonly
				local mse_h			= r(mean)
				local N_h			= r(N)	
				tempname mse_h_folds N_h_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `hres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_h_folds' = (nullmat(`mse_h_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `hres_sq'<.
					mat `N_h_folds' = (nullmat(`N_h_folds'), r(N))
				}
				mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
				mat `N_h_list'			= (nullmat(`N_h_list') \ `N_h')
				mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ `mse_h_folds')
				mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')\ `N_h_folds')
				
				mata: add_result_item(`ename',"`vtilde'","N_h",         "`m'", `N_h')
				mata: add_result_item(`ename',"`vtilde'","N_h_folds",   "`m'", st_matrix("`N_h_folds'"))
				mata: add_result_item(`ename',"`vtilde'","MSE_h",	    "`m'", `mse_h')
				mata: add_result_item(`ename',"`vtilde'","MSE_h_folds", "`m'", st_matrix("`mse_h_folds'"))
				
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
		
		// add shortstack results
		if `ssflag' {
			if ~`tvflag' & ~`lieflag' {	// case 1
		
				// calculate and return mspe and sample size
				tempvar vres_sq
				// shortstack macro has the residuals
				qui gen double `vres_sq' = (`shortstack'_ss_`m')^2 if `touse'
			
				// additive-type model
				qui sum `vres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `vres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
			
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')
				mata: add_result_item(`ename',"`shortstack'_ss","N",            "`m'", `N')
				mata: add_result_item(`ename',"`shortstack'_ss","N_folds",      "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","MSE",          "`m'", `mse')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_folds",    "`m'", st_matrix("`mse_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights",   "`m'", st_matrix("`ssw'"))
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
				
			}
			else if `tvflag' & ~`lieflag' {	// case 2
			
				// calculate and return mspe and sample size
				tempvar vres0_sq vres1_sq
				// shortstack macros have fitted values
				qui gen double `vres0_sq' = (`shortstack'_ss0_`m')^2 if `treatvar' == 0 & `touse'
				qui gen double `vres1_sq' = (`shortstack'_ss1_`m')^2 if `treatvar' == 1 & `touse'
		
				// interactive-type model, return mse separately for treatvar =0 and =1
				qui sum `vres0_sq' if `treatvar' == 0 & `touse', meanonly
				local mse0			= r(mean)
				local N0			= r(N)
				qui sum `vres1_sq' if `treatvar' == 1 & `touse', meanonly
				local mse1			= r(mean)
				local N1			= r(N)
				local N				= `N0'+`N1'
				tempname mse0_folds N0_folds mse1_folds N1_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `vres0_sq' if `treatvar' == 0 & `touse' & `fid'==`k', meanonly
					mat `mse0_folds' = (nullmat(`mse0_folds'), r(mean))
					qui sum `vres1_sq' if `treatvar' == 1 & `touse' & `fid'==`k', meanonly
					mat `mse1_folds' = (nullmat(`mse1_folds'), r(mean))
					qui count if `treatvar' == 0 & `touse' & `fid'==`k' & `vres1_sq'<.
					mat `N0_folds' = (nullmat(`N0_folds'), r(N))
					qui count if `treatvar' == 1 & `touse' & `fid'==`k' & `vres0_sq'<.
					mat `N1_folds' = (nullmat(`N1_folds'), r(N))
				}
	
				mat `mse0_list'			= (nullmat(`mse0_list') \ `mse0')
				mat `N0_list'			= (nullmat(`N0_list') \ `N0')
				mat `mse0_folds_list'	= (nullmat(`mse0_folds_list') \ `mse0_folds')
				mat `N0_folds_list'		= (nullmat(`N0_folds_list')\ `N0_folds')
				
				mat `mse1_list'			= (nullmat(`mse1_list') \ `mse1')
				mat `N1_list'			= (nullmat(`N1_list') \ `N1')
				mat `mse1_folds_list'	= (nullmat(`mse1_folds_list') \ `mse1_folds')
				mat `N1_folds_list'		= (nullmat(`N1_folds_list')\ `N1_folds')
				
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`shortstack'_ss","N`t'",          "`m'", `N`t'')
					mata: add_result_item(`ename',"`shortstack'_ss","N`t'_folds",    "`m'", st_matrix("`N`t'_folds'"))
					mata: add_result_item(`ename',"`shortstack'_ss","MSE`t'",        "`m'", `mse`t'')
					mata: add_result_item(`ename',"`shortstack'_ss","MSE`t'_folds",  "`m'", st_matrix("`mse`t'_folds'"))
					mata: add_result_item(`ename',"`shortstack'_ss","ss_weights`t'", "`m'", st_matrix("`ssw`t''"))
				}
				// final estimator used to stack is a learner item
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
			}
			else if `lieflag' {	// case 3
	
				// calculate and return mspe and sample size
				tempvar hres dres hres_sq dres_sq
				// vtilde has fitted values
				qui gen double `dres_sq' = (`vname' - `shortstack'_ss_`m')^2 if `touse'
				qui gen double `hres_sq' = (`vname' - `shortstack'_h_ss_`m')^2 if `touse'
	
				qui sum `dres_sq' if `touse', meanonly
				local mse			= r(mean)
				local N				= r(N)
				tempname mse_folds N_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `dres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `dres_sq'<.
					mat `N_folds' = (nullmat(`N_folds'), r(N))
				}
				mat `mse_list'			= (nullmat(`mse_list') \ `mse')
				mat `N_list'			= (nullmat(`N_list') \ `N')
				mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
				mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')

				mata: add_result_item(`ename',"`shortstack'_ss","N",         "`m'", `N')
				mata: add_result_item(`ename',"`shortstack'_ss","N_folds",   "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","MSE",       "`m'", `mse')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_folds", "`m'", st_matrix("`mse_folds'"))
				// final estimator used to stack is a learner item; same in "_h" estimation
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
	
				qui sum `hres_sq' if `touse', meanonly
				local mse_h			= r(mean)
				local N_h			= r(N)	
				tempname mse_h_folds N_h_folds
				forvalues k = 1(1)`kfolds' {
					qui sum `hres_sq' if `touse' & `fid'==`k', meanonly
					mat `mse_h_folds' = (nullmat(`mse_h_folds'), r(mean))
					qui count if `touse' & `fid'==`k' & `hres_sq'<.
					mat `N_h_folds' = (nullmat(`N_h_folds'), r(N))
				}
				mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
				mat `N_h_list'			= (nullmat(`N_h_list') \ `N_h')
				mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ `mse_h_folds')
				mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')\ `N_h_folds')

				mata: add_result_item(`ename',"`shortstack'_ss","N_h",          "`m'", `N_h')
				mata: add_result_item(`ename',"`shortstack'_ss","N_h_folds",    "`m'", st_matrix("`N_h_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_h",        "`m'", `mse_h')
				mata: add_result_item(`ename',"`shortstack'_ss","MSE_h_folds",  "`m'", st_matrix("`mse_h_folds'"))
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights",   "`m'", st_matrix("`ssw'"))
				mata: add_result_item(`ename',"`shortstack'_ss","ss_weights_h", "`m'", st_matrix("`ssw_h'"))
			}
		}
	
	}		// end of resampling loop
	
		
	******************************** RETURN RESULTS ************************************

	if ~`tvflag' {
		foreach mname in mse_list N_list mse_folds_list N_folds_list {
			mat rownames ``mname''	= `rnames'
			return mat `mname'		= ``mname''
		}
		return mat mse_folds		= `mse_folds'
		return mat N_folds			= `N_folds'
		return scalar mse			= `mse'
	}
	else {
		foreach mname in mse0_list N0_list mse0_folds_list N0_folds_list mse1_list N1_list mse1_folds_list N1_folds_list {
			mat rownames ``mname''	= `rnames'
			return mat `mname'		= ``mname''
		}
		return scalar mse0			= `mse0'
		return scalar N0			= `N0'
		return scalar mse1			= `mse1'
		return scalar N1			= `N1'
		return mat mse0_folds		= `mse0_folds'
		return mat mse1_folds		= `mse1_folds'
		return mat N0_folds			= `N0_folds'
		return mat N1_folds			= `N1_folds'
	}
	
	if `lieflag' {
		foreach mname in mse_h_list N_h_list mse_h_folds_list N_h_folds_list {
			mat rownames ``mname''	= `rnames'
			return mat `mname'		= ``mname''
		}
		return scalar N_h			= `N_h'
		return mat mse_h_folds		= `mse_h_folds'
		return mat N_h_folds		= `N_h_folds'
		return scalar mse_h			= `mse_h'
	}

	return scalar N			= `N'
	return local cmd_list	`cmd_list'
	return local cmd_h_list	`cmd_h_list'
	
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
							vtype_h(string)			///
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
		parse_estring, vtlist(`vtlist') ename(`ename') estring(`estringh') h vtype_h(`vtype_h') predopt(`predopt') `noisily'
		mata: `ename'.lieflag = 1
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
							vtype_h(string)			///
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
			mata: add_learner_item(`ename',"`vtilde'","vtype_h","`vtype_h'")
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


program get_stack_weights, eclass sortpreserve

	version 16
	syntax varlist(numeric min=2)					///
		[if] [in] [aw fw iw] [,						/// 
									NOIsily			///
									finalest(name)	///
									stype(name)		///
									*				///
									]

	* defaults
	if "`noisily'"==""		local qui qui
	if "`finalest'"==""		local finalest nnls1
	if "`stype'"==""		local stype reg

	marksample touse
	qui count if `touse'
	local N = r(N)
	
	local yvar : word 1 of `varlist'
	local xvars : list varlist - yvar
	
    tempvar wt
    if "`weight'"~="" {
        qui gen double `wt' `exp'
    }
    else {
    	qui gen byte `wt' = 1
    }

	`qui' python: py_get_stack_weights("`yvar'","`xvars'","`touse'","`wt'","`finalest'","`stype'")
	tempname bhat
	// python returns column vectors
	mat `bhat' = r(b)'
	matrix colnames `bhat' = `xvars'
	matrix rownames `bhat' = `yvar'
	ereturn clear
	ereturn post `bhat', depname(`yvar') obs(`N') esample(`touse')
	ereturn scalar N = `N'
	ereturn display

end

version 16.0
python:

import sfi
from sfi import Data,Matrix,Scalar,SFIToolkit
from sklearn.linear_model import LinearRegression,RidgeCV
from sklearn.base import TransformerMixin,BaseEstimator
from sklearn.utils import check_X_y,check_array
from sklearn.utils.validation import _check_sample_weight
import numpy as np
from scipy.optimize import minimize 
from scipy.optimize import nnls 

def py_get_stack_weights(yvar,xvars,touse,wvar,finalest,stype):

    X = Data.get(xvars,selectvar=touse)
    y = Data.get(yvar,selectvar=touse)
    w = Data.get(wvar,selectvar=touse)

    if finalest == "nnls0" and stype == "class": 
        fin_est = LinearRegressionClassifier(fit_intercept=False,positive=True)
    elif finalest == "nnls_sk" and stype == "class": 
        fin_est = LinearRegressionClassifier(fit_intercept=False,positive=True)
    elif finalest == "nnls1" and stype == "class": 
        fin_est = ConstrLSClassifier()
    elif finalest == "ridge" and stype == "class": 
        fin_est = LogisticRegression()
    elif finalest == "nnls0" and stype == "reg": 
        fin_est = LinearRegression(fit_intercept=False,positive=True)
    elif finalest == "nnls_sk" and stype == "reg": 
        fin_est = LinearRegression(fit_intercept=False,positive=True)
    elif finalest == "nnls1" and stype == "reg": 
        fin_est = ConstrLS()
    elif finalest == "ridge" and stype == "reg": 
        fin_est = RidgeCV()
    elif finalest == "singlebest" and stype == "reg": 
        fin_est = SingleBest()
    elif finalest == "ols" and stype == "class": 
        fin_est = LinearRegressionClassifier()    
    elif finalest == "ols" and stype == "reg": 
        fin_est = LinearRegression()    
    elif finalest == "ls1" and stype == "reg":
        fin_est = ConstrLS(unit_interval=False)    
    elif finalest == "ls1" and stype == "class":
        fin_est = ConstrLSClassifier(unit_interval=False)    
    else:
        sfi.SFIToolkit.stata('di as err "specified final estimator not supported"')
        #"
        sfi.SFIToolkit.error(198)
    fin_est.fit(X, y, w)
    b = fin_est.coef_
    Matrix.store("r(b)", b)

class ConstrLS(BaseEstimator):
    _estimator_type="regressor"
    def fit(self, X, y, w):

        X,y = check_X_y(X,y, accept_sparse=True)
        xdim = X.shape[1]

        #Use nnls to get initial guess
        #coef0, rnorm = nnls(X,y)
        #Use LinearRegression to get initial guess
        initial_est = LinearRegression(positive=True,fit_intercept=False)
        initial_est.fit(X, y, w)
        coef0 = initial_est.coef_

        #Define minimisation function
        def fn(coef, X, y):
            return np.linalg.norm(X.dot(coef) - y)
        
        #Constraints and bounds
        cons = {'type': 'eq', 'fun': lambda coef: np.sum(coef)-1}
        if self.unit_interval==True:
            bounds = [[0.0,1.0] for i in range(xdim)] 
        else:
            bounds = None

        w = _check_sample_weight(w, X)
        #If weights vector=1, no weighting needed
        if np.all(w==1):
            #Do minimisation
            fit = minimize(fn,coef0,args=(X, y),method='SLSQP',bounds=bounds,constraints=cons)
        else:
            #Use additional precision
            Xw = np.multiply(X,w,dtype=np.longdouble)
            yw = np.multiply(y,w,dtype=np.longdouble)
            #Do minimisation
            fit = minimize(fn,coef0,args=(Xw, yw),method='SLSQP',bounds=bounds,constraints=cons)
        
        self.coef_ = fit.x
        if abs(sum(self.coef_)-1)>1e-12:
            #Renormalize if coefs still don't sum to 1
            self.coef_ = self.coef_ * 1/sum(self.coef_)
        self.is_fitted_ = True
        return self
        
    def predict(self, X):
        X = check_array(X, accept_sparse=True)
        check_is_fitted(self, 'is_fitted_')
        return np.matmul(X,self.coef_)

    def __init__(self, unit_interval=True):
        self.unit_interval = unit_interval

class SingleBest(BaseEstimator):
    _estimator_type="regressor"
    def fit(self, X, y, w):
        X, y = check_X_y(X, y, accept_sparse=True)
        self.is_fitted_ = True
        ncols = X.shape[1]
        lowest_mse = np.Inf
        for i in range(ncols):
            this_mse=np.mean((y-X[:, i]) ** 2)
            if this_mse < lowest_mse:
                lowest_mse = this_mse
                best = i
        self.best = best
        coef = np.zeros(ncols)
        coef[best] = 1
        self.coef_ = coef
        self.cvalid=X
        return self
    def predict(self, X):
        X = check_array(X, accept_sparse=True)
        check_is_fitted(self, 'is_fitted_')
        return X[:,self.best]

class ConstrLSClassifier(ConstrLS):
    _estimator_type="classifier"
    def predict_proba(self, X):
        return self.predict(X)


end

