*! crossfit v0.6
*! last edited: 9 july 2023
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
		// default is to use pystacked-specific code unless user says "no" or doesn't apply
		if `pystackedmulti' {
			// user says "yes if possible" so set to whatever eqn says
			mata: st_local("pystackedmulti", strofreal(`eqn_info'.pystackedmulti))
		}
	}
	
	if `pystackedmulti' {
		// enter single learner = pystacked and number of pystacked base learners > 1
		_crossfit_pystacked	`if' `in',	///
			ename(`eqn_info')			///
			`noisily'					///
			`cfoptions'			
	}
	else {
		// multiple ddml learners or single learner = pystacked with a single base learner
		_crossfit_other	`if' `in',		///
			ename(`eqn_info')			///
			`noisily'					///
			`cfoptions'			
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
							NOIsily					///
							allowallzero			/// in the LATE model: allow D
							stdfinalest(name)		/// final estimator for standard stacking
							ssfinalest(name)		/// final estimator for short-stacking
							psfinalest(name)		/// final estimator for pooled-stacking
							finalest(name)			/// final estimator for all
							NOSTDstack				/// no standard stacking - use psytacked+voting to get learners only
							*						/// ignored options
							]

	// used throughout	
	local cmd pystacked

	// unless specified, finalest sets the final estimator for all stacking methods
	// if empty, finalest will be the pystacked/ddml default
	if "`ssfinalest'"==""	local ssfinalest `finalest'
	if "`psfinalest'"==""	local psfinalest `finalest'
	if "`stdfinalest'"==""	local stdfinalest `finalest'
	
	mata: st_local("vname", `ename'.vname)
	mata: st_local("vtlist", invtokens(`ename'.vtlist))
	mata: st_local("shortstack", `ename'.shortstack)
	mata: st_local("poolstack", `ename'.poolstack)
	mata: st_local("lieflag", strofreal(`ename'.lieflag))
	
	** error check
	if `lieflag' {
		di as err "internal crossfit error: _crossfit_pystacked does not support model=fiv"
		exit 198
	}

	** indicator for interactive model
	local tvflag	= "`treatvar'"~=""
	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
	** indicator for pooled-stacking
	local psflag	= "`poolstack'"~=""
	** indicator for standard stacking
	local stdflag	= "`nostdstack'"==""
	if `psflag' & ~`stdflag' {
		di as res "pooled stacking requires standard stacking; poolstack option ignored"
		local psflag = 0
	}
	if ~`stdflag' & ~`ssflag' {
		di as err "error - pystacked integration requires using either standard and/or short-stacking"
		exit 198
	}
	// pystacked option for voting instead of stacking
	// votetype is ignored if type=reg; relevant only for type=class
	if ~`stdflag' {
		local nostdstackopt voting votetype(soft)
	}
	
	if "`noisily'"=="" {
		local qui quietly
	}
	
	*** debugging message
	`qui' di as text "entering _crossfit_pystacked..."

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
	
	// stored results; initialize here
	tempname mse_list N_list mse_folds_list N_folds_list
	tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
	tempname mse0_list N0_list mse0_folds_list N0_folds_list
	tempname mse1_list N1_list mse1_folds_list N1_folds_list
	
	// will save base learner CV predictions in mata
	tempname y_stacking_cv y_stacking_cv0 y_stacking_cv1
	
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
		// create list of names of variables created by learners for rep m
		// need to reset for each m
		if ~`tvflag' {
			// PLM
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
		
		******************************** CROSSFITTING ************************************
		
		// cross-fitting fundamentally depends on three cases: plm, interactive, and LIE
		// but pystacked-specific code doesn't support LIE, hence only two cases in this code
		// case 1: `tvflag'==0
		// case 2: `tvflag'==1
		
		if `lastrep'>1 {
			di as text "Resample `m'..."
		}

		// list of fold IDs refers to current set of reps (not full set including appended)
		local thisrep = `m'-`firstrep'+1
		local fid : word `thisrep' of `fidlist'

		// create blank fitted variable(s), other initializations for resample m
		if ~`tvflag'{ // case 1
			tempvar vhat vres
			qui gen `vtype' `vhat'=.
			// base learner predicted values
			// not tempvars; name based on `vtilde' macro; used for poolstacking
			forvalues j=1/`nlearners' {
				cap drop `vtilde'_L`j'_`m'
				qui gen `vtype' `vtilde'_L`j'_`m' = .
			}
		}  
		else { // case 2
			tempvar vhat0
			tempvar vhat1
			qui gen `vtype' `vhat0'=.
			qui gen `vtype' `vhat1'=.
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
	
		// will save pystacked weights and MSEs
		tempname pysw pysw0 pysw1
		tempname pysw_temp pysw1_temp pysw0_temp
		tempname pysm pysm0 pysm1
		tempname pysm_temp pysm1_temp pysm0_temp
		tempvar tomata fidtouse fidtouse1 fidtouse0
		qui gen byte `tomata'=.
		qui gen int `fidtouse'=.
		qui gen int `fidtouse0'=.
		qui gen int `fidtouse1'=.
		
		// crossfit
		di as text "Cross-fitting fold " _c
		forvalues k = 1(1)`kfolds' {
	
			di as text "`k' " _c
			
			// needed for each fold loop
			// predicted values for fold k based on estimation for not-k
			tempvar vhat_k vhat1_k vhat0_k
			
			if ~`tvflag' { // case 1
				
				if (`k'==1) & `stdflag' {
					// initialize mata object to hold dep var and in-sample crossfit CV predictions
					mata: `y_stacking_cv' = J(0,`nlearners'+3,0)
				}
				
				// estimate excluding kth fold
				`qui' `est_main' if `fid'!=`k' & `touse', `est_options' `nostdstackopt' finalest(`stdfinalest')
				`qui' di as text "N=" as res e(N)
				local base_est `e(base_est)'
				local stack_final_est `e(finalest)'
				// check
				assert "`e(cmd)'"=="pystacked"
				assert e(mcount)>1 & e(mcount)<.
				
				// pystacked complains if data change after estimation and before predict
				// so make sure all predict commands are done first
				
				// always save crossfit base learner predicted values
				// OOS (crossfit) base learner predicted values
				tempvar stacking_p
				qui predict double `stacking_p' if `fid'==`k' & `touse', basexb
				
				if `stdflag' {
					// get fitted values for kth fold	
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse'
					// in-sample CV base learner predicted values
					tempvar stacking_p_cv
					qui predict double `stacking_p_cv', basexb cv
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
				}
				
				// rename OOS (crossfit) base learner predicted values
				forvalues j=1/`nlearners' {
					qui replace `vtilde'_L`j'_`m' = `stacking_p'`j' if `fid'==`k' & `touse'
				}
				if `stdflag' {
					// get stacked out-of-sample predicted values
					qui replace `vhat' = `vhat_k' if `fid'==`k' & `touse'
				}
				if `stdflag' {				
					// in-sample CV base learner predicted values
					// accumulated in mata along with corresponding values of dep var
					qui replace `tomata' = `fid'!=`k' & `touse'
					qui replace `fidtouse' = (`k'*`touse')
					fvexpand `stacking_p_cv'*
					mata: `y_stacking_cv' = `y_stacking_cv' \ st_data(., "`fid' `fidtouse' `vname' `r(varlist)'", "`tomata'")
				}
			}
	
			else {	// case 2: interactive models
	
				// pystacked learners
				if (`k'==1) & `stdflag' {
					// initialize mata object to hold dep var and in-sample crossfit CV predictions
					mata: `y_stacking_cv0' = J(0,`nlearners'+3,0)
					mata: `y_stacking_cv1' = J(0,`nlearners'+3,0)
				}
	
				// outcome equation so estimate separately for treated and untreated
	
				// for treatvar = 1
				// estimate excluding kth fold
				`qui' `est_main' if `fid'!=`k' & `treatvar' == 1 & `touse', `est_options' finalest(`stdfinalest')
				`qui' di as text "N=" as res e(N)
				local base_est `e(base_est)'
				local stack_final_est `e(finalest)'
				// check
				assert "`e(cmd)'"=="pystacked"
				assert e(mcount)>1 & e(mcount)<.
				
				// pystacked complains if data change after estimation and before predict
				// so make sure all predict commands are done first
				
				// always save crossfit base learner predicted values
				// OOS (crossfit) base learner predicted values
				tempvar stacking_p
				qui predict double `stacking_p' if `fid'==`k' & `touse', basexb
				
				if `stdflag' {
					// get fitted values for kth fold	
					qui predict `vtype' `vhat1_k' if `fid'==`k' & `touse'
					// in-sample CV base learner predicted values
					tempvar stacking_p_cv
					qui predict double `stacking_p_cv', basexb cv
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
				}
				
				// rename OOS (crossfit) base learner predicted values
				forvalues j=1/`nlearners' {
					qui replace `vtilde'1_L`j'_`m' = `stacking_p'`j' if `fid'==`k' & `touse'
				}
				if `stdflag' {
					// get stacked out-of-sample predicted values
					qui replace `vhat1' = `vhat1_k' if `fid'==`k' & `touse'
				}
				if `stdflag' {
					// in-sample CV base learner predicted values
					// accumulated in mata along with corresponding values of dep var
					qui replace `tomata' = `fid'!=`k' & `touse'
					qui replace `fidtouse1' = (`k'*`touse')
					fvexpand `stacking_p_cv'*
					mata: `y_stacking_cv1' = `y_stacking_cv1' \ st_data(., "`fid' `fidtouse1' `vname' `r(varlist)'", "`tomata'")
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
					`qui' `est_main' if `fid'!=`k' & `treatvar' == 0 & `touse', `est_options' finalest(`stdfinalest')
					`qui' di as text "N=" as res e(N)
					local base_est `e(base_est)'
					local stack_final_est `e(finalest)'
		
					// pystacked complains if data change after estimation and before predict
					// so make sure all predict commands are done first
				
					// always save crossfit base learner predicted values
					// OOS (crossfit) base learner predicted values
					tempvar stacking_p
					qui predict double `stacking_p' if `fid'==`k' & `touse', basexb

					if `stdflag' {
						// get fitted values for kth fold	
						qui predict `vtype' `vhat0_k' if `fid'==`k' & `touse'
						// in-sample CV base learner predicted values
						tempvar stacking_p_cv
						qui predict double `stacking_p_cv', basexb cv
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
					}

					// rename OOS (crossfit) base learner predicted values
					forvalues j=1/`nlearners' {
						qui replace `vtilde'0_L`j'_`m' = `stacking_p'`j' if `fid'==`k' & `touse'
					}
				}
				
				if `stdflag' {
					// get stacked out-of-sample predicted values
					qui replace `vhat0' = `vhat0_k' if `fid'==`k' & `touse'
				}
				// in-sample CV base learner predicted values
				// accumulated in mata along with corresponding values of dep var
				if `stdflag' {
					qui replace `tomata' = `fid'!=`k' & `touse'
					qui replace `fidtouse0' = (`k'*`touse')
					fvexpand `stacking_p_cv'*
					mata: `y_stacking_cv0' = `y_stacking_cv0' \ st_data(., "`fid' `fidtouse0' `vname' `r(varlist)'", "`tomata'")
				}
			}
		}
	
		// last fold
		di as text "...completed cross-fitting" _c
		// if noisily, print new line
		`qui' di

		******************************** SHORTSTACKING ************************************

		// shortstacking. we need to distinguish between 3 cases again. 
		if `ssflag' {
	
			if ~`tvflag' { // case 1
				// stack dep var against all OOS (cross-fit) learner predicted values
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (additive model):"
				`qui' _ddml_nnls `vname' `vtilde_list' if `touse', finalest(`ssfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw
				mat `ssw' = e(b)
				cap drop `shortstack'_ss_`m'
				mat score double `shortstack'_ss_`m' = `ssw' if `touse'
			}
			else {	// case 2: interactive models
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=1):"
				`qui' _ddml_nnls `vname' `vtilde1_list' if `touse' & `treatvar'==1, finalest(`ssfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw1
				mat `ssw1' = e(b)
				cap drop `shortstack'_ss1_`m'
				mat score double `shortstack'_ss1_`m' = `ssw1' if `touse'
				// treatvar == 0
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=0):"
				`qui' _ddml_nnls `vname' `vtilde0_list' if `touse' & `treatvar'==0, finalest(`ssfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw0
				mat `ssw0' = e(b)
				cap drop `shortstack'_ss0_`m'
				mat score double `shortstack'_ss0_`m' = `ssw0' if `touse'
			}
		}
	
		if `ssflag' {
			di as text "...completed short-stacking" _c
		}

		******************************** POOLSTACKING *************************************
		
		if `psflag' {
			// mata object y_stacking_cv has y and predicted yhats of all learners in all crossfits
			if ~`tvflag' { // case 1
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
				tempname psw
				mata: `psw' = st_matrix("e(b)")
				frame change `cframe'
				frame drop `tframe'
				mata: st_matrix("`psw'",`psw')
				mat colnames `psw' =  `vtilde_list'
				cap drop `poolstack'_ps_`m'
				mat score double `poolstack'_ps_`m' = `psw' if `touse'
			}
			else {	// case 2: interactive models
			
				// treatvar=1
				tempname tframe
				qui frame pwf
				local cframe `r(currentframe)'
				frame create `tframe'
				frame change `tframe'
				getmata (`fid' `fidtouse' `vname' `vtilde1_list')=`y_stacking_cv1', force replace
				`qui' di
				`qui' di as text "Pooled-stacking, finalest=`psfinalest' (additive model):"
				`qui' _ddml_nnls `vname' `vtilde1_list', finalest(`psfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname psw1
				mata: `psw1' = st_matrix("e(b)")
				frame change `cframe'
				frame drop `tframe'
				mata: st_matrix("`psw1'",`psw1')
				mat colnames `psw1' = `vtilde1_list'
				cap drop `poolstack'_ps1_`m'
				mat score double `poolstack'_ps1_`m' = `psw1' if `touse'

				// treatvar=0
				tempname tframe
				qui frame pwf
				local cframe `r(currentframe)'
				frame create `tframe'
				frame change `tframe'
				getmata (`fid' `fidtouse' `vname' `vtilde0_list')=`y_stacking_cv0', force replace
				`qui' di
				`qui' di as text "Pooled-stacking, finalest=`psfinalest' (additive model):"
				`qui' _ddml_nnls `vname' `vtilde0_list', finalest(`psfinalest') stype(`stype') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname psw0
				mata: `psw0' = st_matrix("e(b)")
				frame change `cframe'
				frame drop `tframe'
				mata: st_matrix("`psw0'",`psw0')
				mat colnames `psw0' =  `vtilde0_list'
				cap drop `poolstack'_ps0_`m'
				mat score double `poolstack'_ps0_`m' = `psw0' if `touse'
			}
		}
		if `psflag' {
			di as text "...completed pooled-stacking" _c
		}
		
		************************************************************************************

		***************************** ESTIMATION COMPLETE *********************************
		// estimation done, insert newline
		di
		******************************** STORE RESULTS ************************************

		// vtilde, mspe, etc.
		if ~`tvflag' {
	
			// alway label learner predicted values
			forvalues j=1/`nlearners' {
				qui label var `vtilde'_L`j'_`m' "Pred. values E[`vname'|X] using base learner `j', rep `m'"
			}
			
			// save results relating to stacked learner if it exists
			if `stdflag' {
			
				// vtilde has fitted values
				cap drop `vtilde'_`m'
				qui gen `vtype' `vtilde'_`m' = `vhat'
				qui label var `vtilde'_`m' "Pred. values E[`vname'|X] using `cmd', rep `m'"
				tempvar vres_sq
				qui gen double `vres_sq' = (`vname' - `vhat')^2 if `touse'
	
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
				mata: add_learner_item(`ename',"`vtilde'","stack_type","`stype'")
			}
			
			// only one learner so it's the opt; set opt="" if no std stacking
			if `stdflag'	mata: add_learner_item(`ename',"opt","`m'","`vtilde'")
			else			mata: add_learner_item(`ename',"opt","`m'","")
			
		}
		else {
		
			// always label learner predicted values
			forvalues j=1/`nlearners' {
				qui label var `vtilde'0_L`j'_`m' "Pred. values E[`vname'|X] given `treatvar'==0 using base learner `j', rep `m'"
				qui label var `vtilde'1_L`j'_`m' "Pred. values E[`vname'|X] given `treatvar'==1 using base learner `j', rep `m'"
			}
			
			// save results relating to stacked learner if it exists
			if `stdflag' {
			
				cap drop `vtilde'0_`m'
				cap drop `vtilde'1_`m'
				// vtilde is predicted values
				qui gen `vtype' `vtilde'0_`m' = `vhat0'
				qui gen `vtype' `vtilde'1_`m' = `vhat1'
				qui label var `vtilde'0_`m' "Pred. values E[`vname'|X] given `treatvar'==0 using `cmd', rep `m'"
				qui label var `vtilde'1_`m' "Pred. values E[`vname'|X] given `treatvar'==1 using `cmd', rep `m'"
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
				mata: add_learner_item(`ename',"`vtilde'","stack_type","`stype'")
			}
			
			// only one learner so it's the opt; set opt="" if no std stacking
			forvalues t=0/1 {
				if `stdflag'	mata: add_learner_item(`ename',"opt`t'","`m'","`vtilde'")
				else			mata: add_learner_item(`ename',"opt`t'","`m'","")
			}
		}

		// add standard stacking base learner CV predictions
		if `stdflag' {
			if ~`tvflag'{	// case 1
				mata: add_result_item(`ename',"`vtilde'","y_stacking_cv", "`m'", `y_stacking_cv')
			}
			else {
				forvalues t=0/1 {
					mata: add_result_item(`ename',"`vtilde'","y_stacking_cv`t'", "`m'", `y_stacking_cv`t'')
				}
			}
		}
		
		// add shortstack results
		if `ssflag' {
			if ~`tvflag'{	// case 1
		
				label var `shortstack'_ss_`m' "Pred. values E[`vname'|X] using shortstacking, rep `m'"
				// calculate and return mspe and sample size
				tempvar vres_sq
				// shortstack macros have fitted values
				qui gen double `vres_sq' = (`vname' - `shortstack'_ss_`m')^2 if `touse'
			
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
				// save base estimator list with rest of shortstack results
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_base_est","`base_est'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_type","`stype'")
			}
			else {	// case 2
			
				label var `shortstack'_ss0_`m'  "Pred. values E[`vname'|X] given `treatvar'==0 using shortstacking, rep `m'"
				label var `shortstack'_ss1_`m'  "Pred. values E[`vname'|X] given `treatvar'==1 using shortstacking, rep `m'"
				// calculate and return mspe and sample size
				tempvar vres0_sq vres1_sq
				// shortstack macros have fitted values
				qui gen double `vres0_sq' = (`vname' - `shortstack'_ss0_`m')^2 if `treatvar' == 0 & `touse'
				qui gen double `vres1_sq' = (`vname' - `shortstack'_ss1_`m')^2 if `treatvar' == 1 & `touse'
		
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
				// save base estimator list with rest of shortstack results
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_base_est","`base_est'")
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`shortstack'_ss","ss_final_est", "`ssfinalest'")
				mata: add_learner_item(`ename',"`shortstack'_ss","stack_type","`stype'")
			}
		}
	
		// add poolstack results
		if `psflag' {
			if ~`tvflag' {	// case 1
		
				label var `poolstack'_ps_`m' "Pred. values E[`vname'|X] using poolstacking, rep `m'"
				// calculate and return mspe and sample size
				tempvar vres_sq
				// poolstack macros have fitted values
				qui gen double `vres_sq' = (`vname' - `poolstack'_ps_`m')^2 if `touse'
			
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
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`poolstack'_ps","ps_final_est", "`psfinalest'")
				mata: add_learner_item(`ename',"`poolstack'_ps","stack_type","`stype'")
				
			}
			else {	// case 2
			
				label var `poolstack'_ps0_`m'  "Pred. values E[`vname'|X] given `treatvar'==0 using poolstacking, rep `m'"
				label var `poolstack'_ps1_`m'  "Pred. values E[`vname'|X] given `treatvar'==1 using poolstacking, rep `m'"
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
					mata: add_result_item(`ename',"`poolstack'_ps","N`t'",             "`m'", `N`t'')
					mata: add_result_item(`ename',"`poolstack'_ps","N`t'_folds",       "`m'", st_matrix("`N`t'_folds'"))
					mata: add_result_item(`ename',"`poolstack'_ps","MSE`t'",           "`m'", `mse`t'')
					mata: add_result_item(`ename',"`poolstack'_ps","MSE`t'_folds",     "`m'", st_matrix("`mse`t'_folds'"))
					mata: add_result_item(`ename',"`poolstack'_ps","ps_weights`t'",    "`m'", st_matrix("`psw`t''"))
				}
				// final estimator used to stack and stack type are learner items
				mata: add_learner_item(`ename',"`poolstack'_ps","ps_final_est", "`psfinalest'")
				mata: add_learner_item(`ename',"`poolstack'_ps","stack_type","`stype'")
			}
		}
		
		// final clean up
		cap mata: mata drop `y_stacking_cv'
		cap mata: mata drop `y_stacking_cv0' `y_stacking_cv1'
		cap mata: mata drop `psw'
		cap mata: mata drop `psw1' `psw0'
		cap mata: mata drop `ssw'

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
	
	return scalar N			= `N'
	return local cmd_list	pystacked
	
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
							NOIsily					///
							allowallzero			/// in the LATE model: allow D
							ssfinalest(name)		/// final estimator for short-stacking
							finalest(name)			/// final estimator for both
							*						/// ignored options
							]

	// final estimator choice; default is NNLS + coefs sum to 1
	if "`finalest'"==""		local finalest nnls1
	if "`ssfinalest'"==""	local ssfinalest `finalest'
	
	mata: st_local("vname", `ename'.vname)
	mata: st_local("vtlist", invtokens(`ename'.vtlist))
	mata: st_local("nlearners", strofreal(`ename'.nlearners))
	mata: st_local("shortstack", `ename'.shortstack)
	mata: st_local("poolstack", `ename'.poolstack)
	mata: st_local("lieflag", strofreal(`ename'.lieflag))
	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
	** indicator for pooled-stacking
	local psflag	= "`poolstack'"~=""
	** indicator for interactive model
	local tvflag	= "`treatvar'"~=""
	
	if "`noisily'"=="" {
		local qui quietly
	}
	
	*** debugging message
	`qui' di as text "entering _crossfit_other..."

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
		di as err "error - shortstack option relevant only for multiple learners"
		exit 198
	}
	if `psflag' {
		di as err "error - poolstack option available only with pystacked; also unavailable for fiv model"
		exit 198
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
				qui gen double `shortstack'_ss_`m'=.
			}
			forvalues j=1/`nlearners' {
				tempvar vhat`j' vres`j'
				qui gen double `vhat`j''=.
				qui gen double `vres`j''=.
			}
		}  
		else if `tvflag' & ~`lieflag' { // case 2
			if `ssflag' {
				cap drop `shortstack'_ss0_`m'
				cap drop `shortstack'_ss1_`m'
				qui gen double `shortstack'_ss0_`m'=.
				qui gen double `shortstack'_ss1_`m'=.
			}
			forvalues j=1/`nlearners' {
				tempvar vhat0`j' vres0`j'
				tempvar vhat1`j' vres1`j'
				qui gen double `vhat0`j''=.
				qui gen double `vres0`j''=.	
				qui gen double `vhat1`j''=.
				qui gen double `vres1`j''=.	
			}
		}
		else if `lieflag' { // case 3
			// out-of-sample predicted values for E[D|XZ] 
			forvalues j=1/`nlearners' {
				tempvar dhat`j'
				qui gen double `dhat`j''=.  // using learner i
			}
			// out-of-sample predicted values for E[D|X] 
			forvalues j=1/`nlearners' {
				tempvar hhat`j'
				qui gen double `hhat`j''=.  // using learner i in both steps
			}
			// predicted values for E[D|ZX] for each k & learner i
			forvalues k=1/`kfolds' {
				forvalues j=1/`nlearners' {
					tempvar dhat_`j'_`k'
					qui gen double `dhat_`j'_`k''=.
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
					tempvar dhat_isSS_`k' 
					qui gen double `dhat_isSS_`k''=.
				}
				// out-of-sample predicted values for E[D|ZX] from short-stacking by k
				// NB: this is different from "dhatSS", which are from applying constrained regression to full sample
				tempvar dhat_oosSS
				qui gen double `dhat_oosSS'=. 
				// out-of-sample predicted values for E[D|X] 
				forvalues j=1/`nlearners' {
					tempvar hhatSS`j'
					qui gen double `hhatSS`j''=. // using short-stacking for E[D|XZ] & then using learner i for E[D|X]
				}
				// short-stacking predicted values
				tempvar dhatSS
				qui gen double `dhatSS'=. // this will become `shortstack'
				tempvar hhatSS
				qui gen double `hhatSS'=. // this will become `shortstack'_h_ss_
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
					
					// get fitted values for kth fold	
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
						`qui' gen double `vhat_k'=10e-12  if `fid'==`k' & `touse' 	
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
				`qui' _ddml_nnls `vname' `vhats', finalest(`ssfinalest') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw
				mat `ssw' = e(b)
				tempvar vtemp
				qui predict double `vtemp'
				qui replace `shortstack'_ss_`m' = `vtemp'
					
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
				`qui' _ddml_nnls `vname' `vhats1' if `treatvar'==1, finalest(`ssfinalest') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw1
				mat `ssw1' = e(b)
				qui predict double `vtemp'
				qui replace `shortstack'_ss1_`m'=`vtemp'
					
				// treatvar == 0
				forvalues j=1/`nlearners' {
					local vhats0 `vhats0' `vhat0`j''
				}
				tempvar vtemp
				`qui' di
				`qui' di as text as text "Short-stacking, finalest=`ssfinalest' (interactive model, treatvar=0):"
				`qui' _ddml_nnls `vname' `vhats0' if `treatvar'==0, finalest(`ssfinalest') `noisily'
				`qui' di as text "N=" as res e(N)
				tempname ssw0
				mat `ssw0' = e(b)
				qui predict double `vtemp'
				qui replace `shortstack'_ss0_`m'=`vtemp'
		
			}
			else if `lieflag' {
				// apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
				local dhats
				forvalues j=1/`nlearners' {
					local dhats `dhats' `dhat`j''
				}
				`qui' di
				`qui' di as text "Short-stacking, finalest=`ssfinalest' (LIE, OOS E[D|XZ]):"
				`qui' _ddml_nnls `vname' `dhats' if `touse'
				tempname ssw
				mat `ssw'= e(b)
				tempvar vtemp
				qui predict double `vtemp' if `touse'
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
					`qui' _ddml_nnls `vname' `dhats_is' if `fid'!=`k' & `touse' 
					qui predict double `vtemp'
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
						mata: st_local("vtype",return_learner_item(`ename',"`vtilde'","vtype"))				
	
						// replace {D}-placeholder in estimation string with variable name
						local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_isSS_`k''",1)
			
						// estimation	
						`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
						local cmd_h `e(cmd)'
					
						// get fitted values  
						tempvar vtemp
						qui predict double `vtemp' if `touse', `predopt_h'
			
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
					`qui' _ddml_nnls `dhat_oosSS' `hhatSS_list'
					if (`k'==1) {
						mat `ssw_h' = e(b)
					}
					else {
						mat `ssw_h_temp' = e(b)
						mat  `ssw_h_temp' = (`ssw_h' \ `ssw_h_temp' )
					}
					tempvar vtemp
					qui predict double `vtemp'
					qui replace `hhatSS'=`vtemp'
				}
				qui replace `shortstack'_ss_`m'=`dhatSS'
				qui replace `shortstack'_h_ss_`m'=`hhatSS'
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
			}
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
				qui replace `shortstack'_ss1_`m'=`vhat11'
				qui replace `shortstack'_ss0_`m'=`vhat01'
			}
			else if `lieflag' {
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
		
		forvalues j=1/`nlearners' {
			
			local vtilde	: word `j' of `vtlist'
			local cmd		: word `j' of `cmd_list'
			local cmd_h		: word `j' of `cmd_h_list'

			// vtilde, mspe, etc.
			if ~`tvflag' & ~`lieflag' {
		
				cap drop `vtilde'_`m'
				// vtilde is predicted values
				mata: st_local("vtype", return_learner_item(`ename',"`vtilde'","vtype"))
				qui gen `vtype' `vtilde'_`m' = `vhat`j''
				qui label var `vtilde'_`m' "Pred. values E[`vname'|X] using `cmd', rep `m'"
		
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
				// vtilde is predicted values
				mata: st_local("vtype", return_learner_item(`ename',"`vtilde'","vtype"))
				qui gen `vtype' `vtilde'0_`m' = `vhat0`j''
				qui gen `vtype' `vtilde'1_`m' = `vhat1`j''
				qui label var `vtilde'0_`m' "Pred. values E[`vname'|X] given `treatvar'==0 using `cmd', rep `m'"
				qui label var `vtilde'1_`m' "Pred. values E[`vname'|X] given `treatvar'==1 using `cmd', rep `m'"
		
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
				
				mata: st_local("vtype", return_learner_item(`ename',"`vtilde'","vtype"))
				qui gen `vtype' `vtilde'_`m' = `dhat`j''
				qui label var `vtilde'_`m' "Pred. values E[`vname'|X,Z], rep `m'"
				qui gen `vtype' `vtilde'_h_`m' = `hhat`j''
				qui label var `vtilde'_h_`m' "Pred. values E[`vtilde'|X], rep `m'"

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
		
				label var `shortstack'_ss_`m' "Pred. values E[`vname'|X] using shortstacking, rep `m'"
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
			
				label var `shortstack'_ss0_`m'  "Pred. values E[`vname'|X] given `treatvar'==0 using shortstacking, rep `m'"
				label var `shortstack'_ss1_`m'  "Pred. values E[`vname'|X] given `treatvar'==1 using shortstacking, rep `m'"
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
	
				label var `shortstack'_ss_`m' "Pred. values E[`vname'|Z,X] using shortstacking, rep `m'"
				label var `shortstack'_h_ss_`m' "Pred. values E[Dhat|X] of `vname' using shortstacking, rep `m'"
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
