* notes
* vtype was used for predicting values by fold, now applies to vtilde as well
* eqntype replaced by resid option (default=fitted)
* default is additive-type crossfitting; treatvar option triggers interactive-type crossfitting

program define crossfit, rclass sortpreserve

	syntax [anything] [if] [in] ,					/// [anything] is renamed to vname below; currently undocumented option
							[ estring(string asis)	/// estimation string
													/// need asis option in case it includes strings
													///
							ename(name)				/// name for Mata struct; default is "crossfit"
							NOREPLACE				///
							foldvar(varlist)		/// one fold var per resample
							kfolds(integer 5)		/// ignored if foldvars provided
							reps(integer 1)			/// ignored if foldvars provided
							NORANDOM				/// first fold ID uses obs in existing order
							resid					/// return residuals (prediction errors); default is predicted values
							vtilde(namelist)		/// name(s) of fitted variable(s)
							Generate(namelist)		/// synonym for vtilde
							vtype(string)			/// datatype of fitted variable; default=double
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							shortstack(name)		///
													/// 
													/// options specific to LIE/DDML-IV
							estringh(string asis)	/// est string for E[D^|XZ]		
													/// 
							predopt(string asis)	/// undocumented
							NOIsily					///
							allowallzero			/// in the LATE model: allow D
							]

	// renaming for clarity
	local vtlist `vtilde' `generate'
	local vname `anything'
	// clear the local macro
	local vtilde
	local generate

	** indicator for LIE/optimal-IV model
	local lieflag	= "`estringh'"~=""
	** indicator for interactive model
	local tvflag	= "`treatvar'"~=""
	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
	** indicator for initializing eqn struct
	local initflag	= "`noreplace'"==""

	// LIE => we want predicted values not resids
	if `lieflag' & "`resid'"~="" {
		di as res "resid option ignored"
		local resid
	}
	
	if "`noisily'"=="" {
		local qui quietly
	}

	if `initflag' {
		// start with a blank setup
		
		*** set up struct with equation/learner info
		
		if "`ename'"== "" {
			local eqn_info crossfit
		}
		else {
			local eqn_info `ename'
		}
		
		*** variable type
		if "`vtype'"=="" {
			local vtype double
		}
		if "`vtype'"=="none" {
			local vtype
		}
	
		mata: `eqn_info' = init_eStruct()
		initialize_eqn_info,										///
							ename(`eqn_info')						///
							vname(`vname')							///
							vtlist(`vtlist')						///
							shortstack(`shortstack')				///
							estring(`estring') estringh(`estringh')	///
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

	}
	else {
		// ename is a pre-populated eqn struct
		local eqn_info `ename'
		mata: st_local("vname", `eqn_info'.vname)
		mata: st_local("vtlist", invtokens(`eqn_info'.vtlist))
		mata: st_local("nlearners", strofreal(`eqn_info'.nlearners))
		mata: st_local("shortstack", `eqn_info'.shortstack)
		mata: st_local("lieflag", strofreal(`eqn_info'.lieflag))
		local ssflag	= "`shortstack'"~=""
	}
	
	// check
	if `ssflag' & ~(`nlearners' > 1) {
		di as err "warning - shortstack option ignored - #learners must be > 1"
		// clear macro, reset flag
		local shortstack
		local ssflag = 0
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
	}
	else {
		// fold variables not provided, so generate
		// if neither foldvar nor reps provided, set reps to default
		forvalues m=1/`reps' {
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
	
	// datatype for fitted values/residuals
	if "`vtype'"=="" {
		local vtype double
	}

	// stored results; initialize here
	tempname mse_list N_list mse_folds_list N_folds_list
	tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
	tempname mse0_list N0_list mse0_folds_list N0_folds_list
	tempname mse1_list N1_list mse1_folds_list N1_folds_list
		
	// loop over resamples (crossfitting, shortstacking, store results)
	forvalues m=1/`reps' {
	
	// create rnames for renaming the rows of saved matrices
	local rnames `rnames' resample_`m'

	******************************** CROSSFITTING ************************************
	
	// cross-fitting fundamentally depends on three cases 
	// case 1: `tvflag'==0 & `lieflag'==0
	// case 2: `tvflag'==1 & `lieflag'==0
	// case 3: `lieflag'==1
	// nb: `tvflag'==1 & `lieflag'==1 is impossible
	
		if `reps'>1 {
			di as text "Resample `m'..."
		}

		local fid : word `m' of `fidlist'

		// create blank fitted variable(s)
		if ~`tvflag' & ~`lieflag' { // case 1
			if `ssflag' {
				cap drop `shortstack'_`m'
				qui gen `vtype' `shortstack'_`m'=.
			}
			forvalues i=1/`nlearners' {
				tempvar vhat`i' vres`i'
				qui gen `vtype' `vhat`i''=.
				qui gen `vtype' `vres`i''=.
			}
		}  
		else if `tvflag' & ~`lieflag' { // case 2
			if `ssflag' {
				cap drop `shortstack'0_`m'
				cap drop `shortstack'1_`m'
				qui gen `vtype' `shortstack'0_`m'=.
				qui gen `vtype' `shortstack'1_`m'=.
			}
			forvalues i=1/`nlearners' {
				tempvar vhat0`i' vres0`i'
				tempvar vhat1`i' vres1`i'
				qui gen `vtype' `vhat0`i''=.
				qui gen `vtype' `vres0`i''=.	
				qui gen `vtype' `vhat1`i''=.
				qui gen `vtype' `vres1`i''=.	
			}
		}
		else if `lieflag' { // case 3
			// out-of-sample predicted values for E[D|XZ] 
			forvalues i=1/`nlearners' {
				tempvar dhat`i'
				qui gen `vtype' `dhat`i''=.  // using learner i
			}
			// out-of-sample predicted values for E[D|X] 
			forvalues i=1/`nlearners' {
				tempvar hhat`i'
				qui gen `vtype' `hhat`i''=.  // using learner i in both steps
			}
			// predicted values for E[D|ZX] for each k & learner i
			forvalues k=1/`kfolds' {
				forvalues i=1/`nlearners' {
					tempvar dhat_`i'_`k'
					qui gen `vtype' `dhat_`i'_`k''=.
				}
			}
			if `ssflag' { // with short-stacking
				// for final results
				cap drop `shortstack'_`m'
				cap drop `shortstack'_h_`m'
				qui gen `vtype' `shortstack'_`m'=.
				qui gen `vtype' `shortstack'_h_`m'=.			
				// in-sample predicted values for E[D|ZX] for each k: short-stacked
				forvalues k=1/`kfolds' {		
					tempvar dhat_isSS_`k' 
					qui gen `vtype' `dhat_isSS_`k''=.
				}
				// out-of-sample predicted values for E[D|ZX] from short-stacking by k
				// NB: this is different from "dhatSS", which are from applying contrained regression to full sample
				tempvar dhat_oosSS 
				qui gen `vtype' `dhat_oosSS'=. 
				// out-of-sample predicted values for E[D|X] 
				forvalues i=1/`nlearners' {
					tempvar hhatSS`i'
					qui gen `vtype' `hhatSS`i''=. // using short-stacking for E[D|XZ] & then using learner i for E[D|X]
				}
				// short-stacking predicted values
				tempvar dhatSS
				qui gen `vtype' `dhatSS'=. // this will become `shortstack'
				tempvar hhatSS
				qui gen `vtype' `hhatSS'=. // this will become `shortstack'_h
			} 
		}
		else {
			di as err "internal crossfit error"
			exit 198
		}
	
		// will save pystacked weights
		forvalues i=1/`nlearners' {
			tempname pysw_`i' pysw0_`i' pysw1_`i' pyswh_`i'
			tempname pysw_temp_`i' pysw1_temp_`i' pysw0_temp_`i' pyswh_temp_`i'
		}
	
		// crossfit
		di as text "Cross-fitting fold " _c
		forvalues k = 1(1)`kfolds' {
	
			di as text "`k' " _c
	
			forvalues i=1/`nlearners' {
				local vtilde : word `i' of `vtlist'
				mata: st_local("est_main", return_learner_item(`eqn_info',"`vtilde'","est_main"))
				mata: st_local("est_options", return_learner_item(`eqn_info',"`vtilde'","est_options"))
				mata: st_local("predopt",return_learner_item(`eqn_info',"`vtilde'","predopt"))
				mata: st_local("vtype",return_learner_item(`eqn_info',"`vtilde'","vtype"))
				if `lieflag' {
					// LIE locals
					local hhat `vtilde'_h
					mata: st_local("est_main_h", return_learner_item(`eqn_info',"`vtilde'","est_main_h"))
					mata: st_local("est_options_h", return_learner_item(`eqn_info',"`vtilde'","est_options_h"))
					mata: st_local("predopt_h",return_learner_item(`eqn_info',"`vtilde'","predopt_h"))
					mata: st_local("vtype_h",return_learner_item(`eqn_info',"`vtilde'","vtype_h"))			
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
					
					// save pystacked weights
					if ("`cmd'"=="pystacked") {
						if (`k'==1) {
							mat `pysw_`i'' = e(weights)
						}
						else {
							mat `pysw_temp_`i'' = e(weights)
							mat `pysw_`i'' = (`pysw_`i'',`pysw_temp_`i'')
						}
					}
		
					// get fitted values and residuals for kth fold	
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'
		
					// get predicted values
					qui replace `vhat`i'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres`i'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
					// pystacked learners
					if (("`cmd'"=="pystacked") & (`k'==1) & (`m'==1)) {
							// holds for all reps and folds
							local base_est_`i' `e(base_est)'
					}
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
		
					// save pystacked weights
					if ("`cmd'"=="pystacked") {
						if (`k'==1) {
							mat `pysw1_`i'' = e(weights)
						}
						else {
							mat `pysw1_temp_`i'' = e(weights)
							mat `pysw1_`i'' = (`pysw1_`i'',`pysw1_temp_`i'')
						}
					}
		
					// get fitted values for kth fold	
					tempvar vhat_k
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'
					qui replace `vhat1`i'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres1`i'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
		
					// for treatvar = 0

					// first, we need to account for D always = 0 if Z=0
					// check if dvar is always 0
					`qui' count if `treatvar'==0 & `vname'!=0 & `fid'==`k' & `touse' 

					if (`r(N)'==0 & "`allowallzero'"!="") {
						// get fitted values for kth fold	
						tempvar vhat_k	
						`qui' gen `vhat_k'=10e-12  if `fid'==`k' & `touse' 	
					} 
					else {
						// estimate excluding kth fold
						`qui' `est_main' if `fid'!=`k' & `treatvar' == 0 & `touse', `est_options'
			
						// save pystacked weights
						if ("`cmd'"=="pystacked") {
							if (`k'==1) {
								mat `pysw0_`i'' = e(weights)
							}
							else {
								mat `pysw0_temp_`i'' = e(weights)
								mat `pysw0_`i'' = (`pysw0_`i'',`pysw0_temp_`i'')
							}
						}

						// get fitted values for kth fold	
						tempvar vhat_k
						qui predict `vtype' `vhat_k' if `fid'==`k' & `touse', `predopt'					

					}
					qui replace `vhat0`i'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres0`i'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
					
					// pystacked learners
					if (("`cmd'"=="pystacked") & (`k'==1) & (`m'==1)) {
							// holds for all reps and folds
							local base_est_`i' `e(base_est)'
					}
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
		
					// get pystacked weights
					if ("`cmd'"=="pystacked") {
						if ((`k'==1) & (`m'==1)) {
							// holds for all reps and folds
							local base_est_`i' `e(base_est)'
						}
						if (`k'==1) {
							mat `pysw_`i'' = e(weights)
						}
						else {
							mat `pysw_temp_`i'' = e(weights)
							mat `pysw_`i'' = (`pysw_`i'',`pysw_temp_`i'')
						}
					}
					
					// get fitted values (in and out of sample)
					qui predict `vtype' `vhat_k' if `touse', `predopt'
		
					// get *combined* out-of-sample predicted values
					qui replace `dhat`i'' = `vhat_k' if `fid'==`k' & `touse'
		
					// get predicted values in wide format for each i and k
					qui replace `dhat_`i'_`k'' = `vhat_k' if `touse'
		
					// Step II: estimation of E[D^|X]
		
					// replace {D}-placeholder in estimation string with variable name
					local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_`i'_`k''",1)
		
					// estimation	
					`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
					local cmd_h `e(cmd)'
		
					// get pystacked weights
					if ("`cmd_h'"=="pystacked") {
						if ((`k'==1) & (`m'==1)) {
							// holds for all reps and folds
							local base_est_h_`i' `e(base_est)'
						}
						if (`k'==1) {
							// initialize
							mat `pyswh_`i'' = e(weights)
						}
						else {
							mat `pyswh_temp_`i'' = e(weights)
							mat `pyswh_`i'' = (`pyswh_`i'',`pyswh_temp_`i'')
						}
					}
		
					// get fitted values  
					qui predict `vtype_h' `vtil_k' if `touse', `predopt_h'
		
					// get *combined* out-of-sample predicted values
					qui replace `hhat`i'' = `vtil_k' if `fid'==`k' & `touse'
		
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
				forvalues i=1/`nlearners' {
					local vhats `vhats' `vhat`i''
				}
				tempvar vss
				`qui' di
				`qui' di as text "Stacking NNLS (additive model):"
				`qui' _ddml_nnls `vname' `vhats'
				tempname ssw
				mat `ssw' = e(b)
				tempvar vtemp
				qui predict `vtype' `vtemp'
				qui replace `shortstack'_`m' = `vtemp'
					
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'_`m' = `vname' - `shortstack'_`m'
				}
			}
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
				local vhats1 
				local vhats0
				// treatvar == 1
				forvalues i=1/`nlearners' {
					local vhats1 `vhats1' `vhat1`i''
				}
				tempvar vtemp
				`qui' di
				`qui' di as text "Stacking NNLS (interactive model, treatvar=1):"
				`qui' _ddml_nnls `vname' `vhats1' if `treatvar'==1
				tempname ssw0
				mat `ssw0' = e(b)
				qui predict `vtype' `vtemp'
				qui replace `shortstack'1_`m'=`vtemp'
					
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'1_`m' = `vname' - `shortstack'1_`m'
				}
				// treatvar == 0
				forvalues i=1/`nlearners' {
					local vhats0 `vhats0' `vhat0`i''
				}
				tempvar vtemp
				`qui' di
				`qui' di as text "Stacking NNLS (interactive model, treatvar=0):"
				`qui' _ddml_nnls `vname' `vhats0' if `treatvar'==0
				tempname ssw1
				mat `ssw1' = e(b)
				qui predict `vtype' `vtemp'
				qui replace `shortstack'0_`m'=`vtemp'
		
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'0_`m' = `vname' - `shortstack'0_`m'
				}
			}
			else if `lieflag' {
				// apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
				local dhats
				forvalues i=1/`nlearners' {
					local dhats `dhats' `dhat`i''
				}
				`qui' di
				`qui' di as text "Stacking NNLS (LIE, OOS E[D|XZ]):"
				`qui' _ddml_nnls `vname' `dhats' if `touse'
				tempname ssw
				mat `ssw'= e(b)
				tempvar vtemp
				qui predict `vtype' `vtemp' if `touse'
				qui replace `dhatSS'=`vtemp' 
	
				// apply short-stacking to in-sample predicted values of E[D|XZ] *for each k*
				forvalues k = 1(1)`kfolds' {
					local dhats_is
					forvalues i=1/`nlearners' {
						local dhats_is `dhats_is' `dhat_`i'_`k''
					}
					tempvar vtemp
					`qui' di
					`qui' di as text "Stacking NNLS (LIE, in-sample E[D|XZ] fold `k':"
					`qui' _ddml_nnls `vname' `dhats_is' if `fid'!=`k' & `touse' 
					qui predict `vtype' `vtemp'
					qui replace `dhat_isSS_`k'' = `vtemp' if `fid'!=`k' & `touse'
					qui replace `dhat_oosSS' = `vtemp' if `fid'==`k' & `touse'
				}
	
				// need to cross-fit stacked in-sample predicted values against X
				forvalues k = 1(1)`kfolds' {
					forvalues i=1/`nlearners' {
						local vtilde : word `i' of `vtlist'
						mata: st_local("est_main_h", return_learner_item(`eqn_info',"`vtilde'","est_main_h"))
						mata: st_local("est_options_h", return_learner_item(`eqn_info',"`vtilde'","est_options_h"))
						mata: st_local("predopt_h",return_learner_item(`eqn_info',"`vtilde'","predopt_h"))
						mata: st_local("vtype_h",return_learner_item(`eqn_info',"`vtilde'","vtype_h"))				
	
						// replace {D}-placeholder in estimation string with variable name
						local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_isSS_`k''",1)
			
						// estimation	
						`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
						local cmd_h `e(cmd)'
					
						// get fitted values  
						tempvar vtemp
						qui predict `vtype_h' `vtemp' if `touse', `predopt_h'
			
						// get out-of-sample predicted values
						qui replace `hhatSS`i'' = `vtemp' if `fid'==`k' & `touse'
					}
				}
	
				// final stacking for E[D|X]
				tempname sswh sswh_temp
				forvalues k = 1(1)`kfolds' {
					local hhatSS_list
					forvalues i=1/`nlearners' {
						local hhatSS_list `hhatSS_list' `hhatSS`i''
					}
					`qui' di
					`qui' di as text "Stacking NNLS (LIE, E[D|X]):"
					`qui' _ddml_nnls `dhat_oosSS' `hhatSS_list'
					if (`k'==1) {
						mat `sswh' = e(b)
					}
					else {
						mat `sswh_temp' = e(b)
						mat  `sswh_temp' = (`sswh' \ `sswh_temp' )
					}
					tempvar vtemp
					qui predict `vtype' `vtemp'
					qui replace `hhatSS'=`vtemp'
				}
				qui replace `shortstack'_`m'=`dhatSS'
				label var `shortstack'_`m' "short-stacking cross-fitted E[D|Z,X]"
				qui replace `shortstack'_h_`m'=`hhatSS'
				label var `shortstack'_h_`m' "short-stacking cross-fitted E[D|X]"
			}
			else {
				di as err "internal crossfit error"
				exit 198
			}
		}
		else if `ssflag' {
			// single learner case, so shortstack vars are just copies of learner vars
			
			if ~`tvflag' & ~`lieflag' { // case 1
				qui replace `shortstack'_`m' = `vhat1'
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'_`m' = `vname' - `shortstack'_`m'
				}
			}
			else if `tvflag' & ~`lieflag' {	// case 2: interactive models
				qui replace `shortstack'1_`m'=`vhat11'
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'1_`m' = `vname' - `shortstack'1_`m' if `treatvar'==1
				}
				qui replace `shortstack'0_`m'=`vhat01'
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'0_`m' = `vname' - `shortstack'0_`m' if `treatvar'==0
				}
			}
			else if `lieflag' {
				qui replace `shortstack'_`m'=`dhat1'
				label var `shortstack'_`m' "short-stacking cross-fitted E[D|Z,X]"
				qui replace `shortstack'_h_`m'=`hhat1'
				label var `shortstack'_h_`m' "short-stacking cross-fitted E[D|X]"
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
		
		forvalues i=1/`nlearners' {
			
			local vtilde	: word `i' of `vtlist'
			local cmd		: word `i' of `cmd_list'
			local cmd_h		: word `i' of `cmd_h_list'

			// vtilde, mspe, etc.
			if ~`tvflag' & ~`lieflag' {
		
				cap drop `vtilde'_`m'
				if "`resid'"=="" {
					// vtilde is predicted values
					qui gen `vtilde'_`m' = `vhat`i''
					qui label var `vtilde'_`m' "Predicted values cond. exp. of `vname' using `cmd'"
				}
				else {
					// vtilde is residuals
					qui gen `vtilde'_`m' = `vres`i''
					qui label var `vtilde'_`m' "Residuals cond. exp. of `vname' using `cmd'"
				}
		
				// calculate and return mspe and sample size
				tempvar vres_sq
				qui gen double `vres_sq' = `vres`i''^2 if `touse'

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
			
				mata: add_result_item(`eqn_info',"`vtilde'","N",         "`m'", `N')
				mata: add_result_item(`eqn_info',"`vtilde'","N_folds",   "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`eqn_info',"`vtilde'","MSE",       "`m'", `mse')
				mata: add_result_item(`eqn_info',"`vtilde'","MSE_folds", "`m'", st_matrix("`mse_folds'"))
				
				if "`cmd'"=="pystacked" {
					mata: add_result_item(`eqn_info',"`vtilde'","stack_weights","`m'", st_matrix("`pysw_`i''"))
					mata: add_learner_item(`eqn_info',"`vtilde'","stack_base_est","`base_est_`i''")
				}
				
				if `i'==1 {
					local mse_opt		= `mse'
					mata: add_learner_item(`eqn_info',"opt","`m'","`vtilde'")
				}
				else if `mse' < `mse_opt' {
					// overwrite with new opt
					local mse_opt		= `mse'
					mata: add_learner_item(`eqn_info',"opt","`m'","`vtilde'")
				}
				
			}
			else if `tvflag' & ~`lieflag' {
			
				cap drop `vtilde'0_`m'
				cap drop `vtilde'1_`m'
				if "`resid'"=="" {
					// vtilde is predicted values
					qui gen `vtilde'0_`m' = `vhat0`i''
					qui gen `vtilde'1_`m' = `vhat1`i''
					qui label var `vtilde'0_`m' "Predicted values cond. exp. of `vname' given `treatvar'==0 using `cmd'"
					qui label var `vtilde'1_`m' "Predicted values cond. exp. of `vname' given `treatvar'==1 using `cmd'"
				}
				else {
					// vtilde is residuals
					qui gen `vtilde'0_`m' = `vres0`i''
					qui gen `vtilde'1_`m' = `vres1`i''
					qui label var `vtilde'0_`m' "Residuals cond. exp. of `vname' given `treatvar'==0 using `cmd'"
					qui label var `vtilde'1_`m' "Residuals cond. exp. of `vname' given `treatvar'==1 using `cmd'"
				}
		
				// calculate and return mspe and sample size
				tempvar vres0_sq vres1_sq
				// vtilde has fitted values
				qui gen double `vres0_sq' = (`vname' - `vhat0`i'')^2 if `treatvar' == 0 & `touse'
				qui gen double `vres1_sq' = (`vname' - `vhat1`i'')^2 if `treatvar' == 1 & `touse'
		
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
					mata: add_result_item(`eqn_info',"`vtilde'","N`t'",         "`m'", `N`t'')
					mata: add_result_item(`eqn_info',"`vtilde'","N`t'_folds",   "`m'", st_matrix("`N`t'_folds'"))
					mata: add_result_item(`eqn_info',"`vtilde'","MSE`t'",       "`m'", `mse`t'')
					mata: add_result_item(`eqn_info',"`vtilde'","MSE`t'_folds", "`m'", st_matrix("`mse`t'_folds'"))
				}
				
				if "`cmd'"=="pystacked" {
					mata: add_learner_item(`eqn_info',"`vtilde'","stack_base_est","`base_est_`i''")
					mata: add_result_item(`eqn_info',"`vtilde'","stack_weights0","`m'", st_matrix("`pysw0_`i''"))
					mata: add_result_item(`eqn_info',"`vtilde'","stack_weights1","`m'", st_matrix("`pysw1_`i''"))
				}
				
				forvalues t=0/1 {
					if `i'==1 {
						local mse`t'_opt		= `mse`t''
						mata: add_learner_item(`eqn_info',"opt`t'","`m'","`vtilde'")
					}
					else if `mse`t'' < `mse`t'_opt' {
						// overwrite with new opt
						local mse`t'_opt		= `mse`t''
						mata: add_learner_item(`eqn_info',"opt`t'","`m'","`vtilde'")
					}
				}
			}
			else if `lieflag' {
	
				cap drop `vtilde'_`m'
				cap drop `vtilde'_h_`m'
				
				qui gen `vtilde'_`m' = `dhat`i''
				qui label var `vtilde'_`m' "Predicted values E[`vname'|X,Z]"
				qui gen `vtilde'_h_`m' = `hhat`i''
				qui label var `vtilde'_h_`m' "Predicted values E[`vtilde'_`m'|X]"

				// calculate and return mspe and sample size
				tempvar hres dres hres_sq dres_sq
				// vtilde has fitted values
				qui gen double `dres_sq' = (`vname' - `dhat`i'')^2 if `touse'
				qui gen double `hres_sq' = (`vname' - `hhat`i'')^2 if `touse'
	
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
				
				mata: add_result_item(`eqn_info',"`vtilde'","N",         "`m'", `N')
				mata: add_result_item(`eqn_info',"`vtilde'","N_folds",   "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`eqn_info',"`vtilde'","MSE",       "`m'", `mse')
				mata: add_result_item(`eqn_info',"`vtilde'","MSE_folds", "`m'", st_matrix("`mse_folds'"))
	
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
				
				mata: add_result_item(`eqn_info',"`vtilde'","N_h",         "`m'", `N_h')
				mata: add_result_item(`eqn_info',"`vtilde'","N_h_folds",   "`m'", st_matrix("`N_h_folds'"))
				mata: add_result_item(`eqn_info',"`vtilde'","MSE_h",       "`m'", `mse_h')
				mata: add_result_item(`eqn_info',"`vtilde'","MSE_h_folds", "`m'", st_matrix("`mse_h_folds'"))
				
				// MS: this line fails if cmd hasn't been stored by ddml init; will happen if crossfit called directly 
				// mata: st_local("cmd",return_learner_item(`eqn_info',"`vtilde'","cmd"))							
				// instead use cmd selected from cmd_list above
				if "`cmd'"=="pystacked" {
					mata: add_result_item(`eqn_info',"`vtilde'","stack_weights","`m'", st_matrix("`pysw_`i''"))
					mata: add_learner_item(`eqn_info',"`vtilde'","stack_base_est","`base_est_`i''")					
				}
				if "`cmd_h'"=="pystacked" {
					mata: add_result_item(`eqn_info',"`vtilde'","stack_weights_h","`m'", st_matrix("`pyswh_`i''"))
					mata: add_learner_item(`eqn_info',"`vtilde'","stack_base_est_h","`base_est_h_`i''")					
				}
				
				// optimal D and H
				if `i'==1 {
					local mse_opt		= `mse'
					mata: add_learner_item(`eqn_info',"opt","`m'","`vtilde'")
				}
				else if `mse' < `mse_opt' {
					// overwrite with new opt
					local mse_opt		= `mse'
					mata: add_learner_item(`eqn_info',"opt","`m'","`vtilde'")
				}
				if `i'==1 {
					local mse_h_opt		= `mse_h'
					mata: add_learner_item(`eqn_info',"opt_h","`m'","`vtilde'_h")
				}
				else if `mse_h' < `mse_h_opt' {
					// overwrite with new opt
					local mse_h_opt		= `mse_h'
					mata: add_learner_item(`eqn_info',"opt_H","`m'","`vtilde'_h")
				}
				
			}

		}
		
		// add shortstack results
		if `ssflag' {
			if ~`tvflag' & ~`lieflag' {	// case 1
		
				// calculate and return mspe and sample size
				tempvar vres_sq
				// shortstack macro has the residuals
				qui gen double `vres_sq' = (`shortstack'_`m')^2 if `touse'
			
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

				mata: add_result_item(`eqn_info',"`shortstack'","N",         "`m'", `N')
				mata: add_result_item(`eqn_info',"`shortstack'","N_folds",   "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`eqn_info',"`shortstack'","MSE",       "`m'", `mse')
				mata: add_result_item(`eqn_info',"`shortstack'","MSE_folds", "`m'", st_matrix("`mse_folds'"))
				mata: add_result_item(`eqn_info',"`shortstack'","ss_flag", "`m'", 1)
				mata: add_result_item(`eqn_info',"`shortstack'","ss_weights", "`m'", st_matrix("`ssw'"))
				
			}
			else if `tvflag' & ~`lieflag' {	// case 2
			
				// calculate and return mspe and sample size
				tempvar vres0_sq vres1_sq
				// shortstack macros have fitted values
				qui gen double `vres0_sq' = (`shortstack'0_`m')^2 if `treatvar' == 0 & `touse'
				qui gen double `vres1_sq' = (`shortstack'1_`m')^2 if `treatvar' == 1 & `touse'
		
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
				
				mata: add_result_item(`eqn_info',"`shortstack'","ss_flag", "`m'", 1)
				forvalues t=0/1 {
					mata: add_result_item(`eqn_info',"`shortstack'","N`t'",         "`m'", `N`t'')
					mata: add_result_item(`eqn_info',"`shortstack'","N`t'_folds",   "`m'", st_matrix("`N`t'_folds'"))
					mata: add_result_item(`eqn_info',"`shortstack'","MSE`t'",       "`m'", `mse`t'')
					mata: add_result_item(`eqn_info',"`shortstack'","MSE`t'_folds", "`m'", st_matrix("`mse`t'_folds'"))
					mata: add_result_item(`eqn_info',"`shortstack'","ss_weights`t'", "`m'", st_matrix("`ssw`t''"))
				}
			}
			else if `lieflag' {	// case 3
	
				// calculate and return mspe and sample size
				tempvar hres dres hres_sq dres_sq
				// vtilde has fitted values
				qui gen double `dres_sq' = (`vname' - `shortstack'_`m')^2 if `touse'
				qui gen double `hres_sq' = (`vname' - `shortstack'_h_`m')^2 if `touse'
	
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

				mata: add_result_item(`eqn_info',"`shortstack'","N",         "`m'", `N')
				mata: add_result_item(`eqn_info',"`shortstack'","N_folds",   "`m'", st_matrix("`N_folds'"))
				mata: add_result_item(`eqn_info',"`shortstack'","MSE",       "`m'", `mse')
				mata: add_result_item(`eqn_info',"`shortstack'","MSE_folds", "`m'", st_matrix("`mse_folds'"))
	
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

				mata: add_result_item(`eqn_info',"`shortstack'","N_h",         "`m'", `N_h')
				mata: add_result_item(`eqn_info',"`shortstack'","N_h_folds",   "`m'", st_matrix("`N_h_folds'"))
				mata: add_result_item(`eqn_info',"`shortstack'","MSE_h",       "`m'", `mse_h')
				mata: add_result_item(`eqn_info',"`shortstack'","MSE_h_folds", "`m'", st_matrix("`mse_h_folds'"))
				mata: add_result_item(`eqn_info',"`shortstack'","ss_flag", "`m'", 1)
				mata: add_result_item(`eqn_info',"`shortstack'","ss_weights", "`m'", st_matrix("`ssw'"))
				mata: add_result_item(`eqn_info',"`shortstack'","ss_weights_h", "`m'", st_matrix("`sswh'"))
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
	forvalues i=1/`nlearners' {
		local vtilde : word `i' of `vtlist'
		cap return `vtilde'_pysw = `pysw_`i''
		cap return `vtilde'_pysw = `pysw0_`i''
		cap return `vtilde'_pysw = `pysw1_`i''
		cap return `vtilde'_pysw = `pyswh_`i''
	}
	return local cmd_list	`cmd_list'
	return local cmd_h_list	`cmd_h_list'
	
end

program define check_foldvar, rclass
	syntax [anything], fidlist(varlist) touse(varname)
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
		qui tab `fid'
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

program define initialize_eqn_info, rclass

	syntax [anything] [if] [in] ,					/// 
							[						///
							ename(name)				/// name of mata struct
							vname(varname)			/// name of dep var to be orthogonalized
							vtlist(string)			/// names of corresponding tilde variables
							shortstack(name)		/// name of shortstacked variable; may be empty
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
			mata: add_learner_item(`ename',"`vtilde'","vtype_h","`vtype'")
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
