* notes
* vtype was used for predicting values by fold, now applies to vtilde as well
* eqntype replaced by resid option (default=fitted)
* default is additive-type crossfitting; treatvar option triggers interactive-type crossfitting



* single equation version returns scalars:
*   r(mse)
*   r(N)
* and matrices:
*   r(N_folds)
*   r(mse_folds)

* multiple equation version returns matrices:
*   r(N_folds_list)
*   r(mse_folds_list)
*   r(N_list)
*   r(mse_list)
*   r(N_folds)
*   r(mse_folds)


mata:

struct eStruct {
	real matrix			vtlist
	real matrix			vtlist0
	real matrix			vtlist1
	real matrix			elist
	real matrix			emlist
	real matrix			eolist
	real matrix			vtlisth
	real matrix			elisth
	real matrix			emlisth
	real matrix			eolisth
	real scalar			nreps			// number of resamplings
	pointer rowvector	rlist
}

struct eStruct init_eStruct(real scalar reps)
{
	struct eStruct scalar	d

	d.vtlist	= J(1,0,"")
	d.vtlist0	= J(1,0,"")
	d.vtlist1	= J(1,0,"")
	d.elist		= J(1,0,"")
	d.emlist	= J(1,0,"")
	d.eolist	= J(1,0,"")
	d.vtlisth	= J(1,0,"")
	d.elisth	= J(1,0,"")
	d.emlisth	= J(1,0,"")
	d.eolisth	= J(1,0,"")
	d.nreps		= reps
	return(d)
}

struct rStruct {
	real colvector		N
	real colvector		N0
	real colvector		N1
	real colvector		N_h				// (intended for LIE)
	real matrix			N_folds			// sample size by fold; col=fold, row=learner
	real matrix			N0_folds		// sample size by fold; col=fold, row=learner
	real matrix			N1_folds		// sample size by fold; col=fold, row=learner
	real matrix         N_h_folds		// (intended for LIE)
	real colvector		MSE
	real colvector 		MSE0
	real colvector 		MSE1
	real matrix			MSE_folds		// MSE by fold; col=fold, row=learner
	real matrix			MSE0_folds		// MSE by fold; col=fold, row=learner
	real matrix			MSE1_folds		// MSE by fold; col=fold, row=learner
	real colvector      MSE_h 			// (intended for LIE)
	real matrix			MSE_h_folds		// (intended for LIE)
	real matrix 		stack_weights	// weights from use of pystacked
	real matrix 		stack_weights0	// weights from use of pystacked
	real matrix 		stack_weights1	// weights from use of pystacked
	real matrix 		stack_weights_h	// weights from use of pystacked
}

struct rStruct init_rStruct(real scalar k)
{
	// rows are resamplings
	struct rStruct scalar		r
	r.N				= J(0,1,0)
	r.N0			= J(0,1,0)
	r.N1			= J(0,1,0)
	r.N_h			= J(0,1,0)
	r.N_folds		= J(0,k,0)
	r.N0_folds		= J(0,k,0)
	r.N1_folds		= J(0,k,0)
	r.N_h_folds		= J(0,k,0)
	r.MSE			= J(0,1,0)
	r.MSE0			= J(0,1,0)
	r.MSE1			= J(0,1,0)
	r.MSE_h			= J(0,1,0)
	r.MSE_folds		= J(0,k,0)
	r.MSE0_folds	= J(0,k,0)
	r.MSE1_folds	= J(0,k,0)
	r.MSE_h_folds	= J(0,k,0)
	return(r)
}

end

program define initialize_eqn_info, rclass

	syntax [anything] [if] [in] ,					/// 
							[						///
							sname(name)				/// name of mata struct
							reps(integer 1)			/// if not supplied, default is 1
							vtlist(string)			/// names of corresponding tilde variables
							vtlist0(string)			/// names of corresponding tilde variables
							vtlist1(string)			/// names of corresponding tilde variables
							estring(string asis)	/// names of estimation strings
													/// need asis option in case it includes strings
							vtlisth(string)			/// intended for lists of E[D^|X] where D^=E[D|XZ]=vtilde()
							estringh(string asis)	/// names of LIE estimation strings
													/// need asis option in case it includes strings
							NOIsily					///
							]
	
	if "`noisily'"=="" {
		local qui quietly
	}

	tempname t
	
	mata: `sname'			= init_eStruct(`reps')
	// in two steps, to accommodate singleton lists (which are otherwise string scalars and not matrices
	mata: `t' = tokens("`vtlist'")
	mata: `sname'.vtlist	= `t'
	mata: `t' = tokens("`vtlist0'")
	mata: `sname'.vtlist0	= `t'
	mata: `t' = tokens("`vtlist1'")
	mata: `sname'.vtlist1	= `t'
	parse_estring, sname(`sname') estring(`estring') `noisily'

	if "`vtlisth'"~="" {
		mata: `t' = tokens("`vtlisth'")
		mata: `sname'.vtlisth	= `t'
		parse_estring, sname(`sname') estring(`estringh') h `noisily'
	}
	
	mata: st_local("numlearners",strofreal(cols(`sname'.elist)))
	return scalar numlearners = `numlearners'

end

program define parse_estring, rclass

	syntax [anything] [if] [in] ,					/// 
							[						///
							sname(name)				/// name of mata struct
							estring(string asis)	/// names of estimation strings
													/// need asis option in case it includes strings
							h						/// indicates LIE eqn
							NOIsily					///
							]
	
	if "`noisily'"=="" {
		local qui quietly
	}

	// set struct fields
	if "`h'"=="" {
		local elist		elist
		local emlist	emlist
		local eolist	eolist
	}
	else {
		local elist		elisth
		local emlist	emlisth
		local eolist	eolisth
	}
	
	tempname t
	
	local doparse = 1
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
		
		mata: `sname'.`elist'	= (`sname'.`elist', `t')
		syntax [anything] , [*]
		local est_main `anything'
		local est_options `options'
		mata: `t' = "`est_main'"
		mata: `sname'.`emlist'	= (`sname'.`emlist', `t')
		mata: `t' = "`est_options'"
		mata: `sname'.`eolist'	= (`sname'.`eolist', `t')
		

		if "`2'"~="|" & "`3'"~="|" {
			// done parsing
			local doparse = 0
		}

		mac shift 3
		
		local estring `*'
	}

end

program define crossfit, rclass sortpreserve

	syntax [anything] [if] [in] ,					/// 
							[						///
							foldvar(varlist)		/// one fold var per resample
							kfolds(integer 5)		/// ignored if foldvars provided
							reps(integer 1)			/// ignored if foldvars provided
							NORANDOM				/// first fold ID uses obs in existing order
							resid					///
							vtilde(namelist)		/// name(s) of fitted variable(s)
							vname(varname)			/// name of original variable
							estring(string asis)	/// estimation string
													/// need asis option in case it includes strings
							vtype(string)			/// datatype of fitted variable; default=double
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							shortstack(name)		///
													/// 
													/// options specific to LIE/DDML-IV
							estringh(string asis)	/// est string for E[D^|XZ]		
													/// 
							NOIsily					///
							]

	// renaming for clarity
	local vtlist `vtilde'
	// clear the local macro
	local vtilde

	** indicator for LIE/optimal-IV model
	local lieflag	= "`estringh'"~=""
	** indicator for interactive model
	local tvflag	= "`treatvar'"~=""
	** indicator for short-stacking
	local ssflag	= "`shortstack'"~=""
 
	foreach vv in `vtlist' {
		if `lieflag' {
			local vtlisth `vtlisth' `vv'_h
		}
		if `tvflag' {
			local vtlist0 `vtlist0' `vv'0
			local vtlist1 `vtlist1' `vv'1
		}
	}

	// LIE => we want predicted values not resids
	if `lieflag' & "`resid'"~="" {
		di as res "resid option ignored"
		local resid
	}

	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `vname'
	
	if "`noisily'"=="" {
		local qui quietly
	}

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

	*** set up struct with equation/learner info
	
	tempname eqn_info
	initialize_eqn_info,													///
						sname(`eqn_info') reps(`reps')						///
						vtlist(`vtlist') estring(`estring')					///
						vtlist0(`vtlist0') vtlist1(`vtlist1')				///
						vtlisth(`vtlisth') estringh(`estringh')				///
						`noisily'
	local numlearners = r(numlearners)
	
	*** syntax
	if `ssflag' & `numlearners'==1 {
		di as err "error - shortstack option available only for multiple learners"
		exit 198
	}
	
	// datatype for fitted values/residuals
	if "`vtype'"=="" {
		local vtype double
	}
	
	// loop over resamples (crossfitting, shortstacking, store results)
	forvalues m=1/`reps' {

	******************************** CROSSFITTING ************************************
	
	// cross-fitting fundamentally depends on three cases 
	// case 1: `tvflag'==0 & `lieflag'==0
	// case 2: `tvflag'==1 & `lieflag'==0
	// case 3: `lieflag'==1
	// nb: `tvflag'==1 & `lieflag'==1 is impossible

		local fid : word `m' of `fidlist'
		// each resampling is stored on an rStruct
		tempname results_list
		mata: `results_list' = init_rStruct(`kfolds')

		// create blank fitted variable(s)
		if ~`tvflag' & ~`lieflag' { // case 1
			if `ssflag' {
				cap drop `shortstack'_`m'
				qui gen `vtype' `shortstack'_`m'=.
			}
			forvalues i=1/`numlearners' {
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
			forvalues i=1/`numlearners' {
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
			forvalues i=1/`numlearners' {
				tempvar dhat`i'
				qui gen `vtype' `dhat`i''=.  // using learner i
			}
			// out-of-sample predicted values for E[D|X] 
			forvalues i=1/`numlearners' {
				tempvar hhat`i'
				qui gen `vtype' `hhat`i''=.  // using learner i in both steps
			}
			// predicted values for E[D|ZX] for each k & learner i
			forvalues k=1/`kfolds' {
				forvalues i=1/`numlearners' {
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
				forvalues i=1/`numlearners' {
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
		tempname pysw pysw0 pysw1 pysw_t pysw_h
	
		// crossfit
		di
		di as text "Cross-fitting fold " _c
		forvalues k = 1(1)`kfolds' {
	
			di as text "`k' " _c
	
			forvalues i=1/`numlearners' {
				mata: st_local("est_main",`eqn_info'.emlist[`i'])
				mata: st_local("est_options",`eqn_info'.eolist[`i'])
				mata: st_local("vtilde",`eqn_info'.vtlist[`i'])
				if `lieflag' {
					// LIE locals
					mata: st_local("est_main_h",`eqn_info'.emlisth[`i'])
					mata: st_local("est_options_h",`eqn_info'.eolisth[`i'])
					mata: st_local("hhat",`eqn_info'.vtlisth[`i'])
				}
				
				if ~`tvflag' & ~`lieflag' { // case 1
				
					tempvar vhat_k
					
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `touse', `est_options'
					local cmd `e(cmd)'
		
					// save pystacked weights
					if ("`cmd'"=="pystacked") {
						if (`k'==1) {
							mat `pysw' = e(weights)
						}
						else {
							mat `pysw_t' = e(weights)
							mat `pysw' = (`pysw',`pysw_t')
						}
					}
		
					// get fitted values and residuals for kth fold	
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse'
		
					// get predicted values
					qui replace `vhat`i'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres`i'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
		
				}
		
				else if `tvflag' & ~`lieflag' {	// case 2: interactive models
		
					// outcome equation so estimate separately
		
					// for treatvar = 1
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `treatvar' == 1 & `touse', `est_options'
					local cmd `e(cmd)'
		
					// save pystacked weights
					if ("`cmd'"=="pystacked") {
						if (`k'==1) {
							mat `pysw0' = e(weights)
						}
						else {
							mat `pysw_t' = e(weights)
							mat `pysw0' = (`pysw0',`pysw_t')
						}
					}
		
					// get fitted values for kth fold	
					tempvar vhat_k
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse'
					qui replace `vhat1`i'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres1`i'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
		
					// for treatvar = 0
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `treatvar' == 0 & `touse', `est_options'
		
					// save pystacked weights
					if ("`cmd'"=="pystacked") {
						if (`k'==1) {
							mat `pysw1' = e(weights)
						}
						else {
							mat `pysw_t' = e(weights)
							mat `pysw1' = (`pysw1',`pysw_t')
						}
					}
		
					// get fitted values for kth fold	
					tempvar vhat_k
					qui predict `vtype' `vhat_k' if `fid'==`k' & `touse'
					qui replace `vhat0`i'' = `vhat_k' if `fid'==`k' & `touse'
					qui replace `vres0`i'' = `vname' - `vhat_k' if `fid'==`k' & `touse'
	
				}
		
				else if `lieflag' { // case 3
		
					tempvar vhat_k // stores predicted values for E[D|ZX] temporarily
					tempvar vtil_k // stores predicted values for E[D^|X] temporarily
		
					// Step I: estimation of E[D|XZ]=D^
					// estimate excluding kth fold
					`qui' `est_main' if `fid'!=`k' & `touse', `est_options'
					local cmd `e(cmd)'
		
					// get pystacked weights
					if ("`cmd'"=="pystacked") {
						if (`k'==1) {
							mat `pysw' = e(weights)
						}
						else {
							mat `pysw_t' = e(weights)
							mat `pysw' = (`pysw',`pysw_t')
						}
					}
					
					// get fitted values (in and out of sample)
					qui predict `vtype' `vhat_k' if `touse'
		
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
						if (`k'==1) {
							mat `pysw_h' = e(weights)
						}
						else {
							mat `pysw_t' = e(weights)
							mat `pysw_h' = (`pysw_h',`pysw_t')
						}
					}
		
					// get fitted values  
					qui predict `vtype' `vtil_k' if `touse'
		
					// get *combined* out-of-sample predicted values
					qui replace `hhat`i'' = `vtil_k' if `fid'==`k' & `touse'
		
				}
				
				if `k'==1 & `m'==1 {
					local cmd_list `cmd_list' `cmd'
				}
			}
		}
	
		// last fold, insert new line
		di "...completed cross-fitting"
		

		******************************** SHORTSTACKING ************************************
	
		// shortstacking. we need to distinguish between 3 cases again. 
		if `ssflag' & `numlearners'>1 {
	
			if ~`tvflag' & ~`lieflag' { // case 1
				local vhats
				forvalues i=1/`numlearners' {
					local vhats `vhats' `vhat`i''
				}
				tempvar vss
				`qui' di as text "Stacking NNLS (additive model):"
				`qui' _ddml_nnls `vname' `vhats'
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
				forvalues i=1/`numlearners' {
					local vhats1 `vhats1' `vhat1`i''
				}
				tempvar vtemp
				`qui' di as text "Stacking NNLS (interactive model, treatvar=1):"
				`qui' _ddml_nnls `vname' `vhats1' if `treatvar'==1
				qui predict `vtype' `vtemp'
				qui replace `shortstack'1_`m'=`vtemp'
					
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'1_`m' = `vname' - `shortstack'1_`m' if `treatvar'==1
				}
				// treatvar == 0
				forvalues i=1/`numlearners' {
					local vhats0 `vhats0' `vhat0`i''
				}
				tempvar vtemp
				`qui' di as text "Stacking NNLS (interactive model, treatvar=0):"
				`qui' _ddml_nnls `vname' `vhats0'
				qui predict `vtype' `vtemp'
				qui replace `shortstack'0_`m'=`vtemp' 
		
				if "`resid'"~="" {
					// vtilde is the residual
					qui replace `shortstack'0_`m' = `vname' - `shortstack'0_`m' if `treatvar'==0
				}
			}
			else if `lieflag' {
				// apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
				local dhats
				forvalues i=1/`numlearners' {
					local dhats `dhats' `dhat`i''
				}
				`qui' di as text "Stacking NNLS (LIE, OOS E[D|XZ]):"
				`qui' _ddml_nnls `vname' `dhats' if `touse'
				tempvar vtemp
				qui predict `vtype' `vtemp' if `touse'
				qui replace `dhatSS'=`vtemp' 
	
				// apply short-stacking to in-sample predicted values of E[D|XZ] *for each k*
				forvalues k = 1(1)`kfolds' {
					local dhats_is
					forvalues i=1/`numlearners' {
						local dhats_is `dhats_is' `dhat_`i'_`k''
					}
					tempvar vtemp
					`qui' di as text "Stacking NNLS (LIE, in-sample E[D|XZ] fold `k':"
					`qui' _ddml_nnls `vname' `dhats_is' if `fid'!=`k' & `touse' 
					qui predict `vtype' `vtemp'
					qui replace `dhat_isSS_`k'' = `vtemp' if `fid'!=`k' & `touse'
					qui replace `dhat_oosSS' = `vtemp' if `fid'==`k' & `touse'
				}
	
				// need to cross-fit stacked in-sample predicted values against X
				forvalues k = 1(1)`kfolds' {
					forvalues i=1/`numlearners' {
	
						mata: st_local("est_main_h",`eqn_info'.emlisth[`i'])
						mata: st_local("est_options_h",`eqn_info'.eolisth[`i'])
	
						// replace {D}-placeholder in estimation string with variable name
						//list `fid' `dhat_isSS_`k''
						local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_isSS_`k''",1)
			
						// estimation	
						`qui' `est_main_h_k' if `fid'!=`k' & `touse', `est_options_h'
						local cmd_h `e(cmd)'
					
						// get fitted values  
						tempvar vtemp
						qui predict `vtype' `vtemp' if `touse'
			
						// get out-of-sample predicted values
						qui replace `hhatSS`i'' = `vtemp' if `fid'==`k' & `touse'
					}
				}
	
				// final stacking for E[D|X]
				forvalues k = 1(1)`kfolds' {
					local hhatSS_list
					forvalues i=1/`numlearners' {
						local hhatSS_list `hhatSS_list' `hhatSS`i''
					}
					`qui' di as text "Stacking NNLS (LIE, E[D|X]):"
					// sum `dhat_oosSS' `hhatSS_list'
					`qui' _ddml_nnls `dhat_oosSS' `hhatSS_list'
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
	
		if `ssflag' & `numlearners'>1 {
			// last fold, insert new line
			di "...completed short-stacking"
		}
	
		
		******************************** STORE RESULTS ************************************
		
		tempname mse_list N_list mse_folds_list N_folds_list
		tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
		tempname mse0_list N0_list mse0_folds_list N0_folds_list
		tempname mse1_list N1_list mse1_folds_list N1_folds_list
		
		// AA: if we loop over "1 2 3 SS" this code would also calculate MSE for short-stacking in case of LIE
		// MS: this meant some hacky recoding so i gave up; now in a standalone section below
		
		forvalues i=1/`numlearners' {
	
			// vtilde, mspe, etc.
			if ~`tvflag' & ~`lieflag' {
		
				mata: st_local("vtilde",`eqn_info'.vtlist[`i'])
				cap drop `vtilde'_`m'
				if "`resid'"=="" {
					// vtilde is predicted values
					qui gen `vtilde'_`m' = `vhat`i''
				}
				else {
					// vtilde is residuals
					qui gen `vtilde'_`m' = `vres`i''
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
				
				mata: `results_list'.N			= (`results_list'.N			\ `N')
				mata: `results_list'.N_folds	= (`results_list'.N_folds	\ st_matrix("`N_folds'"))
				mata: `results_list'.MSE		= (`results_list'.MSE		\ `mse')
				mata: `results_list'.MSE_folds	= (`results_list'.MSE_folds	\ st_matrix("`mse_folds'"))
				
				if "`cmd'"=="pystacked" {
					// number of learners not known in advance so dimension of matrix not initialised as with other fields
					if `i'==1 {
						mata: `results_list'.stack_weights	= `pysw'
					}
					else {
						mata: `results_list'.stack_weights	= (`results_list'.stack_weights	\ `pysw')
					}
				}
			
			}
			else if `tvflag' & ~`lieflag' {
			
				mata: st_local("vtilde0",`eqn_info'.vtlist0[`i'])
				mata: st_local("vtilde1",`eqn_info'.vtlist1[`i'])
				cap drop `vtilde0'_`m'
				cap drop `vtilde1'_`m'
				if "`resid'"=="" {
					// vtilde is predicted values
					qui gen `vtilde0'_`m' = `vhat0`i''
					qui gen `vtilde1'_`m' = `vhat1`i''
				}
				else {
					// vtilde is residuals
					qui gen `vtilde0'_`m' = `vres0`i''
					qui gen `vtilde1'_`m' = `vres1`i''
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
				
				mata: `results_list'.N0			= (`results_list'.N0			\ `N0')
				mata: `results_list'.N0_folds	= (`results_list'.N0_folds		\ st_matrix("`N0_folds'"))
				mata: `results_list'.MSE0		= (`results_list'.MSE0			\ `mse0')
				mata: `results_list'.MSE0_folds	= (`results_list'.MSE0_folds	\ st_matrix("`mse0_folds'"))
				
				mata: `results_list'.N1			= (`results_list'.N1			\ `N1')
				mata: `results_list'.N1_folds	= (`results_list'.N1_folds		\ st_matrix("`N1_folds'"))
				mata: `results_list'.MSE1		= (`results_list'.MSE1			\ `mse1')
				mata: `results_list'.MSE1_folds	= (`results_list'.MSE1_folds	\ st_matrix("`mse1_folds'"))
				
				if "`cmd'"=="pystacked" {
					// number of learners not known in advance so dimension of matrix not initialised as with other fields
					if `i'==1 {
						mata: `results_list'.stack_weights0	= `pysw0'
						mata: `results_list'.stack_weights1	= `pysw1'
					}
					else {
						mata: `results_list'.stack_weights0	= (`results_list'.stack_weights0 \ `pysw0')
						mata: `results_list'.stack_weights1	= (`results_list'.stack_weights1 \ `pysw1')
					}
				}
				
			}
			else if `lieflag' {
	
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
				
				mata: `results_list'.N_h		= (`results_list'.N_h			\ `N_h')
				mata: `results_list'.N_h_folds	= (`results_list'.N_h_folds		\ st_matrix("`N_h_folds'"))
				mata: `results_list'.MSE_h		= (`results_list'.MSE_h			\ `mse_h')
				mata: `results_list'.MSE_h_folds= (`results_list'.MSE_h_folds	\ st_matrix("`mse_h_folds'"))
				
				if "`cmd'"=="pystacked" {
					// number of learners not known in advance so dimension of matrix not initialised as with other fields
					if `i'==1 {
						mata: `results_list'.stack_weights_h	= `pysw_h'
					}
					else {
						mata: `results_list'.stack_weights_h	= (`results_list'.stack_weights_h \ `pysw_h')
					}
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
				
				mata: `results_list'.N			= (`results_list'.N			\ `N')
				mata: `results_list'.N_folds	= (`results_list'.N_folds	\ st_matrix("`N_folds'"))
				mata: `results_list'.MSE		= (`results_list'.MSE		\ `mse')
				mata: `results_list'.MSE_folds	= (`results_list'.MSE_folds	\ st_matrix("`mse_folds'"))
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
				
				mata: `results_list'.N0			= (`results_list'.N0			\ `N0')
				mata: `results_list'.N0_folds	= (`results_list'.N0_folds		\ st_matrix("`N0_folds'"))
				mata: `results_list'.MSE0		= (`results_list'.MSE0			\ `mse0')
				mata: `results_list'.MSE0_folds	= (`results_list'.MSE0_folds	\ st_matrix("`mse0_folds'"))
				
				mata: `results_list'.N1			= (`results_list'.N1			\ `N1')
				mata: `results_list'.N1_folds	= (`results_list'.N1_folds		\ st_matrix("`N1_folds'"))
				mata: `results_list'.MSE1		= (`results_list'.MSE1			\ `mse1')
				mata: `results_list'.MSE1_folds	= (`results_list'.MSE1_folds	\ st_matrix("`mse1_folds'"))
			}
			else if `lieflag' {	// case 3
	
				// calculate and return mspe and sample size
				tempvar hres dres hres_sq dres_sq
				// vtilde has fitted values
				qui gen double `dres_sq' = (`shortstack'_`m')^2 if `touse'
				qui gen double `hres_sq' = (`shortstack'_h_`m')^2 if `touse'
	
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
				
				mata: `results_list'.N_h		= (`results_list'.N_h			\ `N_h')
				mata: `results_list'.N_h_folds	= (`results_list'.N_h_folds		\ st_matrix("`N_h_folds'"))
				mata: `results_list'.MSE_h		= (`results_list'.MSE_h			\ `mse_h')
				mata: `results_list'.MSE_h_folds= (`results_list'.MSE_h_folds	\ st_matrix("`mse_h_folds'"))
			}
		}
	
		// add struct with results for this resampling to list of results
		mata: `eqn_info'.rlist	= (`eqn_info'.rlist, &`results_list')
	
	}		// end of resampling loop
	
	// check: loop through results
	tempname r
	mata: `r' = init_rStruct(`kfolds')
	di "CHECK: results saved in rStruct?"
	forvalues m=1/`reps' {
		di "rep=`m'"
		mata: `r' = (*(`eqn_info'.rlist[`m']))
		mata: `r'.MSE
		mata: `r'.MSE0
		mata: `r'.MSE1
		mata: `r'.MSE_h
	}
	
	
	******************************** RETURN RESULTS ************************************

	if ~`tvflag' {
		
		return mat mse_folds		= `mse_folds'
		return mat N_folds			= `N_folds'
		return scalar mse			= `mse'

		return mat mse_list			= `mse_list'
		return mat N_list			= `N_list'
		return mat mse_folds_list	= `mse_folds_list'
		return mat N_folds_list		= `N_folds_list'
	}
	else {
		return scalar mse0			= `mse0'
		return scalar N0			= `N0'
		return scalar mse1			= `mse1'
		return scalar N1			= `N1'
		return mat mse0_folds		= `mse0_folds'
		return mat mse1_folds		= `mse1_folds'
		return mat N0_folds			= `N0_folds'
		return mat N1_folds			= `N1_folds'
	
		return mat mse0_list		= `mse0_list'
		return mat N0_list			= `N0_list'
		return mat mse0_folds_list	= `mse0_folds_list'
		return mat N0_folds_list	= `N0_folds_list'
		
		return mat mse1_list		= `mse1_list'
		return mat N1_list			= `N1_list'
		return mat mse1_folds_list	= `mse1_folds_list'
		return mat N1_folds_list	= `N1_folds_list'
	
	}
	if `lieflag' {
		return scalar N_h			= `N_h'
	
		return mat mse_h_folds		= `mse_h_folds'
		return mat N_h_folds		= `N_h_folds'
		return scalar mse_h			= `mse_h'
		
		return mat mse_h_list		= `mse_h_list'
		return mat N_h_list			= `N_h_list'
		return mat mse_h_folds_list	= `mse_h_folds_list'
		return mat N_h_folds_list	= `N_h_folds_list'
	}

	return scalar N			= `N'
	return local cmd		`cmd'
	return local cmd_h		`cmd_h'
	if ("`cmd'"=="pystacked" & `lieflag'==0 & `tvflag'==0) return mat pysw = `pysw' // pystacked weights
	if ("`cmd'"=="pystacked" & `lieflag'==0 & `tvflag'==1) return mat pysw0 = `pysw0'
	if ("`cmd'"=="pystacked" & `lieflag'==0 & `tvflag'==1) return mat pysw1 = `pysw1'
	if ("`cmd'"=="pystacked" & `lieflag'==1) return mat pysw_h 		= `pysw_h'
 
	return local cmd_list	`cmd_list'
	
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