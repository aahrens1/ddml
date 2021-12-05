* notes
* vtype was used for predicting values by fold, now applies to vtilde as well
* eqntype replaced by resid option (default=fitted)
* default is additive-type crossfitting; treatvar option triggers interactive-type crossfitting

* should recode so that 0 and 1 vtilde names not needed
* when looping will vname always be the same? (dep var) treatvar as well?

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
	real matrix vtlist
	real matrix vtlist0
	real matrix vtlist1
	real matrix elist
	real matrix emlist
	real matrix eolist
	real matrix vtlisth
	real matrix elisth
	real matrix emlisth
	real matrix eolisth
}

struct eStruct init_eStruct()
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
	return(d)
}

end

program define initialize_eqn_info, rclass

	syntax [anything] [if] [in] ,					/// 
							[						///
							sname(name)				/// name of mata struct
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
	
	mata: `sname'			= init_eStruct()
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
	
	mata: st_local("numeqns",strofreal(cols(`sname'.elist)))
	return scalar numeqns = `numeqns'

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
		`qui' di "est: `0'"
		syntax [anything] , [*]
		local est_main `anything'
		local est_options `options'
		`qui' di "est_main: `est_main'"
		`qui' di "est_options: `est_options'"
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
							kfolds(integer 0)		/// if not supplied, calculate
							NOIsily					///
							foldvar(name)			/// must be numbered 1...K where K=#folds
							resid					///
							vtilde(namelist)		/// name(s) of fitted variable
							vtilde0(namelist)		/// name of fitted variable
							vtilde1(namelist)		/// name of fitted variable
							vname(varname)			/// name of original variable
							eststring(string asis)	/// estimation string
													/// need asis option in case it includes strings
							vtype(string)			/// datatype of fitted variable; default=double
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							shortstack(name)		///
							shortstackh(name)		///
													/// 
													/// options specific to LIE/DDML-IV
							vtildeh(namelist)		/// intended for E[D^|X] where D^=E[D|XZ]=vtilde()	
							eststringh(string asis)	/// est string for E[D^|XZ]		
							]

	// temporary renaming until options above renamed
	local vtlist `vtilde'
	local vtlisth `vtildeh'
	local vtlist0 `vtilde0'
	local vtlist1 `vtilde1'

	// LIE => we want predicted values not resids
	if "`vtlisth'"~="" {
		local resid
	}

	marksample touse
	
	if "`noisily'"=="" {
		local qui quietly
	}

	*** setup
	
	tempname eqn_info
	initialize_eqn_info, sname(`eqn_info') vtlist(`vtlist') estring(`eststring')	///
						vtlist0(`vtlist0') vtlist1(`vtlist1')						///
						vtlisth(`vtlisth') estringh(`eststringh')					///
						`noisily'
	local numeqns = r(numeqns)
	
	*** syntax
	if "`shortstack'"~="" & `numeqns'==1 {
		di as err "error - shortstack option available only for multiple learners"
		exit 198
	}

	** indicator for LIE/optimal-IV model
	if "`eststringh'"!="" local lie lie
	
	// datatype for fitted values/residuals
	if "`vtype'"=="" {
		local vtype double
	}

	// if kfolds=0 then find number of folds
	// note this requires that foldvar takes values 1..K
	if `kfolds'==0 {
		qui sum `foldvar', meanonly
		local kfolds = r(max)
	}
	
	// cross-fitting fundamentally depends on three cases 
	// case 1: "`treatvar'"=="" & "`lie'"==""
	// case 2: "`treatvar'"!="" & "`lie'"==""
	// case 3: "`lie'"!=""

	// create blank fitted variable(s)
	if "`treatvar'"=="" & "`lie'"=="" { // case 1
		forvalues i=1/`numeqns' {
			tempvar vhat`i' vres`i'
			qui gen `vtype' `vhat`i''=.
			qui gen `vtype' `vres`i''=.
		}
	}  
	else if "`treatvar'"!="" & "`lie'"=="" { // case 2
		forvalues i=1/`numeqns' {
			tempvar vhat0`i' vres0`i'
			tempvar vhat1`i' vres1`i'
			qui gen `vtype' `vhat0`i''=.
			qui gen `vtype' `vres0`i''=.	
			qui gen `vtype' `vhat1`i''=.
			qui gen `vtype' `vres1`i''=.	
		}
	}
	else if "`lie'"!="" { // case 3
		// out-of-sample predicted values for E[D|XZ] 
		forvalues i=1/`numeqns' {
			tempvar dhat`i'
			qui gen `vtype' `dhat`i''=.  // using learner i
		}
		// predicted values for E[D|ZX] for each k & learner i
		forvalues k=1/`kfolds' {
			forvalues i=1/`numeqns' {
				tempvar dhat_`i'_`k'
				qui gen `vtype' `dhat_`i'_`k''=.
			}
		}
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
		forvalues i=1/`numeqns' {
			tempvar dtil`i' dtilSS`i'
			qui gen `vtype' `dtil`i''=.  // using learner i in both steps
			qui gen `vtype' `dtilSS`i''=. // using short-stacking for E[D|XZ] & then using learner i for E[D|X]
		}
		// short-stacking predicted values
		tempvar dhatSS
		qui gen `vtype' `dhatSS'=. 
		tempvar dtilSS
		gen `vtype' `dtilSS'=.
	}

	// save pystacked weights
	tempname pysw pysw0 pysw1 pysw_t pysw_h

	// crossfit
	di
	di as text "Cross-fitting fold " _c
	forvalues k = 1(1)`kfolds' {

		di as text "`k' " _c

		forvalues i=1/`numeqns' {
			mata: st_local("est_main",`eqn_info'.emlist[`i'])
			mata: st_local("est_options",`eqn_info'.eolist[`i'])
			mata: st_local("vtilde",`eqn_info'.vtlist[`i'])
			if "`lie'"!="" {
				// LIE locals
				mata: st_local("est_main_h",`eqn_info'.emlisth[`i'])
				mata: st_local("est_options_h",`eqn_info'.eolisth[`i'])
				mata: st_local("hhat",`eqn_info'.vtlisth[`i'])
			}
			
			if "`treatvar'"=="" & "`lie'"=="" { // case 1
			
				tempvar vhat_k
				
				// estimate excluding kth fold
				`qui' `est_main' if `foldvar'!=`k' & `touse', `est_options'
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
				qui predict `vtype' `vhat_k' if `foldvar'==`k' & `touse'
	
				// get predicted values
				qui replace `vhat`i'' = `vhat_k' if `foldvar'==`k' & `touse'
				qui replace `vres`i'' = `vname' - `vhat_k' if `foldvar'==`k' & `touse'
	
			}
	
			else if "`treatvar'"!="" & "`lie'"=="" {		// case 2: interactive models
	
				// outcome equation so estimate separately
	
				// for treatvar = 1
				// estimate excluding kth fold
				`qui' `est_main' if `foldvar'!=`k' & `treatvar' == 1 & `touse', `est_options'
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
				qui predict `vtype' `vhat_k' if `foldvar'==`k' & `touse'
				qui replace `vhat1`i'' = `vhat_k' if `foldvar'==`k' & `touse'
				qui replace `vres1`i'' = `vname' - `vhat_k' if `foldvar'==`k' & `touse'
	
				// for treatvar = 0
				// estimate excluding kth fold
				`qui' `est_main' if `foldvar'!=`k' & `treatvar' == 0 & `touse', `est_options'
	
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
				qui predict `vtype' `vhat_k' if `foldvar'==`k' & `touse'
				qui replace `vhat0`i'' = `vhat_k' if `foldvar'==`k' & `touse'
				qui replace `vres0`i'' = `vname' - `vhat_k' if `foldvar'==`k' & `touse'

			}
	
			else if "`lie'"!="" { // case 3
	
				tempvar vhat_k // stores predicted values for E[D|ZX] temporarily
				tempvar vtil_k // stores predicted values for E[D^|X] temporarily
	
				// Step I: estimation of E[D|XZ]=D^
				// estimate excluding kth fold
				`qui' `est_main' if `foldvar'!=`k' & `touse', `est_options'
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
				
				// get fitted values (in and out of sample)
				qui predict `vtype' `vhat_k' if `touse'
	
				// get out-of-sample predicted values
				qui replace `dhat`i'' = `vhat_k' if `foldvar'==`k' & `touse'
	
				// save predicted values in wide format for each i and k
				qui replace `dhat_`i'_`k'' = `vhat_k' if `touse'
	
				// Step II: estimation of E[D^|X]
	
				// replace {D}-placeholder in estimation string with variable name
				local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_is`i'_ `k''",1)
	
				// estimation	
				`qui' `est_main_h_k' if `foldvar'!=`k' & `touse', `est_options_h'
				local cmd_h `e(cmd)'
	
				// save pystacked weights
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
	
				// get out-of-sample predicted values
				qui replace `dtil`i'' = `vtil_k' if `foldvar'==`k' & `touse'
	
			}
			
			// should perhaps do this elsewhere
			if `k'==1 {
				local cmd_list `cmd_list' `cmd'
			}
		}
	}

	// last fold, insert new line
	di "...completed cross-fitting"

	// shortstacking. we need to distinguish between 3 cases again. 
	if "`shortstack'"~="" & `numeqns'>1 {

		if "`treatvar'"=="" & "`lie'"=="" { // case 1
			local vhats
			forvalues i=1/`numeqns' {
				local vhats `vhats' `vhat`i''
			}
			tempvar vss
			`qui' di as text "Stacking NNLS:"
			`qui' _ddml_nnls `vname' `vhats'
			tempvar vtemp
			predict `vtype' `vtemp'
			replace `shortstack' = `vtemp'
				
			if "`resid'"~="" {
				// vtilde is the residual
				qui replace `shortstack' = `vname' - `shortstack'
			}
		}
		else if "`treatvar'"!="" & "`lie'"=="" {		// case 2: interactive models
			local vhats1 
			local vhats0
			// treatvar == 1
			forvalues i=1/`numeqns' {
				local vhats1 `vhats' `vhat1`i''
			}
			tempvar vss
			`qui' di as text "Stacking NNLS:"
			`qui' _ddml_nnls `vname' `vhats1' if `treatvar'==1
			predict `vtype' `shortstack'1 
				
			if "`resid'"~="" {
				// vtilde is the residual
				qui replace `shortstack'1 = `vname' - `shortstack'1 if treatvar == 1
			}
			// treatvar == 0
			forvalues i=1/`numeqns' {
				local vhats0 `vhats0' `vhat0`i''
			}
			tempvar vss
			`qui' di as text "Stacking NNLS:"
			`qui' _ddml_nnls `vname' `vhats0'
			predict `vtype' `shortstack'1 
	
			if "`resid'"~="" {
				// vtilde is the residual
				qui replace `shortstack'0 = `vname' - `shortstack'0
			}
		}
		else if "`lie'"!="" {
			// apply short-stacking to cross-fitted (out-of-sample) predicted values of E[D|XZ]
			local dhats
			forvalues i=1/`numeqns' {
				local dhats `dhats' `dhat`i''
			}
			`qui' di as text "Stacking NNLS:"
			`qui' _ddml_nnls `vname' `dhats'
			tempvar vtemp
			predict `vtype' `vtemp'
			replace `dhatSS'=`vtemp' 

			// apply short-stacking to in-sample predicted values of E[D|XZ] *for each k*
			forvalues k = 1(1)`kfolds' {
				local dhats_is
				forvalues i=1/`numeqns' {
					local dhats_is `dhats_is' `dhat_`i'_`k''
				}
				tempvar vtemp
				`qui' di as text "Stacking NNLS:"
				`qui' _ddml_nnls `vname' `dhats_is' if `foldvar'!=`k' & `touse' 
				predict `vtype' `vtemp'
				sum `vtemp'
				replace `dhat_isSS_`k'' = `vtemp' if `foldvar'!=`k'
				replace `dhat_oosSS' = `vtemp' if `foldvar'==`k'
			}

			// need to cross-fit stacked in-sample predicted values against X
			forvalues k = 1(1)`kfolds' {
				forvalues i=1/`numeqns' {

					mata: st_local("est_main_h",`eqn_info'.emlisth[`i'])
					mata: st_local("est_options_h",`eqn_info'.eolisth[`i'])
					mata: st_local("hhat",`eqn_info'.vtlisth[`i'])

					// replace {D}-placeholder in estimation string with variable name
					//list `foldvar' `dhat_isSS_`k''
					local est_main_h_k = subinstr("`est_main_h'","{D}","`dhat_isSS_`k''",1)
		
					// estimation	
					`qui' `est_main_h_k' if `foldvar'!=`k' & `touse', `est_options_h'
					local cmd_h `e(cmd)'
				
					// get fitted values  
					tempvar vtemp
					qui predict `vtype' `vtemp' if `touse'
		
					// get out-of-sample predicted values
					qui replace `dtilSS`i'' = `vtemp' if `foldvar'==`k' & `touse'
				}
			}

			// final stacking for E[D|X]
			forvalues k = 1(1)`kfolds' {
				local dtilSS_list
				forvalues i=1/`numeqns' {
					local dtilSS_list `dtilSS_list' `dtilSS`i''
				}
				`qui' di as text "Stacking NNLS:"
				sum `dhat_oosSS' `dtilSS_list'
				`qui' _ddml_nnls `dhat_oosSS' `dtilSS_list'
				tempvar vtemp
				predict `vtype' `vtemp'
				replace `dtilSS'=`vtemp'
			}
			gen `vtype' `shortstack'=`dhatSS'
			label var `shortstack' "short-stacking cross-fitted E[D|Z,X]"
			gen `vtype' `shortstackh'=`dtilSS'
			label var `shortstackh' "short-stacking cross-fitted E[D|X]"
		}
	}

	// last fold, insert new line
	di "...completed short-stacking"

	tempname mse_list N_list mse_folds_list N_folds_list
	tempname mse_h_list N_h_list mse_h_folds_list N_h_folds_list
	tempname mse0_list N0_list mse0_folds_list N0_folds_list
	tempname mse1_list N1_list mse1_folds_list N1_folds_list
	
	forvalues i=1/`numeqns' {

		// maybe move below into LIE block ... unless this is supposed to be used in the first block below?
		if "`lie'"!="" {
			// vtilde has fitted values
			tempvar hres_sq
			qui gen double `hres_sq' = (`vhat`i'' - `hhat`i'')^2 if `touse'		
		}
	
		// vtilde, mspe, etc.
		if "`treatvar'"=="" {
	
			mata: st_local("vtilde",`eqn_info'.vtlist[`i'])
			if "`resid'"=="" {
				// vtilde is predicted values
				qui gen `vtilde' = `vhat`i''
			}
			else {
				// vtilde is residuals
				qui gen `vtilde' = `vres`i''
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
				qui sum `vres_sq' if `touse' & `foldvar'==`k', meanonly
				mat `mse_folds' = (nullmat(`mse_folds'), r(mean))
				qui count if `touse' & `foldvar'==`k' & `vres_sq'<.
				mat `N_folds' = (nullmat(`N_folds'), r(N))
			}
		
			mat `mse_list'			= (nullmat(`mse_list') \ `mse')
			mat `N_list'			= (nullmat(`N_list') \ `N')
			mat `mse_folds_list'	= (nullmat(`mse_folds_list') \ `mse_folds')
			mat `N_folds_list'		= (nullmat(`N_folds_list')\ `N_folds')
		}
		else {
		
			mata: st_local("vtilde0",`eqn_info'.vtlist0[`i'])
			mata: st_local("vtilde1",`eqn_info'.vtlist1[`i'])
			if "`resid'"=="" {
				// vtilde is predicted values
				qui gen `vtilde0' = `vhat0`i''
				qui gen `vtilde1' = `vhat1`i''
			}
			else {
				// vtilde is residuals
				qui gen `vtilde0' = `vres0`i''
				qui gen `vtilde1' = `vres1`i''
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
				qui sum `vres0_sq' if `treatvar' == 0 & `touse' & `foldvar'==`k', meanonly
				mat `mse0_folds' = (nullmat(`mse0_folds'), r(mean))
				qui sum `vres1_sq' if `treatvar' == 1 & `touse' & `foldvar'==`k', meanonly
				mat `mse1_folds' = (nullmat(`mse1_folds'), r(mean))
				qui count if `treatvar' == 0 & `touse' & `foldvar'==`k' & `vres1_sq'<.
				mat `N0_folds' = (nullmat(`N0_folds'), r(N))
				qui count if `treatvar' == 1 & `touse' & `foldvar'==`k' & `vres0_sq'<.
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
			
		}
		
		if "`lie'"!="" {
			qui sum `hres_sq' if `touse', meanonly
			local mse_h			= r(mean)
			local N_h			= r(N)	
			tempname mse_h_folds N_h_folds
			forvalues k = 1(1)`kfolds' {
				qui sum `hres_sq' if `touse' & `foldvar'==`k', meanonly
				mat `mse_h_folds' = (nullmat(`mse_h_folds'), r(mean))
				qui count if `touse' & `foldvar'==`k' & `hres_sq'<.
				mat `N_h_folds' = (nullmat(`N_h_folds'), r(N))
			}
			mat `mse_h_list'		= (nullmat(`mse_h_list') \ `mse_h')
			mat `N_h_list'			= (nullmat(`N_h_list') \ `N_h')
			mat `mse_h_folds_list'	= (nullmat(`mse_h_folds_list') \ `mse_h_folds')
			mat `N_h_folds_list'	= (nullmat(`N_h_folds_list')\ `N_h_folds')
		}
	}
	
	if "`treatvar'"=="" {
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
	if "`lie'"~="" {
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
	if ("`cmd'"=="pystacked" & "`lie'"=="" & "`treatvar'"=="") return mat pysw = `pysw' // pystacked weights
	if ("`cmd'"=="pystacked" & "`lie'"=="" & "`treatvar'"!="") return mat pysw0 = `pysw0'
	if ("`cmd'"=="pystacked" & "`lie'"=="" & "`treatvar'"!="") return mat pysw1 = `pysw1'
	if ("`cmd'"=="pystacked" & "`lie'"!="") return mat pysw_h 		= `pysw_h'
 
	return local cmd_list	`cmd_list'
	
 end
 