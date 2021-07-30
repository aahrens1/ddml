* notes
* vtype was used for predicting values by fold, now applies to vtilde as well
* eqntype replaced by resid option (default=fitted)
* default is additive-type crossfitting; treatvar option triggers interactive-type crossfitting

program define crossfit, rclass sortpreserve

	syntax [anything] [if] [in] ,					/// 
							[						///
							kfolds(integer 0)		/// if not supplied, calculate
							NOIsily					///
							foldvar(name)			/// must be numbered 1...K where K=#folds
							resid					///
							vtilde(name)			/// name of fitted variable
							vname(varname)			/// name of original variable
							eststring(string asis)	/// estimation string
													/// need asis option in case it includes strings
							vtype(string)			/// datatype of fitted variable; default=double
							treatvar(varname)		/// 1 or 0 RHS variable; relevant for interactive model only
													/// if omitted then default is additive model
							]

	marksample touse
	
	if "`noisily'"=="" {
		local qui quietly
	}

	*** setup
	
	// datatype for fitted values/residuals
	if "`vtype'"=="" {
		local vtype double
	}
	
	// blank fitted variable
	cap gen `vtype' `vtilde'=.
	
	// if kfolds=0 then find number of folds
	// note this requires that foldvar takes values 1..K
	if `kfolds'==0 {
		qui sum `foldvar', meanonly
		local kfolds = r(max)
	}

	// parse estimation string
	// compound double quotes in case it includes quoted strings
	local 0 `"`eststring'"'
	syntax [anything] , [*]
	local est_main `anything'
	local est_options `options'
	`qui' di "est_main: `est_main'"
	`qui' di "est_options: `est_options'"
	
	// crossfit
	di
	di as text "Cross-fitting fold " _c
	forvalues k = 1(1)`kfolds' {

		di as text "`k' " _c

		if "`treatvar'"=="" {
		
			tempvar vtilde_k

			// estimate excluding kth fold
			`qui' `est_main' if `foldvar'!=`k' & `touse', `est_options'
			
			// get fitted values and residuals for kth fold	
			qui predict `vtype' `vtilde_k' if `foldvar'==`k' & `touse'

			if "`resid'"=="" {
				// get predicted values
				qui replace `vtilde' = `vtilde_k' if `foldvar'==`k' & `touse'
			} 
			else {
				// get residuals
				qui replace `vtilde' = `vname' - `vtilde_k' if `foldvar'==`k' & `touse'
			}
		}
		
		else {		// interactive model

			// outcome equation so estimate separately

			// for treatvar = 1
			// estimate excluding kth fold
			`qui' `est_main' if `foldvar'!=`k' & `treatvar' == 1 & `touse', `est_options'
			// get fitted values for kth fold	
			tempvar vtilde_k
			qui predict `vtype' `vtilde_k' if `foldvar'==`k' & `treatvar' == 1 & `touse'
			qui replace `vtilde' = `vtilde_k' if `foldvar'==`k' & `treatvar' == 1 & `touse'

			// for treatvar = 0
			// estimate excluding kth fold
			`qui' `est_main' if `foldvar'!=`k' & `treatvar' == 0 & `touse', `est_options'
			// get fitted values for kth fold	
			tempvar vtilde_k
			qui predict `vtype' `vtilde_k' if `foldvar'==`k' & `treatvar' == 0 & `touse'
			qui replace `vtilde' = `vtilde_k' if `foldvar'==`k' & `treatvar' == 0 & `touse'
			
			// replace with residuals if requested
			if "`resid'"~="" {
				qui replace `vtilde' = `vname' - `vtilde' if `foldvar'==`k' & `touse'
			}
		}
	}
	
	// last fold, insert new line
	di "...complete"

	// calculate and return mspe and sample size
	tempvar vtilde_sq
	if "`resid'"~="" {
		// vtilde has residuals
		qui gen double `vtilde_sq' = `vtilde'^2 if `touse'
	}
	else {
		// vtilde has fitted values
		qui gen double `vtilde_sq' = (`vname' - `vtilde')^2 if `touse'
	}
	
	// mspe
	if "`treatvar'"=="" {
		// additive-type model
		qui sum `vtilde_sq' if `touse', meanonly
		return scalar mse	= r(mean)
		local N				= r(N)
	}
	else {
		// interactive-type model, return mse separately for treatvar =0 and =1
		qui sum `vtilde_sq' if `treatvar' == 0 & `touse', meanonly
		return scalar mse0	= r(mean)
		local N				= r(N)
		qui sum `vtilde_sq' if `treatvar' == 1 & `touse', meanonly
		return scalar mse1	= r(mean)
		local N				= `N' + r(N)
	}
	
	return scalar N			= `N'
 
 end
 