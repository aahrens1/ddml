program define stacking, eclass

	syntax varname [if] [in] ,						/// 
							[						///
							foldvar(varlist)		/// one fold var per resample
							kfolds(integer 5)		/// ignored if foldvars provided
							estring(string asis)	/// estimation string
							vtilde(namelist)		/// name(s) of fitted variable(s)
							stack(name)				///
							eqn(name)				///
							NOIsily					///
							]

	// renaming for clarity
	local vtlist `vtilde'
	// clear the local macro
	local vtilde
	
	if "`noisily'"=="" {
		local qui quietly
	}
	
	marksample touse
	// if dep var is missing, automatically not in estimation sample
	markout `touse' `vname'

	// blank eqn - declare this way so that it's a struct and not transmorphic
	if "`eqn'" == "" {
		local eqn stacking_eqn
	}
	mata: `eqn' = init_eStruct()
	
	tempname mname
	ddml init partial, mname(`mname')
	_ddml_sample if `touse' , mname(`mname')
	// sample, id, fid not needed
	cap drop `mname'*
	
	local doparse = 1
	local vnum = 1
	local hasvtlist = "`vtlist'"~=""
	while `doparse' {
		tokenize `"`estring'"', parse("||")
		
		// catch special case - a single | appears inside the estimation string
		if "`2'"=="|" & "`3'"~="|" {
			local firstpart `1' `2' `3'
			mac shift 2
			local estring `*'
			tokenize `"`estring'"', parse("||")
			local 1 `firstpart'
		}
		
		if `hasvtlist' {
			local vtilde : word `vnum' of `vtlist'
		}
		else {
			local vtilde : word 1 of `1'
			local vtilde Y`vnum'_`vtilde'
			local vtlist `vtlist' `vtilde'
		}
		
		// put in a yeq
		`qui' ddml yeq, learner(`vtilde') mname(`mname') vname(`varlist') : `1'
	
		if "`2'"~="|" & "`3'"~="|" {
			// done parsing
			local doparse = 0
		}

		mac shift 3
		
		local estring `*'
		
		local ++vnum
	}
	
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: `eqn'.shortstack = "`stack'"
	// model struct no longer needed
	mata: mata drop `mname'
	
	`qui' crossfit if `touse',		///
		ename(`eqn') noreplace		///
		`noisily'
	
	mata: st_local("nlearners", strofreal(`eqn'.nlearners))
	
	// stacking weights
	tempname ss_weights
	mata: st_matrix("`ss_weights'",return_result_item(`eqn',"`stack'","ss_weights", "1"))

	foreach var of varlist `vtilde' {
		qui replace `touse' = 0 if `var'==.
	}
	qui count if `touse'
	local N = r(N)
	
	// housekeeping
	// drop created predicted values
	foreach vtilde in `vtlist' {
		cap drop `vtilde'_1
	}
	cap drop `stack'_1
	
	// display results
	di
	di as text "{hline 17}{c TT}{hline 21}"
	di as text "  Method" _c
	di as text _col(18) "{c |}      Weight"
	di as text "{hline 17}{c +}{hline 21}"
	forvalues j=1/`nlearners' {
		local b : word `j' of `vtlist'
		di as text "  `b'" _c
		di as text _col(18) "{c |}" _c
		di as res %15.7f el(`ss_weights',1,`j')
	}
		
	ereturn clear
	ereturn post, depname(`varlist') obs(`N') esample(`touse')
	ereturn local cmd stacking
	ereturn local stack `stack'
	ereturn local predict stacking_p
	ereturn local eqn `eqn'
	ereturn local vtlist `vtlist'
	ereturn scalar nlearners = `nlearners'

end
