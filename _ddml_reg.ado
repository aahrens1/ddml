// does OLS/IV and reports with substitute yname and dnames
program define _ddml_reg, eclass
	syntax [anything] [if] [in] , [ y(name) yname(name) d(namelist) dnames(namelist) z(namelist) znames(namelist) mname(name) rep(integer 1) robust * ]

	if ~replay() {
		
		marksample touse
	
		// add resample suffixes and estimate
		local y_m `y'_`rep'
		add_suffix `d', suffix("_`rep'")
		local d_m `s(vnames)'
		add_suffix `z', suffix("_`rep'")
		local z_m `s(vnames)'
		
		if "`z_m'"=="" {
			qui reg `y_m' `d_m'         if `touse', `robust' `options'
		}
		else {
			// old-style regress syntax: put IVs in parentheses
			qui reg `y_m' `d_m' (`z_m') if `touse', `robust' `options'
		}
		tempname b
		tempname V
		mat `b' = e(b)
		mat `V' = e(V)
		matrix colnames `b' = `dnames'
		matrix rownames `b' = `yname'
	 	matrix colnames `V' = `dnames'
		matrix rownames `V' = `dnames'
		local N = e(N)
		ereturn clear
		ereturn post `b' `V', depname(`yname') obs(`N') esample(`touse')
		if "`robust'"=="robust" {
			ereturn local vce		robust
			ereturn local vctype	Robust
		}
		ereturn local cmd _ddml_reg
		ereturn local y `y_m'
		ereturn local d `d_m'
		ereturn local z `z_m'
		
		// additional estimation results
		ereturn scalar resample = `rep'
		tempname eqn
		mata: `eqn' = init_eStruct()
		// Y eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`yname'")
		mata: st_numscalar("e(`y'_mse)",return_result_item(`eqn',"`y'","MSE","`rep'"))
		mata: st_matrix("e(`y'_mse_folds)",return_result_item(`eqn',"`y'","MSE_folds","`rep'"))
		noi mata: return_learner_item(`eqn',"`y'","cmd")
		mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`y'","cmd")=="pystacked"))
		if `pyswflag' {
			mata: st_matrix("e(`y'_pysw)", mean(return_result_item(`eqn',"`y'","stack_weights","`rep'")'))
		}
		// D eqn results
		local numeqnD	: word count `dnames'
		forvalues i=1/`numeqnD' {
			local dname : word `i' of `dnames'
			local vtilde : word `i' of `d'
			mata: `eqn' = (`mname'.eqnAA).get("`dname'")
			mata: st_numscalar("e(`vtilde'_mse)",return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
			mata: st_matrix("e(`vtilde'_mse_folds)",return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
			mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`vtilde'","cmd")=="pystacked"))
			if `pyswflag' {
				mata: st_matrix("e(`vtilde'_pysw)", mean(return_result_item(`eqn',"`vtilde'","stack_weights","`rep'")'))
			}
		}
	
		// no longer needed
		cap mata: mata drop `eqn'
	}
	
	di as text "E[y|X] = " as res "`e(y)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	di as text "E[D|X] = " as res "`e(d)'"
	if "`e(z)'" ~= "" {
		di as text "E[Z|X] = " as res "`e(z)'"
	}
	ereturn display

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
