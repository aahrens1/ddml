// does OLS/IV and reports with substitute yname and dnames
program define _ddml_reg, eclass
	syntax [anything] [if] [in] , [								///
				y(name) yname(name)								///
				d(namelist) dnames(namelist) dvtnames(namelist)	///
				z(namelist) znames(namelist) zvtnames(namelist)	///
				mname(name) rep(integer 1)						///
				robust											///
				NOREP											///
				*												///
				]

	
	if ~replay() {
		
		marksample touse
		
		mata: st_local("model",`mname'.model)
		local ivhdflag	= "`model'"=="ivhd"
				
		// default vtilde names
		if "`dvtnames'"=="" {
			local dvtnames `d'
		}
		if "`zvtnames'"=="" {
			local zvtnames `z'
		}

		// add resample suffixes and estimate
		// y always gets a suffix
		local y_m `y'_`rep'
		if "`norep'"=="" {	
			add_suffix `d', suffix("_`rep'")
			local d_m `s(vnames)'
			add_suffix `z', suffix("_`rep'")
			local z_m `s(vnames)'
		}
		else {
			local d_m `d'
			local z_m `z'
		}
		
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
		ereturn local model `model'
		ereturn local y `y_m'
		if `ivhdflag'==0 {
			ereturn local d `d_m'
			ereturn local z `z_m'
		}
		else {
			add_suffix `dvtnames', suffix("_`rep'")
			ereturn local d `s(vnames)'
			add_suffix `zvtnames', suffix("_`rep'")
			ereturn local dh `s(vnames)'
		}
		// additional estimation results
		ereturn scalar resample = `rep'
		tempname eqn
		mata: `eqn' = init_eStruct()
		// Y eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`yname'")
		mata: st_numscalar("e(`y'_mse)",return_result_item(`eqn',"`y'","MSE","`rep'"))
		mata: st_matrix("e(`y'_mse_folds)",return_result_item(`eqn',"`y'","MSE_folds","`rep'"))
		mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`y'","cmd")=="pystacked"))
		if `pyswflag' {
			mata: st_matrix("e(`y'_pysw)", mean(return_result_item(`eqn',"`y'","stack_weights","`rep'")'))
		}
		// D eqn results - uses vtnames
		local numeqnD	: word count `dnames'
		forvalues i=1/`numeqnD' {
			local dname : word `i' of `dnames'
			local vtilde : word `i' of `dvtnames'
			local vtilde_h : word `i' of `zvtnames'
			// remove the trailing "_h" so that the AA lookup uses the learner name
			local vtilde_h = substr("`vtilde_h'",1,strlen("`vtilde_h'")-2)
			mata: `eqn' = (`mname'.eqnAA).get("`dname'")
			mata: st_numscalar("e(`vtilde'_mse)",return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
			mata: st_matrix("e(`vtilde'_mse_folds)",return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
			mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`vtilde'","cmd")=="pystacked"))
			if `pyswflag' {
				mata: st_matrix("e(`vtilde'_pysw)", mean(return_result_item(`eqn',"`vtilde'","stack_weights","`rep'")'))
			}
			if `ivhdflag' {
				mata: st_numscalar("e(`vtilde_h'_mse_h)",return_result_item(`eqn',"`vtilde_h'","MSE_h","`rep'"))
				mata: st_matrix("e(`vtilde_h'_mse_h_folds)",return_result_item(`eqn',"`vtilde_h'","MSE_h_folds","`rep'"))
			}
		}
		if `ivhdflag'==0 {
			// Z eqn results; ivhd won't enter
			local numeqnZ	: word count `znames'
			forvalues i=1/`numeqnZ' {
				local zname : word `i' of `znames'
				local vtilde : word `i' of `z'
				mata: `eqn' = (`mname'.eqnAA).get("`zname'")
				mata: st_numscalar("e(`vtilde'_mse)",return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
				mata: st_matrix("e(`vtilde'_mse_folds)",return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
				mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`vtilde'","cmd")=="pystacked"))
				if `pyswflag' {
					mata: st_matrix("e(`vtilde'_pysw)", mean(return_result_item(`eqn',"`vtilde'","stack_weights","`rep'")'))
				}
			}
		}
		// no longer needed
		cap mata: mata drop `eqn'
	}
	
	di as text "y-E[y|X]" _col(11) "= " as res "`e(y)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	if "`e(model)'"~="ivhd" {
		di as text "D-" _c
	}
	di as text "E[D|X]" _col(11)  "= " as res "`e(d)'"
	if "`e(model)'" == "iv" {
		di as text "Z-E[Z|X]" _col(11) "= " as res "`e(z)'"
	}
	else if "`e(model)'" == "ivhd" {
		di as text "E[D^|X,Z] = " as res "`e(dh)'"
	}
	if "`e(model)'" == "ivhd" {
		di as text "Orthogonalised D = D - E[D|X]; optimal IV = E[D^|X,Z] - E[D|X]."
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
