// does OLS/IV and reports with substitute yname and dnames
program define _ddml_reg, eclass
	syntax [anything] [if] [in] , [								///
				y(name) yname(name)								///
				d(namelist) dnames(namelist) dvtnames(namelist)	///
				z(namelist) znames(namelist) zvtnames(namelist)	///
				mname(name)										///
				spec(string) rep(string)						///
				vce(string)										///
				title(string)									///
				medmean(string)									///
				NOREP											///
				replay											///
				*												/// can be e.g. nocons
				]

	mata: st_local("model",`mname'.model)
	local ivhdflag	= "`model'"=="ivhd"
	
	if "`replay'"=="" & "`medmean'"=="" {	// estimate from scratch
		
		marksample touse
		
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		// default vtilde names
		if "`dvtnames'"=="" {
			local dvtnames `d'
		}
		if "`zvtnames'"=="" {
			local zvtnames `z'
		}

		// add resample suffixes
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
		
		// estimate
		if "`z_m'"=="" {
			qui reg `y_m' `d_m'         if `touse', vce(`vce') `options'
		}
		else {
			// old-style regress syntax: put IVs in parentheses
			qui reg `y_m' `d_m' (`z_m') if `touse', vce(`vce') `options'
		}
		tempname b V
		mat `b' = e(b)
		mat `V' = e(V)
		local N = e(N)

		local vce		`e(vce)'
		local vcetype	`e(vcetype)'
		local clustvar	`e(clustvar)'
		local N_clust	=e(N_clust)
		if `ivhdflag' {
			local d		`dvtnames'
			add_suffix	`dvtnames', suffix("_`rep'")
			local d_m	`s(vnames)'
			local dh	`zvtnames'
			add_suffix	`zvtnames', suffix("_`rep'")
			local dh_m	`s(vnames)'
			// clear locals (IV info is stored as dh not z)
			local z
			local z_m
		}
		
		// store post objects
		mata: `A'.put(("N","post"),`N')
		mata: `A'.put(("b","post"),st_matrix("`b'"))
		mata: `A'.put(("V","post"),st_matrix("`V'"))
		mata: `A'.put(("depvar","post"),"`yname'")
		
		// for calling program
		ereturn clear
		mata: st_matrix("e(bmat)",st_matrix("`b'"))
		mata: st_matrix("e(semat)",sqrt(diagonal(st_matrix("`V'"))'))
		
		// store locals
		local list_local title y y_m d d_m dh dh_m z z_m yname dnames vce vcetype
		if "`clustvar'"~=""		local list_local `list_local' clustvar
		foreach obj in `list_local' {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		// store scalars
		local list_scalar
		if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
		foreach obj in `list_scalar' {
			mata: `A'.put(("`obj'","scalar"),``obj'')
		}
		
		// additional estimation results
		ereturn scalar resample = `rep'
		tempname eqn
		mata: `eqn' = init_eStruct()
		// Y eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`yname'")
		// MSE
		mata: `A'.put(("`y'_mse","scalar"),return_result_item(`eqn',"`y'","MSE","`rep'"))
		// MSE folds
		mata: `A'.put(("`y'_mse_folds","matrix"),return_result_item(`eqn',"`y'","MSE_folds","`rep'"))
		// ss weights
		mata: st_local("shortstack_vname", `eqn'.shortstack)
		if "`shortstack_vname'"!="" {
			mata: `A'.put(("`y'_ssw","matrix"), return_result_item(`eqn',"`shortstack_vname'","ss_weights","`rep'"))
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
			// MSE
			mata: `A'.put(("`vtilde'_mse","scalar"),return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
			// MSE folds
			mata: `A'.put(("`vtilde'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
			// ss weights
			mata: st_local("shortstack_vname", `eqn'.shortstack)
			if "`shortstack_vname'"!="" {
				mata: `A'.put(("`d'_ssw","matrix"), return_result_item(`eqn',"`shortstack_vname'","ss_weights","`rep'"))
			}
			if `ivhdflag' {
				// MSE
				mata: `A'.put(("`vtilde_h'_mse","scalar"),return_result_item(`eqn',"`vtilde_h'","MSE_h","`rep'"))
				// MSE folds
				mata: `A'.put(("`vtilde_h'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde_h'","MSE_h_folds","`rep'"))
				// ss weights
				mata: st_local("shortstack_vname", `eqn'.shortstack)
				if "`shortstack_vname'"!="" {
						mata: `A'.put(("`y'_ssw","matrix"), return_result_item(`eqn',"`shortstack_vname'","ss_weights_h","`rep'"))
				}
			}
		}
		if `ivhdflag'==0 {
			// Z eqn results; ivhd won't enter
			local numeqnZ	: word count `znames'
			forvalues i=1/`numeqnZ' {
				local zname : word `i' of `znames'
				local vtilde : word `i' of `z'
				mata: `eqn' = (`mname'.eqnAA).get("`zname'")
				// MSE
				mata: `A'.put(("`vtilde'_mse","scalar"),return_result_item(`eqn',"`vtilde'","MSE","`rep'"))
				// MSE folds
				mata: `A'.put(("`vtilde'_mse_folds","matrix"),return_result_item(`eqn',"`vtilde'","MSE_folds","`rep'"))
				// ss weights
				mata: st_local("shortstack_vname", `eqn'.shortstack)
				if "`shortstack_vname'"!="" {
					mata: `A'.put(("`z'_ssw","matrix"), return_result_item(`eqn',"`shortstack_vname'","ss_weights","`rep'"))
				}
			}
		}
		
		mata: (`mname'.estAA).put(("`spec'","`rep'"),`A')
		
		// no longer needed
		foreach obj in `A' `eqn' {
			cap mata: mata drop `obj'
		}
		
	}
	else if "`replay'"=="" & "`medmean'"~="" {	// aggregate over resamples
		
		tempname b V bagg Vagg Vi
		tempname bvec sbvec bmed Vvec sVvec Vmed
		tempvar esample
		tempname B
		
		// initialize
		mata: st_local("nameD",invtokens(`mname'.nameD))
		local K : word count `nameD'
		mata: st_local("nreps",strofreal(`mname'.nreps))
		mata: `B' = AssociativeArray()
		local isodd = mod(`nreps',2)
		local medrow = ceil(`nreps'/2)
		local N = 0
		
		// bvec a misnomer - usually a vector, but can be a matrix if multiple D variables
		mata: `bvec' = J(`nreps',`K',0)
		mata: `bagg' = J(1,`K',0)
		forvalues m=1/`nreps' {
			mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
			mata: `bvec'[`m',.] = `B'.get(("b","post"))
			// row/colnames etc. - need to do this only once
			if `m'==1 {
				mata: st_local("depvar",`B'.get(("depvar","post")))
				// retrieve locals; if empty, will be ""
				local list_local y d dh z yname dnames vce vcetype clustvar
				foreach obj in `list_local' {
					mata: st_local("`obj'",`B'.get(("`obj'","local")))
				}
				// retrieve scalars (as locals)
				local list_scalar
				if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
				foreach obj in `list_scalar' {
					mata: st_local("`obj'",strofreal(`B'.get(("`obj'","scalar"))))
				}
			}
			// possible that different estimation samples have different #obs
			qui count if `mname'_sample_`m'==1
			local N = `N' + r(N)
		}
		local N = round(`N'/`nreps')
		
		if "`medmean'"=="mn" {
			// mean beta
			mata: `bagg' = mean(`bvec')
			mata: st_matrix("`bagg'",`bagg')
		}
		else if "`medmean'"=="md" {
			// median beta
			forvalues k=1/`K' {
				mata: `sbvec' = sort(`bvec',`k')
				if `isodd' {
					mata: `bagg'[1,`k'] = `sbvec'[`medrow',`k']
				}
				else {
					mata: `bagg'[1,`k'] = (`sbvec'[`medrow',`k'] + `sbvec'[`medrow'+1,`k'])/2
				}
			}
			mata: st_matrix("`bagg'",`bagg')
		}
		else {
			di as err "_ddml_reg error - unrecognized option `medmean'"
			exit 198
		}
		
		mata: `Vagg' = J(`K',`K',0)
		mata: `Vvec' = J(`nreps',1,0)
		if "`medmean'"=="mn" {
			// harmonic mean
			// inefficient - does off-diagonals twice
			forvalues m=1/`nreps' {
				mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
				mata: `Vi' = `B'.get(("V","post"))
				forvalues j=1/`K' {
					forvalues k=1/`K' {
						// abs(.) needed?
						mata: `Vi'[`j',`k'] = `Vi'[`j',`k'] + abs((`bvec'[`m',`j'] - `bagg'[1,`j'])*(`bvec'[`m',`k'] - `bagg'[1,`k']))
					}
				}
				mata: `Vagg' = `Vagg' + 1:/`Vi'
			}
			mata: `Vagg' = `nreps' :/ `Vagg'
			mata: st_matrix("`Vagg'",`Vagg')
		}
		else if "`medmean'"=="md" {
			// median VCV
			// inefficient - does off-diagonals twice
			forvalues j=1/`K' {
				forvalues k=1/`K' {
					forvalues m=1/`nreps' {
						mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
						mata: `Vi' = `B'.get(("V","post"))
						mata: `Vvec'[`m'] = `Vi'[`j',`k']
					}
					// adjustment as per
					// https://docs.doubleml.org/stable/guide/resampling.html#repeated-cross-fitting-with-k-folds-and-m-repetition
					// (generalized to multiple D variables)
					mata: `Vvec' = `Vvec' + abs((`bvec'[.,`j'] :- `bagg'[1,`j']):*(`bvec'[.,`k'] :- `bagg'[1,`k']))
					mata: `sVvec' = sort(`Vvec',1)
					if `isodd' {
						mata: `Vagg'[`j',`k'] = `sVvec'[`medrow',1]
					}
					else {
						mata: `Vagg'[`j',`k'] = (`sVvec'[`medrow',1] + `sVvec'[`medrow'+1,1])/2
					}
				}
			}
			mata: st_matrix("`Vagg'",`Vagg')
		}
		else {
			di as err "_ddml_reg error - unrecognized option `medmean'"
			exit 198
		}
	
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		mata: `A'.put(("N","post"),`N')
		mata: `A'.put(("b","post"),`bagg')
		mata: `A'.put(("V","post"),`Vagg')
		mata: `A'.put(("depvar","post"),"`depvar'")
		
		// store locals
		local list_local title y d dh z yname dnames vce vcetype
		if "`clustvar'"~=""		local list_local `list_local' clustvar
		foreach obj in `list_local' {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		// store scalars
		local list_scalar
		if "`clustvar'"~=""		local list_scalar `list_scalar' N_clust
		foreach obj in `list_scalar' {
			mata: `A'.put(("`obj'","scalar"),``obj'')
		}
		// special case - "_m" subscript doesn't apply to mean/median over resamplings
		// so store without resample subscript
		foreach obj in y d dh z {
			mata: `A'.put(("`obj'_m","local"),"``obj''")
		}
		// special case - vector of betas available only for mean/median
		mata: `A'.put(("b_resamples","matrix"),`bvec')
		
		// store AA with median/mean results
		mata: (`mname'.estAA).put(("`spec'","`medmean'"),`A')
		
		// no longer needed
		foreach obj in `A' `B' `bagg' `bvec' `sbvec' `Vagg' `Vvec' `sVvec' `Vi' {
			cap mata: mata drop `obj'
		}
		
	}
	else {
		// replay
				
		tempname B keys isscalar islocal ismatrix

		mata: `B' = AssociativeArray()
		mata: `B' = (`mname'.estAA).get(("`spec'","`rep'"))
		mata: `keys' = `B'.keys()
		mata: st_local("nentries",strofreal(rows(`keys')))
		mata: `isscalar'	= (`keys'[.,2] :== "scalar")
		mata: `islocal'		= (`keys'[.,2] :== "local")
		mata: `ismatrix'	= (`keys'[.,2] :== "matrix")
		
		tempname b V
		mata: st_matrix("`b'",`B'.get(("b","post")))
		mata: st_matrix("`V'",`B'.get(("V","post")))
		mata: st_local("N",strofreal(`B'.get(("N","post"))))
		mata: st_local("depvar",`B'.get(("depvar","post")))
		
		mata: st_local("yname",`B'.get(("yname","local")))
		mata: st_local("dnames",`B'.get(("dnames","local")))
		
		matrix rownames `b' = `depvar'
		matrix colnames `b' = `dnames'
	 	matrix colnames `V' = `dnames'
		matrix rownames `V' = `dnames'
		
		tempvar esample
		cap gen `esample' = `mname'_sample_`rep'
		if _rc>0 {
			// sample variable doesn't exist; ignore
			local esample
		}
		
		ereturn clear
		ereturn post `b' `V', depname(`depvar') obs(`N') esample(`esample')
		
		ereturn local cmd ddml
		ereturn local model `model'
		ereturn local rep `rep'
		ereturn local spec `spec'
		ereturn local tmname `mname'	// temporary mname
		
		// extract and post scalars, locals, matrices
		forvalues i=1/`nentries' {
			mata: st_local("topost",strofreal(`isscalar'[`i']))
			if `topost' {
				mata: st_local("sname",substr(`keys'[`i',1],1,32))
				mata: st_numscalar("e(`sname')",`B'.get(`keys'[`i',.]))
			}
		}
		forvalues i=1/`nentries' {
			mata: st_local("topost",strofreal(`islocal'[`i']))
			if `topost' {
				mata: st_local("lname",substr(`keys'[`i',1],1,32))
				mata: st_global("e(`lname')",`B'.get(`keys'[`i',.]))
			}
		}
		forvalues i=1/`nentries' {
			mata: st_local("topost",strofreal(`ismatrix'[`i']))
			if `topost' {
				mata: st_local("tmname",substr(`keys'[`i',1],1,32))
				mata: st_matrix("e(`tmname')",`B'.get(`keys'[`i',.]))
			}
		}
		
		// no longer needed
		foreach obj in `B' `keys' `isscalar' `islocal' `ismatrix' {
			cap mata: mata drop `obj'
		}
		
		// display results
		di as text "`e(title)'"
		di as text "y-E[y|X]" _col(11) "= " as res "`e(y_m)'" _c
		di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
		if "`e(model)'"~="ivhd" {
			di as text "D-" _c
		}
		di as text "E[D|X,Z]" _col(11)  "= " as res "`e(d_m)'"
		if "`e(model)'" == "iv" {
			di as text "Z-E[Z|X]" _col(11) "= " as res "`e(z_m)'"
		}
		else if "`e(model)'" == "ivhd" {
			di as text "E[D|X] = " as res "`e(dh_m)'"
		}
		if "`e(model)'" == "ivhd" {
			di as text "Orthogonalised D = D - E[D|X]; optimal IV = E[D|X,Z] - E[D|X]."
		}
		ereturn display
		
		// report warning if clustered SEs requested but doesn't match clustered crossfitting
		mata: st_local("fclustvar",`mname'.fclustvar)
		if "`e(clustvar)'"~="" {
			if "`fclustvar'"=="" {
				di as res "Warning" as text ": crossfit folds do not respect cluster structure used for VCE."
			}
			else if "`fclustvar'"~="`e(clustvar)'" {
				di as res "Warning" as text ": cluster variable for VCE does not match cluster variable for crossfit folds."
			}
		}
	}

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
