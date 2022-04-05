
program _ddml_ate_late, eclass

	version 13

	syntax [anything] [if] [in] ,				///
						[						///
							yvar(varname)		///
							dvar(varname)		///
							zvar(varname)		///
							y0tilde(name)		///
							y1tilde(name)		///
							dtilde(name)		///
							d0tilde(name)		///
							d1tilde(name)		///
							ztilde(name)		///
							mname(name)			///
							spec(string)		///
							rep(string)			///
							replay				///
							title(string)		///
							medmean(string)		///
							clustervar(varname) ///
							ATET ///
						]
		
	mata: st_local("model",`mname'.model)
	
	marksample touse
	
	if "`replay'"=="" & "`medmean'"=="" {	// estimate from scratch
		
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		// 0/1 etc
		local y0		`y0tilde'
		local y1		`y1tilde'
		local d			`dtilde'
		local d0		`d0tilde'
		local d1		`d1tilde'
		local z			`ztilde'
		// add suffixes
		local y0_m		`y0tilde'0_`rep'
		local y1_m		`y1tilde'1_`rep'
		local d_m		`dtilde'_`rep'
		local d0_m		`d0tilde'0_`rep'
		local d1_m		`d1tilde'1_`rep'
		local z_m		`ztilde'_`rep'
		
		if "`model'"=="interactive" {
			mata: ATE("`atet'","`yvar'","`dvar'","`y0_m'", "`y1_m'", "`d_m'","`touse'","`b'","`V'","`clustervar'")
		}
		else {
			mata: LATE("`yvar'","`dvar'","`zvar'","`y0_m'", "`y1_m'", "`d0_m'","`d1_m'","`z_m'","`touse'","`b'","`V'","`clustervar'")
		}
		
		// store post objects
		mata: `A'.put(("N","post"),`r(N)')
		mata: `A'.put(("b","post"),st_matrix("r(b)"))
		mata: `A'.put(("V","post"),st_matrix("r(V)"))
		mata: `A'.put(("depvar","post"),"`yvar'")
		
		// for calling program
		ereturn clear
		mata: st_matrix("e(bmat)",st_matrix("r(b)"))
		mata: st_matrix("e(semat)",sqrt(diagonal(st_matrix("r(V)"))'))
		
		// store locals
		foreach obj in title y0 y0_m y1 y1_m d d_m d0 d0_m d1 d1_m z z_m yvar dvar {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		
		// additional estimation results
		tempname eqn
		mata: `eqn' = init_eStruct()
		// Y eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`yvar'")
		// MSE
		mata: `A'.put(("`y0tilde'_mse","scalar"),return_result_item(`eqn',"`y0tilde'","MSE0","`rep'"))
		mata: `A'.put(("`y1tilde'_mse","scalar"),return_result_item(`eqn',"`y1tilde'","MSE1","`rep'"))
		// MSE folds
		mata: `A'.put(("`y0tilde'_mse_folds","matrix"),return_result_item(`eqn',"`y0tilde'","MSE0_folds","`rep'"))
		mata: `A'.put(("`y1tilde'_mse_folds","matrix"),return_result_item(`eqn',"`y1tilde'","MSE1_folds","`rep'"))
		// pystacked weightss - don't use, mean across folds doesn't mean anything
		// mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`y0tilde'","cmd")=="pystacked"))
		// if `pyswflag' {
		//	mata: `A'.put(("`y0tilde'_pysw","matrix"),mean(return_result_item(`eqn',"`y0tilde'","stack_weights0","`rep'")'))
		// }
		// mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`y1tilde'","cmd")=="pystacked"))
		// if `pyswflag' {
		//	mata: `A'.put(("`y1tilde'_pysw","matrix"),mean(return_result_item(`eqn',"`y1tilde'","stack_weights1","`rep'")'))
		// }
		if "`model'"=="interactive" {
			// D eqn results
			mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
			// MSE
			mata: `A'.put(("`dtilde'_mse","scalar"),return_result_item(`eqn',"`dtilde'","MSE","`rep'"))
			// MSE folds
			mata: `A'.put(("`dtilde'_mse_folds","matrix"),return_result_item(`eqn',"`dtilde'","MSE_folds","`rep'"))
			// pystacked weightss - don't use, mean across folds doesn't mean anything
			// mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`dtilde'","cmd")=="pystacked"))
			// if `pyswflag' {
			// 	mata: `A'.put(("`dtilde'_pysw","matrix"), mean(return_result_item(`eqn',"`dtilde'","stack_weights","`rep'")'))
			// }
		}
		else {
			mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
			// MSE, D
			mata: `A'.put(("`d0tilde'_mse","scalar"),return_result_item(`eqn',"`d0tilde'","MSE0","`rep'"))
			mata: `A'.put(("`d1tilde'_mse","scalar"),return_result_item(`eqn',"`d1tilde'","MSE1","`rep'"))
			// MSE folds, D
			mata: `A'.put(("`d0tilde'_mse_folds","matrix"),return_result_item(`eqn',"`d0tilde'","MSE0_folds","`rep'"))
			mata: `A'.put(("`d1tilde'_mse_folds","matrix"),return_result_item(`eqn',"`d1tilde'","MSE1_folds","`rep'"))
			// pystacked weights, Ds - don't use, mean across folds doesn't mean anything
			// mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`d0tilde'","cmd")=="pystacked"))
			// if `pyswflag' {
			//	mata: `A'.put(("`d0tilde'_pysw","matrix"), mean(return_result_item(`eqn',"`d0tilde'","stack_weights0","`rep'")'))
			// }
			// mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`d1tilde'","cmd")=="pystacked"))
			// if `pyswflag' {
			//	mata: `A'.put(("`d1tilde'_pysw","matrix"), mean(return_result_item(`eqn',"`d1tilde'","stack_weights1","`rep'")'))
			// }
			mata: `eqn' = (`mname'.eqnAA).get("`zvar'")
			// MSE, Z
			mata: `A'.put(("`ztilde'_mse","scalar"),return_result_item(`eqn',"`ztilde'","MSE","`rep'"))
			// MSE folds, Z
			mata: `A'.put(("`ztilde'_mse_folds","matrix"),return_result_item(`eqn',"`ztilde'","MSE_folds","`rep'"))
			// pystacked weights, Zs - don't use, mean across folds doesn't mean anything
			// mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`ztilde'","cmd")=="pystacked"))
			// if `pyswflag' {
			//	mata: `A'.put(("`ztilde'_pysw","matrix"), mean(return_result_item(`eqn',"`ztilde'","stack_weights","`rep'")'))
			// }
		}
	
		mata: (`mname'.estAA).put(("`spec'","`rep'"),`A')
		
		// no longer needed
		foreach obj in `A' `eqn' {
			cap mata: mata drop `obj'
		}
	
	}
	else if "`replay'"=="" & "`medmean'"~="" {	// aggregate over resamples
		
		tempname b V bagg Vagg Vi
		tempname bvec bmed Vvec Vmed
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
		
		mata: `bagg' = J(1,`K',0)
		mata: `bvec' = J(`nreps',`K',0)
		forvalues m=1/`nreps' {
			mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
			mata: `bvec'[`m',.] = `B'.get(("b","post"))
			// row/colnames etc. - need to do this only once
			if `m'==1 {
				mata: st_local("depvar",`B'.get(("depvar","post")))
				// retrieve locals
				foreach obj in y0 y0_m y1 y1_m d d_m d0 d0_m d1 d1_m z z_m yvar dvar {
					mata: st_local("`obj'",`B'.get(("`obj'","local")))
				}
			// possible that different estimation samples have different #obs
			qui count if `mname'_sample_`m'==1
			local N = `N' + r(N)
			}
		}
		local N = round(`N'/`nreps')
		
		if "`medmean'"=="mn" {
			// default is mean
			mata: st_matrix("`bagg'",mean(`bvec'))
		}
		else if "`medmean'"=="md" {
			// median beta
			forvalues k=1/`K' {
				mata: _sort(`bvec',`k')
				if `isodd' {
					mata: `bagg'[1,`k'] = `bvec'[`medrow',`k']
				}
				else {
					mata: `bagg'[1,`k'] = (`bvec'[`medrow',`k'] + `bvec'[`medrow'+1,`k'])/2
				}
			}
			mata: st_matrix("`bagg'",`bagg')
		}
		else {
			di as err "_ddml_ate_late error - unrecognized option `medmean'"
			exit 198
		}
		
		mata: `Vagg' = J(`K',`K',0)
		mata: `Vvec' = J(`nreps',1,0)
		if "`medmean'"=="mn" {
			forvalues m=1/`nreps' {
				mata: `B' = (`mname'.estAA).get(("`spec'","`m'"))
				mata: `Vagg' = `Vagg' + 1/`nreps' * `B'.get(("V","post"))
			}
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
					mata: _sort(`Vvec',1)
					if `isodd' {
						mata: `Vagg'[`j',`k'] = `Vvec'[`medrow',1]
					}
					else {
						mata: `Vagg'[`j',`k'] = (`Vvec'[`medrow',1] + `Vvec'[`medrow'+1,1])/2
					}
				}
			}
			mata: st_matrix("`Vagg'",`Vagg')
		}
		else {
			di as err "_ddml_ate_late error - unrecognized option `medmean'"
			exit 198
		}
	
		tempname A
		mata: `A' = AssociativeArray()
		mata: `A'.reinit("string",2)
		mata: `A'.notfound("")				// so that if a local isn't found, it's an empty string
		
		mata: `A'.put(("N","post"),`N')
		mata: `A'.put(("b","post"),st_matrix("`bagg'"))
		mata: `A'.put(("V","post"),st_matrix("`Vagg'"))
		mata: `A'.put(("depvar","post"),"`depvar'")
		
		// store locals
		foreach obj in title y0 y1 d d0 d1 z yvar dvar {
			mata: `A'.put(("`obj'","local"),"``obj''")
		}
		// special case - "_m" subscript doesn't apply to mean/median over resamplings
		// so store without resample subscript
		foreach obj in title y0 y1 d d0 d1 z {
			mata: `A'.put(("`obj'_m","local"),"``obj''")
		}
		
		mata: (`mname'.estAA).put(("`spec'","`medmean'"),`A')
		
		// no longer needed
		foreach obj in `A' `B' `bagg' `bvec' `Vagg' `Vvec' `Vi' {
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
		
		mata: st_local("yvar",`B'.get(("yvar","local")))
		mata: st_local("dvar",`B'.get(("dvar","local")))
		
		matrix rownames `b' = `depvar'
		matrix colnames `b' = `dvar'
	 	matrix colnames `V' = `dvar'
		matrix rownames `V' = `dvar'
		
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
		ereturn local mname `mname'
		
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
				mata: st_local("mname",substr(`keys'[`i',1],1,32))
				mata: st_matrix("e(`mname')",`B'.get(`keys'[`i',.]))
			}
		}
		
		// no longer needed
		foreach obj in `B' `keys' `isscalar' `islocal' `ismatrix' {
			cap mata: mata drop `obj'
		}
		
		// display results
		di as text "`e(title)'"
		di as text "y-E[y|X,D=0]" _col(14) "= " as res "`e(y0_m)'" _c
		di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
		di as text "y-E[y|X,D=1]" _col(14) "= " as res "`e(y1_m)'"
		if "`e(model)'"=="interactive" {
			di as text "D-E[D|X]" _col(14)  "= " as res "`e(d_m)'"
		}
		else {
			di as text "D-E[D|X,Z=0]" _col(14)  "= " as res "`e(d0_m)'"
			di as text "D-E[D|X,Z=1]" _col(14)  "= " as res "`e(d1_m)'"
			di as text "Z-E[Z|X]" _col(14)  "= " as res "`e(z_m)'"
		}
		ereturn display
	}
	
end


********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

void ATE(   
			string scalar atet,
			string scalar yvar,	 // Y
			string scalar dvar,	 // D
			string scalar y0tilde,  // E[Y|X,D=0]
			string scalar y1tilde,  // E[Y|X,D=1]
			string scalar dtilde,   // E[D|X]
			string scalar sample,   // sample
			string scalar outate,   // output: name of matrix to store b
			string scalar outatese,  // output: name of matrix to store V
			string scalar clustervar
			)
{
	st_view(my_d0x,.,y0tilde,sample)
	st_view(my_d1x,.,y1tilde,sample)
	st_view(md_x,.,dtilde,sample)
	st_view(d,.,dvar,sample)
	st_view(y,.,yvar,sample)
	if (clustervar!="") {
		st_view(clustid,.,clustervar,sample)
		clusterid_uni=uniqrows(clustid)
		nclust = rows(clusterid_uni)
	}

	n = rows(y)

	// psi = psi_b - psi_a*theta
	if (atet=="") {
		psi_b  = (d :* (y :- my_d1x) :/ md_x) :-  ((1 :- d) :* (y :- my_d0x) :/ (1 :- md_x)) :+ my_d1x :- my_d0x 
		psi_a  = J(n,1,-1) 
	}
	else {
		p_hat = mean(d)
		psi_b = (d :* (y :- my_d0x) :/ p_hat) :-  md_x :* (1 :- d) :* (y :- my_d0x) :/ (p_hat :*(1 :- md_x)) 
		psi_a = -d :/ p_hat 
	}
	theta = -mean(psi_b) / mean(psi_a)
	psi = psi_a :* theta :- psi_b

	if (clustervar=="") {
		V =  mean(psi:^2) :/ mean(psi_a):^2 :/ n	
	}
	else {
		gamma = 0
		jhat = 0
		for (i=1;i<=nclust;i++) {
			psi_c = select(psi,clustid:==clusterid_uni[i,1])
			psi_a_c = select(psi_a,clustid:==clusterid_uni[i,1])
			gamma = gamma :+ 1/nclust :* sum(psi_c*psi_c')
			jhat = jhat :+  1/nclust :* sum(psi_a_c)
		}
		V = gamma / jhat:^2 / nclust
	}

	st_numscalar("r(N)",n)
	st_matrix("r(b)",theta)
	st_matrix("r(V)",V)
}

void LATE(  string scalar yvar,      // Y
            string scalar dvar,      // D
            string scalar zvar,      // Z
            string scalar y0tilde,   // E[Y|X,Z=0]
            string scalar y1tilde,   // E[Y|X,Z=1]
            string scalar d0tilde,   // E[D|X,Z=0]
            string scalar d1tilde,   // E[D|X,Z=1]
            string scalar ztilde,    // E[Z|X]
            string scalar sample,    // sample
            string scalar outlate,   // output: name of matrix to store b
            string scalar outlatese,  // output: name of matrix to store V
            string scalar clustervar
            )
{
    st_view(my_z0x,.,y0tilde,sample)
    st_view(my_z1x,.,y1tilde,sample)
    st_view(md_z0x,.,d0tilde,sample)
    st_view(md_z1x,.,d1tilde,sample)
    st_view(mz_x,.,ztilde,sample)
    st_view(d,.,dvar,sample)
    st_view(y,.,yvar,sample)
    st_view(z,.,zvar,sample)
    if (clustervar!="") {
			st_view(clustid,.,clustervar,sample)
			clusterid_uni=uniqrows(clustid)
			nclust = rows(clusterid_uni)
		}

    n = rows(y)

    psi_b =  mean( z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x )
    psi_a =  mean( z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x ) 

		theta = -mean(psi_b) / mean(psi_a)
		psi = psi_a :* theta :- psi_b

		if (clustervar=="") {
			V =  mean(psi:^2) :/ mean(psi_a):^2 :/ n	
		}
		else {
			gamma = 0
			jhat = 0
			for (i=1;i<=nclust;i++) {
				psi_c = select(psi,clustid:==clusterid_uni[i,1])
				psi_a_c = select(psi_a,clustid:==clusterid_uni[i,1])
				gamma = gamma :+ 1/nclust :* sum(psi_c*psi_c')
				jhat = jhat :+  1/nclust :* sum(psi_a_c)
			}
			V = gamma / jhat:^2 / nclust
		}

    st_numscalar("r(N)",n)
    st_matrix("r(b)",theta)
    st_matrix("r(V)",V)
}

end