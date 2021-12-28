program _ddml_medmean, eclass

	syntax [anything] ,						/// 
						[					///
						mname(name)			///
						spec(string)		/// specification
						nreps(integer 0)	/// number of resamples
						median				///
						cmd(string)			/// replay command name
						*					///
						]

	tempname b V bagg Vagg Vi
	tempname bvec bmed Vvec Vmed
	tempvar esample
	
	// initialize
	local isodd = mod(`nreps',2)
	local medrow = ceil(`nreps'/2)
	qui estimates restore `mname'_`spec'_1
	local K = colsof(e(b))
	// assume same for all estimates
	qui gen `esample' = e(sample)
	local N = e(N)
	// row/colnames inherited
	mat `b' = e(b) * 0
	mat `V' = e(V) * 0
	
	// inherited macros
	local yname		`e(depvar)'
	local robust	`e(robust)'
	local Robust	`e(Robust)'
	local model		`e(model)'
	// e(.) varlists are wrong because they have the rep subscript on them
	// fix each list
	foreach vl in y y0 y1 d d0 d1 z z0 z1 dh {
		local vl_m		`e(`vl')'
		foreach v_m in `vl_m' {
			local uloc	=strrpos("`v_m'","_")
			local v		=substr("`v_m'",1,`uloc'-1)
			local `vl'	``vl'' `v'
		}
	}
	
	
	mata: `bagg' = J(1,`K',0)
	mata: `bvec' = J(`nreps',`K',0)
	forvalues m=1/`nreps' {
		qui estimates restore `mname'_`spec'_`m'
		mata: `bvec'[`m',.] = st_matrix("e(b)")
	}
	if "`median'"=="" {
		// default is mean
		mata: st_matrix("`bagg'",mean(`bvec'))
	}
	else {
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
	
	mata: `Vagg' = J(`K',`K',0)
	mata: `Vvec' = J(`nreps',1,0)
	if "`median'"=="" {
		// default is mean
		forvalues m=1/`nreps' {
			qui estimates restore `mname'_`spec'_`m'
			mata: `Vagg' = `Vagg' + 1/`nreps' * st_matrix("e(V)")
		}
		mata: st_matrix("`Vagg'",`Vagg')
	}
	else {
		// median VCV
		// inefficient - does off-diagonals twice
		forvalues j=1/`K' {
			forvalues k=1/`K' {
				forvalues m=1/`nreps' {
					qui estimates restore `mname'_`spec'_`m'
					mata: `Vi' = st_matrix("e(V)")
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
	
	// b and V will have zeros and correct row/colnames
	mat `b' = `bagg' + `b'
	mat `V' = `Vagg' + `V'
	
	ereturn clear
	ereturn post `b' `V', depname(`yname') obs(`N') esample(`esample')

	ereturn local vce		`robust'
	ereturn local vctype	`Robust'
	ereturn local cmd		`cmd'
	ereturn local model		`model'
	ereturn local rep		mean
	foreach vl in y y0 y1 d d0 d1 z z0 z1 dh {
		ereturn local `vl'	``vl''
	}
	

end
