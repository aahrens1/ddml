*! ddml v1.5.0
*! last edited: 18dec2025
*! authors: aa/ms
* utility to take list of coefficient vectors and VCVs and return median/mean aggregation
* can be used interactively
* syntax: either (1) provide Stata matrix of stacked coefficient vectors and stacked VCVs
*         or     (2) provide a namelist of stored Stata estimation results
* b and V are posted and returned in e(b) and e(V)

program _ddml_medmean, eclass
	version 16
	syntax [anything] [if] [in] ,				///
						[						///
							bmat(name)			/// matrix of b (row) Stata vectors
							vmat(name)			/// matrix of (stacked) V Stata matricies
							elist(string)		/// alt - provide list of Stata stored results
							cnames(string)		/// column/row names for b/V
							ROBust				///
							CLUster(varname)	///
							DROPCONStant		/// drop last col/row (=_cons)
							mean				/// default is empty = median
						]

	// assumes no constant in b or V

	// syntax checks
	if "`bmat'"~="" & "`elist'"~="" {
		di as err "error - incompatible options bmat(.) and elist(.)
		exit 198
	}
	if "`vmat'"~="" & "`elist'"~="" {
		di as err "error - incompatible options vmat(.) and elist(.)
		exit 198
	}
	if "`bmat'"=="" & "`elist'"=="" {
		di as err "error - missing option bmat(.) or elist(.)"
		exit 198
	}
	if "`vmat'"=="" & "`elist'"=="" {
		di as err "error - missing option vmat(.) or elist(.)"
		exit 198
	}
	
	// if estimates provided as a list of stored results, append into matrices
	if "`elist'"~="" {
		tempname bmat vmat
		foreach ename in `elist' {
			qui est restore `ename'
			mat `bmat' = nullmat(`bmat') \ e(b)
			mat `vmat' = nullmat(`vmat') \ e(V)
		}
	}
	
	// drop constant (last column) if requested
	if "`dropconstant'"~="" {
		// stacked vertically, so just drop last column
		local numcols = colsof(`bmat')
		local numrows = rowsof(`vmat')
		mat `bmat' = `bmat'[1...,1..(`numcols'-1)]
		tempname vmat2
		mat `vmat2' = `vmat'[1...,1..(`numcols'-1)]
		// drop last row of each V submatrix
		mata: st_matrix("`vmat'",select(st_matrix("`vmat2'"),mod(runningsum(J(`numrows',1,1)),`numcols') :~= 0))
	}

	// if cnames not provided, take from bmat
	if "`cnames'"=="" {
		local cnames : colnames `bmat'
	}
			
	// prep
	tempname bagg Vagg sbmat V_i Vvec sVvec
	local nreps = rowsof(`bmat')
	local K = colsof(`bmat')
	local isodd = mod(`nreps',2)
	local medrow = ceil(`nreps'/2)	
	mata: `bmat' = st_matrix("`bmat'")
	mata: `vmat' = st_matrix("`vmat'")
	mata: `Vagg' = J(`K',`K',0)
	mata: `Vvec' = J(`nreps',1,0)
	
	// mean/median of b
	if "`mean'"=="mean" {
		mata: `bagg' = mean(`bmat')
	}
	else {
		// initialize
		mata: `bagg' = J(1,`K',.)
		forvalues k=1/`K' {
			// leave order of bmat unchanged; sbmat = sorted bmat
			mata: `sbmat' = sort(`bmat',`k')
			if `isodd' {
				mata: `bagg'[1,`k'] = `sbmat'[`medrow',`k']
			}
			else {
				mata: `bagg'[1,`k'] = (`sbmat'[`medrow',`k'] + `sbmat'[`medrow'+1,`k'])/2
			}
		}
	}
	mata: st_matrix("`bagg'",`bagg')
	
	// mean/median of V
	if "`mean'"=="mean" {
		// harmonic mean
		// inefficient - does off-diagonals twice
		forvalues m=1/`nreps' {
			mata: `V_i' = `vmat'[((`m'-1)*`K'+1)::(`m'*`K'),(1..`K')]
			forvalues j=1/`K' {
				forvalues k=1/`K' {
					// abs(.) needed?
					mata: `V_i'[`j',`k'] = `V_i'[`j',`k'] + abs((`bmat'[`m',`j'] - `bagg'[1,`j'])*(`bmat'[`m',`k'] - `bagg'[1,`k']))
				}
			}
			mata: `Vagg' = `Vagg' + 1:/`V_i'
		}
		mata: `Vagg' = `nreps' :/ `Vagg'
	}
	else {
		// median VCV
		// inefficient - does off-diagonals twice
		forvalues j=1/`K' {
			forvalues k=1/`K' {
				forvalues m=1/`nreps' {
					mata: `V_i' = `vmat'[((`m'-1)*`K'+1)::(`m'*`K'),(1..`K')]
					mata: `Vvec'[`m'] = `V_i'[`j',`k']
				}
				// adjustment as per
				// https://docs.doubleml.org/stable/guide/resampling.html#repeated-cross-fitting-with-k-folds-and-m-repetition
				// (generalized to multiple D variables)
				mata: `Vvec' = `Vvec' + abs((`bmat'[.,`j'] :- `bagg'[1,`j']):*(`bmat'[.,`k'] :- `bagg'[1,`k']))
				mata: `sVvec' = sort(`Vvec',1)
				if `isodd' {
					mata: `Vagg'[`j',`k'] = `sVvec'[`medrow',1]
				}
				else {
					mata: `Vagg'[`j',`k'] = (`sVvec'[`medrow',1] + `sVvec'[`medrow'+1,1])/2
				}
			}
		}
	}
	mata: st_matrix("`Vagg'",`Vagg')
	mat colnames `bagg'		= `cnames'
	mat colnames `Vagg' 	= `cnames'
	mat rownames `Vagg' 	= `cnames'

	if "`robust'"~="" & "`cluster'"=="" {
		local vce		robust
		local vcetype	Robust
	}
	else if "`cluster'"~="" {
		local vce		cluster
		local vcetype	Robust
		local clustvar	`cluster'
	}
	else {
		local vce		ols
	}
		
	ereturn post `bagg' `Vagg'
	ereturn local vce		`vce'
	ereturn local vcetype	`vcetype'
	ereturn local clustvar	`clustvar'

	di	
	ereturn display
	
	// clean up mata
	foreach obj in `bmat' `vmat' `bagg' `Vagg' `sbmat' `Vvec' `sVvec' `V_i' {
		cap mata: mata drop `obj'
	}

end
