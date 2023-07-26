*! ddml v1.4
*! last edited: 26july2023
*! authors: aa/ms

program _ddml_export, rclass
	version 16

	syntax , mname(name) fname(string) [ addvars(varlist) * ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	*** extract details of estimation
	mata: model_chars(`mname')
	local nreps		= r(nreps)
	local numeqnD	= r(numeqnD)
	local numeqnZ	= r(numeqnZ)

	// fold IDs
	forvalues m=1/`nreps' {
		local fidlist `fidlist' `mname'_fid_`m'
	}
	
	// check
	if r(crossfitted)==0 {
		di as err "error - model not yet crossfitted, no variables to export"
		exit 198
	}
	
	// collect names of Y variables
	local vlist `r(Y)' `r(Y_L)'

	// collect names of D variables
	forvalues i=1/`numeqnD' {
		local vlist `vlist' `r(D`i')' `r(D`i'_L)' `r(D`i'_h)'
	}
	// collect names of Z variables
	forvalues i=1/`numeqnZ' {
		local vlist `vlist' `r(Z`i')' `r(Z`i'_L)'
	}

	// add rep numbers
	foreach vn in `vlist' {
		forvalues m=1/`nreps' {
			local vreplist `vreplist' `vn'_`m'
		}
	}

	// preserve, drop unneeded vars, rename, export, restore
	preserve

	keep `addvars' `mname'_sample* `fidlist' `vreplist'
	order `addvars' `mname'_sample* `fidlist' `vreplist'
	export delimited using "`fname'", `options'
	
	restore
	
	fvexpand `addvars' `mname'_sample* `fidlist' `vreplist'

	local numvars : word count `r(varlist)'
	return scalar numvars	= `numvars'
	return local vlist		`vlist'
	return local vreplist	`vreplist'
	
end

********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:


end
