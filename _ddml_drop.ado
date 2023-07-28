*! ddml v1.4
*! last edited: 28july2023
*! authors: aa/ms

program _ddml_drop, eclass
	version 16

	syntax , mname(name)		// will already have verified that mname is a valid ddml mStruct
	
	*** extract details of estimation
	mata: model_chars(`mname')
	local nreps		= r(nreps)
	local numeqnD	= r(numeqnD)
	local numeqnZ	= r(numeqnZ)

	// fold IDs
	forvalues m=1/`nreps' {
		local fidlist `fidlist' `mname'_fid_`m'
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

	// drop vars may not exist, so use capture
	foreach vn in `vreplist' {
		cap confirm variable `vn', exact
		if _rc==0	drop `vn'
	}
	
	*** drop id, fold id, sample var
	cap drop `mname'_id
	cap drop `mname'_sample*
	cap drop `mname'_fid*		// multiple folds
	
	*** drop model struct
	cap mata: mata drop `mname'

end
