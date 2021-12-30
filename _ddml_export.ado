
program _ddml_export
	version 13

	syntax , mname(name) fname(string) [ * ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	*** extract details of estimation
	
	// locals used below
	mata: st_local("model",`mname'.model)
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	mata: st_local("kfolds",strofreal(`mname'.kfolds))
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	local numeqnD	: word count `nameD'
	local numeqnZ	: word count `nameZ'
	// fold IDs
	forvalues m=1/`nreps' {
		local fidlist `fidlist' `mname'_fid_`m'
	}
	
	// check
	if ~`crossfitted' {
		di as err "error - model not yet crossfitted, no variables to export"
		exit 198
	}
	
	// collect names of Y variables
	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	mata: st_local("shortstack",invtokens(`eqn'.shortstack))
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	foreach vn in `vtlistY' {
		local vlistY `vlistY' `nameY'
	}
	if "`shortstack'"~="" {
		local vlistY `vlistY' `nameY'
		local vtlistY `vtlistY' `shortstack'
	}
	
	// collect names of D variables
	if `numeqnD' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("shortstack",invtokens(`eqn'.shortstack))
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistDi",invtokens(`eqn'.vtlist))
			local vtlistD `vtlistD' `vtlistDi'
			foreach vn in `vtlistDi' {
				local vlistD `vlistD' `var'
			}
			if "`shortstack'"~="" {
				local vlistD `vlistD' `var'
				local vtlistD `vtlistD' `shortstack'
			}
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			if `lieflag' {
				foreach vn in `vtlistDi' {
					local vlistDH `vlistDH' `var'
					local vtlistDH `vtlistDH' `vn'_h
				}
				if "`shortstack'"~="" {
					local vlistDH `vlistDH' `var'
					local vtlistDH `vtlistDH' `shortstack'_h
				}
			}
		}
	}
	
	// collect names of Z variables
	if `numeqnZ' {
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("shortstack",invtokens(`eqn'.shortstack))
			mata: st_local("vtlistZi",invtokens(`eqn'.vtlist))
			local vtlistZ `vtlistZ' `vtlistZi'
			foreach vn in `vtlistZi' {
				local vlistZ `vlistZ' `var'
			}
			if "`shortstack'"~="" {
				local vlistZ `vlistZ' `var'
				local vtlistZ `vtlistZ' `shortstack'
			}
		}
	}
	
	// treatment effect models have vt variables that come in pairs
	// Y
	if "`model'"=="interactive" | "`model'"=="late" {
		foreach vn in `vtlistY' {
			local vttvlistY `vttvlistY' `vn'0 `vn'1
		}
		local vtlistY `vttvlistY'
		local vlistY `vlistY' `vlistY'
	}
	// D and LATE
	if "`model'"=="late" {
		foreach vn in `vtlistD' {
			local vttvlistD `vttvlistD' `vn'0 `vn'1
		}
		local vtlistD `vttvlistD'
		local vlistD `vlistD' `vlistD'
	}

	// add rep numbers
	foreach vn in `vtlistY' {
		forvalues m=1/`nreps' {
			local vtreplistY `vtreplistY' `vn'_`m'
		}
	}
	local vtlistY `vtreplistY'
	
	foreach vn in `vtlistD' {
		forvalues m=1/`nreps' {
			local vtreplistD `vtreplistD' `vn'_`m'
		}
	}
	local vtlistD `vtreplistD'
	
	foreach vn in `vtlistDH' {
		forvalues m=1/`nreps' {
			local vtreplistDH `vtreplistDH' `vn'_`m'
		}
	}
	local vtlistDH `vtreplistDH'
	
	foreach vn in `vtlistZ' {
		forvalues m=1/`nreps' {
			local vtreplistZ `vtreplistZ' `vn'_`m'
		}
	}
	local vtlistZ `vtreplistZ'
	
	forvalues m=1/`nreps' {
		local vreplistY `vreplistY' `vlistY'
		local vreplistD `vreplistD' `vlistD'
		local vreplistDH `vreplistDH' `vlistDH'
		local vreplistZ `vreplistZ' `vlistZ'
	}
	local vlistY `vreplistY'
	local vlistD `vreplistD'
	local vlistDH `vreplistDH'
	local vlistZ `vreplistZ'
	
	// preserve, drop unneeded vars, rename, export, restore
	preserve
	keep `mname'_id `mname'_sample* `fidlist' `vtlistY' `vtlistD' `vtlistDH' `vtlistZ'

	local numvt	: word count `vtlistY'
	forvalues i=1/`numvt' {
		local vn	: word `i' of `vlistY'
		local vt	: word `i' of `vtlistY'
		local va	= strtoname("`vt'_Y_`vn'")
		rename `vt' `va'
	}
	local numvt	: word count `vtlistD'
	forvalues i=1/`numvt' {
		local vn	: word `i' of `vlistD'
		local vt	: word `i' of `vtlistD'
		local va	= strtoname("`vt'_D_`vn'")
		rename `vt' `va'
	}
	local numvt	: word count `vtlistDH'
	forvalues i=1/`numvt' {
		local vn	: word `i' of `vlistDH'
		local vt	: word `i' of `vtlistDH'
		local va	= strtoname("`vt'_DH_`vn'")
		rename `vt' `va'
	}
	local numvt	: word count `vtlistZ'
	forvalues i=1/`numvt' {
		local vn	: word `i' of `vlistZ'
		local vt	: word `i' of `vtlistZ'
		local va	= strtoname("`vt'_Z_`vn'")
		rename `vt' `va'
	}

	export delimited using "`fname'", `options'
	
	restore

	
end

********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:


end
