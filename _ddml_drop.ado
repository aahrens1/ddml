*! ddml v1.0
*! last edited: 13 dec 2022
*! authors: aa/ms

program _ddml_drop, eclass
	version 13

	syntax , mname(name)		// will already have verified that mname is a valid ddml mStruct

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// locals used below
	mata: st_local("model",`mname'.model)
	
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	local numeqnY	: word count `nameY'				// can be zero if no variables created yet
	local numeqnD	: word count `nameD'
	local numeqnZ	: word count `nameZ'

	if `numeqnY' {
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		// equation struct naming is the model name + dep variable name
		local eqnlist `mname'_`nameY'
		mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
		mata: st_local("ssvname",invtokens(`eqn'.shortstack))
		local vtlist `vtlistY' `ssvname'
	}
	
	if `numeqnD' {
		foreach var in `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			local eqnlist `eqnlist' `mname'_`var'
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			mata: st_local("ssvname",invtokens(`eqn'.shortstack))
			local vtlist `vtlist' `vtlistD' `ssvname'
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			if `lieflag' {
				foreach vn in `vtlistD' {
					local vtlistD_h `vtlistD_h' `vn'_h
				}
				local vtlist `vtlist' `vtlistD_h'			
			}
		}
	}
	
	if `numeqnZ' {
		foreach var in `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			local eqnlist `eqnlist' `mname'_`var'
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			mata: st_local("ssvname",invtokens(`eqn'.shortstack))
			local vtlist `vtlistZ' `ssvname'
		}
	}
	
	*** Add wildcards, then unabbreviate
	foreach var in `vtlist' {
		// variables may not exist yet
		cap unab evtlist : `var'*
		if _rc==0 {
			// variables exist so drop them
			foreach var of varlist `evtlist' {
				cap drop `var'
			}
		}
	}
	
	*** drop id, fold id, sample var
	cap drop `mname'_id
	cap drop `mname'_sample*
	cap drop `mname'_fid*		// multiple folds
	
	*** drop eqn structs
	foreach estruct in `eqnlist' {
		// eqn structs may not exist yet
		cap mata: mata drop `estruct'
	}
	mata: mata drop `eqn'		// temp eqn

	*** drop model structs
	mata: mata drop `mname'

end
