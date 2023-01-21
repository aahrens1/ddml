*! ddml v1.2
*! last edited: 21 jan 2023
*! authors: aa/ms

program _ddml_save
	version 13

	syntax , mname(name) fname(string) [ replace ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()

	// locals used below
	mata: st_local("model",`mname'.model)
	
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens(`mname'.nameZ))
	local numeqnD	: word count `nameD'
	local numeqnZ	: word count `nameZ'
	
	mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	local vtlist `vtlistY'
	
	if `numeqnD' {
		foreach var of varlist `nameD' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("lieflag",strofreal(`eqn'.lieflag))
			mata: st_local("vtlistD",invtokens(`eqn'.vtlist))
			local vtlist `vtlist' `vtlistD'
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
		foreach var of varlist `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`var'")
			mata: st_local("vtlistZ",invtokens(`eqn'.vtlist))
			local vtlist `vtlistZ'
		}
	}
	
	// Add wildcards, add prefixed variables, then unabbreviate
	foreach var in `vtlist' {
		local evtlist `evtlist' `var'*
	}
	local evtlist `mname'_* `evtlist'
	unab evtlist : `evtlist'
	
	// insert onto model struct
	mata: `mname'.strDatavars = "`evtlist'"
	mata: `mname'.matDatavars = st_data(., "`evtlist'")
	
	if "`replace'"~="" {
		// Mata function to delete file
		mata: unlink("`fname'")
	}
	mata: save_model("`fname'",`mname')
	
	// clear from model struct
	mata: `mname'.strDatavars = ""
	mata: `mname'.matDatavars = J(0,0,.)
	
end

mata:

void save_model(					string scalar fname,
									struct mStruct m)
{
	fh = fopen(fname,"w")
	fputmatrix(fh,m)
	fclose(fh)
}

end
