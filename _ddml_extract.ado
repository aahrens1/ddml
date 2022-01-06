program define _ddml_extract, rclass

	syntax [name] , [				///
				mname(name)			///
				vname(name)			///
				keys				///
				key1(string)		///
				key2(string)		///
				key3(string)		///
				subkey1(string)		///
				subkey2(string)		///
				*					///
				]
	
	if "`namelist'" ~= "" {
		// Mata assignment
		local assign `namelist'=
	}
	
	tempname eqn A
	mata: `eqn' = init_eStruct()
	mata: `A' = AssociativeArray()
	
	if "`keys'"~="" {
		di as text "AA keys for `mname'.eqnAA:"
		mata: (`mname'.eqnAA).keys()
		di as text "AA keys for `mname'.estAA:"
		mata: (`mname'.estAA).keys()
		tempname keymat
		mata: `keymat' = (`mname'.estAA).keys()
		mata: st_local("nentries",strofreal(rows(`keymat')))
		forvalues i=1/`nentries' {
			mata: st_local("classname",classname((`mname'.estAA).get((`keymat'[`i',.]))))
			if "`classname'"=="AssociativeArray" {
				mata: st_local("k1",`keymat'[`i',1])
				mata: st_local("k2",`keymat'[`i',2])
				di as text "AA keys for `mname'.estAA, key 1 = `k1', key 2 = `k2':"
				mata: `A' = (`mname'.estAA).get((`keymat'[`i',.]))
				mata: `A'.keys()
			}
		}
		mata: mata drop `keymat'
	}
	
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	if "`keys'"~="" {
		foreach e in `nameY' `nameD' `nameZ' {
			mata: `eqn' = (`mname'.eqnAA).get("`e'")
			di as text "AA keys for eqn `e'.lrnAA"
			mata: (`eqn'.lrnAA).keys()
			di as text "AA keys for eqn `e'.resAA"
			mata: (`eqn'.resAA).keys()
		}
	}
	else if "`vname'"=="" {
		// no vname means extract from the model struct estAA. 2 keys.
		mata: st_local("classname",classname((`mname'.estAA).get(("`key1'","`key2'"))))
		if "`classname'"=="AssociativeArray" {
			// AA with estimation results, 2 keys
			mata: `A' = (`mname'.estAA).get(("`key1'","`key2'"))
			mata: `assign' `A'.get(("`subkey1'","`subkey2'"))
		}
		else {
			// not an AA
			mata: `assign' (`mname'.estAA).get(("`key1'","`key2'"))
		}
	}
	else {
		// vname is in either nameY, nameD or nameZ
		mata: `eqn' = (`mname'.eqnAA).get("`vname'")
		if "`key3'"=="" {
			// only 2 keys, it's lrnAA
			mata: `assign' (`eqn'.lrnAA).get(("`key1'","`key2'"))
		}
		else {
			// 3 keys, it's resAA		
			mata: `assign' (`eqn'.resAA).get(("`key1'","`key2'","`key3'"))
		}
	}
	
	mata: mata drop `eqn'
	mata: mata drop `A'
	
end