program define make_varlists2, rclass

	syntax [anything], mname(name)

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	// locals used below
	mata: st_local("numD",strofreal(cols(`mname'.nameD)))
	mata: st_local("numZ",strofreal(cols(`mname'.nameZ)))
 	mata: st_local("Ytilde",invtokens(`mname'.nameYtilde))
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlist)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))

	// vno_list has list of all orthogonalized variables
	// vn_list has list of all corresponding variables, original names
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_local("vname",`eqn'.Vname)
		mata: st_local("vtilde",`eqn'.Vtilde)
		mata: st_local("eqntype",`eqn'.eqntype)
		if "`eqntype'"=="deq" {
			local Dvn_list `Dvn_list' `vname'
			local Dvno_list `Dvno_list' `vtilde'
		}
		if "`eqntype'"=="zeq" {
			local Zvn_list `Zvn_list' `vname'
			local Zvno_list `Zvno_list' `vtilde'
		}
	}

	// get list of unique original varnames
	local Dvn_uniq	: list uniq Dvn_list
	local Zvn_uniq	: list uniq Zvn_list

	// loop through unique original varnames and create varlists
	// varlists are stored in Stata macros
	// each of these macros is a list of all orthogonalized vars for a given original var
	// macro orthog_lists is a macro with these macro names in it
	tempname Dv_list Do_list
	mata: `Dv_list' = tokens("`Dvn_list'")
	mata: `Do_list' = tokens("`Dvno_list'")
	foreach vn in `Dvn_uniq' {
		tempname Do_index
		mata: `Do_index' = `Dv_list' :== "`vn'"
		// select all orthogonalized variables for unique original varname vn
		mata: st_local("Do_sublist", invtokens(select(`Do_list',`Do_index')))
		tempname Dt_list
		local `Dt_list' `Do_sublist'
		local Dorthog_lists `Dorthog_lists' `Dt_list'
	}

	// do the same for Z
	if (`numZ'>0) {
		tempname Zv_list Zo_list
		mata: `Zv_list' = tokens("`Zvn_list'")
		mata: `Zo_list' = tokens("`Zvno_list'")
		foreach vn in `Zvn_uniq' {
			tempname Zo_index
			mata: `Zo_index' = `Zv_list' :== "`vn'"
			// select all orthogonalized variables for unique original varname vn
			mata: st_local("Zo_sublist", invtokens(select(`Zo_list',`Zo_index')))
			tempname Zt_list
			local `Zt_list' `Zo_sublist'
			local Zorthog_lists `Zorthog_lists' `Zt_list'
		}
	}

	// di "numD: `numD'"
	// di "numeqnsD: `numeqnsD'"
	// di "unique original varnames: `vn_uniq'"
	// di "orthog_lists: `orthog_lists'"
	foreach vl in `Dorthog_lists' {
		// di "`vl' = ``vl''"
		local Dtilde `Dtilde' - ``vl''
	}
	foreach vl in `Zorthog_lists' {
		// di "`vl' = ``vl''"
		local Ztilde `Ztilde' - ``vl''
	}



	// clear from Mata
	mata: mata drop `eqn'
	mata: mata drop `Dv_list' `Do_list' `Do_index' `Zv_list' `Zo_list' `Zo_index'

	return local eq `Ytilde' `Dtilde' `Ztilde'
	return scalar dpos_end = `numD' + 1
	return scalar dpos_start = 2
	if (`numZ'>0) {
		return scalar zpos_start = `numD' +2
		return scalar zpos_end = `numD' + `numZ' + 1
	}
end



mata: 
struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}
end
