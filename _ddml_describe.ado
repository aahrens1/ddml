
program define _ddml_describe

	syntax name(name=mname), [NOCMD all]
	
	local showcmd	= ("`nocmd'"=="")
	local showall	= ("`all'"~="")

	// Will use for flagging minimized MSEs
	mata: st_local("Yopt",`mname'.nameYopt)
	mata: st_local("Dopt",invtokens(`mname'.nameDopt))
	mata: st_local("Zopt",invtokens(`mname'.nameZopt))


	mata: printf("{res}Model: %s\n", `mname'.model)
	di as res "ID: `mname'_id"
	di as res "Fold ID: `mname'_fid"
	mata: printf("{res}Dependent variable (Y): %s\n", `mname'.nameY)
	mata: printf("{res}Dependent variable (orthogonalized): %s\n", invtokens(`mname'.nameYtilde))
	di "Minimum MSE orthogonalized dep var: `Yopt'"
	mata: printf("{res}Causal variable(s) (D): %s\n", invtokens(`mname'.nameD))
	mata: printf("{res}Causal variable(s) (orthogonalized): %s\n", invtokens(`mname'.nameDtilde))
	di "Minimum MSE orthogonalized causal var: `Dopt'"

	if ("`model'"=="iv") {
		mata: printf("{res}Excluded instrumental variable(s): %s\n", `mname'.nameZ)
		mata: printf("{res}Excluded instrumental variable(s) (orthogonalized): %s\n", `mname'.nameZtilde)
		di "Minimum MSE orthogonalized IVs: `Zopt'"
	}
	
	// List equations
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlistY))
	local numeqnsY	= `r(numeqns)'
	di
	di "Number of Y estimating equations: `numeqnsY'"
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
		di "Estimating equation `i': " _c
		mata: printf("{res}N = %6.0f     MSE = %10.6f", `eqn'.N, `eqn'.MSE)
		mata: st_local("vtilde",`eqn'.vtilde)
		local minMSE : list vtilde in Yopt
		if `minMSE' {
			di "*" _c
		}
		mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
		if `r(crossfit)'==0 {
			di " (no crossfit)"
		}
		else {
			di
		}
		mata: printf("{res}  Variable: %s{col 30}Orthogonalized: %s\n", `eqn'.vname, `eqn'.vtilde)
		if `showcmd' {
			mata: printf("{res}  Command: %s\n", `eqn'.eststring)
		}
	}
	mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlistD))
	local numeqnsD	= `r(numeqns)'
	di
	di "Number of D estimating equations: `numeqnsD'"
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		di "Estimating equation `i': " _c
		mata: printf("{res}N = %6.0f     MSE = %10.6f", `eqn'.N, `eqn'.MSE)
		mata: st_local("vtilde",`eqn'.vtilde)
		local minMSE : list vtilde in Dopt
		if `minMSE' {
			di "*" _c
		}
		mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
		if `r(crossfit)'==0 {
			di " (no crossfit)"
		}
		else {
			di
		}
		mata: printf("{res}  Variable: %s{col 30}Orthogonalized: %s\n", `eqn'.vname, `eqn'.vtilde)
		if `showcmd' {
			mata: printf("{res}  Command: %s\n", `eqn'.eststring)
		}
	}
	if ("`model'"=="iv") {
		mata: st_numscalar("r(numeqns)",cols(`mname'.eqnlistZ))
		local numeqnsZ	= `r(numeqns)'
		di
		di "Number of Z estimating equations: `numeqnsZ'"
		forvalues i=1/`numeqnsZ' {
			mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
			mata: printf("{res}N = %6.0f     MSE = %10.6f", `eqn'.N, `eqn'.MSE)
			mata: printf("{res}MSE = %10.6f", `eqn'.MSE)
			mata: st_local("vtilde",`eqn'.vtilde)
			local minMSE : list vtilde in Zopt
			if `minMSE' {
				di "*" _c
			}
			mata: st_numscalar("r(crossfit)",`eqn'.crossfit)
			if `r(crossfit)'==0 {
				di " (no crossfit)"
			}
			else {
				di
			}
			mata: printf("{res}  Variable: %s{col 30}Orthogonalized: %s\n", `eqn'.vname, `eqn'.vtilde)
			if `showcmd' {
				mata: printf("{res}  Command: %s\n", `eqn'.eststring)
			}
		}
	}
	
	if "`Yopt'`Dopt'`Zopt'"~="" {
		di
		di "* indicates minimim MSE estimation"
	}

	if `showall' {
		di
		di as res "Other:"
		di
		di as res "liststruct(.):"
		mata: liststruct(`mname')
		di
		di as res "Y equation pointers:"
		mata: `mname'.eqnlistY
		di as res "D equation pointers:"
		mata: `mname'.eqnlistD
		di as res "Z equation pointers:"
		mata: `mname'.eqnlistZ
	}

	// clear this global from Mata
	mata: mata drop `eqn'

end

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end