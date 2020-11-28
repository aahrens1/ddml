
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
	di as res "Sample indicator: `mname'_esample" _c
	cap count if `mname'_esample
	if _rc==0 {
		di as res " (N=`r(N)')"
	} 
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

	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlist)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))
	di
	di "Number of Y estimating equations: `numeqnsY'"
	desc_equation `mname', eqntype(yeq) optlist(`Yopt') showcmd(`showcmd')
	di
	di "Number of D estimating equations: `numeqnsD'"
	desc_equation `mname', eqntype(deq) optlist(`Dopt') showcmd(`showcmd')
	if ("`model'"=="iv") {
		di
		di "Number of Z estimating equations: `numeqnsZ'"
		desc_equation `mname', eqntype(zeq) optlist(`Zopt') showcmd(`showcmd')
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
		di as res "Equation pointers:"
		mata: `mname'.eqnlist
		di as res "Corresponding tilde names:"
		mata: `mname'.eqnlistNames
	}

end

prog define desc_equation

	syntax name(name=mname), eqntype(string) [ optlist(string) showcmd(integer 0) ]

	tempname eqn
	mata: `eqn' = init_eqnStruct()

	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlist)))

	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_global("r(eqntype)",`eqn'.eqntype)
		if "`eqntype'"==r(eqntype) {
			di "Estimating equation `i': " _c
			mata: printf("{res}N = %6.0f     MSE = %10.6f", `eqn'.N, `eqn'.MSE)
			mata: st_local("vtilde",`eqn'.Vtilde)
			local minMSE : list vtilde in optlist
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
			mata: printf("{res}  Variable: %s{col 30}Orthogonalized: %s\n", `eqn'.Vname, `eqn'.Vtilde)
			if `showcmd' {
				mata: printf("{res}  Command: %s\n", `eqn'.eststring)
			}
		}
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