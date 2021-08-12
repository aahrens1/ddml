* notes:


program define _ddml_describe

	syntax name(name=mname), [NOCMD all byfold]
	
	local showcmd	= ("`nocmd'"=="")
	local showall	= ("`all'"~="")

	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))

	mata: printf("{res}Model: %s\n", `mname'.model)
	di as res "ID: `mname'_id"
	di as res "Fold ID: `mname'_fid"
	di as res "Sample indicator: `mname'_sample" _c
	qui count if `mname'_sample
	di as res " (N=`r(N)')"

	mata: printf("{res}Dependent variable (Y): %s\n", `mname'.nameY)
	mata: printf("{res}Dependent variable (orthogonalized): %s\n", invtokens(`mname'.nameYtilde))
	mata: printf("{res}Causal variable(s) (D): %s\n", invtokens(`mname'.nameD))
	mata: printf("{res}Causal variable(s) (orthogonalized): %s\n", invtokens(`mname'.nameDtilde))
	if ("`model'"=="iv") {
		mata: printf("{res}Excluded instrumental variable(s): %s\n", `mname'.nameZ)
		mata: printf("{res}Excluded instrumental variable(s) (orthogonalized): %s\n", `mname'.nameZtilde)
	}
	
	// List equations

	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlist)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsDH",strofreal(cols(`mname'.nameDHtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))
	if `numeqnsY'>0 {
		di
		di "Number of Y estimating equations: `numeqnsY'"
		desc_equation `mname', eqntype(yeq) optlist(`Yopt') showcmd(`showcmd')
	}
	if `numeqnsD'>0 {
		di
		di "Number of D estimating equations: `numeqnsD'"
		desc_equation `mname', eqntype(deq) optlist(`Dopt') showcmd(`showcmd')
	}
	if `numeqnsZ'>0 {
		di
		di "Number of Z estimating equations: `numeqnsZ'"
		desc_equation `mname', eqntype(zeq) optlist(`Zopt') showcmd(`showcmd')
	}
	if `numeqnsDH'>0 {
		di
		di "Number of DH estimating equations: `numeqnsDH'"
		desc_equation `mname', eqntype(dheq) optlist(`DHopt') showcmd(`showcmd')
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
	
	if `crossfitted' {
		_ddml_crossfit_report `mname', `byfold'
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
			di "Estimating equation `i': "
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