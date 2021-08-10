program define _ddml_crossfit_update_optlist
	syntax name(name=mname), etype(string) [ vlist(string) zett(string) m(integer 1) model(string) ]

	// may be called with empty list (e.g. if no endog regressors)
	if "`vlist'"=="" exit

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	// set struct field name
	if "`etype'"=="yeq" {
		if "`zett'"=="" {
			local optname nameYopt
		}
		else if "`zett'"=="0" {
			local optname nameY0opt
		}
		else if "`zett'"=="1" {
			local optname nameY1opt
		}
	}
	else if "`etype'"=="deq" {
		if "`zett'"=="" {
			local optname nameDopt
		}
		else if "`zett'"=="0" {
			local optname nameD0opt
		}
		else if "`zett'"=="1" {
			local optname nameD1opt
		}
	}
	else if "`etype'"=="dheq" {
		local optname nameDHopt
		** special case optimaliv: need to take MSE & N from deq structure
		local vtildeh
		if "`model'"=="optimaliv" {
			local etype deq
			local zett _h 
			local vtildeh _h
		}
	}
	else if "`etype'"=="zeq" {
		local optname nameZopt
	}
	else {
		di as err "internal error - unknown equation type `etype'"
		exit 198
	}

	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	foreach var of varlist `vlist' {
		// initialize
		local minmse = .
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_global("r(vname)",`eqn'.Vname)
			mata: st_global("r(eqntype)",`eqn'.eqntype)
			mata: st_local("vtilde",`eqn'.Vtilde`vtildeh')
			if "`var'"==r(vname) & "`etype'"==r(eqntype) {
				mata: st_local("command",`eqn'.command)
				// m is the rep number
				mata: st_local("MSE",strofreal(`eqn'.MSE`zett'[`m']))
				mata: st_local("N",strofreal(`eqn'.N`zett'))
				if `MSE' < `minmse' {
					local optvar `vtilde'
					local minmse `MSE'
				}
			}
		}
		local optlist `optlist' `optvar'
	}

	// assumes opt lists have been initialized as void lists
	mata: `mname'.`optname' = (`mname'.`optname' \ tokens("`optlist'"))
	
	mata: mata drop `eqn'

end


mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end
