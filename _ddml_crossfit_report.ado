program define _ddml_crossfit_report
	syntax name(name=mname)

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listD",invtokens(`mname'.nameD))
	mata: st_local("listDH",invtokens(`mname'.nameDH))
	mata: st_local("listZ",invtokens(`mname'.nameZ))
	mata: st_local("reps",strofreal(cols(`mname'.idFold)))
	// first column is ID number so subtract
	local reps = `reps'-1

	if "`model'"=="interactive" | "`model'"=="late" {
		// y equation - all models
		forvalues m=1/`reps' {
			report_equation `mname', etype(yeq) vlist(`nameY') m(`m') model(`model') zett(0)
		}
		forvalues m=1/`reps' {
			report_equation `mname', etype(yeq) vlist(`nameY') m(`m') model(`model') zett(1)
		}
		// interactive model - D eqn only
		if "`model'"=="interactive" {
			forvalues m=1/`reps' {
				report_equation `mname', etype(deq) vlist(`listD') m(`m') model(`model')
			}
		}
		// LATE - D and Z eqns
		else if "`model'"=="late" {
			forvalues m=1/`reps' {
				report_equation `mname', etype(deq) vlist(`listD') m(`m') model(`model') zett(0)
			}
			forvalues m=1/`reps' {
				report_equation `mname', etype(deq) vlist(`listD') m(`m') model(`model') zett(1)
			}
			forvalues m=1/`reps' {
				report_equation `mname', etype(zeq) vlist(`listZ') m(`m') model(`model')
			}
		}
	}		// end interactive/LATE block
	else {
		forvalues m=1/`reps' {
			report_equation `mname', etype(yeq) vlist(`nameY') m(`m') model(`model')
		}
		// use foreach ... in ... to accommodate empty lists
		foreach var in `listD' {
			forvalues m=1/`reps' {
				report_equation `mname', etype(deq) vlist(`var') m(`m') model(`model')
			}
		}
		foreach var in `listDH' {
			forvalues m=1/`reps' {
				report_equation `mname', etype(dheq) vlist(`var') m(`m') model(`model')
			}
		}
		foreach var in `listZ' {
			forvalues m=1/`reps' {
				report_equation `mname', etype(zeq) vlist(`var') m(`m') model(`model')
			}
		}
	}		// end partial linear/IV block

end

program define report_equation
	syntax name(name=mname), etype(string) [ vlist(string) zett(string) m(integer 1) model(string) ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	mata: st_local("model",`mname'.model)

	// set struct field name
	if "`etype'"=="yeq" {
		if "`zett'"=="" {
			local optname nameYopt
			local estring E[y|X]
		}
		else {
		// Either interactive (D=0,1) or LATE (Z=0,1)
			if "`model'"=="interactive" {
				local DZ D
			}
			else {
				local DZ Z
			}
			if "`zett'"=="0" {
				local optname nameY0opt
				local estring E[y|X,`DZ'=0]
			}
			else if "`zett'"=="1" {
				local optname nameY1opt
				local estring E[y|X,`DZ'=1]
			}
		}
	}
	else if "`etype'"=="deq" {
		if "`zett'"=="" {
			local optname nameDopt
			local estring E[D|X]
		}
		else if "`zett'"=="0" {
			local optname nameD0opt
			local estring E[D|X,Z=0]
		}
		else if "`zett'"=="1" {
			local optname nameD1opt
			local estring E[D|X,Z=1]
		}
	}
	else if "`etype'"=="dheq" {
		local optname nameDHopt
		if "`model'"=="optimaliv" {
			local estring E[D^|X,Z] with D^=E^[D|X]
		} 
		else {
			local estring E[D|X,Z]
		}
	}
	else if "`etype'"=="zeq" {
		local optname nameZopt
		local estring E[Z|X]
	}
	else {
		di as err "internal error - unknown equation type `etype'"
		exit 198
	}

	// may be called with empty list (e.g. if no endog regressors)
	local numeqns	: word count `vlist'
	if `numeqns' > 0 {
		if `m'==1 {
			di
			di as res "Mean-squared error for `estring':"
			di _col(2) "Name" _c
			di _col(20) "Orthogonalized" _c
			di _col(40) "Command" _c
			di _col(54) "N" _c
			di _col(65) "MSPE" _c
			di _col(75) "m"
			di "{hline 75}"
		}
		foreach var of varlist `vlist' {
			// m is the rep number
			display_mspe `mname', vname(`var') etype(`etype') m(`m') zett(`zett') model(`model') optname(`optname')
		}
		di "{hline 75}"
	}
	
	mata: mata drop `eqn'

end


* with resampling can N vary? currently assumed to be constant and only MSPE varies
program define display_mspe, rclass

	syntax name(name=mname), vname(varname)		///
								[				///
								etype(string)	///
								zett(string)	///
								m(integer 1)	/// resample number
								model(string) 	///
								optname(string)	///
								]

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("optlist",invtokens(`mname'.`optname'[`m',.]))

	** special case optimaliv: need to take MSE & N from deq structure
	local vtildeh
	if "`model'"=="optimaliv"&"`etype'"=="dheq" {
		local etype deq
		local zett _h 
		local vtildeh _h
	}

	// initialize
	local minmse = .
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_global("r(vname)",`eqn'.Vname)
		mata: st_global("r(eqntype)",`eqn'.eqntype)
		mata: st_local("vtilde",`eqn'.Vtilde`vtildeh')
		local optflag : list optlist & vtilde
		if "`optflag'"~="" {
			local optflag *
		}

		if "`vname'"==r(vname) & "`etype'"==r(eqntype) {
			mata: st_local("command",`eqn'.command)
			mata: st_local("MSE",strofreal(`eqn'.MSE`zett'[`m']))
			mata: st_local("N",strofreal(`eqn'.N`zett'))
			di _col(2) "`vname'" _c
			di _col(20) "`vtilde'" _c
			di _col(40) "`command'" _c
			di _col(50) %6.0f `N' _c
			di _col(60) %10.4f `MSE' "`optflag'" _c
			di _col(74) %2.0f `m'
		}
	}
	
	mata: mata drop `eqn'
	
end

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end
