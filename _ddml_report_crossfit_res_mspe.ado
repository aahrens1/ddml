program define _ddml_report_crossfit_res_mspe
	syntax name(name=mname), etype(string) [ vlist(string) zett(string) m(integer 1) ]

	// set struct field name
	if "`etype'"=="yeq" {
		if "`zett'"=="" {
			local optname nameYopt
			local estring y|X
		}
		else if "`zett'"=="0" {
			local optname nameY0opt
			local estring y|X,D=0
		}
		else if "`zett'"=="1" {
			local optname nameY1opt
			local estring y|X,D=1
		}
	}
	else if "`etype'"=="deq" {
		if "`zett'"=="" {
			local optname nameDopt
			local estring D|X
		}
		else if "`zett'"=="0" {
			local optname nameD0opt
			local estring D|X,Z=0
		}
		else if "`zett'"=="1" {
			local optname nameD1opt
			local estring D|X,Z=1
		}
	}
	else if "`etype'"=="dheq" {
		local optname nameDHopt
		local estring D|X,Z
	}
	else if "`etype'"=="zeq" {
		local optname nameZopt
		local estring Z|X
	}
	else {
		di as err "internal error - unknown equation type `etype'"
		exit 198
	}

	// may be called with empty list (e.g. if no endog regressors)
	local numeqns	: word count `vlist'
	if `numeqns' > 0 {
		
		di
		di as res "Mean-squared error for `estring':"
		di _col(2) "Name" _c
		di _col(20) "Orthogonalized" _c
		di _col(40) "Command" _c
		di _col(54) "N" _c
		di _col(65) "MSPE"
		di "{hline 75}"
		// clear opt list
		// mata: `mname'.`optname' = J(1,0,"")
		foreach var of varlist `vlist' {
			// m is the rep number
			display_mspe `mname', vname(`var') etype(`etype') m(`m') zett(`zett')
			local optlist `optlist' `r(optname)'
		}

		if `m'==1 {
			mata: `mname'.`optname' = tokens("`optlist'")
		}
		else {
			mata: `mname'.`optname' = (`mname'.`optname' \ tokens("`optlist'"))
		}

	}
end


* with resampling can N vary? currently assumed to be constant and only MSPE varies
program define display_mspe, rclass

	syntax name(name=mname), vname(varname)		///
								[				///
								etype(string)	///
								zett(string)	///
								m(integer 1)	/// resample number
								]

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))

	// initialize
	local minmse = .
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_global("r(vname)",`eqn'.Vname)
		mata: st_global("r(eqntype)",`eqn'.eqntype)
		mata: st_local("vtilde",`eqn'.Vtilde)
		if "`vname'"==r(vname) & "`etype'"==r(eqntype) {
			mata: st_local("command",`eqn'.command)
			mata: st_local("MSE",strofreal(`eqn'.MSE`zett'[`m']))
			mata: st_local("N",strofreal(`eqn'.N`zett'))
			if `MSE' < `minmse' {
				local optname `vtilde'
				local minmse `MSE'
			}
			di _col(2) "`vname'" _c
			di _col(20) "`vtilde'" _c
			di _col(40) "`command'" _c
			di _col(50) %6.0f `N' _c
			di _col(60) %10.6f `MSE'
		}
	}
	
	return local optname	`optname'
	return scalar minmse	=`minmse'
	
	mata: mata drop `eqn'
	
end

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end
