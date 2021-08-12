program define _ddml_crossfit_report
	syntax name(name=mname), [byfold]

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	mata: st_local("model",`mname'.model)
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listD",invtokens(`mname'.nameD))
	mata: st_local("listDH",invtokens(`mname'.nameDH))
	mata: st_local("listZ",invtokens(`mname'.nameZ))
	mata: st_local("numreps",strofreal(`mname'.nreps))

	if "`model'"=="interactive" | "`model'"=="late" {
		// y equation - all models
		forvalues m=1/`numreps' {
			report_equation `mname', etype(yeq) vlist(`nameY') m(`m') model(`model') zett(0) `byfold'
		}
		forvalues m=1/`numreps' {
			report_equation `mname', etype(yeq) vlist(`nameY') m(`m') model(`model') zett(1) `byfold'
		}
		// interactive model - D eqn only
		if "`model'"=="interactive" {
			forvalues m=1/`numreps' {
				report_equation `mname', etype(deq) vlist(`listD') m(`m') model(`model') `byfold'
			}
		}
		// LATE - D and Z eqns
		else if "`model'"=="late" {
			forvalues m=1/`numreps' {
				report_equation `mname', etype(deq) vlist(`listD') m(`m') model(`model') zett(0) `byfold'
			}
			forvalues m=1/`numreps' {
				report_equation `mname', etype(deq) vlist(`listD') m(`m') model(`model') zett(1) `byfold'
			}
			forvalues m=1/`numreps' {
				report_equation `mname', etype(zeq) vlist(`listZ') m(`m') model(`model') `byfold'
			}
		}
	}		// end interactive/LATE block
	else {
		forvalues m=1/`numreps' {
			report_equation `mname', etype(yeq) vlist(`nameY') m(`m') model(`model') `byfold'
		}
		// use foreach ... in ... to accommodate empty lists
		foreach var in `listD' {
			if "`model'"=="optimaliv" {
				local fitted fitted
			}
			forvalues m=1/`numreps' {
				report_equation `mname', etype(deq) vlist(`var') m(`m') model(`model') `byfold' `fitted'
			}
		}
		foreach var in `listDH' {
			if "`model'"=="optimaliv" {
				local fitted fitted
			}
			forvalues m=1/`numreps' {
				report_equation `mname', etype(dheq) vlist(`var') m(`m') model(`model') `byfold' `fitted'
			}
		}
		foreach var in `listZ' {
			forvalues m=1/`numreps' {
				report_equation `mname', etype(zeq) vlist(`var') m(`m') model(`model') `byfold'
			}
		}
	}		// end partial linear/IV block

end

program define report_equation
	syntax name(name=mname), etype(string) [ vlist(string) zett(string) m(integer 1) model(string) byfold fitted ]
	
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
			display_mspe `mname', vname(`var') etype(`etype') m(`m') zett(`zett') model(`model') optname(`optname') `byfold' `fitted'
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
								byfold			///
								fitted			/// is vtilde a fitted value? (default=residual)
								]

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("optlist",invtokens(`mname'.`optname'[`m',.]))

	if "`byfold'"~="" {
		local foldvar `mname'_fid_`m'
		qui sum `foldvar', meanonly
		local kfolds = r(max)
	}

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
			if "`byfold'"~="" {
				tempvar spe
				if "`fitted'"=="" {
					qui gen double `spe' = `vtilde'_`m'^2
				}
				else {
					qui gen double `spe' = (`vname'-`vtilde'_`m')^2
				}
				if "`model'"=="interactive" & "`etype'"=="yeq" {
					mata: st_local("Dvar",`mname'.nameD)
					local and_zett & `Dvar'==`zett'
				}
				else if "`model'"=="late" & ("`etype'"=="yeq" | "`etype'"=="deq") {
					mata: st_local("Zvar",`mname'.nameZ)
					local and_zett & `Zvar'==`zett'
				}
				forvalues j=1/`kfolds' {
					di _col(4) "fold=" _col(7) %2.0f `j' _c
					qui count if `j'==`foldvar' `and_zett'
					di _col(50) %6.0f `r(N)' _c
					qui sum `spe' if `j'==`foldvar' `and_zett', meanonly
					di _col(60) %10.4f `r(mean)'
				}
			}
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
