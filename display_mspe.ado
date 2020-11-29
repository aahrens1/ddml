
program define display_mspe, rclass

	syntax name(name=mname), vname(varname) ///
								[ ///
								zett(string) ///
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
		mata: st_local("vtilde",`eqn'.Vtilde)
		if "`vname'"==r(vname) {
			mata: st_local("command",`eqn'.command)
			mata: st_local("MSE",strofreal(`eqn'.MSE`zett'))
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