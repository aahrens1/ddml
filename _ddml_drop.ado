
program _ddml_drop, eclass
	version 13

	syntax , mname(name)

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	mata: st_local("numeqnsY",strofreal(cols(`mname'.eqnlistY)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.eqnlistD)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.eqnlistZ)))

	*** drop id and fold id
	cap drop `mname'_id
	cap drop `mname'_fid

	*** loop through equations and drop Stata variables
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		cap drop `mname'_`vtilde'
	}
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		cap drop `mname'_`vtilde'
	}
	if ("`model'"=="iv") {
		forvalues i=1/`numeqnsZ' {
			mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			cap drop `mname'_`vtilde'
		}
	}
	
	mata: mata drop `mname'

end

********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

end
