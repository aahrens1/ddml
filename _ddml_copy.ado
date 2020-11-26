
program _ddml_copy
	version 13

	syntax , mname(name) newmname(name)

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	tempfile fname
	mata: save_model("`fname'",`mname')

	mata: `newmname' = use_model("`fname'")

	mata: st_local("numeqnsY",strofreal(cols(`mname'.eqnlistY)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.eqnlistD)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.eqnlistZ)))

	*** create id and fold id
	cap drop `mname'_id
	cap drop `mname'_fid
	mata: st_numscalar("r(nobs)",rows(`mname'.id))
	if r(nobs) > _N {
		set obs `r(nobs)'
	}
	qui gen double `newmname'_id = .
	qui gen double `newmname'_fid = .
	mata: st_store( ., ("`newmname'_id", "`newmname'_fid"), (`newmname'.idFold))

	*** loop through equations and create Stata variables
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`newmname'.eqnlistY[1,`i'])
		mata: st_local("vtilde",`eqn'.vtilde)
		cap drop `newmname'_`vtilde'
		qui gen double `newmname'_`vtilde' = .
		mata: st_store( ., ("`newmname'_`vtilde'"), (`eqn'.idVtilde)[.,2])
	}
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`newmname'.eqnlistD[1,`i'])
		mata: st_local("vtilde",`eqn'.vtilde)
		cap drop `newmname'_`vtilde'
		qui gen double `newmname'_`vtilde' = .
		mata: st_store( ., ("`newmname'_`vtilde'"), (`eqn'.idVtilde)[.,2])
	}
	if ("`model'"=="iv") {
		forvalues i=1/`numeqnsZ' {
			mata: `eqn'=*(`newmname'.eqnlistZ[1,`i'])
			mata: st_local("vtilde",`eqn'.vtilde)
			cap drop `newmname'_`vtilde'
			qui gen double `newmname'_`vtilde' = .
			mata: st_store( ., ("`newmname'_`vtilde'"), (`eqn'.idVtilde)[.,2])
		}
	}


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

void save_model(					string scalar fname,
									struct ddmlStruct m)
{
	fh = fopen(fname,"w")
	fputmatrix(fh,m)
	fclose(fh)
}

struct ddmlStruct use_model(		string scalar fname)
{
	struct ddmlStruct scalar	m
	fh = fopen(fname,"r")
	m = fgetmatrix(fh,1)	// nonzero second argument required for "strict"
	fclose(fh)
	return(m)
}

end