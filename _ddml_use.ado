*! ddml v1.4
*! last edited: 25july2023
*! authors: aa/ms

program _ddml_use
	version 13

	syntax , mname(name) fname(string) [ replace ]

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eStruct()

	// does model already exist?
	mata: st_local("isnull",strofreal(findexternal("`mname'")==NULL))

	if `isnull' | "`replace'"=="replace" {
		mata: `mname' = use_model("`fname'")
	}
	else {
		di as err "error - `mname' already exists in Mata memory; use -replace- option"
		exit 198
	}

	/*
	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	di "Model: `model'"
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlist)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("listYtilde",invtokens(`mname'.nameYtilde))
	di "Number of Y estimating equations: `numeqnsY'"
	if `numeqnsD' {
		mata: st_local("listD",invtokens(`mname'.nameD))
		mata: st_local("listDtilde",invtokens(`mname'.nameDtilde))
		di "Number of D estimating equations: `numeqnsD'"
	}
	if `numeqnsZ' {
		mata: st_local("listZ",invtokens(`mname'.nameZ))
		mata: st_local("listZtilde",invtokens(`mname'.nameZtilde))
		di "Number of Z estimating equations: `numeqnsZ'"
	}
	
	*** create id and fold id
	cap drop `mname'_id
	cap drop `mname'_fid
	mata: st_numscalar("r(nobs)",rows(`mname'.id))
	if r(nobs) > _N {
		set obs `r(nobs)'
	}
	qui gen double `mname'_id = .
	// id variable always exists, fold ID may not
	mata: st_numscalar("r(ncols)",cols(`mname'.idFold))
	if r(ncols) > 0 {
		qui gen double `mname'_fid = .
		mata: st_store( ., ("`mname'_id", "`mname'_fid"), (`mname'.idFold))
	}
	else {
		mata: st_store( ., ("`mname'_id"), (`mname'.id))
	}

	*** loop through equations and create Stata variables
	// note that variables may not exist
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		cap drop `mname'_`vtilde'
		mata: st_numscalar("r(ncols)",cols(`eqn'.idVtilde))
		if r(ncols) > 0 {
			qui gen double `mname'_`vtilde' = .
			mata: st_store( ., ("`mname'_`vtilde'"), (`eqn'.idVtilde)[.,2])
		}
	}
	*/

end

********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

struct mStruct use_model(string scalar fname)
{
	struct mStruct scalar	m
	
	if (!fileexists(fname)) {
		errprintf("file %s not found\n", fname)
		exit(601)
	}
	
	fh = fopen(fname,"r")
	m = fgetmatrix(fh,1)	// nonzero second argument required for "strict"
	fclose(fh)
	return(m)
}

end
