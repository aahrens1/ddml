*! ddml v1.1
*! last edited: 28 dec 2022
*! authors: aa/ms

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

	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlist)))

	*** create id and fold id
	cap drop `newmname'_id
	cap drop `newmname'_fid
	mata: st_numscalar("r(nobs)",rows(`mname'.id))
	if r(nobs) > _N {
		set obs `r(nobs)'
	}
	qui gen double `newmname'_id = .
	qui gen byte `newmname'_sample = .
	mata: st_store( ., ("`newmname'_id", "`newmname'_sample"), (`newmname'.idSample))
	// id and sample variables always exist, fold ID may not
	mata: st_numscalar("r(ncols)",cols(`newmname'.idFold))
	if r(ncols) > 0 {
		qui gen double `newmname'_fid = .
		mata: st_store( ., ("`newmname'_fid"), (`newmname'.idFold[.,2]))
	}

	*** loop through equations and create Stata variables
	// note that variables may not exist
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`newmname'.eqnlist[1,`i'])
		mata: st_local("vtilde",`eqn'.Vtilde)
		cap drop `newmname'_`vtilde'
		mata: st_numscalar("r(ncols)",cols(`eqn'.idVtilde))
		if r(ncols) > 0 {
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
