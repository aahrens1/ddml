
program _ddml_export
	version 13

	syntax , mname(name) fname(string) [ * ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	
	*** extract details of estimation
	
	// model
	mata: st_local("numeqnsY",strofreal(cols(`mname'.eqnlistY)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.eqnlistD)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.eqnlistZ)))
	
	forvalues i=1/`numeqnsY' {
		mata: `eqn'=*(`mname'.eqnlistY[1,`i'])
		mata: st_local("vname",`eqn'.Vname)
		mata: st_local("vtilde",`eqn'.Vtilde)
		local vlist		`vlist' `vname'
		local tildelist	`tildelist' `vtilde'
	}
	forvalues i=1/`numeqnsD' {
		mata: `eqn'=*(`mname'.eqnlistD[1,`i'])
		mata: st_local("vname",`eqn'.Vname)
		mata: st_local("vtilde",`eqn'.Vtilde)
		local vlist		`vlist' `vname'
		local tildelist	`tildelist' `vtilde'
	}
	forvalues i=1/`numeqnsZ' {
		mata: `eqn'=*(`mname'.eqnlistZ[1,`i'])
		mata: st_local("vname",`eqn'.Vname)
		mata: st_local("vtilde",`eqn'.Vtilde)
		local vlist		`vlist' `vname'
		local tildelist	`tildelist' `vtilde'
	}
	
	preserve

	foreach vname in `tildelist' {
		cap drop `vname'
		gen `vname' = strofreal(`mname'_`vname',"%26.23e")
	}

	cap drop id
	qui gen id = strofreal(`mname'_id, "%20.0g")
	cap drop fid
	qui gen fid = strofreal(`mname'_fid, "%20.0g")
	keep id fid `tildelist'

	tempfile tfile
	qui save `tfile'
	drop _all
	qui set obs 1
	qui gen id = "id"
	qui gen fid = "fid"
	local i 1
	tokenize `vlist'
	foreach vname in `tildelist' {
		qui gen `vname' = "``i''"
		local ++i
	}
	qui append using `tfile'

	export delimited using `fname', `options'
	
	restore	

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

struct ddmlStruct use_model(		string scalar fname)
{
	struct ddmlStruct scalar	m
	
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
