// to add - option to save the raw data with varnames in `mname'.strDatavars

program _ddml_export
	version 13

	syntax , mname(name) fname(string) [ * ]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()
	
	*** extract details of estimation
	
	// model
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlist)))
	
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_local("vname",`eqn'.Vname)
		mata: st_local("vtilde",`eqn'.Vtilde)
		local vlist		`vlist' `vname'
		local tildelist	`tildelist' `vtilde'
	}
	
	preserve
	
	// remove equation name prefix from tilde variable names
	foreach vname in `tildelist' {
		local newvname : subinstr local vname "`mname'_" ""
		rename `vname' `newvname'
		local newtildelist `newtildelist' `newvname'
	}

	cap drop id
	qui gen id = strofreal(`mname'_id, "%20.0g")
	cap drop fid
	qui gen fid = strofreal(`mname'_fid, "%20.0g")
	keep id fid `newtildelist'
	order id fid, first

	// data
	export delimited using `fname', `options'

	// dictionary
	drop _all
	qui set obs 1
	local i 1
	tokenize `vlist'
	foreach vname in `newtildelist' {
		qui gen `vname' = "``i''"
		local ++i
	}
	export delimited using `fname'_dict, `options'
	
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
