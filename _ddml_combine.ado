** combines tilde variables from various machine learners to one final estimator using nnls

program _ddml_combine, eclass sortpreserve

	syntax [anything] ,								/// 
							[						///
							mname(name)				///
							]	

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	mata: st_local("model",`mname'.model)
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
	mata: st_local("numeqnsY",strofreal(cols(`mname'.nameYtilde)))
	mata: st_local("numeqnsD",strofreal(cols(`mname'.nameDtilde)))
	mata: st_local("numeqnsDH",strofreal(cols(`mname'.nameDHtilde)))
	mata: st_local("numeqnsZ",strofreal(cols(`mname'.nameZtilde)))

	*** NNLS for Y
	mata: st_local("nameY",`mname'.nameY)
	local listYtilde
	forvalues i=1/`numeqns' {
		mata: `eqn'=*(`mname'.eqnlist[1,`i'])
		mata: st_local("vname",`eqn'.Vname)
		mata: st_local("vtilde",`eqn'.Vtilde)
		mata: st_local("eqntype",`eqn'.eqntype)
		mata: st_local("resid",`eqn'.resid)
		if "`eqntype'"=="yeq" & "`vname'"=="`nameY'" {
			// generate fitted values 
			if (`resid'==1) {
				tempvar vhat`j' 
				gen double `vhat`j''=`dvar'-`vtilde' 
				local listYtilde `listYtilde' `vhat`j''
			}
			else {
				local listYtilde `listYtilde' `vtilde'						
			}
		}
	}
	di "_ddml_nnls `nameY' `listYtilde'"
	// need to create new equation
	mata: add_eqn(`mname', "yeq", "`nameY'", "`gen'","_ddml_nnls `nameY' `listYtilde'","double","`prefix'","","")
	// run NNLS
	_ddml_nnls `nameY' `listYtilde'
	// create residuals
	predict double, r
	// need to calculate MSE etc
	// add to equation
	mata: add_to_eqn()

	*** NNLS for D
	if `numeqnsD' {
		mata: st_local("listD",invtokens(`mname'.nameD))
		local j =1
		foreach dvar of varlist `listD' {
		// obtain all tilde variable specific to D1, D2, ...
			local listDtilde
			forvalues i=1/`numeqns' {
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vname",`eqn'.Vname)
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("eqntype",`eqn'.eqntype)
				mata: st_local("resid",`eqn'.resid)
				if "`eqntype'"=="deq" & "`vname'"=="`dvar'" {
					// generate fitted values 
					if (`resid'==1) {
						tempvar vhat`j' 
						gen double `vhat`j''=`dvar'-`vtilde' 
						local listDtilde `listDtilde' `vhat`j''
					}
					else {
						local listDtilde `listDtilde' `vtilde'						
					}
				}
			}
			di "_ddml_nnls `dvar' `listDtilde'"
			list `dvar' `listDtilde'
			mata: add_eqn(`mname', "deq", "`dvar'", "`dvar't`j'","_ddml_nnls `dvar' `listDtilde'","double","`mname'_","","")
			_ddml_nnls `dvar' `listDtilde'
			// create residuals
			predict double, r
			// need to calculate MSE etc
			// add to equation
			mata: add_to_eqn()
			local j = `j'+1
		}
	}

	*** NNLS for Z
	if `numeqnsZ' {
		mata: st_local("listZ",invtokens(`mname'.nameZ))
		local j =1
		foreach zvar of varlist `listZ' {
			local listZtilde
			forvalues i=1/`numeqns' {
				mata: `eqn'=*(`mname'.eqnlist[1,`i'])
				mata: st_local("vname",`eqn'.Vname)
				mata: st_local("vtilde",`eqn'.Vtilde)
				mata: st_local("eqntype",`eqn'.eqntype)
				mata: st_local("resid",`eqn'.resid)
				if "`eqntype'"=="zeq" & "`vname'"=="`zvar'" {
					if (`resid'==1) {
						tempvar vhat`j' 
						gen double `vhat`j''=`zvar'-`vtilde' 
						local listZtilde `listZtilde' `vhat`j''
					}
					else {
						local listZtilde `listZtilde' `vtilde'						
					}
				}
			}
			di "_ddml_nnls `zvar' `listZtilde'"
			list `zvar' `listZtilde'
			mata: add_eqn(`mname', "zeq", "`zvar'", "`gen'","_ddml_nnls `zvar' `listZtilde'","double","`prefix'","","")
			_ddml_nnls `zvar' `listZtilde'
			// create residuals
			predict double, r
			// need to calculate MSE etc
			// add to equation
			mata: add_to_eqn()
			local j = `j'+1
		}
	}
end

mata:
struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}


void add_to_eqn(					struct ddmlStruct m,
									real scalar eqnumber, 
									string scalar vtilde)

{
	pointer(struct eqnStruct) scalar p

	if vtilde!="" {
		eqnumber = selectindex(m.eqnlistNames:==vtilde)
	}

	cmd 			= st_global("r(cmd)")
	mse				= st_numscalar("r(mse)")
	mse_folds		= st_matrix("r(mse_folds)")
	n				= st_numscalar("r(N)")
	n_folds			= st_matrix("r(N_folds)")
	p				= m.eqnlist[1,eqnumber]
	(*p).MSE		= ((*p).MSE \ mse)
	(*p).N			= ((*p).N \ n)
	(*p).command	= cmd

	if (cmd == "pystacked") {
		(*p).stack_weights = st_matrix("r(pysw)")		 
	}

	// MSE by fold list should be initialized to void 0-by-k matrix
	// (otherwise concat fails because of conformability)
	(*p).MSE_folds	= ((*p).MSE_folds \ mse_folds)
	(*p).N_folds	= ((*p).N_folds \ n_folds)
	
	// set crossfitted flag = 1
	(*p).crossfitted	= 1

}

end