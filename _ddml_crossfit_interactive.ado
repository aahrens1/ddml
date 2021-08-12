* last edited: 18 jun 2021

* notes:
* why is crossfitted field set in additive code but not here?
* check it's correct that interactive-type estimation always goes with reporting mse0 and mse1
* are multiple Ds allowed?

*** ddml cross-fitting for the interactive model & LATE
program _ddml_crossfit_interactive, eclass sortpreserve

	syntax [anything] [if] [in] , /// 
							[ kfolds(integer 5)		///
							NOIsily					///
							debug					/// 
							Robust					///
							TABFold					///
							foldlist(numlist)		///
							mname(name)				///
							reps(integer 0)			/// if reps specified then overwrite any existing fold vars
							]

	// no checks included yet
	// no marksample yet

	local debugflag		= "`debug'"~=""
	if ("`noisily'"=="") local qui qui
		
	*** extract details of estimation
	
	// model
	mata: st_local("model",`mname'.model)
	di "Model: `model'"
	mata: st_local("numeqns",strofreal(cols(`mname'.eqnlistNames)))
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
	// initialize opt lists with void matrices of correct dimension
	mata: st_local("numvarsD",strofreal(cols(`mname'.nameD)))
	mata: st_local("numvarsZ",strofreal(cols(`mname'.nameZ)))
	mata: `mname'.nameY0opt = J(0,1,"")
	mata: `mname'.nameY1opt = J(0,1,"")
	mata: `mname'.nameDopt = J(0,`numvarsD',"")
	mata: `mname'.nameD0opt = J(0,`numvarsD',"")
	mata: `mname'.nameD1opt = J(0,`numvarsD',"")
	mata: `mname'.nameZopt = J(0,`numvarsZ',"")
	
	// folds and fold IDs
	mata: st_local("hasfoldvars",strofreal(cols(`mname'.idFold)))
	// if reps is specified then overwrite any existing fold vars.
	if `reps'>0 {
		local hasfoldvars = 0
	}
	else {
		// default reps = 1
		local reps = 1
	}
	// store number of resamplings
	mata: `mname'.nreps = `reps'

	// if empty:
	// add fold IDs to model struct (col 1 = id, col 2 = fold id 1, col 3 = fold id 2 etc.)
	// first initialize with id
	
	if `hasfoldvars'==0 {
		mata: `mname'.idFold = st_data(., ("`mname'_id"))
	}

	forvalues m=1/`reps' {
	
		if `hasfoldvars'==0 {
			*** gen folds
			cap drop `mname'_fid_`m'
			tempvar uni cuni
			qui gen double `uni' = runiform() if `mname'_sample
			qui cumul `uni' if `mname'_sample, gen(`cuni')
			qui gen int `mname'_fid_`m' = ceil(`kfolds'*`cuni') if `mname'_sample
			// add fold id to model struct (col 1 = id, col 2 = fold id)
			mata: `mname'.idFold = (`mname'.idFold , st_data(., ("`mname'_fid_`m'")))
		}
		if ("`tabfold'"!="") {
			di
			di "Overview of frequencies by fold (sample `m'):"
			tab `mname'_fid_`m' if `mname'_sample
			di
		}
	
	}

	// blank eqn - declare this way so that it's a struct and not transmorphic
	// used multiple times below
	tempname eqn
	mata: `eqn' = init_eqnStruct()

	forvalues m=1/`reps' {
	
		di
		di as text "Starting cross-fitting (sample = `m')"

		*** initialize tilde variables
		forvalues i=1/`numeqns' {
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			/// macro m is resampling counter
			cap drop `vtilde'_`m'
			qui gen double `vtilde'_`m'=.
		}
			
		*** do cross-fitting
		
		forvalues i=1/`numeqns' {
		
			// initialize prior to calling crossfit
			mata: `eqn'=*(`mname'.eqnlist[1,`i'])
			mata: st_local("vtilde",`eqn'.Vtilde)
			mata: st_local("vname",`eqn'.Vname)
			mata: st_local("eststring",`eqn'.eststring)
			mata: st_local("eqntype",`eqn'.eqntype)
			// seems to be unused
			// mata: st_local("vtype",`eqn'.vtype)
			local touse `mname'_sample
			// always request fitted values
			local resid  
	
			di as text "Cross-fitting equation `i' (`vname', `vtilde')" _c
	
			/* why not used here but used in interactive code?
			// has the equation already been crossfitted?
			mata: st_numscalar("cvdone",`eqn'.crossfitted)
			if ("`cvdone'"=="1") continue
			*/
			
			if ("`eqntype'"=="yeq") & ("`model'"=="interactive") {
				local treatvar	`listD'
			}
			else if ("`eqntype'"=="yeq") & ("`model'"=="late") {
				local treatvar	`listZ'
			}
			else if ("`eqntype'"=="deq") & ("`model'"=="interactive") {
				local treatvar
			}
			else if ("`eqntype'"=="deq") & ("`model'"=="late") {
				local treatvar	`listZ'
			}
			else if ("`eqntype'"=="zeq") {
				local treatvar
			}
			else {
				di as err "Unknown equation type `eqntype'"
				exit 198
			}
	
			crossfit if `touse',					///
				eststring(`eststring')				///
				kfolds(`kfolds')					///
				foldvar(`mname'_fid_`m')			///
				vtilde(`vtilde'_`m')				///
				vname(`vname')						///
				treatvar(`treatvar')				///
				`resid'
			
			// store MSE and sample size
			if ("`eqntype'"=="yeq") {
				mata: add_to_eqn01(`mname',`i',"`mname'_id `vtilde'", `r(mse0)',`r(N0)',0)
				mata: add_to_eqn01(`mname',`i',"`mname'_id `vtilde'", `r(mse1)',`r(N1)',1)
			}
			else if ("`eqntype'"=="deq") & ("`model'"=="interactive") {
				mata: add_to_eqn(`mname',`i',"`mname'_id `vtilde'", `r(mse)',`r(N)')
			}
			else if ("`eqntype'"=="deq") {
				mata: add_to_eqn01(`mname',`i',"`mname'_id `vtilde'", `r(mse0)',`r(N0)',0)
				mata: add_to_eqn01(`mname',`i',"`mname'_id `vtilde'", `r(mse1)',`r(N1)',1)
			}
			else if ("`eqntype'"=="zeq") {
				mata: add_to_eqn(`mname',`i',"`mname'_id `vtilde'", `r(mse)',`r(N)')
			}
	
		}

		// for each equation: save names of tilde vars with smallest MSE
		// in each resample m and for each variable in nameY/listD/etc.
		// subroutine will handle case of empty lists
		if "`model'"=="interactive" {
			// dependent variable
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(0)
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(1)
			// D variable
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m')
		}
		// late model
		else if "`model'"=="late" {
			// dependent variable
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(0)
			_ddml_crossfit_update_optlist `mname', etype(yeq) vlist(`nameY') m(`m') zett(1)
			// D variable
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m') zett(0)
			_ddml_crossfit_update_optlist `mname', etype(deq) vlist(`listD') m(`m') zett(1)
			// Z variable
			_ddml_crossfit_update_optlist `mname', etype(zeq) vlist(`listZ') m(`m')
		}
		else {
			di as err "internal error - unknown model `model'"
			exit 198
		}
	}

	// set crossfitted field to 1	
	mata: `mname'.crossfitted = 1

	// report results by equation type with resamplings grouped together
	di
	di as res "Reporting crossfitting results:"
	_ddml_crossfit_report `mname'
	
	/*
	// y equation - all models
	forvalues m=1/`reps' {
		_ddml_crossfit_report `mname', etype(yeq) vlist(`nameY') m(`m') model(`model') zett(0)
	}
	forvalues m=1/`reps' {
		_ddml_crossfit_report `mname', etype(yeq) vlist(`nameY') m(`m') model(`model') zett(1)
	}
	// interactive model - D eqn only
	if "`model'"=="interactive" {
		forvalues m=1/`reps' {
			_ddml_crossfit_report `mname', etype(deq) vlist(`listD') m(`m') model(`model')
		}
	}
	// LATE - D and Z eqns
	else if "`model'"=="late" {
		forvalues m=1/`reps' {
			_ddml_crossfit_report `mname', etype(deq) vlist(`listD') m(`m') model(`model') zett(0)
		}
		forvalues m=1/`reps' {
			_ddml_crossfit_report `mname', etype(deq) vlist(`listD') m(`m') model(`model') zett(1)
		}
		forvalues m=1/`reps' {
			_ddml_crossfit_report `mname', etype(zeq) vlist(`listZ') m(`m') model(`model')
		}
	}
	else {
		di as err "internal error - unknown model `model'"
		exit 198
	}
	*/


end


mata:

struct eqnStruct init_eqnStruct()
{
	struct eqnStruct scalar		e
	return(e)
}

void add_to_eqn(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	//(*p).idVtilde	= st_data(., tokens(vnames))
	(*p).MSE		= ((*p).MSE \ mse)
	(*p).N			= n
}

void add_to_eqn01(					struct ddmlStruct m,
									real scalar eqnumber,
									string scalar vnames,
									real scalar mse,
									real scalar n, 
									real scalar Z)
{
	pointer(struct eqnStruct) scalar p
	p				= m.eqnlist[1,eqnumber]
	//(*p).idVtilde	= st_data(., tokens(vnames))
	if (Z==0) {
		(*p).MSE0		= ((*p).MSE0 \ mse)
		(*p).N0			= n
	}
	else {
		(*p).MSE1		= ((*p).MSE1 \ mse)
		(*p).N1			= n
	}

}

end
