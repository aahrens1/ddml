*** ddml estimation: interactive model
* notes:
* add check somewhere that only a single D is allowed.

program _ddml_estimate_ate_late, eclass sortpreserve

	syntax namelist(name=mname) [if] [in] ,		/// 
								[				///
								ROBust			///
								show(string)	/// dertermines which to post
								clear			/// deletes all tilde-variables (to be implemented)
								post(string)	/// specification to post/display
								REP(integer 1)	/// resampling iteration to post/display
								avplot			///
								debug			///
								* ]
	
	// if post not specified, post optimal model
	if "`post'"=="" {
		local post "opt"
	}

	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	// base sample for estimation - determined by if/in
	marksample touse
	// also exclude obs already excluded by ddml sample
	qui replace `touse' = 0 if `mname'_sample==0

	// locals used below
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	mata: st_local("nreps",strofreal(`mname'.nreps))

	if (`crossfitted'==0) {
		di as err "ddml model not cross-fitted; call `ddml crossfit` first"
		exit 198
	}
	
	// ssflag is a model characteristic but is well-defined only if every equation has multiple learners.
	// will need to check this...
	mata: st_local("ssflag",strofreal(`mname'.ssflag))
	
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameY)
	mata: st_local("vtlistY",invtokens(`eqn'.vtlist))
	foreach var in `vtlistY' {
		local Y0tilde `Y0tilde' `var'0
		local Y1tilde `Y1tilde' `var'1
	}
	mata: `eqn' = (`mname'.eqnAA).get(`mname'.nameD)
	mata: st_local("Dtilde",invtokens(`eqn'.vtlist))

	_ddml_allcombos `Y0tilde'- `Y1tilde' - `Dtilde',		///
		addprefix("")										///
		`debug'  

	local ncombos = r(ncombos)
	local tokenlen = `ncombos'*2 -1
	local y0list `r(colstr1)'
	local y1list `r(colstr2)'
	local Dlist `r(colstr3)'
	
	tempname nmat bmat semat
	// should only ever be 1 D eqn
	mata: `nmat' = J(`ncombos',3,"")
	mata: `bmat' = J(`ncombos'*`nreps',`numeqnD',.)
	mata: `semat' = J(`ncombos'*`nreps',`numeqnD',.)
	
	// simplest if put into a Mata string matrix
	tokenize `y0list' , parse("-")
	forvalues i=1/`ncombos' {
		local idx = 2*`i'-1
		mata: `nmat'[`i',1] = strtrim("``idx''")
	}
	tokenize `y1list' , parse("-")
	forvalues i=1/`ncombos' {
		local idx = 2*`i'-1
		mata: `nmat'[`i',2] = strtrim("``idx''")
	}
	tokenize `Dlist' , parse("-")
	forvalues i=1/`ncombos' {
		local idx = 2*`i'-1
		mata: `nmat'[`i',3] = strtrim("``idx''")
	}
	
	*** shortstack names
	if `ssflag' {
		local Y0ss `nameY'_ss0
		local Y1ss `nameY'_ss1
		local Dss `nameD'_ss
	}
	
	forvalues m=1/`nreps' {
		
		// reset locals
		local Y0opt
		local Y1opt
		local Dopt
		
		*** retrieve best model
		mata: `eqn' = (`mname'.eqnAA).get("`nameY'")
		mata: st_local("Y0opt",return_learner_item(`eqn',"opt0","`m'"))
		mata: st_local("Y1opt",return_learner_item(`eqn',"opt1","`m'"))
		local Y0opt `Y0opt'0
		local Y1opt `Y1opt'1
		mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
		mata: st_local("Dopt",return_learner_item(`eqn',"opt","`m'"))
		
		// text used in output below
		if `nreps'>1 {
			local stext " (sample=`m')"
		}
		
		forvalues i = 1/`ncombos' {
			mata: st_local("y0",`nmat'[`i',1])
			mata: st_local("y1",`nmat'[`i',2])
			mata: st_local("d",`nmat'[`i',3])
			// check if opt for this resample; note === for D so order doesn't matter
			local isopt
			local isY0opt : list Y0opt == y0
			local isY1opt : list Y1opt == y1
			local isDopt : list Dopt === d
			if `isY0opt' & `isY1opt' & `isDopt' {
				local optspec`m' = `i'
				local isopt *
				local DMLtitle "Optimal DDML model"
			}
			else {
				local DMLtitle "DDML"
			}
			qui _ddml_ate_late if `touse',	///
				yvar(`nameY')				///
				dvar(`nameD')				///
				y0tilde(`y0')				///
				y1tilde(`y1')				///
				dtilde(`d')					///
				rep(`m')					///
				mname(`mname')				///
				touse(`touse')
			
			estimates store `mname'_`i'_`m', title("Model `mname', specification `i' resample `m'`isopt'")
			mata: `bmat'[(`m'-1)*`ncombos'+`i',.] = st_matrix("e(b)")
			mata: `semat'[(`m'-1)*`ncombos'+`i',.] = sqrt(diagonal(st_matrix("e(V)"))')
			
		}
		
		if `ssflag' {
			
			qui _ddml_ate_late if `touse',	///
				yvar(`nameY')				///
				dvar(`nameD')				///
				y0tilde(`Y0ss')				///
				y1tilde(`Y1ss')				///
				dtilde(`Dss')				///
				rep(`m')					///
				mname(`mname')				///
				touse(`touse')
			estimates store `mname'_ss_`m', title("Model `mname', shorstack resample `m'")
		
		}
	
	}

	if "`show'"=="all" {
		forvalues m=1/`nreps' {
			// text used in output below
			if `nreps'>1 {
				local stext " (sample=`m')"
			}
			// all combos including optimal model
			forvalues i=1/`ncombos' {
				qui estimates restore `mname'_`i'_`m'
				if "`optspec`m''"=="`i'" {
					local DMLtitle "Optimal DDML model"
				}
				else {
					local DMLtitle "DDML"
				}
				di
				di as text "`DMLtitle'`stext':"
				_ddml_ate_late	// uses replay
				di
			}
			// shortstack
			qui estimates restore `mname'_ss_`m'
			di
			di as text "Shortstack DDML model`stext':"
			_ddml_ate_late	// uses replay
			di
		}
	}
		
	if "`show'"=="optimal" | "`show'"=="all" {
		forvalues m=1/`nreps' {
			qui estimates restore `mname'_`optspec`m''_`m'
			di
			di as text "Optimal DDML specification, resample `m': `optspec`m''"
			_ddml_ate_late
			di
		}
	}
		
	if "`show'"=="shortstack" | "`show'"=="all" {
		forvalues m=1/`nreps' {
			qui estimates restore `mname'_ss_`m'
			di
			di as text "Shortstack DDML model, resample `m':"
			_ddml_ate_late
			di
		}
	}

	di
	di as text "Summary DDML estimation results:"
	di as text "spec  r" %14s "Y0 learner" _c
	di as text           %14s "Y1 learner" _c
	di as text           %14s "D learner" _c
	di as text %10s "b" %10s "SE" _c
	di
	forvalues m=1/`nreps' {
		forvalues i=1/`ncombos' {
			mata: st_local("yt0",`nmat'[`i',1])
			mata: st_local("yt1",`nmat'[`i',2])
			mata: st_local("dt",`nmat'[`i',3])
			if "`optspec`m''"=="`i'" {
				di "*" _c
			}
			else {
				di " " _c
			}
			local specrep `: di %3.0f `i' %3.0f `m''
			// pad out to 6 spaces
			local specrep = (6-length("`specrep'"))*" " + "`specrep'"
			di %6s "{stata estimates replay `mname'_`i'_`m':`specrep'}" _c
			di %14s "`yt0'" _c
			di %14s "`yt1'" _c
			di %14s "`dt'" _c
			mata: st_local("b",strofreal(`bmat'[(`m'-1)*`ncombos'+`i',`j']))
			mata: st_local("se",strofreal(`semat'[(`m'-1)*`ncombos'+`i',`j']))
			di %10.3f `b' _c
			local pse (`: di %6.3f `se'')
			di %10s "`pse'" _c
			di
		}

		if `ssflag' {
			qui estimates restore `mname'_ss_`m'
			local specrep `: di "ss" %3.0f `m''
			// pad out to 6 spaces
			local specrep = "  " + "`specrep'"
			di %6s "{stata estimates replay `mname'_ss_`m':`specrep'}" _c			
			di %14s "[shortstack]" _c
			di %14s "[ss]" _c
			di %14s "[ss]" _c
			di %10.3f el(e(b),1,1) _c
			local pse (`: di %6.3f sqrt(el(e(V),1,1))')
			di %10s "`pse'" _c
			di
		}

	}
	
	// post selected estimates; rep is the resample number (default=1)
	if "`post'"=="opt" {
		qui estimates restore `mname'_`optspec`rep''_`rep'
		di
		di as text "Optimal DDML specification, resample `rep': `optspec`rep''"
		_ddml_ate_late
		di
	}
	else if "`post'"=="shortstack" {
		qui estimates restore `mname'_ss_`rep'
		di
		di as text "Shortstack DDML model, resample `rep':"
		_ddml_ate_late
		di	
	}
	else {
		// post macro should be an integer denoting the specification
		qui estimates restore `mname'_`post'_`rep'
		di
		di as text "DDML specification `post', resample `rep':"
		_ddml_ate_late
		di
	}
	
	// temp Mata object no longer needed
	cap mata: mata drop `eqn' `nmat' `bmat' `semat'


end

// adds rep number suffixes to list of varnames
program define add_suffix, sclass
	syntax [anything] , suffix(name)

	// anything is a list of to-be-varnames that need suffix added to them
	foreach vn in `anything' {
		local vnames `vnames' `vn'`suffix'
	}
	
	sreturn local vnames `vnames'
end

