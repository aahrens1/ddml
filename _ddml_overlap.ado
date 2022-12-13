*! ddml v1.0
*! last edited: 13 dec 2022
*! authors: aa/ms

program define _ddml_overlap

	syntax 		, [									///
				mname(name)							///
				replist(numlist integer min=1)		/// list of resamples
				pslist(namelist)					/// list of propensity scores excl resample 
				n(integer 0)						/// number of points (default = N)
				kernel(name)						/// default = triangle
				lopt0(string)						/// line options for d=0
				lopt1(string)						/// line options for d=1
				title(string)						/// title for combined graph
				subtitle(string)					/// subtitle for combined graph
				name(string)						/// name of combined graph; can include ", replace"
				*]
	
	// blank eqn - declare this way so that it's a struct and not transmorphic
	tempname eqn
	mata: `eqn' = init_eStruct()
	
	mata: st_local("model",`mname'.model)
	mata: st_local("crossfitted",strofreal(`mname'.crossfitted))	// flag for crossfitting results available
	mata: st_local("nreps",strofreal(`mname'.nreps))
	mata: st_local("nameY",`mname'.nameY)
	mata: st_local("nameD",invtokens(`mname'.nameD))
	mata: st_local("nameZ",invtokens((`mname'.nameZ)))
	local numeqnD : word count `nameD'
	local numeqnZ : word count `nameZ'
	
	if "`model'"~="interactive" & "`model'"~="late" {
		di as err "error - overlap supported only for interactive or late models"
		exit 198
	}
	if "`model'"=="interactive" & `numeqnD'>1 {
		di as err "error - only one treatment variable currently supported in interactive model"
		exit 198
	}
	if `crossfitted'==0 {
		di as err "error - model not crossfitted"
		exit 198
	}
	
	// default title
	if "`title'"=="" {
		local title "Propensity scores by treatment group"
	}
	// default replist, graph subtitle
	if "`replist'"=="" {
		local replist 1/`nreps'
		if `nreps'>1 & "`subtitle'"=="" {
			local subtitle "all crossfit samples"
		}
	}
	else if "`subtitle'"=="" {
		local subtitle "reps=`replist'"
	}
	
	// default list of propensity scores (prefixes)
	if "`pslist'"=="" {
		// eqn has info about learners
		if "`model'"=="interactive" {
			mata: `eqn' = (`mname'.eqnAA).get("`nameD'")
		}
		else {
			mata: `eqn' = (`mname'.eqnAA).get("`nameZ'")
		}
		mata: st_local("pslist",invtokens(`eqn'.vtlist))
	}
	// labels for propensity scores
	if "`model'"=="interactive" {
		local vlab0 "D=0"
		local vlab1 "D=1"
	}
	else {
		local vlab0 "Z=0"
		local vlab1 "Z=1"
	}
		
	// default number of points
	if `n'==0 {
		qui count if `mname'_sample
		local n=r(N)
	}
	// default kernel
	if "`kernel'"=="" {
		local kernel triangle
	}
	// default line options
	if "`lopt0'"=="" {
		local lopt0 lpattern(solid) lcolor(navy)
	}
	if "`lopt1'"=="" {
		local lopt1 lpattern(shortdash) lcolor(dkorange)
	}
	
	// loop through propensity scores
	foreach dtilde in `pslist' {
		// gname is individual dtilde graph
		local gname `dtilde'
		// loop through resamples
		foreach r of numlist `replist' {
			tempvar x0`r' x1`r' ps0`r' ps1`r' ps`r'
			qui gen `ps`r'' = `dtilde'_`r'
			kdensity `ps`r'' if `nameD'==0, kernel(`kernel') n(`n') nograph gen(`x0`r'' `ps0`r'')
			kdensity `ps`r'' if `nameD'==1, kernel(`kernel') n(`n') nograph gen(`x1`r'' `ps1`r'')
			local gcmd `gcmd'												///
				(line `ps0`r'' `x0`r'', `lopt0')							///
				(line `ps1`r'' `x1`r'', `lopt1')
		}
		label var `ps01' "`vlab0'"
		label var `ps11' "`vlab1'"
		twoway `gcmd',														///
			title("`dtilde'")												///
			xtitle("Propensity score")										///
			ytitle("Density")												///
			legend(order(1 2))												///
			nodraw															///
			name(`gname', replace)
	}
	graph combine `pslist', title("`title'") subtitle("`subtitle'") name(`name')
	// drop separate graphs
	cap graph drop `pslist'
		
end
