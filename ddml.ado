*** feb 1, 2021
* ddml v0.3

* notes:
* e.command		= tokens(estcmd)[1,1] fails if command string starts with a prefix e.g. capture
* check for incompatible y variables disabled - can't accommodate prefixes e.g. capture
* spin off init code into a subroutine?
* init code current calls _ddml_sample to set fold var, kfolds, etc. Allow options with init?

// (no)prefix option not implemented; prefixes not added (where prefix = `model'_)

program ddml, eclass

	version 13
	
	local allargs `0'
	
	// split into before/after :
	tokenize "`allargs'", parse(":")
	local maincmd `1'
	macro shift 2
	local eqn `*'
	
	// parse first part using syntax
	local 0 "`maincmd'"
	// options used here already parsed out
	syntax [anything(name=mainargs)]		///
			 `if' `in'						/// if/in sent to _ddml_sample
			 	 , [						///
					mname(name)				///
					newmname(name)			///
					gen(name)				///
					vname(name)				///
					vtype(string)			///  "double", "float" etc
					REPlace					///
					/* NOPrefix */ 			/// don't add model name as prefix (disabled - interferes with save/use option)
					*						///
					]
	// now parse main args; first element is subcmd
	tokenize "`mainargs'"
	local subcmd `1'
	macro shift
	// restmainargs is main args minus the subcommand
	local restmainargs `*'
	
	local allsubcmds	update describe save export use drop copy sample init yeq deq dheq zeq crossfit estimate
	if strpos("`allsubcmds'","`subcmd'")==0 {
		di as err "error - unknown subcommand `subcmd'"
		exit 198
	}

	*** get latest version
	if "`subcmd'"=="update" {
		net install ddml, from(https://raw.githubusercontent.com/aahrens1/ddml/master/)
	} 
	
	*** describe model
	if substr("`subcmd'",1,4)=="desc" {
		local 0 "`restargs'"
		syntax , [ mname(name)  * ]
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_describe `mname', `options'
	}

	*** save model
	if "`subcmd'"=="save" {
		local fname: word 2 of `mainargs'
		local 0 "`restargs'"
		syntax , [ mname(name) * ]
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_save, mname(`mname') fname(`fname') `options'

	}
	
	*** export model
	if "`subcmd'"=="export" {
		local using: word 2 of `mainargs'
		if "`using'"~="using" {
			di as err "invalid syntax - missing destination filename"
			exit 198
		}
		local fname: word 3 of `mainargs'
		local 0 "`restargs'"
		syntax ,[ mname(name) * ]
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_export, mname(`mname') fname(`fname') `options'

	}
	
	*** use model
	if "`subcmd'"=="use" {
		local fname: word 2 of `mainargs'
		local 0 "`restargs'"
		syntax , [ mname(name) * ]
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		// no need to check name
		_ddml_use, mname(`mname') fname(`fname') `options'
	}

	*** drop model
	if "`subcmd'"=="drop" {
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_drop, mname(`mname')
	}

	*** copy model
	if "`subcmd'"=="copy" {
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		if "`newmname'"=="" {
			di as err "error - newmname(.) option required"
			exit 198
		}
		check_mname "`mname'"
		_ddml_copy, mname(`mname') newmname(`newmname')
	}

	*** initialize new estimation
	if "`subcmd'"=="init" {
		local model: word 1 of `restmainargs'
		local allmodels		partial iv interactive late ivhd
		if strpos("`allmodels'","`model'")==0 {
			di as err "no or wrong model specified." 
			exit 198
		}

		if "`mname'"=="" {
			local mname m0 // sets the default name
		}

		// distinct model: no-LIE optimal IV
		// if "`model'"=="ivhd"&"`nolie'"!="" local model ivhd_nolie
		
		mata: `mname'=init_mStruct()
		cap drop `mname'_id
		qui gen double `mname'_id	= _n
		mata: `mname'.id			= st_data(., "`mname'_id")
		// create and store sample indicator; initialized so all obs are used
		cap drop `mname'_sample
		qui gen byte `mname'_sample = 1
		// fill by hand
		mata: `mname'.model			= "`model'"
		// initialize with default fold var, kfolds, number of resamplings, shortstack
		_ddml_sample `if' `in' , mname(`mname') `options'

	}
	
	*** set sample, foldvar, etc.
	if "`subcmd'"=="sample" {
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_sample `if' `in' , mname(`mname') `options'
	}

	*** add equation  
	if "`subcmd'"=="yeq"|"`subcmd'"=="deq"|"`subcmd'"=="dheq"|"`subcmd'"=="zeq" {

		** check that ddml has been initialized
		// to add
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"

		** vtilde: use 2nd and 1st words of eq (estimator) + "_hat" as the default
		if "`gen'"=="" {
			tokenize `"`eqn'"'
			local gen `2'_`1'_hat
		}

		** vname: use 2nd word of eq (dep var) as the default 
		if "`vname'"=="" {
			tokenize `"`eqn'"'
			local vname `2'
		}

		** check that equation is consistent with model
		mata: st_local("model",`mname'.model)
		if ("`subcmd'"=="zeq"&("`model'"=="ivhd"|"`model'"=="ivhd_nolie"|"`model'"=="partial"|"`model'"=="interactive")) {
			di as err "not allowed; zeq not allowed with `model'"
		}
		
		** check that dep var in eqn isn't already used for some other eqn
		** also set flag for whether dep var is new
		mata: st_local("yvar",`mname'.nameY)
		mata: st_local("dvlist",invtokens(`mname'.nameD))
		mata: st_local("zvlist",invtokens(`mname'.nameZ))
		local posof_y : list posof "`vname'" in yvar
		local posof_d : list posof "`vname'" in dvlist
		local posof_z : list posof "`vname'" in zvlist
		if ("`subcmd'"=="yeq" | "`subcmd'"=="deq") & `posof_z' {
			di as err "not allowed - `vname' already in use in Z eqn"
			exit 198
		}
		if ("`subcmd'"=="deq" | "`subcmd'"=="zeq") & `posof_y' {
			di as err "not allowed - `vname' already in use in Y eqn"
			exit 198
		}
		if ("`subcmd'"=="yeq" | "`subcmd'"=="zeq") & `posof_d' {
			di as err "not allowed - `vname' already in use in D eqn"
			exit 198
		}
		// parsimonious way of getting posof that doesn't require checking eqn type
		local posof = `posof_y' + `posof_d' + `posof_z'
		
		// check syntax of D-eq with LIE
		if "`subcmd'"=="dheq" {
			if regexm(`"`eqn'"',"{D}")==0 {
				di as err "placeholder {D} for E[D^|X] is missing"
				exit 198				
			}
		}
				
		add_eqn_to_model,						///
							mname(`mname')		///
							vname(`vname')		///
							vtilde(`gen')		///
							subcmd(`subcmd')	///
							posof(`posof')		///
							estring(`eqn')
	}

	*** cross-fitting
	if "`subcmd'" =="crossfit" {

		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		
		// why is this needed?
		mata: st_global("r(model)",`mname'.model)
		
		// crossfit
		_ddml_crossfit, `options' mname(`mname') 
		
	}

	*** estimate
	if "`subcmd'" =="estimate" {
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		
		check_mname "`mname'"
				
		mata: st_global("r(model)",`mname'.model)

		if ("`r(model)'"=="partial") {
			_ddml_estimate_partial `mname', `options'
		}
		if ("`r(model)'"=="iv") {
			_ddml_estimate_iv `mname', `options'
		}
		if ("`r(model)'"=="interactive") {
			_ddml_estimate_interactive `mname', `options'
		}
		if ("`r(model)'"=="late") {
			_ddml_estimate_late `mname', `options'
		}
		if ("`r(model)'"=="ivhd") {
			_ddml_estimate_ivhd `mname', `options'
		}
		
		// not yet updated
		/*
		_ddml_ereturn, mname(`mname')
		*/

	}
end

prog define check_mname

	args mname

	mata: st_local("isnull",strofreal(findexternal("`mname'")==NULL))
	if `isnull' {
		di as err "model `mname' not found"
		exit 3259
	}
	
	mata: st_local("eltype",eltype(`mname'))
	if "`eltype'"~="struct" {
		di as err "model `mname' is not a struct"
		exit 3259
	}

	mata: st_local("structname",structname(`mname'))
	if "`structname'"~="mStruct" {
		di as err "model `mname' is not an mStruct"
		exit 3000
	}

end

program define add_eqn_to_model, rclass

	syntax [anything],								/// 
							[						///
							mname(name)				/// name of mata struct with model
							vname(varname)			/// name of dep var in equation (to be orthogonalized)
							vtilde(name)			/// names of tilde variable
							subcmd(string)			/// yeq, deq, dheq or zeq
							estring(string asis)	/// names of estimation strings
													/// need asis option in case it includes strings
							posof(integer 0)		/// position of vname in name list; =0 if a new vname (new eqn)
							NOIsily					///
							]
	
	// used for temporary Mata object
	tempname t
	
	// `mname'_vname is name of eqn object in Mata
	local ename `mname'_`vname'
	if `posof'==0 {
		// vname new to model so need a new eqn struct for it
		mata: `ename' = init_eStruct()
		mata: `ename'.vname = "`vname'"
		// if shortstacking, add (default) shorstack varname
		mata: st_local("ssflag",strofreal(`mname'.ssflag))
		if `ssflag' {
			mata: `ename'.shortstack = "`vname'_ss"
		}
	}
	
	// add vtilde to vtlist if not already there
	mata: st_local("vtlist",invtokens(`ename'.vtlist))
	local vtlist `vtlist' `vtilde'
	local vtlist : list uniq vtlist
	// in two steps, to accommodate singleton lists (which are otherwise string scalars and not matrices
	mata: `t' = tokens("`vtlist'")
	mata: `ename'.vtlist	= `t'
	
	// used below with syntax command
	local 0 `"`estring'"'
	// parse estimation string into main command and options; if and in will be stripped out
	syntax [anything] [if] [in] , [*]
	local est_main `anything'
	local est_options `options'
	if "`subcmd'"=="dheq" {
		mata: add_learner_item(`ename',"`vtilde'","estring_h","`0'")
		mata: add_learner_item(`ename',"`vtilde'","est_main_h","`est_main'")
		mata: add_learner_item(`ename',"`vtilde'","est_options_h","`est_options'")
		mata: `ename'.lieflag = 1
	}
	else {
		mata: add_learner_item(`ename',"`vtilde'","estring","`0'")
		mata: add_learner_item(`ename',"`vtilde'","est_main","`est_main'")
		mata: add_learner_item(`ename',"`vtilde'","est_options","`est_options'")
		// update nlearners - counts deq and dheq as a single learner
		mata: `ename'.nlearners = cols(`ename'.vtlist)
	}

	// insert eqn struct into model struct
	mata: (`mname'.eqnAA).put("`vname'",`ename')
	
	// update rest of model struct
	if `posof'==0 {
		mata: `t' = tokens("`vname'")
		if "`subcmd'"=="yeq" {
			// only ever 1 y eqn
			mata: `mname'.nameY = "`vname'"
			}
		else if "`subcmd'"=="deq" | "`subcmd'"=="dheq" {
			mata: `mname'.nameD = (`mname'.nameD, `t')			
		}
		else if "`subcmd'"=="zeq" {
			mata: `mname'.nameZ = (`mname'.nameZ, `t')			
		}
	}
	
	// no longer needed so clear from Mata
	cap mata: mata drop `t'

end