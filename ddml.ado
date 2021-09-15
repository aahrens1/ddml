*** feb 1, 2021
* ddml v0.3

* notes:
* e.command		= tokens(estcmd)[1,1] fails if command string starts with a prefix e.g. capture
* check for incompatible y variables disabled - can't accommodate prefixes e.g. capture
* spin off init code into a subroutine?
* init code current calls _ddml_sample to set fold var, kfolds, etc. Allow options with init?

program ddml, eclass

	version 13
	
	local allargs `0'
	tokenize "`allargs'", parse(",")
	local mainargs `1'
	macro shift
	local restargs `*'

	// local subcmd : word 1 of `mainargs'
	tokenize "`mainargs'"
	local subcmd `1'
	macro shift
	// restmainargs is main args minus the subcommand
	local restmainargs `*'
	
	local allsubcmds	update describe save export use drop copy sample init yeq deq zeq crossfit estimate dheq
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
		local 0 "`restargs'"
		syntax , [mname(name)]
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_drop, mname(`mname')
	}

	*** copy model
	if "`subcmd'"=="copy" {
		local 0 "`restargs'"
		syntax , mname(name) newmname(name)
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_copy, mname(`mname') newmname(`newmname')
	}

	*** initialize new estimation
	if "`subcmd'"=="init" {
		local model: word 2 of `mainargs'
		if ("`model'"!="partial"&"`model'"!="iv"&"`model'"!="interactive"&"`model'"!="late"&"`model'"!="optimaliv") {
			di as err "no or wrong model specified." 
			exit 1
		}
		local 0 "`restargs'"
		// fold variable is option; default is ddmlfold
		syntax `if' `in', [mname(name) NOLie *]

		if "`mname'"=="" {
			local mname m0 // sets the default name
		}

		// distinct model: no-LIE optimal IV
		if "`model'"=="optimaliv"&"`nolie'"!="" local model optimaliv_nolie

		mata: `mname'=init_ddmlStruct()
		// create and store id variable
		cap drop `mname'_id
		qui gen double `mname'_id	= _n
		mata: `mname'.id			= st_data(., "`mname'_id")
		// create and store sample indicator; initialized so all obs are used
		cap drop `mname'_sample
		qui gen byte `mname'_sample = 1
		// add sample indicator to model struct (col 1 = id, col 2 = fold id)
		mata: `mname'.idSample		= st_data(., ("`mname'_id", "`mname'_sample"))
		// fill by hand
		mata: `mname'.model			= "`model'"
		mata: `mname'.crossfitted	= 0
		
		// initialize with default fold var, kfolds, number of resamplings
		_ddml_sample `if' `in' , mname(`mname') `options'
	}
	
	*** set sample, foldvar, etc.
	if "`subcmd'"=="sample" {
		local 0 "`restmainargs' `restargs'"
		syntax [if] [in] , [mname(name)  * ]
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"
		_ddml_sample `if' `in' , mname(`mname') `options'
	}

	*** add equation  
	if "`subcmd'"=="yeq"|"`subcmd'"=="deq"|"`subcmd'"=="zeq"|"`subcmd'"=="dheq" {

		** parsing
		// macro options has eqn to be estimated set off from the reset by a :
		tokenize `" `restargs' "', parse(":")
		// parse character is in macro `2'
		local eqn `3'
		local 0 "`1'"
		syntax ,	///			
					[				///
					gen(name)		///
					vname(name)		///
					genh(name)		/// (intended for LIE)
					mname(name)		///
					vtype(string)   ///  "double", "float" etc
					REPlace         ///
					NOPrefix 		/// don't add model name as prefix
					]

		** check that ddml has been initialized
		// to add
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"

		** vname: use 2nd word of eq as the default 
		di "`eqn'"
		if "`vname'"=="" {
			local vname : word 2 of `eqn'
		}

		** 
		if "`gen'"=="" {
			local gen `vname'_t
		}

		** drop gen var if it already exists
		if "`replace'"!="" {
			cap drop `gen'
		}

		** check that equation is consistent with model
		mata: st_local("model",`mname'.model)
		if ("`subcmd'"=="zeq"&("`model'"=="optimaliv"|"`model'"=="optimaliv_nolie"|"`model'"=="partial"|"`model'"=="interactive")) {
			di as err "not allowed; zeq not allowed with `model'"
		}
		if ("`subcmd'"=="dheq"&"`model'"=="optimaliv") {
			di as err "not allowed; dheq not allowed with `model' and nolie"
			exit 198
		}
		if ("`subcmd'"=="dheq"&("`model'"!="optimaliv_nolie")) {
			di as err "not allowed; dheq not allowed with `model'"
			exit 198
		}

		** add prefix to vtilde
		if "`noprefix'"=="" {
			local prefix `mname'_
		}
		else {
			local prefix
		}

		** split equation -- only required for D-eq with LIE
		if "`subcmd'"=="deq"&"`model'"=="optimaliv" {
			tokenize `" `eqn' "', parse("|")
			// parse character is in macro `2'
			local eqn `1'
			local eqn_h `3'
			if "`2'`3'"=="" {
				di as err "estimation command for E[D^|X] missing"
				exit 198
			}
			if regexm("`eqn_h'","{D}")==0 {
				di as err "placeholder {D} missing in E[D^|X] command (2nd part)"
				exit 198				
			}
		}

		if "`model'"=="optimaliv"&"`genh'"=="" {
			local genh `gen'_h
		}

		// subcmd macro tells add_eqn(.) which list to add it to
		mata: add_eqn(`mname', "`subcmd'", "`vname'", "`gen'","`eqn'","`vtype'","`prefix'","`eqn_h'","`genh'")
		local newentry `r(newentry)'
		if "`subcmd'"=="yeq" {
			// check if nameY is already there; if it is, must be identical to vname here
			mata: st_global("r(vname)",`mname'.nameY)
			if "`r(vname)'"=="" {
				mata: `mname'.nameY		= "`vname'"
			}
			/*
			*** disabled - doesn't work with cap or other prefixes ***
			else if "`r(vname)'"~="`vname'" {
				di as err "error - incompatible y variables"
				exit 198
			}
			*/
		}
		if "`subcmd'"=="deq" {
			// check if nameD already has vname; if not, add it to the list
			mata: st_global("r(vname)",invtokens(`mname'.nameD))
			if "`r(vname)'"=="" {
				mata: `mname'.nameD		= "`vname'"
			}
			else {
				local dlist `r(vname)' `vname'
				local dlist : list uniq dlist
				mata: `mname'.nameD		= tokens("`dlist'")
			}
		}
		if "`subcmd'"=="deq"&"`model'"=="optimaliv" {
			// check if nameD already has vname; if not, add it to the list
			mata: st_global("r(vname)",invtokens(`mname'.nameD))
			if "`r(vname)'"=="" {
				mata: `mname'.nameDH		= "`vname'"
			}
			else {
				local dlist `r(vname)' `vname'
				local dlist : list uniq dlist
				mata: `mname'.nameDH		= tokens("`dlist'")
			}
		}		
		if "`subcmd'"=="dheq" {
			// check if nameDH already has vname; if not, add it to the list
			mata: st_global("r(vname)",invtokens(`mname'.nameDH))
			if "`r(vname)'"=="" {
				mata: `mname'.nameDH	= "`vname'"
			}
			else {
				local dhlist `r(vname)' `vname'
				local dhlist : list uniq dhlist
				mata: `mname'.nameDH	= tokens("`dhlist'")
			}
		}
		if "`subcmd'"=="zeq" {
			// check if nameZ already has vname; if not, add it to the list
			mata: st_global("r(vname)",invtokens(`mname'.nameZ))
			if "`r(vname)'"=="" {
				mata: `mname'.nameZ		= "`vname'"
			}
			else {
				local zlist `r(vname)' `vname'
				local zlist : list uniq zlist
				mata: `mname'.nameZ		= tokens("`zlist'")
			}
		}
		if `newentry' {
			di as text "Equation successfully added."
		}
		else {
			di as text "Equation successfully replaced."
		}
	}

	*** cross-fitting
	if "`subcmd'" =="crossfit" {

		local 0 "`restargs'"
		syntax , [ mname(name) nnls * ]
		if "`mname'"=="" {
			local mname m0 // sets the default name
		}
		check_mname "`mname'"

		// clear previous results or initialize if no prev results
		_ddml_reset_model_results, mname(`mname')

		mata: st_global("r(model)",`mname'.model)

		if ("`r(model)'"=="partial") {
		_ddml_crossfit_additive , `options' mname(`mname') 
		}
		if ("`r(model)'"=="iv") {
		_ddml_crossfit_additive , `options' mname(`mname') 
		}
		if ("`r(model)'"=="interactive") {
		_ddml_crossfit_interactive , `options' mname(`mname') 
		}
		if ("`r(model)'"=="late") {
		_ddml_crossfit_interactive , `options' mname(`mname') 
		}
		if ("`r(model)'"=="optimaliv") {
		_ddml_crossfit_additive , `options' mname(`mname') 
		}
		if ("`r(model)'"=="optimaliv_nolie") {
		_ddml_crossfit_additive , `options' mname(`mname') 
		}

		if ("`nnls'"!="") _ddml_combine, mname(`mname')

		// set model crossfitted flag = 1
		mata: `mname'.crossfitted	= 1

	}

	*** estimate
	if "`subcmd'" =="estimate" {
		local 0 "`restargs'"
		// mname is required; could make optional with a default name
		syntax , [mname(name) * resample(integer 1)]

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
		if ("`r(model)'"=="optimaliv") {
			_ddml_estimate_optimaliv `mname', `options'
		}
		if ("`r(model)'"=="optimaliv_nolie") {
			_ddml_estimate_optimaliv `mname', `options'
		}

		_ddml_ereturn, mname(`mname')


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
	if "`structname'"~="ddmlStruct" {
		di as err "model `mname' is not a ddmlStruct"
		exit 3000
	}

end
 