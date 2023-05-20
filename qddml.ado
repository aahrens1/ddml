*! ddml v1.2
*! last edited: 21 feb 2023
*! authors: aa/ms

program define qddml, eclass					//  sortpreserve handled in _ivlasso
	syntax [anything] [if] [in] [aw pw],		/// note no "/" after pw
		Model(name)								///
		[										///
		VERBose VVERBose						///
		mname(name)								///
		kfolds(integer 5)						///
		foldvar(varlist)						///
		TABFold									///
		ROBust									///
		CLUster(varname)						///
		vce(string)								///
		debug 									///
		seed(int 0)								///
		cmd(name)								///
		YCMDOPTions(string asis)				///
		DCMDOPTions(string asis)				///
		ZCMDOPTions(string asis)				///
		YCMD(string asis)						///	
		DCMD(string asis)						///
		ZCMD(string asis)						///
		predopt(string asis)					///
		ypredopt(string asis)					///
		dpredopt(string asis)					///
		zpredopt(string asis)					///
		vtype(string)							///
		yvtype(string)							///  "double", "float" etc
		dvtype(string)							///  "double", "float" etc
		zvtype(string)							///  "double", "float" etc
		CMDOPTions(string asis)					///  will be added on to all learners
		pystacked(string asis)					///
		pystacked_y(string asis)				///
		pystacked_d(string asis)				///
		pystacked_z(string asis)				///
		NOIsily 								///
		REPs(integer 0)							///
		shortstack 								///
		poolstack								///
		stdstack								///
		ssfinalest(name)						///
		psfinalest(name)						///
		atet 									///
		ateu									///
		]

	mata: s_ivparse("`anything'")
	** indicators for pystacked and stacking methods
	local pyflag	= "`pystacked'`pystacked_y'`pystacked_d'`pystacked_z'"~="" | "`cmd'`ycmd'`dcmd'`zcmd'"==""
	local ssflag	= "`shortstack'"~=""
	local psflag	= "`poolstack'"~=""
	local stdflag	= "`stdstack'"~=""
	// if no stacking specified, shortstack
	if ~`ssflag' & ~`psflag' & ~`stdflag' {
		local ssflag		=1
		local shortstack	shortstack
	}

	local depvar	`s(depvar)'
	local dendog	`s(dendog)'
	local dexog		`s(dexog)'
	local xctrl		`s(xctrl)'
	local exexog	`s(exexog)'
	if "`debug'"!="" {
		di "dep var = `depvar'"
		di "treatment endog = `dendog'"
		di "treatment exog = `dexog'"
		di "controls = `xctrl'"
		di "excluded instr = `exexog'"		
	}

	** check if pystacked is available
	cap pystacked
	if _rc == 199 {
		local pystacked_avail = 0
	}

	if "`robust'"!=""	local vce robust
	if "`cluster'"~=""	local vce cluster `cluster'
	if "`verbose'"=="" local qui qui
	if "`yvtype'"=="" local yvtype `vtype'
	if "`dvtype'"=="" local dvtype `vtype'
	if "`zvtype'"=="" local zvtype `vtype'
	if "`ypredopt'"=="" local ypredopt `predopt'
	if "`dpredopt'"=="" local dpredopt `predopt'
	if "`zpredopt'"=="" local zpredopt `predopt'
	
	if "`foldvar'"!="" {
		local foldvar foldvar(`foldvar')
		local kfolds
	}
	else if "`kfolds'"!="" {
		local kfolds kfolds(`kfolds')
	}

	if `pyflag' {
		local ycmd		pystacked
		local dcmd		pystacked
		local zcmd		pystacked
		local dhcmd		pystacked
		if "`pystacked_y'"==""	local ycmdoptions	`pystacked'
		else					local ycmdoptions	`pystacked_y'
		if "`pystacked_d'"==""	local dcmdoptions	`pystacked'
		else					local dcmdoptions	`pystacked_d'
		if "`pystacked_z'"==""	local zcmdoptions	`pystacked'
		else					local zcmdoptions	`pystacked_z'
		local dhcmdoptions		`dcmdoptions'
		foreach opt in ycmdoptions dcmdoptions zcmdoptions dhcmdoptions {
			local `opt' : subinstr local `opt' "||" "||", all count(local doublebarsyntax)
			if `doublebarsyntax'==0 {
				// no || syntax so add a comma to the start of the options unless one is there already
				local `opt'		=trim("``opt''")
				if substr("``opt''",1,1) ~= "," {
					local `opt'		, ``opt''
				}
			}
			else {
				// || syntax so check if there is a comma to be followed by options; add a comma if not
				local `opt' : subinstr local `opt' "," ",", all count(local hascomma)
				if !`hascomma'	local `opt' ``opt'' ,
			}
			// if no standard or pooled stacking, use voting to avoid unnecessary stacking CV steps
			// and use nostdstack option so that voting predicted values aren't created
			if `stdflag'==0 & `psflag'==0 {
				// votetype is ignored if type=reg; relevant only for type=class
				local `opt' ``opt'' voting votetype(soft)
				local nostdstack nostdstack
			}
		}
	}
	else {
		if "`cmd'"=="" local cmd pystacked
		if "`ycmd'"=="" local ycmd `cmd'
		if "`dcmd'"=="" local dcmd `cmd'
		if "`zcmd'"=="" local zcmd `cmd'
		local dhcmd `dcmd'
		// include comma for pystacked compatibility
		local ycmdoptions	, `ycmdoptions'
		local dcmdoptions	, `dcmdoptions'
		local zcmdoptions	, `zcmdoptions'
		local dhcmdoptions	, `dcmdoptions'
	}

	**** syntax checks
	if ("`model'"=="fiv") {
		if "`dexog'"!="" {
			di as error "no exogenous treatments allowed"
			exit 198
		}
		if "`xctrl'"=="" {
			local ycmd regress
			local ycmdoptions 
			local dhcmd regress
			local dhcmdoptions 
		}
	} 
	else if ("`model'"=="iv") {
		if "`dexog'"!="" {
			di as error "no exogenous treatments allowed"
			exit 198
		}
		if ("`xctrl'"=="") {
			di as error "no (high-dimensional) controls specified"
			exit 198
		}
	}
	else if ("`model'"=="late") {
		if "`dexog'"!="" {
			di as error "no exogenous treatments allowed"
			exit 198
		}
		if ("`xctrl'"=="") {
			di as error "no (high-dimensional) controls specified"
			exit 198
		}
	}
	else if ("`model'"=="partial") {
		if "`dendog'"!="" {
			di as error "no endogenous treatments allowed"
			exit 198
		}
		if "`exexog'"!="" {
			di as error "no excluded instruments allowed"
			exit 198
		}
		if ("`xctrl'"=="") {
			di as error "no (high-dimensional) controls specified"
			exit 198
		}
	}
	else if ("`model'"=="interactive") {
		if "`dendog'"!="" {
			di as error "no endogenous treatments allowed"
			exit 198
		}
		if "`exexog'"!="" {
			di as error "no excluded instruments allowed"
			exit 198
		}
		if ("`xctrl'"=="") {
			di as error "no (high-dimensional) controls specified"
			exit 198
		}
	}

	*** model name
	if "`mname'"=="" local mname m0		

	*** estimation
	ddml init `model', `kfolds' reps(`reps') cluster(`cluster') `tabfold' `foldvar'

	*** IV-HD
	if ("`model'"=="fiv") & `pyflag' {
		// special treatment for pystacked - split into separate pystacked calls
		// Y eqn spec
		`ycmd' `depvar' `xctrl' `ycmdoptions' `cmdoptions' noestimate
		forvalues m=1/`e(mcount)' {
			di "Y learner `m':"
			local globalopt `e(globalopt)'
			local globalremove noestimate
			local globalopt : list globalopt - globalremove
			ddml E[Y|X], mname(`mname') vname(`depvar') learner(Y`m'_`e(method`m')') predopt(`ypredopt') vtype(`yvtype') `nostdstack':		///
				pystacked `e(depvar)' `e(xvars_o`m')', method(`e(method`m')') pipe1(`e(pipe`m')') cmdopt1(`e(opt`m')') `globalopt'
		}
		// D eqn spec
		`dcmd' `dendog' `xctrl' `exexog' `dcmdoptions' `cmdoptions' noestimate
		forvalues m=1/`e(mcount)' {
			di "D learner `m':"
			local globalopt `e(globalopt)'
			local globalremove noestimate
			local globalopt : list globalopt - globalremove
			ddml E[D|X,Z], mname(`mname') vname(`dendog') learner(D`m'_`e(method`m')') predopt(`dpredopt') vtype(`dvtype') `nostdstack':		///
				pystacked `e(depvar)' `e(xvars_o`m')', method(`e(method`m')') pipe1(`e(pipe`m')') cmdopt1(`e(opt`m')') `globalopt'
		}
		// DH eqn spec
		`dcmd' `dendog' `xctrl' `dhcmdoptions' `cmdoptions' noestimate
		forvalues m=1/`e(mcount)' {
			di "D learner `m':"
			local globalopt `e(globalopt)'
			local globalremove noestimate
			local globalopt : list globalopt - globalremove
			ddml E[D|X], mname(`mname') vname(`dendog') learner(D`m'_`e(method`m')') predopt(`dpredopt') vtype(`dvtype') `nostdstack':		///
				pystacked {D} `e(xvars_o`m')', method(`e(method`m')') pipe1(`e(pipe`m')') cmdopt1(`e(opt`m')') `globalopt'
		}
	}
	else if ("`model'"=="fiv") {
	// non-pystacked
		ddml E[Y|X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype') `nostdstack':		///
			`ycmd' `depvar' `xctrl' `ycmdoptions' `cmdoptions'
		ddml E[D|X,Z], mname(`mname') vname(`dendog') learner(D1_`dcmd') predopt(`dpredopt') vtype(`dvtype') `nostdstack':	///
			`dcmd' `dendog' `xctrl' `exexog' `dcmdoptions' `cmdoptions'
		ddml E[D|X], mname(`mname') vname(`dendog') learner(D1_`dcmd') predopt(`dpredopt') vtype(`dvtype') `nostdstack':	///
			`dhcmd' {D} `xctrl' `dhcmdoptions' `cmdoptions'
	} 

	*** IV 
	else if ("`model'"=="iv") {
		ddml E[Y|X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype') `nostdstack':	///
			`ycmd' `depvar' `xctrl' `ycmdoptions' `cmdoptions'
		local j = 1
		foreach d of varlist `dendog' {
			ddml E[D|X], mname(`mname') vname(`d') predopt(`dpredopt') vtype(`dvtype') `nostdstack':	///
				`dcmd' `d' `xctrl' `dcmdoptions' `cmdoptions'
			local j = `j' + 1
		}
		local j = 1
		foreach z of varlist `exexog' {
			ddml E[Z|X], mname(`mname') vname(`z') predopt(`zpredopt') vtype(`zvtype') `nostdstack':	///
				`zcmd' `z' `xctrl' `zcmdoptions' `cmdoptions'
			local j = `j' + 1
		}
	}

	*** late / interactive IV
	else if ("`model'"=="late"|"`model'"=="interactiveiv") {
		ddml E[Y|Z,X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype') `nostdstack':	///
			`ycmd' `depvar' `xctrl' `ycmdoptions' `cmdoptions'
		ddml E[D|Z,X], mname(`mname') vname(`dendog') predopt(`dpredopt') vtype(`dvtype') `nostdstack':	///
			`dcmd' `dendog' `xctrl' `dcmdoptions' `cmdoptions'
		ddml E[Z|X], mname(`mname') vname(`exexog') predopt(`zpredopt') vtype(`zvtype') `nostdstack':		///
			`zcmd' `exexog' `xctrl' `zcmdoptions' `cmdoptions'
	}

	*** partial linear model
	else if ("`model'"=="partial") {
		ddml E[Y|X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype') `nostdstack':		///
			`ycmd' `depvar' `xctrl' `ycmdoptions' `cmdoptions'
		local j = 1
		foreach d of varlist `dexog' {
			ddml E[D|X], mname(`mname') vname(`d') predopt(`dpredopt') vtype(`dvtype') `nostdstack':		///
				`dcmd' `d' `xctrl' `dcmdoptions' `cmdoptions'
			local j = `j' + 1
		}
	}

	*** interactive model
	else if ("`model'"=="interactive") {
		ddml E[Y|D,X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype') `nostdstack':	///
			`ycmd' `depvar' `xctrl' `ycmdoptions' `cmdoptions'
		ddml E[D|X], mname(`mname') vname(`dexog') predopt(`dpredopt') vtype(`dvtype') `nostdstack':		///
			`dcmd' `dexog' `xctrl' `dcmdoptions' `cmdoptions'
	}	
		
	ddml crossfit, `noisily' `shortstack' `poolstack' ssfinalest(`ssfinalest') psfinalest(`psfinalest')
	if "`verbose'"!="" ddml desc
	ddml estimate, vce(`vce') `atet' `ateu'

end 

mata:
// basic IV parser
// taken from: https://github.com/statalasso/pdslasso/blob/master/ivlasso.ado
void s_ivparse(string scalar cmdline)
{
	// dep var is first word on line
	depvar	= ""
	i 		= 1
	c		= substr(cmdline,1,1)
	while (!((c==" ") | (c=="") | (c=="(") )) {
		depvar	= depvar + c
		i		= i + 1
		c		= substr(cmdline,i,1)
	}
	// loop through remainder, separating into dexog and parenthetical lists
	dexog	= ""
	dendog	= ""
	xctrl	= ""
	exexog	= ""
	pclauselist = J(0, 1, "")
	while (!(c=="")) {
		
		if (c!="(") {
			while (!((c=="(") | (c==""))) {
				dexog = dexog + c
				i	= i + 1
				c	= substr(cmdline,i,1)
			}
		}
		else {
			// start counter and initialize pclause
			pbal=1
			pclause = ""
			i	= i + 1
			c	= substr(cmdline,i,1)
			while (!((pbal==0) | (c==""))) {
				if (c=="(") {
					pbal=pbal+1
				}
				if (c==")") {
					pbal=pbal-1
				}
				if (pbal!=0) {
					pclause = pclause + c
					i	= i + 1
					c	= substr(cmdline,i,1)
				}
			}
			pclauselist = (pclauselist \ pclause)
			i	= i + 1
			c	= substr(cmdline,i,1)
			
			if (strpos(pclause,"=")) {
				epos	= strpos(pclause,"=")
				// strrpos unavailable in Stata 13
				// epos2	= strrpos(pclause,"=")
				epos2	= strlen(pclause) - strpos(strreverse(pclause),"=") + 1
				if (epos==epos2) {
					dendog = dendog + substr(pclause, 1, epos-1)
					exexog = exexog + substr(pclause, epos+1, .)
				}
				else {
					errprintf("\nsyntax error - too many =s\n")
					exit(198)
				}
			}
			else {
				xctrl = xctrl + pclause + " "
			}
			
		}
	}
	depvar	= strtrim(stritrim(depvar))
	dexog	= strtrim(stritrim(dexog))
	dendog	= strtrim(stritrim(dendog))
	xctrl	= strtrim(stritrim(xctrl))
	exexog	= strtrim(stritrim(exexog))
	st_global("s(exexog)",exexog)
	st_global("s(xctrl)",xctrl)
	st_global("s(dendog)",dendog)
	st_global("s(dexog)",dexog)
	st_global("s(depvar)",depvar)
}
end
