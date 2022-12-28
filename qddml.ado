*! ddml v1.1
*! last edited: 18 dec 2022
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
		CMDOPTions(string asis)					///
		NOIsily 								///
		REPs(integer 0)							///
		shortstack 								///
		atet 									///
		NOREG 									///
		]

	mata: s_ivparse("`anything'")

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

	local doreg = "`noreg'"==""

	if "`robust'"!=""	local vce robust
	if "`cluster'"~=""	local vce cluster `cluster'
	if "`vverbose'"=="" local qui qui
	if "`cmd'"=="" local cmd pystacked
	if "`ycmd'"=="" local ycmd `cmd'
	if "`dcmd'"=="" local dcmd `cmd'
	if "`zcmd'"=="" local zcmd `cmd'
	if "`yvtype'"=="" local yvtype `vtype'
	if "`dvtype'"=="" local dvtype `vtype'
	if "`zvtype'"=="" local zvtype `vtype'
	if "`ypredopt'"=="" local ypredopt `predopt'
	if "`dpredopt'"=="" local dpredopt `predopt'
	if "`zpredopt'"=="" local zpredopt `predopt'
	local ycmdoptions `ycmdoptions' `cmdoptions'
	local dcmdoptions `dcmdoptions' `cmdoptions'
	local zcmdoptions `zcmdoptions' `cmdoptions'
	local dhcmd `dcmd'
	local dhcmdoptions `dcmdoptions'

	**** syntax checks
	if ("`model'"=="fiv") {
		if "`dexog'"!="" {
			di as error "no exogenous treatments allowed"
			exit 198
		}
		if "`xctrl'"=="" {
			local ycmd regress
			local ycmdoptions 
			local doreg = 0
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
	`qui' ddml init `model', kfolds(`kfolds') reps(`reps') cluster(`cluster') `tabfold' foldvar(`foldvar')

	*** IV-HD
	if ("`model'"=="fiv") {
		**
		`qui' if (`doreg') ddml E[Y|X], mname(`mname') vname(`depvar') learner(Y0_reg): reg `depvar' `xctrl' 
		`qui' if (`doreg') ddml E[D|X,Z], mname(`mname') vname(`dendog') learner(D0_reg): reg `dendog' `xctrl' `exexog' 
		`qui' if (`doreg') ddml E[D|X], mname(`mname') vname(`dendog') learner(D0_reg): reg {D} `xctrl' 
		**
		`qui' ddml E[Y|X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype'): `ycmd' `depvar' `xctrl', `ycmdoptions'  
		`qui' ddml E[D|X,Z], mname(`mname') vname(`dendog') learner(D1_`dcmd') predopt(`dpredopt') vtype(`dvtype'): `dcmd' `dendog' `xctrl' `exexog', `dcmdoptions' 
		`qui' ddml E[D|X], mname(`mname') vname(`dendog') learner(D1_`dcmd') predopt(`dpredopt') vtype(`dvtype'): `dhcmd' {D} `xctrl', `dhcmdoptions' 
	} 

	*** IV 
	else if ("`model'"=="iv") {
		`qui' if (`doreg') ddml E[Y|X], mname(`mname') vname(`depvar'): reg `depvar' `xctrl'  
		`qui' ddml E[Y|X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype'): `ycmd' `depvar' `xctrl', `ycmdoptions' 
		local j = 1
		foreach d of varlist `dendog' {
			`qui' if (`doreg') ddml E[D|X], mname(`mname') vname(`d'): reg `d' `xctrl'  
			`qui' ddml E[D|X], mname(`mname') vname(`d') predopt(`dpredopt') vtype(`dvtype'): `dcmd' `d' `xctrl', `dcmdoptions' 
			local j = `j' + 1
		}
		local j = 1
		foreach z of varlist `exexog' {
			`qui' if (`doreg') ddml E[Z|X], mname(`mname') vname(`z'): reg `z' `xctrl' 
			`qui' ddml E[Z|X], mname(`mname') vname(`z') predopt(`zpredopt') vtype(`zvtype'): `zcmd' `z' `xctrl', `zcmdoptions' 
			local j = `j' + 1
		}
	}

	*** late / interactive IV
	else if ("`model'"=="late"|"`model'"=="interactiveiv") {
		**
		`qui' if (`doreg') ddml E[Y|Z,X], mname(`mname') vname(`depvar'): reg `depvar' `xctrl' 
		`qui' if (`doreg') ddml E[D|Z,X], mname(`mname') vname(`dendog'): reg `dendog' `xctrl' 
		`qui' if (`doreg') ddml E[Z|X], mname(`mname') vname(`exexog'): reg `exexog' `xctrl' 
		**
		`qui' ddml E[Y|Z,X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype'): `ycmd' `depvar' `xctrl', `ycmdoptions' 
		`qui' ddml E[D|Z,X], mname(`mname') vname(`dendog') predopt(`dpredopt') vtype(`dvtype'): `dcmd' `dendog' `xctrl', `dcmdoptions' 
		`qui' ddml E[Z|X], mname(`mname') vname(`exexog') predopt(`zpredopt') vtype(`zvtype'): `zcmd' `exexog' `xctrl', `zcmdoptions' 
	}

	*** partial linear model
	else if ("`model'"=="partial") {
		`qui' if (`doreg') ddml E[Y|X], mname(`mname') vname(`depvar'): reg `depvar' `xctrl'  
		`qui' ddml E[Y|X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype'): `ycmd' `depvar' `xctrl', `ycmdoptions' 
		local j = 1
		foreach d of varlist `dexog' {
			`qui' if (`doreg') ddml E[D|X], mname(`mname') vname(`d'): reg `d' `xctrl'
			`qui' ddml E[D|X], mname(`mname') vname(`d') predopt(`dpredopt') vtype(`dvtype'): `dcmd' `d' `xctrl', `dcmdoptions' 
			local j = `j' + 1
		}
	}

	*** interactive model
	else if ("`model'"=="interactive") {
		**
		`qui' if (`doreg') ddml E[Y|D,X], mname(`mname') vname(`depvar'): reg `depvar' `xctrl' 
		`qui' if (`doreg') ddml E[D|X], mname(`mname') vname(`dexog'): reg `dexog' `xctrl' 
		**
		`qui' ddml E[Y|D,X], mname(`mname') vname(`depvar') predopt(`ypredopt') vtype(`yvtype'): `ycmd' `depvar' `xctrl', `ycmdoptions' 
		`qui' ddml E[D|X], mname(`mname') vname(`dexog') predopt(`dpredopt') vtype(`dvtype'): `dcmd' `dexog' `xctrl', `dcmdoptions' 
	}	
		
	`qui' ddml crossfit, `noisily' `shortstack'
	if "`verbose'"!="" ddml desc
	ddml estimate, vce(`vce') `atet'

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
