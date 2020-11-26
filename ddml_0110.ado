*! ddml v0.1.10 (1 nov 2020)

*** features to be added 
* - absorb() to partial out fixed effects
* - markout
* - optimalIV and LATE model

program ddml_0110, eclass

    version 13

    syntax [anything] , [*]

    local subcmd: word 1 of `anything'

    *** get latest version
    if "`subcmd'"=="update" {
        net install ddml, from(https://raw.githubusercontent.com/aahrens1/ddml/master/)
    } 
    
    *** initialize new estimation
    if "`subcmd'"=="init" {
        local model: word 2 of `anything'
        if ("`model'"!="partial"&"`model'"!="iv"&"`model'"!="interactive"&"`model'"!="late"&"`model'"!="optimaliv") {
            di as err "no or wrong model specified." 
            exit 1
        }
        ereturn clear
        ereturn local cmd "ddml_init"
        ereturn local model "`model'"
        ereturn scalar yest = 0    // make these hidden later
        ereturn scalar dest = 0  
        ereturn scalar zest = 0
        ereturn scalar crossfit = 0
    } 

    *** add equation  
    if "`subcmd'"=="yeq"|"`subcmd'"=="deq"|"`subcmd'"=="zeq" {

        ** check that equation is consistent with model
        if ("`subcmd'"=="yeq"&"`model'"=="optimaliv") {
            di as err "not allowed; yeq not allowed with `model'"
        }
        if ("`subcmd'"=="zeq"&("`model'"=="optimaliv"|"`model'"=="partial"|"`model'"=="interactive")) {
            di as err "not allowed; deq not allowed with `model'"
        }

        ** check that ddml has been initialized
        if ("`e(cmd)'"!="ddml_init") {
            di as err "no ddml init object found." _c
            di as err "you need to initialize using 'ddml init'"
            exit 1
        }

        ** parsing
        local estname: word 2 of `anything'
        gettoken estname: estname, parse(":")
        gettoken left theeq: anything, parse(":")
        gettoken left theeq: theeq, parse(":") 
        gettoken theeq: theeq, parse("(") match(paren)

        local 0 `theeq'
        syntax [anything] [if] [in], [*]
        if "`if'`in'"!="" {
            di as err "if and in not allowed in equation"
            error 198
        }
        local cmdline `anything'
        local cmd: word 1 of `cmdline'
        local vars: list cmdline - cmd
        local depvar: word 1 of `vars'
        local xvars: list vars - depvar  
        //qui ds `xvars'
        //local xvars = r(varlist)
        local yopts`i' `options'

        ** ereturn
        if "`subcmd'"=="yeq" {
            ereturn local yname`=`e(yest)'+1' "`estname'"
            ereturn local ycmd`=`e(yest)'+1' "`cmd'"
            ereturn local yopts`=`e(yest)'+1' "`options'"
            ereturn local yvar`=`e(yest)'+1' "`depvar'"
            ereturn local yxvars`=`e(yest)'+1' "`xvars'"
            ereturn scalar yest = `e(yest)'+1
            qui gen `estname' = .
        }
        if "`subcmd'"=="deq" {
            ereturn local dname`=`e(dest)'+1' "`estname'"
            ereturn local dcmd`=`e(dest)'+1' "`cmd'"
            ereturn local dopts`=`e(dest)'+1' "`options'"
            ereturn local dvar`=`e(dest)'+1' "`depvar'"
            ereturn local dxvars`=`e(dest)'+1' "`xvars'"
            ereturn scalar dest = `e(dest)'+1
            qui gen `estname' = .
        }
        if "`subcmd'"=="zeq" {
            ereturn local zname`=`e(zest)'+1' "`estname'"
            ereturn local zcmd`=`e(zest)'+1' "`cmd'"
            ereturn local zopts`=`e(zest)'+1' "`options'"
            ereturn local zvar`=`e(zest)'+1' "`depvar'"
            ereturn local zxvars`=`e(zest)'+1' "`xvars'"
            ereturn scalar zest = `e(zest)'+1
            qui gen `estname' = .
        }
    }

    *** cross-fitting
    if "`subcmd'" =="crossfit" {
        _ddml_crossfit, `options'
    }

    *** estimate
    if "`subcmd'" =="estimate" {
        if ("`e(cmd)'"=="ddml_init") {
            di as err "you first need to call 'ddml crossfit'"
            exit 1
        }
        if ("`e(cmd)'"!="ddml_crossfit") {
            di as err "unrecognized command."
            exit 1
        }
        if ("`e(model)'"=="partial") {
            _ddml_estimate_partial, `options'
        }
        if ("`e(model)'"=="iv") {
            _ddml_estimate_iv, `options'
        }
        if ("`e(model)'"=="interactive") {
            _ddml_estimate_interactive, `options'
        }
        if ("`e(model)'"=="late") {
            _ddml_estimate_late, `options'
        }
        if ("`e(model)'"=="optimaliv") {
            _ddml_estimate_optimaliv, `options'
        }
    }
end

*** ddml cross-fitting
program _ddml_crossfit, eclass sortpreserve

    syntax [anything] [if] [in] , /// 
							[ kfolds(integer 2)  ///
							NOIsily ///
							debug /// 
							Robust ///
							TABFold ///
							yrclass ///
							drclass /// 
                            ]

    *** check there are enough equations
    local yestn = e(yest)
    local destn = e(dest) 
    local zestn = e(zest)
    local model = e(model)
    if ("`model'"=="partial"|"`model'"=="interactive") {
        if (`yestn'==0 | `destn'==0) {
            di as err "insufficient equations specified." _c 
            di as err "we need at least ML command for y and one for d."
            exit 1
        }
        if `zestn'>0 {
            di as err "zeq not allowed with `model'"
            exit 1
        }
    }
    if ("`model'"=="iv"|"`model'"=="late") {
        if (`yestn'==0 | `destn'==0 | `zestn'==0) {
            di as err "insufficient equations specified." _c 
            di as err "we need at least one ML command for y, one for d, one for z."
            exit 1
        }
    }
    if ("`model'"=="optimaliv") {
        if (`destn'==0) {
            di as err "insufficient equations specified." _c 
            di as err "we need at least one ML command for d."
            exit 1
        }
    }
    
    *** check that dependent variable in y and d-equation is always the same
    local yvar = e(yvar1)
    local dvar = e(dvar1)
    if (`yestn'>1) {
        forvalues i = 1(1)`yestn' {
            local yvar`i' = e(yvar`i')
            if "`yvar`i''" != "`yvar'" {
                di as err "inconsistent d-variables: `yvar', `yvar`i''"
                exit 1
            }
        }
    }
    if (`destn'>1) {
        forvalues i = 1(1)`destn' {
            local dvar`i' = e(dvar`i')
            if "`dvar`i''" != "`dvar'" {
                di as err "inconsistent d-variables: `dvar', `dvar`i''"
                exit 1
            }
        }
    }
    if (`zestn'>0) {
        local zvar = e(zvar1)
        forvalues i = 1(1)`zestn' {
            local zvar`i' = e(zvar`i')
            if "`zvar`i''" != "`zvar'" {
                di as err "inconsistent z-variables: `zvar', `zvar`i''"
                exit 1
            }
        }
    }

    *** check that variables are binary if required
    if ("`model'"=="late") {
        cap assert `zvar'==1 | `zvar' == 0 | missing(`zvar') // if `sample'
        if _rc==9 {
            di as err "`zvar' should be binary"
        }
    }
    if ("`model'"=="interactive"|"`model'"=="late") {
        cap assert `dvar'==1 | `dvar' == 0 | missing(`dvar') // if `sample'
        if _rc==9 {
            di as err "`dvar' should be binary"
        }
        cap assert `yvar'==1 | `yvar' == 0 | missing(`yvar') // if `sample'
        if _rc==9 {
            di as err "`yvar' should be binary"
        }
    }

    *** gen folds
	tempvar kid uni cuni
	gen double `uni' = runiform()
	cumul `uni', gen(`cuni')
	gen `kid' =ceil(`kfolds'*`cuni')
	if ("`tabfold'"!="") {
		di ""
		di "Overview of frequencies by fold:"
		tab `kid'
		di ""
	}
	//

    *** retrieve y macros
    forvalues i = 1(1)`yestn' {
        local yname`i' = e(yname`i')
        local ycmd`i' = e(ycmd`i')
        local yxvar`i' = e(yxvars`i')
        local yopts`i' = e(yopts`i')
        if "`yopts`i''"=="." { // in case local is empty
            local yopts`i'
        }
    }

	*** retrieve d macros
    forvalues i = 1(1)`destn' {
        local dname`i' = e(dname`i')
        local dcmd`i' = e(dcmd`i')
        local dxvar`i' = e(dxvars`i')	
        local dopts`i' = e(dopts`i')
        local dopts`i' = e(dopts`i')
        if "`dopts`i''"=="." { // in case local is empty
            local dopts`i'
        }
    }

    *** retrieve z macros
    if (`zestn'>0) {
        forvalues i = 1(1)`zestn' {
            local zname`i' = e(zname`i')
            local zcmd`i' = e(zcmd`i')
            local zxvar`i' = e(zxvars`i')	
            local zopts`i' = e(zopts`i')
            local zopts`i' = e(zopts`i')
            if "`zopts`i''"=="." { // in case local is empty
                local zopts`i'
            }
        }
    }

    *** initialize ytilde / dtilde temp vars
    forvalues i = 1(1)`yestn' {
        tempvar ytilde`i'
        tempvar ytilde_temp`i'
        qui gen `ytilde`i'' = .
    }
    forvalues i = 1(1)`destn' {
        tempvar dtilde`i'
        tempvar dtilde_temp`i'
        qui gen `dtilde`i'' = .
    }
    if ("`model'"=="iv"|"`model'"=="late") {
        forvalues i = 1(1)`zestn' {
            tempvar ztilde`i'
            tempvar ztilde_temp`i'
            qui gen `ztilde`i'' = .
        }
    }

    *** show locals (debug-mode)
    if ("`debug'"!="") {
        di "---debug mode----show locals---------------------"
        di "yvar: `yvar'"
        di "dvar: `dvar'"
        di "-------------------------------------------------"
        forvalues i = 1(1)`destn' {
            di "dname`i': `dname`i''"
            di "dcmd`i': `dcmd`i''"
            //di "dcmdline`i': `dcmdline`i''"
            di "dxvar`i': `dxvar`i''"
            di "dopts`i': `dopts`i''"
        }
        di "-----------------------------------"
        forvalues i = 1(1)`yestn' {
            di "yname`i': `yname`i''"
            di "ycmd`i': `ycmd`i''"
            //di "ycmdline`i': `ycmdline`i''"
            di "yxvar`i': `yxvar`i''"
            di "yopts`i': `yopts`i''"
        }
        di "-----------------------------------"
        if ("`model'"=="iv"|"`model'"=="late") {
            forvalues i = 1(1)`zestn' {
                di "zname`i': `zname`i''"
                di "zcmd`i': `zcmd`i''"
                di "zxvar`i': `zxvar`i''"
                di "zopts`i': `zopts`i''"
            }
        }
        di "-----------------------------------"
    }

    *** do cross-fitting
    di "Cross-fitting fold " _c
	forvalues k = 1(1)`kfolds' {
	
		if (`k'==`kfolds') {
			di as text "`k'"
		}
		else {
			di as text "`k' " _c
		}
		
		// ML is applied to I^c sample (all data ex partition k)
		qui {
		
			** y equation(s)
            if ("`model'"=="partial"|"`model'"=="iv") {
                forvalues i = 1(1)`yestn' {
                    cap drop `ytilde_temp`i''
                    `ycmd`i'' `yvar' `yxvar`i'' if `kid'!=`k', `yopts`i''
                    predict `ytilde_temp`i'' if `kid'==`k' 
                    replace `ytilde`i'' = `yvar' - `ytilde_temp`i'' if `kid'==`k'
                }
            }
			
			** y0 and y1 equation(s)
            if ("`model'"=="late"|"`model'"=="interactive") {
                forvalues i = 1(1)`yestn' {
                    cap drop `ytilde_temp`i''
                    `ycmd`i'' `yvar' `yxvar`i'' if `kid'!=`k' & `dvar'==0, `yopts`i''
                    predict `ytilde_temp`i'' if `kid'==`k' & `dvar'==0
                    replace `ytilde`i'' = `yvar' - `ytilde_temp`i'' if `kid'==`k' & `dvar'==0
                }
                forvalues i = 1(1)`yestn' {
                    cap drop `ytilde_temp`i''
                    `ycmd`i'' `yvar' `yxvar`i'' if `kid'!=`k' & `dvar'==1, `yopts`i''
                    predict `ytilde_temp`i'' if `kid'==`k' & `dvar'==1
                    replace `ytilde`i'' = `yvar' - `ytilde_temp`i'' if `kid'==`k' & `dvar'==1
                }
            }

			** d equation(s)
            if ("`model'"=="partial"|"`model'"=="iv"|"`model'"=="interactive"|"`model'"=="optimaliv") {
                forvalues i = 1(1)`destn' {
                    cap drop `dtilde_temp`i''
                    `dcmd`i'' `dvar' `dxvar`i'' if `kid'!=`k', `dopts`i''
                    predict `dtilde_temp`i'' if `kid'==`k'   
                    replace `dtilde`i'' = `dvar' - `dtilde_temp`i'' if `kid'==`k'
                }
            }

			** d0 and d1 equation(s)
            if ("`model'"=="late") {
                forvalues i = 1(1)`destn' {
                    cap drop `dtilde_temp`i''
                    `dcmd`i'' `dvar' `dxvar`i'' if `kid'!=`k' & `zvar'==0, `dopts`i''
                    predict `dtilde_temp`i'' if `kid'==`k' & `zvar'==0  
                    replace `dtilde`i'' = `dvar' - `dtilde_temp`i'' if `kid'==`k' & `zvar'==0
                }
                forvalues i = 1(1)`destn' {
                    cap drop `dtilde_temp`i''
                    `dcmd`i'' `dvar' `dxvar`i'' if `kid'!=`k' & `zvar'==1, `dopts`i''
                    predict `dtilde_temp`i'' if `kid'==`k' & `zvar'==1
                    replace `dtilde`i'' = `dvar' - `dtilde_temp`i'' if `kid'==`k' & `zvar'==1
                }
            }

            ** z equation(s)
            if ("`model'"=="iv"|"`model'"=="late") {
                forvalues i = 1(1)`zestn' {
                    cap drop `ztilde_temp`i''
                    `zcmd`i'' `zvar' `zxvar`i'' if `kid'!=`k', `zopts`i''
                    predict `ztilde_temp`i'' if `kid'==`k'   
                    replace `ztilde`i'' = `zvar' - `ztilde_temp`i'' if `kid'==`k'
                }
            }
		}
		
        // DML 1: currently not supported
		//`qui' reg `ytilde' `dtilde' if `kid'==`k', `robust' nocons
		//qui mat `bhat_k'[`k',1]=_b[`dtilde']
		//qui mat `se_k'[`k',1]=_se[`dtilde']
			
	}	
	//mat list `bhat_k'
	//`qui' mat list `se_k'

    *** calculate and report MSE 
    if ("`model'"=="partial"|"`model'"=="iv") {
        tempname MSEy
        local yminmse = .
        local yminmseid = 1
        mat `MSEy' = J(1,`yestn',.)
        forvalues i = 1(1)`yestn' {
            qui replace `yname`i'' = `ytilde`i''
            tempname ytilde_sq`i'
            qui gen double `ytilde_sq`i'' = (`ytilde`i'')^2
            qui sum `ytilde_sq`i'' , meanonly
            mat `MSEy'[1,`i'] = r(mean)
            local newyminmse = r(mean)
            if (`newyminmse'<`yminmse') {
                local yminmse = `newyminmse'
                local yminmseid = `i'
            }
        }
        di " "
        di as res "Mean-squared error for y|X:"
        di _col(2) _c
        di _col(10) "Name" _c
        di _col(30) "Command" _c
        di _col(50) "MSPE"
        di "{hline 65}"
        forvalues i = 1(1)`yestn' {
            di _col(2) "`i'" _c
            di _col(10) "`yname`i''" _c
            di _col(30) "`ycmd`i''" _c
            if (`yminmseid'==`i') {
                di _col(50) `MSEy'[1,`i'] _c
                di "*"
            } 
            else {
                di _col(50) `MSEy'[1,`i'] 
            }
        }
    }
    if ("`model'"=="late"|"`model'"=="interactive") {
        tempname MSEy0 
        tempname MSEy1
        local yminmse0 = .
        local yminmse1 = .
        local yminmseid1 = 1
        local yminmseid0 = 1
        mat `MSEy0' = J(1,`yestn',.)
        mat `MSEy1' = J(1,`yestn',.)
        forvalues i = 1(1)`yestn' {
            qui replace `yname`i'' = `ytilde`i''
            tempname ytilde_sq`i'
            qui gen double `ytilde_sq`i'' = (`ytilde`i'')^2
            // for D==0
            qui sum `ytilde_sq`i'' if `dvar'==0, meanonly
            mat `MSEy0'[1,`i'] = r(mean)
            local newyminmse0 = r(mean)
            if (`newyminmse0'<`yminmse0') {
                local yminmse0 = `newyminmse0'
                local yminmseid0 = `i'
            }
            // for D==1
            qui sum `ytilde_sq`i'' if `dvar'==1, meanonly
            mat `MSEy1'[1,`i'] = r(mean)
            local newyminmse1 = r(mean)
            if (`newyminmse1'<`yminmse1') {
                local yminmse1 = `newyminmse1'
                local yminmseid1 = `i'
            }
        }
        di " "
        if "`model'"=="interactive" di as res "Mean-squared error for y|D=0,X:"
        else  di as res "Mean-squared error for y|Z=0,X:"
        di _col(2) _c
        di _col(10) "Name" _c
        di _col(30) "Command" _c
        di _col(50) "MSPE"
        di "{hline 65}"
        forvalues i = 1(1)`yestn' {
            di _col(2) "`i'" _c
            di _col(10) "`yname`i''" _c
            di _col(30) "`ycmd`i''" _c
            if (`yminmseid0'==`i') {
                di _col(50) `MSEy0'[1,`i'] _c
                di "*"
            } 
            else {
                di _col(50) `MSEy0'[1,`i'] 
            }
        }
        di " "
        if "`model'"=="interactive" di as res "Mean-squared error for y|D=1,X:"
        else di as res "Mean-squared error for y|Z=1,X:"
        di _col(2) _c
        di _col(10) "Name" _c
        di _col(30) "Command" _c
        di _col(50) "MSPE"
        di "{hline 65}"
        forvalues i = 1(1)`yestn' {
            di _col(2) "`i'" _c
            di _col(10) "`yname`i''" _c
            di _col(30) "`ycmd`i''" _c
            if (`yminmseid1'==`i') {
                di _col(50) `MSEy1'[1,`i'] _c
                di "*"
            } 
            else {
                di _col(50) `MSEy1'[1,`i'] 
            }
        }
    }
    if ("`model'"=="partial"|"`model'"=="iv"|"`model'"=="interactive") {
        tempname MSEd
        mat `MSEd' = J(1,`destn',.)
        local dminmse = .
        local dminmseid = 1
        forvalues i = 1(1)`destn' {
            qui replace `dname`i'' = `dtilde`i''
            tempname dtilde_sq`i'
            qui gen double `dtilde_sq`i'' = (`dtilde`i'')^2
            qui sum `dtilde_sq`i'' , meanonly
            mat `MSEd'[1,`i'] = r(mean)
            local newdminmse = r(mean)
            if (`newdminmse'<`dminmse') {
                local dminmse = `newdminmse'
                local dminmseid = `i'
            }   
        }
        di " "
        di as res "Mean-squared error for d|X:"
        di _col(10) "Name" _c
        di _col(30) "Command" _c
        di _col(50) "MSPE"
        di "{hline 65}"
        forvalues i = 1(1)`destn' {
            di _col(2) "`i'" _c
            di _col(10) "`dname`i''" _c
            di _col(30) "`dcmd`i''" _c
            if (`dminmseid'==`i') {
                di _col(50) `MSEd'[1,`i'] _c
                di "*"
            }
            else {
                di _col(50) `MSEd'[1,`i']
            }
        }
    }
    if ("`model'"=="late") {
        tempname MSEd0 
        tempname MSEd1
        local dminmse0 = .
        local dminmse1 = .
        local dminmseid1 = 1
        local dminmseid0 = 1
        mat `MSEd0' = J(1,`destn',.)
        mat `MSEd1' = J(1,`destn',.)
        forvalues i = 1(1)`destn' {
            qui replace `dname`i'' = `dtilde`i''
            tempname dtilde_sq`i'
            qui gen double `dtilde_sq`i'' = (`dtilde`i'')^2
            // for Z==0
            qui sum `dtilde_sq`i'' if `zvar'==0, meanonly
            mat `MSEd0'[1,`i'] = r(mean)
            local newdminmse0 = r(mean)
            if (`newdminmse0'<`dminmse0') {
                local dminmse0 = `newdminmse0'
                local dminmseid0 = `i'
            }
            // for Z==1
            qui sum `dtilde_sq`i'' if `zvar'==0, meanonly
            mat `MSEd1'[1,`i'] = r(mean)
            local newdminmse1 = r(mean)
            if (`newdminmse1'<`dminmse0') {
                local dminmse1 = `newdminmse1'
                local dminmseid1 = `i'
            }
        }
        di " "
        di as res "Mean-squared error for d|Z=0,X:"
        di _col(2) _c
        di _col(10) "Name" _c
        di _col(30) "Command" _c
        di _col(50) "MSPE"
        di "{hline 65}"
        forvalues i = 1(1)`destn' {
            di _col(2) "`i'" _c
            di _col(10) "`dname`i''" _c
            di _col(30) "`dcmd`i''" _c
            if (`dminmseid0'==`i') {
                di _col(50) `MSEd0'[1,`i'] _c
                di "*"
            } 
            else {
                di _col(50) `MSEd0'[1,`i'] 
            }
        }
        di " "
        di as res "Mean-squared error for d|Z=1,X:"
        di _col(2) _c
        di _col(10) "Name" _c
        di _col(30) "Command" _c
        di _col(50) "MSPE"
        di "{hline 65}"
        forvalues i = 1(1)`destn' {
            di _col(2) "`i'" _c
            di _col(10) "`dname`i''" _c
            di _col(30) "`dcmd`i''" _c
            if (`dminmseid1'==`i') {
                di _col(50) `MSEd1'[1,`i'] _c
                di "*"
            } 
            else {
                di _col(50) `MSEd1'[1,`i'] 
            }
        }
    }
    if ("`model'"=="iv"|"`model'"=="late") {
        tempname MSEz
        mat `MSEz' = J(1,`zestn',.)
        local zminmse = .
        local zminmseid = 1
        forvalues i = 1(1)`zestn' {
            qui replace `zname`i'' = `ztilde`i''
            tempname ztilde_sq`i'
            qui gen double `ztilde_sq`i'' = (`ztilde`i'')^2
            qui sum `ztilde_sq`i'' , meanonly
            mat `MSEz'[1,`i'] = r(mean)
            local newzminmse = r(mean)
            if (`newzminmse'<`zminmse') {
                local zminmse = `newzminmse'
                local zminmseid = `i'
            }   
        }
        di " "
        di as res "Mean-squared error for z|X:"
        di _col(10) "Name" _c
        di _col(30) "Command" _c
        di _col(50) "MSPE"
        di "{hline 65}"
        forvalues i = 1(1)`zestn' {
            di _col(2) "`i'" _c
            di _col(10) "`zname`i''" _c
            di _col(30) "`zcmd`i''" _c
            if (`zminmseid'==`i') {
                di _col(50) `MSEz'[1,`i'] _c
                di "*"
            }
            else {
                di _col(50) `MSEz'[1,`i']
            }
        }
    }
    di as text "* indicates model with minimum MSE."

    *** save all 
    ereturn clear
    ereturn scalar crossfit = 1
    ereturn scalar yest = `yestn'
    ereturn scalar dest = `destn'
    ereturn local cmd ddml_crossfit
    ereturn local depvar `yvar'
    ereturn local dvar `dvar'
    ereturn local model "`model'"
    ereturn scalar crossfit = 1

    * return variable and command names
    forvalues i = 1(1)`yestn' {
        ereturn local y`i' `yname`i''
        ereturn local ycmd`i' `ycmd`i''
    }
    forvalues i = 1(1)`destn' {
        ereturn local d`i' `dname`i''
        ereturn local dcmd`i' `dcmd`i''
    }
    if ("`model'"=="iv"|"`model'"=="late") {
        ereturn scalar zest = `zestn'
        ereturn local zvar `zvar'
        forvalues i = 1(1)`zestn' {
            ereturn local z`i' `zname`i''
            ereturn local zcmd`i' `zcmd`i''
        }
    }

    * return MSE and opt-ID for Y
    if ("`model'"=="partial"|"`model'"=="iv") {
        ereturn scalar yoptid = `yminmseid'
        ereturn matrix ymse = `MSEy'
    }
    if ("`model'"=="late"|"`model'"=="interactive") {
        ereturn scalar y0optid = `yminmseid0'
        ereturn scalar y1optid = `yminmseid1'
        ereturn matrix y0mse = `MSEy0'
        ereturn matrix y1mse = `MSEy1'
    }
    *** return MSE and opt-ID for D
    if ("`model'"=="partial"|"`model'"=="iv"|"`model'"=="interactive") {
        ereturn scalar doptid = `dminmseid'
        ereturn matrix dmse = `MSEd'
    } 
    if ("`model'"=="late") {
        ereturn scalar d1optid = `dminmseid1'
        ereturn scalar d0optid = `dminmseid0'
        ereturn matrix d1mse = `MSEd1'
        ereturn matrix d0mse = `MSEd0'
    }
    *** return MSE and opt-ID for Z
    if ("`model'"=="late"|"`model'"=="iv") {
        ereturn scalar zoptid = `zminmseid'
        ereturn matrix zmse = `MSEz'
    }
end