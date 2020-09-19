*! ddml v0.1.4 (19 sep 2020)

program ddml, eclass

    version 13

    syntax [anything] , [*]

    local subcmd: word 1 of `anything'
    
    *** initialize new estimation
    if "`subcmd'"=="init" {
        local model: word 2 of `anything'
        if ("`model'" != "partial") {
            di as err "no or wrong model specified." _c
            di as err "currently only 'ddml init partial' is supported"
            exit 1
        }
        ereturn clear
        ereturn local cmd "ddml_init"
        ereturn local model "`model'"
        ereturn scalar yest = 0 // number of y-equation estimations
        ereturn scalar dest = 0 // number of d-equation estimations
        ereturn scalar crossfit = 0
    } 
    else {
        if ("`e(cmd)'"!="ddml_init") {
            di as err "no ddml object found." _c
            di as err "you first need to initialize using 'ddml init'"
            exit 1
        }
    }

    *** add y-equation estimation
    if "`subcmd'"=="yeq" {
        local estname: word 2 of `anything'
        gettoken estname: estname, parse(":")
        gettoken left theeq: anything, parse(":")
        gettoken left theeq: theeq, parse(":") 
        gettoken theeq: theeq, parse("(") match(paren)
        ereturn local y`=`e(yest)'+1' "y_`estname'"
        ereturn local ycmd`=`e(yest)'+1' "`theeq'"
        ereturn scalar yest = `e(yest)'+1
    }

    *** add d-equation estimation
    if "`subcmd'"=="deq" {
        local estname: word 2 of `anything'
        gettoken estname: estname, parse(":")
        gettoken left theeq: anything, parse(":")
        gettoken left theeq: theeq, parse(":") 
        gettoken theeq: theeq, parse("(") match(paren)
        ereturn local d`=`e(dest)'+1' "d_`estname'"
        ereturn local dcmd`=`e(dest)'+1' "`theeq'"
        ereturn scalar dest = `e(dest)'+1
    }

    *** cross-fitting
    if "`subcmd'" =="crossfit" {
        _ddml_crossfit, `options'
    }

    *** estimate
    if "`subcmd'" =="estimate" {
        if (`e(crossfit)'==0) {
            di as err "you first need to do 'ddml crossfit'"
            exit 1
        }
        _ddml_estimate, `options'
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
							drclass ]

    *** check there are enough equations
    local yestn = `e(yest)'
    local destn = `e(dest)' 
    if (`yestn'==0 | `destn'==0) {
        di as err "insufficient equations/models specified." _c 
        di as err "we need at least one equation for y and one for d."
        exit 1
    }
    
    *** extract yvar, dvar, xvars and cmds; sense-checking
    local ycmd1 = e(ycmd1)
    local dcmd1 = e(dcmd1)
    local yvar: word 2 of `ycmd1'
    local dvar: word 2 of `dcmd1'
    forvalues i = 1(1)`yestn' {
        local ycmd`i' = e(ycmd`i')
        local yname`i' = e(y`i')
        local yvar`i': word 2 of `ycmd`i'' 
        if "`yvar`i''" != "`yvar'" {
            di as err "inconsistent d-variables: `yvar', `yvar`i''"
            exit 1
        }
    }
    forvalues i = 1(1)`destn' {
        local dcmd`i' = e(dcmd`i')
        local dname`i' = e(d`i')
        local dvar`i': word 2 of `dcmd`i'' 
        if "`dvar`i''" != "`dvar'" {
            di as err "inconsistent d-variables: `dvar', `dvar`i''"
            exit 1
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

    *** parse y command(s)
    forvalues i = 1(1)`yestn' {
        local 0 `ycmd`i''
        syntax [anything] [if] [in], [*]
        if "`if'`in'"!="" {
            di as err "if and in not allowed in response equation"
            error 198
        }
        local ycmdline`i' `anything'
        local ycmd`i': word 1 of `ycmdline`i''
        local yxvar`i': list ycmdline`i' - ycmd`i' 
        local yxvar`i': list yxvar`i' - yvar  
        qui ds `yxvar`i''
        local yxvar`i' = r(varlist)
        local yopts`i' `options'
    }

	*** parse D command(s)
    forvalues i = 1(1)`destn' {
        local 0 `dcmd`i''
        syntax [anything] [if] [in], [*]
        if "`if'`in'"!="" {
            di as err "if and in not allowed in treatment equation"
            error 198
        }
        local dcmdline`i' `anything'
        local dcmd`i': word 1 of `dcmdline`i''
        local dxvar`i': list dcmdline`i' - dcmd`i' 
        local dxvar`i': list dxvar`i' - dvar  	
        qui ds `dxvar`i'' 
        local dxvar`i' = r(varlist)
        local dopts`i' `options'
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

    *** show locals (debug-mode)
    if ("`debug'"!="") {
        di "---debug mode----show locals---------------------"
        di "yvar: `yvar'"
        di "dvar: `dvar'"
        di "-------------------------------------------------"
        forvalues i = 1(1)`destn' {
            di "dcmd`i': `dcmd`i''"
            di "dcmdline`i': `dcmdline`i''"
            di "dxvar`i': `dxvar`i''"
            di "dopts`i': `dopts`i''"
        }
        di "-----------------------------------"
        forvalues i = 1(1)`yestn' {
            di "ycmd`i': `ycmd`i''"
            di "ycmdline`i': `ycmdline`i''"
            di "yxvar`i': `yxvar`i''"
            di "yopts`i': `yopts`i''"
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
			forvalues i = 1(1)`yestn' {
                cap drop `ytilde_temp`i''
				`ycmdline`i'' if `kid'!=`k', `yopts`i''
				predict `ytilde_temp`i'' if `kid'==`k' 
				replace `ytilde`i'' = `yvar' - `ytilde_temp`i'' if `kid'==`k'
			}
			// yrclass-option. currently disabled.
            // this is for commands that don't support post-est predict.
            //else {
			//	cap drop `outsample'
			//	cap drop `estsample'
			//	gen `estsample' = (`kid'!=`k')
			//	gen `outsample' = 1-`estsample'	
			//	`ycmdline', outsample(`outsample') estsample(`estsample') rname(`ytilde_temp')
			//	replace `ytilde' = `ytilde_temp' if `kid'==`k'
			//}
			
			** d equation(s)
            forvalues i = 1(1)`destn' {
                cap drop `dtilde_temp`i''
                `dcmdline`i'' if `kid'!=`k', `dopts`i''
                predict `dtilde_temp`i'' if `kid'==`k'   
                replace `dtilde`i'' = `dvar' - `dtilde_temp`i'' if `kid'==`k'
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
    * y equation(s)
    tempname MSEy
    tempname MSEd
    mat `MSEy' = J(1,`yestn',.)
    mat `MSEd' = J(1,`destn',.)
    local yminmse = .
    local dminmse = .
	forvalues i = 1(1)`yestn' {
        gen double _`yname`i'' = `ytilde`i''
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
    forvalues i = 1(1)`destn' {
        gen double _`dname`i'' = `dtilde`i''
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

    *** display mse
    di " "
    di as res "Mean-squared error for y:"
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
    di " "
    di as res "Mean-squared error for d:"
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
    di as text "* indicates model with minimum MSE."

    *** save all 
    ereturn clear
    ereturn scalar crossfit = 1
    ereturn scalar yest = `yestn'
    ereturn scalar dest = `destn'
    forvalues i = 1(1)`yestn' {
        ereturn local y`i' `yname`i''
        ereturn local ycmd`i' `ycmdline`i''
    }
    forvalues i = 1(1)`destn' {
        ereturn local d`i' `dname`i''
        ereturn local dcmd`i' `dcmdline`i''
    }
    ereturn local cmd ddml_init
    ereturn local model "partial"
    ereturn scalar crossfit = 1
    ereturn scalar yoptid = `yminmseid'
    ereturn scalar doptid = `dminmseid'
    ereturn matrix ymse = `MSEy'
    ereturn matrix dmse = `MSEd'
    ereturn local depvar `yvar'
    ereturn local dvar `dvar'


end

*** ddml estimation
program _ddml_estimate, eclass sortpreserve

    syntax [anything] [if] [in] , /// 
								[  ///
                                robust ///
                                show(string) /// dertermined which to post
                                * ]

    *** set defaults
    //if ("`show'"=="") {
    //    local show opt
    //}
	
    *** save everything that is needed in locals
    local yestn = e(yest)
    local destn = e(dest)
    local yoptid = e(yoptid)
    local doptid = e(doptid)
    local yvar = e(depvar)
    local dvar = e(dvar)
    
    *** do estimation
    forvalues i = 1(1)`yestn' {
        local ytilde`i' _`e(y`i')'
    }
    forvalues i = 1(1)`destn' {
        local dtilde`i' _`e(d`i')'
    }
    if ("`show'"=="all") {
        forvalues i = 1(1)`yestn' {
            forvalues j = 1(1)`destn' {
                if (`i'!=`yoptid' | `j'==`doptid') {
                    di as text "DML with `ytilde`i'' and `dtilde`j'':"
                    qui reg `ytilde`i'' `dtilde`j'', nocons `robust' noheader
                    // display
                    tempname b
                    tempname V 
                    mat `b' = e(b)
                    mat `V' = e(V)
                    matrix colnames `b' = "`dvar'"
                    matrix rownames `b' = "`yvar'"
                    matrix colnames `V' = "`dvar'"
                    matrix rownames `V' = "`dvar'"
                    ereturn clear
                    ereturn post `b' `V' 
                    ereturn display

                }
            }
        }
    }
    *** estimate best model
    di as res "Optimal model: DML with `ytilde`yoptid'' and `dtilde`doptid'':"
    qui reg `ytilde`yoptid'' `dtilde`doptid'', nocons `robust' noheader
    // display
	tempname b
	tempname V 
	mat `b' = e(b)
	mat `V' = e(V)
	matrix colnames `b' = "`dvar'"
	matrix rownames `b' = "`yvar'"
 	matrix colnames `V' = "`dvar'"
	matrix rownames `V' = "`dvar'"
	ereturn clear
	ereturn post `b' `V' 
	ereturn display

end