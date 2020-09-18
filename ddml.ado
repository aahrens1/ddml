** version 5 july 2019

program ddml, eclass sortpreserve

	version 13
	
	syntax [anything] [if] [in] , /// 
								[ kfolds(integer 2)  ///
								NOIsily ///
								debug /// 
								Robust ///
								TABFold ///
								yrclass drclass ]
	
	if ("`noisily'"=="") {
		local qui quietly
	}
	//
	
	dmlparse `anything'
	local ycmd = r(ycmd)
	local dcmd = r(dcmd)
	if ("`debug'"!="") {
		di "y equation: `ycmd'"
		di "d equation: `dcmd'"
	}
	//
	
	*** parse y cmd
	local 0 `ycmd'
	syntax [anything] [if] [in], [*]
	if "`if'`in'"!="" {
		di as err "if and in not allowed in response equation"
		error 198
	}
	local ycmdline `anything'	
	local ycmd: word 1 of `ycmdline'
	local yvar: word 2 of `ycmdline'
	local yxvar: list ycmdline - ycmd  
	local yxvar: list yxvar - yvar  
	qui ds `yxvar'
	local yxvar = r(varlist)
	local yopts `options'
	
	*** parse D cmd 
	local 0 `dcmd'
	syntax [anything] [if] [in], [*]
	if "`if'`in'"!="" {
		di as err "if and in not allowed in treatment equation"
		error 198
	}
	local dcmdline `anything'
	local dcmd: word 1 of `dcmdline'
	local dvar: word 2 of `dcmdline'
	local dxvar: list dcmdline - dcmd  
	local dxvar: list dxvar - dvar  	
	qui ds `dxvar' 
	local dxvar = r(varlist)
	local dopts `options'
	
	*** make sure list of controls is the same
	local notind : list dxvar - yxvar
	local notiny : list yxvar - dxvar
	if ("`notind'"!="" | "`notiny'"!="") {
		di as err "list of controls in response and treatment equation not the same"
		error 198
	}
	//
	
	*** gen folds
	tempvar kid uni cuni
	gen double `uni' = runiform()
	cumul `uni', gen(`cuni')
	gen `kid' =ceil(`kfolds'*`cuni')
	if ("`tabfold'"!="") {
		tab `kid'
	}
	//
	
	tempvar ytilde dtilde ytilde_temp dtilde_temp
	tempname bhat_k se_k
	mat `bhat_k' = J(`kfolds',1,.)
	mat `se_k' = J(`kfolds',1,.)
	qui gen double `ytilde'=.
	qui gen double `dtilde'=.
	
	tempname outsample estsample
		
	di "Cross-fitting fold " _c
	forvalues k = 1(1)`kfolds' {
	
		if (`k'==`kfolds') {
			di as text "`k'"
		}
		else {
			di as text "`k' " _c
		}

		cap drop `ytilde_temp' `dtilde_temp'
		
		// ML is applied to I^c sample (all data ex partition k)
	
		qui {
		
			** y equation
			if ("`yrclass'"=="") {
				`ycmdline' if `kid'!=`k', `yopts'
				predict double `ytilde_temp' if `kid'==`k', r
				replace `ytilde' = `ytilde_temp' if `kid'==`k'
			}
			else {
				cap drop `outsample'
				cap drop `estsample'
				gen `estsample' = (`kid'!=`k')
				gen `outsample' = 1-`estsample'	
				`ycmdline', outsample(`outsample') estsample(`estsample') rname(`ytilde_temp')
				replace `ytilde' = `ytilde_temp' if `kid'==`k'
			}
			
			** d equation
			`dcmdline' if `kid'!=`k', `dopts'
			predict double `dtilde_temp' if `kid'==`k', r
			replace `dtilde' = `dtilde_temp' if `kid'==`k'
		}
		
		`qui' reg `ytilde' `dtilde' if `kid'==`k', `robust' nocons
		qui mat `bhat_k'[`k',1]=_b[`dtilde']
		qui mat `se_k'[`k',1]=_se[`dtilde']
			
	}	
	mat list `bhat_k'
	`qui' mat list `se_k'
	
	*** DML 1: calculate point estimates by fold,
	*** then take average across K estimations
	tempname bhat
	mata: `bhat' = st_matrix("`bhat_k'")
	mata: `bhat' = mean(`bhat')
	
	tempname se
	mata: `se' = st_matrix("`se_k'")
	mata: `se' = sqrt(sum((`se'):^2)/(`kfolds'^2))
	
	
	di as text "point estimate DML1="
	mata: `bhat'
	di as text "standard error DML1="
	mata: `se'
	
	*** DML 2 (recommended method): joint OLS estimation on full sample
	*** see Remark 3.1
	di as text "DML 2:"
	reg `ytilde' `dtilde', nocons `robust'
	
	***
	tempname ytilde_sq
	qui gen double `ytilde_sq' = (`ytilde')^2
	qui sum `ytilde_sq' , meanonly
	di as text "MSE of ytilde"
	di r(mean)
	
	tempname dtilde_sq
	qui gen double `dtilde_sq' = (`dtilde')^2
	qui sum `dtilde_sq' , meanonly
	di as text "MSE of dtilde"
	di r(mean)
	
end
	

