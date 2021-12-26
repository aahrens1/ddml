
program _ddml_ate_late, eclass

	version 13

	syntax [anything] [if] [in] ,				///
						[						///
							yvar(varname)		///
							dvar(varname)		///
							y0tilde(name)		///
							y1tilde(name)		///
							dtilde(name)		///
							mname(name)			///
							rep(integer 1)		///
							touse(varname)		///
						]
		
	if ~replay() {
		
		mata: st_local("model",`mname'.model)
		
		tempname b
		tempname V 
		
		// add resample suffixes and estimate
		local y0_m		`y0tilde'_`rep'
		local y1_m		`y1tilde'_`rep'
		local d_m		`dtilde'_`rep'
		
		markout `touse' `yvar' `dvar' `y0_m' `y1_m' `d_m'
		
		mata: ATE("`yvar'","`dvar'","`y0_m'", "`y1_m'", "`d_m'","`touse'","`b'","`V'")
		
		local N = `r(N)'
		matrix colnames `b' = "`dvar'"
		matrix rownames `b' = "`yvar'"
		matrix colnames `V' = "`dvar'"
		matrix rownames `V' = "`dvar'"
		ereturn clear
		ereturn post `b' `V', depname(`yvar') obs(`N')
		
		ereturn local cmd		_ddml_ate_late
		ereturn local model		`model'
		ereturn local y0		`y0_m'
		ereturn local y1		`y1_m'
		ereturn local d			`d_m'
	}
	
	di as text "y-E[y|X,D=0]" _col(14) "= " as res "`e(y0)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	di as text "y-E[y|X,D=1]" _col(14) "= " as res "`e(y1)'"
	di as text "D-E[D|X]" _col(14)  "= " as res "`e(d)'"
	ereturn display
end


********************************************************************************
*** Mata section															 ***
********************************************************************************

mata:

void ATE(   string scalar yvar,	 // Y
			string scalar dvar,	 // D
			string scalar y0tilde,  // E[Y|X,D=0]
			string scalar y1tilde,  // E[Y|X,D=1]
			string scalar dtilde,   // E[D|X]
			string scalar sample,   // sample
			string scalar outate,   // output: name of matrix to store b
			string scalar outatese  // output: name of matrix to store V
			)
{
	st_view(my_d0x,.,y0tilde,sample)
	st_view(my_d1x,.,y1tilde,sample)
	st_view(md_x,.,dtilde,sample)
	st_view(d,.,dvar,sample)
	st_view(y,.,yvar,sample)

	n = rows(y)

	te
	te  = (d :* (y :- my_d1x) :/ md_x) :-  ((1 :- d) :* (y :- my_d0x) :/ (1 :- md_x)) :+ my_d1x :- my_d0x  
	ate = mean(te)
	ate_V =  variance(te)/n

	st_numscalar("r(N)",n)
	st_matrix(outate,ate)
	st_matrix(outatese,ate_V)
}

end

/*

ATE <- function(y, d, my_d1x, my_d0x, md_x)
{
  return( mean( (d * (y - my_d1x) / md_x) -  ((1 - d) * (y - my_d0x) / (1 - md_x)) + my_d1x - my_d0x ) );
}

SE.ATE <- function(y, d, my_d1x, my_d0x, md_x)
{
  return( sd( (d * (y - my_d1x) / md_x) -  ((1 - d) * (y - my_d0x) / (1 - md_x)) + my_d1x - my_d0x )/sqrt(length(y)) );
}

