
program _ddml_ate_late, eclass

	version 13

	syntax [anything] [if] [in] ,				///
						[						///
							yvar(varname)		///
							dvar(varname)		///
							zvar(varname)		///
							y0tilde(name)		///
							y1tilde(name)		///
							dtilde(name)		///
							d0tilde(name)		///
							d1tilde(name)		///
							ztilde(name)		///
							mname(name)			///
							rep(integer 1)		///
							touse(varname)		///
						]
		
	if ~replay() {
		
		mata: st_local("model",`mname'.model)
		
		tempname b
		tempname V 
		
		// add suffixes
		local y0_m		`y0tilde'0_`rep'
		local y1_m		`y1tilde'1_`rep'
		local d_m		`dtilde'_`rep'
		local d0_m		`d0tilde'0_`rep'
		local d1_m		`d1tilde'1_`rep'
		local z_m		`ztilde'_`rep'
		
		if "`model'"=="interactive" {
			markout `touse' `yvar' `dvar' `y0_m' `y1_m' `d_m'
			mata: ATE("`yvar'","`dvar'","`y0_m'", "`y1_m'", "`d_m'","`touse'","`b'","`V'")
		}
		else {
			markout `touse' `yvar' `dvar' `zvar' `y0_m' `y1_m' `d0_m' `d1_m' `z_m'
			mata: LATE("`yvar'","`dvar'","`zvar'","`y0_m'", "`y1_m'", "`d0_m'","`d1_m'","`z_m'","`touse'","`b'","`V'")
		}
		
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
		if "`model'"=="interactive" {
			ereturn local d		`d_m'
		}
		else {
			ereturn local d0	`d0_m'
			ereturn local d1	`d1_m'
			ereturn local z		`z_m'
		}
		
		// additional estimation results
		ereturn scalar resample = `rep'
		tempname eqn
		mata: `eqn' = init_eStruct()
		// Y eqn results
		mata: `eqn' = (`mname'.eqnAA).get("`yvar'")
		mata: st_numscalar("e(`y0tilde'_mse0)",return_result_item(`eqn',"`y0tilde'","MSE0","`rep'"))
		mata: st_numscalar("e(`y1tilde'_mse1)",return_result_item(`eqn',"`y1tilde'","MSE1","`rep'"))
		mata: st_matrix("e(`y0tilde'_mse0_folds)",return_result_item(`eqn',"`y0tilde'","MSE0_folds","`rep'"))
		mata: st_matrix("e(`y1tilde'_mse1_folds)",return_result_item(`eqn',"`y1tilde'","MSE1_folds","`rep'"))
		mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`y0tilde'","cmd")=="pystacked"))
		if `pyswflag' {
			mata: st_matrix("e(`y0tilde'_pysw)", mean(return_result_item(`eqn',"`y0tilde'","stack_weights0","`rep'")'))
		}
		mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`y1tilde'","cmd")=="pystacked"))
		if `pyswflag' {
			mata: st_matrix("e(`y1tilde'_pysw)", mean(return_result_item(`eqn',"`y1tilde'","stack_weights1","`rep'")'))
		}
		if "`model'"=="interactive" {
			mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
			mata: st_numscalar("e(`dtilde'_mse)",return_result_item(`eqn',"`dtilde'","MSE","`rep'"))
			mata: st_matrix("e(`dtilde'_mse_folds)",return_result_item(`eqn',"`dtilde'","MSE_folds","`rep'"))
			mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`dtilde'","cmd")=="pystacked"))
			if `pyswflag' {
				mata: st_matrix("e(`dtilde'_pysw)", mean(return_result_item(`eqn',"`dtilde'","stack_weights","`rep'")'))
			}
		}
		else {
			mata: `eqn' = (`mname'.eqnAA).get("`dvar'")
			mata: st_numscalar("e(`d0tilde'_mse0)",return_result_item(`eqn',"`d0tilde'","MSE0","`rep'"))
			mata: st_numscalar("e(`d1tilde'_mse1)",return_result_item(`eqn',"`d1tilde'","MSE1","`rep'"))
			mata: st_matrix("e(`d0tilde'_mse0_folds)",return_result_item(`eqn',"`d0tilde'","MSE0_folds","`rep'"))
			mata: st_matrix("e(`d1tilde'_mse1_folds)",return_result_item(`eqn',"`d1tilde'","MSE1_folds","`rep'"))
			mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`d0tilde'","cmd")=="pystacked"))
			if `pyswflag' {
				mata: st_matrix("e(`d0tilde'_pysw)", mean(return_result_item(`eqn',"`d0tilde'","stack_weights0","`rep'")'))
			}
			mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`d1tilde'","cmd")=="pystacked"))
			if `pyswflag' {
				mata: st_matrix("e(`d0tilde'_pysw)", mean(return_result_item(`eqn',"`d1tilde'","stack_weights1","`rep'")'))
			}
			mata: `eqn' = (`mname'.eqnAA).get("`zvar'")
			mata: st_numscalar("e(`ztilde'_mse)",return_result_item(`eqn',"`ztilde'","MSE","`rep'"))
			mata: st_matrix("e(`ztilde'_mse_folds)",return_result_item(`eqn',"`ztilde'","MSE_folds","`rep'"))
			mata: st_local("pyswflag",strofreal(return_learner_item(`eqn',"`ztilde'","cmd")=="pystacked"))
			if `pyswflag' {
				mata: st_matrix("e(`ztilde'_pysw)", mean(return_result_item(`eqn',"`ztilde'","stack_weights","`rep'")'))
			}
		}
	}
	
	di as text "y-E[y|X,D=0]" _col(14) "= " as res "`e(y0)'" _c
	di as text _col(52) "Number of obs   =" _col(70) as res %9.0f `e(N)'
	di as text "y-E[y|X,D=1]" _col(14) "= " as res "`e(y1)'"
	if "`e(model)'"=="interactive" {
		di as text "D-E[D|X]" _col(14)  "= " as res "`e(d)'"
	}
	else {
		di as text "D-E[D|X,Z=0]" _col(14)  "= " as res "`e(d0)'"
		di as text "D-E[D|X,Z=1]" _col(14)  "= " as res "`e(d1)'"
		di as text "Z-E[Z|X]" _col(14)  "= " as res "`e(z)'"
	}
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

void LATE(  string scalar yvar,      // Y
            string scalar dvar,      // D
            string scalar zvar,      // Z
            string scalar y0tilde,   // E[Y|X,Z=0]
            string scalar y1tilde,   // E[Y|X,Z=1]
            string scalar d0tilde,   // E[D|X,Z=0]
            string scalar d1tilde,   // E[D|X,Z=1]
            string scalar ztilde,    // E[Z|X]
            string scalar sample,    // sample
            string scalar outlate,   // output: name of matrix to store b
            string scalar outlatese  // output: name of matrix to store V
            )
{
    st_view(my_z0x,.,y0tilde,sample)
    st_view(my_z1x,.,y1tilde,sample)
    st_view(md_z0x,.,d0tilde,sample)
    st_view(md_z1x,.,d1tilde,sample)
    st_view(mz_x,.,ztilde,sample)
    st_view(d,.,dvar,sample)
    st_view(y,.,yvar,sample)
    st_view(z,.,zvar,sample)
    n = rows(y)
    late =  mean( z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x ) /
                mean( z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x ) 
    late_se = variance(( z :* (y :- my_z1x) :/ mz_x :-  ((1 :- z) :* (y :- my_z0x) :/ (1 :- mz_x)) :+ my_z1x :- my_z0x ) :/ 
               mean( z :* (d :- md_z1x) :/ mz_x :-  ((1 :- z) :* (d :- md_z0x) :/ (1 :- mz_x)) :+ md_z1x :- md_z0x )) / n
    st_numscalar("r(N)",n)
    st_matrix(outlate,late)
    st_matrix(outlatese,late_se)
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

LATE <- function(y, d, z, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)
{
  return( mean( z * (y - my_z1x) / mz_x -  ((1 - z) * (y - my_z0x) / (1 - mz_x)) + my_z1x - my_z0x ) / 
            mean( z * (d - md_z1x) / mz_x -  ((1 - z) * (d - md_z0x) / (1 - mz_x)) + md_z1x - md_z0x ) );
}
SE.LATE <- function(y, d, z, my_z1x, my_z0x, mz_x, md_z1x, md_z0x)
{
  return( sd(( z * (y - my_z1x) / mz_x -  ((1 - z) * (y - my_z0x) / (1 - mz_x)) + my_z1x - my_z0x ) / 
               mean( z * (d - md_z1x) / mz_x -  ((1 - z) * (d - md_z0x) / (1 - mz_x)) + md_z1x - md_z0x )) / sqrt(length(y)) );
}

*/