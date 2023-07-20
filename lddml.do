* do file to create ddml mata library lddml.mlib.
* modify code below to control destination.
* mlib includes the mata utility whichddml().

* Find the right destination folder for lddml.mlib.
* First, save current folder location (will return to it later).
local pwd `c(pwd)'
* Now change to the right folder. Will depend on user.
* This is where the mlib will be created.
* ms:
cap cd "C:\LocalStore\ecomes\Documents\GitHub\ddml"

* Locals used in whichddml; set when compiled
local stata_version `c(stata_version)'
local born_date `c(born_date)'
local current_date `c(current_date)'

version 16.0
mata:
mata clear

void whichddml()
{
""
"Mata library for ddml and related programs,"
"compiled `current_date' under Stata " + "`stata_version'" + " born " + "`born_date'."
"authors AA/MS"
st_sclear()
st_global("r(stata_born_date)","`born_date'")
st_global("r(stata_version)","`stata_version'")
st_global("r(stata_compiled_date)","`current_date")
}

////////////////////////////////////////////////////////////////////////////////

// uniquely identified by vname = dep var in the equation
// equation structure: one for E[y|x], E[d|x] etc
struct eStruct {
	string scalar					vname			// name of variable to be orthogonalized
	string scalar					etype			// will be Y, D or Z
	real matrix						vtlist			// list of orthogonalized (learner) variables
	string scalar					shortstack		// name of shortstack variable
	string scalar					poolstack		// name of poolstack variable
	real scalar						nlearners		// number of learners
	real scalar						lieflag			// =1 if LIE spec with two estimation strings etc.
	real scalar						ateflag			// =1 if treatment variable in ATE/LATE
	real scalar						pystackedmulti	// =#learners if pystacked with multiple learners
	class AssociativeArray scalar	lrnAA			// AssociativeArray with all learners //
													// (keys=vtilde,object)
	class AssociativeArray scalar	resAA			// AssociativeArray with all learner results //
													// (keys=vtilde,object,rep)

}

// should perhaps make vname a required argument for this function
struct eStruct init_eStruct()
{
	struct eStruct scalar	m

	m.vname				= ""
	m.etype				= ""
	m.vtlist			= J(1,0,"")
	m.shortstack		= ""
	m.poolstack			= ""
	m.nlearners			= 0
	m.lieflag			= 0
	m.ateflag			= 0
	m.pystackedmulti	= 0
	
	(m.lrnAA).reinit("string",2)
	(m.resAA).reinit("string",3)
	(m.lrnAA).notfound(NULL)
	(m.resAA).notfound(NULL)
	
	return(m)
}

// clear results from eStruct
void clear_equation_results(struct eStruct e)
{
	(e.resAA).reinit("string",3)
	e.shortstack	= ""
	e.poolstack		= ""
}

// ddml model structure
struct mStruct {
	string scalar					model			// model; partial, iv, late, etc
	real scalar						nreps			// number of resamplings
	real scalar						ncombos			// number of possible specifications (=0 if not yet estimated)
	real scalar						kfolds			// number of crossfitting folds
	class AssociativeArray scalar	eqnAA			// AssociativeArray with all equations
	class AssociativeArray scalar	estAA			// AA wth all estimation results
	string scalar					nameY			// dependent variable 
	string colvector				nameD			// treatment variable(s)
	string colvector				nameZ			// instrument(s)
	string scalar					fclustvar		// name of fold cluster variable (="" if none)
	real scalar						stdflag			// flag for standard stacking (default=1)
	real scalar						ssflag			// flag for short-stacking
	real scalar						psflag			// flag for pooled-stacking
	string scalar					strDatavars		// string with expanded names of Stata variables
	real matrix						matDatavars		// matrix with values of Stata variables
	real scalar						crossfitted   	// =number of reps for which crossfitting; 0 if not
	real scalar						estimated		// =1 if estimation has been done; 0 if not
	real scalar						ycounter		// counter for default y learners
	real scalar						dcounter		// counter for default d learners
	real scalar						zcounter		// counter for default z learners
	real scalar						allpystackedmulti	// =1 if all eqns to use pystacked-specific code/features
}

struct mStruct init_mStruct()
{

	struct mStruct scalar	m
	
	m.model				= ""
	m.nreps				= 0
	m.ncombos			= 0
	m.kfolds			= 0
	m.nameY				= ""
	m.nameD				= J(1,0,"")
	m.nameZ				= J(1,0,"")
	m.stdflag			= 1			// will be set =0 if pystacked uses voting instead
	m.ssflag			= 0
	m.psflag			= 0
	m.fclustvar			= ""
	m.strDatavars		= ""
	m.matDatavars		= J(0,0,.)
	m.crossfitted		= 0
	m.estimated			= 0
	m.allpystackedmulti	= 0
	
	(m.eqnAA).reinit("string",1)
	(m.eqnAA).notfound(NULL)
	
	(m.estAA).reinit("string",2)
	(m.estAA).notfound(NULL)
	
	// initialize counters used for default learner names
	m.ycounter		= 1
	m.dcounter		= 1
	m.zcounter		= 1
	
	return(m)
}

// clear crossfitting and estimation results from all eStructs in mStruct
void clear_model_results(struct mStruct m)
{
	// if not crossfitted, nothing to clear
	if (m.crossfitted>0) {
		struct eStruct scalar			e
		
		// clear results for each equation
		eqnlist = (m.nameY, m.nameD, m.nameZ)
		for (i=1; i<=cols(eqnlist); i++) {
			clear_equation_results((m.eqnAA).get(eqnlist[i]))
		}
		
		// also clear:
		m.ncombos		= 0
		m.crossfitted	= 0
		m.estimated		= 0
		m.ssflag		= 0
		m.psflag		= 0
		m.strDatavars	= ""
		m.matDatavars	= J(0,0,.)
		(m.estAA).reinit("string",2)
		(m.estAA).notfound(NULL)
	}
}

// clear only estimation results in mStruct
void clear_model_estimation(struct mStruct m)
{
	m.estimated		= 0
	m.ncombos		= 0
	(m.estAA).reinit("string",2)
	(m.estAA).notfound(NULL)
}

// add item about learner to eStruct
void add_learner_item(				struct eStruct e,
									string scalar key1,
									string scalar key2,
									transmorphic scalar s
									)
{
	(e.lrnAA).put((key1,key2),s)
}

// retrieve item about learner from eStruct
transmorphic return_learner_item(	struct eStruct e,
									string scalar key1,
									string scalar key2
									)
{
	r = (e.lrnAA).get((key1,key2))
	if (r==NULL) {
		s = sprintf("\nerror - learner item key1(%s), key2(%s) not found\n",key1,key2)
		_error(310,s)
	}
	else {
		return(r)
	}
}

// add result from learner/resample to eStruct
void add_result_item(				struct eStruct e,
									string scalar key1,
									string scalar key2,
									string scalar rep,
									s						//  s can be anything
									)
{
	(e.resAA).put((key1,key2,rep),s)
}

// retrieve result from learner/resample from eStruct
transmorphic return_result_item(	struct eStruct e,
									string scalar key1,
									string scalar key2,
									string scalar rep
									)
{
	r = (e.resAA).get((key1,key2,rep))
	if (r==NULL) {
		s = printf("\nerror - learner item key1(%s), key2(%s), rep(%s) not found\n",key1,key2,rep)
		_error(310,s)
	}
	else {
		return(r)
	}
}


// return number of possible combinations of specs
transmorphic return_ncombos(struct mStruct m)
{
   	struct eStruct scalar	eqn

	eqn = (m.eqnAA).get(m.nameY)
	comboY = cols(eqn.vtlist)

	comboD = 1
	for (i=1;i<=cols(m.nameD);i++) {
		eqn = (m.eqnAA).get(m.nameD[1,i])
		numlrnD = cols(eqn.vtlist)
		comboD = comboD * numlrnD
	}
	comboZ = 1
	for (i=1;i<=cols(m.nameZ);i++) {
		eqn = (m.eqnAA).get(m.nameZ[1,i])
		numlrnZ = cols(eqn.vtlist)
		comboZ = comboZ * numlrnZ
	}
	
	if (m.model=="interactive" | m.model=="interactiveiv") {
	    comboY = comboY^2
	}
	if (m.model=="interactiveiv" | m.model=="fiv") {
	    comboD = comboD^2
	}
	
	// possible specs:
	return(comboY * comboD * comboZ)
}


transmorphic check_spec(							///
							struct mStruct m,		///
							string scalar spec,		///
							string scalar rep		///
						)
{
	B = (m.estAA).get((spec,rep))
	if (B==NULL) {
		s = sprintf("\nerror - spec(%s), rep(%s) not found\n",spec,rep)
		_error(310,s)
	}

}

transmorphic model_chars(struct mStruct m)
{
   	struct eStruct scalar	e
	
	// will be set to zero if LIE or any eqn is not pystackedmulti
	allpystackedmulti				= 1

	numeqnD							= cols(m.nameD)
	numeqnZ							= cols(m.nameZ)
	
	st_global("r(model)",			m.model)
	st_numscalar("r(crossfitted)",	m.crossfitted)
	st_numscalar("r(estimated)",	m.estimated)
	st_numscalar("r(ncombos)",		m.ncombos)
	st_numscalar("r(kfolds)",		m.kfolds)
	st_numscalar("r(nreps)",		m.nreps)
	st_global("r(nameY)",			m.nameY)
	st_global("r(nameD)",			invtokens(m.nameD))
	st_global("r(nameZ)",			invtokens(m.nameZ))
	st_numscalar("r(numeqnD)",		numeqnD)
	st_numscalar("r(numeqnZ)",		numeqnZ)
	
	// LIE model not supported by pystacked-specific code
	if (m.model=="fiv") {
		allpystackedmulti			= 0
	}

	// Y eqns
	e = (m.eqnAA).get(m.nameY)
	shortstack						= e.shortstack
	poolstack						= e.poolstack
	numlrnY							= cols(e.vtlist)
	vtlistY							= invtokens(e.vtlist)
	if (shortstack~="") {
	    vtlistY = vtlistY + " " + shortstack + "_ss"
	}
	if (poolstack~="") {
	    vtlistY = vtlistY + " " + poolstack + "_ps"
	}
	if ((m.model=="interactive") | (m.model=="interactiveiv")) {
	    vtlistY01 = ""
	    vtlistY_tk = tokens(vtlistY)
		numvt = cols(vtlistY_tk)
		for (j=1;j<=numvt;j++) {
			vtlistY01 = vtlistY01 + " " + vtlistY_tk[j] + "0"
			vtlistY01 = vtlistY01 + " " + vtlistY_tk[j] + "1"
		}
		vtlistY = strtrim(vtlistY01)
	}
	st_global("r(Y)",vtlistY)
	st_numscalar("r(numlrnY)",		numlrnY)
	// pystacked as only learner
	if (numlrnY==1) {
		est_main					= return_learner_item(e,e.vtlist,"est_main")
		est_options					= return_learner_item(e,e.vtlist,"est_options")
		if ((tokens(est_main)[1]=="pystacked") & (e.pystackedmulti>0)) {
			st_numscalar("r(numpslrnY)", e.pystackedmulti)
		}
		else if (tokens(est_main)[1]=="pystacked") {
			st_numscalar("r(numpslearnY)", 1)
			// if not pystacked, goes to general code
			// if pystacked with one learner, goes to general code
			allpystackedmulti			= 0
		}
		else {
			st_numscalar("r(numpslrnY)", 0)
			// if not pystacked, goes to general code
			// if pystacked with one learner, goes to general code
			allpystackedmulti			= 0
		}
		if ((e.pystackedmulti>0) & (m.model~="fiv")) {
		    vtlistY_L = ""
		    for (j=1;j<=e.pystackedmulti;j++) {
				if ((m.model=="interactive") | (m.model=="interactiveiv")) {
					vtlistY_L = vtlistY_L + " " + invtokens(e.vtlist) + "0_L" + strofreal(j)
					vtlistY_L = vtlistY_L + " " + invtokens(e.vtlist) + "1_L" + strofreal(j)
				}
				else {
					vtlistY_L = vtlistY_L + " " + invtokens(e.vtlist) + "_L" + strofreal(j)
				}
			}
			st_global("r(Y_L)", strtrim(vtlistY_L))
		}
	}
	else {
		st_numscalar("r(numpslrnY)", 0)
		// if not pystacked, goes to general code
		// if pystacked with one learner, goes to general code
		allpystackedmulti			= 0
	}

	// D eqns
	for (i=1;i<=cols(m.nameD);i++) {
		e = (m.eqnAA).get(m.nameD[1,i])
		shortstack						= e.shortstack
		poolstack						= e.poolstack
		numlrnD							= cols(e.vtlist)
		vtlistD							= invtokens(e.vtlist)
		if (shortstack~="") {
			vtlistD = vtlistD + " " + shortstack + "_ss"
		}
		if (poolstack~="") {
			vtlistD = vtlistD + " " + poolstack + "_ps"
		}
		if (m.model=="interactiveiv") {
			vtlistD01 = ""
			vtlistD_tk = tokens(vtlistD)
			numvt = cols(vtlistD_tk)
			for (j=1;j<=numvt;j++) {
				vtlistD01 = vtlistD01 + " " + vtlistD_tk[j] + "0"
				vtlistD01 = vtlistD01 + " " + vtlistD_tk[j] + "1"
			}
			vtlistD = strtrim(vtlistD01)
		}
		st_global("r(D"+strofreal(i)+")",vtlistD)
		if (m.model=="fiv") {
		    vtlistDh = ""
		    vtlistD_tk = tokens(vtlistD)
			numvt = cols(vtlistD_tk)
			for (j=1;j<=numvt;j++) {
				vtlistDh = vtlistDh + " " + vtlistD_tk[j] + "_h"
			}
			st_global("r(D"+strofreal(i)+"_h)",strtrim(vtlistDh))
		}
		st_numscalar("r(numlrnD"+strofreal(i)+")",	numlrnD)
		// pystacked as only learner
		if (numlrnD==1) {
			est_main					= return_learner_item(e,e.vtlist,"est_main")
			est_options					= return_learner_item(e,e.vtlist,"est_options")
			if ((tokens(est_main)[1]=="pystacked") & (e.pystackedmulti>0)) {
				st_numscalar("r(numpslrnD"+strofreal(i)+")", e.pystackedmulti)
			}
			else if (tokens(est_main)[1]=="pystacked") {
				st_numscalar("r(numpslearnD"+strofreal(i)+")", 1)
				// single pystacked learner, goes to general code
				allpystackedmulti		= 0
			}
			else {
				st_numscalar("r(numpslrnD"+strofreal(i)+")", 0)
				// non-pystacked learner, goes to general code
				allpystackedmulti		= 0
			}
			if ((e.pystackedmulti>0) & (m.model~="fiv")) {
				vtlistD_L = ""
				for (j=1;j<=e.pystackedmulti;j++) {
					if (m.model=="interactiveiv") {
						vtlistD_L = vtlistD_L + " " + invtokens(e.vtlist) + "0_L" + strofreal(j)
						vtlistD_L = vtlistD_L + " " + invtokens(e.vtlist) + "1_L" + strofreal(j)
					}
					else {
						vtlistD_L = vtlistD_L + " " + invtokens(e.vtlist) + "_L" + strofreal(j)
					}
				}
				st_global("r(D"+strofreal(i)+"_L)", strtrim(vtlistD_L))
			}
		}
		else {
			// multiple learners, goes to general code
			allpystackedmulti			= 0
		}
	}
	
	// Z eqns; if none, nothing returned
	for (i=1;i<=numeqnZ;i++) {
		e = (m.eqnAA).get(m.nameZ[1,i])
		shortstack						= e.shortstack
		poolstack						= e.poolstack
		numlrnZ										= cols(e.vtlist)
		vtlistZ							= invtokens(e.vtlist)
		if (shortstack~="") {
			vtlistZ = vtlistZ + " " + shortstack + "_ss"
		}
		if (poolstack~="") {
			vtlistZ = vtlistZ + " " + poolstack + "_ps"
		}
		st_global("r(Z"+strofreal(i)+")",vtlistZ)
		st_numscalar("r(numlrnZ"+strofreal(i)+")",	numlrnZ)
		// pystacked as only learner
		if (numlrnZ==1) {
			est_main					= return_learner_item(e,e.vtlist,"est_main")
			est_options					= return_learner_item(e,e.vtlist,"est_options")
			if ((tokens(est_main)[1]=="pystacked") & (e.pystackedmulti>0)) {
				st_numscalar("r(numpslrnZ"+strofreal(i)+")", e.pystackedmulti)
			}
			else if (tokens(est_main)[1]=="pystacked") {
				st_numscalar("r(numpslearnZ"+strofreal(i)+")", 1)
				// single pystacked learner, goes to general code
				allpystackedmulti		= 0
			}
			else {
				st_numscalar("r(numpslrnZ"+strofreal(i)+")", 0)
				// non-pystacked learner, goes to general code
				allpystackedmulti		= 0
			}
			if ((e.pystackedmulti>0) & (m.model~="fiv")) {
				vtlistZ_L = ""
				for (j=1;j<=e.pystackedmulti;j++) {
					vtlistZ_L = vtlistZ_L + " " + invtokens(e.vtlist) + "_L" + strofreal(j)
				}
				st_global("r(Z"+strofreal(i)+"_L)", strtrim(vtlistZ_L))
			}
		}
		else {
			// multiple learners, goes to general code
			allpystackedmulti			= 0
		}
	}
	
	st_numscalar("r(allpystackedmulti)",	allpystackedmulti)

}

////////////////////////////////////////////////////////////////////////////////

mata mlib create lddml, replace
mata mlib add lddml *()
mata mlib index
mata describe using lddml

end

* finally, return to original folder location
cd "`pwd'"
