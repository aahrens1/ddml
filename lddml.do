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

version 14
mata:
mata clear

void whichddml()
{
""
"Mata library for ddml and related programs,"
"compiled `current_date' under Stata " + "`stata_version'" + " born " + "`born_date'."
"authors AA/MS"
st_sclear()
st_global("s(stata_born_date)","`born_date'")
st_global("s(stata_version)","`stata_version'")
st_global("s(compiled_date)","`current_date")
}

// uniquely identified by vname = dep var in the equation
// equation structure: one for E[y|x], E[d|x] etc
struct eStruct {
	string scalar					vname			// name of variable to be orthogonalized
	real matrix						vtlist			// list of orthogonalized (learner) variables
	string scalar					shortstack		// name of shortstack variable
	string scalar					poolstack		// name of poolstack variable
	real scalar						nlearners		// number of learners
	real scalar						lieflag			// =1 if LIE spec with two estimation strings etc.
	real scalar						ateflag			// =1 if treatment variable in ATE/LATE
	real scalar						pystackedmulti	// =#learner if pystacked with multiple learners
	class AssociativeArray scalar	lrnAA			// AssociativeArray with all learners //
													// (keys=vtilde,object)
	class AssociativeArray scalar	resAA			// AssociativeArray with all learner results //
													// (keys=vtilde,object,rep)

}

// should perhaps make vname a required argument for this function
struct eStruct init_eStruct()
{
	struct eStruct scalar	m
//	class AssociativeArray scalar	A2, A3

	m.vname				= ""
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
	real colvector					id				// id variable (name in Stata will be modelname_id)
	real scalar						nreps			// number of resamplings
	real scalar						ncombos			// number of possible specifications (=0 if not yet estimated)
	real scalar						kfolds			// number of crossfitting folds
	class AssociativeArray scalar	eqnAA			// AssociativeArray with all equations
	class AssociativeArray scalar	estAA			// AA wth all estimation results
	string scalar					nameY			// dependent variable 
	string colvector				nameD			// treatment variable(s)
	string colvector				nameZ			// instrument(s)
	string scalar					fclustvar		// name of fold cluster variable (="" if none)
	real scalar						ssflag			// flag for short-stacking
	real scalar						psflag			// flag for pooled-stacking
	string scalar					strDatavars		// string with expanded names of Stata variables
	real matrix						matDatavars		// matrix with values of Stata variables
	real scalar						crossfitted   	// =number of reps for which crossfitting; 0 if not
	real scalar						estimated		// =1 if estimation has been done; 0 if not
	real scalar						ycounter		// counter for default y learners
	real scalar						dcounter		// counter for default d learners
	real scalar						zcounter		// counter for default z learners
}

struct mStruct init_mStruct()
{

	struct mStruct scalar	m
	
	m.model			= ""
	m.id			= J(0,1,.)
	m.nreps			= 0
	m.ncombos		= 0
	m.kfolds		= 0
	m.nameY			= ""
	m.nameD			= J(1,0,"")
	m.nameZ			= J(1,0,"")
	m.ssflag		= 0
	m.psflag		= 0
	m.fclustvar		= ""
	m.strDatavars	= ""
	m.matDatavars	= J(0,0,.)
	m.crossfitted	= 0
	m.estimated		= 0
	
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

// clear results from all eStructs in mStruct
void clear_model_results(struct mStruct m)
{
	// if not crossfitted, nothing to clear
	if (m.crossfitted>0) {
		struct eStruct scalar			e
		class AssociativeArray scalar	A3
		
		// clear results for each equation
		eqnlist = (m.nameY, m.nameD, m.nameZ)
		for (i=1; i<=cols(eqnlist); i++) {
			clear_equation_results((m.eqnAA).get(eqnlist[i]))
		}
		
		// also clear:
		m.ncombos		= 0
		m.crossfitted	= 0
		m.ssflag		= 0
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
	
	if (m.model=="interactive" | m.model=="late") {
	    comboY = comboY^2
	}
	if (m.model=="late" | m.model=="fiv") {
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


mata mlib create lddml, /* dir(PERSONAL) */ replace
mata mlib add lddml *()
mata mlib index
mata describe using lddml

end

* finally, return to original folder location
cd "`pwd'"
