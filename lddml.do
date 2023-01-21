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
	string scalar					vname		// name of variable to be orthogonalized
	real matrix						vtlist		// list of orthogonalized (learner) variables
	string scalar					shortstack	// name of shortstack variable
	real scalar						nlearners	// number of learners
	real scalar						lieflag		// =1 if LIE spec with two estimation strings etc.
	real scalar						ateflag		// =1 if treatment variable in ATE/LATE
	class AssociativeArray scalar	lrnAA		// AssociativeArray with all learners //
												// (keys=vtilde,object)
	class AssociativeArray scalar	resAA		// AssociativeArray with all learner results //
												// (keys=vtilde,object,rep)

}

// should perhaps make vname a required argument for this function
struct eStruct init_eStruct()
{
	struct eStruct scalar			d
//	class AssociativeArray scalar	A2, A3

	d.vname			= ""
	d.vtlist		= J(1,0,"")
	d.shortstack	= ""
	d.nlearners		= 0
	d.lieflag		= 0
	d.ateflag		= 0
	
	(d.lrnAA).reinit("string",2)
	(d.resAA).reinit("string",3)
	(d.lrnAA).notfound(NULL)
	(d.resAA).notfound(NULL)
	
	return(d)
}

// clear results from eStruct
void clear_equation_results(struct eStruct e)
{
	(e.resAA).reinit("string",3)
	e.shortstack	= ""
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
	real scalar						ssflag			// flag for shortstacking
	string scalar					strDatavars		// string with expanded names of Stata variables
	real matrix						matDatavars		// matrix with values of Stata variables
	real scalar						crossfitted   	// =1 if crossvalidation has been done; 0 if not
	real scalar						estimated		// =1 if estimation has been done; 0 if not
	real scalar						ycounter		// counter for default y learners
	real scalar						dcounter		// counter for default d learners
	real scalar						zcounter		// counter for default z learners
}

struct mStruct init_mStruct()
{

	struct mStruct scalar			d
	
	d.model			= ""
	d.id			= J(0,1,.)
	d.nreps			= 0
	d.ncombos		= 0
	d.kfolds		= 0
	d.nameY			= ""
	d.nameD			= J(1,0,"")
	d.nameZ			= J(1,0,"")
	d.ssflag		= 0
	d.fclustvar		= ""
	d.strDatavars	= ""
	d.matDatavars	= J(0,0,.)
	d.crossfitted	= 0
	d.estimated		= 0
	
	(d.eqnAA).reinit("string",1)
	(d.eqnAA).notfound(NULL)
	
	(d.estAA).reinit("string",2)
	(d.estAA).notfound(NULL)
	
	// initialize counters used for default learner names
	d.ycounter		= 1
	d.dcounter		= 1
	d.zcounter		= 1
	
	return(d)
}

// clear results from all eStructs in mStruct
void clear_model_results(struct mStruct d)
{
	// if not crossfitted, nothing to clear
	if (d.crossfitted==1) {
		struct eStruct scalar			e
		class AssociativeArray scalar	A3
		
		eqnlist = (d.nameY, d.nameD, d.nameZ)
		for (i=1; i<=cols(eqnlist); i++) {
			clear_equation_results((d.eqnAA).get(eqnlist[i]))
		}
		
		// also clear:
		d.ncombos		= 0
		d.crossfitted	= 0
		d.ssflag		= 0
		d.strDatavars	= ""
		d.matDatavars	= J(0,0,.)
		(d.estAA).reinit("string",2)
		(d.estAA).notfound(NULL)
	}
}

// clear only estimation results in mStruct
void clear_model_estimation(struct mStruct d)
{
    d.estimated		= 0
	d.ncombos		= 0
	(d.estAA).reinit("string",2)
	(d.estAA).notfound(NULL)
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
	return((e.lrnAA).get((key1,key2)))
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
	return((e.resAA).get((key1,key2,rep)))
}


mata mlib create lddml, /* dir(PERSONAL) */ replace
mata mlib add lddml *()
mata mlib index
mata describe using lddml

end

* finally, return to original folder location
cd "`pwd'"
