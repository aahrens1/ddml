
cap prog drop ddml2
cap prog drop allcombos
cap prog drop _ddml_estimate_partial
cap prog drop _ddml_crossfit_partial

ddml2 init partial, mname(myest)
ddml2 yeq, gen(lassoy) mname(myest) vname(logpgp95) : lasso2 logpgp95 edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml2 yeq, gen(rigy) mname(myest) vname(logpgp95) : rlasso logpgp95 edes1975 avelf temp* humid* steplow-oilres
ddml2 deq, gen(lassod) mname(myest) vname(avexpr) : lasso2 avexpr edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml2 deq, gen(rigd) mname(myest) vname(lat_abst) : rlasso lat_abst edes1975 avelf temp* humid* steplow-oilres
ddml2 desc, mname(myest)
ddml2 crossfit, mname(myest)
ddml2 desc, mname(myest)
ddml2 estimate, mname(myest)

*************************************************************************

// multiple D variables example
// includes nocrossfit

ddml2 init partial, mname(myest)

ddml2 yeq, gen(lassoy) mname(myest) vname(logpgp95) : lasso2 logpgp95 edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml2 yeq, gen(rigy) mname(myest) vname(logpgp95) nocrossfit: rlasso logpgp95 edes1975 avelf temp* humid* steplow-oilres

ddml2 deq, gen(lassod1) mname(myest) vname(avexpr) : lasso2 avexpr edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml2 deq, gen(rigd1) mname(myest) vname(avexpr) : rlasso avexpr edes1975 avelf temp* humid* steplow-oilres
ddml2 deq, gen(lassod2) mname(myest) vname(lat_abst) : lasso2 lat_abst edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml2 deq, gen(rigd2) mname(myest) vname(lat_abst) nocrossfit : rlasso lat_abst edes1975 avelf temp* humid* steplow-oilres

ddml2 desc, mname(myest)

ddml2 crossfit, mname(myest)
ddml2 desc, mname(myest)

ddml2 estimate, mname(myest)
