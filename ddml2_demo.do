
cap prog drop ddml2
ddml2 init partial, mname(myest)
ddml2 yeq, gen(lassoy) mname(myest) vname(logpgp95) : lasso2 logpgp95 edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml2 yeq, gen(rigy) mname(myest) vname(logpgp95) : rlasso logpgp95 edes1975 avelf temp* humid* steplow-oilres
ddml2 deq, gen(lassod) mname(myest) vname(avexpr) : lasso2 avexpr edes1975 avelf temp* humid* steplow-oilres, lic(aicc) postres
ddml2 deq, gen(rigd) mname(myest) vname(lat_abst) : rlasso lat_abst edes1975 avelf temp* humid* steplow-oilres
cap prog drop ddml2
ddml2 desc, mname(myest)
cap prog drop ddml2
ddml2 crossfit, mname(myest)
cap prog drop ddml2
ddml2 desc, mname(myest)
cap prog drop ddml2
ddml2 estimate, mname(myest)
