
program define _ddml_display_header, rclass

	syntax [anything] , [	///
							str(string) ///
							]

		di
		di as res "Mean-squared error for `str':"
		di _col(2) "Name" _c
		di _col(20) "Orthogonalized" _c
		di _col(40) "Command" _c
		di _col(54) "N" _c
		di _col(65) "MSPE"
		di "{hline 75}"

end
