
program _ddml_save
	version 13

	syntax , mname(name) fname(string) [ replace ]

	if "`replace'"~="" {
		// Mata function to delete file
		mata: unlink("`fname'")
	}
	mata: save_model("`fname'",`mname')

end

mata:

void save_model(					string scalar fname,
									struct ddmlStruct m)
{
	fh = fopen(fname,"w")
	fputmatrix(fh,m)
	fclose(fh)
}

end