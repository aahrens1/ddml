*! loghelp v1.0
*! last edited: 24 july 2023
*! author: ms

program define loghelp, rclass
    version 16.0
	syntax using, inputname(string) [smcl text name append replace msg vsquish]
	
	// default is no log message
	if "`msg'"=="" local nomsg nomsg
	// no newlines between Stata outputs
	local vsquishflag = ("`vsquish'"~="")

	log `using', name(`name') `smcl' `text' `append' `replace' `nomsg'
	// blank line follwing opening message if requested
	if "`msg'"~="" di
	mata: m_loghelp("`inputname'",`vsquishflag')
	log close
end

mata:

void m_loghelp(string scalar inputname,
				real scalar vsquishflag) {

	fh_in  = fopen(inputname, "r")
	while ((line=fget(fh_in))!=J(0,0,"")) {
		if (strmatch(line,"*{stata*")) {
			// line is a stata command
			// remove start of line
			spos = strpos(line,"{stata")
			line = substr(line,spos+6,.)
			// look for " (ascii 34) that starts the string
			spos = strpos(line,char(34))
			// and remove it and any preceding chars (spaces)
			line = substr(line,spos+1,.)
			// look for " (ascii 34) that ends the string
			spos = strrpos(line,char(34))
			// and remove it any any following chars (smcl etc.)
			slen = strlen(line)
			line = substr(line,1,spos-1)
			// echo line and execute it
			printf("{input}. %s\n",line)
			// print blank line unless vsquish specified
			stata(line)
			if (vsquishflag==0) {
				printf("\n")
			}
		}
		else {
			// line is not a stata command
			printf(line)
			printf("\n")
		}
	}
	fclose(fh_in)
}

end

