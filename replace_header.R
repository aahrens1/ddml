
setwd("/Users/kahrens/MyProjects/ddml")

library("stringr")

all.files <- dir()
all.files <- all.files[str_detect(all.files,"ado")]
 
for (i in 1:length(all.files)) {

  txt <- readLines(all.files[i])
  txt[1] <- "*! ddml v1.3"
  txt[2] <- "*! last edited: 11july2023"
  writeLines(txt,con=all.files[i])

}

all.files <- dir()
all.files <- all.files[str_detect(all.files,"sthlp")]

for (i in 1:length(all.files)) {
  
  txt <- readLines(all.files[i])
  txt[2] <- "{* *! version 11july2023}{...}"
  txt[4] <- str_replace(txt[4],"\\{right: v1.2\\}","\\{right: v1.3\\}")
  writeLines(txt,con=all.files[i])
  
}