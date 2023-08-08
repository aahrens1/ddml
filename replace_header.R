setwd("/Users/kahrens/MyProjects/ddml")

library("stringr")

all.files <- dir()
all.files <- all.files[str_detect(all.files,"ado")]

for (i in 1:length(all.files)) {
  
  txt <- readLines(all.files[i])
  txt[1] <- "*! ddml v1.4.2"
  txt[2] <- "*! last edited: 8aug2023"
  writeLines(txt,con=all.files[i])
  
}

all.files <- dir()
all.files <- all.files[str_detect(all.files,"sthlp")]

for (i in 1:length(all.files)) {
  
  txt <- readLines(all.files[i])
  for (j in 1:20) {
    if(!is.na(txt[j])) txt[j] <- str_replace(txt[j],"version 28july2023\\}","version 8aug2023\\}")
    if(!is.na(txt[j])) txt[j] <- str_replace(txt[j],"\\{right: v1.4.1\\}","\\{right: v1.4.2\\}")
  }
  writeLines(txt,con=all.files[i])
  
}