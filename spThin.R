# Load libraries
library(spThin)
library(tidyverse)

#1. set the working directory and keep the infile
setwd("D:/PennState_USA/esdm/pseudoabs/ps")
#2. load excel data
data <-read.csv ("verstr.csv", row.names=1)
attach (data)
#3. thinning of datasets
thinned_datasets_full<-thin (loc.data=data, lat.col="LAT", 
                             long.col="LONG", spec.col="SPEC",
                             thin.par=10, reps=10000, locs.thinned.list.return=TRUE, 
                             write.files=TRUE, max.files=1, 
                             out.dir="D:/PennState_USA/esdm/pseudoabs/ps", out.base="verstr10_ps", 
                             write.log.file=TRUE, log.file="verstr20.txt")

#4.remove duplicates
setwd("D:/PennState_USA/esdm/biomod2/Functional_group/aster")

practice <- read.table(file = "clipboard", 
                       sep = "\t", header=TRUE)
attach (practice)

# Remove duplicated rows based on X_WGS84
data1 <- practice %>% distinct(X_WGS84, .keep_all = TRUE)

##save file
write.csv(data1, file = "aster_remove.csv", row.names = TRUE)
