# MLogitViz testing
# Kelsey Gonzalez
# 2019-12-18
# last update: 2019-12-05



rm(list = ls()) # clean up the environment
require(foreign)
require(mlogit)

ml <- read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")

runmodel("prog", c("ses","write","math", "science", "female","schtyp","awards"), ml)
runmodel("prog", c("ses","write","female","schtyp","awards"), ml)
runmodel("prog", c("ses","write","math","female","schtyp","awards"), ml)

# Let's test it with some real data...
require(readxl)
Pew2013 <- read_xls("Pew2013Latino.xls", col_names = TRUE)
Pew2013 <-na.omit (Pew2013) #Listwise deletion to get the data running quickly
Pew2013$REGION <- as.factor(Pew2013$REGION)
Pew2013$RACE <- as.factor(Pew2013$RACE)
Pew2013$fivecat <- as.factor(Pew2013$fivecat)
Pew2013$PRIMARY_LANGUAGE <- as.factor(Pew2013$PRIMARY_LANGUAGE)
Pew2013$AGE <- as.factor(Pew2013$AGE)
Pew2013$PPARTY <- as.factor(Pew2013$PPARTY)
Pew2013$EDUCCAT2 <- as.factor(Pew2013$EDUCCAT2)
Pew2013$INCOMECAT <- as.factor(Pew2013$INCOMECAT)
Pew2013$ORIGIN <- as.factor(Pew2013$ORIGIN)
Pew2013$CITIZEN <- as.factor(Pew2013$CITIZEN)
Pew2013$RELI <- as.factor(Pew2013$RELI)
Pew2013$OP_IMM <- as.factor(Pew2013$OP_IMM)
Pew2013$OP_TYPICALAMERICAN <- as.factor(Pew2013$OP_TYPICALAMERICAN)
Pew2013$commondich <- as.factor(Pew2013$commondich)
Pew2013$generation <- as.factor(Pew2013$generation)
summary(Pew2013)
Pew2013 <- as.data.frame(Pew2013)

runmodel("fivecat", c("RACE","AGE","ORIGIN","generation", "PPARTY","EDUCCAT2","INCOMECAT"), Pew2013)


indig <- read.csv("Tribal Citizenship Database_import data analysis_FINAL.csv", stringsAsFactors = FALSE)
indig$bia.region <- as.factor(indig$bia.region)
indig$gaming <- as.factor(indig$gaming)
indig$criteria <- as.factor(indig$criteria)
indig$size <- as.numeric(indig$size)

runmodel("criteria", c("gaming", "bia.region", "size"), indig)
class(indig$size)
