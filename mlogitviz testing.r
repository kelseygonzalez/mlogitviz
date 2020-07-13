# MLogitViz testing
# Kelsey Gonzalez
# 2019-12-18
# last update: 2020-07-12



ml <- foreign::read.dta("https://stats.idre.ucla.edu/stat/data/hsbdemo.dta")

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
model = multinom(fivecat ~ RACE + AGE + ORIGIN + generation +  PPARTY + EDUCCAT2 + INCOMECAT, data = Pew2013)



require(tidyverse)
happy <- read_csv("happy.csv")   %>% 
  mutate(year = year,
         happiness = as.factor(happiness),
         realrinc = as.numeric(ifelse(realrinc == "Not applicable", 0, realrinc)),
         sex = as.factor(sex),
         race = as.factor(race),
         racedif = as.factor(ifelse(racedif == 1, "Structural Racisms", "No Structural Racism")),
         edad = as.numeric(ifelse(edad == "89 or older", 89, edad)),
         big_region = as.factor(big_region),
         urban = as.factor(ifelse(urban == 1, "Urban", "Non-Urban")),
         foreignborn = as.factor(ifelse(foreignborn == 1, "Foreign Born", "US Born")),
         pparty2 = as.factor(pparty2),
         married = as.factor(ifelse(married == 1, "Married", "Non-Married")),
         children = as.factor(ifelse(children == 1, "Has Children", "No Children")),
         college = as.factor(ifelse(college == 1, "College Grad", "Non-College Grad")),
         employed = as.factor(ifelse(employed == 1, "Employed", "Non-Employed")),
         satisfied_fin = as.factor(satisfied_fin)) 
         
# this still won't work. There's still a bug somewhere. 
runmodel("happiness", c("realrinc","sex","race", "racedif","edad","big_region"), happy)