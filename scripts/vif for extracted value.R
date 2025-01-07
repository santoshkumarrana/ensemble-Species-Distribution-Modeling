# Variance Inflation factor (VIF) analysis for multicollinearity check
#1. Load libraries
library(car)

#2. Load datasets with extracted values for all bioclimatic variables (based on occurrence data only)
practice <- read.table(file = "clipboard", 
                   sep = "\t", header=TRUE)
attach(practice)
View(practice)

#3.fit regression model (elevation is the response variables, and all bioclimatic variables are predictor variables; make sure to extract the eleveation value for the occurrence points)
model <- lm(elevation~bio1+bio2+bio3+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio12+bio13+bio14+bio15+bio16+bio17+bio18+bio19,data=practice)
vif(model)
summary(model)
alias(model) #if there are alias in model due to sampling size

#Note: The model run should be single step removing a variables having high vif value.
#for instance, if bio12 has the highest vif value, the next model run should be with excluding bio12. These model run should be performed unless vif values are less than 10.
model <- lm(elevation~bio1+bio2+bio3+bio4+bio5+bio6+bio7+bio8+bio9+bio10+bio11+bio13+bio14+bio15+bio16+bio17+bio18+bio19,data=practice)
vif(model)