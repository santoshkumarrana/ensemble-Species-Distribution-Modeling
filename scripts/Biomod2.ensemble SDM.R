###Ensemble species distribution modeling using Biomod2

# load the library
library(biomod2)

#1. Data preparation

#1a. Load species data and define species name, presence/absence data, and XY coordinate
DataSpecies <- read.csv(system.file("external/species/infile/helmax.csv", package="biomod2"), row.names = 1)
head(DataSpecies)

# define the name of studied species
myRespName <- 'helmax'

# the presence/absences data for our species
myResp <- as.numeric(DataSpecies[,myRespName])

# the XY coordinates of species data
myRespXY <- DataSpecies[,c("X_WGS84","Y_WGS84")]

#visualize species data through level plot
level.plot(data.in=myResp, XY=myRespXY, color.gradient = "red",
           cex = 0.5,
           level.range = c(min(myResp), max(myResp)),
           show.scale = TRUE,
           title = "level plot",
           SRC=FALSE,
           save.file="no",
           ImageSize="small",
           AddPresAbs=NULL,
           PresAbsSymbol=c(cex*0.2,10,2))

#multiple plot if infile has multispecies format
myRespName <- c("helmax", "anecyl", "helpau", "artfri")
myResp <- DataSpecies[,myRespName]
multiple.plot(Data = myResp,
              coor = myRespXY, color.gradient='red',
              plots.per.window=9,
              cex=1,
              save.file="no",
              name="multiple plot",
              ImageSize = "small",
              AddPresAbs = NULL,
              PresAbsSymbol = c(cex * 0.2, 10, 2))

##1b.Load the environmental raster layers (could be .img, ArcGIS rasters or any supported format by the raster package)

#load package
library(raster)

#Making stack of predictor variables from hard drive
#ensemble SDM for current (note: myExpl should be defined with the selected variables from VIF and used in final ensemble SDM)
myExpl = raster::stack( system.file("external/bioclim/current/bio1.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio2.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio3.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio4.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio5.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio6.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio7.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio8.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio9.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio10.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio11.grd",package="biomod2"),
						system.file("external/bioclim/current/bio12.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio13.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio14.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio15.grd",package="biomod2"),
                        system.file("external/bioclim/current/bio16.grd",package="biomod2"),
						system.file("external/bioclim/current/bio17.grd",package="biomod2"),
						system.file("external/bioclim/current/bio18.grd",package="biomod2"),
						system.file("external/bioclim/current/bio19.grd",package="biomod2"))


#1c. Formatting data
myBiomodData <- BIOMOD_FormatingData(resp.var = myResp,
                                     expl.var = myExpl,
                                     resp.xy = myRespXY,
                                     resp.name = myRespName,
                                     PA.nb.rep = 0,
                                     PA.nb.absences = 10000,
                                     PA.strategy = 'random')                                     
#1d.Check whether the data are correctly formatted by printing and plotting the created object##
myBiomodData

#2. Model building
#2a. Defining  mododelling options using default option (except, maximumiterations = 5000)
myBiomodOption <- BIOMOD_ModelingOptions(
  MAXENT.Phillips = list( path_to_maxent.jar = '/Users/sqr5923/esdm',
                          maximumiterations = 5000,
                          visible = FALSE,
                          linear = TRUE,
                          quadratic = TRUE,
                          product = TRUE,
                          threshold = TRUE,
                          hinge = TRUE,
                          lq2lqptthreshold = 80,
                          l2lqthreshold = 10,
                          hingethreshold = 15,
                          beta_threshold = -1,
                          beta_categorical = -1,
                          beta_lqp = -1,
                          beta_hinge = -1,
                          defaultprevalence = 0.5))  

Print_Default_ModelingOptions()

# if not use default option from Biomod2
myBiomodOption <- BIOMOD_ModelingOptions()


# 3. Running the models
myBiomodModelOut <- BIOMOD_Modeling(myBiomodData,
                                    models = c('GLM','GAM','GBM','CTA','SRE','RF','FDA','ANN','MARS','MAXENT.Phillips','MAXENT.Phillips.2'),
                                    models.options = myBiomodOption,
                                    NbRunEval=50,
                                    DataSplit=75,
                                    Yweights=NULL,
                                    VarImport=20,
                                    models.eval.meth = c('TSS','ROC'),
                                    SaveObj = TRUE,
                                    rescal.all.models = TRUE,
                                    do.full.models = TRUE,
                                    modeling.id = paste(myRespName,"FirstModeling",sep=""))


#3a. model summary
myBiomodModelOut

#3b. get all models evaluation
myBiomodModelEval <- get_evaluations(myBiomodModelOut)

#3c. print the dimnames of this object
dimnames(myBiomodModelEval)

#3d. let's print the TSS scores of all selected models
myBiomodModelEval["TSS","Testing.data",,,]

#3e. let's print the ROC scores of all selected models
myBiomodModelEval["ROC","Testing.data",,,]

# 4. Plot response curves
# 4a. Load the models for which we want to extract the predicted response curves
myGLMs <- BIOMOD_LoadModels(myBiomodModelOut, models='GLM')
myGAMs <- BIOMOD_LoadModels(myBiomodModelOut, models='GAM')
myGBMs <- BIOMOD_LoadModels(myBiomodModelOut, models='GBM')
myCTAs <- BIOMOD_LoadModels(myBiomodModelOut, models='CTA')
mySREs <- BIOMOD_LoadModels(myBiomodModelOut, models='SRE')
myRFs <- BIOMOD_LoadModels(myBiomodModelOut, models='RF')
myFDAs <- BIOMOD_LoadModels(myBiomodModelOut, models='FDA')
myANNs <- BIOMOD_LoadModels(myBiomodModelOut, models='ANN')
myMARSs <- BIOMOD_LoadModels(myBiomodModelOut, models='MARS')
myMaxents <- BIOMOD_LoadModels(myBiomodModelOut, models='MAXENT.Phillips')
myMaxNets <- BIOMOD_LoadModels(myBiomodModelOut, models='MAXENT.Phillips.2')

# 4b. plot 2D response plots (one by one for each sub models) (note: models  = myGLMs[51] indicates the plot for average model run if NbRunEval=50)
myRespPlot2D <- response.plot2(models  = myGLMs[51],
                               Data = get_formal_data(myBiomodModelOut,'expl.var'), 
                               show.variables= get_formal_data(myBiomodModelOut,'expl.var.names'),
                               do.bivariate = FALSE,
                               fixed.var.metric = 'median',
                               legend = TRUE, 
                               col = c("red","blue"),
                               data_species = get_formal_data(myBiomodModelOut,'resp.var'))

#5. Doing Ensemble Modelling (note: 1 - 4 should be performed twice; once we obtain selected variables from VIF analysis and selecting the suub-models based on the accuracy indices, need to start from 1b for selected variables and model building 3 with selected sub models)
myBiomodEM <- BIOMOD_EnsembleModeling( modeling.output = myBiomodModelOut,
                                       chosen.models = 'all',
                                       em.by = 'all',
                                       eval.metric = c('TSS'),
                                       eval.metric.quality.threshold = c(0.55),
                                       models.eval.meth = c('TSS'),
                                       prob.mean = TRUE,
                                       prob.cv = FALSE,
                                       prob.ci = FALSE,
                                       prob.ci.alpha = 0.05,
                                       prob.median = FALSE,
                                       committee.averaging = FALSE,
                                       prob.mean.weight = TRUE,
                                       prob.mean.weight.decay = 'proportional' )

#5a. print summary
myBiomodEM

#5b. get evaluation scores
get_evaluations(myBiomodEM)

#5c. print variable relative importances (note: change model = 'RF' for other sub models)
variables_importance(model='RF', data=myBiomodData, method="full_rand", nb_rand=21)

#6. projection over the globe under current conditions
myBiomodProj <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExpl,
  proj.name = 'current',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = FALSE,
  build.clamping.mask = FALSE,
  output.format = '.grd')

#6a. summary of created object
myBiomodProj

#6b. if you want to make custom plots, you can also get the projected map
myCurrentProj <- get_predictions(myBiomodProj)
myCurrentProj

#7. ensemble forecasting under current bioclimatic scenario
myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProj)

#7a.meta data
myBiomodEF

#7b. reduce layer names for plotting convegences
plot(myBiomodEF)



#8. ensemble SDM for future bioclimatic scenario

#8a. Load environmental variables for the future bioclimatic scenario (note: use selected variables corresponding to the current bioclimatic variables)
myExplFuture = raster::stack( system.file("external/bioclim/future/bio1.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio2.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio3.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio4.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio5.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio6.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio7.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio8.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio9.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio10.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio11.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio12.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio13.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio14.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio15.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio16.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio17.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio18.grd",package="biomod2"),
                              system.file("external/bioclim/future/bio19.grd",package="biomod2"))
##8b. Future projection
myBiomodProjFuture <- BIOMOD_Projection(
  modeling.output = myBiomodModelOut,
  new.env = myExplFuture,
  proj.name = 'Future',
  selected.models = 'all',
  binary.meth = 'TSS',
  compress = FALSE,
  clamping.mask = FALSE,
  output.format = '.grd')

#8c. files created on hard drive
list.files("glandulossisima/proj_Future/")

#8d. summary of crated oject
myBiomodProjFuture

#8d. if you want to make custom plots, you can also get the projected map
myfutureProj <- get_predictions(myBiomodProjFuture)
myfutureProj


#9. ensemble forecasting future
myBiomodEF <- BIOMOD_EnsembleForecasting(
  EM.output = myBiomodEM,
  projection.output = myBiomodProjFuture)

#9a. meta data
myBiomodEF

#9b. reduce layer names for plotting convegences
plot(myBiomodEF)