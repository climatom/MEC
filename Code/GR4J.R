library(airGR)
library(lubridate)

# # # # # # # # # # # # # # 
# File names / Parse input
# # # # # # # # # # # # # # 
Lat<-51.75
ca<-9948
#if (!grepl("Data",getwd())){
 # setwd("Data") 
#}
setwd("/Users/tommatthews/MEC/Data/") # Change this to the absolute 
# path *IF* running on your own machine

obs_name<-paste("Thames.csv",sep="")
obs<-read.csv(obs_name)
jd=yday(obs$date)

# Compute potential evapotranspiration
obs$pe<-PE_Oudin(jd, obs$temp, Lat, LatUnit = 'deg', 
             TimeStepIn = "daily", TimeStepOut = "daily", RunFortran = FALSE)

# Transform cumec to mm
obs$qmmd<-obs$q/(ca*1e6)*1e3*60^2*24

# Init model
InputsModel <- CreateInputsModel(FUN_MOD = RunModel_GR4J, 
                                DatesR = as.POSIXlt(obs$date),
                                Precip = obs$pr, PotEvap = obs$pe)
# Select period to run
Ind_Run <- seq(which(format(obs$date, format = "%Y-%m-%d") == "1961-01-01"), 
               which(format(obs$date, format = "%Y-%m-%d") == "1990-12-31"))

RunOptions <- CreateRunOptions(FUN_MOD = RunModel_GR4J,
                               InputsModel = InputsModel, IndPeriod_Run = Ind_Run,
                               IniStates = NULL, IniResLevels = NULL, 
                               IndPeriod_WarmUp = NULL)

InputsCrit <- CreateInputsCrit(FUN_CRIT = ErrorCrit_NSE, InputsModel = InputsModel, 
                               RunOptions = RunOptions, Obs = obs$qmmd[Ind_Run])

CalibOptions <- CreateCalibOptions(FUN_MOD = RunModel_GR4J, 
                                   FUN_CALIB = Calibration_Michel)

OutputsCalib <- Calibration_Michel(InputsModel = InputsModel, RunOptions = RunOptions,
                                   InputsCrit = InputsCrit, CalibOptions = CalibOptions,
                                   FUN_MOD = RunModel_GR4J)

Param<-c(699.244,	0.212,	44.701,	2.257)

OutputsModel <- RunModel_GR4J(InputsModel = InputsModel, RunOptions = RunOptions, 
                              Param = Param)

plot(OutputsModel, Qobs = obs$qmmd[Ind_Run])
