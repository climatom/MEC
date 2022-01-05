# Import relevant module(s)
library(lubridate) # For easier date-time handling

# # # # # # # # # # # # # # 
# Constants
# # # # # # # # # # # # # # 
c<-2*pi # Radians in circle
r_mean<-149.6 # average earth-sun distance (million km)
Lf<-3.34*10^5 # latent heat of fusion (J/kg)
water_density<-1000. # kg/m^3
ice_density<-900 # kg/m^3
sheat_ice<-2090 # specific heat capacity of ice (J/kg/K)
target_lat<-67.903 # ~centroid of target glacier (deg N)
target_lon<-18.569 #  ~centroid of target glacier (deg E)
elev_e<-1386.299 # elevation of met data (m a.s.l)
secday=60*60*24 # seconds in day

# # # # # # # # # # # # # # 
# FILENAMES
# *** 
# These will need to be adjusted 
# to match the local file system
# *** 
# # # # # # # # # # # # # # # # # # # # # # # 
din<-"C:/Users/Admin/Desktop/Modelling/Data/"
hyps_name<-paste(din,"stor_hyps.csv",sep="")
met_name<-paste(din,"met.csv",sep="")
val_name<-paste(din,"AnnMB.csv",sep="")
# # # # # # # # # # # # # # # # # # # # # # # 

# # # # # # # # # # # # # # 
# PARAMETERS TO TUNE
# # # # # # # # # # # # # # 
tlapse<--6.5 # C/km
c<--0 # Bias correction for temp (C)
cp<-1.6 # Bias correction for precipitation (scalar [0->inf])
plapse<-100 # %/km
snow_trans<-1.5 # Precip falls as snow if > this (C)
snow_albedo<-0.8 # Dimensionless 
firn_albedo<-0.55 # Dimensionless
ice_albedo<-0.4 # Dimensionless
alb_time_scale<-21.9 # Days
alb_depth_scale<-0.001 # m w.e.
layer_depth <- 2.0 # m
t_tip<-1 # C
t_sens<- 12 # W/m^2/C
t_constant <--25 #W/m^2
trans<-0.5 # Transmissivity for insolation, dimensionless

# # # # # # # # # # # # # # # # # # # # # # # #
# FUNCTIONS
# # # # # # # # # # # # # # # # # # # # # # # #
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
# # # # # # # # # # # # # # # # # # # # # # # #

# T correction / distribution
tdist <- function(t,c,lapse,z){
  tout <- (t+c)+(z-elev_e)*(lapse/1000.0)
  tout
}

# Snowfall distribution
sdist <- function(tiz,p,z=z,snow=rep(0,times=ne)){
  idx<-tiz<snow_trans
  snow[idx]=(p*cp)*exp((z[idx]-elev_e)/1000.0*
                         plapse/(100))/1000.0 # mm -> m l.e.
  snow
}

days_since<-function(snow,dtsnowiz){
  dtsnowiz[snow>0]<-0
  dtsnowiz[snow==0]<-dtsnowiz[snow==0]+1
  dtsnowiz
}

# Insolation
decl <- function(lat,doy){
  dec=deg2rad(23.44)*cos(c*(doy-172.)/365.25)
  dec
}
sin_elev <- function(dec,hour,lat,lon=0){
  lat<-deg2rad(lat)
  out<-sin(lat)*sin(dec) - cos(lat)*cos(dec) * 
    cos(c*hour/24.+deg2rad(lon))
  out
}

sun_dist <- function(doy){
  m<-c*(doy-4.)/365.25
  v<-m+0.0333988*sin(m)+0.0003486*sin(2.*m)+0.0000050*sin(3.*m)
  r<-149.457*((1.-0.0167**2)/(1+0.0167*cos(v)))
  r
}

sin_toa <- function(doy,hour=12,lat,lon=0){
  dec<-decl(lat,doy)
  sin_elev<-sin_elev(dec,hour,lat,lon)
  r<-sun_dist(doy)
  s<-1366.*((r_mean/r)^2)
  toa<-sin_elev*s
  toa[toa<0]<-0.
  toa
}

# TEMP-DEP HEAT FLUX
temp_def <-function(temps,tip,scalar,constant){
  qt <- rep(NA, times=length(temps))
  qt[temps>=tip]<-temps[temps>=tip]*scalar+constant
  qt[temps<tip]<-constant
  qt
}

# ALBEDO
albedo <- function(dtlast,snow_depth){
  
  alb_snow<-firn_albedo+(snow_albedo-firn_albedo)*
    exp(-dtlast/alb_time_scale)
  
  alb<-alb_snow+(ice_albedo-alb_snow)*
    exp(-snow_depth/alb_depth_scale)
  
  alb
}
  
# SEB
qsum<-function(tdep,sin,trans,albedo){
  q<-(tdep+(1-albedo)*sin*trans)*secday # NB Joules/day
  q[q<0]=0
  q
}

# SUB-SURFACE TEMP
tsub<-function(tsublast,fraction,q){
  dt<-fraction*q/
    (ice_density*sheat_ice*layer_depth)
  tsub<-tsublast+dt 
  tsub
}

# REFREEZE FRACTION
refreeze<-function(snowdepth,tsub){
  r<-1-exp(tsub)
  r[snowdepth==0]=0
  r
}

# RUNOFF
runoff<-function(q,r){
  melt<-(1-r)*q/(Lf*1000.0) # m l.e.
}

# # # # # # # # # # # # # #
# PARSE INPUT
# # # # # # # # # # # # # #

# Storglaciaren hypsometry: z_m | area_m2
hyps<-read.csv(hyps_name)
z<-hyps$ï..z_m
z<-z[hyps$area_m2>0]
ne<-length(z)
totArea<-sum(hyps$area_m2)

# E-OBS met data for Storglaciaren: Date | t (C) | p (mm/day)
met<-read.csv(met_name)
met$date<-as.Date(met$date,format="%d/%m/%Y")
nt<-length(met$date)

# Validation
val<-read.csv(val_name)

# Prepare intermediate and output grids
tiz<-matrix(nrow=nt,ncol=ne) # Air T (C) at time i and elevation z
tdep<-matrix(nrow=nt,ncol=ne) # Temperature-dependent heat flux at time i, elevation z
snowiz<-matrix(0,nrow=nt,ncol=ne) # Snowfall at time i, elevation z
sdiz<-matrix(0,nrow=nt,ncol=ne) # Snow depth at time i and elevation z
albiz<-matrix(nrow=nt,ncol=ne) # Albedo at time i and elevation z
albedoiz<-matrix(nrow=nt,ncol=ne) # Albedo at time t, elevation z
qiz<-matrix(nrow=nt,ncol=ne) # Albedo at time t, elevation z
riz<-matrix(nrow=nt,ncol=ne) # Frac meltwater that refrezzes at time t, elevation z
runoffiz<-matrix(nrow=nt,ncol=ne) # Runoff at time t, elevation z
mbiz<-matrix(nrow=nt,ncol=ne)# Mass balance at time i, elevation z
SMB<-matrix(nrow=nt) # Specific MB at time i (m w.e.)
annMB<-rep(NA,70) # Annual specific mass balance
yrs<-rep(NA,70) # Years (useful for subsetting)

# # # # # # # # # # # # # #
# MAIN BEGINS
# # # # # # # # # # # # # #

# Compute sin as the toa for the
# doy in the met df
doy<-yday(met$date)
toa<-sin_toa(doy=doy,lat=target_lat)

# Loop time 
for (i in 1:nt){
  
  # Distribute T
  tiz[i,]<-tdist(met$t[i],c=c,lapse=tlapse,z=z)
  
  # Init on first run
  if (i==1){
    ltm<-mean(met$t)
    tsubiz_last<-tdist(ltm,c=c,lapse=tlapse,z)
    tsubiz_last[tsubiz_last>0]<-0
    dtsnowiz_last<-rep(0,ne)
  
    # Distribute winter balance over the 
    # elevation bands (using regression based on 
    # sample of mean Bw profiles)
    sdiz_last<-z*0.00442-4.797
  }

  # T-dep heat flux
  tdep[i,]<-temp_def(tiz[i,],tip=t_tip,scalar=t_sens,constant=t_constant)  
  
  # Albedo
  albedoiz[i,]<-albedo(dtlast=dtsnowiz_last,snow_depth=sdiz_last)
  
  # SEB 
  qiz[i,]<-qsum(tdep[i,],toa[i],trans,albedoiz[i,])

  # Refreeze fraction
  riz[i,] <- refreeze(sdiz_last,tsubiz_last)
  
  # Runoff
  runoffiz[i,]<-runoff(qiz[i,],riz[i,])
  
  # Distribute new snowfall
  snowiz[i,]<-sdist(tiz[i,],p=met$p[i],z=z)
  
  # Evaluate SMB at this time step
  mbiz[i,]<-snowiz[i,]-runoffiz[i,]
  
  # Evaluate new snow depth
  sdiz[i,]<-sdiz_last+mbiz[i,]
  sdiz[i,sdiz[i,]<0]<-0
  
  # Cache snow depth for next iteration
  sdiz_last<-sdiz[i,]*1.
  
  # New sub-surface temp
  tsubiz_last<-tsub(tsubiz_last,riz[i,],qiz[i,])

  # Update days since last snow
  dtsnowiz_last<-days_since(snowiz[i,],dtsnowiz_last)
  
  # If today is end of summer, reset tsub to annual mean tiz
  if (doy[i]==306){
    tsubiz_last<-tdist(ltm,c=1,lapse=tlapse,z)
    tsubiz_last[tsubiz_last<0]<-0
  }
  SMB[i]<-sum(hyps$area_m2[hyps$area_m2>0]*mbiz[i,])/totArea
}

# Compute annual MB (Sept 01 each year)
i<-1
yri<-1950
for(y in seq(19500901,20190901,10000)){
  idx<-met$date>ymd(y)&met$date<ymd((y+10000))
  annMB[i]<-sum(SMB[idx])
  yrs[i]<-yri+i
  i<-i+1
}

# 1951-2018 correlation with obs?
obs<-val$MBs[val$ï..year>1950]
sim<-annMB[yrs<2019]
plot(sim,obs)
r<-cor(sim,obs) 

# % error in cumulative mass balance - obs and sim
cum_err_pc<-abs(sum(sim)-sum(obs))/abs(sum(obs))*100.0
print(sprintf("Obs CMB is %.2f m.w.e.",sum(obs)))
print(sprintf("Sim CMB is %.2f m.w.e.",sum(sim)))
print(sprintf("[Error is %.1f%%]",cum_err_pc))
print(sprintf("Annual correlation is %.2f",r))        

