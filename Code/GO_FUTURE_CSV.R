library(lubridate)
library(dplyr)

# # # # # # # # # # # # # # # # # # # # 
# Constants (or 'hidden' paramaters)
# DO NOT ADJUST 
# # # # # # # # # # # # # # # # # # # # 
c<-2*pi # Radians in circle
r_mean<-149.6 # average earth-sun distance (million km)
Lf<-3.34*10^5 # latent heat of fusion (J/kg)
water_density<-1000. # kg/m^3
ice_density<-900 # kg/m^3
sheat_ice<-2090 # specific heat capacity of ice (J/kg/K)
target_lat<-67.903 # ~centroid of target glacier (deg N)
target_lon<-18.569 #  ~centroid of target glacier (deg E)
elev_e<-1386.299 # elevation of met data (m a.s.l)
secday<-60*60*24 # seconds in day
va<-1.36 # Coefficient in volume-area scaling; treat as a constant here
max_inc <-1.5 # When glacier is 'growing', insist that area cannot increase by 
# more than this much between adjacent bands
temp_corr<-0 # Bias correction for temp (C)
cp<-2 # Bias correction for precipitation (scalar [0->inf])
snow_trans<-1.5 # Precip falls as snow if > this (C)
snow_albedo<-0.9 # Dimensionless 
firn_albedo<-0.55 # Dimensionless
ice_albedo<-0.35 # Dimensionless
alb_time_scale<-21.9 # Days
alb_depth_scale<-0.001 # m w.e.
layer_depth <- 2.0 # m
future_start<-"2014-12-31"
hist_start<-"1981-01-01"
hist_stop<-"2010-12-31"
proj_start<-"1989-03-31"
spin_up_end<-"1990-08-30"
proj_end<-"2100-12-31"
mean_start<-as.Date(spin_up_end) %m-%months(12) + 2 # Start date for mean comps

# # # # # # # # # # # # # # 
# File names / Parse input
# # # # # # # # # # # # # # 
if (!grepl("Data",getwd())){
  setwd("Data") 
}
#setwd("C:/Users/Admin/Documents/MEC/Data/") # Change this to the absolute 
# path --  e.g., /Users/tommatthews/Documents/MEC/Data *IF* running on your own 
# machine. 
hyps_name<-paste("stor_hyps.csv",sep="")
met_name<-paste("projections.csv",sep="")

# Storglaciaren hypsometry: z_m | area_m2
hyps<-read.csv(hyps_name)
hyps_old<-hyps
z<-hyps$z_m
met_all<-read.csv(met_name)
met_all$date<-as.Date(met_all$date,format="%Y-%m-%d")
#z<-z[hyps$area_m2>0]
ne<-length(z)
totArea<-sum(hyps$area_m2)

# # # # # # # # # # # # # # 
# PARAMETERS TO ** SET **
# DO CHANGE THESE!! 
# # # # # # # # # # # # # # 
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#
tlapse<--6.5 # C/km
plapse<-100 # %/km
t_sens<-10 # W/m^2/C
t_constant <--25 #W/m^2
trans<-0.2 # Transmissivity for insolation, dimensionless
t_tip<-1 # Temp-dep fluxes increase with T above this threshold (C)
#------------------------------------------------------------------------------#
#------------------------------------------------------------------------------#

# # # # # # # # # # # # # # # # # # # # # # # #
# FUNCTIONS
# # # # # # # # # # # # # # # # # # # # # # # #
rad2deg <- function(rad) {(rad * 180) / (pi)}
deg2rad <- function(deg) {(deg * pi) / (180)}
# # # # # # # # # # # # # # # # # # # # # # # #

# T correction / distribution
tdist <- function(t,temp_corr,lapse,z){
  tout <- (t+temp_corr)+(z-elev_e)*(lapse/1000.0)
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
  melt
}

# Prepare output grids
yrs<-year(met_all$date)
nt<-length(met_all$date)
ny<-year(proj_end)-year(spin_up_end)+1

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
annMB<-rep(NA,ny) # Annual specific mass balance
annArea<-rep(NA,ny) # Area at the end of each year (m**2)
annT<-rep(NA,ny) # Annual mean air temperature (m**2)
annP<-rep(NA,ny) # Annual mean precipitation
yrs<-rep(NA,ny) # Years (useful for subsetting)
k<-1 # Annual counter
ref_idx<-seq(length(hyps[,1]))

# # # # # # # # # # # # # #
# MAIN BEGINS
# # # # # # # # # # # # # #

# Compute sin as the toa for the doy in the met df
date_hourly<-seq(as.POSIXct(met_all$date[1]), as.POSIXct(met_all$date[nt]), 
                 by="hour")
doy_hour<-yday(date_hourly)
hour_hour<-hour(date_hourly)
toa<-sin_toa(doy=doy_hour,hour=hour_hour,lat=target_lat)
# Now average to daily
toa<-data.frame(toa)
toa$date<-date_hourly
toa$daydate<-as.Date(toa$date, format = "%Y-%m-%d")
toa <- toa %>% 
  group_by(daydate) %>%
  dplyr::summarize(mean = mean(toa))

# Loop time 
nd<-0
totT<-0
totP<-0

for (i in 1:nt){
  
  # Distribute T
  tiz[i,]<-tdist(met_all$t[i],temp_corr=temp_corr,lapse=tlapse,z=z)
  
  # Init on first run
  if (i==1){
    ltm<-mean(met_all$t[met_all$date<hist_stop]) # Note we use LTM from hist 
    tsubiz_last<-tdist(ltm,temp_corr=temp_corr,lapse=tlapse,z)
    tsubiz_last[tsubiz_last>0]<-0
    dtsnowiz_last<-rep(0,ne)
    
    # Distribute winter balance over the 
    # elevation bands (using regression based on 
    # sample of mean Bw profiles)
    sdiz_last<-z*0.00442-4.797
    sdiz_last[sdiz_last<0]<-0
    
    # Set cmb_ann to be the same as the mean winter balance
    cmb_ann<-sum(hyps$area_m2[hyps$area_m2>0]*sdiz_last[hyps$area_m2 > 0])/totArea
  }
  
  # T-dep heat flux
  tdep[i,]<-temp_def(tiz[i,],tip=t_tip,scalar=t_sens,constant=t_constant) 
  
  # Albedo
  albedoiz[i,]<-albedo(dtlast=dtsnowiz_last,snow_depth=sdiz_last)
  
  # SEB 
  qiz[i,]<-qsum(tdep[i,],toa$mean[i],trans,albedoiz[i,])
  
  # Refreeze fraction
  riz[i,] <- refreeze(sdiz_last,tsubiz_last)
  
  # Runoff
  runoffiz[i,]<-runoff(qiz[i,],riz[i,])
  
  # Distribute new snowfall
  snowiz[i,]<-sdist(tiz[i,],p=met_all$p[i],z=z)
  
  # Evaluate SMB at this time step
  mbiz[i,]<-snowiz[i,]-runoffiz[i,]
  
  # Compute the annual SMB
  SMB[i]<-sum(hyps$area_m2[hyps$area_m2>0]*mbiz[i,hyps$area_m2 > 0])/totArea
  
  # Accumulate SMB
  cmb_ann<-cmb_ann + SMB[i]
  
  # Evaluate new snow depth
  sdiz[i,]<-sdiz_last+mbiz[i,]
  sdiz[i,sdiz[i,]<0]<-0
  
  # Cache snow depth for next iteration
  sdiz_last<-sdiz[i,]*1.
  
  # New sub-surface temp
  tsubiz_last<-tsub(tsubiz_last,riz[i,],qiz[i,])
  
  # Update days since last snow
  dtsnowiz_last<-days_since(snowiz[i,],dtsnowiz_last)
  
  # Accumulate T and P, plus increment the counters
  if (met_all$date[i]>=mean_start){
    totT<-totT+met_all$t[i]
    totP<-totP+met_all$p[i]
    nd<-nd+1 # Increment number of days here. 
    # [1] "2024-08-01"
  }
  
  # If today is end of summer (and it's not the spin-up period)
  if (yday(met_all$date[i])==243 & met_all$date[i]>spin_up_end){
    
    # reset tsub to annual mean tiz
    tsubiz_last<-tdist(ltm,temp_corr=temp_corr,lapse=tlapse,z)
    tsubiz_last[tsubiz_last<0]<-0
    
    # Cache SMB and clear
    annMB[k]<-cmb_ann
    yrs[k]<-year(met_all$date[i])
    
    # Also adjust glacier area and hypsometry
    dA<-(abs(cmb_ann)*totArea)^(1.0/va) # m^2
    
    # Identify lowest elevation with glacier area
    idx<-ref_idx[hyps[,1]==min(hyps[hyps[,2]>0,1])]
    
    if (cmb_ann < 0){
      
      # Subtract area from lowest bands -- and possibly above 
      while(dA>0){
        
        if (dA>hyps[idx,2]){
          
          resid<-dA-hyps[idx,2]
          hyps[idx,2]<-0
          idx<-idx+1
          dA<-resid*1.0
          
        } else {
          
          hyps[idx,2]<-hyps[idx,2]-dA
          dA<-0
          
        }
        
      }
      
      if (sum(hyps[,2])==0){
        print("Glacier disappeared!")
        break
      }  
      
    } else {
      
      # Add area below lowest band
      while(dA>0){
        
        if (hyps[idx,2]+dA>(hyps[(idx+1),2]*max_inc)){ 
          
          resid<-hyps[idx,2]+dA-(hyps[(idx+1),2]*max_inc)
          hyps[idx,2]<-hyps[(idx+1),2]*max_inc
          idx<-idx-1
          dA<-resid*1.0
          
        }
        
        else {
          hyps[idx,2]<-hyps[idx,2]+dA
          dA<-0
        }
        
      }
      
    }
    
    # Reset cumulative mass balance
    cmb_ann<-0
    
    # Store annual stats
    totArea<-sum(hyps[,2])
    annArea[k]<-totArea
    annT[k]<-totT/nd
    annP[k]<-totP/nd 
    
    # Increment counter
    k<-k+1
    
    # Reset nd and running counters
    nd<-0
    totT<-0
    totP<-0
    
    #Provide update
    print(sprintf("Year = %.0f, doy = %.0f --> Area = %.2f km**2",
                  year(met_all$date[i]),yday(met_all$date[i]),totArea/1e6))
    
  }
}
annframe<-data.frame("year"=yrs,"area"=annArea,"MB"=annMB,"temp"=annT,"precip"=annP)
print("--------------------MODEL RUN COMPLETE--------------------")
print(sprintf("Glacier area in year %.0f: %.2f km**2",yrs[ny],annArea[ny]/1e6))
print("----------------------------------------------------------")
