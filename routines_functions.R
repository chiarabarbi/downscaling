library(base)
library(dplyr)
library(fda)
library(geofd)
library(geoR)
library(gstat)
library(ggplot2)
library(xtable)
library(Metrics)
library(tidyr)
library(robCompositions)
library(splines)
library(fda.usc)
library(rgdal)
library(ggfortify)
library(sf)
library(atakrig)
library(readxl)
library(pracma)
require(rgdal)
library(raster)
library(rgeos)
library(viridis)
library(ggpubr)
library(latex2exp)
library(maptools)

setwd("~/Desktop/git_repos/Terna/clustering/downscaling_r")

# if work_with_cond == 1 -> analysis is done in the restricted domain 
# if work_with_cond == 0 -> analysis is done in the whole domain [0,200]

work_with_cond <- 1

step = 1   #passo di discretizzazione 


# dominio per il vento
min_dom = 0
max_dom = 200


##########  utility functions for interpolation

# function that given a dataframe (t,val) returns the linear interpolant evaluated at point xx

interp <- function(xx, tmp){
  approx(tmp$t, tmp$val, xout = xx)$y
}

smoothing <- function(xx, tmp, m){
  sum = 0 
  for( k in 0:m){
    F_v1 = interp(k/m,tmp)
    b_km = bernsteinb(k, m, xx)
    sum = sum + F_v1*b_km
  }
  sum
}

conditioning = function(density, start, stop)
{ 
  index = which(density$t >= start & density$t <= stop)
  data = density$val[index]
  step = density$t[2] - density$t[1]
  integral = trapzc(step,data)
  return(data/integral)
}

pdf_from_smoothing <- function(xx, tmp, m){
  sum = 0 
  for( k in 0:(m-1)){
    F_v1 = interp((k+1)/m,tmp)
    F_v1b =  interp(k/m,tmp)
    b_km1 = bernsteinb(k, m-1, xx)
    sum = sum + (F_v1-F_v1b)*b_km1
  }
  sum*m
  
}

squeeze <- function(x, domain){
  b <- max(domain)
  a <- min(domain)
  (x - a)/(b-a)
}

unsqueeze <- function(x, domain){
  b <- max(domain)
  a <- min(domain)
  x*(b-a)+a
}

trapzc = function(t_step,y)
{ # compute the integral of y with step "t_step"
  return(t_step*(0.5*y[1]+sum(y[2:(length(y)-1)]) +0.5*y[length(y)]))
}

normalize_density <- function(input_df, step){
  N_samples <- dim(input_df)[2]
  n <- dim(input_df)[1]
  ret <- matrix(NA, ncol = N_samples, nrow=n)
  
  for(i in 1:N_samples){
    ret[,i] = input_df[,i]/trapzc(step, input_df[,i])
  }
  ret
}



######### utility functions for fpca

clr = function(density, z, z_step)
{ # transform a density to a clr
  return(log(density)-trapzc(z_step,log(density))/(max(z)-min(z)))
}

clr.mv = function(density.df, z, z_step)
{ # transform a dataset of densities to clr
  N_samples = dim(density.df)[2]
  n=dim(density.df)[1]
  if(length(z)!=n)
  {
    N_samples = dim(density.df)[1]
    n=dim(density)[2]
    density=t(density)
  }
  res=matrix(NA, ncol = N_samples, nrow=n)
  for(i in 1:N_samples)
    res[,i]=log(density.df[,i])-trapzc(z_step,log(density.df[,i]))/(max(z)-min(z))
  return(res)
}


#clr.mv = function(density.df, z, z_step)
#{ # transform a dataset of densities to clr
#  N_samples = dim(density.df)[2]
#  n=dim(density.df)[1]
#  if(length(z)!=n)
#  {
#    N_samples = dim(density.df)[1]
#    n=dim(density)[2]
#    density=t(density)
#  }
#  res=matrix(NA, ncol = N_samples, nrow=n)
#  for(i in 1:N_samples)
#    if(length(which(density.df[,i] == 0)) == 0){
#      res[,i]=log(density.df[,i])-trapzc(z_step,log(density.df[,i]))/(max(z)-min(z))}
#    else{
#      idx = which(density.df[,i] == 0)
#      find_min = min(density.df[,i][which(density.df[,i]>0)])
#      density.df[idx,i] = find_min
#      res[,i]=log(density.df[,i])-trapzc(z_step,log(density.df[,i]))/(max(z)-min(z))
#      }
#  return(res)
#}

clr2density <- function(clr, z, z_step)
{ # back-transform a clr to a density
  if(is.fd(clr))
    return(exp(eval.fd(z,clr))/trapzc(z_step,exp(eval.fd(z,clr))))
  if(!is.fd(clr))
    return(exp(clr)/trapzc(z_step,exp(clr)))
}

clr2density.mv <- function(clr.df, z,z_step)
{ # transform a dataset of clr into densities
  N_samples <- dim(clr.df)[2]
  dens.df <- matrix(0, nrow=length(z), ncol=N_samples)
  for(j in 1:N_samples)
  {
    dens.df[,j]=clr2density(clr = clr.df[,j], z = z, z_step = z_step)
  }
  return(dens.df)
}





