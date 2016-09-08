pars4year<-function(year)
{
  
  if(year==1965)
{parameters<-c(
    brate=0.00126,
    lc=0.125, 
    lj=0.143, 
    ly=0.25, 
    sigma=0.5, 
    q=(0.03798537), 
    ub=0.02559,
    uc=0.000423,
    uj=0.000423, 
    uy=0.001091, 
    ua=0.06551, 
    phica=runif(1,0.1,1), 
    phicj=runif(1,0.1,0.5))
}
 
  
  if(year==1970)
  {parameters<-c(
          brate=0.00117,
          lc=0.125, 
          lj=0.143, 
          ly=0.25, 
          sigma=0.5, 
          q=(0.04454015), 
          ub= 0.022269,
          uc = 0.000414, 
          uj = 0.00121, 
          uy = 0.00137, 
          ua = 0.0083, 
          phica=runif(1,0.1,1), 
          phicj=runif(1,0.1,0.5))
  }

  
  if(year==1975)
{parameters<-c(
         brate=0.0099,
         lc=0.125, 
         lj=0.143, 
         ly=0.25, 
         sigma=0.5, 
         q= , 
         uc = 0.000349, 
         uj = 0.00116, 
         uy = 0.0013, 
         ua = 0.0136,
         phica=runif(1,0.1,1), 
         phicj=runif(1,0.1,0.5)
)
  }
  
  
 
  
if(year==1980)
{parameters<-
  c(brate=0.0096,
         lc=0.125, 
         lj=0.143, 
         ly=0.25, 
         sigma=0.5, 
         q=(0.05357624), 
        uc = 0.0002817785, 
        uj = 0.001088055, 
        uy = 0.001134554, 
        ua = 0.007095612, 
         phica=runif(1,0.1,1), 
         phicj=runif(1,0.1,0.5))
}

 
  
if(year==1985)
{parameters<-c(
         brate=0.0089,
         lc=0.125, 
         lj=0.143, 
         ly=0.25, 
         sigma=0.5, 
         q=(0.05585817), 
         uc = 0.000262,
         uj = 0.000983, 
         uy = 0.001,
         ua = 0.00739,
         phica=runif(1,0.1,1), 
         phicj=runif(1,0.1,0.5))
}
 
if(year==1990)
{parameters<-c(
        brate=0.0114,
        lc=0.125, 
        lj=0.143, 
        ly=0.25, 
        sigma=0.5, 
        q=(0.05746801), 
        uc = 0.00024,
        uj = 0.000999,
        uy = 0.00111,
        ua = 0.00619,
        phica=runif(1,0.1,1), 
        phicj=runif(1,0.1,0.5))

}
  
  return(parameters)
}


