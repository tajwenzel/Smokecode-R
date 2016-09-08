#DO NOT EDIT ANYMORE

#init<-initial.pop
#pars<-a85$trace[length(a85$trace[,1]),]
#solutions<-ode.solution
#year<-1980

llikelihood.wrap<-function(pars)
{
llikelihood.s <- function(pars,init,year,...)
  {
    ode.time <- seq(0, 3000, by = 1)
    
ode.eqs<-dget(file="/Users/Natasha/Dropbox/Smoking Results/Smokecode-R/FUNC_smoke_ode_solver.R", keep.source=TRUE)

combo<<-all.pars(pars,year)

###############################################################
    ## REMEBER PARMS IS NOT THE SAME AS PARS. PARS=UNKNOWN PARAMETERS, PARMS IS ALL PARAMETERS AS GIVEN BY    COMBO
###############################################################
   #init=initial.pop
   #combo<<-all.pars(beta.start,1980)

   #memory allocation
    ode.out<-as.data.frame(rk4(y=init, times= ode.time, func=ode.eqs, parms=combo, rtol=1e-4,atol=1e-4))
    
    #matplot(ode.out,type='l')
  
    
    solutions<-ode.out[2999,];
    
    with(as.list(c(combo,solutions)), 
         {
    pop<-Sc+Ec+Sj+Ej+Sy+Ey+Non+Ha+Na+B;
    
    prevalenceC<- Ec/(Ec + Sc)
    prevalenceJ<- Ej/(Ej + Sj)
    prevalenceY<- Ey/(Ey + Sy)
    prevalenceA<- (Ha + Na)/(Ha + Na + Non + Qa) 
    truequit<- Qa/(Na + Ha + Qa) 
    fic <-(betacc*Ec+betacc*phicj*Ej+betacc*phica*(Ey+Ha+Na))/pop 
    fij <-(betajj*Ej+betaja*(Ey+Ha+Na))/pop 
    fiy <-(betayy*Ey+betaya*(Na+Ha))/pop
    
            
            
        if(
          (any(combo[1:length(combo)] < 0) || any(combo[1:length(combo)] > 1) || prevalenceA <= 0 
            || prevalenceA >= 1 || prevalenceY <= 0 || prevalenceY >= 1 || prevalenceJ <= 0 
            || prevalenceJ >= 1 || prevalenceC <=0 || prevalenceC >= 1 || truequit <= 0 
            || truequit >= 1 || fic <= 0 || fic>=1 || fij <= 0 || fij>=1 || fiy<=0 || fiy>=1)
          ==TRUE)
          
        {LL<<-(-Inf)}
       
     else   #return calculated likelihood for the appropriate year
  
     
     {#bracket is part of the else if, not part of year division
       if(year==1965) #fi checked as of 08/13/2016, prev checked on 08/14
            {LL<<-dbeta(prevalenceC,150, 8159, log=TRUE) + 
                  dbeta(prevalenceJ,1625, 3924,log=TRUE) + 
                  dbeta(prevalenceY,1281, 1147,log=TRUE) + 
                  dbeta(prevalenceA,7260, 9848,log=TRUE) + 
                  dbeta(truequit,93, 422,log=TRUE) + 
                  dbeta(fic,205, 20790,log=TRUE) + 
                  dbeta(fij,612, 14026,log=TRUE) + 
                  dbeta(fiy,80, 7188,log=TRUE);
            }
      
     if(year==1970) #fi checked as of 08/13/2016 prev checked on 08/14
    {LL<<-dbeta(prevalenceC,166,8129 , log=TRUE) + 
      dbeta(prevalenceJ,1951,4980 ,log=TRUE) + 
      dbeta(prevalenceY,1496,1633 ,log=TRUE) + 
      dbeta(prevalenceA,8188,11971 ,log=TRUE) + 
      dbeta(truequit,57, 365,log=TRUE) + 
      dbeta(fic,160, 21811,log=TRUE) + 
      dbeta(fij,689, 14676,log=TRUE) + 
      dbeta(fiy,82, 7484,log=TRUE);
    }
    
      
    if(year==1975) ##fi checked as of 08/13/2016 prev checked on 08/14
    {LL<<-dbeta(prevalenceC,220, 6729, log=TRUE) + 
      dbeta(prevalenceJ,2406,4998 ,log=TRUE) + 
      dbeta(prevalenceY,1612,2299 ,log=TRUE) + 
      dbeta(prevalenceA,9421,14693 ,log=TRUE) + 
      dbeta(truequit,2516, 8192,log=TRUE) + 
      dbeta(fic, 167, 22792,log=TRUE) + 
      dbeta(fij, 735,15299,log=TRUE) + 
      dbeta(fiy,91,7795,log=TRUE);
    }
      
    if(year==1980) ##fi checked as of 08/13/2016 prev checked on 08/14
    {LL<<-dbeta(prevalenceC,117,4937, log=TRUE) + 
      dbeta(prevalenceJ,1276,3568,log=TRUE) + 
      dbeta(prevalenceY,1015,1637,log=TRUE) + 
      dbeta(prevalenceA,4961,9959  ,log=TRUE) + 
      dbeta(truequit,61, 493,log=TRUE) + 
      dbeta(fic, 206,24022,log=TRUE) + 
      dbeta(fij,804,15099,log=TRUE) + 
      dbeta(fiy,69,8090,log=TRUE);
    }
       
    if(year==1985) ##fi checked as of 08/13/2016 prev checked on 08/14
    {LL<<-dbeta(prevalenceC,80,4506 , log=TRUE) + 
      dbeta(prevalenceJ, 1058,3574 ,log=TRUE) + 
      dbeta(prevalenceY,978,1797 ,log=TRUE) + 
      dbeta(prevalenceA,5732,12496 ,log=TRUE) + 
      dbeta(truequit,61, 493,log=TRUE) + 
      dbeta(fic,155, 23445,log=TRUE) + 
      dbeta(fij,673, 15656,log=TRUE) + 
      dbeta(fiy,70, 8389,log=TRUE);
    }
       
    if(year==1990) ##fi checked as of 08/13/2016
    {LL<<-dbeta(prevalenceC,85,3377, log=TRUE) + 
      dbeta(prevalenceJ,910,3381 ,log=TRUE) + 
      dbeta(prevalenceY,782,1832 ,log=TRUE) + 
      dbeta(prevalenceA,6257,15449 ,log=TRUE) + 
      dbeta(truequit,31, 222,log=TRUE) + 
      dbeta(fic,81,18755,log=TRUE) + 
      dbeta(fij,433, 13230,log=TRUE) + 
      dbeta(fiy,6,7458,log=TRUE);}
    
     }     
       })
         
   return(LL)
  }
  return(llikelihood.s)
}


#####MLE for starting parameters

