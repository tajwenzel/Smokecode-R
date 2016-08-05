solver.function <- function(t,init,pars,...) {
  B<-init[1];
  Sc<-init[2];
  Ec<-init[3];
  Sj<-init[4];
  Ej<-init[5];
  Sy<-init[6];
  Ey<-init[7];
  Non<-init[8];
  Ha<-init[9];
  Na<-init[10];
  Qa<-init[11];
  
  Pop<-B+Sc+Ec+Sj+Ej+Sy+Ey+Non+Ha+Na+Qa;
  
  
  with(as.list(combo), {
    
    dB<- ((Sc+Ec)*uc+(Sj+Ej)*uj+(Sy+Ey)*uy+(Non+Ha+Na+Qa)*ua+B*uc)-B*uc-B*lj
    
    dSc <- B*lj-(Sc*betacc*Ec)/Pop-(Sc*betacc*phicj*Ej)/Pop-(Sc*betacc*phica*(Ey+Ha+Na))/Pop-Sc*lc-Sc*uc;
    
    dEc <- (Sc*betacc*Ec)/Pop+(Sc*betacc*phicj*Ej)/Pop+(Sc*betacc*phica*(Ey+Ha+Na))/Pop-Ec*lc-Ec*uc;
    
    dSj <-(-Sj*betajj*Ej)/Pop-(Sj*betaja*(Ey+Ha+Na))/Pop-Sj*lj+Sc*lc-Sj*uj;
    
    dEj <- (Sj*betajj*Ej)/Pop+(Sj*betaja*(Ey+Ha+Na))/Pop-Ej*lj+Ec*lc-Ej*uj;
    
    dSy <- (-Sy*betayy*Ey)/Pop-(Sy*betaya*(Na+Ha))/Pop-Sy*ly-Sy*uy+Sj*lj;
    
    dEy <- (Sy*betayy*Ey)/Pop+(Sy*betaya*(Na+Ha))/Pop-Ey*ly+Ej*lj-Ey*uy;
    
    dNon<- Sy*ly-Non*ua;
    dHa <- Ey*ly*0.1-Ha*ua;
    dNa <- Ey*ly*0.9-q*Na-Na*ua;
    dQa <- q*Na-Qa*ua;
    
    #####ODE list
    
    stored <-list(c(dB,dSc,dEc,dSj,dEj,dSy,dEy,dNon,dHa,dNa,dQa))
    
  })
}