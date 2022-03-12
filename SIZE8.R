#  Use of this program and definitions of the parameters are described in
# the paper:
#   
#   Sample Size Calculation for Complex Clinical Trials
# with Survival Endpoints by J. Shih (1995), Controlled Clinical Trials,
# Vol. 16, pp. 395-407.



#   This SAS program has been tested and checked for accuracy.  The
# author of this program, however, claims no responsibility for any errors
# that may arise.

  
####################################################

ptm1<-proc.time()

fun.seq<-function(char,num){
  return(eval(as.symbol(paste0(char,num))))
}

for (i in 2:4){
  assign(paste0("u",i),i^2)
  print(fun.seq("u",i))
}

par2<-c()

# Define two Markov process functions
#################################
markov<-function(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,r1,r2,logr){
  if (length(pc>1)){pc<-matrix(pc,nrow=1)}
  loss<-matrix(loss,nrow=1)
  noncmpl<-matrix(noncmpl,nrow=1)
  dropin<-matrix(dropin,nrow=1)
  k<-matrix(k,ncol=1)
  if (prop == 1 | ((prop == 0 )& (lag == 0))){
    distr_e=matrix(c(0,0,1,0),ncol=1); distr_c=matrix(c(0,0,0,1),ncol=1)
    # rm(dstr_e, dstr_c)
    trans = diag(4)
    dstr_e=c()
    dstr_c=c()
    for (year in 1:periods){
      ls=1-(1-loss[,year])**(1/n_intrvl)
      dro=1-(1-noncmpl[,year])**(1/n_intrvl)
      dri=1-(1-dropin[,year])**(1/n_intrvl)
      pc1=1-(1-pc[,year])**(1/n_intrvl)
      if (prop == 1){
        pe1 = 1-exp(k*log(1-pc1))
      }else{
        pe1 = 1-exp(k[year,]*log(1-pc1))
      }
      for(ii in 1:n_intrvl){
        trans[,3]=rbind(ls,pe1,(1-(ls+pe1+dro)),dro)
        trans[,4]=rbind(ls,pc1,dri,(1-(ls+pc1+dri)))
        distr_e=trans%*%distr_e; distr_c=trans%*%distr_c
        
        # Adjust for staggered entry
        temp_e=distr_e[c(3, 4),1]*(1-ad_cens[,ii+(year-1)*n_intrvl])
        distr_e[1,]=distr_e[1,]+apply((distr_e[c(3, 4),,drop=F]-temp_e),2,sum)
        distr_e[c(3, 4),]=temp_e
        temp_c=distr_c[c(3, 4),1]*(1-ad_cens[,ii+(year-1)*n_intrvl])
        distr_c[1,]=distr_c[1,]+apply((distr_c[c(3, 4),,drop=F]-temp_c),2,sum)
        distr_c[c(3, 4),]=temp_c
        #end of transition matrix loop;   
      }
      
      dstr_e=cbind(dstr_e,distr_e) 
      dstr_c=cbind(dstr_c,distr_c)
      #end of years loop;   
    }
    # print dstr_e; print dstr_c; 
    
    #end if clause for prop hazards or nonlag models; 
  }else{# /* nonproportional hazards and lag is specified */
    nactv = lag*n_intrvl
    nstates=2*nactv+2
    distr_e = rbind(0,0,1,matrix(0,nactv-1,1), matrix(0,nactv,1))
    distr_c = rbind(0,0,matrix(0,nactv,1), 1, matrix(0,nactv-1,1))
    # rm(dstr_e, dstr_c)
    dstr_e=c()
    dstr_c=c()
    for (year in 1:periods){
      ls=1-(1-loss[,year])**(1/n_intrvl)*matrix(1,1,nactv)
      dro=1-(1-noncmpl[,year])**(1/n_intrvl)*matrix(1,1,nactv)
      dri=1-(1-dropin[,year])**(1/n_intrvl)*matrix(1,1,nactv)
      pc1=1-(1-pc[,year])**(1/n_intrvl)
      pe1 = matrix(0,1,nactv)
      for (l in 1:lag){
        pe1[,((l-1)*n_intrvl+1):(l*n_intrvl)] = 1-(exp(k[l,]*log(1-pc[,year])*
                                                         matrix(1,1,n_intrvl)))**(1/n_intrvl)
      }
      actv = cbind(pc1,pe1)
      a = matrix(0,nactv,nactv); b=a; d=a
      #Start transition creationn and multiplication loop;
      if (lagdout == 1) {
        c = diag(as.vector(dro))
      }else {
        c = rbind(dro,matrix(0,nactv-1,nactv))
      }
      b= diag(as.vector(dri))
      a[2:nactv,1:(nactv-1)] =
        diag(1-((ls+dro+actv[,2:(nactv+1)])[,1:(nactv-1)]))
      d[1:(nactv-1),2:nactv] =
        diag(1-(((ls+dri)[,1:(nactv-1)]+actv[,2:nactv])))
      a[nactv,nactv] = 1-(ls[,1]+dro[,1]+actv[,nactv+1])
      d[1,1] = 1-((ls+dri)[,1]+actv[,1])
      for (ii in 1:n_intrvl){
        trans=rbind(cbind(diag(2),rbind(cbind(ls,ls),cbind(actv[,2:(nactv+1),drop=F],actv[,1:(nactv),drop=F])))
                    ,cbind(matrix(0,2*nactv,2),rbind(cbind(a,b),cbind(c,d))))
        distr_e=trans%*%distr_e
        distr_c=trans%*%distr_c
        #Adjust for staggered entry;
        temp_e=distr_e[3:nstates,1]*(1-ad_cens[,ii+(year-1)*n_intrvl])
        distr_e[1,]=distr_e[1,]+apply(distr_e[3:nstates,,drop=F]-temp_e,2,sum)
        distr_e[3:nstates,]=temp_e
        temp_c=distr_c[3:nstates,1]*(1-ad_cens[,ii+(year-1)*n_intrvl])
        distr_c[1,]=distr_c[1,]+apply(distr_c[3:nstates,,drop=F]-temp_c,2,sum)
        distr_c[3:nstates,]=temp_c
        #end of transition matrix loop;   
      }
      dstr_e=cbind(dstr_e,distr_e)  
      dstr_c=cbind(dstr_c,distr_c)
      #end of years loop;   
    }
    collaps1=3:(nactv+2); collaps2=(nactv+3):(2*nactv+2)
    dstr_e=rbind(dstr_e[1:2,],apply(dstr_e[collaps1,,drop=F],2,sum),apply(dstr_e[collaps2,,drop=F],2,sum))
    dstr_c=rbind(dstr_c[1:2,],apply(dstr_c[collaps1,,drop=F],2,sum),apply(dstr_c[collaps2,,drop=F],2,sum))
    #  print dstr_e; print dstr_c
  }# /*end nonproportional lag model */
  
  if (periods == 1){
    event_c=dstr_c[2,]
    event_e=dstr_e[2,]
    loss_c=dstr_c[1,]
    loss_e=dstr_e[1,]
  }else{
    event_c=dstr_c[2,]-c(0,dstr_c[2,1:(ncol(dstr_c)-1)])
    event_e=dstr_e[2,]-c(0,dstr_e[2,1:(ncol(dstr_e)-1)])
    loss_c=dstr_c[1,]-c(0,dstr_c[1,1:(ncol(dstr_c)-1)])
    loss_e=dstr_e[1,]-c(0,dstr_e[1,1:(ncol(dstr_e)-1)])
  }
  atrisk_c=apply(dstr_c[3:4,,drop=F],2,sum)+loss_c+event_c
  atrisk_e=apply(dstr_e[3:4,,drop=F],2,sum)+loss_e+event_e
  phi=atrisk_c*r1/(atrisk_e*r2)
  theta=log(1-event_c/atrisk_c)/log(1-event_e/atrisk_e)
  
  # theta=(event_c/atrisk_c)/(event_e/atrisk_e)
  # rho=(event_c*r1+event_e*r2)/(apply(event_c*r1+event_e*r2,1,sum))
  rho=(event_c*r1+event_e*r2)/(sum(event_c*r1+event_e*r2))
  gamma=abs(phi*theta/(1+phi*theta)-phi/(1+phi))
  eta=phi/((1+phi)**2)
  
  # if (((length(logr)==1) & (logr[1] != 0))|(length(logr)>1))
    if(sum(logr==0)==0){
    # sig=sqrt(apply((logr**2)*rho*eta,1,sum))
    sig=sqrt(sum((as.vector(logr)**2)*rho*eta))
    sig2 = sum(as.vector(logr)*rho*gamma)
    ed = sig2/sig
  }else{ed=NA}
  
  
  pe2=distr_e[2,]; pc2=distr_c[2,];  pbar=(r1*pc2+r2*pe2)/(r1+r2)
  dropin = dstr_c[3,periods]
  noncmpl = dstr_e[4,periods]
  loss = (r1*dstr_c[1,periods]+r2*dstr_e[1,periods])/(r1+r2)
  z_alpha=qnorm(1-alpha/2)
  
  result=list(z_alpha, pe2, pc2, pbar,ed, dropin,noncmpl,loss)
  names(result)<-c("z_alpha", "pe2", "pc2", "pbar","ed", "dropin","noncmpl", "loss")
  return(result)
}

markov2<-function(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,nptmp2,nptmp,nn,r1,r2,logr){# /* markov chain model with unknown duration */
  # rm(dstr_e, dstr_c)
  if (length(pc>1)){pc<-matrix(pc,nrow=1)}
  loss<-matrix(loss,nrow=1)
  noncmpl<-matrix(noncmpl,nrow=1)
  dropin<-matrix(dropin,nrow=1)
  k<-matrix(k,ncol=1)
  dstr_e=c()
  dstr_c=c()
  
  if ((prop == 1) | ((prop == 0) & (lag == 0))){
    distr_e=matrix(c(0,0,1,0),ncol=1)        
    distr_c=matrix(c(0,0,0,1),ncol=1)
    trans = diag(4)
    for (year in 1:nptmp2){
      ls=1-(1-loss[,year])**(1/n_intrvl)
      dro=1-(1-noncmpl[,year])**(1/n_intrvl)
      dri=1-(1-dropin[,year])**(1/n_intrvl)
      pc1=1-(1-pc[,year])**(1/n_intrvl)
      if (prop == 1){
        pe1 = 1-exp(k*log(1-pc1))
      }else{
        pe1 = 1-exp(k[year,]*log(1-pc1))
      }
      
      if ((year == nptmp2) & (nptmp2 != nptmp)){
        ni = nn%%n_intrvl
        if (ni == 0){ni = n_intrvl}
      }else{
        ni = n_intrvl
      }
      for (ii in 1:ni){
        trans[,3]=rbind(ls,pe1,(1-(ls+pe1+dro)),dro)
        trans[,4]=rbind(ls,pc1,dri,(1-(ls+pc1+dri)))
        distr_e=trans%*%distr_e; distr_c=trans%*%distr_c
        #Adjust for staggered entry;
        temp_e=distr_e[3:4,1]*(1-ad_cens[,ii+(year-1)*ni])
        distr_e[1,]=distr_e[1,]+apply(distr_e[3:4,,drop=F]-temp_e,2,sum)
        distr_e[3:4,]=temp_e
        temp_c=distr_c[3:4,1]*(1-ad_cens[,ii+(year-1)*ni])
        distr_c[1,]=distr_c[1,]+apply(distr_c[3:4,,drop=F]-temp_c,2,sum)
        distr_c[3:4,]=temp_c
        #end of transition matrix loop;   
      }
      dstr_e=cbind(dstr_e,distr_e)
      dstr_c=cbind(dstr_c,distr_c)
      #end of years loop;   
    }
    print (cbind(dstr_e)); print (cbind(dstr_c))
  }else{ #/* nonproportional hazards and lag is specified */
    nactv = lag*n_intrvl
    nstates=2*nactv+2
    distr_e = rbind(0,0,1/matrix(0,nactv-1,1),matrix(0,nactv,1))
    distr_c = rbind(0,0,matrix(0,nactv,1),1,matrix(0,nactv-1,1))
    for (year in 1:nptmp2){
      ls=1-(1-loss[,year])**(1/n_intrvl)*matrix(1,1,nactv)
      dro=1-(1-noncmpl[,year])**(1/n_intrvl)*matrix(1,1,nactv)
      dri=1-(1-dropin[,year])**(1/n_intrvl)*matrix(1,1,nactv)
      pc1=1-(1-pc[,year])**(1/n_intrvl)*matrix(1,1,n_intrvl)
      pe1 = matrix(0,1,nactv)
      for (l in 1:lag){
        pe1[,((l-1)*n_intrvl+1):(l*n_intrvl)] = 1-(exp(k[l,]*log(1-pc[,year])*
                                                         matrix(1,1,n_intrvl)))**(1/n_intrvl)
      }
      actv = cbind(pc1[,1],pe1)
      a = matrix(0,nactv,nactv); b=a; d=a
      # Start transition creationn and multiplication loop;
      if (lagdout == 1) {
        c = diag(as.vector(dro))
      }else{ c = rbind(dro,matrix(0,nactv-1,nactv))}
      b= diag(as.vector(dri))
      a[2:nactv,1:(nactv-1)] =
        diag(1-((ls+dro+actv[,2:(nactv+1)])[,1:(nactv-1)]))
      d[1:(nactv-1),2:nactv] =
        diag(1-(((ls+dri)[,1:(nactv-1)]+actv[,2:nactv])))
      a[nactv,nactv] = 1-(ls[,1]+dro[,1]+actv[,nactv+1])
      d[1,1] = 1-((ls+dri)[,1]+actv[,1])
      
      if ((year == nptmp2) & (nptmp2 != nptmp)){
        ni = nn%%n_intrvl
        if (ni == 0) {ni = n_intrvl}
      }else{
        ni = n_intrvl
      }
      for (ii in 1:ni){
        trans=rbind(cbind(diag(2),rbind(cbind(ls,ls),cbind(actv[,2:(nactv+1),drop=F],actv[,1:(nactv),drop=F])))
                    ,cbind(matrix(0,2*nactv,2),rbind(cbind(a,b),cbind(c,d))))
        # trans=rbind(cbind(diag(2),rbind(cbind(ls,ls),cbind(actv[,2:(nactv+1),drop=F],actv[,1:(nactv),drop=F])))
        #             ,cbind(matrix(0,2*nactv,2),rbind(cbind(a,b),cbind(c,d))))
        distr_e=trans%*%distr_e; distr_c=trans%*%distr_c
        # Adjust for staggered entry;
        temp_e=distr_e[3:nstates,1]*(1-ad_cens[,ii+(year-1)*n_intrvl])
        distr_e[1,]=distr_e[1,]+apply(distr_e[3:nstates,,drop=F]-temp_e,2,sum)
        distr_e[3:nstates,]=temp_e
        temp_c=distr_c[3:nstates,1]*(1-ad_cens[,ii+(year-1)*n_intrvl])
        distr_c[1,]=distr_c[1,]+apply(distr_c[3:nstates,,drop=F]-temp_c,2,sum)
        distr_c[3:nstates,]=temp_c
        # end of transition matrix loop;   
      }
      dstr_e=cbind(dstr_e,distr_e)
      dstr_c=cbind(dstr_c,distr_c)
      # end of years loop;   
    }
    collaps1=3:(nactv+2); collaps2=(nactv+3):(2*nactv+2);
    dstr_e=rbind(dstr_e[1:2,],apply(dstr_e[collaps1,,drop=F],2,sum),apply(dstr_e[collaps2,,drop=F],2,sum))
    dstr_c=rbind(dstr_c[1:2,],apply(dstr_c[collaps1,,drop=F],2,sum),apply(dstr_c[collaps2,,drop=F],2,sum))
    # print dstr_e; print dstr_c; 
  }# /*end nonproportional lag model */
  
  if (nptmp2 == 1){
    event_c=dstr_c[2,]
    event_e=dstr_e[2,]
    loss_c=dstr_c[1,]-c(0,dstr_c[1,1:(ncol(dstr_c)-1)])
    loss_e=dstr_e[1,]-c(0,dstr_e[1,1:(ncol(dstr_e)-1)])
  }else{
    event_c=dstr_c[2,]-c(0,dstr_c[2,1:(ncol(dstr_c)-1)])
    event_e=dstr_e[2,]-c(0,dstr_e[2,1:(ncol(dstr_e)-1)])
    loss_c=dstr_c[1,]-c(0,dstr_c[1,1:(ncol(dstr_c)-1)])
    loss_e=dstr_e[1,]-c(0,dstr_e[1,1:(ncol(dstr_e)-1)])
  }
  atrisk_c=apply(dstr_c[3:4,],2,sum)+loss_c+event_c
  atrisk_e=apply(dstr_e[3:4,],2,sum)+loss_e+event_e
  phi=atrisk_c*r1/(atrisk_e*r2)
  # theta=(event_c/atrisk_c)/(event_e/atrisk_e);
  
  theta=log(1-event_c/atrisk_c)/log(1-event_e/atrisk_e)
  rho=(event_c*r1+event_e*r2)/(sum(event_c*r1+event_e*r2))
  gamma=abs(phi*theta/(1+phi*theta)-phi/(1+phi))
  eta=phi/((1+phi)**2)
  # if ((length(logr)==1 & logr[1] != 0)|(length(logr)>1))
    if(sum(logr==0)==0){
    sig=sqrt(sum((as.vector(logr)**2)*rho*eta))
    sig2 = sum(as.vector(logr)*rho*gamma)
    ed = sig2/sig
  }else{ed=NA}
  pe2=distr_e[2,]; pc2=distr_c[2,];  pbar=(r1*pc2+r2*pe2)/(r1+r2)
  dropin = dstr_c[3,nptmp2]
  noncmpl = dstr_e[4,nptmp2]
  loss = (r1*dstr_c[1,nptmp2]+r2*dstr_e[1,nptmp2])/(r1+r2)
  z_alpha=qnorm(1-alpha/2)
  result=list(z_alpha, pe2, pc2, pbar,ed, dropin,noncmpl,loss)
  names(result)<-c("z_alpha", "pe2", "pc2", "pbar","ed", "dropin","noncmpl", "loss")
  return(result)
}

# Three core functions for power, size and duration
######################################                
power.fun<-function(m1,m2,m3,m4,m5,
                    mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                    nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                    random1,random2,random3,random4,random5,
                    od1,od2,od3,od4,od5,
                    odi1,odi2,odi3,odi4,odi5,
                    alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,r1,r2,logr
){
  fun.seq<-function(char,num){
    return(eval(as.symbol(paste0(char,num))))
  }
  par2<-c()
  for  (iv6 in 1:nomp6){
    n=mp6[iv6,]
    meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
    for (iv1 in 1:fun.seq("nomp",odi1)){
      meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
      for (iv2 in 1:fun.seq("nomp",odi2)){
        meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
        for (iv3 in 1:fun.seq("nomp",odi3)){
          meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
          for (iv4 in 1:fun.seq("nomp",odi4)){
            meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
            for (iv5 in 1:fun.seq("nomp",odi5)){
              pc = mp2[fun.seq("iv",od2),]
              if (prop == 1){
                k = mp1[fun.seq("iv",od1),]
              }else{
                k = mp1[((fun.seq("iv",od1)-1)*periods+1):(fun.seq("iv",od1)*periods),]
              }
              print(mp1)
              print(((fun.seq("iv",od1)-1)*periods+1):(fun.seq("iv",od1)*periods))
              print(k)
              dropin = matrix(mp3[fun.seq("iv",od3),],nrow=1)
              noncmpl = matrix(mp4[fun.seq("iv",od4),],nrow=1)
              loss = matrix(mp5[fun.seq("iv",od5),],nrow=1)
              markov.result1<-markov(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,r1,r2,logr)
              # /* invoke markov chain model */
              ed=markov.result1$ed
              pc2=markov.result1$pc2
              pe2=markov.result1$pe2
              z_alpha=markov.result1$z_alpha
              dropin=markov.result1$dropin
              noncmpl=markov.result1$noncmpl
              loss=markov.result1$loss
              temp1 = mp1[fun.seq("iv",od1),]; temp3=dropin; temp4=noncmpl
              temp5 = loss
              rm(pc, dropin,  noncmpl, loss, k)
              p = (r1*pc2+r2*pe2)/(r1+r2)
              if (((length(logr)==1) & (logr[1] != 0))|(length(logr)>1)){
                dth = n*p
                z_beta = (dth**0.5)*ed-z_alpha
              }else{
                sd1 = ((r1+r2)*p*(1-p)/(r1*r2))**0.5;
                sd2 = (pc2*(1-pc2)/r1+pe2*(1-pe2)/r2)**0.5;
                z_beta = (((n**0.5)*(abs(pc2-pe2))/((r1+r2)**0.5)-z_alpha*sd1)/sd2) ;
                dth = n*p ;
              }
              dth2 = ceiling(dth)
              power = pnorm(z_beta)
              
              if (prop == 0) { temp1 =NA}
              if (random5 == 0){
                par = cbind(alpha,power,temp1,temp3,temp4,temp5,pc2,pe2,
                            n,dth2,periods)
                par2 = rbind(par2,par)
              }
              meanpow1 = meanpow1+power*fun.seq("m",odi5)[iv5,];
              meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,];
              meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,];
              meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,];
            }# /*v5 */
            
            if ((random4 == 0) & (random5 == 1)){
              assign(paste0("temp",odi5),NA)
              # fun.seq("temp",odi5)
              dth2 = ceiling(meandth1)
              par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                          n,dth2,periods)
              par2 = rbind(par2,par)
            }
            meanpow2 = meanpow2+meanpow1*fun.seq("m",odi4)[iv4,]
            meandth2 = meandth2+meandth1*fun.seq("m",odi4)[iv4,]
            meanpc2 = meanpc2+meanpc1*fun.seq("m",odi4)[iv4,]
            meanpe2 = meanpe2+meanpe1*fun.seq("m",odi4)[iv4,]
          }# /*v4*/
          
          if ((random3 == 0) & (random4 == 1)){
            assign(paste0("temp",odi5),NA); assign(paste0("temp",odi4), NA)
            dth2 = ceiling(meandth2)
            par = cbind(alpha,meanpow2,temp1,temp3,temp4,temp5,pc2,pe2,
                        n,dth2,periods)
            par2 = rbind(par2,par)
          }
          meanpow3 = meanpow3+meanpow2*fun.seq("m",odi3)[iv3,]
          meandth3 = meandth3+meandth2*fun.seq("m",odi3)[iv3,]
          meanpc3 = meanpc3+meanpc2*fun.seq("m",odi3)[iv3,]
          meanpe3 = meanpe3+meanpe2*fun.seq("m",odi3)[iv3,]
        }# /*v3*/
        
        if ((random2 == 0) & (random3 == 1)){
          assign(paste0("temp",odi5), NA)
          assign(paste0("temp",odi4), NA)
          assign(paste0("temp",odi3), NA)
          dth2 = ceiling(meandth3)
          par = cbind(alpha,meanpow3,temp1,temp3,temp4,temp5,pc2,pe2,
                      n,dth2,periods)
          par2 = rbind(par2,par)
        }
        meanpow4 = meanpow4 +meanpow3*fun.seq("m",odi2)[iv2,]
        meandth4 = meandth4 +meandth3*fun.seq("m",odi2)[iv2,]
        meanpc4 = meanpc4 +meanpc3*fun.seq("m",odi2)[iv2,]
        meanpe4 = meanpe4 +meanpe3*fun.seq("m",odi2)[iv2,]
      }# /*v2*/
      
      if ((random1 == 0) & (random2 == 1)){
        assign(paste0("temp",odi5), NA)
        assign(paste0("temp",odi4), NA)
        assign(paste0("temp",odi3), NA)
        assign(paste0("temp",odi2), NA)
        dth2 = ceiling(meandth4)
        par = cbind(alpha,meanpow4,temp1,temp3,temp4,temp5,pc2,pe2,
                    n,dth2,periods)
        par2 = rbind(par2,par)
      }
      meanpow5 = meanpow5 +meanpow4*fun.seq("m",odi1)[iv1,]
      meandth5 = meandth5 +meandth4*fun.seq("m",odi1)[iv1,]
      meanpc5 = meanpc5 +meanpc4*fun.seq("m",odi1)[iv1,]
      meanpe5 = meanpe5 +meanpe4*fun.seq("m",odi1)[iv1,]
    }# /*v1*/
    if (random1 == 1){
      assign(paste0("temp",odi5), NA)
      assign(paste0("temp",odi4), NA)
      assign(paste0("temp",odi3), NA)
      assign(paste0("temp",odi2), NA)
      assign(paste0("temp",odi1), NA)
      dth2 = as.integer(meandth5+0.5)
      par = cbind(alpha,meanpow5,temp1,temp3,temp4,temp5,meanpc5,
                  meanpe5,n,dth2,periods)
      par2 = rbind(par2,par)
    }
  }# /*v6*/
  return(par2)
}

sample<-function(m1,m2,m3,m4,m5,
                 mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                 nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                 random1,random2,random3,random4,random5,
                 od1,od2,od3,od4,od5,
                 odi1,odi2,odi3,odi4,odi5,
                 alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,
                 nmax.g,nmin.g,ntol.g,r1,r2,logr){
  fun.seq<-function(char,num){
    return(eval(as.symbol(paste0(char,num))))
  }
  sum.rand = random1+random2+random3+random4+random5
  par2<-c()
  if (sum.rand == 0){
    # all the parameters are fixed, not random variables
    for (iv7 in 1:nomp7){# /*power*/
      power = mp7[iv7,]
      for (iv1 in 1:fun.seq("nomp",odi1)){
        for (iv2 in 1:fun.seq("nomp",odi2)){
          for (iv3 in 1:fun.seq("nomp",odi3)){
            for (iv4 in 1:fun.seq("nomp",odi4)){
              for (iv5 in 1:fun.seq("nomp",odi5)){
                if (prop == 1){
                  k = mp1[fun.seq("iv",od1),]
                }else{
                  k = mp1[((fun.seq("iv",od1)-1)*periods+1):(fun.seq("iv",od1)*periods),]
                }
                pc = matrix(mp2[fun.seq("iv",od2),],nrow=1)
                dropin = matrix(mp3[fun.seq("iv",od3),],nrow=1)
                noncmpl = matrix(mp4[fun.seq("iv",od4),],nrow=1)
                loss = matrix(mp5[fun.seq("iv",od5),],nrow=1)
                
                markov.result2<-markov(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,r1,r2,logr)
                ed=markov.result2$ed
                pc2=markov.result2$pc2
                pe2=markov.result2$pe2
                z_alpha=markov.result2$z_alpha
                dropin=markov.result2$dropin
                noncmpl=markov.result2$noncmpl
                loss=markov.result2$loss
                temp1 = mp1[fun.seq("iv",od1),]; temp3 = dropin; temp4 = noncmpl
                temp5 = loss
                rm(pc, dropin, noncmpl, loss, k)
                p = (r1*pc2+r2*pe2)/(r1+r2)
                if (((length(logr)==1) & (logr[1] != 0))|(length(logr)>1)){
                  z_beta=qnorm(power)
                  dth=((z_alpha+z_beta)/ed)**2
                  n=dth/p
                }else{
                  z_beta=qnorm(power)
                  sd1 = ((r1+r2)*p*(1-p)/(r1*r2))**0.5
                  sd2 = (pc2*(1-pc2)/r1+pe2*(1-pe2)/r2)**0.5
                  n = ((z_alpha*sd1+z_beta*sd2)/(pc2-pe2))**2 
                  n = (r1+r2)*n
                  dth = n*p 
                }
                dth2 = ceiling(dth)
                n2 = ceiling(n)
                
                if (prop == 0) {temp1 =NA}
                par = cbind(alpha,power,temp1,temp3,temp4,temp5,pc2,pe2,
                            n2,dth2,periods)
                par2 = rbind(par2,par)
              }# /*iv5*/
            }# /*iv4*/
          }# /*iv3*/
        }# /*iv2*/
      }# /*iv1*/
    }#/*iv7*/
  }
  
  if (sum.rand != 0){# /*there is heterogeneity in the parameters */
    # Use interval halving method to find n s.t.
    # predicted power >= specified one. 
    for (iv7 in 1:nomp7){# /* power */
      power = mp7[iv7,]
      nmax = nmax.g; nmin = nmin.g; tol = ntol.g
      it = 0; n = (nmin+nmax)/2
      
      fun.iter5<-function(m1,m2,m3,m4,m5,
                          mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                          nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                          random1,random2,random3,random4,random5,
                          od1,od2,od3,od4,od5,
                          odi1,odi2,odi3,odi4,odi5,
                          alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                          nmax.g,nmin.g,ntol.g,
                          nmax,nmin,tol,power,n,iv1,iv2,iv3,iv4){
        # iter5:
        meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        for (iv5 in 1:fun.seq("nomp",odi5)){
          if (prop == 1) {
            k = mp1[fun.seq("iv",od1),]
          }else{ 
            # k = mp1
            k = mp1[((fun.seq("iv",od1)-1)*periods+1):(fun.seq("iv",od1)*periods),]
            }
          pc = mp2[fun.seq("iv",od2),]
          dropin = mp3[fun.seq("iv",od3),]
          noncmpl = mp4[fun.seq("iv",od4),]
          loss = mp5[fun.seq("iv",od5),]
          markov.result3<-markov(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,r1,r2,logr)
          ed=markov.result3$ed
          pc2=markov.result3$pc2
          pe2=markov.result3$pe2
          z_alpha=markov.result3$z_alpha
          dropin=markov.result3$dropin
          noncmpl=markov.result3$noncmpl
          loss=markov.result3$loss
          temp1 = mp1[fun.seq("iv",od1),]; temp3=dropin; temp4=noncmpl
          temp5 = loss
          p = (r1*pc2+r2*pe2)/(r1+r2)
          if (logr > 0){
            dth = n*p
            z_beta = (dth**0.5)*ed-z_alpha
          }else{
            sd1 = ((r1+r2)*p*(1-p)/(r1*r2))**0.5
            sd2 = (pc2*(1-pc2)/r1+pe2*(1-pe2)/r2)**0.5
            z_beta = (((n**0.5)*abs((pc2-pe2))/((r1+r2)**0.5)-z_alpha*sd1)/sd2)
            dth = n*p
          }
          newpow = pnorm(z_beta)
          meanpow1 = meanpow1+newpow*fun.seq("m",odi5)[iv5,]
          meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,]
          meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,]
          meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,]
        }# /*v5*/
        result<-list(meanpow1,meandth1,meanpc1,meanpe1,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow1","meandth1","meanpc1","meanpe1","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter4<-function(m1,m2,m3,m4,m5,
                          mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                          nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                          random1,random2,random3,random4,random5,
                          od1,od2,od3,od4,od5,
                          odi1,odi2,odi3,odi4,odi5,
                          alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                          nmax.g,nmin.g,ntol.g,
                          nmax,nmin,tol,power,n,iv1,iv2,iv3){
        # iter4:
        meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        it=0
        par2<-c()
        for (iv4 in 1:fun.seq("nomp",odi4)){
          result.v5<-fun.iter5(m1,m2,m3,m4,m5,
                               mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                               nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                               random1,random2,random3,random4,random5,
                               od1,od2,od3,od4,od5,
                               odi1,odi2,odi3,odi4,odi5,
                               alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                               nmax.g,nmin.g,ntol.g,
                               nmax,nmin,tol,power,n,iv1,iv2,iv3,iv4)
          meanpow1<-result.v5$meanpow1
          meandth1<-result.v5$meandth1
          meanpc1<-result.v5$meanpc1
          meanpe1<-result.v5$meanpe1
          temp1<-result.v5$temp1
          temp3<-result.v5$temp3
          temp4<-result.v5$temp4
          temp5<-result.v5$temp5
          pc2<-result.v5$pc2
          pe2<-result.v5$pe2
          
          if (random4 == 0){
            # if
            while ((abs(meanpow1-power)/power > tol) | (meanpow1 < power)){
              it = it+1
              if (it > itmax){
                stop(" Number of iterations has exceeded the maximum number.")
              }else{
                if (meanpow1 > power){
                  nmax = n
                }else {nmin = n}
                n = (nmin+nmax)/2
                
                # goto iter5
                result.v5l<-fun.iter5(m1,m2,m3,m4,m5,
                                      mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                      nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                      random1,random2,random3,random4,random5,
                                      od1,od2,od3,od4,od5,
                                      odi1,odi2,odi3,odi4,odi5,
                                      alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                                      nmax.g,nmin.g,ntol.g,
                                      nmax,nmin,tol,power,n,iv1,iv2,iv3,iv4)
                
              }# it <= itmax
              meanpow1<-result.v5l$meanpow1
              meandth1<-result.v5l$meandth1
              meanpc1<-result.v5l$meanpc1
              meanpe1<-result.v5l$meanpe1  
              temp1<-result.v5l$temp1
              temp3<-result.v5l$temp3
              temp4<-result.v5l$temp4
              temp5<-result.v5l$temp5
              pc2<-result.v5l$pc2
              pe2<-result.v5l$pe2
            } # while
            
            n2 = ceiling(n); nmin=nmin.g; nmax = nmax.g
            assign(paste0("temp",odi5),NA)
            dth2 = ceiling(meandth1)
            par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                        n2,dth2,periods)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow2 = meanpow2+meanpow1*fun.seq("m",odi4)[iv4,]
          meandth2 = meandth2+meandth1*fun.seq("m",odi4)[iv4,]
          meanpc2 = meanpc2+meanpc1*fun.seq("m",odi4)[iv4,]
          meanpe2 = meanpe2+meanpe1*fun.seq("m",odi4)[iv4,]
        }# /*v4*/ 
        result<-list(meanpow2,meandth2,meanpc2,meanpe2,par2,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow2","meandth2","meanpc2","meanpe2","par2","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter3<-function(m1,m2,m3,m4,m5,
                          mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                          nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                          random1,random2,random3,random4,random5,
                          od1,od2,od3,od4,od5,
                          odi1,odi2,odi3,odi4,odi5,
                          alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                          nmax.g,nmin.g,ntol.g,
                          nmax,nmin,tol,power,n,iv1,iv2){
        # iter3:
        meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
        it=0
        par2<-c()
        for (iv3 in 1:fun.seq("nomp",odi3)){
          result.v4<-fun.iter4(m1,m2,m3,m4,m5,
                               mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                               nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                               random1,random2,random3,random4,random5,
                               od1,od2,od3,od4,od5,
                               odi1,odi2,odi3,odi4,odi5,
                               alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                               nmax.g,nmin.g,ntol.g,
                               nmax,nmin,tol,power,n,iv1,iv2,iv3)
          meanpow2<-result.v4$meanpow2
          meandth2<-result.v4$meandth2
          meanpc2<-result.v4$meanpc2
          meanpe2<-result.v4$meanpe2
          par2<-result.v4$par2
          temp1<-result.v4$temp1
          temp3<-result.v4$temp3
          temp4<-result.v4$temp4
          temp5<-result.v4$temp5
          pc2<-result.v4$pc2
          pe2<-result.v4$pe2           
          
          if ((random3 == 0) & (random4 == 1)){
            # if 
            while ((abs(meanpow2-power)/power > tol) | (meanpow2 < power)) {
              it = it+1
              if (it > itmax){
                stop(" Number of iterations has exceeded the maximum number.")
              }
              if (meanpow2 > power){
                nmax=n
              }else{nmin = n}
              n = (nmin+nmax)/2
              # goto iter4
              result.v4l<-fun.iter4(m1,m2,m3,m4,m5,
                                    mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                    nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                    random1,random2,random3,random4,random5,
                                    od1,od2,od3,od4,od5,
                                    odi1,odi2,odi3,odi4,odi5,
                                    alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                                    nmax.g,nmin.g,ntol.g,
                                    nmax,nmin,tol,power,n,iv1,iv2,iv3)
              meanpow2<-result.v4l$meanpow2
              meandth2<-result.v4l$meandth2
              meanpc2<-result.v4l$meanpc2
              meanpe2<-result.v4l$meanpe2
              par2<-result.v4l$par2
              temp1<-result.v4l$temp1
              temp3<-result.v4l$temp3
              temp4<-result.v4l$temp4
              temp5<-result.v4l$temp5
              pc2<-result.v4l$pc2
              pe2<-result.v4l$pe2               
            }
            nmin=nmin.g;nmax=nmax.g
            assign(paste0("temp",odi5),NA)
            assign(paste0("temp",odi4) ,NA)
            n2 = ceiling(n); dth2 = ceiling(meandth2)
            par = cbind(alpha,meanpow2,temp1,temp3,temp4,temp5,pc2,pe2,
                        n2,dth2,periods)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow3 = meanpow3+meanpow2*fun.seq("m",odi3)[iv3,]
          meandth3 = meandth3+meandth2*fun.seq("m",odi3)[iv3,]
          meanpc3 = meanpc3+meanpc2*fun.seq("m",odi3)[iv3,]
          meanpe3 = meanpe3+meanpe2*fun.seq("m",odi3)[iv3,]
        }# /*v3*/
        result<-list(meanpow3,meandth3,meanpc3,meanpe3,par2,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow3","meandth3","meanpc3","meanpe3","par2","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter2<-function(m1,m2,m3,m4,m5,
                          mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                          nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                          random1,random2,random3,random4,random5,
                          od1,od2,od3,od4,od5,
                          odi1,odi2,odi3,odi4,odi5,
                          alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                          nmax.g,nmin.g,ntol.g,
                          nmax,nmin,tol,power,n,iv1){
        # iter2:
        meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
        it=0
        par2<-c()
        for (iv2 in 1:fun.seq("nomp",odi2)){
          result.v3<-fun.iter3(m1,m2,m3,m4,m5,
                               mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                               nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                               random1,random2,random3,random4,random5,
                               od1,od2,od3,od4,od5,
                               odi1,odi2,odi3,odi4,odi5,
                               alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                               nmax.g,nmin.g,ntol.g,
                               nmax,nmin,tol,power,n,iv1,iv2)
          meanpow3<-result.v3$meanpow3
          meandth3<-result.v3$meandth3
          meanpc3<-result.v3$meanpc3
          meanpe3<-result.v3$meanpe3
          par2<-result.v3$par2
          temp1<-result.v3$temp1
          temp3<-result.v3$temp3
          temp4<-result.v3$temp4
          temp5<-result.v3$temp5
          pc2<-result.v3$pc2
          pe2<-result.v3$pe2
          
          if ((random2 == 0) & (random3 == 1)){
            while ((abs(meanpow3-power)/power > tol) | (meanpow3 < power)){
              it = it+1
              if (it > itmax){
                stop(" Number of iterations has exceeded the maximum number.")
              }
              if (meanpow3 > power) {
                nmax=n
              }else {nmin = n}
              n = (nmin+nmax)/2
              # goto iter3
              result.v3l<-fun.iter3(m1,m2,m3,m4,m5,
                                    mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                    nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                    random1,random2,random3,random4,random5,
                                    od1,od2,od3,od4,od5,
                                    odi1,odi2,odi3,odi4,odi5,
                                    alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                                    nmax.g,nmin.g,ntol.g,
                                    nmax,nmin,tol,power,n,iv1,iv2)
              meanpow3<-result.v3l$meanpow3
              meandth3<-result.v3l$meandth3
              meanpc3<-result.v3l$meanpc3
              meanpe3<-result.v3l$meanpe3
              par2<-result.v3l$par2
              temp1<-result.v3l$temp1
              temp3<-result.v3l$temp3
              temp4<-result.v3l$temp4
              temp5<-result.v3l$temp5
              pc2<-result.v3l$pc2
              pe2<-result.v3l$pe2               
            }
            nmin=nmin.g; nmax=nmax.g
            assign(paste0("temp",odi5),NA)
            assign(paste0("temp",odi4),NA)
            assign(paste0("temp",odi3),NA)
            n2 = ceiling(n); dth2 = ceiling(meandth3)
            par = cbind(alpha,meanpow3,temp1,temp3,temp4,temp5,pc2,pe2,
                        n2,dth2,periods)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow4 = meanpow4+meanpow3*fun.seq("m",odi2)[iv2,]
          meandth4 = meandth4+meandth3*fun.seq("m",odi2)[iv2,]
          meanpc4 = meanpc4+meanpc3*fun.seq("m",odi2)[iv2,]
          meanpe4 = meanpe4+meanpe3*fun.seq("m",odi2)[iv2,]
        }# /*v2*/
        result<-list(meanpow4,meandth4,meanpc4,meanpe4,par2,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow4","meandth4","meanpc4","meanpe4","par2","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter1<-function(m1,m2,m3,m4,m5,
                          mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                          nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                          random1,random2,random3,random4,random5,
                          od1,od2,od3,od4,od5,
                          odi1,odi2,odi3,odi4,odi5,
                          alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                          nmax.g,nmin.g,ntol.g,
                          nmax,nmin,tol,power,n){
        # iter1:
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        it=0
        par2<-c()
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter2:
          result.v2<-fun.iter2(m1,m2,m3,m4,m5,
                               mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                               nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                               random1,random2,random3,random4,random5,
                               od1,od2,od3,od4,od5,
                               odi1,odi2,odi3,odi4,odi5,
                               alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                               nmax.g,nmin.g,ntol.g,
                               nmax,nmin,tol,power,n,iv1)
          meanpow4<-result.v2$meanpow4
          meandth4<-result.v2$meandth4
          meanpc4<-result.v2$meanpc4
          meanpe4<-result.v2$meanpe4
          par2<-result.v2$par2
          temp1<-result.v2$temp1
          temp3<-result.v2$temp3
          temp4<-result.v2$temp4
          temp5<-result.v2$temp5
          pc2<-result.v2$pc2
          pe2<-result.v2$pe2     
          
          if ((random1 == 0) & (random2 == 1)){
            while ((abs(meanpow4-power)/power > tol) | (meanpow4 < power)){
              it = it+1
              if (it > itmax){
                stop(" Number of iterations has exceeded the maximum number.")
              }
              if (meanpow4 > power){
                nmax=n;
              }else{
                nmin = n
              }
              n = (nmin+nmax)/2
              # goto iter2
              result.v2l<-fun.iter2(m1,m2,m3,m4,m5,
                                    mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                    nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                    random1,random2,random3,random4,random5,
                                    od1,od2,od3,od4,od5,
                                    odi1,odi2,odi3,odi4,odi5,
                                    alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                                    nmax.g,nmin.g,ntol.g,
                                    nmax,nmin,tol,power,n,iv1)
              meanpow4<-result.v2l$meanpow4
              meandth4<-result.v2l$meandth4
              meanpc4<-result.v2l$meanpc4
              meanpe4<-result.v2l$meanpe4
              par2<-result.v2l$par2
              temp1<-result.v2l$temp1
              temp3<-result.v2l$temp3
              temp4<-result.v2l$temp4
              temp5<-result.v2l$temp5
              pc2<-result.v2l$pc2
              pe2<-result.v2l$pe2               
            }
            nmin=nmin.g; nmax=nmax.g
            
            assign(paste0("temp",odi5), NA)
            assign(paste0("temp",odi4), NA)
            assign(paste0("temp",odi3), NA)
            assign(paste0("temp",odi2), NA)
            n2 = ceiling(n); dth2 = ceiling(meandth4);
            par = cbind(alpha,meanpow4,temp1,temp3,temp4,temp5,pc2,pe2,
                        n2,dth2,periods)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow5 = meanpow5+meanpow4*fun.seq("m",odi1)[iv1,]
          meandth5 = meandth5+meandth4*fun.seq("m",odi1)[iv1,]
          meanpc5 = meanpc5+meanpc4*fun.seq("m",odi1)[iv1,]
          meanpe5 = meanpe5+meanpe4*fun.seq("m",odi1)[iv1,]
        }# /*v1*/
        result<-list(meanpow5,meandth5,meanpc5,meanpe5,par2,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow5","meandth5","meanpc5","meanpe5","par2","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      #iter1:
      result.v1<-fun.iter1(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                           nmax.g,nmin.g,ntol.g,
                           nmax,nmin,tol,power,n)
      meanpow5<-result.v1$meanpow5
      meandth5<-result.v1$meandth5
      meanpc5<-result.v1$meanpc5
      meanpe5<-result.v1$meanpe5
      par2<-result.v1$par2
      temp1<-result.v1$temp1
      temp3<-result.v1$temp3
      temp4<-result.v1$temp4
      temp5<-result.v1$temp5
      pc2<-result.v1$pc2
      pe2<-result.v1$pe2       
      
      
      if (random1 == 1){
        while ((abs(meanpow5-power)/power > tol) | (meanpow5 < power)){
          it = it+1
          if (it > itmax){
            stop(" Number of iterations has exceeded the maximum number.")
          }
          if (meanpow5 > power) {
            nmax = n
          }else{
            nmin = n
          }
          n = (nmin+nmax)/2
          # goto iter1
          result.v1l<-fun.iter1(m1,m2,m3,m4,m5,
                                mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                random1,random2,random3,random4,random5,
                                od1,od2,od3,od4,od5,
                                odi1,odi2,odi3,odi4,odi5,
                                alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,it,
                                nmax.g,nmin.g,ntol.g,
                                nmax,nmin,tol,power,n)
          meanpow5<-result.v1l$meanpow5
          meandth5<-result.v1l$meandth5
          meanpc5<-result.v1l$meanpc5
          meanpe5<-result.v1l$meanpe5
          par2<-result.v1l$par2
          temp1<-result.v1l$temp1
          temp3<-result.v1l$temp3
          temp4<-result.v1l$temp4
          temp5<-result.v1l$temp5
          pc2<-result.v1l$pc2
          pe2<-result.v1l$pe2
        }
        nmin=nmin.g; nmax=nmax.g
        
        assign(paste0("temp",odi5), NA)
        assign(paste0("temp",odi4), NA)
        assign(paste0("temp",odi3), NA)
        assign(paste0("temp",odi2), NA)
        assign(paste0("temp",odi1), NA)
        n2 = ceiling(n); dth2 = ceiling(meandth5)
        par = cbind(alpha,meanpow5,temp1,temp3,temp4,temp5,pc2,pe2,
                    n2,dth2,periods)
        par2 = rbind(par2,par)
        it = 0
      }
    }# /*power*/
  }
  return(par2)
}

leng<-function(m1,m2,m3,m4,m5,
               mp1,mp2,mp3,mp4,mp5,mp6,mp7,
               nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
               random1,random2,random3,random4,random5,
               od1,od2,od3,od4,od5,
               odi1,odi2,odi3,odi4,odi5,
               alpha,periods,prop,lag,lagdout,recrate,itmax,
               npmax.g,npmin.g,nptol.g,r1,r2,logr,sbdv){
  fun.seq<-function(char,num){
    return(eval(as.symbol(paste0(char,num))))
  }
  if (recrate == 0){
    stop(" Recruitment rate was not specified.")
  }
  
  # Use interval halving method to find the duration time s.t.
  # mean power >= specified one. 
  for (iv7 in 1:nomp7){# /*power*/
    for (iv6 in 1:nomp6){# /* n */
      power = mp7[iv7,]
      n = mp6[iv6,]
      npmin = npmin.g; npmax = npmax.g; tol = nptol.g
      it = 0
      
      fun.iter62<-function(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                           npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5){
        # iter62:
        # rm(con)
        con<-c()
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        nptmp = (npmin+npmax)/2; nptmp2 = ceiling(nptmp)
        if (nptmp2 > periods){
          con= matrix(mp2[1,periods],1,(nptmp2-periods)) #/* pc */
          if (nomp2 > 1){
            for (k in 2:nomp2){
              con2= matrix(mp2[k,periods],1,(nptmp2-periods)) #/* pc */
              con = rbind(con,con2)
            }
          }
          mp2 =cbind(mp2[1:nomp2,1:periods,drop=F],con) #/*expand the vector of pc rates */
        }
        
        if (nptmp2 > periods){
          con= matrix(mp3[1,periods],1,(nptmp2-periods)) #/* dropin */
          if (nomp3 > 1){
            for (k in 2:nomp3){
              con2= matrix(mp3[k,periods],1,(nptmp2-periods)); #/* dropin */
              con = rbind(con,con2)
            }
          }
          mp3 =cbind(mp3[1:nomp3,1:periods,drop=F],con)
        }
        
        if (nptmp2 > periods){
          con= matrix(mp4[1,periods],1,(nptmp2-periods)) #/* dropout */
          if (nomp4 > 1){
            for (k in 2:nomp4){
              con2= matrix(mp4[k,periods],1,(nptmp2-periods)) #/* dropout */
              con = rbind(con,con2)
            }
          }
          mp4 =cbind(mp4[1:nomp4,1:periods,drop=F],con)
        }
        
        if (nptmp2 > periods){
          con= matrix(mp5[1,periods],1,(nptmp2-periods)) #/* loss of follow-up */
          if (nomp5 > 1){
            for (k in 2:nomp5){
              con2= matrix(mp5[k,periods],1,(nptmp2-periods)) #/* loss of follow-up */
              con = rbind(con,con2)
            }
          }
          mp5 =cbind(mp5[1:nomp5,1:periods,drop=F],con)
        }
        
        # Adjust for staggered entry. 
        n_intrvl=sbdv
        nn = ceiling(nptmp*n_intrvl) #/* no. of subperiods */
        ad_cens=matrix(0,1,nn)
        rcrt_sum=0
        nprec = n/recrate
        minrec = nprec
        if (nptmp < nprec) {minrec = nptmp}
        ntmp = minrec*recrate
        nprec = ceiling(minrec*n_intrvl)
        rcrt = matrix(0,1,nn)
        rcrt[,1:nprec] = matrix(recrate,1,nprec)
        for (ii in 1:nn){
          rcrt_sum = rcrt_sum+rcrt[,ii]
          ad_cens[,nn+1-ii] = rcrt[,ii]/rcrt_sum
        }
        if (prop == 1){
          k = mp1[fun.seq("iv",od1),]
        }else{
          # k = mp1
          k = mp1[((fun.seq("iv",od1)-1)*periods+1):(fun.seq("iv",od1)*periods),]
        }
        pc = mp2[fun.seq("iv",od2),]
        dropin = mp3[fun.seq("iv",od3),]
        noncmpl = mp4[fun.seq("iv",od4),]
        loss = mp5[fun.seq("iv",od5),]
        print(ad_cens)
        markov.result4<-markov2(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,nptmp2,nptmp,nn,r1,r2,logr)
        #/* invoke markov chain model for unknown duration */
        ed=markov.result4$ed
        pc2=markov.result4$pc2
        pe2=markov.result4$pe2
        z_alpha=markov.result4$z_alpha
        dropin=markov.result4$dropin
        noncmpl=markov.result4$noncmpl
        loss=markov.result4$loss
        temp1 = mp1[fun.seq("iv",od1),]; temp3=dropin; temp4=noncmpl
        temp5 = loss
        p = (r1*pc2+r2*pe2)/(r1+r2)
        if (((length(logr)==1) & (logr[1] != 0))|(length(logr)>1)){
          dth = ntmp*p
          z_beta = (dth**0.5)*ed-z_alpha
        }else{
          sd1 = ((r1+r2)*p*(1-p)/(r1*r2))**0.5
          sd2 = (pc2*(1-pc2)/r1+pe2*(1-pe2)/r2)**0.5
          z_beta = (((ntmp**0.5)*(abs(pc2-pe2))/((r1+r2)**0.5)-z_alpha*sd1)/sd2)
          dth = ntmp*p
        }
        newpow = pnorm(z_beta)
        result=list(newpow,nptmp,ntmp,minrec,nprec,dth,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("newpow","nptmp","ntmp","minrec","nprec","dth","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter52<-function(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                           npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4){
        # iter52:
        meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
        it=0
        par2<-c()
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        # newpow.temp<-c()
        for (iv5 in 1:fun.seq("nomp",odi5)){
          result62<-fun.iter62(m1,m2,m3,m4,m5,
                               mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                               nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                               random1,random2,random3,random4,random5,
                               od1,od2,od3,od4,od5,
                               odi1,odi2,odi3,odi4,odi5,
                               alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                               npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
          newpow<-result62$newpow
          nptmp<-result62$nptmp
          ntmp<-result62$ntmp
          minrec<-result62$minrec
          nprec<-result62$nprec
          dth<-result62$dth
          temp1<-result62$temp1
          temp3<-result62$temp3
          temp4<-result62$temp4
          temp5<-result62$temp5
          pc2<-result62$pc2
          pe2<-result62$pe2
          # newpow.temp<-rbind(newpow.temp,c(newpow,nptmp,ntmp,minrec,nprec))
          # if (prop == 0){temp1=NA}
          if (random5 == 0){# /* not random */
            while ((abs(newpow-power)/power > tol) | (newpow < power)){
              it = it+1
              if (it > itmax){
                n2 = ceiling(ntmp); dth2 = ceiling(dth)
                par = cbind(alpha,newpow,temp1,temp3,temp4,temp5,pc2,pe2,
                            n2,dth2,nptmp,it)
                colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
                print(par)
                print(it)
                stop(" no. of interations has exceeded the maximum number.")
              }
              if (newpow > power){
                npmax = nptmp
              }else{
                npmin = nptmp
              }
              # goto iter62
              result62l<-fun.iter62(m1,m2,m3,m4,m5,
                                    mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                    nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                    random1,random2,random3,random4,random5,
                                    od1,od2,od3,od4,od5,
                                    odi1,odi2,odi3,odi4,odi5,
                                    alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                    npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
              newpow<-result62l$newpow
              nptmp<-result62l$nptmp
              ntmp<-result62l$ntmp
              minrec<-result62l$minrec
              nprec<-result62l$nprec
              dth<-result62l$dth
              temp1<-result62l$temp1
              temp3<-result62l$temp3
              temp4<-result62l$temp4
              temp5<-result62l$temp5
              pc2<-result62l$pc2
              pe2<-result62l$pe2
              # newpow.temp<-rbind(newpow.temp,c(newpow,nptmp,ntmp,minrec,nprec))
            }
            n2 = ceiling(ntmp); dth2 = ceiling(dth)
            par = cbind(alpha,newpow,temp1,temp3,temp4,temp5,pc2,pe2,
                        n2,dth2,nptmp)
            par2 = rbind(par2,par)
            it = 0
            npmin = npmin.g; npmax = npmax.g
          }
          meanpow1 = meanpow1+newpow*fun.seq("m",odi5)[iv5,]
          meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,]
          meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,]
          meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,]
        }# /*v5*/
        result=list(meanpow1, meandth1, meanpc1, meanpe1, par2,nptmp,ntmp,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow1", "meandth1", "meanpc1", "meanpe1", "par2","nptmp","ntmp","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter42<-function(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                           npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3){
        # iter42:
        meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
        it=0
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        for (iv4 in 1:fun.seq("nomp",odi4)){
          # iter52:
          result.v52<-fun.iter52(m1,m2,m3,m4,m5,
                                 mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                 nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                 random1,random2,random3,random4,random5,
                                 od1,od2,od3,od4,od5,
                                 odi1,odi2,odi3,odi4,odi5,
                                 alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                 npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4)
          meanpow1<-result.v52$meanpow1
          meandth1<-result.v52$meandth1
          meanpc1<-result.v52$meanpc1
          meanpe1<-result.v52$meanpe1
          par2<-result.v52$par2
          nptmp<-result.v52$nptmp
          ntmp<-result.v52$ntmp
          temp1<-result.v52$temp1
          temp3<-result.v52$temp3
          temp4<-result.v52$temp4
          temp5<-result.v52$temp5
          pc2<-result.v52$pc2
          pe2<-result.v52$pe2   
          
          if ((random4 == 0) & (random5 == 1)){
            while ((abs(meanpow1-power)/power > tol) | (meanpow1 < power)){
              it = it+1
              if (it > itmax){
                stop(" no. of iterations has exceeded the maximum number.")
              }
              if (meanpow1 > power){
                npmax = nptmp
              }else{
                npmin = nptmp
              }
              # goto iter52
              result.v52l<-fun.iter52(m1,m2,m3,m4,m5,
                                      mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                      nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                      random1,random2,random3,random4,random5,
                                      od1,od2,od3,od4,od5,
                                      odi1,odi2,odi3,odi4,odi5,
                                      alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                      npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4)
              meanpow1<-result.v52l$meanpow1
              meandth1<-result.v52l$meandth1
              meanpc1<-result.v52l$meanpc1
              meanpe1<-result.v52l$meanpe1
              par2<-result.v52l$par2
              nptmp<-result.v52l$nptmp
              ntmp<-result.v52l$ntmp
              temp1<-result.v52l$temp1
              temp3<-result.v52l$temp3
              temp4<-result.v52l$temp4
              temp5<-result.v52l$temp5
              pc2<-result.v52l$pc2
              pe2<-result.v52l$pe2        
            }
            n2 = ceiling(ntmp); npmin=npmin.g; npmax = npmax.g
            assign(paste0("temp",odi5),NA)
            dth2 = ceiling(meandth1)
            par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,meanpc1,
                        meanpe1,n2,dth2,nptmp)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow2 = meanpow2+meanpow1*fun.seq("m",odi4)[iv4,]
          meandth2 = meandth2+meandth1*fun.seq("m",odi4)[iv4,]
          meanpc2 = meanpc2+meanpc1*fun.seq("m",odi4)[iv4,]
          meanpe2 = meanpe2+meanpe1*fun.seq("m",odi4)[iv4,]
        }# /*v4*/
        result=list(meanpow2, meandth2, meanpc2, meanpe2, par2,nptmp,ntmp,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow2", "meandth2", "meanpc2", "meanpe2", "par2","nptmp","ntmp","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter32<-function(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                           npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2){
        # iter32:
        meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
        it=0
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        for (iv3 in 1:fun.seq("nomp",odi3)){
          # iter42:
          result.v42<-fun.iter42(m1,m2,m3,m4,m5,
                                 mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                 nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                 random1,random2,random3,random4,random5,
                                 od1,od2,od3,od4,od5,
                                 odi1,odi2,odi3,odi4,odi5,
                                 alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                 npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3)
          meanpow2<-result.v42$meanpow2
          meandth2<-result.v42$meandth2
          meanpc2<-result.v42$meanpc2
          meanpe2<-result.v42$meanpe2
          par2<-result.v42$par2
          nptmp<-result.v42$nptmp
          ntmp<-result.v42$ntmp
          temp1<-result.v42$temp1
          temp3<-result.v42$temp3
          temp4<-result.v42$temp4
          temp5<-result.v42$temp5
          pc2<-result.v42$pc2
          pe2<-result.v42$pe2    
          
          if ((random3 == 0) & (random4 == 1)){
            while ((abs(meanpow2-power)/power > tol) | (meanpow2 < power)){
              it = it+1
              if (it > itmax){
                stop(" no. of iterations has exceeded the maximum number.")
              }
              if (meanpow2 > power){
                npmax = nptmp
              }else{
                npmin = nptmp  
              }
              # goto iter42;
              result.v42l<-fun.iter42(m1,m2,m3,m4,m5,
                                      mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                      nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                      random1,random2,random3,random4,random5,
                                      od1,od2,od3,od4,od5,
                                      odi1,odi2,odi3,odi4,odi5,
                                      alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                      npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3)
              meanpow2<-result.v42l$meanpow2
              meandth2<-result.v42l$meandth2
              meanpc2<-result.v42l$meanpc2
              meanpe2<-result.v42l$meanpe2
              par2<-result.v42l$par2
              nptmp<-result.v42l$nptmp
              ntmp<-result.v42l$ntmp
              temp1<-result.v42l$temp1
              temp3<-result.v42l$temp3
              temp4<-result.v42l$temp4
              temp5<-result.v42l$temp5
              pc2<-result.v42l$pc2
              pe2<-result.v42l$pe2        
            }
            n2 = ceiling(ntmp); npmin=npmin.g;npmax=npmax.g
            assign(paste0("temp",odi5), NA)
            assign(paste0("temp",odi4), NA)
            dth2 = ceiling(meandth2)
            if (prop == 0){temp1 =NA}
            par = cbind(alpha,meanpow2,temp1,temp3,temp4,temp5,meanpc2,
                        meanpe2,n2,dth2,nptmp)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow3 = meanpow3+meanpow2*fun.seq("m",odi3)[iv3,]
          meandth3 = meandth3+meandth2*fun.seq("m",odi3)[iv3,]
          meanpc3 = meanpc3+meanpc2*fun.seq("m",odi3)[iv3,]
          meanpe3 = meanpe3+meanpe2*fun.seq("m",odi3)[iv3,]
        }# /*v3*/
        result=list(meanpow3, meandth3, meanpc3, meanpe3, par2,nptmp,ntmp,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow3", "meandth3", "meanpc3", "meanpe3", "par2","nptmp","ntmp","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter22<-function(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                           npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1){
        # iter22:
        meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
        it=0
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        for (iv2 in 1:fun.seq("nomp",odi2)){
          # iter32:
          result.v32<-fun.iter32(m1,m2,m3,m4,m5,
                                 mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                 nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                 random1,random2,random3,random4,random5,
                                 od1,od2,od3,od4,od5,
                                 odi1,odi2,odi3,odi4,odi5,
                                 alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                 npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2)
          meanpow3<-result.v32$meanpow3
          meandth3<-result.v32$meandth3
          meanpc3<-result.v32$meanpc3
          meanpe3<-result.v32$meanpe3
          par2<-result.v32$par2
          nptmp<-result.v32$nptmp
          ntmp<-result.v32$ntmp
          temp1<-result.v32$temp1
          temp3<-result.v32$temp3
          temp4<-result.v32$temp4
          temp5<-result.v32$temp5
          pc2<-result.v32$pc2
          pe2<-result.v32$pe2
          
          if ((random2 == 0) & (random3 == 1)){
            while ((abs(meanpow3-power)/power > tol) | (meanpow3 < power)){
              it = it+1
              if (it > itmax){
                stop(" no. of iterations has exceeded the maximum number.")
              }
              if (meanpow3 > power){
                npmax = nptmp
              }else{
                npmin = nptmp
              }
              # goto iter32
              result.v32l<-fun.iter32(m1,m2,m3,m4,m5,
                                      mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                      nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                      random1,random2,random3,random4,random5,
                                      od1,od2,od3,od4,od5,
                                      odi1,odi2,odi3,odi4,odi5,
                                      alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                      npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2)
              meanpow3<-result.v32l$meanpow3
              meandth3<-result.v32l$meandth3
              meanpc3<-result.v32l$meanpc3
              meanpe3<-result.v32l$meanpe3
              par2<-result.v32l$par2
              nptmp<-result.v32l$nptmp
              ntmp<-result.v32l$ntmp
              temp1<-result.v32l$temp1
              temp3<-result.v32l$temp3
              temp4<-result.v32l$temp4
              temp5<-result.v32l$temp5
              pc2<-result.v32l$pc2
              pe2<-result.v32l$pe2        
            }
            n2 = ceiling(ntmp); npmin=npmin.g; npmax=npmax.g
            assign(paste0("temp",odi5), NA)
            assign(paste0("temp",odi4), NA)
            assign(paste0("temp",odi3), NA)
            dth2 = ceiling(meandth3)
            if (prop == 0) {temp1 = NA}
            par = cbind(alpha,meanpow3,temp1,temp3,temp4,temp5,meanpc3,
                        meanpe3,n2,dth2,nptmp)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow4 = meanpow4+meanpow3*fun.seq("m",odi2)[iv2,]
          meandth4 = meandth4+meandth3*fun.seq("m",odi2)[iv2,]
          meanpc4 = meanpc4+meanpc3*fun.seq("m",odi2)[iv2,]
          meanpe4 = meanpe4+meanpe3*fun.seq("m",odi2)[iv2,]
        }# /*v2*/
        result=list(meanpow4, meandth4, meanpc4, meanpe4, par2,nptmp,ntmp,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow4", "meandth4", "meanpc4", "meanpe4", "par2","nptmp","ntmp","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      fun.iter12<-function(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                           npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n){
        # iter12:
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        it=0
        fun.seq<-function(char,num){
          return(eval(as.symbol(paste0(char,num))))
        }
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter22:
          result.v22<-fun.iter22(m1,m2,m3,m4,m5,
                                 mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                 nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                 random1,random2,random3,random4,random5,
                                 od1,od2,od3,od4,od5,
                                 odi1,odi2,odi3,odi4,odi5,
                                 alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                 npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1)
          meanpow4<-result.v22$meanpow4
          meandth4<-result.v22$meandth4
          meanpc4<-result.v22$meanpc4
          meanpe4<-result.v22$meanpe4
          par2<-result.v22$par2
          nptmp<-result.v22$nptmp
          ntmp<-result.v22$ntmp
          temp1<-result.v22$temp1
          temp3<-result.v22$temp3
          temp4<-result.v22$temp4
          temp5<-result.v22$temp5
          pc2<-result.v22$pc2
          pe2<-result.v22$pe2    
          
          if ((random1 == 0) & (random2 == 1)){
            while ((abs(meanpow4-power)/power > tol) | (meanpow4 < power)){
              it = it+1
              if (it > itmax){
                stop(" no. of interations has exceeded the maximum number.")
              }
              if (meanpow4 > power){
                npmax = nptmp;
              }else{
                npmin = nptmp
              }
              # goto iter22
              result.v22l<-fun.iter22(m1,m2,m3,m4,m5,
                                      mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                      nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                      random1,random2,random3,random4,random5,
                                      od1,od2,od3,od4,od5,
                                      odi1,odi2,odi3,odi4,odi5,
                                      alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                      npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1)
              meanpow4<-result.v22l$meanpow4
              meandth4<-result.v22l$meandth4
              meanpc4<-result.v22l$meanpc4
              meanpe4<-result.v22l$meanpe4
              par2<-result.v22l$par2
              nptmp<-result.v22l$nptmp
              ntmp<-result.v22l$ntmp
              temp1<-result.v22l$temp1
              temp3<-result.v22l$temp3
              temp4<-result.v22l$temp4
              temp5<-result.v22l$temp5
              pc2<-result.v22l$pc2
              pe2<-result.v22l$pe2 
            }
            n2 = ceiling(ntmp); npmin=npmin.g; npmax=npmax.g
            assign(paste0("temp",odi5), NA) 
            assign(paste0("temp",odi4), NA)
            assign(paste0("temp",odi3), NA)
            assign(paste0("temp",odi2), NA)
            dth2 = ceiling(meandth4)
            if (prop == 0){temp1 = NA}
            par = cbind(alpha,meanpow4,temp1,temp3,temp4,temp5,meanpc4,
                        meanpe4,n2,dth2,nptmp)
            par2 = rbind(par2,par)
            it = 0
          }
          meanpow5 = meanpow5+meanpow4*fun.seq("m",odi1)[iv1,]
          meandth5 = meandth5+meandth4*fun.seq("m",odi1)[iv1,]
          meanpc5 = meanpc5+meanpc4*fun.seq("m",odi1)[iv1,]
          meanpe5 = meanpe5+meanpe4*fun.seq("m",odi1)[iv1,]
        }# /*v1*/
        result=list(meanpow5, meandth5, meanpc5, meanpe5, par2,nptmp,ntmp,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("meanpow5", "meandth5", "meanpc5", "meanpe5", "par2","nptmp","ntmp","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      
      # iter12
      result.v12<-fun.iter12(m1,m2,m3,m4,m5,
                             mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                             nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                             random1,random2,random3,random4,random5,
                             od1,od2,od3,od4,od5,
                             odi1,odi2,odi3,odi4,odi5,
                             alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                             npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n)
      meanpow5<-result.v12$meanpow5
      meandth5<-result.v12$meandth5
      meanpc5<-result.v12$meanpc5
      meanpe5<-result.v12$meanpe5
      par2<-result.v12$par2
      nptmp<-result.v12$nptmp
      ntmp<-result.v12$ntmp
      temp1<-result.v12$temp1
      temp3<-result.v12$temp3
      temp4<-result.v12$temp4
      temp5<-result.v12$temp5
      pc2<-result.v12$pc2
      pe2<-result.v12$pe2
      
      if (random1 == 1){
        while ((abs(meanpow5-power)/power > tol) | (meanpow5 < power)){
          it = it+1
          if (it > itmax){
            stop(" no. of interations has exceeded the maximum number.")
          }
          if (meanpow5 > power){
            npmax = nptmp;
          }else{
            npmin = nptmp
          }
          # goto iter12
          result.v12l<-fun.iter12(m1,m2,m3,m4,m5,
                                  mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                  nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                  random1,random2,random3,random4,random5,
                                  od1,od2,od3,od4,od5,
                                  odi1,odi2,odi3,odi4,odi5,
                                  alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                  npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n)
          meanpow5<-result.v12l$meanpow5
          meandth5<-result.v12l$meandth5
          meanpc5<-result.v12l$meanpc5
          meanpe5<-result.v12l$meanpe5
          par2<-result.v12l$par2
          nptmp<-result.v12l$nptmp
          ntmp<-result.v12l$ntmp
          temp1<-result.v12l$temp1
          temp3<-result.v12l$temp3
          temp4<-result.v12l$temp4
          temp5<-result.v12l$temp5
          pc2<-result.v12l$pc2
          pe2<-result.v12l$pe2 
        }
        n2 = ceiling(ntmp); npmin=npmin.g; npmax=npmax.g
        assign(paste0("temp",odi5),NA)
        assign(paste0("temp",odi4),NA)
        assign(paste0("temp",odi3),NA)
        assign(paste0("temp",odi2),NA)
        assign(paste0("temp",odi1),NA)
        dth2 = ceiling(meandth5)
        if (prop == 0) {temp1 = NA}
        par = cbind(alpha,meanpow5,temp1,temp3,temp4,temp5,meanpc5,
                    meanpe5,n2,dth2,nptmp)
        par2 = rbind(par2,par)
        it = 0
      }
    }# /* n */
  }# /* power */
  return(par2)
}

leng2<-function(m1,m2,m3,m4,m5,
                mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                random1,random2,random3,random4,random5,
                od1,od2,od3,od4,od5,
                odi1,odi2,odi3,odi4,odi5,
                alpha,periods,prop,lag,lagdout,recrate,itmax,
                npmax.g,npmin.g,nptol.g,r1,r2,logr,sbdv){
  if (recrate == 0){
    stop(" Recruitment rate was not specified.")
  }
  fun.seq<-function(char,num){
    return(eval(as.symbol(paste0(char,num))))
  }
  # Use interval halving method to find the duration time s.t.
  # mean power >= specified one. 
  for (iv7 in 1:nomp7){# /*power*/
    for (iv6 in 1:nomp6){# /* n */
      power = mp7[iv7,]
      n = mp6[iv6,]
      npmin = npmin.g; npmax = npmax.g; tol = nptol.g
      it = 0
      fun.iter62<-function(m1,m2,m3,m4,m5,
                           mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                           nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                           random1,random2,random3,random4,random5,
                           od1,od2,od3,od4,od5,
                           odi1,odi2,odi3,odi4,odi5,
                           alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                           npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5){
        # iter62:
        con<-c()
        rm(con)
        par2<-c()
        nptmp = (npmin+npmax)/2; nptmp2 = ceiling(nptmp)
        if (nptmp2 > periods){
          con= matrix(mp2[1,periods],1,(nptmp2-periods)) #/* pc */
          if (nomp2 > 1){
            for (k in 2:nomp2){
              con2= matrix(mp2[k,periods],1,(nptmp2-periods)) #/* pc */
              con = rbind(con,con2)
            }
          }
          mp2 =cbind(mp2[1:nomp2,1:periods,drop=F],con) #/*expand the vector of pc rates */
        }
        
        if (nptmp2 > periods){
          con= matrix(mp3[1,periods],1,(nptmp2-periods)) #/* dropin */
          if (nomp3 > 1){
            for (k in 2:nomp3){
              con2= matrix(mp3[k,periods],1,(nptmp2-periods)); #/* dropin */
              con = rbind(con,con2)
            }
          }
          mp3 =cbind(mp3[1:nomp3,1:periods,drop=F],con)
        }
        
        if (nptmp2 > periods){
          con= matrix(mp4[1,periods],1,(nptmp2-periods)) #/* dropout */
          if (nomp4 > 1){
            for (k in 2:nomp4){
              con2= matrix(mp4[k,periods],1,(nptmp2-periods)) #/* dropout */
              con = rbind(con,con2)
            }
          }
          mp4 =cbind(mp4[1:nomp4,1:periods,drop=F],con)
        }
        
        if (nptmp2 > periods){
          con= matrix(mp5[1,periods],1,(nptmp2-periods)) #/* loss of follow-up */
          if (nomp5 > 1){
            for (k in 2:nomp5){
              con2= matrix(mp5[k,periods],1,(nptmp2-periods)) #/* loss of follow-up */
              con = rbind(con,con2)
            }
          }
          mp5 =cbind(mp5[1:nomp5,1:periods,drop=F],con)
        }
        
        # Adjust for staggered entry. 
        n_intrvl=sbdv
        nn = ceiling(nptmp*n_intrvl) #/* no. of subperiods */
        ad_cens=matrix(0,1,nn)
        rcrt_sum=0
        nprec = n/recrate
        minrec = nprec
        if (nptmp < nprec) {minrec = nptmp}
        ntmp = minrec*recrate
        nprec = ceiling(minrec*n_intrvl)
        rcrt = matrix(0,1,nn)
        rcrt[,1:nprec] = matrix(recrate,1,nprec)
        for (ii in 1:nn){
          rcrt_sum = rcrt_sum+rcrt[,ii]
          ad_cens[,nn+1-ii] = rcrt[,ii]/rcrt_sum
        }
        if (prop == 1){
          k = mp1[fun.seq("iv",od1),]
        }else{
          # k = mp1
          k = mp1[((fun.seq("iv",od1)-1)*periods+1):(fun.seq("iv",od1)*periods),]
        }
        pc = mp2[fun.seq("iv",od2),]
        dropin = mp3[fun.seq("iv",od3),]
        noncmpl = mp4[fun.seq("iv",od4),]
        loss = mp5[fun.seq("iv",od5),]
        print(cbind(k))
        print(cbind(pc))
        print(cbind(dropin))
        print(cbind(noncmpl))
        print(cbind(loss))
        print(cbind(n_intrvl))
        # print(ad_cens)
        print(cbind(alpha,prop,lag,lagdout,nptmp2,nptmp,nn))
        markov.result4<-markov2(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,nptmp2,nptmp,nn,r1,r2,logr)
        #/* invoke markov chain model for unknown duration */
        ed=markov.result4$ed
        pc2=markov.result4$pc2
        pe2=markov.result4$pe2
        z_alpha=markov.result4$z_alpha
        dropin=markov.result4$dropin
        noncmpl=markov.result4$noncmpl
        loss=markov.result4$loss
        temp1 = mp1[fun.seq("iv",od1),]; temp3=dropin; temp4=noncmpl
        temp5 = loss
        
        # print(cbind(nn,nptmp,newpow))
        p = (r1*pc2+r2*pe2)/(r1+r2)
        if (((length(logr)==1) & (logr[1] != 0))|(length(logr)>1)){
          dth = ntmp*p
          z_beta = (dth**0.5)*ed-z_alpha
        }else{
          sd1 = ((r1+r2)*p*(1-p)/(r1*r2))**0.5
          sd2 = (pc2*(1-pc2)/r1+pe2*(1-pe2)/r2)**0.5
          z_beta = (((ntmp**0.5)*(abs(pc2-pe2))/((r1+r2)**0.5)-z_alpha*sd1)/sd2)
          dth = ntmp*p
        }
        newpow = pnorm(z_beta)
        print(cbind(temp1,temp3,temp4,temp5))
        print(cbind(pc2,pe2,ed))        
        print(cbind(nn,nptmp,newpow))
        
        result=list(newpow,nptmp,ntmp,minrec,nprec,dth,temp1,temp3,temp4,temp5,pc2,pe2)
        names(result)<-c("newpow","nptmp","ntmp","minrec","nprec","dth","temp1","temp3","temp4","temp5","pc2","pe2")
        return(result)
      }
      if (random5==0){
        # iter12:
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter22:
          meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
          for (iv2 in 1:fun.seq("nomp",odi2)){
            # iter32:
            meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
            for (iv3 in 1:fun.seq("nomp",odi3)){
              # iter42:
              meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
              for (iv4 in 1:fun.seq("nomp",odi4)){
                # iter52:
                meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
                for (iv5 in 1:fun.seq("nomp",odi5)){
                  # iter62:
                  newpow=0
                  while (((abs(newpow-power)/power > tol) | (newpow < power))){
                    result62<-fun.iter62(m1,m2,m3,m4,m5,
                                         mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                         nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                         random1,random2,random3,random4,random5,
                                         od1,od2,od3,od4,od5,
                                         odi1,odi2,odi3,odi4,odi5,
                                         alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                         npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
                    newpow<-result62$newpow
                    nptmp<-result62$nptmp
                    ntmp<-result62$ntmp
                    minrec<-result62$minrec
                    nprec<-result62$nprec
                    dth<-result62$dth
                    temp1<-result62$temp1
                    temp3<-result62$temp3
                    temp4<-result62$temp4
                    temp5<-result62$temp5
                    pc2<-result62$pc2
                    pe2<-result62$pe2
                    
                    if (prop == 0){temp1=NA}
                    
                    it = it+1
                    if (it > itmax){
                      n2 = ceiling(ntmp); dth2 = ceiling(dth)
                      par = cbind(alpha,newpow,temp1,temp3,temp4,temp5,pc2,pe2,
                                  n2,dth2,nptmp,it)
                      colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
                      print(par)
                      stop(" no. of interations has exceeded the maximum number.")
                    }
                    if (newpow > power){
                      npmax = nptmp
                    }else{
                      npmin = nptmp
                    }
                    # goto iter62
                  }
                  n2 = ceiling(ntmp); dth2 = ceiling(dth)
                  par = cbind(alpha,newpow,temp1,temp3,temp4,temp5,pc2,pe2,
                              n2,dth2,nptmp)
                  par2 = rbind(par2,par)
                  
                }# /*v5*/
              }# /*v4*/
            }# /*v3*/
          }# /*v2*/
        }# /*v1*/
      } # random5=0
      if (random4==0 & random5==1){
        # iter12:
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter22:
          meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
          for (iv2 in 1:fun.seq("nomp",odi2)){
            # iter32:
            meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
            for (iv3 in 1:fun.seq("nomp",odi3)){
              # iter42:
              meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
              for (iv4 in 1:fun.seq("nomp",odi4)){
                # iter52:
                meanpow1=0
                while (((abs(meanpow1-power)/power > tol) | (meanpow1 < power))){
                  par2<-c()
                  meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
                  for (iv5 in 1:fun.seq("nomp",odi5)){
                    # iter62:
                    result62<-fun.iter62(m1,m2,m3,m4,m5,
                                         mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                         nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                         random1,random2,random3,random4,random5,
                                         od1,od2,od3,od4,od5,
                                         odi1,odi2,odi3,odi4,odi5,
                                         alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                         npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
                    newpow<-result62$newpow
                    nptmp<-result62$nptmp
                    ntmp<-result62$ntmp
                    minrec<-result62$minrec
                    nprec<-result62$nprec
                    dth<-result62$dth
                    temp1<-result62$temp1
                    temp3<-result62$temp3
                    temp4<-result62$temp4
                    temp5<-result62$temp5
                    pc2<-result62$pc2
                    pe2<-result62$pe2
                    
                    if (prop == 0){temp1=NA}
                    
                    meanpow1 = meanpow1+newpow*fun.seq("m",odi5)[iv5,]
                    meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,]
                    meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,]
                    meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,]
                  }# /*v5*/
                  it = it+1
                  
                  if (it > itmax){
                    n2 = ceiling(ntmp); dth2 = ceiling(meandth1)
                    assign(paste0("temp",odi5),NA)
                    par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                                n2,dth2,nptmp,it)
                    colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
                    print(par)
                    stop(" no. of iterations has exceeded the maximum number.")
                  }
                  if (meanpow1 > power){
                    npmax = nptmp
                  }else{
                    npmin = nptmp
                  }
                  # goto iter52
                }
                n2 = ceiling(ntmp); npmin=npmin; npmax = npmax
                assign(paste0("temp",odi5),NA)
                dth2 = ceiling(meandth1)
                par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,meanpc1,
                            meanpe1,n2,dth2,nptmp)
                par2 = rbind(par2,par)
              }# /*v4*/
            }# /*v3*/
          }# /*v2*/
        }# /*v1*/
      } # random4=0 and random5=1
      if (random3==0 & random4==1){
        # iter12:
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter22:
          meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
          for (iv2 in 1:fun.seq("nomp",odi2)){
            # iter32:
            meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
            for (iv3 in 1:fun.seq("nomp",odi3)){
              # iter42:
              meanpow2=0
              while (((abs(meanpow2-power)/power > tol) | (meanpow2 < power))){
                par2<-c()
              meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
              for (iv4 in 1:fun.seq("nomp",odi4)){
                # iter52:
                  meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
                  for (iv5 in 1:fun.seq("nomp",odi5)){
                    # iter62:
                    result62<-fun.iter62(m1,m2,m3,m4,m5,
                                         mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                         nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                         random1,random2,random3,random4,random5,
                                         od1,od2,od3,od4,od5,
                                         odi1,odi2,odi3,odi4,odi5,
                                         alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                         npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
                    newpow<-result62$newpow
                    nptmp<-result62$nptmp
                    ntmp<-result62$ntmp
                    minrec<-result62$minrec
                    nprec<-result62$nprec
                    dth<-result62$dth
                    temp1<-result62$temp1
                    temp3<-result62$temp3
                    temp4<-result62$temp4
                    temp5<-result62$temp5
                    pc2<-result62$pc2
                    pe2<-result62$pe2
                    
                    if (prop == 0){temp1=NA}
                    
                    meanpow1 = meanpow1+newpow*fun.seq("m",odi5)[iv5,]
                    meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,]
                    meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,]
                    meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,]
                  }# /*v5*/
                meanpow2 = meanpow2+meanpow1*fun.seq("m",odi4)[iv4,]
                meandth2 = meandth2+meandth1*fun.seq("m",odi4)[iv4,]
                meanpc2 = meanpc2+meanpc1*fun.seq("m",odi4)[iv4,]
                meanpe2 = meanpe2+meanpe1*fun.seq("m",odi4)[iv4,]
              }# /*v4*/
                  it = it+1
                  if (it > itmax){
                    n2 = ceiling(ntmp); dth2 = ceiling(meandth2)
                    assign(paste0("temp",odi5), NA)
                    assign(paste0("temp",odi4), NA)
                    par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                                n2,dth2,nptmp,it)
                    colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
                    print(par)
                    stop(" no. of iterations has exceeded the maximum number.")
                  }
                  if (meanpow2 > power){
                    npmax = nptmp
                  }else{
                    npmin = nptmp  
                  }
                  # goto iter42;
                }
                n2 = ceiling(ntmp); npmin=npmin;npmax=npmax
                assign(paste0("temp",odi5), NA)
                assign(paste0("temp",odi4), NA)
                dth2 = ceiling(meandth2)
                if (prop == 0){temp1 =NA}
                par = cbind(alpha,meanpow2,temp1,temp3,temp4,temp5,meanpc2,
                            meanpe2,n2,dth2,nptmp)
                par2 = rbind(par2,par)
            }# /*v3*/
          }# /*v2*/
        }# /*v1*/
      } # random3=0 and random4=1  
      if (random2==0 & random3==1){
        # iter12:
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter22:
          meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
          for (iv2 in 1:fun.seq("nomp",odi2)){
            # iter32:
            meanpow3=0
            while (((abs(meanpow3-power)/power > tol) | (meanpow3 < power))){
              par2<-c()
            meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
            for (iv3 in 1:fun.seq("nomp",odi3)){
              # iter42:

                meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
                for (iv4 in 1:fun.seq("nomp",odi4)){
                  # iter52:
                  meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
                  for (iv5 in 1:fun.seq("nomp",odi5)){
                    # iter62:
                    result62<-fun.iter62(m1,m2,m3,m4,m5,
                                         mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                         nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                         random1,random2,random3,random4,random5,
                                         od1,od2,od3,od4,od5,
                                         odi1,odi2,odi3,odi4,odi5,
                                         alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                         npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
                    newpow<-result62$newpow
                    nptmp<-result62$nptmp
                    ntmp<-result62$ntmp
                    minrec<-result62$minrec
                    nprec<-result62$nprec
                    dth<-result62$dth
                    temp1<-result62$temp1
                    temp3<-result62$temp3
                    temp4<-result62$temp4
                    temp5<-result62$temp5
                    pc2<-result62$pc2
                    pe2<-result62$pe2
                    
                    if (prop == 0){temp1=NA}
                    
                    meanpow1 = meanpow1+newpow*fun.seq("m",odi5)[iv5,]
                    meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,]
                    meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,]
                    meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,]
                  }# /*v5*/
                  meanpow2 = meanpow2+meanpow1*fun.seq("m",odi4)[iv4,]
                  meandth2 = meandth2+meandth1*fun.seq("m",odi4)[iv4,]
                  meanpc2 = meanpc2+meanpc1*fun.seq("m",odi4)[iv4,]
                  meanpe2 = meanpe2+meanpe1*fun.seq("m",odi4)[iv4,]
                }# /*v4*/
              meanpow3 = meanpow3+meanpow2*fun.seq("m",odi3)[iv3,]
              meandth3 = meandth3+meandth2*fun.seq("m",odi3)[iv3,]
              meanpc3 = meanpc3+meanpc2*fun.seq("m",odi3)[iv3,]
              meanpe3 = meanpe3+meanpe2*fun.seq("m",odi3)[iv3,]
            }# /*v3*/
                it = it+1
                if (it > itmax){
                  n2 = ceiling(ntmp); dth2 = ceiling(meandth3)
                  assign(paste0("temp",odi5), NA)
                  assign(paste0("temp",odi4), NA)
                  assign(paste0("temp",odi3), NA)
                  par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                              n2,dth2,nptmp,it)
                  colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
                  print(par)
                  stop(" no. of iterations has exceeded the maximum number.")
                }
                if (meanpow3 > power){
                  npmax = nptmp
                }else{
                  npmin = nptmp
                }
                # goto iter32
              }
              n2 = ceiling(ntmp); npmin=npmin; npmax=npmax
              assign(paste0("temp",odi5), NA)
              assign(paste0("temp",odi4), NA)
              assign(paste0("temp",odi3), NA)
              dth2 = ceiling(meandth3)
              if (prop == 0) {temp1 = NA}
              par = cbind(alpha,meanpow3,temp1,temp3,temp4,temp5,meanpc3,
                          meanpe3,n2,dth2,nptmp)
              par2 = rbind(par2,par)
          }# /*v2*/
        }# /*v1*/
      } # random2=0 and random3=1  
      if (random1==0 & random2==1){
        # iter12:
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter22:
          meanpow4=0
          while (((abs(meanpow4-power)/power > tol) | (meanpow4 < power))){
            par2<-c()
          meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
          for (iv2 in 1:fun.seq("nomp",odi2)){
            # iter32:
              meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
              for (iv3 in 1:fun.seq("nomp",odi3)){
                # iter42:
                
                meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
                for (iv4 in 1:fun.seq("nomp",odi4)){
                  # iter52:
                  meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
                  for (iv5 in 1:fun.seq("nomp",odi5)){
                    # iter62:
                    result62<-fun.iter62(m1,m2,m3,m4,m5,
                                         mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                         nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                         random1,random2,random3,random4,random5,
                                         od1,od2,od3,od4,od5,
                                         odi1,odi2,odi3,odi4,odi5,
                                         alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                         npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
                    newpow<-result62$newpow
                    nptmp<-result62$nptmp
                    ntmp<-result62$ntmp
                    minrec<-result62$minrec
                    nprec<-result62$nprec
                    dth<-result62$dth
                    temp1<-result62$temp1
                    temp3<-result62$temp3
                    temp4<-result62$temp4
                    temp5<-result62$temp5
                    pc2<-result62$pc2
                    pe2<-result62$pe2
                    
                    if (prop == 0){temp1=NA}
                    
                    meanpow1 = meanpow1+newpow*fun.seq("m",odi5)[iv5,]
                    meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,]
                    meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,]
                    meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,]
                  }# /*v5*/
                  meanpow2 = meanpow2+meanpow1*fun.seq("m",odi4)[iv4,]
                  meandth2 = meandth2+meandth1*fun.seq("m",odi4)[iv4,]
                  meanpc2 = meanpc2+meanpc1*fun.seq("m",odi4)[iv4,]
                  meanpe2 = meanpe2+meanpe1*fun.seq("m",odi4)[iv4,]
                }# /*v4*/
                meanpow3 = meanpow3+meanpow2*fun.seq("m",odi3)[iv3,]
                meandth3 = meandth3+meandth2*fun.seq("m",odi3)[iv3,]
                meanpc3 = meanpc3+meanpc2*fun.seq("m",odi3)[iv3,]
                meanpe3 = meanpe3+meanpe2*fun.seq("m",odi3)[iv3,]
              }# /*v3*/
            meanpow4 = meanpow4+meanpow3*fun.seq("m",odi2)[iv2,]
            meandth4 = meandth4+meandth3*fun.seq("m",odi2)[iv2,]
            meanpc4 = meanpc4+meanpc3*fun.seq("m",odi2)[iv2,]
            meanpe4 = meanpe4+meanpe3*fun.seq("m",odi2)[iv2,]
          }# /*v2*/
              it = it+1
              if (it > itmax){
                n2 = ceiling(ntmp); dth2 = ceiling(meandth4)
                assign(paste0("temp",odi5), NA) 
                assign(paste0("temp",odi4), NA)
                assign(paste0("temp",odi3), NA)
                assign(paste0("temp",odi2), NA)
                par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                            n2,dth2,nptmp,it)
                colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
                print(par)
                stop(" no. of interations has exceeded the maximum number.")
              }
              if (meanpow4 > power){
                npmax = nptmp;
              }else{
                npmin = nptmp
              }
              # goto iter22
            }
            n2 = ceiling(ntmp); npmin=npmin; npmax=npmax
            assign(paste0("temp",odi5), NA) 
            assign(paste0("temp",odi4), NA)
            assign(paste0("temp",odi3), NA)
            assign(paste0("temp",odi2), NA)
            dth2 = ceiling(meandth4)
            if (prop == 0) {temp1 = NA}
            par = cbind(alpha,meanpow4,temp1,temp3,temp4,temp5,meanpc4,
                        meanpe4,n2,dth2,nptmp)
            par2 = rbind(par2,par)
        }# /*v1*/
      } # random1=0 and random2=1  
      if (random1==1){
        # iter12:
        meanpow5=0
        while (((abs(meanpow5-power)/power > tol) | (meanpow5 < power))){
          par2<-c()
        meanpow5 = 0; meandth5 = 0; meanpc5 = 0; meanpe5 = 0
        for (iv1 in 1:fun.seq("nomp",odi1)){
          # iter22:
            meanpow4 = 0; meandth4 = 0; meanpc4 = 0; meanpe4 = 0
            for (iv2 in 1:fun.seq("nomp",odi2)){
              # iter32:
              meanpow3 = 0; meandth3 = 0; meanpc3 = 0; meanpe3 = 0
              for (iv3 in 1:fun.seq("nomp",odi3)){
                # iter42:
                
                meanpow2 = 0; meandth2 = 0; meanpc2 = 0; meanpe2 = 0
                for (iv4 in 1:fun.seq("nomp",odi4)){
                  # iter52:
                  meanpow1 = 0; meandth1 = 0; meanpc1 = 0; meanpe1 = 0
                  for (iv5 in 1:fun.seq("nomp",odi5)){
                    # iter62:
                    result62<-fun.iter62(m1,m2,m3,m4,m5,
                                         mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                                         nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                                         random1,random2,random3,random4,random5,
                                         od1,od2,od3,od4,od5,
                                         odi1,odi2,odi3,odi4,odi5,
                                         alpha,periods,prop,lag,lagdout,recrate,itmax,it,
                                         npmax.g,npmin.g,nptol.g,npmax,npmin,tol,power,n,iv1,iv2,iv3,iv4,iv5)
                    newpow<-result62$newpow
                    nptmp<-result62$nptmp
                    ntmp<-result62$ntmp
                    minrec<-result62$minrec
                    nprec<-result62$nprec
                    dth<-result62$dth
                    temp1<-result62$temp1
                    temp3<-result62$temp3
                    temp4<-result62$temp4
                    temp5<-result62$temp5
                    pc2<-result62$pc2
                    pe2<-result62$pe2
                    
                  if (prop == 0){temp1=NA}
                    
                    meanpow1 = meanpow1+newpow*fun.seq("m",odi5)[iv5,]
                    meandth1 = meandth1+dth*fun.seq("m",odi5)[iv5,]
                    meanpc1 = meanpc1+pc2*fun.seq("m",odi5)[iv5,]
                    meanpe1 = meanpe1+pe2*fun.seq("m",odi5)[iv5,]
                  }# /*v5*/
                  meanpow2 = meanpow2+meanpow1*fun.seq("m",odi4)[iv4,]
                  meandth2 = meandth2+meandth1*fun.seq("m",odi4)[iv4,]
                  meanpc2 = meanpc2+meanpc1*fun.seq("m",odi4)[iv4,]
                  meanpe2 = meanpe2+meanpe1*fun.seq("m",odi4)[iv4,]
                }# /*v4*/
                meanpow3 = meanpow3+meanpow2*fun.seq("m",odi3)[iv3,]
                meandth3 = meandth3+meandth2*fun.seq("m",odi3)[iv3,]
                meanpc3 = meanpc3+meanpc2*fun.seq("m",odi3)[iv3,]
                meanpe3 = meanpe3+meanpe2*fun.seq("m",odi3)[iv3,]
              }# /*v3*/
              meanpow4 = meanpow4+meanpow3*fun.seq("m",odi2)[iv2,]
              meandth4 = meandth4+meandth3*fun.seq("m",odi2)[iv2,]
              meanpc4 = meanpc4+meanpc3*fun.seq("m",odi2)[iv2,]
              meanpe4 = meanpe4+meanpe3*fun.seq("m",odi2)[iv2,]
            }# /*v2*/
          meanpow5 = meanpow5+meanpow4*fun.seq("m",odi1)[iv1,]
          meandth5 = meandth5+meandth4*fun.seq("m",odi1)[iv1,]
          meanpc5 = meanpc5+meanpc4*fun.seq("m",odi1)[iv1,]
          meanpe5 = meanpe5+meanpe4*fun.seq("m",odi1)[iv1,]
        }# /*v1*/
            it = it+1
            if (it > itmax){
              n2 = ceiling(ntmp); dth2 = ceiling(meandth5)
              assign(paste0("temp",odi5),NA)
              assign(paste0("temp",odi4),NA)
              assign(paste0("temp",odi3),NA)
              assign(paste0("temp",odi2),NA)
              assign(paste0("temp",odi1),NA)
              par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                          n2,dth2,nptmp,it)
              colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
              print(par)
              stop(" no. of interations has exceeded the maximum number.")
            }
            if (meanpow5 > power){
              npmax = nptmp;
            }else{
              npmin = nptmp
            }
            # goto iter12
          }
        
        n2 = ceiling(ntmp); npmin=npmin; npmax=npmax
        assign(paste0("temp",odi5),NA)
        assign(paste0("temp",odi4),NA)
        assign(paste0("temp",odi3),NA)
        assign(paste0("temp",odi2),NA)
        assign(paste0("temp",odi1),NA)
        dth2 = ceiling(meandth5)
        if (prop == 0) {temp1 = NA}
        par = cbind(alpha,meanpow5,temp1,temp3,temp4,temp5,meanpc5,
          meanpe5,n2,dth2,nptmp)
        par2 = rbind(par2,par)
        it = 0
      } # random1=1
    }# /* n */
  }# /* power */
  return(par2)
}

write<-function(par2){
  # name = c("alpha", "power", "k", "pcno", "din", "dout", "loss", "pc", "pe",
  # "n", "noevent", "duration")
  alpha = round(par2[,1],2)
  power = round(par2[,2],3)
  k = round(par2[,3],4)
  din = round(par2[,4],3)
  dout = round(par2[,5],3)
  loss = round(par2[,6],3)
  pc = round(par2[,7],4)
  pe = round(par2[,8],4)
  n = round(par2[,9],0)
  event = round(par2[,10],0)
  length= round(par2[,11],2)
  
  result=cbind(alpha, power, k, din, dout, loss, pc, pe, n, event, length)
  row.names(result)<-NULL
  print(result)
  # return(result)
  # options ls = 80;
  
  # proc print data=result;
  # format alpha 4.2 power 5.3 k 6.4 din 5.3 dout 5.3 loss 5.3
  # pc 6.4 pe 6.4 n 5.0 event 4.0 length 5.2;
}

# Read in parameters and compute the output
###############################
ptm2<-proc.time()
lakprog<-function(np=0,
                  pc=0,
                  pcratio =0,
                  k=NA,
                  n=NA,
                  power=NA,
                  din=0,
                  dout=0,
                  loss=0,
                  sbdv=20,
                  alpha = .05,
                  simult=1,
                  rectime=0,
                  recratio=0,
                  recrate=0,
                  ratio=1,
                  prop = 1,
                  lag = 0,
                  lagdout = 1,
                  distpc = 0,
                  distk = 0,
                  distdin = 0,
                  distdout = 0,
                  distloss = 0,
                  diratio=0,
                  doratio=0,
                  loratio=0,
                  nmin = 0,
                  nmax = 20000,
                  ntol = .001,
                  npmin = 0,
                  npmax = 200,
                  nptol = .001,
                  itmax = 999,
                  logr = 1,
                  output = 0){
  
  # proc iml worksize=300; 
  # if there is not sufficient working memory,
  # then the worksize has to be adjusted.  
  #   reset noprint;
  # reset log;
  # reset fw=7;
  fun.seq<-function(char,num){
    return(eval(as.symbol(paste0(char,num))))
  }
  
  #########
  # Save the input parameters at first
  var1 = np
  var2 = pc
  var3 = pcratio
  var4 = k
  var5 = logr
  var6 = prop
  var7 = ratio
  var8 = din
  var9 = diratio
  var10 = dout
  var11 = doratio
  var12 = loss
  var13 = loratio
  var14 = alpha
  var15 = power
  var16 = n
  var17 = simult
  var18 = rectime
  var19 = recratio
  var20 = lag
  var21 = lagdout
  var22 = recrate
  var23 = distk
  var24 = distpc
  var25 = distdin
  var26 = distdout
  var27 = distloss
  nmax.g=nmax
  nmin.g=nmin
  ntol.g=ntol
  npmax.g=npmax
  npmin.g=npmin
  nptol.g=nptol
  
  # Check the validity of the parameter values entered
  periods = np
  sbdv = sbdv
  alpha = alpha
  recrate = matrix(recrate,nrow=1)
  ratio = matrix(ratio,nrow=1)
  simult = simult
  itmax = itmax
  lag = lag
  prop = prop
  
  if(output == 1){
    print(paste0("np = ",np))
    print(paste0("pc = ",pc))
    print(paste0("pcratio = ",pcratio))
    
    print(paste0("k = ",k))
    print(paste0("logr = ",logr))
    print(paste0("prop = ",prop))
    
    print(paste0("ratio = ",ratio))
    print(paste0("din = ",din))
    print(paste0("diratio = ",diratio))
    
    print(paste0("dout = ",dout))
    print(paste0("doratio = ",doratio))
    print(paste0("loss = ",loss))
    
    print(paste0("loratio = ",loratio))
    print(paste0("alpha = ",alpha))
    print(paste0("power = ",power))
    
    print(paste0("n = ",n))
    print(paste0("simult = ",simult))
    print(paste0("rectime = ",rectime))
    
    print(paste0("recratio = ",recratio))
    print(paste0("lag = ",lag))
    print(paste0("lagdout = ",lagdout))
    
    print(paste0("recrate = ",recrate))
    print(paste0("distk = ",distk))
    print(paste0("distpc = ",distpc))
    
    print(paste0("distdin = ",distdin))
    print(paste0("distdout = ",distdout))
    print(paste0("distloss = ",distloss))
  }
  
  if (periods == 0)stop(" No. of periods were not specified.")
  
  if (length(ratio) == 1){
    r1=1
    r2=1
  }else{
    r1=ratio[1,1] 
    r2=ratio[1,2]
  }
  
  if ((!is.na(power))&(nchar(power)> 0)){
    power.loc = 999
  }else{
    power.loc = NA
  }
  
  
  if ((!is.na(k)[1])&(nchar(k)[1]> 0)){
    k.loc = 999
  }else{
    k.loc=NA
    stop("k was not specified.")
  } 
  
  
  if ((!is.na(n))&(nchar(n)> 0)){
    ntot = 999
  }else{
    ntot =NA
  }
  
  if (is.na(power.loc) & is.na(ntot)) stop("Neither power nor sample size are unspecified.")
  
  # rm(power.loc, ntot)
  rm(power.loc,ntot)
  
  # If the hazards are nonproportional and there are lagged noncompliance
  # and dropin, then the number of subintervals in each interval is set at 4
  # by default.
  
  if ((sbdv == 20) & (prop!= 1 & lag!= 0)) {sbdv = 4}
  # check k, pc, din, dout, loss are deterministic or from a distribution
  # d1 = distn for k; d2 = distn for pc; d3 = distn for dropins
  # d4 = distn for dropouts; d5 = distn for loss to follow-up
  # start;
  d1 = distk
  d2 = distpc
  d3 = distdin
  d4 = distdout
  d5 = distloss
  # dseq=c(d1,d2,d3,d4,d5)
  
  # If number of masses equals 1, then the parameter is fixed (deterministic).
  # Otherwise sample size or power is calculated to account for heterogeneity
  # max no. of masses 
  m1max = 30; m2max = 30; m3max = 30; m4max = 30; m5max = 30
  # masses
  m1 = matrix(0,m1max,1)
  m2 = matrix(0,m2max,1)
  m3 = matrix(0,m3max,1)
  m4 = matrix(0,m4max,1)
  m5 = matrix(0,m5max,1)
  # mlist<-vector('list',length=5)
  # mlist[[1]]<-m1
  # mlist[[2]]<-m2 
  # mlist[[3]]<-m3
  # mlist[[4]]<-m4
  # mlist[[5]]<-m5
  # no. of masses
  nom1 = 1; nom2 = 1; nom3 = 1; nom4 = 1; nom5 = 1
  # nomseq<-c(nom1,nom2,nom3,nom4,nom5)
  # syml = "("
  # symr = ")"
  # symcom = ","
  
  
  for(ii in 1:5){
    if (nchar(fun.seq("d",ii))[1] > 1) {
      k=1
      assign(paste0("nom",ii),0)
      tmp = matrix(fun.seq("d",ii), nrow=1); tmp = t(tmp)
      notmp = nrow(tmp)
      # fun.seq("m",ii)[k:(k+notmp-1),] = tmp
      assign(paste0("m",ii),rbind(tmp,fun.seq("m",ii)[-(k:(k+notmp-1)),,drop=FALSE]))
      # fun.seq("nom",ii) = fun.seq("nom",ii)+notmp
      assign(paste0("nom",ii), fun.seq("nom",ii)+notmp)
      k = k+notmp  
    }
  } #ii
  
  
  ####################################################
  # Scan k, n, power, pc, dropins, dropouts & loss to
  #   follow-up rates 
  # start;
  # v1=k
  v1 = var4; v2 = pc; v3 = din; v4 = dout
  v5 = loss; v6 = n; v7 = power
  # vseq=c(v1,v2,v3,v4,v5,v6,v7)
  
  u2 = pcratio; u3 = diratio; u4 = doratio
  u5 = loratio
  # useq<-c(0,u2,u3,u4,u5) 
  
  # max no. of mass points 
  mp1max = 30*periods; mp2max = 30; mp3max = 30; mp4max=30
  mp5max=30; mp6max = 30; mp7max = 30
  # mass points 
  mp1 = matrix(0,mp1max,1); mp2=matrix(0,mp2max,periods); mp3=matrix(0,mp3max,periods)
  mp4 = matrix(0,mp4max,periods); mp5=matrix(0,mp5max,periods); mp6=matrix(0,mp6max,1)
  mp7 = matrix(0,mp7max,1)
  
  #   mplist<-vector('list',length=7)
  # mplist[[1]]<-mp1
  # mplist[[2]]<-mp2 
  # mplist[[3]]<-mp3
  # mplist[[4]]<-mp4
  # mplist[[5]]<-mp5
  # mplist[[6]]<-mp6
  # mplist[[7]]<-mp7
  
  # initialize no. of mass points 
  nomp1 = 1; nomp2 = 1; nomp3=1; nomp4=1
  nomp5=1; nomp6=1; nomp7=1
  # nompseq<-c(nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7)
  
  for (ii in 1:3){
    # 1 = k, 2= n, 3 = power
    kk = ii
    do.judge=1
    if (ii > 1) {kk = ii+4}
    if ((kk == 6) & (is.na(v6))){
      # goto (mplp)
      do.judge=0
    } # n is not specified.
    if ((kk == 7) & (is.na(v7)))  {
      # goto (mplp)
      do.judge=0
    } # power is not specified. 
    # do the following block if do.judge=1
    if (do.judge==1){
      k=1
      assign(paste0("nomp",kk),0)
      tmp = matrix(fun.seq("v",kk),nrow=1); tmp = t(tmp);
      notmp = nrow(tmp);
      # fun.seq("mp",kk)[k:k+notmp-1,] = tmp;
      assign(paste0("mp",kk),rbind(tmp,fun.seq("mp",kk)[-(k:(k+notmp-1)),,drop=FALSE]))
      # fun.seq("mp",kk) = fun.seq("nomp",kk)+notmp;
      assign(paste0("nomp",kk), fun.seq("nomp",kk)+notmp)
      k = k+notmp
      
      #####################################  
      
      # %mplp:
    }
  } #ii
  
  if (prop == 0){
    if ((nomp1%%periods)!= 0) {stop(paste0("Number of hazard ratios and number of periods are not equal: nomp1= ",nomp1,", periods= ",periods))
    }else {
      nomp1 = nomp1/periods
    }
  }
  
  if (((prop == 1) & (nomp1 != nom1) & (nom1 > 1)) |
      ((prop == 0) & (nomp1 != nom1) & (nom1 > 1))) {
    stop(paste0("Number of masses and number of mass points for K are not equal: nom1= ",nom1,", nomp1= ",nomp1))
  }
  
  
  # start;
  # 2:pc, 3:din, 4:dout, 5:loss to follow-up 
  # u2=1;u3=2;u4=matrix(3,3,3)    
  # eval(as.symbol(paste0("u",4)))
  
  for (ii in 2:5){
    
    # tempv = {&&v&ii};
    # tempu = {&&u&ii};
    
    # if (((nchar(fun.seq("v",ii))) <= 1)|is.na(fun.seq("v",ii)))
    if ((length(fun.seq("v",ii))) <= 1) {
      assign(paste0("nomp",ii), 1)
      assign(paste0("mp",ii), matrix(0,1,periods))
    }else{
      tmp = mp2max*periods
      tempv = matrix(0,tmp,1)
      
      # if (((nchar(fun.seq("u",ii))) <= 1)|is.na(fun.seq("u",ii)))   
      if ((length(fun.seq("u",ii))) <= 1) {  # ratio is not specified 
        k=1
        assign(paste0("nomp",ii),0)
        tmp = matrix(fun.seq("v",ii),nrow=1); tmp = t(tmp)
        notmp = nrow(tmp)
        tempv[k:(k+notmp-1),] = tmp
        assign(paste0("nomp",ii),fun.seq("nomp",ii)+notmp)
        k = k+notmp
        ############################
        
        if (((fun.seq("nomp",ii))%%periods)!= 0){
          if (ii == 2){
            stop(paste0("The length of pc does not match the number of periods",fun.seq("nomp",ii)))
          }
          if (ii == 3) {
            stop(paste0("The length of din does not match the number of periods", fun.seq("nomp",ii)))
          }
          if (ii == 4){
            stop(paste0("The length of dout does not match the number of periods",fun.seq("nomp",ii)))     
          }
          
          if (ii == 5){
            stop(paste0("The length of loss does not match the number of periods", fun.seq("nomp",ii)))
          }
        }
        assign(paste0("nomp",ii),fun.seq("nomp",ii)/periods)
        assign(paste0("mp",ii), matrix(0,fun.seq("nomp",ii),periods))
        tempmp<-fun.seq("mp",ii)
        for (jj in 1:fun.seq("nomp",ii)){
          tempmp[jj,] = t(tempv[((jj-1)*periods+1):(jj*periods),])
        }
        assign(paste0("mp",ii), tempmp)
        rm(tempmp)
        #/* %length(&&u&ii) = 1 */
      }else{ # /* ratio is specified */
        tempv = matrix(fun.seq("v",ii),nrow=1)
        no1 = ncol(tempv)/2
        tmp = mp2max*periods
        tempu = matrix(0,tmp,1)
        ###########################
        k = 1; notmp2 = 0
        tmp = matrix(fun.seq("u",ii),nrow=1); tmp = t(tmp)
        notmp = nrow(tmp)
        tempu[k:(k+notmp-1),] = tmp
        notmp2 = notmp2+notmp
        k = k+notmp
        ###########################
        
        
        if((notmp2%%periods) != 0) {
          if (ii == 2){
            stop(paste0("The length of pc does not match the number of periods",fun.seq("nomp",ii)))
          }
          if (ii == 3){
            stop(paste0("The length of din does not match the number of periods",fun.seq("nomp",ii)))
          }
          if (ii == 4){
            stop(paste0("The length of dout does not match the number of periods",fun.seq("nomp",ii)))
          }
          if (ii == 5){
            stop(paste0("The length of loss does not match the number of periods",fun.seq("nomp",ii)))
          }
        }
        no2= notmp2/periods
        assign(paste0("nomp",ii), no1*no2)
        assign(paste0("mp",ii), matrix(0,fun.seq("nomp",ii),periods))
        nn = 1
        tempmp<-fun.seq("mp",ii)
        for (kk in 1:no1){
          for (ll in 1:no2){
            sum= 0
            for (mm in 1:tempv[,2*kk]){
              sum= sum+tempu[(ll-1)*periods+mm,]
            }
            rate = -log(1-tempv[,2*(kk-1)+1])/sum
            for (mm in 1:periods){
              # fun.seq("mp",ii)[nn,mm] = 1-exp(-rate*tempu[(ll-1)*periods+mm,])
              tempmp[nn,mm]= 1-exp(-rate*tempu[(ll-1)*periods+mm,])
            }
            nn = nn + 1
          }# ll
        }# kk
        assign(paste0("mp",ii),tempmp)
        rm(tempmp)
      }# /* %length(&&u&ii) ^= 1 */
    }# /* &&v&ii ^= 0 */
    
    if ((fun.seq("nomp",ii) != fun.seq("nom",ii)) & (fun.seq("nom",ii) > 1)){
      if (ii == 2){
        stop(paste0(" number of masses & mass points for pc are not equal: ", fun.seq("nom",ii), ", ", fun.seq("nomp",ii)))
      }
      if (ii == 3){
        stop(paste0(" number of masses & mass points for din are not equal: ", fun.seq("nom",ii), ", ", fun.seq("nomp",ii)))
      }
      if (ii == 4){
        stop(paste0(" number of masses & mass points for dout are not equal: ",fun.seq("nom",ii), ", ", fun.seq("nomp",ii)))
      }
      if (ii == 5){
        stop(paste0(" number of masses & mass points for loss are not equal: ", fun.seq("nom",ii), ", ", fun.seq("nomp",ii)))
      }
    }
    
  } #/* ii */
  
  if (sum(mp1 != 0)==0){
    stop(" No. of periods were not specified.")
  }
  
  
  # start;
  # /* check logr */
  # logr.loc = matrix(0,1,periods)
  # k = 1; nowt = 0
  #############################
  logr = matrix(var5, nrow=1)
  nowt = ncol(logr)
  # logr.loc[,k:(k+notmp-1)] = tmp
  # nowt = nowt+notmp
  # k = k+notmp
  #############################
  
  # logr = logr.loc[,1:nowt]
  if ((nowt != periods) & (nowt !=1)){
    stop("The lenght of the weight vector and the study length are not equal")
  }
  
  # check rectime 
  # rectime.loc = matrix(0,1,periods)
  rectime = matrix(rectime,nrow=1)
  norect = ncol(rectime)
  ##################################
  
  # check recratio
  # recratio.loc = matrix(0,1,periods)
  # k = 1; norecr = 0
  # tmp = matrix(recratio,nrow=1)
  # notmp = ncol(tmp)
  recratio = matrix(recratio,nrow=1)
  norecr = ncol(recratio)
  # k = k+notmp
  ###################################
  
  
  if ((norect != norecr) & ((length(rectime)>1)|(length(rectime)==1 & rectime[1]!=0))){
    stop( "The length of rectime and the legth of recratio are not equal")
  }
  tempu=tempv=no1=no2=sum<-c()
  rm(tempu, tempv, no1, no2, sum)
  rm(m1max, m2max, m3max, m4max, m5max)
  rm(mp1max, mp2max, mp3max, mp4max, mp5max, mp6max, mp7max)
  
  
  # Check the order of the DO loops to be executed below. Events with
  # heterogenity are of lower orders & are executed sooner in the
  # DO loops than those without heterogeneity.
  high = 0; low = 6
  for (ii in 1:5){ #/* 1=k 2=pc 3=din 4=dout 5=loss */
    if ((is.na(fun.seq("d",ii))|(nchar(fun.seq("d",ii)) <= 1))[1]){ #/*no p.d.f. */
      assign(paste0("od",ii),high+1)
      high = fun.seq("od",ii)
      assign(paste0("odi",high),ii)
      assign(paste0("random",high),0)
    }else{# /* has a p.d.f. */
      assign(paste0("od",ii),low-1)
      low = fun.seq("od",ii)
      assign(paste0("odi",low),ii)
      assign(paste0("random",low),1)
    }
  }
  
  
  ######################################################
  
  #Adjust for staggered entry.
  if ((is.na(power)) | (is.na(n))){
    n_intrvl=sbdv  # of subinterval in each period 
    nn=periods*n_intrvl; ad_cens=matrix(0,1,nn); rcrt_sum=0
    if (simult==0){
      
      # rectime = matrix(rectime,nrow=1)
      # j = unlist(gregexpr("C",toupper(recratio)))[1]
      # if (j > 0){
      # recratio = recratio;
      # }else{
      # recratio = matrix(recratio,nrow=1)
      # }
      
      wk = matrix(c(0,as.integer(rectime*n_intrvl+0.5)),nrow=1)
      rcrt=matrix(0,1,nn)
      for(ii in 1:ncol(rectime)){
        rcrt[1,(wk[1,ii]+1):wk[1,ii+1]]=matrix(recratio[1,ii],1,wk[1,ii+1]-wk[1,ii])
      }
      for(ii in 1:nn){
        rcrt_sum=rcrt_sum+rcrt[1,ii];
        ad_cens[1,nn+1-ii]=rcrt[1,ii]/rcrt_sum;
      }
      rm(wk, rcrt)
    } # end staggered entry.
    rm(rcrt_sum)   
  }
  
  ##################################################
  # call macro %power to calculate power given sample size
  if (is.na(power)){
    result<-power.fun(m1,m2,m3,m4,m5,
                      mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                      nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                      random1,random2,random3,random4,random5,
                      od1,od2,od3,od4,od5,
                      odi1,odi2,odi3,odi4,odi5,
                      alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,r1,r2,logr)
  }
  
  # call macro %sample to calculate sample size with given power 
  
  if (is.na(n)){
    result<-sample(m1,m2,m3,m4,m5,
                   mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                   nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                   random1,random2,random3,random4,random5,
                   od1,od2,od3,od4,od5,
                   odi1,odi2,odi3,odi4,odi5,
                   alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,itmax,
                   nmax.g,nmin.g,ntol.g,r1,r2,logr)
  }
  
  # calcualte duration of the study
  if ((!is.na(power)) & (!is.na(n))){
    result<-leng2(m1,m2,m3,m4,m5,
                  mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                  nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                  random1,random2,random3,random4,random5,
                  od1,od2,od3,od4,od5,
                  odi1,odi2,odi3,odi4,odi5,
                  alpha,periods,prop,lag,lagdout,recrate,itmax,
                  npmax.g,npmin.g,nptol.g,r1,r2,logr,sbdv)
  }
  
  
  # result=par2
  # print(result)
  write(result)
  # return(result)
  
}

ptm3<-proc.time()
lakprog(
  np=4,
  pc=c(0.30,4),
  pcratio=c(1,2,2,2),
  k=0.65,
  power=0.8
)

lakprog(
  np=4,
  pc=c(0.30, 4),
  pcratio=c(1, 2, 2, 2),
  k=0.65,
  power=0.8,
  logr=c(1,2,2,2)
)

lakprog(
  np=4,
  pc=c(0.30, 4),
  pcratio=c(1, 2, 2, 2),
  k=0.65,
  power=0.8,
  din=c(.10, 4),
  diratio=c(1, 1, 1, 1),
  dout=c(.10, 4),
  doratio=c(1, 1, 2, 2),
  loss=c(.10, 4),
  loratio=c(1, 1, 1, 1),
  logr=1
)

lakprog(
  np=5,
  pc=c(0.30, 4),
  pcratio=c(1, 2, 2, 2, 2),
  k=0.65,
  power=0.8,
  din=c(.10, 4),
  diratio=c(1, 1, 1, 1, 1),
  dout=c(.10, 4),
  doratio=c(1, 1, 2, 2, 2),
  loss=c(.10, 4),
  loratio=c(1, 1, 1, 1, 1),
  logr=1,
  simult=0,
  rectime=c(1, 2),
  recratio=c(1, 1)
)

# Example 2a
lakprog(
  np=8,
  pc=c(0.30, 8),
  pcratio=rep(1,8),
  prop=0,
  k=c(1, 0.82, rep(0.65,6)),
  power=0.8
)

lakprog(
  np=24,
  pc=c(0.30, 24),
  pcratio=rep(1,24),
  prop=0,
  k=c(1, 0.9417, 0.8834, 0.825, 0.7667, 0.7084, rep(0.65,18)),
  power=0.8
)

lakprog(
  np=24,
  pc=c(0.30, 24),
  pcratio=rep(1,24),
  prop=0,
  k=c(1, 0.9417, 0.8834, 0.825, 0.7667, 0.7084, rep(0.65,18)),
  power=0.8,
  din=c(.10, 24),
  diratio=rep(1,24),
  dout=c(.10, 24),
  doratio=c(rep(1,12), rep(2,12)),
  lag=7,
  lagdout=1
)

lakprog(
  np=4,
  pc=c(0.30,4),
  pcratio=c(1,2,2,2),
  k=0.65,
  power=0.8,
  n=800,
  din=c(0.10,4),
  diratio=c(1,1,1,1),
  dout=c(0.10,4),
  doratio=c(1,1,2,2),
  loss=c(0.10,4),
  loratio=c(1,1,1,1),
  npmin=1,
  npmax=20,
  recrate=200
)

# Example 4
lakprog(
  np=4,
  pc=c(0.20,4),
  pcratio=c(1,1,1,1),
  k=c(0.65,0.35),
  distk=c(0.5,0.5),
  n=600,
  din=c(0.10,4),
  diratio=c(1,1,1,1),
  dout=c(0.10,4),
  doratio=c(1,1,1,1),
  loss=c(0.10,4),
  loratio=c(1,1,1,1)
)

lakprog( 
  np=5,
  pc=rep(0.016,5),
  prop=1,
  k=0.5981,
  power=0.9,
  logr=0,
  din=c(0.09, 0.045, 0.05, 0.055, 0.06),
  dout=c(0.07, rep(0.035,4)),
  loss=seq(.03,.038,.002)
)

lakprog(
  np=6,
  pc=rep(0.016,6),
  prop=1,
  k=0.5981,
  n=4680,
  logr=0,
  din=c(0.09, 0.045, 0.05, 0.055, 0.06, 0.065),
  dout=c(0.07, rep(0.035,5)),
  loss=seq(.03, .04, .002),
  simult=0,
  rectime=c(.25, .50, 2),
  recratio=c(20, 40, 50)
)

lakprog(
  np=24,
  pc=rep(0.004,24),
  prop=0,
  k=c(1, 0.8660, 0.7321, rep(.5981,21)),
  power=0.9,
  logr=0,
  din=c(rep(.0225,4), rep(.0113,4), rep(.0125,4), rep(.0125,4), rep(.015,4), rep(.0163,4)),
  dout= c(rep(.0175,4), rep(0.0088,20)),
  loss= c(rep(.0075,4), rep(.008,4), rep(.0085,4), rep(.009,4), rep(.0095,4), rep(.01,4)),
  simult=0,
  rectime=seq(1, 8, 1),
  recratio=rep(1,8),
  lag=4,
  lagdout=1
)

lakprog(
  np=5,
  pc=c(0.45, 5),
  pcratio=rep(1,5),
  prop=1,
  k=c(1.4313, 1.2452, 1.0078, .9256, .7862, .6574, .5379),
  distk=c(.05, .118, .188, .248, .238, .108, .005),
  n=400,
  logr=1
)

ptm4<-proc.time()
time1<-ptm4-ptm3
as.numeric(time1[3])

# write(par2)
# 
# size<-function(np=0,
#                pc=0,
#                pcratio =0,
#                k=NA,
#                n=NA,
#                power=NA,
#                din=0,
#                dout=0,
#                loss=0,
#                sbdv=20,
#                alpha = .05,
#                simult=1,
#                rectime=0,
#                recratio=0,
#                recrate=0,
#                ratio=1,
#                prop = 1,
#                lag = 0,
#                lagdout = 1,
#                distpc = 0,
#                distk = 0,
#                distdin = 0,
#                distdout = 0,
#                distloss = 0,
#                diratio=0,
#                doratio=0,
#                loratio=0,
#                nmin = 0,
#                nmax = 20000,
#                ntol = .001,
#                npmin = 0,
#                npmax = 200,
#                nptol = .001,
#                itmax = 999,
#                logr = 1,
#                output = 0){
#   par2<-lakprog(np,
#                 pc,
#                 pcratio,
#                 k,
#                 n,
#                 power,
#                 din,
#                 dout,
#                 loss,
#                 sbdv,
#                 alpha,
#                 simult,
#                 rectime,
#                 recratio,
#                 recrate,
#                 ratio,
#                 prop,
#                 lag,
#                 lagdout,
#                 distpc,
#                 distk,
#                 distdin,
#                 distdout,
#                 distloss,
#                 diratio,
#                 doratio,
#                 loratio,
#                 nmin,
#                 nmax,
#                 ntol,
#                 npmin,
#                 npmax,
#                 nptol,
#                 itmax,
#                 logr,
#                 output)
#   write(par2)
# }
# 
# size(
#   np=4,
#   pc=c(0.30,4),
#   pcratio=c(1,2,2,2),
#   k=0.65,
#   power=0.8
# )
# 

# Example k 
## Sample size
lakprog(
  np=8,
  pc=c(0.30, 8),
  pcratio=rep(1,8),
  prop=0,
  k=c(1, 0.82, rep(0.65,6),1, 0.75, rep(0.55,6)),
  distk=c(0.6,0.4),
  power=0.8
)

## Power
lakprog(
  np=8,
  pc=c(0.30, 8),
  pcratio=rep(1,8),
  prop=0,
  k=c(1, 0.82, rep(0.65,6),1, 0.75, rep(0.55,6)),
  distk=c(0.6,0.4),
  # power=0.8
  n=800
)

## Duration time
lakprog(
  np=4,
  pc=c(0.30,4),
  pcratio=c(1,2,2,2),
  # k=0.65,
  # prop=0,
  # k=c(1,0.82,0.65,0.65),
  k=c(0.65,0.6),
  distk=c(0.6,0.4),
  power=0.8,
  n=800,
  din=c(0.10,4),
  diratio=c(1,1,1,1),
  dout=c(0.10,4),
  doratio=c(1,1,2,2),
  loss=c(0.10,4),
  loratio=c(1,1,1,1),
  npmin=1,
  npmax=20,
  recrate=200
)

# Example pc 
## Sample size
lakprog(
  np=8,
  pcratio=rep(1,8),
  pc=c(0.30, 8,0.25,8),
  distpc=c(0.6,0.4),
  prop=0,
  # k=c(1, 0.82, rep(0.65,6)),
  k=c(1, 0.82, rep(0.65,6),1, 0.75, rep(0.55,6)),
  distk=c(0.6,0.4),
  power=0.8
)

## Power
lakprog(
  np=8,
  # pc=c(0.30, 8),
  pc=c(0.30, 8,0.25,8),
  distpc=c(0.6,0.4),
  pcratio=rep(1,8),
  prop=0,
  k=c(1, 0.82, rep(0.65,6),1, 0.75, rep(0.55,6)),
  distk=c(0.6,0.4),
  # power=0.8
  n=800
)

## Duration time
lakprog(
  np=4,
  # pc=c(0.25,4),
  # pc=c(0.3,4),
  pc=c(0.30,4,0.25,4),
  distpc = c(0.6,0.4),
  pcratio=c(1,2,2,2),
  # k=0.65,
  # prop=0,
  # k=c(1,0.82,0.65,0.65),
  k=0.65,
  # k=c(0.5,0.4),
  # distk=c(0.6,0.4),
  power=0.8,
  n=800,
  din=c(0.10,4),
  diratio=c(1,1,1,1),
  dout=c(0.10,4),
  doratio=c(1,1,2,2),
  loss=c(0.10,4),
  loratio=c(1,1,1,1),
  npmin=1,
  npmax=20,
  recrate=200
)

# Example din 
## Sample size
lakprog(
  np=4, 
  pc=c(0.30, 4, 0.25, 4), 
  distpc=c(0.6, 0.4),
  # pc=c(0.30, 4), 
  pcratio=c(1, 2, 2, 2), 
  # k=0.65,  
  k=c(0.65, 0.55), 
  distk=c(0.6, 0.4),
  power=0.8, 
  # ratio=2 1, 
  # din=c(.10, 4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1, 1, 1, 1),
  dout=c(.10, 4),
  doratio=c(1, 1, 2, 2),
  loss=c(.10, 4),
  loratio=c(1, 1, 1, 1),
  logr=1
)

## Power
lakprog(
  np=4, 
  pc=c(0.30, 4, 0.25, 4), 
  distpc=c(0.6, 0.4),
  # pc=c(0.30, 4), 
  pcratio=c(1, 2, 2, 2), 
  # k=0.65,  
  k=c(0.65, 0.55), 
  distk=c(0.6, 0.4),
  # power=0.8, 
  n=800,
  # ratio=2 1, 
  # din=c(.10, 4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1, 1, 1, 1),
  dout=c(.10, 4),
  doratio=c(1, 1, 2, 2),
  loss=c(.10, 4),
  loratio=c(1, 1, 1, 1),
  logr=1
)

## Duration time
lakprog(
  np=4,
  # pc=c(0.25,4),
  # pc=c(0.3,4),
  pc=c(0.30,4,0.25,4),
  distpc = c(0.6,0.4),
  pcratio=c(1,2,2,2),
  # k=0.65,
  # prop=0,
  # k=c(1,0.82,0.65,0.65),
  k=0.65,
  # k=c(0.5,0.4),
  # distk=c(0.6,0.4),
  power=0.8,
  n=800,
  # din=c(0.10,4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1,1,1,1),
  dout=c(0.10,4),
  doratio=c(1,1,2,2),
  loss=c(0.10,4),
  loratio=c(1,1,1,1),
  npmin=1,
  npmax=20,
  recrate=200
)

# Example dout
## Sample size
lakprog(
  np=4, 
  pc=c(0.30, 4, 0.25, 4), 
  distpc=c(0.6, 0.4),
  # pc=c(0.30, 4), 
  pcratio=c(1, 2, 2, 2), 
  # k=0.65,  
  k=c(0.65, 0.55), 
  distk=c(0.6, 0.4),
  power=0.8, 
  # ratio=c(2, 1),
  # din=c(.10, 4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1, 1, 1, 1),
  # dout=c(.10, 4),
  dout=c(.10, 4, 0.12,4),
  distdout=c(0.6, 0.4),
  doratio=c(1, 1, 2, 2),
  loss=c(.10, 4),
  loratio=c(1, 1, 1, 1),
  logr=1
)

## Power
lakprog(
  np=4, 
  pc=c(0.30, 4, 0.25, 4), 
  distpc=c(0.6, 0.4),
  # pc=c(0.30, 4), 
  pcratio=c(1, 2, 2, 2), 
  # k=0.65,  
  k=c(0.65, 0.55), 
  distk=c(0.6, 0.4),
  # power=0.8, 
  n=800,
  # ratio=2 1, 
  # din=c(.10, 4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1, 1, 1, 1),
  # dout=c(.10, 4),
  dout=c(.10, 4,0.12,4),
  distdout=c(0.6, 0.4),
  doratio=c(1, 1, 2, 2),
  loss=c(.10, 4),
  loratio=c(1, 1, 1, 1),
  logr=1
)

## Duration time
lakprog(
  np=4,
  # pc=c(0.25,4),
  # pc=c(0.3,4),
  pc=c(0.30,4,0.25,4),
  distpc = c(0.6,0.4),
  pcratio=c(1,2,2,2),
  # k=0.65,
  # prop=0,
  # k=c(1,0.82,0.65,0.65),
  # k=0.65,
  k=c(0.5,0.4),
  distk=c(0.6,0.4),
  power=0.8,
  n=800,
  # din=c(0.10,4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1,1,1,1),
  # dout=c(0.10,4),
  dout=c(0.10,4,0.12,4),
  distdout=c(0.6, 0.4),
  doratio=c(1,1,2,2),
  loss=c(0.10,4),
  loratio=c(1,1,1,1),
  npmin=1,
  npmax=20,
  recrate=200
)

# Example loss
## Sample size
lakprog(
  np=4, 
  pc=c(0.30, 4, 0.25, 4), 
  distpc=c(0.6, 0.4),
  # pc=c(0.30, 4), 
  pcratio=c(1, 2, 2, 2), 
  # k=0.65,  
  k=c(0.65, 0.55), 
  distk=c(0.6, 0.4),
  power=0.8, 
  # ratio=c(2, 1),
  # din=c(.10, 4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1, 1, 1, 1),
  # dout=c(.10, 4),
  dout=c(.10, 4, 0.12,4),
  distdout=c(0.6, 0.4),
  doratio=c(1, 1, 2, 2),
  # loss=c(.10, 4),
  loss=c(.10, 4, 0.12,4),
  distloss=c(0.6, 0.4),
  loratio=c(1, 1, 1, 1),
  logr=1
)

## Power
lakprog(
  np=4, 
  pc=c(0.30, 4, 0.25, 4), 
  distpc=c(0.6, 0.4),
  # pc=c(0.30, 4), 
  pcratio=c(1, 2, 2, 2), 
  # k=0.65,  
  k=c(0.65, 0.55), 
  distk=c(0.6, 0.4),
  # power=0.8, 
  n=800,
  # ratio=2 1, 
  # din=c(.10, 4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1, 1, 1, 1),
  # dout=c(.10, 4),
  dout=c(.10, 4,0.12,4),
  distdout=c(0.6, 0.4),
  doratio=c(1, 1, 2, 2),
  # loss=c(.10, 4),
  loss=c(.10, 4, 0.12,4),
  distloss=c(0.6, 0.4),
  loratio=c(1, 1, 1, 1),
  logr=1
)

## Duration time
lakprog(
  np=4,
  # pc=c(0.25,4),
  # pc=c(0.3,4),
  pc=c(0.30,4,0.25,4),
  distpc = c(0.6,0.4),
  pcratio=c(1,2,2,2),
  # k=0.65,
  # prop=0,
  # k=c(1,0.82,0.65,0.65),
  # k=0.65,
  k=c(0.5,0.4),
  distk=c(0.6,0.4),
  power=0.8,
  n=800,
  # din=c(0.10,4),
  din=c(0.10, 4, 0.12, 4),
  distdin=c(0.6, 0.4),
  diratio=c(1,1,1,1),
  # dout=c(0.10,4),
  dout=c(0.10,4,0.12,4),
  distdout=c(0.6, 0.4),
  doratio=c(1,1,2,2),
  # loss=c(0.10,4),
  loss=c(.10, 4, 0.12,4),
  distloss=c(0.6, 0.4),
  loratio=c(1,1,1,1),
  npmin=1,
  npmax=20,
  recrate=200
)
