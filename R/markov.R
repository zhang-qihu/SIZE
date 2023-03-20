markov<-function(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,r1,r2,logr){
  if (length(pc>1)){pc<-matrix(pc,nrow=1)}
  loss<-matrix(loss,nrow=1)
  noncmpl<-matrix(noncmpl,nrow=1)
  dropin<-matrix(dropin,nrow=1)
  k<-matrix(k,ncol=1)
  if (prop == 1 | ((prop == 0 )& (lag == 0))){
    distr_e=matrix(c(0,0,1,0),ncol=1); distr_c=matrix(c(0,0,0,1),ncol=1)

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


    #end if clause for prop hazards or nonlag models;
  }else{# /* nonproportional hazards and lag is specified */
    nactv = lag*n_intrvl
    nstates=2*nactv+2
    distr_e = rbind(0,0,1,matrix(0,nactv-1,1), matrix(0,nactv,1))
    distr_c = rbind(0,0,matrix(0,nactv,1), 1, matrix(0,nactv-1,1))

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

  rho=(event_c*r1+event_e*r2)/(sum(event_c*r1+event_e*r2))
  gamma=abs(phi*theta/(1+phi*theta)-phi/(1+phi))
  eta=phi/((1+phi)**2)

  if(sum(logr==0)==0){
    sig=sqrt(sum((as.vector(logr)**2)*rho*eta))
    sig2 = sum(as.vector(logr)*rho*gamma)
    ed = sig2/sig
  }else{ed=NA}


  pe2=distr_e[2,]; pc2=distr_c[2,];  pbar=(r1*pc2+r2*pe2)/(r1+r2)
  dropin = dstr_c[3,periods]
  noncmpl = dstr_e[4,periods]
  loss = (r1*dstr_c[1,periods]+r2*dstr_e[1,periods])/(r1+r2)
  z_alpha=qnorm(1-alpha/2)

  result=list(z_alpha, pe2, pc2, pbar,ed, dropin,noncmpl,loss,dstr_e[3,],dstr_c[4,])
  names(result)<-c("z_alpha", "pe2", "pc2", "pbar","ed", "dropin","noncmpl", "loss","surv_e","surv_c")
  return(result)
}
