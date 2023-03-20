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
  surv_e<-c()
  surv_c<-c()
  # Use interval halving method to find the duration time s.t.
  # mean power >= specified one.
  for (iv7 in 1:nomp7){# /*power*/
    for (iv6 in 1:nomp6){# /* n */
      power = mp7[iv7,]
      n = mp6[iv6,]
      npmin = npmin.g; npmax = npmax.g; tol = nptol.g
      it = 0
      gr.itmax=0
      temp.result<-c()
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
        surv_e<-c()
        surv_c<-c()
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
          periods.new=npmax.g
          # k = mp1[1:periods.new,]
          k = mp1[((fun.seq("iv",od1)-1)*periods.new+1):(fun.seq("iv",od1)*periods.new),]
          # k = mp1[((fun.seq("iv",od1)-1)*periods+1):(fun.seq("iv",od1)*periods),]
        }
        pc = mp2[fun.seq("iv",od2),]
        dropin = mp3[fun.seq("iv",od3),]
        noncmpl = mp4[fun.seq("iv",od4),]
        loss = mp5[fun.seq("iv",od5),]
        # print(cbind(k))
        # print(cbind(pc))
        # print(cbind(dropin))
        # print(cbind(noncmpl))
        # print(cbind(loss))
        # print(cbind(n_intrvl))
        # print(ad_cens)
        # print(cbind(alpha,prop,lag,lagdout,nptmp2,nptmp,nn))
        markov.result4<-markov2(k,loss,noncmpl,dropin,pc,n_intrvl,ad_cens,alpha,periods,prop,lag,lagdout,nptmp2,nptmp,nn,r1,r2,logr)
        #/* invoke markov chain model for unknown duration */
        ed=markov.result4$ed
        pc2=markov.result4$pc2
        pe2=markov.result4$pe2
        z_alpha=markov.result4$z_alpha
        dropin=markov.result4$dropin
        noncmpl=markov.result4$noncmpl
        loss=markov.result4$loss
        # surv_e=markov.result4$surv_e
        # surv_c=markov.result4$surv_c
        surv_e=rbind(surv_e,markov.result4$surv_e)
        surv_c=rbind(surv_c,markov.result4$surv_c)

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
        # print(cbind(temp1,temp3,temp4,temp5))
        # print(cbind(pc2,pe2,ed))
        # print(cbind(nn,nptmp,newpow))

        result=list(newpow,nptmp,ntmp,minrec,nprec,dth,temp1,temp3,temp4,temp5,pc2,pe2,surv_e,surv_c)
        names(result)<-c("newpow","nptmp","ntmp","minrec","nprec","dth","temp1","temp3","temp4","temp5","pc2","pe2","surv_e","surv_c")
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
                  while (((abs(newpow-power) > tol) | (newpow < power))&(it<itmax)){
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

                    temp.result<-rbind(temp.result,cbind(alpha,newpow,temp1,temp3,temp4,temp5,pc2,pe2,
                                                         ceiling(ntmp),ceiling(dth),nptmp,it))

                    if (it >= itmax){

                      colnames(temp.result) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")
                      if (sum(temp.result[,2]>=power)>0){
                        temp.subset<-temp.result[which(temp.result[,2]>=power),,drop=FALSE]
                        temp.subset<-temp.subset[which.min(temp.subset[,2])[1],]
                        gr.itmax=1
                      }else{
                        temp.subset<-temp.result[which.max(temp.result[,2])[1],]
                        gr.itmax=-1
                      }

                      names(temp.subset) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")


                    }
                    if (newpow > power){
                      npmax = nptmp
                    }else{
                      npmin = nptmp
                    }
                    # goto iter62
                  }
                  surv_e=rbind(surv_e,result62$surv_e)
                  surv_c=rbind(surv_c,result62$surv_c)

                  n2 = ceiling(ntmp); dth2 = ceiling(dth)
                  par = cbind(alpha,newpow,temp1,temp3,temp4,temp5,pc2,pe2,
                              n2,dth2,nptmp)
                  if(gr.itmax!=0){
                    for (i in 1:length(par)){
                      par[i] = temp.subset[i]
                    }
                  }
                  par2<-c()
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
                while (((abs(meanpow1-power) > tol) | (meanpow1 < power))&(it<itmax)){
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
                  temp.result<-rbind(temp.result,cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,meanpc1,
                                                       meanpe1,ceiling(ntmp),ceiling(meandth1),nptmp,it))

                  if (it >= itmax){

                    colnames(temp.result) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

                    if (sum(temp.result[,2]>=power)>0){
                      temp.subset<-temp.result[which(temp.result[,2]>=power),,drop=FALSE]
                      temp.subset<-temp.subset[which.min(temp.subset[,2])[1],]
                      gr.itmax=1
                    }else{
                      temp.subset<-temp.result[which.max(temp.result[,2])[1],]
                      gr.itmax=-1
                    }
                    rm(temp.result)

                    names(temp.subset) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

                  }
                  if (meanpow1 > power){
                    npmax = nptmp
                  }else{
                    npmin = nptmp
                  }
                  # goto iter52
                }
                surv_e=rbind(surv_e,result62$surv_e)
                surv_c=rbind(surv_c,result62$surv_c)

                n2 = ceiling(ntmp); npmin=npmin; npmax = npmax
                assign(paste0("temp",odi5),NA)
                dth2 = ceiling(meandth1)
                par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,meanpc1,
                            meanpe1,n2,dth2,nptmp)
                if(gr.itmax!=0){
                  for (i in 1:length(par)){
                    if(!is.na(par[i])){
                        par[i] = temp.subset[i]
                    }

                  }
                }
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
              while (((abs(meanpow2-power)> tol) | (meanpow2 < power))&(it<itmax)){
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
                temp.result<-rbind(temp.result,cbind(alpha,meanpow2,temp1,temp3,temp4,temp5,meanpc2,
                                                     meanpe2,ceiling(ntmp),ceiling(meandth2),nptmp,it))

                if (it >= itmax){

                  colnames(temp.result) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

                  if (sum(temp.result[,2]>=power)>0){
                    temp.subset<-temp.result[which(temp.result[,2]>=power),,drop=FALSE]
                    temp.subset<-temp.subset[which.min(temp.subset[,2])[1],]
                    gr.itmax=1
                  }else{
                    temp.subset<-temp.result[which.max(temp.result[,2])[1],]
                    gr.itmax=-1
                  }
                  rm(temp.result)

                  names(temp.subset) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")


                }
                if (meanpow2 > power){
                  npmax = nptmp
                }else{
                  npmin = nptmp
                }
                # goto iter42;
              }
              surv_e=rbind(surv_e,result62$surv_e)
              surv_c=rbind(surv_c,result62$surv_c)

              n2 = ceiling(ntmp); npmin=npmin;npmax=npmax
              assign(paste0("temp",odi5), NA)
              assign(paste0("temp",odi4), NA)
              dth2 = ceiling(meandth2)
              if (prop == 0){temp1 =NA}
              par = cbind(alpha,meanpow2,temp1,temp3,temp4,temp5,meanpc2,
                          meanpe2,n2,dth2,nptmp)
              if(gr.itmax!=0){
                for (i in 1:length(par)){
                  if(!is.na(par[i])){
                  par[i] = temp.subset[i]
                  }
                }
              }
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
            while (((abs(meanpow3-power) > tol) | (meanpow3 < power))&(it<itmax)){
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
              temp.result<-rbind(temp.result,cbind(alpha,meanpow3,temp1,temp3,temp4,temp5,meanpc3,
                                                   meanpe3,ceiling(ntmp),ceiling(meandth3),nptmp,it))

              if (it >= itmax){

                colnames(temp.result) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

                if (sum(temp.result[,2]>=power)>0){
                  temp.subset<-temp.result[which(temp.result[,2]>=power),,drop=FALSE]
                  temp.subset<-temp.subset[which.min(temp.subset[,2])[1],]
                  gr.itmax=1
                }else{
                  temp.subset<-temp.result[which.max(temp.result[,2])[1],]
                  gr.itmax=-1
                }
                rm(temp.result)

                names(temp.subset) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")


              }
              if (meanpow3 > power){
                npmax = nptmp
              }else{
                npmin = nptmp
              }
              # goto iter32
            }
            surv_e=rbind(surv_e,result62$surv_e)
            surv_c=rbind(surv_c,result62$surv_c)

            n2 = ceiling(ntmp); npmin=npmin; npmax=npmax
            assign(paste0("temp",odi5), NA)
            assign(paste0("temp",odi4), NA)
            assign(paste0("temp",odi3), NA)
            dth2 = ceiling(meandth3)
            if (prop == 0) {temp1 = NA}
            par = cbind(alpha,meanpow3,temp1,temp3,temp4,temp5,meanpc3,
                        meanpe3,n2,dth2,nptmp)
            if(gr.itmax!=0){
              for (i in 1:length(par)){
                if(!is.na(par[i])){
                par[i] = temp.subset[i]
                }
              }
            }
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
          while (((abs(meanpow4-power) > tol) | (meanpow4 < power))&(it<itmax)){
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
            temp.result<-rbind(temp.result,cbind(alpha,meanpow4,temp1,temp3,temp4,temp5,meanpc4,
                                                 meanpe4,ceiling(ntmp),ceiling(meandth4),nptmp,it))

            if (it >= itmax){

              colnames(temp.result) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

              if (sum(temp.result[,2]>=power)>0){
                temp.subset<-temp.result[which(temp.result[,2]>=power),,drop=FALSE]
                temp.subset<-temp.subset[which.min(temp.subset[,2])[1],]
                gr.itmax=1
              }else{
                temp.subset<-temp.result[which.max(temp.result[,2])[1],]
                gr.itmax=-1
              }
              rm(temp.result)

              names(temp.subset) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

            }
            if (meanpow4 > power){
              npmax = nptmp;
            }else{
              npmin = nptmp
            }
            # goto iter22
          }
          surv_e=rbind(surv_e,result62$surv_e)
          surv_c=rbind(surv_c,result62$surv_c)

          n2 = ceiling(ntmp); npmin=npmin; npmax=npmax
          assign(paste0("temp",odi5), NA)
          assign(paste0("temp",odi4), NA)
          assign(paste0("temp",odi3), NA)
          assign(paste0("temp",odi2), NA)
          dth2 = ceiling(meandth4)
          if (prop == 0) {temp1 = NA}
          par = cbind(alpha,meanpow4,temp1,temp3,temp4,temp5,meanpc4,
                      meanpe4,n2,dth2,nptmp)
          if(gr.itmax!=0){
            for (i in 1:length(par)){
              if(!is.na(par[i])){
              par[i] = temp.subset[i]
              }
            }
          }
          par2 = rbind(par2,par)
        }# /*v1*/
      } # random1=0 and random2=1
      if (random1==1){
        # iter12:
        meanpow5=0
        while (((abs(meanpow5-power) > tol) | (meanpow5 < power))&(it<itmax)){
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
          temp.result<-rbind(temp.result,cbind(alpha,meanpow5,temp1,temp3,temp4,temp5,meanpc5,
                                               meanpe5,ceiling(ntmp),ceiling(meandth5),nptmp,it))

          if (it >= itmax){

            colnames(temp.result) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

            if (sum(temp.result[,2]>=power)>0){
              temp.subset<-temp.result[which(temp.result[,2]>=power),,drop=FALSE]
              temp.subset<-temp.subset[which.min(temp.subset[,2])[1],]
              gr.itmax=1
            }else{
              temp.subset<-temp.result[which.max(temp.result[,2])[1],]
              gr.itmax=-1
            }

            names(temp.subset) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration","it")

          }

          if (meanpow5 > power){
            npmax = nptmp;
          }else{
            npmin = nptmp
          }
          # goto iter12
        }
        surv_e=rbind(surv_e,result62$surv_e)
        surv_c=rbind(surv_c,result62$surv_c)

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
        if(gr.itmax!=0){
          for (i in 1:length(par)){
            if(!is.na(par[i])){
            par[i] = temp.subset[i]
            }
          }
        }
        par2 = rbind(par2,par)
        it = 0
      } # random1=1
    }# /* n */
  }# /* power */
  result<-list(par2,surv_e,surv_c,gr.itmax)
  names(result)<-c("par2","surv_e","surv_c","gr.itmax")
  return(result)
}
