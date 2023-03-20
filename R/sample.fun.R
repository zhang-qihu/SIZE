sample.fun<-function(m1,m2,m3,m4,m5,
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
  surv_e<-c()
  surv_c<-c()
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

                surv_e<-rbind(surv_e,markov.result2$surv_e)
                surv_c<-rbind(surv_c,markov.result2$surv_c)

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
        surv_e<-c()
        surv_c<-c()
        for (iv5 in 1:fun.seq("nomp",odi5)){
          if (prop == 1) {
            k = mp1[fun.seq("iv",od1),]
          }else{

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
          surv_e<-rbind(surv_e,markov.result3$surv_e)
          surv_c<-rbind(surv_c,markov.result3$surv_c)

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
        result<-list(meanpow1,meandth1,meanpc1,meanpe1,temp1,temp3,temp4,temp5,pc2,pe2,surv_e,surv_c)
        names(result)<-c("meanpow1","meandth1","meanpc1","meanpe1","temp1","temp3","temp4","temp5","pc2","pe2","surv_e","surv_c")
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
        surv_e<-c()
        surv_c<-c()
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
            while ((abs(meanpow1-power)/power > tol) | (meanpow1 < power)){
              it = it+1
              if (it > itmax){
                n2 = ceiling(n)
                assign(paste0("temp",odi5),NA)
                dth2 = ceiling(meandth1)
                par = cbind(alpha,meanpow1,temp1,temp3,temp4,temp5,pc2,pe2,
                            n2,dth2,periods)
                colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration")
                print(par)
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

              }
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
            surv_e<-rbind(surv_e,result.v5l$surv_e)
            surv_c<-rbind(surv_c,result.v5l$surv_c)

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
        result<-list(meanpow2,meandth2,meanpc2,meanpe2,par2,temp1,temp3,temp4,temp5,pc2,pe2,surv_e,surv_c)
        names(result)<-c("meanpow2","meandth2","meanpc2","meanpe2","par2","temp1","temp3","temp4","temp5","pc2","pe2","surv_e","surv_c")
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
        surv_e<-c()
        surv_c<-c()
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
            while ((abs(meanpow2-power)/power > tol) | (meanpow2 < power)) {
              it = it+1
              if (it > itmax){
                assign(paste0("temp",odi5),NA)
                assign(paste0("temp",odi4) ,NA)
                n2 = ceiling(n); dth2 = ceiling(meandth2)
                par = cbind(alpha,meanpow2,temp1,temp3,temp4,temp5,pc2,pe2,
                            n2,dth2,periods)
                colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration")
                print(par)
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
            surv_e<-rbind(surv_e,result.v4l$surv_e)
            surv_c<-rbind(surv_e,result.v4l$surv_c)

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
        result<-list(meanpow3,meandth3,meanpc3,meanpe3,par2,temp1,temp3,temp4,temp5,pc2,pe2,surv_e,surv_c)
        names(result)<-c("meanpow3","meandth3","meanpc3","meanpe3","par2","temp1","temp3","temp4","temp5","pc2","pe2","surv_e","surv_c")
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
        surv_e<-c()
        surv_c<-c()
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
                assign(paste0("temp",odi5),NA)
                assign(paste0("temp",odi4),NA)
                assign(paste0("temp",odi3),NA)
                n2 = ceiling(n); dth2 = ceiling(meandth3)
                par = cbind(alpha,meanpow3,temp1,temp3,temp4,temp5,pc2,pe2,
                            n2,dth2,periods)
                colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration")
                print(par)
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
            surv_e<-rbind(surv_e,result.v3l$surv_e)
            surv_c<-rbind(surv_c,result.v3l$surv_c)

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
        result<-list(meanpow4,meandth4,meanpc4,meanpe4,par2,temp1,temp3,temp4,temp5,pc2,pe2,surv_e,surv_c)
        names(result)<-c("meanpow4","meandth4","meanpc4","meanpe4","par2","temp1","temp3","temp4","temp5","pc2","pe2","surv_e","surv_c")
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
        surv_e<-c()
        surv_c<-c()
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
                assign(paste0("temp",odi5), NA)
                assign(paste0("temp",odi4), NA)
                assign(paste0("temp",odi3), NA)
                assign(paste0("temp",odi2), NA)
                n2 = ceiling(n); dth2 = ceiling(meandth4);
                par = cbind(alpha,meanpow4,temp1,temp3,temp4,temp5,pc2,pe2,
                            n2,dth2,periods)
                colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration")
                print(par)
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
            surv_e<-rbind(surv_e,result.v2l$surv_e)
            surv_c<-rbind(surv_c,result.v2l$surv_c)
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
        result<-list(meanpow5,meandth5,meanpc5,meanpe5,par2,temp1,temp3,temp4,temp5,pc2,pe2,surv_e,surv_c)
        names(result)<-c("meanpow5","meandth5","meanpc5","meanpe5","par2","temp1","temp3","temp4","temp5","pc2","pe2","surv_e","surv_c")
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
            assign(paste0("temp",odi5), NA)
            assign(paste0("temp",odi4), NA)
            assign(paste0("temp",odi3), NA)
            assign(paste0("temp",odi2), NA)
            assign(paste0("temp",odi1), NA)
            n2 = ceiling(n); dth2 = ceiling(meandth5)
            par = cbind(alpha,meanpow5,temp1,temp3,temp4,temp5,pc2,pe2,
                        n2,dth2,periods)
            colnames(par) = c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration")
            print(par)
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
        surv_e<-rbind(surv_e,result.v1l$surv_e)
        surv_c<-rbind(surv_c,result.v1l$surv_c)
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

  result<-list(par2,surv_e,surv_c)
  names(result)<-c("par2","surv_e","surv_c")
  return(result)
}
