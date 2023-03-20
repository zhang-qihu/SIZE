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
  surv_e<-c()
  surv_c<-c()
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
              surv_e<-rbind(surv_e,markov.result1$surv_e)
              surv_c<-rbind(surv_c,markov.result1$surv_c)

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
  result<-list(par2,surv_e,surv_c)
  names(result)<-c("par2","surv_e","surv_c")
  return(result)
}
