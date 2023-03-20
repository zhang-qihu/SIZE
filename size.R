#' Sample size, power, and duration time calculation for two-arm
#' clinical trials with survival endpoints.
#' @description Sample size, power or duration of the study are
#' calculated for survival endpoints with time-dependent rates of
#' event, crossover, and loss to follow-up. The calculation allows
#' for the specification of a wide range of complexities commonly
#' occurring in clinical trials including nonproportional hazards,
#' lag in treatment effect, and uncertainties in treatment benefit.
#' It was based on the method of Lakatos (1986, 1988) and previously
#' implemented in SAS.

#' @param np number of periods. The hazard rate is constant within each period.
#' @param pc probability of an event in the control group. \code{pc}
#' may be entered in two ways. First, \code{pc} may be entered as a
#' vector of \code{np} event probabilities. Each entry in the vector
#' is a conditional probability of having an event, given having
#' survived the previous periods. Second, \code{pc} may contain a
#' vector of two values. The first value is the probability of an
#' event and the second value is the number of periods to reach that
#' probability. SIZE allows multiple vectors of \code{pc} and a sample
#'  size or power is calculated for each vector.
#' @param pcratio ratio of the hazard rates for the study endpoint
#' over the \code{np} periods in the control group. If \code{pc} is
#' specified in the second way, \code{pcratio} must be specified.
#' @param k hazard ratio. Values of \code{k} may be entered in two ways:
#' (i) \code{k} = a vector of numbers if prop=0.
#' (ii) \code{k} = scalar if \code{prop}=1.
#' @param n total sample size. If the goal of the analysis is to
#' calculate power, \code{n} must be specified. SIZE allows multiple
#' values of \code{n}. Power is calculated for each \code{n}.
#' @param power If the goal of the analysis is to find sample sizes,
#' power must be specified. SIZE allows multiple values of power.
#' A sample size is calculated for each.
#' @param din the proportion of participants in the control group who
#' switch to the treatment group during the course of follow-up.
#' The default value is 0.
#' @param dout the proportion of participants in the treatment group
#' who switch to the control group during the course of follow-up.
#' The default value is 0.
#' @param loss the proportion of participants lost to follow-up because
#'  of censoring, loss of contact, or competing risks.
#' @param sbdv number of intervals to be divided in each period.
#' Within a given period, the probability of an event is assumed the
#' same across intervals of  equal length and equals
#'  \eqn{1 - (1 - x)^{1/sbdv}}, where \eqn{x} is the probability of
#'  an event in a period. If \code{lag} is not specified, then the
#'  default of \code{sbdv} is 20.  If \code{lag} is specified,
#'  the default is 4. The product of \code{lag} and
#' \code{sbdv} determines the number of intermediate states for the
#' event rate in the treatment group. For example, if \code{lag} = 5 and
#' \code{sbdv} = 4, then there are 5 X 4 = 20 intermediate states
#' between the control and treatment groups. The intermediate states
#' account for the lag in treatment effect.
#' @param alpha type I error for two-tailed alternatives. The default value is 0.05.
#' @param simult 0 or 1. If \code{simult} = 1, then participants
#' enter the study simultaneously. If \code{simult} = 0, participants
#' enter the study according to the accrual pattern specified in
#' \code{rectime} and \code{recratio}.
#' @param rectime  a vector of accrual time. If \code{simult} equals 0,
#' \code{rectime} must be specified.
#' @param recratio a vector of relative accrual rates over the
#' accrual times. The length of $recratio$ and the length of
#' \code{rectime} must be equal.
#' @param recrate the recruitment rate in each period. If both sample
#' size and power are specified and the goal is to solve for the
#' duration time, then \code{recrate} must be specified.
#' @param ratio two numbers referring to the
#' ratio of the sample size in the control group to that in the
#' treatment group. For example, if the allocation ratio is 2:3,
#' then \code{ratio}  = c(2, 3).
#' @param prop 0 or 1. If \code{prop} = 1, then the hazards are
#' proportional and \code{k} remains the same for the duration of the study.
#' If \code{prop} = 0, then the hazard are nonproportional and
#' \code{k} must be a vector of np entries.
#' @param lag lag time, the number of periods to achieve the full
#' treatment effect. It adjusts the event rate for crossovers under
#' the non-proportional hazards model. If \code{lag} is not specified,
#' then participants return to the efficacy level comparable to the
#' level in the opposite arm immediately after crossover.
#' @param lagdout lagged effect of dropout. It equals 1 if crossovers
#'  from the treatment group return to the control group level in the
#'  same fashion as they reach the treatment effect level and equals
#'  0 if the treatment effect vanishes immediately after crossover.
#'  The default value is 1.
#' @param distpc a discrete prior probability distribution for $pc$.
#' Each vector of \code{pc} is considered a mass point and is assigned a
#' probability mass. The number of masses in \code{distpc} must be equal
#' to the number of vectors in \code{pc}. The sum of values of \code{distpc} is 1.
#' @param distk a discrete prior probability distribution for the
#' hazard ratio, \code{k}. Each \code{k} value is considered a mass point and
#' is assigned a probability mass from \code{distk}. The number of masses
#' in \code{distk} must be equal to the number of entries in \code{k}.
#' The sum of values of \code{distk} is 1.
#' @param distdin a discrete prior probability distribution for \code{din}.
#' The sum of values of \code{distdin} is 1.
#' @param distdout a discrete prior probability distribution for \code{dout}.
#' The sum of values of \code{distdout} is 1.
#' @param distloss a discrete prior distribution for \code{loss}.
#' The sum of values of \code{distloss} is 1.
#' @param diratio ratio of the hazard rates for crossovers from the
#' control to the treatment group.
#' @param doratio ratio of the hazard rates for crossovers from the
#' treatment to the control group.
#' @param loratio ratio of the hazard rates for loss.
#' @param nmin,nmax,ntol minimum and maximum of the sample size to be used
#' in the iterative bisection method to search for the sample size
#' such that the calculated power is equal to or grater than the
#' specified power within the relative tolerance, \code{ntol}.
#' Default  = 0, 20000, 0.001.
#' @param npmin,npmax,nptol minimum and maximum of the duration time
#' to be used in the iterative bisection method to search for the
#' duration time such that the calculated power is equal to or
#' greater than the specified power within the relative tolerance,
#' \code{nptol}. Default  = 0, 200, 0.001.
#' @param itmax The maximum number of iterations.
#' @param logr indicator of the test statistic. If \code{logr} = 0,
#' then the binomial formula is used. If \code{logr} = 1, then the log-rank
#' test is used. If \code{logr} = a vector of \code{np} weights, then the
#' weighted log-rank test is used. The default value is 1.
#' @param output 0 or 1. If \code{output} = 1, then all the parameter values
#'  are printed. The default value is 0.
#' @return The output list includes the following elements:
#' \itemize{
#'   \item \code{alpha} - type I error for two-tailed alternatives. The default value is 0.05.
#'   \item \code{power} - the probability of making a correct decision (to reject the null hypothesis) when the null hypothesis is false.
#'   \item \code{k} - hazard ratio.
#'   \item \code{din} - the proportion of participants in the control group who switch to the treatment group during the course of follow-up.
#'   \item \code{dout} - the proportion of participants in the treatment group who switch to the control group during the course of follow-up.
#'   \item \code{loss} - the proportion of participants lost to follow-up because of censoring, loss of contact, or competing risks.
#'   \item \code{pc} - probability of an event in the control group.
#'   \item \code{pe} - probability of an event in the experimental group.
#'   \item \code{n} - total sample size.
#'   \item \code{event} - number of event.
#'   \item \code{duration} - duration time.
#'   }

#' @examples
#' ## Examples with parameter entered with two values
#'
#' ## Sample size
#' size(
#' np=8,
#' pc=c(0.30, 8),
#' pcratio=rep(1,8),
#' k=c(0.65,0.35),
#' distk=c(0.6,0.4),
#' power=0.8)
#'
#'## Power
#'size(
#'np=8,
#'pc=c(0.30, 8),
#'pcratio=rep(1,8),
#'k=c(0.65,0.35),
#'distk=c(0.6,0.4),
#'n=800)
#'
#'## Duration time
#'size(
#'np=4,
#'pc=c(0.30,4),
#'pcratio=c(1,2,2,2),
#'k=0.65,
#'power=0.8,
#'n=800,
#'din=c(0.10,4),
#'diratio=c(1,1,1,1),
#'dout=c(0.10,4),
#'doratio=c(1,1,2,2),
#'loss=c(0.10,4),
#'loratio=c(1,1,1,1),
#'npmin=1,
#'npmax=20,
#'recrate=200)
#'
#' ## Examples with parameter entered with vector
#'
#' ## Sample size
#' size(np=4,
#' pc=c(c(0.35,0.3,0.25,0.2),c(0.3,0.25,0.2,0.15)),
#' distpc=c(0.6, 0.4),
#' k=c(1,0.82,0.65,0.65,1,0.72,0.75,0.45),
#' prop=0,
#' distk=c(0.6, 0.4),
#' power=0.8,
#' din=c(rep(0.005, 4), rep(0.002, 4)),
#' distdin=c(0.6, 0.4),
#' dout=c(rep(0.005, 4), rep(0.002, 4)),
#' distdout=c(0.6, 0.4),
#' loss=c(rep(0.15, 4), rep(0.05, 4)),
#' distloss=c(0.6, 0.4),
#' logr=1)
#'

#'
#' ## Power
#' size(
#' np=4,
#' pc=c(c(0.35,0.3,0.25,0.2),c(0.3,0.25,0.2,0.15)),
#' distpc=c(0.6, 0.4),
#' k=c(1,0.82,0.65,0.65,1,0.72,0.75,0.45),
#' prop=0,
#' distk=c(0.6, 0.4),
#' n=800,
#' din=c(rep(0.005, 4), rep(0.002, 4)),
#' distdin=c(0.6, 0.4),
#' dout=c(rep(0.005, 4), rep(0.002, 4)),
#' distdout=c(0.6, 0.4),
#' loss=c(rep(0.15, 4), rep(0.05, 4)),
#' distloss=c(0.6, 0.4),
#' logr=1)
#'

#'
#' ## Duration time
#' size(
#' np=4,
#' pc=c(0.35,0.3,0.25,0.2),
#' k=c(1,0.82,0.65,0.65),
#' prop=0,
#' power=0.8,
#' n=800,
#' din=rep(0.005, 4),
#' dout=rep(0.005, 4),
#' loss=rep(0.15, 4),
#' npmin=1,
#' npmax=20,
#' recrate=200)
#'


#' @author Qihu Zhang & Joanna H. Shih
#' @importFrom stats pnorm qnorm
#' @references
#' Lakatos, Edward. "Sample size determination in clinical trials with
#' time-dependent rates of losses an noncomplinance."
#' Controlled Clinical Trials 7:189-199,1986.
#'
#' Lakatos, Edward. "Sample sizes based on the log-rank statistic in
#' complex clinical trials." Biometrics 44:229-241, 1988.
#'
#' Shih, Joanna H."Sample Size Calculation for Complex Clinical Trials
#' with Survival Endpoints", Controlled Clinical Trials, 16:395-407,
#' 1995.

#' @export
size<-function(np=0,
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
                  itmax = 99,
                  logr = 1,
                  output = 0){


  fun.seq<-function(char,num){
    return(eval(as.symbol(paste0(char,num))))
  }

  if((length(distk)>1) & (sum(distk)!=1)){
    stop("sum(distk) should be 1.")
  }

  if((length(distpc)>1) & (sum(distpc)!=1)){
    stop("sum(distpc) should be 1.")
  }

  if((length(distdin)>1) & (sum(distdin)!=1)){
    stop("sum(distdin) should be 1.")
  }

  if((length(distdout)>1) & (sum(distdout)!=1)){
    stop("sum(distdout) should be 1.")
  }

  if((length(distloss)>1) & (sum(distloss)!=1)){
    stop("sum(distloss) should be 1.")
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


  rm(power.loc,ntot)

  # If the hazards are nonproportional and there are lagged noncompliance
  # and dropin, then the number of subintervals in each interval is set at 4
  # by default.

  if ((sbdv == 20) & (prop!= 1 & lag!= 0)) {sbdv = 4}
  # check k, pc, din, dout, loss are deterministic or from a distribution
  # d1 = distn for k; d2 = distn for pc; d3 = distn for dropins
  # d4 = distn for dropouts; d5 = distn for loss to follow-up

  d1 = distk
  d2 = distpc
  d3 = distdin
  d4 = distdout
  d5 = distloss


  # If number of masses equals 1, then the parameter is fixed (deterministic).
  # Otherwise sample size or power is calculated to account for heterogeneity
  # max no. of masses
  m1max = 30; m2max = 30; m3max = 30; m4max = 30; m5max = 30

  m1 = matrix(0,m1max,1)
  m2 = matrix(0,m2max,1)
  m3 = matrix(0,m3max,1)
  m4 = matrix(0,m4max,1)
  m5 = matrix(0,m5max,1)

  nom1 = 1; nom2 = 1; nom3 = 1; nom4 = 1; nom5 = 1



  for(ii in 1:5){
    if (nchar(fun.seq("d",ii))[1] > 1) {
      k=1
      assign(paste0("nom",ii),0)
      tmp = matrix(fun.seq("d",ii), nrow=1); tmp = t(tmp)
      notmp = nrow(tmp)

      assign(paste0("m",ii),rbind(tmp,fun.seq("m",ii)[-(k:(k+notmp-1)),,drop=FALSE]))

      assign(paste0("nom",ii), fun.seq("nom",ii)+notmp)
      k = k+notmp
    }
  } #ii


  ####################################################

  # v1=k
  v1 = var4; v2 = pc; v3 = din; v4 = dout
  v5 = loss; v6 = n; v7 = power


  u2 = pcratio; u3 = diratio; u4 = doratio
  u5 = loratio


  # max no. of mass points
  mp1max = 30*periods; mp2max = 30; mp3max = 30; mp4max=30
  mp5max=30; mp6max = 30; mp7max = 30
  # mass points
  mp1 = matrix(0,mp1max,1); mp2=matrix(0,mp2max,periods); mp3=matrix(0,mp3max,periods)
  mp4 = matrix(0,mp4max,periods); mp5=matrix(0,mp5max,periods); mp6=matrix(0,mp6max,1)
  mp7 = matrix(0,mp7max,1)


  # initialize no. of mass points
  nomp1 = 1; nomp2 = 1; nomp3=1; nomp4=1
  nomp5=1; nomp6=1; nomp7=1


  for (ii in 1:3){
    # 1 = k, 2= n, 3 = power
    kk = ii
    do.judge=1
    if (ii > 1) {kk = ii+4}
    if ((kk == 6) & (is.na(v6))){

      do.judge=0
    } # n is not specified.
    if ((kk == 7) & (is.na(v7)))  {

      do.judge=0
    } # power is not specified.
    # do the following block if do.judge=1
    if (do.judge==1){
      k=1
      assign(paste0("nomp",kk),0)
      tmp = matrix(fun.seq("v",kk),nrow=1); tmp = t(tmp);
      notmp = nrow(tmp);

      assign(paste0("mp",kk),rbind(tmp,fun.seq("mp",kk)[-(k:(k+notmp-1)),,drop=FALSE]))

      assign(paste0("nomp",kk), fun.seq("nomp",kk)+notmp)
      k = k+notmp

      #####################################


    }
  } #ii



  if ((prop == 0)&(is.na(n)|is.na(power))){
    if ((nomp1%%periods)!= 0) {stop(paste0("Number of hazard ratios and number of periods are not equal: nomp1= ",nomp1,", periods= ",periods))
    }else {
      nomp1 = nomp1/periods
    }
  }


  if ((prop == 0)&((!is.na(power)) & (!is.na(n)))){
    if ((nomp1%%periods)!= 0) {stop(paste0("Number of hazard ratios and number of periods are not equal: nomp1= ",nomp1,", periods= ",periods))
    }else {
      nomp1 = nomp1/periods
    }
  }

  if (((prop == 1) & (nomp1 != nom1) & (nom1 > 1)) |
      ((prop == 0) & (nomp1 != nom1) & (nom1 > 1))) {
    stop(paste0("Number of masses and number of mass points for K are not equal: nom1= ",nom1,", nomp1= ",nomp1))
  }



  # 2:pc, 3:din, 4:dout, 5:loss to follow-up
  # u2=1;u3=2;u4=matrix(3,3,3)
  # eval(as.symbol(paste0("u",4)))

  for (ii in 2:5){
    if ((length(fun.seq("v",ii))) <= 1) {
      assign(paste0("nomp",ii), 1)
      assign(paste0("mp",ii), matrix(0,1,periods))
    }else{
      tmp = mp2max*periods
      tempv = matrix(0,tmp,1)


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



  # /* check logr */
  #############################
  logr = matrix(var5, nrow=1)
  nowt = ncol(logr)

  #############################


  if ((nowt != periods) & (nowt !=1)){
    stop("The lenght of the weight vector and the study length are not equal")
  }

  # check rectime

  rectime = matrix(rectime,nrow=1)
  norect = ncol(rectime)
  ##################################

  # check recratio

  recratio = matrix(recratio,nrow=1)
  norecr = ncol(recratio)

  ###################################


  if ((norect != norecr) & ((length(rectime)>1)|(length(rectime)==1 & rectime[1]!=0))){
    stop( "The length of rectime and the legth of recratio are not equal")
  }
  tempu=tempv=no1=no2=sum<-c()
  rm(tempu, tempv, no1, no2, sum)
  rm(m1max, m2max, m3max, m4max, m5max)
  rm(mp1max, mp2max, mp3max, mp4max, mp5max, mp6max, mp7max)


  # Check the order of the loops to be executed below. Events with
  # heterogeneity are of lower orders & are executed sooner in the
  # loops than those without heterogeneity.
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
  # calculate power given sample size
  if (is.na(power)){
    result<-power.fun(m1,m2,m3,m4,m5,
                      mp1,mp2,mp3,mp4,mp5,mp6,mp7,
                      nomp1,nomp2,nomp3,nomp4,nomp5,nomp6,nomp7,
                      random1,random2,random3,random4,random5,
                      od1,od2,od3,od4,od5,
                      odi1,odi2,odi3,odi4,odi5,
                      alpha,periods,n_intrvl,ad_cens,prop,lag,lagdout,r1,r2,logr)
  }

  # calculate sample size with given power

  if (is.na(n)){
    result<-sample.fun(m1,m2,m3,m4,m5,
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


  Designed_Quantities0<-result$par2

  Designed_Quantities<-opformat(Designed_Quantities0)
  vect<-paste0("opq",1:length(Designed_Quantities))
  for (i in 1:ncol(Designed_Quantities)) {assign(vect[i],as.vector(Designed_Quantities[,i]))}
  if(sum(is.na(opq3))>0){opq3=var4}
  if(sum(is.na(opq4))>0){opq4=var8}
  if(sum(is.na(opq5))>0){opq5=var10}
  if(sum(is.na(opq6))>0){opq6=var12}
  if(sum(is.na(opq7))>0){opq7=var2}
  result.final<-list(opq1, opq2, opq3, opq4, opq5, opq6, opq7, opq8, opq9, opq10, opq11)

  names(result.final)<-c("alpha", "power", "k", "din", "dout", "loss", "pc", "pe", "n", "event", "duration")
  # print(result.final)
  if (sum(result$gr.itmax)>0){
    warning(paste0("Calculated duration yields ", as.vector(Designed_Quantities0)[2]," power exceeding the desired power (", power,") by more than the specified tolerance level (",nptol,")."))
    cat("\n\n\n")
  }
  if (sum(result$gr.itmax)<0){
    warning(paste0("Calculated duration yields ", as.vector(Designed_Quantities0)[2]," power smaller than the desired power (", power,")."))
    cat("\n\n\n")
  }
  return(result.final)

}
