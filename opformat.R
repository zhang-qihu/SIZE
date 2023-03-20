opformat<-function(par2){

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
  duration= round(par2[,11],2)

  result=cbind(alpha, power, k, din, dout, loss, pc, pe, n, event, duration)
  row.names(result)<-NULL

  return(result)

}
