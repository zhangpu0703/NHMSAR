valid_all <-
function(data,simu,root.filename=" ",path=NULL,title="",id=1,alpha=.05,save=FALSE){

dev.new()
N.samples = dim(data)[2]
N.sim = dim(simu)[2]
Bsim = N.sim/N.samples
id = 1 ; 
qqplot(data[,,id],simu[,,id],pch=20,xlab='Observations',ylab='Simulations',cex=.6)
title(title)
abline(a=0,b=1)
q = matrix(0,Bsim,length(data[,,1]))
for (k in 1:Bsim) {
	q[k,] = sort(simu[,((k-1)*N.samples+1):(k*N.samples),id])
}
IC = matrix(0,length(data[,,id]),2)
for (k in 1:length(data[,,id])) {
	IC[k,] = quantile(q[,k],probs=c(alpha/2,1-alpha/2))
}
lines(sort(data[,,id]),IC[,1],lty=2)
lines(sort(data[,,id]),IC[,2],lty=2)
if (save) {
	filename = paste(path,"qqplot-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=4,height=4)
}
dev.new()
C = cor.MSAR(data,simu,lag=15)
plot(0:14,C$C.data,typ="l",ylab="Correlation",xlab="Time (days)",lwd=2)
title(title)
lines(0:14,C$C.sim,col="red")
lines(0:14,C$CI.sim[1,],col="red",lty=3)
lines(0:14,C$CI.sim[2,],col="red",lty=3)
if (save) {
	filename = paste(path,"Cor-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=4,height=4)
}
dev.new()
u = seq(min(data),max(data),by=.3)
gr.d = ENu_graph(data,u)
gr = ENu_graph(simu,u,add=TRUE,col=2,CI = TRUE,N.s.data=dim(data)[2])
abline(v=.5,lty=2)
if (save) {
	filename = paste(path,"ENu-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=4,height=4)
}
dev.new()
MDO = MeanDurOver(data,simu,u)
if (save) {
	filename = paste(path,"MeanDurOver-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=4,height=4)
}
dev.new()
MDU = MeanDurUnder(data,simu,u)
if (save) {
	filename = paste(path,"MeanDurUnder-",root.filename,".eps",sep="")
	dev.copy2eps(file=filename,width=4,height=4)
}
}
