MeanDurOver <-
function(data,data.sim,u,alpha=.05,col="red",plot=TRUE){
N.samples = dim(data)[2]
Bsim =  dim(data.sim)[2]/N.samples
tps.data = NULL
tps.hh = NULL
F = NULL
Fhh = NULL
for (k in 1:length(u))
{
	F[k] = sum(data<u[k],na.rm=TRUE)
	tps.data[k] = mean(tps_sejour_sup(data,s=u[k]))
	Fhh[k] = sum(data.sim<u[k])
	tps.hh[k] = mean(tps_sejour_sup(data.sim,s=u[k]) )
}
F = F/length(data)
Fhh = Fhh/length(data.sim)

q.tps.nh = matrix(0,Bsim,length(u))
q.tps.hh = matrix(0,Bsim,length(u))
for (ks in 1:Bsim) {
	xx.hh = data.sim[,((ks-1)*N.samples+1):(ks*N.samples),]
	for (k in 1:length(u)){
		q.tps.hh[ks,k] = mean(tps_sejour_sup(xx.hh,u[k]))
	}
}
tps.hh= apply(q.tps.hh,2,mean)

ICinf.hh = NULL
ICsup.hh = NULL
for (k in 1:length(u)) {
	ICinf.hh[k] = quantile(q.tps.hh[,k],alpha/2)
	ICsup.hh[k] = quantile(q.tps.hh[,k],1-alpha/2)
}
CI = matrix(0,length(u),2)
CI[,1] = ICinf.hh
CI[,2] = ICsup.hh

if (plot) {
	plot(F,tps.data,typ="l",ylim=range(c(tps.data,tps.hh)),log="y",lwd=3,xlab = "P(Y<u)",ylab="Mean sojourn duration over u (log scale)")
	lines(Fhh,tps.hh,col=col)
	lines(Fhh,ICinf.hh,col=col,lty=2)
	lines(Fhh,ICsup.hh,col=col,lty=2)

}
return(list(F=F,mdo.data = tps.data,F.sim=Fhh,mdo.sim = q.tps.hh,CI=CI,mdo.sim.all=q.tps.hh))
}
