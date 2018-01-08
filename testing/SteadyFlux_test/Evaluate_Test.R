setwd("/Users/lauracondon/Documents/Parflow/EcoSlim/Steady_Flux_Test")
runname="SLIM_flux_test"
fluxrate=0.00125 	#flux rate in m/hr

dx=1
dy=1
dz =0.05
pfdt=1  #parflow timestep

icnp=10 #particles per cell for ic from slimin.txt
pcnp=10 #particles per evap-trans ic from slimin.txt

por = 0.25
denh20=1000

#Calucate Total ET mass out based on PF flux
ETout=fluxrate*10*pfdt*dx*dy*denh20

#Calcualte IC and Precip particle masses
satvol=dx*dy*dz*por
rainvol=fluxrate*pfdt*dx*dy
rainvol/satvol  #Rain and ET flux rate as a fraction of saturated cell volume
ic_mass=satvol/icnp
p_mass=rainvol/pcnp
p_mass/ic_mass   #ratio of mass per precip particle to mass per ic particle


#Read in ET outputs
dir=paste("EcoSLIM", fluxrate, sep="")
fin=paste("./", dir, "/", runname, "_ET_output.txt", sep="")
ET=matrix(scan(fin , skip=1), ncol=7, byrow=T)
colnames(ET)=c("Time", "ET_age", "ET_comp", "ET_mass1", "ET_mass2", "ET_mass3", "ET_np")
nstep=nrow(ET)

#Get a moving average for smoothing
window=501 #doing a centered movign average so it must be odd
half=floor(window/2)
ET_smooth=rep(NA,nstep)
for(i in ceiling(window/2):(nstep-ceiling(window/2))){
	ET_smooth[i]=mean(ET[(i-half):(i+half),5])
}

#plot
#par(mfrow=c(2,1))
plot(ET[,1], ET[,5], type='l', pch=18, xlab="Time Step", ylab="ET Mass out", col='darkgrey')
abline(h=ETout, col=2, lwd=2)
abline(h=mean(ET[,5]), col='green', lwd=2)
lines(ET[,1], ET_smooth, col="blue", lwd=2)
legend('topright', legend=c("PF ET Mass", "Mean Ecoslim ET Mass", "Smoothed ET"), col=c('red', 'green', 'blue'), lwd=c(2,2))

quartz()
plot(ET[,1], ET[,3], type='l', xlab="Time Step", ylab="Composition")

