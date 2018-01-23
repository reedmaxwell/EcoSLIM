rm(list=ls())
setwd("/Users/lauracondon/Documents/Parflow/EcoSlim/SteadyFlux_test_1_23")
#setwd("/Users/lauracondon/Documents/ParFlow/EcoSlim_Testing/SteadyFlux_test")

###########################################################
# Setup
#######################################################
runname="SLIM_flux_test"
#fluxrate=0.00125 	#flux rate in m/hr
fluxrate=0.0025 	#flux rate in m/hr

#parflow grid
dx=1
dy=1
dz =0.05
pfdt=1  #parflow timestep
por = 0.25

#EcoSLIM inputs
t1=1
t2=3000
nt=t2-t1+1
icnp=10 #particles per cell for ic from slimin.txt
pcnp=10 #particles per evap-trans ic from slimin.txt
denh20=1000

###########################################################
# Total ET analysis
#######################################################
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
#dir="test_ETCenter"
fin=paste("./", dir, "/", runname, "_ET_output.txt", sep="")
ET=matrix(scan(fin , skip=1), ncol=7, byrow=T)
colnames(ET)=c("Time", "ET_age", "ET_comp1", "ET_comp2", "ET_comp3", "ET_mass",  "ET_np")
nstep=nrow(ET)
#lines(ET[,7], col='blue')
sum(ET[,7])
#test=ET[,7]

#Get a moving average for smoothing
window=501 #doing a centered movign average so it must be odd
half=floor(window/2)
ET_smooth=comp_smooth=rep(NA,nstep)

for(i in ceiling(window/2):(nstep-ceiling(window/2))){
	ET_smooth[i]=mean(ET[(i-half):(i+half),6])
	comp_smooth[i]=mean(ET[(i-half):(i+half),3])
}

#plot
par(mfrow=c(2,1))
plot(ET[,1], ET[,6], type='l', pch=18, xlab="Time Step", ylab="ET Mass out", col='darkgrey',  main="ET Mass") #, ylim=c(7,18))
abline(h=ETout, col=2, lwd=2)
abline(h=mean(ET[,6]), col='green', lwd=2)
lines(ET[,1], ET_smooth, col="blue", lwd=2)
legend('topright', legend=c(paste("PF (", ETout, ")", sep=""), paste("Mean Ecoslim (",round(mean(ET[,6]),2), ")", sep="" ), "Smoothed ET"), col=c('red', 'green', 'blue'), lwd=c(2,2))
temp=round(mean(ET[,6]),2)
print(mean(ET[,6])/ETout)

#quartz()
compplot=ET[,3]*1+ET[,4]*2
#plot(ET[,1], ET[,3], type='l', xlab="Time Step", ylab="Composition",col='darkgrey')
plot(ET[,1],compplot, type='l', xlab="Time Step", ylab="Composition",col='darkgrey')
lines(ET[,1], comp_smooth, col="blue", lwd=2)


###########################################################
# Gridded ET analysis
#######################################################
#ET_loc=1:10*30  #grabbing just the 10 cells with an ET flux 
#nloc=length(ET_loc)
nloc=10

#make a matrix with the ET timeseries for each grid cell
ET_grid=array(0, dim=c(nloc, 10, nt))
colnames(ET_grid)=c('X', 'Y', 'Z', 'ET_npart', ' ET_mass', ' ET_age', 'ET_comp', 'EvapTrans_Rate', 'Saturation', 'Porosity')
#loop over time and read files
for(i in 1:nt){
	t=i+t1-1
	print(c(i,t))
	#dir=paste("EcoSLIM", fluxrate, sep="")
	fin=sprintf("./%s/%s_ET_summary.%05d.txt", dir, runname, t)
	temp=matrix(scan(fin, skip=1), ncol=10, byrow=T)
	#ET_grid[,,i]=temp[ET_loc,]	
	ET_grid[,,i]=temp
}

#desired mass per grid
ET_PF=fluxrate*dx*dy*denh20  #ET from ParFlow EvapTrans

#Calculate the ratio of ET mas flux from ParFlow to the water mass in the cell
rat=matrix(0, ncol=nt, nrow=nloc) 
for(i in 1:nt){
 rat[, i]=(ET_PF/(ET_grid[,9,i]*ET_grid[,10,i]*dx*dy*dz*denh20)) 
}


#Make smothed timeseris for plotting
window=501 #doing a centered movign average so it must be odd
half=floor(window/2)
massg_smooth=compg_smooth=ageg_smooth=matrix(NA, nrow=nloc, ncol=nstep)
for(i in 1:nloc){
	for(j in ceiling(window/2):(nstep-ceiling(window/2))){
		massg_smooth[i,j]=mean(ET_grid[i, 5, (j-half):(j+half)])
		compg_smooth[i,j]=mean(ET_grid[i, 7, (j-half):(j+half)])
		ageg_smooth[i,j]=mean(ET_grid[i, 6, (j-half):(j+half)])
	}
}


#plotting mass and storage ratios
#Mass
par(mfrow=c(5,2), mar=c(4,4,2,2))
for(i in 1:10){
	ET_ES=mean(ET_grid[i,5,])
	title=paste("Cell, x=", i, " Mean Mass=", round(ET_ES,2), " Ratio=", round(mean(rat[i]),2), sep="")
	plot(t1:t2, ET_grid[i,5,], type='l', col='darkgrey', xlab='time', ylab='mass', main=title, axes=F, ylim=c(0,3))
	axis(1)
	axis(2)
	abline(h=ET_PF, col='red', lwd=2)
	abline(h=ET_ES, col='green')
	lines(ET[,1], massg_smooth[i,], col="blue", lwd=1)
	par(new=T)
	plot(t1:t2, rat[i,], type='l', col='purple',lty=2, axes=F, xlab="", ylab="", ylim=c(0,5))
	axis(4, col='purple')
	box()
}

#Composition
quartz()
par(mfrow=c(5,2), mar=c(4,4,2,2))
for(i in 1:10){
	plot(t1:t2, ET_grid[i,7,], type='l', col='darkgrey', xlab='time', ylab='Composition', main=paste(" Composition Cell, x=", i), ylim=c(0,2))
	lines(ET[,1], compg_smooth[i,], col="blue", lwd=1)
	abline(h=1.9, col='black', lty=3)
}

#Age
quartz()
par(mfrow=c(5,2), mar=c(4,4,2,2))
for(i in 1:10){
plot(t1:t2, ET_grid[i,6,], type='l', col='darkgrey', xlab='time', ylab='age', main=paste(" Age Cell, x=", i), ylim=c(0,3000))
}


#Combined plot using only the smoothed lines
#collist=c('#f7fbff','#deebf7','#c6dbef','#9ecae1','#6baed6','#4292c6','#2171b5','#08519c','#08306b')
collist=rainbow(10)

par(mfrow=c(1,3))
#mass
plot(t1:t2, massg_smooth[1,], type='l', col=collist[1], xlab='time', ylab='mass', main="Time averaged Mass flux", ylim=c(1.1,1.4), xlim=c(half, (t2-half)))
abline(h=ET_PF, lwd=2)
for(i in 1:nloc){
	lines(t1:t2, massg_smooth[i,], col=collist[i])
}
cellnames=paste("x=", 1:nloc)
legend('topright', cellnames, lty=rep(1, nloc), col=collist, ncol=2)

#composition
plot(t1:t2, compg_smooth[1,], type='l', col=collist[1], xlab='time', ylab='Composition', main="Composition", ylim=c(0.5,2), xlim=c(half, (t2-half)))
abline(h=2, lwd=2)
abline(h=1)
for(i in 1:nloc){
	lines(t1:t2, compg_smooth[i,], col=collist[i])
}
cellnames=paste("x=", 1:nloc)
#legend('bottomright', cellnames, lty=rep(1, nloc), col=collist)

#Age
plot(t1:t2, ageg_smooth[1,], type='l', col=collist[1], xlab='time', ylab='mass', main="Age", ylim=c(0,2500), xlim=c(half, (t2-half)))
for(i in 1:nloc){
	lines(t1:t2, ageg_smooth[i,], col=collist[i])
}
cellnames=paste("x=", 1:nloc)
#legend('topright', cellnames, lty=rep(1, nloc), col=collist)


###########################################################
# Gridded ET analysis  from point data 
#######################################################
#Read in the binary file with the exited particles
fin=paste("./", dir, "/SLIM_flux_test_exited_particles.bin" , sep="")
length=file.info(fin)$size/8
ncol=8
print(length)/ncol
to.read = file(fin,'rb')
part=matrix(0, nrow=(length/ncol), ncol=ncol)
for(i in 1:(length/ncol)){
	part[i,]= readBin(to.read, double(), endian='little', size=8,n=ncol)
}
close(to.read)
colnames(part)=c("Time", "x", "y", "z", "Ptime", "Mass", "Comp", "ExitStatus")
npart=nrow(part)


#Map the particles to the grid
ET_check=matrix(0, nrow=nt, ncol=10)
for(i in 1:npart){
	if(part[i,8]==2){
		ttemp=part[i,1]
		xtemp=ceiling(part[i,2])
		#xtemp=floor(part[i,2])+1
		#if(xtemp>10 & xtemp<=11){xtemp=10}
		ET_check[ttemp,xtemp]=ET_check[ttemp,xtemp]+part[i,6]
	}
}

#Check for differences with the gridded outputs
for(i in 1:10){
	#print(i)
	test=ET_grid[i,5,]
	test2=ET_check[,i]
	dif=test-test2
	print(paste("Grid_cell:", i, " Range_Mass_Dif:", round(min(dif),3), ",", round(max(dif),3), " Count_Mass_Dif:", length(which(abs(dif)>0.01)), sep=""))
	print(which(abs(dif)>0.01))
}

#plot the total flux gridded and not
plot(ET[,1], ET[,6], type='l', pch=18, xlab="Time Step", ylab="ET Mass out", col='darkgrey',  main="ET Mass", ylim=c(7,18))
abline(h=ETout, col=2, lwd=2)
lines(ET[,1],apply(ET_check, 1, sum),col='blue' )
range(ET[,6]-apply(ET_check,1,sum))
diftest=ET[,6]-apply(ET_check,1,sum)

#look at the differences in cell 10
i=10
plot(ET_grid[i,5,], ET_check[,i])
lines(c(0,100), c(0,100))
plot(t1:t2, ET_grid[i,5,], type='l', col='darkgrey', xlab='time', ylab='mass', main=title, axes=F, ylim=c(0,3))
mean(ET_grid[i,5,])
mean(ET_check[,i])	




##############


abline(h=mean(ET[,6]), col='green', lwd=2)
lines(ET[,1], ET_smooth, col="blue", lwd=2)
legend('topright', legend=c(paste("PF (", ETout, ")", sep=""), paste("Mean Ecoslim (",round(mean(ET[,6]),2), ")", sep="" ), "Smoothed ET"), col=c('red', 'green', 'blue'), lwd=c(2,2))
temp=round(mean(ET[,6]),2)
print(mean(ET[,6])/ETout)



# plot
par(mfrow=c(5,2), mar=c(4,4,2,2))
for(i in 1:10){
	#plot(ET_grid[i,5,], ET_check[,i])
	#lines(c(0,100), c(0,100))

	ET_ES=mean(ET_grid[i,5,])
	title=paste("Cell, x=", i, " Mean Mass=", round(ET_ES,2), " Ratio=", round(mean(rat[i]),2), sep="")
	plot(t1:t2, ET_grid[i,5,], type='l', col='darkgrey', xlab='time', ylab='mass', main=title, axes=F, ylim=c(0,3))
	axis(1)
	axis(2)
	abline(h=ET_PF, col='red', lwd=2)
	abline(h=ET_ES, col='green')
	lines(ET[,1], massg_smooth[i,], col="blue", lwd=1)
	par(new=T)
	plot(t1:t2, rat[i,], type='l', col='purple',lty=2, axes=F, xlab="", ylab="", ylim=c(0,5))
	axis(4, col='purple')
	box()
}









t=2034
ilist=which(part[,1]==t)
part[ilist, ]
ET_grid[,5,t]
ET_check[t,]
round(part[ilist, 2],2)


#testing the rounding
test1=ceiling(part[,2])
test2=floor(part[,2])+1
ilist=which(diftest!=0)
part[ilist,2]


