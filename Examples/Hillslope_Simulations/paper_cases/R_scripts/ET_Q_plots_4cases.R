### plot four panels with ET and Q for all four cases
### RMM 25-May-18, (edited 7-Sep-18)

##This script loads in EcoSLIM output and makes Figures 2 and 3 from Maxwell et al Ecohydrology 2018

## you may need to edit the working directory, uncomment and adjust as needed
## setwd('~/EcoSLIM/Examples/Hillslope_Simulations/paper_cases/R_scripts')

## load shrubs with ER forcing
filename='../ER_hillslope_shrub/SLIM_hillslope_ER_shrub_ET_output.txt'
ET_er_shrub = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../ER_hillslope_shrub/SLIM_hillslope_ER_shrub_surface_outflow.txt'
Outflow_er_shrub = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../ER_hillslope_shrub/SLIM_hillslope_ER_shrub_water_balance.txt'
PET_er_shrub = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=3,byrow=T)


ET_day_er_shrub = matrix(0,nrow=20*365,ncol=7,byrow=T)
Outflow_day_er_shrub =  matrix(0,nrow=20*365,ncol=7,byrow=T)
PET_day_er_shrub =  matrix(0,nrow=20*365,ncol=3,byrow=T)

for (i in 1:(20*365)) {
	for (j in 1:7) {
	ET_day_er_shrub[i,j] = mean(ET_er_shrub[((i-1)*24+1):(i*24),j])
	Outflow_day_er_shrub[i,j] = mean(Outflow_er_shrub[((i-1)*24+1):(i*24),j])
	}
}

for (i in 1:(20*365)) {
	for (j in 1:3) {
	PET_day_er_shrub[i,j] = mean(PET_er_shrub[((i-1)*24+1):(i*24),j])
	}
}
options(scipen=5)

## load trees with ER forcing
filename='../ER_hillslope_trees/SLIM_hillslope_ER_trees_ET_output.txt'
ET_er_trees = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../ER_hillslope_trees/SLIM_hillslope_ER_trees_surface_outflow.txt'
Outflow_er_trees = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../ER_hillslope_trees/SLIM_hillslope_ER_trees_water_balance.txt'
PET_er_trees = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=3,byrow=T)


ET_day_er_trees = matrix(0,nrow=20*365,ncol=7,byrow=T)
Outflow_day_er_trees =  matrix(0,nrow=20*365,ncol=7,byrow=T)
PET_day_er_trees =  matrix(0,nrow=20*365,ncol=3,byrow=T)

for (i in 1:(20*365)) {
	for (j in 1:7) {
	ET_day_er_trees[i,j] = mean(ET_er_trees[((i-1)*24+1):(i*24),j])
	Outflow_day_er_trees[i,j] = mean(Outflow_er_trees[((i-1)*24+1):(i*24),j])
	}
}

for (i in 1:(20*365)) {
	for (j in 1:3) {
	PET_day_er_trees[i,j] = mean(PET_er_trees[((i-1)*24+1):(i*24),j])
	}
}
options(scipen=5)

## load shrubs with LW forcing
filename='../LW_hillslope_shrub/SLIM_hillslope_LW_shrub_ET_output.txt'
ET_lw_shrub = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../LW_hillslope_shrub/SLIM_hillslope_LW_shrub_surface_outflow.txt'
Outflow_lw_shrub = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../LW_hillslope_shrub/SLIM_hillslope_LW_shrub_water_balance.txt'
PET_lw_shrub = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=3,byrow=T)


ET_day_lw_shrub = matrix(0,nrow=20*365,ncol=7,byrow=T)
Outflow_day_lw_shrub =  matrix(0,nrow=20*365,ncol=7,byrow=T)
PET_day_lw_shrub =  matrix(0,nrow=20*365,ncol=3,byrow=T)

for (i in 1:(20*365)) {
	for (j in 1:7) {
	ET_day_lw_shrub[i,j] = mean(ET_lw_shrub[((i-1)*24+1):(i*24),j])
	Outflow_day_lw_shrub[i,j] = mean(Outflow_lw_shrub[((i-1)*24+1):(i*24),j])
	}
}

for (i in 1:(20*365)) {
	for (j in 1:3) {
	PET_day_lw_shrub[i,j] = mean(PET_lw_shrub[((i-1)*24+1):(i*24),j])
	}
}

## load trees with LW forcing
filename='../LW_hillslope_trees/SLIM_hillslope_LW_trees_ET_output.txt'
ET_lw_trees = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../LW_hillslope_trees/SLIM_hillslope_LW_trees_surface_outflow.txt'
Outflow_lw_trees = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=7,byrow=T)

filename='../LW_hillslope_trees/SLIM_hillslope_LW_trees_water_balance.txt'
PET_lw_trees = matrix(scan(file=filename,skip=1),nrow=8760*20,ncol=3,byrow=T)


ET_day_lw_trees = matrix(0,nrow=20*365,ncol=7,byrow=T)
Outflow_day_lw_trees =  matrix(0,nrow=20*365,ncol=7,byrow=T)
PET_day_lw_trees =  matrix(0,nrow=20*365,ncol=3,byrow=T)

for (i in 1:(20*365)) {
	for (j in 1:7) {
	ET_day_lw_trees[i,j] = mean(ET_lw_trees[((i-1)*24+1):(i*24),j])
	Outflow_day_lw_trees[i,j] = mean(Outflow_lw_trees[((i-1)*24+1):(i*24),j])
	}
}

for (i in 1:(20*365)) {
	for (j in 1:3) {
	PET_day_lw_trees[i,j] = mean(PET_lw_trees[((i-1)*24+1):(i*24),j])
	}
}




## normalize ET and Q to mm/day
# as total amount, convert to mm/d
norm = (1000*200*1)/(24*1000)

#plot flow and ET residence time comp Y20
#
fout = 'Figure2_Q-y20-4-panel.pdf'
pdf(file=fout)
par(mfrow=c(2,2))
plot(Outflow_day_er_shrub[,1]/24/365,(Outflow_day_er_shrub[,3]+Outflow_day_er_shrub[,4]+Outflow_day_er_shrub[,5])*Outflow_day_er_shrub[,2]/24/365,cex=0.5,col="black",type='l',xlab='', ylab='Outflow [mm/d] and Residence Time [y]',xlim=c(19,20.),ylim=c(0,5.5),main="A) ER Shrubs",xaxt='n')
axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))

par(new=T)
plot(Outflow_day_er_shrub[,1]/24/365,Outflow_day_er_shrub[,6]/norm,cex=0.5,col='blue',type='p',axes=F,xlim=c(19,20.),ylim=c(0,5.5),xlab="",ylab="") 
legend('topright',c('Outflow Residence Time'),col=c('black'),lty=1,box.lwd = 0,cex=0.7)
legend('topleft',c('Outflow'),col=c('blue'),pch=21,box.lwd = 0, cex=0.7)

plot(Outflow_day_er_trees[,1]/24/365,(Outflow_day_er_trees[,3]+Outflow_day_er_trees[,4]+Outflow_day_er_trees[,5])*Outflow_day_er_trees[,2]/24/365,cex=0.5,col="black",type='l',xlab='', ylab='Outflow [mm/d] and Residence Time [y]',xlim=c(19,20.),ylim=c(0,5.5),main="B) ER Trees",xaxt='n')
axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))

par(new=T)
plot(Outflow_day_er_trees[,1]/24/365,Outflow_day_er_trees[,6]/norm,cex=0.5,col='blue',type='p',axes=F,xlim=c(19,20.),ylim=c(0,5.5),xlab="",ylab="") 

plot(Outflow_day_lw_shrub[,1]/24/365,(Outflow_day_lw_shrub[,3]+Outflow_day_lw_shrub[,4]+Outflow_day_lw_shrub[,5])*Outflow_day_lw_shrub[,2]/24/365,cex=0.5,col="black",type='l',xlab='', ylab='Outflow [mm/d] and Residence Time [y]',xlim=c(19,20.),ylim=c(0,5.5),main="C) LW Shrubs",xaxt='n')
axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))

par(new=T)
plot(Outflow_day_lw_shrub[,1]/24/365,Outflow_day_lw_shrub[,6]/norm,cex=0.5,col='blue',type='p',axes=F,xlim=c(19,20.),ylim=c(0,5.5),xlab="",ylab="") 

plot(Outflow_day_lw_trees[,1]/24/365,(Outflow_day_lw_trees[,3]+Outflow_day_lw_trees[,4]+Outflow_day_lw_trees[,5])*Outflow_day_lw_trees[,2]/24/365,cex=0.5,col="black",type='l',xlab="",xaxt='n', ylab='',xlim=c(19,20.),ylim=c(0,5.5),main="D) LW Trees")
axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))

par(new=T)
plot(Outflow_day_lw_trees[,1]/24/365,Outflow_day_lw_trees[,6]/norm,cex=0.5,col='blue',type='p',axes=F,xlim=c(19,20.),ylim=c(0,5.5),xlab="",ylab="Outflow [mm/d] and Residence Time [y]") 

dev.off()

## normalize ET and Q to mm/day
# as total amount, convert to mm/d
norm = (1000*200*1)/(24*1000)

#plot  ET residence time comp Y20
fout = 'Figure3_ET-y20-4-panel.pdf'
pdf(file=fout)

par(mfrow=c(2,2))
plot(ET_day_er_shrub[,1]/24/365,(ET_day_er_shrub[,3]+ET_day_er_shrub[,4]+ET_day_er_shrub[,5])*ET_day_er_shrub[,2]/24/365,cex=0.5,col="black",type='l',xlab='', ylab='ET [mm/d], ET Residence Time [y]',xlim=c(19,20.),ylim=c(0,1.5),main="A) ER Shrubs",xaxt='n')
#points(Outflow_day_er_shrub[,1]/24/365,(Outflow_day_er_shrub[,3]+Outflow_day_er_shrub[,4])*Outflow_day_er_shrub[,2]/24/365,cex=0.5,col="blue",type='l')
par(new=T)
plot(ET_day_er_shrub[,1]/24/365,ET_day_er_shrub[,6]/norm,cex=0.5,col="blue",type='p',axes=F,xlim=c(19,20.),ylim=c(0,1.5),xlab="",ylab="") 
axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))
#axis(4, ylim=c(0,15/norm))
#mtext(4,text="ET [mm/d]",level=2)
legend('topright',c('ET Residence Time'),col=c('black'),lty=1,box.lwd = 0)
legend('topleft',c('ET'),col=c('blue'),pch=21,box.lwd = 0)

plot(ET_day_er_trees[,1]/24/365,(ET_day_er_trees[,3]+ET_day_er_trees[,4]+ET_day_er_trees[,5])*ET_day_er_trees[,2]/24/365,cex=0.5,col="black",type='l',xlab='', ylab='ET [mm/d], ET Residence Time [y]',xlim=c(19,20.),ylim=c(0,1.5),main="B) ER Trees",xaxt='n')
#points(Outflow_day_er_trees[,1]/24/365,(Outflow_day_er_trees[,3]+Outflow_day_er_trees[,4])*ET_day_er_trees[,2]/24/365,cex=0.5,col="blue",type='l')
par(new=T)
plot(ET_day_er_trees[,1]/24/365,ET_day_er_trees[,6]/norm,cex=0.5,col='blue',type='p',axes=F,xlim=c(19,20.),ylim=c(0,1.5),xlab="",ylab="") 
axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))

plot(ET_day_lw_shrub[,1]/24/365,(ET_day_lw_shrub[,3]+ET_day_lw_shrub[,4]+ET_day_lw_shrub[,5])*ET_day_lw_shrub[,2]/24/365,cex=0.5,col="black",type='l',xlab='', ylab='ET [mm/d], ET Residence Time [y]',xlim=c(19,20.),ylim=c(0,1.5),main="C) LW Shrubs",xaxt='n')
axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))

par(new=T)
plot(ET_day_lw_shrub[,1]/24/365,ET_day_lw_shrub[,6]/norm,cex=0.5,col='blue',type='p',axes=F,xlim=c(19,20.),ylim=c(0,1.5),xlab="",ylab="") 

plot(ET_day_lw_trees[,1]/24/365,(ET_day_lw_trees[,3]+ET_day_lw_trees[,4]+ET_day_lw_trees[,5])*ET_day_lw_trees[,2]/24/365,cex=0.5,col="black",type='l',xlab="",xaxt='n', ylab='ET [mm/d], ET Residence Time [y]',xlim=c(19,20.),ylim=c(0,1.5),main="D) LW Trees")
par(new=T)
plot(ET_day_lw_trees[,1]/24/365,ET_day_lw_trees[,6]/norm,cex=0.5,col='blue',type='p',axes=F,xlim=c(19,20.),ylim=c(0,1.5),xlab="",ylab="") 

axis(1, at=c(19+1/12,19+3/12,19+5/12,19+7/12,19+9/12,19+11/12), labels=c("Oct","Dec","Feb","Apr","Jun","Aug"))

dev.off()



