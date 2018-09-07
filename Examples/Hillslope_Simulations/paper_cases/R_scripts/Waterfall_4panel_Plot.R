### Waterfall and ET Hillslope plots
### LEC, RMM 25-May-18, (edited 7-Sep-18)

##This script loads in EcoSLIM output and makes Figures 4 and 5 from Maxwell et al Ecohydrology 2018

library("MASS")
library(fields)
rm(list=ls())

## you may need to edit the working directory  uncomment and adjust as needed
##setwd('~/EcoSLIM/Examples/Hillslope_Simulations/paper_cases/R_scripts')
source("./slicedens.R")

months=c("OCT","NOV","DEC","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP")

## make waterfall four panel (Fig 4)
fout = 'Figure4_waterfall_4-panel.pdf'
pdf(file=fout,width=9, height=9)

par(mfrow=c(2,2))

#Read EcoSlim output files and put in independent matrices for plotting
## load trees with ER forcing
fin='./ER_hillslope_trees/SLIM_hillslope_ER_trees_exited_particles.bin'
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

#3D Histograms of particle ages
out = which( (part[,8] == 1) & (part[,1]>=(19*8760)))
et = which((part[,8] == 2) & (part[,1]>=(19*8760)))
den3d_out_er_trees <- kde2d( (part[out,5])/24/365,part[out,1]/24/365, n = 50)

den3d_et_er_trees_time <- kde2d( (part[et,5])/24/365,part[et,1]/24/365, n = 50)
den3d_et_er_trees_space <- kde2d( (part[et,5])/24/365,part[et,2], n = 100)

#Outflow
x=part[out,5]/24/365
y=part[out,1]/24/365
z <- y
trans=0.7
fcol <- rbind(c(0,.1,.5,trans), c(.3,.8,.8,trans), c(1,1,0, trans))
fcol=fcol[3:1,]
lcol <- rbind(c(0,.3,.3,.8), c(.1,.1,.2,.7), c(0,0,1,.65))


slicedens(x,y,z,fcol=fcol, bcol='white', lcol=lcol,gboost=.9, yinc=1, slices=12, xlab="Outflow Residence Time [YR] ",ylab="",main="A) ER Trees Outflow",yaxlab=months)
        
#ET
x=part[et,5]/24/365
y=part[et,1]/24/365
z <- y

slicedens(x,y,z,fcol=fcol, bcol='white', lcol=lcol,gboost=.90, yinc=7, slices=12, xlab="ET Residence Time [YR]",ylab="", main="B) ER Trees ET",yaxlab=months)


## load shrubs with LW forcing
fin='../LW_hillslope_shrub/SLIM_hillslope_LW_shrub_exited_particles.bin'
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

#3D Histograms of particle ages
out = which( (part[,8] == 1) & (part[,1]>=(19*8760)))
et = which((part[,8] == 2) & (part[,1]>=(19*8760)))
den3d_out_lw_shrubs <- kde2d( (part[out,5])/24/365,part[out,1]/24/365, n = 50)
den3d_et_lw_shrubs_time <- kde2d( (part[et,5])/24/365,part[et,1]/24/365, n = 50)
den3d_et_lw_shrubs_space <- kde2d( (part[et,5])/24/365,part[et,2], n = 100)


## load trees with LW forcing
fin='../LW_hillslope_trees/SLIM_hillslope_LW_trees_exited_particles.bin'
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

#3D Histograms of particle ages
out = which( (part[,8] == 1) & (part[,1]>=(19*8760)))
et = which((part[,8] == 2) & (part[,1]>=(19*8760)))
den3d_out_lw_trees <- kde2d( (part[out,5])/24/365,part[out,1]/24/365, n = 50)
den3d_et_lw_trees_time <- kde2d( (part[et,5])/24/365,part[et,1]/24/365, n = 50)
den3d_et_lw_trees_space <- kde2d( (part[et,5])/24/365,part[et,2], n = 100)

#Outflow
x=part[out,5]/24/365
y=part[out,1]/24/365
z <- y
trans=0.7
fcol <- rbind(c(0,.1,.5,trans), c(.3,.8,.8,trans), c(1,1,0, trans))
fcol=fcol[3:1,]
lcol <- rbind(c(0,.3,.3,.8), c(.1,.1,.2,.7), c(0,0,1,.65))


slicedens(x,y,z,fcol=fcol, bcol='white', lcol=lcol,gboost=.9, yinc=1, slices=12, xlab="Outflow Residence Time [YR]",main="C) LW Trees Outflow",ylab="",yaxlab=months)
          
#ET
x=part[et,5]/24/365
y=part[et,1]/24/365
z <- y

slicedens(x,y,z,fcol=fcol, bcol='white', lcol=lcol,gboost=.90, yinc=7, slices=12, xlab="ET Residence Time [YR]",main="D) LW Trees ET",ylab="",yaxlab=months)

dev.off()

## load shrubs with ER forcing
fin='./ER_hillslope_shrub/SLIM_hillslope_ER_shrub_exited_particles.bin'
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

#3D Histograms of particle ages
out = which( (part[,8] == 1) & (part[,1]>=(19*8760)))
et = which((part[,8] == 2) & (part[,1]>=(19*8760)))
den3d_out_er_shrubs <- kde2d( (part[out,5])/24/365,part[out,1]/24/365, n = 50)
den3d_et_er_shrubs_time <- kde2d( (part[et,5])/24/365,part[et,1]/24/365, n = 50)
den3d_et_er_shrubs_space <- kde2d( (part[et,5])/24/365,part[et,2], n = 100)

## make ET distribution four panel (Fig 5)
fout = 'Figure5_ET_hillsllope_4-panel.pdf'
pdf(file=fout)

par(mfrow=c(2,2))

image(den3d_et_er_shrubs_space,xlab = "ET residence time [yr]", ylab = "X Location Hillslope [m]",col = heat.colors(9),zlim=c(0.0001,0.05),xlim=c(0,10),breaks=c(0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05),main="A) ER Shrubs")

image(den3d_et_er_trees_space,xlab = "ET residence time [yr]", ylab = "X Location Hillslope [m]",col = heat.colors(9),zlim=c(0.0001,0.05),xlim=c(0,10),breaks=c(0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05),main="B) ER Trees")

image(den3d_et_lw_shrubs_space,xlab = "ET residence time [yr]", ylab = "X Location Hillslope [m]",col = heat.colors(9),zlim=c(0.0001,0.05),xlim=c(0,10),breaks=c(0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05),main="C) LW Shrubs")

image(den3d_et_lw_trees_space,xlab = "ET residence time [yr]", ylab = "X Location Hillslope [m]",col = heat.colors(9),zlim=c(0.0001,0.05),xlim=c(0,10),breaks=c(0.000001,0.000005,0.00001,0.00005,0.0001,0.0005,0.001,0.005,0.01,0.05),main="D) LW Trees")

image.plot(den3d_et_lw_trees_space, zlim=c(0.0001,0.05), cex=0.7,col=(heat.colors(9)), legend.only=T, legend.width=0.5, legend.shrink=0.5, legend.args=list(text='            Density [-]',cex=0.85,side=3,line=0.1), smallplot=c(.65,.68, .55,.75))
dev.off()


