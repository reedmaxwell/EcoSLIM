
# x, y, z: data
# slices: number of horizontal slices through the data
# lboost: coefficient to increase the height of the lines
# gboost: coefficient to increase the height of the graph (ylim)
# xinc: horizontal offset for each succesive slice
# (typically something like 1/80)
# yinc: vertical offset for each succesive slice
# bcol: background color
# fcol: fill color for each slice (polygon)
# lcol: line color for each slice
# lwidth: line width
# cutprop:
# Vary transarency
# transprop=FALSE
# tmax = .9
# tmin = .2

slicedens<-function(x,y,z=NULL,slices=50,lboost=1,gboost=1,
                    xinc=0,yinc=0.01, bcol="black", fcol="black",
                    lcol="white",lwidth=1, cutprop=FALSE,
                    transprop=FALSE, tmax = .8, tmin = .05,
                    heightprop=FALSE, xlab=FALSE, ylab=FALSE, height=100,main=FALSE,yaxlab=FALSE) {

  # This function takes a matrix of one or more rgb sets
  # as well as a degree [0,1] and returns a combined
  # color.
  color.mix <- function(cols, degree=0) {
    if (is.null(nrow(cols))) {
      if (class(cols)=="numeric")
        return(rgb(cols[1],cols[2],cols[3],cols[4]))
      return(cols)
    }
    # Define a function to find elementwise minimum
    (deg <- degree*(nrow(cols)-1)+1)
    emin <- function(x, y=0) apply(cbind(x, y), 1, min)
    (r <- 1-emin(abs(1:nrow(cols)-deg),1))
    (comb <- apply(cols*r,2,sum))
    mm <- function(x) max(min(x,1),0)
    rgb(mm(comb[1]),
        mm(comb[2]),
        mm(comb[3]),
        mm(comb[4]))
  }

  ycut<-min(y)+((0:(slices))*(max(y)-min(y))/slices)

  height<-gboost*((slices*yinc)+max(density(x)$y))

  #plot( c(min(x)-((max(x)-min(x))/10),max(x)+((max(x)-min(x))/10)),
   #     c(0,height),
    #    xaxt="n",yaxt="n",ylab="",xlab="")
  plot( c(min(x)-((max(x)-min(x))/10),max(x)+((max(x)-min(x))/10)),
        c(0,height),
        yaxt="n",ylab=ylab,xlab=xlab,yaxt='n',xlim=c(0,10),main=main)
        
#        axis(2, at=(0:(slices-1)*yinc), labels=yaxlab)
        
#        axis(2, at=c(0,1,2,3,4,5,6,7,8,9,10,11), #labels=c("OCT","NOV","DEC","JAN","FEB","MAR","APR","MAY","JUN","JUL","AUG","SEP"))
axis(2, at=(0:(slices-1)*yinc), labels=F)  labseq=seq(0, (slices-1), 2)  axis(2, at=(labseq*yinc), labels=yaxlab[(labseq+1)])

  rect(par("usr")[1],par("usr")[3],par("usr")[2],par("usr")[4],col=bcol)

  # Calcuate the 'degree' for each z value which will be used to
  # choose the color of each slice.
  if (length(z)==length(y)) {
    zmin <- min(z)
    zmax <- max(z)
    zrange <- max(zmax-zmin)
  }

  # Define ifcol and ilfol for later reference.
  # Unless noted otherwise, degree=0
  # Meaning the first color will be selected from
  # the rgb color matrix.
  ifcol <- fcol; ilcol <- lcol;  zdeg <- 0

  # Define the degree which is the color degree that each slice will
  # contain
  if (length(z)==length(y)){
    meanz <- NULL
    for(i in slices:1)
      meanz[i]<- mean(z[y>=min(ycut[i]) & y<max(ycut[i+1])])
    zdegree<-(meanz-min(meanz, na.rm=TRUE))/
      (max(meanz, na.rm=TRUE)-min(meanz, na.rm=TRUE))
  }

  # Loop through and plot slices
  for(i in slices:1) {
    miny<-ycut[i]
    maxy<-ycut[i+1]


    gx<-(i-1)*(max(x)-min(x))*xinc

    if (cutprop) {
      slLength <- slices*sum(y>=miny & y<maxy)/length(y)
      if (i==slices) gy <- (i-1)*(height)*yinc
      if (i<slices)  gy <- gyLast-(height)*yinc*slLength
      gyLast <- gy
    } else {
    		#gy<-(i-1)*(height)*yinc
    		gy<-(i-1)*yinc
    	}
    	
    	
    if (transprop) {
      trange <- tmax-tmin
      if (is.null(nrow(ifcol)))
        ifcol[4] <- min(trange*slices*sum(y>=miny & y<maxy)/length(y)+tmin,tmax)
      if (!is.null(nrow(ifcol)))
        ifcol[,4] <- min(trange*slices*sum(y>=miny & y<maxy)/length(y)+tmin,tmax)
    }



    # If z is an input vector then grab the color degree from it
    if (length(z)==length(y)) zdeg<-zdegree[i]

    # Added the try because distributions without defined upper
    # and lower bounds can give the density function trouble.
    try({
      # Use the color.mixer function to select a color
      fcol<-color.mix(ifcol, zdeg);
      lcol<-color.mix(ilcol, zdeg);
      # Calculte density curves and plot them
      dd<-density(x[y>=miny & y<maxy]);
      if (heightprop) vscale <- lboost*slices*sum(y>=miny & y<maxy)/length(y)
      if (!heightprop) vscale <- lboost
      if(dd$y[1]>dd$y[length(dd$y)]){dd$y[1]=dd$y[length(dd$y)]}
      polygon(dd$x+gx,vscale*dd$y+gy,col=fcol, border=fcol);
      lines(dd$x+gx,vscale*dd$y+gy,col=lcol,lwd=lwidth);
      #print(paste(gx, gy));
    })
  }
}

