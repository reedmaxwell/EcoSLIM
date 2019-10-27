#function to read EcoSLIM binary file OUTPUTS
# Inputs -
# filename - name of .bin file to reads
# filetype = 'exited' or 'restart'
# nind = The number of indicators specified in the EcloSLIM driver file

ES_read <- function(filename, type = "exited", nind=0){
if(nind>0){
  ind_names=c(paste("IndAge", 1:nind, sep=""), paste("IndLen", 1:nind, sep=""))
}
  if(type == "restart"){
    to.read = file(filename,"rb")
    npart=readBin(to.read, integer(), endian="little",size=4,n=1)
    pid=readBin(to.read, integer(), endian="little",size=4,n=1)
    print(npart)
    print(pid)

    n.col=17+nind*2

    #NOTE: These are written out transposed from the exited particles file see below
    data = matrix(0,ncol=n.col,nrow=npart,byrow=F)
    for (i in 1:n.col) {
      data[,i] = readBin(to.read, double(), endian="little",size=8,n=npart)
    }
    close(to.read)
    data[1,]
    particle_restart <- data.frame(data)
    #X, Y, Z, Age, Saturated_Age, Mass, Source, Status, Particle_Concentration, Exit_Status, ParticleID,
    #InitialX, InitialY, InitialZ, InsertTime, Path_Length, Saturated_Path_Length
    col.names=c("X","Y","Z","age","sat_age","mass","source","status", "conc","exit_status",
              "ID","init_X","init_Y","init_Z","ins_time","path_len","spath_len")
    if(nind>0){col.names=c(col.names,ind_names)}
    colnames(particle_restart) <- col.names


    return(particle_restart)

  } else if(type == "exited"){
    #lines = file.info(filename)$size
    n.col=16+nind*2
    lines = (file.info(filename)$size/8)/(n.col)

    data = matrix(0,ncol=n.col,nrow=lines,byrow=F)
    to.read = file(filename,"rb")
    for (i in 1:lines) {
      data[i,1:n.col] = readBin(to.read, single(), endian="little",size=8,n=n.col)
    }
    readBin(to.read, double(), endian="little",size=8,n=n.col)
    close(to.read)
    #1 = Time #2 = X #3 = Y  #4 = Z #5 = age #6 = mass #7 = source #8 = out(1), ET(2)

    exited_particles <- data.frame(data)
    # Time, X, Y, Z, Age, Saturated_Age, Mass, Source, Exit_Status, ParticleID,
    #InitialX, InitialY, InitialZ, InsertTime, Path_Length, Saturated_Path_Length
    col.names= c("Time","X","Y","Z","age","sat_age","mass","source","out_as","ID","init_X",
                                    "init_Y","init_Z","ins_time","path_len","spath_len")

    if(nind>0){col.names=c(col.names,ind_names)}
    colnames( exited_particles) <- col.names

    return(exited_particles)

  } else {
    print("error")
  }
}
