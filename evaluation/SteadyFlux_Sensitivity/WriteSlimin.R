#This script writes an input file for EcoSlim given command line arguments
#UNDER CONSTRUCTION

args <- commandArgs(trailingOnly = TRUE);
day=args[1]
print(day)


print("test")
fout="slimin.txt"
fout2=paste(fout, day)
wd=getwd()
print(wd)
print(fout2)
write.table(c(0,0), fout2, row.names=F, col.names=F, append=F)