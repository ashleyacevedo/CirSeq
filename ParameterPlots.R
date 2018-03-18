args = commandArgs(T)
data = read.table(paste(args[1],"/QualityMetrics.txt",sep=""), sep="\t", header=TRUE)
attach(data)

#mutation frequency
pdf(paste(args[1], "/MutationFrequency.pdf", sep=""))
plot(Q, Mismatches/TotalBases, type = "p",log="y", xlab="Average Quality Score",ylab="Mutation Frequency", axes=F, pch=16, col="black", ylim=c(.0001,1), xlim=c(0,40),cex=1.5)
axis(1, at = c(0,10,20,30,40), labels = expression(0,10,20,30,40))
axis(2, at = 10^c(-4,-3,-2,-1,0), labels = expression(10^-4,10^-3, 10^-2,10^-1, 1))
axis(2, at = c(.0002,.0003,.0004,.0005,.0006,.0007,.0008,.0009), labels = c("","","","", "","","",""), tck=-.01, lwd=0, lwd.ticks=1)
axis(2, at = c(.002,.003,.004,.005,.006,.007,.008,.009), labels = c("","","","", "","","",""), tck=-.01, lwd=0, lwd.ticks=1)
axis(2, at = c(.02,.03,.04,.05,.06,.07,.08,.09), labels = c("","","","", "","","",""), tck=-.01, lwd=0, lwd.ticks=1)
axis(2, at = c(.2,.3,.4,.5,.6,.7,.8,.9), labels = c("","","","", "","","",""), tck=-.01, lwd=0, lwd.ticks=1)

#tstv
pdf(paste(args[1], "/Transition:Transversion.pdf", sep=""))
plot(Q, Ts/Tv,type ="p", xlab = "Average Quality Score", ylab= "Transition:Transversion", pch=16, axes=F, col="black",ylim=c(0,40), xlim=c(0,40), cex=1.5)
axis(1, at = c(0,10,20,30,40), labels = expression(0,10,20,30,40))
axis(2, at = c(0,10,20,30,40), labels = expression(0,10,20,30,40))

detach(data)