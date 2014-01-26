#!/usr/bin/env Rscript

subsample <- function(the.data, key.col, min, max) {
    the.data[sapply(n[,key.col], function(x){x > min && x < max}),]
}

args <- commandArgs()
infile = args[6]

n <- read.table(infile)
colnames(n) <- c("gi", "identity", "coverage")

pdf(paste(infile,".seqquery.plot.pdf",sep=""))

plot(n[,2:3], axes=0)
axis(side=1,at=seq(0.0,10,0.1))
axis(side=2,at=seq(0.0,10,0.1))

for(i in seq(0, 10, 0.2)) {
    abline(h=i, col=rgb(1,0,0,0.5))
    abline(v=i, col=rgb(0,0,1,0.5))
}

for(i in seq(0.1, 10, 0.2)) {
    abline(h=i, col=rgb(1,0,0,0.25))
    abline(v=i, col=rgb(0,0,1,0.25))
}

dev.off()
