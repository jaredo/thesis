library(RColorBrewer);cols <- paste(brewer.pal(8,"Set1"),'CC',sep='')[-6];library(xtable)
black <- rgb(0,0,0,alpha=0.1)
cohorts <- c('clean-carl-flip_b37','orkney889_b37','valborbera_b37','clean-fvg-flip_b37','KORCULA_for1000G','SPLIT_for1000G','VIS_for1000G')
slrpIBDmatrix <- function(fname,ids,minibd=1,mx=2*183.11) {
  slrp <- read.table(fname,as.is=TRUE)
  slrp$V3[slrp$V3<minibd] <- 0
#  ids <- sort(unique(c(slrp$V1,slrp$V2)))
  N <- length(ids)*2
  ibdmat <- matrix(0,nrow=N,ncol=N)
  for(i in 1:nrow(slrp)) {
    id1 <- slrp[i,1]+1
    id2 <- slrp[i,2]+1
    size <- slrp[i,3]
    ibdmat[id1,id2] <-ibdmat[id1,id2] + size
    ibdmat[id2,id1] <-ibdmat[id1,id2]
  }
  return((ibdmat[seq(1,N,2),seq(1,N,2)]+ibdmat[seq(2,N,2),seq(2,N,2)]  )/mx)
}

slrp.propibd <- function(fname,ids,markerfile,threshold,rrm) {
  markers <- read.table(markerfile,as.is=TRUE)[,4]
  slrp <- read.table(fname,as.is=TRUE)
#  slrp$V3[slrp$V3<minibd] <- 0
  N <- length(ids)*2
  ibdprop <- matrix(0,nrow=N,ncol=length(markers))
  for(i in 1:nrow(slrp)) {
    id1 <- slrp[i,1]+1
    id2 <- slrp[i,2]+1
    if(!(missing(threshold) & missing(rrm))) {
      if (rrm[ceiling(id1/2),ceiling(id2/2)]<threshold) {
        ii <- markers>= slrp[i,6] & markers<= slrp[i,7]
        ibdprop[id1,ii] <- ibdprop[id1,ii]+1
        ibdprop[id2,ii] <- ibdprop[id2,ii]+1
      }
    }
  }
  return(ibdprop)
}


readrrm <- function(rrmfilename,famfilename,ids) {
  fam <- read.table(famfilename,as.is=T)
  rrm <- matrix(scan(rrmfilename),nrow=nrow(fam),ncol=nrow(fam))
  if(!missing(ids)) {
    ii <- match(ids,fam[,2])
    rrm <- rrm[ii,ii]
  }
#  diag(rrm) <-diag(rrm)
  diag(rrm) <- NA
  rrm
}


slrp <- list()
for(fname in cohorts) {
    slrp[[fname]] <- read.table(paste("~/Dropbox/dphil/isolates/unrelated/",fname,"-chr10-unrelated-slrp_yield.txt",sep=""),header=TRUE)
      slrp[[fname]]$yield <- 100 * slrp[[fname]]$NPHASED / slrp[[fname]]$NHET
      fam.all <- read.table(paste("~/Dropbox/dphil/isolates/genotypes/",fname,".fam",sep=""),as.is=TRUE)
      fam.unrelated <- read.table(paste("~/Dropbox/dphil/isolates/genotypes/",fname,"-chr10-unrelated.fam",sep=""),as.is=TRUE)
      rrm <- matrix(scan(paste("~/Dropbox/dphil/isolates/rrm/",fname,".rrm.gz",sep="")),nrow=nrow(fam.all))
      diag(rrm) <- NA
      ii <- match(slrp[[fname]][,1],fam.all[,2])
      slrp[[fname]]$max <- r <- apply(rrm[ii,ii],1,max,na.rm=TRUE)
  }

#slrp.yield <- do.call("rbind",lapply(switch.error,function(x) c(x[[1]][1,4],x[[2]][1,4])))




pibd <- numeric(length(data))
pibd.unrelated <- numeric(length(data))

mibd <- numeric(length(data))
mibd.unrelated <- numeric(length(data))

data <- list()
for(i in 1:length(cohorts)) {
  f <- cohorts[i]
  print(f)
  data[[f]] <- list(ids=read.table(paste("~/Dropbox/dphil/isolates/genotypes/",f,"-chr10-unrelated.fam",sep=""),as.is=T)[,2])
  data[[f]]$rrm <- readrrm(paste("~/Dropbox/dphil/isolates/rrm/",f,".rrm.gz",sep=''),paste("~/Dropbox/dphil/isolates/genotypes/",f,".fam",sep=''),data[[f]]$ids)
  data[[f]]$n <- nrow(data[[f]]$rrm)
  bimfile <- paste("~/Dropbox/dphil/isolates/genotypes/",f,"-chr10-unrelated.bim",sep="")
  fname <- paste("~/Dropbox/dphil/isolates/ibd/",f,"-chr10-unrelated-slrp.ibd.gz",sep="")

  data[[f]]$chunks <- read.table(fname,as.is=T)$V3
  
  data[[f]]$r_ibd <- slrpIBDmatrix(fname,data[[f]]$ids)
  sumibd <- slrp.propibd(fname,data[[f]]$ids,bimfile,1.0,data[[f]]$rrm)
  sumibd.unrelated <- slrp.propibd(fname,data[[f]]$ids,bimfile,.125,data[[f]]$rrm)
  propibd <- sumibd >0
  propibd.unrelated <- sumibd.unrelated >0
  data[[f]]$pibd <- 100*rowMeans(propibd)
  data[[f]]$pibd.unrelated <- 100*rowMeans(propibd.unrelated)
  data[[f]]$mibd <- rowMeans(sumibd)
  data[[f]]$mibd.unrelated <- rowMeans(sumibd.unrelated)
  
}


for(i in 1:length(cohorts)) {
    f <- cohorts[i]
  fname <- paste("~/Dropbox/dphil/isolates/ibd/",f,"-chr10-unrelated-slrp.ibd.gz",sep="")
  data[[f]]$chunks <- read.table(fname,as.is=T)$V3
}


f1 <- function(d)   100*rowSums(d$ibd)
f2 <- function(d,thresh) {
  rrm <- d$rrm
  diag(rrm) <- 0
  ibd <- d$ibd
  ibd[rrm>thresh] <- 0
  100*rowSums(ibd)
}
f3 <- function(d) 100*rowMeans(d$probibd)
f4 <- function(d,thresh) {
  rrm <- d$rrm
  diag(rrm) <- 0

  100*rowMeans(d$probibd)
}


#nms <- c("CARL","ORKNEY","VALBORBERA","FVG","KORCULA","SPLIT","VIS")
nms <- c("CARL","Orkney","Valborbera","FVG","Korcula","Split","VIS")



slrp.yield <- 100*do.call("rbind",lapply(slrp,function(x) c(sum(x[,3])/sum(x[,2]),sum(x[x$max<.125,3])/sum(x[x$max<.125,2]))))



## pdf("percentage_chromosome_IBD.pdf",width=12,height=6)

## pibd <- lapply(data,function(x) x$pibd)
## par(mar=c(5,4,1,1)+.1,mfrow=c(1,2))
## ii <- order(unlist(lapply(data,function(x) mean(x$mibd.unrelated))))
## boxplot(pibd[ii],col='light blue',ylab="% Chromosome 10 that is IBD1",names=NA,axes=0,ylim=c(0,100),main="All individuals");axis(2);box()

## pibd.unrelated <- lapply(data,function(x) x$pibd.unrelated)
## text(x=1:length(cohorts),y=-6,nms[ii],srt=45,xpd=1,adj=1)
## boxplot(pibd.unrelated[ii],col='light blue',ylab="",names=NA,axes=0,ylim=c(0,100),main="Close relatives removed");axis(2);box()
## text(x=1:length(cohorts),y=-7,nms[ii],srt=45,xpd=1,adj=1)

## dev.off()

mibd <- lapply(data,function(d) d$mibd)
mibd.unrelated <- lapply(data,function(d) d$mibd.unrelated)

nms <- c("CARL","ORKNEY","VALBORBERA","FVG","KORCULA","SPLIT","VIS")

## pdf("meanibd_chromosome_IBD.pdf",width=12,height=6)

## par(mar=c(5,4,1,1)+.1,mfrow=c(1,2))
## ii <- order(unlist(lapply(mibd,mean)))
## yl <- range(mibd)
## boxplot(mibd[ii],col='light blue',ylab="Mean number of individuals that are IBD1",names=NA,axes=0,ylim=yl,main="All individuals");axis(2);box()
## text(x=1:length(cohorts),y=-1,nms[ii],srt=45,xpd=1,adj=1)
## boxplot(mibd.unrelated[ii],col='light blue',ylab="",names=NA,axes=0,ylim=yl,main="Close relatives removed");axis(2);box()
## text(x=1:length(cohorts),y=-1,nms[ii],srt=45,xpd=1,adj=1)

## dev.off()


#cohort_avg_totsharing <- lapply(data,function(x) 100*x$r_ibd[upper.tri(x$r_ibd)])

#cohort_avg_sharing.unrelated <- lapply(data,function(x) 100*colSums(x$r_ibd)[apply(x$rrm,1,max,na.rm=T)<.125] / (sum(apply(x$rrm,1,max,na.rm=T)<.125)-1))

pdf("meanibd_chromosome_IBD_unrelated.pdf",width=4,height=4)

par(mar=c(4,4,1,1)+.1,cex=1)
ii <- order(unlist(lapply(data,function(x) mean(x$mibd.unrelated))))
mibd.unrelated <- (lapply(data,function(x) x$mibd.unrelated)[ii])
yl <- range(mibd.unrelated)
boxplot(mibd.unrelated,col='light blue',ylab="Average #individuals IBD1",names=NA,axes=0,ylim=yl,main="");axis(2);box()
text(x=1:length(cohorts),y=-1,nms[ii],srt=45,xpd=1,adj=1)

dev.off()


## boxplot(lapply(data,function(d) d$rrm[d$rrm<.2]),col='light blue',ylab="Relatedness coefficient")

## reg1 <- c(55.50165,56.56438)*1e6
## reg2 <- c(88.18071,88.18071+diff(reg1))

## info <- subset(read.table("~/Dropbox/dphil/isolates/ibd/orkney889_b37-10-imputed_info.gz",as.is=TRUE,header=TRUE,comment.char=""),exp_freq_a1>0)

## info1 <- subset(info,position>reg1[1] & position<reg1[2] & snp_id=="---" & exp_freq_a1>0)
## info2 <- subset(info,position>reg2[1] & position<reg2[2] & snp_id=="---"& exp_freq_a1>0)

## plot(density(info1$info,from=0,to=1),main="",xlab="Impute2 INFO score",col=cols[1],lwd=2)
## lines(density(info2$info,from=0,to=1),col=cols[2],lwd=2)
## legend("topleft",legend=c("Low IBD","High IBD"),col=cols[1:2],lwd=2)
## map <- read.table("~/Dropbox/dphil/isolates/orkney889_b37-chr10-unrelated.bim",as.is=TRUE)
## ibd <- colMeans(orkney.propibd)
## info$ibd <- approx(map$V4,ibd,xout=info$position)$y

## est <- bkde2D(info[1000:(nrow(info)-1000),c("info","ibd")],bandwidth=c(0.01,.1),gridsize=c(256,256))
## image(est$x1,est$x2,est$fhat,col=heat.colors(1000))
## contour(est$x1,est$x2,est$fhat,add=TRUE)

## plot(density(info$info[info$ibd<=median(ibd)],from=0,to=1),main="",xlab="Impute2 INFO score",col=cols[1],lwd=2)
## lines(density(info$info[info$ibd>10],from=0,to=1),col=cols[1],lwd=2)



pdf("ibd_summary.pdf",height=4,width=8)
cols <- c("dark grey", "light grey")
cols <- brewer.pal(2,"Pastel1")[1:2]
  
par(mar=c(4.5,4.5,1,1)+.1,mfrow=c(1,2))
ii <- order(slrp.yield[,1])
xp <- colMeans(barplot(t(slrp.yield[ii,]),beside=T,col=cols,ylab='SLRP % heterozygotes phased',ylim=c(0,100),las=2,axisnames=F));
box();abline(h=axTicks(2),col='light grey')
text(x=xp,y=-2,nms[ii],srt=45,xpd=1,adj=1)
barplot(t(slrp.yield[ii,]),beside=T,col=cols,ylab='',ylim=c(0,100),legend.text=c('All founders','Distantly related'),
                args.legend=list(x="topleft",bg="white"),las=2,add=T,axisnames=F)
box()

## pairwise_sharing <- lapply(data,function(x) 100*x$r_ibd[upper.tri(x$r_ibd)])
## boxplot(pairwise_sharing[ii],col='light blue',ylab="Pairwise average sharing",names=NA,axes=0,main="",ylim=c(0,5));axis(2);box()
## text(x=1:length(cohorts),y=-.1,nms[ii],srt=45,xpd=1,adj=1)

## xhat <- seq(100*1/(2*183.11),50,.01)
## plot(NA,xlim=c(0,5),ylim=c(0.82,1),ylab="CDF",xlab="Pairwise IBD sharing (%)");grid()
## for(i in 1:length(nms)) lines(xhat,ecdf(pairwise_sharing[[i]])(xhat),col=cols[i],lwd=2)
## #for(i in 1:length(nms)) lines(density(pairwise_sharing[[i]],from=0,bw=.01),col=cols[i],lwd=3)
## legend("bottomright",legend=nms,col=cols,lwd=2,bg='white')

yl <- range(mibd)
boxplot(mibd.unrelated[ii],col=cols[2],ylab="Mean number of individuals that are IBD1",names=NA,axes=0,ylim=yl,main="");axis(2);box();
text(x=1:length(cohorts),y=-1,nms[ii],srt=45,xpd=1,adj=1)

## cohort_avg_sharing <- lapply(data,function(x) 100*colSums(x$r_ibd) / (ncol(x$r_ibd)-1))
## boxplot(cohort_avg_sharing[ii],col='light blue',ylab="Cohort average sharing",names=NA,axes=0,main="");axis(2);box()
## abline(h=axTicks(2),col='light grey')
## boxplot(cohort_avg_sharing[ii],col='light blue',ylab="Cohort average sharing",names=NA,axes=0,main="",add=T)
## text(x=1:length(cohorts),y=-.075,nms[ii],srt=45,xpd=1,adj=1)

dev.off()


## chunks <- lapply(data,function(x) x$chunks[x$chunks>0])
## boxplot(chunks[ii],col='light blue',ylab="IBD tract size",names=NA,axes=0,ylim=c(0,20),main="");axis(2);box()
## text(x=1:length(cohorts),y=-1,nms[ii],srt=45,xpd=1,adj=1)




plot(density(chunks[[1]],from=0),xlim=c(0,25),col=brewer.pal(8,"Set1")[1],lwd=2,ylim=c(0,.32))
for(i in 2:7) lines(density(chunks[[i]],from=0),xlim=c(0,25),col=brewer.pal(8,"Set1")[-6][i],lwd=2)
legend("topright",legend=nms[ii],col=brewer.pal(8,"Set1")[-6][1:7],lwd=2)
