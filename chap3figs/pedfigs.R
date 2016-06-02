map <- read.table("nuclearfam.map",as.is=T)
r <- read.table("nuclearfam.recombinations",as.is=T)
r <- r[!duplicated(r[,1:2]),]
match(r$V3,map$V4)


nsnp <- nrow(map)*2
n <- 5
tmp <- scan("nuclearfam.ped",what="character")
ii <- (rep(1:6,n) + rep(seq(0,(nsnp+6)*(n-1),6+nsnp),each=6))
fam <- as.data.frame(matrix(tmp[ii],ncol=6,nrow=n,byrow=T))

g0 <- matrix(tmp[-ii],ncol=n,nrow=nsnp)[,order(fam[,3])]
ref <- rep(g[seq(1,nsnp,2),1],each=2)
G <- (g0==ref)*1
fam <- fam[order(fam[,3]),]

i1 <- seq(1,nsnp,2)
i2 <- i1+1

H <- matrix(nrow=nsnp,ncol=2*n)
for(i in 0:(n-1)) {
  H[,i*2+1] <- G[i1,i+1]
  H[,i*2+2] <- G[i2,i+1] 
}
    
H[145:155,]


