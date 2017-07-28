fst<-c(0.05, 0.2)
nind<-200
nloci<-5000

### --- simulate ancestral allele frequencies
anc.pi<-rbeta(nloci, 0.8, 0.8)

### --- simulate cluster allele frequencies
p.k1.fst1<-data.frame(pop1=numeric(nloci))
p.k2.fst1<-data.frame(pop1=numeric(nloci), pop2=numeric(nloci))
p.k3.fst1<-data.frame(pop1=numeric(nloci), pop2=numeric(nloci), pop3=numeric(nloci))

for(i in 1:nloci){
  p.k1.fst1[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))

  p.k2.fst1[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))
  p.k2.fst1[i,2]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))

  p.k3.fst1[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))
  p.k3.fst1[i,2]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))
  p.k3.fst1[i,3]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))
}

### --- simulate cluster allele frequencies
p.k1.fst2<-data.frame(pop1=numeric(nloci))
p.k2.fst2<-data.frame(pop1=numeric(nloci), pop2=numeric(nloci))
p.k3.fst2<-data.frame(pop1=numeric(nloci), pop2=numeric(nloci), pop3=numeric(nloci))

for(i in 1:nloci){
  p.k1.fst2[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))

  p.k2.fst2[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))
  p.k2.fst2[i,2]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))

  p.k3.fst2[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))
  p.k3.fst2[i,2]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))
  p.k3.fst2[i,3]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))
}

### --- simulate ancestry combinations (Q matrix lower triangle plus
### diagonal) for sampled individuals, choose from different hybrid
### classes.  Want to model different ancestry classes, not just
### admixture coefficients

### these will parameters for a multinomial draw for locus specific ancestries

Q.k1 <- rep(1, nind)
Q.k2 <- rbind(
  matrix(c(1,0,0), ncol=3, nrow=50, byrow=T),  ## 50 parentals type 1
  matrix(c(0,0,1), ncol=3, nrow=50, byrow=T),  ## 50 parentals type 2
  matrix(c(0,1,0), ncol=3, nrow=20, byrow=T),  ## 20 F1s
  matrix(c(0.25,0.5,0.25), ncol=3, nrow=20, byrow=T),  ## 20 F2s
  matrix(c(0.5,0.5,0), ncol=3, nrow=20, byrow=T),  ## 20 BC1s to type 1
  matrix(c(0,0.5,0.5), ncol=3, nrow=20, byrow=T),  ## 20 BC1s to type 2
  matrix(c(0.375,0.25,0.375), ncol=3, nrow=20, byrow=T)  ## 20 F3s
  )
Q.k3 <- rbind(  ### NB: there are no three way hybrids
              matrix(c(1,0,0,0,0,0), ncol=6, nrow=30, byrow=T),  ## 30 parentals type 1
              matrix(c(0,0,1,0,0,0), ncol=6, nrow=30, byrow=T),  ## 30 parentals type 2
              matrix(c(0,0,0,0,0,1), ncol=6, nrow=30, byrow=T),  ## 30 parentals type 3
              matrix(c(0,1,0,0,0,0), ncol=6, nrow=10, byrow=T),  ## 10 F1s of 1/2
              matrix(c(0,0,0,1,0,0), ncol=6, nrow=10, byrow=T),  ## 10 F1s of 1/3
              matrix(c(0,0,0,0,1,0), ncol=6, nrow=10, byrow=T),  ## 10 F1s of 2/3
              matrix(c(0.25,0.5,0.25,0,0,0), ncol=6, nrow=10, byrow=T),  ## 10 F2s between type 1 and 2
              matrix(c(0.25,0,0,0.5,0,0.25), ncol=6, nrow=10, byrow=T),  ## 10 F2s between type 1 and 3
              matrix(c(0,0,0.25,0,0.5,0.25), ncol=6, nrow=10, byrow=T),  ## 10 F2s between type 2 and 3
              matrix(c(0.5,0.5,0,0,0,0), ncol=6, nrow=10, byrow=T),  ## 10 BC1s between 1/2 F1 and type 1
              matrix(c(0.5,0,0,0.5,0,0), ncol=6, nrow=10, byrow=T),  ## 10 BC1s between 1/3 F1 and type 1
              matrix(c(0,0,0.5,0,0.5,0), ncol=6, nrow=10, byrow=T),  ## 10 BC1s between 2/3 F1 and type 2
              matrix(c(0.375,0.25,0.375,0,0,0), ncol=6, nrow=10, byrow=T),  ## 10 F3s of 1/2
              matrix(c(0,0,0,0.375,0.25,0.375), ncol=6, nrow=10, byrow=T)   ## 10 F3s of 2/3
              )

### sample locus specific ancestries for all individuals from a
### multinomial with Q parameters ... translate these into allele
### copies so we can do Bernoulli trials to build up genotypes.  Store
### only genotypes for now


## k1
g.k1.fst1<-matrix(0, nrow=nind, ncol=nloci)
for ( locus in 1:nloci ){
  g.k1.fst1[,locus]<-rbinom(nind, 2, p.k1.fst1$pop1[locus])
}
g.k1.fst2<-matrix(0, nrow=nind, ncol=nloci)
for ( locus in 1:nloci ){
  g.k1.fst2[,locus]<-rbinom(nind, 2, p.k1.fst2$pop1[locus])
}

## k2
z<-numeric(nloci)
g.k2.fst1<-matrix(0, nrow=nind, ncol=nloci)
for (ind in 1:nrow(Q.k2)){
  z<-rmultinom(nloci, 1, Q.k2[ind,])  ## ancestry vector
  for ( locus in 1:nloci ){
    zz<-which(z[,locus] == 1)  ## ancestry for two allele copies in an ind
    if(zz == 1){
     g.k2.fst1[ind, locus] <- rbinom(1,2, prob=p.k2.fst1$pop1[locus]) 
    }
    else if (zz == 2){
      g.k2.fst1[ind, locus] <- rbinom(1,1, prob=p.k2.fst1$pop1[locus]) +
        rbinom(1,1, prob=p.k2.fst1$pop2[locus])
    }
    else{
      g.k2.fst1[ind, locus] <- rbinom(1,2, prob=p.k2.fst1$pop2[locus])
    }
  }
}

g.k2.fst2<-matrix(0, nrow=nind, ncol=nloci)
for (ind in 1:nrow(Q.k2)){
  z<-rmultinom(nloci, 1, Q.k2[ind,])
  for ( locus in 1:nloci ){
    zz<-which(z[,locus] == 1)
    if(zz == 1){
     g.k2.fst2[ind, locus] <- rbinom(1,2, prob=p.k2.fst2$pop1[locus]) 
    }
    else if (zz == 2){
      g.k2.fst2[ind, locus] <- rbinom(1,1, prob=p.k2.fst2$pop1[locus]) +
        rbinom(1,1, prob=p.k2.fst2$pop2[locus])
    }
    else{
      g.k2.fst2[ind, locus] <- rbinom(1,2, prob=p.k2.fst2$pop2[locus]) 
    }
  }
}



#k3
g.k3.fst1<-matrix(0, nrow=nind, ncol=nloci)
for (ind in 1:nrow(Q.k3)){
  z<-rmultinom(nloci, 1, Q.k3[ind,])  ## ancestry vector
  for ( locus in 1:nloci ){
    zz<-which(z[,locus] == 1)  ## ancestry for two allele copies in an ind
    if(zz == 1){
     g.k3.fst1[ind, locus] <- rbinom(1,2, prob=p.k3.fst1$pop1[locus]) 
    }
    else if (zz == 2){
      g.k3.fst1[ind, locus] <- rbinom(1,1, prob=p.k3.fst1$pop1[locus]) +
        rbinom(1,1, prob=p.k3.fst1$pop2[locus])
    }
    else if (zz == 3){
      g.k3.fst1[ind, locus] <- rbinom(1,2, prob=p.k3.fst1$pop2[locus]) 
    }
    else if (zz == 4){
      g.k3.fst1[ind, locus] <- rbinom(1,1, prob=p.k3.fst1$pop1[locus]) +
        rbinom(1,1, prob=p.k3.fst1$pop3[locus])
    }
    else if (zz == 5){
      g.k3.fst1[ind, locus] <- rbinom(1,1, prob=p.k3.fst1$pop2[locus]) +
        rbinom(1,1, prob=p.k3.fst1$pop3[locus])
    }
    else {
      g.k3.fst1[ind, locus] <- rbinom(1,2, prob=p.k3.fst1$pop3[locus])
    }
  }
}

g.k3.fst2<-matrix(0, nrow=nind, ncol=nloci)
for (ind in 1:nrow(Q.k3)){
  z<-rmultinom(nloci, 1, Q.k3[ind,])  ## ancestry vector
  for ( locus in 1:nloci ){
    zz<-which(z[,locus] == 1)  ## ancestry for two allele copies in an ind
    if(zz == 1){
     g.k3.fst2[ind, locus] <- rbinom(1,2, prob=p.k3.fst2$pop1[locus]) 
    }
    else if (zz == 2){
      g.k3.fst2[ind, locus] <- rbinom(1,1, prob=p.k3.fst2$pop1[locus]) +
        rbinom(1,1, prob=p.k3.fst2$pop2[locus])
    }
    else if (zz == 3){
      g.k3.fst2[ind, locus] <- rbinom(1,2, prob=p.k3.fst2$pop2[locus]) 
    }
    else if (zz == 4){
      g.k3.fst2[ind, locus] <- rbinom(1,1, prob=p.k3.fst2$pop1[locus]) +
        rbinom(1,1, prob=p.k3.fst2$pop3[locus])
    }
    else if (zz == 5){
      g.k3.fst2[ind, locus] <- rbinom(1,1, prob=p.k3.fst2$pop2[locus]) +
        rbinom(1,1, prob=p.k3.fst2$pop3[locus])
    }
    else {
      g.k3.fst2[ind, locus] <- rbinom(1,2, prob=p.k3.fst2$pop3[locus])
    }
  }
}


### write input files for entropy ... do not model variance in coverage for now
### 
### a) SNP counts    ... model 10x coverage to begin
### b) genotype likelihoods  (make sure we get the same as for SNP counts)

## make SNP count matrix from genotype and append to file
for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k1.fst1.txt", append=T)
  write.table(make.snpcount(g.k1.fst1[,locus]), file="k1.fst1.txt", append=T,
              col.names=F, row.names=F)
}
for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k1.fst2.txt", append=T)
  write.table(make.snpcount(g.k1.fst2[,locus]), file="k1.fst2.txt", append=T,
              col.names=F, row.names=F)
}

for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k2.fst1.txt", append=T)
  write.table(make.snpcount(g.k2.fst1[,locus]), file="k2.fst1.txt", append=T,
              col.names=F, row.names=F)
}
for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k2.fst2.txt", append=T)
  write.table(make.snpcount(g.k2.fst2[,locus]), file="k2.fst2.txt", append=T,
              col.names=F, row.names=F)
}

for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k3.fst1.txt", append=T)
  write.table(make.snpcount(g.k3.fst1[,locus]), file="k3.fst1.txt", append=T,
              col.names=F, row.names=F)
}
for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k3.fst2.txt", append=T)
  write.table(make.snpcount(g.k3.fst2[,locus]), file="k3.fst2.txt", append=T,
              col.names=F, row.names=F)
}


## How ancestry informative are the markers?
## plot(p.k2.fst2[,1], p.k2.fst2[,2])
## hist(abs(p.k2.fst2[,1] - p.k2.fst2[,2]))

### how informative are the markers about Fst?
## hist(abs(p.k1.fst1[,1] - anc.pi))

make.snpcount<-function(x){
  y<-matrix(0, nrow=length(x), 2)
  for(i in 1:length(x)){
    if(x[i]==2){
      y[i,]<-c(0,10)
    }
    else if (x[i]==1){
      y[i,]<-c(5,5)
    }
    else{
      y[i,]<-c(10,0)
    }
  }
  return(y)
}
    
