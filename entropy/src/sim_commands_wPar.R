### include family data (offspring; only pure, F1, BC1 or F2)
## for k=2 only

fst<-c(0.05, 0.2)
nind<-200
nloci<-5000

nchr<-20
start.chr<-seq(1,nloci,nloci/nchr)
end.chr<-seq(nloci/nchr,nloci,nloci/nchr)
locusids<- paste(rep(1:nchr,each=nloci/nchr),rep(1:(nloci/nchr),times=nchr),sep=":")


### --- simulate ancestral allele frequencies
anc.pi<-rbeta(nloci, 0.8, 0.8)

### --- simulate cluster allele frequencies
p.k2.fst1<-data.frame(pop1=numeric(nloci), pop2=numeric(nloci))

for(i in 1:nloci){
  p.k2.fst1[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))
  p.k2.fst1[i,2]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[1]), (1-anc.pi[i]) * (-1 + 1/fst[1]))
}

### --- simulate cluster allele frequencies
p.k2.fst2<-data.frame(pop1=numeric(nloci), pop2=numeric(nloci))

for(i in 1:nloci){
  p.k2.fst2[i,1]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))
  p.k2.fst2[i,2]<-rbeta(1,anc.pi[i] * (-1 + 1/fst[2]), (1-anc.pi[i]) * (-1 + 1/fst[2]))
}

### --- simulate ancestry combinations (Q matrix lower triangle plus
### diagonal) for sampled individuals, choose from different hybrid
### classes.  Want to model different ancestry classes, not just
### admixture coefficients

### these will parameters for a multinomial draw for locus specific ancestries

Q.k2 <- rbind(
  matrix(c(1,0,0), ncol=3, nrow=50, byrow=T),  ## 50 parentals type 1
  matrix(c(0,0,1), ncol=3, nrow=50, byrow=T),  ## 50 parentals type 2
  matrix(c(0,1,0), ncol=3, nrow=20, byrow=T),  ## 20 F1s
  matrix(c(0.25,0.5,0.25), ncol=3, nrow=20, byrow=T),  ## 20 F2s
  matrix(c(0.5,0.5,0), ncol=3, nrow=20, byrow=T),  ## 20 BC1s to type 1
  matrix(c(0,0.5,0.5), ncol=3, nrow=20, byrow=T),  ## 20 BC1s to type 2
  matrix(c(0.375,0.25,0.375), ncol=3, nrow=20, byrow=T)  ## 20 F3s
  )

### sample locus specific ancestries for all individuals from a
### multinomial with Q parameters ... translate these into allele
### copies so we can do Bernoulli trials to build up genotypes.  Store
### only genotypes for now


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


## add parents and offspring to k2, for fst1
##--------------------------------------------------------------------------
## make parental ancestries
moms.z1<- rbind(
    matrix(0,nrow=5,ncol=nloci), # 5 parental type 1
    matrix(0,nrow=5,ncol=nloci), # 5 F1
    matrix(1,nrow=5,ncol=nloci) # 5 parental type 2
    )
moms.z2<- rbind(
    matrix(0,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci)
    )

dads.z1<- rbind(
    matrix(0,nrow=5,ncol=nloci), # 5 parental type 1
    matrix(0,nrow=5,ncol=nloci), # 5 F1
    matrix(1,nrow=5,ncol=nloci) # 5 parental type 2
    )
dads.z2<- rbind(
    matrix(0,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci)
    )

## sample parental genotypes
nmoms<-dim(moms.z1)[1]
ndads<-dim(dads.z1)[1]

moms.g<-array(dim=c(nmoms,nloci,2),dimnames=c("mom","locus","copy"))
dads.g<-array(dim=c(ndads,nloci,2),dimnames=c("dad","locus","copy"))

for (mom in 1:nmoms){
    for ( locus in 1:nloci ){
        if(moms.z1[mom,locus] == 0){
            moms.g[mom, locus,1] <- rbinom(1,1, prob=p.k2.fst1$pop1[locus]) 
        }
        else if (moms.z1[mom,locus] == 1){
            moms.g[mom, locus,1] <- rbinom(1,1, prob=p.k2.fst1$pop2[locus]) 
        }

        if(moms.z2[mom,locus] == 0){
            moms.g[mom, locus,2] <- rbinom(1,1, prob=p.k2.fst1$pop1[locus]) 
        }
        else if (moms.z2[mom,locus] == 1){
            moms.g[mom, locus,2] <- rbinom(1,1, prob=p.k2.fst1$pop2[locus]) 
        }
    }
}

for (dad in 1:ndads){
    for ( locus in 1:nloci ){
        if(dads.z1[dad,locus] == 0){
            dads.g[dad, locus,1] <- rbinom(1,1, prob=p.k2.fst1$pop1[locus]) 
        }
        else if (dads.z1[dad,locus] == 1){
            dads.g[dad, locus,1] <- rbinom(1,1, prob=p.k2.fst1$pop2[locus]) 
        }

        if(dads.z2[dad,locus] == 0){
            dads.g[dad, locus,2] <- rbinom(1,1, prob=p.k2.fst1$pop1[locus]) 
        }
        else if (dads.z2[dad,locus] == 1){
            dads.g[dad, locus,2] <- rbinom(1,1, prob=p.k2.fst1$pop2[locus]) 
        }
    }
}

## make gametes
## recombi: 1 crossover per chromosome, fixed
noffperfam<-ndads
numbering<- c(paste("0",1:9,sep=""),10:max(nmoms,ndads))

moms.gamete<-matrix(NA,nrow=nmoms*noffperfam,ncol=nloci)
offids<-character(nmoms*noffperfam)

for (mom in 1:nmoms){
    for (fam in 1:noffperfam){
        chrindex<-numeric(nloci)
        brake<-start.chr + floor(runif(n=nchr,min=0,max=(nloci/nchr)-1))

        for (chr in 1:nchr){
            copy<-sample(c(1,2),1)
            if (copy == 1) {
                chrindex[start.chr[chr]:brake[chr]]<- 1
                chrindex[(brake[chr]+1):end.chr[chr]]<- 2
            }
            else if (copy == 2){
                chrindex[start.chr[chr]:brake[chr]]<- 2
                chrindex[(brake[chr]+1):end.chr[chr]]<- 1
            }
        }


        moms.gamete[noffperfam*(mom-1)+fam,which(chrindex==1)] <- moms.g[mom,which(chrindex==1),1]
        moms.gamete[noffperfam*(mom-1)+fam,which(chrindex==2)] <- moms.g[mom,which(chrindex==2),2]
        offids[noffperfam*(mom-1)+fam]<- paste("off",numbering[mom],numbering[fam],sep=".")
    }
}

dads.gamete<-matrix(NA,nrow=ndads*noffperfam,ncol=nloci)
offids2<-character(ndads*noffperfam)

for (fam in 1:noffperfam){
    for (dad in 1:ndads){
        chrindex<-numeric(nloci)
        brake<-start.chr + floor(runif(n=nchr,min=0,max=(nloci/nchr)-1))

        for (chr in 1:nchr){
            copy<-sample(c(1,2),1)
            if (copy == 1) {
                chrindex[start.chr[chr]:brake[chr]]<- 1
                chrindex[(brake[chr]+1):end.chr[chr]]<- 2
            }
            else if (copy == 2){
                chrindex[start.chr[chr]:brake[chr]]<- 2
                chrindex[(brake[chr]+1):end.chr[chr]]<- 1
            }
        }


        dads.gamete[noffperfam*(fam-1)+dad,which(chrindex==1)] <- dads.g[dad,which(chrindex==1),1]
        dads.gamete[noffperfam*(fam-1)+dad,which(chrindex==2)] <- dads.g[dad,which(chrindex==2),2]
        offids2[noffperfam*(fam-1)+dad]<- paste("off",numbering[fam],numbering[dad],sep=".")
    }
}

### make offspring
## (each mom mates once with each of the dads)
g.offspring<- moms.gamete + dads.gamete
rownames(g.offspring)<- offids
g.mom<-moms.g[,,1] + moms.g[,,2]
rownames(g.mom)<- paste("mom",numbering[1:nmoms],sep=".")
g.dad<-dads.g[,,1] + dads.g[,,2]
rownames(g.dad)<- paste("dad",numbering[1:ndads],sep=".")


## bind to population reference samples
g.wPar.k2.fst1<- rbind(g.k2.fst1,g.mom,g.dad,g.offspring)



###############################################################################
### write input files for entropy ... do not model variance in coverage for now
### 
### a) SNP counts    ... model 10x coverage to begin
### b) genotype likelihoods  (make sure we get the same as for SNP counts)

## a) make SNP count matrix from genotype and append to file

for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k2.wPar.fst1.txt", append=T)
  write.table(make.snpcount(g.wPar.k2.fst1[,locus]), file="k2.wPar.fst1.txt", append=T,
              col.names=F, row.names=F)
}

## b) make genotype likelihoods from genotype data
for(locus in 1:nloci){
  write.table(t(make.genolikes(g.wPar.k2.fst1[,locus],locusids[locus])), file="k2.wPar.fst1.like.txt", append=T,
               col.names=F, row.names=F, quote=F)
}
  
    

## make header file (add number of loci and inds afterwards) dim(g.wPar.k2.fst1); length(offids)
refids<-paste("ref",c(paste("00",1:9,sep=""),paste("0",10:99,sep=""),100:nind),sep=".")
ids<- c(refids,rownames(g.mom),rownames(g.dad),rownames(g.offspring))

classids<- c(rep("r00",nind),paste("m",numbering[1:nmoms],sep=""),rep("r00",ndads),rep(paste("p",numbering[1:nmoms],sep=""),each=noffperfam))

write.table(rbind(ids,classids),file="header.txt",quote=F,row.names=F,col.names=F)



## write hybrid class for each ind
refclass<- rep(c("par1","par2","F1","F2","BC1.1","BC1.2","F3"),times=c(50,50,20,20,20,20,20)) # reference samples
parclass<- c(rep(c("par1","F1","par2"),each=5),rep(c("par1","F1","par2"),each=5)) # parentals
offclass<- c(rep(rep(c("par1","BC1.1","F1"),each=5),times=5), # offspring
             rep(rep(c("BC1.1","F2","BC1.2"),each=5),times=5),
             rep(rep(c("F1","BC1.2","par2"),each=5),times=5))

hybridclasses<- c(refclass, parclass, offclass)

write.table(hybridclasses,file="hybridclasses.txt",quote=F,row.names=F,col.names=F)


## How ancestry informative are the markers?
## plot(p.k2.fst1[,1], p.k2.fst1[,2])
## hist(abs(p.k2.fst1[,1] - p.k2.fst1[,2]))



## c) for low coverage
#---#---#---#---#---#

## sample reads according to coverage
mycoverage<- 3

draw.reads<- function(genotype,coverage=mycoverage){
    rbinom(1,coverage,(genotype/2))
}

observed<- apply(g.wPar.k2.fst1,c(1,2),draw.reads)

# make SNP countmatrix
countmat<- cbind(as.vector(rbind(rep("locus",nloci),observed)),
                 as.vector(rbind(1:nloci,mycoverage-observed)))

# genotype likelihoods P(X | genotype, nreads)
all0<- c(1,dbinom(0,mycoverage,0.5),0.00001)
all1<- c(0.00001,dbinom(0,mycoverage,0.5),1)
all0<- round(-10*log10(all0/sum(all0)),digits=0)
all1<- round(-10*log10(all1/sum(all1)),digits=0)
het<- c(0.025,0.95,0.025)
het<- round(-10*log10(het),digits=0)


all0<- paste(all0[1],all0[2],all0[3],sep=" ")
all1<- paste(all1[1],all1[2],all1[3],sep=" ")
het<- paste(het[1],het[2],het[3],sep=" ")


calc.like<- function(observed.counts,coverage=mycoverage){
    if (observed.counts==0) {genolike<- all0}
    else if (observed.counts==coverage) {genolike<- all1}
    else {genolike<- het}
    return(genolike)
}

likelimat<-apply(observed,c(1,2),calc.like)

likelimat<- cbind(locusids,t(likelimat))



write.table(countmat, file="k2.wPar.fst1.cov3.txt", quote=F, col.names=F, row.names=F)
write.table(likelimat, file="k2.wPar.fst1.cov3.like.txt", col.names=F, row.names=F, quote=F)
#---#---#---#---#---#








## add parents and offspring to k2, for fst2
##--------------------------------------------------------------------------
## make parental ancestries
moms.z1<- rbind(
    matrix(0,nrow=5,ncol=nloci), # 5 parental type 1
    matrix(0,nrow=5,ncol=nloci), # 5 F1
    matrix(1,nrow=5,ncol=nloci) # 5 parental type 2
    )
moms.z2<- rbind(
    matrix(0,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci)
    )

dads.z1<- rbind(
    matrix(0,nrow=5,ncol=nloci), # 5 parental type 1
    matrix(0,nrow=5,ncol=nloci), # 5 F1
    matrix(1,nrow=5,ncol=nloci) # 5 parental type 2
    )
dads.z2<- rbind(
    matrix(0,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci),
    matrix(1,nrow=5,ncol=nloci)
    )

## sample parental genotypes
nmoms<-dim(moms.z1)[1]
ndads<-dim(dads.z1)[1]

moms.g<-array(dim=c(nmoms,nloci,2),dimnames=c("mom","locus","copy"))
dads.g<-array(dim=c(ndads,nloci,2),dimnames=c("dad","locus","copy"))

for (mom in 1:nmoms){
    for ( locus in 1:nloci ){
        if(moms.z1[mom,locus] == 0){
            moms.g[mom, locus,1] <- rbinom(1,1, prob=p.k2.fst2$pop1[locus]) 
        }
        else if (moms.z1[mom,locus] == 1){
            moms.g[mom, locus,1] <- rbinom(1,1, prob=p.k2.fst2$pop2[locus]) 
        }

        if(moms.z2[mom,locus] == 0){
            moms.g[mom, locus,2] <- rbinom(1,1, prob=p.k2.fst2$pop1[locus]) 
        }
        else if (moms.z2[mom,locus] == 1){
            moms.g[mom, locus,2] <- rbinom(1,1, prob=p.k2.fst2$pop2[locus]) 
        }
    }
}

for (dad in 1:ndads){
    for ( locus in 1:nloci ){
        if(dads.z1[dad,locus] == 0){
            dads.g[dad, locus,1] <- rbinom(1,1, prob=p.k2.fst2$pop1[locus]) 
        }
        else if (dads.z1[dad,locus] == 1){
            dads.g[dad, locus,1] <- rbinom(1,1, prob=p.k2.fst2$pop2[locus]) 
        }

        if(dads.z2[dad,locus] == 0){
            dads.g[dad, locus,2] <- rbinom(1,1, prob=p.k2.fst2$pop1[locus]) 
        }
        else if (dads.z2[dad,locus] == 1){
            dads.g[dad, locus,2] <- rbinom(1,1, prob=p.k2.fst2$pop2[locus]) 
        }
    }
}

## make gametes
## recombi: 1 crossover per chromosome, fixed
noffperfam<-ndads
numbering<- c(paste("0",1:9,sep=""),10:max(nmoms,ndads))

moms.gamete<-matrix(NA,nrow=nmoms*noffperfam,ncol=nloci)
offids<-character(nmoms*noffperfam)

for (mom in 1:nmoms){
    for (fam in 1:noffperfam){
        chrindex<-numeric(nloci)
        brake<-start.chr + floor(runif(n=nchr,min=0,max=(nloci/nchr)-1))

        for (chr in 1:nchr){
            copy<-sample(c(1,2),1)
            if (copy == 1) {
                chrindex[start.chr[chr]:brake[chr]]<- 1
                chrindex[(brake[chr]+1):end.chr[chr]]<- 2
            }
            else if (copy == 2){
                chrindex[start.chr[chr]:brake[chr]]<- 2
                chrindex[(brake[chr]+1):end.chr[chr]]<- 1
            }
        }


        moms.gamete[noffperfam*(mom-1)+fam,which(chrindex==1)] <- moms.g[mom,which(chrindex==1),1]
        moms.gamete[noffperfam*(mom-1)+fam,which(chrindex==2)] <- moms.g[mom,which(chrindex==2),2]
        offids[noffperfam*(mom-1)+fam]<- paste("off",numbering[mom],numbering[fam],sep=".")
    }
}

dads.gamete<-matrix(NA,nrow=ndads*noffperfam,ncol=nloci)
offids2<-character(ndads*noffperfam)

for (fam in 1:noffperfam){
    for (dad in 1:ndads){
        chrindex<-numeric(nloci)
        brake<-start.chr + floor(runif(n=nchr,min=0,max=(nloci/nchr)-1))

        for (chr in 1:nchr){
            copy<-sample(c(1,2),1)
            if (copy == 1) {
                chrindex[start.chr[chr]:brake[chr]]<- 1
                chrindex[(brake[chr]+1):end.chr[chr]]<- 2
            }
            else if (copy == 2){
                chrindex[start.chr[chr]:brake[chr]]<- 2
                chrindex[(brake[chr]+1):end.chr[chr]]<- 1
            }
        }


        dads.gamete[noffperfam*(fam-1)+dad,which(chrindex==1)] <- dads.g[dad,which(chrindex==1),1]
        dads.gamete[noffperfam*(fam-1)+dad,which(chrindex==2)] <- dads.g[dad,which(chrindex==2),2]
        offids2[noffperfam*(fam-1)+dad]<- paste("off",numbering[fam],numbering[dad],sep=".")
    }
}

### make offspring
## (each mom mates once with each of the dads)
g.offspring<- moms.gamete + dads.gamete
rownames(g.offspring)<- offids
g.mom<-moms.g[,,1] + moms.g[,,2]
rownames(g.mom)<- paste("mom",numbering[1:nmoms],sep=".")
g.dad<-dads.g[,,1] + dads.g[,,2]
rownames(g.dad)<- paste("dad",numbering[1:ndads],sep=".")


## bind to population reference samples
g.wPar.k2.fst2<- rbind(g.k2.fst2,g.mom,g.dad,g.offspring)



###############################################################################
### write input files for entropy ... do not model variance in coverage for now
### 
### a) SNP counts    ... model 10x coverage to begin
### b) genotype likelihoods  (make sure we get the same as for SNP counts)

## a) make SNP count matrix from genotype and append to file

for(locus in 1:nloci){
  cat(paste("locus", locus, "\n"), file="k2.wPar.fst2.txt", append=T)
  write.table(make.snpcount(g.wPar.k2.fst2[,locus]), file="k2.wPar.fst2.txt", append=T,
              col.names=F, row.names=F)
}

## b) make genotype likelihoods from genotype data
for(locus in 1:nloci){
  write.table(t(make.genolikes(g.wPar.k2.fst2[,locus],locusids[locus])), file="k2.wPar.fst2.like.txt", append=T,
               col.names=F, row.names=F, quote=F)
}
  
    
## How ancestry informative are the markers?
## plot(p.k2.fst2[,1], p.k2.fst2[,2])
## hist(abs(p.k2.fst2[,1] - p.k2.fst2[,2]))

## c) for low coverage
#---#---#---#---#---#

## sample reads according to coverage
mycoverage<- 3

draw.reads<- function(genotype,coverage=mycoverage){
    rbinom(1,coverage,(genotype/2))
}

observed<- apply(g.wPar.k2.fst2,c(1,2),draw.reads)

# make SNP countmatrix
countmat<- cbind(as.vector(rbind(rep("locus",nloci),observed)),
                 as.vector(rbind(1:nloci,mycoverage-observed)))

# genotype likelihoods
all0<- c(1,dbinom(0,mycoverage,0.5),0.00001)
all1<- c(0.00001,dbinom(0,mycoverage,0.5),1)
all0<- round(-10*log10(all0/sum(all0)),digits=0)
all1<- round(-10*log10(all1/sum(all1)),digits=0)
het<- c(0.025,0.95,0.025)
het<- round(-10*log10(het),digits=0)


all0<- paste(all0[1],all0[2],all0[3],sep=" ")
all1<- paste(all1[1],all1[2],all1[3],sep=" ")
het<- paste(het[1],het[2],het[3],sep=" ")


calc.like<- function(observed.counts,coverage=mycoverage){
    if (observed.counts==0) {genolike<- all0}
    else if (observed.counts==coverage) {genolike<- all1}
    else {genolike<- het}
    return(genolike)
}

likelimat<-apply(observed,c(1,2),calc.like)

likelimat<- cbind(locusids,t(likelimat))



write.table(countmat, file="k2.wPar.fst2.cov3.txt", quote=F, col.names=F, row.names=F)
write.table(likelimat, file="k2.wPar.fst2.cov3.like.txt", col.names=F, row.names=F, quote=F)
#---#---#---#---#---#



########################################################################################
## functions ###
#######################################################################################

## a)
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


## b)
##### genotype likelihoods #####
## (makes only some sense for relatively high coverage)

## some set ups
coverage<- 10

## probability of observing only 1 allele when heterozygous 
phom<- dbinom(0,coverage,0.5)
phom<- round(-10*log10(phom)) # -10*log10(p-value) output by bcftools
phet<- 1 - (2* dbinom(0,coverage,0.5))
phet<- round(-10*log10(phet))


make.genolikes<- function(x,locusname="locus"){
    y<- numeric(length(x)*3 + 1)
    pos1<- seq(2,length(y),3) # position for genotype prob allele1/allele1
    pos2<- seq(3,length(y),3) # position for genotype prob allele1/allele2
    pos3<- seq(4,length(y),3) # position for genotype prob allele2/allele2

    y[1]<- locusname
    for (i in 1:length(x)){
        if(x[i]==2){
            y[pos1[i]]<- round(-10*log10(0.00001)) # arbitrary
            y[pos2[i]]<- phom # if all reads are type 1 there is a chance that ind is heterozygous with only one allele state observed    
            y[pos3[i]]<- phet # prob that ind is not a heterozygote if all reads are type 1
        }
        else if (x[i]==1){
            y[pos1[i]]<- round(-10*log10(0.005))  # arbitrary
            y[pos2[i]]<- round(-10*log10(0.99))   # arbitrary 
            y[pos3[i]]<- round(-10*log10(0.005))  # arbitrary     
        }
        else{
            y[pos1[i]]<- phet
            y[pos2[i]]<- phom       
            y[pos3[i]]<- round(-10*log10(0.00001))          
        }
    }
    return(y)
}





