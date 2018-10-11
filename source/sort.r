rm(list=ls())
setwd("/Users/vinicius/opt/genomic_module/source/")


(ped <- read.table('unsort.ped'))

n = nrow(ped)
pedr = matrix(data = 0, nrow = n, ncol = 3)
gen = tgen = perm = iperm = rep(1,n)
miss=1

## Gen count
while(TRUE){
  flag1 = sum(gen)
  for(i in 1:n){
    indv = ped[i,1]
    sire = ped[i,2]
    dam  = ped[i,3]

    # sire
    if (sire > 0) {
      if( gen[indv] >= gen[sire] ) gen[sire]=gen[sire]+1
    }

    # dam
    if (dam > 0) {
      if( gen[indv] >= gen[dam] ) gen[dam]=gen[dam]+1
    }
  }
  
  flag2=sum(gen)
  if (flag1==flag2) break
  miss=miss+1
  if(miss>50) stop("circular pedigree!!")
}
for(i in 1:n){
  tgen[i]=gen[ped[i,1]]
}

idx = order(tgen, decreasing = T)
ped = ped[idx,]

## recod
for(i in 1:n){
  id=ped[i,1]
  iperm[i]=id
  perm[id]=i
}

