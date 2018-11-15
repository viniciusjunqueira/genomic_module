rm(list=ls())
setwd("/Users/vinicius/opt/genomic_module/source/")

#
# Loading unsorted pedigree
#
(ped <- read.table('unsort.ped'))

#
# Initiatlizating needed variables
#
n=nrow(ped)
gen=tgen=perm=iperm=rep(1,n)
miss=1

#
# Step that computes 'generic' indices for each animal based on progeny. This is useful
# to sort individuals
#
while(TRUE){
  flag1 = sum(gen)
  for(i in 1:n){
    indv = ped[i,1]
    sire = ped[i,2]
    dam  = ped[i,3]
    # sire
    if (sire > 0){
      if( gen[indv] >= gen[sire] ) gen[sire]=gen[sire]+1
    }
    # dam
    if (dam > 0){
      if( gen[indv] >= gen[dam] ) gen[dam]=gen[dam]+1
    }
  }
  flag2=sum(gen)
  if (flag1==flag2) break
  miss=miss+1
  if(miss>50) stop(" circular pedigree!! ")
}
tgen=gen[ped[,1]]

#
# Organize indices to be in needed order
#
bool=rep(F,n)
Beg=max(tgen)
curr=1
indices=rep(0,n)
while(TRUE){
  if(all(bool)) break
  idx = which(tgen == Beg)
  bool[idx] = TRUE
  indices[idx] = curr:(curr+length(idx)-1)
  curr=max(indices)+1
  Beg=Beg-1
}
indices=order(tgen, decreasing = TRUE)
ped=ped[indices,]

#
# Recode pedigree based on pre-defined order
#
for(i in 1:n){
  id=ped[i,1]
  iperm[i]=id # Original identification
  perm[id]=i  # Recoded identification
}

