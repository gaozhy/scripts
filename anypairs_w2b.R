# script to calculate pairwise semantic similarities for words, using word2vec semantic space
# P Hoffman, 31/10/2018

require(lsa)

load("C:/Sync/Data/R/word vectors/word2vec_vectors.RData")
svd=goo
# These vectors were downloaded from: https://code.google.com/archive/p/word2vec/
# They were trained on the Google News dataset (see website for details)


#### now read in yourword pairs
# next line assumes you have copied them into the clipboard
pairs= as.matrix(read.table("clipboard"))
# alternatively could read from a file, e.g.
# pairs=as.matrix(read.csv('myfile.csv'))

npair=dim(pairs)[1]
res=0

for (i in 1:npair) {
  
  x=tolower(pairs[i,1])
  for (j in 1:npair){
    z=tolower(pairs[i,j])
  
  
  #check if both words present in database
    if(x %in% rownames(svd)==FALSE | z %in% rownames(svd)==FALSE) {
      #skip
      res[i,j]=NA
    } else
    
      res[i,j]=cosine(svd[x,],svd[z,])  
  }
}

# this copies your cosine similarities to clipboard
write.table(res,file="clipboard",r=F,c=F)
# (of course you could write to a file instead)
