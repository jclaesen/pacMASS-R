library(tictoc)
setwd("C:\\Users\\jurgen\\Documents\\Research projects\\postDoc\\proteomics\\predict atomic composition with mass\\pacMass_prep\\transferred to Eva")

tic()
library("BRAIN")
#library("plyr")

source("functions.R")
AC <<-  read.table("AC_matrix_2.txt", header=TRUE, row.names=NULL, sep="\t")
isoDistr <<- read.table("iso_matrix_2.txt", header=TRUE, row.names=NULL, sep="\t")

data <- read.csv("test1.csv")

###part 0: convert m0/z to uncharged m0: data is a file which contains the monoisotopic mass of the selected peaks/peptides
monoMass <- data$m.z*data$z - data$z*1.00794
Aaseq <- data$Sequence
idAC <- data.frame(matrix(0,nrow=length(monoMass), ncol=6),stringsAsFactors=FALSE)
colnames(idAC) <- c("mass","C", "H", "N", "O", "S")

for(l in 1:length(monoMass)){
    tmpO <- Aaseq[l]
    tmp <- getAtomsFromSeq(tmpO)
    idAC$mass[l] <- monoMass[l]
    idAC[l,2:6] <- unlist(tmp)
}

resultsMONO_noS <- mapply(calculateAC2, nbS=0, totalWeight=monoMass, ppm=10)
resultsMONO_oneS <- mapply(calculateAC2, nbS=1, totalWeight=monoMass, ppm=10)
resultsMONO_twoS <- mapply(calculateAC2, nbS=2, totalWeight=monoMass, ppm=10)
resultsMONO_threeS <- mapply(calculateAC2, nbS=3, totalWeight=monoMass, ppm=10)

toc()

monoMass <- c(1045.53451, 1061.6782, 1046.53451, 1047.53451, 1048.53451)

tic()
resultsMONO_noS <- mapply(calculateAC2, nbS=0, totalWeight=monoMass, ppm=5)
resultsMONO_oneS <- mapply(calculateAC2, nbS=1, totalWeight=monoMass, ppm=5)
resultsMONO_twoS <- mapply(calculateAC2, nbS=2, totalWeight=monoMass, ppm=5)
resultsMONO_threeS <- mapply(calculateAC2, nbS=3, totalWeight=monoMass, ppm=5)
toc()

results <- microbenchmark(mapply(calculateAC2, nbS=0, totalWeight=monoMass, ppm=5))

monoMass = 1045.53451
tic()
resultsMONO_noS <- calculateAC2(nbS=0, totalWeight=monoMass,ppm=5)
resultsMONO_oneS <- calculateAC2(nbS=1, totalWeight=monoMass,ppm=5)
resultsMONO_twoS <- calculateAC2(nbS=2, totalWeight=monoMass,ppm=5)
resultsMONO_threeS <- calculateAC2(nbS=3, totalWeight=monoMass,ppm=5)
toc()



#check if AC predicted by pacMASS match the AC of the peptide sequences
candidates <- rep(0,length(monoMass))
matches_noS <- rep(0,length(monoMass))
matches_oneS <- rep(0,length(monoMass))
matches_twoS <- rep(0,length(monoMass))
matches_threeS <- rep(0,length(monoMass))

idxFound_noS <- rep(0,length(monoMass))
idxFound_oneS <- rep(0,length(monoMass))
idxFound_twoS <- rep(0,length(monoMass))
idxFound_threeS <- rep(0,length(monoMass))

idAC2 <- idAC[,c(2:6)]

for(k in 1:length(monoMass)){
    candidates[k] <- sum(dim(resultsMONO_noS[[k]])[1],dim(resultsMONO_oneS[[k]])[1],dim(resultsMONO_twoS[[k]])[1],dim(resultsMONO_threeS[[k]])[1])
    toBeCompared <- idAC2[k,]
    if(!class(dim(resultsMONO_noS[[k]])[1])=="NULL"){
        tmp_noS <- which(apply(resultsMONO_noS[[k]][,1:5],1,function(x) identical(as.integer(x),as.integer(toBeCompared))))
        if(length(tmp_noS)!=0){
            idxFound_noS[k] <- as.numeric(tmp_noS)
            matches_noS[k] <- length(tmp_noS)
        }
    }
    if(!class(dim(resultsMONO_oneS[[k]])[1])=="NULL"){
        tmp_oneS <- which(apply(resultsMONO_oneS[[k]][,1:5],1,function(x) identical(as.integer(x),as.integer(toBeCompared))))
        if(length(tmp_oneS)!=0){
            idxFound_oneS[k] <- as.numeric(tmp_oneS)
            matches_oneS[k] <- length(tmp_oneS)
        }
    }
    if(!class(dim(resultsMONO_twoS[[k]])[1])=="NULL"){
        tmp_twoS <- which(apply(resultsMONO_twoS[[k]][,1:5],1,function(x) identical(as.integer(x),as.integer(toBeCompared))))
        if(length(tmp_twoS)!=0){
            idxFound_twoS[k] <- as.numeric(tmp_twoS)
            matches_twoS[k] <- length(tmp_twoS)
        }
    }
    if(!class(dim(resultsMONO_threeS[[k]])[1])=="NULL"){
        tmp_threeS <- which(apply(resultsMONO_threeS[[k]][,1:5],1,function(x) identical(as.integer(x),as.integer(toBeCompared))))
        if(length(tmp_threeS)!=0){
            idxFound_threeS[k] <- as.numeric(tmp_threeS)
            matches_threeS[k] <- length(tmp_threeS)
        }
    }
}

found_noS <- sum(matches_noS)
found_oneS <- sum(matches_oneS)
found_twoS <- sum(matches_twoS)
found_threeS <- sum(matches_threeS)

