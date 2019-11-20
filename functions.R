calculateAC2 <- function(nbS, totalWeight, ppm, alpha=.05){

    weight <- c(12,1.0078250321,14.0030740052,15.9949146,31.97207070)

    idx <- which(AC$S==nbS)
    AC2 <- AC[idx, ]
    RR12 <- isoDistr[idx,2]/isoDistr[idx,1]
    RR23 <- isoDistr[idx,3]/isoDistr[idx,2]
    RR34 <- isoDistr[idx,4]/isoDistr[idx,3]
    RR45 <- isoDistr[idx,5]/isoDistr[idx,4]

    tolerance <- ppm*totalWeight/10^6

    test <- calculateRR2(nbS=nbS, monoMass=totalWeight, alpha=alpha)
    idx2 <- which(RR12>=test[1,2] & RR12<=test[1,3])
    idx3 <- which(RR23>=test[2,2] & RR23<=test[2,3])
    idx4 <- which(RR34>=test[3,2] & RR34<=test[3,3])
    idx5 <- which(RR45>=test[4,2] & RR45<=test[4,3])

    idxR <- which(RR45>=test[4,2] & RR45<=test[4,3] & RR12>=test[1,2] & RR12<=test[1,3] & RR23>=test[2,2] & RR23<=test[2,3] & RR34>=test[3,2] & RR34<=test[3,3])
    if(length(idxR)){
        nbCombined <- c(min(AC2$C[idxR]),min(AC2$H[idxR]),min(AC2$N[idxR]),min(AC2$O[idxR]),nbS)
        nbCombinedUP <- c(max(AC2$C[idxR]),max(AC2$H[idxR]),max(AC2$N[idxR]),max(AC2$O[idxR]),nbS)

        tAble <- lapply(nbCombinedUP-nbCombined, seq.int,from=0,by=1)

        nominalMass <- round(calculateNM2(nbS=nbS, monoMass=totalWeight, alpha=alpha))
        nM <-"NA"

        if(length(unique(as.vector(nominalMass)))==1){
            nM <- unique(as.vector(nominalMass))

            if(nbCombined[3]%%2==0){

	            if(nM%%2!=0){
	
		            nbCombined[3] <- nbCombined[3]-1
	            }

            }else{
	
	            if(nM%%2==0){
	
		            nbCombined[3] <- nbCombined[3]-1
	            }
            }
    
            tAble[[3]]<-seq.int(0,nbCombinedUP[3]-nbCombined[3],by=2)

            #new rule for Hydrogen: see conclusions refiningNominalMass.R (pacMASS folder)
            remainderH <- nbCombined[2]%%4
            remainderNM <- nM%%4

            if(remainderNM==0){#remainderH has to be 0

	            if(remainderH!=0){
	
		            nbCombined[2] <- nbCombined[2]-remainderH
	            }

            }else if(remainderNM==1){#remainderH has to be 3
	
	            if(remainderH==0){
	
		            nbCombined[2] <- nbCombined[2]-1 #ok
	            }
                else if(remainderH==1){

                    nbCombined[2] <- nbCombined[2]-2  #ok

                }
                else if(remainderH==2){

                    nbCombined[2] <- nbCombined[2]-3  #ok
                }
            
            }else if(remainderNM==2){#remainderH has to be 2
	
	            if(remainderH==0){
	
		            nbCombined[2] <- nbCombined[2]-2 #ok
	            }
                else if(remainderH==1){

                    nbCombined[2] <- nbCombined[2]-3 #ok

                }
                else if(remainderH==3){

                    nbCombined[2] <- nbCombined[2]-1 #ok
                }
            
            }else if(remainderNM==3){#remainderH has to be 1
	
	            if(remainderH==0){
	
		            nbCombined[2] <- nbCombined[2]-3 #ok
	            }
                else if(remainderH==2){

                    nbCombined[2] <- nbCombined[2]-1 #ok

                }
                else if(remainderH==3){

                    nbCombined[2] <- nbCombined[2]-2 #ok
                }
            }

            tAble[[2]]<-seq.int(0,nbCombinedUP[2]-nbCombined[2],by=4)
        }


        massC <- tAble[[1]]*weight[1]
        massH <- tAble[[2]]*weight[2]
        massN <- tAble[[3]]*weight[3]
        massO <- tAble[[4]]*weight[4]
        massS <- tAble[[5]]*weight[5]

        combinationsAC <- expand.grid(tAble[[1]],tAble[[2]],tAble[[3]],tAble[[4]],tAble[[5]])
        combinationsMASS <- expand.grid(massC,massH,massN,massO,massS)
        colnames(combinationsAC) <- c("C","H","N","O","S")

        masses <- rowSums(combinationsMASS)
        maSS <- masses+sum(nbCombined*weight)
        idX <- which(maSS<=(totalWeight+tolerance)& maSS>=(totalWeight-tolerance))
        res <- t(apply(combinationsAC[idX,],1,function(x,nbCombined) x+nbCombined,nbCombined))
        out2 <- cbind(res,maSS[idX])

        if(dim(out2)[2]!=0){
            if(is.na(nM)){
                valences <- c(4,1,5,6,6,0)
                out3 <-out2[which((out2%*%valences)%%2==0),]
                return(out3[order(out3[,1]),])
            }else{
                return(out2[order(out2[,1]),])
            }
        }else{

            return(NA)
        }
    }else{
        return(NA)
    }
}

calculateNM2 <- function(nbS, monoMass, alpha=.05){

    if(nbS==0){

        beta <- c(-0.0281454811,  0.9995094428)
        MvNM <- 0.00339356789
        corrFact  <- 1.00000380
        meanMass <- 1197.76598
        diffMass <- 113449218239

    }

    if(nbS==1){

        beta <- c(0.00962820343, 0.99950448286)
        MvNM <- 0.00441192839
        corrFact  <- 1.00000821
        meanMass <- 1524.47582
        diffMass <- 78595101626
    }

    if(nbS==2){

        beta <- c(0.0551256003, 0.9994982808)
        MvNM <- 0.00599835281
        corrFact  <- 1.00002548
        meanMass <- 1948.05935
        diffMass <- 29653059529

    }

    if(nbS==3){

        beta <- c(0.104524965, 0.999492285)
        MvNM <- 0.00812510070
        corrFact  <- 1.00008876
        meanMass <- 2376.95397
        diffMass <- 7870618506

    }

    if(nbS==4){

        beta <- c(0.147943249, 0.999492098)
        MvNM <- 0.01110538580
        corrFact  <- 1.00031260
        meanMass <- 2708.08370
        diffMass <- 1808231040

    }

    if(nbS==5){

        beta <- c(0.209461865, 0.999489582)
        MvNM <- 0.01232500256
        corrFact  <- 1.00098912
        meanMass <- 2934.94833
        diffMass <- 451568554

    }

    if(nbS==6){

        beta <- c(0.175224999, 0.999522901)
        MvNM <- 0.01537179455
        corrFact  <- 1.00289855
        meanMass <- 3101.84379
        diffMass <- 162241578


    }

    if(nbS>6){

        beta <- c(0.188238156, 0.999551461)
        MvNM <- 0.01871637
        corrFact  <- 1.00671141
        meanMass <- 3201.99966
        diffMass <- 58949427

    }
    NM <- beta[1] + beta[2]*monoMass
    if(alpha==0.05){
       out <- c(NM, NM-1.96*sqrt(MvNM)*sqrt(corrFact+(monoMass-meanMass)^2/diffMass), NM+1.96*sqrt(MvNM)*sqrt(corrFact+(monoMass--meanMass)^2/diffMass))
     
    }else if(alpha==0.01){
        out <- c(NM, NM-3.29*sqrt(MvNM)*sqrt(corrFact+(monoMass-meanMass)^2/diffMass), NM+3.29*sqrt(MvNM)*sqrt(corrFact+(monoMass--meanMass)^2/diffMass))
    }else{
        out <- c(NM, NM-1.96*sqrt(MvNM)*sqrt(corrFact+(monoMass-meanMass)^2/diffMass), NM+1.96*sqrt(MvNM)*sqrt(corrFact+(monoMass--meanMass)^2/diffMass))
        print("Changed alpha to 0.05")
    }
    return(out)
}

calculateRR2 <- function(nbS, monoMass, alpha=0.05){

	out<-matrix(0,nrow=4,ncol=3)
	
	if(nbS==0){

		beta0 <- c(-0.0189816820,  0.060423108,  0.0305081511,  0.0301270267)
		beta1 <- c(0.5674546622,  0.235662205,  0.2217594086,  0.1735174427)
		beta2 <- c(-0.0216932234,  0.029637619, -0.0238256966, -0.0174466220)
		beta3 <- c(0.0075673616, -0.009869113,  0.0065015479,  0.0039611966)
		beta4 <- c(-0.0009059258,  0.001128834, -0.0006612119, -0.0003479714)

		MvRR  <- c(0.001161310, 0.0001495041, 3.082393e-05, 2.037520e-05)
		corrFact  <- 1.000004
		meanMass  <- 1.197766
		diffMass  <- 113449.21824
		
	}

	if(nbS==1){

		beta0 <- c(-0.025138947,  0.42040388,  0.053970091,  0.0847557940)
		beta1 <- c(0.565173102, -0.29035784,  0.321820468,  0.1669168496)
		beta2 <- c(-0.022024778,  0.35711591, -0.115852624, -0.0196523913)
		beta3 <- c(0.008973144, -0.10020531,  0.035812461,  0.0045485236)
		beta4 <- c(-0.001156060,  0.01020321, -0.003818716, -0.0003689985
)

		MvRR <- c(0.001502627, 0.0002023816, 4.707465e-05, 2.014207e-05)
		corrFact  <- 1.000008
		meanMass  <- 1.524476
		diffMass  <- 78595.10163	


	}

	if(nbS==2){

		beta0 <- c(-0.033893697,  0.72349524,  0.032262337,  0.229835870)
		beta1 <- c(0.565384632, -0.64468507,  0.429087497, -0.021658747)
		beta2 <- c(-0.026136440,  0.53059249, -0.181136449,  0.097994130)
		beta3 <- c(0.012468914, -0.13753241,  0.051527272, -0.027555700)
		beta4 <- c(-0.001777315,  0.01309831, -0.005173147,  0.002782747)

		MvRR <- c(0.002005110, 0.0002673778, 7.106671e-05, 3.152397e-05)
		corrFact  <- 1.000025
		meanMass  <- 1.948059
		diffMass  <- 29653.05953
		


	}

	if(nbS==3){

		beta0 <- c(-0.043228065,  0.97393090,  0.014513355,  0.360287349)
		beta1 <- c(0.565105658, -0.87427089,  0.490702874, -0.163575517)
		beta2 <- c(-0.025561801,  0.61317791, -0.204433409,  0.171747422)
		beta3 <- c(0.011648239, -0.14907758,  0.053527995, -0.044669951)
		beta4 <- c(-0.001532171,  0.01349895, -0.005008776,  0.004245513)

		MvRR <- c(0.002688456, 0.0003123977, 9.609067e-05, 4.079007e-05)
		corrFact  <- 1.000089
		meanMass  <- 2.376954
		diffMass  <- 7870.61851


	}

	if(nbS==4){

		beta0 <- c(-0.091521425,  1.10056929, -0.001487005,  0.43455886)
		beta1 <- c(0.659964385, -0.87136236,  0.538845132, -0.19351672)
		beta2 <- c(-0.097607118,  0.54700360, -0.220722981,  0.16787935)
		beta3 <- c(0.032790006, -0.12107593,  0.054711103, -0.03969725)
		beta4 <- c(-0.003661898,  0.01012674, -0.004898160,  0.00348629)

		MvRR <- c(0.003429367, 0.0003440896, 1.218722e-04, 4.696665e-05)
		corrFact  <- 1.000313
		meanMass  <- 2.708084
		diffMass  <- 1808.23104


	}

	if(nbS==5){

		beta0 <- c(-0.051732786,  1.37515374, -0.038979107,  0.562275988)
		beta1 <- c(0.543672943, -1.15038799,  0.605477833, -0.324066003)
		beta2 <- c(-0.010471838,  0.67223458, -0.243958456,  0.230450372)
		beta3 <- c(0.005775558, -0.14776536,  0.057156580, -0.053389629)
		beta4 <- c(-0.000745030,  0.01231632, -0.004884386,  0.004608342)
		
		MvRR <- c(0.003735255, 0.0003593458, 1.352894e-04, 5.205649e-05)
		corrFact  <- 1.000989
		meanMass  <- 2.934948
		diffMass  <- 451.56855	

	}

	if(nbS==6){

		beta0 <- c(-0.235720198,  1.271460962, -0.074691912,  0.532268779)
		beta1 <- c(0.764850633, -0.795410368,  0.644406611, -0.197992545)
		beta2 <- c(-0.109997374,  0.420266433, -0.243158716,  0.143211133)
		beta3 <- c(0.024599540, -0.078510897,  0.052561006, -0.029599963)
		beta4 <- c(-0.002115823,  0.005616518, -0.004201548,  0.002308974)

		MvRR <- c(0.004005769, 0.0003745374, 1.476534e-04, 5.550083e-05)
		corrFact  <- 1.002899
		meanMass  <- 3.101844
		diffMass  <- 162.24158

	}

	if(nbS>6){

		
		beta0 <- c(0.275916431,  1.43225979,  0.315863733,  0.741687800)
		beta1 <- c(0.084096306, -1.03919326,  -0.046926127, -0.553722694)
		beta2 <- c(0.177867705,  0.62458122, 0.201073973,  0.394132137)
		beta3 <- c(-0.023207940,  -0.14517056, -0.064888740, -0.101692499)
		beta4 <- c(0.000249031,  0.01285476, 0.006830557,  0.009554752)
		
		MvRR <- c(0.005878334, 0.00153863, 0.0006399718, 0.0008042911)
		corrFact  <- 1.006711
		meanMass  <- 3.202
		diffMass  <- 58.7158
		
	}

	RR <- beta0 + beta1*monoMass/1000 + beta2*(monoMass/1000)^2 + beta3*(monoMass/1000)^3 + beta4*(monoMass/1000)^4
		
	out[,1] <- RR
    if(alpha==0.05){
        out[,2] <- RR-1.96*sqrt(MvRR)*sqrt(corrFact+(monoMass/1000-meanMass)^2/diffMass)
	    out[,3] <- RR+1.96*sqrt(MvRR)*sqrt(corrFact+(monoMass/1000--meanMass)^2/diffMass)     
    }else if(alpha==0.01){
        out[,2] <- RR-3.29*sqrt(MvRR)*sqrt(corrFact+(monoMass/1000-meanMass)^2/diffMass)
        out[,3] <- RR+3.29*sqrt(MvRR)*sqrt(corrFact+(monoMass/1000--meanMass)^2/diffMass)
    }else{
        print("Changed alpha to 0.05")
        out[,2] <- RR-1.96*sqrt(MvRR)*sqrt(corrFact+(monoMass/1000-meanMass)^2/diffMass)
	    out[,3] <- RR+1.96*sqrt(MvRR)*sqrt(corrFact+(monoMass/1000--meanMass)^2/diffMass)         
    }
	colnames(out) <- c("fit","lwb","upb")
	return(out)
}
