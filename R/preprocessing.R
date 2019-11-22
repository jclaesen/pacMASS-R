filterMonoMass <- function(monoMass, lowerLimit=0, upperLimit=4000){
  
  if((monoMass >= lowerLimit) & (monoMass <= upperLimit))
    return(monoMass)
  
}

convertNeutralMass <- function(mz, z){
  
  neutralMass <- (mz-1.00794)*z
  return(neutralMass)
}
