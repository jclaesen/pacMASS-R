pacmass <- function(mz=NULL, z=NULL, fileNAME=NULL, nbS, ppm=5, alpha=0.05){

  #' Predict the elemental composition of proteins and peptides
  #'
  #' @param mz monoisotopic mass-to-charge value
  #' @param z charge
  #' @param fileNAME name of a file containing mz and z.
  #' @param nbS number of S-atoms
  #' @param ppm mass error tolerance
  #' @param alpha significance level
  #' @details The default value of mz, and fileNAME is NULL. However one of these parameters has to be defined in order for pacmass to work.
  #'
  #' When mz is defined and z=NULL, it is assumed that mz is the mass of a neutral molecule.
  #'
  #' When fileNAME is defined, the file has to be a comma-separated (csv) file, a tab-separated txt file or a tab-separated tsv file. The mass values have to be in a column with as name 'mz', the charges have to be in column with as name 'z'.
  #' @return A list with the predicted elemental compositions
  #' @examples
  #' pacmass(mz=1709.897, nbS=0, ppm=10, alpha=0.05)
  #' pacmass(mz=855.956, z=2, nbS=0)
  #' pacmass(mz=855.956, z=2, nbS=c(0,1))
  #' pacmass(mz=c(855.956,735.693), z=c(2,3), nbS=0)
  #' @export
  #' @import utils
  #' @import stats


  #globalVariables(names(".sysdata.rda"), package="pacMASS", add=FALSE)


  if (is.null(fileNAME)) {
    if (is.null(z)) {
      monoMass <- mz
    } else if (is.numeric(z)) {
      monoMass <- convertNeutralMass(mz, z)
    } else {
      stop(paste0("the defined charge: ",z, " is not valid"))
    }
  } else if (is.character(fileNAME)) {
    if (endsWith(fileNAME,'.txt') | endsWith(fileNAME, '.tsv')) {
      mzFile <- read.table(fileNAME, header=TRUE, sep="\t")
      if (all(c("mz","z")%in%colnames(mzFile))) {
        monoMass <- convertNeutralMass(mzFile$mz, mzFile$z)
      } else {
        stop("The file does not have mz and/or z as column names")
      }
    } else if (endsWith(fileNAME,'.csv')) {
      mzFile <- read.csv(fileNAME, header=TRUE)
      if (all(c("mz","z")%in%colnames(mzFile))) {
        monoMass <- convertNeutralMass(mzFile$mz, mzFile$z)
      } else {
        stop("The file does not have mz and/or z as column names")
      }
    } else stop(paste0("file: ", fileNAME,"is not supported"))
  }

  lenS <- length(nbS)
  results <- setNames(vector("list", lenS), as.character(nbS))

  for(k in 1:lenS){
    results[[k]] <- lapply(monoMass, calculateAC, nbS=nbS[k], ppm=ppm, alpha=alpha)
  }

  return(results)
}
