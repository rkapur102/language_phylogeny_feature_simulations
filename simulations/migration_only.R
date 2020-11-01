# NOTE: This is a simulation used solely to examine birth, death, and migration. Feature values are assigned but do not change at all throughout

# Set a seed.
set.seed(20)

# phylogeny-related parameters
numlangs <- 300
numtimesteps <- 500 # number of time steps in entire simulation, usually 500 or 1000
chlargeq <- 530 # ch for large quadrilateral
startingsizelargeq <- 1
startingsizesmallq <- 5
degincrease <- 1 # how may times the quadrilaterals increase each time ch isn't met, also the initial quadrilateral dimension. decreases as ch is neared
tsuntilnf <- 300 # number of time steps that must pass until languages are considered new families
br <- 0.02 # birth rate
dr <- 0.02 # death rate
mr <- 0.11 # migration rate
nhr <- 0.85 # stability rate

# feature-related parameters
num_feature_values <- 3
feature_val_names <- c("no tones", "simple tone system", "complex tone system")
values <- 1:num_feature_values # for a three-valued feature. Example used is TONE

featurestability <- 0.3
diffusionprob <- 0.45
internalchangeprob <- 0.45

# populated places. files compiled by Wichmann (2017)
# MODIFY FILEPATH
populatedplaces <- read.table(file="YOUR_FILEPATH_TO/COLING_supplementary_material_language_feature_simulations/populated_places/EURASIA.txt")
populatedplaces <- populatedplaces[sample(nrow(populatedplaces),50000),]
populatedplaces <- populatedplaces[,-3]

# set up languages dataframe
setuplanguages <- function(nlangs) {
  languages <- data.frame(matrix(NA_real_, nrow = nlangs, ncol = 10))
  colnames(languages) <- c("number", "latitude", "longitude", "family", "ancestors", "descendants", "splitfrom", "timeofsplit", "valbefore", "valafter")
  
  languages$number <- 1:nrow(languages)
  languages$timeofsplit <- 1
  languages$family <- 1:nrow(languages)
  languages$valbefore <- unlist(lapply(languages$valbefore, function(x) sample(values, size=1, replace=TRUE)))
  
  for(i in 1:nrow(languages)) {
    place <- sample(1:nrow(populatedplaces), size=1) 
    languages$latitude[i] <- populatedplaces$V1[place]
    languages$longitude[i] <- populatedplaces$V2[place]
    populatedplaces <- populatedplaces[-place,]
  }
  
  return(languages)
}

# migration code
migration <- function( populatedplaces, coor, startingsizelargeq, startingsizesmallq, chlargeq, degincrease, languages ) {
  # coor = language's current coordinates
  
  add <- startingsizelargeq # big enough so that the language can move
  rg <- c()
  
  count <- 0
  
  while ( length(rg) < 10000 ) {
    # which migration points have a latitude that is greater than the ( latitude of the language that will move - add ) and also less than the latitude of the language that will move + add)
    i1 <- intersect(which(populatedplaces$V1 < coor[1] + add), which(populatedplaces$V1 > coor[1] - add))
    # same for longitude
    i2 <- intersect(which(populatedplaces$V2 < coor[2] + add), which(populatedplaces$V2 > coor[2] - add))
    # rg has the sorted row vector of the common values in vectors i1 and i2 - which ones fit into both bads
    # chooses the ones that fits both with intersect()
    
    rg <- intersect(i1,i2)
    
    if( degincrease != .1 ) {
      if( length(rg) > ( chlargeq - 250 ) && count < 1 ) { # this is a pretty populated place
        degincrease <- .1
      } else if( length(rg) > ( chlargeq - 100 ) ) {
        degincrease <- .1
      } else if( length(rg) > ( chlargeq - 200 ) ) {
        degincrease <- .5
      }
    }
    
    add <- add + degincrease
    count <- count + 1
    
    if ( length(rg) > chlargeq ) {
      break
    }
  }
  
  firstquad <- rg # quadrilateral returns index of populated places in the populated places vector that were in the quadrilateral
  
  locations <- sample(firstquad, 10)
  
  locationquadrilaterals <- data.frame(matrix(NA_real_, nrow = 10, ncol = 2))
  colnames(locationquadrilaterals) <- c("place", "count")
  
  locationquadrilaterals$place <- locations
  locationquadrilaterals$count <- locations
  
  locationquadrilaterals$count <- unlist(lapply(locationquadrilaterals$count, function(x) {
    coordinates <- c( populatedplaces$V1[x], populatedplaces$V2[x] )
    x <- quadrilateral( populatedplaces, coordinates, startingsizesmallq, languages )
  }))
  
  finalquad <- unlist(locationquadrilaterals$place[which.min(locationquadrilaterals$count)] , use.names = FALSE)
  
  info <- c(populatedplaces$V1[finalquad], populatedplaces$V2[finalquad], finalquad) 
  return(info)
}

# below is also for migration
quadrilateral <- function( populatedplaces, coor, startingsizesmallq, languages ) { # for migration
  # drawing just 1 quadrilateral
  add <- startingsizesmallq
  
  l1 <- intersect(which(languages$latitude < coor[1] + add), which(languages$latitude > coor[1] - add))
  l2 <- intersect(which(languages$longitude < coor[2] + add), which(languages$longitude > coor[2] - add))
  
  return(length(intersect(l1, l2))) # number of languages in smaller quadrilateral
}

# simulation
simulate <- function(languages, largestfamnum, largestlangnum, numtimesteps, chlargeq, degincrease, values, featurestability, diffusionprob, internalchangeprob) { 
  # copy starting data to the end data column. 
  languages$valafter <- languages$valbefore
  
  # power law dataframe
  powerlaw <- data.frame(matrix(NA_real_, nrow = 0, ncol = 3))
  colnames(powerlaw) <- c("family", "nlangs", "timestep")
  
  cat("Beginning simulation for a ", num_feature_values, "-valued feature. Number of time steps: ", numtimesteps, ".", "\n", sep="")
  
  # where new feature values should be placed
  featurecol <- 10
  
  # for each time step
  
  # cannot vectorize the actual simulation because each time step is dependent on the results of the previous. parallel operations with the apply family of functions will not work.
  for(i in 1:numtimesteps) {
    cat("Step ", i, "\n", sep="")
    
    # go through each language at each time step
    a <- 1
    all_langs <- nrow(languages) # current number of languages
    while( a <= all_langs ) {
      # check if the language needs to be considered a new family
      if( !is.na(languages$timeofsplit[a]) && ( i > languages$timeofsplit[a] + tsuntilnf ) ) {
        largestfamnum <- largestfamnum + 1
        languages$family[a] <- largestfamnum
        otherlang <- which(languages$number==languages$splitfrom[a])
        languages$splitfrom[otherlang] <- NA_real_
        languages$timeofsplit[otherlang] <- NA_real_
        languages$timeofsplit[a] <- NA_real_
        languages$splitfrom[a] <- NA_real_
        
        descendants <- languages$descendants[[a]]
        for( convert in 1:length(descendants) ) {
          languages$family[which(languages$number==descendants[convert])] <- largestfamnum
        }
      }
      
      # NON-FEATURE RELATED DYNAMICS
      
      pmoptions <- c(br, dr, mr, nhr) # language birth rate, language death rate, migration rate, stability rate
      pmchoice <- sample(c("b","d","m", "s"), size=1, replace=TRUE, prob=pmoptions)
      
      # this accounts for birth, death, and then migration
      if( pmchoice=="b" ) { # language birth
        newlang <- languages[a,] # creating space for the new language
        
        newcoordinates <- migration(populatedplaces, c(languages$latitude[a], languages$longitude[a]), startingsizelargeq, startingsizesmallq, chlargeq, degincrease, languages) # the newborn language's starting location will be a migration point from the language it originated from
        
        largestlangnum <- largestlangnum + 1
        newlang$number[1] <- largestlangnum # lang number
        newlang$latitude[1] <- as.numeric(newcoordinates[1]) # lat
        newlang$longitude[1] <- as.numeric(newcoordinates[2]) # long
        newlang$splitfrom[1] <- languages$number[a] # partner in the language split
        newlang$timeofsplit[1] <- i # time of split
        
        languages$splitfrom[a] <- newlang$number[1] # partner in the language split
        languages$timeofsplit[a] <- i # time of split
        
        # new languages inherit feature values, so no need to change
        
        # adding the parent language to the list of ancestors
        newlang$ancestors[1] <- ifelse(is.na(newlang$ancestors[1]), list(c(languages$number[a])), list(c(as.vector(newlang$ancestors[[1]]), languages$number[a])))
        
        # adding this language to the list of descendants 
        languages$descendants[a] <- ifelse(is.na(languages$descendants[a]), list(c(newlang$number[1])), list(c(as.vector(languages$descendants[[a]]), newlang$number[1])))
        
        # adding the new language to the list of descendants for all ancestors the parent may have
        if( !is.na(languages$ancestors[a]) ) {
          for( b in 1:length(languages$ancestors[[a]]) ) {
            index <- which(languages$number==languages$ancestors[[a]][b])
            if( is.na(languages$descendants[index]) ) {
              languages$descendants[index] <- list(c(newlang$number[1]))
            } else {
              languages$descendants[index] <- list(c(as.vector(languages$descendants[[index]]), newlang$number[1]))
            }
          }
        }
        
        languages <- rbind( languages, newlang ) # adding new language to end of the master dataframe
        
        a <- a + 1
        all_langs <- nrow(languages)
      } else if( pmchoice == "d" ) { # death
        # removes the dead language from the lists of descendants and ancestors for all other languages
        languages$descendants <- lapply(languages$descendants, function(x){ y <- x[x != languages$number[a]]; if(length(y) < 1) { y <- NA_real_ }; y })
        languages$ancestors <- lapply(languages$ancestors, function(x){ y <- x[x != languages$number[a]]; if(length(y) < 1) { y <- NA_real_ }; y })
        
        languages <- languages[-a,] 
        
        if( length(languages) == 0 ) {
          break
        }
        all_langs <- nrow(languages)
        next # because the language has been removed from the simulation, no feature change manipulation necessary
      } else if( pmchoice == "m" ) { # migration
        newcoordinates <- migration(populatedplaces, c(languages$latitude[a], languages$longitude[a]), startingsizelargeq, startingsizesmallq, chlargeq, degincrease, languages )
        
        # replacing the new coordinates with the old coordinates so they can free up again
        populatedplaces$V1[as.numeric(newcoordinates[3])] <- languages$latitude[a] 
        populatedplaces$V2[as.numeric(newcoordinates[3])] <- languages$longitude[a]      
        
        languages$latitude[a] <- as.numeric(newcoordinates[1])
        languages$longitude[a] <- as.numeric(newcoordinates[2])
        
        a <- a + 1
        all_langs <- nrow(languages)
      } else {
        a <- a + 1
      }
      
      if( all_langs == 0 ) {
        break
      }
    }
    
    if ( all_langs == 0 ) {
      break
    }
    
    # POWER LAW
    
    families = as.factor(languages$family)
    
    temp <- data.frame(matrix(NA_real_, nrow = nlevels(families), ncol = 3) )
    colnames(temp) <- c("family", "nlangs", "timestep")
    
    levels <- as.numeric(levels(families))
    temp$family <- 1:nrow(temp)
    temp$nlangs <- 1:nrow(temp)
    temp$timestep <- i
    
    temp$nlangs <- unlist(lapply(temp$nlangs, function(x) {length(which(languages$family == levels[x]))}))
    
    powerlaw <- rbind(powerlaw, temp)
  }

  return(list(languages, powerlaw))
}

languages_before <- setuplanguages(numlangs)
output <- simulate(languages_before, numlangs, numlangs, numtimesteps, chlargeq, degincrease, values, featurestability, diffusionprob, internalchangeprob)
languages_after <- as.data.frame(output[[1]])
power_law <- as.data.frame(output[[2]])
