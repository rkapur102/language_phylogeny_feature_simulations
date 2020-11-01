# This simulation examines evolutionary dynamics of unstructured features and calculates descriptive statistics on the feature values throughout the simulation.

set.seed(1)

# phylogeny-related parameters
numlangs <- 300
numtimesteps <- 1000 # number of time steps in entire simulation, usually 500 or 1000
chlargeq <- 530 # ch for large quadrilateral
startingsizelargeq <- 1
startingsizesmallq <- 5
# degincrease changes to 0.01 later
degincrease <- 1 # how may times the quadrilaterals increase each time ch isn't met, also the initial quadrilateral dimension
tsuntilnf <- 300 # number of time steps that must pass until languages are considered new families
br <- 0.02 # birth rate
dr <- 0.02 # death rate
mr <- 0.10 # migration rate
nhr <- 0.86 # nothing happening rate

# feature-related parameters
num_feature_values <- 3
feature_val_names <- c("no tones", "simple tone system", "complex tone system")
values <- 1:num_feature_values # for a three-valued feature

featurestability <- 0.3
diffusionprob <- 0.45
internalchangeprob <- 0.45

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

diffusion <- function( langnum, languages, values ) {
  add <- 2
  rg <- c()
  while ( length(rg) < 10000 ) {
    # which migration points have a latitude that is greater than the ( latitude of the language that will move - add ) and also less than the latitude of the language that will move + add)
    i1 <- intersect(which(languages$latitude[-langnum] < languages$latitude[langnum] + add), which(languages$latitude[-langnum] > languages$latitude[langnum] - add))
    # same for longitude
    i2 <- intersect(which(languages$longitude[-langnum] < languages$longitude[langnum] + add), which(languages$longitude[-langnum] > languages$longitude[langnum] - add))
    # rg has the sorted row vector of the common values in vectors i1 and i2 - which ones fit into both bads
    # chooses the ones that fits both i1 and i2 with intersect()
    rg <- intersect(i1,i2)
    add <- add + 2
    if ( length(rg) > 10 ) { # at least 10 nearby languages
      break
    }
  }
  
  # feature 1, feature 2, ...
  probabilities <- numeric(num_feature_values)
  
  rg <- ifelse(rg >= langnum,rg+1,rg)
  
  # get frequency of each feature value as fraction
  probabilities <- unlist(lapply(values, function(x) sum(languages$valafter[rg] == x)))/length(rg)
  
  outcome <- sample(values, size=1, replace=TRUE, prob=probabilities)
  return(outcome)
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

# below is for migration/diffusion
quadrilateral <- function( populatedplaces, coor, startingsizesmallq, languages ) { # for migration
  # drawing just 1 quadrilateral
  add <- startingsizesmallq
  
  l1 <- intersect(which(languages$latitude < coor[1] + add), which(languages$latitude > coor[1] - add))
  l2 <- intersect(which(languages$longitude < coor[2] + add), which(languages$longitude > coor[2] - add))
  
  return(length(intersect(l1, l2))) # number of languages in smaller quadrilateral
}

# descriptive stats

# genealogical homogeneity
parkvall_gh <- function(languages, values, hworld) {
  familyhgs <- data.frame(matrix(NA_real_, nrow = 0, ncol = 2))
  colnames(familyhgs) <- c("family", "hg")
  
  for( i in unique(languages$family)) {
    if( length(which(languages$family==i)) <=1 ) {
      next
    } else {
      currenthg <- data.frame(matrix(NA_real_, nrow = 1, ncol = 2))
      colnames(currenthg) <- c("family", "hg")
      langsinfamily <- which(languages$family==i)
      
      # get frequency of each value for each family
      feature_val_counts <- unlist(lapply(values, function(x) sum(languages$valafter[langsinfamily] == x)))
      
      if( sum(feature_val_counts==0) >= 2 ) {
        next
      }
      
      feature_val_counts <- feature_val_counts/length(langsinfamily)
      
      dg <- 1 - sum(feature_val_counts^2)
      currenthg$family[1] <- i
      currenthg$hg <- 1/dg
      
      familyhgs <- rbind(familyhgs,currenthg)
    }
  }
  
  hfam <- mean(familyhgs$hg)
  cfam <- hfam/hworld
  
  return(cfam)
}
wh_metricc <- function(languages, values) {
  related <- data.frame(matrix(NA_real_, nrow = 0, ncol = 2))
  colnames(related) <- c("family", "R")
  
  # gets weighted related proportions for all families
  for( i in 1:nrow(languages) ) {
    if( (languages$family[i] %in% related$family)==TRUE || length(which(languages$family==languages$family[i]))<=1 ) {
      next
    } else {
      currentrelated <- data.frame(matrix(NA_real_, nrow = 1, ncol = 2))
      colnames(currentrelated) <- c("family", "R")
      langsinfamily <- which(languages$family==languages$family[i])
      
      notonepairs <- 0
      simpletonepairs <- 0
      complextonepairs <- 0
      
      paircount <- 0
      
      n <- length(langsinfamily)
      
      for( a in 1:length(langsinfamily) ) {
        for( b in (a+1):length(langsinfamily) ) {
          if( a==length(langsinfamily) ) {
            break
          }
          paircount <- paircount + 1
          if( languages$valafter[langsinfamily[a]]==languages$valafter[langsinfamily[b]] ) {
            if( languages$valafter[langsinfamily[a]] == values[1] ) {
              notonepairs <- notonepairs + 1
            } else if( languages$valafter[langsinfamily[a]] == values[2] ) {
              simpletonepairs <- simpletonepairs + 1
            } else {
              complextonepairs <- complextonepairs + 1
            }
          }
        }
      }
      
      if( paircount==0 ) {
        next
      }
      
      R <- notonepairs + simpletonepairs + complextonepairs
      R <- R/paircount
      
      R <- R/sqrt(n)
      
      currentrelated$family[1] <- i
      currentrelated$R[1] <- R
      
      related <- rbind(related,currentrelated)
    }
  }
  
  # gets unrelated proportion
  numunrelatedpairs <- 0
  unrelatednotonepairs <- 0
  unrelatedsimpletonepairs <- 0
  unrelatedcomplextonepairs <- 0
  
  for( x in 1:nrow(languages) ) {
    for( y in (x+1):nrow(languages) ) {
      if( languages$family[x]==languages$family[y] || x==nrow(languages) ) {
        next
      }
      
      numunrelatedpairs <- numunrelatedpairs + 1
      
      if( languages$valafter[x]==languages$valafter[y] ) {
        if( languages$valafter[x] == values[1] ) {
          unrelatednotonepairs <- unrelatednotonepairs + 1
        } else if( languages$valafter[x] == values[2] ) {
          unrelatedsimpletonepairs <- unrelatedsimpletonepairs + 1
        } else {
          unrelatedcomplextonepairs <- unrelatedcomplextonepairs + 1
        }
      }
    }
  }
  
  Uf <- unrelatednotonepairs + unrelatedsimpletonepairs + unrelatedcomplextonepairs
  Uf <- Uf/numunrelatedpairs
  
  Rf <- sum(related$R)
  
  Sf <- (Rf-Uf)/(1-Uf)
  
  return(Sf)
}

# areal homogeneity
parkvall_ah <- function(languages, values, hworld, add) {
  # add is the "radius" of the quadrilateral
  
  areahgs <- data.frame(matrix(NA_real_, nrow = 0, ncol = 2))
  colnames(areahgs) <- c("area", "hg")
  
  # draw quadrilaterals around 15 languages
  
  for(i in 1:15) {
    currenthg <- data.frame(matrix(NA_real_, nrow = 1, ncol = 2))
    colnames(currenthg) <- c("area", "hg")
    
    random_lang <- sample(1:nrow(languages),size=1)
    
    lat <- languages$latitude[random_lang]
    long <- languages$longitude[random_lang]
    
    i1 <- intersect(which(languages$latitude[-random_lang] < lat + add), which(languages$latitude[-random_lang] > lat - add))
    
    i2 <- intersect(which(languages$longitude[-random_lang] < long + add), which(languages$longitude[-random_lang] > long - add))
    
    # basically rg
    langs_in_area <- intersect(i1,i2)
    
    if( length(langs_in_area)<=1 ) {
      next
    }
    
    # get frequency of each value for each family
    feature_val_counts <- unlist(lapply(values, function(x) sum(languages$valafter[langs_in_area] == x)))
    
    if( sum(feature_val_counts==0) >= 2 ) {
      next
    }
    
    feature_val_counts <- feature_val_counts/length(langs_in_area)
    
    dg <- 1 - sum(feature_val_counts^2)
    currenthg$area[1] <- i
    currenthg$hg <- 1/dg
    areahgs <- rbind(areahgs,currenthg)
  }
  
  afam <- mean(areahgs$hg)
  cfam <- afam/hworld
  return(cfam)
}

# simulation
simulate <- function(languages, largestfamnum, largestlangnum, numtimesteps, chlargeq, degincrease, values, featurestability, diffusionprob, internalchangeprob, br, dr, mr, nhr) { 
  # copy starting data to the end data column. end data column will be constantly updated to account for feature change
  languages$valafter <- languages$valbefore
  
  # DESCRIPTIVE STATS
  descriptivestats <- data.frame(matrix(NA_real_, nrow = numtimesteps, ncol = 7))
  colnames(descriptivestats) <- c("gh","ah_2", "ah_5", "ah_10", "ah_25", "ah_50", "c")
  
  # power law dataframe
  powerlaw <- data.frame(matrix(NA_real_, nrow = 0, ncol = 3))
  colnames(powerlaw) <- c("family", "nlangs", "timestep")
  
  cat("Beginning simulation for a ", num_feature_values, "-valued feature. Number of time steps: ", numtimesteps, ".", "\n", sep="")
  
  # where new feature values should be placed
  featurecol <- 10
  
  # for each time step
  
  # cannot vectorize the actual simulation because each time step is dependent on the results of the previous. parallel operations will not work.
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
      
      # FEATURE RELATED DYNAMICS
      
      # if the language didn't die off. if it did, no feature work can be done
      if(pmchoice != "d") {
        # all the values for the feature in question, resets for each language
        # values <- c("1 No tones", "2 Simple tone system", "3 Complex tone system") 
        probsofoptions <- c(featurestability, internalchangeprob, diffusionprob)
        
        # choose what type of feature change
        option <- sample(c(1, 2, 3), size=1, replace=TRUE, prob=probsofoptions)
        
        # if option 1, no change, so do nothing
        
        if( option==2 ) { # internal change
          # UNSTRUCTURED FEATURES
          temp_values <- values[!values==languages[a-1,featurecol]]
          languages[a-1,featurecol] <- sample(temp_values, size=1, replace=TRUE, prob=rep(1/length(temp_values), length(temp_values)))
        } else if( option == 3 ) { # diffusion
          # UNSTRUCTURED FEATURES
          languages[a-1,featurecol] <- diffusion(a-1, languages, values)
        }
      }
    }
    if ( all_langs == 0 ) {
      break
    }
    
    # DESCRIPTIVE STATS

    # calculate world homogeneity
    feature_freqs_world <- unlist(lapply(values, function(x) sum(languages$valafter[1:nrow(languages)] == x)))
    feature_freqs_world <- feature_freqs_world/nrow(languages)
    dgworld <- 1 - sum(feature_freqs_world^2)
    hworld <- 1/dgworld
    
    # calculate descriptive stats
    descriptivestats$gh[i] <- parkvall_gh(languages,values,hworld)
    descriptivestats$ah_2[i] <- parkvall_ah(languages,values,hworld,2)
    descriptivestats$ah_5[i] <- parkvall_ah(languages,values,hworld,5)
    descriptivestats$ah_10[i] <- parkvall_ah(languages,values,hworld,10)
    descriptivestats$ah_25[i] <- parkvall_ah(languages,values,hworld,25)
    descriptivestats$ah_50[i] <- parkvall_ah(languages,values,hworld,50)
    descriptivestats$c[i] <- wh_metricc(languages,values)
    
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

  return(list(languages, descriptivestats, powerlaw))
}

languages_before <- setuplanguages(numlangs)
output <- simulate(languages_before, numlangs, numlangs, numtimesteps, chlargeq, degincrease, values, featurestability, diffusionprob, internalchangeprob)
languages_after <- as.data.frame(output[[1]])
power_law <- as.data.frame(output[[2]])
