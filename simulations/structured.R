# This simulation examines evolutionary dynamics of structured features. Supports transitional probabilities. 

# Set seed for reproducibility.
set.seed(20)

require(openxlsx)

# phylogeny-related parameters
# phylogeny-related parameters
numlangs <- 300
numtimesteps <- 100 # number of time steps in entire simulation, usually 500 or 1000
chlargeq <- 320 # ch for large quadrilateral
startingsizelargeq <- 1
startingsizesmallq <- 5
# degincrease changes to 0.01 later
degincrease <- 1 # how may times the quadrilaterals increase each time ch isn't met, also the initial quadrilateral dimension
tsuntilnf <- 300 # number of time steps that must pass until languages are considered new families
br <- 0.02 # birth rate
dr <- 0.02 # death rate
mr <- 0.05 # migration rate
nhr <- 0.91 # nothing happening rate

# feature-related parameters
num_feature_values <- 2
feature_val_names <- c("simple", "complex")
values <- 1:num_feature_values # for a two-valued feature

# probabilities of changing from feature to feature are different
# adjacency matrix to represent it, where weights on edges are probabilities
# the graph will be directed for structural features, undirected for unstructured ones
# if we want to have the case where there is no specific probability for changing from feature to feature, then set all probabilities to be equal

# samples a two-feature. set your desired transitional probabilities here
transitionals <- c(.8, .2)
matrix_vals <- c(.8,.8,.2,.2)
probabilities <- data.frame(matrix(matrix_vals,num_feature_values,num_feature_values))

# feature stability, diffusion, and internal change probabilities are now defined in feature section of simulation code -> they depend on given transitional probabilities

# populated places
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
    # chooses the ones that fits both i1 and i2
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
  return( outcome )
}

# migration code
migration <- function( populatedplaces, coor, startingsizelargeq, startingsizesmallq, chlargeq, degincrease, languages) {
  # coor = language's current coordinates
  
  add <- startingsizelargeq # make it big enough so that the language can move
  rg <- c()
  
  count <- 0
  
  while ( length(rg) < 10000 ) {
    # which migration points have a latitude that is greater than the ( latitude of the language that will move - add ) and also less than the latitude of the language that will move + add)
    i1 <- intersect(which(populatedplaces$V1 < coor[1] + add), which(populatedplaces$V1 > coor[1] - add))
    # same for longitude
    i2 <- intersect(which(populatedplaces$V2 < coor[2] + add), which(populatedplaces$V2 > coor[2] - add))
    # rg has the sorted row vector of the common values in vectors i1 and i2 - which ones fit into both bads
    # chooses the ones that fits both
    
    rg <- intersect(i1,i2)
    
    if( degincrease != .1 ) {
      if( length(rg) > ( chlargeq - 250 ) && count < 1 ) { # this is a highly populated place
        # browser()
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
simulate <- function(languages, largestfamnum, largestlangnum, numtimesteps, chlargeq, degincrease, values, featurestability, diffusionprob, internalchangeprob, transitionals, probabilities) { 
  # copy starting data to the end data column
  languages$valafter <- languages$valbefore
  
  # PERCENTAGE OF DOMINANT FEATURE VALUE
  dom_feature_percent <- data.frame(matrix(NA_real_, nrow = numtimesteps, ncol = 1))
  colnames(dom_feature_percent) <- c("percent")
  
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
      
      pmoptions <- c(br, dr, mr, nhr) # language birth rate, language death rate, migration rate, nothing happening rate
      pmchoice <- sample(c("b","d","m", "s"), size=1, replace=TRUE, prob=pmoptions)
      
      # this accounts for birth, death, and then migration
      if( pmchoice=="b" ) { # language birth
        newlang <- languages[a,] # creating space for the new language
        
        newcoordinates <- migration(populatedplaces, c(languages$latitude[a], languages$longitude[a]), startingsizelargeq, startingsizesmallq, chlargeq, degincrease, languages) # the newborn language's starting
        # location will be a migration point from the language it originated from
        
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
        current_val <- languages[a-1,featurecol][[1]]
        
        featurestability <- probabilities[current_val,current_val]
        diffusionprob <- .75 * (1-featurestability)
        internalchangeprob <- .25 * (1-featurestability)
        
        # these are dependent on change probabilities defined above
        probsofoptions <- c(featurestability, internalchangeprob, diffusionprob) # 50% chance of staying the same, 40% chance of internal change, 10% chance of diffusion
        
        # choose what type of feature change
        option <- sample(c(1, 2, 3), size=1, replace=TRUE, prob=probsofoptions)
        
        # if option 1, no change, so do nothing
        
        if( option==2 ) { # internal change
          # FOR 2-VALUED STRUCTURED FEATURE
          languages[a-1,featurecol] <- values[values != current_val] # with internal change, it always changes. diffusion, can stay the same.
          # FOR 2+-VALUED STRUCTURED FEATURE, USE THIS:
          # temp_values <- values[values != current_val]
          # internalchangeoutcome <- sample(temp_values, size=1, replace=TRUE, prob=rep((1/length(temp_values)),length(temp_values))) 
          # languages[a-1,featurecol] <- ifelse(internalchangeoutcome > current_val, current_val + 1, ifelse(internalchangeoutcome < current_val, current_val - 1, current_val))
          
        } else if( option == 3 ) { # diffusion
          # FOR STRUCTURED FEATURES 
          diffusion_result <- diffusion(a-1, languages, values)
          languages[a-1,featurecol] <- ifelse(diffusion_result > current_val, current_val + 1, ifelse(diffusion_result < current_val, current_val - 1, current_val))
        }
      }
    }
    if ( all_langs == 0 ) {
      break
    }
    
    # PERCENTAGE OF DOMINANT FEATURE VALUE
    proportions <- table(languages[,featurecol])/length(languages[,featurecol])
    dom_feature_percent$percent[i] = proportions[1]
  }
  
  return(list(languages, dom_feature_percent))
}

languages_before <- setuplanguages(numlangs)
output <- simulate(languages_before, numlangs, numlangs, numtimesteps, chlargeq, degincrease, values, featurestability, diffusionprob, internalchangeprob, transitionals, probabilities)
languages_after <- as.data.frame(output[[1]])
dom_feature_percent <- as.data.frame(output[[2]])
