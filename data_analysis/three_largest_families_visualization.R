# Code to plot the three largest families after a simulation is run
# languages_after must exist (i.e. a simulation must be run and results must be in the RStudio Global Environment). 

require(maps)
require(mapdata)
require(ggplot2)
require(ggrepel)

families <- unique(languages_after$family)
langs_per_f <- unlist(lapply(families, function(x) {
  length(which(languages_after$family == x))
}))

total <- data.frame(families, langs_per_f)
biggest_families <- total[order(-langs_per_f),]
biggest_families <- biggest_families[c(1,2,3),]$families

langs_three_biggest <- languages_after[which(languages_after$family %in% biggest_families),]

# This function can be found in map_plots.R, make sure it is in the Global Environment before running.
plot_families(langs_three_biggest, "3 Biggest Families, 500 ts, INSERT PARAMS") 
