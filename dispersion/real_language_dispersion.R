require(GeoRange)
require(proj4)
require(maps)
require(mapdata)
require(ggplot2)
require(geosphere)

set.seed(5) 
options(scipen=5) # no scientific notation displayed anywhere

# read data 
# make sure to expand the zip file first
glotto <- read.csv(file="YOUR_FILEPATH_TO/COLING_supplementary_material_language_feature_simulations/dispersion/GlottologData/glotto_languoid.csv", header=TRUE, sep=",")

# remove unnecessary columns
drops <- c("bookkeeping","description", "markup_description", "child_family_count", "child_language_count", "child_dialect_count", "country_ids", "iso639P3code", "parent_id")
glotto <- glotto[ , !(names(glotto) %in% drops)]

# save just family info, for names if needed
glotto_families <- glotto[which(glotto$level == "family"),1:4]
# glotto_families

# remove rows with na/missing values
glotto[glotto==""] <- NA
glotto <- glotto[complete.cases(glotto), ] 

# remove all pseudofamilies - these are spread across multiple continents and not classified phylogenetically
pseudo_families <- read.csv(file="YOUR_FILEPATH_TO/COLING_supplementary_material_language_feature_simulations/dispersion/GlottologData/glotto_families.csv", header=TRUE, sep=",")
pseudo_families <- pseudo_families[which(pseudo_families$category == "Pseudo Family"), c("id","name", "macroareas")]
pfms <- as.vector(pseudo_families$id)
glotto <- glotto[!(glotto$id %in% pfms | glotto$family_id %in% pfms),]

# create list of dataframes
# each dataframe = languages/dialects in a family with their name, latitude, and longitude as columns
# dataframes will be in a list, order of list matches up with order of family names in family_ids
# ^ how to tell which family dataframe to access

family_ids <- unique(glotto$family_id)

create_family_df <- function(name) {
  return(glotto[which(glotto$family_id == name),][,-2])
}

family_df_list <- lapply(family_ids, create_family_df)

# indices are same as family_ids
family_sizes <- unlist(lapply(family_df_list, nrow))

calculate_dispersion <- function(df) {
  longs <- as.vector(df$longitude)
  lats <- as.vector(df$latitude)
  return(CHullArea(longs, lats))
}

dispersions <- unlist(lapply(family_df_list, calculate_dispersion))

# of the families
name_size_dispersion_df <- data.frame(family_ids, family_sizes, dispersions)
colnames(name_size_dispersion_df) <- c("names", "sizes", "dispersions")

# can sort according to family size and according to dispersion
sorted_size <- name_size_dispersion_df[order(name_size_dispersion_df$sizes, decreasing = TRUE),]
sorted_dispersion <- name_size_dispersion_df[order(name_size_dispersion_df$dispersions, decreasing = TRUE),]


linearMod <- lm(dispersions ~ 0 + sizes, data=sorted_size)
print(unname(coef(linearMod)["sizes"]))
# based on real language data from Glottolog, real dispersion distribution slope is 8.84981 (with y intercept set to 0)

lm_eqn <- function(df){
  m <- lm(dispersions ~ sizes, df)
  a <- unname(coef(m)[1])
  b <- unname(coef(m)[2])
  r2 <- summary(m)$r.squared
  eq <- substitute(italic(y) == a + b %.% italic(x), 
                   list(a = format(a, digits = 3),
                        b = format(b, digits = 3)))
  return(list(as.character(as.expression(eq)), c(a, b, r2)))
}

plot_size_vs_dispersion <- function(name_size_dispersion_df, info) {
  data_for_eqn <- name_size_dispersion_df[!(name_size_dispersion_df$sizes %in% c(1,2)),]
  first_part_piecewise <- name_size_dispersion_df[(name_size_dispersion_df$sizes %in% c(1,2)),] 
  
  plot1 <- ggplot(name_size_dispersion_df, aes(x = sizes, y = dispersions)) + 
    geom_point(shape=16, color="blue") +
    geom_smooth(method="lm", formula = y ~ x, data = data_for_eqn, color='blue') + 
    geom_text(x = 350, y = 17000, label = lm_eqn(data_for_eqn)[[1]], parse = TRUE, size = 15) + 
    labs(title=info)
  
  plot1
}

plot_size_vs_dispersion(sorted_size[,2:3], "Dispersion vs. Family Size Distribution of Real Languages")

# regression eqn is only based on families with 3 or more languages - cannot draw convex hull polygon with only two vertices
reg_eqn <- lm_eqn(sorted_size[!(sorted_size$sizes %in% c(1,2)),])[[2]] # c(a, b, r2). eqn follows form y = a + bx

plot_dispersion_on_map <- function(df) {
  mapWorld <- borders("world", fill="gray50", colour="gray50")
  values <- ggplot(df, aes(x=longitude, y=latitude)) + 
    mapWorld +
    geom_point(size=3.5, show.legend = FALSE) # + 
  # labs(title=info) 
  values
}
