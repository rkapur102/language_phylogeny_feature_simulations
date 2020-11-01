# IMPORTANT: some functions are in real_language_dispersion.R. load into environment before running this
# Must run a simulation and generate the languages_after dataframe before running this script.

require(GeoRange)

family_nums_s <- unique(languages_after$family)

create_family_df_simulation <- function(num) {
  return(languages_after[which(languages_after$family == num),][,1:3])
}

family_df_list_s <- lapply(family_nums_s, create_family_df_simulation)

family_sizes_s <- unlist(lapply(family_df_list_s, nrow))

calculate_dispersion <- function(df) {
  longs <- as.vector(df$longitude)
  lats <- as.vector(df$latitude)
  return(CHullArea(longs, lats))
}

dispersions_s <- unlist(lapply(family_df_list_s, calculate_dispersion))

num_size_dispersion_df_s <- data.frame(family_nums_s, family_sizes_s, dispersions_s)
colnames(num_size_dispersion_df_s) <- c("num", "sizes", "dispersions")

sorted_size_s <- num_size_dispersion_df_s[order(num_size_dispersion_df_s$sizes, decreasing = TRUE),]
sorted_dispersion_s <- num_size_dispersion_df_s[order(num_size_dispersion_df_s$dispersions, decreasing = TRUE),]

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

plot_size_vs_dispersion_s <- function(name_size_dispersion_df, info) {
  data_for_eqn <- name_size_dispersion_df[!(name_size_dispersion_df$sizes %in% c(1,2)),]
  first_part_piecewise <- name_size_dispersion_df[(name_size_dispersion_df$sizes %in% c(1,2)),] 
  
  plot1 <- ggplot(name_size_dispersion_df, aes(x = sizes, y = dispersions)) + 
    geom_point(shape=16, color="red") +
    geom_smooth(method="lm", formula = y ~ x, data = data_for_eqn, color='red') + 
    geom_text(x = 10, y = 1000, label = lm_eqn(data_for_eqn)[[1]], parse = TRUE, size = 15) + 
    labs(title=info)
  
  plot1
}

plot_title <- "Dispersion vs. Family Size for Simulation, 500 ts, B/D-2% M-%11"
plot_size_vs_dispersion_s(sorted_size_s[,2:3], plot_title)
reg_eqn <- lm_eqn(sorted_size_s[!(sorted_size_s$sizes %in% c(1,2)),])[[2]] # c(a, b, r2). eqn follows form y = a + bx
