# First, run simulation. Ensure languages_before, languages_after, and descriptive_stats are stored in RStudio Global Environment

require(maps)
require(mapdata)
require(ggplot2)
require(gridExtra)
require(ggformula)
require(cobs)
require(ggrepel)

# FAMILIES

# to plot simulated languages, colored by family, on map
plot_families <- function(languages, info) {
  languages$family <- as.factor(languages$family)
  mapWorld <- borders("world", fill="gray50", colour="gray50", xlim=c(-10,200), ylim=c(0,100))
  values <- ggplot(languages, aes(x=longitude, y=latitude, color=family)) + # shape=family)) + 
    scale_shape_manual(values=1:nlevels(languages$family)) +
    mapWorld +
    geom_point(size=6, show.legend = FALSE) + 
    labs(title=info) 
  values
}
plot_families(languages_before, "Before, 500 ts, B/D-2% M-%11")
plot_families(languages_after, "After, 500 ts, B/D-2% M-%11")

# to plot populated places on map
populatedplaces <- populatedplaces[which(populatedplaces$V2 > -15),]
colnames(populatedplaces) <- c("latitude", "longitude")
plot_places<- function(populatedplaces) {
  mapWorld <- borders("world", fill="gray50", colour="gray50", xlim=c(-10,200), ylim=c(0,100))
  values <- ggplot(populatedplaces, aes(x=longitude, y=latitude)) + 
    mapWorld +
    geom_point(size=0.1, color="black") +
    labs(xlab="longitude", ylab="latitude")
  values
}
plot_places(populatedplaces)

# to plot largest families
fams <- levels(as.factor(languages_after$family))
num_fams <- length(fams)

gg_color_hue <- function(n) {
  hues = seq(15, 375, length = n + 1)
  hcl(h = hues, l = 65, c = 100)[1:n]
}

all_colors <- gg_color_hue(num_fams)
families_ordered <- families[order(families)]
my_colors <- all_colors[which(families_ordered %in% biggest_families)]

plot_largest_families <- function(languages, my_colors) {
  mapWorld <- borders("world", fill="gray50", colour="gray50", xlim=c(-10,200), ylim=c(0,100))
  values <- ggplot(languages, aes(x=longitude, y=latitude, color=as.factor(family))) + 
    mapWorld +
    geom_point(size=3.5) + 
    scale_color_manual(name="families", labels=c("first", "second", "third"), values=my_colors)
  values
}

plot_largest_families(langs_three_biggest, my_colors)

# FEATURE VALUE MAPS

before <- languages_before[,-10]
colnames(before)[9] <- "value"
after <- languages_after[,-9]
colnames(after)[9] <- "value"

plot_values <- function(languages) {
  # View(languages)
  # map <- NULL
  mapWorld <- borders("world", fill="gray50", colour="gray50", xlim=c(-10,200), ylim=c(0,100))
  values <- ggplot(languages, aes(x=longitude, y=latitude)) + 
    mapWorld +
    geom_point(size=3.5, mapping = aes(color=as.factor(value))) +
    scale_color_discrete(name = "value", labels = feature_val_names) # to reformat legend
  values
}
plot_values(before)
plot_values(after)

# DESCRIPTIVE STATS

descriptive_stats <- cbind(1:numtimesteps, descriptive_stats)
colnames(descriptive_stats)[1] <- "ts"

plot_stats <- function(info,metric,colname, y_axis_lim, descriptivestats, master_ds_line) {
  # SIMULATION 1
  x1 <- descriptivestats[descriptivestats$sim==1,]$ts
  y1 <- descriptivestats[descriptivestats$sim==1,][,colname]
  constraints <- rbind(
    c(0,1,descriptivestats[descriptivestats$sim==1,][descriptivestats$ts==1,][,colname][1]),
    c(0,750,descriptivestats[descriptivestats$sim==1,][descriptivestats$ts==750,][,colname][1])
  )
  wt1 <- cobs(x1, y1, nknots = 7, pointwise=constraints)
  fit1 <- predict(wt1, x1)[, 'fit']
  df1 <- data.frame(x1, y1, fit1)
  
  # SIMULATION 2
  x2 <- descriptivestats[descriptivestats$sim==2,]$ts
  y2 <- descriptivestats[descriptivestats$sim==2,][,colname]
  constraints <- rbind(
    c(0,1,descriptivestats[descriptivestats$sim==2,][descriptivestats$ts==1,][,colname][1]),
    c(0,750,descriptivestats[descriptivestats$sim==2,][descriptivestats$ts==750,][,colname][1])
  )
  wt2 <- cobs(x2, y2, nknots = 7, pointwise=constraints)
  fit2 <- predict(wt2, x2)[, 'fit']
  df2 <- data.frame(x2, y2, fit2)
  
  # SIMULATION 3
  x3 <- descriptivestats[descriptivestats$sim==3,]$ts
  y3 <- descriptivestats[descriptivestats$sim==3,][,colname]
  constraints <- rbind(
    c(0,1,descriptivestats[descriptivestats$sim==3,][descriptivestats$ts==1,][,colname][1]),
    c(0,750,descriptivestats[descriptivestats$sim==3,][descriptivestats$ts==750,][,colname][1])
  )
  wt3 <- cobs(x3, y3, nknots = 7, pointwise=constraints)
  fit3 <- predict(wt3, x3)[, 'fit']
  df3 <- data.frame(x3, y3, fit3)
  
  # SIMULATION 4
  x4 <- descriptivestats[descriptivestats$sim==4,]$ts
  y4 <- descriptivestats[descriptivestats$sim==4,][,colname]
  constraints <- rbind(
    c(0,1,descriptivestats[descriptivestats$sim==4,][descriptivestats$ts==1,][,colname][1]),
    c(0,750,descriptivestats[descriptivestats$sim==4,][descriptivestats$ts==750,][,colname][1])
  )
  wt4 <- cobs(x4, y4, nknots = 7, pointwise=constraints)
  fit4 <- predict(wt4, x4)[, 'fit']
  df4 <- data.frame(x4, y4, fit3)
  
  stats <- ggplot(descriptivestats, aes(x=ts, y=get(colname))) + 
    labs(title=info,x="time step",y=metric) +
    ylim(y_axis_lim[1], y_axis_lim[2]) +
    scale_x_continuous(breaks = seq(from = 0, to = 750, by = 150),
                       labels = seq(from = 0, to = 750, by = 150)) + 
    geom_point(size = 1, alpha=0.05, aes(color = as.factor(sim))) +
    geom_point(data = master_ds_line, size = 1.7, aes(color = as.factor(sim))) + 
    # NOW, INDIVIDUAL
    geom_line(data = df1, aes(x1, fit1), color = "#777777", size = 0.9) +
    geom_line(data = df2, aes(x2, fit2), color = "#111111", size = 0.9) +
    geom_line(data = df3, aes(x3, fit3), color = "#555555", size = 0.9) +
    geom_line(data = df4, aes(x4, fit4), color = "#999999", size = 0.9) +
    scale_color_manual(name = "simulation", labels = c("L-D/L-IC", "H-D/L-IC", "L-D/H-IC", "H-D/H-IC"), values = c("#777777", "#111111", "#555555", "#999999"))  +
    theme(legend.position="none") 
  stats
}

# genealogical homogeneity
master_ds_line <- master_ds[which(master_ds$ts %in% c(1,750)),] # to include endpoints on graphs
p0 <- plot_stats("Comparison Across Simulations, Genealogical Homogeneity (Parkvall)", "parvall gh", "gh", c(1.1,1.4), master_ds[complete.cases(master_ds),], master_ds_line)
p6 <- plot_stats("Comparison Across Simulations, Genealogical Homogeneity (Wichmann & Holman)", "wh metric c gh", "c", c(0,75), master_ds, master_ds_line)
p0
p6

# areal homogeneity
p1 <- plot_stats("size = 2°", "areal homogeneity", "ah_2", c(1,1.4), master_ds, master_ds_line)
p2 <- plot_stats("size = 5°", "areal homogeneity", "ah_5", c(1,1.4), master_ds, master_ds_line)
p3 <- plot_stats("size = 10°", "areal homogeneity", "ah_10", c(1,1.4), master_ds, master_ds_line)
p4 <- plot_stats("size = 25°", "areal homogeneity", "ah_25", c(1,1.4), master_ds, master_ds_line)
p5 <- plot_stats("size = 50°", "areal homogeneity", "ah_50", c(1,1.4), master_ds, master_ds_line)

grid.arrange(p1,p2,p3,p4,p5, ncol=5, top="Comparison Across Simulations and Area Section Sizes, Areal Homogeneity (Parkvall)")
