require(ggplot2)

# Need to run simulation and generate power_law dataframe first.

# powerlaw @ each time step
pleachts <- function(powerlaw, numtimesteps, info) {
  # see https://stackoverflow.com/questions/53979358/how-to-override-an-aes-color-controlled-by-a-variable-based-on-a-condition
  library(RColorBrewer)
  
  createPalette <- function(n, colors = 'Greens') {
    max_colors <- brewer.pal.info[colors, ]$maxcolors # Get maximum colors in palette
    palette <- brewer.pal(min(max_colors, n), colors) # Get RColorBrewer palette
    if (n > max_colors) {
      palette <- colorRampPalette(palette)(n) # make it longer i n > max_colros
    }
    
    # assume that  n-th color should be black
    palette[n] <- "#000000"
    
    # return palette
    palette[1:n]
  }
  
  mypalette <- createPalette(length(levels(as.factor(powerlaw[as.numeric(powerlaw$timestep) %% (numtimesteps/10) == 0, ]$timestep))), 'Spectral') #palettes from RColorBrewer
  
  # see https://stackoverflow.com/questions/18305852/power-regression-in-r-similar-to-excel
  power_eqn = function(df, start = list(a=300,b=1)) {
    m = nls(as.numeric(reorder(family,-nlangs)) ~ a*nlangs^b, start = start, data = df)
    eq <- substitute(italic(y) == a  ~italic(x)^b*","~~italic('se')~"="~se*","~~italic(p)~"="~pvalue,
                     list(a = format(coef(m)[[1]], digits = 6), # a
                          b = format(coef(m)[[2]], digits = 6), # b
                          # r2 = format(summary(m)$r.squared, digits = 3), 
                          se = format(summary(m)$parameters[2,'Std. Error'], digits = 6), # standard error
                          pvalue = format(summary(m)$coefficients[2,'Pr(>|t|)'], digits=6) )) # p value (based on t statistic)
    as.character(as.expression(eq))                 
  }
  
  powerlaw$timestep <- as.factor(powerlaw$timestep)
  powerlaw <- powerlaw[as.numeric(powerlaw$timestep) %% (numtimesteps/10) == 0, ]
  
  # nonlinear least squares regression
  plot1 <- ggplot(powerlaw, aes(x = as.numeric(reorder(family,-nlangs)), y = nlangs, color=as.factor(timestep) ) ) + 
    geom_point(aes(x = as.numeric(reorder(family,-nlangs)), y = nlangs ), data=powerlaw[powerlaw$timestep == numtimesteps, ], color="black", shape=1 ) + 
    stat_smooth( method = 'nls', formula = 'y~a*x^b', method.args = list(start= c(a =1,b=1)),se=FALSE, fullrange=TRUE) +
    geom_text(aes(x = as.numeric(reorder(family,-nlangs)), y = nlangs), x = quantile(powerlaw$family)[4], y = max(powerlaw[powerlaw$timestep == numtimesteps, ]$nlangs), label = power_eqn(powerlaw[powerlaw$timestep == numtimesteps, ]), parse = TRUE, size=5, color="black") + 
    theme(axis.ticks.x = element_blank() ) +  
    labs( x = "family", y = "number of languages", title=info ) + 
    scale_color_manual(name = "timestep", values = mypalette) 
  
  plot1
}
pleachts(power_law, 500, "Power Law Distribution in Increments, 500 ts")

# powerlaw @ a given time step
plendofs <- function(powerlaw, whichts, plottitle) { 
  # see https://stackoverflow.com/questions/18305852/power-regression-in-r-similar-to-excel
  
  power_eqn = function(df, start = list(a=300,b=1)) {
    m = nls(family ~ a*nlangs^b, start = start, data = df)
    eq <- substitute(italic(y) == a  ~italic(x)^b, #*","~~italic('se')~"="~se,
                     list(a = format(coef(m)[[1]], digits = 6), # a
                          b = format(coef(m)[[2]], digits = 6)))# , # b
                          # r2 = format(summary(m)$r.squared, digits = 3), 
                          # se = format(summary(m)$parameters[2,'Std. Error'], digits = 6))) # standatd error
    as.character(as.expression(eq))                 
  }
  
  powerlaw <- powerlaw[powerlaw$timestep==whichts,] # plot the desired timestep
  powerlaw <- powerlaw[order(powerlaw$nlangs, decreasing = TRUE),]
  powerlaw$family <- c(1:nrow(powerlaw))
  
  # nonlinear least squares regression
  plot1 <- ggplot(powerlaw, aes(x = family, y = nlangs)) + 
    geom_bar(stat = "identity") + 
    geom_point(shape=1, color="black") + 
    stat_smooth(method = 'nls', formula = 'y~a*x^b', method.args = list(start= c(a=20,b=-.5)),se=FALSE, color="black") +
    geom_text(x = 15, y = 40, label = power_eqn(powerlaw), parse = TRUE, size = 7) + 
    theme(legend.position = "none", axis.ticks.x = element_blank(), plot.title=element_text(hjust=0.5)) + 
    labs( x = "rank of family", y = "number of languages", title=plottitle ) + 
  plot1
}
plendofs(power_law, 500, "Power Law Distribution @ 500th Timestep, 500 ts")

# to plot one of the curves corresponding to a timestep where a family has the highest number of languages out of the entire simulation
highestnlangs <- max(powerlaw[as.numeric(powerlaw$timestep) %% 10 == 0, ]$nlangs)
tswithhighestnlang <- powerlaw[which(powerlaw$nlangs==highestnlangs),]$timestep
sample(tswithhighestnlang, 1)
plendofs(powerlaw, sample(tswithhighestnlang, 1))