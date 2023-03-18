# Source Data and Functions -----------------------------------------------
library(scales) 
library(ggplot2)
library(tidyverse)
library(hrbrthemes)
library(viridis)
library(patchwork)
library(plotly)
library(MASS)
library(gridExtra)
library(grid)
library(testthat)

#WORKING_DIR has to contain data_ready_to_analyze.RData
WORKING_DIR <-  "/Users/d/git_for_cv/phase_analysis_public"

file2load <-
  list.files(
    path = WORKING_DIR,
    pattern = "data_ready_to_analyze.RData",
    full.names = TRUE,
    recursive = TRUE
  )

load(file2load)
source("/Users/d/git_for_cv/phase_analysis_public/phase_plot_functions.R")
source("/Users/d/git_for_cv/phase_analysis_public/functions_delta_value.R")

PLOT_PATH = paste(WORKING_DIR,"Rplots",sep = "/")

if( ! file.exists ( PLOT_PATH ) ){
  dir.create(PLOT_PATH, showWarnings = TRUE)
} 


# Prepare Data for Visualization ------------------------------------------


#Select only those cells in which foci have formed
dat.w_foci <- subset(df.for.analysis, GFP_foci != 0 & RFP_foci != 0)

#Create a list in which each element is one xmer
l.w.foci <- split(x = dat.w_foci, f = dat.w_foci$xmer)

#For each xmer calculate delta and add this information
l.w.foci <- f.add.binodal(l.w.foci)

#Create a list where each element is one biological repeat (i.e. colony) of one xmer
l.w.foci.by.col <- split(dat.w_foci,paste(dat.w_foci$xmer,dat.w_foci$biol_replica, sep = "_"))

#Calculate delta for each colony so that later the average and sd per xmer can be calculated
l.w.foci.by.col <- f.add.binodal(l.w.foci.by.col)

#Retrieve data frame containing only information useful for delta calculation for each colony
#x_peak = x value of peak on density plot, y_peak is corresponding y value
#y_at_half is y_peak/2 and x_at_half is the x value of an actual data point which y value is closest to y_at_half
#all data with the suffix _conc are the concentrations in nM 
bycolony_peaks <- as.data.frame(t(unlist(sapply(l.w.foci.by.col,'[[',"peak_and_half"))))

#Retrieve sorted vector of xmer names
xmers <- sort(unique(dat.w_foci$xmer))

bycolony_peaks[,"xmer"] <- as.factor(rep(xmers,each=6))

bycolonies_mean <- aggregate(x = bycolony_peaks[,c("x_peak","x_peak_conc")], by = list(XMER = bycolony_peaks$xmer), FUN = mean)
bycolonies_mean <- bycolonies_mean[c(1,4:15,2,3),]
bycolonies_sd <- aggregate(x = bycolony_peaks[,c("x_peak","x_peak_conc")], by = list(XMER = bycolony_peaks$xmer), FUN = sd)
bycolonies_sd <- bycolonies_sd[c(1,4:15,2,3),]
bycolonies_se <- bycolonies_sd
bycolonies_se[,2:3] <- bycolonies_se[,2:3]/sqrt(6)


# Data Visualization ------------------------------------------------------

# Check correct determination of delta by plotting density and phase diagram incl delta for every colony
pdf(file= paste(PLOT_PATH,"/check_delta_determination_all_xmer_with_foci.pdf",sep="/"), useDingbats = F)
par(mfrow=c(2,2))

tick_loc = log2(c(1e0,1e1,1e2,1e3,1e4,1e5))
ax_labels = expression(10^0,10^1, 10^2, 10^3, 10^4,10^5)

for(entry in names(l.w.foci.by.col)[order(names(l.w.foci.by.col))]){
  
  plot(l.w.foci.by.col[[entry]]$density,main = "")
  
  x_peak <- l.w.foci.by.col[[entry]]$peak_and_half["x_peak"]
  y_peak <- l.w.foci.by.col[[entry]]$peak_and_half["y_peak"]
  x_peak_half <- l.w.foci.by.col[[entry]]$peak_and_half["x_at_half"]
  y_peak_half <- l.w.foci.by.col[[entry]]$peak_and_half["y_at_half"]
  
  points(x = x_peak, y = y_peak, col = "red")
  points(x = x_peak_half, y = y_peak_half, col = "blue")
  
  mtext(paste("x =",round(x_peak,2),"x_half =",round(x_peak_half,2)), side = 3)
  
  plot(l.w.foci.by.col[[entry]]$x,
       l.w.foci.by.col[[entry]]$y,
       main = paste(entry,"/ n =", length(l.w.foci.by.col[[entry]]$x)),
       pch = 20, cex = 0.4, xlab = "", ylab = "", xlim = c(0,17), ylim = c(0,17),
       col = rgb(0,0,0,0.7), axes = F
  )
  
  axis(1, at = tick_loc, labels = ax_labels)
  axis(2, at = tick_loc, labels = ax_labels)
  
  points(x= x_peak, y= x_peak, bg="red", pch=22, cex = 2)
}

dev.off()

# Bar Plot Delta Values----
# Here the Delta Value for each xmer is summarized in a bar plot including sd values

pdf(file= paste(PLOT_PATH,"/barplot_delta_values_with_foci.pdf",sep="/"), useDingbats = F)
par(mfrow=c(1,1), mar=c(7,5,4,4))

barcolors <-  c(rep("orange", times = 8),rep("red", times = 3),rep("purple", times = 2),"green", "blue")
barlabel <- bycolonies_mean[,"XMER"]

#Plot delta concentrations values (in nM)
delta_nM <- bycolonies_mean[,"x_peak_conc"]
sd_nM <- bycolonies_sd[,"x_peak_conc"]

delta_nM.barplot <- barplot(
  height = delta_nM,
  names.arg = barlabel,
  col = barcolors,
  las = 2, 
  ylab = expression(paste(Delta," (nM)")),
  ylim = c(0,200)
)
arrows(
  delta_nM.barplot,
  delta_nM + sd_nM,
  delta_nM.barplot,
  delta_nM - sd_nM,
  angle = 90,
  code = 3,
  length = .05
)

#Plot log2 transformed delta values (i.e. delta before transforming back to nM)
delta_log <- bycolonies_mean[,"x_peak"]
sd_log <- bycolonies_sd[,"x_peak"]

delta_log.barplot <- barplot(
  height = delta_log,
  names.arg = barlabel,
  col = barcolors,
  las = 2,
  ylab = expression(paste("log2(",Delta,")")),
  ylim = c(0,9)
)
arrows(
  delta_log.barplot,
  delta_log + sd_log,
  delta_log.barplot,
  delta_log - sd_log,
  angle = 90,
  code = 3,
  length = .05
)

dev.off()





# Plot one phase diagram per valence group --------------------------------
tick_loc = log2(c(1e0,1e1,1e2,1e3,1e4,1e5))
ax_labels = expression(10^0,10^1, 10^2, 10^3, 10^4, 10^5)

scaffolds2plot <- c("1olg_326-356","v406","v601","v801","v1201","v2401")
pdf(file= paste(PLOT_PATH,"/one_phase_diagram_per_valence_w_foci.pdf",sep="/"), width = 6, height = 4)
#R default for mar is c(5.1, 4.1, 4.1, 2.1) c(bottom, left, top, right)
par(mfrow=c(2,3), mar = c(2,2.5,2,1), oma = c(2,2,0,0))


for(entry in scaffolds2plot){
  x_peak <- l.w.foci[[entry]]$peak_and_half["x_peak"]
  
  plot(l.w.foci[[entry]]$x, 
       l.w.foci[[entry]]$y, 
       pch = 16, cex = 0.3, xlab="", ylab="",xlim = c(0,17), ylim = c(0,17),
       col = rgb(0,0,0,0.3), axes = F,
       
       #plot grid before data with panel.first
       panel.first = c(abline(h = tick_loc, lty = 2, col = 'lightgrey'),
                       abline(v = tick_loc, lty = 2, col = 'lightgrey'),
                       lines(x= c(0,17),y = c(0,17), col = "grey"))
       
  )
  
  mtext(paste(entry,"/ n =", length(l.w.foci[[entry]]$x)), cex = .75, font = 2)
  
  axis(1, at = tick_loc, labels = ax_labels)
  axis(2, at = tick_loc, labels = ax_labels)
  
  points(x= x_peak, y= x_peak, bg="red", pch=22, cex = 2)
}

mtext(text = "[Dimeric Component] (nM)", side = 1, line = .6, outer = T)
mtext(text = "[Multimeric Component] (nM)", side = 2, line = .5, outer = T)

dev.off()



# Plot all phase diagrams in individual .tiff files (with foci) -----------
if( ! file.exists ( paste0(PLOT_PATH,"/phase_diagrams_sep_files_with_foci") ) ){
  dir.create(paste0(PLOT_PATH,"/phase_diagrams_sep_files_with_foci"), showWarnings = TRUE)
} 

tick_loc = log2(c(1e0,1e1,1e2,1e3,1e4,1e5))
ax_labels = expression(10^0,10^1, 10^2, 10^3, 10^4,10^5)

for(entry in names(l.w.foci)[order(names(l.w.foci))]){
  
  
  tiff(paste0(PLOT_PATH,"/phase_diagrams_sep_files_with_foci/",entry,"_with_foci.tiff"), 
       height = 4, width = 4, units = 'in', res=300, compression = "lzw")
  par(mfrow=c(1,1), mar=c(4.2,4,2,.5)) #c(bottom, left, top, right)
  
  x_peak <- l.w.foci[[entry]]$peak_and_half["x_peak"]
  
  plot(l.w.foci[[entry]]$x, 
       l.w.foci[[entry]]$y, 
       pch = 16, cex = 0.4, xlab = "[Dimer Component] (nM)", ylab = "[Multimer Component] (nM)", xlim = c(0,17), ylim = c(0,17),
       main = paste(entry,"/ n =", length(l.w.foci[[entry]]$x)),
       col = rgb(0,0,0,0.25), axes = F,
       panel.first = c(abline(h = tick_loc, lty = 2, col = 'lightgrey'),
                       abline(v = tick_loc, lty = 2, col = 'lightgrey'),
                       lines(x= c(0,17),y = c(0,17), col = "grey"))
  )
  
  
  axis(1, at = tick_loc, labels = ax_labels)
  axis(2, at = tick_loc, labels = ax_labels)
  
  points(x= x_peak, y= x_peak, bg="red", pch=22, cex = 2)
  
  dev.off()
}



# Plot 2D density ---------------------------------------------------------
# Often used to deal with overplotting (apart from reducing sample size)
# However, here we want to test whether there is hidden information
# Specifically, we want to see whether the phase boundary is markedly different and 
# whether we can use the density plot to deduce the lowest point on the phase boundary
# Therefore, the more sample points the better as we hypothesize that the density in 
# our point of interest increases and might become easier to quantify

xmer_labels_sorted <- names(l.w.foci)[c(1,4:15,2,3)]
tick_loc = log2(c(1e1,1e2,1e3,1e4))
ax_labels = expression(10^1, 10^2, 10^3, 10^4)
ax_limits = c(min(tick_loc)-2,max(tick_loc)+2)

p <- list()

for(xmer in xmer_labels_sorted){
  test <- data.frame(l.w.foci[[xmer]]$x,l.w.foci[[xmer]]$y)
  colnames(test) <- c("x","y")
  
  p[[xmer]] <- ggplot(test, aes(x=x, y=y) ) +
    stat_density_2d(aes(fill = ..density..), geom = "raster", contour = FALSE) +
    scale_x_continuous(breaks = tick_loc, labels = ax_labels, limits = ax_limits, expand = c(0, 0)) +
    scale_y_continuous(breaks = tick_loc, labels = ax_labels, limits = ax_limits, expand = c(0, 0)) +
    scale_fill_viridis() +
    theme(legend.position='none')+
    labs(x = NULL, y=NULL,title = xmer)
}


ml <-
  marrangeGrob(
    p,
    layout_matrix = matrix(
      1:16,
      nrow = 4,
      ncol = 4,
      byrow = T
    ),
    bottom = "[Dimer] (nM)",
    left = "[Multimer] (nM)",
    top = "2D Density of the Phase Diagram for Cells with Condensates"
  )
#use layout_matrix to fill plots by row

ggsave(
  filename = paste(PLOT_PATH,"2d_density.pdf",sep="/"), 
  ml,
  width = 10, height = 10, dpi = 300, units = "in", device='pdf'
)




























































