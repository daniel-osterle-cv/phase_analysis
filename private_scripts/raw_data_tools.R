#' This does the initial processing of the txt files that are usually in images_output after 
#' running uscope_analyze_FLUO
#' 
#'
#' @param design.file 
#' @param min.cell.size 
#' @param max.cell.size 
#' @param brightfield.cutoff 
#' @param normalize.concentration 
#' @param chroma.val.gfp.new 
#' @param chroma.val.rfp.new 
#'
#' @return
#' @export
#'
#' @examples
f.process.raw.data <-
  function(design.file,
           min.cell.size = 1000,
           max.cell.size = 2500,
           brightfield.cutoff = 0.8,
           normalize.concentration = FALSE,
           chroma.val.gfp.new = 38455.519,
           chroma.val.rfp.new = 3842.237) {
    
    #Load Raw Data
    l.data.raw <- microscope.load.data(design.file)
    l.data.raw <- uscope.process.add.ncells(data = l.data.raw)
    
    lsos() ## check how much memory the object takes
    design.file <- uscope.process.estimate.background(l.data.raw, design.file)
    
    l.data.temp    <- uscope.process.reorder(l.data.raw, design = design.file)
    l.data.temp    <- uscope.process.remove.first.pic(l.data.temp)
    l.data.temp    <- uscope.process.remove.background(l.data.temp, design.file)
    
    print("Number of cells in raw data:")
    uscope.count.cells(l.data.temp)
    
    l.data.temp    <- uscope.process.remove.small(l.data.temp, MIN.size = min.cell.size, MAX.size = max.cell.size)
    l.data.temp   <- uscope.process.BF(l.data.temp)
    l.data.temp   <- uscope.process.remove.BF.outliers(l.data.temp, cutoff = brightfield.cutoff)
    
    l.data.processed   <- uscope.process.add.ncells(data = l.data.temp)
    
    print("Number of cells after processing raw data:")
    uscope.count.cells(l.data.processed)
    
    # Tidy Up Data ------------------------------------------------------------
    # Create an empty list called 'data'
    l.data.for.analysis <- vector("list")
    
    # Select wells with actual data, i.e. fluo channels are non-NA
    well <- design$PDEF$well[which(design$PDEF$GFP != "NA" & design$PDEF$RFP != "NA")]
    
    # Loop through each well ID
    for (w in well) {
      
      # Define sel and len to use in rep() function below
      sel <- design$PDEF$well == w
      len <- length(l.data.processed$comp[[w]]$cell)
      
      # Create a data frame for the current well ID with the following columns
      l.data.for.analysis[[w]] <- data.frame (
        l.data.processed$comp[[w]]$cell,  # Cell ID
        l.data.processed$comp[[w]]$area,  # Cell area
        l.data.processed$comp[[w]]$pic,   # Picture ID
        
        l.data.processed$comp[[w]]$RFP_int_b5,  # RFP intensity (bin 5 = median)
        l.data.processed$comp[[w]]$GFP_int_b5,  # GFP intensity (bin 5 = median)
        
        l.data.processed$comp[[w]]$f1_inRFP_toRFPmed,  # RFP foci
        l.data.processed$comp[[w]]$f1_inGFP_toGFPmed,  # GFP foci
        
        l.data.processed$comp[[w]]$f1_inGFParea,  # RFP foci size
        
        l.data.processed$comp[[w]]$inGFPnfoci,  # Number of GFP foci per cell
        
        l.data.processed$comp[[w]]$x,  # x-coordinate of cell centroid
        l.data.processed$comp[[w]]$y,  # y-coordinate of cell centroid
        
        
        rep(design$PDEF$tech_replica[sel], times = len),   # Replicate number technical repeats (i.e. including duplications of biol replica) (same value for all cells in current well)
        rep(design$PDEF$bio_replica[sel], times = len),   # Replicate number biological repeats (i.e. real repeats) (same value for all cells in current well)
        
        rep(design$PDEF$GFP[sel], times = len),   # GFP kinase type (same value for all cells in current well)
        rep(design$PDEF$RFP[sel], times = len),   # Replicate number (same value for all cells in current well)
        rep(design$PDEF$well[sel], times = len)   # Well ID (same value for all cells in current well)
      )
      
      #Set column names for the current data frame
      colnames(l.data.for.analysis[[w]]) <- c(
        "cell",                  # Cell ID
        "cellsize",              # Cell area (in pixels)
        "pic",                   # Picture ID
        
        "RFP",                # Median RFP intensity 
        "GFP",                # Maximum GFP intensity
        
        "RFP_foci",              # RFP foci intensity
        "GFP_foci",              # GFP foci intensity
        
        "GFP_foci_size",         # GFP foci size (in pixels)
        
        "nfoci",             # Number of GFP foci per cell
        
        "x",                     # X coordinate of cell center (in pixels)
        "y",                     # Y coordinate of cell center (in pixels)
        
        "tech_replica",          # Replicate number technical repeats (i.e. including duplications of biol replica) (same value for all cells in current well)
        "biol_replica",          # Replicate number biological repeats (i.e. real repeats) (same value for all cells in current well)
        
        "xmer",                  # Replicate ID (same value for all cells in current well)
        "dimer",                 # Name of kinase being measured (same value for all cells in current well)
        "well"                   # Well ID (same value for all cells in current well)
      )
    }
    
    # Combine data frames for each well into a single data frame
    df.for.analyis <- do.call(rbind, l.data.for.analysis)
    

    # NORMALIZATION -----------------------------------------------------------
    if (normalize.concentration == TRUE) {
    
      #library(ggplot2)
      #library(scales) 
      ## create linear model
      
      int.summary = read.csv("/Users/d/git_for_cv/phase_analysis_public/private_scripts/Results-summary.csv")
      int.summary = int.summary[2:8,]
      
      int.summary$GFP = int.summary$GFP-100
      int.summary$RFP = int.summary$RFP-101
      
      
      # ggplot()+
      #   geom_point(data=int.summary, aes(x=GFP, y=GFP.conc), color="green", size=3)+
      #   geom_point(data=int.summary, aes(x=RFP, y=RFP.conc), color="red", size=3)+
      #   geom_smooth(method="lm", data=int.summary, aes(x=GFP, y=GFP.conc), color="green", se = F)+
      #   geom_smooth(method="lm", data=int.summary, aes(x=RFP, y=RFP.conc), color="red", se = F)+
      #   theme(axis.line = element_blank(),panel.grid.major = element_blank(), panel.grid.minor = element_blank(), panel.background = element_blank())+ 
      #   scale_y_continuous(trans=log10_trans(),limits = c(1,20000), breaks = c(1, 10, 100, 1000, 10000))+
      #   scale_x_continuous(trans=log10_trans(),limits = c(1,20000), breaks = c(1, 10, 100, 1000, 10000))
      # 
      # 
      
      lm_fit.R = lm( log10(RFP.conc) ~log10(RFP), data= int.summary)
      summary(lm_fit.R)
      
      
      lm_fit.G = lm(log10(GFP.conc)~log10(GFP ), data= int.summary)
      summary(lm_fit.G)
      
      #get data of normalizing strain to predict concentrations
      
      datnorm = read.csv("/Users/d/git_for_cv/phase_analysis_public/private_scripts/datnorm.csv")
      
      datnorm.l = subset(datnorm, id=="ctrl-log")
      datnorm.s = subset(datnorm, id=="ctrl-sat")
      
      datnorm.l $predicted.GFP = 10^(predict(lm_fit.G, datnorm.l))
      datnorm.l $predicted.RFP = 10^(predict(lm_fit.R, datnorm.l))
      
      datnorm.s $predicted.GFP = 10^(predict(lm_fit.G, datnorm.s))
      datnorm.s $predicted.RFP = 10^(predict(lm_fit.R, datnorm.s))
      
      
      ##################################
      
      
      
    
      
    #These are the reference values of GFP and RFP as measured with 
    #the chromaslide when the regression was performed
    chroma.val.gfp.ref = 49930.831
    chroma.val.rfp.ref = 5900.757
    
    #Save raw intensity values just in case
    df.for.analyis$GFPraw <- df.for.analyis$GFP
    df.for.analyis$RFPraw <- df.for.analyis$RFP
    
    #Not sure why this is still here... delete?
    df.for.analyis$GFP <- df.for.analyis$GFPraw
    df.for.analyis$RFP <- df.for.analyis$RFPraw
    
    #Calculate the factor for calibration
    gfp.calibration.factor <- chroma.val.gfp.ref/chroma.val.gfp.new
    rfp.calibration.factor <- chroma.val.rfp.ref/chroma.val.rfp.new
    
    #Normalize GFP & RFP values by multiplying with calibration factor
    df.for.analyis$GFP <- df.for.analyis$GFP * gfp.calibration.factor
    df.for.analyis$RFP <- df.for.analyis$RFP * rfp.calibration.factor
    
    #Calculate the Concentration using the intensity-concentration regression model
    df.for.analyis$GFPconc <- 10^(predict(lm_fit.G, df.for.analyis))
    df.for.analyis$RFPconc <- 10^(predict(lm_fit.R, df.for.analyis))

    #Rename p53-E9 to pdbCode_AAsequence
    df.for.analyis$xmer[which(df.for.analyis$xmer == "p53-E9")] <- "1olg_326-356"
    }
    
    return(df.for.analyis)
    
  }














































