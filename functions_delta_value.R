#' Function to calculate Delta Value
#'
#' @param l.in Input list, each element contains one dataframe per variant
#' @param x_data which data to use on x-axis for phase diagram plot, binodal analysis
#' @param y_data which data to use on y-axis for phase diagram plot, binodal analysis
#' @param diagonal_distance How much distance to allow from diagonal, strong influence on binodal value
#' @param bandwidth Which bw in density() function to use
#' @param peak_range_lim Given prior knowledge, where could peak be located. Specify range limits as vector.
#'
#' @return Each element contains original dataframe as variant and additional parameters
#' @export
#'
#' @examples
f.add.binodal <-
  function(l.in = list(),
           x_data = "RFPconc",
           y_data = "GFPconc",
           diagonal_distance = 1,
           bandwidth = 0.5,
           peak_range_lim,
           peak_method = "largest") {
    for (name in names(l.in)) {
      l.in[[name]] <-
        list("df" = l.in[[name]],
             "x" = log2(l.in[[name]][[x_data]]),
             "y" = log2(l.in[[name]][[y_data]]))
      
      #TRUE/FALSE vector: Is the point less than 'diagonal_distance' away from diagonal (orthogonal distance)
      l.in[[name]][["close_to_diag"]]  <-
        abs(l.in[[name]]$x - l.in[[name]]$y) < diagonal_distance
      
      l.in[[name]][["avg_dist"]]  <-
        (l.in[[name]]$x + l.in[[name]]$y) / 2
      
      #logic: density(avg_dist[close_to_diag], from = , to = , bw = )
      l.in[[name]][["density"]]  <-
        density(l.in[[name]][["avg_dist"]][l.in[[name]][["close_to_diag"]]], bw = bandwidth)
      
      l.in[[name]][["peak_and_half"]]  <-
        get_peak(d = l.in[[name]]$density,
                 peak_range = peak_range_lim,
                 method = peak_method)
    }
    l.in
  }


#then returns the x and y coordinates of the rightmost maximum
#input has to be a density

#' Get x,y values of a density's peak and half-peak (half the y, with x bigger than x_peak)
#'
#'
#'
#'
#' @param d Input Density
#' @param peak_range Given prior knowledge, in which range is peak of interest. Specify vector
#' @param method "largest" (default) or "last". "largest" gets highest peak, "last" gets peak with highest x-value
#'
#' @return Returns vector with x_peak, y_peak, x_at_half, y_at_half
#' @export
#'
#' @examples
get_peak <- function(d, peak_range, method = "largest"){
  
  #Find all local maxima by checking that second derivative is negative
  maxima_indices <- which(diff(sign(diff(d$y)))==-2)+1
  
  if (method == "largest") { #Default: get largest peak, i.e. highest y-value
    
    if (missing(peak_range)) {
      peak_y <- max(d$y[maxima_indices])
      peak_x <- d$x[which(d$y == peak_y)]
      
      peak_x_conc <- 2^(peak_x)
      peak_y_conc <- 2^(peak_y)
      
    } else {
      maxima_x <- d$x[maxima_indices]
      x_in_peak_range <- maxima_x[which(maxima_x > peak_range[1] & maxima_x < peak_range[2])]
      
      peak_y <- max(d$y[x_in_peak_range])
      peak_x <- d$x[which(d$y == peak_y)]
      
      peak_x_conc <- 2^(peak_x)
      peak_y_conc <- 2^(peak_y)
      
    }
    
  } else if(method == "last"){ #get last (rightmost) peak, i.e. highest x-value
    
    if (missing(peak_range)) {
      last_index <- tail(maxima_indices, n=1)
      peak_x <- d$x[last_index]
      peak_y <- d$y[last_index]
      
      peak_x_conc <- 2^(peak_x)
      peak_y_conc <- 2^(peak_y)
      
    } else {
      maxima_x <- d$x[maxima_indices]
      x_in_peak_range <- maxima_x[which(maxima_x > peak_range[1] & maxima_x < peak_range[2])]
      peak_x <- tail(x_in_peak_range, n=1)
      peak_y <- d$y[which(d$x == peak_x)]
      
      peak_x_conc <- 2^(peak_x)
      peak_y_conc <- 2^(peak_y)
    }
  }
  
  
  x_after <- d$x > peak_x
  y2 <- d$y
  y2[!x_after] = -1
  
  #This is how you search which entry in vector x is closest to your.number
  #If an entry is exactly your.number then the result is 0, i.e. the minimum
  #which.min(abs(x-your.number))
  half_index = which.min(abs(y2 - ( peak_y/ 2)))
  
  x_peak_half <- d$x[half_index]
  y_peak_half <- d$y[half_index]
  
  x_peak_half_conc <- 2^(x_peak_half)
  y_peak_half_conc <- 2^(y_peak_half)
  
    return(
      c(
        "x_peak" = peak_x,
        "y_peak" = peak_y,
        "x_at_half" = x_peak_half,
        "y_at_half" = y_peak_half,
        "x_peak_conc" = peak_x_conc,
        "y_peak_conc" = peak_y_conc,
        "x_at_half_conc" = x_peak_half_conc,
        "y_at_half_conc" = y_peak_half_conc
      )
    )
}

























