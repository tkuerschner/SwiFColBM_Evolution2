# title: "source_fun"
# author: "Tobias Kuerschner"
# date: "30 June 2021"

###############################################################################################################################################
# Function 1
#
# "PhaseDiff" - This function takes two timeseries (TSerA, TSerB) and the time (timeObj) and calculates the phasedifference of those two series 
#
###############################################################################################################################################


PhaseDiff <- function(TSerA, TSerB, timeObj)
{
  library(psd) #required for spectral analyis
  t <- timeObj #time
  y1 <- TSerA #timeseries A
  y2 <- TSerB #timeseries B
  
  if (base::sum(y1) == 0 ||
      (base::sum(y2) == 0))
    #catching failed runs
  {
    return(NA) #throwing na for failed runs
  } else{
    # spectral analysis
    out1 <- psd::pspectrum(y1)
    out2 <- psd::pspectrum(y2)
    # frequency with the highest peak in the spectrum
    f1 <- out1$freq[base::which.max(out1$spec)]
    f2 <- out2$freq[base::which.max(out2$spec)]
    if (f1 >= f2)
    {
      f <- f1
    } else{
      f <- f2
    }
    # fitting procedure:
    fit1 <-
      base::lm(y1 ~ base::sin(2 * pi * f * t) + base::cos(2 * pi * f * t))
    fit2 <-
      base::lm(y2 ~ base::sin(2 * pi * f * t) + base::cos(2 * pi * f * t))
    #calculation of phase of y1:
    a1 <- fit1$coefficients[2]
    b1 <- fit1$coefficients[3]
    ph1 <- base::atan(b1 / a1)
    #calculation of phase of y2:
    a2 <- fit2$coefficients[2]
    b2 <- fit2$coefficients[3]
    ph2 <- base::atan(b2 / a2)
    phase_difference <- base::as.numeric((ph2 - ph1) / pi)
    # result:
    phase_difference
    return(phase_difference)
  }
}

#########################################
#             End function 1            #
#########################################


###############################################################################################################################################
# Function 2
#
# "MassDecompTrend" -  helper function to return only the "trend" part of a decomposed time series 
#
###############################################################################################################################################

MassDecompTrend<-function(ts){
  
  temp1<-stats::decompose(ts)
  return(temp1$trend)
}

#########################################
#             End function 2            #
#########################################

