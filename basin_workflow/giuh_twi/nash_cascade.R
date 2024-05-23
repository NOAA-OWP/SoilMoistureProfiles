# @author Ahmad Jan Khattak
# @email ahmad.jan@noaa.gov
# @date  March 21, 2024

# The function explicitly simulates Nash cascades runoff for parameters n and k
# n : number of reservoirs
# ksub : rate constant [1/hour]; determines the amount of water that can flow from one reservoir to another
# time_h : time in hours (equal to the length of giuh array) 
Nash_Cascade_Runoff <- function(rain_mm, nsubsteps, time_h, factor, N) {
  
  n <- N   
  k <- 100  # this is fixed and determine the ksub, which is calibrated by using `factor`
  
  ksub = (k+1.0)^(1.0/nsubsteps)-1.0  #modify k per hour to k per substep.  (substep time constant)
  ksub = ksub * factor
  dt = 1.0   #  hour
  dtsub = dt/nsubsteps      # hour
  S <- rep(0,n)
  
  total_sub_ts = as.integer(time_h+1.0e-06)/dtsub; # total subtimesteps
  total_sub_ts = as.integer(total_sub_ts)
  
  S[1] <- rain_mm
  volin <- 0.0
  volin <- volin + S[1]
  
  # loop over timesteps
  volout <- 0
  Q_hr <- vector(mode = "numeric")
  Q_ts <- vector(mode = "numeric")
  qts <- 0.0
  
  for (ts in 1:total_sub_ts) {
    for (j in 1:n) {
      Q_r <- S[j]*ksub  # this is Q from res. j to (j+1).  If j==n, this is Qout
      dS <- Q_r*dtsub
      S[j] = S[j] - dS
      
      if (j < n) {
        S[j+1] = S[j+1] + dS
      } else {
        Qout=Q_r
      }
    }
    
    qts = qts + Q_r
    volout <-  volout + dS
    
    if (ts %% nsubsteps == 0) {
      Q_hr <- c(Q_hr, Qout)
      Q_ts <- c(Q_ts, qts*dtsub)
      qts = 0.0
    }
      
  }
    
  volend <- sum(S)
  volout <- sum(Q_ts)
  
  return(list(Q_ts = Q_ts, ksub = ksub, volend = round(volend, 4), volout = volout))
  
}

# The 
# if calib_n_k is set to TRUE, the code calibarates both Nash cascade parameters N and K 
# default calib_n_k is FALSE, i.e. only calibrate K (one parameter)
get_nash_params <- function(giuh_dat_values, calib_n_k = FALSE) {
  
  rainfall_input <- 15    ## mm rainfall (good point, make this to 1.0 and recover params)
  n_cats = length(giuh_dat_values$divide_id)
  
  N_vector <- vector(mode = "numeric", length = n_cats)
  K_vector <- vector(mode = "numeric", length = n_cats)
  
  # loop over all catchments
  for (ncat in 1:n_cats) {
    peak_error_old <- 1000
    substeps <- 10     # hardcoded
      
    d = fromJSON(giuh_dat_values$giuh[ncat]) 
    cat_giuh_ordinates = d$frequency
    
    runoff_giuh <- vector(mode = "numeric", length = length(cat_giuh_ordinates))
    time_h = length(cat_giuh_ordinates) # simulation time in hours for the total runoff
      
    # using giuh ordinates and rainfall, compute runoff per hour
    for (ix in 1:length(cat_giuh_ordinates)) {
      runoff_giuh[ix] =  rainfall_input * cat_giuh_ordinates[ix]  #this should sum up to 1.0
    }
    
    # find index of the peak runoff
    peak_index <- which.max(cat_giuh_ordinates)
    
    # total number of iteration needed to get the convergence; max is 15
    itr <- 1
    if (calib_n_k == TRUE) {
      itr <- 15
    }
    
    # Choose N based on peak index
    if (calib_n_k == FALSE & peak_index  < 3 ) {
      N <- 2
    } else if (calib_n_k == FALSE) {
      N <- 5
    }
      
    for (i in 1:itr) {
      
      if (calib_n_k == TRUE) {
        N <- i+1
        }
        
      # these are local parameters used in the convergence scheme
      error <- 1000 
      error_old <- 1000
      factor <- 1.0
      flag <- FALSE
      
      # the loop alters Nash K (or subk) through a parameter `factor`
      for (j in 1:20) {
        
        factor <- factor * 1.1
        storage <- rep(0, N)
        
        # call to Nash cascade function that routes runoff and returns N, K, and other parameters
        nash_dat <- Nash_Cascade_Runoff(rainfall_input, substeps, time_h, factor, N)
        
        
        # store data in temporary locations; if N and K values are accepted
        # based on the mass balance error and peak location, they will be used 
        # to update the true ksub (which is K) and N 
        runoff_temp_ts = nash_dat$Q_ts # runoff stored in a temporary vector
        ksub_temp = nash_dat$ksub
        volend_temp = nash_dat$volend
        volout_temp = nash_dat$volout
        
        # compute Root Mean Squared Error
        error <- rmse(runoff_temp_ts, runoff_giuh)
        
        # if the suggested N and K increased error in this iteration than error in the
        # previous iteration, then N and K from previous iteration are considered
        # optimized values, we accept those values and stop
        if (error_old < error) {
          flag <- TRUE
          break
        } else {
          runoff_old <- runoff_temp_ts
          ksub_old <-  ksub_temp
          error_old <- error
          volend_old <- volend_temp
          volout_old <- volout_temp
        }
      } # END of K loop optimizer
      
    } #END of N and K loop optimizer (note: this loop runs only once if N is pre-determined based on the peak runoff)
      
    peak_err <- abs(runoff_old[peak_index] - runoff_giuh[peak_index])
      
    if (peak_err > peak_error_old) {
      cat("Rejected... ", N, peak_err, peak_error_old, "\n")
      break
    } else {
      #cat("Accepted... ", N, ksub_old, sum(runoff_old), volend_old, volout_old, volend_old + volout_old, sum(runoff_old) + volend_old, "\n")
      
      peak_error_old <- peak_err
      total_direct_runoff <- runoff_old
      volend <- volend_old
      volout <- volout_old
      ksub <- ksub_old
    }
    
    N_vector[ncat] = N
    K_vector[ncat] = ksub
    
    #cat ("N & K: ", ncat, N, ksub, "\n")
    #var_d <- list(N = N, K = ksub)
    #cat("Final: ", sum(total_direct_runoff), volend, sum(total_direct_runoff) + volend, "\n")
    #return(list(total_direct_runoff = total_direct_runoff, var_d = var_d))
    } # END of iteration over all catchments loop
  
  df_N_K <- data.frame("divide_id" = giuh_dat_values$divide_id, "N_nash" = N_vector, "K_nash" = K_vector)
  
  return (df_N_K)
}
