# conditions for simulation - supplement #
# initiated: 14-may-2021 #

# (supplement: just 20 more reps per condition of the simulation 11-may-2021)

# 1. conditions ####
conditions <- list(cv = c("mse"),
                   dimension = c("low", "high"),
                   signal_level =  c(0.2, 0.5, 0.8),
                   py_pattern = c(1,2),
                   reps = c(1:20))

condition_df <- data.frame(cv = NA,
                           dimension = NA,
                           signal_level = NA,
                           py_pattern = NA,
                           reps = NA)

counts <- 0
  for (cvz in 1:1){
    for (dimensionz in 1:2){
      for (signal_levelz in 1:3){
          for (patternz in 1:2){
            for (repsz in 1:20){
              
              counts <- counts + 1
              
              cv_now <- conditions$cv[cvz]
              dimension_now <- conditions$dimension[dimensionz]
              signal_level_now <- conditions$signal_level[signal_levelz]
              py_pattern_now <- conditions$py_pattern[patternz]
              reps_now <- conditions$reps[repsz]
              
              condition_df[counts,] <- c(cv_now,
                                         dimension_now,
                                         signal_level_now,
                                         py_pattern_now,
                                         reps_now)
              
              print(counts) 
            }
          }
      }
    }
  }

