#helper function: clean up name of two compounds
clean_gluc <- function(df){
  df <- df |> 
    mutate(compound=gsub('GLUC', 'gluc',gsub("_","-",toupper(compound))),
           compound=gsub('THCOH', '11-OH-THC', compound))
  return(df)
}



#helper function: calculate sensitivity & specificity
make_calculations <- function(dataset, dataset_removedups, split, compound, 
                              start = start, stop = stop, tpt_use = tpt_use){
  ## remove NAs
  df <- dataset_removedups %>% 
    dplyr::select(treatment, compound, timepoint_use) %>%
    rename(compound = 2) %>%
    filter(complete.cases(.))
  if(nrow(df)>0){
    if(stop <= 0){
      output <- df %>% 
        summarise(TP = 0,
                  FN = 0,
                  FP = sum(compound >= split),
                  TN = sum(compound < split)) 
    }else{
      if(split == 0){
        output_pre <- df %>% 
          filter(tpt_use == "pre-smoking") %>%
          summarise(TP = 0,
                    FN = 0,
                    FP = sum(compound >= split),
                    TN = sum(compound < split)) 
        
        output <- df %>% 
          filter(tpt_use != "pre-smoking") %>%
          summarise(TP = sum(treatment != "Placebo" & compound > split),
                    FN = sum(treatment != "Placebo" & compound <= split),
                    FP = sum(treatment == "Placebo" & compound > split),
                    TN = sum(treatment == "Placebo" & compound < split))
        
        output <- output + output_pre
      }else{
        ## calculate values if pre-smoking
        output_pre <- df %>% 
          filter(tpt_use == "pre-smoking") %>%
          summarise(TP = 0,
                    FN = 0,
                    FP = sum(compound >= split),
                    TN = sum(compound < split)) 
        
        output <- df %>% 
          filter(tpt_use != "pre-smoking") %>%
          summarise(TP = sum(treatment != "Placebo" & compound >= split),
                    FN = sum(treatment != "Placebo" & compound < split),
                    FP = sum(treatment == "Placebo" & compound >= split),
                    TN = sum(treatment == "Placebo" & compound < split))
        
        output <- output + output_pre
      }
    }
  }
  # clean things up; make calculations on above values
  output <- output %>%
    mutate(detection_limit = split,
           compound = compound,
           time_start = start,
           time_stop = stop,
           time_window = tpt_use,
           NAs = nrow(dataset) - nrow(df),
           N = nrow(dataset_removedups),
           N_removed = nrow(dataset) - nrow(dataset_removedups),
           Sensitivity = (TP/(TP + FN)), 
           Specificity = (TN /(TN + FP)),
           PPV = (TP/(TP+FP)),
           NPV = (TN/(TN + FN)),
           Efficiency = ((TP + TN)/(TP + TN + FP + FN))*100
    )
  
  return(output)
}


#helper function: apply make_calculations() based on specific cutoffs
sens_spec <- function(dataset, compound, start, stop, tpt_use, 
                      lowest_value = 0.5, splits = NULL, ...){
  # if it's not all NAs...
  if(sum(is.na(dataset[,compound])) != nrow(dataset)){
    # specify what splits should be used for calculations
    if(is.null(splits)){
      limits <- dataset[is.finite(rowSums(dataset[,compound])),compound] #sums up values in the specified compound col, unless its infinite
      ## define lower and upper limits
      lower = min(limits, na.rm=TRUE)
      upper = max(limits, na.rm=TRUE)
      ## determine splits to use for calculations
      tosplit = pull(limits[,1])[limits[,1]>0]
      ## only split if there are detectable limits:
      if(length(tosplit)>=1){
        splits = c(lowest_value, quantile(tosplit, probs=seq(0, 1, by = 0.01), na.rm=TRUE))
        splits = unique(splits)
      #splits: a vector containing lowest_value, and values of 100 percentiles
      }else{
        splits = 0
      }
    }else{
      splits = splits
    }
    # filter to include timepoint of interest
    dataset <- dataset %>% 
      filter(time_from_start > start & time_from_start <= stop & !is.na(timepoint_use))
    dataset_removedups <- dataset %>%
      filter(!is.na(timepoint_use)) %>% 
      group_by(timepoint_use) %>% 
      distinct(id, .keep_all = TRUE) %>% 
      ungroup()
    
    ## create empty output variable which we'll fill in
    ## iterate through each possible dose and calculate
    output <- map_dfr(as.list(splits), ~make_calculations(dataset, 
                                                          dataset_removedups, 
                                                          split = .x,
                                                          compound,
                                                          start = start,
                                                          stop = stop, 
                                                          tpt_use = tpt_use))
  }
  
  return(output)
}




#passes parameters into sens_spec() and return a data frame
#params:
  #dataset:WB/OF/BR
  #cpd: string specifying compound
  #timepoints: timepoints data frame
  #splits: list of cutoffs
#returns data frame with columns:
  #detection_limit = cutoff,
  #compound = compound,
  #time_start = start of timepoint,
  #time_stop = end of timepoint,
  ##time_window = timepoint,
  #NAs = nrow(dataset) - nrow(df),
  #N = nrow(dataset_removedups),  
  #N_removed = nrow(dataset) - nrow(dataset_removedups),
  #Sensitivity = (TP/(TP + FN)), 
  #Specificity = (TN /(TN + FP)),
  #PPV = (TP/(TP+FP)),
  #NPV = (TN/(TN + FN)),
  #Efficiency = ((TP + TN)/(TP + TN + FP + FN))*100

sens_spec_cpd <- function(dataset, cpd, timepoints, splits = NULL){
  args2 <- list(start = timepoints$start, 
                stop = timepoints$stop, 
                tpt_use = timepoints$timepoint)
  out <- args2 %>% 
    pmap_dfr(sens_spec, dataset, compound = cpd, splits = splits)
  return(out)
}


#parameter:
#output: output_WB/OF/BR
average_sens_spec = function(output) {
  output |>
    group_by(compound, detection_limit) |>
    summarize(average_sensitivity = mean(Sensitivity, na.rm = TRUE), average_specificity = mean(Specificity, na.rm = TRUE)) |>
    ungroup(compound) 
}






#sens_spec() but with limits on the splits. only used once. 
sens_spec_OFTHC <- function(dataset, compound, start, stop, tpt_use, 
                            lowest_value = 0.5, splits = NULL, ...){
  # if it's not all NAs...
  if(sum(is.na(dataset[,compound])) != nrow(dataset)){
    # specify what splits should be used for calculations
    if(is.null(splits)){
      limits <- dataset[is.finite(rowSums(dataset[,compound])),compound] #sums up values in the specified compound col, unless its infinite
      ## define lower and upper limits
      lower = 0
      upper = 2
      ## determine splits to use for calculations
      tosplit = pull(limits[,1])[limits[,1]>0]
      ## only split if there are detectable limits:
      if(length(tosplit)>=1){
        splits = c(lowest_value, quantile(tosplit, probs=seq(0, 1, by = 0.01), na.rm=TRUE))
        splits = unique(splits)
        #splits: a vector containing lowest_value, and values of 100 percentiles
      }else{
        splits = 0
      }
    }else{
      splits = splits
    }
    # filter to include timepoint of interest
    dataset <- dataset %>% 
      filter(time_from_start > start & time_from_start <= stop & !is.na(timepoint_use))
    dataset_removedups <- dataset %>%
      filter(!is.na(timepoint_use)) %>% 
      group_by(timepoint_use) %>% 
      distinct(id, .keep_all = TRUE) %>% 
      ungroup()
    
    ## create empty output variable which we'll fill in
    ## iterate through each possible dose and calculate
    output <- map_dfr(as.list(splits), ~make_calculations(dataset, 
                                                          dataset_removedups, 
                                                          split = .x,
                                                          compound,
                                                          start = start,
                                                          stop = stop, 
                                                          tpt_use = tpt_use))
  }
  
  return(output)
}

#again, only used once
sens_spec_cpd_OFTHC <- function(dataset, cpd, timepoints, splits = NULL){
  args2 <- list(start = timepoints$start, 
                stop = timepoints$stop, 
                tpt_use = timepoints$timepoint)
  out <- args2 %>% 
    pmap_dfr(sens_spec, dataset, compound = cpd, splits = seq(0, 2, by = 0.02))
  return(out)
}


