#remove all duplicate id
drop_dups <- function(dataset){
  out <- dataset |> 
    filter(!is.na(timepoint_use)) |> 
    group_by(timepoint_use) |> 
    distinct(id, .keep_all = TRUE) |> 
    ungroup()
  return(out)
} 



#plot value of compound across time, coloured by treatment
#signature:
  #map(compounds_WB/compounds_BR/compounds_OF, ~ compound_scatterplot_group( 
  #dataset=BR/OF/BR, 
  #compound=.x, 
  #timepoints=timepoints_WB/timepoints_BR/timepoints_OF))
compound_scatterplot_group_by_treatment <- function(dataset, compound, timepoints){
  if(max(dataset[,compound],na.rm=TRUE)==0){
    print(
      dataset |> 
        filter(!is.na(time_from_start)) |>
        ggplot(aes_string(x="time_from_start", 
                          y=compound,
                          color="treatment")) + 
        geom_point() +
        geom_vline(data=timepoints, aes(xintercept=as.numeric(stop)), 
                   linetype="dashed", 
                   color="gray28") +
        scale_color_manual(values=c("#19831C", "#A27FC9","#A27FC9")) +
        scale_y_continuous(limits=c(0,3)) +
        theme_classic() +
        theme(legend.position="bottom",
              legend.title=element_blank()) +
        labs(x='Time From Start (min)',
             y=gsub('GLUC', 'gluc',gsub("_", "-", toupper(compound))))
    )}else{
      print(
        dataset |> 
          filter(!is.na(time_from_start)) |>
          ggplot(aes_string(x="time_from_start", 
                            y=compound,
                            color="treatment")) + 
          geom_point() +
          geom_vline(data=timepoints, aes(xintercept=as.numeric(stop)), 
                     linetype="dashed", 
                     color="gray28")  +
          scale_color_manual(values=c("#19831C", "#A27FC9","#A27FC9")) +
          theme_classic() +
          theme(legend.position="bottom",
                legend.title=element_blank()) +
          labs(x='Time From Start (min)',
               y=gsub('GLUC', 'gluc', gsub("_", "-", toupper(compound))))
      )
    }
}


#plot value of compound across time colour by treatment, log transformed
#signature:
#map(compounds_WB/compounds_BR/compounds_OF, ~ compound_scatterplot_group( 
#dataset=BR/OF/BR, 
#compound=.x, 
#timepoints=timepoints_WB/timepoints_BR/timepoints_OF))
compound_scatterplot_group_by_treatment_log <- function(dataset, compound, timepoints){
  if(max(dataset[,compound],na.rm=TRUE)==0){
    print(
      dataset |> 
        filter(!is.na(time_from_start)) |>
        mutate(log_transformed_compound = log(get(compound)))|>
        ggplot(aes_string(x="time_from_start", 
                          y="log_transformed_compound",
                          color="treatment")) + 
        geom_point() +
        geom_vline(data=timepoints, aes(xintercept=as.numeric(stop)), 
                   linetype="dashed", 
                   color="gray28") +
        scale_color_manual(values=c("#19831C", "#A27FC9","#A27FC9")) +
        scale_y_continuous(limits=c(0,3)) +
        theme_classic() +
        theme(legend.position="bottom",
              legend.title=element_blank()) +
        labs(x='Time From Start (min)',
             y=gsub('GLUC', 'gluc',gsub("_", "-", toupper(compound))))
    )}else{
      print(
        dataset |> 
          filter(!is.na(time_from_start)) |>
          mutate(log_transformed_compound = log(get(compound)))|>
          ggplot(aes_string(x="time_from_start", 
                            y="log_transformed_compound",
                            color="treatment")) + 
          geom_point() +
          geom_vline(data=timepoints, aes(xintercept=as.numeric(stop)), 
                     linetype="dashed", 
                     color="gray28")  +
          scale_color_manual(values=c("#19831C", "#A27FC9","#A27FC9")) +
          theme_classic() +
          theme(legend.position="bottom",
                legend.title=element_blank()) +
          labs(x='Time From Start (min)',
               y=gsub('GLUC', 'gluc', gsub("_", "-", toupper(compound))))
      )
    }
}




