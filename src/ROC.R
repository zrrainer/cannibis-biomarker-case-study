#plot detection limit(cutoff) vs. sensitivity and specificity
#parameters:
  #output: output_OF/WB/BR from a map_fdr sens_spec_cpd() function
  #tpts: length(unique(output_WB/OF/BR$time_start))
  #tissue: string specifying which matrix
#returns: 
  #two plots - dettection limit (cutoff) vs. sensitivity and specificity
ss_plot <- function(output, tpts=8, tissue){
  to_include = output %>%
    group_by(compound) %>% 
    summarize(mean_detection = mean(detection_limit)) %>% 
    filter(mean_detection > 0)
  
  output <-  output %>% 
    mutate(iszero = ifelse(time_start<0,TRUE,FALSE),
           Sensitivity = round(Sensitivity*100,0),
           Specificity = round(Specificity*100,0)) %>%
    filter(compound %in% to_include$compound,
           time_window != "pre-smoking") %>%
    clean_gluc() %>% 
    mutate(compound = fct_relevel(as.factor(compound), "THC"))
  
  output <- output %>%  mutate(
    legend = paste0(time_window,' (N=', N,')'))
  
  blue_colors = c('#C2F8FF', '#A2DDED', '#86BEDC', '#6C9FCA', 
                  '#547EB9', '#3F5EA8', '#2D4096', '#1E2385',
                  '#181173', '#180762', '#180051')
  values = c(blue_colors[1:tpts])
  
  print(ggplot(output, aes(x = detection_limit, y = Sensitivity, group = fct_inorder(legend))) + 
          geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
          geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
          facet_grid(~compound, scales = "free_x") +
          labs(x = 'Detection Limit',
               y = 'Sensitivity') +
          ylim(0,100) +
          scale_color_manual(values = values, name = 'Time Window') +
          theme_classic(base_size = 12) + 
          theme(axis.title = element_text(size=16), 
                panel.grid = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_text(size = 12))  
  )
  print(
    ggplot(output, aes(x = detection_limit, y = Specificity, group = fct_inorder(legend))) + 
      geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
      geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
      facet_grid(~compound, scales = "free_x") +
      ylim(0,100) +
      labs(title = tissue,
           x = 'Detection Limit',
           y = 'Specificity') +
      scale_color_manual(values = values, name = 'Time Window') +
      theme_classic(base_size = 12) + 
      theme(axis.title = element_text(size=16),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12))
  )
}


#exact same function as ss_plot() but has xmin, xmax parameter
ss_plot_xscaled <- function(output, tpts=8, tissue, xmin, xmax){
  to_include = output %>%
    group_by(compound) %>% 
    summarize(mean_detection = mean(detection_limit)) %>% 
    filter(mean_detection > 0)
  
  output <-  output %>% 
    mutate(iszero = ifelse(time_start<0,TRUE,FALSE),
           Sensitivity = round(Sensitivity*100,0),
           Specificity = round(Specificity*100,0)) %>%
    filter(compound %in% to_include$compound,
           time_window != "pre-smoking") %>%
    clean_gluc() %>% 
    mutate(compound = fct_relevel(as.factor(compound), "THC"))
  
  output <- output %>%  mutate(
    legend = paste0(time_window,' (N=', N,')'))
  
  blue_colors = c('#C2F8FF', '#A2DDED', '#86BEDC', '#6C9FCA', 
                  '#547EB9', '#3F5EA8', '#2D4096', '#1E2385',
                  '#181173', '#180762', '#180051')
  values = c(blue_colors[1:tpts])
  
  print(ggplot(output, aes(x = detection_limit, y = Sensitivity, group = fct_inorder(legend))) + 
          geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
          geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
          facet_grid(~compound, scales = "free_x") +
          labs(x = 'Detection Limit',
               y = 'Sensitivity') +
          ylim(0,1) +
          scale_color_manual(values = values, name = 'Time Window') +
          theme_classic(base_size = 12) + 
          theme(axis.title = element_text(size=16), 
                panel.grid = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_text(size = 12)) +
          scale_x_continuous(limits = c(xmin,xmax))
  )
  print(
    ggplot(output, aes(x = detection_limit, y = Specificity, group = fct_inorder(legend))) + 
      geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
      geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
      facet_grid(~compound, scales = "free_x") +
      ylim(0,100) +
      labs(title = tissue,
           x = 'Detection Limit',
           y = 'Specificity') +
      scale_color_manual(values = values, name = 'Time Window') +
      theme_classic(base_size = 12) + 
      theme(axis.title = element_text(size=16),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12))+
      scale_x_continuous(limits = c(xmin,xmax))
  )
}

#plot detection limit vs. average sens/spec

#ss_plot() repurposed to plot average sens/spec
ss_plot_avg <- function(output_avg, tpts=8, tissue){
  to_include = output_avg %>%
    group_by(compound) %>% 
    summarize(mean_detection = mean(detection_limit)) %>% 
    filter(mean_detection > 0)
  
  output_avg <-  output_avg %>% 
    mutate(average_sensitivity = round(average_sensitivity*100,0),
           average_specificity = round(average_specificity*100,0)) %>%
    clean_gluc() %>% 
    mutate(compound = fct_relevel(as.factor(compound), "THC"))
  
  blue_colors = c('#C2F8FF', '#A2DDED', '#86BEDC', '#6C9FCA', 
                  '#547EB9', '#3F5EA8', '#2D4096', '#1E2385',
                  '#181173', '#180762', '#180051')
  values = c(blue_colors[1:tpts])
  #####
  print(ggplot(output_avg, aes(x = detection_limit, y = average_sensitivity)) + 
          geom_point(size = 0.9, show.legend = FALSE) +
          geom_path(size=1.2)+
          facet_grid(~compound, scales = "free_x") +
          labs(x = 'Detection Limit',
               y = 'Average Sensitivity') +
          ylim(0,100) +
          scale_color_manual(values = values, name = 'Time Window') +
          theme_classic(base_size = 12) + 
          theme(axis.title = element_text(size=16), 
                panel.grid = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_text(size = 12))  
  )
  print(ggplot(output_avg, aes(x = detection_limit, y = average_specificity)) + 
          geom_point(size = 0.9, show.legend = FALSE) +
          geom_path(size=1.2)+
          facet_grid(~compound, scales = "free_x") +
          labs(x = 'Detection Limit',
               y = 'Average Sensitivity') +
          ylim(0,100) +
          scale_color_manual(values = values, name = 'Time Window') +
          theme_classic(base_size = 12) + 
          theme(axis.title = element_text(size=16), 
                panel.grid = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_text(size = 12))  
  )
}

#plot detection limit vs. average sens/spec on the same graph
#parameter
  #output_avg: df from average_sens_spec()
  #tpts: length(unique(output_WB$time_start))
  #xmin, xmax: limit on the x axis
ss_plot_avg_together <- function(output_avg, tpts=8, tissue){
  to_include = output_avg %>%
    group_by(compound) %>% 
    summarize(mean_detection = mean(detection_limit)) %>% 
    filter(mean_detection > 0)
  
  output_avg <-  output_avg %>% 
    mutate(average_sensitivity = round(average_sensitivity*100,0),
           average_specificity = round(average_specificity*100,0)) %>%
    clean_gluc() %>% 
    mutate(compound = fct_relevel(as.factor(compound), "THC"))
  
  blue_colors = c('#ac4e57', '#7496a4')
  values = c(blue_colors[1:tpts])
  #####
  print(ggplot(output_avg, aes(x = detection_limit, y = average_sensitivity)) + 
          geom_point(size = 0.9, show.legend = FALSE) +
          geom_path(aes(x = detection_limit, y = average_sensitivity, color = "average sensitivity"), size=1.2) + 
          geom_path(aes(x = detection_limit, y = average_specificity, color = "average specificity"), size=1.2) + 
          scale_y_continuous(name = "average sensitivity", sec.axis = sec_axis(~./1, name = "average specificity")) +
          facet_grid(~compound, scales = "free_x") +
          labs(x = 'Detection Limit',
               y = 'Average Sensitivity/Specificity') +
          ylim(0,100) +
          scale_color_manual(values = values, name = 'Time Window') +
          theme_classic(base_size = 12) + 
          theme(axis.title = element_text(size=16), 
                panel.grid = element_blank(),
                strip.background = element_blank(),
                strip.text.x = element_text(size = 12)) 
  )
}






#parameters:
  #output: output_OF/WB/BR from a map_fdr sens_spec_cpd() function
  #tpts: length(unique(output_WB/OF/BR$time_start))
  #tissue: string specifying which matrix
#returns: plot sensitivity vs. specificity
roc_plot <- function(output, tpts=8, tissue){
  to_include = output %>%
    group_by(compound) %>% 
    summarize(mean_detection = mean(detection_limit)) %>% 
    filter(mean_detection > 0)
  
  output <-  output %>% 
    mutate(iszero = ifelse(time_start<0,TRUE,FALSE),
           Sensitivity = round(Sensitivity*100,0),
           Specificity = round(Specificity*100,0)) %>%
    filter(compound %in% to_include$compound,
           time_window != "pre-smoking") %>%
    clean_gluc() %>% 
    mutate(compound = fct_relevel(as.factor(compound), "THC"))
  
  output <- output %>% mutate(
    legend = paste0(time_window,' (N=', N,')'))
  
  blue_colors = c('#C2F8FF', '#86BEDC', 
                  '#547EB9', '#2D4096',
                  '#181173', '#180051')
  values = c(blue_colors[1:tpts])
  print(
    ggplot(output, aes(x=(100-Specificity), y = Sensitivity, group = fct_inorder(legend))) +
      geom_point(aes(color=fct_inorder(legend)), size = 0.9, show.legend = FALSE) +
      geom_path(aes(color=fct_inorder(legend)), size=1.2) + 
      facet_grid(~compound) +
      xlim(0, 100) +
      ylim(0, 100) +
      labs(title = tissue,
           x = '(100-Specificity)',
           y = 'Sensitivity') +
      scale_color_manual(values = values, name = 'Time Window') +
      theme_classic(base_size = 12) + 
      theme(axis.title = element_text(size=16),
            panel.grid = element_blank(),
            strip.background = element_blank(),
            strip.text.x = element_text(size = 12),
            axis.text = element_text(size=12) )
  )
}