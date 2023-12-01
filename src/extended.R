dosage_plot = function(df,compound) {

  ggplot(data = df,
         mapping = aes(y = thc, x = compound, color = treatment)) +
    geom_point() +
    geom_smooth(method = "lm", formula = y ~ x)
  

}