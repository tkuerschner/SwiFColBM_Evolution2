
# title: "source_lib"
# author: "Tobias Kuerschner"
# date: "2 December 2019"


### Packages



for (pckg in
     c
     (
       'tidyverse',
       'nlme',
       'viridis',
       'scico',
       'ggeffects',
       'psd',
       'TSstudio',
       'zoo',
       'data.table',
       'ggpubr',
       'TSstudio',
       'MullerPlot',
       'ggmuller',
       'wsyn',
       'png',
       'gris',
       'patchwork',
       'lme4'
       
     ))
{
  if (!require(pckg, character.only = T))
    install.packages(pckg, dependencies = TRUE)
  require(pckg, character.only = T)
}
