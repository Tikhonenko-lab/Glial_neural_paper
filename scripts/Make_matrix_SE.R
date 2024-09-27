library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(viridis)
library(purrr)
library(ggrepel)


#Find 20% of the samples containing microexons
#Input : Dataframe for NBL/HGG/DMG etc.
find_percent_samples <- function(df){
  uniq<- as.numeric(length(unique(df$Kids_First_Biospecimen_ID)))
  df<- df %>% mutate(Percent_of_samples = (Num_of_samples/uniq) *100  )
  return(df)

}


#add size of the microexon
#Input: Dataframe for NBL/HGG/DMG etc 
add_size_microexon<- function(df_micro){

  df_micro <-df_micro %>% mutate(Size_micro = exonEnd - exonStart_0base)
  return (df_micro)
}

