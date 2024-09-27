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


#Find common microexons given a list oSf microexon
#Input : List containing microexon ids
find_common_microexons <- function(list1){
  common_list <- list1 %>% purrr::reduce(intersect)
  return(common_list)
}



#Function to return exon1 and number of the tries == 1
return_filtered_df <-function (df){
  micro <-  df %>% filter(Size_micro <=30) %>% filter(Num_times_microexon == 1 ) %>% filter(Size_micro %% 3 == 0)
  return (micro)

}

#Find the com/mon list in all df
find_common_elements<- function(list1, df){
  new_df<- df %>% filter(Microexon %in% list1)
  return (new_df)


#Function to read the data from OpenPedCan and extract only SE 
read_only_SE <- function(file){

}

#Return a pivoted df after selecting a few columsn
return_pivot <- function(df){
  df_new <-df %>% select(IncLevel1, Kids_First_Biospecimen_ID, Microexon)
  j1 <-  as.data.frame(pivot_wider(df_new, names_from = Kids_First_Biospecimen_ID,  values_from = IncLevel1))
  return(j1)
}

