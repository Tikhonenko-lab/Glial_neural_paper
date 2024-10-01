### Created by: Pamela Mishra/Sept/2024/
## Function: To run regular repetitive tasks on OpenPedCan V13 data

library(dplyr)
library(tidyr)
library(data.table)
library(pheatmap)
library(viridis)
library(purrr)
library(ggrepel)


#Find 20% of the samples containing microexons
#Input : Dataframe for NBL/HGG/DMG etc.
find_percent_samples_after_filter <- function(df){
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
}

#Function to read the data from OpenPedCan and extract only SE 
#Extract only RNA-Seq samples and non-NA cancer types
read_only_SE <- function(){

      SE<- fread("../data/OpenPedCan-analysis/data/v13/histologies.tsv")
      
      #make a new df with a reduced number of columns
      cols<- c("splicing_case", "geneSymbol", "sample_id", "chr", "strand", "exonStart_0base", "exonEnd"  , "upstreamES",
         "upstreamEE",  "downstreamES", "downstreamEE", "IJC_SAMPLE_1", "SJC_SAMPLE_1","IncLevel1")

      SE_trunc<- SE %>% select(all_of(cols))

      #Take histology data in
      hist <- fread("../data/OpenPedCan-analysis/data/v13/histologies.tsv")
      SE_trunc <- SE_trunc %>% filter(splicing_case == "SE")

      hist<- hist %>% filter(experimental_strategy =="RNA-Seq")
      #hist<- hist %>% filter(cancer_group !="NA")

      SE_trunc_hist<- left_join(SE_trunc, hist, by= "Kids_First_Biospecimen_ID" )

      SE_trunc_final <- SE_trunc_hist %>% select(all_of(cols), cancer_group)
      #Then check for NAs
      any(is.na(SE_trunc_final))

      #Some cancer_groups have "NA" in their histology
      #Make sure to remove these lines
      SE_trunc_final<- SE_trunc_final %>% filter(! if_any(c("cancer_group"), is.na))
      saveRDS(SE_trunc_final, file="SE_noNA_data.RData")

      return (SE_trunc_final)
}

#Return a pivoted df after selecting a few columsn
return_pivot <- function(df){
  df_new <-df %>% select(IncLevel1, Kids_First_Biospecimen_ID, Microexon)
  j1 <-  as.data.frame(pivot_wider(df_new, names_from = Kids_First_Biospecimen_ID,  values_from = IncLevel1))
  return(j1)
}

#Make for each cancer type
#Specify cancer_type ==  NBL/HGG/LGG etc
make_cancer_matrix <- function(SE_matrix, cancer_type){
         cancer_matrix <- SE_matrix %>% filter(cancer_group= cancer_type)

         cancer_matrix <- cancer_matrix %>% mutate(Microexon_full =paste(geneSymbol,chr,strand,exonStart_0base,exonEnd,upstreamES, upstreamEE,downstreamES, downstreamEE, sep="_"))
         
	 cancer_matrix <- cancer_matrix %>% mutate(Microexon_part =paste(geneSymbol,chr,strand,exonStart_0base,exonEnd, sep="_"))
         
	 #Microexon size
	 cancer_matrix  = add_size_microexon(cancer_matrix)
	 #Filter by microexon size and fdivisible by 3
	 cancer_matrix  = return_filtered_df(cancer_matrix)
	 return (cancer_matrix)
}



