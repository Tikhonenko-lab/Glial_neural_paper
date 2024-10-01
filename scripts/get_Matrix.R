#Written by Pamela Mishra/Sept 2024
library(data.table)
library(dplyr)

source ("Modules_make_matrix.R")

#AFter making the cancer matrix, filter and count percentage of samples have the microexon AFTER FILTERING.
#FIlters the cancer_matrix again
#df is the cancer matrix. Can be NBL, HGG, LGG etc
run_all<-function(df){

  df <- df[df$IJC_SAMPLE_1 + df$SJC_SAMPLE_1 >=50,]
  df<-df %>% group_by(Microexon_full) %>% mutate(Num_times_microexon_full_after_filter= n())
  df<- find_percent_samples_after_filter(df)
  return (df)
}

#START HERE: 
#Read the original splicing file. THis takes time and memroy.
#It will not the read the file if, there is not enough RAM
#Need to add columns and rows iun future
SE_matrix <- read_only_SE()

#Make matrix for the cancer type
NBL <- make_cancer_matrix(SE_matrix, "NBL")
#Filters the cancer matrix by numnber of total reads and 
NBL <- run_all(NBL)
#Select 90% of NBL samples have the microexon_full.
NBL_90 <- NBL[NBL$Percent_of_samples >=90,]
