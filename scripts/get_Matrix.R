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

#fOR THE HEATMAP
h1<-NBL_90 %>% select(Kids_First_Biospecimen_ID, IncLevel1, Microexon_full)
h2<- return_pivot(h1)
h2<- h2 %>%  tibble::column_to_rownames("Microexon_full")
h2<-as.matrix(h2)
#if samples have more NA's then remove those NAs (diffrent module)

pheatmap(h2, color= viridis(50),show_colnames = FALSE, 
         na_col = "grey" ,fontsize_row = 2.1,
          treeheight_col = 0 , fontsize=8, 
         main= "NBL microexons <=30,div by 3,90% samples have row microexon")
