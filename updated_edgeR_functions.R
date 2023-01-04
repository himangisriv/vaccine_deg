
library(tidyverse)
library(RColorBrewer)  
library(limma)
library(edgeR) 
library(RColorBrewer)
library(dplyr)
library('variancePartition')
library('edgeR')
library('BiocParallel')







data_file<-'X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_counts.csv'
treatment_file<-'trt_pubid_2022-DEC-19.csv' # Full Path Please
sample_meta_data<-'X:/fast/gilbert_p/fg_data/SCRI/TBVPX-203/RNA/2019Dec/salmon_metadata.csv'
day_subset<-c("0", "3")
treatment_gp<- c("2 µg ID93 + 5 µg GLA-SE (2 Vaccine Injections)", 
                 "2 µg ID93 + 5 µg GLA-SE (3 Vaccine Injections)")




analysis_type=c("limaa_voom","limma_voom_duplicate_correlation","dream")
boolen<-c(TRUE,FALSE,FALSE)
final_analyis_type<-analysis_type[which(unlist(boolen)==TRUE)]



inputs <- input_argument_list %>% 
  
  .[c("data_file","sample_file","treatment_file","treatment_group_1",
      "treatment_group_2","day_subset_1","day_subset_2",
      "final_analysis_type")]



input_argument_list<-c(data_file,
                       sample_meta_data,
                       treatment_file,treatment_gp,day_subset,final_analyis_type)

names(input_argument_list) <- c("data_file","sample_file","treatment_file",
                                "treatment_group_1","treatment_group_2",
                                "day_subset_1","day_subset_2",
                                "final_analysis_type")




# Create directories for results. 
write_output_directory_name<-function(input_argument_list)
{
  output_analysis_dir =input_argument_list %>% 
    .[c("final_analysis_type","day_subset_2", "day_subset_1")] %>%
    stringr::str_c(.,collapse = "_")
  
  

  write_results_dir = map(output_analysis_dir, 
                          dir.create,
                          recursive = TRUE, 
                          showWarnings = TRUE)
}






# Create filenames for CSV for results.
write_output_file_name <- function(input_argument_list)
{
  filename<-c()
  filename<-input_argument_list %>% 
    .[c("final_analysis_type","day_subset_2", "day_subset_1")] %>%
    stringr::str_c(.,collapse = "_") %>%
    paste(.,".csv",sep = "")
}





data_filtering <- function(count_matrix)
{
  cpm <- cpm(count_matrix)
  lcpm <- cpm(count_matrix, log=TRUE)
  # fitering  genes
  keep <- rowSums(cpm(count_matrix)>1) >= 15
  d0 <- count_matrix[keep,]
  nsamples <- ncol(d0)
  col <- brewer.pal(nsamples, "Paired")
  par(mfrow=c(1,2))
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="A. Raw data", xlab="Log-cpm")
  lcpm <- cpm(d0, log=TRUE)
  plot(density(lcpm[,1]), col=col[1], lwd=2, ylim=c(0,0.26), las=2, main="", xlab="")
  title(main="B. Filtered data", xlab="Log-cpm")
  
  return(d0)
}

data_normalization_edgeR <- function(d0)
{
  # creating a DGE object
  d0 <- DGEList(d0)
  d0 <- calcNormFactors(d0,method = "TMM")
  return(d0)
  
}






data_preprocessing <- function(input_argument_list) 
  
{
  
  inputs <- input_argument_list %>% 
    
    .[c("data_file","sample_file","treatment_file","treatment_group_1",
        "treatment_group_2","day_subset_1","day_subset_2",
        "final_analysis_type")]
  
  
  
  
  data <- readr::read_csv(inputs["data_file"],
                          col_names = TRUE) %>% 
    tibble::column_to_rownames(var = "gene_id")
  
  trt <- readr::read_csv(treatment_file)
  # Original treatment file: column name is Treatment.Group. 
  # When reading readr::read_csv() it reads `treatment group` 
  # See comments about changing input
  
  trt <- trt[trt$Treatment_Group %in% inputs[c("treatment_group_1",
                                               "treatment_group_2")], ]
  # will be `subject_id` after fixing input file
  trt$Subject_ID<-as.character(trt$Subject_ID)
  
  
  
  sample <- readr::read_csv(sample_meta_data,
                            col_names = TRUE) %>% 
    dplyr::select(c("samplename","ptid","day")) %>% 
    dplyr::filter(day %in% inputs[c("day_subset_1","day_subset_2")]) %>%
    dplyr::mutate(Subject_ID = stringr::str_remove_all(ptid,"-")) %>%
    dplyr::mutate(Subject_ID = stringr::str_remove_all(Subject_ID,"203"))
  
  
  
  
  n_row_trt = dim(trt)[1]
  
  meta_data <- dplyr::left_join(trt, 
                                sample,
                                by = "Subject_ID") %>%
    tibble::column_to_rownames(var = "samplename")
  
  
  n_row_meta_data = dim(meta_data)[1]
  # stopifnot(n_row_trt == n_row_meta_data)
  # might check sum(is.na(meta_data['columname_that_was_in_sample_df']) == 0
  
  # Add a comment: <idx2> is ...
  idx2 <- match(rownames(meta_data), colnames(data))
  count_matrix<- data[,idx2]
  count_matrix<-as.matrix(count_matrix[,colSums(is.na(count_matrix))==0])
  
  
  result_data<-list(data,sample,trt,meta_data,count_matrix)
  names(result_data) <- c("original_count_data",
                          "original_sample_meta_data",
                          "subset_treatment_data",
                          "subset_meta_data",
                          "final_count_matrix")
  
  
  return(result_data)
}







#constructing a design matrix, voom transformation and fitting the data
fitting_Edge_R_data<- function(input_argument_list)
{
  inputs <- input_argument_list %>% 
    
    .[c("data_file","sample_file","treatment_file","treatment_group_1",
        "treatment_group_2","day_subset_1","day_subset_2",
        "final_analysis_type")]
  
  
  ty<-data_preprocessing(input_argument_list)
  data1<-data_filtering(ty$final_count_matrix)
  d_ss<-data_normalization_edgeR(data1)
  day_s=as.character(ty$subset_meta_data$day)
  meta_data<-ty$subset_meta_data
  pub_id<-as.character(ty$subset_meta_data$pubid)
  
  
  
  if(inputs["final_analysis_type"]== "limaa_voom" )
  {
    
    design <- model.matrix(~day_s+pub_id, d_ss$samples)
    v <- voom(d_ss, design, plot=TRUE)
    fit <- lmFit(v, design)
    fit_test <- eBayes(fit)
    dt <- decideTests(fit_test)
    result<-topTreat(fit_test,coef='day_s3')
  }
  else if (inputs["final_analysis_type"]== "limma_voom_duplicate_correlation" )
  {
    
    design <- model.matrix(~day_s, d_ss$samples)
    v <- voom(d_ss, design, plot=TRUE)
    corfit <- duplicateCorrelation(v, design, block=meta_data$pubid)
    v <- voom(d_ss, design, block=meta_data$pubid, correlation = corfit$consensus)
    fit <- lmFit(v, design, block=meta_data$pubid, correlation = corfit$consensus)
    fit_test <- eBayes(fit)
    dt <- decideTests(fit_test)
    result<-topTreat(fit_test,coef='day_s3', number=3 )
  }
  return(result)
}











