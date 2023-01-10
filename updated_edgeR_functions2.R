arguments_file<-"X:/fast/gilbert_p/hsrivast/vaccine_1/arguments_file.csv"

library(tidyverse)
library(RColorBrewer)  
library(limma)
library(edgeR) 
library(RColorBrewer)
library(dplyr)
library('variancePartition')
library('edgeR')
library('BiocParallel')
library(DESeq2)


#' Title
#'
#' @param arguments_file 
#'
#' @return
#' @export
#'
#' @examples
extract_pipeline_input_from_input_file <- function(arguments_file)
  {
  
  # Remove non-existant configurations. 
  input_file <- arguments_file %>% 
    read_csv(col_types = cols(),locale = readr::locale(encoding = "latin1")) %>% 
    filter(!(argument %>% is.na & value1 %>% is.na)) 
  
  # Extract argument types. 
  arg_type <- input_file %>% 
    select(argument) %>% 
    unlist(use.names = F)
  
  # Link argument to argument type in a named list data structure.
  pipeline_input <- input_file %>% 
                   select(value1) %>% 
                   unlist(use.names = F) %>% 
                   as.list() %>% 
                   set_names(arg_type)
  
  
  return(pipeline_input)
}




#' Title
#'
#' @param path 
#' @param key 
#'
#' @return
#' @export
#'
#' @examples
write_output_directory_name<-function(path,key)
{
  newpath <- file.path(dirname(path), key)
  
  if (!dir.exists(newpath))
  {
    write_results_dir = map(newpath, 
                            dir.create,
                            recursive = TRUE, 
                            showWarnings = TRUE)
  }else{
    print("directory already exists")
  }
}




#' Title
#'
#' @param pipeline_input 
#'
#' @return
#' @export
#'
#' @examples
analysis_type_to_run <- function(pipeline_input)
  {
  path<-pipeline_input %>% magrittr::extract2('output_directory')
  s<-tibble(limma_voom = pipeline_input %>% magrittr::extract2('limma_voom'),
         limma_voom_duplicate_cor = pipeline_input %>% 
                      magrittr::extract2('limma_voom_duplicate_correlation'))%>% 
    
    gather() %>% 
    filter(value == 'TRUE')
    
    map2(path,s$key,write_output_directory_name)

    
}




#' Title
#'
#' @param arguments_file 
#'
#' @return
#' @export
#'
#' @examples
data_preprocessing <- function(arguments_file) 
  
{
  
  pipeline_input <- extract_pipeline_input_from_input_file(arguments_file)
  
  
  
  
  data <- readr::read_csv(pipeline_input[["data_file"]],
                          col_names = TRUE) %>% 
    tibble::column_to_rownames(var = "gene_id")
  
  trt <- readr::read_csv(pipeline_input[["treatment_file"]])
  # Original treatment file: column name is Treatment.Group. 
  # When reading readr::read_csv() it reads `treatment group` 
  # See comments about changing input
  
  trt <- trt[trt$Treatment_Group %in% strsplit(pipeline_input[["treatment_gp"]],",")[[1]], ]
  # will be `subject_id` after fixing input file
  trt$Subject_ID<-as.character(trt$Subject_ID)
  
  
  
  sample <- readr::read_csv(pipeline_input[['sample_meta_data']],
                            col_names = TRUE) %>% 
    dplyr::select(c("samplename","ptid","day")) %>% 
    dplyr::filter(day %in% strsplit(pipeline_input[['day_subset']],",")[[1]]) %>%
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



#' Title
#'
#' @param count_matrix 
#'
#' @return
#' @export
#'
#' @examples
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

#' Title
#'
#' @param d0 
#'
#' @return
#' @export
#'
#' @examples
data_normalization_edgeR <- function(d0)
{
  # creating a DGE object
  d0 <- DGEList(d0)
  d0 <- calcNormFactors(d0,method = "TMM")
  return(d0)
  
}


#' Title
#'
#' @param res 
#' @param lfcthresh 
#' @param sigthresh 
#' @param main 
#' @param legendpos 
#' @param labelsig 
#' @param textcx 
#' @param ... 
#'
#' @return
#' @export
#'
#' @examples
volcanoplot_edge_R <- function (res, lfcthresh=2, sigthresh=0.05, main="Volcano Plot for 63_vs_3", legendpos="bottomright", labelsig=TRUE, textcx=1, ...) {
  Gene=rownames(res)
  with(res, plot(logFC, -log10(P.Value), pch=20, main=main, ...))
  with(subset(res, adj.P.Val<sigthresh ), points(logFC, -log10(P.Value), pch=20, col="red", ...))
  with(subset(res, abs(logFC)>lfcthresh), points(logFC, -log10(P.Value), pch=20, col="orange", ...))
  with(subset(res, adj.P.Val<sigthresh & abs(logFC)>lfcthresh), points(logFC, -log10(P.Value), pch=20, col="green", ...))
  if (labelsig) {
    require(calibrate)
    with(subset(res, adj.P.Val<sigthresh & abs(logFC)>lfcthresh), points(logFC, -log10(P.Value), labs=Gene, cex=textcx, ...))
  }
  legend(legendpos, xjust=1, yjust=1, legend=c(paste("FDR<",sigthresh,sep=""), paste("|logFC|>",lfcthresh,sep=""), "both"), pch=20, col=c("red","orange","green"))
}




#constructing a design matrix, voom transformation and fitting the data
#' Title
#'
#' @param input_argument_list 
#'
#' @return
#' @export
#'
#' @examples
fitting_data<- function(arguments_file)
{
  pipeline_input <- extract_pipeline_input_from_input_file(arguments_file)
  
  
  ty<-data_preprocessing(arguments_file)
  data1<-data_filtering(ty$final_count_matrix)
  d_ss<-data_normalization_edgeR(data1)
  day_s=as.character(ty$subset_meta_data$day)
  meta_data<-ty$subset_meta_data
  pub_id<-as.character(ty$subset_meta_data$pubid)
  
  
  
  if(pipeline_input[["limma_voom"]]== "True" )
  {
    
    design <- model.matrix(~day_s+pub_id, d_ss$samples)
    v <- voom(d_ss, design, plot=TRUE)
    fit <- lmFit(v, design)
    fit_test <- eBayes(fit)
    dt <- decideTests(fit_test)
    plotMD(fit_test, column=1, status=dt[,2], main=colnames(fit_test)[2], xlim=c(-8,13))
    result<-topTreat(fit_test,coef='day_s3',number=Inf)
    #find out significant genes
    resSigind = result[ which(result$adj.P.Val < 0.1 & result$logFC > 0), ]
    resultSigrep = result[ which(result$padj < 0.1 & result$logFC < 0), ]
    resultSig = rbind(resultSigind, resultSigrep)
    # order the resultult from deseq based on pvalues
    result <- result[order(result$padj), ]
    
    ## Merge with normalized count data
    resultdata <- merge(as.data.frame(result), 
                        as.data.frame(d_ss$counts), 
                        by = "row.names", 
                        sort = FALSE)
    names(resultdata)[1] <- "Gene"
    ## Examine plot of p-values
    hist(result$P.Value, breaks=50, col="grey")
  }
  else if (pipeline_input[["limma_voom_duplicate_correlation"]]== "True" )
  {
    
    design <- model.matrix(~day_s, d_ss$samples)
    v <- voom(d_ss, design, plot=TRUE)
    corfit <- duplicateCorrelation(v, design, block=meta_data$pubid)
    v <- voom(d_ss, design, block=meta_data$pubid, correlation = corfit$consensus)
    fit <- lmFit(v, design, block=meta_data$pubid, correlation = corfit$consensus)
    fit_test <- eBayes(fit)
    dt <- decideTests(fit_test)
    plotMD(fit_test, column=1, status=dt[,2], main=colnames(fit_test)[2], xlim=c(-8,13))
    result<-topTreat(fit_test,coef='day_s3',number=Inf)
    
    #find out significant genes
    resSigind = result[ which(result$padj < 0.1 & result$logFC > 0), ]
    resultSigrep = result[ which(result$adj.P.Val < 0.1 & result$logFC < 0), ]
    resultSig = rbind(resultSigind, resultSigrep)
    # order the resultult from deseq based on pvalues
    result <- result[order(result$adj.P.Val), ]
    
    ## Merge with normalized count data
    resultdata <- merge(as.data.frame(result), 
                        as.data.frame(d_ss$counts), 
                        by = "row.names", 
                        sort = FALSE)
    names(resultdata)[1] <- "Gene"
    ## Examine plot of p-values
    hist(result$P.Value, breaks=50, col="grey")
  }
  else if (pipeline_input[["Deseq"]]== "True" )
  {
    
    meta_data$day <- factor(meta_data$day)
    meta_data$samplename <- rownames(meta_data)
    data_ds<-as.data.frame(ty$final_count_matrix)
    idx <- match(rownames(meta_data), colnames(data_ds))
    data_ds<- data_ds[,idx]
    all(rownames(meta_data) == colnames(data_ds))
    
    dds <- DESeqDataSetFromMatrix(countData = round(data_ds),
                                  colData = meta_data,
                                  design= ~ pubid+day)
    
    #fiter genes
    keep <- rowSums(counts(dds) >= 10) >= 10
    sum(keep)
    dds <- dds[keep,]
    nrow(dds)
    
    
    
    #estimate dispersion in DEseq obect
    dds <- estimateSizeFactors(dds)
    sizeFactors(dds)
    dds = estimateDispersions( dds )
    
    
    
    #plot dispersion estimate
    
    plotDispEsts(dds)
    cds <- DESeq(dds)
    res <- results(cds)
    table(res$padj<0.1)
    
    
    #find out significant genes
    resSigind = res[ which(res$padj < 0.1 & res$log2FoldChange > 0), ]
    resSigrep = res[ which(res$padj < 0.1 & res$log2FoldChange < 0), ]
    resSig = rbind(resSigind, resSigrep)
    rownames(resSigind)
    
    
    # order the result from deseq based on pvalues
    res <- res[order(res$padj), ]
    
    ## Merge with normalized count data
    resdata <- merge(as.data.frame(res), 
                     as.data.frame(counts(dds, normalized = TRUE)), 
                     by = "row.names", 
                     sort = FALSE)
    names(resdata)[1] <- "Gene"
    head(resdata)
    
    
    ## Examine plot of p-values
    hist(res$pvalue, breaks=50, col="grey")
    
  }
  return(result)
}


