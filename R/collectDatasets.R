
#' @export
fetchData <- function(accession, output_dir){
  require("BiocManager")
  require("devtools")
  if (!requireNamespace("BiocManager", quietly = TRUE))
    suppressWarnings(install.packages("BiocManager"))
  require(BiocManager)
  if (!require(GEOquery)) suppressWarnings(BiocManager::install("GEOquery"))
  suppressWarnings(library(GEOquery))

  if (!require(GEOmetadb)) suppressWarnings(BiocManager::install("GEOmetadb"))
  suppressWarnings(library(GEOmetadb))

  if (!require(AnnotationDbi)) suppressWarnings(BiocManager::install("AnnotationDbi"))
  suppressWarnings(library(AnnotationDbi))

  if (!require(stringr)) suppressWarnings(install.packages("stringr"))
  suppressWarnings(library(stringr))

  if(!require(affy)) suppressWarnings(BiocManager::install("affy"))
  suppressWarnings(library(affy))
  if(!require(annotate)) suppressWarnings(BiocManager::install("annotate"))
  suppressWarnings(library(annotate))

  if(!require(org.Hs.eg.db)) suppressWarnings(BiocManager::install("org.Hs.eg.db"))
  suppressWarnings(library(org.Hs.eg.db))

  setwd(output_dir)
  gse <- suppressMessages(getGEO(paste("GSE",accession,sep=""), GSEMatrix = TRUE))
  message("Connected to GSE. Attempting to map the GPL number to the microarray platform.")
  finalData <- map_accession_to_platform(gse, accession)
  return(finalData)

  unlink(paste(output_dir,"/",id,sep=""), recursive = T)

}

map_accession_to_platform <- function(data, accession_id){
  total <- list()
  for(i in 1:length(data)){
    platform_id<-data[[i]]@experimentData@other[["platform_id"]]

    if(length(data)==1){
      message("Only one GPL found in dataset.")
    } else{
      if(i==1){
        message("Multiple GPLs found in dataset. ")
      }
      platform_id <- strsplit(platform_id, "\n")[[1]]
    }
    table <- read.table("https://gist.githubusercontent.com/seandavi/bc6b1b82dc65c47510c7/raw/b2027d7938896dce6145a1ebbcea75826813f6e1/platformMap.txt", header=T)
    if(platform_id[i] %in% table$gpl){
      platform <- table[which(table$gpl==platform_id[i]),"bioc_package"]
      message(paste("Platform used for",platform_id[i],":", platform))
    }
    else{
      message(paste("Platform title not found from the GPL number:", platform_id[i]))
      platform <- "temporary"
    }
    val <- map_platform_to_cdf(data, platform, accession_id, platform_id[i])
    total <- append(total, val)
  }
  return(total)
}

map_platform_to_cdf <- function(data, platform, accession_id, platform_id){
  tryCatch({
    if(length(data)>1){

      species <- unique(data[[paste("GSE", accession_id,"-",platform_id, "_series_matrix.txt.gz",sep = "")]]@phenoData@data[["organism_ch1"]])


    } else{

      species <- unique(data[[paste("GSE", accession_id, "_series_matrix.txt.gz",sep = "")]]@phenoData@data[["organism_ch1"]])
    }

    print(paste("SPECIES: ",species))

    if(length(species)>1){
      message(paste("More than one species detected. CDF to be used:", species[[1]]))
      species<-species[[1]]
    }
    species1 <- strsplit(as.character(species), " ")

    species1 <- species1[[1]]

    new_string_species <- paste(str_sub(species1[1], 0, 1),str_sub(species1[[2]], 0, 1),sep="")
    new_string_species <- str_to_lower(new_string_species, locale="en")
    if(platform=="temporary"){
      if(length(data)>1){
        some_data <-  data[[paste("GSE", accession_id,"-",platform_id, "_series_matrix.txt.gz",sep = "")]]@phenoData@data[["instrument_model"]]

        message(paste("Platform used: ", unique(some_data)))

      } else{
        some_data <- data[[paste("GSE", accession_id, "_series_matrix.txt.gz",sep = "")]]@phenoData@data[["instrument_model"]]

        message(paste("Platform used: ", unique(some_data)))
      }


      cdf_url <- readline(prompt=paste("Could not find cdf name from GPL of dataset, ", platform_id, ".\nPlease provide the link to the CDF file for dataset, or type Q to skip dataset. \nAlternatively type N to proceed without using CDF file.", sep=""))
      if(cdf_url=="Q"){
        return()
      }
    } else{
      cdf_url <- paste("http://mbni.org/customcdf/23.0.0/entrezg.download/",platform,new_string_species,"entrezgcdf_23.0.0.tar.gz",sep="")
    }
    message("Installing CDF: ")
    if(cdf_url!="N"){
      suppressMessages(devtools::install_url(cdf_url))
      cdf_library <- paste(platform,new_string_species,"entrezgcdf",sep="")
      suppressMessages(require(cdf_library, character.only=TRUE))

      message(paste("CDF installed: ",cdf_library, sep=""))
    } else{
      message("CDF not utilized for dataset.")
      cdf_library="N"
    }

    files <- getGEOSuppFiles(paste("GSE",accession_id,sep=""))
    wd <- paste(getwd(),"/GSE",accession_id,sep="")
    setwd(wd)
    return(runCdfCode(data, wd,accession_id,cdf_library, platform_id))
  }
  )
}


runCdfCode <- function(data, wd,accession,cdfname, platform_id){
  message("Unzipping TAR file into individual CEL files...")
  untar(paste("GSE",accession,"_RAW.tar",sep=""), list = FALSE, exdir = paste(wd,"/GSE",accession,"_EXTRACTED",sep=""), extras = NULL, verbose = FALSE, restore_times =  TRUE, support_old_tars = Sys.getenv("R_SUPPORT_OLD_TARS", FALSE), tar = Sys.getenv("TAR"))
  cel_files <- untar(paste("GSE",accession,"_RAW.tar",sep=""), list = TRUE, exdir = paste(wd,"/GSE",accession,"_EXTRACTED",sep=""), extras = NULL, verbose = FALSE, restore_times =  TRUE, support_old_tars = Sys.getenv("R_SUPPORT_OLD_TARS", FALSE), tar = Sys.getenv("TAR"))
  setwd(paste(wd,"/GSE",accession,"_EXTRACTED",sep=""))
  if(accession!="32646"){

    for(filename in cel_files){
      tryCatch({
        gunzip(filename, destname = gsub("[.]gz$", "", filename), overwrite = FALSE,
               remove = TRUE, BFR.SIZE = 1e+07)
      },
      error = function(cond){


      })
    }
  }
  cels_in_dir <- Sys.glob("*.CEL")
  message("Unzipping complete.")
  if(cdfname!="N"){
    message("Converting CEL files based on custom CDF...")
    if(length(cels_in_dir)==0){
      message("No CEL Files found!")
      return("")
    } else{
      return(transform_data(cdfname, accession, data, platform_id, cels_in_dir))
    }
  } else{
    if(length(cels_in_dir)==0){
      message("No CEL Files found!")
      return("")
    } else{
    return(transform_data("NOCDF", accession, data, platform_id, cels_in_dir))
    }
  }
}



transform_data <- function(cdfname, accession, gse, platform_id, filenamesV1){

  #install.packages("devtools")
  #library(devtools)
  #install.packages("richierocks/pathological")
  #library(pathological)


  #BiocManager::install("annotate")
  #require(cdfname, character.only=TRUE)
  full_gse_name <- paste("GSE",accession,sep="")
  if(length(gse)==1){
    filenames <- colnames(gse[[paste(full_gse_name,"_series_matrix.txt.gz", sep="")]]@assayData$exprs)
  } else{
    filenames <- colnames(gse[[paste(full_gse_name,"-",platform_id,"_series_matrix.txt.gz", sep="")]]@assayData$exprs)
  }
 # filenames <- paste(filenames,".CEL",sep="")
  matches <- as.data.frame(sapply(filenames, function (y) sapply(filenamesV1, function (x) grepl(y, x))))
  matches <- matches[which(apply(matches, 1, any)),]
  # for(row in length(rownames(matches))){
  #   if(sum(matches[row,])==0){
  #     matches <- matches[which(rownames(matches))!=rownames[row],]
  #   }
  # }
  filenames <- rownames(matches)
  allFilesInDir <- list.files(getwd())
  newNames<-str_to_upper(allFilesInDir)
  file.rename(allFilesInDir,newNames)
  if(cdfname=="NOCDF"){
    Data<-ReadAffy(filenames=filenames)
  } else{
  Data<-ReadAffy(cdfname = cdfname, filenames=filenames)
  }
  #eset<-rma(Data)
  eset<-mas5(Data)
  ID<-featureNames(eset)
  ID2<-sub("_at","",ID)
  #BiocManager::install("org.Hs.eg.db")

  # require(org.Hs.eg.db)
  GS <- as.matrix(getSYMBOL(ID2, 'org.Hs.eg'))
  ematrix<-exprs(eset)
  rows <- GS
  cols =  c("GeneSymbol",colnames(ematrix))
  ematrix <- cbind(rows,ematrix)
  ematrix <- ematrix[which(ematrix[,1] != "NA"),] #remove NAs
  ematrix <- ematrix[order(ematrix[,1]),] #sort by gene name
  ematrix <- rbind(cols, ematrix)
  dir.create(paste("/Users/anjanbharadwaj/GEOpipeline/GSE",accession,"_",platform_id,sep=""))
  setwd(paste("/Users/anjanbharadwaj/GEOpipeline/GSE",accession,"_",platform_id,sep=""))
  write.table(ematrix,file=paste("NormalizedExpressionArray.customCDF.mas5.txt",sep=""),sep="\t", col.names=F, row.names=F,quote=FALSE)
  exprs <- as.data.frame(ematrix)
  rownames(exprs) <- exprs[,1]
  exprs <- exprs[,2:length(colnames(exprs))]
  exprs <- exprs[2:length(rownames(exprs)),]

  colnames(exprs) <- str_replace_all(colnames(exprs), ".CEL", "")
  colnames(exprs) <- str_replace_all(colnames(exprs), ".cel", "")


  if(length(gse)>1){
  listOfData <- list(exprs, gse[[paste(full_gse_name,"-",platform_id,"_series_matrix.txt.gz", sep="")]]@phenoData@data)
  } else{
    listOfData <- list(exprs,gse[[paste(full_gse_name,"_series_matrix.txt.gz", sep="")]]@phenoData@data)
  }

  indices<- sapply(listOfData[[1]], is.factor)
  listOfData[[1]][indices] <- lapply(listOfData[[1]][indices], function(x) as.numeric(as.character(x)))
  saveRDS(listOfData,"GSEData.rds")
  #combine(paste("GSE",accession,"NormalizedExpressionArray.customCDF.mas5.txt",sep=""), gse)
  message(paste("Your data has been downloaded and is available at: ", getwd(),"/GSEData.rds",sep=""))


  return(listOfData)

}

