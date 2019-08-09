# GEOpipeline
Pipeline for automatically processing Gene Expression Omnibus data by retrieving it from online, finding the appropriate CDF file, normalizing the data, and grouping together gene expression data and phenotype information into one simple object. 

In the GEO pipeline, the user enters their GEO accession number and the location of the output directory where files will be created. First, the pipeline queries the GEO database and downloads the specific GSE dataset from online. Depending on the GPL platform number associated with the dataset, it will search through an online mapping file I curated (https://raw.githubusercontent.com/anjanbharadwaj/GPLtoBiocAnnotations/master/gplToBioc.txt) to find the corresponding annotation package. Once this is done, the pipeline will search through http://brainarray.mbni.med.umich.edu/Brainarray/Database/CustomCDF/23.0.0/entrezg.asp to find the corresponding CDF (chip-definition-file) for the given annotation package. 

The CDF file is used to normalize the dataset in a standard fashion. Since data on GEO can be submitted publicly by any researcher, and they may have used different normalization schemes for their data, we instead download the raw microarray data, which comes in the format of several hundreds-thousands of .CEL files). These .CEL files are processed using the CDF file, and allows us to convert the raw data into a normalized gene expression table. This allows us to compare gene expression data from different microarray datasets, as all of the datasets will be in the same space if we use these CDF files.

Finally, phenotype data for each of the samples/cells is retrieved from the GEO database and is included in the output object/file, which is a list comprised of: (1) the gene expression table and (2) the phenotype table. 

------



## Installation

1. install.packages("devtools")
2. devtools::install_github("anjanbharadwaj/GEOpipeline", dependencies = T)


## Usage

### GEOPipeline::fetchData(accession, location)

#### Parameters

- accession = GEO accession number
  - Example: for https://www.ncbi.nlm.nih.gov/gds/?term=GSE11121, accession = "11121"
- location = filepath of output directory
  - Example: "/Users/johndoe/Documents"
- Example: GEOPipeline::fetchData("11121", "/Users/Name/Documents")

#### Values

- Returns a list with the following data
  - list[[1]] = gene expression table
  - list[[2]] = phenotype information pulled from GEO database
- An RDS file with both of these objects combined will be availabe at [location]/[GSE Accession]-[GPL platform]/GSEData.rds



## Questions

Any questions can be directed to anjanbharadwaj02@gmail.com. Please view R documentation files provided in the package. Vignettes will be added in the near future.
