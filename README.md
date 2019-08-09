# GEOpipeline
Pipeline for automatically processing Gene Expression Omnibus data by retrieving it from online, finding the appropriate CDF file, normalizing the data, and grouping together gene expression data and phenotype information into one simple object. 

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
