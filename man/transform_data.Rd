\name{transform_data}
\alias{transform_data}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Finish CDF Normalization and Package Data into 1 List
%%  ~~function to do ... ~~
}
\description{
  This is a helper function that finishes CDF normalization and packages all data (gene expression and phenotypic) into one simple object.
  %%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
transform_data(cdfname, accession, gse, platform_id, filenamesV1)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
	      \item{cdfname}{Name of CDF File}
  \item{accession}{Accession ID corresponding to GSE object}
	    \item{gse}{GSE object from GEOquery}
    \item{platform_id}{GPL platform number}
    	     \item{filenamesV1}{List of CEL file names}


}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
  \item{expression}{The normalized gene expression data}
  \item{phenotype_data}{The corresponding phenotype metadata for each of the samples in the gene expression data}
%%  ~Describe the value returned
%%  If it is a LIST, use
%%  \item{comp1 }{Description of 'comp1'}
%%  \item{comp2 }{Description of 'comp2'}
%% ...
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Anjan Bharadwaj
%%  ~~who you are~~
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{
##transform_data ("hgu133ahsentrezgcdf", "11121",
## gse, "GPL570", ##c("GSM12890.CEL","GSM138290.CEL"))


}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~database }% __ONLY ONE__ keyword per line
