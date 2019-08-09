\name{runCdfCode}
\alias{runCdfCode}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  CDF Normalize Data
%%  ~~function to do ... ~~
}
\description{
  This is a helper function that uses the ReadAffy function based on the CDF found and the gene expression data provided (in the format of CEL files)
  .%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
runCdfCode(data, wd,accession,cdfname, platform_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
    \item{data}{GSE object from GEOquery}
     \item{wd}{Working directory filepath}
  \item{accession}{Accession ID corresponding to GSE object}
      \item{cdfname}{Name of CDF File}
    \item{platform_id}{GPL platform number}

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

##runCdfCode (data, "/Users/johndoe/Documents/GEODATA/##GSE11121",
##"11121","hgu133ahsentrezgcdf", "GPL570")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~database }% __ONLY ONE__ keyword per line
