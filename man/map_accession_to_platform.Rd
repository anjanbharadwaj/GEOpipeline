\name{map_accession_to_platform}
\alias{map_accession_to_platform}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Retrieve GPL Platform
%%  ~~function to do ... ~~
}
\description{
  This is an internal function used to retrieve the GPL platform number based on the GSE accession number.%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
map_accession_to_platform(data, accession_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{GSE object from GEOquery}
  \item{accession_id}{Accession ID corresponding to GSE object}
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
##map_accession_to_platform(data, "11121")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~database }% __ONLY ONE__ keyword per line
