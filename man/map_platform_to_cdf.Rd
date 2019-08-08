\name{map_platform_to_cdf}
\alias{map_platform_to_cdf}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Retrieve CDF Platform Name
%%  ~~function to do ... ~~
}
\description{
  This is a helper function that maps the GPl number to a CDF platform name on the BrainArray website.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
map_platform_to_cdf(data, platform, accession_id, platform_id)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{data}{GSE object from GEOquery}
  \item{platform}{GPL Platform Name}
  \item{accession_id}{Accession ID corresponding to GSE object}
  \item{platform_id}{GPL Platform Number}
  
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
##map_platform_to_cdf (data,"hgu133a", "11121","GPL570")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~database }% __ONLY ONE__ keyword per line
