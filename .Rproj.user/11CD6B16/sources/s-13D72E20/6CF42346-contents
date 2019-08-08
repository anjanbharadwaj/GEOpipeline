\name{fetchData}
\alias{fetchData}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
  Retrieve GEO Data
%%  ~~function to do ... ~~
}
\description{
  The function will utilize the other helper functions in the package to download data from the online Gene Expression Omnibus, find the proper CDF file to use from the BrainArray website (alternatively, if it cannot be found online, the user will be prompted to enter the link to a CDF file of their choice), and normalize the data based on the CDF file. If the dataset from online is not in the .CEL format, the tool will not perform normalization.
%%  ~~ A concise (1-5 lines) description of what the function does. ~~
}
\usage{
fetchData(accession, output_dir)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{accession}{The numerical portion of a GSE accession ID, in character format. Example"" for accession GSE11121, the user should provide the argument as "11121"}
  \item{output_dir}{A character describing the output file directory. Example: "/Users/johndoe/Documents/GEODATA"}
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
##fetchData("11121", "/Users/johndoe/Documents/GEODATA")

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~utilities }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~database }% __ONLY ONE__ keyword per line
