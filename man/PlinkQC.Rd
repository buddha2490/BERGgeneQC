\name{PlinkQC}
\alias{PlinkQC}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
R interface for running a Plink QC workflow.
}
\description{
This is an interface between R and PLINK so we don???t have to directly run PLINK code.  It will run a bunch of QC statistics for samples and SNPs and generate a report of exclusions.
}
\usage{
plinkQC(home)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{home}{
Home file path for the projects.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
An output RDS file will be created for each batch and saved to the ~/PLINK directory
}
\references{
%% ~put references to the literature/web site here ~
}
\author{
Brian Carter
}
\note{
%%  ~~further notes~~
}

%% ~Make other sections like Warning with \section{Warning }{....} ~

\seealso{
%% ~~objects to See Also as \code{\link{help}}, ~~~
}
\examples{

src <- file.path("s:/diet validation genotyping")
PlinkQC(src)
}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
