\name{finalReport}
\alias{finalReport}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Summarizes all of the output files used in the GenomeStudio processing workflow and returns a QC Report.
}
\description{
This function summarizes much of the GenomeStudio output and the various QC statistics.  It will return a summary table for each sample (useful for converting sample IDs), and a SNP summary table for each marker.  Finally, a set of exclusion tables will be returned.

}
\usage{
finalReport(home)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{home}{
Home directory location
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
The function will save a list of three data frames to the Summary Reports subfolder

  SummaryReport - Summary exclusion tables

  SampleSummary - All QC metrics for each sample in the batch

  SNPSummary - All QC metrics for the SNPs in the batch.

The function will also return the SummaryReport exclusion table for each batch

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

datasrc <- file.path("s:/diet validation genotyping")
Report_12Feb2020 <- finalReport(datasrc)

}

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
