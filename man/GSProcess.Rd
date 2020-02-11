\name{GSProcess}
\alias{GSProcess}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Process the tables produced by GenomeStudio and saved in the "Output Tables" folder
}
\description{
GenomeStudio will produce three sets of tables:  The Samples table, SNP table, and FullData tables.  Each of these has a piece of data we will need.  This function will process all the output files.
}
\usage{
GSProcess(home)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{home}{
Home file path for the projects.
}
}
\details{
This function can take quite a while to run, probably about half an hour to an hour for a batch of 96 samples.  In order to reduce the total run time, the function will search for previously completed projects in the directories so that they are not rerun.
}
\value{
Three sets of RDS files are saved in the "Processed Files" folder:

Samples_...RDS - sample summary file for each of the samples in the batch

SNPTable_...RDS - SNP summary statistics for the batch.

FullDataTable...RDS - Additional SNP statistics by batch.

This also saves the BAF/LRR statistics for each sample.  Each sample recieves a single file and is saved at "~/BAF/Files"
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
GSProcess(src)
}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
