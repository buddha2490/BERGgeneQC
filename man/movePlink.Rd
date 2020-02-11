\name{movePlink}
\alias{movePlink}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
The PLINK files created by GenomeStudio are named and organized badly.  This function will rename everything based on the batch name and prepare them for the QC analysis.
}
\description{
The function will recognize only the folders that are need to be processed and will ignore the folders already completed.
}
\usage{
movePlink(home)
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
movePlink(src)

}
% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
