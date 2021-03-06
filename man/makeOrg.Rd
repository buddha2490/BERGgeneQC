\name{makeOrg}
\alias{makeOrg}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{
Creates a set of output directories for the GenomeStudio workflow
}
\description{
The user creates a home directory where all the files will be located.  This home directory will include folders for the GenomeStudio projects and the "Output Tables" directory where all the GenomeStudio output was saved.  makeOrg() will then check the contents of the home directory and create any and all folders necessary for the subsequent analyses.
}
\usage{
makeOrg(home)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{e}{
home directory filepath object.  This should be the parent folder that includes the GenomeStudio projects and the "Output Tables" directories.
}
}
\details{
%%  ~~ If necessary, more details than the description above ~~
}
\value{
No values are returned from this function.
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
makeOrg(src)

}
