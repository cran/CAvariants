\name{asbestos}
\alias{asbestos}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Selikoff's data, two-way contingency table.
}
\description{
%%  ~~ A concise (1-5 lines) description of the dataset. ~~
The data set consists of 4 rows and 5 columns. 
The rows represent the degree of severity of asbestosis and the 
columns concern the time of exposure to asbestos of workers
}
\usage{data(asbestos)}
\format{
  The format is:
row names [1:4] "None" "Grade1" "Grade2" "Grade3"    \cr
col names [1:5] "0-9" "10-19" "20-29" "30-39" "40+" 
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Selikoff IJ 1981 Household risks with inorganic fibers. Bulletin of the New York Academy of Medicine, 57, 947 -- 961. \cr
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. John Wiley & Sons.
}
\examples{
asbestos <-
structure(c(310, 36, 0, 0, 212, 158, 9, 0, 21, 35, 17, 4, 25, 
102, 49, 18, 7, 35, 51, 28), .Dim = 4:5, .Dimnames = list(c("none", 
"grade1", "grade2", "grade3"), c("0-9", "10-19", "20-29", "30-39", 
"40+")))
dim(asbestos)
dimnames(asbestos) 
}
\keyword{datasets}
