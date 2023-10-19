\name{shopdataM}
\alias{shopdataM}
\docType{data}
\title{
%%   ~~ data name/kind ... ~~
Two-way contingency table of Dutch shoplifting (1977-1978) 
}
\description{
This two-way contingency table summarises, in part, the results of a survey 
of the Dutch Central Bureau of Statistics (Israels, 1987). The table considers a 
sample of 20819 men who were suspected of shoplifting in stores of the Netherlands between 1977 and 1978.}
\usage{data(shopdataM)}
\format{
The format is:\cr
row names [1:13] "clothing" "accessories" "tobacco" "stationary" ...\cr
col names [1:9] "M12<" "M13" "M16" "M19" ...
}
%\details{
%%  ~~ If necessary, more details than the __description__ above ~~
%}
%\source{
%%  ~~ reference to a publication or URL from which the data were obtained ~~
%}
\references{
Beh EJ and Lombardo R 2014 Correspondence Analysis: Theory, Practice and New Strategies. Wiley. \cr
Israels A 1987 Eigenvalue Techniques for Qualitative Data. DSWO Press, Leiden.\cr
Lombardo R Beh EJ 2016 Variants of Simple Correspondence Analysis. The R Journal, 8 (2), 167--184.
}
\examples{
shopdataM <- structure(c(81, 66, 150, 667, 67, 24, 47, 430, 743, 132, 32, 
197, 209, 138, 204, 340, 1409, 259, 272, 117, 637, 684, 408, 
57, 547, 550, 304, 193, 229, 527, 258, 368, 98, 246, 116, 298, 
61, 402, 454, 384, 149, 151, 84, 146, 141, 61, 40, 13, 71, 52, 
138, 252, 942, 297, 313, 92, 251, 167, 193, 30, 16, 130, 111, 
280, 624, 359, 109, 136, 36, 96, 67, 75, 11, 16, 31, 54, 200, 
195, 178, 53, 121, 36, 48, 29, 50, 5, 6, 14, 41, 152, 88, 137, 
68, 171, 37, 56, 27, 55, 17, 3, 11, 50, 211, 90, 45, 28, 145, 
17, 41, 7, 29, 28, 8, 10, 28, 111, 34), .Dim =  c(13L,9L), .Dimnames = list(
c("clothing", "accessories", "tobacco", "stationary", "books", 
"records", "household", "candy", "toys", "jewelry", "perfumes", 
"hobby", "other"), c("M12<", "M13", "M16", "M19", "M25", 
"M35", "M45", "M57", "M65+")))
dim(shopdataM)
}
\keyword{datasets}
