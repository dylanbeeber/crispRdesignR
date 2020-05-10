\name{crispRdesignRUI}
\alias{crispRdesignRUI}
\title{UI caller for crispRdesignR}
\description{
Activates the shiny UI for the crispRdesignR package
}
\usage{
crispRdesignRUI()
}
\value{
No return value, called to initiate user interface.
}
\author{
Dylan Beeber
}
\examples{
requireNamespace("gbm", quietly = TRUE)
requireNamespace("Biostrings", quietly = TRUE)
requireNamespace("stringr", quietly = TRUE)
if (interactive()) {
  crispRdesignRUI()
  }
}
