\name{sgRNA_design}
\alias{sgRNA_design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sgRNA Target Design}
\description{
sgRNA_design returns information to design sgRNA sequences based on a given target sequence, a genome to annotate off-target information, and a genome annoation file (.gtf), to annotate the off-target findings.
}
\usage{
sgRNA_design(usersequence, genomename, gtfname)
}
%- maybe also 'usage' for other objects documented here.
\arguments{
  \item{usersequence}{
    The target sequence to generate sgRNA guides for. Can either be a character sequence containing DNA bases or the name of a fasta file in the working directory.
}
  \item{genomename}{
    The name of a geneome (from the BSgenome package) to check off-targets for.
  }
  \item{gtfname}{The name of a genome annotation file (.gtf) in the working directory to check off-target sequences against.
  }
}
\details{
Important Note: When designing sgRNA for large genomes (billions of base pairs), use short query DNA sequences (under 500 bp). Depending on your hardware checking for off-targets can be quite computationally intensive and may take several hours if not limited to smaller query sequences.
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
alldata <- sgRNA_design("GGCAGAGCTTCGTATGTCGGCGATTCATCTCAAGTAGAAGATCCTGGTGCAGTAGGCCTATGTGAGTTTTTGAAGGGGGTTCAAAGCGCCTTGTAA", BSgenome.Scerevisiae.UCSC.sacCer2, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz")

alldata <- sgRNA_design("DAK1.fasta", BSgenome.Scerevisiae.UCSC.sacCer2, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz")

% Add one or more standard keywords, see file 'KEYWORDS' in the
% R documentation directory.
\keyword{ ~kwd1 }% use one of  RShowDoc("KEYWORDS")
\keyword{ ~kwd2 }% __ONLY ONE__ keyword per line
