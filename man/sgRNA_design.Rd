\name{sgRNA_design}
\alias{sgRNA_design}
%- Also NEED an '\alias' for EACH other topic documented here.
\title{sgRNA Target Design}
\description{
sgRNA_design returns information to design sgRNA sequences based on a given target sequence, a genome to annotate off-target information, and a genome annoation file (.gtf), to annotate the off-target findings.
}
\usage{
sgRNA_design(userseq, genomename, gtfname, calloffs = TRUE, annotateoffs = TRUE)
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
  \item{userPAM}{An optional argument used to set a custom PAM for the sgRNA. If not set, the function will default to the "NGG" PAM. Warning: Doench efficieny scores are only accurate for the "NGG" PAM.
  }
  \item{calloffs}{If TRUE, the function will search for off-targets in the genome chosen specified by the genomename argument. If FALSE, off-target calling will be skipped.
  }
  \item{annotateoffs}{If TRUE, the function will provide annotations for the off-targets called using the genome annotation file specified by the gtfname argument. If FALSE, off-target annotation will be skipped.
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
alldata <- sgRNA_design("GGCAGAGCTTCGTATGTCGGCGATTCATCTCAAGTAGAAGATCCTGGTGCAGTAGG", BSgenome.Scerevisiae.UCSC.sacCer2, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz")

alldata <- sgRNA_design("DAK1.fasta", BSgenome.Scerevisiae.UCSC.sacCer2, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz")

alldata <- sgRNA_design("DAK1.fasta", BSgenome.Scerevisiae.UCSC.sacCer2, "Saccharomyces_cerevisiae.R64-1-1.92.gtf.gz", calloffs = TRUE, annotateoffs = FALSE)
