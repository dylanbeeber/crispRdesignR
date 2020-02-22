\name{FindsgRNAfunction}
\alias{sgRNA_design_function}
\title{
sgRNA target design for Shiny App}
\description{
Warning: This function should not be directly called by the user - it must be called though RunShiny.R

Designs sgRNA based on inputs provided in the Shiny App.
}
\usage{
sgRNA_design_function(userseq, genomename, gtf,
designprogress, userPAM, calloffs, annotateoffs)
}
\arguments{
  \item{userseq}{
    The target sequence to generate sgRNA guides for. Can either be a character sequence containing DNA bases or the name of a fasta file in the working directory.
}
  \item{genomename}{
    The name of a geneome (from the BSgenome package) to check off-targets for.
  }
  \item{gtf}{The name of a genome annotation file (.gtf) in the working directory to check off-target sequences against.
  }
  \item{designprogress}{Assists in communicating the progress of the sgRNA design to the Shiny App.
  }
  \item{userPAM}{An optional argument used to set a custom PAM for the sgRNA. If not set, the function will default to the "NGG" PAM. Warning: Doench efficieny scores are only accurate for the "NGG" PAM.
  }
  \item{calloffs}{If TRUE, the function will search for off-targets in the genome chosen specified by the genomename argument. If FALSE, off-target calling will be skipped.
  }
  \item{annotateoffs}{If TRUE, the function will provide annotations for the off-targets called using the genome annotation file specified by the gtfname argument. If FALSE, off-target annotation will be skipped.
  }
}
\author{
Dylan Beeber
}