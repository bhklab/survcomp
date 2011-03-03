\name{upp7g}
\alias{upp7g}
\docType{data}
\title{Subset of UPP dataset containing gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from the UPP datasets (see section details). The subset contains the seven genes introduced by Desmedt et al. 2008}
\format{
  ExpressionSet with 7 features and 251 samples, containing:
  \itemize{
    \item \code{exprs(upp7g)}: Matrix containing gene expressions as measured by Affymetrix hgu133a and hgu133b technology (single-channel, oligonucleotides).
    \item \code{fData(upp7g)}: AnnotatedDataFrame containing annotations of Affy microarray platform hgu133a and hgu133b.
    \item \code{pData(upp7g)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(upp7g)}: MIAME object containing information about the dataset.
    \item \code{annotation(upp7g)}: Name of the affy chip.
  }
}
\details{
This dataset represents a subset of the study published by Miller et al. 2005. The subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as UPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
}
\source{
\bold{upp:}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3494}
}
\references{
Lance D. Miller, Johanna Smeds, Joshy George, Vinsensius B. Vega, Liza Vergara, Alexander Ploner, Yudi Pawitan, Per Hall, Sigrid Klaar, Edison T. Liu and Jonas Bergh (2005) "An expression signature for p53 status in human breast cancer predicts mutation status, transcriptional effects, and patient survival" \emph{Proceedings of the National Academy of Sciences of the United States of America} \bold{102}(38):13550-13555
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(breastCancerData)
## show the first 5 columns of the expression data
exprs(upp7g)[ ,1:5]
## show the first 6 rows of the phenotype data
head(pData(upp7g))
## show first 20 feature names
featureNames(upp7g)
## show the experiment data summary
experimentData(upp7g)
## show the used platform
annotation(upp7g)
## show the abstract for this dataset
abstract(upp7g)
}
\keyword{datasets}