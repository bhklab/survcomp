\name{mainz7g}
\alias{mainz7g}
\docType{data}
\title{Subset of MAINZ dataset containing gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from the MAINZ datasets (see section details). The subset contains the seven genes introduced by Desmedt et al. 2008}
\format{
  ExpressionSet with 7 features and 200 samples, containing:
  \itemize{
    \item \code{exprs(mainz7g)}: Matrix containing gene expressions as measured by Affymetrix hgu133a technology (single-channel, oligonucleotides).
    \item \code{fData(mainz7g)}: AnnotatedDataFrame containing annotations of Affy microarray platform hgu133a.
    \item \code{pData(mainz7g)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(mainz7g)}: MIAME object containing information about the dataset.
    \item \code{annotation(mainz7g)}: Name of the affy chip.
  }
}
\details{
This dataset represents a subset of the study published by Schmidt et al. 2008. The subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as UPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
}
\source{
\bold{mainz:}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121}
}
\references{
Marcus Schmidt, Daniel Boehm, Christian von Toerne, Eric Steiner, Alexander Puhl, Heryk Pilch, Hans-Anton Lehr, Jan G. Hengstler, Hainz Koelbl and Mathias Gehrmann (2008) "The Humoral Immune System Has a Key Prognostic Impact in Node-Negative Breast Cancer", \emph{Cancer Research}, \bold{68}(13):5405-5413
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(breastCancerData)
## show the first 5 columns of the expression data
exprs(mainz7g)[ ,1:5]
## show the first 6 rows of the phenotype data
head(pData(mainz7g))
## show first 20 feature names
featureNames(mainz7g)
## show the experiment data summary
experimentData(mainz7g)
## show the used platform
annotation(mainz7g)
## show the abstract for this dataset
abstract(mainz7g)
}
\keyword{datasets}