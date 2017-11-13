\name{unt7g}
\alias{unt7g}
\docType{data}
\title{Subset of UNT dataset containing gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from the UNT datasets (see section details). The subset contains the seven genes introduced by Desmedt et al. 2008}
\format{
  ExpressionSet with 7 features and 137 samples, containing:
  \itemize{
    \item \code{exprs(unt7g)}: Matrix containing gene expressions as measured by Affymetrix hgu133a and hgu133b technology (single-channel, oligonucleotides).
    \item \code{fData(unt7g)}: AnnotatedDataFrame containing annotations of Affy microarray platform hgu133a and hgu133b.
    \item \code{pData(unt7g)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(unt7g)}: MIAME object containing information about the dataset.
    \item \code{annotation(unt7g)}: Name of the affy chip.
  }
}
\details{
This dataset represents a subset of the study published by Sotiriou et al. 2006. The subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as UPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
}
\source{
\bold{unt:}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2990}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532}
}
\references{
Christos Sotiriou, Pratyaksha Wirapati, Sherene Loi, Adrian Harris, Steve Fox, Johanna Smeds, Hans Nordgren, Pierre Farmer, Viviane Praz, Benjamin Haibe-Kains, Christine Desmedt, Denis Larsimont, Fatima Cardoso, Hans Peterse, Dimitry Nuyten, Marc Buyse, Marc J. Van de Vijver, Jonas Bergh, Martine Piccart and Mauro Delorenzi (2006) "Gene Expression Profiling in Breast Cancer: Understanding the Molecular Basis of Histologic Grade To Improve Prognosis", \emph{Journal of the National Cancer Institute}, \bold{98}(4):262-272
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(breastCancerData)
## show the first 5 columns of the expression data
exprs(unt7g)[ ,1:5]
## show the first 6 rows of the phenotype data
head(pData(unt7g))
## show first 20 feature names
featureNames(unt7g)
## show the experiment data summary
experimentData(unt7g)
## show the used platform
annotation(unt7g)
## show the abstract for this dataset
abstract(unt7g)
}
\keyword{datasets}