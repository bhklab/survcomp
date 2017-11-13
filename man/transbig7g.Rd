\name{transbig7g}
\alias{transbig7g}
\docType{data}
\title{Subset of the TRANSBIG dataset containing gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from the TRANSBIG dataset (see section details). The subset contains the seven genes introduced by Desmedt et al. 2008}
\format{
  ExpressionSet with 7 features and 198 samples, containing:
  \itemize{
    \item \code{exprs(transbig7g)}: Matrix containing gene expressions as measured by Affymetrix hgu133a technology (single-channel, oligonucleotides).
    \item \code{fData(transbig7g)}: AnnotatedDataFrame containing annotations of Affy microarray platform hgu133a.
    \item \code{pData(transbig7g)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(transbig7g)}: MIAME object containing information about the dataset.
    \item \code{annotation(transbig7g)}: Name of the affy chip.
  }
}
\details{
This dataset represents a subset of the study published by Desmedt et al. 2007. The subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as UPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
}
\source{
\bold{transbig:}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7390}
}
\references{
Christine Desmedt, Fanny Piette, Sherene Loi, Yixin Wang, Francoise Lallemand, Benjamin Haibe-Kains, Giuseppe Viale, Mauro Delorenzi, Yi Zhang, Mahasti Saghatchian d Assignies, Jonas Bergh, Rosette Lidereau, Paul Ellis, Adrian L. Harris, Jan G. M. Klijn, John A. Foekens, Fatima Cardoso, Martine J. Piccart, Marc Buyse and Christos Sotiriou (2007) "Strong Time Dependence of the 76-Gene Prognostic Signature for Node-Negative Breast Cancer Patients in the TRANSBIG Multicenter Independent Validation Series", \emph{Clinical Cancer Research}, \bold{13}(11):3207-3214
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(breastCancerData)
## show the first 5 columns of the expression data
exprs(transbig7g)[ ,1:5]
## show the first 6 rows of the phenotype data
head(pData(transbig7g))
## show first 20 feature names
featureNames(transbig7g)
## show the experiment data summary
experimentData(transbig7g)
## show the used platform
annotation(transbig7g)
## show the abstract for this dataset
abstract(transbig7g)
}
\keyword{datasets}