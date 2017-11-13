\name{vdx7g}
\alias{vdx7g}
\docType{data}
\title{Subset of VDX dataset containing gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from the VDX datasets (see section details). The subset contains the seven genes introduced by Desmedt et al. 2008}
\format{
  ExpressionSet with 7 features and 344 samples, containing:
  \itemize{
    \item \code{exprs(vdx7g)}: Matrix containing gene expressions as measured by Affymetrix hgu133a technology (single-channel, oligonucleotides).
    \item \code{fData(vdx7g)}: AnnotatedDataFrame containing annotations of Affy microarray platform hgu133a.
    \item \code{pData(vdx7g)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(vdx7g)}: MIAME object containing information about the dataset.
    \item \code{annotation(vdx7g)}: Name of the affy chip.
  }
}
\details{
This dataset represents a subset of the study published by Wang. et al. 2005 and Minn et al. 2007 . The subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as UPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
}
\source{
\bold{vdx:}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034}
}
\references{
Y. Wang, J. G. Klijn, Y. Zhang, A. M. Sieuwerts, M. P. Look, F. Yang, D. Talantov, M. Timmermans, M. E. Meijer-van Gelder, J. Yu, T. Jatkoe, E. M. Berns, D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", \emph{Lancet}, \bold{365}:671-679

Andy J. Minn, Gaorav P. Gupta, David Padua, Paula Bos, Don X. Nguyen, Dimitry Nuyten, Bas Kreike, Yi Zhang, Yixin Wang, Hemant Ishwaran, John A. Foekens, Marc van de Vijver and Joan Massague (2007) "Lung metastasis genes couple breast tumor size and metastatic spread", \emph{Proceedings of the National Academy of Sciences}, \bold{104}(16):6740-6745
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(breastCancerData)
## show the first 5 columns of the expression data
exprs(vdx7g)[ ,1:5]
## show the first 6 rows of the phenotype data
head(pData(vdx7g))
## show first 20 feature names
featureNames(vdx7g)
## show the experiment data summary
experimentData(vdx7g)
## show the used platform
annotation(vdx7g)
## show the abstract for this dataset
abstract(vdx7g)
}
\keyword{datasets}