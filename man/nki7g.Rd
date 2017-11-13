\name{nki7g}
\alias{nki7g}
\docType{data}
\title{Subset of NKI dataset containing gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from the NKI datasets (see section details). The subset contains the seven genes introduced by Desmedt et al. 2008}
\format{
  ExpressionSet with 7 features and 337 samples, containing:
  \itemize{
    \item \code{exprs(nki7g)}: Matrix containing gene expressions as measured by Agilent technology (dual-channel, oligonucleotides).
    \item \code{fData(nki7g)}: AnnotatedDataFrame containing annotations of Agilent microarray platform.
    \item \code{pData(nki7g)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(nki7g)}: MIAME object containing information about the dataset.
    \item \code{annotation(nki7g)}: Name of the agilent chip.
  }
}
\details{
This dataset represents a subset of the study published by van't Veer et al. 2002 and van de Vijver et al. 2002. The subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as UPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
}
\source{
\bold{nki:}

\url{http://www.rii.com/publications/2002/vantveer.html}
}
\references{
Laura J. van't Veer, Hongyue Dai, Marc J. van de Vijver, Yudong D. He, Augustinus A.M. Hart, Mao Mao, Hans L. Peterse, Karin van der Kooy, Matthew J. Marton, Anke T. Witteveen, George J. Schreiber, Ron M. Kerkhoven, Chris Roberts, Peter S. Linsley, Rene Bernards and Stephen H. Friend (2002) "Gene expression profiling predicts clinical outcome of breast cancer", \emph{Nature}, \bold{415}:530-536

M. J. van de Vijver, Y. D. He, L. van't Veer, H. Dai, A. M. Hart, D. W. Voskuil, G. J. Schreiber, J. L. Peterse, C. Roberts, M. J. Marton, M. Parrish, D. Atsma, A. Witteveen, A. Glas, L. Delahaye, T. van der Velde, H. Bartelink, S. Rodenhuis, E. T. Rutgers, S. H. Friend and R. Bernards (2002) "A Gene Expression Signature as a Predictor of Survival in Breast Cancer", \emph{New England Journal of Medicine}, \bold{347}(25):1999-2009
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(breastCancerData)
## show the first 5 columns of the expression data
exprs(nki7g)[ ,1:5]
## show the first 6 rows of the phenotype data
head(pData(nki7g))
## show first 20 feature names
featureNames(nki7g)
## show the experiment data summary
experimentData(nki7g)
## show the used platform
annotation(nki7g)
## show the abstract for this dataset
abstract(nki7g)
}
\keyword{datasets}