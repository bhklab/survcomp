\name{breastCancerData}
\alias{breastCancerData}
\docType{data}
\title{Sample data containing six datasets for gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from 6 different datasets (see section details). The subsets contain the seven gene signature introduced by Desmedt et al. 2008.}
\usage{data(breastCancerData)}
\format{
  Six ExpressionSets. Example for 'mainz7g': eSet with 7 features and 200 samples, containing:
  \itemize{
    \item \code{exprs(mainz7g)}: Matrix containing gene expressions as measured by Affymetrix hgu133a technology (single-channel, oligonucleotides).
    \item \code{fData(mainz7g)}: AnnotatedDataFrame containing annotations of Affy microarray platform hgu133a.
    \item \code{pData(mainz7g)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(mainz7g)}: MIAME object containing information about the dataset.
    \item \code{annotation(mainz7g)}: Name of the affy chip.
  }
}
\details{
This dataset represents subsets of the studies published by Schmidt et al. 2008 [mainz7g], Wang. et al. 2005 and Minn et al. 2007 [vdx7g], Miller et al. 2005 [upp7g], Sotiriou et al. 2006 [unt7g], Desmedt et al. 2007 and TRANSBIG [transbig7g], van't Veer et al. 2002 and van de Vijver et al. 2002 [nki7g]. Each subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as uPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
}
\source{
\bold{mainz:}
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE11121}

\bold{transbig:}
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE7390}

\bold{upp:}
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE3494}

\bold{unt:}
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2990}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE6532}

\bold{vdx:}
\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE2034}

\url{http://www.ncbi.nlm.nih.gov/geo/query/acc.cgi?acc=GSE5327}

\bold{nki:}
\url{http://www.rii.com/publications/2002/vantveer.html}
}
\references{
Marcus Schmidt, Daniel Boehm, Christian von Toerne, Eric Steiner, Alexander Puhl, Heryk Pilch, Hans-Anton Lehr, Jan G. Hengstler, Hainz Koelbl and Mathias Gehrmann (2008) "The Humoral Immune System Has a Key Prognostic Impact in Node-Negative Breast Cancer", \emph{Cancer Research}, \bold{68}(13):5405-5413

Christine Desmedt, Fanny Piette, Sherene Loi, Yixin Wang, Francoise Lallemand, Benjamin Haibe-Kains, Giuseppe Viale, Mauro Delorenzi, Yi Zhang, Mahasti Saghatchian d Assignies, Jonas Bergh, Rosette Lidereau, Paul Ellis, Adrian L. Harris, Jan G. M. Klijn, John A. Foekens, Fatima Cardoso, Martine J. Piccart, Marc Buyse and Christos Sotiriou (2007) "Strong Time Dependence of the 76-Gene Prognostic Signature for Node-Negative Breast Cancer Patients in the TRANSBIG Multicenter Independent Validation Series", \emph{Clinical Cancer Research}, \bold{13}(11):3207-3214

Lance D. Miller, Johanna Smeds, Joshy George, Vinsensius B. Vega, Liza Vergara, Alexander Ploner, Yudi Pawitan, Per Hall, Sigrid Klaar, Edison T. Liu and Jonas Bergh (2005) "An expression signature for p53 status in human breast cancer predicts mutation status, transcriptional effects, and patient survival" \emph{Proceedings of the National Academy of Sciences of the United States of America} \bold{102}(38):13550-13555

Christos Sotiriou, Pratyaksha Wirapati, Sherene Loi, Adrian Harris, Steve Fox, Johanna Smeds, Hans Nordgren, Pierre Farmer, Viviane Praz, Benjamin Haibe-Kains, Christine Desmedt, Denis Larsimont, Fatima Cardoso, Hans Peterse, Dimitry Nuyten, Marc Buyse, Marc J. Van de Vijver, Jonas Bergh, Martine Piccart and Mauro Delorenzi (2006) "Gene Expression Profiling in Breast Cancer: Understanding the Molecular Basis of Histologic Grade To Improve Prognosis", \emph{Journal of the National Cancer Institute}, \bold{98}(4):262-272

Y. Wang, J. G. Klijn, Y. Zhang, A. M. Sieuwerts, M. P. Look, F. Yang, D. Talantov, M. Timmermans, M. E. Meijer-van Gelder, J. Yu, T. Jatkoe, E. M. Berns, D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", \emph{Lancet}, \bold{365}:671-679

Andy J. Minn, Gaorav P. Gupta, David Padua, Paula Bos, Don X. Nguyen, Dimitry Nuyten, Bas Kreike, Yi Zhang, Yixin Wang, Hemant Ishwaran, John A. Foekens, Marc van de Vijver and Joan Massague (2007) "Lung metastasis genes couple breast tumor size and metastatic spread", \emph{Proceedings of the National Academy of Sciences}, \bold{104}(16):6740-6745

Laura J. van't Veer, Hongyue Dai, Marc J. van de Vijver, Yudong D. He, Augustinus A.M. Hart, Mao Mao, Hans L. Peterse, Karin van der Kooy, Matthew J. Marton, Anke T. Witteveen, George J. Schreiber, Ron M. Kerkhoven, Chris Roberts, Peter S. Linsley, Rene Bernards and Stephen H. Friend (2002) "Gene expression profiling predicts clinical outcome of breast cancer", \emph{Nature}, \bold{415}:530-536

M. J. van de Vijver, Y. D. He, L. van't Veer, H. Dai, A. M. Hart, D. W. Voskuil, G. J. Schreiber, J. L. Peterse, C. Roberts, M. J. Marton, M. Parrish, D. Atsma, A. Witteveen, A. Glas, L. Delahaye, T. van der Velde, H. Bartelink, S. Rodenhuis, E. T. Rutgers, S. H. Friend and R. Bernards (2002) "A Gene Expression Signature as a Predictor of Survival in Breast Cancer", \emph{New England Journal of Medicine}, \bold{347}(25):1999-2009
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(breastCancerData)
########
## Example for the mainz7g dataset
########
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