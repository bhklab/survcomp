\name{sampleData}
\alias{sampleData}
\docType{data}
\title{Sample data containing six datasets for gene expression, annotations and clinical data.}
\description{
This dataset contains a subset of the gene expression, annotations and clinical data from 6 different datasets (see section details). The subsets contain the seven gene signature introduced by Desmedt et al. 2008/}
\usage{data(sampleData)}
\format{
  Six ExpressionSets. Example for 'mainzSample': eSet with 7 features and 200 samples, containing:
  \itemize{
    \item \code{exprs(mainzSmall)}: Matrix containing gene expressions as measured by Affymetrix hgu133a technology (single-channel, oligonucleotides).
    \item \code{fData(mainzSmall)}: AnnotatedDataFrame containing annotations of Affy microarray platform hgu133a.
    \item \code{pData(mainzSmall)}: AnnotatedDataFrame containing Clinical information of the breast cancer patients whose tumors were hybridized.
    \item \code{experimentalData(mainzSmall)}: MIAME object containing information about the dataset.
    \item \code{annotation(mainzSmall)}: Name of the affy chip.
  }
}
\details{
This dataset represents subsets of the studies published by Schmidt et al. 2008 [mainz], Wang. et al. 2005 and Minn et al. 2007 [vdx], Miller et al. 2005 [upp], Sotiriou et al. 2006 [unt], Desmedt et al. 2007 and TRANSBIG [transbig], van't Veer et al. 2002 and van de Vijver et al. 2002 [nki]. Each subset contains the genes AURKA (also known as STK6, STK7, or STK15), PLAU (also known as uPA), STAT1, VEGF, CASP3, ESR1, and ERBB2, as introduced by Desmedt et al. 2008. The seven genes represent the proliferation, tumor invasion/metastasis, immune response, angiogenesis, apoptosis phenotypes, and the ER and HER2 signaling, respectively.
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
Marcus Schmidt and Daniel Boehm and Christian von Toerne and Eric Steiner and Alexander Puhl and Heryk Pilch and Hans-Anton Lehr and Jan G. Hengstler and Hainz Koelbl and Mathias Gehrmann (2008)"The Humoral Immune System Has a Key Prognostic Impact in Node-Negative Breast Cancer", \emph{Cancer Research}, \bold{68}(13):5405-5413

Christine Desmedt and Fanny Piette and Sherene Loi and Yixin Wang and Francoise Lallemand and Benjamin Haibe-Kains and Giuseppe Viale and Mauro Delorenzi and Yi Zhang and Mahasti Saghatchian d Assignies and  Jonas Bergh and Rosette Lidereau and Paul Ellis and Adrian L. Harris and Jan G. M. Klijn and John A. Foekens and Fatima Cardoso and Martine J. Piccart and Marc Buyse and Christos Sotiriou (2007) "Strong Time Dependence of the 76-Gene Prognostic Signature for Node-Negative Breast Cancer Patients in the TRANSBIG Multicenter Independent Validation Series", \emph{Clinical Cancer Research}, \bold{13}(11):3207-3214

Lance D. Miller and Johanna Smeds and Joshy George and Vinsensius B. Vega and Liza Vergara and Alexander Ploner and Yudi Pawitan and Per Hall and Sigrid Klaar and Edison T. Liu and Jonas Bergh (2005) "An expression signature for p53 status in human breast cancer predicts mutation status, transcriptional effects, and patient survival" \emph{Proceedings of the National Academy of Sciences of the United States of America} \bold{102}(38):13550-13555

Christos Sotiriou and Pratyaksha Wirapati and Sherene Loi and Adrian Harris and Steve Fox and Johanna Smeds and Hans Nordgren and Pierre Farmer and Viviane Praz and Benjamin Haibe-Kains and Christine Desmedt and Denis Larsimont and Fatima Cardoso and Hans Peterse and Dimitry Nuyten and Marc Buyse and Marc J. Van de Vijver and Jonas Bergh and Martine Piccart and Mauro Delorenzi (2006) "Gene Expression Profiling in Breast Cancer: Understanding the Molecular Basis of Histologic Grade To Improve Prognosis", \emph{Journal of the National Cancer Institute}, \bold{98}(4):262-272

Y. Wang and J. G. Klijn and Y. Zhang and A. M. Sieuwerts and M. P. Look and F. Yang and D. Talantov and M. Timmermans and M. E. Meijer-van Gelder and J. Yu and T. Jatkoe and E. M. Berns and D. Atkins and J. A. Foekens (2005) "Gene-Expression Profiles to Predict Distant Metastasis of Lymph-Node-Negative Primary Breast Cancer", \emph{Lancet}, \bold{365}:671-679

Minn, Andy J. and Gupta, Gaorav P. and Padua, David and Bos, Paula and Nguyen, Don X. and Nuyten, Dimitry and Kreike, Bas and Zhang, Yi and Wang, Yixin and Ishwaran, Hemant and Foekens, John A. and van de Vijver, Marc and Massague, Joan (2007) "Lung metastasis genes couple breast tumor size and metastatic spread", \emph{Proceedings of the National Academy of Sciences}, \bold{104}(16):6740-6745

Laura J. van't Veer and Hongyue Dai and Marc J. van de Vijver and Yudong D. He and Augustinus A.M. Hart and Mao Mao and Hans L. Peterse and Karin van der Kooy and Matthew J. Marton and Anke T. Witteveen and George J. Schreiber and Ron M. Kerkhoven and Chris Roberts and Peter S. Linsley an Rene Bernards and Stephen H. Friend (2002) "Gene expression profiling predicts clinical outcome of breast cancer", \emph{Nature}, \bold{415}:530-536

M. J. van de Vijver and Y. D. He and L. van't Veer and H. Dai and A. M. Hart and D. W. Voskuil and G. J. Schreiber and J. L. Peterse and C. Roberts and M. J. Marton and M. Parrish and D. Atsma and A. Witteveen and A. Glas and L. Delahaye and T. van der Velde and H. Bartelink and S. Rodenhuis and E. T. Rutgers and S. H. Friend and R. Bernards (2002) "A Gene Expression Signature as a Predictor of Survival in Breast Cancer", \emph{New England Journal of Medicine}, \bold{347}(25):1999-2009
}
\examples{
## load Biobase package
library(Biobase)
## load the dataset
data(sampleData)
########
## Example for the mainzSample dataset
########
## show the first 5 rows and columns of the expression data
exprs(mainzSample)[1:5,1:5]
## show the first 6 rows of the phenotype data
head(pData(mainzSample))
## show first 20 feature names
featureNames(mainzSample)[1:20]
## show the experiment data summary
experimentData(mainzSample)
## show the used platform
annotation(mainzSample)
## show the abstract for this dataset
abstract(mainzSample)
}
\keyword{datasets}