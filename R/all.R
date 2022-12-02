#'
#' ASRgenomics: A package with complementary genomic functions
#'
#' @description
#'
#' This package presents a series of molecular and genetic routines in the R
#' environment with the aim of assisting in analytical pipelines before and after
#' use of ASReml-R or another library to perform analyses such as Genomic Selection
#' (GS) or Genome-Wide Association Analyses (GWAS).
#'
#' The main tasks considered are:
#' \itemize{
#' \item Preparing and exploring pedigree, phenotypic and genomic data.
#' \item Calculating and evaluating genomic matrices and their inverse.
#' \item Complementing and expanding results from genomic analyses.
#' }
#'
#' The functions implemented consider aspects such as: filtering SNP data for
#' quality control; assessing a kinship matrix \eqn{\boldsymbol{K}} by
#' reporting diagnostics (statistics and plots); performing Principal
#' Component Analyses (PCA) based on
#' kinship or SNP matrices for understanding population structure;
#' calculating the genomic matrix \eqn{\boldsymbol{G}} and its inverse
#' and assessing their quality; matching pedigree- against genomic-based matrices;
#' tuning up a genomic matrix (bend, blend or align); and obtaining the hybrid
#' matrix \eqn{\boldsymbol{H}} as required with single-step GBLUP (ssGBLUP).
#'
#' For the full manual please visit: \url{https://asreml.kb.vsni.co.uk/wp-content/uploads/sites/3/ASRgenomics_Manual.pdf}
#'
#' @docType package
#' @name ASRgenomics
#'
#' @import data.table
#' @import factoextra
#' @import ggplot2
#' @import scattermore
#' @importFrom utils head
#' @importFrom methods new getFunction
#' @importFrom stats prcomp var cov2cor aggregate cor cov na.omit
#' @importFrom Matrix nearPD
#' @importFrom AGHmatrix Gmatrix
#' @importFrom crayon blue
#' @importFrom cowplot plot_grid
#' @importFrom ellipse ellipse
#' @importFrom superheat superheat

NULL

# \if{html}{
#   \out{<div style="cursor: pointer; text-align: center;" onclick="window.location='https://vsni.co.uk/';">}\figure{VSNi_Logo.png}{options: style = "max-width:10\%;"}\out{</div>}
# }

# Add global variables to avoid NOTES.
utils::globalVariables(
  c("value", "Value", "density",
    "PC1", "PC2", "Ind", "SNP", "MAF",
    "het", "Row", "Col", "n.states",
    "..map", "state2", "state1", "code0",
    "code1A", "code1B", "code2"))
#' Obtains the inverse of the genomic relationship matrix G
#'
#' @description
#' Generates the inverse of a genomic relationship matrix \eqn{\boldsymbol{G}} that is provided.
#' This input matrix should be of the full form (\eqn{n \times n})
#' with individual names assigned to
#' \code{rownames} and \code{colnames}. Several checks for the stability of the matrix are
#' presented based on the reciprocal conditional number.
#'
#' In case of an ill-conditioned matrix,
#' options of blending, bending or aligning before inverting are available.
#' These options will be deprecated (discontinued)
#' in future versions of \code{ASRgenomics}
#' as they can be better implemented in the function \code{G.tuneup()}.
#'
#' Based on procedures published by Nazarian and Gezan \emph{et al.} (2016).
#'
#' @param G Input of the symmetric genomic relationship matrix \eqn{\boldsymbol{G}} in
#' full form (\eqn{n \times n}), to obtain its inverse (default = \code{NULL}).
#' @param A Input of the pedigree relationship matrix \eqn{\boldsymbol{A}}
#' to perform blending or aligning in full form.
#' It should be of the same dimension as the \eqn{\boldsymbol{G}} matrix
#'  (\eqn{n \times n}) (default = \code{NULL}) (to be deprecated).
#' @param rcn.thr A threshold for identifying  the \eqn{\boldsymbol{G}} matrix as an ill-conditioned matrix.
#' Based on the reciprocal conditional number (default = \code{1e-12}).
#' @param blend If \code{TRUE} a "blending" with identity matrix \eqn{\boldsymbol{I}} or pedigree relationship matrix
#' \eqn{\boldsymbol{A}} (if provided) is performed (default = \code{FALSE}) (to be deprecated).
#' @param pblend If blending is requested this is the proportion of the identity matrix \eqn{\boldsymbol{I}} or
#' pedigree relationship matrix \eqn{\boldsymbol{A}} to blend for (default = \code{0.02}) (to be deprecated).
#' @param bend If \code{TRUE} a "bending" is performed by making the matrix near positive definite (default = \code{FALSE}) (to be deprecated).
#' @param eig.tol Defines relative positiveness (\emph{i.e.}, non-zero) of eigenvalues compared to the largest one.
#' It determines which threshold of eigenvalues will be treated as zero (default = \code{NULL}) (to be deprecated).
#' @param align If \code{TRUE} the genomic relationship matrix \eqn{\boldsymbol{G}} is aligned to the
#' pedigree relationship matrix \eqn{\boldsymbol{A}} (default = \code{FALSE}) (to be deprecated).
#' @param digits Set up the number of digits in used to round the output matrix (default = \code{8}).
#' @param sparseform If \code{TRUE} it generates an inverse matrix in sparse form to be used directly in \pkg{asreml} with
#' required attributes (default = \code{FALSE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with three of the following elements:
#' \itemize{
#' \item \code{Ginv}: the inverse of \eqn{\boldsymbol{G}} matrix in full form (only if \code{sparseform = FALSE}).
#' \item \code{Ginv.sparse}: the inverse of \eqn{\boldsymbol{G}} matrix in sparse form (only if \code{sparseform = TRUE}).
#' \item \code{status}: the status (\code{ill-conditioned} or \code{well-conditioned}) of
#' the inverse of \eqn{\boldsymbol{G}} matrix.
#' \item \code{rcn}: the reciprocal conditional number of the inverse of \eqn{\boldsymbol{G}} matrix.
#' }
#'
#' @references
#' Nazarian A., Gezan S.A. 2016. GenoMatrix: A software package for pedigree-based
#' and genomic prediction analyses on complex traits. Journal of Heredity 107:372-379.
#'
#' @export
#'
#' @examples
#' # Example: An ill-conditioned matrix.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Get the inverse of G.
#' GINV <- G.inverse(G = G, bend = FALSE, blend = FALSE, align = FALSE)
#' GINV$Ginv[1:5, 1:5]
#' GINV$status
#'

G.inverse <- function(G = NULL, A = NULL, rcn.thr = 1e-12,
                      blend = FALSE, pblend = 0.02,
                      bend = FALSE, eig.tol = NULL, align = FALSE,
                      digits = 8, sparseform = FALSE, message = TRUE){

  # TODO check why eig.tol is not used in the function.

  # Deprecation traps ---------------------------------------------------------------------------

  if (!is.null(A) | blend | bend | align | !is.null(eig.tol)){
    warning("The arguments \'A', \'blend', \'pblend', \'bend', \'align', \'eig.tol' are still active",
            " but will be discontinued in future versions. Please use \'G.tuneup'",
            " to perform the respective actions.")
  }

  # Traps ---------------------------------------------------------------------------------------

  # Check if the class of G is matrix
  if (is.null(G) || !inherits(G, "matrix")) {
    stop("G should be a valid object of class matrix.")
  }
  # Check the rownames/colnames
  if (is.null(rownames(G))){
    stop("Individual names not assigned to rows of matrix G.")
  }
  if (is.null(colnames(G))){
    stop('Individual names not assigned to columns of matrix G.')
  }
  if ((identical(rownames(G), colnames(G))) == FALSE){
    stop("Rownames and colnames of matrix G do not match.")
  }

  # Checks for other input.
  if (pblend < 0 | pblend > 1) {
    stop("Specification of pblend must be between 0 and 1.")
  }

  if (rcn.thr <= 0) {
    stop("Value for rcn.thr must be positive.")
  }

  if(!is.null(eig.tol)){
    if (eig.tol <= 0) {
      stop("Value for eig.tol must be positive.")
    }
  }

  # Get dimensions.
  n <- dim(G)[1]
  p <- dim(G)[2]

  # Reciprocal Condition Number RCN
  rcn <- rcond(G)
  if (message){
    message('Reciprocal conditional number for original matrix is: ', rcn)
  }

  # Check the RCN default value and inverting the matrix
  if (rcn > rcn.thr){

    # Try to get the inverse of G.
    ginverse <- tryCatch(
      expr = {chol2inv(chol(G))},
      error = function(holder){return(NULL)}
    )

    if (is.null(ginverse)){
      stop("Matrix G is not positive definite.")
    }

    rownames(ginverse) <- rownames(G) # Add
    colnames(ginverse) <- colnames(G) # Add
    ginverse <- round(ginverse, digits)
  }
  if (rcn < rcn.thr & message){
    message('Reciprocal conditional number of G is: ',rcn, ', which is lower than the threshold of ',rcn.thr)
    message("G seems to be an ill-conditioned matrix.")
  }
  if (rcn < rcn.thr & !isTRUE(blend) & !isTRUE(bend) & !isTRUE(align)){
    stop('Consider bending, blending or aligning before invertion to make
         this matrix more stable.')
  }

  if (isTRUE(blend) || isTRUE(bend) || isTRUE(align)){
    Gb <- G.tuneup(G=G, A=A, blend=blend, pblend=pblend, bend=bend, align=align,
                   rcn=FALSE,sparseform=FALSE, determinant=FALSE, message=FALSE)$Gb

    # Try to get the inverse of Gb (tuned).
    ginverse <- tryCatch(
      expr = {chol2inv(chol(Gb))},
      error = function(holder){return(NULL)}
    )

    if (is.null(ginverse)){
      stop("Matrix G (tuned) is not positive definite.")
    }

    rownames(ginverse) <- rownames(Gb)
    colnames(ginverse) <- colnames(Gb)
  }

  ginverse <- round(ginverse, digits)

  # Reciprocal Condition Number RCN
  rcninv<-rcond(ginverse)
  if (message){
    message('Reciprocal conditional number for inverted matrix is: ', rcninv)
  }

  # Evaluating Matrix Empirically
  # Note: these are simple comparisons to check stability of matrix
  n <- ncol(ginverse)
  CN.1 <- ginverse[1,2]/sqrt(ginverse[1,1]*ginverse[1,1])
  CN.N <- ginverse[(n-1),n]/sqrt(ginverse[(n-1),(n-1)]*ginverse[n,n])
  max_diag <- abs(max(diag(ginverse)))
  max_off_diag <- max(abs(ginverse-diag(ginverse)))
  if(abs(CN.1)>0.99 | abs(CN.N)>0.99 | max_diag>1000 |  max_off_diag>1000) {
    if (message){
      message("Inverse of matrix G appears to be ill-conditioned.")
    }
    status <- 'ill-conditioned'
  } else {
    if (message){
      message("Inverse of matrix G does not appear to be ill-conditioned.")
    }
    status <- 'well-conditioned'
  }

  # Obtaining Matrix in Sparse Form if requested (ready for ASReml-R v.4)
  if (isTRUE(sparseform)) {
    ginverse.sparse <- full2sparse(ginverse)
    attr(ginverse.sparse, "INVERSE") <- TRUE
    return(list(Ginv.sparse = ginverse.sparse, rcn=rcninv, status=status))
  } else{
    attr(ginverse, "rowNames") <- rownames(ginverse)
    attr(ginverse, "colNames") <- colnames(ginverse)
    attr(ginverse, "INVERSE") <- TRUE
    return(list(Ginv = ginverse, rcn=rcninv, status=status))
  }

}
#' Obtains the genomic matrix from SNP data for additive or dominant relationships
#'
#' Generates the genomic numerator relationship matrix for
#' additive (VanRaden or Yang) or dominant (Su or Vitezica) relationships.
#' Matrix provided \eqn{\boldsymbol{M}} is of the form \eqn{n \times p}, with \eqn{n} individuals and \eqn{p} markers.
#' Individual and
#' marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as 0, 1, 2 (integers or decimal numbers).
#' Missing values, if present, need to be specified.
#'
#' @param M A matrix with SNP data of form \eqn{n \times p}, with \eqn{n} individuals and \eqn{p} markers.
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as 0, 1, 2 (integers or decimal numbers) (default = \code{NULL}).
#' @param method The method considered for calculation of genomic matrix.
#' Options are: \code{"VanRaden"} and \code{"Yang"} for additive matrix,
#' and \code{"Su"} and \code{"Vitezica"} for dominant matrix (default =  \code{"VanRaden"}).
#' @param na.string A character that is interpreted as missing values (default = \code{"NA"}).
#' @param sparseform If \code{TRUE} it generates a matrix in sparse form to be used directly in \pkg{asreml} with
#' required attributes (default = \code{FALSE}).
#' @param digits Set up the number of digits used to round the output matrix (default = 8).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details
#' Note: If data is provided with missing values, it will process calculations
#' of relationships on pairwise non-missing data.
#'
#' It uses function \code{Gmatrix} for calculations
#' from R package  \pkg{AGHmatrix} (Amadeu \emph{et al.} 2019).
#'
#' @return A list with one of these two elements:
#' \itemize{
#' \item \code{G}: the \eqn{\boldsymbol{G}} matrix in full form (only if \code{sparseform = FALSE}).
#' \item \code{G.sparse}: the \eqn{\boldsymbol{G}} matrix in sparse form (only if \code{sparseform = TRUE}).
#' }
#'
#' @references
#' Amadeu, R.R., Cellon, C., Olmstead, J.W., Garcia, A.A.F, Resende, M.F.R. and P.R. Munoz. 2016.
#' AGHmatrix: R package to construct relationship matrices for autotetraploid and diploid species:
#' A blueberry example. The Plant Genome 9(3). doi: 10.3835/plantgenome2016.01.0009
#'
#' @export
#'
#' @examples
#' # Example: Requesting a full matrix by VanRanden.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' \donttest{
#' # Example: Requesting a sparse form by VanRanden.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = TRUE)$G.sparse
#' head(G)
#' head(attr(G, "rowNames"))
#' }
#'

G.matrix <- function(M = NULL, method = "VanRaden", na.string = "NA",
                     sparseform = FALSE, digits = 8, message = TRUE){

  use.dll <- FALSE
  # Check if the class of M is matrix
  if (is.null(M) || !inherits(M, "matrix")) {
    stop("M should be a valid object of class matrix.")
  }
  M <- as.matrix(M)
  n <- dim(M)[1]
  p <- dim(M)[2]
  if (p < n) {warning('There are more individuals than markers in this data!')}

  # Check the rownames/colnames of M
  if(is.null(colnames(M))) {
    stop('Marker names not assigned to columns of matrix M.')
  }
  if(is.null(rownames(M))) {
    stop('Individual names not assigned to rows of matrix M.')
  }

  # Call libraries to obtain G matrix
  # Using AGHmatrix
  if (!use.dll) {
    #options(echo=FALSE,message=FALSE)
    if (message){
      GHAT <- Gmatrix(M, method=method, missingValue=na.string, maf=0)
    } else {
      GHAT <- silent_(Gmatrix(M, method=method, missingValue=na.string, maf=0))
    }

    #options(echo=TRUE)
  # Using our own DLL
  } else {
    # # Preparing things for the DLL
    # for (i in 1:p) { M[,i] <- as.double(M[,i]) }  # Coersing M to be double precision
    # M <- as.matrix(M)
    # if (method == 'VanRaden') {mm = 2}
    # if (method == 'Yang') {mm = 1}
    # if (method == 'Su') {stop('Method not available in DLL.')}
    # if (method == 'Vitizica') {stop('Method not available in DLL.')}
    # # dyn.load("ghat4.dll", type='Fortran')
    # dyn.load("src/ghat4.dll", type='Fortran')  # Not sure this is the best way!
    # GHAT <- base::.Fortran("ghat", molM=M, ghatM=matrix(as.double(0),n,n),
    #                        mm=as.integer(2), n=as.integer(n), m=as.integer(p))$ghatM
    # #dyn.unload("ghat4.dll")
  }

  # round GHAT
  GHAT <- round(GHAT, digits)

  # Making sure GHAT has the rownames and colnames
  #GHAT <- as.matrix(GHAT)
  if (length(rownames(M)) == 0) {
    rownames(GHAT) <- as.character(1:n)
    colnames(GHAT) <- as.character(1:n)
  } else {
    rownames(GHAT) <- rownames(M)
    colnames(GHAT) <- rownames(M)
  }

  if (isTRUE(sparseform)) {
    GHAT.sparse <- full2sparse(GHAT)
    return(list(G.sparse = GHAT.sparse))
    } else{
    return(list(G = GHAT))
    }


}
#' Generates the conditional predictions of random effects (BLUPs)
#'
#' @description
#' Predicts random effects values for individuals with unobserved responses (here called \code{x},
#' a vector of length \eqn{nx}) based on known random effect values for individuals with
#' observed responses (here called \code{y}, a vector of length \eqn{ny}). This is done using the
#' common genomic relationship matrix  \eqn{\boldsymbol{G}} for all
#' individuals (full matrix of dimension \eqn{(nx + ny) \times (nx + ny)}).
#'
#' The prediction of unobserved responses will be performed through the
#' multivariante Normal conditional distribution. These predictions are identical to
#' what would be obtained if the entire set of individuals (\eqn{nx + ny}) were included into a
#' GBLUP animal model fit with individuals in the set \code{x} coded as missing.
#'
#' The user needs to provide the matrix \eqn{\boldsymbol{G}} in full form.
#' Individual names (\eqn{nx + ny}) should be assigned to \code{rownames} and \code{colnames}, and these
#' can be in any order. If the variance-covariance matrix of the set \code{y} is provided,
#' standard errors of random effects in set \code{x} are calculated.
#'
#' @param G Input of the genomic relationship matrix \eqn{\boldsymbol{G}} in full form
#' (of dimension \eqn{(nx + ny) \times (nx + ny)}) (default = \code{NULL}).
#' @param gy Input of random effects (\emph{e.g.} breeding values) for individuals with known values.
#' Individual names should be assigned to \code{rownames} of this vector
#' and be found on the matrix \eqn{\boldsymbol{G}} (default = \code{NULL}).
#' @param vcov.gy The variance-covariance matrix associated with the random effects from the
#' individuals with known values (set \code{y}, of dimension \eqn{ny \times ny})
#' (default = \code{NULL}).
#'
#' @return A data frame with the predicted random effect values for individuals with
#' unobserved responses in the set \code{x}. If the variance-covariance matrix is provided,
#' standard errors are included in an additional column.
#'
#' @export
#'
#' @examples
#' \dontrun{
#' library(asreml) # Load asreml.
#'
#' # Example: Apple data creating 100 missing observations.
#'
#' # Prepare G (nx + ny).
#' G <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = FALSE)$G
#' dim(G)
#'
#' # Prepare subset of data.
#' # Select only 147 from 247 individuals from pheno.apple and geno.apple.
#' Gy <- G[1:147, 1:147]
#' phenoy <- pheno.apple[1:147, ]
#'
#' # Obtain the BLUPs for the 147 individuals using ASReml-R.
#'
#' # Blend Gy.
#' Gyb <- G.tuneup(G = Gy, blend = TRUE, pblend = 0.02)$Gb
#'
#' # Get the Gy inverse
#' Gyinv  <- G.inverse(G = Gyb, sparseform = TRUE)$Ginv.sparse
#'
#' # Fit a GBLUP model
#' phenoy$INDIV <- as.factor(phenoy$INDIV)
#' modelGBLUP <-
#'  asreml(
#'   fixed = JUI_MOT ~ 1,
#'   random = ~vm(INDIV, Gyinv),
#'   workspace = 128e06,
#'   data = phenoy)
#'
#' # Obtain Predictions - BLUP (set y).
#' BLUP <- summary(modelGBLUP,coef = TRUE)$coef.random
#' head(BLUP)
#' gy <- as.matrix(BLUP[, 1])
#' rownames(gy) <- phenoy$INDIV
#'
#' # Ready to make conditional predictions.
#' g.cond <- G.predict(G = G, gy = gy)
#' head(g.cond)
#' }
#'

G.predict <- function(G = NULL, gy = NULL, vcov.gy = NULL){

  # Check if the class of G is matrix
  if (is.null(G) || !inherits(G, "matrix")) {
    stop("G should be a valid object of class matrix.")
  }
  if (is.null(rownames(G))){
    stop("Individual names not assigned to rows of matrix G.")
  }
  if (is.null(colnames(G))){
    stop("Individual names not assigned to columns of matrix G.")
  }
  if ((identical(rownames(G), colnames(G))) == FALSE){
    stop("Rownames and colnames of matrix G do not match.")
  }

  # Check if the class of gy is matrix
  if (is.null(gy) || !inherits(gy, "matrix")) {
    stop("gy should be a valid object of class matrix.")
  }
  if (is.null(rownames(gy))){
    stop("Individual names not assigned to rows of gy.")
  }

  # Check for consistency between G and resp
  resp <- gy
  respind <- rownames(resp)
  Gind <- rownames(G)

  yidx <-  which(Gind %in% respind == TRUE)
  xidx <- which(Gind %in% respind == FALSE)
  if (length(xidx) == 0) {
    stop("Individuals in G are the same as those in the vector gy, nothing to predict.")
  }

  # Creating indexes
  G_yy <- G[yidx, yidx]
  G_xy <- G[xidx, yidx]

  # Obtain inverse of G_yy
  G_yy.inv <- G.inverse(G=G_yy, bend=FALSE, blend=FALSE, sparseform=FALSE,
                        message=FALSE)
  if (G_yy.inv$status == 'ill-conditioned') {
    stop("Inverse of matrix G is ill-conditioned, use G.tuneup and provide G.")
  } else {
    G_yy.inv  <- G_yy.inv$Ginv
  }

  # Obtain predictions
  beta <- G_xy %*% G_yy.inv
  pred <- beta %*% resp
  pred <- data.frame(predicted.value=pred)

  # Obtaining vcov(preds)
  if (!is.null(vcov.gy)) {
    # Check if the class of Kinv is matrix
    if (is.null(vcov.gy) || !inherits(vcov.gy, "matrix")) {
      stop("Input vcov.resp should be a valid object of class matrix.")
    }
    if (is.null(rownames(vcov.gy))){
      stop("Individual names not assigned to rows of matrix vcov.gy.")
    }
    if (is.null(colnames(vcov.gy))){
      stop("Individual names not assigned to columns of matrix vcov.gy.")
    }
    if ((identical(row.names(vcov.gy), colnames(vcov.gy))) == FALSE){
      stop("Rownames and colnames of matrix vcov.gy do not match.")
    }
    vcov.pred <- beta %*% vcov.gy %*% t(beta)
    st.pred <- sqrt(diag(vcov.pred))
    pred <- data.frame(predicted.value=pred, std.error=st.pred)
  }

  return(pred=pred)
}
#' Tune-up the the genomic relationship matrix G
#'
#' @description
#' Generates a new matrix that can be \emph{blended}, \emph{bended} or \emph{aligned}
#' in order to make it stable for future use or matrix inversion. The input matrix should
#' be of the full form (\eqn{n \times n}) with individual names assigned to
#' \code{rownames} and \code{colnames}.
#'
#' This routine provides three options of tune-up:
#' \itemize{
#'  \item{\emph{Blend}}{. The \eqn{\boldsymbol{G}} matrix is blended (or averaged) with another matrix.
#'  \eqn{\boldsymbol{G}^\ast=(1-p) \boldsymbol{G} + p \boldsymbol{A}}, where
#'  \eqn{\boldsymbol{G}^\ast} is the blended matrix, \eqn{\boldsymbol{G}}
#'  is the original matrix,  \eqn{p} is the
#'  proportion of the identity matrix or pedigree-based relationship matrix to consider,
#'  and \eqn{\boldsymbol{A}} is the matrix to blend.
#'  Ideally, the pedigree-based relationship matrix should be used, but if this is not available (or it is
#'  of poor quality), then it is replaced by an identity matrix \eqn{\boldsymbol{I}}.}
#'  \item{\emph{Bend}}{. It consists on adjusting the original \eqn{\boldsymbol{G}} matrix to obtain a
#'  near positive definite matrix, which is done by making all negative or very small eigenvalues
#'  slightly positive.}
#'  \item{\emph{Align}}{. The original \eqn{\boldsymbol{G}} matrix is aligned to the
#'  matching pedigree relationship
#'  matrix \eqn{\boldsymbol{A}} where an \eqn{\alpha} and \eqn{\beta} parameters are obtained.
#'  More information can be found in the manual or in Christensen \emph{et al.} (2012).}
#' }
#'
#' The user should provide the matrices \eqn{\boldsymbol{G}} and
#' \eqn{\boldsymbol{A}} in full form (\eqn{n \times n})
#' and matching individual names should be
#' assigned to the \code{rownames} and \code{colnames} of the matrices.
#'
#' Based on procedures published by Nazarian and Gezan \emph{et al.} (2016).
#'
#' @param G Input of the genomic matrix \eqn{\boldsymbol{G}} to tune-up
#'  in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param A Input of the matching pedigree relationship matrix \eqn{\boldsymbol{A}}
#' in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param blend If \code{TRUE} a \emph{blending} with identity matrix \eqn{\boldsymbol{I}} or pedigree relationship matrix
#' \eqn{\boldsymbol{A}} (if provided) is performed (default = \code{FALSE}).
#' @param pblend If blending is requested this is the proportion of the identity matrix \eqn{\boldsymbol{I}} or
#' pedigree-based relationship matrix \eqn{\boldsymbol{A}} to blend for (default = \code{0.02}).
#' @param bend If \code{TRUE} a \emph{bending} is performed by making the matrix near positive definite (default = \code{FALSE}).
#' @param eig.tol Defines relative positiveness (\emph{i.e.}, non-zero) of eigenvalues compared to the largest one.
#' It determines which threshold of eigenvalues will be treated as zero (default = \code{1e-06}).
#' @param align If \code{TRUE} the genomic matrix \eqn{\boldsymbol{G}} is \emph{aligned} to the
#' matching pedigree relationship matrix \eqn{\boldsymbol{A}} (default = \code{FALSE}).
#' @param rcn If \code{TRUE} the reciprocal conditional number of the original and the bended,
#' blended or aligned matrix will be calculated (default = \code{TRUE}).
#' @param digits Set up the number of digits used to round the output matrix (default = \code{8}).
#' @param sparseform If \code{TRUE} it generates an inverse matrix in sparse form to be used directly in \pkg{asreml} with
#' required attributes (default = \code{FALSE}).
#' @param determinant If \code{TRUE} the determinant will be calculated, otherwise, this is obtained for
#' matrices of a dimension of less than 1,500 \eqn{\times} 1,500 (default = \code{TRUE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with six of the following elements:
#' \itemize{
#' \item \code{Gb}: the inverse of \eqn{\boldsymbol{G}} matrix in full form
#' (only if sparseform = \code{FALSE}).
#' \item \code{Gb.sparse}: if requested, the inverse of \eqn{\boldsymbol{G}} matrix in
#' sparse form (only if sparseform = \code{TRUE}).
#' \item \code{rcn0}: the reciprocal conditional number of the original matrix.
#' Values near zero are associated with an ill-conditioned matrix.
#' \item \code{rcnb}: the reciprocal conditional number of the blended, bended or aligned
#' matrix. Values near zero are associated with an ill-conditioned matrix.
#' \item \code{det0}: if requested, the determinant of the original matrix.
#' \item \code{blend}: if the matrix was \emph{blended}.
#' \item \code{bend}: if the matrix was \emph{bended}.
#' \item \code{align}: if the matrix was \emph{aligned}.
#' }
#'
#' @references
#' Christensen, O.F., Madsen, P., Nielsen, B., Ostersen, T. and Su, G. 2012.
#' Single-step methods for genomic evaluation in pigs. Animal 6:1565-1571.
#' doi:10.1017/S1751731112000742.
#'
#' Nazarian A., Gezan S.A. 2016. GenoMatrix: A software package for pedigree-based
#' and genomic prediction analyses on complex traits. Journal of Heredity 107:372-379.
#'
#' @export
#'
#' @examples
#' # Example: Apple dataset.
#'
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Blend G matrix.
#' G_blended <- G.tuneup(G = G, blend = TRUE, pblend = 0.05)
#' G_blended$Gb[1:5, 1:5]
#'
#' # Bend G matrix.
#' G_bended <- G.tuneup(G = G, bend = TRUE, eig.tol = 1e-03)
#' G_bended$Gb[1:5, 1:5]
#'
#' \donttest{
#' # Example: Loblolly Pine dataset with pedigree - Aligned G matrix.
#'
#' A <- AGHmatrix::Amatrix(ped.pine)
#' dim(A)
#'
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.2, ind.callrate = 0.20,
#'  na.string = "-9")$M.clean
#'
#' # Get G matrix.
#' G <- G.matrix(M = M.clean, method = "VanRaden", na.string = "-9")$G
#' G[1:5, 1:5]
#' dim(G)
#'
#' # Match G and A.
#' Aclean <- match.G2A(A = A, G = G, clean = TRUE, ord = TRUE, mism = TRUE)$Aclean
#'
#' # Align G with A.
#' G_align <- G.tuneup(G = G, A = Aclean, align = TRUE)
#' G_align$Gb[1:5, 1:5]
#' }
#'

G.tuneup <- function(G = NULL, A = NULL, blend = FALSE, pblend = 0.02,
                     bend = FALSE, eig.tol = 1e-06, align = FALSE, rcn = TRUE,
                     digits = 8, sparseform = FALSE, determinant = TRUE, message = TRUE){

  # Check if the class of G is matrix
  if (is.null(G) || !inherits(G, "matrix")) {
    stop('G should be a valid object of class matrix.')
  }
  # Check the rownames/colnames of G
  if (is.null(rownames(G))){
    stop('Individual names not assigned to rows of matrix G.')
  }
  if (is.null(colnames(G))){
    stop('Individual names not assigned to columns of matrix G.')
  }
  if ((identical(rownames(G), colnames(G))) == FALSE){
    stop("Rownames and colnames of matrix G do not match.")
  }
  # Checks on other parameters
  if (pblend < 0 | pblend > 1) {
    stop("Specification of pblend must be between 0 and 1.")
  }
  if (eig.tol <= 0) {
    stop("Value for eig.tol must be positive.")
  }

  # Check option parameters
  if (isTRUE(blend) & isTRUE(bend) & isTRUE(align)){
    stop('More than one option (blend, bend or align) was requested. Choose only one of them')
  }
  if (!isTRUE(blend) & !isTRUE(bend) & !isTRUE(align)){
    stop('No option (blend, bend or aling) were requested.')
  }
  n <- dim(G)[1]
  p <- dim(G)[2]

  # #G <- as.matrix(forceSymmetric(G))    #Force a square matrix to be a symmetric Matrix
  # if (p != n) { stop('Matrix G is not a square matrix.') }

  if (!is.null(A)) {
    # Check if the class of A is matrix
    if (!inherits(A, "matrix")) {
      stop('A should be a valid object of class matrix.')
    }
    if (is.null(rownames(A))){
      stop('Individual names not assigned to rows of matrix A.')
    }
    if (is.null(colnames(A))){
      stop('Individual names not assigned to columns of matrix A.')
    }
    nA <- dim(A)[1]
    pA <- dim(A)[2]
    #A <- as.matrix(forceSymmetric(A))    #Force a square matrix to be a symmetric Matrix
    if (pA != nA) { stop('Matrix A is not a square matrix.') }
    if (p != pA | n != nA) { warning('Matrix A and G are not of the same dimensions.
                                     Use match_G2A() to obtain subset of matched genotypes.') }
    }

  # Reciprocal Condition Number RCN
  if (isTRUE(rcn)){
  rcn0 <- rcond(G)
  if (message) {
    message('Reciprocal conditional number for original matrix is: ', rcn0)
  }
  } else { rcn0 <- NULL}

  # Checking for determinant
  if (n > 1500) {
    if (!determinant) {
      if (message) {
        message('Determinant is not calculated as matrix is too large.')
      }
      detG <- NULL
    } else {
      detG <- det(G)
      if (message) {
        message('Determinant for original matrix is: ', detG)
      }
    }
  } else {
    detG <- det(G)
    if (message) {
      message('Determinant for original matrix is: ', detG)
    }
  }



  # Performing blend with I if requested
  if (isTRUE(blend)){
    if (is.null(A)) {
      Gb <- (1-pblend)*G + pblend*diag(x=1, nrow=n, ncol=n)
      if (message) {
        message('Matrix was BLENDED using an identity matrix.')
      }
    }
    # Performing blend with A if requested
    if (!is.null(A)) {
      if (p != pA | n != nA) { stop('Matrix A and G are not of the same dimension.') }
      if ( sum(rownames(G) == rownames(A)) < n) {
        warning('Names of rows/columns do not match between G and A matrix.')
        stop('You should use the function match.G2A to match these two matrices.')
      }
      if ( sum(rownames(G) == rownames(A)) == n) {
        Gb <- (1-pblend)*G + pblend*A
        if (message) {
          message('Matrix was BLENDED using the provided A matrix.')
        }
      }
    }
  }

  # Obtaining the near positive definite matrix (bend)
  if (isTRUE(bend)) {
    Gb <- as.matrix(Matrix::nearPD(G, corr=FALSE, keepDiag=FALSE, do2eigen=TRUE,
                           doSym=TRUE, doDykstra=TRUE, only.values=FALSE,
                           eig.tol=eig.tol, conv.tol=1e-07, posd.tol=1e-02,
                           maxit=100, conv.norm.type="I", trace=FALSE)$mat)
    if (message) {
      message('Matrix was BENDED.')
    }
  }

  # Performing alignment with A matrix (align)
  if (isTRUE(align)) {
    if (is.null(A)) {
    stop('A matrix needs to be provided for align be performed.')}
    Xm <- cbind(c(1,1),c(mean(diag(G)),mean(G)))
    y <- as.matrix(c(mean(diag(A)),mean(A)))
    beta <- solve(Xm) %*% y
    Gb <- beta[1] + beta[2]*G
    if (message) {
      message('Matrix was ALIGNED.')
    }
  }

  # Rounding Gb
  Gb <- round(Gb, digits)
  ### Need to make sure rownames(Gb) and colnames(Gb) are correct

  # Reciprocal Condition Number RCN
  if (isTRUE(rcn)){
  rcnb <- rcond(Gb)
  if (message) {
    message('Reciprocal conditional number for tune-up matrix is: ', rcnb)
  }
  } else { rcnb <- NULL}

  # Obtaining Matrix in Sparse Form if requested (ready for ASReml-R v.4)
  if (isTRUE(sparseform)) {
    Gb.sparse <- full2sparse(Gb)
    return(list(Gb.sparse=Gb.sparse, rcn0=rcn, det0=detG, rcnb=rcnb, blend=blend, bend=bend, align=align))
  } else {
    return(list(Gb=Gb, rcn0=rcn0, det0=detG, rcnb=rcnb, blend=blend, bend=bend, align=align))
  }


}
#' Generates the inverse of the hybrid H matrix
#'
#' The single-step GBLUP approach combines the information from the pedigree relationship matrix
#' \eqn{\boldsymbol{A}} and the genomic relationship matrix \eqn{\boldsymbol{G}} in one
#' hybrid relationship matrix called \eqn{\boldsymbol{H}}.
#' This function will calculate directly the inverse of this matrix \eqn{\boldsymbol{H}}.
#' The user should provide the matrices \eqn{\boldsymbol{A}} or
#' its inverse (only one of these is required) and the
#' inverse of the matrix \eqn{\boldsymbol{G}} (\eqn{\boldsymbol{G_{inv}}}) in its full form. Individual names should
#' be assigned to  \code{rownames} and \code{colnames}, and individuals from
#' \eqn{\boldsymbol{G_{inv}}} are verified to be all a subset within individuals from
#' \eqn{\boldsymbol{A}} (or \eqn{\boldsymbol{A_{inv}}}).
#'
#' @param A Input of the pedigree relationship matrix \eqn{\boldsymbol{A}}
#' in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param Ainv Input of the inverse of the pedigree relationship matrix
#' \eqn{\boldsymbol{A}^{-1}} in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param Ginv Input of the inverse of the genomic relationship matrix
#' \eqn{\boldsymbol{G}^{-1}} in full form (\eqn{ng \times ng}) (default = \code{NULL}).
#' @param lambda The scaling factor for \eqn{(\boldsymbol{G}^{-1}-\boldsymbol{A}^{-1}_{22})} (default = \code{NULL}).
#' @param tau The scaling factor for \eqn{\boldsymbol{G}^{-1}} (default = \code{1}).
#' @param omega The scaling factor for \eqn{\boldsymbol{A}^{-1}_{22}} (default = \code{1}).
#' @param sparseform If \code{TRUE} it generates the requested matrix in sparse form to be used
#' directly in \pkg{asreml} with required attributes (default = \code{FALSE}).
#' @param keep.order If \code{TRUE} the original order of the individuals from the
#' \eqn{\boldsymbol{A}} or \eqn{\boldsymbol{A_{inv}}} matrix is kept.
#' Otherwise the non-genotyped individuals are placed first and
#' then genotyped individuals (default = \code{TRUE}).
#' @param digits Set up the number of digits used to round the output matrix (default = \code{8}).
#' @param inverse If \code{TRUE} it generates the inverse of \eqn{\boldsymbol{H}} matrix (default = \code{TRUE})
#' (to be deprecated).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return
#' The inverse of the hybrid matrix \eqn{\boldsymbol{H}} matrix, in full or sparse form with
#' required attributes to be used in \pkg{asreml}.
#'
#' @details
#' The generation of the \eqn{\boldsymbol{H^{-1}}} matrix contains a few scaling factors
#' to help with the calculation of this inverse and
#' to allow further exploration of the combination of the
#' information from the \eqn{\boldsymbol{A^{-1}}} and \eqn{\boldsymbol{G^{-1}}}.
#' We follow the specifications described by Martini \emph{et. al} (2018),
#' which is done by specifying the parameters \eqn{\lambda}, or the pair
#' \eqn{\tau} and \eqn{\omega}.
#'
#' \if{html}{
#' The general expression used is:
#'   \deqn{\boldsymbol{H^{-1}}=\boldsymbol{A^{-1}}+\begin{bmatrix}\boldsymbol{0}&\boldsymbol{0}\\\boldsymbol{0}&(\tau\boldsymbol{G^{-1}}-\omega\boldsymbol{{A_{22}^{-1}}})\end{bmatrix}}
#' }
#'
#' \if{html}{
#' and a more common representation of the above expression is found when \eqn{\tau = \omega = \lambda}, as shown below:
#'   \deqn{\boldsymbol{H^{-1}}=\boldsymbol{A^{-1}}+\begin{bmatrix}\boldsymbol{0}&\boldsymbol{0}\\\boldsymbol{0}&\lambda(\boldsymbol{G^{-1}}-\boldsymbol{{A_{22}^{-1}}})\end{bmatrix}}
#' }
#'
#' If \code{inverse = FALSE} the \eqn{\boldsymbol{H}}
#' matrix is provided instead of its inverse. This option will be deprecated and
#' it is better to use the function \link{H.matrix}.
#'
#' \if{html}{
#' The \eqn{\boldsymbol{H}} matrix is obtained with the following equations:
#'   \deqn{\boldsymbol{H}=\boldsymbol{A}+\begin{bmatrix}\boldsymbol{A}_{12}\boldsymbol{A}_{22}^{-1}(\boldsymbol{G}-\boldsymbol{A}_{22})\boldsymbol{A}_{22}^{-1}\boldsymbol{A}_{21}&\boldsymbol{A}_{12}\boldsymbol{A}_{22}^{-1}(\boldsymbol{G}-\boldsymbol{A}_{22})\\(\boldsymbol{G}-\boldsymbol{A}_{22})\boldsymbol{A}_{22}^{-1}\boldsymbol{A}_{21}&(\boldsymbol{G}-\boldsymbol{A}_{22})\end{bmatrix}}
#' }
#'
#' @references
#' Christensen, O.F., Lund, M.S. 2010. Genomic prediction matrix when some animals
#' are not genotyped. Gen. Sel. Evol. 42(2):1–8.
#'
#' Christensen, O., Madsen, P., Nielsen, B., Ostersen, T., and Su, G. 2012. Single-step methods
#' for genomic evaluation in pigs. Animal 6(10):1565–1571.
#'
#' Legarra, A., Aguilar, I., and Misztal, I. 2009. A relationship matrix including full
#' pedigree and genomic information. J. Dairy Sci. 92:4656-4663.
#'
#' Martini, J.W.R., Schrauf, M.F., Garcia-Baccino, C.A., Pimentel, E.C.G., Munilla, S.,
#' Rogberg-Muñoz, A., Cantet, R.J.C., Reimer, C., Gao, N., Wimmer, V., and Simianer, H. 2018.
#' The effect of the \eqn{H^{-1}} scaling factors \eqn{\tau} and \eqn{\omega}
#' on the structure of \eqn{H} in the single-step procedure.
#' Genet. Sel. Evol. 50:1-9.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Get A matrix.
#' A <- AGHmatrix::Amatrix(data = ped.pine)
#' A[1:5,1:5]
#' dim(A)
#'
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.2, ind.callrate = 0.20,
#'  na.string = "-9",
#'  plots = FALSE)$M.clean
#'
#' # Get G matrix.
#' G <- G.matrix(M.clean, method = "VanRaden", na.string = "-9")$G
#' G[1:5, 1:5]
#' dim(G)
#'
#' # Match G and A.
#' check <- match.G2A(
#'  A = A, G = G,
#'  clean = TRUE, ord = TRUE, mism = TRUE, RMdiff = TRUE)
#'
#' # Align G matrix with A.
#' G_align <- G.tuneup(G = check$Gclean, A = check$Aclean, align = TRUE, sparseform = FALSE)$Gb
#'
#' # Get Ginverse using the G aligned.
#' Ginv <- G.inverse(G = G_align, sparseform = FALSE)$Ginv
#' Ginv[1:5, 1:5]
#' dim(Ginv)
#'
#' # Obtain Hinv.
#' Hinv <- H.inverse(A = A, Ginv = Ginv, lambda = 0.90, sparseform = TRUE)
#' head(Hinv)
#' attr(Hinv, "INVERSE")
#' }
#'

H.inverse <- function(A = NULL, Ainv = NULL, Ginv = NULL,
                      lambda = NULL,tau = 1, omega = 1,
                      sparseform = FALSE, keep.order = TRUE,
                      digits = 8, inverse = TRUE,
                      message = TRUE){

  # Check if we have lambda and make some checks
  if(!is.null(lambda)){
    tau <- lambda
    omega <- lambda
    if(message){
      message("A lambda value was provided and it will be used instead of tau and omega.")
    }
    if (lambda < 0 | lambda > 1) {
      stop("Value of lambda must be between 0 and 1.")
    }
  }else{
    if(message){
      message("No lambda value was provided, tau and omega scaling factors will be considered.")
    }
  }


  # Checks on some input parameters
  if (tau < 0 | tau > 2) {
    stop("Value of tau must be between 0 and 1.")
  }
  if (omega < -1 | omega > 1) {
    stop("Value of omega must be between -1 and 1.")
  }
  if (tau == 0 & omega == 0) {
    stop("Values of tau and omega can not be both equal to zero.")
  }


  # We need A or Ainv
  if (is.null(A) & is.null(Ainv)) {
    stop('A or Ainv needs to be provided.')
  }

  if (!is.null(Ainv)) { # We have Ainv and get A
    if (!inherits(Ainv, "matrix")) {
      stop("Ainv should be a valid object of class matrix.")
    }
    if (is.null(rownames(Ainv))){
      stop("Individual names not assigned to rows of matrix Ainv.")
    }
    if (is.null(colnames(Ainv))){
      stop("Individual names not assigned to columns of matrix Ainv.")
    }
    if ((identical(rownames(Ainv), colnames(Ainv))) == FALSE){
      stop("Rownames and colnames of matrix Ainv do not match.")
    }
    crit<-sum(diag(Ainv)>=2)
    if (crit == 0) {
      warning('It seems like an A matrix was provided not an Ainv matrix.')
    }
    A <- G.inverse(G=Ainv, blend=FALSE, bend=FALSE, message = FALSE)$Ginv
    attributes(A)$INVERSE <- NULL
  } else { # Then we have A and we get Ainv to use
    if (!inherits(A, "matrix")) {
      stop("A should be a valid object of class matrix.")
    }
    if (is.null(rownames(A))){
      stop("Individual names not assigned to rows of matrix A.")
    }
    if (is.null(colnames(A))){
      stop("Individual names not assigned to columns of matrix A.")
    }
    if ((identical(rownames(A), colnames(A))) == FALSE){
      stop("Rownames and colnames of matrix A do not match.")
    }
    crit<-sum(diag(A)>2)
    if (crit >0) {
      warning('It seems like an Ainv matrix was provided not an A matrix.')
    }
    Ainv <- G.inverse(G=A, blend=FALSE, bend=FALSE, message = FALSE)$Ginv

  }


  # We need Ginv
  # Check components of Ginv
  if (is.null(Ginv) || !inherits(Ginv, "matrix")) {
    stop("Ginv should be a valid object of class matrix.")
  }
  if (is.null(rownames(Ginv))){
    stop("Individual names not assigned to rows of matrix Ginv.")
  }
  if (is.null(colnames(Ginv))){
    stop("Individual names not assigned to columns of matrix Ginv.")
  }
  if ((identical(rownames(Ginv), colnames(Ginv))) == FALSE){
    stop("Rownames and colnames of matrix Ginv do not match.")
  }

  # Check for consistency between Ainv and Ginv
  Aind <- rownames(Ainv)
  Gind <- rownames(Ginv)

  notGenotyped <- which((Aind %in% Gind) == FALSE) # which ind has Ainv but not Ginv
  notPedigree <- which(Gind %in% Aind == FALSE)    # which ind has Ginv but no Ainv

  # If Ginv has different ind than Ainv, stop here.
  if (length(notPedigree) > 0) {
    stop("Matrix Ginv has ", length(notGenotyped), " individuals that are not present in A or Ainv.")
  }
  # Check if Ainv has different ind than Ginv
  if (length(notGenotyped) > 0 & message){
    message("Matrix A or Ainv has ", length(notGenotyped), " individuals that are not present in Ginv.")
  }
  if (length(notGenotyped) < 0 & message){
    message("Matrix A or Ainv has all individuals that are present in Ginv.")
  }


  # Creating indexes (sorted matrices)
  idA <- rownames(Ainv)
  idG <- rownames(Ginv)
  idH <- unique(c(idG,idA))
  idH <- rev(idH)
  A <- A[idH,idH]
  Ainv <- Ainv[idH,idH]
  genotyped <- idH %in% idG == TRUE
  A11inv <- Ainv[!genotyped, !genotyped]
  A12inv <- Ainv[!genotyped, genotyped]
  A22inv <- solve(A[genotyped,genotyped])
  Ginv <- Ginv[idH[genotyped],idH[genotyped]]

  # Making sure both Ginv and A22inv agree on sorted individuals
  if (all(rownames(Ginv) == row.names(A22inv)) == FALSE) {  # Then we have a problem
    stop('Order of matrix Ginv does not match subset of A or Ainv.')
  }

  # Traditional method using Ginv & Ainv to obtain Hinv directly.
  if (isTRUE(inverse)){
    # Obtaining Hinv
    H11inv <- matrix(0, nrow=dim(A11inv)[1], ncol=dim(A11inv)[2] )
    H12inv <- matrix(0, nrow=dim(A12inv)[1], ncol=dim(A12inv)[2] )
    H22inv <- tau*Ginv - omega*A22inv
    Hinv <- Ainv + cbind(rbind(H11inv, t(H12inv)), rbind(H12inv, H22inv))

    # round Hinv
    Hinv <- round(Hinv, digits)

    # Same order as Ainverse (if requested)
    if (keep.order) {
      Hinv <- Hinv[match(Aind, rownames(Hinv)), match(Aind, rownames(Hinv))]
    }
    attr(Hinv, "INVERSE") <- TRUE

    # Generating sparseform matrix
    if (sparseform) {
      Hinv <- full2sparse(Hinv)
      attr(Hinv, "INVERSE") <- TRUE
    }
    return(Hinv=Hinv)
  } else {
    A11 <- A[!genotyped, !genotyped]
    A12 <- A[!genotyped, genotyped]
    A21 <- t(A12)
    A22 <- A[genotyped,genotyped]

    #H22 = solve((tau * Ginv + (1 - omega) * A22inv))
    H22 <- G.inverse((tau * Ginv + (1 - omega) * A22inv), message = FALSE)
    if (H22$status == 'ill-conditioned') {
      stop("Matrix H is ill-conditioned, try different scaling factors.")
    } else {
      H22  <- H22$Ginv
    }

    # Old repetitive code.
    # H11 = A12 %*% A22inv %*% (H22 - A22) %*% A22inv %*% A21
    # H12 = A12 %*% A22inv %*% (H22 - A22)
    # H21 = (H22 - A22) %*% A22inv %*% A21
    # H22 = (H22 - A22)

    # Pre-calculated matrices.
    A22A21 <- A22inv %*% A21

    # Finalize multiplications.
    H22 <- (H22 - A22)
    H12 <- A12 %*% A22inv %*% H22
    H11 <- H12 %*% A22A21
    H21 <- H22 %*% A22A21

    # Bind matrices.
    H <- A + cbind(rbind(H11, H21), rbind(H12, H22))

    # round H
    H <- round(H, digits)

    # Same order as A (if requested)
    if (keep.order) {
      H <- H[match(Aind, rownames(H)), match(Aind, rownames(H))]
    }

    # Generating sparseform matrix
    if (sparseform) {
      H <- full2sparse(H)
      attr(H, "INVERSE") <- FALSE
    }
    return(H=H)
  }

}
#' Generates the hybrid \eqn{H} matrix
#'
#' The single-step GBLUP approach combines the information from the pedigree relationship matrix
#' \eqn{\boldsymbol{A}} and the genomic relationship matrix \eqn{\boldsymbol{G}} in one
#' hybrid relationship matrix called \eqn{\boldsymbol{H}}.
#' This function will calculate directly this matrix \eqn{\boldsymbol{H}}.
#' The user should provide the matrices \eqn{\boldsymbol{A}} or
#' its inverse (only one of these is required) and the
#' inverse of the matrix \eqn{\boldsymbol{G}} (\eqn{\boldsymbol{G_{inv}}}) in its full form.
#' Individual names should
#' be assigned to  \code{rownames} and \code{colnames}, and individuals from
#' \eqn{\boldsymbol{G_{inv}}} are verified to be all a subset within individuals from
#' \eqn{\boldsymbol{A}} (or \eqn{\boldsymbol{A_{inv}}}).
#' This function is a wrapper of  the \link{H.inverse} function.
#'
#' @param A Input of the pedigree relationship matrix \eqn{\boldsymbol{A}}
#' in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param Ainv Input of the inverse of the pedigree relationship matrix
#' \eqn{\boldsymbol{A}^{-1}} in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param Ginv Input of the inverse of the genomic relationship matrix
#' \eqn{\boldsymbol{G}^{-1}} in full form (\eqn{ng \times ng}) (default = \code{NULL}).
#' @param lambda The scaling factor for \eqn{(\boldsymbol{G}^{-1}-\boldsymbol{A}^{-1}_{22})} (default = \code{NULL}).
#' @param tau The scaling factor for \eqn{\boldsymbol{G}^{-1}} (default = \code{1}).
#' @param omega The scaling factor for \eqn{\boldsymbol{A}^{-1}_{22}} (default = \code{1}).
#' @param sparseform If \code{TRUE} it generates the requested matrix in sparse form to be used
#' directly in \pkg{asreml} with required attributes (default = \code{FALSE}).
#' @param keep.order If \code{TRUE} the original order of the individuals from the
#' \eqn{\boldsymbol{A}} or \eqn{\boldsymbol{A_{inv}}} matrix is kept.
#' Otherwise the non-genotyped individuals are placed first and
#' then genotyped individuals (default = \code{TRUE}).
#' @param digits Set up the number of digits used to round the output matrix (default = \code{8}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return
#' The hybrid matrix \eqn{\boldsymbol{H}} matrix, in full or sparse form.
#'
#' @md
#' @details
#' This function is currently equivalent to using \link{H.inverse} with (\code{inverse = FALSE}).
#'
#' \if{html}{
#' The \eqn{\boldsymbol{H}} matrix is obtained with the following equations:
#'   \deqn{\boldsymbol{H}=\boldsymbol{A}+\begin{bmatrix}\boldsymbol{A}_{12}\boldsymbol{A}_{22}^{-1}(\boldsymbol{G}-\boldsymbol{A}_{22})\boldsymbol{A}_{22}^{-1}\boldsymbol{A}_{21}&\boldsymbol{A}_{12}\boldsymbol{A}_{22}^{-1}(\boldsymbol{G}-\boldsymbol{A}_{22})\\(\boldsymbol{G}-\boldsymbol{A}_{22})\boldsymbol{A}_{22}^{-1}\boldsymbol{A}_{21}&(\boldsymbol{G}-\boldsymbol{A}_{22})\end{bmatrix}}
#' }
#'
#'
#' @references
#' Christensen, O.F., Lund, M.S. 2010. Genomic prediction matrix when some animals
#' are not genotyped. Gen. Sel. Evol. 42(2):1–8.
#'
#' Christensen, O., Madsen, P., Nielsen, B., Ostersen, T., and Su, G. 2012. Single-step methods
#' for genomic evaluation in pigs. Animal 6(10):1565–1571.
#'
#' Legarra, A., Aguilar, I., and Misztal, I. 2009. A relationship matrix including full
#' pedigree and genomic information. J. Dairy Sci. 92:4656-4663.
#'
#' Martini, J.W.R., Schrauf, M.F., Garcia-Baccino, C.A., Pimentel, E.C.G., Munilla, S.,
#' Rogberg-Muñoz, A., Cantet, R.J.C., Reimer, C., Gao, N., Wimmer, V., and Simianer, H. 2018.
#' The effect of the \eqn{H^{-1}} scaling factors \eqn{\tau} and \eqn{\omega}
#' on the structure of \eqn{H} in the single-step procedure.
#' Genet. Sel. Evol. 50:1-9.
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Get A matrix.
#' A <- AGHmatrix::Amatrix(data = ped.pine)
#' A[1:5,1:5]
#' dim(A)
#'
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.2, ind.callrate = 0.20,
#'  na.string = "-9",
#'  plots = FALSE)$M.clean
#'
#' # Get G matrix.
#' G <- G.matrix(M = M.clean, method = "VanRaden", na.string = "-9")$G
#' G[1:5, 1:5]
#' dim(G)
#'
#' # Match G2A.
#' check <- match.G2A(
#'  A = A, G = G,
#'  clean = TRUE, ord = TRUE, mism = TRUE, RMdiff = TRUE)
#'
#' # Align G matrix with A.
#' G_align <- G.tuneup(G = check$Gclean, A = check$Aclean, align = TRUE, sparseform = FALSE)$Gb
#'
#' # Get Ginverse using the aligned G.
#' Ginv <- G.inverse(G = G_align, sparseform = FALSE)$Ginv
#' Ginv[1:5, 1:5]
#' dim(Ginv)
#'
#' # Obtaining H.
#' H <- H.matrix(A = A, G = Ginv, lambda = 0.90, sparseform = FALSE)
#' H[1:5, 1:5]
#' }
#'

H.matrix <- function(A = NULL, Ainv = NULL, Ginv = NULL,
                     lambda = NULL,tau=1, omega = 1,
                     sparseform = FALSE, keep.order = TRUE,
                     digits = 8, message = TRUE){

  return(
    H.inverse(A = A, Ainv = Ainv, Ginv = Ginv,
              lambda = lambda,tau = tau, omega = omega,
              sparseform = sparseform, keep.order = keep.order,
              digits = digits, message = message,
              inverse = FALSE # This is the only one locked.
    )
  )
}


#' Enhanced heatmap plot for a kinship matrix K
#'
#' Generates a heatmap with dendrogram based on a provided kinship matrix.
#' This matrix can be a pedigree relationship matrix \eqn{\boldsymbol{A}}, a
#' genomic relationship matrix \eqn{\boldsymbol{G}} or a hybrid relationship
#' matrix \eqn{\boldsymbol{H}}.
#' Individual names should be assigned to \code{rownames} and \code{colnames}.
#' It sorts individuals according to dendrogram in both columns and rows.
#'
#' Uses the library \code{superheat} from Barter and Yu (2018) to generate plots.
#'
#' @param K Input of a kinship matrix in full format (\eqn{n \times n}) (default = \code{NULL}).
#' @param dendrogram If \code{TRUE} a dendrogram is added to the columns based on the
#' kinship matrix (default = \code{TRUE}).
#' @param clustering.method The clustering method considered for the dendrogram.
#' Options are: \code{"hierarchical"} and \code{"kmeans"} (default = \code{"hierarchical"}).
#' @param dist.method The method considered to calculate the distance matrix between
#' individuals used for hierarchical clustering. Options are: \code{"euclidean"},
#' \code{"maximum"}, \code{"manhattan"}, \code{"canberra"}, \code{"binary"} and
#' \code{"minkowski"} (default = \code{"euclidean"}).
#' @param row.label If \code{TRUE} the individual names (\code{rownames}) are added as labels to
#' the left of the heatmap (default = \code{TRUE}).
#' @param col.label If \code{TRUE} the individual names (\code{colnames}) are added as labels to
#' the bottom of the heatmap (default = \code{FALSE}).
#'
#' @return
#' A plot with the properties specified by the above arguments.
#'
#' @references
#' Barter, R.L. and Yu, B. 2018. Superheat: An R package for creating beautiful
#' and extendable heatmaps for visualizing complex data.
#' J. Comput. Graph. Stat. 27(4):910-922.
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Plot a subset of the individuals.
#' kinship.heatmap(K = G[1:10, 1:10], dendrogram = TRUE, row.label = TRUE, col.label = TRUE)
#'


kinship.heatmap <- function(K = NULL, dendrogram = TRUE,
                            clustering.method = c("hierarchical", "kmeans"),
                            dist.method = c("euclidean", "maximum", "manhattan",
                                            "canberra", "binary", "minkowski"),
                            row.label = TRUE, col.label = FALSE){

  # Check if the class of K is matrix
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix.")
  }
  # Check the rownames/colnames
  if (is.null(rownames(K))){
    stop("Individual names not assigned to rows of matrix K.")
  }
  if (is.null(colnames(K))){
    stop('Individual names not assigned to columns of matrix K.')
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }
  # Check if the are missing values
  if (any(is.na(K))){
    stop("Matrix K contains some missing data.")
  }
  clustering.method <- match.arg(clustering.method)
  dist.method <- match.arg(dist.method)

  if (isTRUE(col.label)) {
    labCol <- "variable"
  } else {
    labCol <- "none"
  }
  if (isTRUE(row.label)) {
    labRow <- "variable"
  } else {
    labRow <- "none"
  }

  pp <- superheat(K,
            pretty.order.rows = TRUE,
            pretty.order.cols = TRUE,
            col.dendrogram = dendrogram,
            dist.method = dist.method,
            clustering.method = clustering.method,
            scale = FALSE,
            left.label = labRow,
            left.label.col = "white",
            force.left.label = TRUE,
            bottom.label = labCol,
            bottom.label.col = "white",
            bottom.label.text.angle = 90,
            bottom.label.text.size = 2.5,
            force.bottom.label = TRUE,
            print.plot = TRUE,
            legend.text.size = 10, # default 12
            legend.width = 2.0,
            left.label.text.size = 2.5)

  return(pp)
}


#' Performs a Principal Component Analysis (PCA) based on a kinship matrix K
#'
#' Generates a PCA and summary statistics from a given kinship matrix for
#' population structure. This matrix
#' can be a pedigree-based relationship matrix \eqn{\boldsymbol{A}}, a genomic
#' relationship matrix \eqn{\boldsymbol{G}} or a hybrid relationship matrix
#' \eqn{\boldsymbol{H}}. Individual names should be assigned to \code{rownames} and
#' \code{colnames}. There is additional output such as plots and other data frames
#' to be used on other downstream analyses (such as GWAS).
#'
#' It calls function \code{eigen()} to obtain eigenvalues and later generate the PCA and the
#' \code{factoextra} R package to extract and visualize results.
#'
#' @param K Input of a kinship matrix in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param scale If \code{TRUE} the PCA analysis will scale the kinship matrix, otherwise
#' it is used in its original scale (default = \code{TRUE}).
#' @param label If \code{TRUE} then includes in output individuals names (default = \code{FALSE}).
#' @param ncp The number of PC dimensions to be shown in the screeplot, and to provide
#' in the output data frame (default = \code{10}).
#' @param groups Specifies a vector of class factor that will be used to define different
#' colors for individuals in the PCA plot. It must be presented in the same order as the individuals
#' in the kinship matrix (default = \code{NULL}).
#' @param ellipses If \code{TRUE}, ellipses will will be drawn around each of the define levels in
#' \code{groups} (default = \code{FALSE}).
#'
#' @return A list with the following four elements:
#' \itemize{
#' \item \code{eigenvalues}: a data frame with the eigenvalues and its variances associated with each dimension
#' including only the first \code{ncp} dimensions.
#' \item \code{pca.scores}: a data frame with scores (rotated observations on the new components) including
#' only the first \code{ncp} dimensions.
#' \item \code{plot.pca}: a scatterplot with the first two-dimensions (PC1 and PC2) and their scores.
#' \item \code{plot.scree}: a barchart with the percentage of variances explained by the \code{ncp} dimensions.
#' }
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#' G[1:5, 1:5]
#'
#' # Perform the PCA.
#' G_pca <- kinship.pca(K = G, ncp = 10)
#' ls(G_pca)
#' G_pca$eigenvalues
#' head(G_pca$pca.scores)
#' G_pca$plot.pca
#' G_pca$plot.scree
#'
#' # PCA plot by family (17 groups).
#' grp <- as.factor(pheno.apple$Family)
#' G_pca_grp <- kinship.pca(K = G, groups = grp, label = FALSE, ellipses = FALSE)
#' G_pca_grp$plot.pca
#'

kinship.pca <- function(K=NULL, scale = TRUE, label = FALSE, ncp = 10,
                        groups = NULL, ellipses = FALSE){

  # Check if the class of K is matrix
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix.")
  }
  # Check the rownames/colnames
  if (is.null(rownames(K))){
    stop("Individual names not assigned to rows of matrix K.")
  }
  if (is.null(colnames(K))){
    stop('Individual names not assigned to columns of matrix K.')
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }
  # Check if the number of individuals is greater than 3
  if (nrow(K) < 3 | ncol(K) < 3){
      stop("Matrix K needs at least 3 individuals.")
  }
  # Check if the are missing values
  if (any(is.na(K))){
    stop("K matrix contains some missing data.")
  }
  if (ncp < 0 | ncp > nrow(K)) {
    stop("Value ncp must be positive and smaller than the number of rows in matrix K.")
  }

  # Generating the pca
  if (scale) {
    K <- cov2cor(K)
  }

  Peig <- eigen(K) # PCA-eigenvalues
  loadings <- Peig$vectors  # Loadings
  PCAcoord <- as.matrix(K) %*% as.matrix(loadings)  # New PCA coordinates
  colnames(PCAcoord) <- paste('PC',c(1:ncol(K)),sep='')
  colnames(loadings) <- paste('PC',c(1:ncol(K)),sep='')
  rownames(loadings) <- rownames(K)

  # Preparing the prcomp class object
  sdev <- sqrt(Peig$values)
  rotation <- loadings
  x <- PCAcoord
  pca <- list(sdev=sdev, rotation=rotation, center=scale, scale=scale, x=x)
  class(pca) <- "prcomp"

  # Percentage of variances explained by each principal component
  scree_plot <- fviz_eig(pca, addlabels=TRUE, ncp=ncp,
                                     barfill = "#0072B2",
                                     barcolor = "#0072B2",
                                     ggtheme = theme_classic())
  # Extract the eigenvalues/variances of the principal dimensions
  eig_var <- get_eig(pca)

  # Plot PCA
  if (isTRUE(label)) {
    if(is.null(groups)) {
       pca_plot <- fviz_pca_ind(pca, geom=c("point","text"),
                                            repel=TRUE,
                                            col.ind = "#0072B2",
                                            ggtheme = theme_classic())
    } else {
       pca_plot <- fviz_pca_ind(pca,
                        geom = c("point","text"),
                        repel = TRUE,
                        col.ind = groups, # color by groups
                        mean.point = FALSE,
                        legend.title = "Groups",
                        ggtheme = theme_classic())
    }
  }
  if (isFALSE(label)) {
    if(is.null(groups)) {
       pca_plot <- fviz_pca_ind(pca, geom="point",
                                            col.ind = "#0072B2",
                                            ggtheme = theme_classic())
    } else {
       pca_plot <- fviz_pca_ind(pca,
                        geom = "point",
                        col.ind = groups, # color by groups,
                        mean.point = FALSE,
                        legend.title = "Groups",
                        ggtheme = theme_classic())
    }
  }

  # Process ellipses if requested.
  if (!is.null(groups) & ellipses){

    # TODO this can be more memory efficient.

    group.comps <- cbind.data.frame(pca$x[, 1:2], groups)

    # Get centroids.
    centroids <- aggregate(cbind(PC1, PC2) ~ groups, data = group.comps, FUN = mean)

    # Get ellipses.
    ellipses.data  <- do.call(
      rbind,
      lapply(unique(group.comps$groups), function(t) {

        data.frame(
          groups = as.character(t),
          ellipse(
            cov(group.comps[group.comps$groups == t, 1:2]),
            centre = as.matrix(centroids[t, 2:3]), level = 0.95),
          stringsAsFactors=FALSE)
      }
      )
    )

    # Add ellipses to plot.
    pca_plot <- pca_plot +
      geom_path(data = ellipses.data, linewidth = .5, inherit.aes = F,
                aes(x = PC1, y = PC2, color = groups))
  }

  # Scores (rotated X observations on the new components) for ncp components
  scores <- pca$x[,c(1:ncp)]
  eigenvalues <- eig_var[c(1:ncp),]

  return(list(pca.scores=scores, eigenvalues=eigenvalues, plot.scree=scree_plot, plot.pca=pca_plot))

  }
#' Function that silences everything (e.g., \code{cat()}, \code{print()}, \code{message()}, ...)
#'
#' @param code code to be silenced.
#'
#' @return None.
#'
#' @keywords internal

silent_ <- function(code) {
  sink(tempfile()) ; on.exit(sink()) ; invisible(force(code))
}

#' Creates a dummy map if not provided
#'
#' @param marker.id vector with names of markers to compose dummy map.
#' @param message logical value indicating whether diagnostic messages should be printed on screen (default = \code{TRUE}).
#'
#' @return Data frame with dummy map. A single chromosome/linkage group is created and marker
#' distances are a sequence from one to the number of markers.
#'
#' @keywords internal

dummy.map_ <- function(marker.id =  NULL, message = TRUE) {

  # Report if required.
  if (message) message("Creating dummy map.")

  # Get dummy map.
  map <- data.frame(marker = marker.id,
                    chrom = 1,
                    pos = seq_along(marker.id))
}


#' Assess condition of the inverse of \strong{K}
#'
#' @param Kinv An inverse relationship matrix in full or sparse form.
#'
#' @return An object of class character with the condition of the matrix:
#' well-conditioned or ill-conditioned.
#'
#' @keywords internal

# TODO check with can be further improved.
Kinv.condition <- function(Kinv){

  if(nrow(Kinv) != ncol(Kinv)){
    Kinv <- sparse2full(Kinv)
  }

  n <- ncol(Kinv)
  CN.1 <- Kinv[1, 2]/sqrt(Kinv[1, 1] * Kinv[1, 1])
  CN.N <- Kinv[(n - 1), n]/sqrt(Kinv[(n - 1), (n - 1)] * Kinv[n, n])
  max_diag <- abs(max(diag(Kinv)))
  max_off_diag <- max(abs(Kinv - diag(Kinv)))
  if (abs(CN.1) > 0.99 | abs(CN.N) > 0.99 | max_diag > 1000 | max_off_diag > 1000) {
    status <- "ill-conditioned"
  }
  else {
    status <- "well-conditioned"
  }

  return(status)
}
#' Check data class
#'
#' @param data_ Data to be checked.
#' @param class_ The expected class of data_.
#'
#' @keywords internal

check.data_ <- function(data_ = NULL,
                        class_ = NULL){

  object_ <- get(data_, envir = parent.frame())

  # Test if data class is compliant.
  if( !any(class(object_) %in% class_) ) # If data has different classes: break
    stop(paste0("The ", data_, " argument should be of class(es) ",
                paste0(class_, collapse = " or ") , "."))

}

#' Check data mode
#'
#' @param data_ Data to be checked.
#' @param mode_ The expected mode of data_.
#'
#' @keywords internal

check.data.mode_ <- function(data_ = NULL,
                        mode_ = NULL){

  object_ <- get(data_, envir = parent.frame())

  # Test if data mode is compliant.
  if( !any(mode(object_) %in% mode_) ) # If data has different classes: break
    stop(paste0("The ", data_, " argument should be of mode(es) \'",
                paste0(mode_, collapse = " or ") , "'."))
}

# #' Check class of objects inside list
# #'
# #' @param data_ Data to be checked.
# #' @param class_ The expected class of data_.
# #'
# #' @keywords internal
#
# check.objects.list_ <- function(data_ = NULL,
#                         class_ = NULL){
#
#   object_ <- get(data_, envir = parent.frame())
#
#   # Test if data class is compliant.
#   if( !any(class(object_) %in% class_) ) # If data has different classes: break
#     stop(paste0("Objects inside list \'", data_, "' should be of class(es) ",
#                 paste0(class_, collapse = " or ") , "."))
#
# }

#' Check logical arguments
#'
#' @param arg_ The boolean argument to be checked.
#'
#' @keywords internal

check.logical_ <- function(arg_ = NULL){

  # Check if logical and stop if not.
  if( !is.logical( get(arg_, envir = parent.frame()) ) )
    stop(paste0("The value of \'", arg_, "' argument should be of class \'logical' (TRUE or FALSE)."))
}

#' Check string arguments
#'
#' @param data_ Data to be checked.
#' @param mandatory_ If the argument is mandatory for the analysis.
#' @param arg_ The string with the name of the function argument (e.g., \code{"gen"}).
#' @param rename_ If the respective column should be renamed to the argument name (e.g., \code{"genotype"} to \code{"gen"}).
#' @param class_ The expected class of the variable in data.
#' @param class.action_ The action to be taken if the variable has the wrong class.
#' Options are: \code{"message"}, \code{"warning"}, \code{"stop"}.
#' @param message_ If \code{class.action_ == "message"}, write \code{message = "message"} to capture upstream message command.
#' @param arg.update_ If the value passed to the argument should be updated to the argument name (e.g., if \code{gen = "geno"}, then \code{gen == "gen"}).
#'
#' @details This functions uses the \code{get} and \code{assign} which are need access to
#' objects that are one environment up on the hierarchy. The \code{envir} is set to
#' \code{parent.frame}. If the function is looking for something two or more environments up,
#' the arguments of \code{parent.frame} have to be changed.
#'
#' @keywords internal

check.args_ <- function(data_ = NULL,
                        mandatory_ = FALSE,
                        arg_ = NULL,
                        class_ = NULL,
                        class.action_ = NULL,
                        message_ = message){

  # Capture relevant info.
  data.frame_ <- get(data_, envir = parent.frame())
  real.var.value_ <- get(arg_, envir = parent.frame())
  class.fun_ <- paste0("is.", class_)

  # Evaluate if arg is not null.
  if( !is.null(real.var.value_) ){

    # Check real.var.value_ is in data (mandatory stop).
    if( !real.var.value_ %in% names(data.frame_) )
      stop(paste0("\'", real.var.value_,
                  "' does not correspond to a variable name of \'pheno.data'."))

    # Check class of variable in data.
    if( !getFunction(class.fun_)(data.frame_[[real.var.value_]]) )

      # Do action unless action is message and message is FALSE.
      if( !(class.action_ == "message" & isFALSE(message_)) )
        getFunction(class.action_)(
          paste0("Variable \'", real.var.value_,
                 "' should be of class \'", class_, "'."))

  }

  # Evaluate if arg is null.
  if( is.null(real.var.value_) & mandatory_ )
      stop(paste0("The argument \'", arg_, "' is mandatory."))
}
#' Generates a sparse form matrix from a full form matrix
#'
#' Modifies the input square matrix into its sparse form with three columns per line,
#' corresponding to the set: \code{Row, Col, Value}. This matrix defines the lower triangle
#' row-wise of the original matrix and it is sorted as columns within row.
#' Values of zero on the matrix are dropped by default. Individual
#' names should be assigned to \code{rownames} and \code{colnames}.
#'
#' Based on the function published by Borgognone \emph{et al.} (2016).
#'
#' @param K A square matrix in full form (\eqn{n \times n}) (default = \code{NULL}).
#' @param drop.zero If \code{TRUE} observations equal to zero are dropped
#' (default = \code{TRUE}).
#'
#' @return A matrix in sparse form with columns: \code{Row, Col, Value}
#' together with the attributes \code{rowNames} and \code{colNames}.
#' If attribute \code{INVERSE} is found this is also passed to the sparse matrix.
#'
#' @references
#' Borgognone, M.G., Butler, D.G., Ogbonnaya, F.C and Dreccer, M.F. 2016.
#' Molecular marker information in the analysis of multi-environment trial helps
#' differentiate superior genotypes from promising parents. Crop Science 56:2612-2628.
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = FALSE)$G
#' G[1:5, 1:5]
#'
#' # Transform matrix into sparse.
#' G.sparse <- full2sparse(K = G)
#' head(G.sparse)
#' head(attr(G.sparse, "rowNames"))
#' head(attr(G.sparse, "colNames"))
#'


full2sparse <- function (K = NULL, drop.zero = TRUE) {

  if (is.null(K) || !inherits(K, "matrix")) {
    stop('K should be a valid object of class matrix.')
  }
  # Check the attributes of K
  if (is.null(rownames(K))){
    stop('Individual names not assigned to rows of matrix K.')
  }
  if (is.null(colnames(K))){
    stop('Individual names not assigned to columns of matrix K.')
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }

  INVERSE <- attr(K, "INVERSE")

  if(drop.zero) {
    which <- (K != 0 & lower.tri(K, diag = TRUE))
  } else {
    which <- lower.tri(K, diag = TRUE)
  }

  sparse.frame <- data.frame(Row = t(row(K))[t(which)], Col = t(col(K))[t(which)],
                   Value = t(K)[t(which)])

  sparse.frame <- as.matrix(sparse.frame)

  # Add attributes.
  attr(sparse.frame, "rowNames") <- rownames(K)
  attr(sparse.frame, "colNames") <- colnames(K)

  # Pass collected attributes.
  if (!is.null(INVERSE)) {attr(sparse.frame, "INVERSE") <- INVERSE}

  return(sparse.frame)

}
#' Genotypic data for apple dataset
#'
#' Genotypic data on 247 apple clones (\emph{i.e.}, genotypes) with a total of
#' 2,828 SNP markers (coded as 0, 1, 2 and there are no missing records).
#' Dataset obtained from supplementary material in Kumar \emph{et al.} (2015).
#'
#' @docType data
#'
#' @usage geno.apple
#'
#' @format matrix
#'
#' @keywords datasets
#'
#' @references
#' Kumar S., Molloy C., Muñoz P., Daetwyler H., Chagné D., and Volz R. 2015.
#' Genome-enabled estimates of additive and nonadditive genetic variances and prediction
#' of apple phenotypes across environments. G3 Genes, Genomes, Genetics 5:2711-2718.
#'
#' @examples
#' geno.apple[1:5, 1:5]
#'
#' @name geno.apple
NULL
#' Genotypic data of 655 genotypes for loblolly pine dataset
#'
#' Genotypic data for a total of 4,853 SNPs (coded as 0, 1, 2 and -9 for missing)
#' on 655 genotypes of Loblolly Pine (\emph{Pinus taeda} L.).
#' Dataset modified from supplementary material from Resende \emph{et al.} (2012).
#' This dataset differs from the original as some genotypes were made artificially
#' missing by full-sib family.
#'
#' @docType data
#'
#' @usage geno.pine655
#'
#' @format matrix
#'
#' @keywords datasets
#'
#' @references
#' Resende, M.F.R., Munoz, P. Resende, M.D.V., Garrick, D.J., Fernando, R.L., Davis, J.M.,
#' Jokela, E.J., Martin, T.A., Peter, G.F., and Kirst, M. 2012. Accuracy of genomic
#' selection methods in a standard data set of loblolly pine (\emph{Pinus taeda} L.).
#' Genetics 190:1503-1510.
#'
#' @examples
#' geno.pine655[1:5, 1:5]
#'
#' @name geno.pine655
NULL
#' Genotypic data of 926 genotypes for loblolly pine dataset
#'
#' Genotypic data for a total of 4,853 SNPs (coded as 0, 1, 2 and -9 for missing)
#' on 926 genotypes of Loblolly Pine (\emph{Pinus taeda} L.).
#' Dataset obtained from supplementary material in Resende \emph{et al.} (2012).
#'
#' @docType data
#'
#' @usage geno.pine926
#'
#' @format matrix
#'
#' @keywords datasets
#'
#' @references
#' Resende, M.F.R., Munoz, P. Resende, M.D.V., Garrick, D.J., Fernando, R.L., Davis, J.M.,
#' Jokela, E.J., Martin, T.A., Peter, G.F., and Kirst, M. 2012. Accuracy of genomic
#' selection methods in a standard data set of loblolly pine (\emph{Pinus taeda} L.).
#' Genetics 190:1503-1510.
#'
#' @examples
#' geno.pine926[1:5, 1:5]
#'
#' @name geno.pine926
NULL
#' Genotypic data for Atlantic salmon dataset
#'
#' Genotypic data on 1,481 Atlantic salmon samples. A total of 17,156
#' SNP markers (coded as 0, 1, 2 and \code{NA} for missing) are included in this dataset.
#' Dataset obtained from supplementary material in Robledo \emph{et al.} (2018).
#'
#' @docType data
#'
#' @usage geno.salmon
#'
#' @format matrix
#'
#' @keywords datasets
#'
#' @references
#' Robledo D., Matika O., Hamilton A., and Houston R.D. 2018.
#' Genome-wide association and genomic selection for resistance
#' to amoebic gill disease in Atlantic salmon.
#' G3 Genes, Genomes, Genetics 8:1195-1203.
#'
#' @examples
#' geno.salmon[1:5, 1:5]
#'
#' @name geno.salmon
NULL
#' Estimates minor allele frequency (MAF)
#'
#' @param M The additive \eqn{n \times p} matrix coded with 0, 1, 2, or \code{NA}
#' (default = \code{NULL}).
#'
#' @return A vector with the MAF of each molecular marker
#'
#' @keywords internal

maf <- function(M = NULL){

  # Get the frequency of the markers.
  maf <- colMeans(M, na.rm = TRUE)/2

  # Identify the MAF.
  maf <- apply(cbind(maf, 1 - maf), 1, min)

  return(maf)
}

#' Estimates observed and expected heterozygosity
#'
#' @param M The additive \eqn{n \times p} matrix coded with 0, 1, 2 coding (default = \code{NULL}).
#'
#' @return A list with vectors containing the observed (ho) and expected (he) heterozygosity.
#'
#' @keywords internal

heterozygosity <- function(M = NULL){

  # Get q.
  q <- maf(M)

  # Get p.
  p <- 1 - q

  # Get the expected heterozygosity.
  he <- 2 * p * q

  # Get the obseved heterozygosity.
  ho <- colMeans(M == 1, na.rm = TRUE)

  return(data.frame(ho = ho, he = he))
}

#' Estimates call rate
#'
#' @param M The additive \eqn{n \times p} matrix with any coding (default = \code{NULL}).
#' @param margin A character indicating the margin for call rate calculations.
#' Options are: \code{row} and \code{col} (default = \code{row}).
#'
#' @return A vector containing the call rate.
#'
#' @keywords internal

callrate <- function(M = NULL, margin = c("row", "col")){

  # Collect input.
  margin <- match.arg(margin)

  # CR by row.
  if (margin == "row") cr <- 100 - rowSums(is.na(M))/ncol(M) * 100

  # CR by col.
  if (margin == "col") cr <- 100 - colSums(is.na(M))/nrow(M) * 100

  return(cr)

}

#' Estimates the population level inbreeding (Fis) by marker
#'
#' @param M The additive \eqn{n \times p} matrix with any coding (default = \code{NULL}).
#' @param margin A character indicating the margin for call rate calculations.
#' Options are: \code{row} and \code{col} (default = \code{col}).
#'
#' @return A vector containing the Fis for the markers.
#'
#' @keywords internal

Fis <- function(M = NULL, margin = c("col", "row")){

  # Collect input.
  margin <- match.arg(margin)

  # H by row.
  if (margin == "col") H <- heterozygosity(M = M)

  # H by col.
  if (margin == "row") H <- heterozygosity(M = t(M))

  # Calculate Fis.
  Fis <- ifelse(test = H[, "he"] == 0,
                yes = 0,
                no = 1 - (H[,"ho"] / H[,"he"]))

  return(Fis)

}
#' Reports summary statistics, plots and filter options for a given kinship matrix K
#'
#' It reports summary statistics, plots and allows for some filter options
#' for diagonal and off-diagonal elements for a given kinship matrix.
#' The input matrix can be a pedigree-based
#' relationship matrix \eqn{\boldsymbol{A}}, a genomic relationship matrix \eqn{\boldsymbol{G}} or a
#' hybrid relationship matrix \eqn{\boldsymbol{H}}.
#' Individual names should be assigned to \code{rownames} and \code{colnames}.
#'
#' @param K Input of a kinship matrix in full format (\eqn{n \times n}) (default = \code{NULL}).
#' @param diagonal.thr.large A threshold value to flag large diagonal values (default = \code{1.2}).
#' @param diagonal.thr.small A threshold value to flag small diagonal values (default = \code{0.8}).
#' @param duplicate.thr A threshold value to flag possible duplicates. Any calculation larger than the
#' threshold based on
#' \eqn{\boldsymbol{k}_{i,i}\mathbin{/}\sqrt{\boldsymbol{k}_{i,i} \times \boldsymbol{k}_{j,j}}}
#' is identified as a duplicate (default = \code{0.95}).
#' @param clean.diagonal If \code{TRUE} returns a kinship matrix filtered by values smaller than
#' \code{diagonal.thr.large} and larger than \code{diagonal.thr.small} (default = \code{FALSE}).
#' @param clean.duplicate If \code{TRUE} return a kinship matrix without the flagged duplicate individuals.
#' All individuals involved are removed (default = \code{FALSE}).
#' @param plots If \code{TRUE} generates graphical output of the diagonal and off-diagonal
#' values of the kinship matrix (default = \code{TRUE}).
#' @param sample.plot A numeric value between 0 and 1 indicating the proportion
#' of the data points to be sampled for fast plotting of off-diagonal values.
#' Note that for proportions other than 1, the method is not exact and low
#' proportions are only recommended for large kinship matrices (default = \code{1}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{list.diagonal}: a data frame with the list of flagged large or small diagonal values.
#' \item \code{list.duplicate}: a data frame with the list of possible duplicates.
#' \item \code{clean.kinship}: output of kinship matrix filtered without the flagged diagonal
#'  and/or duplicate individuals.
#' \item \code{plot.diagonal}: histogram with the distribution of diagonal values from the kinship matrix.
#' \item \code{plot.offdiag}: histogram with the distribution of off-diagonal values from kinship matrix.
#' }
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.apple, method = "VanRaden")$G
#'
#' # Diagnose G.
#' G_summary <- kinship.diagnostics(
#'  K = G,
#'  diagonal.thr.large = 1.3, diagonal.thr.small = 0.7, clean.diagonal = TRUE,
#'  duplicate.thr = 0.8, clean.duplicate = TRUE,
#'  sample.plot = 0.50)
#' ls(G_summary)
#' dim(G_summary$clean.kinship)
#' G_summary$clean.kinship[1:5, 1:5]
#' G_summary$list.duplicate
#' G_summary$list.diagonal
#' G_summary$plot.diag
#' G_summary$plot.offdiag
#'

kinship.diagnostics <- function(K = NULL,
                                diagonal.thr.large = 1.2, diagonal.thr.small = 0.8,
                                duplicate.thr = 0.95, clean.diagonal = FALSE,
                                clean.duplicate = FALSE,
                                plots = TRUE, sample.plot = 1, message = TRUE){

  # Check if the class of K is matrix.
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix")
  }
  # Check the attributes of K
  if (is.null(rownames(K))){
    stop('Individual names not assigned to rows of matrix K.')
  }
  if (is.null(colnames(K))){
    stop('Individual names not assigned to columns of matrix K.')
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }
  # Check on other input
  if (duplicate.thr < 0 | duplicate.thr > 1) {
    stop("Specification of duplicate.thr must be between 0 and 1.")
  }
  if (diagonal.thr.large < diagonal.thr.small) {
    stop("Value of diagonal.thr.large has to be equal or larger than diagonal.thr.small.")
  }
  if (diagonal.thr.large < 0 | diagonal.thr.small < 0) {
    stop("Values of diagonal.thr.large and diagonal.thr.small have the be both positive.")
  }
  if (sample.plot <= 0 | sample.plot > 1) {
    stop("Values of sample.plot must be between 0 and 1.")
  }

  # Preparing submatrices
  n <- nrow(K)
  #indNames <- rownames(K)
  diagK <- diag(K)  # Vector of diagonal
  #offdiag <- K[lower.tri(K, diag=FALSE)]
  Kcorr <- cov2cor(K)
  K.sparse <- full2sparse(K)
  corrS <- full2sparse(Kcorr)
  K.sparse <- data.frame(K.sparse, Corr=corrS[,3])
  offK <- K.sparse[K.sparse$Row != K.sparse$Col,]
  rm(K.sparse, Kcorr, corrS)

  # SOME STATISTICS
  # Some general statistics and report in 'rp'
  if (message){
    message("Matrix dimension is: ", n, "x", n)
  }

  rank <- try(qr(K)$rank, silent=TRUE)

  # Commented out
  #if (message){
  #  if(class(rank) == "try-error"){
  #    message("Rank cannot be obtained due to missing values!")
  #  } else {
  #    message("Rank of matrix is: ",rank)
  #  }
  #}

  range.diagonal <- c(min=min(diagK, na.rm=TRUE), max=max(diagK, na.rm=TRUE))
  if (message){
    message("Range diagonal values: ", round(range.diagonal[1], 5),
            " to ",round(range.diagonal[2], 5))
  }
  mean.diag <- mean(diagK, na.rm=TRUE)
  if (message){
    message("Mean diagonal values: ", round(mean.diag, 5))
  }
  range.off.diagonal <- c(min=min(offK$Value, na.rm=TRUE), max=max(offK$Value, na.rm=TRUE))
  if (message){
    message("Range off-diagonal values: ", round(range.off.diagonal[1], 5),
            " to ",round(range.off.diagonal[2], 5))
  }
  mean.off.diag <- mean(offK$Value, na.rm=TRUE)
  if (message){
    message("Mean off-diagonal values: ", round(mean.off.diag, 5))
  }

  ##########################
  # DEALING with diagonal

  df.list.diag <- data.frame(
    value = sort(diagK[diagK > diagonal.thr.large
                       | diagK < diagonal.thr.small], decreasing=TRUE))

  # Generating list of flagged potential duplicates
  df.list.duplicate <- offK[offK$Corr > duplicate.thr,]
  df.list.duplicate <- data.frame(df.list.duplicate, Indiv.A=rownames(K)[df.list.duplicate$Row],
                                  Indiv.B=colnames(K)[df.list.duplicate$Col])
  rownames(df.list.duplicate) <- NULL
  df.list.duplicate <- df.list.duplicate[,c(5,6,3,4)]
  df.list.duplicate <- df.list.duplicate[order(df.list.duplicate$Corr, decreasing=TRUE),]
  rownames(df.list.duplicate) <- NULL

  # Generating new K if requested with clean.diagonal.
  if (clean.diagonal){
    if (nrow(df.list.diag) > 0){
      Kclean <- K[-which(rownames(K) %in% row.names(df.list.diag)),
                  -which(rownames(K) %in% row.names(df.list.diag))]
    } else {
      if (message){
        message("No individuals filtered out by the diagonal thresholds, as none were found.")
      }

      # I added this line because it will be NULL if nothing has changed.
      Kclean <- NULL
    }
  } else {
    Kclean <- NULL
  }

  # Generating new K if requested with clean.duplicate.
  if (clean.duplicate){
    if (nrow(df.list.duplicate) > 0){
      idx.offdiag <- unique(c(df.list.duplicate$Indiv.A, df.list.duplicate$Indiv.B))
      # idx.offdiag <- unique(df.list.duplicate$Indiv.A)
      if (is.null(Kclean)) {
        Kclean <- K[-which(rownames(K) %in% idx.offdiag),
                    -which(rownames(K) %in% idx.offdiag)]
      } else {
        Kclean <- Kclean[-which(rownames(Kclean) %in% idx.offdiag),
                         -which(rownames(Kclean) %in% idx.offdiag)]
      }
    } else {
      if (message){
        message("No individuals filtered by duplicate.thr = ", duplicate.thr,", as none were found.")
      }
    }
  }

  # Removing K
  rm(K)

  # SOME STATISTICS
  # Report extreme cases
  count <- length(diagK[diagK > diagonal.thr.large | diagK < diagonal.thr.small])
  if (message){
    message("There are ", count, " extreme diagonal values, outside < ", diagonal.thr.small,
            " and > ", diagonal.thr.large)
  }
  count <- nrow(df.list.duplicate)
  if (message){
    message("There are ", count, " records of possible duplicates, based on: k(i,j)/sqrt[k(i,i)*k(j,j)] >  ", duplicate.thr)
  }

  # Obtain histogram of diagonal ----------------------------------------------------------------

  if (plots){

    # Get data.frame for plotting.
    diagK <- as.data.frame(diagK)
    names(diagK)  <- "value"

    # Generate main plot.
    p1 <- ggplot(diagK, aes(x=value)) +
      geom_histogram(aes(y = after_stat(density)), fill='#0072B2', bins=40) +
      geom_density(alpha=0.3, fill="grey", position='identity') +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            plot.title = element_text(face = "bold")) +
      ggtitle("Diagonal Values")

    # Generate boxplot.
    p2 <- ggplot(aes(value), data = diagK) +
      geom_boxplot() +
      theme_void()

    # Combine plots.
    plot.diag <- plot_grid(p1, p2, rel_heights = c(1, 0.2), ncol = 1, nrow = 2, align = "hv")

    # Obtain histogram of off-diagonal ------------------------------------------------------------

    # Sampling for plot if requested.
    if (sample.plot != 1){
      offK <- offK[sample(
        x = 1:nrow(offK),
        size = floor(sample.plot * nrow(offK)),
        replace = FALSE), , drop = FALSE]
    }

    p1 <- ggplot(offK, aes(x = Value)) +
      geom_histogram(aes(y = after_stat(density)), fill = '#0072B2', bins = 40) +
      geom_density(alpha = 0.3, fill = "grey", position = 'identity') +
      theme_classic() +
      theme(axis.title.x=element_blank(),
            plot.title = element_text(face = "bold")) +
      ggtitle("Off-diagonal Values")

    # Boxplot
    p2 <- ggplot(aes(Value), data = offK) +
      geom_boxplot() +
      theme_void()

    # Combine plots.
    plot.offdiag <- plot_grid(p1, p2, rel_heights = c(1, 0.2), ncol = 1, nrow = 2, align = "hv")

  } else {
    # Nulify plots if not requested.
    plot.diag <- plot.offdiag <- NULL
  }


  # Finalize ------------------------------------------------------------------------------------

  if (nrow(df.list.duplicate) == 0) {df.list.duplicate <- NULL}
  if (nrow(df.list.diag) == 0) {df.list.diag <- NULL}

  return(list(list.diagonal=df.list.diag, list.duplicate=df.list.duplicate,
              clean.kinship=Kclean, plot.diag=plot.diag, plot.offdiag=plot.offdiag))

}
#' Check the genomic relationship matrix G against
#' the pedigree relationship matrix A or vice versa
#'
#' Assesses a given genomic relationship matrix \eqn{\boldsymbol{G}} against the
#' pedigree relationship matrix \eqn{\boldsymbol{A}}, or vice versa,
#' to determine the matched and mismatched individuals.
#' If requested, it provides the cleaned versions containing only the matched individuals
#' between both matrices. The user should provide the matrices \eqn{\boldsymbol{G}}and
#' \eqn{\boldsymbol{A}} in full form (\eqn{ng \times ng} and \eqn{na \times na}, respectively).
#' Individual names should be assigned to \code{rownames} and \code{colnames} for both matrices.
#'
#' @param G Input of the genomic relationship matrix \eqn{\boldsymbol{G}} in full form (\eqn{ng \times ng}) (default = \code{NULL}).
#' @param A Input of the pedigree relationship matrix \eqn{\boldsymbol{A}} in full form (\eqn{na \times na}) (default = \code{NULL}).
#' @param clean If \code{TRUE} generates new clean \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}}
#' matrices in full form containing only matched individuals (default = \code{TRUE}).
#' @param ord If \code{TRUE} it will order by ascending order of individual names
#' both of the clean \eqn{\boldsymbol{A}} and \eqn{\boldsymbol{G}} matrices (default = \code{TRUE}).
#' @param mism If \code{TRUE} generates two data frames with mismatched individual names
#' from the \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}} matrices (default = \code{FALSE}).
#' @param RMdiff If \code{TRUE} it generates the matrix (in lower diagonal row-wise sparse form) of matched
#' observations from both the \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}} matrices.
#' This matrix can be used to identify inconsistent values between matched matrices, but it can be very large
#' (default = \code{FALSE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{Gclean}: a matrix with the portion of \eqn{\boldsymbol{G}} containing only matched individuals.
#' \item \code{Aclean}: a matrix with the portion of \eqn{\boldsymbol{A}} containing only matched individuals.
#' \item \code{mismG}: a vector containing the names of the individuals from matrix \eqn{\boldsymbol{G}} that are
#' missing in matrix \eqn{\boldsymbol{A}}.
#' \item \code{mismA}: a vector containing the names of the individuals from matrix \eqn{\boldsymbol{A}} that are
#' missing in matrix \eqn{\boldsymbol{G}}.
#' \item \code{RM}: a data frame with the observations from both the \eqn{\boldsymbol{G}} and \eqn{\boldsymbol{A}}
#' matched matrices, together with their absolute relationship difference.
#' \item \code{plotG2A}: scatterplot with the pairing of matched pedigree- against genomic-based
#' relationship values.
#' This graph might take a long to plot with large datasets.
#' }
#'
#' @export
#'
#' @examples
#' \donttest{
#' # Get A matrix.
#' A <- AGHmatrix::Amatrix(data = ped.pine)
#' A[1:5,1:5]
#' dim(A)
#'
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.2, ind.callrate = 0.20,
#'  na.string = "-9",
#'  plots = FALSE)$M.clean
#'
#' # Get G matrix.
#' G <- G.matrix(M = M.clean, method = "VanRaden", na.string = "-9")$G
#' G[1:5, 1:5]
#' dim(G)
#'
#' # Match G2A.
#' check <- match.G2A(
#'  A = A, G = G,
#'  clean = TRUE, ord = TRUE, mism = TRUE, RMdiff = TRUE)
#' ls(check)
#' dim(check$Aclean)
#' dim(check$Gclean)
#' check$Aclean[1:5, 1:5]
#' check$Gclean[1:5, 1:5]
#' head(check$mismG)
#' head(check$mismA)
#' check$plotG2A
#' head(check$RM)
#' }
#'

match.G2A <- function(A = NULL, G = NULL, clean = TRUE, ord = TRUE,
                      mism = FALSE, RMdiff = FALSE, message = TRUE){

  if (is.null(A) || !inherits(A, "matrix")) {
    stop("A should be a valid object of class matrix.")
  }
  if (is.null(G) || !inherits(G, "matrix")) {
    stop("G should be a valid object of class matrix.")
  }
  # Check the rownames/colnames
  if (is.null(rownames(A))){
    stop("Individual names not assigned to rows of matrix A.")
  }
  if (is.null(colnames(A))){
    stop('Individual names not assigned to columns of matrix A.')
  }
  if ((identical(rownames(A), colnames(A))) == FALSE){
    stop("Rownames and colnames of matrix A do not match.")
  }
  if (is.null(rownames(G))){
    stop("Individual names not assigned to rows of matrix G.")
  }
  if (is.null(colnames(G))){
    stop("Individual names not assigned to columns of matrix G.")
  }
  if ((identical(rownames(G), colnames(G))) == FALSE){
    stop("Rownames and colnames of matrix G do not match.")
  }

  # Check for consistency between A and G
  Aind <- row.names(A)
  Gind <- row.names(G)

  if (all(Aind %in% Gind) & message){
    message("All ", ncol(A), " individuals from matrix A match those individuals from matrix G.")
  }
  if (all(Gind %in% Aind) & message){
    message("All ", ncol(G), " individuals from matrix G match those individuals from matrix A.")
  }

  notGenotyped <- which((Aind %in% Gind) == FALSE)
  notPedigree <- which((Gind %in% Aind) == FALSE)

  # If G has different individuals than A, remove these individuals from G
  if (length(notPedigree) > 0) {
    if (message){
      message("Matrix G has ", length(notPedigree), " individuals (out of ", ncol(G), ") NOT present on matrix A.")
    }
    if (clean) {
      Gclean <- G[-notPedigree,-notPedigree]
    } else {
      Gclean <- NULL
    }
    if (mism) {
      rpG <- rownames(G[notPedigree,])
    } else {
      rpG <- NULL
    }
  } else {
    Gclean <- G
    rpG <- NULL
  }

  # If A has different ind than G, remove these ind from A
  if (length(notGenotyped) > 0) {
    if (message){
      message("Matrix A has ", length(notGenotyped), " individuals (out of ", ncol(A), ") NOT present on matrix G.")
    }
    if (clean) {
      Aclean <- A[-notGenotyped, -notGenotyped]
    } else {
      Aclean <- NULL
    }
    if (mism) {
      rpP <- rownames(A[notGenotyped,])
    } else {
      rpP <- NULL
    }
  } else {
    Aclean <- A
    rpP <- NULL
  }

  # Check order of A and G
  if (all(row.names(Aclean) == row.names(Gclean)) == FALSE){
    if (!ord) {
      if (message){
        message("Order of individual names from matched matrices A and G DO NOT agree.")
      }
    }
    if (ord){
      Gclean <- Gclean[order(rownames(Gclean), decreasing=FALSE),
                       order(colnames(Gclean), decreasing=FALSE)]
      Aclean <- Aclean[order(rownames(Aclean), decreasing=FALSE),
                       order(colnames(Aclean), decreasing=FALSE)]
    }
  }

  # TODO this section has to be improved since A.sparse can also be NULL!

  # Sparse form matrices (faster for plots useful for RM matrix)
  G.sparse <- full2sparse(K=Gclean, drop.zero=FALSE)
  A.sparse <- full2sparse(K=Aclean, drop.zero=FALSE)

  # Generating plot of Aclean vs Gclean
  LL <- min(A.sparse[,3], G.sparse[,3])
  UL <- max(A.sparse[,3], G.sparse[,3])

  # Improved version of plot (faster). GG
  p <- ggplot(data.frame(AValue = A.sparse[,3], GValue = G.sparse[,3]),
              aes(x = AValue, y = GValue, color = 'black')) +
         geom_scattermost(xy = cbind(A.sparse[,3], G.sparse[,3]), color = "#0072B2", pointsize = 2) +
         geom_abline(linetype = "dashed") +
         xlim(LL, UL) + ylim(LL, UL)+
         labs(x="Pedigree Relationship (A matrix)",
         y = "Genomic Relationship (G matrix)") +
         theme_classic()

    # RM matrix for diagnostics
  if (isTRUE(RMdiff)) {
    if (!isTRUE(clean)) {
      stop("Option clean must be TRUE to produce RM data frame.")
    }
    RM <- data.frame(A.sparse[,1:2], AValue=round(A.sparse[,3],6), GValue=round(G.sparse[,3],6))
    RM$absdiff <- round(abs(RM$AValue - RM$GValue),6)
    RM$IDRow <- rownames(Aclean)[RM$Row]
    RM$IDCol <- rownames(Aclean)[RM$Col]
    RM$Diag <- 0
    RM$Diag[RM$Row == RM$Col] <- 1
    RM <- RM[c(1,2,6,7,3,4,5,8)]
  } else {
    RM <- NULL
  }

  AValue <- GValue <-  NULL

  return(list(Aclean=Aclean, Gclean=Gclean, mismG=rpG, mismA=rpP, RM=RM, plotG2A=p))


}

# A.sparse <- rbind(A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse, A.sparse)
# A.sparse <- rbind(A.sparse, A.sparse)
# G.sparse <- rbind(G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse, G.sparse)
# G.sparse <- rbind(G.sparse, G.sparse)

#' Check any kinship matrix K against phenotypic data
#'
#' Assesses a given kinship matrix against the provided phenotypic data to determine
#' if all genotypes are in the kinship matrix or not. It also reports which individuals
#' match or are missing from one set or another.
#' If requested, a reduced kinship matrix is generated that has only the matched individuals.
#' The input kinship matrix can be a pedigree-based relationship matrix \eqn{\boldsymbol{A}},
#' a genomic-based relationship matrix \eqn{\boldsymbol{G}}, or a hybrid
#' relationship matrix \eqn{\boldsymbol{H}}.
#' Individual names should be assigned to \code{rownames} and \code{colnames} of input matrix.
#'
#' @param K Input of a kinship matrix in full form (\eqn{n \times n}) (default = \code{NULL});
#' @param pheno.data A data fame with the phenotypic data to assess (for \code{n} individuals)
#' (default = \code{NULL}).
#' @param indiv The string for the column name for genotypes/individuals in the phenotypic data (default = \code{NULL}).
#' @param clean If \code{TRUE}, generates a new clean kinship matrix containing only the matched
#' phenotyped individuals (default = \code{FALSE}).
#' @param ord If \code{TRUE}, it will order the kinship matrix as in the phenotypic data, which is
#' recommended for some downstream genomic analyses (default = \code{TRUE}).
#' @param mism If \code{TRUE}, it generates data frames with matched and mismatched individual's names
#' from the kinship matrix and the phenotypic data (default = \code{FALSE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#'  \itemize{
#' \item \code{mismatchesK}: a vector containing the names of the individuals from the provided kinship matrix
#' that \emph{mismatch} with the phenotypic data.
#' \item \code{matchesK}: a vector containing the names of the individuals from the provided kinship matrix
#' that \emph{match} with the phenotypic data.
#' \item \code{mismatchesP}: a vector containing the names of phenotyped individuals
#' that \emph{mismatch} with those from the kinship matrix.
#' \item \code{matchesP}: a vector containing the names of phenotyped individuals
#' that \emph{match} with those from the kinship matrix.
#' \item \code{Kclean}: a clean kinship matrix containing only the matched phenotyped individuals.
#' }
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' G <- G.matrix(M = geno.pine655, method = "VanRaden", na.string = "-9", sparseform = FALSE)$G
#' dim(G)
#'
#' # Match G and the phenotype.
#' check <-
#'  match.kinship2pheno(
#'   K = G, pheno.data = pheno.pine,
#'   indiv = "Genotype",
#'   clean = TRUE, mism = TRUE)
#' ls(check)
#' length(check$matchesK)
#' length(check$mismatchesK)
#' length(check$matchesP)
#' length(check$mismatchesP)
#' dim(check$Kclean)
#'


match.kinship2pheno <- function(K = NULL,
                                pheno.data = NULL, indiv = NULL,
                                clean = FALSE, ord = TRUE, mism = FALSE,
                                message = TRUE){

  # Check if the class of K is matrix
  if (is.null(K) || !inherits(K, "matrix")) {
    stop("K should be a valid object of class matrix")
  }
  if (is.null(rownames(K))){
    stop("Individual names not assigned to rows of matrix K")
  }
  if (is.null(colnames(K))){
    stop("Individual names not assigned to columns of matrix K")
  }
  if ((identical(rownames(K), colnames(K))) == FALSE){
    stop("Rownames and colnames of matrix K do not match.")
  }

  # Checks
  kinship_names <- rownames(K)
  pheno_names <- pheno.data[[indiv]]

  if( all(kinship_names %in% pheno_names) == TRUE){
    if (message){
      message("All individuals within the kinship matrix match the phenotyped individuals.")
    }
    mismatchesK <- NULL
    matchesK <- which(kinship_names %in% pheno_names == TRUE)
  } else {
    mismatchesK <- which(kinship_names %in% pheno_names != TRUE)
    if (message){
      message("Kinship matrix contains ", length(mismatchesK),
              " individuals that DO NOT match the phenotyped individuals.")
    }
    matchesK <- which(kinship_names %in% pheno_names == TRUE)
    if (message){
      message("Kinship matrix contains ", length(matchesK),
              " individuals that match the phenotyped individuals.")
    }
  }

  if( all(pheno_names %in% kinship_names) == TRUE){
    if (message){
      message("All phenotyped individuals match the individuals within the kinship matrix.")
    }
    mismatchesP <- NULL
    matchesP <- which(pheno_names %in% kinship_names == TRUE)
  } else {
    mismatchesP <- which(pheno_names %in% kinship_names != TRUE)
    if (message){
      message("Phenotypic data contains ", length(mismatchesP),
              " individuals that DO NOT match the kinship matrix individuals.")
    }
    matchesP <- which(pheno_names %in% kinship_names == TRUE)
    if (message){
    message("Phenotypic data contains ", length(matchesP),
            " individuals that match the kinship matrix individuals.")
    }
  }

  if(length(mismatchesK) != 0 & clean){
    if (message){
      message("Individuals within the kinship matrix that do not match those in the phenotypic data will be removed from this matrix.")
    }
    #Kclean <- K[matchesP, matchesP]
    Kclean <- K[matchesK, matchesK]
  } else {
    Kclean <- NULL
  }

  if(ord) {
    #list.order <- match(rownames(Kclean), pheno_names)
    list.order <- match(pheno_names, rownames(Kclean))
    Kclean <- Kclean[list.order,list.order]
  }
  if(!mism) {
    mismatchesK <- NULL; matchesK <- NULL
    mismatchesP <- NULL; matchesP <- NULL
  }

  return(list(mismatchesK=mismatchesK, matchesK=matchesK,
              mismatchesP=mismatchesP, matchesP=matchesP,
              Kclean=Kclean))

}
#' Pedigree data for loblolly pine dataset
#'
#' Individual pedigree data for a total of 2,034 records of loblolly pine
#' (\emph{Pinus taeda} L.).
#' Missing parental information coded as 0.
#' Dataset obtained from supplementary material in Resende \emph{et al.} (2012).
#'
#' @docType data
#'
#' @usage ped.pine
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references
#' Resende, M.F.R., Munoz, P. Resende, M.D.V., Garrick, D.J., Fernando, R.L., Davis, J.M.,
#' Jokela, E.J., Martin, T.A., Peter, G.F., and Kirst, M. 2012. Accuracy of genomic
#' selection methods in a standard data set of loblolly pine (\emph{Pinus taeda} L.).
#' Genetics 190:1503-1510.
#'
#' @examples
#' ped.pine |> head()
#'
#' @name ped.pine
NULL

#' Pedigree data for Atlantic salmon dataset
#'
#' Pedigree data of 1,481 Atlantic salmon samples.
#' Missing parental information coded as 0.
#' Dataset obtained from supplementary material in Robledo \emph{et al.} (2018).
#'
#' @docType data
#'
#' @usage ped.salmon
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references
#' Robledo D., Matika O., Hamilton A., and Houston R.D. 2018.
#' Genome-wide association and genomic selection for resistance
#' to amoebic gill disease in Atlantic salmon.
#' G3 Genes, Genomes, Genetics 8:1195-1203.
#'
#' @examples
#' ped.salmon |> head()
#'
#' @name ped.salmon
NULL
#' Phenotypic data for apple dataset
#'
#' Phenotypic data on 247 apple clones (\emph{i.e.}, genotypes) evaluated for
#' several fruit quality traits at two New Zealand sites, Motueka (MOT) and Hawkes Bay (HB).
#' Dataset obtained from supplementary material in Kumar \emph{et al.} (2015).
#'
#' @docType data
#'
#' @usage pheno.apple
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references
#' Kumar S., Molloy C., Muñoz P., Daetwyler H., Chagné D., and Volz R. 2015.
#' Genome-enabled estimates of additive and nonadditive genetic variances and prediction
#' of apple phenotypes across environments. G3 Genes, Genomes, Genetics 5:2711–2718.
#'
#' @examples
#' pheno.apple |> head()
#'
#' @name pheno.apple
NULL
#' Phenotypic data for loblolly pine dataset
#'
#' Deregressed estimated breeding values (DEBV) for the trait diameter at breast height (DBH)
#' at 6 years of age from trees grown at site Nassau.
#' The dataset contains a total of 861 genotypes of loblolly pine
#' (\emph{Pinus taeda} L.).
#' Dataset obtained from supplementary material in Resende \emph{et al.} (2012).
#'
#' @docType data
#'
#' @usage pheno.pine
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references
#' Resende, M.F.R., Munoz, P. Resende, M.D.V., Garrick, D.J., Fernando, R.L., Davis, J.M.,
#' Jokela, E.J., Martin, T.A., Peter, G.F., and Kirst, M. 2012. Accuracy of genomic
#' selection methods in a standard data set of loblolly pine (\emph{Pinus taeda} L.).
#' Genetics 190:1503-1510.
#'
#' @examples
#' pheno.pine |> head()
#'
#' @name pheno.pine
NULL
#' Phenotypic data for Atlantic salmon dataset
#'
#' Phenotypic data on 1,481 Atlantic salmon individuals.
#' All fish were phenotyped for mean gill score (mean of the left
#' gill and right gill scores) and amoebic load (qPCR values using
#' \emph{Neoparamoeba perurans} specific primers, amplified from one of the gills).
#' Dataset obtained from supplementary material in Robledo \emph{et al.} (2018).
#'
#' @docType data
#'
#' @usage pheno.salmon
#'
#' @format data.frame
#'
#' @keywords datasets
#'
#' @references
#' Robledo D., Matika O., Hamilton A., and Houston R.D. 2018.
#' Genome-wide association and genomic selection for resistance
#' to amoebic gill disease in Atlantic salmon.
#' G3 Genes, Genomes, Genetics 8:1195-1203.
#'
#' @examples
#' pheno.salmon |> head()
#'
#' @name pheno.salmon
NULL
#' Quality control filtering of molecular matrix M for downstream analyses
#'
#' Reads molecular data in the format 0, 1, 2 and performs some basic quality control
#' filters and simple imputation.
#' Matrix provided is of the full form (\eqn{n \times p}), with \eqn{n} individuals and \eqn{p} markers.
#' Individual and marker names are assigned to \code{rownames} and \code{colnames},
#' respectively. Filtering can be done with the some of the following options by
#' specifying thresholds for:
#' missing values on individuals, missing values on markers, minor allele frequency,
#' inbreeding Fis value (of markers), and observed heterozygosity (of markers).
#' String used for identifying missing values can be specified.
#' If requested, missing values will be imputed based on the mean of each SNP.
#'
#' @param M A matrix with SNP data of full form (\eqn{n \times p}), with \eqn{n} individuals and \eqn{p} markers
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as 0, 1, 2 (integer or numeric) (default = \code{NULL}).
#' @param base If \code{TRUE} matrix \eqn{\boldsymbol{M}} is considered as bi-allele SNP data format (character)
#' and the SNPs are recoded to numerical values before performing the quality control filters
#' (default = \code{FALSE}) (currently deprecated).
#' @param na.string A character that will be interpreted as \code{NA} values (default = \code{"NA"}).
#' @param map (Optional) A data frame with the map information with \eqn{p} rows (default = \code{NULL}).
#' @param marker A character indicating the name of the column in data frame \code{map} with the identification
#' of markers. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param chrom A character indicating the name of the column in data frame \code{map} with the identification
#' of chromosomes (default = \code{NULL}).
#' @param pos A character indicating the name of the column in data frame \code{map} with the identification
#' of marker positions (default = \code{NULL}).
#' @param ref A character indicating the name of the column in the map containing the reference allele for
#' recoding. If absent, then conversion will be based on the major allele (most frequent).
#' The marker information of a given individuals with two of the specified major alleles
#' in \code{ref} will be coded as 2 (default = \code{NULL}).
#' @param marker.callrate A numerical value between 0 and 1 used to remove SNPs with a rate
#' of missing values equal or larger than this value (default = 1, \emph{i.e.} no removing).
#' @param ind.callrate A numerical value between 0 and 1 used to remove individuals with a
#' rate of missing values equal or larger than this value (default = 1, \emph{i.e.} no removing).
#' @param maf A numerical value between 0 and 1 used to remove SNPs with a Minor Allele Frequency
#' (MAF) below this value (default = 0, \emph{i.e.} no removing).
#' @param heterozygosity A numeric value indicating the maximum value of accepted observed heterozygosity (Ho)
#' (default = 1, \emph{i.e.} no removing).
#' @param Fis A numeric value indicating the maximum value of accepted inbreeding (Fis) following
#' the equation \eqn{|1 - (Ho/He)|} (default = 1, \emph{i.e.} no removing).
#' @param impute If \code{TRUE} imputation of missing values is done using the mean of each SNP
#' (default = \code{FALSE}).
#' @param Mrecode If \code{TRUE} it provides the recoded \eqn{\boldsymbol{M}} matrix from the bi-allelic to numeric SNP
#' (default = \code{FALSE}) (currently deprecated).
#' @param plots If \code{TRUE} generates graphical output of the quality control based on the
#' original input matrix (default = \code{TRUE}).
#' @param digits Set up the number of digits used to round the output matrix (default = 2).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following elements:
#' \itemize{
#' \item \code{M.clean}: the cleaned \eqn{\boldsymbol{M}} matrix after the quality control filters have been applied.
#' \item \code{map}: if provided, a cleaned \code{map} data frame after the quality control filters have been applied.
#' \item \code{plot.missing.ind}: a plot of missing data per individual (original marker matrix).
#' \item \code{plot.missing.SNP}: a plot of missing data per SNP (original marker matrix).
#' \item \code{plot.heteroz}: a plot of observed heterozygocity per SNP (original marker matrix).
#' \item \code{plot.Fis}: a plot of Fis per SNP (original marker matrix).
#' \item \code{plot.maf}: a plot of the minor allele frequency (original marker matrix).
#' }
#'
#' @md
#' @details
#' \strong{Warning}: The arguments \code{base}, \code{ref}, and \code{Mrecode}
#' currently are deprecated and will
#' be removed on the next version of \code{ASRgenomics}.
#' Use function \link{snp.recode} to recode the matrix prior to using \code{qc.filtering}.
#'
#' The filtering process is carried out as expressed in the following simplified pseudo-code
#' that consists on a loop repeated twice:
#'
#' \strong{for i in 1 to 2}
#'
#' &nbsp; &nbsp; Filter markers based on call rate.
#'
#' &nbsp; &nbsp; Filter individuals based on call rate.
#'
#' &nbsp; &nbsp; Filter markers based on minor allele frequency.
#'
#' &nbsp; &nbsp; Filter markers based on observed heterozygosity.
#'
#' &nbsp; &nbsp; Filter markers based on inbreeding.
#'
#' \strong{end for}
#'
#' @export
#'
#' @examples
#' # Example: Pine dataset (coded as 0,1,2 with missing as -9).
#'
#' M.clean <- qc.filtering(
#'  M = geno.pine926,
#'  maf = 0.05,
#'  marker.callrate = 0.9, ind.callrate = 0.9,
#'  heterozygosity = 0.9, Fis = 0.6,
#'  na.string = "-9")
#' ls(M.clean)
#' M.clean$M.clean[1:5, 1:5]
#' dim(M.clean$M.clean)
#' head(M.clean$map)
#' M.clean$plot.maf
#' M.clean$plot.missing.ind
#' M.clean$plot.missing.SNP
#' M.clean$plot.heteroz
#' M.clean$plot.Fis
#'
#' \donttest{
#' # Example: Salmon dataset (coded as 0,1,2 with missing as NA).
#'
#' M.clean <- qc.filtering(
#'  M = geno.salmon,
#'  maf = 0.02,
#'  marker.callrate = 0.10, ind.callrate = 0.20,
#'  heterozygosity = 0.9, Fis = 0.4)
#' M.clean$M.clean[1:5, 1:5]
#' dim(M.clean$M.clean)
#' head(M.clean$map)
#' M.clean$plot.maf
#' M.clean$plot.missing.ind
#' M.clean$plot.missing.SNP
#' M.clean$plot.heteroz
#' M.clean$plot.Fis
#' }
#'


qc.filtering <- function(M = NULL, base = FALSE, na.string = NA,
                         map = NULL, marker = NULL, chrom = NULL, pos = NULL, ref = NULL,
                         marker.callrate = 1, ind.callrate = 1, maf = 0,
                         heterozygosity = 1, Fis = 1,
                         impute = FALSE, Mrecode = FALSE,
                         plots = TRUE, digits = 2, message = TRUE) {

  # Deprecation traps ---------------------------------------------------------------------------

  if (!is.null(ref) | base | Mrecode){
    stop("The recoding has been deprecated in \'qc.filtering()', please use \'snp.recode()' to perform this task.")
  }

  # Traps ---------------------------------------------------------------------------------------

  # Check if the class of M is matrix.
  if (is.null(M) || !inherits(M, "matrix"))
    stop("M should be a valid object of class matrix.")

  if (is.null(colnames(M)))
    stop("Marker names not assigned to columns of matrix M.")

  if (is.null(rownames(M)))
    stop("Individuals names not assigned to rows of matrix M.")

  # Other input checks.
  if (marker.callrate < 0 | marker.callrate > 1)
    stop("Specification of marker.callrate must be be between 0 and 1.")

  if (ind.callrate < 0 | ind.callrate > 1)
    stop("Specification of ind.callrate must be between 0 and 1.")

  if (maf < 0 | maf > 1)
    stop("Specification of maf must be between 0 and 1.")

  if (Fis < 0 | Fis > 1)
    stop("Specification of Fis must be between 0 and 1.")

  if (heterozygosity < 0 | heterozygosity > 1)
    stop("Specification of heterozygosity must be between 0 and 1.")

  # Check map if provided.
  if (!is.null(map)) {

    # Check map class.
    check.data_(data_ = "map", class_ = "data.frame")

    # Check map names.
    # Check mandatory variables in map.
    if(is.null(marker)){stop("\'marker' must be provided if \'map' is provided.")}

    # Check if they are present in the map.
    map.name.hit <- c(marker, chrom, pos) %in% names(map)
    if (!all(map.name.hit)){
      stop("Value provided to argument \'", c("marker", "chrom", "pos")[!map.name.hit],
           "' does not correspond to a variable in \'map'.")
    }

    # Match map and M.
    if (!identical(as.character(map[[marker]]), colnames(M))){
      stop("map[[marker]] and colnames(M) must be identical.")
    }
  }

  # # This is a slow test. Maybe not worth it. It is not necessary here.
  # if(!all(unique(c(M)) %in% c(0, 1, 2, na.string)) & message){
  #   message("Some of the values in M are not one of the following: ",
  #           paste0(c(0, 1, 2, na.string), collapse = ", "), ".")
  # }

  # Body ----------------------------------------------------------------------------------------

  # Initial info about the matrix.
  if (message) {
    message("Initial marker matrix M contains ", nrow(M),
            " individuals and ", ncol(M), " markers.")
  }

  # Replace na.string by NA.
  if (!is.na(na.string)) {
    if (na.string == "NA") { na.string <- NA }
  }
  if (!is.na(na.string)) {
    if (message){
     message('A total of ', sum(M %in% na.string),
             " values were identified as missing with the string ",
     na.string, " and were replaced by NA.")
    }
    M[M %in% na.string] <- NA
  }

  # Check if all are compliant.
  if (!all(M %in% c(0, 1, 2, NA))) {
    stop("Data must be in numeric format: 0, 1, 2 and NA.")
  }

  # Remove markers with no valid information.
  miss.all <- colSums(is.na(M)) == nrow(M)

  if (any(miss.all)) {

    # Apply the removal.
    M <- M[, !miss.all]

    # Report.
    if (message){
      message("A total of ", sum(miss.all), " markers were removed for only having missing data.")
    }
  } ; rm(miss.all)


  # Generating some plots ---------------------------------------------------------------------------------

  if (plots){

    # Missing of individuals.
    # missingInd_DF <- data.frame(Ind =rowMeans(is.na(M)))
    missingInd_DF <- data.frame(Ind = (100 - callrate(M = M, margin = "row"))/100)
    missingInd_plot <- ggplot(missingInd_DF, aes(x=Ind)) +
      geom_histogram(fill='#0072B2', bins=40) +
      theme_classic() +
      xlab("Missing data per Individual")+
      ylab("Count")
    rm(missingInd_DF)

    # Missing of markers.
    # missingSNP_DF <- data.frame(SNP=colMeans(is.na(M)))
    missingSNP_DF <- data.frame(SNP = (100 - callrate(M = M, margin = "col"))/100)
    missingSNP_plot <- ggplot(missingSNP_DF, aes(x=SNP)) +
      geom_histogram(fill='#0072B2', bins=40) +
      theme_classic() +
      xlab("Missing data per SNP") +
      ylab("Count")

    # Histogram of MAF.
    qDF <- data.frame(MAF = maf(M = M))
    maf_plot <- ggplot(qDF, aes(x=MAF)) +
      geom_histogram(fill = '#0072B2', bins = 40) +
      theme_classic() +
      xlab("Minor Allele Frequency (MAF)")+
      ylab("Count")

    # Histogram of heterozygotes.
    het_DF <- data.frame(het = heterozygosity(M = M)[, "ho"])
    het_plot <- ggplot(het_DF, aes(x = het)) +
      geom_histogram(fill = '#0072B2', bins = 40) +
      theme_classic() +
      xlab("Heterozygotes")+
      ylab("Count")

    # Histogram of Fis.
    fis_DF <- data.frame(fis = abs(Fis(M = M)))
    fis_plot <- ggplot(fis_DF, aes(x = fis)) +
      geom_histogram(fill = '#0072B2', bins = 40) +
      theme_classic() +
      xlab("Fis")+
      ylab("Count")

  } else {
    missingInd_plot <- NULL
    missingSNP_plot <- NULL
    maf_plot <- NULL
    het_plot <- NULL
    fis_plot <- NULL
  }

  # Filtering markers -------------------------------------------------------------------------------------

  # Filtering process - 2 rounds (initializing objects).
  cr_mk_out <- 0 ; cr_id_out <- 0
  maf_out <- 0 ; fis_out <- 0
  h_out <- 0
  cr_mk_filter <- TRUE ; cr_id_filter <- TRUE
  maf_filter <- TRUE ; fis_filter <- TRUE
  h_filter <- TRUE

  for (blank_ in 1:2) {

    # Filtering markers by CR.
    if (marker.callrate < 1){
      cr_mk_filter <- 1 - callrate(M = M, margin = "col")/100 <= marker.callrate
      cr_mk_out <- cr_mk_out + sum(!cr_mk_filter)
      M <- M[, cr_mk_filter, drop = FALSE]
      rm(cr_mk_filter)
    }

    # Filtering callrate of individuals.
    if (ind.callrate < 1){
      cr_id_filter <- 1 - callrate(M = M, margin = "row")/100 <= marker.callrate
      cr_id_filter <- rowMeans(is.na(M)) <= ind.callrate
      cr_id_out <- cr_id_out + sum(!(cr_id_filter))
      M <- M[cr_id_filter, , drop = FALSE]
      rm(cr_id_filter)
    }

    # Filtering markers by MAF.
    if (maf > 0){
      q <- maf(M = M)
      maf_filter <- q > maf - .Machine$double.eps
      maf_out <- maf_out + sum(!maf_filter, na.rm = TRUE)
      M <- M[, maf_filter, drop = FALSE]
      rm(maf_filter, q)
    }

    # Filtering markers by heterozygosity.
    if (heterozygosity < 1){

      # Get observed heterozygosity.
      h <- heterozygosity(M = M)[, "ho"]

      # Get incidence vector.
      h_filter <- h <= heterozygosity

      # Add current run to the sum.
      h_out <- h_out + sum(!h_filter, na.rm = TRUE)

      # Apply filter.
      M <- M[, h_filter, drop = FALSE]

      # Remove objects.
      rm(h_filter, h)
    }

    if (Fis < 1){

      # Get Fis.
      fis <- Fis(M = M, margin = "col")

      # Get incidence vector.
      fis_filter <- abs(fis) <= Fis

      # Add current run to the sum.
      fis_out <- fis_out + sum(!fis_filter, na.rm = TRUE)

      # Apply filter.
      M <- M[, fis_filter, drop = FALSE]

      # Remove objects.
      rm(fis_filter, fis)
    }
  }

  # Some intermediate reporting.
  if (message){
    message("A total of ", cr_mk_out, " markers were removed because ",
            "their proportion of missing values was equal or larger than ",
            marker.callrate, ".")

    message("A total of ", cr_id_out, " individuals were removed because ",
            "their proportion of missing values was equal or larger than ",
            ind.callrate, ".")

    message("A total of ", maf_out, " markers were removed because ",
                 "their MAF was smaller than ", maf, ".")

    message("A total of ", h_out, " markers were removed because ",
                 "their heterozygosity was larger than ", heterozygosity, ".")

    message("A total of ", fis_out, " markers were removed because ",
                 "their |F| was larger than ", Fis, ".")

    missing.SNP <- sum(is.na(M))
    prop.miss <- 100*missing.SNP/(ncol(M)*nrow(M))
    message("Final cleaned marker matrix M contains ", round(prop.miss,2),
            "% of missing SNPs.")

    message("Final cleaned marker matrix M contains ", nrow(M),
            " individuals and ", ncol(M), " markers.")
    }

  # Simple mean imputation.
  if (impute){
    missing.SNP <- sum(is.na(M))
    if (missing.SNP == 0 & isTRUE(message)) {
      message('No imputation was performed as there are no missing marker data.')
    } else {

      # Loop through markers and impute with mean.
      for(i in 1:ncol(M)){
        M[is.na(M[,i]), i] <- mean(M[,i], na.rm=TRUE)
      }

      # Polishing the dataset.
      M <- round(M, digits)
      if (isTRUE(message)){
        prop.miss <- 100*missing.SNP/(ncol(M)*nrow(M))
        message("A total of ", missing.SNP, " missing values were imputed, ",
                "corresponding to ", round(prop.miss,2), "% of the total number of SNPs.")
      }
    }
  }

  # Finalize ------------------------------------------------------------------------------------

  # Remove eventual cleaned markers from M in the map.
  if (!is.null(map)){

    # Get match index.
    matches <- na.omit(match(colnames(M), as.character(map[[marker]])))

    # Applies to map; This separation is done because of data.table.
    map <- map[matches,]
  }

  return(list(M.clean = M, map = map, plot.missing.ind = missingInd_plot,
         plot.missing.SNP = missingSNP_plot, plot.heteroz = het_plot, plot.Fis = fis_plot,
         plot.maf = maf_plot))

}
#' Performs a Principal Component Analysis (PCA) based on a molecular matrix M
#'
#' Generates a PCA and summary statistics from a given molecular matrix
#' for population structure. Matrix
#' provided is of full form (\eqn{n \times p}), with n individuals and p markers. Individual and
#' marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' SNP data is coded as 0, 1, 2 (integers or decimal numbers). Missing values are
#' not accepted and these need to be imputed (see function \code{qc.filtering()}
#' for implementing mean imputation). There is additional output such as plots and
#' other data frames
#' to be used on other downstream analyses (such as GWAS).
#'
#' It calls function \code{prcomp()} to generate the PCA and the
#' \code{factoextra} R package to extract and visualize results.
#' Methodology uses normalized allele frequencies as proposed by Patterson \emph{et al.} (2006).
#'
#' @param M A matrix with SNP data of full form (\eqn{n \times p}), with \eqn{n}
#' individuals and \eqn{p} markers (default = \code{NULL}).
#' @param label If \code{TRUE} then includes in output individuals names (default = \code{FALSE}).
#' @param ncp The number of PC dimensions to be shown in the screeplot, and to provide
#' in the output data frame (default = \code{10}).
#' @param groups Specifies a vector of class factor that will be used to define different
#' colors for individuals in the PCA plot. It must be presented in the same order as the individuals
#' in the molecular \eqn{\boldsymbol{M}} matrix (default = \code{NULL}).
#' @param ellipses If \code{TRUE}, ellipses will will be drawn around each of the define levels in
#' \code{groups} (default = \code{FALSE}).
#'
#' @return A list with the following four elements:
#' \itemize{
#' \item \code{eigenvalues}: a data frame with the eigenvalues and its variances associated with each dimension
#' including only the first \code{ncp} dimensions.
#' \item \code{pca.scores}: a data frame with scores (rotated observations on the new components) including
#' only the first \code{ncp} dimensions.
#' \item \code{plot.pca}: a scatterplot with the first two-dimensions (PC1 and PC2) and their scores.
#' \item \code{plot.scree}: a barchart with the percentage of variances explained by the \code{ncp} dimensions.
#' }
#'
#' @references
#' Patterson N., Price A.L., and Reich, D. 2006. Population structure and eigenanalysis.
#' PLoS Genet 2(12):e190. doi:10.1371/journal.pgen.0020190
#'
#' @export
#'
#' @examples
#' # Perform the PCA.
#' SNP_pca <- snp.pca(M = geno.apple, ncp = 10)
#' ls(SNP_pca)
#' SNP_pca$eigenvalues
#' head(SNP_pca$pca.scores)
#' SNP_pca$plot.pca
#' SNP_pca$plot.scree
#'
#' # PCA plot by family (17 groups).
#' grp <- as.factor(pheno.apple$Family)
#' SNP_pca_grp <- snp.pca(M = geno.apple, groups = grp, label = FALSE)
#' SNP_pca_grp$plot.pca
#'


snp.pca <- function(M = NULL, label = FALSE, ncp = 10,
                    groups = NULL, ellipses = FALSE){

  # Check if the class of M is matrix
  if (is.null(M) || !inherits(M, "matrix")) {
    stop("M should be a valid object of class matrix.")
  }
  if (is.null(colnames(M))){
    stop("Marker names not assigned to columns of matrix M.")
  }
  if (is.null(rownames(M))){
    stop("Individuals names not assigned to rows of matrix M.")
  }
  # Check if the are missing values
  if (any(is.na(M))){
    stop("M matrix contains some missing data, consider performing some imputation.")
  }
  if (ncp < 0 | ncp > nrow(M)) {
    stop("Value ncp must be positive and smaller than the number of rows in matrix M.")
  }

  ## PCA by Petterson et al. (2006)
  M.mean <- colMeans(M)/2
  M.scale <- sqrt(M.mean * (1 - M.mean))
  M.norm <- matrix(NA, nrow = nrow(M), ncol = ncol(M))
  for (i in 1:ncol(M)) {  # Done by SNP column
    M.norm[,i] <- (M[,i]/2 - M.mean[i]) / M.scale[i]
  }

  # Pass name of individuals along.
  rownames(M.norm) <- rownames(M)
  colnames(M.norm) <- colnames(M)

  # Generating the pca
  pca <- prcomp(var(t(M.norm)), scale.=FALSE)  # Original it takes a lont time

  # Percentage of variances explained by each principal component
  scree_plot <- fviz_eig(pca, addlabels=TRUE, ncp=ncp,
                                     barfill = "#0072B2",
                                     barcolor = "#0072B2",
                                     ggtheme = theme_classic())
  # Extract the eigenvalues/variances of the principal dimensions
  eig_var <- get_eig(pca)

  # Plot PCA
  if (isTRUE(label)) {
    if(is.null(groups)) {
      pca_plot <- fviz_pca_ind(pca, geom=c("point","text"),
                                           repel=TRUE,
                                           col.ind = "#0072B2",
                                           ggtheme = theme_classic())
    } else {
      pca_plot <- fviz_pca_ind(pca,
                                           geom = c("point","text"),
                                           repel = TRUE,
                                           col.ind = groups, # color by groups
                                           mean.point = FALSE,
                                           legend.title = "Groups",
                                           ggtheme = theme_classic())
    }
  }
  if (isFALSE(label)) {
    if(is.null(groups)) {
      pca_plot <- fviz_pca_ind(pca, geom="point",
                                           col.ind = "#0072B2",
                                           ggtheme = theme_classic())
    } else {
      pca_plot <- fviz_pca_ind(pca,
                                           geom = "point",
                                           col.ind = groups, # color by groups,
                                           mean.point = FALSE,
                                           legend.title = "Groups",
                                           ggtheme = theme_classic())
    }
  }

  # Process ellipses if requested.
  if (!is.null(groups) & ellipses){

    # TODO this can be more memory efficient (the data frame is being recreated).

    group.comps <- cbind.data.frame(pca$x[, 1:2], groups)

    # Get centroids.
    centroids <- aggregate(cbind(PC1, PC2) ~ groups, data = group.comps, FUN = mean)

    # Get ellipses.
    ellipses.data  <- do.call(
      rbind,
      lapply(unique(group.comps$groups), function(t) {

        data.frame(
          groups = as.character(t),
          ellipse(
            cov(group.comps[group.comps$groups == t, 1:2]),
                centre = as.matrix(centroids[t, 2:3]), level = 0.95),
          stringsAsFactors=FALSE)
      }
      )
    )

    # Add ellipses to plot.
    pca_plot <- pca_plot +
      geom_path(data = ellipses.data, linewidth = .5, inherit.aes = F,
                aes(x = PC1, y = PC2, color = groups))
  }


  # Scores (rotated X observations on the new components) for ncp components.
  scores <- pca$x[,c(1:ncp)]
  eigenvalues <- eig_var[c(1:ncp),]

  return(list(pca.scores=scores, eigenvalues=eigenvalues, plot.scree=scree_plot, plot.pca=pca_plot))

}
#' Reduces the number of redundant markers on a molecular matrix M by pruning
#'
#' For a given molecular dataset \eqn{\boldsymbol{M}} (in the format 0, 1 and 2)
#' it produces a reduced molecular matrix by eliminating "redundant"
#' markers using pruning techniques. This function finds and drops some of the
#' SNPs in high linkage disequilibrium (LD).
#'
#' Pruning is recommended as redundancies can affect
#' the quality of matrices used for downstream analyses.
#' The algorithm used is based on the Pearson's correlation between markers
#' as a \emph{proxy} for LD. In the event of a pairwise correlation higher
#' than the selected threshold markers will be eliminated as specified by: call rate,
#' minor allele frequency. In case of tie, one marker will be dropped at random.
#'
#' @param M A matrix with marker data of full form (\eqn{n \times p}), with \eqn{n} individuals
#' and \eqn{p} markers. Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as 0, 1, 2 (integer or numeric) (default = \code{NULL}).
#' @param map (Optional) A data frame with the map information with \eqn{p} rows.
#' If \code{NULL} a dummy map is generated considering a single chromosome and sequential positions
#' for markers. A \code{map} is mandatory if \code{by.chrom = TRUE}, where also option \code{chrom}
#' must also be non-null.
#' @param marker A character indicating the name of the column in data frame \code{map}
#' with the identification
#' of markers. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param chrom A character indicating the name of the column in data frame \code{map} with the identification
#' of chromosomes. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param pos A character indicating the name of the column in data frame \code{map} with the identification
#' of marker positions (default = \code{NULL}).
#' @param method A character indicating the method (or algorithm) to be used as reference for
#' identifying redundant markers.
#' The only method currently available is based on correlations (default = \code{"correlation"}).
#' @param criteria A character indicating the criteria to choose which marker to drop
#' from a detected redundant pair.
#' Options are: \code{"callrate"} (the marker with fewer missing values will be kept) and
#' \code{"maf"} (the marker with higher minor allele frequency will be kept) (default = \code{"callrate"}).
#' @param pruning.thr A threshold value to identify redundant markers with Pearson's correlation larger than the
#' value provided (default = \code{0.95}).
#' @param by.chrom If TRUE the pruning is performed independently by chromosome (default = \code{FALSE}).
#' @param window.n A numeric value with number of markers to consider in each
#' window to perform pruning (default = \code{50}).
#' @param overlap.n A numeric value with number of markers to overlap between consecutive windows
#' (default = \code{5}).
#' @param iterations An integer indicating the number of sequential times the pruning procedure
#' should be executed on remaining markers.
#' If no markers are dropped in a given iteration/run, the algorithm will stop (default = \code{10}).
#' @param seed An integer to be used as seed for reproducibility. In case the criteria has the
#' same values for a given pair of markers, one will be dropped at random (default = \code{NULL}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details Filtering markers (\link{qc.filtering}) is of high relevance before pruning.
#' Poor quality markers (\emph{e.g.}, monomorphic markers) may prevent correlations from being
#' calculated and may affect eliminations.
#'
#' @return
#' \itemize{
#'  \item{\code{Mpruned}}{: a matrix containing the pruned marker \emph{M} matrix.}
#'  \item{\code{map}}{: an data frame containing the pruned map.}
#' }
#'
#' @export
#'
#' @examples
#' # Read and filter genotypic data.
#' M.clean <- qc.filtering(
#'  M = geno.pine655,
#'  maf = 0.05,
#'  marker.callrate = 0.20, ind.callrate = 0.20,
#'  Fis = 1, heterozygosity = 0.98,
#'  na.string = "-9",
#'  plots = FALSE)$M.clean
#'
#' # Prune correlations > 0.9.
#' Mpr <- snp.pruning(
#'  M = M.clean, pruning.thr = 0.90,
#'  by.chrom = FALSE, window.n = 40, overlap.n = 10)
#' head(Mpr$map)
#' Mpr$Mpruned[1:5, 1:5]
#'

snp.pruning <- function(
    M = NULL, map = NULL, marker = NULL, chrom = NULL, pos = NULL,
    method = c('correlation'), criteria = c("callrate", "maf"),
    pruning.thr = 0.95, by.chrom = FALSE, window.n = 50, overlap.n = 5,
    iterations = 10,
    # n.cores = 1,
    seed = NULL, message = TRUE) {


  # Traps ---------------------------------------------------------------------------------------

  # Check M class.
  check.data_(data_ = "M", class_ = "matrix")

  # Get maf and callrate.
  # Logic: as we are using correlations, viable maf and callrate are essential.
  maf <- maf(M = M)
  callrate <- callrate(M = M, margin = "col")

  # Check callrate.
  if (any(callrate == 0)){
    stop("There are markers will all samples missing. Please use qc.filtering() before pruning.")
  }

  # Check maf.
  if (any(maf == 0)){
    stop("There are markers with minor allele frequency equal to 0. Please use qc.filtering() before pruning.")
  }

  # Create map if not provided.
  if (is.null(map)) {

    map <- dummy.map_(colnames(M))
    marker <- "marker" ; chrom <- "chrom" ; pos <- "pos"

  } else {
    # Check map class.
    check.data_(data_ = "map", class_ = "data.frame")

    # Check map names.
    # Check mandatory variables in map.
    if(is.null(marker)){stop("The \'marker' option must be specified if \'map' is provided.")}
    if(is.null(chrom)){stop("The \'chrom' option must be specified if \'map' is provided.")}

    # Check if they are present in the map.
    map.name.hit <- c(marker, chrom, pos) %in% names(map)
    if (!all(map.name.hit)){
      stop("Value provided to argument \'", c("marker", "chrom", "pos")[!map.name.hit], "' does not correspond to a variable in
           data frame \'map'.")
    }

    # Match map and M.
    if (!identical(as.character(map[[marker]]), colnames(M))){
      stop("map[[marker]] and colnames(M) must be identical. Please check input.")
    }
  }

  # Check method.
  method <- match.arg(method)

  # Check criteria.
  criteria <- match.arg(criteria)

  # Check threshold.
  if (pruning.thr <= 0  | pruning.thr > 1){
    stop("The condition for pruning.thr is between 0 and 1.")
  }

  # Check by.chrom.
  check.logical_(arg_ = "by.chrom")

  # Check window.
  if(window.n <= 1){
    stop("The \'window.n' argument should have an integer larger than 1.")
  }

  # Check overlap.
  if(overlap.n <= 0){
    stop("The \'overlap.n' argument should have an integer larger or eqaul than 0.")
  }

  if(overlap.n >= window.n){
    stop("The \'overlap.n' argument should be lower than the \'window.n' argument.")
  }

  # Check iterations.
  if(iterations <= 0){
    stop("The \'iterations' argument should have an positive integer.")
  }

  # if(n.cores <= 0){
  #   stop("The \'n.cores' argument should have an integer > 0.")
  # }

  # Check message.
  check.logical_(arg_ = "message")

  # Body ----------------------------------------------------------------------------------------

  # Setting seed.
  if (!is.null(seed)) { set.seed(seed = seed) }

  # Identify tiebraking criteria.
  # Logic: the criteria to select markers are maf or callrate.
  if (criteria == "maf"){
    map$criteria <- maf
  }

  if (criteria == "callrate"){
    map$criteria <- callrate
  }

  # Selection dummy.
  map$sel <- 1

  # Ordering by maf if requested.
  # TODO check how this affects the code below.
  # if (maf.order) {
  #   map <- map[order(map$maf, decreasing = FALSE),]
  #   M <- M[, map[[marker]]]
  # }

  # Collect garbage.
  rm(maf, callrate)

  # Marker drop function ------------------------------------------------------------------------

  # Call function that drops markers based on the correlation.
  marker.drop <- function(curr.set.index = NULL){

    init.set.pos <- sets[curr.set.index]  # Position on map and M.

    # Selecting section of M based on set and overlap.
    if (n.sets == 0){
      window.M <- cur.M
    } else if (curr.set.index == n.sets) {
      window.M <- cur.M[, sets[curr.set.index]:ncol(cur.M)]
    } else {
      window.M <- cur.M[, sets[curr.set.index]:(sets[curr.set.index + 1] + overlap.n - 1)]
    }

    # Generating Corr matrix in sparse (no diagonals).
    C.sparse <- suppressWarnings(cor(window.M, use = 'pairwise.complete.obs'))

    # Replace NAs with 0 to avoid problems with full2sparse.
    # Cause: if a correlation cannot be calculated for some reason,
    # we set it up as 0 and do not remove the marker.
    # This is a conservative approach as we do not have enough info about the marker.
    is.na(C.sparse) <- 0

    # Transform to sparse.
    C.sparse <- as.data.table(full2sparse(C.sparse))
    C.sparse[, Value := abs(Value)]
    C.sparse <- C.sparse[Row != Col,]

    # Order so we check from largest to smaller correlation.
    setorder(C.sparse, -Value)

    # Initiate indices to remove.
    rm.pos <- c()

    while (C.sparse$Value[1] >= pruning.thr) {

      # Identify current Row and Col corrected for set position.
      row.pos <- C.sparse[1, Row] + init.set.pos - 1
      col.pos <- C.sparse[1, Col] + init.set.pos - 1

      # Selecting which marker to keep based on the number of missing.
      # Row and Col equal then random.
      if (cur.map$criteria[row.pos] == cur.map$criteria[col.pos]) {

        # Get a random TRUE or FALSE to select the marker.
        if (sample(x = c(TRUE, FALSE), size = 1)) {
          # Drop marker on col.

          # Update C.sparse.
          C.sparse <- C.sparse[Col !=  Col[1] & Row != Col[1], ]

          # Append position to remove.
          rm.pos <- append(rm.pos, col.pos)

        } else {
          # Drop marker on row.

          # Update C.sparse.
          C.sparse <- C.sparse[Col !=  Row[1] & Row != Row[1], ]

          # Append position to remove.
          rm.pos <- append(rm.pos, row.pos)

        }
      }

      # Row better than Col.
      else if (cur.map$criteria[row.pos] > cur.map$criteria[col.pos]) {

        # Update C.sparse.
        C.sparse <- C.sparse[Col !=  Col[1] & Row != Col[1], ]

        # Append position to remove.
        rm.pos <- append(rm.pos, col.pos)
      }

      # Col better than Row.
      else {

        # Update C.sparse.
        C.sparse <- C.sparse[Col !=  Row[1] & Row != Row[1], ]

        # Append position to remove.
        rm.pos <- append(rm.pos, row.pos)
      }
    }

    return(rm.pos)
  }

  # Collect info for summary --------------------------------------------------------------------

  if (message){

    # Total number of markers.
    original.n.markers <- ncol(M)

    # Number of markers by chromosome.
    if (by.chrom){
      original.n.markers.chrom <- table(map$chrom)
    }
  }

  # Iterate through data ------------------------------------------------------------------------

  if (message){
    message(blue("\nInitiating pruning procedure."))
    message("Initial marker matrix M contains ", nrow(M),
            " individuals and ", ncol(M), " markers.")
  }

  # If by.chrom is requested, the range is obtained.
  if (by.chrom){

    # Get chromosome range.
    chrom.range <- unique(map[[chrom]])

    if (message){
      message("Requesting pruning by chromosome.")
    }

  } else {

    # Create dummy chromosome range.
    chrom.range <- 1

    if (message){
      message("Requesting pruning without chromosome indexing.")
    }
  }

  # Number of times to iterate in each chromosome.
  iter.range <- 1:iterations

  # Loop across chromosomes.
  for (cur.chrom in chrom.range){

    if (length(chrom.range) > 1 & message){
      message(paste0("Chromosome: ", cur.chrom))
    }

    # Split datasets if by.chrom was requested.
    # Get split index.
    if(by.chrom){

      split.index <- map[, chrom] == cur.chrom

      # Get data to for current chromosome.
      cur.map <- map[split.index, ]
      cur.M <- M[, split.index]

      # Save other chromosomes.
      map <- map[!split.index, ]
      M <- M[, !split.index, drop = FALSE]

    } else {

      # Collect map and M for calculations.
      # Logic: this is required if map is passed but by.chrom is FALSE.
      # Original map and M have to be NULL because they are bound later.
      cur.map <- map
      map <- NULL
      cur.M <- M
      M <- NULL
    }

    # Pre-set objects needed in loop.
    drop.pos <- NULL

    # This must be a loop because it is conditional.
    # Logic: this section is conditional to the previous one, so, no parallelization in R.
    for (iter in iter.range) {

      if (message){
        message("  Iteration: ", iter)
      }

      # Tag markers to eliminate.
      cur.map$sel[unlist(drop.pos)] <- 0

      # Stop criteria. If there was no drop on the last run. Break.
      if (iter > 1 & all(cur.map$sel == 1)) { break }

      # Eliminate markers.
      cur.M <- cur.M[, cur.map$sel == 1]
      cur.map <- cur.map[cur.map$sel == 1, ]

      # Defining step size (based on current data).
      step <- window.n - overlap.n

      # Get sets based on step size.
      sets <- seq(1, ncol(cur.M), step)
      n.sets <- length(sets) - 1

      # Get range of sets to loop across.
      sets.range <- 1:(n.sets)

      # Looping across all sets.
      # if (n.cores > 1){
      #   drop.pos <- mclapply(X = sets.range, mc.cores = n.cores,
      #                                  FUN = marker.drop, mc.set.seed = seed)
      # } else {
      drop.pos <- lapply(X = sets.range, FUN = marker.drop)
      # }
    }

    # Bind chromosomes.
    map <- rbind(map, cur.map)
    M <- cbind(M, cur.M)

  }

  # Summary -------------------------------------------------------------------------------------

  if (message){

    message("\nFinal pruned marker matrix M contains ", nrow(M),
            " individuals and ", ncol(M), " markers.")
    #message("A total of ", ncol(M), " markers were kept after pruning.")
    message("A total of ", original.n.markers - ncol(M), " markers were pruned.")

    if (by.chrom){
      # Number of markers by chromosome.
      message(paste0("A total of ", table(map$chrom), " markers were kept in chromosome ",
                     names(table(map$chrom)), ".", collapse = "\n"))
      message(paste0("A total of ", original.n.markers.chrom - table(map$chrom),
                     " markers were pruned from chromosome ", names(table(map$chrom)),
                     ".", collapse = "\n"))
    }

    # maf and call rate report.
    message("Range of minor allele frequency after pruning: ",
            paste0(round(range(maf(M = M)), 2), collapse = " ~ "))
    message("Range of marker call rate after pruning: ",
            paste0(round(range(callrate(M = M, margin = "col")), 2), collapse = " ~ "))
    message("Range of individual call rate after pruning: ",
            paste0(round(range(callrate(M = M, margin = "row")), 2), collapse = " ~ "))

  }

  # Finilize ------------------------------------------------------------------------------------

  map <- map[, !names(map) %in% c("criteria", "sel")]

  # Return pruned map and molecular matrix.
  return(list(map = map, Mpruned = M))

}
#' Recodes the molecular matrix M for downstream analyses
#'
#' Reads molecular data in format of bi-allelic nucleotide bases (AA,
#' AG, GG, CC, etc.) and recodes them as 0, 1, 2 and \code{NA} to be used in other
#' downstream analyses.
#'
#' @param M A character matrix with SNP data of full form (\eqn{n \times p}),
#' with \eqn{n} individuals and \eqn{p} markers
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as AA, AG, GG, CC, etc (default = \code{NULL}).
#' @param recoding A character indicating the recoding option to be performed.
#' Currently, only the nucleotide bases (AA, AG, ...) to allele count is available (\code{"ATGCto012"})
#' (default = \code{"ATGCto012"}).
#' @param map (Optional) A data frame with the map information with \eqn{p} rows.
#' If \code{NULL} a dummy map is generated considering a single chromosome and sequential
#' positions for markers and includes reference allele and alternative allele (default = \code{NULL}).
#' @param marker A character indicating the name of the column in data frame \code{map} with the identification
#' of markers. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param ref A character indicating the name of the column in the map containing the reference allele for
#' recoding. If absent, then conversion will be based on the major allele (most frequent).
#' The marker information of a given individual with two of the specified major alleles
#' in \code{ref} will be coded as 2. This is mandatory if \code{map} is provided (default = \code{NULL}).
#' @param alt A character indicating the name of the column in the map containing the alternative allele for
#' recoding. If absent, then it will be inferred from the data. The marker information of a given individual
#' with two of the specified alleles in \code{alt} will be coded as 0 (default = \code{NULL}).
#' @param na.string A character that is interpreted as missing values (default = \code{"NA"}).
#' @param rename.markers If \code{TRUE} marker names (as provided in \strong{M}) will be expanded
#' to store the reference and alternative alleles. For example, from AX-88234566 to AX-88234566_C_A.
#' In the event of unidentified alleles, 0 will be used (default = \code{TRUE}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @return A list with the following two elements:
#' \itemize{
#' \item \code{Mrecode}: the molecular matrix \eqn{\boldsymbol{M}} recoded to
#' 0, 1, 2 and \code{NA}.
#' \item \code{mapr}: the data frame with the map information including
#' reference allele and alternative allele.
#' }
#'
#' @export
#'
#' @examples
#' # Create bi-allelic base data set.
#' Mnb <- matrix(c(
#'   "A-",  NA, "GG",   "CC",   "AT",   "CC",   "AA",   "AA",
#'   "AAA", NA, "GG",   "AC",   "AT",   "CG",   "AA",   "AT",
#'   "AA",  NA, "GG",   "CC",   "AA",   "CG",   "AA",   "AA",
#'   "AA",  NA, "GG",   "AA",   "AA",    NA,    "AA",   "AA",
#'   "AT",  NA, "GG",   "AA",   "TT",   "CC",   "AT",   "TT",
#'   "AA",  NA,   NA,   "CC",    NA,    "GG",   "AA",   "AA",
#'   "AA",  NA,   NA,   "CC",   "TT",   "CC",   "AA",   "AT",
#'   "TT",  NA, "GG",   "AA",   "AA",   "CC",   "AA",   "AA"),
#'   ncol = 8, byrow = TRUE, dimnames = list(paste0("ind", 1:8),
#'                                        paste0("m", 1:8)))
#' Mnb
#'
#' # Recode without map (but map is created).
#' Mr <- snp.recode(M = Mnb, na.string = NA)
#' Mr$Mrecode
#' Mr$map
#'
#' # Create map.
#' mapnb <- data.frame(
#'  marker = paste0("m", 1:8),
#'  reference = c("A", "T", "G", "C", "T", "C", "A", "T"),
#'  alternative = c("T", "G", "T", "A", "A", "G", "T", "A")
#'  )
#'  mapnb
#'
#' # Recode with map without alternative allele.
#' Mr <- snp.recode(M = Mnb, map = mapnb, marker = "marker", ref = "reference",
#'            na.string = NA, rename.markers = TRUE)
#' Mr$Mrecode
#' Mr$map
#'
#' # Notice that the alternative allele is in the map as a regular variable,
#' # but in the names it is inferred from data (which might be 0 (missing)).
#'
#' # Recode with map with alternative allele.
#' Mr <- snp.recode(M = Mnb, map = mapnb, marker = "marker",
#'  ref = "reference", alt = "alternative",
#'  na.string = NA, rename.markers = TRUE)
#' Mr$Mrecode
#' Mr$map # Now the alternative is also on the names.
#'
#' # We can also recode without renaming the markers.
#' Mr <- snp.recode(M = Mnb, map = mapnb, marker = "marker", ref = "reference",
#'            na.string = NA, rename.markers = FALSE)
#' Mr$Mrecode
#' Mr$map # Now the alternative is also on the names.
#'

snp.recode <- function(M = NULL, map = NULL, marker = NULL, ref = NULL, alt = NULL,
                       recoding = c("ATGCto012"),
                       na.string = NA, rename.markers = TRUE,
                       message = TRUE){

  # Traps ---------------------------------------------------------------------------------------

  # Check recoding.
  # recoding is just a placeholder for now.
  recoding <- match.arg(recoding)

  # Check class of M.
  check.data_(data_ = "M", class_ = "matrix")

  # Check if class of values is character.
  check.data.mode_(data_ = "M", mode_ = "character")

  # Create map if not provided (same as in pruning).
  if (!is.null(map)) {

    # Check map class.
    check.data_(data_ = "map", class_ = "data.frame")

    # Check map names. Check mandatory variables in map.
    if(is.null(marker)){stop("The \'marker' option must be specified if \'map' is provided.")}
    if(is.null(ref)){stop("The \'ref' option must be specified if \'map' is provided.")}

    # Check if they are present in the map.
    map.name.hit <- c(marker, ref, alt) %in% names(map)
    if (!all(map.name.hit)){
      stop("Value provided to argument \'", c("marker", "ref", "alt")[!map.name.hit],
           "' does not correspond to a variable in \'map'.")
    }

    # Match map and M.
    if (!identical(as.character(map[[marker]]), colnames(M))){
      stop("map[[marker]] and colnames(M) must be identical.")
    }
  }

  # Check if all values are composed of 2 letters.
  if (any(nchar(M) != 2, na.rm = TRUE)) {
    warning("Marker(s) not compliant with bi-allelic coding: ",
            paste0(colnames(M)[ceiling(which(nchar(M) != 2) / nrow(M))], collapse = ", "),
            ".\n  The respective datapoints have been replaced with NA.")
    M[nchar(M) != 2] <- NA
  }

  special.char <- apply(X = M, MARGIN = 2, FUN = function(col){
    any(grepl(pattern = "[[:punct:]]", x = col))
  })

  # Removing eventual special characters (e.g. A-).
  if (any(special.char)){
    warning("Special characters identified in marker(s): ",
            paste0(names(special.char)[special.char], collapse = ", "),
            ".\n  The respective datapoints have been replaced with NA.")
    M[grepl(pattern = "[[:punct:]]", x = M)] <- NA
  }

  # Body ----------------------------------------------------------------------------------------

  # Replace na.string by NA.
  if (!is.na(na.string)) {
    if (na.string == "NA") { na.string <- NA }
  }
  if (!is.na(na.string)) {
    if (message){
      message('A total of ', sum(M %in% na.string),
              " values were identified as missing with the string ",
              na.string, " and were replaced by NA.")
    }
    M[M %in% na.string] <- NA
  }

  # Function to get sorted states of a marker.
  get.states_ <- function(m = NULL){
    sort( # TODO sorting is not really required now.
      unique(
        unlist(
          strsplit(x = m, split = ""))))
  }

  # Identify the states.
  states <- apply(X = M, MARGIN = 2, FUN = get.states_, simplify = FALSE)

  # Initiate main frame.
  reference.frame <- data.table()

  # Get number of states.
  reference.frame[, n.states := sapply(states, length) ]

  # Get markers names.
  reference.frame[, marker := colnames(M)]

  # Check for more than 2 states.
  if (any(reference.frame$n.states > 2)) {
    stop("Markers with more than two allelic states: ",
         paste0(colnames(M)[which(reference.frame$n.states > 2)], collapse = ", "),".")
  }

  # Add states to the frame.
  reference.frame[, c("state1", "state2"):=
                    list(sapply(states, function(m) m[1]),
                         sapply(states, function(m) m[2]))]

  # Add reference state to reference frame.
  if(is.null(ref)) {

    # Logic: the following code calculates the MAF based on the bases.
    # It pastes all together, and the separate all and counts the composing letters.
    # It this is too slow we might have to change a bit.
    # The ifelse trick is required if a marker only has missing data.
    reference.frame[, ref :=
                      apply(X = M, MARGIN = 2, FUN = function(m){
                        cur.ref <- names(
                          which.max(
                            table(
                              strsplit(
                                paste0(na.omit(m), collapse = ""),
                                "")
                              )
                            )
                          )
                        return(
                          ifelse(test = is.null(cur.ref),
                                 yes = NA,
                                 no = cur.ref))
                      })
    ]

  } else {
    reference.frame[, ref := ..map[[ref]]]

    if (!is.null(alt)){
      reference.frame[, alt := ..map[[alt]]]
    }
  }

  # Identify the alternative state (based on the reference state).
  if (is.null(alt)){
    reference.frame[, alt := ifelse(test = {state2 == ref}, yes = state1, no = state2)]
  }

  # Find wrong coding in reference allele.
  # Find wrong references.
  wrong.code <- reference.frame[state1 != ref & state2 != ref, marker]

  if (length(wrong.code) > 0) {
    stop("The provided reference (\'ref') missmatches the allele codings in: ", wrong.code, ".")
  }

  # Find wrong coding in alternative allele.
  if (!is.null(alt)){

    # Find wrong references.
    wrong.code <- reference.frame[state1 != alt & state2 != alt, marker]
    if (length(wrong.code) > 0) {
      stop("The provided reference (\'alt') missmatches the allele codings in: ", wrong.code, ".")
    }
  } ; rm(wrong.code)

  # TODO try replacing on a new matrix.
  # Create combinations for comparisons.
  # NAs are checked in the alternative state for the first 3, and on the reference for the 4.
  reference.frame[!is.na(alt), code0 := paste0(alt, alt)]
  reference.frame[!is.na(alt), code1A := paste0(ref, alt)]
  reference.frame[!is.na(alt), code1B := paste0(alt, ref)]
  reference.frame[!is.na(ref), code2 := paste0(ref, ref)]

  # Rename markers if requested.
  if (rename.markers){
    # reference.frame[, marker := paste0(marker, "_", ref, "_", alt)]
    reference.frame[, marker :=
                      paste0(marker, "_",
                             replace(x = ref, list = is.na(ref), values = "0"), "_",
                             replace(x = alt, list = is.na(alt), values = "0"))]

  }

  # Replace letters with numbers.
  M <- sapply(1:ncol(M), FUN = function(index){
    m <- M[, index]
    tmp.ref <- reference.frame[index,]
    m[m %in% na.omit(tmp.ref[["code0"]])] <- 0
    m[m %in% na.omit(tmp.ref[["code1A"]])] <- 1
    m[m %in% na.omit(tmp.ref[["code1B"]])] <- 1
    m[m %in% na.omit(tmp.ref[["code2"]])] <- 2
    return(m)
  })

  # Reassign names to M.
  colnames(M) <- reference.frame[["marker"]]

  # Transform to numeric.
  mode(M) <- "numeric"

  # Finalize ------------------------------------------------------------------------------------

  # Report.
  if (message) {
    message("Matrix M was recoded from bi-allelic nucleotide bases to numeric.")
  }

  # Prepare ref to export.
  if (is.null(map)){
    map <- dummy.map_(marker.id = reference.frame[["marker"]], message = FALSE)

    # Add reference and alternative alleles to map.
    map$ref <- reference.frame$ref
    map$alt <- reference.frame$alt

  } else {

    # If map is not NULL and rename is requested. Collect names from reference frame.
    if(rename.markers){
      map[[marker]] <- reference.frame[["marker"]]
    }
  }


  # Return the output list.
  return(list(Mrecode = M, map = map))

}
#' Generates a full matrix form from a sparse form matrix
#'
#' Modifies the input sparse form matrix into its full form.
#' The sparse form has three columns per line, corresponding to the set:
#' \code{Row, Col, Value}, and is defined by a lower triangle row-wise
#' of the full matrix and is sorted as columns within row.
#' Individual names should be assigned as attributes: \code{attr(K, "rowNames")}
#' and \code{attr(K, "colNames")}. If these are not provided they are considered
#' as 1 to \eqn{n}.
#'
#' Based on a function from ASReml-R 3 library by Butler \emph{et al.} (2009).
#'
#' @param K A square matrix in sparse form (default = \code{NULL}).
#'
#' @return A full square matrix where individual names are assigned to
#' \code{rownames} and \code{colnames}.
#' If attribute \code{INVERSE} is found this is also passed to the full matrix.
#'
#' @references
#' Butler, D.G., Cullis, B.R., Gilmour, A.R., and Gogel, B.J. 2009.
#' ASReml-R reference manual. Version 3. The Department of Primary
#' Industries and Fisheries (DPI&F).
#'
#' @export
#'
#' @examples
#' # Get G matrix.
#' Gsp <- G.matrix(M = geno.apple, method = "VanRaden", sparseform = TRUE)$G.sparse
#' head(Gsp)
#' head(attr(Gsp, "rowNames"))
#'
#' # Transform into full matrix.
#' G <- sparse2full(K = Gsp)
#' G[1:5, 1:5]
#'

sparse2full <- function(K = NULL) {

  if (is.null(K) || !inherits(K, "matrix")) {
    stop('K should be a valid object of class matrix.')
  }

  # Collect inverse attribute if any.
  INVERSE <- attr(K, "INVERSE")

  # Collect dimnames to apply on new matrix.
  rownames <- attr(K, "rowNames")

  # Collect number of rows.
  nrow <- max(K[, 1])

  # Collect number of columns.
  ncol <- max(K[, 2])

  # Get dummy rownames if not provided in attributes.
  if (is.null(rownames)) {
    rownames <- as.character(1:nrow)
  }

  # Assign relationships to relative positions.
  K.full <- rep(0, nrow * ncol)
  K.full[(K[, 2] - 1) * nrow + K[, 1]] <- K[, 3]
  K.full[(K[, 1] - 1) * nrow + K[, 2]] <- K[, 3]

  # Reshape to matrix.
  K.full <- matrix(K.full, nrow = nrow, ncol = ncol, byrow = FALSE)

  # Assign row and colnames.
  attr(K.full, "colNames") <-
    attr(K.full, "rowNames") <-
      colnames(K.full) <-
        rownames(K.full) <-
          rownames

  if (!is.null(INVERSE)) {attr(K.full, "INVERSE") <- INVERSE}

  return(K.full)

}
#' Generates a molecular matrix M for hypothetical crosses based on the
#' genomic information of the parents
#'
#' This function generates (or imputes) a molecular matrix for offspring
#' from hypothetical crosses based on the genomic information from the parents.
#' This is a common procedure in species such as maize, where only the parents
#' (inbred lines) are genotyped, and this information is used to generate/impute
#' the genotypic data of each of the hybrid offspring.
#' This function can be also used for bulked DNA analyses, in order to obtain an
#' bulked molecular matrix for full-sib individuals were only parents are genotyped.
#'
#'
#' @param M A matrix with marker data of full form (\eqn{n \times p}), with \eqn{n} individuals
#' (mothers and fathers) and \eqn{p} markers.
#' Individual and marker names are assigned to \code{rownames} and \code{colnames}, respectively.
#' Data in matrix is coded as 0, 1, 2 (integer or numeric) (default = \code{NULL}).
#' @param ped A data frame with three columns containing only the pedigree of the hypothetical offspring.
#' (not pedigree of parents)
#' It should include the three columns for individual, mother and father (default = \code{NULL}).
#' @param indiv A character indicating the column in \code{ped} data frame containing the identification
#' of the offspring (default = \code{NULL}).
#' @param mother A character indicating the column in \code{ped} data frame containing the identification
#' of the mother (default = \code{NULL}).
#' @param father A character indicating the column in \code{ped} data frame containing the identification
#' of the father (default = \code{NULL}).
#' @param heterozygote.action Indicates the action to take when heterozygotes are found in a marker.
#' Options are: \code{"useNA"}, \code{"exact"}, \code{"fail"}, and \code{"expected"}.
#' See details for more information (default = \code{"useNA"})
#' @param na.action  Indicates the action to take when missing values are found in a marker.
#' Options are: \code{"useNA"} and \code{"expected"}.
#' See details for more information (default = \code{"useNA"}).
#' @param message If \code{TRUE} diagnostic messages are printed on screen (default = \code{TRUE}).
#'
#' @details
#' For double-haploids, almost the totality of the markers (except for genotyping errors)
#' will be homozygotic reads. But in many other cases (including recombinant inbred lines)
#' there will be a proportion of heterozygotic reads. In these case, it is
#' very difficult to infer (impute) the exact genotype of a given offspring individual.
#' For example, if parents are 0 (AA) and 1 (AC) then offsprings will differ given this
#' Mendelian sampling. However, different strategies exist to determine the
#' expected value for that specific cross (if required), which are detailed below
#' using the option \code{heterozygote.action}.
#'
#' \itemize{
#'
#'  \item{If \code{heterozygote.action =  "useNA"},
#' the generated offspring will have, for the heterozygote read, an \code{NA},
#' and no markers are removed.
#' Hence, no attempt will be done to impute/estimate this value.}
#'
#'  \item{If \code{heterozygote.action = "exact"},
#' any marker containing one or more heterozygote reads will be removed.
#' Hence, inconsistent markers are fully removed from the \eqn{\boldsymbol{M}} matrix.}
#'
#'  \item{If \code{heterozygote.action = "fail"},
#' function stops and informs of the presence of heterozygote reads.}
#'
#'  \item{If \code{heterozygote.action = "expected"},
#' then an algorithm is implemented, on the heterozygote read to determine its
#' expected value. For example, if parents are 0 and 1, then the expected value
#' (with equal probability) is 0.5. For a cross between two heterozygotes,
#' the expected value is: \eqn{0(1/4) + 1(1/2) + 2(1/4) = 1}. And for a cross
#' between 1 and 2, the expected value is: \eqn{1(1/2) + 2(1/2) = 1.5}}
#'
#' }
#'
#' Missing value require special treatment, and an imputation strategy is detailed
#' below as indicated using the option \code{na.action}.
#'
#' \itemize{
#'
#'  \item{If \code{na.action = "useNA"},
#'  if at least one of the parental reads is missing values for a given marker then it will be assigned
#'  as missing for the hypothetical cross. Hence, no attempt will be done to impute/estimate this value.
#   }
#'
#'  \item{If \code{na.action = "expected"},
#'  then an algorithm is implemented that will impute the expected read of the cross
#'  if the genotype of \strong{one of the parents is missing} (\emph{e.g.}, cross between 0 and NA).
#'  Calculations are based on parental allelic frequencies \eqn{p} and \eqn{q} for the given marker.
#'  The expressions for expected values are detailed below.
#'  }
#'
#' \itemize{
#'
#' \item{If the genotype of the non-missing parent read is 0.
#'
#' \eqn{q^2} (probability that the missing parent is 0) x 0 (expected value of the offspring from a 0 x 0 cross: \eqn{0(1/1)}) +
#'
#' \eqn{2pq} (probability that the missing parent is 1) x 0.5 (expected offspring from a 0 x 1 cross: \eqn{0(1/2) + 1(1/2)}) +
#'
#' \eqn{q^2} (probability that the missing parent is 2) x 1 (expected offspring from a 0 x 2 cross: \eqn{1(1/1)})}
#'
#' \item{If the genotype of the non-missing parent read is 1.
#'
#' \eqn{q^2} (probability that the missing parent is 0) x 0.5 (offspring: \eqn{0(1/2) + 1(1/2)}) +
#'
#' \eqn{2pq} (probability that the missing parent is 1) x 1 (offspring: \eqn{0(1/4) + 1(1/2) + 2(1/4)}) +
#'
#' \eqn{q^2} (probability that the missing parent is 2) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)})}
#'
#' \item{If the genotype of the non-missing parent read is 2.
#'
#' \eqn{q^2} (probability that the missing parent is 0) x 1 (offspring: \eqn{1(1/1)}) +
#'
#' \eqn{2pq} (probability that the missing parent is 1) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)}) +
#'
#' \eqn{q^2} (probability that the missing parent is 2) x 2 (offspring: \eqn{2(1/1)})}
#' }
#' }
#' }
#'
#' Similarly, the calculation of the expected read of a cross when \strong{both parents are missing} is
#' also based on population allelic frequencies for the given marker.
#' The expressions for expected values are detailed below.
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{q^2 \times q^2} (probability that both parents are 0) x 0 (expected value of the offspring from a 0 x 0 cross: 0(1/1)) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2 \times (q^2 \times 2pq)} (probability that the first parent is 0 and the second is 1; this requires
#' the multiplication by 2 because it is also possible that the first parent is 1 and the second is 0)
#' x 0.5 (offspring: \eqn{0(1/2) + 1(1/2)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2 \times (q^2 \times p^2)} (this could be 0 x 2 or 2 x 0) x 1 (offspring: \eqn{1(1/1)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2pq \times 2pq} (both parents are 1) x 1 (offspring: \eqn{0(1/4) + 1(1/2) + 2(1/4)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{2 \times (2pq \times q2)} (this could be 1 x 2 or 2 x 1) x 1.5 (offspring: \eqn{1(1/2) + 2(1/2)}) +
#'
#'  &nbsp; &nbsp;  &nbsp; &nbsp;  &nbsp; &nbsp; \eqn{p^2 \times p^2} (both parents are 2) x 2 (offspring: \eqn{2(1/1)})
#'
#' Note that the use of \code{na.action = "expected"} is recommended when
#' a large number of offspring will conform the hybrid cross (such as
#' with bulked DNA analyses) for family groups with reasonable number of individuals.
#'
#' \strong{Warning}. If \code{"expected"} is used for \code{heterozygote.action} or \code{na.action},
#' direct transformation of the molecular data to other codings (\emph{e.g.},
#' dominance matrix coded as \code{c(0,1,0)}) is not recommended.
#'
#' @return
#' A molecular matrix \eqn{\boldsymbol{M}} containing the genotypes generated/imputed for the
#' hypothetical cross.
#'
#' @export
#' @md
#'
#' @examples
#' # Create dummy pedigree (using first 10 as parents).
#' ped <- data.frame(
#'  male = rownames(geno.apple)[1:5],
#'  female = rownames(geno.apple)[6:10])
#' ped$offs <- paste(ped$male, ped$female, sep = "_")
#' ped
#'
#' # Select portion of M for parents.
#' Mp <- geno.apple[c(ped$male, ped$female), 1:15]
#'
#' # Get genotype of crosses removing markers with heterozygotes.
#' synthetic.cross(
#'  M = Mp, ped = ped,
#'  indiv = "offs", mother = "female", father = "male",
#'  heterozygote.action = "exact",
#'  na.action = "useNA")
#'
#' # Request the synthetic cross to be NA in the respective samples.
#' synthetic.cross(
#'  M = Mp, ped = ped,
#'  indiv = "offs", mother = "female", father = "male",
#'  heterozygote.action = "useNA",
#'  na.action = "useNA")
#'
#' # Get genotype of crosses and use expected values.
#' synthetic.cross(
#'  M = Mp, ped = ped,
#'  indiv = "offs", mother = "female", father = "male",
#'  heterozygote.action = "expected", na.action = "expected")
#'

synthetic.cross <- function(M = NULL, ped = NULL, indiv = NULL, mother = NULL, father = NULL,
                            heterozygote.action = c("useNA", "exact", "fail", "expected"),
                            na.action = c("useNA", "expected"),
                            # n.cores = 1,
                            message = TRUE){

  # Traps ---------------------------------------------------------------------------------------

  # Check na.action.
  na.action <- match.arg(na.action)

  # Check heterozygote.action.
  heterozygote.action <- match.arg(heterozygote.action)

  # Check M class.
  check.data_(data_ = "M", class_ = "matrix")

  # Check ped class.
  check.data_(data_ = "ped", class_ = "data.frame")

  # Check indiv.
  check.args_(data_ = "ped", mandatory_ = TRUE, arg_ = "indiv", class_ = "character")

  # Check indiv.
  check.args_(data_ = "ped", mandatory_ = TRUE, arg_ = "mother", class_ = "character")

  # Check indiv.
  check.args_(data_ = "ped", mandatory_ = TRUE, arg_ = "father", class_ = "character")

  # # Check ped names.
  # ped.name.hit <- c(indiv, mother, father) %in% names(ped)
  # if (!all(ped.name.hit)){
  #   stop("Value provided to argument \'", c("indiv", "mother", "father")[!ped.name.hit],
  #        "' does not correspond to a variable in \'ped'.")
  # }

  # Check heterozygote action.
  heterozygote.action <- match.arg(heterozygote.action)

  # Check na action.
  na.action <- match.arg(na.action)

  # Check message.
  check.logical_(arg_ = "message")

  # Stop if heterozygous action is fail.
  if (heterozygote.action == "fail"){

    # Identify markers with heterozygotes and remove.
    het.markers <- apply(X = M, MARGIN = 2, FUN = function(m) {any(m == 1, na.rm = TRUE)})
    if (any(het.markers)){
      stop("Stop requested as some of the markers have have heterozygous states (1), e.g., ",
           paste0(names(head(het.markers[het.markers], 5)), collapse = ", "), "...")
    }
  }

  # Check if all parents have genotype in M.
  if(!all(c(ped[[mother]], ped[[father]]) %in% rownames(M))){
    stop("Some parents do not have genotypic information in matrix M.")
  }

  # Data manipulation traps ---------------------------------------------------------------------

  # Get unique combinations.
  unique.comb <- paste(ped[[mother]], ped[[father]], sep = "_")

  # Get reciprocals.
  unique.comb.rec <- paste(ped[[father]], ped[[mother]], sep = "_")
  recs <- unique.comb %in% unique.comb.rec

  # Identify duplicates.
  dups <- duplicated(unique.comb)
  if(any(dups)){

    # Report duplicates.
    message("A total of ", sum(dups), " duplicated rows were found in the supplied pedigree.")
    message("Removing \'ped' lines: ", which(dups), ".")

    # Remove duplicates.
    ped <- ped[!dups, ]

  }

  if(any(recs)){

    # Report reciprocals.
    warning("Reciprocals found in the supplied pedigree. These were not removed.")

    # TODO try to remove reciprocals?
  }

  rm(dups, recs)

  # Prepare data --------------------------------------------------------------------------------

  # Collect initial number of markers.
  initial.n.markers <- ncol(M)

  # Check and remove eventual heterozygotes in data if exact method is chosen.
  if (heterozygote.action == "exact"){

    if (message){
      message(blue("\nExact method selected. Eliminating markers containing one or more heterozygotic read."))
    }

    # Identify markers with heterozygotes and remove.
    het.markers <- apply(X = M, MARGIN = 2, FUN = function(m) {any(m == 1, na.rm = TRUE)})
    M <- M[, !het.markers, drop = FALSE]

    if (message){
      message("Total number of dropped markers: ", initial.n.markers - ncol(M))
      message("Total number of remaining markers: ", ncol(M))
    }

    # Stop if no marker left.

    if (ncol(M) == 0){
      stop("No marker were left after removal of heterozygote reads.")
    }
  }

  # Generate hybrid space -----------------------------------------------------------------------

  # Collect marker names (this is necessary in some cases, e.g., 1 marker left).
  marker.ids <- colnames(M)

  # Get range to loop through.
  range.hybrids <- 1:nrow(ped)

  if (message){
    message(blue("\nGenerating in hypothetical crosses genotypic information."))
  }


  # Call function to get offspring genotype.
  get.off.gen <- function(cur.offspring = NULL){

    # Collect genotype of combination.
    M.mother <- M[ped[cur.offspring, mother],]
    M.father <- M[ped[cur.offspring, father],]
    M.offspring <- rbind(M.mother, M.father)

    if (heterozygote.action == "useNA"){
      M.offspring[M.offspring %in% 1] <- NaN
    }

    # Get genotypes with NA.
    if (na.action == "useNA"){
      return(colMeans(M.offspring))
    }

    # Get genotypes with expected values.
    if (na.action == "expected"){

      # Initiate evo object.
      evo <- NULL

      # Get expected means.
      kid <- colMeans(M.offspring)

      # Identify genotypes missing in the current kid
      missing.markers <- is.na(kid) & !is.nan(kid)

      if (any(missing.markers)){

        # Get the frequency of pseudo-major allele 2.
        freq2 <- colMeans(M, na.rm = T)/2

        # Loop across all missing markers from the current kid.
        evo <- sapply(X = which(missing.markers), function(m) {

          # Get current marker.
          marker <- M.offspring[, m]

          # Get pseudo-p.
          f.pseudo.major <- freq2[m]

          # Get pseudo-q.
          f.pseudo.minor <- (1 - f.pseudo.major)

          # Logic: pseudo-q and p are required so we know which frequency to use on
          # the multiplications below. The MAF is not good here because
          # the major allele might not be represented by a 2 in the molecular matrix as
          # we dont know where this data comes from. Using pseudo-q and p, there is no
          # need to know the reference allele! We have to mach the genotype with
          # the possible offspring.

          # If there is one missing parent.
          if(sum(is.na(marker)) == 1){
            par.gen.no.miss <- sum(marker, na.rm = TRUE)

            # If the genotype of the non-missing parent is 0.
            if (par.gen.no.miss == 0){
              evo <-
                # q2
                f.pseudo.minor^2 * 0 +
                # 2pq
                2 * f.pseudo.minor * f.pseudo.major * 0.5 +
                # q2
                f.pseudo.major^2 * 1
            }

            # If the genotype of the non-missing parent is 1.
            if (par.gen.no.miss == 1){

              evo <-
                # q2
                f.pseudo.minor^2 * 0.5 +
                # 2pq
                2 * f.pseudo.minor * f.pseudo.major * 1 +
                # q2
                f.pseudo.major^2 * 1.5
            }

            # If the genotype of the non-missing parent is 2.
            if (par.gen.no.miss == 2){
              evo <-
                # q2
                f.pseudo.minor^2 * 1 +
                # 2pq
                2 * f.pseudo.minor * f.pseudo.major * 1.5 +
                # q2
                f.pseudo.major^2 * 2
            }
          }

          # If there are two missing parents.
          if(sum(is.na(marker)) == 2){
            f.pseudo.major <- freq2[m]

            # All possible combiations of unknown parents.
            evo <-

              # q2 x q2 (0 x 0)
              f.pseudo.minor^2 * f.pseudo.minor^2   * 0 +

              # 2 * q2 x 2pq (this could be 0 x 1 or 1 x 0)
              2 * (f.pseudo.minor^2 * 2 * f.pseudo.minor * f.pseudo.major)   * 0.5+

              # 2 * q2 x p2 (this could be 0 x 2 or 2 x 0)
              2 * (f.pseudo.minor^2 * f.pseudo.major^2)   * 1 +

              # 2pq x 2pq
              2 * f.pseudo.minor * f.pseudo.major * 2 * f.pseudo.minor * f.pseudo.major   *  1 +

              # 2 * 2pq x q2 (this could be 1 x 2 or 2 x 1)
              2 * (2 * f.pseudo.minor * f.pseudo.major * f.pseudo.major^2)   * 1.5 +

              # p2 x p2
              f.pseudo.major^2 * f.pseudo.major^2   * 2

            # Simplification (not sure if works on all cases - probably yes).
            # evo <- f.pseudo.minor^2 * 0 +
            #   2 * f.pseudo.minor * f.pseudo.major * 1 +
            #   f.pseudo.major^2 * 2

          }

          # Return the expected value of the offspring.
          return(evo)
        })


      }

      # Replace na in kid.
      kid[which(missing.markers)] <- evo

      # Return the imputed kid.
      return(kid)
    }
  }

  # Run function.
  # if (n.cores > 1){
  #   M <- mclapply(X = range.hybrids, mc.cores = n.cores, FUN = get.off.gen)
  # } else {
    M <- lapply(X = range.hybrids, FUN = get.off.gen)
  # }

  # Get hybrid matrix.
  M <- do.call(rbind, M)

  # Replace NaN with NA because of heterozygotes.
  M[is.nan(M)] <- NA

  # Add names of hybrids to new matrix.
  rownames(M) <- ped[[indiv]]

  # Add maker ids.
  colnames(M) <- marker.ids

  # Return.
  return(M)
}

