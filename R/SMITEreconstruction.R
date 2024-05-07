#' SMITE Reconstruction
#'
#' Use a model (x) created by SMITE.calib() to estimate the reconstruction target (b) based on the 'forward' coral variable matrix (A).
#' @param A The 'forward' coral variable matrix (t x p; time in rows, coral parameters in columns).
#' @param x A vector of SMITE model parameters (length p), ideally created from the SMITE.calib() function.
#' @param Amu A vector of length p giving the mean values of each coral variable in the calibration dataset.
#' @param Asd A vector of length p giving the standard deviations of each coral variable in the calibration dataset.
#' @param bmu A single-element atomic vector giving the mean of the reconstruction target in the calibration dataset.
#' @param bsd A single-element atomic vector giving the standard deviation of the reconstruction target in the calibration dataset.
#' @return bhat - Predicted reconstruction target values. The SEP from the SMITE.calib() function can be used to estimate minimum errors (e.g., not associated with age model error).
#' @export
#' @examples
#' # Load data from Hughes et al. (in revision)
#' library(SMITER)
#' data(BMDA_1B_Comp)
#' data(BMDA_1B_EComp)
#'
#' # Designate forward matrix (A), reconstruction target (b), and corresponding errors (Ae, be).
#' proxy_names <- names(bmda_1b_comp[,c(5:ncol(bmda_1b_comp))])
#'
#' A <- bmda_1b_comp[,proxy_names]
#' Ae <- bmda_1b_ecomp[,proxy_names]
#' b <- bmda_1b_comp[,'Temp']
#' be <- rep(0.02, length(b))
#'
#' # Execute SMITE Calibration with no truncation (i.e., all information retained).
#' SMITE <- SMITE.calib(A = A, b = b, Ae = Ae, be = be, eigenclean = ncol(A))
#'
#' # Extract SMITE model parameters
#' x <- SMITE$x$x.mu
#'
#' # Reconstruction (should be nearly identical to result$recon$bhat.mu, only minor differences due to bootstrapping)
#' bhat <- SMITE.recon(
#'   A = A,
#'   x = x,
#'   Amu = colMeans(A),
#'   Asd = apply(A, 2, function(x) sd(x)),
#'   bmu = mean(b),
#'   bsd = sd(b)
#' )
#'
#' cbind(bhat, result$recon$bhat.mu)

SMITE.recon <- function(A, x, Amu, Asd, bmu, bsd) {

  # =============================== #
  # Formatting
  # =============================== #

  A <- as.matrix(A)
  Amu <- as.vector(Amu)
  Asd <- as.vector(Asd)

  x <- as.matrix(x)
  bmu <- as.vector(bmu)
  bsd <- as.vector(bsd)

  # =============================== #
  # Flags
  # =============================== #

  if(ncol(A) != length(x)) {
    stop("The model parameters (x) must be the same length as the number of columns in the forward matrix (A).")
  }

  if(ncol(A) != length(Amu)) {
    stop("The means of the forward matrix from the calibration dataset (Au) must be the same length as the number of columns in the forward matrix (A).")
  }

  if(ncol(A) != length(Asd)) {
    stop("The standard deviations of the forward matrix from the calibration dataset (Asd) must be the same length as the number of columns in the forward matrix (A).")
  }

  if(length(bmu) != 1 | length(bsd) != 1) {
    stop("The mean or standard deviation of the reconstruction target from the calibration dataset (bmu or bsd) should only be length 1.")
  }

  # =============================== #
  # Create normal forward matrix
  # =============================== #

  Anorm <- matrix(NA, nrow = nrow(A), ncol = ncol(A))

  for(i in 1:ncol(A)) {
    Anorm[,i] <- (A[,i] - Amu[i]) / Asd[i]
  }

  # =============================== #
  # Reconstruction
  # =============================== #

  bnorm <- Anorm %*% x
  bhat <- (bnorm * bsd) + bmu

  return(bhat)
}
