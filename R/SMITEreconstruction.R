#' SMITE Reconstruction
#'
#' Use a model (m) created by SMITE.calib() to estimate the reconstruction target (b) based on the 'forward' coral variable matrix (A).
#' @param A The 'forward' coral variable matrix (t x p; time in rows, coral parameters in columns).
#' @param m A vector of SMITE model parameters (length p), ideally created from the SMITE.calib() function.
#' @param Amu A vector of length p giving the mean values of each coral variable in the calibration dataset.
#' @param Asd A vector of length p giving the standard deviations of each coral variable in the calibration dataset.
#' @param bmu A single-element atomic vector giving the mean of the reconstruction target in the calibration dataset.
#' @param bsd A single-element atomic vector giving the standard deviation of the reconstruction target in the calibration dataset.
#' @return bhat - Predicted reconstruction target values. The SEP from the SMITE.calib() function can be used to estimate minimum errors (e.g., not associated with age model error).
#' @export

SMITE.recon <- function(A, m, Amu, Asd, bmu, bsd) {

  # =============================== #
  # Formatting
  # =============================== #

  A <- as.matrix(A)
  Amu <- as.vector(Amu)
  Asd <- as.vector(Asd)

  m <- as.matrix(m)
  bmu <- as.vector(bmu)
  bsd <- as.vector(bsd)

  # =============================== #
  # Flags
  # =============================== #

  if(ncol(A) != length(m)) {
    stop("The model parameters (m) must be the same length as the number of columns in the forward matrix (A).")
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

  bnorm <- Anorm %*% m
  bhat <- (bnorm * bsd) + bmu

  return(bhat)
}
