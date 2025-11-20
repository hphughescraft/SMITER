#' SMITE Reconstruction
#'
#' Use a model (x) created by SMITE.calib() to estimate the reconstruction target (b) based on the 'forward' coral variable matrix (A).
#' @param A The 'forward' coral variable matrix (t x p; time in rows, coral parameters in columns).
#' @param x A matrix of individual bootstrapped SMITE model parameters (p x it), ideally created from the 'x_bootstrap' SMITE.calib() function.
#' @param Amu A vector of length p giving the mean values of each coral variable in the calibration dataset.
#' @param Asd A vector of length p giving the standard deviations of each coral variable in the calibration dataset.
#' @param bmu A single-element atomic vector giving the mean of the reconstruction target in the calibration dataset.
#' @param bsd A single-element atomic vector giving the standard deviation of the reconstruction target in the calibration dataset.
#' @param it The number of bootstrap iterations that will be performed to obtain confidence intervals on reconstructed values.
#' @param alpha The significance level for the confidence interval (i.e., 0.05 = 95-percent confidence).
#' @param res_model An object from SMITE.calib(), bootstrapped iterations of the standard deviation of the residuals (sigma) and the AR1&2 coefficients ('ar1', 'ar2').
#' @param res_model_type select between c('gaussian', 'ar1', and 'ar2') to determine how residual errors will be propagated.
#' @return bhat - Predicted reconstruction target values with alpha confidence intervals.
#' @importFrom stats sd cor.test quantile na.omit acf filter rnorm
#' @export
#' @examples
#'
#' # Load data from Hughes et al. (2024)
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
#'
#' # Reconstruction (errors from model parameters, observations, and calibration are propagated into reconstruction).
#' bhat <- SMITE.recon(
#'   A = A,
#'   x = SMITE$x_bootstrap,
#'   Amu = colMeans(A),
#'   Asd = apply(A, 2, sd),
#'   bmu = mean(b),
#'   bsd = sd(b),
#'   res_model = SMITE$res_model,
#'   res_model_type = 'gaussian', # Can choose between 'gaussian', 'ar1', or 'ar2'.
#'   alpha = 0.05,
#'   it = 10000
#' )
#'
#' cbind(bhat, SMITE$bhat$Mu)

SMITE.recon <- function(A, x, Amu, Asd, bmu, bsd,
                        res_model, res_model_type = c("gaussian", "ar1", "ar2"),
                        it = 10000, alpha = 0.05) {

  # =============================== #
  # Formatting
  # =============================== #

  A <- as.matrix(A)
  Ae <- as.matrix(Ae)
  Amu <- as.vector(Amu)
  Asd <- as.vector(Asd)

  bmu <- as.vector(bmu)
  bsd <- as.vector(bsd)

  x <- as.matrix(x)

  # =============================== #
  # Flags
  # =============================== #

  if(ncol(A) != length(Amu)) {
    stop("The means of the forward matrix from the calibration dataset (Au) must be the same length as the number of columns in the forward matrix (A).")
  }

  if(ncol(A) != length(Asd)) {
    stop("The standard deviations of the forward matrix from the calibration dataset (Asd) must be the same length as the number of columns in the forward matrix (A).")
  }

  if(length(bmu) != 1 | length(bsd) != 1) {
    stop("The mean or standard deviation of the reconstruction target from the calibration dataset (bmu or bsd) should only be length 1.")
  }

  if(length(res_model_type) > 1) {
    stop("Choose a residual model type (gaussian, ar1, ar2).")
  }

  if(is.null(alpha)) {
    stop("Please define the significance level 'alpha'.")
  }

  alpha_symm <- alpha / 2


  # =============================== #
  # Create normal forward matrix
  # =============================== #

  Anorm <- sweep(A, 2, Amu, "-")
  Anorm <- sweep(Anorm, 2, Asd, "/")
  Aenorm <- sweep(Ae, 2, Asd, "/")

  # =============================== #
  # Reconstruction
  # =============================== #

  bhat_mat <- matrix(NA, nrow = nrow(A), ncol = it)

  for(i in 1:it) {
    si <- sample(c(1:it), size = 1)

    # SMITE model parameters
    xi <- x[,si]

    # A-matrix resample
    Ap <- matrix(
      rnorm(length(Anorm), mean = c(Anorm), sd = c(Aenorm)),
      nrow = nrow(Anorm)
    )

    # Residuals
    if(res_model_type == 'gaussian') {
      res <- rnorm(nrow(A), mean = 0, sd = res_model['sigma', si])

    } else if(res_model_type == 'ar1') {
      ar1 <- res_model['ar1', si]
      se <- res_model['sigma', si] * sqrt(pmax(0, 1 - ar1^2))
      wn <- rnorm(nrow(A), mean = 0, sd = se)

      res <- as.numeric(filter(wn, filter = ar1, method = "recursive"))

    } else if(res_model_type == 'ar2') {
      ar1 <- res_model['ar1', si]
      ar2 <- res_model['ar2', si]
      se <- res_model['sigma', si] * sqrt(pmax(0, 1 - ar1^2 - ar2^2))
      wn <- rnorm(nrow(A), mean = 0, sd = se)

      res <- as.numeric(filter(wn, filter = c(ar1, ar2), method = "recursive"))

    }

    bhat_mat[,i] <- (((Ap %*% xi) * bsd) + bmu) + res

  }

  bhat <- as.data.frame(
    t(
      apply(
        bhat_mat, 1, function(x)
          quantile(x, probs = c(alpha_symm, 0.5, 1 - alpha_symm))
        )
      )
    )

  names(bhat) <- c("Low", "Mu", "High")

  return(bhat)
}
