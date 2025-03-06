#' SMITE Calibration
#'
#' Create a model (x) based on the 'forward' coral variable matrix (A) and the reconstruction target (b).
#' @param A The 'forward' coral variable matrix (t x p; time in rows, coral parameters in columns).
#' @param b The reconstruction target (t x 1).
#' @param Ae A t x p matrix containing the errors for the forward matrix.
#' @param be Error estimates for the reconstruction target. It can either be one value or a vector of equal length to the reconstruction target (t x 1).
#' @param it The number of bootstrap Monte Carlo iterations (recommended at 10,000).
#' @param noise Either 'white' or 'red', describing how the noise is propagated into the model.
#' @param acc Autocorrelation Coefficient; if noise is specified as 'red', it describes the degree of autocorrelation in the error term.
#' @param eigenclean Describes how (if at all) the singular values should be truncated to 'clean' the inverse solution. If the value is between 0 and 1, it will remove singular values based on the cumulative variance explained. If the value is greater than 1, it will return that many singular values (from the highest order).
#' @param alpha The significance level for the confidence interval (i.e., 0.05 = 95-percent confidence).
#' @return S - Confidence intervals on the singular values of the forward matrix,
#' @return x - SMITE model parameters, or loadings, for each column of the forward matrix.
#' @return bhat - Predicted target anomalies
#' @return recon - Predicted absolute target values
#' @return e - Relevant error metrics. It enumerates the Standard Error of Prediction (SEP), Root-Mean-Square-Error (RMSE), and the correlation coefficient (r).
#' @return e.xval - Relevant error metrics for the cross-validated data points specifically.
#' @export
#' @examples
#' # Load data from Hughes et al. (in revision)
#' library(SMITER)
#' data(BMDA_1B_Comp)
#' data(BMDA_1B_EComp)
#'
#' # Designate forward matrix (A), reconstruction target (b), and corresponding errors (Ae, be).
#'
#' proxy_names <- names(bmda_1b_comp[,c(5:ncol(bmda_1b_comp))])
#'
#' A <- bmda_1b_comp[,proxy_names]
#' Ae <- bmda_1b_ecomp[,proxy_names]
#' b <- bmda_1b_comp[,'Temp']
#' be <- rep(0.02, length(b))
#'
#' # Execute SMITE Calibration with no truncation (i.e., all information retained).
#' result <- SMITE.calib(A = A, b = b, Ae = Ae, be = be, eigenclean = ncol(A))


SMITE.calib <- function(A, b, Ae = NULL, be = NULL, it = 10000, noise = "white", acc = NULL,
                        eigenclean = NULL, alpha = 0.05, xval = NULL) {


  # =============================== #
  # Flags
  # =============================== #

  if(is.null(A) | is.null(b)) {
    stop("Must define the forward matrix (A) and the reconstruction target (b).")
  } else {
    A <- as.matrix(A)
    b <- as.matrix(b)
  }

  if(is.null(Ae) | is.null(be)) {
    stop("Must define errors in the forward matrix (Ae) and the errors in the reconstruction target (be).")
  } else {
    Ae <- as.matrix(Ae)
    be <- as.matrix(be)
  }

  if(any(dim(Ae) != dim(A))) {
    stop("Error matrix (Ae) must be same dimensions as the forward matrix (A).")
  }


  if(any(dim(be) != dim(b))) {
    stop("Error estimates of the reconstruction target (be) must be either one value or the same dimensions as the reconstruction target (b).")
  }

  if(is.null(noise)) {
    stop("Please specify whether the noise is 'red' (autocorrelated) or 'white' (Gaussian).")
  }

  if(noise == 'red' & is.null(acc)) {
    stop("Please specify the Autocorrelation Coefficient for the red noise.")
  }

  if(is.null(xval)) {
    print("No window size for cross-validation (xval) specified. Using minimum window size of 1.")
    xval <- 1
  }

  if(xval <= 0 | xval > (nrow(A) / 2)) {
    stop("Cross-validation window should be at least 1 and no larger than half the size of the data.")
  }

  # =============================== #
  # Begin bootstrap
  # =============================== #
  k_list <- list()

  for(k in 1:it) {

    # =============================== #
    # Perturb A (Ap) & b (bp)
    # =============================== #

    # Preallocate
    Ap <- matrix(data = NA, nrow = nrow(A), ncol = ncol(A))
    bp <- matrix(data = NA, nrow = nrow(b), ncol = ncol(b))

    # =============================== #
    # Nested For-Loop: Ap allocation
    # =============================== #

    if(noise == 'white') {
      for(j in 1:ncol(Ap)) {
        for(i in 1:nrow(Ap)) {
          # Perturb each cell of G by the error defined in the corresponding cell of Ge
          Ap[i,j] <- rnorm(n = 1, mean = A[i,j], sd = Ae[i,j])
        }
      }
    } else if(noise == 'red') {

      Ar <- matrix(data = NA, nrow = nrow(A), ncol = ncol(A))

      for(j in 1:ncol(Ar)) {
        for(i in 1:nrow(Ar)) {
          if(i == 1) {
            Ar[i,j] <- rnorm(1, 0, Ae[i,j])
          } else {
            Ar[i,j] <- (Ar[i - 1, j] * acc) + rnorm(1, 0, Ae[i,j])
          }
        }
      }

      Ap <- A + Ar
    } else {
      stop("Please specify whether the noise is red or white.")
    }


    for(i in 1:nrow(Ap)) {
      bp[i] <- rnorm(n = 1, mean = b[i], sd = be[i])
    }

    Ap_norm <- apply(Ap, MARGIN = 2, FUN = function(x) (x - mean(x)) /  sd(x))
    bp_norm <- (bp - mean(bp)) / sd(bp)

    # =============================== #
    # Cross-validation
    # =============================== #

    # Set starting index
    si <- sample(1:(nrow(Ap) - xval + 1), 1)

    # Store indexed rows
    rr <- si:(si + xval - 1)
    Apv <- Ap_norm[rr,]
    bpv <- bp_norm[rr,]

    # Remove indexed rows
    Apx <- Ap_norm[-(rr),]
    bpx <- bp_norm[-(rr),]



    # =============================== #
    # SVD of Ap
    # =============================== #
    G <- svd(Apx)

    U <- G$u
    V <- G$v
    S <- G$d

    # =============================== #
    # Eigenclean
    # =============================== #
    if(is.null(eigenclean)) {
      eigenclean <- ncol(A)
    }

    if(eigenclean > 0 & eigenclean < 1) { # Cumulative Variance Explained
      p <- (cumsum(S)/sum(S)) < eigenclean
      S <- S[p]
      U <- U[,p]
      V <- V[,p]
    } else if(eigenclean >= 1) { # Number of singular values returned
      p <- eigenclean
      S <- S[1:p]
      U <- U[,1:p]
      V <- V[,1:p]
    }

    # =============================== #
    # Moore-Penrose Pseudoinverse
    # =============================== #

    x <- V %*% diag(1/S) %*% t(U) %*% bpx
    bhat <- Apx %*% x

    # =============================== #
    # Store S, m, & bhat, & residuals
    # =============================== #

    # Place removed rows back into dhat
    if(si == 1) { # Prevents over-allocating from the beginning
      bhat <- c(rep(NA, xval),
                bhat[si:length(bhat)])
    } else if(si > length(bhat)){ # Prevents over-allocating from the end
      bhat <- c(bhat[1:(si - 1)],
                rep(NA, xval))
    } else { # Original splice code
      bhat <- c(bhat[1:(si - 1)],
                rep(NA, xval),
                bhat[si:length(bhat)])
    }

    res <- bhat - bp_norm
    res_xval <- (Apv %*% x) - bpv

    k_list[[k]] <- list(
      "S" = S,
      "x" = x,
      "bhat" = bhat,
      "res" = res,
      "res_xval" = res_xval
    )
  } # End Ap Allocation Bootstrap

  # =============================== #
  # Extract distributions
  # =============================== #

  S_mat <- matrix(NA, nrow = length(S), ncol = it)
  x_mat <- matrix(NA, nrow = nrow(x), ncol = it)
  bhat_mat <- matrix(NA, nrow = length(bhat), ncol = it)
  res_mat <- matrix(NA, nrow = length(res), ncol = it)
  res_xval_mat <- matrix(NA, nrow = length(res_xval), ncol = it)

  for(i in 1:it) {
    S_mat[,i] <- k_list[[i]]$S
    x_mat[,i] <- k_list[[i]]$x
    bhat_mat[,i] <- k_list[[i]]$bhat
    res_mat[,i] <- k_list[[i]]$res
    res_xval_mat[,i] <- k_list[[i]]$res_xval
  }

  # ======================================= #
  # Calculate means and confidence interval
  # ======================================= #

  if(is.null(alpha)) {
    stop("Please define the significance level 'alpha'.")
  }

  alpha_symm <- alpha / 2
  cint_index <- seq(floor(alpha_symm * it), floor((1 - alpha_symm) * it), 1)

  # ====== #
  # S
  # ====== #
  S_mu <- apply(S_mat, 1, function(x) mean(x))
  S_sd <- apply(S_mat, 1, function(x) sd(x))
  S_cint <- apply(S_mat, 1, function(x) x[order(x)][cint_index])
  S_low <- S_cint[1,]
  S_high <- S_cint[nrow(S_cint),]

  S_df <- data.frame(
    "S.low" = S_low,
    "S.mu" = S_mu,
    "S.high" = S_high
  )

  # ====== #
  # x
  # ====== #
  x_mu <- apply(x_mat, 1, function(f) mean(f))
  x_sd <- apply(x_mat, 1, function(f) sd(f))
  x_cint <- apply(x_mat, 1, function(f) f[order(f)][cint_index])
  x_low <- x_cint[1,]
  x_high <- x_cint[nrow(x_cint),]

  x_df <- data.frame(
    "x.low" = x_low,
    "x.mu" = x_mu,
    "x.high" = x_high
  )

  rownames(x_df) <- colnames(A)

  # ============ #
  # bhat & recon
  # ============ #
  bhat_mu <- apply(bhat_mat, 1, function(f) mean(f, na.rm = TRUE))
  bhat_sd <- apply(bhat_mat, 1, function(f) sd(f, na.rm = TRUE))
  bhat_cint <- apply(bhat_mat, 1, function(f) f[order(f)][cint_index])
  bhat_low <- bhat_high <- rep(NA, length(bhat_mu))

  for(i in 1:ncol(bhat_cint)) {
    m <- na.omit(bhat_cint[,i])
    n <- m[ceiling(alpha_symm * length(m)):ceiling((1 - alpha_symm) * length(m))]

    bhat_low[i] <- min(n)
    bhat_high[i] <- max(n)
  }

  bhat_df <- data.frame(
    "bhat.low" = bhat_low,
    "bhat.mu" = bhat_mu,
    "bhat.high" = bhat_high
  )

  recon_df <- (bhat_df * sd(b)) + mean(b)

  # ============ #
  # Residuals
  # ============ #
  res_mu <- apply(res_mat, 1, function(f) mean(f, na.rm = TRUE))
  res_sd <- apply(res_mat, 1, function(f) sd(f, na.rm = TRUE))
  res_cint <- apply(res_mat, 1, function(f) f[order(f)][cint_index])
  res_low <- res_high <- rep(NA, length(res_mu))

  for(i in 1:ncol(res_cint)) {
    m <- na.omit(res_cint[,i])
    n <- m[ceiling(alpha_symm * length(m)):ceiling((1 - alpha_symm) * length(m))]

    res_low[i] <- min(n)
    res_high[i] <- max(n)
  }

  res_df <- data.frame(
    "res.low" = res_low * sd(b),
    "res.mu" = res_mu * sd(b),
    "res.high" = res_high * sd(b)
  )

  # Cross-validated

  res_xval_mu <- apply(res_xval_mat, 1, function(f) mean(f, na.rm = TRUE))
  res_xval_sd <- apply(res_xval_mat, 1, function(f) sd(f, na.rm = TRUE))
  res_xval_cint <- apply(res_xval_mat, 1, function(f) f[order(f)][cint_index])
  res_xval_low <- res_xval_high <- rep(NA, length(res_xval_mu))

  for(i in 1:ncol(res_xval_cint)) {
    m <- na.omit(res_xval_cint[,i])
    n <- m[ceiling(alpha_symm * length(m)):ceiling((1 - alpha_symm) * length(m))]

    res_xval_low[i] <- min(n)
    res_xval_high[i] <- max(n)
  }

  res_xval_df <- data.frame(
    "res.low" = res_xval_low * sd(b),
    "res.mu" = res_xval_mu * sd(b),
    "res.high" = res_xval_high * sd(b)
  )

  # =============================== @
  # RETURN
  # =============================== @

  SMITE_list <- list(
    "S" = S_df,
    "x" = x_df,
    "bhat" = bhat_df,
    "recon" = recon_df,
    "e" = c(
      "SEP" = (mean(c(bhat_df$bhat.high - bhat_df$bhat.mu,
                      bhat_df$bhat.mu - bhat_df$bhat.low)) / 1.96) * sd(b),
      "RMSE" = (sum(sqrt((b - recon_df$bhat.mu)^2))) / (nrow(recon_df)),
      "r" = cor.test(b, bhat_df$bhat.mu)$estimate
    ),
    "e.xval" = c(
      "SEP" = (mean(c((b - res_df$res.mu) - (b - res_df$res.high),
                      (b - res_df$res.low) - (b - res_df$res.mu))) / 1.96),
      "RMSE" = mean(abs(res_df$res.mu)),
      "r" = cor.test(b, (b - res_df$res.mu))$estimate
    )
  )

  return(SMITE_list)

}
