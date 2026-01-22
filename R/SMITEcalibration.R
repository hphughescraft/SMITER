#' SMITE Calibration
#'
#' Create a model (x) based on the 'forward' coral variable matrix (A) and the reconstruction target (b).
#' @param A The 'forward' coral variable matrix (t x p; time in rows, coral parameters in columns).
#' @param b The reconstruction target (t x 1).
#' @param Ae A t x p matrix containing the errors for the forward matrix.
#' @param be Error estimates for the reconstruction target. It can either be one value or a vector of equal length to the reconstruction target (t x 1).
#' @param it The number of bootstrap Monte Carlo iterations (recommended at 10,000).
#' @param eigenclean Describes how (if at all) the singular values should be truncated to 'clean' the inverse solution. If the value is between 0 and 1, it will remove singular values based on the cumulative variance explained. If the value is greater than 1, it will return that many singular values (from the highest order).
#' @param alpha The significance level for the confidence interval (i.e., 0.05 = 95-percent confidence).
#' @param xval Cross-validation window size.
#' @param weights TRUE or FALSE on whether you want to incorporate uncertainty (Ae and be) into calibration via an iterative Jacobian reweighting loop. Default is TRUE.
#' @param itj The number of iterations for the Jacobian reweighting. Default is 25.
#' @param tol The tolerance for the iterative Jacobian reweighting. Default is 1e-8.
#' @return S - Confidence intervals on the singular values of the forward matrix,
#' @return x - Summary of SMITE model parameters, or loadings, for each column of the forward matrix.
#' @return x_bootstrap - Individual bootstrap iterations of SMITE model parameters, needed for propagating uncertainty in predictions.
#' @return res_model - Standard deviation and autoregressive coefficients (1 & 2) of the residuals of each bootstrap.
#' @return bhat - Predicted target values from in-sample blocks
#' @return bhat_xval - predicted target values from out-of-sample blocks
#' @return e_stats - Confidence intervals for R and RMSE (non-cross-validated and cross-validated). The size of the confidence interval is based on alpha.
#' @importFrom stats sd cor.test quantile na.omit acf filter rnorm
#' @export
#' @examples
#' # Load data from Hughes et al. (2024)
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


SMITE.calib <- function(A, b, Ae = NULL, be = NULL, it = 10000,
                        eigenclean = NULL, alpha = 0.05, xval = NULL,
                        weights = TRUE, itj = 25, tol = 1e-8) {


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

  if(is.null(xval)) {
    print("No window size for cross-validation (xval) specified. Using nrow(A) / 5.")
    xval <- ceiling(nrow(A) / 5)
  }

  if(xval <= 0) {
    stop("Cross-validation window should be at least 1.")
  }

  # =============================== #
  # Begin bootstrap
  # =============================== #
  k_list <- list()

  for(k in 1:it) {

    # =============================== #
    # Perturb A (Ap) & b (bp)
    # =============================== #

    # Initialize
    Ap <- matrix(nrow = 0, ncol = ncol(A))
    Aep <- matrix(nrow = 0, ncol = ncol(Ae))
    bp <- matrix(nrow = 0, ncol = ncol(b))
    bep <- matrix(nrow = 0, ncol = ncol(be))

    # =============================== #
    # Nested For-Loop: Ap allocation
    # =============================== #

    block_track <- c()

    while(nrow(Ap) < nrow(A)) {
      start <- sample(1:(nrow(A) - xval + 1), 1)
      start
      block <- start:(start + xval - 1)

      A_block <- A[block,]
      Ae_block <- Ae[block,]
      b_block <- as.matrix(b[block,])
      be_block <- as.matrix(be[block,])

      Ap <- rbind(Ap, A_block)
      Aep <- rbind(Aep, Ae_block)
      bp <- rbind(bp, b_block)
      bep <- rbind(bep, be_block)

      block_track <- c(block_track, block)
    }

    # Trim to N
    Ap <- Ap[1:nrow(A),]
    Aep <- Aep[1:nrow(Ae),]
    bp <- bp[1:nrow(b),]
    bep <- bep[1:nrow(be),]
    block_track <- block_track[1:nrow(A)]

    Ap_norm <- apply(Ap, 2, function(x) (x - mean(x)) /  sd(x))

    # Divide Aep by sigma_Ap to normalized error
    Aep_norm <- sweep(Aep, 2, apply(Ap, MARGIN = 2, function(x) sd(x)), "/")

    bp_norm <- as.matrix((bp - mean(bp)) / sd(bp))
    bep_norm <- as.matrix(bep / sd(bp))


    # =============================== #
    # Cross-validation
    # =============================== #

    # Set starting index
    si <- sample(1:(nrow(Ap) - xval + 1), 1)

    # Store indexed rows
    rr <- si:(si + xval - 1)
    Apv <- Ap_norm[rr,]
    Aepv <- Aep_norm[rr,]
    bpv <- bp_norm[rr,]
    bepv <- bep_norm[rr,]


    # Remove indexed rows
    Apx <- Ap_norm[-(rr),]
    Aepx <- Aep_norm[-(rr),]
    bpx <- bp_norm[-(rr),]
    bepx <- bep_norm[-(rr),]

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
      print("No regularization (eigenclean) specified. Implementing 0 regularization.")
      eigenclean <- ncol(A)
    }

    if(eigenclean > 0 & eigenclean <= 1) { # Cumulative Variance Explained
      p <- which((cumsum(S^2)/sum(S^2)) <= eigenclean)

      if(length(p) < 2) {
        stop("Eigenclean CVE is too small. Select a larger ")
      }

      S <- S[p]
      U <- U[,p]
      V <- V[,p]
    } else if(eigenclean > 1) { # Number of singular values returned
      p <- eigenclean
      S <- S[1:p]
      U <- U[,1:p]
      V <- V[,1:p]
    }

    # =============================== #
    # Moore-Penrose Pseudoinverse
    # =============================== #

    x <- V %*% diag(1/S) %*% t(U) %*% bpx

    # =============================== #
    # Iterative Jacobian Reweighting
    # =============================== #

    # Setup constants
    # Predictor-error correlation
    r_A <- cov2cor(cov(Aepx))

    # Cross-covariance between Ae and be
    rho_Aebe <- apply(Aepx, 2, function(x) cor(x, bepx))
    rho_Aebe <- ifelse(is.na(rho_Aebe), 0, rho_Aebe) # Guard against NAs

    if(weights) {

      # Store original x
      x_og <- x

      for(j in 1:itj) {

        # Compute bhat from transient x
        bhat_itj <- Apx %*% x

        # Compute orthogonal distance
        r <- (bhat_itj - bpx) / sqrt(2)

        # =============================== #
        # Orthogonal residual uncertainty
        # =============================== #

        n <- nrow(Apx)
        p <- ncol(Apx)

        xv <- as.numeric(x)
        be_sig <- as.numeric(bepx)

        rperp_var <- numeric(n)

        for (i in seq_len(n)) {

          # Ae for ith observation
          Di <- as.numeric(Aepx[i, ])

          # Map Ae to b-space
          v <- Di * xv

          # Covariance of Ae mapped to b-space
          var_dbhat <- as.numeric(t(v) %*% r_A %*% v)

          # Variance of be
          var_db <- be_sig[i]^2

          # Cross-covariance term: x^T c_i
          cov_dbhat_db <- sum(xv * (rho_Aebe * Di * be_sig[i]))

          # Orthogonal residual variance
          rperp_var[i] <- 0.5 * (var_dbhat + var_db - 2 * cov_dbhat_db)
        }

        # Numerical guard
        rperp_var[rperp_var <= 0] <- min(rperp_var[rperp_var > 0], na.rm = TRUE)

        # Weights
        w <- 1 / sqrt(rperp_var)
        W <- diag(w)

        # Weighted SVD
        Apxw <- W %*% Apx
        bpxw <- W %*% bpx

        G <- svd(Apxw)
        U <- G$u
        V <- G$v
        S <- G$d

        # Regularization
        if (is.null(eigenclean)) {
          eigenclean <- ncol(A)
        }

        if (eigenclean > 0 & eigenclean <= 1) {

          p <- which((cumsum(S^2) / sum(S^2)) < eigenclean)
          S <- S[p]
          U <- U[, p, drop = FALSE]
          V <- V[, p, drop = FALSE]

        } else if (eigenclean > 1) {

          p <- eigenclean
          S <- S[1:p]
          U <- U[, 1:p, drop = FALSE]
          V <- V[, 1:p, drop = FALSE]

        }

        # Updated x
        x_j <- V %*% diag(1 / S) %*% t(U) %*% bpxw

        # Convergence check (relative)
        denom <- max(1e-12, sqrt(sum(x_og^2)))
        if (sqrt(sum((x_j - x_og)^2)) / denom < tol) {
          x <- x_j
          break
        }

        x_og <- x_j
        x <- x_j

      }
    }

    bhat <- c((Apx %*% x))
    bhat_xval <- c((Apv %*% x))

    # =============================== #
    # Store S, m, & bhat, & residuals
    # =============================== #

    # Statistics (use universal sd(b) for scaling)
    r <- cor.test(bhat, bpx)$estimate
    rmse <- sqrt(mean((bpx - bhat)^2)) * sd(b)

    r_xval <- cor.test(bpv, bhat_xval)$estimate
    rmse_xval <- sqrt(mean((bpv - bhat_xval)^2)) * sd(b)

    # Residuals
    res <- (bpx * sd(b)) - (bhat * sd(b))
    res_xval <- (bpv * sd(b)) - (bhat_xval * sd(b))


    k_list[[k]] <- list(
      "S" = S,
      "x" = x,
      "bhat" = (bhat * sd(b)) + mean(b),
      "bhat_xval" = (bhat_xval * sd(b)) + mean(b),
      "res" = res,
      "res_xval" = res_xval,
      "block_track" = block_track,
      "id_xval" = rr,
      "r" = r,
      "rmse" = rmse,
      "r_xval" = r_xval,
      "rmse_xval" = rmse_xval
    )
  } # End Ap Allocation Bootstrap

  # =============================== #
  # Extract distributions
  # =============================== #

  S_mat <- matrix(NA, nrow = length(S), ncol = it)
  x_mat <- matrix(NA, nrow = nrow(x), ncol = it)
  bhat_mat <- matrix(NA, nrow = length(b), ncol = it)
  bhat_xval_mat <- matrix(NA, nrow = length(b), ncol = it)
  res_mat <- matrix(NA, nrow = length(b), ncol = it)
  res_xval_mat <- matrix(NA, nrow = length(b), ncol = it)
  stats_mat <- matrix(NA, nrow = 4, ncol = it,
                      dimnames = list(
                        c("r", "rmse", "r_xval", "rmse_xval"),
                        seq(1, it, 1)
                      ))

  for(i in 1:it) {

    block_track <- k_list[[i]]$block_track

    # Validation blocks (mapped back to original data)
    idx <- block_track[k_list[[i]]$id_xval]

    # Training blocks (mapped back to original data)
    id <- block_track[-k_list[[i]]$id_xval]

    S_mat[,i] <- k_list[[i]]$S
    x_mat[,i] <- k_list[[i]]$x

    bhat_mat[id, i] <- k_list[[i]]$bhat
    bhat_xval_mat[idx, i] <- k_list[[i]]$bhat_xval

    res_mat[id, i] <- k_list[[i]]$res
    res_xval_mat[idx, i] <- k_list[[i]]$res_xval

    stats_mat[1,i] <- k_list[[i]]$r
    stats_mat[2,i] <- k_list[[i]]$rmse
    stats_mat[3,i] <- k_list[[i]]$r_xval
    stats_mat[4,i] <- k_list[[i]]$rmse_xval
  }

  # ======================================= #
  # Calculate means and confidence interval
  # ======================================= #

  if(is.null(alpha)) {
    stop("Please define the significance level 'alpha'.")
  }

  alpha_symm <- alpha / 2

  # ====== #
  # S
  # ====== #
  S_df <- as.data.frame(
    t(apply(S_mat, 1, function(x) quantile(na.omit(x), probs = c(alpha_symm, 0.500, 1 - alpha_symm))))
  )

  names(S_df) <- c("Low", "Mu", "High")

  # ====== #
  # x
  # ====== #
  x_df <- as.data.frame(
    t(apply(x_mat, 1, function(x) quantile(na.omit(x), probs = c(alpha_symm, 0.500, 1 - alpha_symm))))
  )

  names(x_df) <- c("Low", "Mu", "High")

  rownames(x_df) <- colnames(A)

  # ============ #
  # bhat
  # ============ #

  bhat_df <- as.data.frame(
    t(apply(bhat_mat, 1,
            function(x) quantile(na.omit(x), probs = c(alpha_symm, 0.500, 1 - alpha_symm))))
  )

  bhat_xval_df <- as.data.frame(
    t(apply(bhat_xval_mat, 1,
            function(x) quantile(na.omit(x), probs = c(alpha_symm, 0.500, 1 - alpha_symm))))
  )

  names(bhat_df) <- names(bhat_xval_df) <- c("Low", "Mu", "High")

  # ============ #
  # Residuals
  # ============ #
  res_model <-   rbind(
    apply(res_xval_mat, 2, function(x) sd(x, na.rm = TRUE)),
    apply(res_xval_mat, 2, function(x) acf(na.omit(x), plot = FALSE)$acf[c(2, 3)])
  )

  rownames(res_model) <- c('sigma', 'ar1', 'ar2')


  # ================ #
  # error statistics
  # ================ #
  stats_df <- as.data.frame(
    t(apply(stats_mat, 1, function(x) quantile(x, probs = c(alpha_symm, 0.500, 1 - alpha_symm))))
    )

  colnames(stats_df) <- c("Low", "Mu", "High")

  # =============================== @
  # RETURN
  # =============================== @

  SMITE_list <- list(
    "S" = S_df,
    "x" = x_df,
    "x_bootstrap" = x_mat,
    "res_model" = res_model,
    "bhat" = bhat_df,
    "bhat_xval" = bhat_xval_df,
    "e_stats" = stats_df
  )

  return(SMITE_list)

}

