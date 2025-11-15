#' SMITE Calibration
#'
#' Create a model (x) based on the 'forward' coral variable matrix (A) and the reconstruction target (b).
#' @param A The 'forward' coral variable matrix (t x p; time in rows, coral parameters in columns).
#' @param b The reconstruction target (t x 1).
#' @param Ae A t x p matrix containing the errors for the forward matrix.
#' @param be Error estimates for the reconstruction target. It can either be one value or a vector of equal length to the reconstruction target (t x 1).
#' @param it The number of bootstrap Monte Carlo iterations (recommended at 10,000).
#' @param acc Autocorrelation Coefficient; if noise is specified as 'red', it describes the degree of autocorrelation in the error term.
#' @param eigenclean Describes how (if at all) the singular values should be truncated to 'clean' the inverse solution. If the value is between 0 and 1, it will remove singular values based on the cumulative variance explained. If the value is greater than 1, it will return that many singular values (from the highest order).
#' @param alpha The significance level for the confidence interval (i.e., 0.05 = 95-percent confidence).
#' @param xval Cross-validation window size.
#' @param weights TRUE or FALSE on whether you want to incorporate uncertainty (Ae and be) into calibration. Default is TRUE.
#' @return S - Confidence intervals on the singular values of the forward matrix,
#' @return x - SMITE model parameters, or loadings, for each column of the forward matrix.
#' @return bhat - Predicted target anomalies from in-sample blocks
#' @return bhat.xval - predicted target anomalies from out-of-sample blocks
#' @return recon - Predicted absolute target values from out-of-sample blocks
#' @return residuals - Residuals from the calibration (bhat - b)
#' @return residuals.xval - Only cross-validated residuals from the calibration
#' @return e - Relevant error metrics. It enumerates the Standard Error of Prediction (SEP), Root-Mean-Square-Error (RMSE), and the correlation coefficient (r).
#' @return e.xval - Relevant error metrics for the cross-validated data points specifically.
#' @return e.stats - Confidence intervals for R and RMSE (non-cross-validated and cross-validated). The size of the confidence interval is based on alpha.
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


SMITE.calib <- function(A, b, Ae = NULL, be = NULL, it = 10000, acc = NULL,
                        eigenclean = NULL, alpha = 0.05, xval = NULL, weights = TRUE) {


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
    print("No window size for cross-validation (xval) specified. Using minimum window size of 1.")
    xval <- 1
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

    # Test
    # plot(bpx, col = 'blue', type = 'l')
    # lines(Apx[,2] * -1, col = 'red')

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

    if(weights) {

      # Compute row-wise total variance: be_i + sum(x_j^2 * Ae_ij)
      tvar <- rep(0, nrow(Apx))
      for (i in 1:nrow(Apx)) {
        Avar <- sum((x^2) * (Aepx[i,]^2))  # variance contribution from A
        tvar[i] <- bepx[i]^2 + Avar     # total variance per row (b + A)
      }

      # Weight matrix: inverse of total variance
      W <- diag(sqrt(1 / tvar))  # sqrt(W) for premultiplication


      # =============================== #
      # Weighted: SVD of Ap
      # =============================== #
      Apxw <- W %*% Apx
      bpxw <- W %*% bpx

      G <- svd(Apxw)

      U <- G$u
      V <- G$v
      S <- G$d

      # =============================== #
      # Weighted: Eigenclean
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
      # Weighted: Moore-Penrose Pseudoinverse
      # =============================== #

      x <- V %*% diag(1/S) %*% t(U) %*% bpxw


    }

    bhat <- c(((Apx %*% x) * sd(b)) + mean(b))
    bhat_xval <- c(((Apv %*% x) * sd(b)) + mean(b))

    # Test
    # plot(x = c(1:length(b))[-rr], y = bpx, col = 'blue', type = 'l')
    # lines(x = rr, y = bhat_xval, col = 'red')
    # lines(x = rr, y = bpv, col = 'purple')
    # lines(x = c(1:length(b))[-rr], y = bhat, col = 'red')

    # =============================== #
    # Store S, m, & bhat, & residuals
    # =============================== #

    # Statistics (use universal sd(b) for scaling)
    r <- cor.test(bhat, bpx)$estimate
    rmse <- sqrt(mean((bpx - bhat)^2)) * sd(b)

    r_xval <- cor.test(bpv, bhat_xval)$estimate
    rmse_xval <- sqrt(mean((bpv - bhat_xval)^2)) * sd(b)

    k_list[[k]] <- list(
      "S" = S,
      "x" = x,
      "bhat" = bhat,
      "bhat_xval" = bhat_xval,
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
  # bhat & recon
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

  # Test
  # plot(b, type = 'l', col = 'blue')
  # lines(bhat_xval_df$Mu, col = 'red')

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
    "bhat" = bhat_df,
    "bhat_xval" = bhat_xval_df,
    "e_stats" = stats_df
  )

  return(SMITE_list)

}

