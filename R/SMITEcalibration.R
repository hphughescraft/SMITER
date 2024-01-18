#' SMITE Calibration
#'
#' Create a model (m) based on the 'forward' coral variable matrix (A) and the reconstruction target (b).
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
#' @return m - SMITE model parameters, or loadings, for each column of the forward matrix.
#' @return bhat - Predicted target anomalies
#' @return recon - Predicted absolute target values
#' @return e - Relevant error metrics. It enumerates (half) the distance of the 95-percent confidence interval (the average distance from the mean to the upper and lower bounds) for the singular values and the model parameters. It also enumerates the Standard Error of Prediction (SEP), Root-Mean-Square-Error (RMSE), and the correlation coefficient (r).
#' @export

SMITE.calib <- function(A, b, Ae = NULL, be = NULL, it = 10000, noise = "white", acc = NULL,
                        eigenclean = NULL, alpha = 0.05) {
  
  # =============================== #
  # Convert A & b to matrices
  # =============================== #
  
  A <- as.matrix(A)
  b <- as.matrix(b)
  
  # =============================== #
  # Perturb the A matrix (Ap) & b (bp)
  # =============================== #
  
  if(is.null(Ae) == FALSE) {
    
    # =============================== #
    # Wrangle
    # =============================== #
    Ae <- as.matrix(Ae)
    
    if(any(dim(Ae) != dim(A))) {
      stop("Error matrix (Ae) must be same dimensions as the forward matrix (A).")
    }
    
    be <- as.matrix(be)
    
    if(any(dim(be) != dim(b))) {
      stop("Error estimates of the reconstruction target (be) must be either one value or the same dimensions as the reconstruction target (b).")
    }
    
    if(is.null(noise)) {
      stop("Please specify whether the noise is 'red' (autocorrelated) or 'white' (Gaussian).")
    }
    
    if(noise == 'red' & is.null(acc)) {
      stop("Please specify the Autocorrelation Coefficient for the red noise.")
    }
    
    # =============================== #
    # Begin bootstrap
    # =============================== #
    k.list <- list()
    
    for(k in 1:it) {
      
      # Preallocate a Gp & dp matrix #
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
      
      
      # =============================== #
      # SVD of Gp
      # =============================== #
      Ap.norm <- apply(Ap, MARGIN = 2, FUN = function(x) (x - mean(x)) /  sd(x))
      bp.norm <- (bp - mean(bp)) / sd(bp)
      
      G <- svd(Ap.norm)
      
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
      
      m <- V %*% diag(1/S) %*% t(U) %*% bp.norm
      bhat <- Ap.norm %*% m
      
      # =============================== #
      # Store S, m, & bhat
      # =============================== #
      
      k.list[[k]] <- list(
        "S" = S,
        "m" = m,
        "bhat" = bhat
      )
    } # Bootstrap
    
    # =============================== #
    # Extract distributions
    # =============================== #
    
    S.mat <- matrix(NA, nrow = length(S), ncol = it)
    m.mat <- matrix(NA, nrow = nrow(m), ncol = it)
    bhat.mat <- matrix(NA, nrow = nrow(bhat), ncol = it)
    
    for(i in 1:it) {
      S.mat[,i] <- k.list[[i]]$S
      m.mat[,i] <- k.list[[i]]$m
      bhat.mat[,i] <- k.list[[i]]$bhat
    }
    
    # ======================================= #
    # Calculate means and confidence interval
    # ======================================= #
    
    if(is.null(alpha)) {
      stop("Please define the significance level 'alpha'.")
    }
    
    alpha.symm <- alpha / 2
    cint.index <- seq(floor(alpha.symm * it), floor((1 - alpha.symm) * it), 1)
    
    # ====== #
    # S
    # ====== #
    S.mu <- apply(S.mat, 1, function(x) mean(x))
    S.sd <- apply(S.mat, 1, function(x) sd(x))
    S.cint <- apply(S.mat, 1, function(x) x[order(x)][cint.index])
    S.low <- S.cint[1,]
    S.high <- S.cint[nrow(S.cint),]
    
    S.df <- data.frame(
      "S.low" = S.low,
      "S.mu" = S.mu,
      "S.high" = S.high
    )
    
    # ====== #
    # m
    # ====== #
    m.mu <- apply(m.mat, 1, function(x) mean(x))
    m.sd <- apply(m.mat, 1, function(x) sd(x))
    m.cint <- apply(m.mat, 1, function(x) x[order(x)][cint.index])
    m.low <- m.cint[1,]
    m.high <- m.cint[nrow(m.cint),]
    
    m.df <- data.frame(
      "m.low" = m.low,
      "m.mu" = m.mu,
      "m.high" = m.high
    )
    
    rownames(m.df) <- colnames(G)
    
    # ====== #
    # bhat
    # ====== #
    bhat.mu <- apply(bhat.mat, 1, function(x) mean(x))
    bhat.sd <- apply(bhat.mat, 1, function(x) sd(x))
    bhat.cint <- apply(bhat.mat, 1, function(x) x[order(x)][cint.index])
    bhat.low <- bhat.cint[1,]
    bhat.high <- bhat.cint[nrow(bhat.cint),]
    
    bhat.df <- data.frame(
      "bhat.low" = bhat.low,
      "bhat.mu" = bhat.mu,
      "bhat.high" = bhat.high
    )
    
    recon.df <- (bhat.df * sd(b)) + mean(b)
    # =============================== @
    # RETURN
    # =============================== @
    
    SMITE.list <- list(
      "S" = S.df,
      "m" = m.df,
      "bhat" = bhat.df,
      "recon" = recon.df,
      "e" = c(
        "S" = mean((S.df$S.high - S.df$S.low) / 2),
        "m" = mean((m.df$m.high - m.df$m.low) / 2),
        "bhat" = mean(c(bhat.df$bhat.high - bhat.df$bhat.mu, bhat.df$bhat.mu - bhat.df$bhat.low)) / 1.96,
        "SEP" = (mean(c(bhat.df$bhat.high - bhat.df$bhat.mu, bhat.df$bhat.mu - bhat.df$bhat.low)) / 1.96) * sd(b),
        "RMSE" = (sum(sqrt((b - recon.df$bhat.mu)^2))) / (nrow(recon.df) - 2),
        "r" = cor.test(b, bhat.df$bhat.mu)$estimate
      )
    )
    
    return(SMITE.list)
    
  } else if(is.null(Ae)) {
    stop("Must define errors on the forward matrix (A).")
  }
}
