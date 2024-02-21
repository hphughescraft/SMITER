\name{SMITE.recon}
\alias{SMITE.recon}
\title{SMITE Reconstruction}
\description{Apply a SMITE model (m) to a 'forward' coral variable matrix (A) to produce predicted values for a reconstruction target (bhat).}
\usage{
SMITE.recon(
  A,
  m,
  Amu,
  Asd,
  bmu,
  bsd
)
}

\arguments{
  \item{A}{The 'forward' coral variable matrix (t x p; time in rows, coral parameters in columns).}
  \item{m}{A vector of SMITE model parameters (length p), ideally created from the SMITE.calib() function.}
  \item{Amu}{A vector of length p giving the mean values of each coral variable in the calibration dataset.}
  \item{Asd}{A vector of length p giving the standard deviations of each coral variable in the calibration dataset.}
  \item{bmu}{A single-element atomic vector giving the mean of the reconstruction target in the calibration dataset.}
  \item{bsd}{A single-element atomic vector giving the standard deviation of the reconstruction target in the calibration dataset.}
}

\value{
Returns a vector of predicted reconstruction target values (bhat). Minimum errors (i.e., not included age model errors) for bhat can be estimated from the SEP calculated in the SMITE.calib() function.
}

\references{
Hughes, H. P., Thompson, D. M., Foster, G. L., Lees, J., Surge, D. (in revision). “Synthetic and Practical Reconstructions of SST and pH Using the Novel Multiproxy SMITE Method.” PLOS ONE.
}

\author{Hunter P. Hughes}

\examples{
# Load data from Hughes et al. (in revision)
library(SMITER)
data(BMDA_1B_Comp)
data(BMDA_1B_EComp)

# Designate forward matrix (A), reconstruction target (b), and corresponding errors (Ae, be).

proxy_names <- names(bmda_1b_comp[,c(5:ncol(bmda_1b_comp))])

A <- bmda_1b_comp[,proxy_names]
Ae <- bmda_1b_ecomp[,proxy_names]
b <- bmda_1b_comp[,'Temp']
be <- rep(0.02, length(b))

# Execute SMITE Calibration with no truncation (i.e., all information retained).

result <- SMITE.calib(A = A, b = b, Ae = Ae, be = be, eigenclean = ncol(A))

# Extract SMITE model parameters #
m <- result$m$m.mu

# Reconstruction (should be nearly to result$recon$bhat.mu, only minor differences due to bootstrapping)
bhat <- SMITE.recon(
  A = A,
  m = m,
  Amu = colMeans(A),
  Asd = apply(A, 2, function(x) sd(x)),
  bmu = mean(b),
  bsd = sd(b)
)

result$recon
bhat

}

