#' BMDA 1B Error Composite Data
#'
#' A dataset from Porites astreoides coral 1B (Hughes et al., 2024) containing error estimates from age modeled geochemical variables along with in situ climate variables.

#
#' @format A data frame with 33 rows and 11 columns:
#' \describe{
#'   \item{Date}{Monthly (decimal dates), ranging from June 2010 to February 2013.}
#'   \item{Temp}{in situ temperature measured from Hogs Reef Buoy.}
#'   \item{pH}{in situ pH measured from Hogs Reef Buoy.}
#'   \item{pHcf}{pH of the calcifying fluid, estimated using boron systematics (e.g., D'Olivo et al., 2017) and in situ temperature.}
#'   \item{BCa}{Boron-to-calcium ratios (umol/mol).}
#'   \item{SrCa}{Strontium-to-calcium ratios (mmol/mol).}
#'   \item{MgCa}{Magnesium-to-calcium ratios (mmol/mol).}
#'   \item{UCa}{Uranium-to-calcium ratios (umol/mol).}
#'   \item{LiCa}{Lithium-to-calcium ratios (umol/mol).}
#'   \item{LiMg}{Lithium-to-magnesium ratios (mmol/mol), calculated by normalizing Li/Ca to Mg/Ca.}
#'   \item{d11B}{delta-11 Boron isotope ratios (permil).}
#' }
#'
#' @source Hughes et al., (2024)
"bmda_1b_ecomp"
