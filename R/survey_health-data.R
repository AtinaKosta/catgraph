#' Synthetic health survey data (categorical variables)
#'
#' A synthetic but realistically structured data frame of 600 respondents
#' with eight categorical health and demographic variables. The data are
#' generated with known, calibrated association strengths so that
#' \code{\link{catgraph}} produces an informative graph: strong associations
#' (smoking ~ lung disease), moderate ones (age group ~ exercise frequency,
#' exercise ~ BMI, BMI ~ diet quality, age group ~ health insurance), and
#' near-zero ones (sex ~ diet quality). Approximately 5\% of values are
#' missing completely at random (MCAR).
#'
#' @format A data frame with 600 rows and 8 columns:
#' \describe{
#'   \item{sex}{Factor with levels \code{female}, \code{male}.}
#'   \item{age_group}{Factor with levels \code{18-34}, \code{35-54},
#'     \code{55+}.}
#'   \item{smoking_status}{Factor with levels \code{never}, \code{former},
#'     \code{current}.}
#'   \item{lung_disease}{Factor with levels \code{no}, \code{yes}.}
#'   \item{exercise_freq}{Factor with levels \code{low}, \code{moderate},
#'     \code{high}.}
#'   \item{bmi_category}{Factor with levels \code{underweight},
#'     \code{normal}, \code{overweight}, \code{obese}.}
#'   \item{diet_quality}{Factor with levels \code{poor}, \code{fair},
#'     \code{good}.}
#'   \item{health_insurance}{Factor with levels \code{no}, \code{yes}.}
#' }
#'
#' @details
#' The data are entirely synthetic. No real individuals are represented.
#' The association structure was designed to serve as a pedagogical example
#' for \pkg{catgraph}: users can verify that the graph recovers the known
#' strong, moderate, and weak associations.
#'
#' Missingness was introduced completely at random (MCAR), independently
#' per column, at a 5\% rate. Pairwise deletion (as implemented throughout
#' \pkg{catgraph}) is therefore unbiased for this dataset.
#'
#' @source Synthetic data generated for package demonstration. No real
#'   individuals are represented.
#'
#' @examples
#' data(survey_health)
#' str(survey_health)
#' cg <- catgraph(survey_health)
#' summary(cg)
#' plot_heatmap(cg, title = "Health survey associations")
"survey_health"
