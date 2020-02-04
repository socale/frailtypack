#' Shiny application for modelisation and prediction of frailty models
#'
#' This function loads the shiny package and runs the application for modelisation and prediction of several frailty models using package frailtypack.
#' 
#' 
#' @usage runShiny()
#' @return No value returned.
#' @export
#' @references Rizopoulos D. (2016)
#'
#' @examples
#' \dontrun{
#' 
#' runShiny()
#' 
#' }
#' 
runShiny <- function () {
  if (requireNamespace("shiny", quietly = TRUE)) {
      shiny::runApp(system.file("shiny_frailtypack", package = "frailtypack"))
  } else {
    cat("\npackage 'shiny' is not installed...")
  }
}