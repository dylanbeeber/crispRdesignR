#' @export

##crispRdesignRUI
crispRdesignRUI <- function() {
  requireNamespace("gbm", quietly = TRUE)
  app_path <- system.file("apps", "RunShiny.R", package = "crispRdesignR")
  shiny::runApp(app_path)
}
