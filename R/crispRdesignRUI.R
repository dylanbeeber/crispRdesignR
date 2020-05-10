#' @export

##crispRdesignRUI
crispRdesignRUI <- function() {
  app_path <- system.file("apps", "RunShiny.R", package = "crispRdesignR")
  shiny::runApp(app_path)
}
