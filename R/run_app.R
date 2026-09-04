#' Launch the FRDA Transcriptomic Atlas app
#'
#' @param data_mode Data access mode. Defaults to "cloud". Use "local" for local access of data for Developers
#'
#' @return A shiny.appobj that runs the app
#' @export
run_app <- function() {

  shiny::shinyApp(
    ui = app_ui(),
    server = function(input, output, session) {
      app_server(input, output, session)
    }
  )
}
