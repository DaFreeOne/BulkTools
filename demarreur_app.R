port <- as.integer(Sys.getenv("SHINY_PORT", "5288"))
host <- Sys.getenv("SHINY_HOST", "0.0.0.0")
app_dir <- Sys.getenv("SHINY_APP_DIR", "/app")

shiny::runApp(
  appDir = app_dir,
  host = host,
  port = port,
  launch.browser = FALSE
)
