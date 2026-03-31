
WITHOUT DOCKER:
Run the app :
1) Run in R terminal :
R > shiny::runApp(host = "0.0.0.0", port = 5288, launch.browser = FALSE)


2) Open in browser :
URL > http://localhost:5288



WITH DOCKER : 
1) Build the Docker : 
Move yourself so that you are in the same directory as the DockerFile
sh > docker build -t shiny_docker -f DockerFile .

2) Run the Docker :
