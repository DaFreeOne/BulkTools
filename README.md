
WITHOUT DOCKER:
Run the app :
1) Run in R terminal :
> shiny::runApp(host = "0.0.0.0", port = 5288, launch.browser = FALSE)


2) Open in browser :
> http://localhost:5288



WITH DOCKER : 
1) Build the Docker : 
Move yourself so that you are in the same directory as the DockerFile
Run this command in the terminal :
> docker build -t shiny_docker .

2) To run the Docker, run this in the terminal :
> docker run --rm \
       -p 5288:5288 \
       -e SHINY_PORT=5288 \
       -e SHINY_ROOT_PATH=/browse \
       -e SHINY_ROOT_NAME=home \
       --mount type=bind,src="$HOME",target=/browse \
       shiny_docker