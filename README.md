
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

2) To run the Docker : 
Edit root path "SHINY_ROOT_PATH" and "target" to the root from which you want to browse your local files.

Run this in the terminal :
FOR LINUX :
> docker run --rm -p 5288:5288 \
    -e SHINY_PORT=5288 \
    -e SHINY_ROOT_PATH=/browse \
    -e SHINY_ROOT_NAME=home \
    -v /home/quentin/data:/browse \
    qtea1/shiny_docker

    
FOR WINDOWS : 
> docker run --rm -p 5288:5288 `
    -e SHINY_PORT=5288 `
    -e SHINY_ROOT_PATH=/browse `
    -e SHINY_ROOT_NAME=home `
    --mount type=bind,source="C:\Users\Quentin\Documents",target=/browse `
    qtea1/shiny_docker

