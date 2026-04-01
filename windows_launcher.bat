@echo off
start "" /b cmd /c "timeout /t 3 /nobreak >nul && start "" http://localhost:5288/" 

docker run --rm -p 5288:5288 ^
  -e SHINY_PORT=5288 ^
  -e SHINY_ROOT_PATH=/browse ^
  -e SHINY_ROOT_NAME=home ^
  --mount type=bind,source="C:\Users",target=/browse ^
  --mount type=bind,source="C:\shiny_out",target=/out ^
  qtea1/bulktools:latest

pause