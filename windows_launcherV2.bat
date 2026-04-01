@echo off
setlocal

set "BROWSE_DIR=%USERPROFILE%"
set "OUT_DIR=%USERPROFILE%\shiny_out"

if not exist "%OUT_DIR%" mkdir "%OUT_DIR%"

start "" /b cmd /c "timeout /t 3 /nobreak >nul && start "" http://localhost:5288/"

docker run --rm -p 5288:5288 ^
  -e SHINY_PORT=5288 ^
  -e SHINY_ROOT_PATH=/browse ^
  -e SHINY_ROOT_NAME=home ^
  --mount type=bind,source="%BROWSE_DIR%",target=/browse ^
  --mount type=bind,source="%OUT_DIR%",target=/out ^
  qtea1/bulktools:latest

pause