@echo off
rem Doble clic para abrir el editor de la simulacion en el navegador.
cd /d "%~dp0"
if not exist ".venv\Scripts\python.exe" (
  echo No se encontro .venv. Sigue los pasos de instalacion del README.
  pause
  exit /b 1
)
".venv\Scripts\python.exe" -m mcplan editor
if errorlevel 1 pause
