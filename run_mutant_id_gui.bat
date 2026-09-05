@echo off
setlocal
cd /d "%~dp0"

python -c "import numpy, matplotlib, PIL, Bio" >nul 2>&1
if errorlevel 1 (
    echo Missing Python packages.
    echo Install them with:
    echo   python -m pip install -r requirements-mutant-id.txt
    pause
    exit /b 1
)

python mutant_id_gui.py
if errorlevel 1 (
    echo.
    echo The mutant ID GUI exited with an error.
)
pause
