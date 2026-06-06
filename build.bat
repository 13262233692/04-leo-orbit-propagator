@echo off
echo ========================================
echo LEO Orbit Propagator - Build Script
echo ========================================

where python >nul 2>nul
if %errorlevel% neq 0 (
    echo Error: Python not found in PATH
    exit /b 1
)

echo [1/5] Checking for CMake...
where cmake >nul 2>nul
if %errorlevel% neq 0 (
    echo Error: CMake not found in PATH
    exit /b 1
)
echo CMake found.

echo.
echo [2/5] Installing Python dependencies...
python -m pip install --upgrade pip
python -m pip install -r requirements.txt
if %errorlevel% neq 0 (
    echo Error: Failed to install Python dependencies
    exit /b 1
)

echo.
echo [3/5] Checking for pybind11...
if not exist "third_party\pybind11\CMakeLists.txt" (
    echo Downloading pybind11...
    if not exist "third_party" mkdir third_party
    cd third_party
    
    where git >nul 2>nul
    if %errorlevel% equ 0 (
        git clone --depth 1 --branch v2.12.0 https://github.com/pybind/pybind11.git
    ) else (
        echo Warning: git not found, please install pybind11 manually
        echo into third_party/pybind11
        cd ..
        exit /b 1
    )
    cd ..
)
echo pybind11 found.

echo.
echo [4/5] Building C++ extension module...
python setup.py build_ext --inplace
if %errorlevel% neq 0 (
    echo Error: Build failed
    exit /b 1
)

echo.
echo [5/5] Installing package in development mode...
python -m pip install -e .

echo.
echo ========================================
echo Build completed successfully!
echo ========================================
