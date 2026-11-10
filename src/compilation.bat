@echo off
REM Compile temporal_scaling (gcc / MinGW-w64) en executable STATIQUE : aucune DLL a
REM distribuer (ni ffmpeg, ni MinGW). Seules des DLL Windows systeme sont utilisees.
REM
REM  - OpenCL : OpenCL.dll est chargee a l'execution (--opencl), rien a lier.
REM  - Fenetre d'affichage (--display) : opengl32 gdi32 user32 winmm (fournis par Windows).
REM  - ffmpeg integre (--source / --vf) : necessite ffmpeg STATIQUE de vcpkg. Une fois :
REM        vcpkg install ffmpeg:x64-mingw-static
REM    Le dossier est cherche dans %VCPKG_ROOT% puis C:\vcpkg. S'il est absent, le programme
REM    est compile sans --source (les pipes -i et --input-cmd continuent de fonctionner).
setlocal
cd /d "%~dp0"

where gcc >nul 2>nul
if errorlevel 1 (
    echo Erreur : gcc introuvable dans le PATH.
    exit /b 1
)

set "VCPKG=%VCPKG_ROOT%"
if "%VCPKG%"=="" set "VCPKG=C:\vcpkg"
set "FF=%VCPKG%\installed\x64-mingw-static"

set "FFOPT="
set "FFLIBS="
if exist "%FF%\include\libavformat\avformat.h" (
    set "FFOPT=-DTS_FFMPEG -I"%FF%\include""
    set "FFLIBS=-L"%FF%\lib" -lavformat -lavfilter -lavcodec -lswscale -lswresample -lavutil -latomic -lsecur32 -lncrypt -lcrypt32 -lws2_32 -lole32 -lbcrypt"
    echo ffmpeg statique trouve : --source et --vf disponibles.
) else (
    echo ffmpeg statique introuvable dans %FF%
    echo Compilation SANS --source. Pour l'activer : vcpkg install ffmpeg:x64-mingw-static
)

echo Compilation de temporal_scaling...
gcc -O3 -w -static -s -Wl,--gc-sections %FFOPT% -o temporal_scaling.exe temporal_scaling.c %FFLIBS% -lm -lpthread -lopengl32 -lgdi32 -luser32 -lwinmm
if errorlevel 1 (
    echo ECHEC de la compilation.
    exit /b 1
)

echo OK : %~dp0temporal_scaling.exe
endlocal

pause
