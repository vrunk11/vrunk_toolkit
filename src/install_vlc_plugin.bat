@echo off
REM Installe le plugin VLC "tscale" dans le dossier plugins de VLC. La copie dans
REM "Program Files" demande les droits administrateur : si le script n'est pas lance en
REM administrateur, il se relance tout seul (une fenetre de confirmation Windows s'affiche).
REM
REM Ensuite : VLC > Outils > Preferences > Video > Sortie > "Temporal scaling (OpenCL GPU)"
REM (ses reglages : Preferences > Tout > Video > Sortie > tscale). Pour desinstaller :
REM supprimer le fichier copie.
setlocal
cd /d "%~dp0"
set "SRC=%~dp0vlc_plugin\video_output\libtscale_plugin.dll"
set "VLCDIR=C:\Program Files\VideoLAN\VLC"
set "DST=%VLCDIR%\plugins\video_output"

net session >nul 2>&1
if errorlevel 1 (
    echo Droits administrateur requis : relance du script en administrateur...
    powershell -NoProfile -Command "Start-Process -FilePath '%~f0' -Verb RunAs"
    exit /b 0
)

if not exist "%SRC%" ( echo Plugin absent : lance d'abord build_vlc_plugin.bat & pause & exit /b 1 )
if not exist "%DST%" ( echo Dossier VLC introuvable : %DST% & pause & exit /b 1 )

REM VLC doit etre ferme : il verrouille ses plugins pendant qu'il tourne
tasklist /fi "imagename eq vlc.exe" | find /i "vlc.exe" >nul
if not errorlevel 1 (
    echo Ferme VLC avant l'installation, puis relance ce script.
    pause
    exit /b 1
)

copy /y "%SRC%" "%DST%\"
if errorlevel 1 (
    echo ECHEC de la copie vers %DST%
    pause
    exit /b 1
)

REM regenere le cache des plugins de VLC pour qu'il voie le nouveau module
"%VLCDIR%\vlc-cache-gen.exe" "%VLCDIR%\plugins"
echo.
echo OK : plugin installe dans %DST%
echo Lance VLC, puis Outils ^> Preferences ^> Video ^> Sortie ^> "Temporal scaling (OpenCL GPU)".
endlocal
pause
