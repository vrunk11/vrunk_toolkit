@echo off
REM Installe le plugin VLC "tscale" dans le dossier plugins de VLC (necessite les droits admin :
REM clic droit > Executer en tant qu'administrateur). Apres cela, la sortie video "Temporal scaling
REM (OpenCL GPU)" est disponible dans VLC : Outils > Preferences > Video > Sortie, et ses reglages
REM dans Preferences > Tout > Video > Sortie > tscale. Lancer VLC une fois avec --reset-plugins-cache
REM si VLC ne la voit pas. Pour desinstaller : supprimer le fichier copie.
setlocal
cd /d "%~dp0"
set "SRC=%~dp0vlc_plugin\video_output\libtscale_plugin.dll"
set "DST=C:\Program Files\VideoLAN\VLC\plugins\video_output"

if not exist "%SRC%" ( echo Plugin absent : lance d'abord build_vlc_plugin.bat & exit /b 1 )
if not exist "%DST%" ( echo Dossier VLC introuvable : %DST% & exit /b 1 )

copy /y "%SRC%" "%DST%\" >nul
if errorlevel 1 (
    echo ECHEC de la copie : relance ce script en administrateur.
    exit /b 1
)
"C:\Program Files\VideoLAN\VLC\vlc-cache-gen.exe" "C:\Program Files\VideoLAN\VLC\plugins" >nul 2>nul
echo OK : plugin installe dans %DST%
endlocal
pause
