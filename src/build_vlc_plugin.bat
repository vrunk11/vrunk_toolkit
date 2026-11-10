@echo off
REM Compile le plugin VLC "tscale" : sortie video qui met l'image a l'echelle sur le GPU
REM (OpenCL, memes algorithmes que temporal_scaling.exe) puis l'affiche.
REM
REM   build_vlc_plugin.bat  [dossier_sdk_vlc]
REM
REM Prerequis : gcc (MinGW-w64) dans le PATH, et le VLC 3.0.x installe dans
REM C:\Program Files\VideoLAN\VLC (on y prend libvlccore.dll pour l'edition de liens).
REM Le SDK de VLC (en-tetes) n'est pas livre avec l'installateur : sans argument, ce script
REM le telecharge (archive officielle vlc-3.0.24-win64.7z, ~39 Mo, download.videolan.org)
REM et n'en extrait que le dossier "sdk" dans vlc_plugin\sdk. Les plugins VLC 3.0.x sont
REM compatibles entre eux : change VLCVER ci-dessous si ton VLC est une autre 3.0.x.
REM
REM Resultat : vlc_plugin\video_output\libtscale_plugin.dll
REM Pour l'utiliser : vlc_tscale.bat (sans droits admin) ou install_vlc_plugin.bat.
setlocal
cd /d "%~dp0"

set "VLCVER=3.0.24"
set "VLCDIR=C:\Program Files\VideoLAN\VLC"
set "ROOT=%~dp0vlc_plugin"
set "SDK=%~1"

where gcc >nul 2>nul
if errorlevel 1 ( echo Erreur : gcc introuvable dans le PATH. & exit /b 1 )
if not exist "%VLCDIR%\libvlccore.dll" ( echo Erreur : %VLCDIR%\libvlccore.dll introuvable ^(VLC 3.0.x installe ?^). & exit /b 1 )

if "%SDK%"=="" set "SDK=%ROOT%\sdk"
if not exist "%SDK%\include\vlc\plugins\vlc_common.h" (
    echo SDK VLC absent : telechargement de vlc-%VLCVER%-win64.7z ...
    if not exist "%ROOT%" mkdir "%ROOT%"
    curl -L --fail -o "%ROOT%\vlc-%VLCVER%-win64.7z" "https://download.videolan.org/pub/videolan/vlc/%VLCVER%/win64/vlc-%VLCVER%-win64.7z"
    if errorlevel 1 ( echo ECHEC du telechargement. & exit /b 1 )
    if not exist "%ROOT%\dl" mkdir "%ROOT%\dl"
    tar -xf "%ROOT%\vlc-%VLCVER%-win64.7z" -C "%ROOT%\dl" "vlc-%VLCVER%/sdk"
    if errorlevel 1 ( echo ECHEC de l'extraction du SDK. & exit /b 1 )
    if exist "%ROOT%\sdk" rmdir /s /q "%ROOT%\sdk"
    move "%ROOT%\dl\vlc-%VLCVER%\sdk" "%ROOT%\sdk" >nul
    rmdir /s /q "%ROOT%\dl"
    del "%ROOT%\vlc-%VLCVER%-win64.7z"
    set "SDK=%ROOT%\sdk"
)
if not exist "%SDK%\include\vlc\plugins\vlc_common.h" ( echo Erreur : SDK VLC introuvable dans %SDK% & exit /b 1 )

if not exist "%ROOT%\obj" mkdir "%ROOT%\obj"
if not exist "%ROOT%\video_output" mkdir "%ROOT%\video_output"

REM bibliotheque d'import de libvlccore.dll (celle du VLC installe)
copy /y "%VLCDIR%\libvlccore.dll" "%ROOT%\obj\libvlccore.dll" >nul
pushd "%ROOT%\obj"
gendef libvlccore.dll >nul 2>nul
dlltool -d libvlccore.def -l libvlccore.dll.a -D libvlccore.dll
popd

echo Compilation du moteur (OpenCL)...
gcc -c -O3 -w -I"%~dp0." -o "%ROOT%\obj\ts_core.o" "%~dp0ts_core.c"
if errorlevel 1 ( echo ECHEC ts_core & exit /b 1 )

echo Compilation du plugin VLC...
gcc -c -O2 -w -D_WIN32_WINNT=0x0601 -include winsock2.h -Dpoll=WSAPoll -D__PLUGIN__ -DMODULE_STRING=\"tscale\" -I"%SDK%\include\vlc\plugins" -I"%SDK%\include" -I"%~dp0." -o "%ROOT%\obj\vlc_tscale.o" "%~dp0vlc_tscale.c"
if errorlevel 1 ( echo ECHEC vlc_tscale & exit /b 1 )

gcc -shared -o "%ROOT%\video_output\libtscale_plugin.dll" "%ROOT%\obj\vlc_tscale.o" "%ROOT%\obj\ts_core.o" -L"%ROOT%\obj" -lvlccore -lopengl32 -lgdi32 -luser32 -lwinmm -Wl,-Bstatic -lpthread -Wl,-Bdynamic -static-libgcc
if errorlevel 1 ( echo ECHEC edition de liens & exit /b 1 )

echo OK : %ROOT%\video_output\libtscale_plugin.dll
endlocal
