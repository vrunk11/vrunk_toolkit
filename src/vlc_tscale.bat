@echo off
REM Lance VLC avec la sortie video "tscale" : mise a l'echelle sur le GPU (OpenCL) avec les
REM algorithmes de temporal_scaling, puis affichage. VLC reste le lecteur (interface, codecs,
REM sous-titres, audio, synchronisation). Pas besoin de droits admin : le plugin est charge
REM depuis vlc_plugin\ via la variable VLC_PLUGIN_PATH (compile-le d'abord : build_vlc_plugin.bat).
REM
REM Usage : vlc_tscale.bat [fichier]     (ou glisser le fichier sur ce script)
REM
REM Desentrelacement : celui de VLC (applique avant notre mise a l'echelle ; aussi reglable dans
REM VLC : menu Video > Desentrelacement). Les modes a double cadence (yadif2x, bwdif, bob...)
REM produisent 2 images par trame : le plugin suit la cadence recue, "frame ratio" s'y ajoute.
REM Choisis l'algorithme/le GPU ici ou dans VLC : Outils > Preferences > Tout > Video > Sortie.

SET /p out_h= output vertical resolution (Entree = automatique : suit la fenetre ; ou 1080, 2160... ; 0 = hauteur de l'ecran) :
SET /p ratio= frame ratio (Entree/0 = automatique : vise la cadence "fps" par pulldown ; ou 1 = normal, 2, 3... : fixe, ignore fps) :
set "fps=0"
if "%ratio%"=="" set "ratio=0"
if "%ratio%"=="0" SET /p fps= target fps in Hz (Entree = frequence de l'ecran ; ou 60, 59.94, 120...) :
if "%fps%"=="" set "fps=0"
SET /p mode= mode (0~7, 5 = kaiser) :
SET /p tmp_mode= temporal mode (0~3, 0 = duplication ; 1+ retarde l'image d'une trame) :
SET /p deint= deinterlacing (Entree = automatique selon le fichier ; 0 = jamais ; ou yadif, yadif2x, bwdif, bob, linear, x, mean, blend, phosphor, ivtc) :
set "deint_args="
if not "%deint%"=="" (
    if "%deint%"=="0" (
        set "deint_args=--deinterlace=0"
    ) else (
        set "deint_args=--deinterlace=1 --deinterlace-mode=%deint%"
    )
)
if "%out_h%"=="" set "out_h=-1"
if "%ratio%"=="" set "ratio=0"
if "%mode%"=="" set "mode=5"
if "%tmp_mode%"=="" set "tmp_mode=0"

set "VLC_PLUGIN_PATH=%~dp0vlc_plugin"
if not exist "%VLC_PLUGIN_PATH%\video_output\libtscale_plugin.dll" (
    echo Plugin introuvable : lance d'abord build_vlc_plugin.bat
    pause
    exit /b 1
)

REM --no-plugins-cache : VLC relit le dossier du plugin a chaque lancement
REM --avcodec-hw=none  : decodage logiciel -> les sous-titres/OSD restent affiches (le decodage
REM                      materiel renvoie des surfaces DX11 sur lesquelles VLC ne sait pas les melanger)
"C:\Program Files\VideoLAN\VLC\vlc.exe" --no-plugins-cache --avcodec-hw=none --vout=tscale ^
  --tscale-height=%out_h% --tscale-ratio=%ratio% --tscale-fps=%fps% --tscale-mode=%mode% --tscale-tmp-mode=%tmp_mode% %deint_args% %*

pause
