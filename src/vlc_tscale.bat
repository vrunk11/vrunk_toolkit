@echo off
REM Lance VLC avec la sortie video "tscale" : mise a l'echelle sur le GPU (OpenCL) avec les
REM algorithmes de temporal_scaling, puis affichage. VLC reste le lecteur (interface, codecs,
REM sous-titres, audio, synchronisation). Pas besoin de droits admin : le plugin est charge
REM depuis vlc_plugin\ via la variable VLC_PLUGIN_PATH (compile-le d'abord : build_vlc_plugin.bat).
REM
REM Usage : vlc_tscale.bat [fichier]     (ou glisser le fichier sur ce script)
REM
REM Desentrelacement : utilise celui de VLC. Pour du contenu entrelace, ajoute par exemple
REM    --deinterlace=1 --deinterlace-mode=yadif2x
REM (yadif2x double la cadence : mets alors "frame ratio" a 1 ou plus selon ta cible).
REM Choisis l'algorithme/le GPU ici ou dans VLC : Outils > Preferences > Tout > Video > Sortie.

SET /p out_h= output vertical resolution (ex 1080 ou 2160 ; 0 = hauteur de l'ecran) :
SET /p ratio= frame ratio (images affichees par image decodee, 1 = normal) :
SET /p mode= mode (0~7, 5 = kaiser) :
SET /p tmp_mode= temporal mode (0~3, 0 = duplication ; 1+ retarde l'image d'une trame) :

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
  --tscale-height=%out_h% --tscale-ratio=%ratio% --tscale-mode=%mode% --tscale-tmp-mode=%tmp_mode% %*

pause
