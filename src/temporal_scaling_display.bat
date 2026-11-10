@echo off
REM Lecture directe, SANS pipe ni processus externe : temporal_scaling ouvre le fichier,
REM le decode (ffmpeg integre), le met a l'echelle en YUV 4:2:0 -> 4:4:4 sur le GPU (OpenCL),
REM convertit en RGB sur le GPU, l'affiche dans sa propre fenetre et joue le son.
REM Tout le reste est detecte tout seul : resolution, rapport d'aspect (DAR), cadence,
REM matrice de couleur, et contenu entrelace (desentrelacement automatique si tu ne choisis rien).
REM Usage : glisser le fichier video sur ce script (ou : temporal_scaling_display.bat fichier)
REM Touches dans la fenetre : Echap = quitter, F ou F11 = plein ecran, Espace = pause.
REM Ce script doit etre a cote de temporal_scaling.exe (compile avec ffmpeg integre :
REM voir compilation.bat).

SET /p out_h= output vertical resolution (ex 1080 ou 2160) :
SET /p ratio= frame ratio :
SET /p mode= mode (0~7, 5 = kaiser) :

echo.
"%~dp0temporal_scaling.exe" --deint-list
echo.
SET /p deint= deinterlacer ALGO:CADENCE (ex 1:2 ; Entree = automatique, ou une chaine de filtres ffmpeg) :

set "deintopt="
if not "%deint%"=="" set deintopt=--deint "%deint%"

"%~dp0temporal_scaling.exe" --source "%~1" --out-height %out_h% --fps-ratio %ratio% --mode %mode% %deintopt% --opencl --display

pause
