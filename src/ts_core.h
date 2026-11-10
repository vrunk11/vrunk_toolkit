// ============================================================
// ts_core.h — API du scaler temporal_scaling pour être embarqué (plugin VLC).
//
// Pas de dépendance à VLC ici : ts_core.c compile le même moteur que
// temporal_scaling.c (cartes de mise à l'échelle, noyaux, OpenCL) et l'expose
// sous la forme de trois fonctions. L'entrée est du YUV planaire 8 bits
// (4:2:0 / 4:2:2 / 4:4:4), la sortie du RGB24 (une image par phase temporelle).
// ============================================================
#ifndef TS_CORE_H
#define TS_CORE_H

#include <stddef.h>
#include <stdint.h>

#ifdef __cplusplus
extern "C" {
#endif

typedef struct TsCtx TsCtx;

enum { TS_CHROMA_420 = 0, TS_CHROMA_422 = 1, TS_CHROMA_444 = 2 };

typedef struct {
	int    in_w, in_h;          // taille visible de l'image source (paire pour 4:2:0 / 4:2:2)
	int    chroma;              // TS_CHROMA_*
	int    out_w, out_h;        // taille de sortie (doit être > à l'entrée)
	int    mode;                // 0..7, comme --mode
	int    ratio;               // multiplicateur de cadence (images produites par image reçue)
	int    tmp_mode;            // 0..3, comme --tmp-mode (1+ : une image de retard)
	double curve_base;          // --curve-base (> 1)
	double kaiser_beta;         // --kaiser-beta
	double aniso_strength;      // --aniso-strength
	int    ocl_device;          // index de --ocl-list, -1 = meilleur GPU
	int    yuv_709;             // matrice YUV->RGB : 1 = BT.709, 0 = BT.601
	int    yuv_full;            // 1 = plage pleine, 0 = limitée
} TsParams;

// Ouvre le scaler (initialise OpenCL, construit cartes et tables). NULL si échec ;
// err reçoit un message lisible.
TsCtx *ts_open(const TsParams *p, char *err, size_t errlen);

// Nom du périphérique OpenCL utilisé.
const char *ts_device_name(void);

// Nombre d'octets d'une image de sortie RGB24 (out_w*out_h*3).
size_t ts_frame_bytes(const TsCtx *c);

// Traite une image source : planes/pitches = plans Y,U,V (avec leur pas de ligne).
// out reçoit `ratio` images RGB24 consécutives (ratio * ts_frame_bytes octets).
// Avec tmp_mode > 0 les images rendues sont celles de l'image *précédente* (retard
// d'une image : le fondu a besoin de la suivante). Renvoie 0 si OK.
int ts_push(TsCtx *c, const uint8_t *const planes[3], const int pitches[3], uint8_t *out);

// Vérification rapide, sans rien initialiser de lourd : paramètres valides et un périphérique
// OpenCL utilisable est présent. Renvoie 0 si ts_open a de bonnes chances de réussir.
int ts_check(const TsParams *p, char *err, size_t errlen);

void ts_close(TsCtx *c);

#ifdef __cplusplus
}
#endif
#endif
