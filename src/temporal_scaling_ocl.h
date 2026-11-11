// ============================================================
// temporal_scaling_ocl.h — backend OpenCL pour temporal_scaling
//
// Objectifs :
//  - tourner sur un GPU Intel HD 520 (OpenCL 1.2/2.0) comme sur des GPU
//    AMD / NVIDIA récents : le code n'utilise que de l'OpenCL C 1.2
//    "full profile" sans extension (pas de fp64, pas d'images, pas
//    d'atomiques, pas de mémoire locale).
//  - ne demander aucun SDK : la bibliothèque OpenCL (OpenCL.dll,
//    libOpenCL.so.1) est chargée dynamiquement et les déclarations
//    nécessaires sont définies ici. Si elle est absente, le programme
//    continue sur le CPU.
//
// Ce fichier est inclus par temporal_scaling.c (pas de compilation séparée).
// Sous Linux, ajouter -ldl si la glibc est ancienne.
//
// Le GPU reproduit l'arithmétique entière du chemin CPU 2D (modes 0..7,
// RGB et planaire). Seul le chemin anisotrope des modes 6/7 utilise du
// flottant et peut différer d'une unité sur de rares pixels.
// ============================================================
#ifndef TEMPORAL_SCALING_OCL_H
#define TEMPORAL_SCALING_OCL_H

#include <stdint.h>
#include <stddef.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifdef _WIN32
	#include <windows.h>
	#define OCL_CALL __stdcall
#else
	#include <dlfcn.h>
	#define OCL_CALL
#endif

// ---------- déclarations OpenCL minimales (ABI stable) ----------
typedef int32_t  ocl_int;
typedef uint32_t ocl_uint;
typedef uint64_t ocl_ulong;
typedef struct _ocl_platform *ocl_platform;
typedef struct _ocl_device   *ocl_device;
typedef struct _ocl_context  *ocl_context;
typedef struct _ocl_queue    *ocl_queue;
typedef struct _ocl_mem      *ocl_mem;
typedef struct _ocl_program  *ocl_program;
typedef struct _ocl_kernel   *ocl_kernel;
typedef struct _ocl_event    *ocl_event;

#define OCL_DEVICE_TYPE_CPU          (1u << 1)
#define OCL_DEVICE_TYPE_GPU          (1u << 2)
#define OCL_DEVICE_TYPE_ALL          0xFFFFFFFFu
#define OCL_DEVICE_TYPE_INFO         0x1000
#define OCL_DEVICE_MAX_COMPUTE_UNITS 0x1002
#define OCL_DEVICE_MAX_CLOCK_FREQ    0x100C
#define OCL_DEVICE_GLOBAL_MEM_SIZE   0x101F
#define OCL_DEVICE_AVAILABLE         0x1027
#define OCL_DEVICE_COMPILER_AVAILABLE 0x1028
#define OCL_DEVICE_NAME              0x102B
#define OCL_PLATFORM_NAME            0x0902
#define OCL_CONTEXT_PLATFORM         0x1084
#define OCL_MEM_READ_WRITE           (1u << 0)
#define OCL_MEM_READ_ONLY            (1u << 2)
#define OCL_MEM_COPY_HOST_PTR        (1u << 5)
#define OCL_PROGRAM_BUILD_LOG        0x1183

typedef struct {
	ocl_int (OCL_CALL *GetPlatformIDs)(ocl_uint, ocl_platform *, ocl_uint *);
	ocl_int (OCL_CALL *GetPlatformInfo)(ocl_platform, ocl_uint, size_t, void *, size_t *);
	ocl_int (OCL_CALL *GetDeviceIDs)(ocl_platform, ocl_ulong, ocl_uint, ocl_device *, ocl_uint *);
	ocl_int (OCL_CALL *GetDeviceInfo)(ocl_device, ocl_uint, size_t, void *, size_t *);
	ocl_context (OCL_CALL *CreateContext)(const intptr_t *, ocl_uint, const ocl_device *, void *, void *, ocl_int *);
	ocl_queue (OCL_CALL *CreateCommandQueue)(ocl_context, ocl_device, ocl_ulong, ocl_int *);
	ocl_program (OCL_CALL *CreateProgramWithSource)(ocl_context, ocl_uint, const char **, const size_t *, ocl_int *);
	ocl_int (OCL_CALL *BuildProgram)(ocl_program, ocl_uint, const ocl_device *, const char *, void *, void *);
	ocl_int (OCL_CALL *GetProgramBuildInfo)(ocl_program, ocl_device, ocl_uint, size_t, void *, size_t *);
	ocl_kernel (OCL_CALL *CreateKernel)(ocl_program, const char *, ocl_int *);
	ocl_int (OCL_CALL *SetKernelArg)(ocl_kernel, ocl_uint, size_t, const void *);
	ocl_int (OCL_CALL *EnqueueNDRangeKernel)(ocl_queue, ocl_kernel, ocl_uint, const size_t *, const size_t *, const size_t *, ocl_uint, const ocl_event *, ocl_event *);
	ocl_mem (OCL_CALL *CreateBuffer)(ocl_context, ocl_ulong, size_t, void *, ocl_int *);
	ocl_int (OCL_CALL *EnqueueReadBuffer)(ocl_queue, ocl_mem, ocl_uint, size_t, size_t, void *, ocl_uint, const ocl_event *, ocl_event *);
	ocl_int (OCL_CALL *EnqueueWriteBuffer)(ocl_queue, ocl_mem, ocl_uint, size_t, size_t, const void *, ocl_uint, const ocl_event *, ocl_event *);
	ocl_int (OCL_CALL *Finish)(ocl_queue);
	ocl_int (OCL_CALL *Flush)(ocl_queue);
	ocl_int (OCL_CALL *WaitForEvents)(ocl_uint, const ocl_event *);
	ocl_int (OCL_CALL *ReleaseEvent)(ocl_event);
	ocl_int (OCL_CALL *ReleaseMemObject)(ocl_mem);
	ocl_int (OCL_CALL *ReleaseKernel)(ocl_kernel);
	ocl_int (OCL_CALL *ReleaseProgram)(ocl_program);
	ocl_int (OCL_CALL *ReleaseCommandQueue)(ocl_queue);
	ocl_int (OCL_CALL *ReleaseContext)(ocl_context);
} OclApi;

// ---------- géométrie fournie par le programme principal ----------
typedef struct {
	int    in_w, in_h, out_w, out_h;
	int    chan;        // 3 = RGB entrelacé, 1 = plan unique
	int    map_idx;     // 0 = carte luma/RGB, 1 = carte chroma
	size_t in_off;      // début du plan dans le buffer d'entrée (octets)
	size_t out_off;     // début du plan dans le buffer de sortie (octets)
	size_t out_step;    // octets d'une frame de sortie de ce plan
} OclPlane;

typedef struct {
	int mode, ratio, bps_in, bps_out, aniso_fp;
	int rgb_out;                            // 1 : convertir la sortie 4:4:4 en RGB24 sur le GPU
	int yuv_709, yuv_full;                  // matrice / plage pour cette conversion
	int nplanes;
	OclPlane planes[3];
	size_t in_len, out_len;                 // octets
	const BresenhamMap *maps[2];
	const int *lutA, *lutB, *lutC;          // radial LUT (B/C = NULL => A)
	int lut_size;
	const int *krad;                        // LUT kaiser (aniso) ou NULL
	int krad_size;
	const int *tmp_weights;
} OclSetup;

// ---------- état global ----------
typedef struct {
	int         loaded;
	void       *lib;
	OclApi      api;
	ocl_platform plats[16];
	ocl_device   devs[64];
	int          dev_plat[64];
	int          ndev;
	ocl_device   dev;
	ocl_platform plat;
	ocl_context  ctx;
	ocl_queue    q;
	ocl_program  prog;
	ocl_kernel   k_blend, k_grad, k_conf, k_dist, k_scale, k_scale4, k_yuv;
	ocl_mem      d_rgb;            // sortie RGB24 (si rgb_out)
	char         dev_name[512];

	OclSetup     s;
	ocl_mem      d_in, d_next, d_blend, d_out;
	ocl_mem      d_map[2][6];     // w_pos w_src w_flag h_pos h_src h_flag
	ocl_mem      d_grad[3], d_dist[3];
	ocl_mem      d_tab[2][3];      // tables compactes par carte (Y / chroma) et par LUT
	ocl_mem      d_dense[2];       // table dense 36 poids (noyau scale4_k)
	int          use_w4;           // 1 = noyau 4 pixels par work-item
	ocl_mem      d_qx[2], d_qy[2]; // indice de noyau par colonne / ligne de sortie
	int          nfx[2];
	ocl_mem      d_krad;
	const unsigned char *last_next;   // frame "next" déjà présente sur le GPU
	int          next_valid;
	int          ready;
} OclState;

static OclState g_ocl;

// ---------- noyaux OpenCL ----------
#define OCL_STR(...) #__VA_ARGS__

static const char ocl_src_pragma[] = "#pragma OPENCL FP_CONTRACT OFF\n";

static const char ocl_src_body[] = OCL_STR(

__constant int S5[25] = {
	-1,  -2, 0,  2,  1,
	-4,  -8, 0,  8,  4,
	-6, -12, 0, 12,  6,
	-4,  -8, 0,  8,  4,
	-1,  -2, 0,  2,  1 };

int rd(__global const uchar *b, uint i, int bps)
{
	return (bps == 1) ? ((int)b[i] << 8) : (int)(((__global const ushort *)b)[i]);
}

void wr(__global uchar *b, uint i, int v, int bps)
{
	v = v < 0 ? 0 : (v > 65535 ? 65535 : v);
	if(bps == 1) b[i] = (uchar)(v >> 8);
	else ((__global ushort *)b)[i] = (ushort)v;
}

// ---- fondu temporel : o = (a*(256-tw) + b*tw) >> 8 ----
__kernel void blend_k(__global const uchar *a, __global const uchar *b,
                      __global uchar *o, uint n, int bps, int tw)
{
	uint i = (uint)get_global_id(0);
	if(i >= n) return;
	int itw = WEIGHT_MAX - tw;
	if(bps == 1) {
		o[i] = (uchar)((a[i] * itw + b[i] * tw) >> WEIGHT_SHIFT);
	} else {
		int va = (int)(((__global const ushort *)a)[i]);
		int vb = (int)(((__global const ushort *)b)[i]);
		((__global ushort *)o)[i] = (ushort)((va * itw + vb * tw) >> WEIGHT_SHIFT);
	}
}

// ---- luma 8 bits d'un pixel (clampé) ----
int lum(__global const uchar *s, int bps, int chan, int x, int y, int w, int h)
{
	x = x < 0 ? 0 : (x >= w ? w - 1 : x);
	y = y < 0 ? 0 : (y >= h ? h - 1 : y);
	uint p = (uint)(y * w + x) * (uint)chan;
	if(chan == 1) {
		return (bps == 1) ? (int)s[p] : ((int)(((__global const ushort *)s)[p]) >> 8);
	}
	int r, g, b;
	if(bps == 1) { r = s[p]; g = s[p + 1]; b = s[p + 2]; }
	else {
		__global const ushort *q = (__global const ushort *)s;
		r = q[p] >> 8; g = q[p + 1] >> 8; b = q[p + 2] >> 8;
	}
	return (77 * r + 150 * g + 29 * b) >> 8;
}

// ---- gradient Sobel 5x5 -> plans G : gx, gy, mag, conf ----
__kernel void grad_k(__global const uchar *src, uint src_off, __global short *G,
                     int w, int h, int bps, int chan)
{
	int x = (int)get_global_id(0), y = (int)get_global_id(1);
	if(x >= w || y >= h) return;
	src += src_off;
	int gx = 0, gy = 0;
	for(int dy = -2; dy <= 2; dy++) {
		for(int dx = -2; dx <= 2; dx++) {
			int v = lum(src, bps, chan, x + dx, y + dy, w, h);
			gx += S5[(dy + 2) * 5 + (dx + 2)] * v;
			gy += S5[(dx + 2) * 5 + (dy + 2)] * v;
		}
	}
	gx >>= 6; gy >>= 6;
	int N = w * h, i = y * w + x;
	G[i]         = (short)gx;
	G[N + i]     = (short)gy;
	G[2 * N + i] = (short)(abs(gx) + abs(gy));
}

int mag_at_(__global const short *mg, int x, int y, int w, int h)
{
	x = x < 0 ? 0 : (x >= w ? w - 1 : x);
	y = y < 0 ? 0 : (y >= h ? h - 1 : y);
	return (int)mg[y * w + x];
}

int is_art(__global const short *mg, int sx, int sy, int mag_c, int w, int h)
{
	int mx_l = mag_at_(mg, sx - 1, sy, w, h);
	int mx_r = mag_at_(mg, sx + 1, sy, w, h);
	int my_u = mag_at_(mg, sx, sy - 1, w, h);
	int my_d = mag_at_(mg, sx, sy + 1, w, h);
	int narrow_x = (mx_l * NARROW_DROP < mag_c) && (mx_r * NARROW_DROP < mag_c);
	int narrow_y = (my_u * NARROW_DROP < mag_c) && (my_d * NARROW_DROP < mag_c);
	int cont_x = (mx_l * CONT_KEEP >= mag_c) && (mx_r * CONT_KEEP >= mag_c);
	int cont_y = (my_u * CONT_KEEP >= mag_c) && (my_d * CONT_KEEP >= mag_c);
	int is_cross = narrow_x && narrow_y;
	int is_axial = (narrow_x && cont_y) || (narrow_y && cont_x);
	return is_cross || is_axial;
}

// ---- carte de confiance (plan 3 de G) ----
__kernel void conf_k(__global short *G, int w, int h, int mode, int aniso_fp)
{
	int x = (int)get_global_id(0), y = (int)get_global_id(1);
	if(x >= w || y >= h) return;
	int N = w * h, i = y * w + x;
	int mag = (int)G[2 * N + i];
	int conf;
	if(mode == 7) {
		if(mag <= ANISO_LOW || is_art(G + 2 * N, x, y, mag, w, h)) conf = 0;
		else if(mag >= ANISO_HIGH) conf = WEIGHT_MAX;
		else conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
	} else {
		if(mag <= ANISO_LOW) conf = 0;
		else if(mag >= ANISO_HIGH) conf = WEIGHT_MAX;
		else conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
	}
	conf = (conf * aniso_fp) >> WEIGHT_SHIFT;
	G[3 * N + i] = (short)conf;
}

// ---- distance à la couture ----
__kernel void dist_k(__global const short *G, __global short *dist, int w, int h)
{
	int x = (int)get_global_id(0), y = (int)get_global_id(1);
	if(x >= w || y >= h) return;
	int N = w * h, i = y * w + x;
	__global const short *cf = G + 3 * N;
	if(cf[i] > 0) { dist[i] = 0; return; }
	const int R = COUTURE_N, R2 = COUTURE_N * COUTURE_N;
	int best2 = R2 + 1;
	for(int dj = -R; dj <= R; dj++) {
		int yy = y + dj;
		if(yy < 0 || yy >= h) continue;
		int dj2 = dj * dj;
		for(int di = -R; di <= R; di++) {
			int d2 = di * di + dj2;
			if(d2 == 0 || d2 > R2) continue;
			if(d2 >= best2) continue;
			int xx = x + di;
			if(xx < 0 || xx >= w) continue;
			if(cf[yy * w + xx] > 0) best2 = d2;
		}
	}
	if(best2 > R2) dist[i] = (short)((COUTURE_N + 1) * WEIGHT_MAX);
	else dist[i] = (short)((int)(sqrt((float)best2) * (float)WEIGHT_MAX + 0.5f));
}

// ---- convolution 6x6 sur la liste compacte des taps non nuls ----
// Noyau (KSTRIDE shorts) : [0]=n, octets 2..37 = indices de tap (dj+2)*6+(di+2),
// shorts 19.. = poids. Seuls les n taps utiles sont lus.
int4 conv(__global const uchar *src, int bps, int chan, int sx, int sy,
          int w, int h, __global const short *tab)
{
	int4 acc = (int4)(0, 0, 0, 0);
	const int n = (int)tab[0];
	__global const uchar *ib = (__global const uchar *)(tab + 1);
	__global const short *wp = tab + 19;
	if(sx >= 2 && sx + 3 < w && sy >= 2 && sy + 3 < h)
	{
		// intérieur : plus aucun test de bord, simple décalage linéaire
		const uint base = (uint)((sy - 2) * w + (sx - 2)) * (uint)chan;
		for(int t = 0; t < n; t++) {
			int id = (int)ib[t];
			int dj = id / 6, di = id - dj * 6;
			int wt = (int)wp[t];
			uint p = base + (uint)(dj * w + di) * (uint)chan;
			acc.x += rd(src, p, bps) * wt;
			if(chan == 3) {
				acc.y += rd(src, p + 1, bps) * wt;
				acc.z += rd(src, p + 2, bps) * wt;
			}
		}
		return acc;
	}
	for(int t = 0; t < n; t++) {
		int id = (int)ib[t];
		int dj = id / 6, di = id - dj * 6;
		int wt = (int)wp[t];
		int ly = sy + dj - 2; ly = ly < 0 ? 0 : (ly >= h ? h - 1 : ly);
		int lx = sx + di - 2; lx = lx < 0 ? 0 : (lx >= w ? w - 1 : lx);
		uint p = (uint)(ly * w + lx) * (uint)chan;
		acc.x += rd(src, p, bps) * wt;
		if(chan == 3) {
			acc.y += rd(src, p + 1, bps) * wt;
			acc.z += rd(src, p + 2, bps) * wt;
		}
	}
	return acc;
}

// ---- noyau anisotrope orienté ; .w = 1 si valide, 0 sinon ----
int4 aniso(__global const uchar *src, int bps, int chan, int sx, int sy,
           int fx, int fy, int gx, int gy, int conf, int w, int h,
           __global const int *krad)
{
	float glen = sqrt((float)gx * (float)gx + (float)gy * (float)gy);
	if(glen < 0.5f) glen = 0.5f;
	float nx = (float)gx / glen, ny = (float)gy / glen;
	float ratio = 1.0f + (ANISO_RATIO - 1.0f) * ((float)conf / (float)WEIGHT_MAX);
	float ffx = (float)fx / (float)WEIGHT_MAX;
	float ffy = (float)fy / (float)WEIGHT_MAX;
	int4 acc = (int4)(0, 0, 0, 0);
	int wsum = 0;
	for(int dj = -2; dj <= 3; dj++) {
		for(int di = -2; di <= 3; di++) {
			float ddx = (float)di - ffx;
			float ddy = (float)dj - ffy;
			float d_across = ddx * nx + ddy * ny;
			float d_along  = -ddx * ny + ddy * nx;
			float da = d_along / ratio;
			float r = sqrt(d_across * d_across + da * da);
			int ri = (int)(r * (float)WEIGHT_MAX + 0.5f);
			int ww = (ri < KAISER_RAD_SIZE) ? krad[ri] : 0;
			if(ww == 0) continue;
			int lx = sx + di, ly = sy + dj;
			lx = lx < 0 ? 0 : (lx >= w ? w - 1 : lx);
			ly = ly < 0 ? 0 : (ly >= h ? h - 1 : ly);
			uint p = (uint)(ly * w + lx) * (uint)chan;
			acc.x += rd(src, p, bps) * ww;
			if(chan == 3) {
				acc.y += rd(src, p + 1, bps) * ww;
				acc.z += rd(src, p + 2, bps) * ww;
			}
			wsum += ww;
		}
	}
	if(wsum > 0) return (int4)(acc.x / wsum, acc.y / wsum, acc.z / wsum, 1);
	return (int4)(0, 0, 0, 0);
}

// interpolation bilinéaire en virgule fixe d'un plan de short
int bil(__global const short *pl, int i00, int i10, int i01, int i11, int fx, int fy)
{
	int t = ((int)pl[i00] * (WEIGHT_MAX - fx) + (int)pl[i10] * fx) >> WEIGHT_SHIFT;
	int b = ((int)pl[i01] * (WEIGHT_MAX - fx) + (int)pl[i11] * fx) >> WEIGHT_SHIFT;
	return (t * (WEIGHT_MAX - fy) + b * fy) >> WEIGHT_SHIFT;
}

// ---- conversion YUV 4:4:4 planaire -> RGB24 (une frame de sortie r) ----
// Matrice BT.601/709 et plage limitée/pleine fixées à la compilation.
__kernel void yuv2rgb_k(__global const uchar *src, __global uchar *dst,
                        uint y_off, uint u_off, uint v_off, uint step,
                        int w, int h, int r)
{
	int x = (int)get_global_id(0), y = (int)get_global_id(1);
	if(x >= w || y >= h) return;
	const uint i = (uint)(y * w + x);
	const uint base = (uint)r * step + i;
	int Y = (int)src[y_off + base];
	int U = (int)src[u_off + base] - 128;
	int V = (int)src[v_off + base] - 128;
	// coefficients 8.8 : [cy, v->r, u->g, v->g, u->b]
	const int full = YUV_FULL, hd = YUV_709;
	const int cy  = full ? 256 : 298;
	const int crv = hd ? (full ? 403 : 459) : (full ? 359 : 409);
	const int cgu = hd ? (full ?  48 :  55) : (full ?  88 : 100);
	const int cgv = hd ? (full ? 120 : 136) : (full ? 183 : 208);
	const int cbu = hd ? (full ? 475 : 541) : (full ? 454 : 516);
	if(!full) Y -= 16;
	int R = (cy * Y + crv * V + 128) >> 8;
	int G = (cy * Y - cgu * U - cgv * V + 128) >> 8;
	int B = (cy * Y + cbu * U + 128) >> 8;
	R = R < 0 ? 0 : (R > 255 ? 255 : R);
	G = G < 0 ? 0 : (G > 255 ? 255 : G);
	B = B < 0 ? 0 : (B > 255 ? 255 : B);
	const uint o = ((uint)r * (uint)(w * h) + i) * 3u;
	dst[o] = (uchar)R; dst[o + 1] = (uchar)G; dst[o + 2] = (uchar)B;
}

// ---- modes 3/4/5, plans seuls : 4 pixels de sortie par work-item ----
// Les pixels voisins partagent presque toujours la même source (sx, sy) : on
// charge la fenêtre 6x6 une seule fois en registres puis on applique les 4
// noyaux denses (36 poids). Mêmes sommes entières que scale_k, donc mêmes
// résultats. Hors cas uniforme / hors intérieur : repli pixel par pixel.
// (même signature que scale_k ; tabA désigne ici la table dense)
__kernel void scale4_k(
	__global const uchar *src, uint src_off,
	__global uchar *dst, uint dst_off,
	__global const int *w_pos, __global const int *w_src, __global const int *w_flag,
	__global const int *h_pos, __global const int *h_src, __global const int *h_flag,
	__global const int *w_qi, __global const int *h_qi, int nfx,
	int wbase, int hbase,
	int in_w, int in_h, int out_w, int out_h,
	__global const short *tabA, __global const short *tabB, __global const short *tabC,
	__global const int *krad, __global const short *G, __global const short *dist,
	int aniso_fp)
{
	const int bps_in = BPS_IN, bps_out = BPS_OUT;
	const int x0 = (int)get_global_id(0) * 4, y = (int)get_global_id(1);
	if(x0 >= out_w || y >= out_h) return;
	src += src_off;
	dst += dst_off;
	const int nx = (out_w - x0) < 4 ? (out_w - x0) : 4;
	const int sy = h_pos[hbase + y] >> WEIGHT_SHIFT;
	const int sx0 = w_pos[wbase + x0] >> WEIGHT_SHIFT;
	const uint hq = (uint)h_qi[hbase + y];

	int uni = (sy >= 2 && sy + 3 < in_h && sx0 >= 2 && sx0 + 3 < in_w);
	for(int i = 1; i < nx; i++)
		uni = uni && ((w_pos[wbase + x0 + i] >> WEIGHT_SHIFT) == sx0);

	if(uni)
	{
		int win[36];
		const uint base = (uint)((sy - 2) * in_w + (sx0 - 2));
		#pragma unroll
		for(int j = 0; j < 6; j++) {
			#pragma unroll
			for(int i = 0; i < 6; i++)
				win[j * 6 + i] = rd(src, base + (uint)(j * in_w + i), bps_in);
		}
		for(int i = 0; i < nx; i++)
		{
			__global const short *tb = tabA + ((hq * (uint)nfx + (uint)w_qi[wbase + x0 + i]) * 36u);
			int acc = 0;
			#pragma unroll
			for(int k = 0; k < 36; k += 4) {
				int4 wv = convert_int4(vload4(0, tb + k));
				acc += win[k] * wv.x + win[k + 1] * wv.y + win[k + 2] * wv.z + win[k + 3] * wv.w;
			}
			wr(dst, (uint)(y * out_w + x0 + i), acc >> WEIGHT_SHIFT, bps_out);
		}
	}
	else
	{
		// repli : table compacte (tabB), noyau par pixel
		for(int i = 0; i < nx; i++)
		{
			int wp = w_pos[wbase + x0 + i];
			int sx = wp >> WEIGHT_SHIFT;
			uint tidx = (hq * (uint)nfx + (uint)w_qi[wbase + x0 + i]) * 56u;
			int4 v = conv(src, bps_in, 1, sx, sy, in_w, in_h, tabB + tidx) >> (int4)(WEIGHT_SHIFT);
			wr(dst, (uint)(y * out_w + x0 + i), v.x, bps_out);
		}
	}
}

// ---- mise à l'échelle : un work-item par pixel de sortie ----
__kernel void scale_k(
	__global const uchar *src, uint src_off,
	__global uchar *dst, uint dst_off,
	__global const int *w_pos, __global const int *w_src, __global const int *w_flag,
	__global const int *h_pos, __global const int *h_src, __global const int *h_flag,
	__global const int *w_qi, __global const int *h_qi, int nfx,
	int wbase, int hbase,
	int in_w, int in_h, int out_w, int out_h,
	__global const short *tabA, __global const short *tabB, __global const short *tabC,
	__global const int *krad, __global const short *G, __global const short *dist,
	int aniso_fp)
{
	const int chan = CHAN, bps_in = BPS_IN, bps_out = BPS_OUT, mode = MODE;
	int x = (int)get_global_id(0), y = (int)get_global_id(1);
	if(x >= out_w || y >= out_h) return;
	src += src_off;
	dst += dst_off;
	uint opix = (uint)(y * out_w + x) * (uint)chan;

	if(mode >= 3)
	{
		int hp = h_pos[hbase + y];
		int sy = hp >> WEIGHT_SHIFT, fy = hp & (WEIGHT_MAX - 1);
		int wp = w_pos[wbase + x];
		int sx = wp >> WEIGHT_SHIFT, fx = wp & (WEIGHT_MAX - 1);
		uint tidx = ((uint)h_qi[hbase + y] * (uint)nfx + (uint)w_qi[wbase + x]) * 56u;
		int4 v;

		if(mode == 6 || mode == 7)
		{
			int N = in_w * in_h;
			int sx1 = sx + 1 < in_w ? sx + 1 : sx;
			int sy1 = sy + 1 < in_h ? sy + 1 : sy;
			int i00 = sy * in_w + sx,  i10 = sy * in_w + sx1;
			int i01 = sy1 * in_w + sx, i11 = sy1 * in_w + sx1;
			int igx = bil(G,         i00, i10, i01, i11, fx, fy);
			int igy = bil(G + N,     i00, i10, i01, i11, fx, fy);
			int mag = bil(G + 2 * N, i00, i10, i01, i11, fx, fy);

			int conf;
			if(mode == 7) {
				if(mag <= ANISO_LOW || is_art(G + 2 * N, sx, sy, mag, in_w, in_h)) conf = 0;
				else if(mag >= ANISO_HIGH) conf = WEIGHT_MAX;
				else conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
			} else {
				if(mag <= ANISO_LOW) conf = 0;
				else if(mag >= ANISO_HIGH) conf = WEIGHT_MAX;
				else conf = (mag - ANISO_LOW) * WEIGHT_MAX / (ANISO_HIGH - ANISO_LOW);
			}
			conf = (conf * aniso_fp) >> WEIGHT_SHIFT;

			if(conf == 0)
			{
				int d_fp = bil(dist, i00, i10, i01, i11, fx, fy);
				int dmax = (COUTURE_N + 1) * WEIGHT_MAX;
				int alpha = (d_fp >= dmax) ? 0 : (WEIGHT_MAX - (d_fp * WEIGHT_MAX / dmax));
				int4 s = conv(src, bps_in, chan, sx, sy, in_w, in_h, tabB + tidx) >> (int4)(WEIGHT_SHIFT);
				if(alpha == 0) {
					v = s;
				} else {
					int4 m = conv(src, bps_in, chan, sx, sy, in_w, in_h, tabC + tidx) >> (int4)(WEIGHT_SHIFT);
					int ia = WEIGHT_MAX - alpha;
					v = (m * alpha + s * ia) >> (int4)(WEIGHT_SHIFT);
				}
			}
			else
			{
				int4 a = aniso(src, bps_in, chan, sx, sy, fx, fy, igx, igy, conf, in_w, in_h, krad);
				if(a.w != 0) v = a;
				else v = conv(src, bps_in, chan, sx, sy, in_w, in_h, tabA + tidx) >> (int4)(WEIGHT_SHIFT);
			}
		}
		else
		{
			v = conv(src, bps_in, chan, sx, sy, in_w, in_h, tabA + tidx) >> (int4)(WEIGHT_SHIFT);
		}

		wr(dst, opix, v.x, bps_out);
		if(chan == 3) {
			wr(dst, opix + 1, v.y, bps_out);
			wr(dst, opix + 2, v.z, bps_out);
		}
	}
	else
	{
		// modes 0/1/2 : plus proche voisin / bilinéaire sur pixels dupliqués
		int sy = h_src[hbase + y];
		int hl = (mode >= 1) ? h_flag[hbase + y] : 0;
		int sy2 = sy;
		if(hl ==  1) sy2 = (sy + 1 < in_h) ? sy + 1 : sy;
		if(hl == -1) sy2 = (sy - 1 >= 0)   ? sy - 1 : sy;
		int sx = w_src[wbase + x];
		int wl = (mode >= 1) ? w_flag[wbase + x] : 0;
		int sxb = sx;
		if(wl ==  1) sxb = (sx + 1 < in_w) ? sx + 1 : sx;
		if(wl == -1) sxb = (sx - 1 >= 0)   ? sx - 1 : sx;

		for(int c = 0; c < chan; c++)
		{
			int a = rd(src, (uint)(sy * in_w + sx) * (uint)chan + (uint)c, bps_in);
			int v;
			if(hl == 0) {
				if(wl == 0) v = a;
				else v = (a + rd(src, (uint)(sy * in_w + sxb) * (uint)chan + (uint)c, bps_in)) >> 1;
			} else {
				int a2 = rd(src, (uint)(sy2 * in_w + sx) * (uint)chan + (uint)c, bps_in);
				if(wl == 0) v = (a + a2) >> 1;
				else {
					int b1 = rd(src, (uint)(sy * in_w + sxb) * (uint)chan + (uint)c, bps_in);
					int b2 = rd(src, (uint)(sy2 * in_w + sxb) * (uint)chan + (uint)c, bps_in);
					v = (a + b1 + a2 + b2) >> 2;
				}
			}
			wr(dst, opix + (uint)c, v, bps_out);
		}
	}
}

);

// ---------- chargement de la bibliothèque ----------
static void *ocl_sym(void *lib, const char *name)
{
#ifdef _WIN32
	return (void *)GetProcAddress((HMODULE)lib, name);
#else
	return dlsym(lib, name);
#endif
}

static int ocl_load_library(void)
{
	if(g_ocl.loaded) return 0;
#ifdef _WIN32
	g_ocl.lib = (void *)LoadLibraryA("OpenCL.dll");
#elif defined(__APPLE__)
	g_ocl.lib = dlopen("/System/Library/Frameworks/OpenCL.framework/OpenCL", RTLD_LAZY);
#else
	g_ocl.lib = dlopen("libOpenCL.so.1", RTLD_LAZY);
	if(!g_ocl.lib) g_ocl.lib = dlopen("libOpenCL.so", RTLD_LAZY);
#endif
	if(!g_ocl.lib) return -1;

	#define OCL_LOAD(field, sym) \
		do { *(void **)&g_ocl.api.field = ocl_sym(g_ocl.lib, sym); \
		     if(!g_ocl.api.field) return -1; } while(0)
	OCL_LOAD(GetPlatformIDs,          "clGetPlatformIDs");
	OCL_LOAD(GetPlatformInfo,         "clGetPlatformInfo");
	OCL_LOAD(GetDeviceIDs,            "clGetDeviceIDs");
	OCL_LOAD(GetDeviceInfo,           "clGetDeviceInfo");
	OCL_LOAD(CreateContext,           "clCreateContext");
	OCL_LOAD(CreateCommandQueue,      "clCreateCommandQueue");
	OCL_LOAD(CreateProgramWithSource, "clCreateProgramWithSource");
	OCL_LOAD(BuildProgram,            "clBuildProgram");
	OCL_LOAD(GetProgramBuildInfo,     "clGetProgramBuildInfo");
	OCL_LOAD(CreateKernel,            "clCreateKernel");
	OCL_LOAD(SetKernelArg,            "clSetKernelArg");
	OCL_LOAD(EnqueueNDRangeKernel,    "clEnqueueNDRangeKernel");
	OCL_LOAD(CreateBuffer,            "clCreateBuffer");
	OCL_LOAD(EnqueueReadBuffer,       "clEnqueueReadBuffer");
	OCL_LOAD(EnqueueWriteBuffer,      "clEnqueueWriteBuffer");
	OCL_LOAD(Finish,                  "clFinish");
	OCL_LOAD(Flush,                   "clFlush");
	OCL_LOAD(WaitForEvents,           "clWaitForEvents");
	OCL_LOAD(ReleaseEvent,            "clReleaseEvent");
	OCL_LOAD(ReleaseMemObject,        "clReleaseMemObject");
	OCL_LOAD(ReleaseKernel,           "clReleaseKernel");
	OCL_LOAD(ReleaseProgram,          "clReleaseProgram");
	OCL_LOAD(ReleaseCommandQueue,     "clReleaseCommandQueue");
	OCL_LOAD(ReleaseContext,          "clReleaseContext");
	#undef OCL_LOAD
	g_ocl.loaded = 1;
	return 0;
}

// ---------- énumération des périphériques ----------
typedef struct {
	char       plat_name[128];
	char       name[256];
	ocl_ulong  type;
	ocl_uint   cu, mhz, avail, compiler;
	ocl_ulong  mem;
} OclDevInfo;

static int ocl_enumerate(void)
{
	OclApi *a = &g_ocl.api;
	ocl_uint np = 0;
	g_ocl.ndev = 0;
	if(a->GetPlatformIDs(16, g_ocl.plats, &np) != 0 || np == 0) return -1;
	for(ocl_uint p = 0; p < np; p++) {
		ocl_uint nd = 0;
		ocl_device tmp[32];
		if(a->GetDeviceIDs(g_ocl.plats[p], OCL_DEVICE_TYPE_ALL, 32, tmp, &nd) != 0) continue;
		for(ocl_uint d = 0; d < nd && g_ocl.ndev < 64; d++) {
			g_ocl.devs[g_ocl.ndev]     = tmp[d];
			g_ocl.dev_plat[g_ocl.ndev] = (int)p;
			g_ocl.ndev++;
		}
	}
	return g_ocl.ndev > 0 ? 0 : -1;
}

static void ocl_dev_info(int idx, OclDevInfo *o)
{
	OclApi *a = &g_ocl.api;
	ocl_device d = g_ocl.devs[idx];
	memset(o, 0, sizeof(*o));
	a->GetPlatformInfo(g_ocl.plats[g_ocl.dev_plat[idx]], OCL_PLATFORM_NAME, sizeof(o->plat_name) - 1, o->plat_name, NULL);
	a->GetDeviceInfo(d, OCL_DEVICE_NAME,              sizeof(o->name) - 1, o->name, NULL);
	a->GetDeviceInfo(d, OCL_DEVICE_TYPE_INFO,         sizeof(o->type),     &o->type, NULL);
	a->GetDeviceInfo(d, OCL_DEVICE_MAX_COMPUTE_UNITS, sizeof(o->cu),       &o->cu, NULL);
	a->GetDeviceInfo(d, OCL_DEVICE_MAX_CLOCK_FREQ,    sizeof(o->mhz),      &o->mhz, NULL);
	a->GetDeviceInfo(d, OCL_DEVICE_AVAILABLE,         sizeof(o->avail),    &o->avail, NULL);
	a->GetDeviceInfo(d, OCL_DEVICE_COMPILER_AVAILABLE,sizeof(o->compiler), &o->compiler, NULL);
	a->GetDeviceInfo(d, OCL_DEVICE_GLOBAL_MEM_SIZE,   sizeof(o->mem),      &o->mem, NULL);
}

// Affiche les périphériques OpenCL (option --ocl-list)
static int ocl_list_devices(void)
{
	if(ocl_load_library() != 0) { fprintf(stderr, "OpenCL: library not found\n"); return -1; }
	if(ocl_enumerate() != 0)    { fprintf(stderr, "OpenCL: no device found\n");   return -1; }
	for(int i = 0; i < g_ocl.ndev; i++) {
		OclDevInfo o; ocl_dev_info(i, &o);
		fprintf(stderr, "  [%d] %s | %s | %s | %u CU @ %u MHz | %u MB%s\n",
			i, o.plat_name, o.name,
			(o.type & OCL_DEVICE_TYPE_GPU) ? "GPU" : (o.type & OCL_DEVICE_TYPE_CPU) ? "CPU" : "other",
			o.cu, o.mhz, (unsigned)(o.mem >> 20),
			(o.avail && o.compiler) ? "" : " (unusable)");
	}
	return 0;
}

// ---------- initialisation : contexte, file, programme ----------
#define OCL_FAIL(msg, e) do { fprintf(stderr, "OpenCL: %s (error %d)\n", msg, (int)(e)); return -1; } while(0)

// want < 0 : choix automatique du meilleur GPU ; sinon index de --ocl-list
static int ocl_init(int want)
{
	OclApi *a = &g_ocl.api;
	ocl_int e;

	if(ocl_load_library() != 0) { fprintf(stderr, "OpenCL: library not found\n"); return -1; }
	if(ocl_enumerate() != 0)    { fprintf(stderr, "OpenCL: no device found\n");   return -1; }

	int pick = -1;
	if(want >= 0) {
		if(want >= g_ocl.ndev) { fprintf(stderr, "OpenCL: device index %d out of range (see --ocl-list)\n", want); return -1; }
		pick = want;
	} else {
		double best = -1.0;
		for(int i = 0; i < g_ocl.ndev; i++) {
			OclDevInfo o; ocl_dev_info(i, &o);
			if(!(o.type & OCL_DEVICE_TYPE_GPU) || !o.avail || !o.compiler) continue;
			double score = (double)(o.cu ? o.cu : 1) * (double)(o.mhz ? o.mhz : 1);
			if(score > best) { best = score; pick = i; }
		}
		if(pick < 0) { fprintf(stderr, "OpenCL: no usable GPU found (see --ocl-list)\n"); return -1; }
	}

	OclDevInfo oi; ocl_dev_info(pick, &oi);
	if(!oi.avail || !oi.compiler) { fprintf(stderr, "OpenCL: device [%d] is unusable\n", pick); return -1; }
	snprintf(g_ocl.dev_name, sizeof(g_ocl.dev_name), "%s (%s)", oi.name, oi.plat_name);

	g_ocl.dev  = g_ocl.devs[pick];
	g_ocl.plat = g_ocl.plats[g_ocl.dev_plat[pick]];

	intptr_t props[3] = { OCL_CONTEXT_PLATFORM, (intptr_t)g_ocl.plat, 0 };
	g_ocl.ctx = a->CreateContext(props, 1, &g_ocl.dev, NULL, NULL, &e);
	if(!g_ocl.ctx) OCL_FAIL("clCreateContext failed", e);
	g_ocl.q = a->CreateCommandQueue(g_ocl.ctx, g_ocl.dev, 0, &e);
	if(!g_ocl.q) OCL_FAIL("clCreateCommandQueue failed", e);

	fprintf(stderr, "OpenCL: using %s, %u CU @ %u MHz, %u MB\n",
		g_ocl.dev_name, oi.cu, oi.mhz, (unsigned)(oi.mem >> 20));
	return 0;
}

// Compile le programme spécialisé pour ce run : le mode, le nombre de canaux
// et les profondeurs de bits sont des constantes de compilation, ce qui
// supprime les branchements et réduit la pression sur les registres (gain
// net sur les petits GPU intégrés).
static int ocl_build_program(const OclSetup *s)
{
	OclApi *a = &g_ocl.api;
	ocl_int e;

	const char *srcs[2] = { ocl_src_pragma, ocl_src_body };
	g_ocl.prog = a->CreateProgramWithSource(g_ocl.ctx, 2, srcs, NULL, &e);
	if(!g_ocl.prog) OCL_FAIL("clCreateProgramWithSource failed", e);

	char opts[640];
	snprintf(opts, sizeof(opts),
		"-DWEIGHT_SHIFT=%d -DWEIGHT_MAX=%d -DNQUANT=%d -DNQUANT_SHIFT=%d "
		"-DKAISER_RAD_SIZE=%d -DANISO_LOW=%d -DANISO_HIGH=%d -DANISO_RATIO=%.6ff "
		"-DNARROW_DROP=%d -DCONT_KEEP=%d -DCOUTURE_N=%d "
		"-DMODE=%d -DCHAN=%d -DBPS_IN=%d -DBPS_OUT=%d -DYUV_709=%d -DYUV_FULL=%d",
		WEIGHT_SHIFT, WEIGHT_MAX, NQUANT, NQUANT_SHIFT,
		KAISER_RAD_SIZE, ANISO_LOW, ANISO_HIGH, (double)ANISO_RATIO,
		NARROW_DROP, CONT_KEEP, COUTURE_N,
		s->mode, s->planes[0].chan, s->bps_in, s->bps_out, s->yuv_709, s->yuv_full);

	e = a->BuildProgram(g_ocl.prog, 1, &g_ocl.dev, opts, NULL, NULL);
	if(e != 0) {
		size_t len = 0;
		a->GetProgramBuildInfo(g_ocl.prog, g_ocl.dev, OCL_PROGRAM_BUILD_LOG, 0, NULL, &len);
		char *log = (char *)malloc(len + 1);
		if(log) {
			log[0] = 0;
			a->GetProgramBuildInfo(g_ocl.prog, g_ocl.dev, OCL_PROGRAM_BUILD_LOG, len, log, NULL);
			log[len] = 0;
			fprintf(stderr, "OpenCL build log:\n%s\n", log);
			free(log);
		}
		OCL_FAIL("kernel compilation failed", e);
	}

	#define OCL_KERNEL(var, name) \
		do { g_ocl.var = a->CreateKernel(g_ocl.prog, name, &e); \
		     if(!g_ocl.var) OCL_FAIL("clCreateKernel failed: " name, e); } while(0)
	OCL_KERNEL(k_blend, "blend_k");
	OCL_KERNEL(k_grad,  "grad_k");
	OCL_KERNEL(k_conf,  "conf_k");
	OCL_KERNEL(k_dist,  "dist_k");
	OCL_KERNEL(k_scale, "scale_k");
	OCL_KERNEL(k_scale4, "scale4_k");
	OCL_KERNEL(k_yuv, "yuv2rgb_k");
	#undef OCL_KERNEL
	return 0;
}

// ---------- buffers & tables ----------
static ocl_mem ocl_buf(size_t bytes, ocl_ulong flags, const void *host)
{
	ocl_int e;
	if(bytes == 0) bytes = 4;
	ocl_mem m = g_ocl.api.CreateBuffer(g_ocl.ctx, flags, bytes, (void *)host, &e);
	if(!m) fprintf(stderr, "OpenCL: cannot allocate %.1f MB (error %d)\n", bytes / 1048576.0, (int)e);
	return m;
}

// Une carte n'utilise que quelques positions fractionnaires distinctes (le
// motif se répète avec la période du rapport d'échelle). Plutôt qu'une table
// de NQUANT*NQUANT noyaux (19 Mo, lue au hasard, qui sature le cache des
// petits GPU), on ne garde que les noyaux réellement utilisés : quelques Ko
// qui tiennent en cache. Les poids sont strictement les mêmes.
//
// Pour chaque colonne/ligne de sortie, calcule l'indice du noyau quantifié
// parmi les positions distinctes ; renvoie leur nombre (liste dans `list`).
static int ocl_quant_axis(const int *pos, int len, int *idx_out, int *list)
{
	int map[NQUANT], seen[NQUANT];
	memset(seen, 0, sizeof(seen));
	for(int i = 0; i < len; i++)
		seen[(pos[i] & (WEIGHT_MAX - 1)) << NQUANT_SHIFT] = 1;
	int n = 0;
	for(int q = 0; q < NQUANT; q++) {
		map[q] = -1;
		if(seen[q]) { map[q] = n; list[n++] = q; }
	}
	for(int i = 0; i < len; i++)
		idx_out[i] = map[(pos[i] & (WEIGHT_MAX - 1)) << NQUANT_SHIFT];
	return n;
}

#define OCL_KSTRIDE 56   // shorts par noyau (voir conv() côté OpenCL)
// Table compacte nfy*nfx noyaux (n, indices des taps, poids), taps non nuls seulement.
// Les poids ne dépendent pas du stride : la même table sert à Y et UV.
static ocl_mem ocl_make_table(const int *lut, int lut_size,
                              const int *ylist, int nfy, const int *xlist, int nfx)
{
	const size_t nk = (size_t)nfy * nfx;
	int16_t *t = (int16_t *)calloc(nk * OCL_KSTRIDE, sizeof(int16_t));
	if(!t) { fprintf(stderr, "OpenCL: out of memory (kernel table)\n"); return NULL; }
	for(int a = 0; a < nfy; a++)
	for(int b = 0; b < nfx; b++)
	{
		int16_t w[KERNEL_MAX_TAPS];
		int8_t  dx[KERNEL_MAX_TAPS], dy[KERNEL_MAX_TAPS];
		const int n = kernel_taps(lut, lut_size, ylist[a], xlist[b], w, dx, dy);
		int16_t *d = t + ((size_t)a * nfx + b) * OCL_KSTRIDE;
		d[0] = (int16_t)n;
		for(int i = 0; i < n; i++) {
			((uint8_t *)(d + 1))[i] = (uint8_t)((dy[i] + 2) * 6 + (dx[i] + 2));
			d[19 + i] = w[i];
		}
	}
	ocl_mem m = ocl_buf(nk * OCL_KSTRIDE * sizeof(int16_t), OCL_MEM_READ_ONLY | OCL_MEM_COPY_HOST_PTR, t);
	free(t);
	return m;
}

// Table dense (36 poids par noyau, ordre (dj+2)*6+(di+2)) pour scale4_k.
static ocl_mem ocl_make_dense(const int *lut, int lut_size,
                              const int *ylist, int nfy, const int *xlist, int nfx)
{
	const size_t nk = (size_t)nfy * nfx;
	int16_t *t = (int16_t *)calloc(nk * 36, sizeof(int16_t));
	if(!t) { fprintf(stderr, "OpenCL: out of memory (kernel table)\n"); return NULL; }
	for(int a = 0; a < nfy; a++)
	for(int b = 0; b < nfx; b++)
	{
		int16_t w[KERNEL_MAX_TAPS];
		int8_t  dx[KERNEL_MAX_TAPS], dy[KERNEL_MAX_TAPS];
		const int n = kernel_taps(lut, lut_size, ylist[a], xlist[b], w, dx, dy);
		int16_t *d = t + ((size_t)a * nfx + b) * 36;
		for(int i = 0; i < n; i++)
			d[(dy[i] + 2) * 6 + (dx[i] + 2)] = w[i];
	}
	ocl_mem m = ocl_buf(nk * 36 * sizeof(int16_t), OCL_MEM_READ_ONLY | OCL_MEM_COPY_HOST_PTR, t);
	free(t);
	return m;
}

// Indices + tables pour une carte (idx 0 = luma/RGB, 1 = chroma).
static int ocl_setup_tables(int idx, const OclSetup *s)
{
	const BresenhamMap *m = s->maps[idx];
	const ocl_ulong f = OCL_MEM_READ_ONLY | OCL_MEM_COPY_HOST_PTR;
	int xlist[NQUANT], ylist[NQUANT];
	int *qx = (int *)malloc((size_t)m->w_length * sizeof(int));
	int *qy = (int *)malloc((size_t)m->h_length * sizeof(int));
	if(!qx || !qy) { free(qx); free(qy); return -1; }
	const int nfx = ocl_quant_axis(m->w_pos, m->w_length, qx, xlist);
	const int nfy = ocl_quant_axis(m->h_pos, m->h_length, qy, ylist);
	g_ocl.nfx[idx] = nfx;
	g_ocl.d_qx[idx] = ocl_buf((size_t)m->w_length * sizeof(int), f, qx);
	g_ocl.d_qy[idx] = ocl_buf((size_t)m->h_length * sizeof(int), f, qy);
	free(qx); free(qy);
	if(!g_ocl.d_qx[idx] || !g_ocl.d_qy[idx]) return -1;

	g_ocl.d_tab[idx][0] = ocl_make_table(s->lutA, s->lut_size, ylist, nfy, xlist, nfx);
	if(!g_ocl.d_tab[idx][0]) return -1;
	g_ocl.d_tab[idx][1] = s->lutB ? ocl_make_table(s->lutB, s->lut_size, ylist, nfy, xlist, nfx) : g_ocl.d_tab[idx][0];
	g_ocl.d_tab[idx][2] = s->lutC ? ocl_make_table(s->lutC, s->lut_size, ylist, nfy, xlist, nfx) : g_ocl.d_tab[idx][0];
	if(!g_ocl.d_tab[idx][1] || !g_ocl.d_tab[idx][2]) return -1;
	if(g_ocl.use_w4) {
		g_ocl.d_dense[idx] = ocl_make_dense(s->lutA, s->lut_size, ylist, nfy, xlist, nfx);
		if(!g_ocl.d_dense[idx]) return -1;
	}
	return 0;
}

static int ocl_upload_map(int idx, const BresenhamMap *m)
{
	const ocl_ulong f = OCL_MEM_READ_ONLY | OCL_MEM_COPY_HOST_PTR;
	const size_t wb = (size_t)m->w_length * sizeof(int);
	const size_t hb = (size_t)m->h_length * sizeof(int);
	g_ocl.d_map[idx][0] = ocl_buf(wb, f, m->w_pos);
	g_ocl.d_map[idx][1] = ocl_buf(wb, f, m->w_src);
	g_ocl.d_map[idx][2] = ocl_buf(wb, f, m->w_flag);
	g_ocl.d_map[idx][3] = ocl_buf(hb, f, m->h_pos);
	g_ocl.d_map[idx][4] = ocl_buf(hb, f, m->h_src);
	g_ocl.d_map[idx][5] = ocl_buf(hb, f, m->h_flag);
	for(int i = 0; i < 6; i++) if(!g_ocl.d_map[idx][i]) return -1;
	return 0;
}

// À appeler une fois les cartes et LUT définitives (après tous les ajustements).
static int ocl_setup(const OclSetup *s)
{
	g_ocl.s = *s;
	if(ocl_build_program(s) != 0) return -1;
	const ocl_ulong rw = OCL_MEM_READ_WRITE;

	g_ocl.d_in    = ocl_buf(s->in_len,  rw, NULL);
	g_ocl.d_next  = ocl_buf(s->in_len,  rw, NULL);
	g_ocl.d_blend = ocl_buf(s->in_len,  rw, NULL);
	g_ocl.d_out   = ocl_buf(s->out_len, rw, NULL);
	if(!g_ocl.d_in || !g_ocl.d_next || !g_ocl.d_blend || !g_ocl.d_out) return -1;
	if(s->rgb_out) {
		g_ocl.d_rgb = ocl_buf(s->out_len, rw, NULL);
		if(!g_ocl.d_rgb) return -1;
	}

	if(ocl_upload_map(0, s->maps[0]) != 0) return -1;
	if(s->maps[1] == s->maps[0]) {
		// même carte pour la chroma : on réutilise les mêmes buffers
		for(int i = 0; i < 6; i++) g_ocl.d_map[1][i] = g_ocl.d_map[0][i];
	} else if(ocl_upload_map(1, s->maps[1]) != 0) return -1;

	g_ocl.use_w4 = (s->mode == 5 && s->planes[0].chan == 1 && !getenv("TS_OCL_NOW4"));
	if(s->mode >= 3) {
		if(ocl_setup_tables(0, s) != 0) return -1;
		if(s->maps[1] == s->maps[0]) {
			g_ocl.nfx[1] = g_ocl.nfx[0];
			g_ocl.d_qx[1] = g_ocl.d_qx[0];
			g_ocl.d_qy[1] = g_ocl.d_qy[0];
			for(int i = 0; i < 3; i++) g_ocl.d_tab[1][i] = g_ocl.d_tab[0][i];
			g_ocl.d_dense[1] = g_ocl.d_dense[0];
		} else if(ocl_setup_tables(1, s) != 0) return -1;
	} else {
		// modes 0..2 : tampons factices (jamais lus)
		for(int m = 0; m < 2; m++) {
			g_ocl.d_qx[m] = ocl_buf(4, rw, NULL);
			g_ocl.d_qy[m] = ocl_buf(4, rw, NULL);
			for(int i = 0; i < 3; i++) g_ocl.d_tab[m][i] = ocl_buf(4, rw, NULL);
			g_ocl.nfx[m] = 1;
		}
	}

	if(s->krad)
		g_ocl.d_krad = ocl_buf((size_t)s->krad_size * sizeof(int), OCL_MEM_READ_ONLY | OCL_MEM_COPY_HOST_PTR, s->krad);
	else
		g_ocl.d_krad = ocl_buf(4, rw, NULL);
	if(!g_ocl.d_krad) return -1;

	for(int p = 0; p < s->nplanes; p++) {
		const size_t n = (size_t)s->planes[p].in_w * s->planes[p].in_h;
		if(s->mode >= 6) {
			g_ocl.d_grad[p] = ocl_buf(n * 4 * sizeof(int16_t), rw, NULL);
			g_ocl.d_dist[p] = ocl_buf(n * sizeof(int16_t),     rw, NULL);
		} else {
			g_ocl.d_grad[p] = ocl_buf(4, rw, NULL);
			g_ocl.d_dist[p] = ocl_buf(4, rw, NULL);
		}
		if(!g_ocl.d_grad[p] || !g_ocl.d_dist[p]) return -1;
	}
	g_ocl.ready = 1;
	return 0;
}

// ---------- traitement d'une frame ----------
#define OCL_ARG(k, i, v) \
	do { ocl_int _e = g_ocl.api.SetKernelArg((k), (ocl_uint)(i), sizeof(v), &(v)); \
	     if(_e != 0) { fprintf(stderr, "OpenCL: clSetKernelArg(%d) failed (%d)\n", (int)(i), (int)_e); return -1; } } while(0)

static int ocl_run(ocl_kernel k, int dim, size_t gx, size_t gy)
{
	// Groupes 2D explicites (16x8 = 128 work-items, accepté par tous les GPU) :
	// des pixels voisins partagent les mêmes lignes source, donc mieux en cache
	// qu'une forme choisie par le pilote. Les noyaux testent déjà les bornes.
	static size_t lw = 16, lh = 8;
	static int lws_init = 0;
	if(!lws_init) {
		const char *env = getenv("TS_OCL_LWS");   // réglage avancé, ex. "32x4" ; "0" = pilote
		if(env) { unsigned w = 0, h = 0; if(sscanf(env, "%ux%u", &w, &h) == 2 && w && h) { lw = w; lh = h; } else lw = lh = 0; }
		lws_init = 1;
	}
	size_t g[2] = { gx, gy };
	size_t l[2] = { lw, lh };
	const size_t *lp = NULL;
	if(dim == 2 && lw && lh) {
		g[0] = (gx + lw - 1) / lw * lw;
		g[1] = (gy + lh - 1) / lh * lh;
		lp = l;
	}
	ocl_int e = g_ocl.api.EnqueueNDRangeKernel(g_ocl.q, k, (ocl_uint)dim, NULL, g, lp, 0, NULL, NULL);
	if(e != 0) { fprintf(stderr, "OpenCL: kernel launch failed (error %d)\n", (int)e); return -1; }
	return 0;
}

// Traite une frame d'entrée (in) avec son successeur (next, == in si aucun)
// et écrit les `ratio` frames de sortie dans out (même disposition que le CPU).
static double ocl_now(void)
{
#ifdef _WIN32
	LARGE_INTEGER f, c; QueryPerformanceFrequency(&f); QueryPerformanceCounter(&c);
	return (double)c.QuadPart / (double)f.QuadPart;
#else
	struct timespec t; clock_gettime(CLOCK_MONOTONIC, &t); return t.tv_sec + t.tv_nsec * 1e-9;
#endif
}
static double g_t_up, g_t_gpu, g_t_down; static int g_prof = -1;

static int ocl_submit_frame(const unsigned char *in, const unsigned char *next,
                             unsigned char *out, ocl_event *ev)
{
	if(!g_ocl.ready) return -1;
	OclApi *a = &g_ocl.api;
	const OclSetup *s = &g_ocl.s;
	const int has_next = (next != in);
	ocl_int e;
	if(g_prof < 0) g_prof = (getenv("TS_OCL_PROFILE") != NULL);
	double t0 = g_prof ? ocl_now() : 0, t1 = 0, t2 = 0;

	// La frame "next" de l'appel précédent est la frame "in" de celui-ci :
	// on échange les tampons au lieu de la renvoyer (en 4K : 12 Mo économisés
	// par image). "next" n'est envoyée que si un fondu temporel l'utilise.
	int need_next = 0;
	for(int r = 0; r < s->ratio; r++) if(s->tmp_weights[r] > 0) need_next = 1;

	if(g_ocl.next_valid && g_ocl.last_next == in) {
		ocl_mem t = g_ocl.d_in; g_ocl.d_in = g_ocl.d_next; g_ocl.d_next = t;
	} else {
		e = a->EnqueueWriteBuffer(g_ocl.q, g_ocl.d_in, 0, 0, s->in_len, in, 0, NULL, NULL);
		if(e != 0) OCL_FAIL("upload failed", e);
	}
	g_ocl.next_valid = 0;
	if(has_next && need_next) {
		e = a->EnqueueWriteBuffer(g_ocl.q, g_ocl.d_next, 0, 0, s->in_len, next, 0, NULL, NULL);
		if(e != 0) OCL_FAIL("upload failed", e);
		g_ocl.last_next = next;
		g_ocl.next_valid = 1;
	}

if(g_prof) { a->Finish(g_ocl.q); t1 = ocl_now(); }
	// gradient / confiance / couture : une fois par plan, sur la frame d'origine
	if(s->mode >= 6)
	{
		for(int p = 0; p < s->nplanes; p++)
		{
			const OclPlane *pl = &s->planes[p];
			ocl_uint off = (ocl_uint)pl->in_off;
			int w = pl->in_w, h = pl->in_h, bps = s->bps_in, ch = pl->chan;
			int mode = s->mode, afp = s->aniso_fp;
			ocl_mem G = g_ocl.d_grad[p], D = g_ocl.d_dist[p];

			OCL_ARG(g_ocl.k_grad, 0, g_ocl.d_in);
			OCL_ARG(g_ocl.k_grad, 1, off);
			OCL_ARG(g_ocl.k_grad, 2, G);
			OCL_ARG(g_ocl.k_grad, 3, w);
			OCL_ARG(g_ocl.k_grad, 4, h);
			OCL_ARG(g_ocl.k_grad, 5, bps);
			OCL_ARG(g_ocl.k_grad, 6, ch);
			if(ocl_run(g_ocl.k_grad, 2, (size_t)w, (size_t)h) != 0) return -1;

			OCL_ARG(g_ocl.k_conf, 0, G);
			OCL_ARG(g_ocl.k_conf, 1, w);
			OCL_ARG(g_ocl.k_conf, 2, h);
			OCL_ARG(g_ocl.k_conf, 3, mode);
			OCL_ARG(g_ocl.k_conf, 4, afp);
			if(ocl_run(g_ocl.k_conf, 2, (size_t)w, (size_t)h) != 0) return -1;

			OCL_ARG(g_ocl.k_dist, 0, G);
			OCL_ARG(g_ocl.k_dist, 1, D);
			OCL_ARG(g_ocl.k_dist, 2, w);
			OCL_ARG(g_ocl.k_dist, 3, h);
			if(ocl_run(g_ocl.k_dist, 2, (size_t)w, (size_t)h) != 0) return -1;
		}
	}

	for(int r = 0; r < s->ratio; r++)
	{
		const int tw = s->tmp_weights[r];
		ocl_mem src = g_ocl.d_in;

		if(tw > 0 && has_next)
		{
			ocl_uint n = (ocl_uint)(s->in_len / (size_t)s->bps_in);
			int bps = s->bps_in;
			OCL_ARG(g_ocl.k_blend, 0, g_ocl.d_in);
			OCL_ARG(g_ocl.k_blend, 1, g_ocl.d_next);
			OCL_ARG(g_ocl.k_blend, 2, g_ocl.d_blend);
			OCL_ARG(g_ocl.k_blend, 3, n);
			OCL_ARG(g_ocl.k_blend, 4, bps);
			OCL_ARG(g_ocl.k_blend, 5, tw);
			if(ocl_run(g_ocl.k_blend, 1, (size_t)n, 1) != 0) return -1;
			src = g_ocl.d_blend;
		}

		for(int p = 0; p < s->nplanes; p++)
		{
			const OclPlane *pl = &s->planes[p];
			ocl_uint soff = (ocl_uint)pl->in_off;
			ocl_uint doff = (ocl_uint)(pl->out_off + (size_t)r * pl->out_step);
			int wbase = r * pl->out_w, hbase = r * pl->out_h;
			int in_w = pl->in_w, in_h = pl->in_h, out_w = pl->out_w, out_h = pl->out_h;
			int ch = pl->chan, bi = s->bps_in, bo = s->bps_out, mode = s->mode, afp = s->aniso_fp;
			ocl_mem *mp = g_ocl.d_map[pl->map_idx];
			ocl_mem G = g_ocl.d_grad[p], D = g_ocl.d_dist[p];

			ocl_kernel K = g_ocl.use_w4 ? g_ocl.k_scale4 : g_ocl.k_scale;
			ocl_mem tA = g_ocl.use_w4 ? g_ocl.d_dense[pl->map_idx] : g_ocl.d_tab[pl->map_idx][0];
			ocl_mem tB = g_ocl.use_w4 ? g_ocl.d_tab[pl->map_idx][0] : g_ocl.d_tab[pl->map_idx][1];
			int i = 0;
			OCL_ARG(K, i, src);          i++;
			OCL_ARG(K, i, soff);         i++;
			OCL_ARG(K, i, g_ocl.d_out);  i++;
			OCL_ARG(K, i, doff);         i++;
			OCL_ARG(K, i, mp[0]);        i++;   // w_pos
			OCL_ARG(K, i, mp[1]);        i++;   // w_src
			OCL_ARG(K, i, mp[2]);        i++;   // w_flag
			OCL_ARG(K, i, mp[3]);        i++;   // h_pos
			OCL_ARG(K, i, mp[4]);        i++;   // h_src
			OCL_ARG(K, i, mp[5]);        i++;   // h_flag
			OCL_ARG(K, i, g_ocl.d_qx[pl->map_idx]); i++;
			OCL_ARG(K, i, g_ocl.d_qy[pl->map_idx]); i++;
			OCL_ARG(K, i, g_ocl.nfx[pl->map_idx]);  i++;
			OCL_ARG(K, i, wbase);        i++;
			OCL_ARG(K, i, hbase);        i++;
			OCL_ARG(K, i, in_w);         i++;
			OCL_ARG(K, i, in_h);         i++;
			OCL_ARG(K, i, out_w);        i++;
			OCL_ARG(K, i, out_h);        i++;
			OCL_ARG(K, i, tA);                       i++;
			OCL_ARG(K, i, tB);                       i++;
			OCL_ARG(K, i, g_ocl.d_tab[pl->map_idx][2]); i++;
			OCL_ARG(K, i, g_ocl.d_krad); i++;
			OCL_ARG(K, i, G);            i++;
			OCL_ARG(K, i, D);            i++;
			OCL_ARG(K, i, afp);          i++;
			if(ocl_run(K, 2, g_ocl.use_w4 ? (size_t)((out_w + 3) / 4) : (size_t)out_w, (size_t)out_h) != 0) return -1;
		}
	}

	// sortie 4:4:4 -> RGB24 sur le GPU (une passe par frame de sortie)
	ocl_mem final_buf = g_ocl.d_out;
	if(s->rgb_out)
	{
		final_buf = g_ocl.d_rgb;
		for(int r = 0; r < s->ratio; r++)
		{
			ocl_uint yo = (ocl_uint)s->planes[0].out_off, uo = (ocl_uint)s->planes[1].out_off;
			ocl_uint vo = (ocl_uint)s->planes[2].out_off, st = (ocl_uint)s->planes[0].out_step;
			int w = s->planes[0].out_w, h = s->planes[0].out_h, rr = r;
			ocl_mem srcb = g_ocl.d_out, dstb = g_ocl.d_rgb;
			OCL_ARG(g_ocl.k_yuv, 0, srcb);
			OCL_ARG(g_ocl.k_yuv, 1, dstb);
			OCL_ARG(g_ocl.k_yuv, 2, yo);
			OCL_ARG(g_ocl.k_yuv, 3, uo);
			OCL_ARG(g_ocl.k_yuv, 4, vo);
			OCL_ARG(g_ocl.k_yuv, 5, st);
			OCL_ARG(g_ocl.k_yuv, 6, w);
			OCL_ARG(g_ocl.k_yuv, 7, h);
			OCL_ARG(g_ocl.k_yuv, 8, rr);
			if(ocl_run(g_ocl.k_yuv, 2, (size_t)w, (size_t)h) != 0) return -1;
		}
	}

if(g_prof) { a->Finish(g_ocl.q); t2 = ocl_now(); }
	// Lecture asynchrone : on rend la main tout de suite pour que l'appelant
	// puisse soumettre la frame suivante avant que celle-ci soit terminée
	// (la file est ordonnée : calculs et copies s'enchaînent sans trou).
	// En mode profil, on reste synchrone pour mesurer chaque étape.
	*ev = NULL;
	if(g_prof) {
		e = a->EnqueueReadBuffer(g_ocl.q, final_buf, 1, 0, s->out_len, out, 0, NULL, NULL);
		if(e != 0) OCL_FAIL("download failed", e);
		g_t_up += t1 - t0; g_t_gpu += t2 - t1; g_t_down += ocl_now() - t2;
	} else {
		e = a->EnqueueReadBuffer(g_ocl.q, final_buf, 0, 0, s->out_len, out, 0, NULL, ev);
		if(e != 0) OCL_FAIL("download failed", e);
		a->Flush(g_ocl.q);
	}
	return 0;
}

// Attend la fin d'une frame soumise (ev peut être NULL : déjà terminée).
static int ocl_wait_frame(ocl_event ev)
{
	if(!ev) return 0;
	ocl_int e = g_ocl.api.WaitForEvents(1, &ev);
	g_ocl.api.ReleaseEvent(ev);
	if(e != 0) OCL_FAIL("wait failed", e);
	return 0;
}

// Version synchrone.
static int ocl_process_frame(const unsigned char *in, const unsigned char *next,
                             unsigned char *out)
{
	ocl_event ev;
	if(ocl_submit_frame(in, next, out, &ev) != 0) return -1;
	return ocl_wait_frame(ev);
}

static void ocl_release(void)
{
	if(!g_ocl.loaded) return;
	if(g_prof) fprintf(stderr, "OpenCL profile: upload %.3fs, compute %.3fs, download %.3fs\n", g_t_up, g_t_gpu, g_t_down);
	OclApi *a = &g_ocl.api;
	if(g_ocl.q) a->Finish(g_ocl.q);
	#define REL(m) do { if(m) { a->ReleaseMemObject(m); (m) = NULL; } } while(0)
	REL(g_ocl.d_in); REL(g_ocl.d_next); REL(g_ocl.d_blend); REL(g_ocl.d_out); REL(g_ocl.d_rgb);
	for(int i = 0; i < 6; i++) {
		if(g_ocl.d_map[1][i] == g_ocl.d_map[0][i]) g_ocl.d_map[1][i] = NULL;
		REL(g_ocl.d_map[0][i]); REL(g_ocl.d_map[1][i]);
	}
	for(int p = 0; p < 3; p++) { REL(g_ocl.d_grad[p]); REL(g_ocl.d_dist[p]); }
	for(int i = 0; i < 3; i++) if(g_ocl.d_tab[1][i] == g_ocl.d_tab[0][i]) g_ocl.d_tab[1][i] = NULL;
	if(g_ocl.d_qx[1] == g_ocl.d_qx[0]) g_ocl.d_qx[1] = NULL;
	if(g_ocl.d_qy[1] == g_ocl.d_qy[0]) g_ocl.d_qy[1] = NULL;
	for(int m = 0; m < 2; m++) {
	
	if(g_ocl.d_tab[m][1] == g_ocl.d_tab[m][0]) g_ocl.d_tab[m][1] = NULL;
	
	if(g_ocl.d_tab[m][2] == g_ocl.d_tab[m][0]) g_ocl.d_tab[m][2] = NULL;
	
	for(int i = 0; i < 3; i++) REL(g_ocl.d_tab[m][i]);
	
	REL(g_ocl.d_qx[m]); REL(g_ocl.d_qy[m]);
	}
	if(g_ocl.d_dense[1] == g_ocl.d_dense[0]) g_ocl.d_dense[1] = NULL;
	REL(g_ocl.d_dense[0]); REL(g_ocl.d_dense[1]);
	REL(g_ocl.d_krad);
	#undef REL
	#define RELK(k) do { if(k) { a->ReleaseKernel(k); (k) = NULL; } } while(0)
	RELK(g_ocl.k_blend); RELK(g_ocl.k_grad); RELK(g_ocl.k_conf); RELK(g_ocl.k_dist);
	RELK(g_ocl.k_scale); RELK(g_ocl.k_scale4); RELK(g_ocl.k_yuv);
	#undef RELK
	if(g_ocl.prog) { a->ReleaseProgram(g_ocl.prog);      g_ocl.prog = NULL; }
	if(g_ocl.q)    { a->ReleaseCommandQueue(g_ocl.q);    g_ocl.q    = NULL; }
	if(g_ocl.ctx)  { a->ReleaseContext(g_ocl.ctx);       g_ocl.ctx  = NULL; }
	// état réinitialisable : un nouvel ocl_init/ocl_setup peut suivre (changement de taille)
	g_ocl.next_valid = 0;
	g_ocl.last_next  = NULL;
	g_ocl.ready = 0;
}

#endif // TEMPORAL_SCALING_OCL_H
