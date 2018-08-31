// Onur G. Guleryuz 1995, 1996, 1997,
// University of Illinois at Urbana-Champaign,
// Princeton University,
// Polytechnic University.

#include <opencv2\opencv.hpp>
#include <stdio.h>
#include <math.h>
#include "macros.h"
#include "util.h"
#include "wav_filters.h"
#include "wav_trf.h"
#include "wav_gen.h"

/***********************************************************/
// Basic operations for thresholding and clipping.

// round and clip the pixels to comply with a grayscale image.
void round_and_clip_to_0_255(float **im, int N1, int N2)

{
	int i, j;

	for (i = 0; i < N1; i++)
		for (j = 0; j < N2; j++) {

			if (im[i][j]>255)
				im[i][j] = 255;
			else if (im[i][j] < 0)
				im[i][j] = 0;
			im[i][j] = (float)((int)(im[i][j] + .5));
		}
}

/***********************************************************/
// Example transform/inverse-transform routines.
/*
Ni:rows
Nj:cols
*/

void wavelet_transform(float **buffer_t, int Ni, int Nj, int noInverse)
{
	int i;
	int shift_arr_r[MAX_ARR_SIZE], shift_arr_c[MAX_ARR_SIZE];

	// 4 levels of wavelets.
	const int levs = 4;
	int Nl, Nh;
	float *lp, *hp;
	float mse = 0;

	// Choose wavelet filter.
	choose_filter('D', 9);

	for (i = levs - 1; i >= 0; i--) {

		shift_arr_r[i] = shift_arr_c[i] = 0;
	}

	if (noInverse){

		// Select the forward bank of filters.
		lp = MFLP; Nl = Nflp;
		hp = MFHP; Nh = Nfhp;

		wavpack2d_inpl(buffer_t, Ni, Nj, levs, lp, Nl, hp, Nh, 1, shift_arr_r, shift_arr_c);
		//wav2d_inpl(buffer_t, Ni, Nj, levs, lp, Nl, hp, Nh, 1, shift_arr_r, shift_arr_c);
	}
	else{
		// Select the inverse bank of filters.
		lp = MILP; Nl = Nilp;
		hp = MIHP; Nh = Nihp;

		// Inverse transform.
		wavpack2d_inpl(buffer_t, Ni, Nj, levs, lp, Nl, hp, Nh, 0, shift_arr_r, shift_arr_c);
		//wav2d_inpl(buffer_t, Ni, Nj, levs, lp, Nl, hp, Nh, 0, shift_arr_r, shift_arr_c);
	}

}

/***********************************************************/