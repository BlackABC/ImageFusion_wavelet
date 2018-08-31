// Onur G. Guleryuz 1995, 1996, 1997,
// University of Illinois at Urbana-Champaign,
// Princeton University,
// Polytechnic University.

#include "wav_basic.h"
#include "wav_filters_extern.h"
#include "alloc.h"
#include "wav_gen.h"

// Macro that sets filter parameters.
// Called via choose_filter().
#define set_filter_param(fname) \
{ \
	MFLP= fname ## _FLP; \
	MFHP= fname ## _FHP; \
	MILP= fname ## _ILP; \
	MIHP= fname ## _IHP; \
	\
	/* Filter size. */ \
	Nflp= fname ## _Nflp; \
	Nfhp= fname ## _Nfhp; \
	Nilp= fname ## _Nilp; \
	Nihp= fname ## _Nihp; \
	\
	/* Filter symmetry */ \
	PS= fname ## _PS; \
	\
	/* Filter starting points. */ \
	begflp=mfl= fname ## _begflp; \
	begfhp=mfh= fname ## _begfhp; \
	begilp=mil= fname ## _begilp; \
	begihp=mih= fname ## _begihp; \
}

// Set wavelet filter parameters depending on name
// and tap. See wav_filters.h and wav_filters_extern.h 
// This routine must be called before doing any transforms
// to select which wavelet bank to use.
void choose_filter(char name,int tap)

{
	switch (name) {
		case 'A':
			switch (tap) {
				default:
					set_filter_param(Antonini);
			}
			break;

		case 'H':
			switch (tap) {
				default:
					set_filter_param(Haar);
			}
			break;

		case 'N':
			switch (tap) {
				case -1:
					set_filter_param(Nick_flip);
					break;

				default:
					set_filter_param(Nick);
			}
			break;

		default:
			switch (tap) {
				case 9:
					set_filter_param(D79);
					break;

				default:
					set_filter_param(D8);
			}
	}
}

// 2D "in place" wavelet transform. Input data is in
// float **image and output is returned inside image.
//
// levs>=1 is the number of wavelet levels.
//
// Ni and Nj are the dimensions of the image.
// Ni and Nj should be divisible by 2^levs.
//
// shift_arr_row and shift_arr_col (each of dimension levs)
// contain the amount of shift (0 or 1)
// that should be applied to filters in each level for transforms over rows and columns.
// The use of shift_arr_row and shift_arr_col is mainly intended for 
// overcomplete filtering/transforms.
//
// forw != 0 for a forward transform (inverse is calculated otherwise).
//
// float *lp of dimension Nl contains the low pass wavelet filter.
// float *hp of dimension Nh contains the high pass wavelet filter.
//
void wav2d_inpl(float **image,int Ni,int Nj,int levs,float *lp,int Nl,
						float *hp,int Nh,char forw,int *shift_arr_row,int *shift_arr_col)

{
	int i,dimi,dimj;

	dimi=(forw)?Ni:(Ni>>(levs-1));
	dimj=(forw)?Nj:(Nj>>(levs-1));

	for(i=0;i<levs;i++) {

		if(forw) {

			// Do the rows.

			// Adjust filter parameters for overcomplete filtering.
			begflp=mfl+shift_arr_row[i];   
			begfhp=mfh-shift_arr_row[i];

			filt_n_dec_all_rows(image,dimi,dimj,lp,Nl,hp,Nh);

			// Now the columns.

			begflp=mfl+shift_arr_col[i];   
			begfhp=mfh-shift_arr_col[i];

			filt_n_dec_all_cols(image,dimi,dimj,lp,Nl,hp,Nh);
		}
		else {

			// Rows.
			begilp=mil+shift_arr_row[levs-1-i];
			begihp=mih-shift_arr_row[levs-1-i];

			ups_n_filt_all_rows(image,dimi,dimj,lp,Nl,hp,Nh,levs-1-i,shift_arr_row);

			// Columns.
			begilp=mil+shift_arr_col[levs-1-i];
			begihp=mih-shift_arr_col[levs-1-i];

			ups_n_filt_all_cols(image,dimi,dimj,lp,Nl,hp,Nh,levs-1-i,shift_arr_col);
		}

		dimi=(forw)?(dimi>>1):(dimi<<1);
		dimj=(forw)?(dimj>>1):(dimj<<1);
	}
}



// 2-D "in place" wavelet packet transform.
// See wav2d_inpl() for an explanation of parameters.
//
// This routine will grow a full packet tree for levels upto LL_ONLY_AFTER_LEV 
// and just a wavelet tree afterward, i.e., if levs>LL_ONLY_AFTER_LEV
// we grow wavelets only over LL band after LL_ONLY_AFTER_LEV.
// LL_ONLY_AFTER_LEV is set in wav_gen.h Just set it to some large number
// if you want packets all the way.
// 
// Ni and Nj should be divisible by 2^levs.
// shift_arr_row and shift_arr_col contain the amount of shift (0 or 1)
// that should be applied to filters in each level.
// This is mainly intended for overcomplete filtering.
// See wav2d_inpl() comments.
//
void wavpack2d_inpl(float **image,int Ni,int Nj,int levs,float *lp,int Nl,
						float *hp,int Nh,char forw,int *shift_arr_row,int *shift_arr_col)
{
	int i,k,l,dimi,dimj;
	int p,q,rcnt,ccnt;
	float **buffer;

	buffer=allocate_2d_float(Ni,Nj,0);

	dimi=(forw)?Ni:(Ni>>(levs-1));
	dimj=(forw)?Nj:(Nj>>(levs-1));

	for(i=0;i<levs;i++) {

		// Number of row bands.
		rcnt=Ni/dimi;

		// Number of column bands.
		ccnt=Nj/dimj;

		if(forw&&(i>=LL_ONLY_AFTER_LEV)) {
			// Grow wavelets only over LL band after LL_ONLY_AFTER_LEV.
			rcnt=ccnt=1;
		}
		else if((!forw)&&(levs-1-i>=LL_ONLY_AFTER_LEV)) {
			// Inverse wavelets only over LL band after LL_ONLY_AFTER_LEV.
			rcnt=ccnt=1;
		}

		// Do the next level of the wavelet tree over all of
		// the bands we are supposed to grow on.
		for(p=0;p<rcnt;p++)
			for(q=0;q<ccnt;q++) {

				for(k=0;k<dimi;k++)
					for(l=0;l<dimj;l++)
						buffer[k][l]=image[p*dimi+k][q*dimj+l];

				if(forw) {

					// Do the rows.

					// Adjust filter parameters for overcomplete filtering.
					begflp=mfl+shift_arr_row[i];   
					begfhp=mfh-shift_arr_row[i];

					filt_n_dec_all_rows(buffer,dimi,dimj,lp,Nl,hp,Nh);

					// Now the columns.
					begflp=mfl+shift_arr_col[i];   
					begfhp=mfh-shift_arr_col[i];

					filt_n_dec_all_cols(buffer,dimi,dimj,lp,Nl,hp,Nh);
				}
				else {

					// Rows.
					begilp=mil+shift_arr_row[levs-1-i];
					begihp=mih-shift_arr_row[levs-1-i];

					ups_n_filt_all_rows(buffer,dimi,dimj,lp,Nl,hp,Nh,levs-1-i,shift_arr_row);

					// Columns.
					begilp=mil+shift_arr_col[levs-1-i];
					begihp=mih-shift_arr_col[levs-1-i];

					ups_n_filt_all_cols(buffer,dimi,dimj,lp,Nl,hp,Nh,levs-1-i,shift_arr_col);
				}

				for(k=0;k<dimi;k++)
					for(l=0;l<dimj;l++)
						image[p*dimi+k][q*dimj+l]=buffer[k][l];
			}

		dimi=(forw)?(dimi>>1):(dimi<<1);
		dimj=(forw)?(dimj>>1):(dimj<<1);
	}

	free_2d_float(buffer,Ni);
}

// Forward complex wavelet transform using Kingsbury's paper.
// Please check with Nick Kingsbury's paper to confirm that
// I am doing things correctly here.
//
// Input is in float **im, and the results are returned in 4
// **float's, trf[0],...,trf[3].
//
// Ni, Nj multiples of 2^levs.
// levs>=1.
// See wav2d_inpl() for an explanation of parameters.
//
void complex_wav_forw(float **im,float ***trf,int Ni,int Nj,int levs)

{
	int i,j,k,l;
	int Nl,Nh,dimi,dimj,hdimi,hdimj;
	float *lp,*hp;
	int shift_arr_r[MAX_ARR_SIZE],shift_arr_c[MAX_ARR_SIZE];
	float tmp;

	// Choose Antonini.
	choose_filter('A',0);

	lp=MFLP;Nl=Nflp;  
	hp=MFHP;Nh=Nfhp;

	// First level, fully overcomplete resulting 4x redundancy.
	// Results are in trf[i].
	for(i=0;i<4;i++) {

		for(k=0;k<Ni;k++)
			for(l=0;l<Nj;l++)
				trf[i][k][l]=im[k][l];

		// Set shifts as specified by Kingsbury's paper.
		if(i/2==0)
			shift_arr_r[0]=1;
		else
			shift_arr_r[0]=0;

		if(i%2==0)
			shift_arr_c[0]=1;
		else
			shift_arr_c[0]=0;

		// Transform.
		wav2d_inpl(trf[i],Ni,Nj,1,lp,Nl,hp,Nh,1,shift_arr_r,shift_arr_c);
	}

	// Now grow complex wavelets over the LL band
	// using Kingsbury's orthogonal bank.
	if(levs>1) {

		for(i=0;i<4;i++) {

			dimi=Ni/2;
			dimj=Nj/2;
			for(j=1;j<levs;j++) {
	
				if(i/2==0) {
					// Regular orth. bank.
					choose_filter('N',0);
				}
				else {
					// Flipped orth. bank (see Kingsbury).
					// Flipped is same as inverse.
					choose_filter('N',-1);
				}
		
				lp=MFLP;Nl=Nflp;  
				hp=MFHP;Nh=Nfhp;

				// Transform over rows.
				filt_n_dec_all_rows(trf[i],dimi,dimj,lp,Nl,hp,Nh);
	
				if(i%2==0) {
					// Regular.
					choose_filter('N',0);
				}
				else {
					// Flipped.
					choose_filter('N',-1);
				}
	
				lp=MFLP;Nl=Nflp;  
				hp=MFHP;Nh=Nfhp;

				// Transform over columns.
				filt_n_dec_all_cols(trf[i],dimi,dimj,lp,Nl,hp,Nh);
		
				dimi>>=1;
				dimj>>=1;
		
			}
		}
	}

	// Generate the complex wavelet coefficients.
	hdimi=Ni/(1<<levs);
	hdimj=Nj/(1<<levs);

	for(k=0;k<Ni;k++)
		for(l=0;l<Nj;l++) {

			// After this computation the complex wavelet coefficients
			// are given by trf[0]+jtrf[2] and trf[3]+jtrf[1];
			tmp=trf[0][k][l]+trf[3][k][l];
			trf[3][k][l]=trf[0][k][l]-trf[3][k][l];
			trf[0][k][l]=tmp;

			tmp=trf[1][k][l]+trf[2][k][l];
			// It is possible that the below is 2 - 1.
			trf[2][k][l]=trf[1][k][l]-trf[2][k][l];
			trf[1][k][l]=tmp;
		}
}

// Inverse complex wavelet transform.
// trf is modified in intermediate computations.
// Ni, Nj multiples of 2^levs
// levs>=1.
//
// See complex_wav_forw() comments.
//
float **complex_wav_inv(float ***trf,int Ni,int Nj,int levs)

{
	int i,j,k,l;
	int Nl,Nh,dimi,dimj,hdimi,hdimj;
	float *lp,*hp;
	int shift_arr_r[MAX_ARR_SIZE],shift_arr_c[MAX_ARR_SIZE];
	float tmp;
	float **im;


	// Clean the filter shift arrays.
	for(i=0;i<1<<levs;i++)
		shift_arr_r[i]=shift_arr_c[i]=0;

	// Get back the intermediate wavelet coefficients from the
	// complex wavelet coefficients.
	hdimi=Ni/(1<<levs);
	hdimj=Nj/(1<<levs);

	for(k=0;k<Ni;k++)
		for(l=0;l<Nj;l++) {

			// Prior to this computation the complex wavelet coefficients
			// are given by trf[0]+jtrf[2] and trf[3]+jtrf[1];
			tmp=(trf[0][k][l]+trf[3][k][l])/2;
			trf[3][k][l]=(trf[0][k][l]-trf[3][k][l])/2;
			trf[0][k][l]=tmp;

			tmp=(trf[1][k][l]+trf[2][k][l])/2;
			trf[2][k][l]=(trf[1][k][l]-trf[2][k][l])/2;
			trf[1][k][l]=tmp;
		}

	if(levs>1) {

		for(i=0;i<4;i++) {
	
			dimi=Ni>>(levs-1);
			dimj=Nj>>(levs-1);
	
				
			for(j=1;j<levs;j++) {
	
				if(i/2==0) {
					// Regular inverse, it just happens that orth. inverses are flipped.
					choose_filter('N',0);
				}
				else {
					// Flipped inverse.
					choose_filter('N',-1);
				}

				lp=MILP;Nl=Nilp;  
				hp=MIHP;Nh=Nihp;

				// Inverse transform over rows.
				ups_n_filt_all_rows(trf[i],dimi,dimj,lp,Nl,hp,Nh,levs-1-j,shift_arr_r);
		
	
				if(i%2==0) {
					// Regular inverse.
					choose_filter('N',0);
				}
				else {
					// Flipped inverse.
					choose_filter('N',-1);
				}
	
				lp=MILP;Nl=Nilp;  
				hp=MIHP;Nh=Nihp;

				// Inverse transform over columns.
				ups_n_filt_all_cols(trf[i],dimi,dimj,lp,Nl,hp,Nh,levs-1-j,shift_arr_c);
		
				dimi<<=1;
				dimj<<=1;
			}
		}
	}

	// Initialize final result with zeros.
	im=allocate_2d_float(Ni,Nj,1);

	// Choose Antonini.
	choose_filter('A',0);

	lp=MILP;Nl=Nilp;  
	hp=MIHP;Nh=Nihp;

	// Now the first level.
	for(i=0;i<4;i++) {

		// Set shifts as specified by Kingsbury's paper.
		if(i/2==0)
			shift_arr_r[0]=1;
		else
			shift_arr_r[0]=0;

		if(i%2==0)
			shift_arr_c[0]=1;
		else
			shift_arr_c[0]=0;

		// Inverse transform.
		wav2d_inpl(trf[i],Ni,Nj,1,lp,Nl,hp,Nh,0,shift_arr_r,shift_arr_c);

		// Accumulate the results of 4x redundancy.
		for(k=0;k<Ni;k++)
			for(l=0;l<Nj;l++)
				im[k][l]+=trf[i][k][l];

	}

	// Normalize results by 4 to account for 4x redundancy accumulation.
	for(k=0;k<Ni;k++)
		for(l=0;l<Nj;l++)
			im[k][l]/=4;

	return(im);
}


// Forward complex wavelet packet transform.
// Ni, Nj multiples of 2^levs.
// levs>=1.
//
// See complex_wav_forw() comments.
// See wavpack2d_inpl() comments for LL_ONLY_AFTER_LEV.
//
// Not fully tested. Please check with literature.
//
void complex_wav_pack_forw(float **im,float ***trf,int Ni,int Nj,int levs)

{
	int i,j,k,l;
	int Nl,Nh,dimi,dimj,hdimi,hdimj;
	float *lp,*hp;
	int shift_arr_r[MAX_ARR_SIZE],shift_arr_c[MAX_ARR_SIZE];
	float tmp,**buffer;
	int p,q,rcnt,ccnt;

	// Choose Antonini.
	choose_filter('A',0);

	lp=MFLP;Nl=Nflp;  
	hp=MFHP;Nh=Nfhp;

	// First level, fully overcomplete resulting 4x redundancy.
	// Results are in trf[i].
	for(i=0;i<4;i++) {

		for(k=0;k<Ni;k++)
			for(l=0;l<Nj;l++)
				trf[i][k][l]=im[k][l];

		// Set shifts as specified by Kingsbury's paper.
		if(i/2==0)
			shift_arr_r[0]=1;
		else
			shift_arr_r[0]=0;

		if(i%2==0)
			shift_arr_c[0]=1;
		else
			shift_arr_c[0]=0;

		// Transform.
		wav2d_inpl(trf[i],Ni,Nj,1,lp,Nl,hp,Nh,1,shift_arr_r,shift_arr_c);
	}

	// Now grow using Kingsbury's orthogonal bank.
	if(levs>1) {

		buffer=allocate_2d_float(Ni/2,Nj/2,0);

		for(i=0;i<4;i++) {

			dimi=Ni/2;
			dimj=Nj/2;
			for(j=1;j<levs;j++) {

				// Number of row/column bands.
				rcnt=Ni/dimi;
				ccnt=Nj/dimj;

				if(j>=LL_ONLY_AFTER_LEV) {
					// Grow wavelets only over LL band after LL_ONLY_AFTER_LEV.
					rcnt=ccnt=1;
				}
	
				for(p=0;p<rcnt;p++)
					for(q=0;q<ccnt;q++) {

						// Copy data to buffer.
						for(k=0;k<dimi;k++)
							for(l=0;l<dimj;l++)
								buffer[k][l]=trf[i][p*dimi+k][q*dimj+l];
						
						if(i/2==0) {
							// Regular orth. bank.
							choose_filter('N',0);
						}
						else {
							// Flipped orth. bank (see Kingsbury).
							// Flipped is same as inverse.
							choose_filter('N',-1);
						}
				
						lp=MFLP;Nl=Nflp;  
						hp=MFHP;Nh=Nfhp;

						// Transform over rows.
						filt_n_dec_all_rows(buffer,dimi,dimj,lp,Nl,hp,Nh);
			
						if(i%2==0) {
							// Regular.
							choose_filter('N',0);
						}
						else {
							// Flipped.
							choose_filter('N',-1);
						}
			
						lp=MFLP;Nl=Nflp;  
						hp=MFHP;Nh=Nfhp;

						// Transform over columns.
						filt_n_dec_all_cols(buffer,dimi,dimj,lp,Nl,hp,Nh);

						// Copy results to trf[i].
						for(k=0;k<dimi;k++)
							for(l=0;l<dimj;l++)
								trf[i][p*dimi+k][q*dimj+l]=buffer[k][l];
					}
		
				dimi>>=1;
				dimj>>=1;
		
			}
		}
		free_2d_float(buffer,Ni/2);
	}

	// Generate the complex wavelet coefficients.
	hdimi=Ni/(1<<levs);
	hdimj=Nj/(1<<levs);

	for(k=0;k<Ni;k++)
		for(l=0;l<Nj;l++) {

			// After this computation the complex wavelet coefficients
			// are given by trf[0]+jtrf[2] and trf[3]+jtrf[1];
			tmp=trf[0][k][l]+trf[3][k][l];
			trf[3][k][l]=trf[0][k][l]-trf[3][k][l];
			trf[0][k][l]=tmp;

			tmp=trf[1][k][l]+trf[2][k][l];
			trf[2][k][l]=trf[1][k][l]-trf[2][k][l];
			trf[1][k][l]=tmp;
		}
}

// Inverse complex wavelet packet transform.
// trf is modified in intermediate computations.
// Ni, Nj multiples of 2^levs
// levs>=1.
//
// See complex_wav_pack_forw() comments.
//
float **complex_wav_pack_inv(float ***trf,int Ni,int Nj,int levs)

{
	int i,j,k,l;
	int Nl,Nh,dimi,dimj,hdimi,hdimj;
	float *lp,*hp;
	int shift_arr_r[MAX_ARR_SIZE],shift_arr_c[MAX_ARR_SIZE];
	float tmp;
	float **im,**buffer;
	int p,q,rcnt,ccnt;


	// Clean the filter shift arrays.
	for(i=0;i<1<<levs;i++)
		shift_arr_r[i]=shift_arr_c[i]=0;

	// Get back the intermediate wavelet coefficients from the
	// complex wavelet coefficients.
	hdimi=Ni/(1<<levs);
	hdimj=Nj/(1<<levs);

	for(k=0;k<Ni;k++)
		for(l=0;l<Nj;l++) {

			// Prior to this computation the complex wavelet coefficients
			// are given by trf[0]+jtrf[2] and trf[3]+jtrf[1];
			tmp=(trf[0][k][l]+trf[3][k][l])/2;
			trf[3][k][l]=(trf[0][k][l]-trf[3][k][l])/2;
			trf[0][k][l]=tmp;

			tmp=(trf[1][k][l]+trf[2][k][l])/2;
			trf[2][k][l]=(trf[1][k][l]-trf[2][k][l])/2;
			trf[1][k][l]=tmp;
		}

	if(levs>1) {

		buffer=allocate_2d_float(Ni/2,Nj/2,0);
		for(i=0;i<4;i++) {
	
			dimi=Ni>>(levs-1);
			dimj=Nj>>(levs-1);
	
				
			for(j=1;j<levs;j++) {
	
				// Number of row/column bands.
				rcnt=Ni/dimi;
				ccnt=Nj/dimj;

				if(j>=LL_ONLY_AFTER_LEV) {
					// Grow wavelets only over LL band after LL_ONLY_AFTER_LEV.
					rcnt=ccnt=1;
				}

				for(p=0;p<rcnt;p++)
					for(q=0;q<ccnt;q++) {

						// Copy data to buffer.
						for(k=0;k<dimi;k++)
							for(l=0;l<dimj;l++)
								buffer[k][l]=trf[i][p*dimi+k][q*dimj+l];

						if(i/2==0) {
							// Regular inverse, it just happens that orth. inverses are flipped.
							choose_filter('N',0);
						}
						else {
							// Flipped inverse.
							choose_filter('N',-1);
						}

						lp=MILP;Nl=Nilp;  
						hp=MIHP;Nh=Nihp;

						// Inverse transform over rows.
						ups_n_filt_all_rows(buffer,dimi,dimj,lp,Nl,hp,Nh,levs-1-j,shift_arr_r);
				
			
						if(i%2==0) {
							// Regular inverse.
							choose_filter('N',0);
						}
						else {
							// Flipped inverse.
							choose_filter('N',-1);
						}
			
						lp=MILP;Nl=Nilp;  
						hp=MIHP;Nh=Nihp;

						// Inverse transform over columns.
						ups_n_filt_all_cols(buffer,dimi,dimj,lp,Nl,hp,Nh,levs-1-j,shift_arr_c);

						// Copy results to trf[i].
						for(k=0;k<dimi;k++)
							for(l=0;l<dimj;l++)
								trf[i][p*dimi+k][q*dimj+l]=buffer[k][l];
					}
		
				dimi<<=1;
				dimj<<=1;
			}
		}

		free_2d_float(buffer,Ni/2);
	}

	// Initialize final result with zeros.
	im=allocate_2d_float(Ni,Nj,1);

	// Choose Antonini.
	choose_filter('A',0);

	lp=MILP;Nl=Nilp;  
	hp=MIHP;Nh=Nihp;

	// Now the first level.
	for(i=0;i<4;i++) {

		// Set shifts as specified by Kingsbury's paper.
		if(i/2==0)
			shift_arr_r[0]=1;
		else
			shift_arr_r[0]=0;

		if(i%2==0)
			shift_arr_c[0]=1;
		else
			shift_arr_c[0]=0;

		// Inverse transform.
		wav2d_inpl(trf[i],Ni,Nj,1,lp,Nl,hp,Nh,0,shift_arr_r,shift_arr_c);

		// Accumulate the results of 4x redundancy.
		for(k=0;k<Ni;k++)
			for(l=0;l<Nj;l++)
				im[k][l]+=trf[i][k][l];

	}

	// Normalize results by 4 to account for 4x redundancy accumulation.
	for(k=0;k<Ni;k++)
		for(l=0;l<Nj;l++)
			im[k][l]/=4;

	return(im);
}
