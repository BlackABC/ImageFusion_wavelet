//// Onur G. Guleryuz 1995, 1996, 1997,
//// University of Illinois at Urbana-Champaign,
//// Princeton University,
//// Polytechnic University.
//
//#include <stdio.h>
//#include <math.h>
//#include "macros.h"
//#include "alloc.h"
//#include "wav_filters.h"
//#include "wav_trf.h"
//#include "wav_gen.h"
//
//// global variable on whether we should grow packets or not.
//int PACKET=0;
//
///***********************************************************/
//// file io related
//
//FILE *myfopen(char *str,char *fmt)
//
//{
//	FILE *res;
//	
//	res=fopen(str,fmt);
//	check_ptr(res,"myfopen");
//	return(res);
//}
//
//float **load_image(char *name,int N1,int N2)
//
//{
//	int i,j;
//	float **res;
//	FILE *mf;
//	
//	res=allocate_2d_float(N1,N2,0);
//	mf=myfopen(name,"rb");
//	for(i=0;i<N1;i++)
//		for(j=0;j<N2;j++) {
//			if(feof(mf)) {
//				printf("eof reached in loading %4d x%4d image\n",N1,N2);
//				printf("at row %d and column %d\n",i,j);
//				printf("perhaps you have the dimensions wrong?\n");
//				exit(1);
//			}
//			else
//				res[i][j]=(float)fgetc(mf);
//		}
//	fclose(mf);
//	return(res);
//}
//
//
///***********************************************************/
//// Basic operations for thresholding and clipping.
//
//// round and clip the pixels to comply with a grayscale image.
//void round_and_clip_to_0_255(float **im,int N1,int N2)
//
//{
//	int i,j;
//
//	for(i=0;i<N1;i++)
//		for(j=0;j<N2;j++) {
//
//			if(im[i][j]>255)
//				im[i][j]=255;
//			else if(im[i][j]<0)
//				im[i][j]=0;
//			im[i][j]=(float)((int)(im[i][j]+.5));
//		}
//}
//
//
///***********************************************************/
//// Example transform/inverse-transform routines.
//
//void complx_wavelet_transform_test(float **original,int Ni,int Nj)
//
//{
//	int i,j;
//	float **buffer_t;
//	float **ctrf[4],**cim;
//	float mse=0;
//	const int levs=3;
//
//	buffer_t=allocate_2d_float(Ni,Nj,0);
//
//	for(i=0;i<4;i++) 
//		ctrf[i]=allocate_2d_float(Ni,Nj,0);
//
//	// Copy "original" input to the buffer.
//	for(i=0;i<Ni;i++)
//		for(j=0;j<Nj;j++)
//			buffer_t[i][j]=original[i][j];
//
//	if(PACKET) {
//		// Forward packet transform.
//		// Allocates components of ctrf.
//		complex_wav_pack_forw(buffer_t,ctrf,Ni,Nj,levs);
//
//		// Inverse transform.
//		// Allocates the temporary result array cim.
//		cim=complex_wav_pack_inv(ctrf,Ni,Nj,levs);
//	}
//	else {
//		// Forward transform.
//		// Allocates components of ctrf.
//		complex_wav_forw(buffer_t,ctrf,Ni,Nj,levs);
//
//		// Inverse transform.
//		// Allocates the temporary result array cim.
//		cim=complex_wav_inv(ctrf,Ni,Nj,levs);
//	}
//
//	mse=0;
//	for(i=0;i<Ni;i++)
//		for(j=0;j<Nj;j++) {
//
//			mse+=(original[i][j]-cim[i][j])*(original[i][j]-cim[i][j]);
//		}
//	mse/=Ni*Nj;
//
//	printf("\nThe mse due to transformation and inverse transformation is: \n");
//	printf("\t%8.4f, ( %8.4f dB)\n ",mse,10*log10(255*255/mse));
//	printf("(Very small amount of mse is OK due to numerical errors.)\n");
//
//	// Now round to see if we are getting perfect reconstruction.
//	round_and_clip_to_0_255(cim,Ni,Nj);
//
//	mse=0;
//	for(i=0;i<Ni;i++) {
//
//		for(j=0;j<Nj;j++) {
//
//			mse+=(original[i][j]-cim[i][j])*(original[i][j]-cim[i][j]);
//		}
//	}
//	
//	if(mse!=0) {
//
//		printf("\nmse after rounding is non-zero!\n");
//	}
//	else {
//
//		printf("\nmse after rounding is zero.\n");
//	}
//
//	// Free ctrf components allocated.
//	for(i=0;i<4;i++)
//		free_2d_float(ctrf[i],Ni);
//
//	// Free buffer.
//	free_2d_float(buffer_t,Ni);
//
//	// Free temporary storage for results.
//	free_2d_float(cim,Ni);
//}
//
//
//void wavelet_transform_test(float **original,int Ni,int Nj)
//
//{
//	int i,j;
//	float **buffer_t;
//	int shift_arr_r[MAX_ARR_SIZE],shift_arr_c[MAX_ARR_SIZE];
//
//	// 4 levels of wavelets.
//	const int levs=4;
//	int Nl,Nh;
//	float *lp,*hp;
//	float mse=0;
//
//	// Choose wavelet filter.
//	choose_filter('D',9);
//
//	// Main buffer for operations.
//	buffer_t=allocate_2d_float(Ni,Nj,0);
//
//	// At every level of a 1-D wavelet transform	there are two transforms possible
//	// obtained by the original bank or its one pixel shifted version.
//	// In L levels there are 2^L possibilities with each of these possibilities
//	// forming a complete transform (i.e., invertible, etc.). If we evaluate
//	// more than one then what we end up with is an overcomplete transform.
//	// In 2-D there are 4 possibilities at each level and hence 4^L total
//	// possibilities.
//	//
//	// At level l, this code selects either the "regular" (shift_arr_r[l]=0)
//	// or the one shifted bank (shift_arr_r[l]=1) over rows and likewise
//	// over columns using shift_arr_c[l].
//	//
//	// In order to generate the fully overcomplete transform we effectively
//	// put the binary representation of numbers from 0 to 2^L-1 (represented
//	// by oc_max) in shift_arr_r and shift_arr_c to traverse all 4^L possibilities.
//	//
//	// Here no shifts necessary since we are computing a single complete transform
//	// (the "regular" case). See wavelet_denoise() in test_denoise.c for the use of shifts in
//	// overcomplete transformations.
//	for(i=levs-1;i>=0;i--) {
//
//		shift_arr_r[i]=shift_arr_c[i]=0;
//	}
//
//			
//
//	for(i=0;i<Ni;i++) {
//		
//		for(j=0;j<Nj;j++) {
//			
//			buffer_t[i][j]=original[i][j];
//		}
//	}
//				
//			
//	// Select the forward bank of filters.
//	lp=MFLP;Nl=Nflp;  
//	hp=MFHP;Nh=Nfhp;
//				
//	if(PACKET) {
//		
//		// Packet transform.
//		wavpack2d_inpl(buffer_t,Ni,Nj,levs,lp,Nl,hp,Nh,1,shift_arr_r,shift_arr_c);
//	}
//	else {
//		
//		// Regular wavelets.
//		wav2d_inpl(buffer_t,Ni,Nj,levs,lp,Nl,hp,Nh,1,shift_arr_r,shift_arr_c);
//	}
//	
//	// Select the inverse bank of filters.
//	lp=MILP;Nl=Nilp;  
//	hp=MIHP;Nh=Nihp;
//				
//	// Inverse transform.
//	if(PACKET)
//		wavpack2d_inpl(buffer_t,Ni,Nj,levs,lp,Nl,hp,Nh,0,shift_arr_r,shift_arr_c);
//	else
//		wav2d_inpl(buffer_t,Ni,Nj,levs,lp,Nl,hp,Nh,0,shift_arr_r,shift_arr_c);
//				
//			
//	mse=0;
//	for(i=0;i<Ni;i++)
//		for(j=0;j<Nj;j++) {
//
//			mse+=(original[i][j]-buffer_t[i][j])*(original[i][j]-buffer_t[i][j]);
//		}
//	mse/=Ni*Nj;
//
//	printf("\nThe mse due to transformation and inverse transformation is: \n");
//	printf("\t%8.4f, ( %8.4f dB)\n ",mse,10*log10(255*255/mse));
//	printf("(Very small amount of mse is OK due to numerical errors.)\n");
//
//	// Now round to see if we are getting perfect reconstruction.
//	round_and_clip_to_0_255(buffer_t,Ni,Nj);
//
//	mse=0;
//	for(i=0;i<Ni;i++) {
//
//		for(j=0;j<Nj;j++) {
//
//			mse+=(original[i][j]-buffer_t[i][j])*(original[i][j]-buffer_t[i][j]);
//		}
//	}
//	
//	if(mse!=0) {
//
//		printf("\nmse after rounding is non-zero!\n");
//	}
//	else {
//
//		printf("\nmse after rounding is zero.\n");
//	}
//
//	// Free buffers.
//	free_2d_float(buffer_t,Ni);
//}
//
///***********************************************************/
//
//int main(int argc,char **argv)
//
//{
//	float **original;
//	const int Ni=512,Nj=512;
//
//	// This is some version of peppers that I lifted off the Internet.
//	// I don't know if it is the "standard" peppers.
//	original=load_image("peppers.raw",Ni,Nj);
//
//	PACKET=0;
//
//	printf("-------------------------------\n");
//	printf("\nTesting wavelet transform:\n");
//	printf("\t(Check source to see which banks are being utilized.)\n");
//	wavelet_transform_test(original,Ni,Nj);
//
//	printf("-------------------------------\n");
//	printf("\nTesting complex wavelet transform:\n");
//	complx_wavelet_transform_test(original,Ni,Nj);
//
//	PACKET=1;
//
//	printf("-------------------------------\n");
//	printf("\nTesting wavelet packet transform:\n");
//	wavelet_transform_test(original,Ni,Nj);
//
//	printf("-------------------------------\n");
//	printf("\nTesting complex wavelet packet transform:\n");
//	// Complex wavelet packets fully tested. Please check with literature.
//	complx_wavelet_transform_test(original,Ni,Nj);
//
//	return(1);
//}
