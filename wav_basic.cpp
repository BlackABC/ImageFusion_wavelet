// Onur G. Guleryuz 1995, 1996, 1997,
// University of Illinois at Urbana-Champaign,
// Princeton University,
// Polytechnic University.

#include <stdio.h>
#include <stdlib.h>
#include "wav_filters_extern.h"
#include "wav_gen.h"
#include "alloc.h"

// Extend data by ofs pixels on both ends using periodic 
// symmetry. data should be "pointing in" to allow up to data[-ofs]
// and data[N-1+ofs]. 
// ofs>=N is handled.

void extend_periodic(float *data,int N,int ofs)

{
	int i,k;

	k=N-1;
	for(i=1;i<=ofs;i++) {

		data[-i]=data[N-i];
		data[k+i]=data[i-1];
	}
}

// Extend data by ofs pixels on both ends using mirror 
// symmetry. data should be "pointing in" to allow up to data[-ofs]
// and data[N-1+ofs]. 
//
// There are two possible "mirrors", i.e., data[-1]=data[1] (mirror on 0) 
// or data[-1]=data[0] (mirror "between" -1 and 0). phase is either 0, 1, or 2
// to select among these cases. See if(phase==... below.	
	
void extend_mirror(float *data,int N,int ofs,int phase)

{
	int i,k;
	int phase1,phase2;

	if((phase<0)||(phase>2)) {

		printf("extend_mirror: illegal phase\n");
		exit(1);
	}

	if(phase==2){
		// At left boundary the mirror is on 0.
		phase1=0;

		// At right boundary the mirror is on N-1.
		phase2=1;
	}
	else {

		// phase==0, at left mirror on 0, and at right it's between.
		// phase==1, at left mirror between, and at right it's on N-1.
		phase1=phase2=phase;
	}

	k=N-1;
		
	if(N==1) {

		for(i=1;i<=ofs;i++) {

			data[-i]=data[0];
			data[k+i]=data[0];
		}
	}
	else
		for(i=1;i<=ofs;i++) {

			data[-i]=data[i-phase1];
			data[k+i]=data[N-phase2-i];
		}
}

// This routine is used in evaluating forward wavelet transforms.
//
// data of dimension N contains the input.
// float *filter contains the Nf filter taps.
// ofs is the decimation factor.
//
// filtering starts at beg and continues until N+beg, yielding N/ofs samples. 
// data needs to be padded suitably at both ends.
// This hoopla about the starting point is in place to ensure 
// correct boundary processing when evaluating wavelet transforms.
//
// coef returns the result of the filtering and decimation.
// 
// Viewed as a transform the basis function that generates the
// coefficient of index 0 is "<- fliplr(filter) -> (beg)", 
// i.e., the scalar product obtained by positioning 
// the flipped filter (formed into a row vector) 
// with the right side positioned on top of sample beg in data.

int filter_n_decimate(float *data,float *coef,int N,float *filter,
						int Nf,int beg,int ofs)

{
	int i,j,k;
	float temp;

	k=0;
	for(i=beg;i<N+beg;i+=ofs) {

		temp=0;
		for(j=i;j>i-Nf;j--)
			temp+=data[j]*filter[i-j];

		coef[k]=temp;
		k++;
	}

	return(k);
}

// This routine is used in evaluating inverse wavelet transforms.
//
// Upsample and filter to implement a synthesis filter bank.
// Basically the inverse of filter_n_decimate() above.
// No actual upsampling etc., to avoid zero multiplies.
//
// coef is the input and contains coefficients resulting from 
// filtering and decimation by ofs. 
//
// Because this routine is typically called twice for inverse wavelet
// transforms, once for the low band and once for high, data is incremented via +=.
// Hence data should initially be set to 0.
//
// Viewed as a transform, filter_n_decimate above starts by putting the 
// "flipped" forward filter at begf, i.e.,  <- fliplr(ffilter) ->(begf)
// The resulting coefficient is at index 0.
// Here we are upsampling by ofs and filtering with the inv. filter
// time advanced by fbeg so the coefficient at index 0, multiplies the basis function
// (-fbeg)<-(ifilter)->(Nf-1-fbeg=beg).
// 
// How is that for confusing? All of this makes sense,
// when you write it down carefully. But who has the time ...:)

void upsample_n_filter(float *coef,float *data,int N,float *filter,
						int Nf,int beg,int ofs)

{
	int i,j,l,n,p,fbeg;
	float temp;

	fbeg=Nf-beg-1;

	for(i=0;i<N;i++) {

		l=0;
		n=(i+fbeg)%ofs;
		p=(i+fbeg)/ofs;

		temp=0;
		for(j=n;j<Nf;j+=ofs) {

			temp+=coef[p-l]*filter[j];
			l++;
		}

		// += for successive calls for low and high bands. 
		data[i]+=temp;
	}
}

// Used in forward wavelet trf.
//
// All rows of the image are run through the dual filterbank
// given by lp and hp. Results are returned in image "in place".
// low freq. coefficients start at 0 and end at N/2, where high freq.
// info starts.
void filt_n_dec_all_rows(float **image,int Ni,int Nj,float *lp,int Nl,
							float *hp,int Nh)

{
	int i,j,ext,ofs;
	float *data;

	ext=max(Nl,Nh);
	data=allocate_1d_float((Nj+2*ext),0);
	// Point in for required extensions.
	data+=ext;

	// offset for the location where high band coefficients
	// should be copied.
	ofs=Nj>>1;

	for(i=0;i<Ni;i++) {
		
		// Copy row.
		for(j=0;j<Nj;j++)
			data[j]=image[i][j];

		// Extend.
		if(PS)
			extend_periodic(data,Nj,ext);
		else {
			// Mirrors at end points.
			extend_mirror(data,Nj,ext,2);
		}

		filter_n_decimate(data,image[i],Nj,lp,Nl,begflp,2);
		filter_n_decimate(data,image[i]+ofs,Nj,hp,Nh,begfhp,2);
		}

	data-=ext;
	free((void *)data);
	}

// Used in forward wavelet trf.
//
// Index swapped version of filt_n_dec_all_rows.
void filt_n_dec_all_cols(float **image,int Ni,int Nj,float *lp,int Nl,
							float *hp,int Nh)

{
	int i,j,ext,ofs;
	float *data,*data2;

	ext=max(Nl,Nh);
	data=allocate_1d_float((Ni+2*ext),0);
	// Point in for required extensions.
	data+=ext;

	// Need second array since cannot access columns directly via pointers.
	data2=allocate_1d_float(Ni,0);

	// offset for the location where high band coefficients
	// should be copied.
	ofs=Ni>>1;

	for(j=0;j<Nj;j++) {
		
		// Copy column.
		for(i=0;i<Ni;i++)
			data[i]=image[i][j];

		// Extend.
		if(PS)
			extend_periodic(data,Ni,ext);
		else {
			// Mirrors at end points.
			extend_mirror(data,Ni,ext,2);
		}

		filter_n_decimate(data,data2,Ni,lp,Nl,begflp,2);
		filter_n_decimate(data,data2+ofs,Ni,hp,Nh,begfhp,2);

		for(i=0;i<Ni;i++)
			image[i][j]=data2[i];
		}

	data-=ext;
	free((void *)data);
	free((void *)data2);
	}

// Used in inverse wavelet trf.
//
// All rows of the image are run through the dual synthesis filterbank
// given by lp and hp having number of taps Nl and Nh respectively. 
// Results are returned in image "in place".
//
// Here we actaully need to know the applied shift to the transform,
// hence the passed parameters lev and shift_arr.

void ups_n_filt_all_rows(float **image,int Ni,int Nj,float *lp,int Nl,
							float *hp,int Nh,int lev,int *shift_arr)

{
	int i,j,k,ext,ofs;
	float *data1,*data2;

	ext=max(Nl,Nh);
	ofs=Nj>>1;

	data1=allocate_1d_float((ofs+2*ext),0);
	data2=allocate_1d_float((ofs+2*ext),0);
	data1+=ext;data2+=ext;

	for(i=0;i<Ni;i++) {	

		for(j=0;j<ofs;j++) {

			k=j+ofs;
			// low pass and high pass coefficients.
			data1[j]=image[i][j];image[i][j]=0;
			data2[j]=image[i][k];image[i][k]=0;
		}

		// Take care of the extension at the boundaries.
		if(PS) {

			extend_periodic(data1,ofs,ext);	
			extend_periodic(data2,ofs,ext);	
		}
		else {

			// shift_arr[lev]=0 OR 1.
			// The symmetric banks used are positioned so that on the left side
			// the forward lowpass filter has its point of symmetry exactly on top of 0
			// and forward high pass filter has its point of symmetry exactly on top of 1. 
			// This is coordinated by calling filter_n_decimate() with the proper value of beg.
			//
			// The above is reversed at the right side assuming decimation by 2. 
			// Hence shift_arr[lev]==0 => at left boundary the mirror is at 0, 
			// at right boundary it's between, see extend_mirror().

			// If shift_arr[lev]==1 then the positions for the low and high are reversed
			// since the filters are shifted so if
			// shift_arr[lev]==1 => at left boundary the mirror is between, 
			// at right boundary it's at 0.
			extend_mirror(data1,ofs,ext,shift_arr[lev]);	
			extend_mirror(data2,ofs,ext,1-shift_arr[lev]);	
		}

		// invert low pass.
		upsample_n_filter(data1,image[i],Nj,lp,Nl,begilp,2);
		// invert high pass.
		upsample_n_filter(data2,image[i],Nj,hp,Nh,begihp,2);
	}

	data1-=ext;data2-=ext;
	free((void *)data1);
	free((void *)data2);
}

// Used in inverse wavelet trf.
//
// Index swapped version of ups_n_filt_all_rows.
void ups_n_filt_all_cols(float **image,int Ni,int Nj,float *lp,int Nl,
							float *hp,int Nh,int lev,int *shift_arr)

{
	int i,j,k,ext,ofs;
	float *data1,*data2,*data3;

	ext=max(Nl,Nh);
	ofs=Ni>>1;

	data1=allocate_1d_float((ofs+2*ext),0);
	data2=allocate_1d_float((ofs+2*ext),0);
	data1+=ext;data2+=ext;

	// Need third array since cannot access columns directly via pointers.
	data3=allocate_1d_float(Ni,0);

	for(j=0;j<Nj;j++) {	

		for(i=0;i<ofs;i++) {

			k=i+ofs;
			// low pass and high pass coefficients.
			data1[i]=image[i][j];
			data2[i]=image[k][j];
			data3[i]=data3[k]=0;
		}

		if(PS) {

			extend_periodic(data1,ofs,ext);	
			extend_periodic(data2,ofs,ext);	
		}
		else {

			extend_mirror(data1,ofs,ext,shift_arr[lev]);	
			extend_mirror(data2,ofs,ext,1-shift_arr[lev]);	
		}

		upsample_n_filter(data1,data3,Ni,lp,Nl,begilp,2);
		upsample_n_filter(data2,data3,Ni,hp,Nh,begihp,2);

		for(i=0;i<Ni;i++)
			image[i][j]=data3[i];
	}

	data1-=ext;data2-=ext;
	free((void *)data1);
	free((void *)data2);
	free((void *)data3);
}

