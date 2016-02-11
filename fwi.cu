#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <cuda_runtime.h>

extern "C" {
#include <rsf.h>
}

#ifndef true
#define true    (1)
#endif
#ifndef false
#define false   (0)
#endif
#ifndef EPS
#define EPS	SF_EPS
#endif

#define PI 	SF_PI
#define Block_Size1 16	/* 1st dim block size */
#define Block_Size2 16	/* 2nd dim block size */
#define Block_Size  512	/* vector computation blocklength */

#include "mod_kernel.cu"

void sf_check_gpu_error (const char *msg) 
/*< check GPU errors >*/
{
    cudaError_t err = cudaGetLastError ();
    if (cudaSuccess != err) { 
	sf_error ("Cuda error: %s: %s", msg, cudaGetErrorString (err)); 
	exit(0);   
    }
}

void expand(float*vv, float *v0, int nz, int nx, int nz1, int nx1)
/*< round up the model size to be multiples of block size >*/
{
	int i1,i2,i11,i22;

	for(i2=0; i2<nx; i2++)
	for(i1=0; i1<nz; i1++)
	{
		i11=(i1<nz1)?i1:(nz1-1);
		i22=(i2<nx1)?i2:(nx1-1);
		vv[i1+i2*nz]=v0[i11+nz1*i22];
	}	
}

void window(float *v0,float *vv, int nz, int nx, int nz1, int nx1)
/*< window the portion to be the same size as initial model >*/
{
	int i1, i2;

	for(i2=0; i2<nx1; i2++)
	for(i1=0; i1<nz1; i1++)
		v0[i1+i2*nz1]=vv[i1+nz*i2];
}

void matrix_transpose(float *matrix, float *trans, int n1, int n2)
/*< matrix transpose: matrix tansposed to be trans >*/
{
	int i1, i2;

	for(i2=0; i2<n2; i2++)
	for(i1=0; i1<n1; i1++)
	    trans[i2+n2*i1]=matrix[i1+n1*i2];
}

int main(int argc, char *argv[])
{
	/* variables on host */
	bool verb, precon, csdgather;
	int is, it, iter, niter, distx, distz, csd, rbell;
	int nz, nx, nz1, nx1, nt, ns, ng;
	int sxbeg, szbeg, gxbeg, gzbeg, jsx, jsz, jgx, jgz;/*  parameters of acquisition geometery */
	float dx, dz, fm, dt, dtx, dtz, mstimer,amp, obj1, obj, beta, epsil, alpha;
	float *v0, *vv, *dobs, *trans, *objval, *ptr=NULL;
	sf_file vinit, shots, vupdates, grads, objs, illums;

	/* variables on device */
	int 	*d_sxz, *d_gxz;			
	float 	*d_wlt, *d_vv, *d_illum, *d_lap, *d_vtmp, *d_sp0, *d_sp1, *d_gp0, *d_gp1,*d_bndr;
	float	*d_dobs, *d_dcal, *d_derr, *d_g0, *d_g1, *d_cg, *d_pars, *d_alpha1, *d_alpha2;

    	/* initialize Madagascar */
    	sf_init(argc,argv);

    	/* set up I/O files */
    	vinit=sf_input ("in");   /* initial velocity model, unit=m/s */
	shots=sf_input("shots"); /* recorded shots from exact velocity model */
    	vupdates=sf_output("out"); /* updated velocity in iterations */ 
    	grads=sf_output("grads");  /* gradient in iterations */ 
	objs=sf_output("objs");/* values of objective function in iterations */
	illums=sf_output("illums");/* source illumination in iterations */

    	/* get parameters from velocity model and recorded shots */
	if (!sf_getbool("verb",&verb)) verb=true;/* vebosity */
    	if (!sf_histint(vinit,"n1",&nz1)) sf_error("no n1");/* n1 */
    	if (!sf_histint(vinit,"n2",&nx1)) sf_error("no n2");/* n2 */
    	if (!sf_histfloat(vinit,"d1",&dz)) sf_error("no d1");/* d1 */
   	if (!sf_histfloat(vinit,"d2",&dx)) sf_error("no d2");/* d2 */
	if (!sf_getbool("precon",&precon)) precon=false;/* precondition or not */
    	if (!sf_getint("niter",&niter))   niter=100;	/* number of iterations */
	if (!sf_getint("rbell",&rbell))	  rbell=2;	/* radius of bell smooth */

   	if (!sf_histint(shots,"n1",&nt)) sf_error("no nt");
	/* total modeling time steps */
   	if (!sf_histint(shots,"n2",&ng)) sf_error("no ng");
	/* total receivers in each shot */
   	if (!sf_histint(shots,"n3",&ns)) sf_error("no ns");
	/* number of shots */
   	if (!sf_histfloat(shots,"d1",&dt)) sf_error("no dt");
	/* time sampling interval */
   	if (!sf_histfloat(shots,"amp",&amp)) sf_error("no amp");
	/* maximum amplitude of ricker */
   	if (!sf_histfloat(shots,"fm",&fm)) sf_error("no fm");
	/* dominant freq of ricker */
   	if (!sf_histint(shots,"sxbeg",&sxbeg)) sf_error("no sxbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"szbeg",&szbeg)) sf_error("no szbeg");
	/* x-begining index of sources, starting from 0 */
   	if (!sf_histint(shots,"gxbeg",&gxbeg)) sf_error("no gxbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"gzbeg",&gzbeg)) sf_error("no gzbeg");
	/* x-begining index of receivers, starting from 0 */
   	if (!sf_histint(shots,"jsx",&jsx)) sf_error("no jsx");
	/* source x-axis  jump interval  */
   	if (!sf_histint(shots,"jsz",&jsz)) sf_error("no jsz");
	/* source z-axis jump interval  */
   	if (!sf_histint(shots,"jgx",&jgx)) sf_error("no jgx");
	/* receiver x-axis jump interval  */
   	if (!sf_histint(shots,"jgz",&jgz)) sf_error("no jgz");
	/* receiver z-axis jump interval  */
   	if (!sf_histint(shots,"csdgather",&csd)) sf_error("csdgather or not required");
	/* default, common shot-gather; if n, record at every point*/

	/* put the labels, legends and parameters for output rsf file */
	sf_putint(vupdates,"n1",nz1);	
	sf_putint(vupdates,"n2",nx1);
	sf_putfloat(vupdates,"d1",dz);
	sf_putfloat(vupdates,"d2",dx);
	sf_putstring(vupdates,"label1","Depth");
	sf_putstring(vupdates,"label2","Distance");
	sf_putstring(vupdates,"label3","Iteration");
	sf_putint(vupdates,"n3",niter);
	sf_putint(vupdates,"d3",1);
	sf_putint(vupdates,"o3",1);
	sf_putint(grads,"n1",nz1);	
	sf_putint(grads,"n2",nx1);
	sf_putint(grads,"n3",niter);
	sf_putfloat(grads,"d1",dz);
	sf_putfloat(grads,"d2",dx);
	sf_putint(grads,"d3",1);
	sf_putint(grads,"o3",1);
	sf_putstring(grads,"label1","Depth");
	sf_putstring(grads,"label2","Distance");
	sf_putstring(grads,"label3","Iteration");
	sf_putint(illums,"n1",nz1);	
	sf_putint(illums,"n2",nx1);
	sf_putfloat(illums,"d1",dz);
	sf_putfloat(illums,"d2",dx);
	sf_putint(illums,"n3",niter);
	sf_putint(illums,"d3",1);
	sf_putint(illums,"o3",1);
	sf_putint(objs,"n1",niter);
	sf_putint(objs,"n2",1);
	sf_putint(objs,"d1",1);
	sf_putint(objs,"o1",1);

	dtx=dt/dx; 
	dtz=dt/dz; 
	csdgather=(csd>0)?true:false;
	/* round the size up to multiples of Block size */
	nx=(int)((nx1+Block_Size1-1)/Block_Size1)*Block_Size1;
	nz=(int)((nz1+Block_Size2-1)/Block_Size2)*Block_Size2; 
	dim3 dimg=dim3(nz/Block_Size1, nx/Block_Size2), dimb=dim3(Block_Size1, Block_Size2); 

	v0=(float*)malloc(nz1*nx1*sizeof(float));/* initial velocity model */
	vv=(float*)malloc(nz*nx*sizeof(float));	 /* extended velocity model, size=multiple of 16x16 block */
	dobs=(float*)malloc(ng*nt*sizeof(float));/* observations, one shot */
	trans=(float*)malloc(ng*nt*sizeof(float));/* transposed one shot */
	objval=(float*)malloc(niter*sizeof(float));/* objective/misfit function */
	sf_floatread(v0, nz1*nx1, vinit);	/* read the initial velcity model, size=nz1*nx1 */
	expand(vv, v0, nz, nx, nz1, nx1);	/* expand the model to be of size nz*nx */
	memset(dobs, 0, ng*nt*sizeof(float));	
	memset(objval, 0, niter*sizeof(float));

    	cudaSetDevice(0);
	sf_check_gpu_error("Failed to initialize device!");
	/* allocate memory for device variables */
	cudaMalloc(&d_vv, nz*nx*sizeof(float));	/* velocity */
	cudaMalloc(&d_sp0, nz*nx*sizeof(float));/* source wavefield p0 */
	cudaMalloc(&d_sp1, nz*nx*sizeof(float));/* source wavefield p1 */
	cudaMalloc(&d_gp0, nz*nx*sizeof(float));/* geophone/receiver wavefield p0 */
	cudaMalloc(&d_gp1, nz*nx*sizeof(float));/* geophone/receiver wavefield p1 */
	cudaMalloc(&d_wlt, nt*sizeof(float));	/* ricker wavelet */
	cudaMalloc(&d_sxz, ns*sizeof(float));	/* source positions */
	cudaMalloc(&d_gxz, ng*sizeof(float));	/* geophone positions */
	cudaMalloc(&d_bndr, nt*(2*nz+nx)*sizeof(float));/* boundaries for wavefield reconstruction */
	cudaMalloc(&d_dobs, ng*nt*sizeof(float));/* observed seismic data */
	cudaMalloc(&d_dcal, ng*sizeof(float));	/* calculated/synthetic seismic data */
	cudaMalloc(&d_derr, ns*ng*nt*sizeof(float));/* residual/error between synthetic and observation */
	cudaMalloc(&d_g0, nz*nx*sizeof(float));	/* gradient at previous step */
	cudaMalloc(&d_g1, nz*nx*sizeof(float));	/* gradient at curret step */
	cudaMalloc(&d_cg, nz*nx*sizeof(float));	/* conjugate gradient */
	cudaMalloc(&d_lap, nz*nx*sizeof(float));/* laplace of the source wavefield */
	cudaMalloc(&d_illum, nz*nx*sizeof(float));/* illumination of the source wavefield */
	cudaMalloc(&d_pars, 4*sizeof(float));	/* d_pars[0]: obj; d_pars[1]: beta; d_pars[2]: epsilon; d_pars[3]: alpha; */
	cudaMalloc(&d_alpha1, ng*sizeof(float));/* d_alpha1[]: numerator of alpha, length=ng */
	cudaMalloc(&d_alpha2, ng*sizeof(float));/* d_alpha2[]: denominator of alpha, length=ng	*/
	cudaMalloc(&d_vtmp, nz*nx*sizeof(float));/* temporary velocity computed with epsil */
	sf_check_gpu_error("Failed to allocate required memory!");

	/* initialize varibles */
	cudaMemcpy(d_vv, vv, nz*nx*sizeof(float), cudaMemcpyHostToDevice);
	cudaMemset(d_sp0, 0, nz*nx*sizeof(float));
	cudaMemset(d_sp1, 0, nz*nx*sizeof(float));
	cudaMemset(d_gp0, 0, nz*nx*sizeof(float));
	cudaMemset(d_gp1, 0, nz*nx*sizeof(float));
	cuda_ricker_wavelet<<<(nt+511)/512,512>>>(d_wlt, amp, fm, dt, nt);
	/* configure the source/geophone geometry */
	if (!(sxbeg>=0 && szbeg>=0 && sxbeg+(ns-1)*jsx<nx1 && szbeg+(ns-1)*jsz<nz1))	
	{ sf_error("sources exceeds the computing zone!\n"); exit(1);}
	cuda_set_sg<<<(ns+511)/512,512>>>(d_sxz, sxbeg, szbeg, jsx, jsz, ns, nz);
	distx=sxbeg-gxbeg;
	distz=szbeg-gzbeg;
	if (!(gxbeg>=0 && gzbeg>=0 && gxbeg+(ng-1)*jgx<nx1 && gzbeg+(ng-1)*jgz<nz1))	
	{ sf_error("geophones exceeds the computing zone!\n"); exit(1);}
	if (csdgather)	{
		if (!(	(sxbeg+(ns-1)*jsx)+(ng-1)*jgx-distx <nx1  && (szbeg+(ns-1)*jsz)+(ng-1)*jgz-distz <nz1))	
		{ sf_error("geophones exceeds the computing zone!\n"); exit(1);}
	}
	cuda_set_sg<<<(ng+511)/512,512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
	cudaMemset(d_bndr, 0, nt*(2*nz+nx)*sizeof(float));
	cudaMemset(d_dobs, 0, ng*nt*sizeof(float));
	cudaMemset(d_dcal, 0, ng*sizeof(float));
	cudaMemset(d_derr, 0, ns*ng*nt*sizeof(float));
	cudaMemset(d_g0, 0, nz*nx*sizeof(float));
	cudaMemset(d_g1, 0, nz*nx*sizeof(float));
	cudaMemset(d_cg, 0, nz*nx*sizeof(float));
	cudaMemset(d_lap, 0, nz*nx*sizeof(float));
	cudaMemset(d_illum, 0, nz*nx*sizeof(float));
	cudaMemset(d_pars, 0, 4*sizeof(float));
	cudaMemset(d_alpha1, 0, ng*sizeof(float));
	cudaMemset(d_alpha2, 0, ng*sizeof(float));
	cudaMemset(d_vtmp, 0, nz*nx*sizeof(float));

	/* creat timing variables on device */
	cudaEvent_t start, stop;
  	cudaEventCreate(&start);	
	cudaEventCreate(&stop);

	/* Our FWI is carried out via the following steps: 
		step 1. do modeling to calculate synthetic seismograms;
		step 2. backpropagate wavefield to obtain the gradient;
		step 3. compute the conjugate gradient;
		step 4. find a good step length;
		step 5. update the current velocity model with estimated step length and conjugate gradient;
		step 6. repeat step 1-5 in every iteration;
	*/
	for(iter=0; iter<niter; iter++)
	{
		cudaEventRecord(start);/* record starting time */

		sf_seek(shots, 0L, SEEK_SET);
		cudaMemcpy(d_g0, d_g1, nz*nx*sizeof(float), cudaMemcpyDeviceToDevice);
		cudaMemset(d_g1, 0, nz*nx*sizeof(float));
		cudaMemset(d_illum, 0, nz*nx*sizeof(float));
		for(is=0;is<ns;is++)
		{
			/* read one observed shot record */
			sf_floatread(dobs, ng*nt, shots);
			matrix_transpose(dobs, trans, nt, ng);
			cudaMemcpy(d_dobs, trans, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			/* configure the source/geophone geometry */
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
			}
			cudaMemset(d_sp0, 0, nz*nx*sizeof(float));
			cudaMemset(d_sp1, 0, nz*nx*sizeof(float));
			/* advance time steps for source wavefield */
			for(it=0; it<nt; it++)
			{
				/* forward modeling: source */
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, true);
				cuda_step_forward<<<dimg,dimb>>>(d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;
				/* record the calculated synthetic seismograms */
				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				cuda_cal_residuals<<<(ng+511)/512, 512>>>(d_dcal, &d_dobs[it*ng], &d_derr[is*ng*nt+it*ng], ng);
				/* saving the boundaries */
				cuda_rw_bndr<<<(2*nz+nx+511)/512,512>>>(&d_bndr[it*(2*nz+nx)], d_sp0, nz, nx, true);
			}

			cudaMemset(d_gp0, 0, nz*nx*sizeof(float));
			cudaMemset(d_gp1, 0, nz*nx*sizeof(float));
			/* backward propagation */
			for(it=nt-1; it>-1; it--)
			{
				/* source backpropagation */
				ptr=d_sp0;d_sp0=d_sp1;d_sp1=ptr;
				cuda_rw_bndr<<<(2*nz+nx+255)/256,256>>>(&d_bndr[it*(2*nz+nx)], d_sp1, nz, nx, false);
				cuda_step_backward<<<dimg,dimb>>>(d_illum, d_lap, d_sp0, d_sp1, d_vv, dtz, dtx, nz, nx);
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, false);

				/* extrapolate residual wavefield */
				cuda_add_source<<<(ng+511)/512, 512>>>(d_gp1, &d_derr[is*ng*nt+it*ng], d_gxz, ng, true);
				cuda_step_forward<<<dimg,dimb>>>(d_gp0, d_gp1, d_vv, dtz, dtx, nz, nx);

				/* do summation for the gradient: equation 9 */
				cuda_cal_gradient<<<dimg,dimb>>>(d_g1, d_lap, d_gp1, nz, nx);
				ptr=d_gp0; d_gp0=d_gp1; d_gp1=ptr;
			}
		}
		/* calculate the value of the objective function */
		cuda_cal_objective<<<1, Block_Size>>>(&d_pars[0], d_derr, ns*ng*nt);
		cudaMemcpy(&obj, &d_pars[0], sizeof(float), cudaMemcpyDeviceToHost);

		/* compute source illumination and save it into snaps */
		cudaMemcpy(vv, d_illum, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, vv, nz, nx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, illums);

		/* compute the gradient of FWI by scaling, precondition incorporated here: equations 9 and 10 */
		cuda_scale_gradient<<<dimg,dimb>>>(d_g1, d_vv, d_illum, nz, nx, precon);
		cudaMemcpy(vv, d_g1, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, vv, nz, nx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, grads);
		/* Gaussian smoothing for the sharp gradient */
		cuda_bell_smoothz<<<dimg,dimb>>>(d_g1, d_illum, rbell, nz, nx);
		cuda_bell_smoothx<<<dimg,dimb>>>(d_illum, d_g1, rbell, nz, nx);

		/* calculate the factor beta in conjugate gradient method: equation 7 */
		if (iter>0) cuda_cal_beta<<<1, Block_Size>>>(&d_pars[1], d_g0, d_g1, d_cg, nz*nx); 
		cudaMemcpy(&beta, &d_pars[1], sizeof(float), cudaMemcpyDeviceToHost);
		/* compute the conjugate gradient */
		cuda_cal_conjgrad<<<dimg, dimb>>>(d_g1, d_cg, beta, nz, nx);
		/* estimate epsilon according to equation 11 */
		cuda_cal_epsilon<<<1, Block_Size>>>(d_vv, d_cg, &d_pars[2], nz*nx);
		cudaMemcpy(&epsil, &d_pars[2], sizeof(float), cudaMemcpyDeviceToHost);

		sf_seek(shots, 0L, SEEK_SET);
		cudaMemset(d_alpha1, 0, ng*sizeof(float));
		cudaMemset(d_alpha2, 0, ng*sizeof(float));
		/* obtain a tentative velocity model to estimate a good stepsize alpha */
		cuda_cal_vtmp<<<dimg, dimb>>>(d_vtmp, d_vv, d_cg, epsil, nz, nx);
		for(is=0;is<ns;is++)
		{
			sf_floatread(dobs, ng*nt, shots);
			matrix_transpose(dobs, trans, nt, ng);
			cudaMemcpy(d_dobs, trans, ng*nt*sizeof(float), cudaMemcpyHostToDevice);
			/* configure the source/geophone geometry */
			if (csdgather)	{
				gxbeg=sxbeg+is*jsx-distx;
				cuda_set_sg<<<(ng+511)/512, 512>>>(d_gxz, gxbeg, gzbeg, jgx, jgz, ng, nz);
			}
			cudaMemset(d_sp0, 0, nz*nx*sizeof(float));
			cudaMemset(d_sp1, 0, nz*nx*sizeof(float));
			/* remodeling with tentative velocity model vtmp */
			for(it=0; it<nt; it++)
			{
				cuda_add_source<<<1,1>>>(d_sp1, &d_wlt[it], &d_sxz[is], 1, true);
				cuda_step_forward<<<dimg,dimb>>>(d_sp0, d_sp1, d_vtmp, dtz, dtx, nz, nx);
				ptr=d_sp0; d_sp0=d_sp1; d_sp1=ptr;

				cuda_record<<<(ng+511)/512, 512>>>(d_sp0, d_dcal, d_gxz, ng);
				/* compute the numerator and the denominator of alpha: equations 5 and 12 */
				cuda_sum_alpha12<<<(ng+511)/512, 512>>>(d_alpha1, d_alpha2, d_dcal, &d_dobs[it*ng], &d_derr[is*ng*nt+it*ng], ng);
			}
		}
		/* find a good stepsize alpha: equation 5*/
		cuda_cal_alpha<<<1,Block_Size>>>(&d_pars[3], d_alpha1, d_alpha2, epsil, ng);
		cudaMemcpy(&alpha, &d_pars[3], sizeof(float), cudaMemcpyDeviceToHost);

		/* update the velocity model according to previous velocity, conjugate gradient and estimated stepsize */
		cuda_update_vel<<<dimg,dimb>>>(d_vv, d_cg, alpha, nz, nx);
		cudaMemcpy(vv, d_vv, nz*nx*sizeof(float), cudaMemcpyDeviceToHost);
		window(v0, vv, nz, nx, nz1, nx1);
		sf_floatwrite(v0, nz1*nx1, vupdates);

		cudaEventRecord(stop);/* record ending time */
  		cudaEventSynchronize(stop);
  		cudaEventElapsedTime(&mstimer, start, stop);

		/* compute the normalized objective function */
		if(iter==0) 	{obj1=obj; objval[iter]=1.0;}
		else		objval[iter]=obj/obj1;

		if(verb) {/* output important information at each FWI iteration */
			sf_warning("obj=%f  beta=%f  epsil=%f  alpha=%f",obj, beta, epsil, alpha);
			sf_warning("iteration %d finished: %f (s)",iter+1, mstimer*1e-3);
		}
	}
	/* destroy timing varibles */
	cudaEventDestroy(start);
	cudaEventDestroy(stop);

	sf_floatwrite(objval,iter,objs);
	sf_fileclose(shots); 

	/* free varibles on device */
	cudaFree(d_vv);
	cudaFree(d_sp0);
	cudaFree(d_sp1);
	cudaFree(d_gp0);
	cudaFree(d_gp1);
	cudaFree(d_wlt);
	cudaFree(d_sxz);
	cudaFree(d_gxz);
	cudaFree(d_bndr);
	cudaFree(d_dobs);
	cudaFree(d_dcal);
	cudaFree(d_derr);
	cudaFree(d_g0);
	cudaFree(d_g1);
	cudaFree(d_cg);
	cudaFree(d_lap);
	cudaFree(d_illum);
	cudaFree(d_pars);
	cudaFree(d_alpha1);
	cudaFree(d_alpha2);
	cudaFree(d_vtmp);
	sf_check_gpu_error("Failed to free the allocated memory!");
	/* free varibles on host */
	free(v0);
	free(vv);
	free(dobs);
	free(trans);
	free(objval);

	exit(0);
}
