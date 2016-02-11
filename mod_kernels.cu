__global__ void cuda_set_sg(int *sxz, int sxbeg, int szbeg, int jsx, int jsz, int ns, int nz)
/*< set the positions of sources/geophones >*/
{
	int id=blockDim.x*blockIdx.x + threadIdx.x;
    	if (id<ns) sxz[id]=(szbeg+id*jsz)+nz*(sxbeg+id*jsx);
}

__global__ void cuda_ricker_wavelet(float *wlt, float amp, float fm, float dt, int nt)
/*< generate ricker wavelet with time deley >*/
{
	int it=blockDim.x*blockIdx.x + threadIdx.x;
    	if (it<nt)
	{
	    	float tmp = PI*fm*(it*dt-1.0/fm);
	    	tmp *=tmp;
		wlt[it]=amp*(1.0-2.0*tmp)*expf(-tmp);
	}
}

__global__ void cuda_add_source(float *p, float *source, int *sxz, int ns, bool add)
/*< add==true, add (inject) the source; add==false, subtract the source >*/
{
	int id=blockDim.x*blockIdx.x + threadIdx.x;
    	if (id<ns)
	{
		if (add)	p[sxz[id]]+=source[id];
		else 		p[sxz[id]]-=source[id];
	}	
}

__global__ void cuda_record(float*p, float *seis, int *gxz, int ng)
/*< record the seismogram at time it >*/
{
	int id=threadIdx.x+blockDim.x*blockIdx.x;
    	if (id<ng) seis[id]=p[gxz[id]];
}

__global__ void cuda_step_forward(float *p0, float *p1, float *vv, float dtz, float dtx, int nz, int nx)
/*< step forward: dtz=dt/dx; dtx=dt/dz; >*/
{
	int i1=threadIdx.x+blockIdx.x*blockDim.x;
	int i2=threadIdx.y+blockIdx.y*blockDim.y;
	int id=i1+i2*nz;

	__shared__ float s_p0[Block_Size2+2][Block_Size1+2];
	__shared__ float s_p1[Block_Size2+2][Block_Size1+2];
	if(threadIdx.x<1)
	{
		s_p0[threadIdx.y+1][threadIdx.x]=(blockIdx.x>0)?p0[id-1]:0.0;	
		s_p1[threadIdx.y+1][threadIdx.x]=(blockIdx.x>0)?p1[id-1]:0.0;
	}
	if(threadIdx.x>=blockDim.x-1)
	{
		s_p0[threadIdx.y+1][threadIdx.x+2]=(blockIdx.x<gridDim.x-1)?p0[id+1]:0.0;
		s_p1[threadIdx.y+1][threadIdx.x+2]=(blockIdx.x<gridDim.x-1)?p1[id+1]:0.0;
	}
	if(threadIdx.y<1)
	{
		s_p0[threadIdx.y][threadIdx.x+1]=(blockIdx.y>0)?p1[id-nz]:0.0;
	 	s_p1[threadIdx.y][threadIdx.x+1]=(blockIdx.y>0)?p1[id-nz]:0.0;
	}
	if(threadIdx.y>=blockDim.y-1)
	{
		s_p0[threadIdx.y+2][threadIdx.x+1]=(blockIdx.y<gridDim.y-1)?p1[id+nz]:0.0;
		s_p1[threadIdx.y+2][threadIdx.x+1]=(blockIdx.y<gridDim.y-1)?p1[id+nz]:0.0;
	}
	s_p0[threadIdx.y+1][threadIdx.x+1]=p0[id];
	s_p1[threadIdx.y+1][threadIdx.x+1]=p1[id];
	__syncthreads();

	if (i1<nz && i2<nx)
	{
		float v1=vv[id]*dtz;
		float v2=vv[id]*dtx; 
		float c1=v1*v1*(s_p1[threadIdx.y+1][threadIdx.x+2]-2.0*s_p1[threadIdx.y+1][threadIdx.x+1]+s_p1[threadIdx.y+1][threadIdx.x]);
		float c2=v2*v2*(s_p1[threadIdx.y+2][threadIdx.x+1]-2.0*s_p1[threadIdx.y+1][threadIdx.x+1]+s_p1[threadIdx.y][threadIdx.x+1]);
/*
		if(i1==0)// top boundary is free surface boundary condition, commentted!!
		{
			c1=v1*(-s_p1[threadIdx.y+1][threadIdx.x+1]+s_p1[threadIdx.y+1][threadIdx.x+2]
						+s_p0[threadIdx.y+1][threadIdx.x+1]-s_p0[threadIdx.y+1][threadIdx.x+2]);
			if(i2>0 && i2<nx-1) c2=0.5*c2;
		}
*/
		if(i1==nz-1) /* bottom boundary */
		{
			c1=v1*(s_p1[threadIdx.y+1][threadIdx.x]-s_p1[threadIdx.y+1][threadIdx.x+1]
						-s_p0[threadIdx.y+1][threadIdx.x]+s_p0[threadIdx.y+1][threadIdx.x+1]);
			if(i2>0 && i2<nx-1) c2=0.5*c2;
		}

		if(i2==0)/* left boundary */
		{
			if(i1>0 && i1<nz-1) c1=0.5*c1;
			c2=v2*(-s_p1[threadIdx.y+1][threadIdx.x+1]+s_p1[threadIdx.y+2][threadIdx.x+1]
						+s_p0[threadIdx.y+1][threadIdx.x+1]-s_p0[threadIdx.y+2][threadIdx.x+1]);

		}

		if(i2==nx-1) /* right boundary */
		{
			if(i1>0 && i1<nz-1) c1=0.5*c1;
			c2=v2*(s_p1[threadIdx.y][threadIdx.x+1]-s_p1[threadIdx.y+1][threadIdx.x+1]
						-s_p0[threadIdx.y][threadIdx.x+1]+s_p0[threadIdx.y+1][threadIdx.x+1]);
		}
		p0[id]=2.0*s_p1[threadIdx.y+1][threadIdx.x+1]-s_p0[threadIdx.y+1][threadIdx.x+1]+c1+c2;
	}
}
