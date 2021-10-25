#include<thrust/device_vector.h>
#include<thrust/host_vector.h>
#include"base.h"
#include<stdio.h>
#include<assert.h>
void fill_n_v(device_vector_double& data , int n, typename device_vector_double::value_type v){
	 data.resize(n,v);
}
void fill_n_v(device_vector_int& data , int n, typename device_vector_int::value_type v){
	 data.resize(n,v);
}
__global__ void NotGetDoubleDataByIndexKernel(double* src, int* index, double* dst, int maxId)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= maxId) return;
	int theIndex = index[id];
	dst[theIndex] = src[id];
}
__global__ void GetDoubleDataByIndexKernel(double* src, int* index, double* dst, int maxId)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= maxId) return;
	int theIndex = index[id];
	src[id] = dst[theIndex];
}

void do_MapDoubleDataByIndex(double* src, int* index, double* dst, int number, bool get)
{
	int blockDim = 128;
	int gridDim = number/blockDim;
	gridDim = gridDim * blockDim < number ? gridDim + 1 : gridDim;
	if(get)
		GetDoubleDataByIndexKernel<<< gridDim, blockDim >>>(src, index, dst, number);
	else
		NotGetDoubleDataByIndexKernel<<< gridDim, blockDim >>>(src, index, dst, number);
	cudaDeviceSynchronize();
}

__global__ void VoltKernel(double* new_volt, double* volt, double* current, int* xpos, int* xneg,
						   int* ypos, int* yneg, int* zpos, int* zneg,
						   double DD, double dt, int maxId)
{
	int id = threadIdx.x + blockIdx.x * blockDim.x;
	if(id >= maxId) return;
	double base = volt[id];
	double baseX2 = 2 * base;
	double vxpos = volt[ xpos[id] ];
	double vxneg = volt[ xneg[id] ];
	double xdiff = vxpos + vxneg - baseX2;
	double vypos = volt[ ypos[id] ];
	double vyneg = volt[ yneg[id] ];
	double ydiff = vypos + vyneg - baseX2;
	double vzpos = volt[ zpos[id] ];
	double vzneg = volt[ zneg[id] ];
	double zdiff = vzpos + vzneg - baseX2;
	double dv = ( DD * (xdiff + ydiff + zdiff) - current[id] ) * dt;
	/*
	if(id == 1282860) {
		printf("%lf %lf %lf %lf %lf %lf\n", vxpos, vxneg, vypos, vyneg, vzpos, vzneg);
		printf("%lf, %lf, %lf\n", xdiff, ydiff, zdiff);
		printf("%lf + %lf =  %lf at VoltKernel\n", base, dv, base + dv);
	}
	*/
	new_volt[ id ] = base + dv;
}

void do_CalcVolt(double* new_volt, double* volt, double* current, int* xpos, int* xneg,
				 int* ypos, int* yneg, int* zpos, int* zneg,
				 double DD, double dt, int number)
{
	int blockDim = 128;
	int gridDim = number/blockDim;
	
	gridDim = gridDim * blockDim < number ? gridDim + 1 : gridDim;
	assert(gridDim * blockDim >= number);
	cudaDeviceSynchronize();
	
	VoltKernel<<< gridDim, blockDim >>>(new_volt, volt, current, xpos, xneg,
										ypos, yneg, zpos, zneg, DD, dt, number);
	cudaDeviceSynchronize();
}
