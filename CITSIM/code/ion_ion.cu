#include <stdio.h>

__global__ void CoulombForce2d(int N, float *DXd, float *DYd, float *Fd, float *Gd) {
  int idx = blockIdx.x * blockDim.x + threadIdx.x;
  float dx,dy,r,r12,r32,rm32;
  dx = DXd[idx]; dy = DYd[idx];
  r = dx*dx + dy*dy; r12 = sqrt(r); r32 = r*r12; rm32 = 1./r32;
  Fd[idx] = dx*rm32; Gd[idx] = dy*rm32;
}

extern "C" int ion_ion2d(int N, float *X, float *Y, float *F, float *G) {
  float *Xd,*Yd,*Fd,*Gd;
  dim3 dG, dB;
  long long int nBytes;
  int nMultiPr = 16;
  int nThreads;
  /* Allocate device memory */
  nBytes = N*sizeof(float);  
  cudaMalloc( (void **) &Xd, nBytes ); cudaMalloc( (void **) &Yd, nBytes );
  cudaMalloc( (void **) &Fd, nBytes ); cudaMalloc( (void **) &Gd, nBytes );
  /* Move ion position data to device */
  cudaMemcpy( Xd, X, nBytes, cudaMemcpyHostToDevice );
  cudaMemcpy( Yd, Y, nBytes, cudaMemcpyHostToDevice );  
  /* Device now has current ion positions */  
  /* Distribute computation of Coulomb force among device multiprocessors, thread processors
     dG dimension and size of grid   expressed in blocks
     dB dimension and size of blocks expressed in threads
     
     Create as many blocks as there are multiprocessors
     Distribute N ions equally among blocks
  */
  nThreads = N / nMultiPr + 1;
  dG.x = nMultiPr; dG.y = 1;
  dB.x = nThreads; dB.y = 1; dB.z = 1;
  /* Integrate SDE forward in time */
  CoulombForce2d<<<dG,dB>>>(N,Xd,Yd,Fd,Gd);
  /* Move computed electrostatic force to host memory */  
  cudaMemcpy( F, Fd, nBytes, cudaMemcpyDeviceToHost );
  cudaMemcpy( G, Gd, nBytes, cudaMemcpyDeviceToHost );
  return 0;
}
//-----------------------------------------------------------------------------------------
//-------------3D _ N-Body ion-ion interaction _ Prins algorithm-------------------------------
//-----------------------------------------------------------------------------------------

__device__ float3 bodyBodyInteraction(float4 bi, float4 bj, float3 ai)  
{  
  float3 r;  
  // r_ij  [3 FLOPS]  
  r.x = bj.x - bi.x;  
  r.y = bj.y - bi.y;  
  r.z = bj.z - bi.z;  
  // distSqr = dot(r_ij, r_ij) + EPS^2  [6 FLOPS]  
   float distSqr = r.x * r.x + r.y * r.y + r.z * r.z + 1e-6;  
  // invDistCube =1/distSqr^(3/2)  [4 FLOPS (2 mul, 1 sqrt, 1 inv)]  
   float distSixth = distSqr * distSqr * distSqr;  
   float invDistCube = 1.0f/sqrtf(distSixth);  
  // s = m_j * invDistCube [1 FLOP]  
   float s = bi.w * bj.w * invDistCube;  // *bi.w (bruno)
  // a_i =  a_i + s * r_ij [6 FLOPS]  
  ai.x += r.x * s;  
  ai.y += r.y * s;  
  ai.z += r.z * s;  
  return ai;  
}  
    
__device__ float3 tile_calculation(float4 myPosition, float3 accel)  
{  
  int i;
  extern __shared__ float4 shPosition[];
  for (i = 0; i < blockDim.x; i++) {
    accel = bodyBodyInteraction(myPosition, shPosition[i], accel);  
  }  
  return accel;  
} 
    
__global__ void calculate_forces(int N, void *devX,void *devA) 
{ 
  extern __shared__ float4 shPosition[];  
  float4 *globalX = (float4 *)devX;  
  float4 *globalA = (float4 *)devA;  
  float4 myPosition;  
  int i, tile;  
  float3 acc = {0.0f, 0.0f, 0.0f};  
  int gtid = blockIdx.x * blockDim.x + threadIdx.x; 
  //printf("Hello from block %d, thread %d\n", blockIdx.x, threadIdx.x);
  myPosition = globalX[gtid];
  /*printf("myPosition (x,y,z,w)=%f,%f,%f,%f\n",myPosition.x,
  myPosition.y,myPosition.z,myPosition.w); */
  for (i = 0, tile = 0; i < N; i += 192, tile++) {  
    int idx = tile * blockDim.x + threadIdx.x;  
    shPosition[threadIdx.x] = globalX[idx];  
    __syncthreads();  
    acc = tile_calculation(myPosition, acc);  
    __syncthreads();  
}  
      // Save the result in global memory for the integration step.  
      float4 acc4 = {acc.x, acc.y, acc.z, 0.0f};  
      globalA[gtid] = acc4;
      /*printf("acc.x, acc.y, acc.z=%f \t %f \t %f \n",acc.x, 
        acc.y, acc.z);*/
    }  

extern "C" int ion_ion3d(int N, float *X, float *F) {
  float4 *hostX = (float4 *) X;
  float4 *hostA = (float4 *) F;
  float4 *devX,*devA; 
  long long int nBytes;
  /*Define GPU grid*/
  dim3 DimGrid((N-1)/192+1,1,1);
  dim3 DimBlock(192,1,1);
  //printf("Hello from block ionion3d");
  cudaSetDevice(1);
  /* Allocate device memory */
  nBytes = N*sizeof(float4);
  cudaMalloc( (void **) &devX, nBytes );
  cudaMalloc( (void **) &devA, nBytes );
  /* Move ion position and ion charge data to device */
  cudaMemcpy( devX, hostX, nBytes, cudaMemcpyHostToDevice );
  cudaMemcpy( devA, hostA, nBytes, cudaMemcpyHostToDevice );  
  /* call Prins algorithm*/
  calculate_forces<<<DimGrid,DimBlock,192*sizeof(float4)>>>(N,devX,devA);
  /* Move computed electrostatic force to host memory */  
  cudaMemcpy( hostA, devA, nBytes, cudaMemcpyDeviceToHost );
  cudaFree(devX);
  cudaFree(devA);
  return 0;
}

//~ extern "C" int ion_ion3d(int N, float *X, float *Y, float *Z, float *W, float *Fx, float *Fy, float *Fz) {
  //~ float *Xd,*Yd,*Zd,*Wd,*Fxd,*Fyd,*Fzd;
  //~ long long int nBytes;
  //~ /*Define GPU grid*/
  //~ dim3 DimGrid((N-1)/192+1,1,1);
  //~ dim3 DimBlock(192,1,1);
  //~ /* Allocate device memory */
  //~ nBytes = N*sizeof(float);  
  //~ cudaMalloc( (void **) &Xd, nBytes ); cudaMalloc( (void **) &Yd, nBytes );
  //~ cudaMalloc( (void **) &Zd, nBytes ); cudaMalloc( (void **) &Wd, nBytes );
  //~ cudaMalloc( (void **) &Fxd, nBytes ); cudaMalloc( (void **) &Fyd, nBytes );
  //~ cudaMalloc( (void **) &Fzd, nBytes );
  //~ /* Move ion position and ion charge data to device */
  //~ cudaMemcpy( Xd, X, nBytes, cudaMemcpyHostToDevice );
  //~ cudaMemcpy( Yd, Y, nBytes, cudaMemcpyHostToDevice );  
  //~ cudaMemcpy( Zd, Z, nBytes, cudaMemcpyHostToDevice );
  //~ cudaMemcpy( Wd, W, nBytes, cudaMemcpyHostToDevice );
  //~ /* Integrate SDE forward in time */
  //~ calculate_forces<<<DimGrid,DimBlock>>>(N,Xd,Yd,Zd,Wd,Fxd,Fyd,Fzd);
  //~ /* Move computed electrostatic force to host memory */  
  //~ cudaMemcpy( Fx, Fxd, nBytes, cudaMemcpyDeviceToHost );
  //~ cudaMemcpy( Fy, Fyd, nBytes, cudaMemcpyDeviceToHost );
  //~ cudaMemcpy( Fz, Fzd, nBytes, cudaMemcpyDeviceToHost );
  //~ return 0;
//~ }
    
