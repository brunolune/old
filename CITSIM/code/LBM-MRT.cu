#include <stdio.h>

#ifndef LBM_H_INCLUDED
#include "LBM-MRT.h"
#endif

//
// velD3Q19[19]={OOO,OOP,OOM,OPO,OPP,OPM,OMO,OMP,OMM,POO,POP,POM,PPO,PMO,MOO,MOP,MOM,MPO,MMO};
//                 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18

__device__ void ComputeEqQ_D3Q19(float *f, float *qeq) {
   qeq[0]=0.; for (int n=0; n<19; n++) qeq[0] += f[n];
   qeq[1] = f[1]-f[2]+f[4]-f[5]+f[7]-f[8]+f[10]-f[11]+f[15]-f[16];
   qeq[2] = f[3]+f[4]+f[5]-f[6]-f[7]-f[8]+f[12]-f[13]+f[17]-f[18];
   qeq[3] = f[9]+f[10]+f[11]+f[12]+f[13]-f[14]-f[15]-f[16]-f[17]-f[18];
}

__device__ void ComputeEqF_D3Q19(float *qeq, float *feq) {
   float p,u,v,w,U32,eu;
   p=qeq[0]; u=qeq[1]; v=qeq[2]; w=qeq[3];
   U32=1.5*(u*u+v*v+w*w);
   eu= 0.;  feq[ 0]=0.333333333*(p - U32);
   eu= u;   feq[ 1]=0.055555556*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-u;   feq[ 2]=0.055555556*(p - U32 + (3. + 4.5*eu)*eu);
   eu= v;   feq[ 3]=0.055555556*(p - U32 + (3. + 4.5*eu)*eu);
   eu= v+u; feq[ 4]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu= v-u; feq[ 5]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-v;   feq[ 6]=0.055555556*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-v+u; feq[ 7]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-v-u; feq[ 8]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu= w;   feq[ 9]=0.055555556*(p - U32 + (3. + 4.5*eu)*eu);
   eu= w+u; feq[10]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu= w-u; feq[11]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu= w+v; feq[12]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu= w-v; feq[13]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-w;   feq[14]=0.055555556*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-w+u; feq[15]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-w-u; feq[16]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-w+v; feq[17]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   eu=-w-v; feq[18]=0.027777778*(p - U32 + (3. + 4.5*eu)*eu);
   for (int n=0; n<19; n++) if (feq[n]<0.) feq[n]=0.;
}

__device__ void ComputeEqM_D3Q19(float *qeq, float *meq) {
   float rho,jx,jy,jz,jx2,jy2,jz2,we,wej,wxx,twothirds,jdotj;
   rho=qeq[0]; jx=qeq[1]; jy=qeq[2]; jz=qeq[3];
   jx2=jx*jx; jy2=jy*jy; jz2=jz*jz;
   jdotj=jx2+jy2+jz2;
   we=0.0; wej=-475.0/63.0; wxx=0.0;
   //we=3.0; wej=-5.5; wxx=-0.5;
   twothirds=2.0/3.0;
   // Compute the equilibrium moments needed in collision operator for MRT LBM
   meq[0] = rho;                  //rho  - density
   meq[1] = -11.0*rho+19.0*jdotj; //e    - part of kinetic energy independent of density
   meq[2] = we*rho+wej*jdotj;     //eps  - part of kinetic energy square independent of density
   meq[3] = jx;                   //jx   - momentum in x direction
   meq[4] = -twothirds*jx;        //qx   - energy flux independent of mass flux
   meq[5] = jy;                   //jy   - momentum in y direction
   meq[6] = -twothirds*jy;        //qy   - energy flux independent of mass flux
   meq[7] = jz;                   //jz   - momentum in z direction
   meq[8] = -twothirds*jz;        //qz   - energy flux independent of mass flux
   meq[9] = 2.0*jx2-jy2-jz2;      //3pxx - symmetric traceless viscous stress tensor
   meq[10]= wxx*meq[9];           //3Pixx- quartic order moment
   meq[11]= jy2-jz2;              //pww  - symmetric traceless viscous stress tensor
   meq[12]= wxx*meq[11];          //Piww - quartic order moment
   meq[13]= jx*jy;                //pxy  - symmetric traceless viscous stress tensor
   meq[14]= jy*jz;                //pyz  - symmetric traceless viscous stress tensor
   meq[15]= jz*jx;                //pzx  - symmetric traceless viscous stress tensor
   meq[16]= 0.0;                  //mx   - cubic order moment
   meq[17]= 0.0;                  //my   - cubic order moment
   meq[18]= 0.0;                  //mz   - cubic order moment
}

// Translate entries 1 index position to right
__global__ void StreamRight(float *Fin, float *Fout, int Lx, int Mx, int My) {
  int i = threadIdx.x;  // Lattice i-index
  int j =  blockIdx.x;  // Lattice j-index
  int k =  blockIdx.y;  // Lattice k-index
  long int loc;
  __shared__ float line[mL];
  loc = i + Mx*(j + My*k);
  if (i<Lx)
     line[i] = Fin[loc];
  __syncthreads();
  if (i>0)
     Fout[loc] = line[i-1];
  else
     Fout[loc] = line[0];
}

// Translate entries 1 index position to left
__global__ void StreamLeft(float *Fin, float *Fout, int Lx, int Mx, int My) {
  int i = threadIdx.x;  // Lattice i-index
  int j =  blockIdx.x;  // Lattice j-index
  int k =  blockIdx.y;  // Lattice k-index
  long int loc;
  __shared__ float line[mL];
  loc = i + Mx*(j + My*k);
  if (i<Lx)
     line[i] = Fin[loc];
  __syncthreads();
  if  (i<Lx-1)
      Fout[loc] = line[i+1];
  else
      Fout[loc] = line[Lx-1];
}

// Translate entries 1 index position up
__global__ void StreamUp(float *Fin, float *Fout, int Lx, int Ly, int Lz, int Mx, int My) {
  int tix=threadIdx.x;
  int tiy=threadIdx.y;
  int i=threadIdx.x + blockIdx.x*blockDim.x; // Lattice i-index
  int j=threadIdx.y + blockIdx.y*blockDim.y; // Lattice j-index
  int k,MxMy;
  long int loc0,loc1,locBC;
  __shared__ float tile[mB][mB+1];
  MxMy = Mx*My;
  for (k=0; k<Lz; k++) {
    locBC = i + MxMy*k; loc0 = locBC + Mx*(j-1); loc1 = loc0 + Mx;
    if ((i<Lx) & (j<Ly)) {
      if (j>0)
        tile[tix][tiy] = Fin[loc0];
      else
        tile[tix][tiy] = Fin[locBC];
    }
    __syncthreads();
    if ((i<Lx) & (j<Ly))
         Fout[loc1] = tile[tix][tiy];
  }
}
// Translate entries 1 index position down
__global__ void StreamDown(float *Fin, float *Fout, int Lx, int Ly, int Lz, int Mx, int My) {
  int tix=threadIdx.x;
  int tiy=threadIdx.y;
  int i=threadIdx.x + blockIdx.x*blockDim.x; // Lattice i-index
  int j=threadIdx.y + blockIdx.y*blockDim.y; // Lattice j-index
  int k,MxMy;
  long int loc0,loc1,locBC;
  __shared__ float tile[mB][mB+1];
  MxMy = Mx*My;
  for (k=0; k<Lz; k++) {
    locBC = i + Mx*j + MxMy*k; loc1 = locBC; loc0 = loc1 + Mx;
    if ((i<Lx) & (j<Ly)) {
      if (j<Ly-1)
        tile[tix][tiy] = Fin[loc0];
      else
        tile[tix][tiy] = Fin[locBC];
    }
    __syncthreads();
    if ((i<Lx) & (j<Ly))
        Fout[loc1] = tile[tix][tiy];
  }
}
// Translate entries 1 index position up
__global__ void StreamFront(float *Fin, float *Fout, int Lx, int Ly, int Lz, int Mx, int My) {
  int tix=threadIdx.x;
  int tiy=threadIdx.y;
  int i=threadIdx.x + blockIdx.x*blockDim.x; // Lattice i-index
  int k=threadIdx.y + blockIdx.y*blockDim.y; // Lattice k-index
  int j,MxMy;
  long int loc0,loc1,locBC;
  __shared__ float tile[mB][mB+1];
  MxMy = Mx*My;
  for (j=0; j<Ly; j++) {
    locBC = i + Mx*j; loc0 = locBC + MxMy*(k-1); loc1 = loc0 + MxMy;
    if ((i<Lx) & (k<Lz)) {
      if (k>0)
        tile[tix][tiy] = Fin[loc0];
      else
        tile[tix][tiy] = Fin[locBC];
    }
    __syncthreads();
    if ((i<Lx) & (k<Lz))
      Fout[loc1] = tile[tix][tiy];
  }
}
// Translate entries 1 index position down
__global__ void StreamBack(float *Fin, float *Fout, int Lx, int Ly, int Lz, int Mx, int My) {
  int tix=threadIdx.x;
  int tiy=threadIdx.y;
  int i=threadIdx.x + blockIdx.x*blockDim.x; // Lattice i-index
  int k=threadIdx.y + blockIdx.y*blockDim.y; // Lattice k-index
  int j,MxMy;
  long int locBC,loc0,loc1;
  __shared__ float tile[mB][mB+1];
  MxMy = Mx*My;
  for (j=0; j<Ly; j++) {
    locBC = i + Mx*j + MxMy*k; loc1 = locBC; loc0 = loc1 + MxMy;
    if ((i<Lx) & (k<Lz)) {
      if (k<Lz-1)
        tile[tix][tiy] = Fin[loc0];
      else
        tile[tix][tiy] = Fin[locBC];
    }
    __syncthreads();
    if ((i<Lx) & (k<Lz))
        Fout[loc1] = tile[tix][tiy];
  }
}

//  Initialize PDFs
__global__ void initF(int nLBMmodel, int nPDF, int mQ, float tnow, int Mx, int My,
                      int Fsize, int Qsize, float* Phi, int *G, float *Q, float *F) {
  int il = threadIdx.x;  // Lattice i-index
  int jl =  blockIdx.x;  // Lattice j-index
  int kl =  blockIdx.y;  // Lattice k-index
  long int loc0,loc;
  int n;
  float f[38],q[8];
  loc0 = il + Mx*(jl + My*kl);
  if (G[loc0]<0) {
    loc = loc0;
    for (n=0; n<nPDF; n++) {
      F[loc]=0.; loc+=Fsize;
    }
    return;
  }
  loc = loc0;
  for (n=0; n<mQ; n++) {
    q[n]=Q[loc]; loc+=Qsize;
  };
  switch (nLBMmodel) {
    case D3Q19:
    case D3Q19MRT: {
      ComputeEqF_D3Q19(q,f);
      break;
    }
  }
  loc = loc0;
  for (n=0; n<nPDF; n++) {
    F[loc]=f[n]; loc+=Fsize;
  }
}

// Set F values for an inactive site
// Modify F values from most recent time step (stored in f0) to impose flow BC. Place result in f1
__device__ void SetF_InactiveSite(int nPDF, float *f) {
  int i;
  f[0]=1.0;
  for (i=1; i<nPDF; i++) f[i]=0.;
}

// Modify F values from most recent time step (stored in f0) to impose flow BC. Place result in f1
__device__ void SetF_FlowBC(int nLBMmodel, float *f0, float *f1, float *feq, int bcnrml, float *qbc) {
  int n,bc,nrml;
  float q[maxQ];
  float rho,rhov,rhow; // ,rhou,deltarhov,deltarhow;
  nrml = bcnrml % 100; bc = bcnrml - nrml;
  switch (nLBMmodel) {
    case D3Q19:
    case D3Q19MRT:   {                          // Fluid flow, impose continuum values
    // velD3Q19[19]={OOO,OOP,OOM,OPO,OPP,OPM,OMO,OMP,OMM,POO,POP,POM,PPO,PMO,MOO,MOP,MOM,MPO,MMO};
    //                 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
      switch (bc) {
        case IMPOSED_INFLOW: {
          for (n=0; n<4; n++) q[n] = qbc[n];
          // q[0] contains desired pressure
          // q[1] contains desired normal velocity. Project on Cartesian components
          // vec[0] - z component, vec[1] - y component, vec[2] - x component
          q[2]=q[1]*vecD3Q19[nrml][1] ; q[3]=q[1]*vecD3Q19[nrml][0]; q[1]=q[1]*vecD3Q19[nrml][2];
          ComputeEqF_D3Q19(q,f1);
          //if (nrml==1)
          //printf("Setting imposed inflow q[0:3]=%f %f %f %f bcnrml=%d bc=%d nrml=%d\n",q[0],q[1],q[2],q[3],bcnrml,bc,nrml);
          break;
        }
        case IMPOSED_OUTFLOW: {
          for (n=0; n<4; n++) q[n] = qbc[n];
          // q[0] contains desired pressure
          // q[1] contains desired normal velocity. Project on Cartesian components
          q[2]=-q[1]*vecD3Q19[nrml][1] ; q[3]=-q[1]*vecD3Q19[nrml][0]; q[1]=-q[1]*vecD3Q19[nrml][2];
          ComputeEqF_D3Q19(q,f1);
          //if (nrml==2)
          //printf("Setting imposed outflow q[0:3]=%f %f %f %f bcnrml=%d bc=%d nrml=%d f1=%6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f %6.4f \n",q[0],q[1],q[2],q[3],bcnrml,bc,nrml,f1[0],f1[1],f1[2],f1[3],f1[4],f1[5],f1[6],f1[7],f1[8],f1[9],f1[10],f1[11],f1[12],f1[13],f1[14],f1[15],f1[16],f1[17],f1[18]);
          break;
        }
        case xINFLOW: {
          for (n=0; n<4; n++) q[n] = qbc[n];
          rho=q[0]; rhov=q[2]; rhow=q[3];
          q[1] = rho - f0[0] - 2*(f0[2]+f0[5]+f0[8]+f0[11]+f0[16]) - f0[3] - f0[6] - f0[9] - f0[12] - f0[13] - f0[14] - f0[17] - f0[18];
          ComputeEqF_D3Q19(q,feq);
          f1[1]  = rho - f0[0] - f0[2] - f0[3] - 2*(f0[5]+f0[8]+f0[11]+f0[16]) - f0[6] - f0[9] - f0[12] - f0[13] - f0[14] - f0[17] - f0[18] - feq[4] + feq[5] - feq[7] + feq[8] - feq[10] + feq[11] - feq[15] + feq[16];
          f1[4]  =(rhov - f0[3] + f0[6] + 2*f0[8] - f0[12] + f0[13] - f0[17] + f0[18] + feq[4] - feq[5] + feq[7] - feq[8])/2.;
          f1[7]  =(-rhov + f0[3] + 2*f0[5] - f0[6] + f0[12] - f0[13] + f0[17] - f0[18] + feq[4] - feq[5] + feq[7] - feq[8])/2.;
          f1[10] =(rhow - f0[9] - f0[12] - f0[13] + f0[14] + 2*f0[16] + f0[17] + f0[18] + feq[10] - feq[11] + feq[15] - feq[16])/2.;
          f1[15] =(-rhow + f0[9] + 2*f0[11] + f0[12] + f0[13] - f0[14] - f0[17] - f0[18] + feq[10] - feq[11] + feq[15] - feq[16])/2.;
          break;
        }
        case xOUTFLOW: {
          for (n=0; n<4; n++) q[n] = qbc[n];
          rho=q[0]; rhov=q[2]; rhow=q[3];
          q[1] = -rho + f0[0] + 2*(f0[1]+f0[4]+f0[7]+f0[10]+f0[15]) + f0[3] + f0[6] + f0[9] + f0[12] + f0[13] + f0[14] + f0[17] + f0[18];
          ComputeEqF_D3Q19(q,feq);
          f1[2] = rho - f0[0] - f0[1] - f0[3] - 2*f0[4] - f0[6] - 2*f0[7] - f0[9] - 2*f0[10] - f0[12] - f0[13] - f0[14] - 2*f0[15] - f0[17] - f0[18] + feq[4] - feq[5] + feq[7] - feq[8] + feq[10] - feq[11] + feq[15] - feq[16];
          f1[5] = (rhov - f0[3] + f0[6] + 2*f0[7] - f0[12] + f0[13] - f0[17] + f0[18] - feq[4] + feq[5] - feq[7] + feq[8])/2.;
          f1[8] = (-rhov + f0[3] + 2*f0[4] - f0[6] + f0[12] - f0[13] + f0[17] - f0[18] - feq[4] + feq[5] - feq[7] + feq[8])/2.;
          f1[11] = (rhow - f0[9] - f0[12] - f0[13] + f0[14] + 2*f0[15] + f0[17] + f0[18] - feq[10] + feq[11] - feq[15] + feq[16])/2.;
          f1[16] = (-rhow + f0[9] + 2*f0[10] + f0[12] + f0[13] - f0[14] - f0[17] - f0[18] - feq[10] + feq[11] - feq[15] + feq[16])/2.;
          break;
        }
        default: {
          printf("Warning: undefined flow BC: bcnrml=%d bc=%d nrml=%d\n",bcnrml,bc,nrml);
          for (int n=0; n<19; n++) f1[n]=f0[n];
        }
      }
      break;
    }
  }
}

// Modify F values from most recent time step (stored in f0) to impose wall BC. Place result in f1
__device__ void SetF_WallBC(int nLBMmodel, float *f0, float *f1, int wBC) {
  switch (nLBMmodel) {
    case D3Q19:
    case D3Q19MRT: {  // Fluid flow, implement bounce-back
// velD3Q19[19]={OOO,OOP,OOM,OPO,OPP,OPM,OMO,OMP,OMM,POO,POP,POM,PPO,PMO,MOO,MOP,MOM,MPO,MMO};
//                 0   1   2   3   4   5   6   7   8   9  10  11  12  13  14  15  16  17  18
      f1[ 0] = f0[ 0];
      f1[ 1] = f0[ 2];
      f1[ 2] = f0[ 1];
      f1[ 3] = f0[ 6];
      f1[ 4] = f0[ 8];
      f1[ 5] = f0[ 7];
      f1[ 6] = f0[ 3];
      f1[ 7] = f0[ 5];
      f1[ 8] = f0[ 4];
      f1[ 9] = f0[14];
      f1[10] = f0[16];
      f1[11] = f0[15];
      f1[12] = f0[18];
      f1[13] = f0[17];
      f1[14] = f0[ 9];
      f1[15] = f0[11];
      f1[16] = f0[10];
      f1[17] = f0[13];
      f1[18] = f0[12];
      break;
    }
  }
}

// Set PDF boundary conditions, and also transfer from newest values in F1 to F0.
// F0,F1 have same content at end of this time ste
__global__ void SetBC(int nLBMmodel, int nPDF, int Fsize, int nQ, int Qsize,
                      int Mx, int My, int Lx, int Ly, int Lz,
                      float *F0, float *F1, int *G, float *Q1) {
  int i = threadIdx.x;  // Lattice i-index
  int j =  blockIdx.x;  // Lattice j-index
  int k =  blockIdx.y;  // Lattice k-index
  long int loc,loc0;
  int n,bc;
  int wall=100;  // 1-99 to encode normal directions
  float f0[maxF],f1[maxF],feq[maxF],qbc[maxQ];
  loc0 = i + Mx*(j + My*k);
  loc=loc0;
  // Place F1 values from most recent time step into f0,f1.
  for (n=0; n<nPDF; n++) {
    f1[n]=f0[n]=F1[loc]; loc+=Fsize;
  }
  bc = G[loc0];
  // Place boundary continuum values from Q into q
  loc=loc0;
  for (n=0; n<nQ; n++)
    { qbc[n]=Q1[loc]; loc+=Qsize;}
  // Second-order at mid-lattice placed wall, full bounceback. All inactive/wall nodes are bounced back
  if ((bc<0) || ((bc>0) && (bc<wall))) {
    SetF_WallBC(nLBMmodel,f0,f1,bc);
  }
  // Flow boundary conditions
  if (bc>=wall) {
    SetF_FlowBC(nLBMmodel,f0,f1,feq,bc,qbc);
  }
  // Interior sites with bc=0 do not require any processing
  loc=loc0;
  for (n=0; n<nPDF; n++)
    // Transfer f1 values (from BC's or most recent time step) to F0 and F1 in preparation for next time step
    {F0[loc]=F1[loc]=f1[n]; loc+=Fsize;}
}

__global__ void CollideD3Q19(int nPDF, int Fsize, float tau,
                             int Mx, int My,
                             float *F1) {
  int i = threadIdx.x;  // Lattice i-index
  int j =  blockIdx.x;  // Lattice j-index
  int k =  blockIdx.y;  // Lattice k-index
  long int loc,loc0;
  int n;
  float f0[19],f1[19],feq[19],qeq[4];
  loc0 = i + Mx*(j + My*k);
  loc=loc0;
  for (n=0; n<nPDF; n++)
    {f0[n]=F1[loc]; loc+=Fsize;}
  ComputeEqQ_D3Q19(f0,qeq);
  ComputeEqF_D3Q19(qeq,feq);
  for (n=0; n<nPDF; n++)
      f1[n] = f0[n] + (feq[n] - f0[n])/tau;
  loc=loc0;
  for (n=0; n<nPDF; n++)
    {F1[loc]=f1[n]; loc+=Fsize;}
}

__global__ void CollideD3Q19MRT(int nPDF, int Fsize,
                                int Mx, int My,
                                float *Phi, float *Tau,
                                float *F1) {
  int il= threadIdx.x;  // Lattice i-index
  int jl=  blockIdx.x;  // Lattice j-index
  int kl=  blockIdx.y;  // Lattice k-index
  long int loc,loc0;
  int i,j,ij,n;
  float f0[19],f1[19],meq[19],m0[19],qeq[4];
  float sum;
  loc0 = il + Mx*(jl + My*kl);
  loc=loc0;
  for (n=0; n<nPDF; n++)
    {f0[n]=F1[loc]; loc+=Fsize;}
  //Compute equilibrium moments
  ComputeEqQ_D3Q19(f0,qeq); ComputeEqM_D3Q19(qeq,meq);
  //Compute current moments by linear transformation of f
  ij=0;
  for (i=0; i<nPDF; i++) {
    sum=0.0;
    for (j=0; j<nPDF; j++) {
      sum+=Phi[ij]*f0[j];
      ij++;
    }
    m0[i]=sum;
  }
  //compute new PDFs
  for (i=0; i<nPDF; i++) {
      sum=0.0;
      for (j=0; j<nPDF; j++) {
        ij=j*19+i;
        sum+=Tau[j]*(meq[j]-m0[j])*Phi[ij];
      }
      f1[i] = f0[i] + sum;
  }
  loc=loc0;
  for (n=0; n<nPDF; n++)
    {F1[loc]=f1[n]; loc+=Fsize;}
}

__global__ void ContVD3Q19(int nPDF, int mQ, int Fsize, int Qsize,
                           int Mx, int My, float *F1, float *Q1) {
  int i = threadIdx.x;  // Lattice i-index
  int j =  blockIdx.x;  // Lattice j-index
  int k =  blockIdx.y;  // Lattice k-index
  long int loc,loc0;
  int m,n;
  float f1[19],qeq[4];
  loc0 = i + Mx*(j + My*k);
  loc=loc0;
  for (n=0; n<nPDF; n++)
    { f1[n]=F1[loc]; loc+=Fsize;}
  ComputeEqQ_D3Q19(f1,qeq);
  loc=loc0;
  for (m=0; m<mQ; m++)
    { Q1[loc]=qeq[m]; loc+=Qsize;}
}
// Host routines

// 1. Initialization routines

// Define velocity codes, number of streaming steps, weights
void LBMvelinit(int nLBMmodel, int nPDF, int *vel, int *s, float *w) {
  int n;
  switch (nLBMmodel) {
    case D3Q19:
    case D3Q19MRT: {
      for (n=0; n<nPDF; n++) {
        vel[n]=velD3Q19[n]; s[n]=StreamD3Q19[n]; w[n]=wD3Q19[n]/WD3Q19;
      }
      break;
    }
  }
}

void LBMinitF(int nLBMmodel, int nPDF, int mQ, float tnow, int *L, int *M, int Fsize, int Qsize,
              float *Phi, int *G, float *Q, float *F) {
  dim3 dG, dB;
  // A CUDA block processes a line of data along x-axis, one thread per lattice site
  dB.x = L[0]; dB.y = 1;    dB.z = 1;
  // The CUDA grid corresponds to dimensions along y and z
  dG.x = L[1]; dG.y = L[2];
  initF<<<dG,dB>>>(nLBMmodel,nPDF,mQ,tnow,M[0],M[1],Fsize,Qsize,Phi,G,Q,F);
}

void LBMinitMRT(float *p, float *Tau, float *Phi) {
  float sum,visc,svsc;
  int i,j,ij;
  float ex2[19],ey2[19],ez2[19],e2[19],s[19],tau[19];
  float phi[19*19];
  //Relaxation matrix eigenvalues
  //   viscosity gives momentum relaxation time
  visc=p[24]; svsc=1.0/(3.0*visc+0.5);
  //   load diagonal matrix of relaxation time scales
  s[ 0]=0.0 ; s[ 1]=1.19; s[ 2]=1.4 ; s[ 3]=0.0 ; s[ 4]=1.2 ;
  s[ 5]=0.0 ; s[ 6]=1.20; s[ 7]=0.0 ; s[ 8]=1.2 ; s[ 9]=svsc;
  s[10]=1.4 ; s[11]=svsc; s[12]=1.4 ; s[13]=svsc; s[14]=svsc;
  s[15]=1.0 ; s[16]=1.98; s[17]=1.98; s[18]=1.98;
  //Initialize Phi matrix
  for (i=0;i<19;i++) {
    ex2[i]=ex[i]*ex[i];
    ey2[i]=ey[i]*ey[i];
    ez2[i]=ez[i]*ez[i];
    e2[i]=ex2[i]+ey2[i]+ez2[i];
  }
  // m = (Phi^T) f. Compute rows of Phi^T=(phi0^T ... phiN^T)
  ij=0;
  //row 0
  for (i=0;i<19;i++) {
    phi[ij]=1.0;
    ij++;
  }
  //row 1
  for (i=0;i<19;i++) {
    phi[ij]=19.0*e2[i]-30.0;
    ij++;
  }
  //row 2
  for (i=0;i<19;i++) {
    phi[ij]=0.5*(21.0*e2[i]*e2[i]-53.0*e2[i]+24.0);
    ij++;
  }
  //row 3
  for (i=0;i<19;i++) {
    phi[ij]=ex[i];
    ij++;
  }
  //row 4
  for (i=0;i<19;i++) {
    phi[ij]=ex[i]*(5.0*e2[i]-9.0);
    ij++;
  }
  //row 5
  for (i=0;i<19;i++) {
    phi[ij]=ey[i];
    ij++;
  }
  //row 6
  for (i=0;i<19;i++) {
    phi[ij]=ey[i]*(5.0*e2[i]-9.0);
    ij++;
  }
  //row 7
  for (i=0;i<19;i++) {
    phi[ij]=ez[i];
    ij++;
  }
  //row 8
  for (i=0;i<19;i++) {
    phi[ij]=ez[i]*(5.0*e2[i]-9.0);
    ij++;
  }
  //row 9
  for (i=0;i<19;i++) {
    phi[ij]=3.0*ex2[i]-e2[i];
    ij++;
  }
  //row 10
  for (i=0;i<19;i++) {
    phi[ij]=(3.0*e2[i]-5.0)*(3.0*ex2[i]-e2[i]);
    ij++;
  }
  //row 11
  for (i=0;i<19;i++) {
    phi[ij]=ey2[i]-ez2[i];
    ij++;
  }
  //row 12
  for (i=0;i<19;i++) {
    phi[ij]=(3.0*e2[i]-5.0)*(ey2[i]-ez2[i]);
    ij++;
  }
  //row 13
  for (i=0;i<19;i++) {
    phi[ij]=ex[i]*ey[i];
    ij++;
  }
  //row 14
  for (i=0;i<19;i++) {
    phi[ij]=ey[i]*ez[i];
    ij++;
  }
  //row 15
  for (i=0;i<19;i++) {
    phi[ij]=ez[i]*ex[i];
    ij++;
  }
  //row 16
  for (i=0;i<19;i++) {
    phi[ij]=ex[i]*(ey2[i]-ez2[i]);
    ij++;
  }
  //row 17
  for (i=0;i<19;i++) {
    phi[ij]=ey[i]*(ez2[i]-ex2[i]);
    ij++;
  }
  //row 18
  for (i=0;i<19;i++) {
    phi[ij]=ez[i]*(ex2[i]-ey2[i]);
    ij++;
  }
  // Scale relaxation matrix eigenvalues by squared norm of phi row vectors to get tau
  for (i=0;i<19;i++) {
    sum=0.0; ij=19*i;
    for (j=0;j<19;j++) {
      sum+=phi[ij]*phi[ij]; ij++;
    }
    tau[i]=s[i]/sum;
  }
  cudaMemcpy(Tau,tau,19*sizeof(float), cudaMemcpyHostToDevice );
  cudaMemcpy(Phi,phi,19*19*sizeof(float), cudaMemcpyHostToDevice );
}

// 2. Boundary value routines

void LBMSetBC(int nLBMmodel, int nPDF, int Fsize, int nQ, int Qsize, int *L, int *M,
              long int *mF, float *F0, float *F1, int *G, float *QBC) {
  dim3 dG, dB;
  // Go through entire volume
  // A CUDA block processes a line of data along x-axis, one thread per lattice site
  dB.x = L[0]; dB.y = 1;    dB.z = 1;
  // The CUDA grid corresponds to dimensions along y and z
  dG.x = L[1]; dG.y = L[2];
  SetBC<<<dG,dB>>>(nLBMmodel,nPDF,Fsize,nQ,Qsize,M[0],M[1],L[0],L[1],L[2],F0,F1,G,QBC);
}

// 3. Streaming routines

void stream(int dir, int *L, int *M, float *F0, float *F1) {
  dim3 dG, dB;
  int nBx,nBy,nBz;
  switch (dir) {
    case 0: {
      dB.x = L[0]; dB.y = 1; dB.z = 1;
      dG.x = L[1]; dG.y = L[2];
      StreamRight<<<dG,dB>>>(F0, F1, L[0], M[0], M[1]); // stream in positive x-direction
      break;
    }
    case 1: {
      dB.x = L[0]; dB.y = 1; dB.z = 1;
      dG.x = L[1]; dG.y = L[2];
      StreamLeft<<<dG,dB>>>(F0, F1, L[0], M[0], M[1]); // stream in negative x-direction
      break;
    }
    case 2: {
      dB.x = mB; dB.y = mB; dB.z = 1;
      nBx = L[0]/mB; if (L[0] % mB > 0) nBx++;
      nBy = L[1]/mB; if (L[1] % mB > 0) nBy++;
      dG.x = nBx; dG.y = nBy;
      StreamUp<<<dG,dB>>>(F0, F1, L[0], L[1], L[2], M[0], M[1]);  // stream in positive y-direction
      break;
    }
    case 3: {
      dB.x = mB; dB.y = mB; dB.z = 1;
      nBx = L[0]/mB; if (L[0] % mB > 0) nBx++;
      nBy = L[1]/mB; if (L[1] % mB > 0) nBy++;
      dG.x = nBx; dG.y = nBy;
      StreamDown<<<dG,dB>>>(F0, F1, L[0], L[1], L[2], M[0], M[1]);  // stream in negative y-direction
      break;
    }
    case 4: {
      dB.x = mB; dB.y = mB; dB.z = 1;
      nBx = L[0]/mB; if (L[0] % mB > 0) nBx++;
      nBz = L[2]/mB; if (L[2] % mB > 0) nBz++;
      dG.x = nBx; dG.y = nBz;
      StreamFront<<<dG,dB>>>(F0, F1, L[0], L[1], L[2], M[0], M[1]);  // stream in positive z-direction
      break;
    }
    case 5: {
      dB.x = mB; dB.y = mB; dB.z = 1;
      nBx = L[0]/mB; if (L[0] % mB > 0) nBx++;
      nBz = L[2]/mB; if (L[2] % mB > 0) nBz++;
      dG.x = nBx; dG.y = nBz;
      StreamBack<<<dG,dB>>>(F0, F1, L[0], L[1], L[2], M[0], M[1]);  // stream in negative z-direction
      break;
    }
  }
}

void LBMstream(int nLBMmodel, int nPDF, int *vel, int *nStream, int *L, int *M, long int *mF, float *F0, float *F1) {
  int d,n;
  int mask;
  int nS[64];
  mask=1;
  for (n=1; n<nPDF; n++) nS[n]=nStream[n];
  for (d=0; d<6; d++) { // Loop over directions
    for (n=1; n<nPDF; n++) { // Loop over PDFs
      if (vel[n] & mask) {
        if (nS[n] % 2 == 1)
          stream(d,L,M,F0+mF[n],F1+mF[n]);
        else
          stream(d,L,M,F1+mF[n],F0+mF[n]);
        nS[n]--;
      }
    }
    mask = mask << 1;
  }
}

// 4. Collision (relaxation) routines

void LBMcollide(int nLBMmodel, int nPDF, int Fsize,
                int *vel, float *w, float tau,
                float *Tau, float *Phi,
                int *L, int *M,
                float *F1) {
  dim3 dG, dB;
  dB.x = L[0]; dB.y = 1; dB.z = 1;
  dG.x = L[1]; dG.y = L[2];
  switch (nLBMmodel) {
    case D3Q19: {
      CollideD3Q19<<<dG,dB>>>(nPDF,Fsize,tau,M[0],M[1],F1);
      break;
    }
    case D3Q19MRT: {
      CollideD3Q19MRT<<<dG,dB>>>(nPDF,Fsize,M[0],M[1],Phi,Tau,F1);
      break;
    }
  }
}

void LBMcontV(int nLBMmodel, int nPDF, int nCont, int Fsize, int Qsize,
              int *L, int *M, float *F1, float *Q1) {
  dim3 dG, dB;
  dB.x = L[0]; dB.y = 1; dB.z = 1;
  dG.x = L[1]; dG.y = L[2];
  switch (nLBMmodel) {
    case (D3Q19):
    case (D3Q19MRT): {
      ContVD3Q19<<<dG,dB>>>(nPDF,nCont,Fsize,Qsize,M[0],M[1],F1,Q1);
      break;
    }
  }
}

// Entry points for external calls
/* Known LBM models:
Code  Name       nPDF  nCont
0     (reserved)
1     D3Q19      19    4
2     D3Q19MRT   19    4
*/
extern "C" int microLBM(float *p, float *q, int *g, float *qBC, float *f, int Mx, int My, int Mz) {
  int nt,nLatticeSites;
  static int Fsize,Qsize;
  long int nBytesGeom,nBytesPDF,nBytesCont,nBytes,nMB,nBytesQBC,nBytesTau;
  int nPDF,nCont,nPDFmodels[4] = {0,19,19,19},nContmodels[4]={0,4,4,4};
  int i,nLBMmodel,nQ,nBound;
  int L[3],M[3];
  static long int offF[27];
  static int LBMvel[64],LBMstr[64];
  static float LBMw[64];
  int InitializePDFs,CopyPDFsToHost;
  static float *F0,*F1,*Q1,*QBC,*Tau,*Phi;
  static int *G;
  int DeviceMemory=5376; // Max device memory in MB (3GB)
  int MegaByte=1024*1024;
  float tnow,tau;
  //----------------------------------------------------------------------------
  /* Load run parameters */
  nLBMmodel = (int) p[0];                    // LBM model
  for (i=0; i<3; i++) L[i] = (int) p[i+1];   // p[1:3] = lattice dimensions
  for (i=0; i<3; i++) M[i] = (int) p[i+4];   // p[4:6] = lattice memory space
  nt = (int) p[7];                           // number of time steps
  InitializePDFs = (int) p[8];               // flag to initialize flow field to this density
  tau = p[9];                                // Collide operator relaxation time
  nBound = (int) p[10];                       // Nr. of boundary condition sets
  nQ = (int) p[11];                          // Nr. of Q components (e.g., 4 for Navier-Stokes)
  CopyPDFsToHost = (int) p[12];              // flag to return PDFs
  tnow = p[25];                              // current time
  //----------------------------------------------------------------------------
  nPDF = nPDFmodels[nLBMmodel]; nCont = nContmodels[nLBMmodel];
  /* Allocate device memory */
  nLatticeSites=1; for (i=0; i<3; i++) nLatticeSites *= M[i];
  nBytesGeom = nLatticeSites      *sizeof(int);
  nBytesPDF  = nLatticeSites*nPDF *sizeof(float);
  nBytesCont = nLatticeSites*nCont*sizeof(float);
  nBytesQBC  = nQ*nBound*sizeof(float);
  nBytesTau  = nPDF*sizeof(float);
  // Total space
  nBytes = nBytesGeom + 2*nBytesPDF + nBytesCont; nMB = nBytes/1024/1024;
  if (nMB>DeviceMemory) {
    printf("Not enough memory on GPU device\n");
    printf("M=(%d,%d,%d) \n",M[0],M[1],M[2]);
    printf("nLatticeSites=%d \n",nLatticeSites);
    printf("nBytesGeom=%d MB\n",(int) (nBytesGeom/MegaByte));
    printf("nBytesPDFs=%d MB\n",(int) (2*nBytesPDF/MegaByte));
    printf("nBytesCont=%d MB\n",(int) (nBytesCont/MegaByte));
    exit(1);
  }
  // for (int n=0; n<24; n++) printf("QBC[%d]=%f\n",n,qBC[n]);
  // Allocate space on GPU device memory:
  cudaMalloc( (void **) &G , nBytesGeom);    // geometry flags
  cudaMalloc( (void **) &F0, nBytesPDF );    // old PDFs
  cudaMalloc( (void **) &F1, nBytesPDF );    // new PDFs
  cudaMalloc( (void **) &Q1, nBytesCont);    // continuum parameters from LBM PDFs
  cudaMalloc( (void **) &QBC,nBytesQBC);
  cudaMalloc( (void **) &Tau,nBytesTau);
  cudaMalloc( (void **) &Phi,nPDF*nBytesTau);
  // Size of F,Q arrays
  Fsize = nLatticeSites; Qsize = nLatticeSites;
  // Offsets to each PDF
  for (i=0; i<nPDF; i++) offF[i] = i*nLatticeSites;
  // Load geometry data to device memory
  cudaMemcpy( G, g, nBytesGeom, cudaMemcpyHostToDevice ); // geometry flags
  // Load boundary conditions to device memory
  cudaMemcpy(QBC,qBC,nBytesQBC, cudaMemcpyHostToDevice );
  // Load continuum field values to device memory (always needed in order to transmit boundary conditions)
  cudaMemcpy(Q1,q,nBytesCont, cudaMemcpyHostToDevice );
  // Load lattice model velocity defintions, weights
  LBMvelinit(nLBMmodel,nPDF,LBMvel,LBMstr,LBMw);
  // Carry out lattice model specific initializations
  switch (nLBMmodel) {
    case D3Q19: {
      break;
    }
    case D3Q19MRT: { // Compute moment transformation matrix and transfer to GPU arrays Tau, Phi
      LBMinitMRT(p,Tau,Phi);
      break;
    }
  }
  if (InitializePDFs) {
    // Set initial PDFs from continuum field variables
    LBMinitF(nLBMmodel,nPDF,nCont,tnow,L,M,Fsize,Qsize,Phi,G,Q1,F1);   // F1 initialized
  } else {
    // Load PDF values to device memory
    cudaMemcpy(F1,f,nBytesPDF, cudaMemcpyHostToDevice );
  }
  /* Evolve LBM forward in time */
  for (int n=0; n<nt; n++) {
    LBMSetBC  (nLBMmodel,nPDF,Fsize,nQ,Qsize,L,M,offF,F0,F1,G,Q1);    // F1 modified by BC's placed in F0. F1,F0 have same content
    LBMstream (nLBMmodel,nPDF,LBMvel,LBMstr,L,M,offF,F0,F1);          // F0 streamed to F1
    LBMcollide(nLBMmodel,nPDF,Fsize,LBMvel,LBMw,tau,Tau,Phi,L,M,F1);  // F1 relaxed to F1
    //                                          SRT MRT
  }
  LBMcontV(nLBMmodel,nPDF,nCont,Fsize,Qsize,L,M,F1,Q1);
  cudaMemcpy( q, Q1, nBytesCont, cudaMemcpyDeviceToHost ); // Copy continuum field values
  //printf("Copying updated continuum values from GPU to CPU\n");
  if (CopyPDFsToHost) {
    //printf("nBytesPDF = %d \n",nBytesPDF);
    cudaMemcpy( f, F1, nBytesPDF,  cudaMemcpyDeviceToHost ); // Copy PDF values
  }
  cudaFree(G); cudaFree(F0); cudaFree(F1); cudaFree(Q1); cudaFree(QBC);
  return 0;
}

// Utility routines:

/* FindWalls: Given the array g with boundary codes:
     -1   = inactive node
      0   = active node
      >63 = flow condition node
   modify g such that any inactive node linked by a lattice direction to an active node becomes a wall
   node with the appropriate inward-pointing normal direction. Lattice direction numbering convention is
   encoded in least-significant 6 bits:
*/
extern "C" int FindNormals(float *p, int *g, int Mx, int My, int Mz) {
  int TRUE=1, INTERIOR=0, INACTIVE=-1; // FALSE=0
  int i,j[3],idx[3],n,d,nLBMmodel,NeighborInLattice,nloc,loc;
  int L[3],M[3],fWallSitePrint,fWallSiteDebug,nWallSites,nFlowSites,BCcategory;
  int FoundInteriorNeighbor,InactiveSite,BCFlowSiteNoNormal;
  int DefineWallSiteNormal,DefineFlowSiteNormal;
  // int nEnlarge,k1Enlarge,k2Enlarge,Enlarge;
  nLBMmodel = (int) p[0];                    // LBM model
  for (i=0; i<3; i++) L[i] = (int) p[i+1];   // p[1:3] = lattice dimensions
  for (i=0; i<3; i++) M[i] = (int) p[i+4];   // p[4:6] = lattice memory space
  fWallSitePrint = (int) p[31]; fWallSiteDebug = (int) p[30];
  // k1Enlarge = pLBM(15); k2Enlarge = pLBM(16); nEnlarge = (int) pLBM(17);
  nWallSites=0; nFlowSites=0;
  if (fWallSiteDebug)
    printf("FindNormals: nLBMmodel=%d, L[0]=%4d, L[1]=%4d, L[2]=%4d, M[0]=%4d, M[1]=%4d, M[2]=%4d\n",
            nLBMmodel,L[0],L[1],L[2],M[0],M[1],M[2]);
  for (j[0]=0; j[0]<L[0]; j[0]++) {
    for (j[1]=0; j[1]<L[1]; j[1]++) {
      for (j[2]=0; j[2]<L[2]; j[2]++) {
        loc = j[0] + M[0]*(j[1] + M[1]*j[2]);
        if (g[loc] != INTERIOR) {
           // See if the inward-pointing normal direction must be defined at this site
           // Is this site is inactive? (but perhaps a wall site for which a normal needs to be defined)
           InactiveSite = g[loc]<0;
           // Is this site on a flow condition boundary for which the normal has not yet been defined?
           BCFlowSiteNoNormal = g[loc] % 100 == 0;
           BCcategory = g[loc]/100;
           if (fWallSiteDebug)  printf("Checking site i,j,k=%d %d %d bc=%4d\n",j[0],j[1],j[2],g[loc]);
           switch (nLBMmodel) {
             case (D3Q19):
             case (D3Q19MRT): {
               for (n=1; n<nDirD3Q19; n++) {
                 // Find indices of neighboring lattice sites along each direction
                 for (d=0; d<3; d++) idx[d] = j[d] + dirD3Q19[n][2-d];
                 NeighborInLattice = TRUE;
                 for (d=0; d<3; d++)
                   NeighborInLattice = NeighborInLattice && (idx[d]>=0) && (idx[d]<L[d]);
                 if (fWallSiteDebug)
                   printf("Along direction nr. %d, idx = %4d %4d %d, NeighborInLattice = %d\n",
                           n,idx[0],idx[1],idx[2],NeighborInLattice);
                 if (NeighborInLattice) {
                   // idx[:] define a valid neighbor site
                   nloc = idx[0] + M[0]*(idx[1] + M[1]*idx[2]);
                   // Is this neighbor lattice site is in the interior?
                   FoundInteriorNeighbor = g[nloc]==INTERIOR;
                   // Is this a wall site, i.e. an inactive site with an interior neighbor?
                   // If so, an interior pointing normal direction must be defined
                   DefineWallSiteNormal = InactiveSite && FoundInteriorNeighbor;
                   // Does this flow BC site have an interior neighbor?
                   // If so, an interior pointing normal direction must be defined
                   DefineFlowSiteNormal = BCFlowSiteNoNormal && FoundInteriorNeighbor;
                   if (fWallSiteDebug)
                     printf("FoundInteriorNeighbor=%d, DefineWallSiteNormal=%d, DefineFlowSiteNormal=%d\n",
                             FoundInteriorNeighbor,DefineWallSiteNormal,DefineFlowSiteNormal);
                   if (DefineWallSiteNormal) {
                     g[loc]=indxD3Q19[n]; nWallSites++;
                     if (fWallSitePrint) {
                       printf("Found wall node at i=%4d, j=%4d, k=%4d, n=%2d normal=%2d G=%4d\n",j[0],j[1],j[2],n,nrmlD3Q19[n],g[loc]);
                       printf("Neighbor at i=%4d j=%4d k=%4d with G=%4d\n",idx[0],idx[1],idx[2],g[nloc]);
                     }
                     break;
                   };
                   if (DefineFlowSiteNormal) {
                     g[loc] = BCcategory*100 + indxD3Q19[n]; nFlowSites++;
                     if (fWallSitePrint) {
                       if (BCcategory==1)
                         printf("Have defined inflow normal at i=%4d, j=%4d, k=%4d, n=%2d normal=%2d G=%4d\n",j[0],j[1],j[2],n,nrmlD3Q19[n],g[loc]);
                       else
                         printf("Have defined outflow normal at i=%4d, j=%4d, k=%4d, n=%2d normal=%2d G=%4d\n",j[0],j[1],j[2],n,nrmlD3Q19[n],g[loc]);
                     }
                     break;
                   };
                 }
               }
               // Flow BC sites not identified as connected to an interior site are set inactive
               if( (g[loc]>=100) && (g[loc] % 100 == 0)) g[loc]=INACTIVE;
               break;
             }
           }
        }
      }
    }
  }
  printf("   FindNormals has defined normals at %6d wall sites, %6d flow sites\n",nWallSites,nFlowSites);
  return 0;
}

/* LBMbc: Given an array of incompressible fluid parameters in SI units, modify it to
          LBM computational units
*/
extern "C" int LBMbc(float *p, float *qBC) {
  int i,n,nBound,nQ,nLBMmodel;
  // Variables with 'Ref' suffix are expressed in SI units
  float dxRef,dtRef,pRef,cSoundRef,viscRef,uRef,tau;
  // Load reference parameters
  nLBMmodel = (int) p[0];
  dxRef = p[20]; pRef = p[22]; cSoundRef = p[23]; viscRef = p[24];
  // Compute additional reference quantities
  uRef = cSoundRef/csLBM[nLBMmodel];        // Reference velocity
  dtRef = dxRef/uRef;                       // Reference time
  p[21] = dtRef;
  tau = 0.5*(1.+6.*viscRef*dtRef/dxRef/dxRef); // BGK relaxation time (non-dimensional)
  p[9] = tau;
  // Transform qBC from physical quantities to LBM non-dimensional quantities
  nBound = (int) p[10];                     // Nr. of boundary condition sets
  nQ     = (int) p[11];                     // Nr. of Q components (e.g., 4 for Navier-Stokes)
  for (i=0; i<nBound; i++) {
    qBC[nQ*i] = qBC[nQ*i]/pRef;
    for (n=1; n<nQ; n++) qBC[nQ*i+n] = qBC[nQ*i+n]/uRef;
  };
  return 0;
}
