#ifndef LBM_H_INCLUDED
#define LBM_H_INCLUDED
// LBM direction encoding. 3 letter variable formed from O=zero, P=plus, M=minus
// encode in 6 bits xxyyzz with less significant bit associated with P direction
// This is not the standard enumeration of LBM directions found in the literature,
// but is efficient for implementation of the streaming step on a GPU (see LBMstream routine)
// LDM model indices                        D3Q19
#define OOO  0  // 00 00 00                 0
#define OOP  1  // 00 00 01                 1
#define OOM  2  // 00 00 10                 2
                // 00 00 11  not used
#define OPO  4  // 00 01 00                 3
#define OPP  5  // 00 01 01                 4
#define OPM  6  // 00 01 10                 5
                // 00 01 11  not used
#define OMO  8  // 00 10 00                 6
#define OMP  9  // 00 10 01                 7
#define OMM 10  // 00 10 10                 8
                // 00 10 11  not used
                // 00 11 00  not used
                // 00 11 01  not used
                // 00 11 10  not used
                // 00 11 11  not used
#define POO 16  // 01 00 00                 9
#define POP 17  // 01 00 01                10
#define POM 18  // 01 00 10                11
                // 01 00 11  not used
#define PPO 20  // 01 01 00                12
#define PPP 21  // 01 01 01                not used
#define PPM 22  // 01 01 10                not used
                // 01 01 11  not used
#define PMO 24  // 01 10 00                13
#define PMP 25  // 01 10 01                not used
#define PMM 26  // 01 10 10                not used
                // 01 10 11  not used
                // 01 11 00  not used
                // 01 11 01  not used
                // 01 11 10  not used
                // 01 11 11  not used
#define MOO 32  // 10 00 00                14
#define MOP 33  // 10 00 01                15
#define MOM 34  // 10 00 10                16
                // 10 00 11  not used
#define MPO 36  // 10 01 00                17
#define MPP 37  // 10 01 01                not used
#define MPM 38  // 10 01 10                not used
                // 10 01 11  not used
#define MMO 40  // 10 10 00                18
#define MMP 41  // 10 10 01                not used
#define MMM 42  // 10 10 10                not used
                // 10 10 11  not used
                // 10 11 00  not used
                // 10 11 01  not used
                // 10 11 10  not used
                // 10 11 11  not used
                // 11 xx xx  not used

// Constants
// Buffer size for x-direction streaming
#define mL 512
// Tile size for y,z-directions streaming
#define mB 16
// sqrt(2.)
#define S2 0.70710678118654746
// Maximum PDF components
#define maxF 24
// Maximum Q components
#define maxQ 12

// LBM model codes
#define D3Q19 1
#define D3Q19MRT 2

// Number of directions in LBM models
#define nDirD3Q19 19

// Boundary condition codes
#define NR_BOUNDARIES 5
#define BCinterior 0
#define BCinactive 1
#define BCwall     2
#define BCinflow   3
#define BCoutflow  4
#define IMPOSED_INFLOW 100
#define IMPOSED_OUTFLOW 200
#define COMPUTED_INFLOW 300
#define COMPUTED_OUTFLOW 400
#define xINFLOW  1100
#define xOUTFLOW 1200

static float    WD3Q19    =36.;
//                          0   1   2   3   4   5   6   7   8   9   10  11  12  13  14  15  16  17  18
static float    wD3Q19[19]={12., 2., 2., 2., 1., 1., 2., 1., 1., 2., 1., 1., 1., 1., 2., 1., 1., 1., 1.};
static int    velD3Q19[19]={OOO,OOP,OOM,OPO,OPP,OPM,OMO,OMP,OMM,POO,POP,POM,PPO,PMO,MOO,MOP,MOM,MPO,MMO};
static int StreamD3Q19[19]={  0,  1,  1,  1,  2,  2,  1,  2,  2,  1,  2,  2,  2,  2,  1,  2,  2,  2,  2};

// Multiple relaxation time direction vectors
static float ez[19]={ 0., 0., 0., 0., 0., 0., 0., 0., 0., 1., 1., 1., 1., 1.,-1.,-1.,-1.,-1.,-1.};
static float ey[19]={ 0., 0., 0., 1., 1., 1.,-1.,-1.,-1., 0., 0., 0., 1.,-1., 0., 0., 0., 1.,-1.};
static float ex[19]={ 0., 1.,-1., 0., 1.,-1., 0., 1.,-1., 0., 1.,-1., 0., 0., 0., 1.,-1., 0., 0.};
//                   OOO,OOP,OOM,OPO,OPP,OPM,OMO,OMP,OMM,POO,POP,POM,PPO,PMO,MOO,MOP,MOM,MPO,MMO

// Direction codes in order of increasing distance from current node as used in automated wall normal identification routine
//                Direction    0       1       2       4       8       16      32        5      10         6        9       20      40
//                           zyx     zyx     zyx      zyx     zyx      zyx     zyx      zyx     zyx       zyx      zyx      zyx     zyx
static int   nrmlD3Q19[19]={ OOO,    OOP,    OOM,     OPO,    OMO,     POO,    MOO,     OPP,    OMM,      OPM,     OMP,     PPO,    MMO,      PMO,     MPO,     POP,    MOM,      POM,     MOP};
static int dirD3Q19[19][3]={{0,0,0},{0,0,1},{0,0,-1},{0,1,0},{0,-1,0},{1,0,0},{-1,0,0},{0,1,1},{0,-1,-1},{0,1,-1},{0,-1,1},{1,1,0},{-1,-1,0},{1,-1,0},{-1,1,0},{1,0,1},{-1,0,-1},{1,0,-1},{-1,0,1}};
static int   indxD3Q19[19]={ 0,      1,      2,       3,      4,       5,      6,       7,      8,        9,       10,      11,     12,       13,      14,      15,     16,       17,      18};

__constant__ __device__ float vecD3Q19[19][3]={{0.,0.,0.},{0.,0.,1.},{0.,0.,-1.},{0.,1.,0.},{0.,-1.,0.},{1.,0.,0.},{-1.,0.,0.},{0.,S2,S2},{0.,-S2,-S2},{0.,S2,-S2},{0.,-S2,S2},{S2,S2,0.},{-S2,-S2,0.},{S2,-S2,0.},{-S2,S2,0.},{S2,0.,S2},{-S2,0.,-S2},{S2,0.,-S2},{-S2,0.,S2}};


// No slip boundary conditions
static int noslpD3Q19[19]={OOO,OOM,OOP,OMO,OMM,OMP,OPO,OPM,OPP,MOO,MOM,MOP,MMO,MPO,POO,POM,POP,PMO,PPO};
// x,y,z component reflection boundary conditions
static int xMrflD3Q19[19]={OOO,OOP,OOP,OPO,OPP,OPP,OMO,OMP,OMP,POO,POP,POP,PPO,PMO,MOO,MOP,MOP,MPO,MMO};
static int xPrflD3Q19[19]={OOO,OOM,OOM,OPO,OPM,OPM,OMO,OMM,OMM,POO,POM,POM,PPO,PMO,MOO,MOM,MOM,MPO,MMO};
static int yMrflD3Q19[19]={OOO,OOP,OOM,OPO,OPP,OPM,OPO,OPP,OPM,POO,POP,POM,PPO,PPO,MOO,MOP,MOM,MPO,MPO};
static int yPrflD3Q19[19]={OOO,OOP,OOM,OMO,OMP,OMM,OMO,OMP,OMM,POO,POP,POM,PMO,PPO,MOO,MOP,MOM,MMO,MMO};
static int zMrflD3Q19[19]={OOO,OOP,OOM,OPO,OPP,OPM,OMO,OMP,OMM,POO,POP,POM,PPO,PMO,POO,POP,POM,PPO,PMO};
static int zPrflD3Q19[19]={OOO,OOP,OOM,OPO,OPP,OPM,OMO,OMP,OMM,MOO,MOP,MOM,MPO,MMO,MOO,MOP,MOM,MPO,MMO};


// LBM model sound speed (in nondimensional lattice units)
//                     N/A,        D3Q19:1/sqrt(3)
static float csLBM[4]={0.577350269,   0.577350269, 0.577350269, 0.577350269};

#endif    // Close initial #ifndef
