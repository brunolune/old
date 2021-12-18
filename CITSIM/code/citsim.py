#!/usr/bin/python
# Makefile - Defines targets for mass spectrum computation stages
# citsim.py: Python script to simulate mass spectrum experiment

from pylab import *    # Import scientific computation environment
import os,sys,time     # Operating system utilities
#===============================================================================
# Experiments start at t=0 and progress through the following stages:
#  1. Ion cooling during time interval dtCool (ring voltage maintained constant)
#  2. Ion ejection during time interval dtEject (linear ring voltage increase)

# Ion trap geometry is specified in IonTrap.edp

# Fundamental units: length in mm, time in ms, voltage in V, mass in amu
# Derived units: frequency in kHz

# Experiment parameters
nCycles  = 1                      # Number of injection/cooling/ejection cycles

# Voltages applied to ion trap
ringRF   = 6000.0                 # Ring electrode RF frequency (kHz)
ringV    = 200.0                  # Ring electrode initial voltage amplitude (V)
ringdVdt = 30.0                   # Ring electrode voltage ramp rate (V/ms)
entV     = 0.0                    # Entry cap voltage
extV     = 1.0                    # Exit  cap voltage
tini     = 0.0                    # Experiment start time (ms)
dtCool   = 1.0                    # Trap cooling time (ms)
dtEject  = 3.0                    # Trap ejection time (ms)

# Initial ion conditions
nIons = 1                         # Number of ions initially injected into trap
Tion  = 500.0                     # Ion temperature (K)
# Radius of sphere containing initial positions (randomly distributed)
rSphere = 1.0e-2                  # (mm)
# Deform initial spherical placement into ellipsoid
axEll = array([1.,1.,1.])         # non-dimensional
# Ion isotope masses (amu)
mIsotope = array([124.,126.,128.,129.,130.,131.,132.,134.,136.])
# Ion isotope charges (in electron-units)
cIsotope = array([-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.,-1.])
# Ion isotope frequencies
pIsotope = array([0.09,0.09,1.92,26.44,4.08,21.18,26.89,10.44,8.87])

# Background conditions
pres     = 0.                     # Background pressure (Torr)
mBuffer  = 4.0                    # Buffer gas mass
                                 
#===============================================================================
WorkDir = os.environ['PWD']
print "Mass spectrum directory: ",WorkDir
os.chdir(WorkDir)
sys.path.insert(0, WorkDir)

# Read values from ion trap geometry definitions
import ITgeom
nDim=ITgeom.nDim; h=ITgeom.Sh; rExit=ITgeom.rExit
nX=ones(nDim,'int'); dX=zeros(nDim,'float'); Xmin=zeros(nDim,'float')
nX[0]  =ITgeom.Snx;   nX[1]  =ITgeom.Sny
dX[0]  =ITgeom.Sdx;   dX[1]  =ITgeom.Sdy
Xmin[0]=ITgeom.Sxmin; Xmin[1]=ITgeom.Symin
nGrdPts=nX[0]*nX[1]
if nDim==3:
  nX[2]=ITgeom.Snz; dX[2]=ITgeom.Sdz; Xmin[2]=ITgeom.Szmin; nGrdPts=nGrdPts*nX[2]

print "Reading trap geometry and electric field forces"
print "nDim =",nDim," nX =",nX

def dist(i,j,n):
  return sqrt( (n[i,0]-n[j,0])**2 + (n[i,1]-n[j,1])**2 )
  
# Read ion trap geometry mesh
IonTrapMeshFile = open("../IonTrapSh.msh",'r')
nline=0; NrNodes=0; NrElems=0; nNode=0; nElem=0;
for line in IonTrapMeshFile:
  line = line.strip()
  cols = line.split()
  if nline==0: # Initial line with number of nodes, elements
    NrNodes = int(cols[0]); NrElems = int(cols[1])
    ITnodes = zeros([NrNodes,nDim+1],'double')
    ITelems = zeros([NrElems,nDim+2],'int')
    ElemEdgeLen = zeros([NrElems,3],'double') # Element data: edge lengths (3), 
  if (0<nline) and (nline<=NrNodes): # Node coordinate, boundary condition line
    ITnodes[nNode] = map(float, cols)
    nNode = nNode + 1
  if (NrNodes<nline) and (nline<=NrNodes+NrElems): # Element node line
    nodes = map(int, cols); ITelems[nElem] = nodes
    ITelems[nElem,0:nDim+1] = ITelems[nElem,0:nDim+1]-1
    n = ITelems[nElem,0:nDim+1]
    ElemEdgeLen[nElem,0] = dist(n[0],n[1],ITnodes)
    ElemEdgeLen[nElem,1] = dist(n[1],n[2],ITnodes)
    ElemEdgeLen[nElem,2] = dist(n[2],n[0],ITnodes)
    nElem = nElem + 1
  nline = nline + 1
IonTrapMeshFile.close()

# Read forces produced by unit voltage on entry, ring, exit electrodes
ient=0; irng=1; iext=2;                  # Electrode indices
nF = zeros(nDim+2,'int'); nF[0:nDim] = nX[0:nDim]; nF[nDim]=nDim; nF[nDim+1]=3
Fgrid=zeros(nF,dtype='double',order='F') # Grid values of electric forces
Slin=zeros(nGrdPts,'double')             # linear ordering of a scalar field

def readF(fname,nX,Sl):
  SolFile = open(fname,"r")
  text = SolFile.read(); SolFile.close()
  text = text.strip(); fields = text.split()
  Sl[:] = map(float, fields[1:])  

# 0. On entry cap electrode
k=ient
readF("../UnitEntryCapFx.sol",nX,Slin); Fgrid[:,:,0,k] = Slin.reshape(nX,order='F')
readF("../UnitEntryCapFy.sol",nX,Slin); Fgrid[:,:,1,k] = Slin.reshape(nX,order='F')
if nDim==3:
  readF("../UnitEntryCapFz.sol",nX,Slin); Fgrid[:,:,2,k] = Slin.reshape(nX,order='F')

# 1. On ring electrode
k=irng
readF("../UnitRingFx.sol",nX,Slin); Fgrid[:,:,0,k] = Slin.reshape(nX,order='F')
readF("../UnitRingFy.sol",nX,Slin); Fgrid[:,:,1,k] = Slin.reshape(nX,order='F')
if nDim==3:
  readF("../UnitRingFz.sol",nX,Slin); Fgrid[:,:,2,k] = Slin.reshape(nX,order='F')

# 2. On exit electrode
k=iext
readF("../UnitExitCapFx.sol",nX,Slin); Fgrid[:,:,0,k] = Slin.reshape(nX,order='F')
readF("../UnitExitCapFy.sol",nX,Slin); Fgrid[:,:,1,k] = Slin.reshape(nX,order='F')
if nDim==3:
  readF("../UnitExitCapFz.sol",nX,Slin); Fgrid[:,:,2,k] = Slin.reshape(nX,order='F')

# Scale forces
Fgrid = 9.64853e1*Fgrid

# Function to switch from trap coordinates to ion trajectory computation coordinates
def Xtraj(Xtrap,h,Xmin):
  return (Xtrap-Xmin)/h

# Function to switch from ion trajectory computation coordinates to trap coordinates
def Xtrap(Xtraj,h,Xmin):
  return h*Xtraj+Xmin

# Load ion trajectory computational kernels
from ion_traj import *

# Normalize isotope distribution


# Carry out experiment
nplt  = 500                                     # Number of time slices at which to return ion positions
trmp  = tini + dtCool                           # Time at which rampup of ring electrode voltage starts
tfin  = trmp + dtEject                          # Time at which trap is cleared (end of voltage rampup)
dt    = 1.0e-3                                  # Verlet integrator time step (ms)
Omega = 2.*pi*ringRF                            # RF pulsation
Xplt  = zeros([nIons,nDim,nplt],'double',order='F')   # Trajectories
Vplt  = zeros([nIons,nDim,nplt],'double',order='F')   # Velocities
# Load parameters into a single array for cleaner call interface
tfin = 8.e-1; dt=2.0e-3
params = zeros(32,'double',order='F')
params[0]=tini; params[1]=trmp; params[2]=tfin; params[3]=dt
params[4]=ringV; params[5]=ringdVdt; params[6]=entV; params[7]=extV
params[8]=Omega
for kcyc in range(nCycles):
  print 'Experiment cycle: ',kcyc
  # Define initial ion positions  
  r=rand(nIons)*rSphere; phi=rand(nIons)*2.*pi
  Xplt[:,0,0] = r[:]*cos(phi[:]); Xplt[:,1,0] = r[:]*sin(phi[:])
  for l in range(nIons):
    Xplt[l,:,0] = Xtraj(Xplt[l,:,0],h,Xmin)     # Switch to ion trajectory computation coordinates (grid based)
  # Define initial ion positions
  
  phi=rand(nIons)*2.*pi
  # Carry out time integration
  trajectory2d(params,Xplt,Vplt,Fgrid)
  figure(); np=int(params[9])
  for n in range(nIons):
    plot(Xplt[n,0,0:np],Xplt[n,1,0:np],'o')
  show()





