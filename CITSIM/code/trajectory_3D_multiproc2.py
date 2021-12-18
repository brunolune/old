#!/usr/bin/python
# Makefile - Defines targets for mass spectrum computation stages
# citsim.py: Python script to simulate mass spectrum experiment

from pylab import *    # Import scientific computation environment
import matplotlib.animation as animation
from random import *
import os,sys,time     # Operating system utilities
WorkDir = os.environ['PWD']
print "Mass spectrum directory: ",WorkDir
os.chdir(WorkDir)
sys.path.insert(0, WorkDir)

# Physical constants
eC = 1.602176565e-19                               # electron charge (C)
amu = 1.66053886e-27                               # atomic mass unit (kg)
kB = 1.3806503e-23                                 # Boltzmann constant (J/K)
eamu = eC/amu                                      # C/amu
kBamu = kB/amu                                     # mm^2/ms^2/amu/K
ke = 8.9875517873681764e9                          # ke=1/(4*pi*eps0) (N/m^2/C^2)
fCoul = ke*eC*eC/amu                               # Coulomb force factor
Torr = 101325./760.                                # Torr in N/m^2
eps0 = 8.854187817620e-12                          # Vacuum permittivity (F/m)
msec = 1.0e-3                                      # milisecond in s
# Langevin collision probability factor
fLanCollP = eC/(2*eps0)*sqrt(4.*pi*eps0/(1.e30*amu))*(Torr/kB)*msec

import trajectory_params as tp
# Make parameters from trajectory_params available to this Python namespace
nCycles  = tp.nCycles

ringV    = tp.ringV
ringU    = tp.ringU          # constant DC voltage on Ring electrode
ringRF   = tp.ringRF
ringdVdt = tp.ringdVdt
tini     = tp.tini
dtCool   = tp.dtCool
dtEject  = tp.dtEject

entV     = tp.entV
entRF    = tp.entRF
entVRF   = tp.entVRF
DelPhiEnt= tp.DelPhiEnt

extV     = tp.extV
extRF    = tp.extRF
extVRF   = tp.extVRF
DelPhiExt= tp.DelPhiExt

Tion  = tp.Tion
nIons = tp.nIons
RSphere = tp.RSphere
rSphere = tp.rSphere
xcSphere = tp.xcSphere
tBirth = tp.tBirth
fixedTOB = tp.fixedTOB
Xinitconds = tp.Xinitconds
Vinitconds = tp.Vinitconds
fixedInitPos = tp.fixedInitPos
fixedInitVel = tp.fixedInitVel
axEll = tp.axEll
mIsotope = tp.mIsotope
cIsotope = tp.cIsotope
pIsotope = tp.pIsotope

tfact  = tp.tfact
ionion = tp.ionion
ionionGPU = tp.ionionGPU

presGas  = tp.presGas
tempGas  = tp.tempGas
mGas     = tp.mGas
CollisionParameters  = tp.CollisionParameters
CollisionModel = tp.CollisionModel

ComputeSpectrum=tp.ComputeSpectrum
dtHistBin = tp.dtHistBin
nSigBin = tp.nSigBin
rDetected = tp.rDetected
zDetected = tp.zDetected
rEndcapHole = tp.rEndcapHole
zEndcapHole = tp.zEndcapHole
DetectedIon = tp.DetectedIon
tMassScale = tp.tMassScale

nplt = tp.nplt
kx=tp.kx; ky=tp.ky                      # Skip ratios for quiver
TrajectoryPlot=tp.TrajectoryPlot
ElectrodePlots = tp.ElectrodePlots
AnimateIons = tp.AnimateIons
nCoolFrames = tp.nCoolFrames
nEjectFrames = tp.nEjectFrames

nDimFly=tp.nDimFly #add nDimFly to specify the dimension of the flight (bruno)

# End of trajectory_params specifications

#kTm    = kBamu*Tion;    (moved in fly() for being able to reset pIsotope... with loops.py prog) (bruno)
#nIsotopes = mIsotope.shape[0]
#CDFIsotope = zeros(mIsotope.shape)
#CDFIsotope[0] = pIsotope[0]
#for l in range(1,nIsotopes):
#  CDFIsotope[l] = CDFIsotope[l-1] + pIsotope[l]

# Read values from ion trap geometry definitions
import ITgeom
nDimField=ITgeom.nDim; h=ITgeom.Sh; rExit=ITgeom.rExit
nX=ones(nDimFly,'int'); dX=zeros(nDimFly,'float'); Xmin=zeros(3,'float') #use nDimFly instead of nDim (bruno)
nX[0]  =ITgeom.Snx;   nX[1]  =ITgeom.Sny
dX[0]  =ITgeom.Sdx;   dX[1]  =ITgeom.Sdy
Xmin[0]=ITgeom.Sxmin; Xmin[1]=ITgeom.Symin
nGrdPts=nX[0]*nX[1]
if nDimFly==3:
  nX[2]=ITgeom.Sny; dX[2]=ITgeom.Sdy; Xmin[2]=ITgeom.Symin; #removed nGrdPts change since we want to use field in 2D (bruno)
if nDimField==3:
   nGrdPts=nGrdPts*nX[2]
nB = ITgeom.nBorder
zB = zeros(nB); rB = zeros(nB); bcB = zeros(nB).astype('int')
for i in range(nB):
  zB[i] = ITgeom.zBorder[i]; rB[i] = ITgeom.rBorder[i]; bcB[i] = ITgeom.bcBorder[i]

flowBC = zeros([nB,nDimFly+3]); #changed to nDimFly to avoid error when calling trajectory3D (bruno)
nBC = size(presGas)   # Length of presGas vectors specifies number of fluid background/boundary conditions
iBC = 1               # iBC=0 is background, iBC = 1 is first boundary (exit)
if nDimField==2:
  for i in range(nB):
    flowBC[i,0] = zB[i]; flowBC[i,1] = rB[i]; flowBC[i,2] = bcB[i]
    if bcB[i] == 1:    # Load fluid conditions
      flowBC[i,3] = presGas[iBC]; flowBC[i,4] = tempGas[iBC]
      iBC = iBC + 1

# Finished reading ion trap geometry

print "Reading trap geometry and electric field forces"
print "nDimField =",nDimField," nX =",nX

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
    ITnodes = zeros([NrNodes,nDimField+1],'double')
    ITelems = zeros([NrElems,nDimField+2],'int')
    ElemEdgeLen = zeros([NrElems,3],'double') # Element data: edge lengths (3), 
  if (0<nline) and (nline<=NrNodes): # Node coordinate, boundary condition line
    ITnodes[nNode] = map(float, cols)
    nNode = nNode + 1
  if (NrNodes<nline) and (nline<=NrNodes+NrElems): # Element node line
    nodes = map(int, cols); ITelems[nElem] = nodes
    ITelems[nElem,0:nDimField+1] = ITelems[nElem,0:nDimField+1]-1
    n = ITelems[nElem,0:nDimField+1]
    ElemEdgeLen[nElem,0] = dist(n[0],n[1],ITnodes)
    ElemEdgeLen[nElem,1] = dist(n[1],n[2],ITnodes)
    ElemEdgeLen[nElem,2] = dist(n[2],n[0],ITnodes)
    nElem = nElem + 1
  nline = nline + 1
IonTrapMeshFile.close()
zc=arange(nX[0])*h+Xmin[0]; rc=arange(nX[1])*h+Xmin[1]; #no change here should work since nX,Xmin dimension nDimFly=3 (bruno)
zqv=arange(Xmin[0],Xmin[0]+nX[0]*h,kx*h); rqv=arange(Xmin[1],Xmin[1]+nX[1]*h,ky*h);

# Read forces produced by unit voltage on entry, ring, exit electrodes
ient=0; irng=1; iext=2;                  # Electrode indices
nF = zeros(nDimField+2,'int'); nF[0:nDimField] = nX[0:nDimField]; nF[nDimField]=nDimField; nF[nDimField+1]=3 #no change here should work since nDimField=2 (bruno)
nV = zeros(nDimField+1,'int'); nV[0:nDimField] = nX[0:nDimField]; nV[nDimField]=3;
aGrid=zeros(nF,dtype='double',order='F') # Grid values of acceleration induced on 1 amu ion by electrode fields
VGrid=zeros(nV,dtype='double',order='F') # Grid values of electrostatic potentials
Slin=zeros(nGrdPts,'double')             # linear ordering of a scalar field

def readF(fname,Sl):
  SolFile = open(fname,"r")
  text = SolFile.read(); SolFile.close()
  text = text.strip(); fields = text.split()
  Sl[:] = map(float, fields[1:])  

# 0. On entry cap electrode
k=ient
readF("../UnitEntryCapPotential.sol",Slin); VGrid[:,:,k] = Slin.reshape(nX[0:2],order='F') #changed nX to nX[0:2] for use fields in 2d (bruno)
readF("../UnitEntryCapFx.sol",Slin); aGrid[:,:,0,k] = Slin.reshape(nX[0:2],order='F')
readF("../UnitEntryCapFy.sol",Slin); aGrid[:,:,1,k] = Slin.reshape(nX[0:2],order='F')
if nDimField==3:
  readF("../UnitEntryCapFz.sol",Slin); aGrid[:,:,2,k] = Slin.reshape(nX,order='F')
  
if ElectrodePlots: # ElectrodePlots should be set to 0 (bruno)
  figure(1); subplot(2,2,2)  
  contour(zc,rc,VGrid[:,:,k].T,25)
  hold(True); axis('equal') 
  title('Electric potential and force from entry cap electrode'); xlabel('z (mm)'); ylabel(' r (mm)')
  (nqx,nqy)=shape(aGrid[::kx,::ky,0,k])
  quiver(zqv,rqv,aGrid[::kx,::ky,0,k].T,aGrid[::kx,::ky,1,k].T)
  
# 1. On ring electrode
k=irng
readF("../UnitRingPotential.sol",Slin); VGrid[:,:,k] = Slin.reshape(nX[0:2],order='F')
readF("../UnitRingFx.sol",Slin); aGrid[:,:,0,k] = Slin.reshape(nX[0:2],order='F')
readF("../UnitRingFy.sol",Slin); aGrid[:,:,1,k] = Slin.reshape(nX[0:2],order='F')
if nDimField==3:
  readF("../UnitRingFz.sol",Slin); aGrid[:,:,2,k] = Slin.reshape(nX,order='F')
if ElectrodePlots: # ElectrodePlots should be set to 0 (bruno)
  figure(1); subplot(2,2,3)
  contour(zc,rc,VGrid[:,:,k].T,25)
  hold(True); axis('equal') 
  title('Electric potential and force from ring electrode'); xlabel('z (mm)'); ylabel(' r (mm)')
  quiver(zqv,rqv,aGrid[::kx,::ky,0,k].T,aGrid[::kx,::ky,1,k].T)

# 2. On exit electrode
k=iext
readF("../UnitExitCapPotential.sol",Slin); VGrid[:,:,k] = Slin.reshape(nX[0:2],order='F')
readF("../UnitExitCapFx.sol",Slin); aGrid[:,:,0,k] = Slin.reshape(nX[0:2],order='F')
readF("../UnitExitCapFy.sol",Slin); aGrid[:,:,1,k] = Slin.reshape(nX[0:2],order='F')
if nDimField==3:
  readF("../UnitExitCapFz.sol",Slin); aGrid[:,:,2,k] = Slin.reshape(nX,order='F')
if ElectrodePlots: # ElectrodePlots should be set to 0 (bruno)
  figure(1); subplot(2,2,4)
  contour(zc,rc,VGrid[:,:,k].T,25)
  hold(True); axis('equal') 
  title('Electric potential and force from exit cap electrode'); xlabel('z (mm)'); ylabel(' r (mm)')
  quiver(zqv,rqv,aGrid[::kx,::ky,0,k].T,aGrid[::kx,::ky,1,k].T) 

# Scale forces
aGrid = -eamu*aGrid   # aGrid now contains acceleration induced by unit voltage on an electrode on
                      # a point mass of one atomic mass unit and one positive charge unit.
                      # Acceleration is in CITSIM units (mm)/(ms)^2 = (m/s^2) * 1.0e3

# Load ion trajectory computational kernels
if nDimFly==2:
  from ion_traj_2D import *
elif nDimFly==3:
  from ion_traj_3D_multiproc import *

#set nDim=2 to nDimFly=3 from here (bruno)
#from here nDim is equivalent to nDimFly
nDim=nDimFly

# Define experiment *****************************************************************************************************************************************************
def fly():
  global aGrid,mIon,cIon,NrIso,nXplt,tplt,Xplt,Vplt,Fplt,Ncolplt,ntEject,factAnim,nIso,IonExit,dtHistBin,\
         tExit,params,AvgFeplt,AvgFcolplt,jobs,p,nIonstarts,nIonstops                     # all variables that are declared in fly() that we want to be able to access in ipython 
  
  kTm    = kBamu*Tion; #moved in fly() for being able to reset pIsotope... with loops.py prog (bruno)
  nIsotopes = mIsotope.shape[0] 
  CDFIsotope = zeros(mIsotope.shape)
  CDFIsotope[0] = pIsotope[0]
  for l in range(1,nIsotopes):
    CDFIsotope[l] = CDFIsotope[l-1] + pIsotope[l]
  
  # Compute time needed by an ion to traverse one grid spacing (move it here to take into account parameters change in loopsx.py, bruno) 
  mmin = min(mIsotope)
  maxV = max(abs(ringV),abs(entV),abs(extV))
  vel = sqrt(2.*maxV/mmin*eamu)
  dtgrd = h/vel
  
  trmp  = tini + dtCool                           # Time at which rampup of ring electrode voltage starts
  tfin  = trmp + dtEject                          # Time at which trap is cleared (end of voltage rampup)
  if ringRF>0:                                    # Verlet integrator time step (ms)
    dt = 1./ringRF/tfact
  else:
    dt = dtgrd
  Omega    = 2.*pi*ringRF                            # RF pulsations
  OmegaEnt = 2.*pi*entRF
  OmegaExt = 2.*pi*extRF
  NrIso  = zeros(nIsotopes,'int')
  nIso   = zeros(nIons,'int');   nI = zeros(nIons,'int')
  mIon  = zeros(nIons,'double'); mz = zeros(nIons,'double')
  cIon  = zeros(nIons,'double'); cz = zeros(nIons,'double')
  IonExit = zeros([nIons,4*nDim+1],'double',order='F')   # Trap exit time, location, velocity  
  tplt   = zeros(nplt,'double')                          # Times
  Xplt   = zeros([nIons,nDim,nplt],'double',order='F')   # Trajectories
  Vplt   = zeros([nIons,nDim,nplt],'double',order='F')   # Velocities
  Fplt   = zeros([nIons,nDim,nplt],'double',order='F')   # Forces
  AvgFeplt  = zeros([nIons,nplt],'double',order='F')        # average electric Force
  AvgFcolplt  = zeros([nIons,nplt],'double',order='F')      # average collision Force
  Xcol   = zeros([nIons,nDim,nplt],'double',order='F')   # Ion position before collision
  Ncolplt = zeros([nIons,nplt],'double',order='F')    # Nb of collision along the traj
  Vcol   = zeros([nIons,nDim,nplt],'double',order='F')   # Ion velocity before collision 
  Fcol   = zeros([nIons,nDim,nplt],'double',order='F')   # Collision force
  nXplt  = zeros(nIons,'int32')                          # Number of plots ion was still in trap
  Xh = zeros(nDim); Vh = zeros(nDim)
  
  # Load parameters into a single array for cleaner call interface
  params = zeros(512,'double',order='F') #to be used in fortran program ion_traj_xD.f90
  params[0]=tini; params[1]=trmp; params[2]=tfin; params[3]=dt
  params[4]=ringV; params[5]=ringdVdt; params[6]=entV; params[7]=extV
  params[8]=Omega; params[9]=OmegaEnt; params[10]=OmegaExt
  params[11]=DelPhiEnt; params[12]=DelPhiExt
  params[13]=entVRF; params[14]=extVRF
  params[15]=ionion; params[16]=ionionGPU  
  params[17]=fCoul
  params[18]=mGas
  params[19]= max(presGas) - min(presGas)
  params[20]= CollisionModel
  params[21]= h
  params[22:25]= Xmin[0:3]
  alphae = CollisionParameters[0]  
  params[26]=fLanCollP*(presGas[0]/tempGas[0])*dt*sqrt(alphae)
  params[27]=kB*tempGas[0]/(mGas*amu)
  params[28]=sqrt(params[27])
  diamHS = CollisionParameters[3]
  # params[29] contains dt*sigma*n=dt*sigma*p/(k*T) constant part of collision factor calculation within HS1 collision model
  params[29]=dt*CollisionParameters[2]*presGas[0]*Torr/(kB*tempGas[0]) / 1.e9 * h 
  print 'cs=',CollisionParameters[2],'cor=',CollisionParameters[1]
  if nDim==2:
    params[30]=sqrt(pi * kBamu*tempGas[0]/(2*mGas))/h # Average (mean) buffer gas speed in 2D (bruno)
    params[31]=sqrt(2 * kBamu*tempGas[0]/mGas)/h    # Median buffer gas speed in 2D (bruno)
  elif nDim==3:
    params[30]=sqrt(8/pi * kBamu*tempGas[0]/mGas)/h # Average (mean) buffer gas speed in 3D
    params[31]=sqrt(3 * kBamu*tempGas[0]/mGas)/h    # Median buffer gas speed in 3D (bruno)
  params[32]=ringU
  params[33]=nDimField
  params[34]=nDimFly
  # Specify rectangle defining trap interior
  params[50]=rEndcapHole; params[51]=zEndcapHole
  params[52]=tBirth # tBirth determine if ions created randomly within 1/ringRF (0) or every 1/ringRF/nIons (1) (bruno)
  params[53]=fixedTOB
  params[54]=fixedInitPos
  params[55]=fixedInitVel
  
  nGasBC = size(presGas)
  if nGasBC>9:
    print 'Maximum of 9 gas boundary conditions may be specified in trajectory_params'
    print 'presGas=',presGas
    print '9<nGasBC=',nGasBC
    sys.exit()
  params[100] = nGasBC
  params[101:101+nGasBC] = presGas[0:nGasBC]
  params[111:111+nGasBC] = tempGas[0:nGasBC]
  params[120:130] = CollisionParameters[0:10]
  params[148]=0
  
  #fig = figure(0); clf();
  #figWM = get_current_fig_manager(); figWM.window.SetSize((1000,800))

  #seed(0)                               # Set seed to ensure reproducible results
  
  #Main loop
  for kcyc in range(nCycles):
    print 'Experiment cycle: ',kcyc
    
    # Define initial ion positions    
    if fixedInitPos:
    #fixed init conds (bruno)
      Xplt[:,0,0]=Xinitconds[0] #fixed init conds (bruno)
      Xplt[:,1,0]=Xinitconds[1]
      if nDim==3:
        Xplt[:,2,0]=Xinitconds[2]
    else:
      # random initial positions
      r = zeros(nIons); phi = zeros(nIons); tht = zeros(nIons)
      dr = RSphere-rSphere
      for l in range(nIons):
        r[l]=random()*dr+rSphere; phi[l]=(random())*pi; tht[l]=random()*2.*pi #changed phi between [0,pi] instead of [-pi/2,pi/2] (bruno)
        if nDim==2:
          Xplt[:,0,0] = axEll[0]*r[:]*cos(tht[:]) + xcSphere[0]
          Xplt[:,1,0] = axEll[1]*r[:]*sin(tht[:]) + xcSphere[1]
        else:
          Xplt[:,0,0] = axEll[0]*r[:]*cos(tht[:])*sin(phi[:]) + xcSphere[0]
          Xplt[:,1,0] = axEll[1]*r[:]*sin(tht[:])*sin(phi[:]) + xcSphere[1]
          Xplt[:,2,0] = axEll[2]*r[:]*cos(phi[:]) + xcSphere[2]

    print 'Initial ion distribution:'
    if nDim==2:
      print ' Ion    m/Z   z(mm)   r(mm)  Vz(mm/ms)  Vr(mm/ms)'
    else:
      print ' Ion    m/Z       x       y       z        Vx       Vy       Vz'
    
    for l in range(nIons):
      r = random(); n = 0;
      while (r>CDFIsotope[n]) & (n<nIsotopes-1):
        n = n + 1
      nI[l] = n; NrIso[n] = NrIso[n] + 1
      mz[l] = mIsotope[n]; cz[l] = cIsotope[n]

    # Sort ions in order of increasing mass and define initial velocities
    perm = argsort(mz)
    for l in range(nIons):
      mIon[l] = mz[perm[l]]; cIon[l] = cz[perm[l]]; nIso[l] = nI[perm[l]]
      # Define initial ion velocities
      if fixedInitVel:
        Vplt[:,0,0] = Vinitconds[0] #fixed init conds (bruno)
        Vplt[:,1,0] = Vinitconds[1]
        if nDim==3:
          Vplt[:,2,0]=Vinitconds[2]
      else:
        #random initial velocities
        sigma = sqrt(kTm/mIon[l])
        for d in range(nDim):
          Vplt[l,d,0] = gauss(0.,sigma)
          
      #print initial conditions
      if nDim==2:
        print '%3d    %4d %7.4f %7.4f  %7.1f  %7.1f' % \
              (l,mIon[l],Xplt[l,0,0],Xplt[l,1,0],Vplt[l,0,0],Vplt[l,1,0])                             
      else:
        print '%3d    %4d %7.3f %7.3f %7.3f  %7.3f  %7.3f  %7.3f' % \
              (l,mIon[l],Xplt[l,0,0],Xplt[l,1,0],Xplt[l,2,0],Vplt[l,0,0],Vplt[l,1,0],Vplt[l,2,0])
    
    nIonstarts=np.arange(1,nIons+1,nIons/nbproc).tolist()
    nIonstops=np.append((np.delete(np.arange(1,nIons+1,nIons/nbproc),0)-1),nIons).tolist()
    
    import multiprocessing as mp
    import logging
    
    # Carry out time integration
    print 'nDim b4 calling ion_traj',nDim
    cput0=time.clock()
    if nDim==2:
        trajectory2d(params,mIon,cIon,IonExit,nXplt,tplt,Xplt,Vplt,Fplt,Ncolplt,Xcol,Vcol,Fcol,aGrid,flowBC)
    elif nDim==3:
        if __name__ == '__main__':
            multiprocessing.log_to_stderr(logging.DEBUG)
            nIonsp=zeros(nbproc)
           
            jobs = []
            for i in range(nbproc):
                nIonsp[i]=nIonstops[i]-nIonstarts[i]+1
                IonExitp = zeros([nIonsp[i],4*nDim+1],'double',order='F')
                Xpltp   = zeros([nIonsp[i],nDim,nplt],'double',order='F')
                Vpltp   = zeros([nbproc,nIonsp[i],nDim,nplt],'double',order='F')
                Fpltp   = zeros([nbproc,nIonsp[i],nDim,nplt],'double',order='F')
                AvgFepltp  = zeros([nbproc,nIonsp[i],nplt],'double',order='F')
                AvgFcolpltp  = zeros([nbproc,nIonsp[i],nplt],'double',order='F')
                Xcolp   = zeros([nbproc,nIonsp[i],nDim,nplt],'double',order='F')
                Ncolpltp = zeros([nbproc,nIonsp[i],nplt],'double',order='F')
                Vcolp   = zeros([nbproc,nIonsp[i],nDim,nplt],'double',order='F')
                Fcolp   = zeros([nbproc,nIonsp[i],nDim,nplt],'double',order='F')
                nXpltp  = zeros([nbproc,nIonsp[i]],'int32')
                
                p = multiprocessing.Process (target=trajectory3d,args=(params,mIon,cIon,IonExit[i,:,:],\
                nXplt[i,:,:],tplt[nIonstarts[i,:,:],Xplt[i,:,:],\
                Vplt[i,:,:],Fplt[i,:,:],Ncolplt[i,:,:],\
                Xcolp[:,:,:],Vcolp[:,:,:,i],Fcolp[:,:,:,i],AvgFepltp[:,:,:,i],AvgFcolpltp[:,:,:,i],aGrid,flowBC))
                
                jobs.append(p)
                p.start()
                
            for j in jobs:
                j.join()
                
            for i in range(nbproc):
                IonExit[nIonstarts[i]:nIonstops[i],:,:] = IonExitp[i,:,:,:,]  
                Xplt[nIonstarts[i]:nIonstops[i],:,:] = Xpltp[i,:,:,:,]
                Vplt[nIonstarts[i]:nIonstops[i],:,:] = Vpltp[i,:,:,:,]
                Fplt[nIonstarts[i]:nIonstops[i],:,:] = Fpltp[i,:,:,:,]
                AvgFeplt[nIonstarts[i]:nIonstops[i],:,:] = AvgFepltp[i,:,:,:,]
                AvgFcolplt[nIonstarts[i]:nIonstops[i],:,:] = AvgFcolpltp[i,:,:,:,]
                Xcol[nIonstarts[i]:nIonstops[i],:,:] = Xcolp[i,:,:,:,]
                Ncolplt[nIonstarts[i]:nIonstops[i],:,:] = Ncolpltp[i,:,:,:,]
                Vcol[nIonstarts[i]:nIonstops[i],:,:] = Vcolp[i,:,:,:,]
                Fcol[nIonstarts[i]:nIonstops[i],:,:] = Fcolp[i,:,:,:,]
                nXplt[nIonstarts[i]:nIonstops[i],:,:] = nXpltp[i,:,:,:,]
        
    cput1=time.clock()
    print 'Trajectory computation time = ',cput1-cput0,' secs'
    
    ntEject = 0
    for nt in range(max(nXplt)):
      if tplt[nt] < tini + dtCool:
        ntEject = nt
    factAnim = max(1,(nplt-ntEject)/nEjectFrames)
    nAnimFrames = nCoolFrames + nEjectFrames

    #print main parameters of the simulation (bruno)
    print 'ringRF=',ringRF,', ringV=',ringV,', dtEject=',dtEject,', presGas[0]=', presGas[0],', entVRF=',entVRF,', entRF=',entRF	
    
    # Plot results
    # Trajectories plot
    if TrajectoryPlot:
      if nDimFly==2:
        figure(0); clf()
        figWM = get_current_fig_manager(); figWM.window.SetSize((800,600))
        clrs = ['blue','green','red','magenta','skyblue','pink','lightgreen','orange','darkblue','darkgreen','darkred','yellow']      
        for k in range(nIsotopes):
          plot([], [], color=clrs[k], marker='o', ms=6, label='M=' + str(mIsotope[k]) + 'amu')
        legend()      
        for l in range(nIons):
          nbpts = abs(nXplt[l])      #changed np to nbpts to avoid confusion with numpy abbreviation (bruno)
          if (nbpts>0):
            plot(Xplt[l,0,0:nbpts],Xplt[l,1,0:nbpts],'-',color=clrs[l%12]) #color=clrs[nIso[l]]
            plot(Xplt[l,0,nbpts-1:nbpts],Xplt[l,1,nbpts-1:nbpts],color=clrs[l%12],marker='o')
          
        hold(True); axis('equal')
        Vsum = ringV*VGrid[:,:,irng].T + entV*VGrid[:,:,ient].T + extV*VGrid[:,:,iext].T
        CPV1=contour(zc, rc,Vsum,25,colors='gray',linestyles='dashed')
        CPV2=contour(zc,-rc,Vsum,25,colors='gray',linestyles='dashed')
        clabel(CPV1, CPV1.levels, inline=True, fmt = '%g', fontsize=10)
        clabel(CPV2, CPV2.levels, inline=True, fmt = '%g', fontsize=10)
        zSeg=zeros(2); rSeg=zeros(2)
        linsty = [':r',':g','-k']
        maxz = max(zB); maxr = max(rB)
        for i in range(nB-1):
          zSeg[0]=zB[i]; zSeg[1]=zB[i+1]
          rSeg[0]=rB[i]; rSeg[1]=rB[i+1]
          plot(zSeg,rSeg,linsty[bcB[i]]); plot(zSeg,-rSeg,linsty[bcB[i]])  
        xlabel('x (mm)'); ylabel('r (mm)'); title('Ion trajectories')
        draw(); pause(0.01)
      elif nDimFly==3:
        import matplotlib as mpl
        from mpl_toolkits.mplot3d import Axes3D
        import matplotlib.pyplot as plt
        mpl.rcParams['legend.fontsize'] = 10
        clrs = ['blue','green','red','magenta','skyblue','pink','lightgreen','orange','darkblue','darkgreen','darkred','yellow']
        fig = plt.figure()
        ax = fig.gca(projection='3d')
        for l in range(nIons):
          plt.hold(True)
          ax.plot(Xplt[l,0,0:nXplt[l]],Xplt[l,1,0:nXplt[l]],Xplt[l,2,0:nXplt[l]],color=clrs[l%12],label='parametric curve')
          ax.set_xlabel('x (mm)')
          ax.set_ylabel('y (mm)')
          ax.set_zlabel('z (mm)')
          plt.show()
        #Cylinders
        #Ring
        xr=np.linspace(-0.390, 0.390, 100)
        zr=np.linspace(-0.5, 0.5, 100)
        Xring, Zring=np.meshgrid(xr, zr)
        Yring = np.sqrt(0.25-Zring**2)
        #left endcap interior hole
        xlecih=np.linspace(-0.640,-0.650, 2)
        zlecih=np.linspace(-0.15, 0.15, 30)
        Xlecih, Zlecih=np.meshgrid(xlecih, zlecih)
        Ylecih = np.sqrt(0.025-Zlecih**2)
        #left endcap exterior hole
        xleceh=np.linspace(-0.790, -0.890, 3)
        zleceh=np.linspace(-0.5, 0.5, 50)
        Xleceh, Zleceh=np.meshgrid(xleceh, zleceh)
        Yleceh = np.sqrt(0.25-Zleceh**2)
        #right endcap interior hole
        xrecih=np.linspace(0.640,0.650, 2)
        zrecih=np.linspace(-0.15, 0.15, 30)
        Xrecih, Zrecih=np.meshgrid(xrecih, zrecih)
        Yrecih = np.sqrt(0.025-Zrecih**2)
        #right endcap exterior hole
        xreceh=np.linspace(0.790, 0.890, 3)
        zreceh=np.linspace(-0.5, 0.5, 50)
        Xreceh, Zreceh=np.meshgrid(xreceh, zreceh)
        Yreceh = np.sqrt(0.25-Zreceh**2)
        # Draw parameters
        rstride = 15
        cstride = 11
        ax.plot_surface(Xring, Yring, Zring, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xring, -Yring, Zring, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xlecih, Ylecih, Zlecih, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xlecih, -Ylecih, Zlecih, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xleceh, Yleceh, Zleceh, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xleceh, -Yleceh, Zleceh, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xrecih, Yrecih, Zrecih, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xrecih, -Yrecih, Zrecih, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xreceh, Yreceh, Zreceh, alpha=0.2, rstride=rstride, cstride=cstride)
        ax.plot_surface(Xreceh, -Yreceh, Zreceh, alpha=0.2, rstride=rstride, cstride=cstride)
        plt.show
        ax.scatter(IonExit[:,1],IonExit[:,2],IonExit[:,3], c='red', alpha=0.5)
        plt.savefig('results/date_' + Date+'_data_'+Data+'_fulltrajs.png')
        
    # Exit times
    if ComputeSpectrum:
      figure(1); clf()
      figWM = get_current_fig_manager(); figWM.window.SetSize((800,600))
      # Form list of ion exit times through exit cap
      tExit = []
      for l in range(nIons):
        if DetectedIon(IonExit[l,0],IonExit[l,1],IonExit[l,2],IonExit[l,3],IonExit[l,4],IonExit[l,5],IonExit[l,6]):
          tExit.append(IonExit[l,0])
      if len(tExit)>5:
       tExit_mean = mean(tExit); tExit_std = std(tExit); nDetected = len(tExit)
       print nDetected,' detected ions out of ',nIons,' introduced in trap ',(100.*nDetected)/nIons,' %'
       print 'Exit time statistics: mean = ',tExit_mean,'(ms) standard deviation = ',tExit_std,' (ms)'
       if dtHistBin<=0:
         dtHistBin = max(1.0e-3,tExit_std/nSigBin)
       nHistBin = (max(tExit)-min(tExit))/dtHistBin
       MassData=[]
       if nHistBin==0: #to avoid error when no ions have been detected (bruno)
         nHistbin=1
       for t in tExit:
         MassData.append(tMassScale(t))    
       hist(MassData,nHistBin,facecolor='green')
       if tMassScale(1)==1:    
         xlabel('t (ms)')
       else:
         xlabel('M (amu)')
       ylabel('n (ions)'); title('Ion detection histogram')
    
    # Electrode plots
    if ElectrodePlots:
      figure(3); subplot(2,2,1)
      for k in range(nIons):
        nbpts = nXplt[k]  #changed np to nbpts to avoid confusion with numpy abbreviation (bruno)
        if (nump>2):
          plot(Xplt[k,0,0:nbpts],Xplt[k,1,0:nbpts],'-',color=clrs[nIso[k]])
          plot(Xplt[k,0,nbpts-1:nbpts],Xplt[k,1,nbpts-1:nbpts],marker='o',color= clrs[nIso[k]])
      hold(True); axis('equal'); contour(zc,rc,VGrid[:,:,irng].T,25,colors='gray',linewidth=0.5);
      xlabel('z (mm)'); ylabel('r (mm)'); title('Ion trajectories')    
      for i in range(nB-1):
        zSeg[0]=zB[i]; zSeg[1]=zB[i+1]
        rSeg[0]=rB[i]; rSeg[1]=rB[i+1]
        plot(zSeg,rSeg,linsty[bcB[i]]); plot(zSeg,-rSeg,linsty[bcB[i]])  
    
    # Animation
    if AnimateIons:
      xtxt = -0.25*maxz; ytxt = 0.8*maxr
      fig = figure(4)
      ax = axes( xlim=(min(zB),max(zB)), ylim=(-max(rB), max(rB)) )
      # Ions is a list plot lines. Each plot holds the locations of one isotope
      Ions = []
      for k in range(nIsotopes):
        l, = ax.plot([], [], 'o', color=clrs[k], ms=6, label='M=' + str(mIsotope[k]) + 'amu')
        Ions.append(l) 
      ttxt = ax.text(xtxt,ytxt,' ')     
    
      def init():
        animWM = get_current_fig_manager(); animWM.window.SetSize((1000,800))
        CPV1=contour(zc, rc,Vsum,25,colors='gray',linestyles='dashed')
        CPV2=contour(zc,-rc,Vsum,25,colors='gray',linestyles='dashed')
        clabel(CPV1, CPV1.levels, inline=True, fmt = '%g', fontsize=10)
        clabel(CPV2, CPV2.levels, inline=True, fmt = '%g', fontsize=10)
        for i in range(nB-1):
          zSeg[0]=zB[i]; zSeg[1]=zB[i+1]
          rSeg[0]=rB[i]; rSeg[1]=rB[i+1]
          plot(zSeg,rSeg,linsty[bcB[i]]); plot(zSeg,-rSeg,linsty[bcB[i]])  
        xlabel('z (mm)'); ylabel('r (mm)'); title('Ion trajectories')
        axis('equal')
        legend()
        for k in range(nIsotopes):
          Ions[k].set_data([], [])
        return 0
        #return Ions,
       
      def anim(i):
        # update pieces of the animation        
        nt = (i-nCoolFrames)*factAnim + ntEject
        if 0<=nt and nt<nplt:
          ttxt.set_text('t=%7.3f ms' % tplt[nt])
          i1=0
          for k in range(nIsotopes):
            i2 = i1 + NrIso[k]
            Ions[k].set_data(Xplt[i1:i2,0,nt],Xplt[i1:i2,1,nt])
            i1 = i2
        return 0

      ani = animation.FuncAnimation(fig, anim, init_func=init, frames=min(nAnimFrames,max(nXplt)), interval=1, blit=False)
      ani.save('IonAnimation.mp4')
    
    show()
