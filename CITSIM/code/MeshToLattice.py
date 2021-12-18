#!/usr/bin/python
# Makefile - Defines targets for mass spectrum computation stages
# MeshToLattice.py: Python script to process unstructured mesh and define geometry for
#                   lattice Boltzmann computation

from pylab import *    # Import scientific computation environment
import os,sys,time     # Operating system utilities
WorkDir = os.environ['PWD']
print "Generating lattice from mesh values in directory: ",WorkDir
os.chdir(WorkDir)
sys.path.insert(0, WorkDir)

def dist(i,j,n):
  return sqrt( (n[i,0]-n[j,0])**2 + (n[i,1]-n[j,1])**2 )

nDim=2; mG=512-1; mB=16; mX=ones(3,'int')      # Default values
import ITgeom
nDim=ITgeom.nDim;
print "nDim =",nDim," mG =",mG," mB =",mB

if nDim==2:
  mX[0]=ITgeom.Snx; mX[1]=ITgeom.Sny

if nDim==3:
  mX[0]=ITgeom.Snx; mX[1]=ITgeom.Sny; mX[2]=ITgeom.Snz
 
# Read ion trap geometry mesh
IonTrapMeshFile = open('IonTrapTh.msh','r')
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

# Process geometry
# 1. Find extents
Xmin=zeros(nDim,'double'); Xmax=zeros(nDim,'double'); dX=zeros(nDim,'double')
for i in range(nDim):
  Xmin[i] = ITnodes[:,i].min()
  Xmax[i] = ITnodes[:,i].max()
  dX[i] = Xmax[i] - Xmin[i]
# 2. Find min/max element edge lengths
minEdgeLen = ElemEdgeLen.min(); maxEdgeLen = ElemEdgeLen.max()
# 3. Construct a uniform grid over domain
imax = dX.argmax()
h  = dX[imax]/mG
mX = zeros(nDim,'int'); nX = zeros(nDim,'int')
for i in range(nDim):
  mX[i]=int(dX[i]/h); nX[i]=mX[i]+1
  if nX[i]%mB != 0:
    nX[i] = mB*(nX[i]/mB+1)
print "Lattice Boltzmann computational lattice with ",nX," nodes"
# 4. Rescale geometry: Xmin->Origin, unit length=h
for i in range(NrNodes):
  ITnodes[i,0:nDim] = (ITnodes[i,0:nDim]-Xmin[0:nDim])/h

nV = zeros(nDim+1,'int'); nV[0:nDim] = nX[0:nDim]; nV[nDim]=3;
nF = zeros(nDim+2,'int'); nF[0:nDim] = nX[0:nDim]; nF[nDim]=3; nF[nDim+1]=3
Ggrid=zeros(nX,dtype='int32',order='F')  # Geometry flags for lattice Boltzmann 
Vgrid=zeros(nV,dtype='double',order='F') # Grid values of electric potentials
Fgrid=zeros(nF,dtype='double',order='F') # Grid values of electric forces

print 'Done.'

