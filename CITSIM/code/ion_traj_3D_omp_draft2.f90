! ion-traj.f90 - Fortran routines to compute ion trajectories
MODULE CITSIMparams3D
  INTEGER, PARAMETER :: nDim=3,mxIons=10000,maxparams=512,nElec=3,nstages=3,nPDF=9,num_threads_max=24 !changed nDim=3 for 3D
  DOUBLE PRECISION, PARAMETER :: dtmin=1.0d-9
END MODULE CITSIMparams3D

SUBROUTINE bilinear3d(ie,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
  USE CITSIMparams3D
  IMPLICIT NONE
  INTEGER ie,i0,i1,j0,j1,nx,ny
  DOUBLE PRECISION :: dhx,dhy
  DOUBLE PRECISION :: aGrid(0:nx-1,0:ny-1,nDim-1,nElec)
  DOUBLE PRECISION :: Fcell(nDim-1,nDim-1,nDim-1),Feu(nDim-1) !changed nDim to ndim-1 for 3D (bruno)
  Fcell(1,1,:) = aGrid(i0,j0,:,ie); Fcell(2,1,:) = aGrid(i1,j0,:,ie)
  Fcell(1,2,:) = aGrid(i0,j1,:,ie); Fcell(2,2,:) = aGrid(i1,j1,:,ie)
  !print *,'Fcell=',Fcell(:,:,ie)
  Feu(:) = (1.d0-dhx)*(1.d0-dhy)*Fcell(1,1,:) + dhx*(1.d0-dhy)*Fcell(2,1,:) + &
           (1.d0-dhx)*dhy*Fcell(1,2,:) + dhx*dhy*Fcell(2,2,:)
END SUBROUTINE bilinear3d

SUBROUTINE ComputeElectrodeForces3d(tnow,params,InTrap,aGrid,X,F,nIons,nx,ny,cIons,chunk_size,cindex)
  USE CITSIMparams3D
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: twopi=6.28318530717959
  LOGICAL InTrap(nIons),OutOfTrap
  INTEGER, INTENT(IN) :: nIons,nx,ny,chunk_size,cindex
  DOUBLE PRECISION, INTENT(IN) :: cIons(mxIons)
  DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1),X(mxIons,nDim)
  DOUBLE PRECISION,INTENT(OUT) :: F(mxIons,nDim)
  INTEGER :: i0,i1,j0,j1,l,k0,k1, axfmode,axfnsteps
  !INTEGER, SAVE :: axfstep = 0. !not thread safe
  DOUBLE PRECISION :: tini,trmp,tfin,tnow,dt,Omega,OmegaEnt,OmegaExt,DelPhiEnt,DelPhiExt
  DOUBLE PRECISION :: Uent,Uext,ringV,entV,extV,entVRF,extVRF,Udc,dUdt,U0
  DOUBLE PRECISION :: dhx,dhy,ysign,hGrid
  DOUBLE PRECISION :: aGrid(0:nx-1,0:ny-1,nDim-1,nElec) !changed nDim to ndim-1 for 3D (bruno) 
  DOUBLE PRECISION :: Req,angle
  DOUBLE PRECISION :: Fl(nDim),Feu(nDim-1)
  DOUBLE PRECISION :: axfreq0,axfreq1,axfreqstep,axfdt,dfdt
  !Load parameters:
  tini=params(0); trmp=params(1); tfin=params(2); dt=params(3)
  Omega=params(8); OmegaEnt=params(9); OmegaExt=params(10)
  U0=params(4); dUdt=params(5); Uent=params(6); Uext=params(7); Udc=params(32)
  DelPhiEnt=params(11); DelPhiExt=params(12)
  entVRF=params(13); extVRF=params(14)
  hGrid=params(21)
  axfmode=params(56)       ! 0 = constant axial RF, 1 = linear incr. in axf,  2= linear incr. steps in axf, 3 = linear incr. steps in m/z
  axfnsteps=params(57)     ! nb of steps for axfreq 
  axfreq0=params(58)
  axfreq1=params(59)
  dfdt=params(60)          !axial freq. linear incr. rate
  !IF (tnow==tini+dt) THEN
  !axfstep=0
  !ENDIF
  !Electrodes voltages are set here:
  ringV = Udc+(U0 + dUdt*MAX(tnow-trmp,0.))*COS(Omega*tnow)
  entV = Uent + entVRF*COS(OmegaEnt*tnow+DelPhiEnt)   
  extV = Uext + extVRF*COS(OmegaExt*tnow+DelPhiExt)
  IF (axfmode==1) THEN
    entV = Uent + entVRF*COS(twopi*(axfreq1-dfdt*tnow)*tnow+DelPhiEnt)   
    extV = 0 !Uext + extVRF*COS(OmegaExt*tnow+DelPhiExt)
  END IF
  IF (axfmode==2) THEN
    axfreqstep=(axfreq1-axfreq0)/axfnsteps
    axfdt=(tfin-tini)/axfnsteps
    !IF (tnow >= (axfstep+1)*axfdt) THEN
    !    axfstep=axfstep+1
    !    print *,'axfstep=',axfstep
    !    print *,'int(tnow/axfdt)=',FLOOR(tnow/axfdt)
    !    print *,'axf =',(axfreq1-FLOOR(tnow/axfdt)*axfreqstep),'kHz'
    !    OPEN(UNIT=17,FILE='results/axf_ejection_log.txt',FORM='FORMATTED',POSITION='APPEND')
    !    WRITE(17,*) 'axf =',(axfreq1-FLOOR(tnow/axfdt)*axfreqstep),'kHz'
    !    CLOSE(17)
    !END IF
    entV = Uent + entVRF*COS(twopi*(axfreq1-FLOOR(tnow/axfdt)*axfreqstep)*tnow+DelPhiEnt)   
    extV = 0 !Uext + extVRF*COS(OmegaExt*tnow+DelPhiExt)
  END IF
  !Calculate the electric force
  F = 0.d0
  DO l=1,chunk_size
    IF (.NOT. InTrap(cindex+l)) RETURN
    Req=SQRT(X(cindex+l,2)**2+X(cindex+l,3)**2)
    angle=ATAN2(X(cindex+l,2),X(cindex+l,3)) !y=X(l,2),z=X(l,3) (bruno)
    i0=FLOOR(X(cindex+l,1)); j0=FLOOR(Req) !calculate corresponding r position for 3D (bruno)
    OutOfTrap = (i0<0) .OR. (i0>nx-1) .OR. (j0>ny-1)
    IF (OutOfTrap) RETURN
    ! Compute position within cell    
    dhx=X(cindex+l,1)-i0;      i1=i0+1
    dhy=Req-j0; j1=j0+1 !Req instead of X(l,2) for 3D (bruno)
    Fl(:) = 0.d0      
    CALL bilinear3d(1,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
    Fl(1) = Fl(1) + entV*Feu(1)                      ! Entry cap electrode
    Fl(2) = Fl(2) + SIN(angle)*entV*Feu(2)
    Fl(3) = Fl(3) + COS(angle)*entV*Feu(2) !Changed for 3D (bruno)
    CALL bilinear3d(2,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
    Fl(1) = Fl(1) + ringV*Feu(1)                     ! Ring electrode
    Fl(2) = Fl(2) + SIN(angle)*ringV*Feu(2)
    Fl(3) = Fl(3) + COS(angle)*ringV*Feu(2)!Changed for 3D (bruno)
    CALL bilinear3d(3,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
    Fl(1) = Fl(1) + extV*Feu(1)                      ! Exit cap electrode
    Fl(2) = Fl(2) + SIN(angle)*extV*Feu(2)
    Fl(3) = Fl(3) + COS(angle)*extV*Feu(2)!Changed for 3D (bruno)
    F(cindex+l,:) = cIons(cindex+l)*Fl(:)/hGrid     ! Scale electrode-induced accelerations to grid units
    !added *cIons(l) in previous line to avoid inconsistency in the way charges were accounted for the different forces in Verlet calc (bruno)
  END DO
END SUBROUTINE ComputeElectrodeForces3d

!~ SUBROUTINE ComputeIonIonForces3d(params,InTrap,cIons,Xion,Fion,nIons)
!~   USE, INTRINSIC :: ISO_C_BINDING  
!~   USE CITSIMparams3D
!~   IMPLICIT NONE
!~   LOGICAL InTrap(nIons)
!~   INTEGER, INTENT(IN) :: nIons
!~   DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1)
!~   DOUBLE PRECISION, INTENT(IN) :: cIons(mxIons)  
!~   DOUBLE PRECISION, INTENT(IN) :: Xion(mxIons,nDim)
!~   DOUBLE PRECISION, INTENT(OUT) :: Fion(mxIons,nDim)
!~   !
!~   LOGICAL IonIonGPU
!~   INTEGER i,j
!~   INTEGER (C_INT) :: N,ReturnCode
!~   REAL (C_FLOAT) :: X(4,mxIons),F(4,mxIons) !C_FLOAT for interoperability between c and fortran
!~   DOUBLE PRECISION :: rmin=0.01,rmin3=1d-6, EPS2=1d-6
!~   DOUBLE PRECISION :: dX(nDim),r12,r,r2,r3,Fij(nDim),fCoul,fCoulh,h
!~   !
!~   INTERFACE
!~     INTEGER (C_INT) FUNCTION ion_ion3d(N,X,F) BIND(C,NAME='ion_ion3d') !modif for 3D (bruno)
!~       USE, INTRINSIC :: ISO_C_BINDING
!~       IMPLICIT NONE
!~       INTEGER, PARAMETER :: nDim=3,mxIons=10000 !changed nDim=3, mxIons=10000 (bruno)
!~       INTEGER (C_INT), VALUE :: N   !value means that N is passed by value
!~       REAL (C_FLOAT) :: X(4,N),F(4,N) !modified for 3D (bruno)
!~     END FUNCTION ion_ion3d
!~   END INTERFACE
!~   !  
!~   IonIonGPU = params(16) > 0
!~   ! Scaling factor for Coulomb force computation
!~   fCoul = params(17); h = params(21); fCoulh= fCoul/h**3 !?? why (10.*h)**3, removed 10 (bruno)
!~   Fion = 0.d0
!~   IF (IonIonGPU) THEN
!~     N = nIons; F=0.
!~     !Assign fortran columns to C rows for compatibility 
!~     X(1,1:N) = Xion(1:N,1); X(2,1:N) = Xion(1:N,2)
!~     X(3,1:N) = Xion(1:N,3); X(4,1:N) = cIons(1:N)
!~     ! GPU call
!~     ReturnCode = ion_ion3d(N,X,F)
!~     !Assign returning forces C arrays rows to Fortran Fion columns 
!~     Fion(1:nIons,1) = -fCoulh*F(1,1:nIons); Fion(1:nIons,2) = -fCoulh*F(2,1:nIons);
!~     Fion(1:nIons,3) = -fCoulh*F(3,1:nIons)
!~   ELSE
!~     DO i=1,nIons
!~     !print *,'cIon(',i,')=',cIons(i)
!~       IF (.NOT. InTrap(i)) CYCLE
!~       DO j=i+1,nIons
!~         IF (.NOT. InTrap(j)) CYCLE
!~         dX(:) = Xion(j,:)-Xion(i,:) 
!~         r2 = dX(1)**2 + dX(2)**2+ dX(3)**2+EPS2; r = SQRT(r2); r3 = r*r2
!~         ! Note: grid units are used. Check for subgrid separation.
!~         !IF (r<rmin) THEN
!~         !  dX(:) = dX(:)*rmin/r; r3 = rmin3
!~         !END IF
!~         Fij(:) = dX(:)*cIons(i)*cIons(j)/r3 ! added *cIons(i)*cIons(j) to account for multiply charged ions (bruno)
!~         Fion(i,:) = Fion(i,:) - Fij(:) !Fij pointing from i to j so its contribution to Fion(i,:) need to be reverted (bruno)
!~         Fion(j,:) = Fion(j,:) + Fij(:) !j feels the opposite force (bruno)
!~       END DO
!~     END DO
!~     Fion(1:nIons,:) = fCoulh*Fion(1:nIons,:)
!~   END IF
!~ END SUBROUTINE ComputeIonIonForces3d
!~ 
!~ SUBROUTINE ComputeGasFlow3d(params,nB,i1,i2,j1,j2,flowBC,fGas)
!~   USE, INTRINSIC :: ISO_C_BINDING  
!~   IMPLICIT NONE
!~   INTEGER, PARAMETER :: nDim=3,maxparams=512,nPDF=9 !changed nDim=3 (bruno)
!~   DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1)
!~   INTEGER, INTENT(IN) :: nB
!~   INTEGER (C_INT), INTENT(IN) :: i1,i2,j1,j2
!~   DOUBLE PRECISION, INTENT(IN) :: flowBC(nB,nDim+3)
!~   INTEGER (C_INT) :: ReturnCode
!~   REAL (C_FLOAT) :: fGas(i1:i2,j1:j2,0:nPDF-1)
!~   !
!~   INTERFACE
!~     INTEGER (C_INT) FUNCTION LBM2d(i1,i2,j1,j2,fGas) BIND(C,NAME='LBM2d')
!~       USE, INTRINSIC :: ISO_C_BINDING
!~       IMPLICIT NONE
!~       INTEGER, PARAMETER :: nDim=3,nPDF=9 !changed nDim=3 (bruno)
!~       INTEGER (C_INT) :: i1,i2,j1,j2
!~       REAL (C_FLOAT) :: fGas(i1:i2,j1:j2,0:nPDF-1)
!~     END FUNCTION LBM2d
!~   END INTERFACE
!~   !  GPU call  
!~   !ReturnCode = LBM2d(i1,i2,j1,j2,fGas)
!~ END SUBROUTINE ComputeGasFlow3d
!~ 
!~ SUBROUTINE IonGasCollisionF3dCM(params,mIon,vxIon,vyIon,vzIon,FxIon,FyIon,FzIon, &
!~                             uGas,vGas,wGas,sigxGas,sigyGas,sigzGas)
!~   USE CITSIMparams3D
!~   IMPLICIT NONE
!~   DOUBLE PRECISION, PARAMETER :: twopi=6.28318530717959
!~   DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1),mIon,vxIon,vyIon,vzIon
!~   DOUBLE PRECISION, INTENT(IN) :: uGas,vGas,wgas,sigxGas,sigyGas,sigzGas
!~   DOUBLE PRECISION, INTENT(OUT) :: FxIon,FyIon,FzIon
!~   DOUBLE PRECISION :: dt,mGas,eRestitution,h
!~   DOUBLE PRECISION :: rU(6),radBM,thtBM,radBMz,thtBMz,z0,z1,z2
!~   DOUBLE PRECISION :: vxGas,vyGas,vzGas,vxIon2,vyIon2,vzIon2,vxGas2,vyGas2,vzGas2
!~   DOUBLE PRECISION :: impactangle,impacttheta,mVrelAz,mVrelEl,mVrelR
!~   DOUBLE PRECISION :: XImpact,YImpact,ZImpact,XImpact2,YImpact2,ZImpact2,XImpact3,YImpact3,ZImpact3,NImpact3
!~   DOUBLE PRECISION :: VIonCMrad,VGasCMrad,XtransVector,YtransVector,ZtransVector
!~   DOUBLE PRECISION :: VIonCMradaftercol,VxIonCMaftercol,VyIonCMaftercol,VzIonCMaftercol,VxIonLABaftercol,VyIonLABaftercol
!~   DOUBLE PRECISION :: VzIonLABaftercol
!~   !This is my collision algorithm where collisions are treated in the center of mass frame
!~   !Load parameters:
!~   dt = params(3);  mGas = params(18); eRestitution = params(121)
!~   h= params(21)
!~   !Generate random uniform numbers
!~   CALL RANDOM_NUMBER(rU)
!~   !Box-Muller transform for gaussian distributed random numbers
!~   radBM = SQRT(-2.*LOG(rU(2))); thtBM = twopi*rU(3)
!~   radBMz = SQRT(-2.*LOG(rU(4))); thtBMz = twopi*rU(5)
!~   z0 = radBM*COS(thtBM); z1 = radBM*SIN(thtBM); z2 = radBMz*COS(thtBMz)
!~   !Draw gas velocity from Maxwell distribution 
!~   !uGas, vGas, wGas are local average gas velocity provided by LBM or DSMC
!~   !as well as the MB standard deviation in (grid units)/ms
!~   vxGas = (sigxGas*z0+uGas); vyGas = (sigyGas*z1+vGas); vzGas = (sigzGas*z2+wgas)
!~   !We place ourselves in the reference frame of the center of mass
!~   vxIon2=vxIon-(mGas*vxGas+mIon*vxIon)/(mIon+mGas)
!~   vyIon2=vyIon-(mGas*vyGas+mIon*vyIon)/(mIon+mGas)
!~   vzIon2=vzIon-(mGas*vzGas+mIon*vzIon)/(mIon+mGas)
!~   vxGas2=vxGas-(mGas*vxGas+mIon*vxIon)/(mIon+mGas)
!~   vyGas2=vyGas-(mGas*vyGas+mIon*vyIon)/(mIon+mGas)
!~   vzGas2=vzGas-(mGas*vzGas+mIon*vzIon)/(mIon+mGas)
!~   !Pick random impact angles
!~   impactangle=asin(sqrt(0.999999999*rU(6)))
!~   impacttheta=twopi*rU(1)  ! Random contact angle, phi is angle along which collision occurs
!~   !Determination of impact position (standard convention for spheric coordinates)
!~   mVrelR=sqrt((vxGas-vxIon)**2+(vyGas-vyIon)**2+(vzGas-vzIon)**2)
!~   mVrelAz=atan2(vyIon2,vxIon2)
!~   mVrelEl=(twopi/4)-acos(vzIon2/mVrelR)
!~   !XImpact defines the position of incidence on a sphere of radius 1
!~   !in a reference frame whose z direction is oriented by mVrel 
!~   XImpact=cos(twopi/4-impactangle)*cos(impacttheta)
!~   YImpact=cos(twopi/4-impactangle)*sin(impacttheta)
!~   ZImpact=sin(twopi/4-impactangle)
!~   !transformation back into the reference frame with initial orientation
!~   !Impact3 gives the coordinate of the Impact in the initial reference frame
!~   XImpact2=ZImpact*sin(twopi/4-mVrelEl)+XImpact*cos(twopi/4-mVrelEl)
!~   YImpact2=YImpact
!~   ZImpact2=ZImpact*cos(twopi/4-mVrelEl)-XImpact*sin(twopi/4-mVrelEl)
!~   XImpact3=XImpact2*cos(mVrelAz)-YImpact2*sin(mVrelAz)
!~   YImpact3=XImpact2*sin(mVrelAz)+YImpact2*cos(mVrelAz)
!~   ZImpact3=ZImpact2
!~   NImpact3=sqrt(XImpact3**2+YImpact3**2+ZImpact3**2)
!~   !At this point we have determined the collision
!~   !We will now calculate the radial force exerted by each ball on each other
!~   !The collision is treated using velocities in the center of mass
!~   !Impact3 is chosen as norm vector along radial direction
!~   !calculating radial part of velocities in CM:
!~   VIonCMrad=-(vxIon2*XImpact3/NImpact3+vyIon2*YImpact3/NImpact3+vzIon2*ZImpact3/NImpact3)
!~   !nvion2*cos(impactangle) gives the same result!
!~   !calculating transverse vector
!~   XtransVector=vxIon2+VIonCMrad*XImpact3
!~   YtransVector=vyIon2+VIonCMrad*YImpact3
!~   ZtransVector=vzIon2+VIonCMrad*ZImpact3
!~   !Calculating resulting velocities after collision using radial components
!~   !of the force only:
!~   VIonCMradaftercol=-eRestitution*VIonCMrad
!~   !Recomposing ion velocities in CM after collision
!~   VxIonCMaftercol=-VIonCMradaftercol*XImpact3+XtransVector
!~   VyIonCMaftercol=-VIonCMradaftercol*YImpact3+YtransVector
!~   VzIonCMaftercol=-VIonCMradaftercol*ZImpact3+ZtransVector
!~   !Transformation back into the LAB reference frame of the velocities after
!~   !collision:
!~   VxIonLABaftercol=VxIonCMaftercol+(mGas*vxGas+mIon*vxIon)/(mIon+mGas)
!~   VyIonLABaftercol=VyIonCMaftercol+(mGas*vyGas+mIon*vyIon)/(mIon+mGas)
!~   VzIonLABaftercol=VzIonCMaftercol+(mGas*vzGas+mIon*vzIon)/(mIon+mGas)
!~   !Find force that applied over time step dt which would produce the same change of ion velocity
!~   !(Note: since velocity was in (grid units)/ms, force is in (amu) (grid units)/ms^2)
!~   FxIon = mIon*(VxIonLABaftercol-vxIon)/dt
!~   FyIon = mIon*(VyIonLABaftercol-vyIon)/dt
!~   FzIon = mIon*(VzIonLABaftercol-vzIon)/dt
!~ END SUBROUTINE IonGasCollisionF3dCM
!~ 
SUBROUTINE IonGasCollisionF3dSIM(params,mIon,vxIon,vyIon,vzIon,FxIon,FyIon,FzIon, &
                            uGas,vGas,wGas,sigxGas,sigyGas,sigzGas)
  USE CITSIMparams3D
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: twopi=6.28318530717959
  DOUBLE PRECISION, INTENT(IN) :: mIon,vxIon,vyIon,vzIon,params(0:maxparams-1)
  DOUBLE PRECISION, INTENT(IN) :: uGas,vGas,wgas,sigxGas,sigyGas,sigzGas
  DOUBLE PRECISION, INTENT(OUT) :: FxIon,FyIon,FzIon
  DOUBLE PRECISION :: dt,mGas,eRestitution,rU(6),radBM,thtBM,radBMz,thtBMz,z0,z1,z2,vxGas,vyGas,vzGas
  DOUBLE PRECISION :: vxion2,vyion2,vzion2,speedionr,azionr,elionr
  DOUBLE PRECISION :: impactangle,impacttheta,vtion,vrion,vrion2
  DOUBLE PRECISION :: v3(3),v4(3),v5(3),v6(3),v7(3)
  DOUBLE PRECISION :: cangle,sangle,ctheta,stheta,celback,selback,cazback,sazback

  !This is SIMION's algorithm retranscribed
  !Load parameters:
  dt = params(3);  mGas = params(18); eRestitution = params(121)
  !Generate uniform random numbers
  CALL RANDOM_NUMBER(rU)
  !Box-Muller transform
  radBM = SQRT(-2.*LOG(rU(2))); thtBM = twopi*rU(3)
  radBMz = SQRT(-2.*LOG(rU(4))); thtBMz = twopi*rU(5)
  z0 = radBM*COS(thtBM); z1 = radBM*SIN(thtBM); z2 = radBMz*COS(thtBMz)
  !Draw velocity from Maxwell distribution for buffer gas and express in (grid units)/ms
  !uGas, vGas, wGas are local average gas velocity provided by LBM or DSMC
  vxGas = (sigxGas*z0+uGas); vyGas = (sigyGas*z1+vGas); vzGas = (sigzGas*z2+wgas)
  !print *,"vxGas =",vxGas,", vyGas =",vyGas,", vzGas =",vzGas
  !We place ourselves in the reference frame of the gas particle
  vxion2=vxIon-vxGas
  vyion2=vyIon-vyGas
  vzion2=vzIon-vzGas
  !Transformation in spheric coordinate (SIMION's axis orientation):
  speedionr=sqrt(vxion2**2+vyion2**2+vzion2**2)
  azionr=-atan2(vzion2,vxion2)
  elionr=(twopi/4)-acos(vyIon2/speedionr)
  !Pick 2 angles for collision:
  impactangle=asin(sqrt(0.999999999*rU(6)))
  impacttheta=twopi*rU(1)
  !Projection of radial component onto line of collision:
  vrion=speedionr*cos(impactangle)
  vtion=speedionr*sin(impactangle)
  !collision
  vrion2=vrion*(mIon-mGas*eRestitution)/(mIon+mGas)
  !determination of line of collision
  !CALL elevation_rotate(twopi/4-impactangle,vrion2,vtion,0.d0,v3)
  cangle=cos(twopi/4-impactangle);sangle=sin(twopi/4-impactangle)
  v3(1)=vrion2*cangle-vtion*sangle
  v3(2)=vrion2*sangle+vtion*cangle
  v3(3)=0
  !CALL azimuth_rotate(impacttheta,v3(1),v3(2),v3(3),v4)
  ctheta = COS(impacttheta); stheta = SIN(impacttheta)
  v4(3)=v3(3)*ctheta-v3(1)*stheta
  v4(1)=v3(3)*stheta+v3(1)*ctheta
  v4(2)=v3(2)
  !CALL elevation_rotate(-twopi/4+elionr,v4(1),v4(2),v4(3),v5)
  celback= COS(-twopi/4+elionr); selback = SIN(-twopi/4+elionr)
  v5(1)=v4(1)*celback-v4(2)*selback
  v5(2)=v4(1)*selback+v4(2)*celback
  v5(3)=v4(3)
  !CALL azimuth_rotate(azionr,v5(1),v5(2),v5(3),v6)
  cazback= COS(azionr); sazback = SIN(azionr)
  v6(3)=v5(3)*cazback-v5(1)*sazback
  v6(1)=v5(3)*sazback+v5(1)*cazback
  v6(2)=v5(2)
  !Transform back into LAB reference frame
  v7(1)=v6(1)+vxGas
  v7(2)=v6(2)+vyGas
  v7(3)=v6(3)+vzGas
  ! Find force that applied over time step dt would produce same change of ion velocity
  ! (Note: since velocity was in (grid units)/ms, force is in (amu) (grid units)/ms^2
  FxIon = mIon*(v7(1)-vxIon)/dt
  FyIon = mIon*(v7(2)-vyIon)/dt
  FzIon = mIon*(v7(3)-vzIon)/dt
  !print *, 'nvIon= ',sqrt(vxIon**2+vyIon**2+vzIon**2),' nvIonaftercol= ',sqrt(v7(1)**2+v7(2)**2+v7(3)**2),&
  !' FxIon=',FxIon,' FyIon=',FyIon,' FzIon=',FzIon
END SUBROUTINE IonGasCollisionF3dSIM
!~ 
!~ SUBROUTINE IonGasCollisionF3dpointwise(params,mIon,vxIon,vyIon,vzIon,FxIon,FyIon,FzIon, &
!~                             uGas,vGas,wGas,sigxGas,sigyGas,sigzGas)
!~   USE CITSIMparams3D
!~   IMPLICIT NONE
!~   DOUBLE PRECISION, PARAMETER :: twopi=6.28318530717959
!~   DOUBLE PRECISION, INTENT(IN) :: mIon,vxIon,vyIon,vzIon,params(0:maxparams-1)
!~   DOUBLE PRECISION, INTENT(IN) :: uGas,vGas,wgas,sigxGas,sigyGas,sigzGas
!~   DOUBLE PRECISION, INTENT(OUT) :: FxIon,FyIon,FzIon
!~   DOUBLE PRECISION :: dt,mGas,eRestitution,rU(6),radBM,thtBM,radBMz,thtBMz,z0,z1,z2,vxGas,vyGas,vzGas,phi,cphi,sphi
!~   DOUBLE PRECISION :: mtot,vxion2,vyion2,vzion2
!~   
!~   !Load parameters:
!~   dt = params(3);  mGas = params(18); eRestitution = params(121)
!~   !Generate uniform random numbers
!~   CALL RANDOM_NUMBER(rU)
!~   !Box-Muller transform
!~   radBM = SQRT(-2.*LOG(rU(2))); thtBM = twopi*rU(3)
!~   radBMz = SQRT(-2.*LOG(rU(4))); thtBMz = twopi*rU(5)
!~   z0 = radBM*COS(thtBM); z1 = radBM*SIN(thtBM); z2 = radBMz*COS(thtBMz)
!~   !Draw velocity from Maxwell distribution for buffer gas and express in (grid units)/ms
!~   !uGas, vGas, wGas are local average gas velocity provided by LBM or DSMC
!~   vxGas = (sigxGas*z0+uGas); vyGas = (sigyGas*z1+vGas); vzGas = (sigzGas*z2+wgas)
!~   !Just expresses conservation of energy and momemtum:
!~   mtot=mIon+mGas
!~   vxion2=(mIon*vxIon+mGas*vxGas-eRestitution*mGas*(vxIon-vxGas))/mtot
!~   vyion2=(mIon*vyIon+mGas*vyGas-eRestitution*mGas*(vyIon-vyGas))/mtot
!~   vzion2=(mIon*vzIon+mGas*vzGas-eRestitution*mGas*(vzIon-vzGas))/mtot
!~   FxIon = mIon*(vxIon2-vxIon)/dt
!~   FyIon = mIon*(vxIon2-vyIon)/dt
!~   FzIon = mIon*(vxIon2-vzIon)/dt
!~ END SUBROUTINE IonGasCollisionF3dpointwise
!~ 
!~ 
!~ SUBROUTINE StoreCollision3d(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)
!~   USE CITSIMparams3D
!~   IMPLICIT NONE
!~   INTEGER, INTENT(IN) :: l,nIons,nplt
!~   DOUBLE PRECISION, INTENT(INOUT) :: Xcol(nIons,nDim,nplt),Vcol(nIons,nDim,nplt),Fcol(nIons,nDim,nplt)
!~   DOUBLE PRECISION, INTENT(IN) :: Xion(mxIons,nDim),Vion(mxIons,nDim),Fion(mxIons,nDim)
!~   !  
!~   INTEGER nxt
!~   !
!~   nxt = Xcol(l,1,nplt)
!~   IF (nxt >= nplt) RETURN  ! No more space to store collisions
!~   Xcol(l,1:nDim,nxt) = Xion(l,1:nDim)
!~   Vcol(l,1:nDim,nxt) = Vion(l,1:nDim)
!~   Fcol(l,1:nDim,nxt) = Fion(l,1:nDim)
!~   Xcol(l,1,nplt) = nxt + 1  
!~ END SUBROUTINE StoreCollision3d
!~ 

SUBROUTINE  ComputeIonGasForces3d(params,GasFlow,InTrap,mIon,cIon,csIon,Xion,Vion,nIons,Fion, &
                                  nplt,Xcol,Vcol,Fcol,flowData3d,mxLBM,myLBM,mzLBM,&
                                  x0LBM,y0LBM,z0LBM,dxLBM,dyLBM,dzLBM,tLBM,dtLBM,lion)
  USE, INTRINSIC :: ISO_C_BINDING
  USE CITSIMparams3D
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: GasFlow,InTrap(nIons) !GasFlow is not used, should be removed?
  INTEGER, INTENT(IN) :: nIons,nplt,mxLBM,myLBM,mzLBM
  DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1),mIon(nIons),cIon(nIons),csIon(nIons)
  DOUBLE PRECISION, INTENT(IN) :: Xion(mxIons,nDim),Vion(mxIons,nDim)
  DOUBLE PRECISION, INTENT(OUT) :: Fion(mxIons,nDim)
  DOUBLE PRECISION, INTENT(INOUT) :: Xcol(nIons,nDim,nplt),Vcol(nIons,nDim,nplt),Fcol(nIons,nDim,nplt)
  DOUBLE PRECISION, INTENT(IN) :: x0LBM,y0LBM,z0LBM,dxLBM,dyLBM,dzLBM,tLBM,dtLBM
  REAL (KIND=C_FLOAT), INTENT(IN) :: flowData3d(mxLBM,myLBM,mzLBM,8)  
  !
  INTEGER, PARAMETER :: Langevin=1, HardSphere=2, LatticeBoltzmann=3, DSMC=4, Hybrid=5,DSMC3d=6
  INTEGER :: CollisionModel,lion,iLBM,jLBM,kLBM
  DOUBLE PRECISION, PARAMETER :: sqrtpi=1.77245385090552,kBamu=8314.47148427469
  DOUBLE PRECISION :: rUnif,h,dt,Pcoll,pres,pres_temp,temp,mGas,mReduced
  DOUBLE PRECISION :: sig,sigma,sigmah,sigma2,sigLBM,sigxGas,sigyGas,sigzGas,uGas,vGas,wGas
  DOUBLE PRECISION :: AvgGasVel,MedGasVel,AvgRelVel,s,VionMag,CF,Pfactor,dist
  DOUBLE PRECISION :: Xmin(nDim),XLBM(nDim),RefLBMVel,LBMVel
  DOUBLE PRECISION :: AvgGasVelThermal,LBMVel2,MedGasVelThermal,invlambda
  
  Fion = 0.d0
  !Load parameters
  CollisionModel = params(20); h= params(21); dt = params(3); mGas = params(18)
  pres = params(101); temp = params(111)
  Xmin(1) = params(22); Xmin(2) = params(23)
  sigma = params(28); sigmah=sigma/h
  !If no collisions or wrong mGas or temp
  IF ((mGas<=0) .OR. (pres<=0) .OR. (temp<=0)) RETURN
  !Select collision model
  SELECT CASE (CollisionModel)
    CASE (Hybrid)
      ! Pressure and temperature defined in presGas[0] and tempGas[0]
      !params[30]=sqrt(8/pi * kBamu*tempGas[0]/mGas)/h # Average (mean) buffer gas speed in 3D in trajectory_3D.py
      !params[31]=sqrt(3 * kBamu*tempGas[0]/mGas)/h    # Median buffer gas speed in 3D (bruno) in trajectory_3D.py
      AvgGasVel = params(30)
      MedGasVel = params(31)
      !DO l=1,nIons
        !IF (.NOT. InTrap(lion)) CYCLE
        mReduced = mIon(lion)*mGas/(mIon(lion)+mGas)
        VionMag = SQRT(Vion(lion,1)**2 + Vion(lion,2)**2+ Vion(lion,3)**2) !added Vion(l,3)**2 for 3D (bruno)
       ! print *,'AvgGasVel=', AvgGasVel
        s = VionMag/MedGasVel
        !print *,'s=', s
        IF (s<1.0e-3) THEN
          AvgRelVel = AvgGasVel
        ELSE
          AvgRelVel = 0.5*AvgGasVel*((s+0.5/s)*sqrtpi*ERF(s) + EXP(-s**2))   
        END IF
        invlambda=csIon(lion)*params(29)/params(122) !invlambda to replace param(29)
        CF = 1. - EXP(-invlambda*AvgRelVel-params(26)*cIon(lion)/SQRT(mReduced))
        
        !print *,'tau/dt=',1/(invlambda*AvgRelVel+params(26)*cIon(l)/SQRT(mReduced))
        !print *,'CF=',CF
        !print *,'dt=',dt
        !print *,'params(29)=,',params(29)
        !print *,'params(122)=,',params(122)
        !print *,'params(26)=',params(26)
        !print *,'csIon(l)=',csIon(l)
        CALL RANDOM_NUMBER(rUnif)
        IF (rUnif < CF) THEN
          CALL IonGasCollisionF3dSIM(params,mIon(lion),Vion(lion,1),Vion(lion,2),Vion(lion,3),Fion(lion,1),&
          Fion(lion,2),Fion(lion,3),0.d0,0.d0,0.d0,sigmah,sigmah,sigmah)
          Xcol(lion,1,nplt) = Xcol(lion,1,nplt) + 1
          !CALL StoreCollision3d(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)        
        END IF
      !END DO
    CASE (DSMC3d)
      ! Extract buffer gas conditions
      !DO l=1,nIons
       ! IF (.NOT. InTrap(lion)) CYCLE   
        mReduced = mIon(lion)*mGas/(mIon(lion)+mGas)
        VionMag = SQRT(Vion(lion,1)**2 + Vion(lion,2)**2+ Vion(lion,3)**2) !added Vion(l,3)**2 for 3D (bruno)
        ! Interpolate LBM data to current ion position
        !
        ! Ion position in LBM reference frame, units (mm)
        XLBM(1:nDim) = h*Xion(lion,1:nDim) + Xmin(1:nDim)      
        iLBM = ABS((XLBM(1)-x0LBM)/dxLBM); iLBM = MAX(1,iLBM); iLBM = MIN(iLBM,mxLBM)
        jLBM = ABS((XLBM(2)-y0LBM)/dyLBM); jLBM = MAX(1,jLBM); jLBM = MIN(jLBM,myLBM)
        kLBM = ABS((XLBM(3)-z0LBM)/dzLBM); kLBM = MAX(1,kLBM); kLBM = MIN(kLBM,mzLBM) !added for 3D (bruno)
        !IF (flowData3d(iLBM,jLBM,kLBM,1)==0) THEN
        !IF (XLBM(2)>0.1 .OR. XLBM(3)>0.1) THEN
        !print *,'--------------------------------'
        !print *,'XLBM(1:nDim)',XLBM(1:nDim)
        !print *,'iLBM,jLBM, kLBM=',iLBM,jLBM,kLBM
        !print *,'P=',flowData3d(iLBM,jLBM,kLBM,1),', T=',flowData3d(iLBM,jLBM,kLBM,8)
        !print *,'--------------------------------'
        !END IF
        AvgGasVelThermal = sqrt(8./3.141592653589793*kBamu*flowData3d(iLBM,jLBM,kLBM,8)/mGas)/h ! Average (mean) buffer gas speed
        MedGasVelThermal = sqrt(3*kBamu*flowData3d(iLBM,jLBM,kLBM,8)/mGas)/h    ! Median buffer gas speed
        LBMVel = SQRT((flowData3d(iLBM,jLBM,kLBM,2)**2 + flowData3d(iLBM,jLBM,kLBM,3)**2+ flowData3d(iLBM,jLBM,kLBM,4)**2)/h**2)
        !print *,'LBMVel=',LBMVel
        AvgGasVel = AvgGasVelThermal + LBMVel
        !print *,'AvgGasVel=', AvgGasVel,'AvgGasVelThermal=', AvgGasVelThermal,'LBMVel=',LBMVel
        MedGasVel = MedGasVelThermal + LBMVel
        s = VionMag/MedGasVel
        !print *,"s=", s
        IF (s<1.0e-3) THEN
          AvgRelVel = AvgGasVel
        ELSE
          AvgRelVel = 0.5*AvgGasVel*((s+0.5/s)*sqrtpi*ERF(s) + EXP(-s**2))   
        END IF
        ! params[29]=dt*CollisionParameters[2]*presGas[0]*Torr/(kB*tempGas[0]) / 1.e9 * h
        ! params[26]=fLanCollP*(presGas[0]/tempGas[0])*dt*sqrt(alphae)
        ! params[130]=fLanCollP*dt*sqrt(alphae) #same as params[26] without presGas[0] and tempGas[0]
        ! params[131]=dt*Torr/kB/1.e9*h #same as params[29] without CollisionParameters[2],presGas[0] and tempGas[0]
        ! params(101)=presGas[0]
        pres_temp = flowData3d(iLBM,jLBM,kLBM,1)*760./101325./flowData3d(iLBM,jLBM,kLBM,8) !pres in Torr
        !print*,'flowData3d(iLBM,jLBM,kLBM,1)=',flowData3d(iLBM,jLBM,kLBM,1)
        !CF = 1. - EXP(-pres_temp*(csIon(l)*params(131)*AvgRelVel))
        CF = 1. - EXP(-pres_temp*(csIon(lion)*params(131)*AvgRelVel+params(130)*cIon(lion)/SQRT(mReduced))) 
        !print *,'CF=',CF
        !print *,'dt=',dt
        !print *,'params(29)=,',params(29)
        !print *,'params(29)~param(122)*params(131)*pres_temp=',params(122)*params(131)*pres_temp
        !print *,'params(122)=',params(122)
        !print *,'csIon(l)=',csIon(l)
        !print *,'params(26)=',params(26)
        !print *,'params(26)~params(130)*pres_temp=',params(130)*pres_temp
        !invlambda=csIon(l)*params(29)/params(122) !invlambda to replace param(29)
        !CF = 1. - EXP(-invlambda/params(101)*pres * AvgRelVel)
        CALL RANDOM_NUMBER(rUnif)
        IF (rUnif < CF) THEN
          ! Following are expressed as ratios and will be rescaled to sigmah in IonGasCollisionF
          uGas = flowData3d(iLBM,jLBM,kLBM,2)/h !non dim. units (lattice units)
          vGas = flowData3d(iLBM,jLBM,kLBM,3)/h !non dim. units (lattice units)
          wGas = flowData3d(iLBM,jLBM,kLBM,4)/h !non dim. units (lattice units)
          !print *,'uGas,vGas,wGas',uGas,vGas,wGas
          sigxGas = flowData3d(iLBM,jLBM,kLBM,5)/h
          sigyGas = flowData3d(iLBM,jLBM,kLBM,6)/h
          sigzGas = flowData3d(iLBM,jLBM,kLBM,7)/h
          CALL IonGasCollisionF3dSIM(params,mIon(lion),Vion(lion,1),Vion(lion,2),Vion(lion,3),Fion(lion,1),Fion(lion,2), &
                              Fion(lion,3),uGas,vGas,wGas,sigxGas,sigyGas,sigzGas)
          Xcol(lion,1,nplt) = Xcol(lion,1,nplt) + 1
          !CALL StoreCollision3d(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion) 
        END IF
  END SELECT  
END SUBROUTINE ComputeIonGasForces3d

SUBROUTINE trajectory3d(params,mIon,cIon,csIon,IonExit,nXplt,tplt,Xplt,Vplt,Fplt,Ncolplt,Xcol,Vcol,Fcol,AvgFeplt,AvgFcolplt &
                        ,Pplt,VelGasplt,aGrid,flowBC,nIons,nplt,nx,ny,nB)
  USE OMP_LIB
  USE, INTRINSIC :: ISO_C_BINDING 
  USE CITSIMparams3D
  IMPLICIT NONE
  ! Interface declarations
  INTEGER, INTENT(IN) :: nIons,nx,ny,nplt,nB
  INTEGER, INTENT(INOUT) :: nXplt(nIons)
  DOUBLE PRECISION, INTENT(INOUT) :: params(0:maxparams-1)
  DOUBLE PRECISION, INTENT(INOUT) :: tplt(nplt),Xplt(nIons,nDim,nplt),Vplt(nIons,nDim,nplt),Fplt(nIons,nDim,nplt)
  DOUBLE PRECISION, INTENT(INOUT) :: AvgFeplt(nIons,nplt),AvgFcolplt(nIons,nDim,nplt),Pplt(nIons,nplt),VelGasplt(nIons,nplt)
  DOUBLE PRECISION, INTENT(INOUT) :: Xcol(nIons,nDim,nplt),Vcol(nIons,nDim,nplt),Fcol(nIons,nDim,nplt),Ncolplt(nIons,nplt)
  DOUBLE PRECISION, INTENT(IN) :: mIon(nIons),cIon(nIons),csIon(nIons),aGrid(0:nx-1,0:ny-1,nDim-1,nElec)
  DOUBLE PRECISION, INTENT(IN) :: flowBC(nB,nDim+3)
  DOUBLE PRECISION, INTENT(INOUT) :: IonExit(nIons,0:4*nDim)
  ! Local variables
  INTEGER, PARAMETER :: Langevin=1, HardSphere=2, LatticeBoltzmann=3,DSMC=4, Hybrid=5, DSMC3d=6 
  DOUBLE PRECISION, PARAMETER :: twopi=6.283185307179586
  CHARACTER :: fname*16,line*80
  LOGICAL InTrap(nIons),Active(nIons),IonIon,IonGas,GasFlow
  INTEGER :: i,j,k,n,l,c,mLBM(nDim),iError,CollisionModel,tBirth,nbIonExit,cindex !inxtplt,iplt,its
  INTEGER (C_INT) :: nActive
  REAL (C_FLOAT), ALLOCATABLE, DIMENSION(:,:,:) :: fGas
  DOUBLE PRECISION :: dt2,tini,trmp,tfin,dt,U0,dUdt,Uent,Uext,Omega !tnow
  DOUBLE PRECISION :: dtplt,RFperiod,dphi,phi,cphi,sphi !tnxtplt
  DOUBLE PRECISION :: X(mxIons,nDim,nstages),V(mxIons,nDim,nstages)
  DOUBLE PRECISION :: XF(mxIons,nDim),VF(mxIons,nDim),F(mxIons,nDim)
  DOUBLE PRECISION :: Fionion(mxIons,nDim), Fiongas(mxIons,nDim), Fe(mxIons,nDim)
  DOUBLE PRECISION :: rExtentMin,rExtentMax,zExtentMin,zExtentMax,hGrid,Xmin(nDim)
  DOUBLE PRECISION :: AvgKESum(nIons),AvgDistSum(nIons),AvgDTCSum(nIons),AvgKEplt(nIons),AvgDistplt(nIons),AvgDTCplt(nIons)
  DOUBLE PRECISION :: AvgFesum(nIons),AvgFcolsum(nIons,nDim),tnow(num_threads_max),tnxtplt(num_threads_max)
  DOUBLE PRECISION :: d2(nDim),tBirths(nIons),dtBirth,test(2)
  INTEGER :: inxtplt(num_threads_max),iplt(num_threads_max),its(num_threads_max),chunk_size(num_threads_max)
  INTEGER :: mxLBM,myLBM,mzLBM,nCols,iLBM,jLBM,kLBM,axfmode,thread_num
  INTEGER :: num_threads
  DOUBLE PRECISION :: x0LBM,y0LBM,z0LBM,dxLBM,dyLBM,dzLBM,tLBM,dtLBM
  REAL (C_FLOAT), ALLOCATABLE, DIMENSION(:,:,:) :: flowData
  REAL (C_FLOAT), ALLOCATABLE, DIMENSION(:,:,:,:) :: flowData3d
  
  print *,"at ion_traj, start of trajectory3d subroutine, nDim=",nDim
  
  ! Return trajectories for n ions at nplt time slices in nDim dimensions
  
  ! Assume all ions are in trap initially
  InTrap = .TRUE.; Active = .FALSE.
  
  ! Load parameters:
  ! times
  axfmode=params(56)
  tini=params(0); trmp=params(1); tfin=params(2); dt=params(3); tBirth=params(52)
  dtplt=(tfin-tini)/nplt !Sorin used (nplt-1) (bruno)
  tnxtplt=tini+dtplt; inxtplt=2; iplt=1
  IF (dtplt<dt) dtplt=dt
  !print *,"(tfin-tini)/dt=",(tfin-tini)/dt !(bruno)
  print *,"1:dt=",dt,"dtplt=",dtplt
  ! RF amplitudes, ramping rate
  U0=params(4); dUdt=params(5); Uent=params(6); Uext=params(7)
  ! RF frequency
  Omega=params(8)
  ! Nb of threads
  num_threads=params(61)
  ! simulation options
  IonIon  = params(15) > 0
  IonGas  = params(20) > 0
  GasFlow = params(19) > 0
  CollisionModel = params(20)
  ! geometry parameters
  hGrid = params(21); Xmin(1) = params(22); Xmin(2) = params(23); Xmin(3) = params(24) !added Xmin(3) in 3D (bruno)
  print *,"Xmin(:)=",Xmin(:)
  rExtentMax = params(50)/hGrid; rExtentMin = -rExtentMax
  !yExtentMax = rExtentMax; xExtentMax = rExtentMax; !!rExtentMax and rExtentMin can be used directly 
  !yExtentMin = rExtentMin; xExtentMin = rExtentMin;
  zExtentMax = (params(51) - Xmin(1))/hGrid; zExtentMin = (-params(51) - Xmin(1))/hGrid
  print *,"rExtentMax,zExtentMax,rExtentMin,zExtentMin=",rExtentMax,zExtentMax,rExtentMin,zExtentMin
  
  ! dtBirth calculation and tBirths generation
  IF (tBirth==0) CALL RANDOM_NUMBER(tBirths) ! Generation of random times of birth (bruno)
  RFperiod=1./(Omega/twopi)                  ! RF period in ms (bruno)
  dtBirth=RFperiod/nIons                     ! Calculation of dtBirth for regular times of Birth
  tBirths=tBirths*RFperiod
  
  ! Xcol(:,1,nplt) initialization
  Xcol(:,1,nplt) = 1  ! Use the last 1,nplt component to identify next available storage for a collision
  
  ! Read results from LBM buffer gas flow simulation
  IF (CollisionModel .EQ. LatticeBoltzmann) THEN
    PRINT *,'Reading flow data ...'
    OPEN(UNIT=12,FILE='../flow/flow.data',STATUS='OLD')
    READ(12,*)mxLBM,myLBM,x0LBM,y0LBM,dxLBM,dyLBM,tLBM,dtLBM
    PRINT *,' Flow computation in z-r plane on ',mxLBM,' by ',myLBM,' lattice'
    ALLOCATE(flowData(mxLBM,myLBM,5),STAT=iError)
    IF (iError /= 0) THEN
      PRINT *,'Cannot allocate space for lattice Boltzmann flow computation'
      STOP
    END IF 
    DO i=1,mxLBM
      DO j=1,myLBM
        ! flowData contains p u v (SI units Pa, m/s, m/s), sigx sigy (normalized unit vector components, non-dimensional)
        READ(12,*)flowData(i,j,1:5)        
      END DO
    END DO
    PRINT *,' ... done'
    CLOSE(12)
  END IF
  IF (CollisionModel .EQ. DSMC3d) THEN                                                                            
    PRINT *,'Reading DSMC flow data ...'                                                                          
    OPEN(UNIT=12,FILE='../flow/DSMCflow3d.data',STATUS='OLD')  
    !OPEN(UNIT=13,FILE='../flow/DSMCflow3dtest.data',STATUS='NEW')                                                   
    READ(12,*)mxLBM,myLBM,mzLBM,x0LBM,y0LBM,z0LBM,dxLBM,dyLBM,dzLBM,tLBM,dtLBM                                    
    PRINT *,' Flow computation in 3d on ',mxLBM,' by ',myLBM,' by ',mzLBM,' lattice'                              
    ALLOCATE(flowData3d(mxLBM,myLBM,mzLBM,8),STAT=iError)                                                         
    IF (iError /= 0) THEN                                                                                         
      PRINT *,'Cannot allocate space for lattice Boltzmann flow computation'                                      
      STOP                                                                                                        
    END IF                                                                                                        
    DO k=1,mzLBM                                                                                                  
      DO j=1,myLBM                                                                                                
        DO i=1,mxLBM                                                                                              
          ! DSMCflowData contains p u v w (SI units Pa, m/s, m/s),                                                
          ! sigx sigy sigz (normalized unit vector components, non-dimensional),and T                             
          READ(12,*)flowData3d(i,j,k,1:8)                                                                         
        END DO                                                                                                    
      END DO                                                                                                      
    END DO                                                                                                        
    PRINT *,' ... done'  
!~     DO k=1,mzLBM                                                                                                  
!~       DO j=1,myLBM                                                                                                
!~         DO i=1,mxLBM                                                                                              
!~           ! DSMCflowData contains p u v w (SI units Pa, m/s, m/s),                                                
!~           ! sigx sigy sigz (normalized unit vector components, non-dimensional),and T                             
!~           WRITE(13,*)flowData3d(i,j,k,1:8)                                                                         
!~         END DO                                                                                                    
!~       END DO                                                                                                      
!~     END DO                                                                                            
    CLOSE(12) 
!~     CLOSE(13)                                                                                                    
  END IF   
  
  !Set counter and sums to zero (bruno)
  !for stability plots
  its=0
  AvgKESum=0.
  AvgDistSum=0.
  AvgDTCSum=0.
  
  ! Initial ion positions  
  DO l=1,nIons
    ! Scale coordinates such that grid spacing = unity (to enable quick identification of grid potential)
    X(l,1:nDim,1) = (Xplt(l,1:nDim,1) - Xmin(1:nDim))/hGrid
    ! Start Verlet algorithm with one Euler time step using given velocities
    X(l,1:nDim,2) = X(l,1:nDim,1) + dt*Vplt(l,1:nDim,1)/hGrid  
  END DO  
  
  ! first time increment after Euler time step
  tnow=tini+dt; dt2=dt**2
  nbIonExit=0
  
!~   PRINT *,'Active at start *****************************************************************************************************'
!~   PRINT *, Active
!~   PRINT *, "tBirth=",tBirth
!~   PRINT *,"before loop++++++++++++++++++++++++++++++++++++++++++++++++++"
!~   PRINT *,"tnow=",tnow
!~   PRINT *,"tfin=",tfin
  
  !*********************************PARALLEL CODE********************************************************************************
  !Define chunk_size(num_threads)
  chunk_size=nIons/num_threads
  IF (MOD(nIons,num_threads)==0) THEN
  chunk_size(num_threads)=nIons/num_threads
  ELSE
  chunk_size(num_threads)=chunk_size(1)+MOD(nIons,num_threads)
  END IF
  !$ call omp_set_num_threads(num_threads)
  !$omp parallel do private(thread_num,cindex) schedule(static)
  DO c=1,num_threads
  
  cindex=(c-1)*chunk_size(1)
  
  !$omp single
  print *, "num_threads=",num_threads
  !$omp end single
  
  
      !IF (.NOT. InTrap(l)) CYCLE

!~       !$omp critical
!~       !$ thread_num = omp_get_thread_num()
!~       !$ print *, "parallel do loop Start !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!"
!~       !$ print *, "This thread = ",thread_num,", #ion=",l
!~       !$ print *, "tnow=",tnow
!~         PRINT *,"******************************************************************************"
!~         PRINT *,Active
!~       !$omp end critical
!~      
      DO WHILE (tnow(c)<tfin)
      
        its(c)=its(c)+1  !time steps counter (bruno)
        !Active ions are created randomly within the first RF cycle or at regular dtBirths/nIons (bruno) 
!~         !$omp critical
!~         !$ print *, "This thread = ",thread_num, ", its=",its,", #ion=",l
!~         !$omp end critical
!~         
        
        IF (tnow(c)<RFperiod) THEN
        !$omp single
        SELECT CASE (tBirth)
        CASE (0) ! random times of birth        
!~         PRINT *,"tnow(c)=",tnow(c)
!~         PRINT *,"tBirths=",tBirths
        WHERE (tnow(c) .GT. tBirths) Active = .TRUE.
!~         PRINT *,"******************************************************************************"
!~         PRINT *,Active
!~         !$ thread_num = omp_get_thread_num()
!~         !$ print *, "This thread = ",thread_num, ", its=",its,", #chunk=",c
        CASE (1) ! regular times of birth
          nActive = MIN(CEILING((tnow(l)-tini)/dtBirth)*nIons,nIons)
          Active(1:nActive) = .TRUE.
        CASE (2) ! fixed times of birth
          IF (tnow(c)>params(53)) THEN
           Active(1:nIons) = .TRUE.
          END IF
        END SELECT
        !$omp end single
        END IF
        
!~         IF (.NOT. Active(l)) THEN
!~         tnow(l)=tnow(l)+dt
!~         CYCLE
!~         END IF

        XF(cindex+1:cindex+chunk_size(c),1:nDim) = X(cindex+1:cindex+chunk_size(c),1:nDim,2)  ! X for forces at midpoint    
        
        ! Determine if ions still inside the trap
        DO l=1,chunk_size(c)
            ! Conditions to determine if ions are still within the trap 
            IF (InTrap(cindex+l) .AND. ((XF(cindex+l,1)<0) .OR. (XF(cindex+l,1)>nx) .OR. &
             ((XF(cindex+l,2)**2+XF(cindex+l,3)**2)>ny**2))) THEN
            InTrap(cindex+l) = .FALSE.
            !nbIonExit=nbIonExit+1
            IonExit(cindex+l,0) = tnow(c)
            IonExit(cindex+l,1:nDim) = hGrid*X(cindex+l,1:nDim,2) + Xmin(1:nDim)
            IonExit(cindex+l,nDim+1:2*nDim) = (X(cindex+l,1:nDim,2)-X(cindex+l,1:nDim,1))/dt*hGrid
            IonExit(cindex+l,2*nDim+1:3*nDim) = hGrid*X(cindex+l,1:nDim,1) + Xmin(1:nDim)
            !$omp critical
            PRINT 1001,cindex+l,mIon(cindex+l),cIon(cindex+l),&
            csIon(cindex+l)*1e13,IonExit(cindex+l,1),IonExit(cindex+l,2)&
            ,IonExit(cindex+l,3),tnow(c),inxtplt(c)+1   !added IonExit(l,3) for 3D
            1001 FORMAT('Ion ',i5,' of mass ',F8.1,' amu, charge ', F4.1,'+ and ccs ',E10.4,&
            ' mm^2 has moved out of trap area at x='&
            ,F8.4,'(mm) y=',F8.4,'(mm) z=',F8.4,'(mm) t=',F10.5,'(ms) iplt=',i7)
            !$omp end critical
            !IF (axfmode==2) THEN
            !    OPEN(UNIT=17,FILE='results/axf_ejection_log.txt',FORM='FORMATTED',POSITION='APPEND')
            !    WRITE(17,*) 1001,l,mIon(l),cIon(l),csIon(l)*1e13,IonExit(l,1),IonExit(l,2),IonExit(l,3),tnow,inxtplt+1
            !    CLOSE(17)
            !END IF
            EXIT
            ELSE IF (InTrap(cindex+l)) THEN
            ! Calculate avg KE, avg flight distance, avg dist. to center (DTC) (bruno)
            d2=(X(cindex+l,1:nDim,2)-X(cindex+l,1:nDim,1))**2
            AvgKESum(cindex+l)=AvgKESum(cindex+l)+d2(1)+d2(2)+d2(3) !+d2(3) one more component in 3D (bruno)
            AvgDistSum(cindex+l)=AvgDistSum(cindex+l)+SQRT(d2(1)+d2(2)+d2(3))!+d2(3) one more component in 3D (bruno)
            AvgDTCSum(cindex+l)=AvgDTCSum(cindex+l)+SQRT((X(cindex+l,1,2)+Xmin(1)/hGrid)**2+&
            (X(cindex+l,2,2)-Xmin(2)/hGrid)**2+(X(cindex+l,3,2)-Xmin(3)/hGrid)**2)
            END IF
        END DO
        
        !  Electrode force
        CALL ComputeElectrodeForces3d(tnow(c),params,InTrap,aGrid,XF,F,nIons,nx,ny,cIon,chunk_size(c),cindex) !changed to ComputeElectrodeForces3d 
        Fe(cindex+1:cindex+chunk_size(c),:)=F(cindex+1:cindex+chunk_size(c),:)
        !$omp critical
        !$ thread_num = omp_get_thread_num()
        !$ print *, "This thread = ",thread_num, ", its=",its(1:num_threads),", #chunk=",c,", Chunk_size(c)=",chunk_size(c)
        print *,"F(cindex+1:cindex+chunk_size(c),1)=",F(cindex+1:cindex+chunk_size(c),1)
        !$omp end critical
!~         
!~ 
!~         !  Ion-ion interaction force
!~         !IF (IonIon) THEN
!~         !  CALL ComputeIonIonForces3d(params,InTrap,cIon,XF,Fionion,nIons)
!~         !  F(1:nIons,:) = F(1:nIons,:) + Fionion(1:nIons,:)
!~           !IF (its==1) THEN
!~           !  OPEN(UNIT=17,FILE='results/Fionion_160421.txt',STATUS='REPLACE',FORM='FORMATTED')
!~           !  DO j=1,nIons
!~           !  WRITE(17,*) j,XF(j,1),XF(j,2),XF(j,3),Fionion(j,1),Fionion(j,2),Fionion(j,3)
!~           !  END DO
!~           !CLOSE(17)
!~           !END IF 
!~         !END IF
!~         
        !  Ion-buffer interaction force
        !  At this stage X(2) is known and V(1.5)=(X(2)-X(1))/dt is O(dt^2) accurate
        !  Collision is assumed to occur from t1.5 to t2.5 and modeled as a force at t2
!~         IF (IonGas) THEN
!~          ! DO l=1,nIons
!~             VF(l,:) = (X(l,:,2) - X(l,:,1))/dt
!~           !END DO
!~           CALL ComputeIonGasForces3d(params,GasFlow,InTrap,mIon,cIon,csIon,XF,VF,nIons,Fiongas, &
!~                                      nplt,Xcol,Vcol,Fcol,flowData3d,mxLBM,myLBM,mzLBM,&
!~                                      x0LBM,y0LBM,z0LBM,dxLBM,dyLBM,dzLBM,tLBM,dtLBM,l)
!~           F(l,:) = F(l,:) + Fiongas(l,:)
!~         END IF
!~         
        ! Verlet update    
        DO l=1,chunk_size(c)
            IF (InTrap(cindex+l) .AND. Active(cindex+l)) THEN
                X(cindex+l,:,3) = 2.d0*X(cindex+l,:,2) - X(cindex+l,:,1) + dt2*F(cindex+l,:)/mIon(cindex+l) ! removed *cIon(l) due to inconsistency in the way charges were accounted for different forces (bruno)
                V(cindex+l,:,2) = (X(cindex+l,:,3) - X(cindex+l,:,1))/(2.d0*dt)
                !calculate avg forces
                AvgFesum(cindex+l)=AvgFesum(cindex+l)+SQRT(SUM(Fe(cindex+l,:)**2))
                AvgFcolsum(cindex+l,:)=AvgFcolsum(cindex+l,:)+Fiongas(cindex+l,:)
            END IF
        END DO
        
        ! Store position at this time if requested (change to maintain latest time steps in last 20% of nplt)
        IF (tnow(c)>=tnxtplt(c)) THEN
          IF (inxtplt(c)>nplt) PRINT *,'Warning current plot = ',inxtplt(c),' > nplt=',nplt,' (overwriting last plot)'!,'nbIonExit=',nbIonExit
          iplt(c)=MIN(inxtplt(c),nplt)
          tplt(iplt(c)) = tnow(c)
          DO l=1,chunk_size(c)
            IF (InTrap(cindex+l)) THEN
              Xplt(cindex+l,1:nDim,iplt(c)) = hGrid*X(cindex+l,1:nDim,2) + Xmin(1:nDim) ! Scaled back from grid coordinates to physical coordinates
              Vplt(cindex+l,1:nDim,iplt(c)) = hGrid*V(cindex+l,1:nDim,2)
              Fplt(cindex+l,1:nDim,iplt(c)) = F(cindex+l,1:nDim)
              Ncolplt(cindex+l,iplt(c))=Xcol(cindex+l,1,nplt)-1 !to save the nb of collisions along a traj
              !$omp critical
              !$ thread_num = omp_get_thread_num()
              print *,"----------------------------inside recording section-------------"
              !$ print *, "This thread = ",thread_num, ", its=",its(1:num_threads),", #chunk=",c,", Chunk_size(c)=",chunk_size(c)
              print *,"iplt(c)=",iplt(c)
              print *,"Xplt(cindex+1:cindex+chunk_size(c),1,iplt(c))=",F(cindex+1:cindex+chunk_size(c),1)
              !$omp end critical
              IF (CollisionModel .EQ. DSMC3d) THEN
                  ! Ion position in LBM reference frame, units (mm)
                  iLBM = ABS((Xplt(cindex+l,1,iplt(c))-x0LBM)/dxLBM); iLBM = MAX(1,iLBM); iLBM = MIN(iLBM,mxLBM)
                  jLBM = ABS((Xplt(cindex+l,2,iplt(c))-y0LBM)/dyLBM); jLBM = MAX(1,jLBM); jLBM = MIN(jLBM,myLBM)
                  kLBM = ABS((Xplt(cindex+l,3,iplt(c))-z0LBM)/dzLBM); kLBM = MAX(1,kLBM); kLBM = MIN(kLBM,mzLBM) !added for 3D (bruno)
                  Pplt(cindex+l,iplt(c))=flowData3d(iLBM,jLBM,kLBM,1)
                  VelGasplt(cindex+l,iplt(c))=SQRT(flowData3d(iLBM,jLBM,kLBM,2)**2 + flowData3d(iLBM,jLBM,kLBM,3)**2+ &
                   flowData3d(iLBM,jLBM,kLBM,4)**2)
              END IF
              IF (iplt(c)>1) THEN
              AvgFcolplt(cindex+l,:,iplt(c))=AvgFcolsum(cindex+l,:)/&
              (Ncolplt(cindex+l,iplt(c))-Ncolplt(cindex+l,iplt(c)-1))
              AvgFeplt(cindex+l,iplt(c))=AvgFesum(cindex+l)/&
              ((tplt(iplt(c))-tplt(iplt(c)-1))/dt)
              AvgFesum(cindex+l)=0
              AvgFcolsum(cindex+l,:)=0
              ELSE
              AvgFcolplt(cindex+l,:,iplt(c))=AvgFcolsum(cindex+l,:)/(Ncolplt(cindex+l,iplt(c)))
              AvgFeplt(cindex+l,iplt(c))=AvgFesum(cindex+l)/(tplt(iplt(c))/dt)
              AvgFesum(cindex+l)=0
              AvgFcolsum(cindex+l,:)=0
              END IF
              nXplt(cindex+l) = iplt(c)
            END IF
          END DO
          tnxtplt(c)=tnxtplt(c)+dtplt; inxtplt(c)=inxtplt(c)+1
        END IF
        
        ! Increment time (moved the time increment after the Xplt update, Bruno)
        ! in this way we store position and velocity of euler time step as well, and be consistent for the rest tof the time steps.
        tnow(c)=tnow(c)+dt
        
        ! Update positions and velocities 
        DO l=1,chunk_size(c)
          IF (InTrap(cindex+l) .AND. Active(cindex+l)) THEN
          X(cindex+l,1:nDim,1)=X(cindex+l,1:nDim,2); X(cindex+l,1:nDim,2)=X(cindex+l,1:nDim,3)
          V(cindex+l,1:nDim,1)=V(cindex+l,1:nDim,2)
          END IF
        END DO
        
    END DO ! End of time Main Loop
  
  END DO ! End nIons main loop
  !$omp end parallel do
  
!~   PRINT *,"Active at end *****************************************************************************************************"
!~   PRINT *, Active
!~   PRINT *, "tBirth=",tBirth
!~   PRINT *,"after loop++++++++++++++++++++++++++++++++++++++++++++++++++"
!~   PRINT *,"tnow=",tnow
!~   PRINT *,"tfin=",tfin
   
  ! Use params to transfer calculated means (bruno)
!  AvgKEplt=0.
!  AvgDistplt=0.
!  AvgDTCplt=0.
!  params(9) = iplt
!~   params(149) = its !(bruno)
!~   IF (nIons<116) THEN ! can only be transfered with params if nIons is small
!~   params(150:150+3*nIons-1)=0
!~     DO l=1,nIons !1.66053892e-27kg/amu / 1.60217657e-19J/eV (=) 1.03642692e-8 = factor to have KE in eV
!~       params(150:150+nIons-1)=AvgKESum(:)*(hGrid/dt)**2*0.5*1.03642692e-8*mIon(l)/its  
!~       params(150+nIons:150+2*nIons-1)=AvgDistSum(:)*hGrid
!~       params(150+2*nIons:150+3*nIons-1)=AvgDTCSum(:)*hGrid/its 
!~       !DO i=1,iplt-1
!~       !  AvgKEplt(l)=AvgKEplt(l)+0.5*mIon(l)*1.03642692e-8*(Vplt(l,1,i)**2+Vplt(l,2,i)**2)
!~       !  AvgDistplt(l)=AvgDistplt(l)+SQRT((Xplt(l,1,i+1)-Xplt(l,1,i))**2+(Xplt(l,2,i+1)-Xplt(l,2,i))**2)
!~       !  AvgDTCplt(l)=AvgDTCplt(l)+SQRT((Xplt(l,1,i))**2+(Xplt(l,2,i))**2)
!~       !END DO
!~     END DO  
!~   END IF
  
  !print *,'params(199:205)',params(199:205)
  !print *,'iplt,...:',iplt,AvgKEplt/(iplt-1),AvgDistplt,AvgDTCplt/(iplt-1)
  !print *, 'Xcol(:,1,nplt)',Xcol(:,1,nplt)
  ! Write one test trajectory to a file (bruno)
  !DO i=1,nIons
  !  OPEN(UNIT=17,FILE='results/trajectory.txt',STATUS='REPLACE',FORM='FORMATTED')
  !    !WRITE(17,*) params(9),params(199:202)
  !    DO j=11,nplt
  !    WRITE(17,*) j-10,tplt(j),Xplt(i,:,j),Vplt(i,:,j)
  !  END DO
  !END DO 
  !CLOSE(17) 
  !dphi=twopi/nIons
  !DO i=1,nIons
  !  phi=i*dphi; cphi=COS(phi); sphi=SIN(phi)
  !  ! Store trajectories in data files
  !  WRITE(fname,'(A4,I4.4,A4)')'traj',i,'.vtk'
  !  OPEN(UNIT=17,FILE=fname,STATUS='REPLACE',FORM='FORMATTED')
  !  WRITE(line,*) '# vtk DataFile Version 3.0'; WRITE(17,'(A80)')ADJUSTL(line)
  !  WRITE(line,*) 'Ion ',i,' trajectory, m=',mIon(i); WRITE(17,'(A80)')ADJUSTL(line)
  !  WRITE(line,*) 'ASCII'; WRITE(17,'(A80)')ADJUSTL(line)
  !  WRITE(line,*) 'DATASET POLYDATA'; WRITE(17,'(A80)')ADJUSTL(line)
  !  WRITE(line,*) 'POINTS ',nXplt(i),' float'; WRITE(17,'(A80)')ADJUSTL(line)  
  !  DO n=1,nXplt(i)
  !    WRITE(line,*)Xplt(i,1,n),Xplt(i,2,n)*cphi,Xplt(i,2,n)*sphi; WRITE(17,'(A80)')ADJUSTL(line)
  !  END DO
  !  WRITE(17,*)
  !  WRITE(line,*)'LINES ',nXplt(i)-1,3*(nXplt(i)-1); WRITE(17,'(A80)')ADJUSTL(line)
  !  DO n=1,nXplt(i)-1
  !   WRITE(line,*)2,n,n+1; WRITE(17,'(A80)')ADJUSTL(line)
  !  END DO
  !  WRITE(17,*)
  !  WRITE(line,*)'POINT_DATA ',nXplt(i); WRITE(17,'(A80)')ADJUSTL(line)
  !  WRITE(line,*)'SCALARS Velocity float'; WRITE(17,'(A80)')ADJUSTL(line)
  !  WRITE(line,*)'LOOKUP_TABLE default'; WRITE(17,'(A80)')ADJUSTL(line)
  !  DO n=1,nXplt(i)
  !    WRITE(line,*) SQRT(SUM(Vplt(i,1:nDim,n)**2));  WRITE(17,'(A80)')ADJUSTL(line)
  !  END DO
  !  CLOSE(17)
  !  ! Store collisions in data files
  !  nCols = Xcol(i,1,nplt)-1
  !  IF ((CollisionModel>0) .AND. (nCols>0)) THEN
  !    WRITE(fname,'(A4,I4.4,A4)')'coll',i,'.vtk'
  !    OPEN(UNIT=17,FILE=fname,STATUS='REPLACE',FORM='FORMATTED')
  !    WRITE(line,*) '# vtk DataFile Version 3.0'; WRITE(17,'(A80)')ADJUSTL(line)
  !    WRITE(line,*) 'Ion ',i,' collisions, m=',mIon(i); WRITE(17,'(A80)')ADJUSTL(line)
  !    WRITE(line,*) 'ASCII'; WRITE(17,'(A80)')ADJUSTL(line)
  !    WRITE(line,*) 'DATASET STRUCTURED_GRID'; WRITE(17,'(A80)')ADJUSTL(line)
  !    WRITE(line,*) 'DIMENSIONS ',nCols,' 1 1'; WRITE(17,'(A80)')ADJUSTL(line)
  !    WRITE(line,*) 'POINTS ',nCols,' float'; WRITE(17,'(A80)')ADJUSTL(line)  
  !    DO n=1,nCols
  !      Xcol(i,1:nDim,n) = hGrid*Xcol(i,1:nDim,n) + Xmin(1:nDim) ! Scaled back from grid coordinates to physical coordinates
 !       !Vcol(i,1:nDim,n) = hGrid*Vcol(i,1:nDim,n)
  !      WRITE(line,*)Xcol(i,1,n),Xcol(i,2,n)*cphi,Xcol(i,2,n)*sphi; WRITE(17,'(A80)')ADJUSTL(line)
  !    END DO
  !    WRITE(17,*)
  !    WRITE(line,*)'POINT_DATA ',nCols; WRITE(17,'(A80)')ADJUSTL(line)
  !    WRITE(line,*)'VECTORS CollForces float'; WRITE(17,'(A80)')ADJUSTL(line)
  !    DO n=1,nCols
  !     WRITE(line,*)Fcol(i,1,n),Fcol(i,2,n)*cphi,Fcol(i,2,n)*sphi; WRITE(17,'(A80)')ADJUSTL(line)
  !    END DO   
  !    CLOSE(17)
  !  END IF
  !END DO
END SUBROUTINE trajectory3d
