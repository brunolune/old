! ion-traj.f90 - Fortran routines to compute ion trajectories
MODULE CITSIMparams2D
  INTEGER, PARAMETER :: nDim=2,mxIons=10000,maxparams=512,nElec=3,nstages=3,nPDF=9
  DOUBLE PRECISION, PARAMETER :: dtmin=1.0d-9
END MODULE CITSIMparams2D

SUBROUTINE bilinear(ie,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
  USE CITSIMparams2D
  IMPLICIT NONE
  INTEGER ie,i0,i1,j0,j1,nx,ny
  DOUBLE PRECISION :: dhx,dhy
  DOUBLE PRECISION :: aGrid(0:nx-1,0:ny-1,nDim,nElec)
  DOUBLE PRECISION :: Fcell(nDim,nDim,nDim),Feu(nDim)
  Fcell(1,1,:) = aGrid(i0,j0,:,ie); Fcell(2,1,:) = aGrid(i1,j0,:,ie)
  Fcell(1,2,:) = aGrid(i0,j1,:,ie); Fcell(2,2,:) = aGrid(i1,j1,:,ie)
  !print *,'Fcell=',Fcell(:,:,ie)
  Feu(:) = (1.d0-dhx)*(1.d0-dhy)*Fcell(1,1,:) + dhx*(1.d0-dhy)*Fcell(2,1,:) + &
           (1.d0-dhx)*dhy*Fcell(1,2,:) + dhx*dhy*Fcell(2,2,:)
END SUBROUTINE bilinear

SUBROUTINE ComputeElectrodeForces2d(tnow,params,InTrap,aGrid,X,F,nIons,nx,ny,cIons)
  USE CITSIMparams2D
  IMPLICIT NONE
  LOGICAL InTrap(nIons),OutOfTrap
  INTEGER, INTENT(IN) :: nIons,nx,ny
  DOUBLE PRECISION, INTENT(IN) :: cIons(mxIons)
  DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1),X(mxIons,nDim)
  DOUBLE PRECISION,INTENT(OUT) :: F(mxIons,nDim)
  INTEGER :: i0,i1,j0,j1,l
  DOUBLE PRECISION :: tini,trmp,tfin,tnow,dt,Omega,OmegaEnt,OmegaExt,DelPhiEnt,DelPhiExt
  DOUBLE PRECISION :: Uent,Uext,ringV,entV,extV,entVRF,extVRF,Udc,dUdt,U0
  DOUBLE PRECISION :: dhx,dhy,ysign,hGrid
  DOUBLE PRECISION :: aGrid(0:nx-1,0:ny-1,nDim,nElec)
  DOUBLE PRECISION :: Fl(nDim),Feu(nDim)
  !Load parameters
  tini=params(0); trmp=params(1); tfin=params(2); dt=params(3)
  U0=params(4); dUdt=params(5); Uent=params(6); Uext=params(7); Udc=params(32)
  Omega=params(8); OmegaEnt=params(9); OmegaExt=params(10)
  DelPhiEnt=params(11); DelPhiExt=params(12)
  entVRF=params(13); extVRF=params(14)
  hGrid=params(21)
  !
  ringV = Udc+(U0 + dUdt*MAX(tnow-trmp,0.))*COS(Omega*tnow)
  entV = Uent + entVRF*COS(OmegaEnt*tnow+DelPhiEnt)   
  extV = Uext + extVRF*COS(OmegaExt*tnow+DelPhiExt)
  F = 0.d0
  DO l=1,nIons
    IF (.NOT. InTrap(l)) CYCLE    
    IF (j0>ny-1) THEN
    END IF
    i0=FLOOR(X(l,1)); j0=FLOOR(X(l,2))
    IF (j0<0) THEN
      ysign = -1.0
      j0 = -j0
    ELSE
      ysign = 1.0
    END IF
    OutOfTrap = (i0<0) .OR. (i0>nx-1) .OR. (j0>ny-1)
    IF (OutOfTrap) CYCLE
    ! Compute position within cell    
    dhx=X(l,1)-i0;      i1=i0+1
    dhy=ABS(X(l,2))-j0; j1=j0+1
    Fl(:) = 0.d0      
    CALL bilinear(1,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
    Fl(1) = Fl(1) + entV*Feu(1)                      ! Entry cap electrode
    Fl(2) = Fl(2) + ysign*entV*Feu(2)
    CALL bilinear(2,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
    Fl(1) = Fl(1) + ringV*Feu(1)                     ! Ring electrode
    Fl(2) = Fl(2) + ysign*ringV*Feu(2) 
    CALL bilinear(3,i0,i1,j0,j1,nx,ny,dhx,dhy,aGrid,Feu)
    Fl(1) = Fl(1) + extV*Feu(1)                      ! Exit cap electrode
    Fl(2) = Fl(2) + ysign*extV*Feu(2)
    F(l,:) = cIons(l)*Fl(:)/hGrid                          ! Scale electrode-induced accelerations to grid units
  END DO
END SUBROUTINE ComputeElectrodeForces2d

SUBROUTINE ComputeIonIonForces2d(params,InTrap,cIons,Xion,Fion,nIons) !to be modified not done in 2D
  USE, INTRINSIC :: ISO_C_BINDING  
  USE CITSIMparams2D
  IMPLICIT NONE
  LOGICAL InTrap(nIons)
  INTEGER, INTENT(IN) :: nIons
  DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1)
  DOUBLE PRECISION, INTENT(IN) :: cIons(mxIons)  
  DOUBLE PRECISION, INTENT(IN) :: Xion(mxIons,nDim)
  DOUBLE PRECISION, INTENT(OUT) :: Fion(mxIons,nDim)
  !
  LOGICAL IonIonGPU
  INTEGER i,j
  INTEGER (C_INT) :: N,ReturnCode
  REAL (C_FLOAT) :: X(mxIons),Y(mxIons),F(mxIons),G(mxIons) !C_FLOAT just tell that these arrays are interoperable between c and fortran
  DOUBLE PRECISION :: rmin=0.5,rmin3=0.125
  DOUBLE PRECISION :: dX(nDim),r12,r,r2,r3,Fij(nDim),fCoul,fCoulh,h
  !
  INTERFACE
    INTEGER (C_INT) FUNCTION ion_ion2d(N,X,Y,F,G) BIND(C,NAME='ion_ion2d')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INTEGER, PARAMETER :: nDim=2,mxIons=1000
      INTEGER (C_INT), VALUE :: N
      REAL (C_FLOAT) :: X(N),Y(N),F(N),G(N)
    END FUNCTION ion_ion2d
  END INTERFACE
  !  
  IonIonGPU = params(16) > 0
  ! Scaling factor for Coulomb force computation
  fCoul = params(17); h = params(21); fCoulh= fCoul/(10.*h)**3  
  Fion = 0.d0
  IF (IonIonGPU) THEN
    N = nIons; F=0.; G=0.
    X(1:N) = Xion(1:N,1); Y(1:N) = Xion(1:N,2)
    ! GPU call
    !ReturnCode = ion_ion2d(N,X,Y,F,G)
    Fion(1:nIons,1) = F(1:nIons); Fion(1:nIons,2) = G(1:nIons)
  ELSE
    DO i=1,nIons
      IF (.NOT. InTrap(i)) CYCLE
      DO j=i+1,nIons
        IF (.NOT. InTrap(j)) CYCLE
        dX(:) = Xion(j,:)-Xion(i,:)
        r2 = dX(1)**2 + dX(2)**2; r = SQRT(r2); r3 = r*r2
        ! Note: grid units are used. Check for subgrid separation.
        IF (r<rmin) THEN
          dX(:) = dX(:)*rmin/r; r3 = rmin3
        END IF        
        Fij(:) = dX(:)/r3
        Fion(i,:) = Fion(i,:) - Fij(:)
        Fion(j,:) = Fion(j,:) + Fij(:)
      END DO
    END DO
    Fion(1:nIons,:) = fCoulh*Fion(1:nIons,:)
  END IF
END SUBROUTINE ComputeIonIonForces2d

SUBROUTINE ComputeGasFlow2d(params,nB,i1,i2,j1,j2,flowBC,fGas)
  USE, INTRINSIC :: ISO_C_BINDING  
  IMPLICIT NONE
  INTEGER, PARAMETER :: nDim=2,maxparams=512,nPDF=9
  DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1)
  INTEGER, INTENT(IN) :: nB
  INTEGER (C_INT), INTENT(IN) :: i1,i2,j1,j2
  DOUBLE PRECISION, INTENT(IN) :: flowBC(nB,nDim+3)
  INTEGER (C_INT) :: ReturnCode
  REAL (C_FLOAT) :: fGas(i1:i2,j1:j2,0:nPDF-1)
  !
  INTERFACE
    INTEGER (C_INT) FUNCTION LBM2d(i1,i2,j1,j2,fGas) BIND(C,NAME='LBM2d')
      USE, INTRINSIC :: ISO_C_BINDING
      IMPLICIT NONE
      INTEGER, PARAMETER :: nDim=2,nPDF=9
      INTEGER (C_INT) :: i1,i2,j1,j2
      REAL (C_FLOAT) :: fGas(i1:i2,j1:j2,0:nPDF-1)
    END FUNCTION LBM2d
  END INTERFACE
  !  GPU call  
  !ReturnCode = LBM2d(i1,i2,j1,j2,fGas)
END SUBROUTINE ComputeGasFlow2d

SUBROUTINE IonGasCollisionF(params,mIon,vxIon,vyIon,FxIon,FyIon, &
                            uGas,vGas,sigxGas,sigyGas)
  USE CITSIMparams2D
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: twopi=6.28318530717959
  DOUBLE PRECISION, INTENT(IN) :: mIon,vxIon,vyIon,params(0:maxparams-1)
  DOUBLE PRECISION, INTENT(IN) :: uGas,vGas,sigxGas,sigyGas
  DOUBLE PRECISION, INTENT(OUT) :: FxIon,FyIon
  DOUBLE PRECISION :: rU(3),phi,cphi,sphi,radBM,thtBM,z0,z1,vxGas,vyGas,dt,mGas
  DOUBLE PRECISION :: eRestitution,vGasC,vIonC,vGasN,vIonN,wIonC,wIonN,FIonC,h
  !
  dt = params(3);  mGas = params(18); eRestitution = params(121)
  h= params(21);
  CALL RANDOM_NUMBER(rU)
  phi = twopi*rU(1)  ! Random contact angle, phi is angle along which collision occurs
  ! Box-Muller transform
  radBM = SQRT(-2.*LOG(rU(2))); thtBM = twopi*rU(3)
  z0 = radBM*COS(thtBM); z1 = radBM*SIN(thtBM)
  ! Draw velocity from Maxwell distribution for buffer gas and express in (grid units)/ms      
  vxGas = (sigxGas*z0+uGas); vyGas = (sigyGas*z1+vGas) !!!!I just suppressed *h which was for all results obtained on 060914 :(!!!
  cphi = COS(phi); sphi = SIN(phi)
  vGasC =  vxGas*cphi + vyGas*sphi         ! gas velocity along collision direction
  vIonC =  vxIon*cphi + vyIon*sphi         ! ion velocity along collision direction
  vGasN = -vxGas*sphi + vyGas*cphi         ! gas velocity normal to collision
  vIonN = -vxIon*sphi + vyIon*cphi         ! ion velocity normal to collision
  wIonC = (mIon*vIonC + mGas*vGasC + mGas*eRestitution*(vGasC-vIonC))/(mIon+mGas)
  wIonN = vIonN
  ! Find force that applied over time step dt would produce same change of ion velocity
  ! (Note: since velocity was in (grid units)/ms, force is in (amu) (grid units)/ms^2
  FIonC = mIon*(wIonC-vIonC)/dt
  FxIon = FIonC*cphi; FyIon = FIonC*sphi
  !print *,'vxG=',vxGas,' vyG=',vyGas,' vCI=',vIonC,' vCG=',vGasC,' wCI=',wIonC,' FIonC=',FIonC          
END SUBROUTINE IonGasCollisionF

SUBROUTINE IonGasCollisionWikipediaF(params,mIon,vxIon,vyIon,FxIon,FyIon, &
                            uGas,vGas,sigxGas,sigyGas)
  USE CITSIMparams2D
  IMPLICIT NONE
  DOUBLE PRECISION, PARAMETER :: twopi=6.28318530717959
  DOUBLE PRECISION, INTENT(IN) :: mIon,vxIon,vyIon,params(0:maxparams-1)
  DOUBLE PRECISION, INTENT(IN) :: uGas,vGas,sigxGas,sigyGas
  DOUBLE PRECISION, INTENT(OUT) :: FxIon,FyIon
  DOUBLE PRECISION :: rU(3),phi,radBM,thtBM,z0,z1,vxGas,vyGas,dt,mGas
  DOUBLE PRECISION :: eRestitution,h
  DOUBLE PRECISION :: theta1,theta2,vIonxr,vIonyr,vGasxr,vGasyr,nvIon,nvGas,vIonfxr,vIonfx,vIonfy
  !
  dt = params(3);  mGas = params(18); eRestitution = params(121)
  h= params(21)
  CALL RANDOM_NUMBER(rU)
  phi=ASIN(SQRT(0.999999999*rU(1))) ! Random contact angle, phi is angle along which collision occurs
  ! Box-Muller transform
  radBM = SQRT(-2.*LOG(rU(2))); thtBM = twopi*rU(3)
  z0 = radBM*COS(thtBM); z1 = radBM*SIN(thtBM)
  ! Draw velocity from Maxwell distribution for buffer gas and express in (grid units)/ms      
  vxGas = (sigxGas*z0+uGas); vyGas = (sigyGas*z1+vGas)
  nvIon=SQRT(vxIon**2+vyIon**2)
  theta1=ATAN2(vyIon,vxIon)
  nvGas=SQRT(vxGas**2+vyGas**2)
  theta2=ATAN2(vyGas,vxGas)
  vIonxr=nvIon*COS(theta1-phi)
  vIonyr=nvIon*SIN(theta1-phi)
  vGasxr=nvGas*COS(theta2-phi)
  vGasyr=nvGas*SIN(theta2-phi)
  vIonfxr=(mIon*vIonxr + mGas*vGasxr + mGas*eRestitution*(vGasxr-vIonxr))/(mIon+mGas)
  vIonfx=vIonfxr*COS(phi)-vIonyr*SIN(phi)
  vIonfy=vIonfxr*SIN(phi)+vIonyr*COS(phi)
  ! Find force that applied over time step dt would produce same change of ion velocity
  ! (Note: since velocity was in (grid units)/ms, force is in (amu) (grid units)/ms^2
  FxIon = mIon*(vIonfx-vxIon)/dt
  FyIon = mIon*(vIonfy-vyIon)/dt
  !print *,'vxG=',vxGas,' vyG=',vyGas,' vCI=',vIonC,' vCG=',vGasC,' wCI=',wIonC,' FIonC=',FIonC          
END SUBROUTINE IonGasCollisionWikipediaF

SUBROUTINE StoreCollision(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)
  USE CITSIMparams2D
  IMPLICIT NONE
  INTEGER, INTENT(IN) :: l,nIons,nplt
  DOUBLE PRECISION, INTENT(INOUT) :: Xcol(nIons,nDim,nplt),Vcol(nIons,nDim,nplt),Fcol(nIons,nDim,nplt)
  DOUBLE PRECISION, INTENT(IN) :: Xion(mxIons,nDim),Vion(mxIons,nDim),Fion(mxIons,nDim)
  !
  INTEGER nxt
  !
  nxt = Xcol(l,1,nplt)
  IF (nxt >= nplt) RETURN  ! No more space to store collisions
  Xcol(l,1:nDim,nxt) = Xion(l,1:nDim)
  Vcol(l,1:nDim,nxt) = Vion(l,1:nDim)
  Fcol(l,1:nDim,nxt) = Fion(l,1:nDim)
  Xcol(l,1,nplt) = nxt + 1  
END SUBROUTINE StoreCollision

SUBROUTINE  ComputeIonGasForces2d(params,GasFlow,InTrap,mIon,cIon,csIon,Xion,Vion,nIons,Fion, &
                                  nplt,Xcol,Vcol,Fcol,                                  &
                                  flowData,mxLBM,myLBM,x0LBM,y0LBM,dxLBM,dyLBM,tLBM,dtLBM)
  USE, INTRINSIC :: ISO_C_BINDING
  USE CITSIMparams2D
  IMPLICIT NONE
  !
  LOGICAL, INTENT(IN) :: GasFlow,InTrap(nIons)
  INTEGER, INTENT(IN) :: nIons,nplt
  DOUBLE PRECISION, INTENT(IN) :: params(0:maxparams-1),mIon(nIons),cIon(nIons),csIon(nIons)
  DOUBLE PRECISION, INTENT(IN) :: Xion(mxIons,nDim),Vion(mxIons,nDim)
  DOUBLE PRECISION, INTENT(OUT) :: Fion(mxIons,nDim)
  DOUBLE PRECISION, INTENT(INOUT) :: Xcol(nIons,nDim,nplt),Vcol(nIons,nDim,nplt),Fcol(nIons,nDim,nplt)
  INTEGER, INTENT(IN) :: mxLBM,myLBM
  DOUBLE PRECISION, INTENT(IN) :: x0LBM,y0LBM,dxLBM,dyLBM,tLBM,dtLBM
  REAL (KIND=C_FLOAT), INTENT(IN) :: flowData(mxLBM,myLBM,6)  
  !
  INTEGER, PARAMETER :: Langevin=1, HardSphere=2, LatticeBoltzmann=3,DSMC=4, Hybrid=5
  INTEGER :: CollisionModel,l,iLBM,jLBM
  DOUBLE PRECISION, PARAMETER :: twopi=6.28318530717959,sqrtpi=1.77245385090552,kBamu=8314.47148427469
  DOUBLE PRECISION :: rUnif,h,dt,Pcoll,pres,temp,mGas,mReduced
  DOUBLE PRECISION :: sig,sigma,sigmah,sigma2,sigLBM,sigxGas,sigyGas,uGas,vGas
  DOUBLE PRECISION :: AvgGasVel,MedGasVel,AvgRelVel,s,VionMag,CF,Pfactor,dist
  DOUBLE PRECISION :: Xmin(nDim),XLBM(nDim),RefLBMVel,LBMVel
  DOUBLE PRECISION :: AvgGasVelThermal,LBMVel2,MedGasVelThermal,invlambda
  !
  Fion = 0.d0
  CollisionModel = params(20); h= params(21); dt = params(3); mGas = params(18)
  pres = params(101); temp = params(111)
  Xmin(1) = params(22); Xmin(2) = params(23)
  sigma = params(28); sigmah=sigma/h
  IF ((mGas<=0) .OR. (pres<=0) .OR. (temp<=0)) RETURN
  SELECT CASE (CollisionModel) 
    CASE (Langevin)
      Pfactor = params(26)
      DO l=1,nIons
        IF (.NOT. InTrap(l)) CYCLE
        mReduced = mIon(l)*mGas/(mIon(l)+mGas)
        CALL RANDOM_NUMBER(rUnif)  ! Get uniform distributed random number
        ! Compute collision probability
        Pcoll = Pfactor*cIon(l)/SQRT(mReduced)
        IF (Pcoll>rUnif) THEN
          CALL IonGasCollisionF(params,mIon(l),Vion(l,1),Vion(l,2),Fion(l,1),Fion(l,2),0.d0,0.d0,sigmah,sigmah)
          Xcol(l,1,nplt) = Xcol(l,1,nplt) + 1
          !CALL StoreCollision(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)        
        END IF
      END DO
    CASE (HardSphere)
      ! Extract buffer gas conditions
      AvgGasVel = params(30)
      MedGasVel = params(31)
      DO l=1,nIons
        IF (.NOT. InTrap(l)) CYCLE
        VionMag = SQRT(Vion(l,1)**2 + Vion(l,2)**2)
        s = VionMag/MedGasVel
        IF (s<1.0e-3) THEN
          AvgRelVel = AvgGasVel
        ELSE
          AvgRelVel = 0.5*AvgGasVel*((s+0.5/s)*sqrtpi*ERF(s) + EXP(-s**2))   
        END IF
        ! params(29)=dt * kB*tempGas[0]/(pi*diamHS**2)/(presGas[0]*101325./760.) * 1.e9/h
        invlambda=csIon(l)*params(29)/params(122) !invlambda to replace param(29)
        CF = 1. - EXP(-invlambda*AvgRelVel)
        CALL RANDOM_NUMBER(rUnif)
        IF (rUnif < CF) THEN          
          CALL IonGasCollisionWikipediaF(params,mIon(l),Vion(l,1),Vion(l,2),Fion(l,1),Fion(l,2),0.d0,0.d0,sigmah,sigmah)
          Xcol(l,1,nplt) = Xcol(l,1,nplt) + 1
          !CALL StoreCollision(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)        
        END IF
      END DO
    CASE (Hybrid)
      ! Pressure and temperature defined in presGas[0] and tempGas[0]
      AvgGasVel = params(30)
      MedGasVel = params(31)
      DO l=1,nIons
        IF (.NOT. InTrap(l)) CYCLE
        mReduced = mIon(l)*mGas/(mIon(l)+mGas)
        VionMag = SQRT(Vion(l,1)**2 + Vion(l,2)**2)
        s = VionMag/MedGasVel
        IF (s<1.0e-3) THEN
          AvgRelVel = AvgGasVel
        ELSE
          AvgRelVel = 0.5*AvgGasVel*((s+0.5/s)*sqrtpi*ERF(s) + EXP(-s**2))   
        END IF
        ! params(29)=dt * kB*tempGas[0]/(pi*diamHS**2)/(presGas[0]*101325./760.) * 1.e9/h
        invlambda=csIon(l)*params(29)/params(122) !invlambda to replace param(29)
        CF = 1. - EXP(-invlambda*AvgRelVel-params(26)*cIon(l)/SQRT(mReduced))
        CALL RANDOM_NUMBER(rUnif)
        IF (rUnif < CF) THEN          
          CALL IonGasCollisionWikipediaF(params,mIon(l),Vion(l,1),Vion(l,2),Fion(l,1),Fion(l,2),0.d0,0.d0,sigmah,sigmah)
          Xcol(l,1,nplt) = Xcol(l,1,nplt) + 1
          !CALL StoreCollision(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)        
        END IF
      END DO
    CASE (LatticeBoltzmann)
      ! Extract buffer gas conditions
      AvgGasVelThermal = params(30)
      MedGasVelThermal = params(31)      
      DO l=1,nIons
        IF (.NOT. InTrap(l)) CYCLE   
        VionMag = SQRT(Vion(l,1)**2 + Vion(l,2)**2)
        ! Interpolate LBM data to current ion position
        !
        ! Ion position in LBM reference frame, units (mm)
        XLBM(1:nDim) = h*Xion(l,1:nDim) + Xmin(1:nDim)
        iLBM = (XLBM(1)-x0LBM)/dxLBM; iLBM = MAX(1,iLBM); iLBM = MIN(iLBM,mxLBM)
        jLBM = (XLBM(2)-y0LBM)/dyLBM; jLBM = MAX(1,jLBM); jLBM = MIN(jLBM,myLBM)
        !print *,'iLBM,jLBM=',iLBM,jLBM
        LBMVel2 = (flowData(iLBM,jLBM,2)**2 + flowData(iLBM,jLBM,3)**2)/h**2
        LBMVel = SQRT(LBMVel2)
        AvgGasVel = AvgGasVelThermal + LBMVel
        MedGasVel = MedGasVelThermal + LBMVel
        s = VionMag/MedGasVel
        IF (s<1.0e-3) THEN
          AvgRelVel = AvgGasVel
        ELSE
          AvgRelVel = 0.5*AvgGasVel*((s+0.5/s)*sqrtpi*ERF(s) + EXP(-s**2))   
        END IF
        ! params(29)=dt * kB*tempGas[0]/(pi*diamHS**2)/(presGas[0]*101325./760.) * 1.e9/h 
        ! params(101)=presGas[0]
        pres = flowdata(iLBM,jLBM,1)*760./101325.
        invlambda=csIon(l)*params(29)/params(122) !invlambda to replace param(29)
        CF = 1. - EXP(-invlambda/params(101)*pres * AvgRelVel)
        CALL RANDOM_NUMBER(rUnif)
        IF (rUnif < CF) THEN
          sigLBM = SQRT(flowData(iLBM,jLBM,4) + flowData(iLBM,jLBM,5))
          ! Following are expressed as ratios and will be rescaled to sigmah in IonGasCollisionF
          uGas = flowData(iLBM,jLBM,2)/h
          vGas = flowData(iLBM,jLBM,3)/h
          sigma = params(28)
          sigxGas = sigma*flowData(iLBM,jLBM,4)/h
          sigyGas = sigma*flowData(iLBM,jLBM,5)/h
          CALL IonGasCollisionF(params,mIon(l),Vion(l,1),Vion(l,2),Fion(l,1),Fion(l,2), &
                                uGas,vGas,sigxGas,sigyGas)
          CALL StoreCollision(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)
        END IF
      END DO
    CASE (DSMC)
      ! Extract buffer gas conditions
      DO l=1,nIons
        IF (.NOT. InTrap(l)) CYCLE   
        VionMag = SQRT(Vion(l,1)**2 + Vion(l,2)**2)
        ! Interpolate LBM data to current ion position
        !
        ! Ion position in LBM reference frame, units (mm)
        XLBM(1:nDim) = h*Xion(l,1:nDim) + Xmin(1:nDim)
        iLBM = (XLBM(1)-x0LBM)/dxLBM; iLBM = MAX(1,iLBM); iLBM = MIN(iLBM,mxLBM)
        jLBM = (XLBM(2)-y0LBM)/dyLBM; jLBM = MAX(1,jLBM); jLBM = MIN(jLBM,myLBM)
        !print *,'iLBM,jLBM=',iLBM,jLBM
        AvgGasVelThermal = sqrt(8/3.141592653589793*kBamu*flowData(iLBM,jLBM,6)/mGas)/h ! Average (mean) buffer gas speed
        MedGasVelThermal = sqrt(2 * kBamu*flowData(iLBM,jLBM,6)/mGas)/h    ! Median buffer gas speed
        LBMVel = SQRT((flowData(iLBM,jLBM,2)**2 + flowData(iLBM,jLBM,3)**2)/h**2)
        AvgGasVel = AvgGasVelThermal + LBMVel
        MedGasVel = MedGasVelThermal + LBMVel
        s = VionMag/MedGasVel
        IF (s<1.0e-3) THEN
          AvgRelVel = AvgGasVel
        ELSE
          AvgRelVel = 0.5*AvgGasVel*((s+0.5/s)*sqrtpi*ERF(s) + EXP(-s**2))   
        END IF
        ! params(29)=dt * kB*tempGas[0]/(pi*diamHS**2)/(presGas[0]*101325./760.) * 1.e9/h 
        ! params(101)=presGas[0]
        pres = flowdata(iLBM,jLBM,1)*760./101325. !pres in Torr
        invlambda=csIon(l)*params(29)/params(122) !invlambda to replace param(29)
        CF = 1. - EXP(-invlambda/params(101)*pres * AvgRelVel)
        CALL RANDOM_NUMBER(rUnif)
        IF (rUnif < CF) THEN
          ! Following are expressed as ratios and will be rescaled to sigmah in IonGasCollisionF
          uGas = flowData(iLBM,jLBM,2)/h !non dim. units (lattice units)
          vGas = flowData(iLBM,jLBM,3)/h !non dim. units (lattice units)
          sigxGas = flowData(iLBM,jLBM,4)/h
          sigyGas = flowData(iLBM,jLBM,5)/h
          CALL IonGasCollisionF(params,mIon(l),Vion(l,1),Vion(l,2),Fion(l,1),Fion(l,2), &
                                uGas,vGas,sigxGas,sigyGas)
          CALL StoreCollision(l,nIons,nplt,Xcol,Vcol,Fcol,Xion,Vion,Fion)
        END IF
      END DO
  END SELECT  
END SUBROUTINE ComputeIonGasForces2d

SUBROUTINE trajectory2d(params,mIon,cIon,csIon,IonExit,nXplt,tplt,Xplt,Vplt,Fplt,Ncolplt,&
                        Xcol,Vcol,Fcol,aGrid,flowBC,nIons,nplt,nx,ny,nB)
  USE, INTRINSIC :: ISO_C_BINDING 
  USE CITSIMparams2D
  IMPLICIT NONE
  ! Interface declarations
  INTEGER, INTENT(IN) :: nIons,nx,ny,nplt,nB
  INTEGER, INTENT(INOUT) :: nXplt(nIons)
  DOUBLE PRECISION, INTENT(INOUT) :: params(0:maxparams-1)
  DOUBLE PRECISION, INTENT(INOUT) :: tplt(nplt),Xplt(nIons,nDim,nplt),Vplt(nIons,nDim,nplt),Fplt(nIons,nDim,nplt)
  DOUBLE PRECISION, INTENT(INOUT) :: Xcol(nIons,nDim,nplt),Vcol(nIons,nDim,nplt),Fcol(nIons,nDim,nplt),Ncolplt(nIons,nplt)
  DOUBLE PRECISION, INTENT(IN) :: mIon(nIons),cIon(nIons),csIon(nIons),aGrid(0:nx-1,0:ny-1,nDim,nElec)
  DOUBLE PRECISION, INTENT(IN) :: flowBC(nB,nDim+3)
  DOUBLE PRECISION, INTENT(INOUT) :: IonExit(nIons,0:4*nDim)
  ! Local variables
  DOUBLE PRECISION, PARAMETER :: twopi=6.283185307179586
  CHARACTER :: fname*16,line*80
  LOGICAL InTrap(nIons),Active(nIons),IonIon,IonGas,GasFlow
  INTEGER, PARAMETER :: Langevin=1, HardSphere=2, LatticeBoltzmann=3, DSMC=4, Hybrid=5
  INTEGER :: i,j,k,n,l,inxtplt,iplt,mLBM(nDim),iError,CollisionModel,tBirth,its
  INTEGER (C_INT) :: nActive
  REAL (C_FLOAT), ALLOCATABLE, DIMENSION(:,:,:) :: fGas 
  DOUBLE PRECISION :: dt2,tini,trmp,tfin,tnow,dt,U0,dUdt,Uent,Uext,Omega
  DOUBLE PRECISION :: tnxtplt,dtplt,RFperiod,dphi,phi,cphi,sphi
  DOUBLE PRECISION :: X(mxIons,nDim,nstages),V(mxIons,nDim,nstages)
  DOUBLE PRECISION :: XF(mxIons,nDim),VF(mxIons,nDim),F(mxIons,nDim)
  DOUBLE PRECISION :: Fionion(mxIons,nDim), Fiongas(mxIons,nDim)
  DOUBLE PRECISION :: rExtentMin,rExtentMax,zExtentMin,zExtentMax,hGrid,Xmin(nDim)
  DOUBLE PRECISION :: AvgKESum(nIons),AvgDistSum(nIons),AvgDTCSum(nIons),AvgKEplt(nIons),AvgDistplt(nIons),AvgDTCplt(nIons)
  DOUBLE PRECISION :: d2(nDim),tBirths(nIons),dtBirth,test(2)
  !
  INTEGER :: mxLBM,myLBM,nCols
  DOUBLE PRECISION :: x0LBM,y0LBM,dxLBM,dyLBM,tLBM,dtLBM
  REAL (C_FLOAT), ALLOCATABLE, DIMENSION(:,:,:) :: flowData
  print *,"at ion_traj, start of trajectory2d subroutine, nDim=",nDim
  
  ! Return trajectories for n ions at nplt time slices in nDim dimensions
  
  ! Assume all ions are in trap initially
  InTrap = .TRUE.; Active = .FALSE.
  
  ! Load parameters:
  ! times
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
  ! simulation options
  IonIon  = params(15) > 0
  IonGas  = params(20) > 0
  GasFlow = params(19) > 0
  CollisionModel = params(20)
  ! geometry parameters
  hGrid = params(21); Xmin(1) = params(22); Xmin(2) = params(23)
  print *,"Xmin(:)=",Xmin(:)
  rExtentMax = params(50)/hGrid; rExtentMin = -rExtentMax
  !yExtentMax = rExtentMax; xExtentMax = rExtentMax; !!rExtentMax and rExtentMin can be used directly 
  !yExtentMin = rExtentMin; xExtentMin = rExtentMin;
  zExtentMax = (params(51) - Xmin(1))/hGrid; zExtentMin = (-params(51) - Xmin(1))/hGrid  
  print *,"rExtentMax,zExtentMax,rExtentMin,zExtentMin=",rExtentMax,zExtentMax,rExtentMin,zExtentMin
  
  ! dtBirth calculation and tBirths generation
  IF (tBirth==0) CALL RANDOM_NUMBER(tBirths) ! Generation of random times of births (bruno)
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
  IF (CollisionModel .EQ. DSMC) THEN
    PRINT *,'Reading DSMC flow data ...'
    OPEN(UNIT=12,FILE='../flow/DSMCflow.data',STATUS='OLD')
    READ(12,*)mxLBM,myLBM,x0LBM,y0LBM,dxLBM,dyLBM,tLBM,dtLBM !I have to keep these variables for consistency
    PRINT *,' Flow computation in z-r plane on ',mxLBM,' by ',myLBM,' lattice'
    ALLOCATE(flowData(mxLBM,myLBM,6),STAT=iError)
    IF (iError /= 0) THEN
      PRINT *,'Cannot allocate space for lattice Boltzmann flow computation'
      STOP
    END IF 
    DO i=1,mxLBM
      DO j=1,myLBM
        ! DSMCflowData contains p u v (SI units Pa, m/s, m/s), sigx sigy (normalized unit vector components, non-dimensional),and T
        READ(12,*)flowData(i,j,1:6)
      END DO
    END DO
    PRINT *,' ... done'
    CLOSE(12)
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

  !Main loop
  DO WHILE (tnow<tfin)
    its=its+1  !time steps counter (bruno)
    
    !Active ions are created randomly within the first RF cycle or at regular dtBirths/nIons (bruno) 
    IF (tnow<RFperiod) THEN
      SELECT CASE (tBirth)
        CASE (0) ! random times of birth
        WHERE (tnow .GT. tBirths) Active = .TRUE.
        CASE (1) ! regular times of birth
        nActive = MIN( CEILING((tnow-tini)/dtBirth)*nIons , nIons)
        Active(1:nActive) = .TRUE.
        CASE (2) ! fixed times of birth
        IF (tnow>params(53)) THEN
          Active(1:nIons) = .TRUE.
        END IF
      END SELECT
    END IF  
    
    XF(1:nIons,1:nDim) = X(1:nIons,1:nDim,2)  ! X for forces at midpoint
    
    ! Determine if ions still inside the trap
    DO l=1,nIons
      
      ! Conditions to determine if ions are still within the trap 
      IF (InTrap(l) .AND. ((XF(l,1)<0) .OR. (XF(l,1)>nx) .OR. &
                           (XF(l,2)<-ny) .OR. (XF(l,2)>ny) .OR. &
                           ((XF(l,1)>zExtentMax).AND.((XF(l,2)>rExtentMax).OR.(XF(l,2)<rExtentMin))) .OR. &
                             ((XF(l,1)<zExtentMin).AND.((XF(l,2)>rExtentMax).OR.(XF(l,2)<rExtentMin))))) THEN
        InTrap(l) = .FALSE.
        IonExit(l,0) = tnow
        IonExit(l,1:nDim) = hGrid*X(l,1:nDim,2) + Xmin(1:nDim)
        IonExit(l,nDim+1:2*nDim) = (X(l,1:nDim,2)-X(l,1:nDim,1))/dt*hGrid
        IonExit(l,2*nDim+1:3*nDim) = hGrid*X(l,1:nDim,1) + Xmin(1:nDim)
        PRINT 1001,l,mIon(l),cIon(l),csIon(l)*1e13,IonExit(l,1),IonExit(l,2),tnow,inxtplt+1    !added IonExit(l,3) for 3D
        1001 FORMAT('Ion ',i5,' of mass ',F8.1,' amu, charge ', F4.1,'+ and ccs ',E10.4,' mm^2 has moved out of trap area at x='&
        ,F8.4,'(mm) y=',F8.4,'(mm) t=',F10.5,'(ms) iplt=',i7)
      ELSE IF (InTrap(l)) THEN
        ! Calculate avg KE, avg flight distance, avg dist. to center (DTC) (bruno)
        d2=(X(l,1:nDim,2)-X(l,1:nDim,1))**2
        AvgKESum(l)=AvgKESum(l)+d2(1)+d2(2)
        AvgDistSum(l)=AvgDistSum(l)+SQRT(d2(1)+d2(2))
        AvgDTCSum(l)=AvgDTCSum(l)+SQRT((X(l,1,2)+Xmin(1)/hGrid)**2+(X(l,2,2)-Xmin(2)/hgrid)**2)
      END IF
    END DO
    
    !  Electrode force
    CALL ComputeElectrodeForces2d(tnow,params,InTrap,aGrid,XF,F,nIons,nx,ny,cIon)
    
    !  Ion-ion interaction force
    IF (IonIon) THEN
      CALL ComputeIonIonForces2d(params,InTrap,cIon,XF,Fionion,nIons)
      F(1:nIons,:) = F(1:nIons,:) + Fionion(1:nIons,:)
    END IF
    
    !  Ion-buffer interaction force
    !  At this stage X(2) is known and V(1.5)=(X(2)-X(1))/dt is O(dt^2) accurate
    !  Collision is assumed to occur from t1.5 to t2.5 and modeled as a force at t2
    IF (IonGas) THEN
      DO l=1,nIons
        VF(l,:) = (X(l,:,2) - X(l,:,1))/dt
      END DO
      CALL ComputeIonGasForces2d(params,GasFlow,InTrap,mIon,cIon,csIon,XF,VF,nIons,Fiongas, &
                                 nplt,Xcol,Vcol,Fcol,                                 &
                                 flowData,mxLBM,myLBM,x0LBM,y0LBM,dxLBM,dyLBM,tLBM,dtLBM)
      F(1:nIons,:) = F(1:nIons,:) + Fiongas(1:nIons,:)
    END IF
    
    ! Verlet update
    DO l=1,nIons
        IF (.NOT. InTrap(l)) CYCLE
        IF (.NOT. Active(l)) CYCLE
        X(l,:,3) = 2.d0*X(l,:,2) - X(l,:,1) + dt2*F(l,:)/mIon(l)
        V(l,:,2) = (X(l,:,3) - X(l,:,1))/(2.d0*dt)
    END DO
    
    ! Store position at this time if requested (change to maintain latest time steps in last 20% of nplt)
    IF (tnow>=tnxtplt) THEN
      IF (inxtplt>nplt) PRINT *,'Warning current plot = ',inxtplt,' > nplt=',nplt,' (overwriting last plot)'
      iplt=MIN(inxtplt,nplt)
      tplt(iplt) = tnow
      DO l=1,nIons
        IF (InTrap(l)) THEN      
          Xplt(l,1:nDim,iplt) = hGrid*X(l,1:nDim,2) + Xmin(1:nDim) ! Scaled back from grid coordinates to physical coordinates
          Vplt(l,1:nDim,iplt) = hGrid*V(l,1:nDim,2)
          Fplt(l,1:nDim,iplt) = F(l,1:nDim)
          Ncolplt(l,iplt)=Xcol(l,1,nplt)-1 !to save the nb of collisions along a traj
          nXplt(l) = iplt
        END IF
      END DO
      tnxtplt=tnxtplt+dtplt; inxtplt=inxtplt+1
    END IF
    
    ! Increment time (moved the time increment after the Xplt update, Bruno)
    ! in this way we store position and velocity of euler time step as well, and be consistent for the rest tof the time steps.
    tnow=tnow+dt
    
    ! Update positions and velocities 
    DO l=1,nIons
      IF (.NOT. InTrap(l)) CYCLE
      IF (.NOT. Active(l)) CYCLE
      X(l,1:nDim,1)=X(l,1:nDim,2); X(l,1:nDim,2)=X(l,1:nDim,3)
      V(l,1:nDim,1)=V(l,1:nDim,2)
    END DO
    
  END DO ! End of Main Loop
  
  ! Use params to transfer calculated means (bruno)
!  AvgKEplt=0.
!  AvgDistplt=0.
!  AvgDTCplt=0.
  params(9) = iplt
  params(149) = its !(bruno)
  IF (nIons<116) THEN ! can only be transfered with params if nIons is small
  params(150:150+3*nIons-1)=0
    DO l=1,nIons !1.66053892e-27kg/amu / 1.60217657e-19J/eV (=) 1.03642692e-8 = factor to have KE in eV
      params(150:150+nIons-1)=AvgKESum(:)*(hGrid/dt)**2*0.5*1.03642692e-8*mIon(l)/its  
      params(150+nIons:150+2*nIons-1)=AvgDistSum(:)*hGrid
      params(150+2*nIons:150+3*nIons-1)=AvgDTCSum(:)*hGrid/its 
      !DO i=1,iplt-1
      !  AvgKEplt(l)=AvgKEplt(l)+0.5*mIon(l)*1.03642692e-8*(Vplt(l,1,i)**2+Vplt(l,2,i)**2)
      !  AvgDistplt(l)=AvgDistplt(l)+SQRT((Xplt(l,1,i+1)-Xplt(l,1,i))**2+(Xplt(l,2,i+1)-Xplt(l,2,i))**2)
      !  AvgDTCplt(l)=AvgDTCplt(l)+SQRT((Xplt(l,1,i))**2+(Xplt(l,2,i))**2)
      !END DO
    END DO  
  END IF
  
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
END SUBROUTINE trajectory2d
