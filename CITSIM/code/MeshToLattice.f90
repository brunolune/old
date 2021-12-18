! VmeshToGrid.f90 - Fortran routines to compute Cartesian grid potential and
!                   geometry values from unstructured mesh results

DOUBLE PRECISION FUNCTION TrngArea(P)
  DOUBLE PRECISION :: P(3,2),A,X(3),Y(3)
  X(:) = P(:,1); Y(:) = P(:,2)
  A = X(2)*Y(3)-X(3)*Y(2) + &
      X(3)*Y(1)-X(1)*Y(3) + &
      X(1)*Y(2)-X(2)*Y(1)
END FUNCTION TrngArea

SUBROUTINE mesh_to_grid_2d(mX,NrNodes,ITnodes,NrElems,ITelems,nX1,nX2, &
                           Vgrid,Ggrid,Fgrid)
! Python call: mesh_to_grid_2d(h,Xmin,mX,ITnodes,ITelems,Vgrid,Ggrid,Fgrid)
  IMPLICIT DOUBLE PRECISION (A-H,O-Z)
  INTEGER, PARAMETER :: nD=2,nN=3,nE=3 ! Nr. of dims, nodes per element, electrodes
  !
  INTEGER, INTENT(IN) :: mX(nD),NrNodes,NrElems,nX1,nX2
  DOUBLE PRECISION, INTENT(IN) :: ITnodes(NrNodes,nD+1)
  INTEGER, INTENT(IN) :: ITelems(NrElems,nN+1)
  DOUBLE PRECISION, INTENT(INOUT) :: Vgrid(nX1,nX2,3)
  INTEGER, INTENT(INOUT) :: Ggrid(nX1,nX2)
  DOUBLE PRECISION, INTENT(INOUT) :: Fgrid(nX1,nX2,3,3)
  !
  INTEGER :: iTmin(nD),iTmax(nD)
  DOUBLE PRECISION :: TrngPts(nN,nD)
  ! Grid has unit cell spacing and extends from (0,0) to (mX(1),mX(2))
  ! Indices: i,j,k - x,y,z; n - nodes; l - elements
  DO l=1,NrElems                       ! Loop over all triangles
    DO m=1,nN                          ! Get triangle vertices
      n = ITelems(l,m)
      TrngPts(m,1:nD) = ITnodes(n,1:nD)
    END DO
    TrngA = TrngArea(TrngPts)
    ! Enclosing index rectangle
    DO m=1,nD
      iTmin(m) =   FLOOR(MINVAL(TrngPts(:,m)))
      iTmax(m) = CEILING(MAXVAL(TrngPts(:,m)))
    END DO
    ! Loop over grid points in enclosing triangle
    DO i=iTmin(1),iTmax(1)
      DO j=iTmin(2),iTmax(2)
        
      END DO
    END DO
  END DO
END SUBROUTINE mesh_to_grid_2d
