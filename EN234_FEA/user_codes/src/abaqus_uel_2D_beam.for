!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 2D linear elastic continuum element
!
!    The file also needs the following subroutines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_1D_integrationpoints           - defines integration points for 1D line integral
!          abq_facenodes_2D                       - returns list of nodes on an element face
!          abq_inverse_LU                         - computes the inverse of an arbitrary matrix by LU decomposition
!          abq_UEL_invert2d                       - computes inverse and determinant of a 2x2 matrix
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL_BEAM(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
     1     PROPS,NPROPS,COORDS,MCRD,NNODE,U,DU,V,A,JTYPE,TIME,DTIME,
     2     KSTEP,KINC,JELEM,PARAMS,NDLOAD,JDLTYP,ADLMAG,PREDEF,NPREDF,
     3     LFLAGS,MLVARX,DDLMAG,MDLOAD,PNEWDT,JPROPS,NJPROP,PERIOD)
    !
      INCLUDE 'ABA_PARAM.INC'
    !
    !
      DIMENSION RHS(MLVARX,*),AMATRX(NDOFEL,NDOFEL),PROPS(*),
     1   SVARS(*),ENERGY(8),COORDS(MCRD,NNODE),U(NDOFEL),
     2   DU(MLVARX,*),V(NDOFEL),A(NDOFEL),TIME(2),PARAMS(*),
     3   JDLTYP(MDLOAD,*),ADLMAG(MDLOAD,*),DDLMAG(MDLOAD,*),
     4   PREDEF(2,NPREDF,NNODE),LFLAGS(*),JPROPS(*)

    !
    !       Variables that must be computed in this routine
    !       RHS(i)                     Right hand side vector.  In EN234_FEA the dimensions are always RHS(MLVARX,1)
    !       AMATRX(i,j)                Stiffness matrix d RHS(i)/ d DU(j)
    !       SVARS(1:NSVARS)            Element state variables.  Must be updated in this routine
    !       ENERGY(1:8)
    !                                  Energy(1) Kinetic Energy
    !                                  Energy(2) Elastic Strain Energy
    !                                  Energy(3) Creep Dissipation
    !                                  Energy(4) Plastic Dissipation
    !                                  Energy(5) Viscous Dissipation
    !                                  Energy(6) Artificial strain energy
    !                                  Energy(7) Electrostatic energy
    !                                  Energy(8) Incremental work done by loads applied to the element
    !       PNEWDT                     Allows user to control ABAQUS time increments.
    !                                  If PNEWDT<1 then time step is abandoned and computation is restarted with
    !                                  a time increment equal to PNEWDT*DTIME
    !                                  If PNEWDT>1 ABAQUS may increase the time increment by a factor PNEWDT
    !
    !       Variables provided for information
    !       NDOFEL                     Total # DOF for the element
    !       NRHS                       Dimension variable
    !       NSVARS                     Total # element state variables
    !       PROPS(1:NPROPS)            User-specified properties of the element
    !       NPROPS                     No. properties
    !       JPROPS(1:NJPROPS)          Integer valued user specified properties for the element
    !       NJPROPS                    No. integer valued properties
    !       COORDS(i,N)                ith coordinate of Nth node on element
    !       MCRD                       Maximum of (# coords,minimum of (3,#DOF)) on any node
    !       U                          Vector of DOF at the end of the increment
    !       DU                         Vector of DOF increments
    !       V                          Vector of velocities (defined only for implicit dynamics)
    !       A                          Vector of accelerations (defined only for implicit dynamics)
    !       JTYPE                      Integer identifying element type (the number n in the Un specification in the input file)
    !       TIME(1:2)                  TIME(1)   Current value of step time
    !                                  TIME(2)   Total time
    !       DTIME                      Time increment
    !       KSTEP                      Current step number (always 1 in EN234_FEA)
    !       KINC                       Increment number
    !       JELEM                      User assigned element number in ABAQUS (internally assigned in EN234_FEA)
    !       PARAMS(1:3)                Time increment parameters alpha, beta, gamma for implicit dynamics
    !       NDLOAD                     Number of user-defined distributed loads defined for this element
    !       JDLTYP(1:NDLOAD)           Integers n defining distributed load types defined as Un or (if negative) UnNU in input file
    !       ADLMAG(1:NDLOAD)           Distributed load magnitudes
    !       DDLMAG(1:NDLOAD)           Increment in distributed load magnitudes
    !       PREDEF(1:2,1:NPREDF,1:NNODE)   Predefined fields.
    !       PREDEF(1,...)              Value of predefined field
    !       PREDEF(2,...)              Increment in predefined field
    !       PREDEF(1:2,1,k)            Value of temperature/temperature increment at kth node
    !       PREDEF(1:2,2:NPREDF,k)     Value of user defined field/field increment at kth node (not used in EN234FEA)
    !       NPREDF                     Number of predefined fields (1 for en234FEA)
    !       LFLAGS                     Control variable
    !       LFLAGS(1)                  Defines procedure type
    !       LFLAGS(2)                  0 => small displacement analysis  1 => Large displacement (NLGEOM option)
    !       LFLAGS(3)                   1 => Subroutine must return both RHS and AMATRX (always true in EN234FEA)
    !                                   2 => Subroutine must return stiffness AMATRX = -dF/du
    !                                   3 => Subroutine must return daming matrix AMATRX = -dF/dudot
    !                                   4 => Subroutine must return mass matrix AMATRX = -dF/duddot
    !                                   5 => Define the RHS only
    !                                   6 => Define the mass matrix for the initial acceleration calculation
    !                                   100 => Define perturbation quantities for output
    !       LFLAGS(4)                   0 => General step   1 => linear perturbation step
    !       LFLAGS(5)                   0 => current approximation to solution based on Newton correction; 1 => based on extrapolation
    !       MLVARX                      Dimension variable (equal to NDOFEL in EN234FEA)
    !       PERIOD                      Time period of the current step
    !
    !
    ! Local Variables
      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ksize,NNODEs
      integer      :: face_node_list(3)                       ! List of nodes on an element face
    !
      double precision  ::  xi(2,5)                           ! Area integration points
      double precision  ::  w(5)                              ! Area integration weights
      double precision  ::  N(5)                              ! 2D shape functions
      double precision  ::  dNdxi(9,2)                        ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                         ! Spatial derivatives
      double precision  ::  dxdxi(2,2)                        ! Derivative of spatial coords wrt normalized coords
    
    !   Variables below are for computing integrals over element faces
      double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
      double precision  ::  xi1(5)                            ! 1D integration points
      double precision  ::  w1(5)                             ! Integration weights
      double precision  ::  N1(3)                             ! 1D shape functions
      double precision  ::  dN1dxi(3)                         ! 1D shape function derivatives
      double precision  ::  norm(2)                           ! Normal to an element face
      double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
    !
      double precision  ::  scoords(2,4)                      ! Coordinates of slave nodes
      double precision  ::  T(8,6)                            ! Relates dofs of master, slave nodes
      double precision  ::  elam1(2), elam2(2)                ! Laminar basis vector components
      double precision  ::  R(2,3)                            ! Relates strain components in global, laminar bases
      double precision  ::  RBT(2,6)                          ! For simplification
      double precision  ::  Tension, Shear, Moment            ! Beam forces to be saved as SVARS
      
      double precision  ::  strain(2)                         ! Strain vector contains [e11, e12]
      double precision  ::  stress(2)                         ! Stress vector contains [s11, s12]
      double precision  ::  D(2,2)                            ! stress = D*(strain)  (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(3,8)                            ! strain = B*(dof_total)
      double precision  ::  dxidx(2,2), determinant, det0     ! Jacobian inverse and determinant
      double precision  ::  costi, sinti                      ! Beam rotation components (same for 1,2)
      double precision  ::  h, wi, L                          ! Beam element geometry
      double precision  ::  E, xnu, G                         ! Material properties

    ! ABAQUS UEL implementing a 2D, small strain, linear elastic, beam element     
      
      NNODEs = 4

      !Establish 1D integration sceheme in xi2 - Assume weights for trapezoidal integration
      n_points = 5
      xi = 0.d0
      xi(2,1) = 1.d0
      xi(2,2) = 0.5d0
      xi(2,4) = -0.5d0
      xi(2,5) = -1.d0
      w(1) = 0.25d0
      w(2) = 0.5d0
      w(3) = 0.5d0
      w(4) = 0.5d0
      w(5) = 0.25d0
 
      if (MLVARX<2*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 2*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif

      RHS(1:MLVARX,1)           = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      Tension                   = 0.d0
      Shear                     = 0.d0
      Moment                    = 0.d0

      ! Extract element properties
      h   = PROPS(1)
      wi  = PROPS(2)
      E   = PROPS(3)
      xnu = PROPS(4)
      
      ! Define constutive matrix
      G = E/(2.d0*(1.d0+xnu))
      D = 0.d0
      D(1,1) = E
      D(2,2) = G
      
      ! Define internal element geometry
      L = sqrt((COORDS(1,2)-COORDS(1,1))**2.d0 + 
     1                                (COORDS(2,2)-COORDS(2,1))**2.d0)
      
      sinti = (COORDS(2,2)-COORDS(2,1))/L
      costi = (COORDS(1,2)-COORDS(1,1))/L
      
      scoords = 0.d0
      scoords(1,1) = COORDS(1,1) + (h/2.d0)*sinti !x11
      scoords(1,2) = COORDS(1,2) + (h/2.d0)*sinti !x12
      scoords(1,3) = COORDS(1,2) - (h/2.d0)*sinti !x13
      scoords(1,4) = COORDS(1,1) - (h/2.d0)*sinti !x14
      scoords(2,1) = COORDS(2,1) - (h/2.d0)*costi !x21
      scoords(2,2) = COORDS(2,2) - (h/2.d0)*costi !x22
      scoords(2,3) = COORDS(2,2) + (h/2.d0)*costi !x23
      scoords(2,4) = COORDS(2,1) + (h/2.d0)*costi !x24
      
      T = 0.d0
      T(1,3) =   COORDS(2,1) - scoords(2,1)
      T(2,3) = -(COORDS(1,1) - scoords(1,1))
      T(3,6) =   COORDS(2,2) - scoords(2,2)
      T(4,6) =   COORDS(1,2) - scoords(1,2)
      T(5,6) =   COORDS(2,2) - scoords(2,3)
      T(6,6) =   COORDS(1,2) - scoords(1,3)
      T(7,3) =   COORDS(2,1) - scoords(2,4)
      T(8,3) = -(COORDS(1,1) - scoords(1,4))
      T(1,1) = 1.d0
      T(2,2) = 1.d0
      T(7,1) = 1.d0
      T(8,2) = 1.d0
      T(3,4) = 1.d0
      T(4,5) = 1.d0
      T(5,4) = 1.d0
      T(6,5) = 1.d0

      ENERGY(1:8) = 0.d0

    ! Loop over 1D trapezoidal integration points
      do kint = 1, n_points
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODEs,N,dNdxi)
        dxdxi = matmul(scoords(1:2,1:NNODEs),dNdxi(1:NNODEs,1:2))
        
        call abq_UEL_invert2d(dxdxi,dxidx,determinant)
        dNdx(1:NNODEs,1:2) = matmul(dNdxi(1:NNODEs,1:2),dxidx)
        
        B = 0.d0
        B(1,1:2*NNODEs-1:2) = dNdx(1:NNODEs,1)
        B(2,2:2*NNODEs:2)   = dNdx(1:NNODEs,2)
        B(3,1:2*NNODEs-1:2) = dNdx(1:NNODEs,2)
        B(3,2:2*NNODEs:2)   = dNdx(1:NNODEs,1)
        
        elam1 =(1.d0/sqrt(dxdxi(1,1)**2.d0+dxdxi(2,1)**2.d0))*dxdxi(:,1)
        elam2 =(1.d0/sqrt(dxdxi(1,2)**2.d0+dxdxi(2,2)**2.d0))*dxdxi(:,2)
        
        R = 0.d0
        R(1,1) = elam1(1)**2.d0
        R(1,2) = elam1(2)**2.d0
        R(1,3) = elam1(1)*elam1(2)
        R(2,1) = 2.d0*elam1(1)*elam2(1)
        R(2,2) = 2.d0*elam1(2)*elam2(2)
        R(2,3) = elam1(1)*elam2(2) + elam1(2)*elam2(1)
        
        RBT = matmul(R,matmul(B,T))
        
        strain = matmul(R,matmul(B,matmul(T,U)))
        stress = matmul(D,strain)
        
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1   - matmul(transpose(RBT),stress)*
     2                                wi*L*(h/2.d0)*w(kint)

        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1  + matmul(transpose(RBT),matmul(D,RBT))*
     2                                 wi*L*(h/2.d0)*w(kint)
        
        Tension = Tension + stress(1)*wi*(h/2.d0)*w(kint)
        
        Shear = Shear + stress(2)*wi*(h/2.d0)*w(kint)
        
        Moment = Moment + stress(1)*xi(2,kint)*wi*(h**2.d0/4.d0)*w(kint)
        
        if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(4*kint-3:4*kint) = stress(1:4)
        endif
      end do

      SVARS(1) = Tension
      SVARS(2) = Shear
      SVARS(3) = Moment

      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
    !
    !   Apply distributed loads
    !
    !   Distributed loads are specified in the input file using the Un option in the input file.
    !   n specifies the face number, following the ABAQUS convention
    !
      
      return

      END SUBROUTINE UEL_BEAM