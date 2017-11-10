!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 3D linear elastic continuum element
!
!    The file also needs the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!=========================== ABAQUS format user element subroutine ===================

      
      
      SUBROUTINE UEL_C(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      integer      :: i,j,n_points,kint
    !
      double precision  ::  xi(2,9)                           ! Area integration points
      double precision  ::  w(9)                              ! Area integration weights
      double precision  ::  N(9)                              ! 2D shape functions
      double precision  ::  dNdxi(9,2)                        ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                         ! Spatial derivatives
      double precision  ::  Nbar(9)
      double precision  ::  dNbardxi(9,2)
      double precision  ::  dNbardx(9,2)
      double precision  ::  dxdxi(2,2)                        ! Derivative of spatial coords wrt normalized coords

      double precision  ::  sol(6), dsol(6)                   ! Sol vector contains [mu, c, dmudx1, dmudx2, dcdx1, dcdx2]
      double precision  ::  q(6)                              ! q vector defined in class
      double precision  ::  D(6,6)                            ! D matrix defined in class
      double precision  ::  B(6,18)                           ! p = B*U
      double precision  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
      double precision  ::  diffusion_coeft,kappa,theta       ! Material properties
      double precision  ::  c                                 ! concentration

    !
    !     Example ABAQUS UEL implementing 2D phase field model

      if (NNODE == 3) n_points = 4
      if (NNODE == 4) n_points = 4
      if (NNODE == 6) n_points = 4
      if (NNODE == 8) n_points = 9
      if (NNODE == 9) n_points = 9

      call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)

      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0

      diffusion_coeft = PROPS(1)
      kappa = PROPS(2)
      theta = PROPS(3)

      D = 0.d0
      D(1,1) = 1.d0
      D(2,2) = 1.d0/DTIME
      D(3,5) = -kappa
      D(4,6) = -kappa
      D(5,3) = theta*diffusion_coeft
      D(6,4) = theta*diffusion_coeft
  
    !     --  Loop over integration points
      do kint = 1, n_points
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
        dxdxi = matmul(coords(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
        determinant = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
        dxidx(1,1:2) =  [ dxdxi(2,2),-dxdxi(1,2)]/determinant
        dxidx(2,1:2) =  [-dxdxi(2,1),dxdxi(1,1) ]/determinant
        dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)
        B = 0.d0
        B(1,1:2*NNODE-1:2) = N(1:NNODE)
        B(2,2:2*NNODE:2) = N(1:NNODE)
        B(3,1:2*NNODE-1:2) = dNdx(1:NNODE,1)
        B(4,1:2*NNODE-1:2) = dNdx(1:NNODE,2)
        B(5,2:2*NNODE:2) = dNdx(1:NNODE,1)
        B(6,2:2*NNODE:2) = dNdx(1:NNODE,2)

        sol = matmul(B(1:6,1:2*NNODE),U(1:2*NNODE))      ! The p vector at the end of the step
        dsol = matmul(B(1:6,1:2*NNODE),DU(1:2*NNODE,1))    ! Increment in the p vector

        c = sol(2)
      
        q(1) = sol(1) - c*(c*c-1.d0)
        q(2) = dsol(2)/DTIME
        q(3) = -kappa*sol(5)
        q(4) = -kappa*sol(6)
        q(5) = diffusion_coeft*(sol(3)+(theta-1.d0)*dsol(3))
        q(6) = diffusion_coeft*(sol(4)+(theta-1.d0)*dsol(4))

        D(1,2) = -3.d0*c*c+1.d0


        RHS(1:2*NNODE,1) = RHS(1:2*NNODE,1) 
     1   - matmul(transpose(B(1:6,1:2*NNODE)),q)*w(kint)*determinant

        AMATRX(1:2*NNODE,1:2*NNODE) = AMATRX(1:2*NNODE,1:2*NNODE) 
     1   + matmul(transpose(B(1:6,1:2*NNODE)),
     2             matmul(D,B(1:6,1:2*NNODE)))*w(kint)*determinant

      end do
  

      return

      END SUBROUTINE UEL_C
