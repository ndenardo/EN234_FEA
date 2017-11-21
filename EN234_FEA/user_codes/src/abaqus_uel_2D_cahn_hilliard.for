!
!    ABAQUS format UEL subroutine
!
!    This file is compatible with both EN234_FEA and ABAQUS/Standard
!
!    The example implements a standard fully integrated 2D linear elastic continuum element
!
!    The file also needs the following subrouines:
!          abq_UEL_2D_integrationpoints           - defines integration points for 2D continuum elements
!          abq_UEL_2D_shapefunctions              - defines shape functions for 2D continuum elements
!          abq_UEL_1D_integrationpoints           - defines integration points for 1D line integral
!          abq_facenodes_2D                       - returns list of nodes on an element face
!=========================== ABAQUS format user element subroutine ===================

      SUBROUTINE UEL(RHS,AMATRX,SVARS,ENERGY,NDOFEL,NRHS,NSVARS,
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
      integer      :: i,j,n_points,kint, nfacenodes, ipoin, ksize
      integer      :: face_node_list(3)                      ! List of nodes on an element face

    ! Define 8-noded shape functions
      double precision  ::  xi(2,9)                           ! Area integration points
      double precision  ::  w(9)                              ! Area integration weights
      double precision  ::  N(9)                              ! Interpolation functions
      double precision  ::  dNdxi(9,2)                        ! 2D shape function derivatives
      double precision  ::  dNdx(9,2)                         ! Spatial derivatives
      double precision  ::  dxdxi(2,2)                        ! Derivative of spatial coords wrt normalized coords
      double precision  ::  dxidx(2,2), determinant           ! Jacobian inverse and determinant
      
    ! Define 4-noded shape functions corresponding to corner nodes
      double precision  ::  Nbar(9)                           ! Interpolation functions
      double precision  ::  dNbardxi(9,2)                     ! 2D shape function derivatives
      double precision  ::  dNbardx(9,2)                      ! Spatial derivatives
      
    ! Variables below are for computing integrals over element faces (not implemented)
      double precision  ::  face_coords(2,3)                  ! Coords of nodes on an element face
      double precision  ::  xi1(6)                            ! 1D integration points
      double precision  ::  w1(6)                             ! Integration weights
      double precision  ::  N1(3)                             ! 1D shape functions
      double precision  ::  dN1dxi(3)                         ! 1D shape function derivatives
      double precision  ::  norm(2)                           ! Normal to an element face
      double precision  ::  dxdxi1(2)                         ! Derivative of 1D spatial coord wrt normalized areal coord
    !
      double precision  ::  strain(4)                         ! Strain vector contains [e11, e22, e33, 2e12, 2e13, 2e23]
      double precision  ::  stress(4)                         ! Stress vector contains [s11, s22, s33, s12, s13, s23]
      double precision  ::  De(4,4)                           ! Elastic stress = De*(strain) (NOTE FACTOR OF 2 in shear strain)
      double precision  ::  B(9,24)                           ! strain = B*(dof_total)
      double precision  ::  E, xnu, De11, De12, De44          ! Elastic moduli and contitutive properties 
      double precision  ::  Omega, We, Kappa, Diff, Theta     ! Diffusion properties
                                          
      double precision  ::  c                                 ! Concentration states at integration points
      double precision  ::  f, dfdc                           ! function of concentration and its derivative
      double precision  ::  P(9), dP(9)                       ! Vectors used in computing Q
      double precision  ::  Q(9)                              ! Vector of elastic, concentration stresses
      double precision  ::  D(9,9)                            ! Constitutive matrix
      double precision  ::  skk                               ! Summation of normal stress values
      
    !     ABAQUS UEL implementing 2D elements for solution of cahn-hilliard equation, elasticity
      
      
      n_points = 4
     
      call abq_UEL_2D_integrationpoints(n_points, NNODE, xi, w)
 
      if (MLVARX<2*NNODE) then
        write(6,*) ' Error in abaqus UEL '
        write(6,*) ' Variable MLVARX must exceed 2*NNODE'
        write(6,*) ' MLVARX = ',MLVARX,' NNODE = ',NNODE
        stop
      endif
      
      RHS(1:MLVARX,1) = 0.d0
      AMATRX(1:NDOFEL,1:NDOFEL) = 0.d0
      stress = 0.d0
      strain = 0.d0
      skk    = 0.d0
      
      E     = PROPS(1)
      xnu   = PROPS(2)
      Omega = PROPS(3)
      We    = PROPS(4)
      Kappa = PROPS(5)
      Diff  = PROPS(6)
      Theta = PROPS(7)
      
      De11 = (1.d0-xnu)             * (E/((1.d0+xnu)*(1.d0-2.d0*xnu)))
      De12 = xnu                    * (E/((1.d0+xnu)*(1.d0-2.d0*xnu)))
      De44 = ((1.d0-2.d0*xnu)/2.d0) * (E/((1.d0+xnu)*(1.d0-2.d0*xnu)))
      
      De = 0.d0
      De(1:3,1:3) = De12
      De(1,1) = De11
      De(2,2) = De11
      De(3,3) = De11
      De(4,4) = De44

      ENERGY(1:8) = 0.d0

    ! -- Loop over integration points
      do kint = 1, n_points
        
      ! define shape functions for different nodal schemes
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),NNODE,N,dNdxi)
        call abq_UEL_2D_shapefunctions(xi(1:2,kint),4,Nbar,dNbardxi)
        
        dxdxi = matmul(COORDS(1:2,1:NNODE),dNdxi(1:NNODE,1:2))
        
      ! invert dxdxi
        determinant = dxdxi(1,1)*dxdxi(2,2)-dxdxi(1,2)*dxdxi(2,1)
        dxidx(1,1:2) =  [ dxdxi(2,2),-dxdxi(1,2)]/determinant
        dxidx(2,1:2) =  [-dxdxi(2,1),dxdxi(1,1) ]/determinant

        dNdx(1:NNODE,1:2) = matmul(dNdxi(1:NNODE,1:2),dxidx)
        dNbardx(1:4,1:2) = matmul(dNbardxi(1:4,1:2),dxidx)
        
      ! Assemble the B matrix
        B = 0.d0
        B(1,1:2*NNODE-3:4) = dNdx(1:NNODE/2,1)
        B(2,2:2*NNODE-2:4) = dNdx(1:NNODE/2,2)
        B(3,1:2*NNODE-3:4) = dNdx(1:NNODE/2,2)
        B(3,2:2*NNODE-2:4) = dNdx(1:NNODE/2,1)
        B(4,3:2*NNODE-1:4) = Nbar(1:4)
        B(5,4:2*NNODE:4)   = Nbar(1:4)
        B(6,3:2*NNODE-1:4) = dNbardx(1:4,1)
        B(7,3:2*NNODE-1:4) = dNbardx(1:4,2)
        B(8,4:2*NNODE:4)   = dNbardx(1:4,1)
        B(9,4:2*NNODE:4)   = dNbardx(1:4,2)
        
        B(1,2*NNODE+1:3*NNODE-1:2) = dNdx(NNODE/2+1:NNODE,1)
        B(2,2*NNODE+2:3*NNODE:2)   = dNdx(NNODE/2+1:NNODE,2)
        B(3,2*NNODE+1:3*NNODE-1:2) = dNdx(NNODE/2+1:NNODE,2)
        B(3,2*NNODE+2:3*NNODE:2)   = dNdx(NNODE/2+1:NNODE,1)
        
      ! Compute the P, dP vectors to get mu, c, dc
        P  = 0.d0
        dP = 0.d0
        P  = matmul(B(1:9,1:3*NNODE),U(1:3*NNODE))
        dP = matmul(B(1:9,1:3*NNODE),DU(1:3*NNODE,1))
        
        c  = 0.d0
        c  = P(5)  !c = c  + dc

      ! Assemble D Matrix
        dfdc = 0.d0
        dfdc = We*(12.d0*(c**2.d0) - 12.d0*c + 2.d0)
        
        D = 0.d0        
        D(1:2,1:2) = De(1:2,1:2)
        D(3,3) = De(4,4)
        D(4,1) = -Omega*(De(1,1)+De(1,2)+De(1,3))/3.d0
        D(1,5) = -Omega*(De(1,1)+De(1,2)+De(1,3))/3.d0
        D(4,2) = -Omega*(De(2,1)+De(2,2)+De(2,3))/3.d0
        D(2,5) = -Omega*(De(2,1)+De(2,2)+De(2,3))/3.d0
        D(4,4) = 1.d0
        D(4,5) = -dfdc - ((Omega**2.d0)/9.d0) * 
     1   (De(1,1)+De(1,2)+De(1,3)+De(2,1)+De(2,2)+De(2,3)+
     2    De(3,1)+De(3,2)+De(3,3))
        D(5,5) = 1.d0/DTIME
        D(8,6) = Theta*Diff
        D(9,7) = Theta*Diff
        D(6,8) = -Kappa
        D(7,9) = -Kappa
        
      ! Compute stress values, skk
        strain(1) = P(1) - Omega*(c)/3.d0
        strain(2) = P(2) - Omega*(c)/3.d0
        strain(3) = 0    - Omega*(c)/3.d0
        strain(4) = P(3) 
        stress(1:4) = matmul(De,strain(1:4))
        skk = stress(1) + stress(2) + stress(3)
        
      ! Assemble Q vector
        f = 0.d0
        f = 2.d0*We*(c)*(c-1.d0)*((2.d0*c)-1.d0) ! = f(c + dc)
        
        Q    = 0.d0
        Q(1) = stress(1)
        Q(2) = stress(2)
        Q(3) = stress(4)
        Q(4) = P(4) - f - Omega*skk/3.d0
        Q(5) = dP(5)/DTIME
        Q(6) = -Kappa*P(8)
        Q(7) = -Kappa*P(9)
        Q(8) = Diff*(P(6)+(Theta-1.d0)*dP(6))
        Q(9) = Diff*(P(7)+(Theta-1.d0)*dP(7))
        
      ! Assemble element stiffness and residual 
        RHS(1:3*NNODE,1) = RHS(1:3*NNODE,1)
     1  - matmul(transpose(B(1:9,1:3*NNODE)),Q)*w(kint)*determinant 

        AMATRX(1:3*NNODE,1:3*NNODE) = AMATRX(1:3*NNODE,1:3*NNODE)
     1  + matmul(transpose(B(1:9,1:3*NNODE)),matmul(D,B(1:9,1:3*NNODE)))
     2                                             *w(kint)*determinant
                
        ENERGY(2) = ENERGY(2)
     1   + 0.5D0*dot_product(Q,P)*w(kint)*determinant           ! Store the elastic strain energy

        if (NSVARS>=n_points*4) then   ! Store stress at each integration point (if space was allocated to do so)
            SVARS(4*kint-3:4*kint) = stress(1:4)
        endif
        
      end do

      PNEWDT = 1.d0          ! This leaves the timestep unchanged (ABAQUS will use its own algorithm to determine DTIME)
      
      return

      END SUBROUTINE UEL