!
!    ABAQUS format user material subroutine for small strain hypoelastic material
!
!

      SUBROUTINE UMAT(STRESS,STATEV,DDSDDE,SSE,SPD,SCD,
     1 RPL,DDSDDT,DRPLDE,DRPLDT,
     2 STRAN,DSTRAN,TIME,DTIME,TEMP,DTEMP,PREDEF,DPRED,CMNAME,
     3 NDI,NSHR,NTENS,NSTATV,PROPS,NPROPS,COORDS,DROT,PNEWDT,
     4 CELENT,DFGRD0,DFGRD1,NOEL,NPT,LAYER,KSPT,KSTEP,KINC)
!
      INCLUDE 'ABA_PARAM.INC'
!     WARNING - the aba_param.inc file declares
!        Implicit double precision (a-h,o-z)
!     This means that, by default, any variables with
!     first letter between a-h or o-z are double precision.
!     The rest are integers.
!     Note that this also means that if you type a variable
!     name incorrectly, the compiler won't catch your typo.
!
      CHARACTER*80 CMNAME
      DIMENSION STRESS(NTENS),STATEV(NSTATV),
     1 DDSDDE(NTENS,NTENS),DDSDDT(NTENS),DRPLDE(NTENS),
     2 STRAN(NTENS),DSTRAN(NTENS),TIME(2),PREDEF(1),DPRED(1),
     3 PROPS(NPROPS),COORDS(3),DROT(3,3),DFGRD0(3,3),DFGRD1(3,3)

!
!      Relevant manual sections:
!
!      Theory guide 1.4.3 discusses ABAQUS strain rates (but doesnt mention
!      User subroutine reference manual 1.1.44 - describes UMAT
!      Theory manual 3.2.2 Solid element library has discussion of Hughes-Winget algorithm
!      and approximation to consistent tangent used by ABAQUS
!
!      DDSDDE(NTENS,NTENS)
!         Jacobian matrix of the constitutive model.
!         DDSDDE(I,J) defines the change in the Ith stress component
!         at the end of the time increment caused by an infinitesimal
!         perturbation of the Jth component of the strain increment array.
!         Unless you invoke the unsymmetric equation solution capability
!         for the user-defined material, ABAQUS/Standard will use only
!         the symmetric part of DDSDDE. The symmetric part of the matrix
!         is calculated by taking one half the sum of the matrix and its transpose.
!         EN234FEA always uses the full matrix (no symmetry is assumed)

!      STRESS(NTENS)
!         This array is passed in as the stress tensor at the beginning
!         of the increment and must be updated in this routine to be the
!         stress tensor at the end of the increment. If you specified
!         initial stresses (“Initial conditions,” Section 19.2.1), this
!         array will contain the initial stresses at the start of the
!         analysis. The size of this array depends on the value of NTENS
!         as defined below. In finite-strain problems the stress tensor
!         has already been rotated to account for rigid body motion in
!         the increment before UMAT is called, so that only the corotational
!         part of the stress integration should be done in UMAT. The
!         measure of stress used is “true” (Cauchy) stress.
!
!   NB      When used in ABAQUS with hybrid elements the stress array has different dimensions
!           and additional variables must be defined.   See the ABAQUS manual for details.
!           This behavior is not implemented in EN234FEA
!
!      STATEV(NSTATV)
!         An array containing the solution-dependent state variables.
!         These are passed in as the values at the beginning of the
!         increment unless they are updated in user subroutines USDFLD
!        (“USDFLD,” Section 25.2.39) or UEXPAN (“UEXPAN,” Section 25.2.20),
!        in which case the updated values are passed in. In all cases
!         STATEV must be returned as the values at the end of the increment.
!         The size of the array is defined as described in
!        “Allocating space” in “User subroutines: overview,” Section 25.1.1.
!
!         In finite-strain problems any vector-valued or tensor-valued
!         state variables must be rotated to account for rigid body
!         motion of the material, in addition to any update in the
!         values associated with constitutive behavior. The rotation
!         increment matrix, DROT, is provided for this purpose.
!
!      SSE, SPD, SCD
!         Specific elastic strain energy, plastic dissipation, and
!         “creep” dissipation, respectively. These are passed in as
!         the values at the start of the increment and should be
!         updated to the corresponding specific energy values at
!         the end of the increment. They have no effect on the solution,
!         except that they are used for energy output.
!
!     Only in a fully coupled thermal-stress analysis
!      RPL
!         Volumetric heat generation per unit time at the end of the increment
!         caused by mechanical working of the material.
!
!     DDSDDT(NTENS)
!          Variation of the stress increments with respect to the temperature.
!
!     DRPLDE(NTENS)
!           Variation of RPL with respect to the strain increments.
!
!     DRPLDT
!           Variation of RPL with respect to the temperature.
!
!     Variables that can be updated
!
!     PNEWDT
!        Ratio of suggested new time increment to the time increment being
!        used (DTIME, see discussion later in this section). This variable
!        allows you to provide input to the automatic time incrementation
!        algorithms in ABAQUS/Standard (if automatic time incrementation is chosen).
!        For a quasi-static procedure the automatic time stepping that ABAQUS/Standard
!        uses, which is based on techniques for integrating standard creep laws
!        (see “Quasi-static analysis,” Section 6.2.5), cannot be controlled from within
!        the UMAT subroutine.
!        PNEWDT is set to a large value before each call to UMAT.
!        If PNEWDT is redefined to be less than 1.0, ABAQUS/Standard must abandon the
!        time increment and attempt it again with a smaller time increment. The
!        suggested new time increment provided to the automatic time integration
!        algorithms is PNEWDT × DTIME, where the PNEWDT used is the minimum value
!        for all calls to user subroutines that allow redefinition of PNEWDT for this
!        iteration.
!        If PNEWDT is given a value that is greater than 1.0 for all calls to user
!        subroutines for this iteration and the increment converges in this iteration,
!        ABAQUS/Standard may increase the time increment. The suggested new time increment
!        provided to the automatic time integration algorithms is PNEWDT × DTIME, where
!        the PNEWDT used is the minimum value for all calls to user subroutines for
!        this iteration.
!        If automatic time incrementation is not selected in the analysis procedure,
!        values of PNEWDT that are greater than 1.0 will be ignored and values of
!        PNEWDT that are less than 1.0 will cause the job to terminate.
!
!    Variables passed in for information
!
!     STRAN(NTENS)
!         An array containing the total strains at the beginning of the increment.
!         If thermal expansion is included in the same material definition, the
!         strains passed into UMAT are the mechanical strains only (that is, the
!         thermal strains computed based upon the thermal expansion coefficient have
!         been subtracted from the total strains). These strains are available for output
!         as the “elastic” strains.
!
!         In finite-strain problems the strain components have been rotated to account for
!         rigid body motion in the increment before UMAT is called and are approximations
!         to logarithmic strain.

!     DSTRAN(NTENS)
!         Array of strain increments. If thermal expansion is included in the same
!         material definition, these are the mechanical strain increments (the total
!         strain increments minus the thermal strain increments).
!
!     TIME(1)
!         Value of step time at the beginning of the current increment.
!
!     TIME(2)
!          Value of total time at the beginning of the current increment.
!
!     DTIME
!        Time increment.
!
!     TEMP
!         Temperature at the start of the increment.
!
!     DTEMP
!         Increment of temperature.
!
!     PREDEF
!        Array of interpolated values of predefined field variables at this point
!        at the start of the increment, based on the values read in at the nodes.
!
!      DPRED
!        Array of increments of predefined field variables.
!
!      CMNAME
!        User-defined material name, left justified. Some internal material models are given names starting with the “ABQ_” character string. To avoid conflict, you should not use “ABQ_” as the leading string for CMNAME.
!
!      NDI
!        Number of direct stress components at this point.
!
!      NSHR
!        Number of engineering shear stress components at this point.
!
!      NTENS
!        Size of the stress or strain component array (NDI + NSHR).
!
!      NSTATV
!         Number of solution-dependent state variables that are associated with
!         this material type (defined as described in “Allocating space” in “User
!         subroutines: overview,” Section 25.1.1).
!
!      PROPS(NPROPS)
!         User-specified array of material constants associated with this user material.
!
!      NPROPS
!         User-defined number of material constants associated with this user material.
!
!      COORDS
!         An array containing the coordinates of this point. These are the current
!         coordinates if geometric nonlinearity is accounted for during the step
!         (see “Procedures: overview,” Section 6.1.1); otherwise, the array contains
!         the original coordinates of the point.
!
!     DROT(3,3)
!          Rotation increment matrix. This matrix represents the increment of rigid
!          body rotation of the basis system in which the components of stress
!          (STRESS) and strain (STRAN) are stored. It is provided so that vector- or
!          tensor-valued state variables can be rotated appropriately in this subroutine:
!          stress and strain components are already rotated by this amount before UMAT
!          is called. This matrix is passed in as a unit matrix for small-displacement
!          analysis and for large-displacement analysis if the basis system for the
!          material point rotates with the material (as in a shell element or when a
!          local orientation is used).
!
!      CELENT
!          Characteristic element length, which is a typical length of a line across
!          an element for a first-order element; it is half of the same typical length
!          for a second-order element. For beams and trusses it is a characteristic length
!          along the element axis. For membranes and shells it is a characteristic length
!          in the reference surface. For axisymmetric elements it is a characteristic length
!          in the  plane only. For cohesive elements it is equal to the constitutive
!          thickness.
!
!      DFGRD0(3,3)
!          Array containing the deformation gradient at the beginning of the increment.
!          See the discussion regarding the availability of the deformation gradient for
!          various element types.
!
!     DFGRD1(3,3)
!            Array containing the deformation gradient at the end of the increment.
!           The components of this array are set to zero if nonlinear geometric effects
!           are not included in the step definition associated with this increment. See
!           the discussion regarding the availability of the deformation gradient for
!           various element types.
!
!      NOEL
!           Element number.
!
!      NPT
!           Integration point number.
!
!      LAYER
!          Layer number (for composite shells and layered solids).
!
!      KSPT
!          Section point number within the current layer.
!
!      KSTEP
!         Step number.
!
!     KINC
!         Increment number.

!      user coding to define DDSDDE, STRESS, STATEV, SSE, SPD, SCD
!      and, if necessary, RPL, DDSDDT, DRPLDE, DRPLDT, PNEWDT
!
!     Local variables

      double precision :: edev(6)
      double precision :: evol
      double precision :: se,ee
      double precision :: s0,e0,n,K,dsedee,Et,Es

      integer :: j

       s0 = PROPS(1)
       e0 = PROPS(2)
       n = PROPS(3)
       K = PROPS(4)

       evol = sum(STRAN(1:3)+DSTRAN(1:3))
       edev(1:3) = STRAN(1:3)+DSTRAN(1:3) - evol/3.d0
       edev(4:6) = 0.5d0*(STRAN(4:6)+DSTRAN(4:6))

       ee = dsqrt(dot_product(edev(1:3),edev(1:3)) +
     1               2.d0*dot_product(edev(4:6),edev(4:6)))/dsqrt(1.5d0)

       DDSDDE(1:6,1:6) = 0.d0
       if (ee<e0) then
          se = s0*( dsqrt( (1.d0+n*n)/((n-1.d0)*(n-1.d0)) 
     1        - (n/(n-1.d0) - ee/e0)**2.d0 ) - 1.d0/(n-1.d0) )
          dsedee = (s0/e0)*(n/(n-1.d0) - ee/e0)/
     1 dsqrt((1.d0+n*n)/((n-1.d0)*(n-1.d0)) - (n/(n-1.d0)-ee/e0)**2.d0)
          if (ee==0.d0) then
             stress = 0.d0
             Es = n*s0/e0
          else
             Es = se/ee
             Et = dsedee
             stress = 2.d0*se*edev/(3.d0*ee)
             stress(1:3) = stress(1:3) + K*evol
             DDSDDE(1:6,1:6) = 
     1          4.d0*(Et-Es)*spread(edev,dim=2,ncopies=6)*
     1                       spread(edev,dim=1,ncopies=6)/(9.d0*ee*ee)
          endif
       else
          se = s0*(ee/e0)**(1.d0/n)
          stress = 2.d0*se*edev/(3.d0*ee)
          stress(1:3) = stress(1:3) + K*evol
          Et = se/(n*ee)
          Es = se/ee
          DDSDDE(1:6,1:6) = 4.d0*(Et-Es)*spread(edev,dim=2,ncopies=6)*
     1                      spread(edev,dim=1,ncopies=6)/(9.d0*ee*ee)
       endif
       forall(j=1:3) DDSDDE(j,j) = DDSDDE(j,j) + 2.d0*Es/3.d0
       forall(j=4:6) DDSDDE(j,j) = DDSDDE(j,j) + Es/3.d0
       DDSDDE(1:3,1:3) =  DDSDDE(1:3,1:3) + (K-2.d0*Es/9.d0)

       return



      RETURN
      END SUBROUTINE UMAT