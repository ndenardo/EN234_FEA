!
!    ABAQUS format user material subroutine for explicit dynamics
!
!

      subroutine vumat(nblock, ndir, nshr, nstatev, nfieldv, nprops,
     1  lanneal, stepTime, totalTime, dt, cmname, coordMp, charLength,
     2  props, density, strainInc, relSpinInc,
     3  tempOld, stretchOld, defgradOld, fieldOld,
     4  stressOld, stateOld, enerInternOld, enerInelasOld,
     5  tempNew, stretchNew, defgradNew, fieldNew,
     6  stressNew, stateNew, enerInternNew, enerInelasNew )
!
!      include 'vaba_param.inc'
!
      implicit double precision (a-h,o-z)

      dimension props(nprops), density(nblock), coordMp(nblock,*),
     1  charLength(nblock), strainInc(nblock,ndir+nshr),
     2  relSpinInc(nblock,nshr), tempOld(nblock),
     3  stretchOld(nblock,ndir+nshr),
     4  defgradOld(nblock,ndir+nshr+nshr),
     5  fieldOld(nblock,nfieldv), stressOld(nblock,ndir+nshr),
     6  stateOld(nblock,nstatev), enerInternOld(nblock),
     7  enerInelasOld(nblock), tempNew(nblock),
     8  stretchNew(nblock,ndir+nshr),
     9  defgradNew(nblock,ndir+nshr+nshr),
     1  fieldNew(nblock,nfieldv),
     2  stressNew(nblock,ndir+nshr), stateNew(nblock,nstatev),
     3  enerInternNew(nblock), enerInelasNew(nblock)

       character*80 cmname
!
!      Local variables
!
       integer k,ntens,nit,maxit

       double precision dedev(ndir+nshr),sdevstar(ndir+nshr)
       double precision sestar,skkstar
       double precision E,xnu,Y,e0,edot0,em
       double precision S,H,td,Omega,alpha
       double precision f,dfde
       double precision eplas,deplas,ta
       double precision dta, dtade, s0, ds0de
       double precision C, dCde
       double precision A, B, A2, B2, A3, A4, B4
       double precision dA2de, dB2de
!
!      Conventions for storing tensors:
!
!      Deformation gradients are provided as components in the global basis
!      (ABAQUS also allows the user to define a fixed local coord system defined with *ORIENTATION)
!
!      Stresses and stretches are stored as components in a basis that 'rotates with the material'
!      These are defined as follows:
!          Let e_i be the initial basis vectors (the global ijk directions, or user-specified basis vectors)
!          Let R be the rotation tensor, defined through the polar decomposition of deformation gradient F=RU=VR
!          Then define rotated basis vectors m_i = R e_i.
!          The co-rotational components of a tensor are defined as A_ij = m_i . A . m_j
!          The components A_ij can also be interpreted as the global components of the tensor  R^T A R
!
!
!      The components of symmetric tensors (stress,strainInc and stretch) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12]
!      The components of unsymmetric tensors (defgrad and relSpinInc) are stored as vectors in the following order
!                      For 3D problems (ndir=3,nshr=3) [11,22,33,12,23,31,21,32,13]
!                      For 2D problems (ndir=3,nshr=1) [11,22,33,12,21]
!
!      The stresses are Cauchy (true) stress
!
!      nblock                   No. integration points in data block (data are provided in 'blocks' of integration points)
!                               EN234FEA always has only 1 int pt in each block
!      ndir                     No. direct tensor components
!      nshr                     No. shear tensor components
!      nstatev                  No. user defined state variables (declared in input file)
!      nprops                   No. material properties
!      lanneal                  Annealing flag - if set to 1, then stresses are zeroed, and state vars should be reset to initial state
!      stepTime                 Elapsed time in this step at end of increment
!      totalTime                Total elapsed time at end of time increment
!      dt                       time step
!      cmname                   Material name
!      coordMp(n,i)             ith coord of nth integration point
!      charLength(i)            Characteristic length of element (in EN234FEA this is (element volume)^1/3)
!      props(k)                 kth material property
!      density                  density
!      strainInc(n,i)           ith strain increment component for nth int pt.  Strain components are in a
!                               basis that rotates with material, stored in order [de11, de22, de33, de12, de23, de13]
!                               NOTE THAT SHEAR STRAIN COMPONENTS ARE NOT DOUBLED IN THE VECTOR
!      relSpinInc(n,i)          ith component of relative spin increment tensor.   The global components of relSpinInc are
!                               dW_ij - dR_ij   where dW is the skew part of displacement gradient increment, and
!                               dR_ij is the rotation increment.
!      tempOld                  Temperature at start of increment
!      stretchOld(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradOld(n,i)          ith component of deformation gradient in fixed global basis at start of increment
!      fieldOld(n,i)            ith user-defined field variable at nth integration point at start of increment
!      stressOld(n,i)           ith component of Cauchy stress in co-rotational basis (equivalent to R^T sigma R in fixed basis)
!      stateOld(n,i)            ith user defined material state variable at start of increment
!      enerInternOld(n)         specific internal energy at nth integration point, start of increment
!      enerInelasOld(n)         specific irreversible dissipation at nth integration point, start of increment
!      tempNew                  Temperature at end of increment
!      stretchNew(n,i)          ith component of right stretch V in co-rotational basis (equivalent to U in global basis) at start of increment
!      defgradNew(n,i)          ith component of deformation gradient in fixed global basis at end of increment
!      fieldNew(n,i)            ith user-defined field variable at nth integration point at end of increment
!      stressNew(n,i)           ith component of Cauchy stress in co-rotational basis (equivale nt to R^T sigma R in fixed basis)
!      stateNew(n,i)            ith user defined material state variable at end of increment
!      enerInternNew(n)         specific internal energy at nth integration point, end of increment
!      enerInelasNew(n)         specific irreversible dissipation at nth integration point, end of increment
!
!       Coded for 3D problems; small stretch assumption (works approximately for finite rotations)
!
      E = props(1)
      xnu = props(2)
      Y = props(3)
      e0 = props(4)
      em  = props(5)
      edot0 = props(6)
      S = props(7)
      H = props(8)
      td = props(9)
      Omega = props(10)
      alpha = props(11)
      
      dta   = 0.d0
      dtade = 0.d0 
      C     = 0.d0 
      dCde  = 0.d0 
      s0    = 0.d0
      ds0de = 0.d0
      
      A  = 0.d0
      A2 = 0.d0
      A3 = 0.d0
      A4 = 0.d0
      B  = 0.d0
      B2 = 0.d0
      B4 = 0.d0
      
      dA2de = 0.d0
      dB2de = 0.d0
      
      ntens = ndir+nshr

      do k = 1,nblock

        eplas = stateOld(k,1)
        deplas = stateOld(k,3)
        ta = stateOld(k,2)

        !Deviatoric strain increment - used to get elastic predictor(s)
        dedev(1:ndir) = strainInc(k,1:ndir)
     1                      - sum(strainInc(k,1:ndir))/3.d0
        dedev(ndir+1:ntens) = strainInc(k,ndir+1:ntens)     !   No factor of 2 in dedev
        
        !Used in calculating (new) cauchy stress sigmaij^n+1 
        skkstar = sum(stressOld(k,1:ndir))                
     1          + E*sum(strainInc(k,1:ndir))/(1.d0-2.d0*xnu)
        
        !Deviatoric stress at the start of the step
        sdevstar(1:ndir) = stressOld(k,1:ndir)
     1                   - sum(stressOld(k,1:ndir))/3.d0
        sdevstar(ndir+1:ntens) = stressOld(k,ndir+1:ntens)
        sdevstar(1:ntens) = sdevstar(1:ntens)
     1                    + E*dedev(1:ntens)/(1.d0+xnu)
      
        !Predictor for effective stress
        sestar = dsqrt(1.5d0)*
     1    dsqrt(dot_product(sdevstar(1:ndir),sdevstar(1:ndir)) +
     2    2.d0*dot_product(sdevstar(ndir+1:ntens),
     3                     sdevstar(ndir+1:ntens)) )


        
      !Check if stress is in the elastic or plastic region.  If stress elastic, update by elastic, if plastic,plastic
      !Elastic increment (either no stress, or no plastic strain rate)
        if (sestar/Y<1.d-09.or.edot0==0.d0) then
           stressNew(k,1:ntens) = sdevstar(1:ntens)
           stressNew(k,1:ndir) = stressNew(k,1:ndir) + skkstar/3.d0
           stateNew(k,1) = eplas
           stateNew(k,2) = 0.d0
           stateNew(k,3) = ta
           cycle
        endif


      !Internal NR Loop - solve for plastic strain increment and increment in ta
       err = 1.d0
       maxit = 100
       nit = 1
       tol = 1.d-6*Y
       if (deplas==0.d0) deplas = 1.d-09/Y
       
       do while (err>tol)

          if (deplas<1.d-88) then
              
              dta = dt - ta*(deplas/Omega)**2.d0 
              
              dtade = -2.d0*ta*deplas/(Omega**2.d0)
              
          else
              
              A = Omega*dt/deplas
              B = exp(-deplas/Omega)
              
              A2 = A + ta*B - ta
              B2 = 1.d0 + Omega*B/deplas
              
              dta = A2/B2
              
              dA2de = -A/deplas - ta*B/Omega
              dB2de = -Omega*B/deplas**2.d0 - B/deplas
              
              dtade = (dA2de*B2 - A2*dB2de)/B2**2.d0 
              
          endif
      
         A3 = 1.d0 + (eplas + deplas)/e0
         s0 = Y*A3**em
         ds0de = (Y*em/e0)*A3**(em-1.d0)
         
         B4 = (ta + dta)/td
         A4 = exp(-B4**alpha)
         C = 1.d0 - A4
         dCde = A4*(alpha/td)*dtade*B4**(alpha-1.d0)
         
         f = sestar/S - 1.5d0*E*deplas/((1.d0+xnu)*S) - s0/S - H*C
     1        - log(deplas/(dt*edot0))
             
         dfde = -1.5d0*E/((1.d0+xnu)*S) - ds0de/S - H*dCde - 1.d0/deplas

         deplas_new = deplas - f/dfde
         
         if (deplas_new<0.d0) then
            deplas = deplas/10.d0
         else
            deplas = deplas_new
         endif

         nit = nit + 1
         err = dabs(f)
         if (nit>maxit) then
            write(6,*) ' Newton iterations in UMAT failed to converge '
            stop
         endif
       end do
    
      !Update stress and state variables
        stressNew(k,1:ntens) =
     1           ( 1.d0 - 1.5d0*deplas*E/(1.d0+xnu)/sestar )
     2                           *sdevstar(1:ntens)
        
        stressNew(k,1:ndir) =
     1              stressNew(k,1:ndir) + skkstar/3.d0

        stateNew(k,1) = eplas + deplas
        stateNew(k,3) = deplas
        stateNew(k,2) = ta + dta

      end do
!
      End Subroutine vumat