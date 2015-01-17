      module integrals_amplitudes
      integer, save :: SpOrbNumb                                   ! Total number of spin-orbitals
      integer, save :: OnePartIntNumb, TwoPartIntNumb              ! Total number of one- and two- particle integrals, respectively
      real*8, allocatable, save :: amplitudes(:)                   ! 1-D array with zero-order Hamiltonian parametrization
      real*8, allocatable, save :: OnePartInt(:), TwoPartInt(:)    ! 1-D arrays with one- and two- particle integrals, respectively
      real*8, save :: PhysVacEner, MeanEner, NuclRepEner           ! Physical Vacuum energy, Mean energy of the lowest energy cluster, nuclear repulsion energy
      real*8, save :: ZeroOrderVacEner
!      real*8, allocatable, save :: xRes(:), yRes(:), zRes(:)

      contains

      real*8 function twoparticle1 (SpOrb_I, SpOrb_J, SpOrb_K, SpOrb_L)
      implicit none
      integer SpOrb_I, SpOrb_J, SpOrb_K, SpOrb_L
      integer OrbCopy_I, OrbCopy_J, OrbCopy_K, OrbCopy_L
      integer SerNumb_IK, SerNumb_JL, SerNumb_IJKL

      if ((SpOrb_I + SpOrb_K)/2*2 .ne. SpOrb_I + SpOrb_K .or. & 
          (SpOrb_J + SpOrb_L)/2*2 .ne. SpOrb_J + SpOrb_L ) then

            twoparticle1 = 0.00D0

       else

          if (SpOrb_I .gt. SpOrb_K) then
             if (SpOrb_I/2*2 .eq. SpOrb_I) then
                OrbCopy_I = (SpOrb_K)/2
                OrbCopy_K = (SpOrb_I)/2
              else
                OrbCopy_I = (SpOrb_K + 1)/2
                OrbCopy_K = (SpOrb_I + 1)/2
             endif   
           else
             if (SpOrb_I/2*2 .eq. SpOrb_I) then
                OrbCopy_I = (SpOrb_I)/2
                OrbCopy_K = (SpOrb_K)/2
              else
                OrbCopy_I = (SpOrb_I + 1)/2
                OrbCopy_K = (SpOrb_K + 1)/2
             endif
          endif

          if (SpOrb_J .gt. SpOrb_L) then
             if (SpOrb_J/2*2 .eq. SpOrb_J) then
                OrbCopy_J = (SpOrb_L)/2
                OrbCopy_L = (SpOrb_J)/2
              else
                OrbCopy_J = (SpOrb_L + 1)/2
                OrbCopy_L = (SpOrb_J + 1)/2
             endif
           else
             if (SpOrb_J/2*2 .eq. SpOrb_J) then
                OrbCopy_J = (SpOrb_J)/2
                OrbCopy_L = (SpOrb_L)/2
              else
                OrbCopy_J = (SpOrb_J + 1)/2
                OrbCopy_L = (SpOrb_L + 1)/2
             endif
          endif

          SerNumb_IK = ((OrbCopy_I - 1) * (SpOrbNumb - OrbCopy_I))/2 + OrbCopy_K
          SerNumb_JL = ((OrbCopy_J - 1) * (SpOrbNumb - OrbCopy_J))/2 + OrbCopy_L

          SerNumb_IJKL = (min(SerNumb_IK, SerNumb_JL) - 1) * (SpOrbNumb/2*(SpOrbNumb/2 + 1) - & 
                          min(SerNumb_IK, SerNumb_JL))/2 + max(SerNumb_IK, SerNumb_JL)
          
          twoparticle1 = TwoPartInt(SerNumb_IJKL)

       endif
       end function

      subroutine ZeroOrderVacuumEnergy(OrbNumb, PhysVac)
      implicit none
      integer i, j
      integer OrbNumb
      integer PhysVac(OrbNumb)

      ZeroOrderVacEner = 0.00D0

        do i = 1, OrbNumb
          if (PhysVac(i) .ne. 0) then
            ZeroOrderVacEner = ZeroOrderVacEner + float(PhysVac(i)) * amplitudes(i)
          endif
        enddo

      end subroutine

      subroutine integral_transformation(OrbNumb, CorePlusValNumb, PhysVac)
      implicit none
      integer i, j
      integer res1, res2
      integer CorePlusValNumb
      integer SpOrb_Left, SpOrb_Right
      integer OrbNumb
      integer PhysVac(OrbNumb)
      integer SpOrbPhysVac(2 * OrbNumb)

      open (10, file = 'Fock Matrix Diagonal Elements')

      res1 = 1
      res2 = 2
      SpOrbPhysVac = 0

        do i = 1, CorePlusValNumb
          if (PhysVac(i) .eq. 1) then
             SpOrbPhysVac(2 * i - 1) = 1
           elseif (PhysVac(i) .eq. -1) then
             SpOrbPhysVac(2 * i) = 1
           elseif (PhysVac(i) .eq. 2) then
             SpOrbPhysVac(2 * i - 1) = 1
             SpOrbPhysVac(2 * i) = 1
          endif
        enddo

        do SpOrb_Left = 1, SpOrbNumb
          do SpOrb_Right = SpOrb_Left, SpOrbNumb

             if (SpOrb_Left/2*2 .ne. SpOrb_Left .and. SpOrb_Right/2*2 .ne. SpOrb_Right) then
                do i = 1, 2 * OrbNumb
                  if (SpOrbPhysVac(i) .eq. 1) then
                    OnePartInt(res1) = OnePartInt(res1) + twoparticle1(SpOrb_Left, i, SpOrb_Right, i) - & 
                                      twoparticle1(SpOrb_Left, i, i, SpOrb_Right)
                  endif
                enddo
                  if (SpOrb_Left .eq. SpOrb_Right) then
                    write(10, 30) SpOrb_Left, OnePartInt(res1)
                  endif
               res1 = res1 + 2

             elseif (SpOrb_Left/2*2 .eq. SpOrb_Left .and. SpOrb_Right/2*2 .eq. SpOrb_Right) then
                do i = 1, 2 * OrbNumb
                  if (SpOrbPhysVac(i) .eq. 1) then
                    OnePartInt(res2) = OnePartInt(res2) + twoparticle1(SpOrb_Left, i, SpOrb_Right, i) - & 
                                      twoparticle1(SpOrb_Left, i, i, SpOrb_Right)
                  endif
                enddo
                  if (SpOrb_Left .eq. SpOrb_Right) then
                    write(10, 30) SpOrb_Left, OnePartInt(res2)
                  endif
               res2 = res2 + 2

             endif
          enddo
        enddo
      close (10)
!
!    TO BE DELETED
!
!
!      open (220, file = 'OneElectron')
!      write (220, *) OnePartIntNumb
!      write (220, *) ' '
!       do i = 1, OnePartIntNumb
!          write (220, 240) OnePartInt(2*i - 1), OnePartInt(2*i)
!       enddo
!      close (220)
!  240 format(2F12.6)
!
!
!
   30 format (I4, F12.6)
      end subroutine

      real*8 function oneparticle1 (SpOrbLeft, SpOrbRight)
      implicit none
      integer i, j, temp
      integer SpOrbLeft, SpOrbRight
      integer SpOrbLeftCopy, SpOrbRightCopy
      integer OrbLeftCopy, OrbRightCopy
      integer SerNumb
      real*8, external :: twoparticle1

      SpOrbLeftCopy = SpOrbLeft
      SpOrbRightCopy = SpOrbRight

        if (SpOrbLeftCopy .gt. SpOrbRightCopy) then 
            temp = SpOrbLeftCopy
            SpOrbLeftCopy = SpOrbRightCopy
            SpOrbRightCopy = temp
        endif

        if ((SpOrbLeftCopy + SpOrbRightCopy)/2*2 .ne. SpOrbLeftCopy + SpOrbRightCopy) then
            oneparticle1 = 0.0D0

        else

            if (SpOrbLeftCopy/2*2 .eq. SpOrbLeftCopy) then 
                 OrbLeftCopy = SpOrbLeftCopy/2
                 OrbRightCopy = SpOrbRightCopy/2
                 SerNumb = ((OrbLeftCopy - 1) * (SpOrbNumb - OrbLeftCopy))/2 + OrbRightCopy
                 oneparticle1 = OnePartInt(2 * SerNumb)

            else
                 OrbLeftCopy = (SpOrbLeftCopy + 1)/2
                 OrbRightCopy = (SpOrbRightCopy + 1)/2
                 SerNumb = ((OrbLeftCopy - 1) * (SpOrbNumb - OrbLeftCopy))/2 + OrbRightCopy
                 oneparticle1 = OnePartInt(2 * SerNumb - 1)
            
            endif
        endif

       end function 

       end module
