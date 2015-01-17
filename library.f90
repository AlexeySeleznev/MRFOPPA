      subroutine integrals_amplitudes_read (OrbNumb)
      use integrals_amplitudes
      implicit none
      integer i
      integer OrbNumb

      allocate (amplitudes(OrbNumb))
      allocate (OnePartInt(2 * OnePartIntNumb), TwoPartInt(TwoPartIntNumb))

      open (10, file = 'One.txt')
      open (20, file = 'Two.txt')
      open (30, file = 'Energy1.txt')
      
        do i = 1, OrbNumb
          read (30, *) amplitudes(i)
        enddo
        do i = 1, OnePartIntNumb
          read (10, *) OnePartInt(2 * i)
          OnePartInt(2 * i - 1) = OnePartInt(2 * i)
        enddo
        do i = 1, TwoPartIntNumb
          read (20, *) TwoPartInt(i)
        enddo

      close (10)
      close (20)
      close (30)
      end subroutine

      subroutine PhysVacEnerCalc (OrbNumb, ValPlusCoreNumb, PhysVac)
      use integrals_amplitudes
      implicit none
      integer i, j
      integer OrbNumb, ValPlusCoreNumb
      integer PhysVac(OrbNumb)
      integer SpOrbPhysVac(2 * OrbNumb)

      PhysVacEner = 0.00D0
      SpOrbPhysVac = 0

        do i = 1, OrbNumb
          if (PhysVac(i) .eq. 1) then
             SpOrbPhysVac(2 * i - 1) = 1
           elseif (PhysVac(i) .eq. -1) then
             SpOrbPhysVac(2 * i) = 1
           elseif (PhysVac(i) .eq. 2) then
             SpOrbPhysVac(2 * i - 1) = 1
             SpOrbPhysVac(2 * i) = 1
          endif
        enddo

        do i = 1, 2 * ValPlusCoreNumb
          if (SpOrbPhysVac(i) .ne. 0) then
            PhysVacEner = PhysVacEner + oneparticle1(i, i)
          endif
        enddo
        do i = 1, 2 * ValPlusCoreNumb - 1
          do j = i + 1, 2 * ValPlusCoreNumb
            if (SpOrbPhysVac(i) .ne. 0 .and. SpOrbPhysVac(j) .ne. 0) then
               PhysVacEner = PhysVacEner + twoparticle1(i, j, i, j) - & 
                             twoparticle1(i, j, j, i)
            endif
          enddo
        enddo

      PhysVacEner = PhysVacEner + NuclRepEner
      end subroutine

      subroutine order (n, L, phase)
      implicit none
      integer i, j, k
      integer n
      integer tem
      integer L(n)
      real*8 phase
      k = 0
      phase = 1.00D0
         do i = 1, n - 1
             do j = i + 1, n
                if (L(i) .eq. L(j) .and. L(i) .ne. 0 .and. L(j) .ne. 0) then
                    k = k + 1
                 elseif (L(i) .gt. L(j) .and. L(i) .ne. 0 .and. L(j) .ne. 0) then
                   phase = phase * (-1.00D0)
                   tem = L(i)
                   L(i) = L(j)
                   L(j) = tem
                 elseif (L(i) .eq. 0 .and. L(i) .lt. L(j)) then
                   tem = L(i)
                   L(i) = L(j)
                   L(j) = tem
                endif
             enddo
         enddo
             if (k .ne. 0) then
                phase = 0.00D0
             endif
      end subroutine

      subroutine interphase_calculator (CreatDim, AnnihDim, OrbNumb, Creation, Annihilation, SpinOrbPhysVacMod, interphase)
      implicit none
      integer i, j
      integer OrbNumb
      integer CreatDim, AnnihDim
      integer Creation(CreatDim), Annihilation(AnnihDim)
      integer SpinOrbPhysVacMod(2 * OrbNumb), Dublicate(2 * OrbNumb)
      integer sum
      integer multiplier/1/
      real*8 interphase

      Dublicate = SpinOrbPhysVacMod
!      print *, '/////////////////\\\\\\\\\\\\\\\\'
!      print *, 'SpinOrbPhysVacMod = ', Dublicate(:)
!      print *, 'Annihilation = ', Annihilation(:)

      if (Annihdim .gt. 0) then
         do i = AnnihDim, 1, -1
            Dublicate(Annihilation(i)) = Dublicate(Annihilation(i)) - 1
              if (Annihilation(i) .gt. 1) then
                sum = 0
                  do j = 1, Annihilation(i) - 1
                    sum = sum + Dublicate(j)
                  enddo
!                print *, '   sum = ', sum
                multiplier = multiplier * ((-1)**(sum))
              endif
         enddo
      endif
!      print *, 'SpinOrbPhysVacMod = ', Dublicate(:)
!      print *, 'multiplier = ', multiplier
!      print *, 'Creation = ', Creation(:)

      if (Creatdim .gt. 0) then
         do i = CreatDim, 1, -1
            Dublicate(Creation(i)) = Dublicate(Creation(i)) + 1
              if (Creation(i) .gt. 1) then
                sum = 0
                  do j = 1, Creation(i) - 1
                    sum = sum + Dublicate(j)
                  enddo
!                print *, '   sum = ', sum
                multiplier = multiplier * ((-1) ** (sum))
              endif
         enddo
      endif
!      print *, 'SpinOrbPhysVacMod = ', Dublicate(:)
!      print *, 'multiplier = ', multiplier
!      print *, '////////////////\\\\\\\\\\\\\\\\\\\\\\\'
      interphase = real(multiplier)
      end subroutine

      subroutine merge_spin_tensor (creation_dim, annihilation_dim, NonZeroAnnihNumb, stack_dim, creation, & 
                                    annihilation, stack, phase)
      implicit none
      integer i, j
      integer creation_dim, annihilation_dim, stack_dim
      integer creation(creation_dim), annihilation(annihilation_dim), stack(stack_dim)
      real*8 phase
      integer counter
      integer multiplier
      integer counter_position
      integer NonZeroAnnihNumb
      integer MergeNumb
      integer NonMergeNumb

!       print *, 'START MERGE'
!       print *, 'creation = ', creation(:)
!       print *, 'annihilation = ', annihilation(:)
!       print *, 'stack = ', stack(:)
       phase = 1.00D0
       MergeNumb = 0
       NonMergeNumb = 0
       counter_position = 1

       do j = 1, stack_dim
        counter = 0
         do i = 1, creation_dim
           if (creation(i) .eq. stack(j)) then
             multiplier = (-1) ** (NonZeroAnnihNumb + creation_dim + NonMergeNumb - i)
             phase = phase * real(multiplier)
!             print *, 'creation_dim, NonZeroAnnihNumb, NonMergeNumb, i, mult = ', creation_dim, NonZeroAnnihNumb, NonMergeNumb, & 
!                                                                                  i, multiplier
             MergeNumb = MergeNumb + 1
             counter = counter + 1
             creation(i) = 0
             cycle
           endif
         enddo
           if (counter .eq. 0) then
            annihilation(counter_position) = stack(j)
!            multiplier = (-1) ** (creation_dim - MergeNumb)
!            phase = phase * real(multiplier)
            NonMergeNumb = NonMergeNumb + 1
!            print *, 'creation_dim - MergeNumb, phase = ', multiplier, phase
            counter_position = counter_position + 1
           endif
       enddo
!       print *, 'AFTER OPERATION'
!       print *, 'creation = ', creation(:)
!       print *, 'annihilation = ', annihilation(:)
!       print *, 'phase = ', phase
!       print *, 'END MERGE'
       end subroutine

!      subroutine sorting (L_dim, L, last_annihilation, phase)
!      integer i, j
!      integer L_dim
!      integer L(L_dim)
!      integer last_annihilation
 !n, s, x, k, L(n), dim
!      k = 0
!      s = 1
!      dim = n - count(L == 0)
!              do i = 1, dim -  1
!                   do j = i + 1, dim
!                        if (L(i) == L(j)) then
!                             k = k + 1
!                         elseif (L(i) > L(j)) then
!                             s = s * (-1)
!                             x = L(i)
!                             L(i) = L(j)
!                             L(j) = x
!                        endif
!                   enddo
!              enddo
!                if (k /= 0) then
!                    s = 0
!                endif
!      end subroutine
