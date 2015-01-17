      subroutine mediator (OrbNumb, ElemOpNumb, MaxSpin, PathNumb, length, S, M, & 
                           NonZeroNumb, StartNumb, Gap, DestNumb, Compliance, path, q_val, k_q_val, & 
                           coef_q_val, coef_q_val_config)
      use spin_tensor

      implicit none
      integer i, j, k
      integer counter_front, counter_back, counter_up, counter_down, counter_gap, counter_k_q_val
      integer OrbNumb, PathNumb, StartNumb, ElemOpNumb, NonZeroNumb, CreatZeroNumb, AnnihZeroNumb
      integer n, S, M, length, MaxSpin
      integer Gap, DestNumb
      integer Compliance(OrbNumb) 
      integer path(MaxSpin)
      integer a(MaxSpin)
      integer SpinOrbPhysVacMod(2*OrbNumb)
      integer stack(Gap), k_q_val(length * PathNumb)
      real*8 phase, phase_merge, phase_order, interphase
      real*8 coef_q_val(length * PathNumb), dublicate_coef_q_val(length)
      real*8 coef_q_val_config(length * PathNumb)

      integer q_val(length * PathNumb, ElemOpNumb)
      integer, target :: dublicate_q_val(length, ElemOpNumb)
      integer, pointer :: semi_dublicate_q_val(:, :)
      integer, pointer :: semi_dublicate_q_val_creation(:), semi_dublicate_q_val_annihilation(:)

!      print *, 'INSIDE'
!      print *, 'MaxSpin = ', MaxSpin

      dublicate_q_val = 0
      dublicate_coef_q_val = 1.00D0
      SpinOrbPhysVacMod = 0

      counter_front = 0
      counter_back = 0
      counter_up = 0
      counter_down = 0
      counter_gap = 0
      NonZeroNumb = 0
      phase = 1.00D0
      do i = 1, OrbNumb
        if (Compliance(i) .eq. 1) then
          dublicate_q_val(:, 1 + counter_front) = 2 * i - 1
          dublicate_q_val(:, 2 + counter_front) = 2 * i
          counter_front = counter_front + 2
          SpinOrbPhysVacMod(2 * i - 1) = 1
          SpinOrbPhysVacMod(2 * i) = 1
         elseif (Compliance(i) .eq. 2) then
          dublicate_q_val(:, ElemOpNumb - counter_back - 1) = 2 * i - 1
          dublicate_q_val(:, ElemOpNumb - counter_back) = 2 * i
          counter_back = counter_back + 2
         elseif (Compliance(i) .eq. 30) then
          a(1 + counter_down) = i
          counter_down = counter_down + 1
          SpinOrbPhysVacMod(2 * i - 1) = 1
          SpinOrbPhysVacMod(2 * i) = 1  
         elseif (Compliance(i) .eq. 3421) then
          dublicate_q_val(:, ElemOpNumb - counter_back) = 2 * i
          counter_back = counter_back + 1
         elseif (Compliance(i) .eq. 3422) then
          dublicate_q_val(:, ElemOpNumb - counter_back) = 2 * i - 1
          counter_back = counter_back + 1
         elseif (Compliance(i) .eq. 3411) then
          dublicate_q_val(:, 1 + counter_front) = 2 * i - 1
          counter_front = counter_front + 1
         elseif (Compliance(i) .eq. 3412) then
          dublicate_q_val(:, 1 + counter_front) = 2 * i
          counter_front = counter_front + 1
         elseif (Compliance(i) .eq. 341) then
          stack(1 + counter_gap) = 2 * i - 1
          counter_gap = counter_gap + 1
          a(1 + counter_up + DestNumb) = i
          counter_up = counter_up + 1
         elseif (Compliance(i) .eq. 342) then
          stack(1 + counter_gap) = 2 * i
          counter_gap = counter_gap + 1
          a(1 + counter_up + DestNumb) = i
          counter_up = counter_up + 1
         elseif (Compliance(i) .eq. 4) then
          a(1 + counter_up + DestNumb) = i
          counter_up = counter_up + 1
         elseif (Compliance(i) .eq. 10) then
          SpinOrbPhysVacMod(2 * i - 1) = 1
          SpinOrbPhysVacMod(2 * i) = 1
        endif
      enddo

!      print *, 'SpinOrbPhysVacMod = ', SpinOrbPhysVacMod(:)

!      print *, '1'
      if (MaxSpin .gt. 0) then
        n = MaxSpin
        semi_dublicate_q_val => dublicate_q_val(:, ElemOpNumb/2 - DestNumb + 1 : ElemOpNumb/2 - DestNumb + 1 + MaxSpin)
        call spin_tensor_up (S, M, MaxSpin, n, a, path, semi_dublicate_q_val, & 
                              dublicate_coef_q_val, length, DestNumb)
        nullify (semi_dublicate_q_val)
      endif
!      print *, '2'

      if (MaxSpin .eq. 0) then

         NonZeroNumb = 1

         semi_dublicate_q_val_annihilation => dublicate_q_val(1, 1 : ElemOpNumb/2)
         semi_dublicate_q_val_creation => dublicate_q_val(1, ElemOpNumb/2 + 1 : ElemOpNumb)
         call order(ElemOpNumb/2, semi_dublicate_q_val_creation, phase_order)
         call order(ElemOpNumb/2, semi_dublicate_q_val_annihilation, phase_order)
!        
         interphase = 1.00D0
!         call interphase_calculator (ElemOpNumb/2, ElemOpNumb/2, OrbNumb, &  
!                               semi_dublicate_q_val_creation(1 : ElemOpNumb/2), &  
!                               semi_dublicate_q_val_annihilation(1 : ElemOpNumb/2), & 
!                               SpinOrbPhysVacMod, interphase)
         nullify (semi_dublicate_q_val_annihilation)
         nullify (semi_dublicate_q_val_creation)

         counter_k_q_val = 0
         do j = 1, ElemOpNumb
           if (dublicate_q_val(1, j) .ne. 0) then
             counter_k_q_val = counter_k_q_val + 1
             q_val(NonZeroNumb + StartNumb, counter_k_q_val) = dublicate_q_val(1, j)
           endif
         enddo

         k_q_val(NonZeroNumb + StartNumb) = counter_k_q_val
         coef_q_val(NonZeroNumb + StartNumb) = dublicate_coef_q_val(1)
         coef_q_val_config(NonZeroNumb + StartNumb) = dublicate_coef_q_val(1) * interphase

      elseif (MaxSpin .gt. 0) then

!     print *, '**************************************************************'
       do i = 1, length

          if (dublicate_q_val(i, ElemOpNumb/2 - DestNumb + 1) .ne. 0 .and. dublicate_q_val(i, ElemOpNumb/2 + 1) .ne. 0) then

            NonZeroNumb = NonZeroNumb + 1

            semi_dublicate_q_val_annihilation => dublicate_q_val(i, 1 : ElemOpNumb/2)
            semi_dublicate_q_val_creation => dublicate_q_val(i, ElemOpNumb/2 + 1 : ElemOpNumb)
!        print *, 'counter_front, counter_back = ', counter_front, counter_back
!        print *, 'S, M, n, length, DestNumb, PathNumb, ElemOpNumb = ', S, M, n, length, DestNumb, PathNumb, ElemOpNumb
!        print *, 'path = ', path(:)
!        print *, 'a = ', a(:)
!        print *, 'q_val = '
!        do k = 1, PathNumb * length
!        write(*, 50) (q_val(k, j), j = 1, ElemOpNumb)
!        enddo
!        print *, 'dublicate_q_val = '
!        do k = 1, length
!        write(*, 50) (dublicate_q_val(k, j), j = 1, ElemOpNumb)
!        enddo
!        print *, 'semi_dublicate_q_val_creation ='
!        write(*, 50) (semi_dublicate_q_val_creation(j), j = 1, ElemOpNumb/2)
!    50  format(50I6)
!        print *, 'dublicate_coef_q_val = ', dublicate_coef_q_val(:)

            call order(ElemOpNumb/2, semi_dublicate_q_val_creation, phase_order)
            call order(ElemOpNumb/2, semi_dublicate_q_val_annihilation, phase_order)
            CreatZeroNumb = count(semi_dublicate_q_val_creation .eq. 0)
            AnnihZeroNumb = count(semi_dublicate_q_val_annihilation .eq. 0)
!
            interphase = 1.00D0
!            call interphase_calculator (ElemOpNumb/2 - CreatZeroNumb, ElemOpNumb/2 - AnnihZeroNumb, OrbNumb, &  
!                                 semi_dublicate_q_val_creation(1 : ElemOpNumb/2 - CreatZeroNumb), &  
!                                 semi_dublicate_q_val_annihilation(1 : ElemOpNumb/2 - AnnihZeroNumb), & 
!                                 SpinOrbPhysVacMod, interphase)
!            print *, 'Interphase = ', interphase
            
              if (Gap .ne. 0) then
!                 print *, 'BEFORE MERGE semi_dublicate_q_val_annihilation = ', semi_dublicate_q_val_annihilation(:)
!                 print *, 'CreatZeroNumb, AnnihZeroNumb, ElemOpNumb/2 - AnnihZeroNumb', CreatZeroNumb, & 
!                           AnnihZeroNumb, ElemOpNumb/2 - AnnihZeroNumb
                 call merge_spin_tensor (ElemOpNumb/2 - CreatZeroNumb, AnnihZeroNumb, ElemOpNumb/2 - AnnihZeroNumb, &  
                                 Gap, semi_dublicate_q_val_creation(1 : ElemOpNumb/2 - CreatZeroNumb), &  
                                 semi_dublicate_q_val_annihilation(ElemOpNumb/2 - AnnihZeroNumb + 1 : ElemOpNumb/2), &  
                                 stack, phase_merge)
                 call order (ElemOpNumb/2, semi_dublicate_q_val_annihilation, phase_order)
                 dublicate_coef_q_val(i) = dublicate_coef_q_val(i) * phase_merge * phase_order
              endif

            nullify (semi_dublicate_q_val_creation)
            nullify (semi_dublicate_q_val_annihilation)

            counter_k_q_val = 0

!            print *, 'dublicate_q_val = ', dublicate_q_val(i, :)
!            print *, 'NonZeroNumb, StartNumb = ', NonZeroNumb, StartNumb
            do j = 1, ElemOpNumb
              if (dublicate_q_val(i, j) .ne. 0) then
                counter_k_q_val = counter_k_q_val + 1
                q_val(NonZeroNumb + StartNumb, counter_k_q_val) = dublicate_q_val(i, j)
              endif
            enddo
!            print *, 'counter_q_val = ', counter_k_q_val
           
            k_q_val(NonZeroNumb + StartNumb) = counter_k_q_val
            coef_q_val(NonZeroNumb + StartNumb) = dublicate_coef_q_val(i)
            coef_q_val_config(NonZeroNumb + StartNumb) = dublicate_coef_q_val(i) * interphase

          endif

       enddo

      endif

      end subroutine

