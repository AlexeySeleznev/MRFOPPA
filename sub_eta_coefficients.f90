      subroutine sub_eta_coefficients (OrbNumb, ElConfNumb, PhysVac, S_final_state, SubExOpNumb)
      use elconfig
      implicit none
      integer i, j, k, l                                  ! dummy variables
      integer k_q_sub
      integer OrbNumb
      integer ElConfNumb
      integer ElemOpNumb
      integer SubExOpNumb
      integer SingElNumb
      integer S_final_state
      integer PhysVac(OrbNumb)
      integer ConfFunc(OrbNumb)
      integer Compliance(OrbNumb)
      integer Spin_Bufer
      integer sum_odd_creation, sum_odd_annihilation
      integer sum_even_creation, sum_even_annihilation
      real*8  phase
      real*8  eta

      integer, allocatable, target, dimension(:, :) :: operators
      integer, pointer, dimension(:, :) :: operators_creation, operators_annihilation
      integer, pointer, dimension(:) :: operators_creation_1D, operators_annihilation_1D
      integer, allocatable, dimension(:) :: operators_array

      if (S_final_state .ne. 0 .and. S_final_state .ne. 1) then
          open (75, file = 'SUB_EL_CONFIG_ALL.temp')
      else
          open (75, file = 'SUB_EL_CONFIG.temp')
      endif
      open (80, file = 'Q_SUB_ALL.temp')
      open (85, file = 'K_Q_SUB_ALL.temp')
      open (90, file = 'ETA_COEFFICIENTS.temp')

      Spin_Bufer = 0
      SubExOpNumb = 0
      do i = 1, OrbNumb
        if (PhysVac(i) .eq. 1) then
           Spin_Bufer = Spin_Bufer + 1
        elseif (PhysVac(i) .eq. -1) then
           Spin_Bufer = Spin_Bufer - 1
        endif
      enddo

         do l = 1, ElConfNumb
            ElemOpNumb = 0
            SingElNumb = 0
            read (75, *) (ConfFunc(j), j = 1, OrbNumb)
            do i = 1, OrbNumb
               if (PhysVac(i) .eq. 2) then
                 if (ConfFunc(i) .eq. 2) then
                    Compliance(i) = 0
                 elseif (ConfFunc(i) .eq. 0) then
                    Compliance(i) = 1
                    ElemOpNumb = ElemOpNumb + 2
                 elseif (ConfFunc(i) .eq. 1) then
                    Compliance(i) = 30
                    ElemOpNumb = ElemOpNumb + 1
                    SingElNumb = SingElNumb + 1
                 endif
               elseif (PhysVac(i) .eq. 1) then
                 if (ConfFunc(i) .eq. 2) then
                    Compliance(i) = 3421
                    ElemOpNumb = ElemOpNumb + 1
                 elseif (ConfFunc(i) .eq. 1) then
                    Compliance(i) = 341
                    ElemOpNumb = ElemOpNumb + 2
                    SingElNumb = SingElNumb + 1
                 elseif (ConfFunc(i) .eq. 0) then
                    Compliance(i) = 3411
                    ElemOpNumb = ElemOpNumb + 1
                 endif
               elseif (PhysVac(i) .eq. -1) then
                 if (ConfFunc(i) .eq. 2) then
                    Compliance(i) = 3422
                    ElemOpNumb = ElemOpNumb + 1
                 elseif (ConfFunc(i) .eq. 1) then
                    Compliance(i) = 342
                    ElemOpNumb = ElemOpNumb + 2
                    SingElNumb = SingElNumb + 1
                 elseif (ConfFunc(i) .eq. 0) then
                    Compliance(i) = 3412
                    ElemOpNumb = ElemOpNumb + 1
                 endif
               elseif (PhysVac(i) .eq. 0) then
                 if (ConfFunc(i) .eq. 2) then
                    Compliance(i) = 2
                    ElemOpNumb = ElemOpNumb + 2
                 elseif (ConfFunc(i) .eq. 1) then
                    Compliance(i) = 4
                    SingElNumb = SingElNumb + 1
                    ElemOpNumb = ElemOpNumb + 1
                 elseif (ConfFunc(i) .eq. 0) then
                    Compliance(i) = 0
                 endif
               endif
            enddo
            allocate (operators(2 ** SingElNumb, ElemOpNumb))
            operators = 0
            operators_annihilation => operators(:, 1 : ElemOpNumb/2)
            operators_creation => operators(:, ElemOpNumb/2 + 1 : ElemOpNumb)
            call operator_constructor(OrbNumb, Compliance, operators_creation, & 
                                      operators_annihilation, ElemOpNumb/2, ElemOpNumb/2, & 
                                      ElemOpNumb/2, 2 ** SingElNumb, OrbNumb)
            nullify (operators_creation, operators_annihilation)
            do i = 1, 2 ** SingElNumb
                operators_annihilation_1D => operators(i, 1 : ElemOpNumb/2)
                operators_creation_1D => operators(i, ElemOpNumb/2 + 1 : ElemOpNumb)
                call order(ElemOpNumb/2, operators_annihilation_1D, phase)
                call order(ElemOpNumb/2, operators_creation_1D, phase)
                sum_odd_creation = 0
                sum_odd_annihilation = 0
                sum_even_creation = 0
                sum_even_annihilation = 0
                do j = 1, ElemOpNumb/2
                  if (operators_creation_1D(j) .ne. 0 .and. (operators_creation_1D(j)/2) * 2 .ne. & 
                                                                            operators_creation_1D(j)) then
                    sum_odd_creation = sum_odd_creation + 1
                  endif
                  if (operators_creation_1D(j) .ne. 0 .and. (operators_creation_1D(j)/2) * 2 .eq. & 
                                                                            operators_creation_1D(j)) then
                    sum_even_creation = sum_even_creation + 1
                  endif
                  if (operators_annihilation_1D(j) .ne. 0 .and. (operators_annihilation_1D(j)/2) * 2 .ne. & 
                                                                              operators_annihilation_1D(j)) then
                    sum_odd_annihilation = sum_odd_annihilation + 1
                  endif
                  if (operators_annihilation_1D(j) .ne. 0 .and. (operators_annihilation_1D(j)/2) * 2 .eq. & 
                                                                              operators_annihilation_1D(j)) then
                    sum_even_annihilation = sum_even_annihilation + 1
                  endif
                enddo
                nullify (operators_creation_1D, operators_annihilation_1D)
                if (sum_odd_creation + sum_even_annihilation + Spin_Bufer .eq. sum_odd_annihilation + & 
                                                                                        sum_even_creation) then
                    allocate (operators_array(ElemOpNumb))
                    k_q_sub = 0
                    operators_array = 0
                    do j = 1, ElemOpNumb
                       if (operators(i, j) .ne. 0) then
                          k_q_sub = k_q_sub + 1
                          operators_array(k_q_sub) = operators(i, j)
                       endif
                    enddo
                    write (85, 100) k_q_sub
                    write (80, 100) (operators_array(j), j = 1, k_q_sub)
                    call eta_coefficients_calculator (ElemOpNumb, k_q_sub, operators_array, eta)
                    write (90, *) eta
                    deallocate (operators_array)
                    SubExOpNumb = SubExOpNumb + 1
                endif    
            enddo
            deallocate (operators)
         enddo     
      close(75)
      close(80)
      close(85)
      close(90)
  100 format(900I4)
      end subroutine

      subroutine eta_coefficients_calculator (ElemOpNumb, k_q_sub, operators, eta)
      use integrals_amplitudes
      use InitialStateData
      implicit none
      integer i, j
      integer ElemOpNumb
      integer MaxOpNumb
      integer k_q_sub
      integer operators(ElemOpNumb)
      real*8 eta
      real*8 Em, ht
      integer, allocatable, dimension(:) :: q_I, q_J
      

      eta = 0.00D0
      MaxOpNumb = max(MaxInitDim, k_q_sub)
      allocate (q_I(MaxOpNumb), q_J(MaxOpNumb))
      q_I = 0
      q_J = 0

      Em = ZeroOrderVacEner
      do j = 1, k_q_sub/2
         Em = Em + amplitudes((operators(k_q_sub/2 + j) + 1)/2) - amplitudes((operators(j) + 1)/2)
      enddo
      do i = 1, k_q_sub
         q_I(i) = operators(i)
      enddo

      do i = 1, InitExOpNumb
         do j = 1, k_q_init(i)
            q_J(j) = q_init(i, j)
         enddo
         call V_element (k_q_sub, k_q_init(i), q_I, q_J, MaxOpNumb, ht)
         eta = eta + ht * ksi(i)
      enddo

      eta = eta / (MeanEner - Em)

      deallocate (q_I, q_J)
      end subroutine





