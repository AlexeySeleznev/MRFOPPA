      subroutine zero_order_EVP (ValConfFuncNumb, InitStateNumb, S_init_final_same, S_curr_final)
!
!      VARIABLES DESCRIPTION
!
      use InitialStateData
      implicit none
      integer i, j, k, l, f                                 !dummy variables
      integer shift_k_I, shift_k_J                          !dummy variables
      integer info
      
      integer ValConfFuncNumb                               ! Total number of valent configuration functions
      integer InitStateNumb                                 ! Serial number of initial state wave function among zero-order EVP solutions
      integer k_q_val_dim, MaxOpNumb                        ! Q_VAL dimensions
      integer SolutNumbToBeShown                            ! Number of solutions to be shown
      logical S_init_final_same                             ! S_init_state = S_final_state --> True
                                                            ! S_init_state != S_final_state --> False
      logical S_curr_final                                  ! S_curr_state = S_final_state --> True
                                                            ! S_curr_state = S_init_state --> False
      integer cf_q_val(ValConfFuncNumb)
      integer, allocatable, dimension(:) :: k_q_val
      real*8, allocatable, dimension(:) :: coef_q_val
      integer, allocatable, dimension(:, :) :: q_val
      real*8, allocatable, dimension(:) :: working_array
      real*8, allocatable, dimension(:) :: EigenVaLModSp    ! Hamiltonian eigenvalues in model space
      real*8 hmod(ValConfFuncNumb, ValConfFuncNumb)         ! Hamiltonian matrix
      real*8 ht                                             ! dummy variables
      real   time1, time2

!
!     INPUT READING
!

      k_q_val_dim = 0
      if (S_curr_final) then
         open (10, file = 'CF_Q_VAL.temp')
      elseif (.not. S_curr_final) then
         open (10, file = 'CF_Q_VAL_INIT.temp')
      endif
        do i = 1, ValConfFuncNumb
          read(10, *) cf_q_val(i)
          k_q_val_dim = k_q_val_dim + cf_q_val(i)
        enddo
      close (10)

      allocate (k_q_val(k_q_val_dim))
      allocate (coef_q_val(k_q_val_dim))

      if (S_curr_final) then
        open (20, file = 'K_Q_VAL.temp')
        open (30, file = 'COEF_Q_VAL.temp')
      elseif (.not. S_curr_final) then
        open (20, file = 'K_Q_VAL_INIT.temp')
        open (30, file = 'COEF_Q_VAL_INIT.temp')
      endif
        do i = 1, k_q_val_dim
          read(20, *) k_q_val(i)
          read(30, *) coef_q_val(i)
        enddo
      close (20)
      close (30)

      MaxOpNumb = MAXVAL(k_q_val)
      allocate (q_val(k_q_val_dim, MaxOpNumb))
      q_val = 0

      if (S_curr_final) then
        open (40, file = 'Q_VAL.temp')
      elseif (.not. S_curr_final) then
        open (40, file = 'Q_VAL_INIT.temp')
      endif
        do i = 1, k_q_val_dim
          if (k_q_val(i) .ne. 0) then
            read(40, *) (q_val(i, j), j = 1, k_q_val(i))
          endif
        enddo
      close (40)

     call CPU_TIME(time1)

      hmod = 0.00D0

!
!     GENERATION OF MATRIX ELEMENTS
!

!
!         external cycles
!
      shift_k_I = 0

        do i = 1, ValConfFuncNumb
          shift_k_J = 0

            do j = 1, i
!
!         inner cycles 
!
               do k = 1 + shift_k_I, shift_k_I + cf_q_val(i)
                  do l = 1 + shift_k_J, shift_k_J + cf_q_val(j)

                     if (k_q_val(k) .eq. 0 .and. k_q_val(l) .eq. 0) then
                       call V_element (k_q_val(k), k_q_val(l), 0, 0, & 
                                         MaxOpNumb, ht)
                       hmod(i, j) = hmod(i, j) + ht * coef_q_val(k) * coef_q_val(l)

                     elseif (k_q_val(k) .ne. 0 .and. k_q_val(l) .eq. 0) then
                       call V_element (k_q_val(k), k_q_val(l), q_val(k, :), 0, & 
                                         MaxOpNumb, ht)
                       hmod(i, j) = hmod(i, j) + ht * coef_q_val(k) * coef_q_val(l)

                     elseif (k_q_val(k) .eq. 0 .and. k_q_val(l) .ne. 0) then
                       call V_element (k_q_val(k), k_q_val(l), 0, q_val(l, :), & 
                                         MaxOpNumb, ht)
                       hmod(i, j) = hmod(i, j) + ht * coef_q_val(k) * coef_q_val(l)

                     elseif (k_q_val(k) .ne. 0 .and. k_q_val(l) .ne. 0) then
                       call V_element (k_q_val(k), k_q_val(l), q_val(k, :), q_val(l, :), & 
                                         MaxOpNumb, ht)
                       hmod(i, j) = hmod(i, j) + ht * coef_q_val(k) * coef_q_val(l)
                     endif

                  enddo
               enddo
!
!         end of inner cycles
!
              shift_k_J = shift_k_J + cf_q_val(j)
            enddo
          shift_k_I = shift_k_I + cf_q_val(i)
        enddo
!
!         end of external cycles
!

      call CPU_TIME(time2)

      write (*, *) ' '
      write (*, 50) time2 - time1

!
!     EVP SOLUTION
!

      call CPU_TIME(time1)

      allocate (EigenValModSp (ValConfFuncNumb))
      allocate (working_array (3 * ValConfFuncNumb))

      EigenValModSp = 0
      working_array = 0
!
!     LAPACK subroutine DSYEV
!
      call dsyev ( 'V', 'L', ValConfFuncNumb, hmod, ValConfFuncNumb, EigenValModSp, & 
                    working_array, 3 * ValConfFuncNumb, info)

      deallocate(working_array)

      call CPU_TIME(time2)
      write (*, *) ' '
      write (*, 60) time2 - time1

!
!     SOLUTION DISPLAY AND INITIAL STATE WAVE FUNCTION FORMATION
!

      if (.not. S_init_final_same .and. .not. S_curr_final .or. S_init_final_same .and. S_curr_final) then

!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
      InitExOpNumb = 0                 ! variable is defined in the InitialStateData module
      allocate (cf_q_init_logical(ValConfFuncNumb))         ! logical array is defined in the InitialStateData module
      cf_q_init_logical = .false.
      open (40, file = 'Zero_Order_Initial_State.temp')
         do j = 1, ValConfFuncNumb
            if (abs(hmod(j, InitStateNumb)) .gt. 0.000001) then
              write (40, *) j, hmod(j, InitStateNumb)
              InitExOpNumb = InitExOpNumb + cf_q_val(j)
              cf_q_init_logical(j) = .true.
            endif
         enddo
      close(40)

      allocate (k_q_init(InitExOpNumb))                    !
      allocate (ksi(InitExOpNumb))                         ! Arrays are defined in the InitialStateData module
      allocate (q_init(InitExOpNumb, MaxOpNumb))           !
      MaxInitDim = MaxOpNumb
      l = 0
      k = 0
         do j = 1, ValConfFuncNumb
            if (cf_q_init_logical(j)) then
               do i = 1 + l, cf_q_val(j) + l
                  k = k + 1
                  k_q_init(k) = k_q_val(i)
                  ksi(k) = coef_q_val(i) * hmod(j, InitStateNumb)
                    if (k_q_val(i) .ne. 0) then
                       do f = 1, k_q_val(i)
                           q_init(k, f) = q_val(i, f)
                       enddo
                    endif
               enddo
               l = l + cf_q_val(j)
            else
               l = l + cf_q_val(j)
            endif
         enddo
!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!

      open (45, file = 'Hmod.temp', FORM = 'UNFORMATTED')
         do i = 1, ValConfFuncNumb
            write (45) (hmod(i, j), j = 1, ValConfFuncNumb)
         enddo
      close (45)

      endif

      if (.not. S_curr_final) then
         if (InitStateNumb .le. 5 .and. ValConfFuncNumb .gt. 5) then
            SolutNumbToBeShown = 5
         else
            SolutNumbToBeShown = InitStateNumb
         endif
      else
            SolutNumbToBeShown = ValConfFuncNumb
      endif

      do i = 1, SolutNumbToBeShown
         write(*, *) ' '
         write(*, 70) i, EigenValModSp(i)
         write(*, *) '--------------------------------------------------------------------'
           do j = 1, ValConfFuncNumb
              if (abs(hmod(j, i)) .ge. 0.1) then
                 write (*, 80) hmod(j, i), j
              endif
           enddo
      enddo
      deallocate (k_q_val)
      deallocate (coef_q_val)
      deallocate (q_val)

   50 format ('      ..EXECUTION TIME OF MATRIX ELEMENT GENERATION......', F9.6)
   60 format ('      ..EXECUTION TIME OF EVP SOLUTION...................', F9.6)
   70 format ('Eigenvector #  ', I4, '    with eigenvalue     ', F15.9)
   80 format ('   ', F9.6, I5)
      end subroutine
