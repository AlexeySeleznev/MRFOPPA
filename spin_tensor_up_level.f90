      subroutine spin_tensor_up_level (OrbNumb, ElConfNumb, target_S, PhysVac, FirstOrderIdent, S_curr_final, ConFuncNumb)

!
!     VARIABLES DESCRIPTION
!

      implicit none
      integer i, j, k         !dummy variables
      integer counter         !dummy variables
      integer coun1, coun2    !dummy variables
      integer coun3, coun4    !dummy variables
      integer CurrConf        !dummy variables

      integer OrbNumb         ! Number of orbitals
      integer ElConfNumb      ! Total number of electron configurations
      integer ConFuncNumb     ! Total number of configuration functions 

      integer target_S        ! 2 * Target S^2 eigenvalue for configuration functions to be built
      integer MaxSpin         ! Maximal possible S^2 eigenvalue for configuration functions corresponding to the electron configuration under consideration
      integer ElNumb          ! Number of electrons
      integer ElemOpNumb      ! Maximal number of elementary operators within one complex operator
      integer MaxPathNumb     ! Maximal number of paths
      integer length          ! Maximal number of complex operators within one configuration function
      integer PathNumb        ! Number of paths
      integer Gap             ! Number of non-coupled electrons within Physical Vacuum
      integer DestNumb        ! Number of coupled electrons within Physical Vacuum to be annihilated (2 --> 1)
      integer SingPathNumb, TripPathNumb, QuinPathNumb, SepPathNumb ! Number of paths
      integer DublPathNumb, QuarPathNumb, SextPathNumb, OctPathNumb ! Number of paths      
      
      logical SingVac ! true - physical vacuum is a singlet S^2 eigenfunction
                      ! false -physical vacuum is a dublet S^2 eigenfunction
      logical FirstOrderIdent ! false - model space
                              ! true - subsidiary space
      logical S_curr_final    ! true if S current is equal to S_final
                              ! false if S current is equal to S_initial
      integer :: PhysVac(OrbNumb) ! 1D-Array with Physical Vacuum
      integer :: ConfFunc(OrbNumb) ! 1D-Array with Target configuration function
      integer :: Compliance(OrbNumb) ! 1D-Array with Compliance

      integer, allocatable, dimension(:, :) :: Paths ! 2D-Array with paths
      integer, allocatable, dimension(:, :) :: Singlet_Paths, Triplet_Paths, Quintet_Paths, Septet_Paths ! 2D-Arrays with paths
      integer, allocatable, dimension(:, :) :: Dublet_Paths, Quartet_Paths, Sextet_Paths, Octet_Paths ! 2D-Arrays with paths
      integer, allocatable, dimension(:) :: CurrSValue, NextSValue ! 1D-Arrays with current and next S^2 eigenvalues for different paths

      integer, allocatable, target, dimension(:, :) :: q_sing_val, q_trip_val, q_quin_val, q_sep_val  ! 2D-Array with elementary operators for constructing a spin-tensor operator
      integer, allocatable, dimension(:) :: k_q_sing_val, k_q_trip_val, k_q_quin_val, k_q_sep_val ! 1D-Array with numbers of elementary operators constituting complex operators within a spin-tensor operator
      integer, allocatable, dimension(:) :: cf_q_sing_val, cf_q_trip_val, cf_q_quin_val, cf_q_sep_val ! 1D-Array with numbers of complex operators constituting a spin-tensor operator
      real*8,  allocatable, target, dimension(:) :: coef_q_sing_val, coef_q_trip_val, coef_q_quin_val, coef_q_sep_val ! 1D-Array with coefficients before complex operators within a spin-tensor operator
      real*8,  allocatable, target, dimension(:) :: coef_q_sing_val_config, coef_q_trip_val_config, & 
                                                    coef_q_quin_val_config, coef_q_sep_val_config ! 1D-Array with coefficients before determinants within a configuration function
      integer, allocatable, target, dimension(:, :) :: q_dubl_val, q_quar_val, q_sext_val, q_oct_val  ! 2D-Array with elementary operators for constructing a spin-tensor operator
      integer, allocatable, dimension(:) :: k_q_dubl_val, k_q_quar_val, k_q_sext_val, k_q_oct_val ! 1D-Array with numbers of elementary operators constituting complex operators within a spin-tensor operator
      integer, allocatable, dimension(:) :: cf_q_dubl_val, cf_q_quar_val, cf_q_sext_val, cf_q_oct_val ! 1D-Array with numbers of complex operators constituting a spin-tensor operator
      real*8,  allocatable, target, dimension(:) :: coef_q_dubl_val, coef_q_quar_val, coef_q_sext_val, coef_q_oct_val ! 1D-Array with coefficients before complex operators within a spin-tensor operator
      real*8,  allocatable, target, dimension(:) :: coef_q_dubl_val_config, coef_q_quar_val_config, & 
                                                    coef_q_sext_val_config, coef_q_oct_val_config ! 1D-Array with coefficients before complex operators within a configuration function

      ConFuncNumb = 0

!
!     Number of electrons
!

      ElNumb = 0
        do i = 1, OrbNumb
          ElNumb = ElNumb + abs(PhysVac(i))
        enddo

!     Physical Vacuum: singlet or dublet

        if (ElNumb/2*2 .eq. ElNumb) then    ! even or odd
          SingVac = .true.
        else
          SingVac = .false.
        endif

!
!     Temporary files
!

     if (.not. FirstOrderIdent .and. S_curr_final) then
      open (100, file = 'CF_Q_VAL.temp') 
      open (110, file = 'K_Q_VAL.temp')
      open (120, file = 'Q_VAL.temp')
      open (130, file = 'COEF_Q_VAL.temp')
      open (140, file = 'CONF_FUNC_STRUCTURE')
      write (140, 20) (PhysVac(i), i = 1, OrbNumb)
      open (150, file = 'VAL_EL_CONFIG.temp')

     elseif (.not. FirstOrderIdent .and. .not. S_curr_final) then
      open (100, file = 'CF_Q_VAL_INIT.temp')
      open (110, file = 'K_Q_VAL_INIT.temp')
      open (120, file = 'Q_VAL_INIT.temp')
      open (130, file = 'COEF_Q_VAL_INIT.temp')
      open (140, file = 'CONF_FUNC_STRUCTURE_INIT')
      write (140, 20) (PhysVac(i), i = 1, OrbNumb)
      open (150, file = 'VAL_EL_CONFIG_INIT.temp')

     elseif (FirstOrderIdent) then
      open (100, file = 'CF_Q_SUB.temp')
      open (110, file = 'K_Q_SUB.temp')
      open (120, file = 'Q_SUB.temp')
      open (130, file = 'COEF_Q_SUB.temp')
      open (140, file = 'CONF_FUNC_STRUCTURE_SUB')
      write (140, 20) (PhysVac(i), i = 1, OrbNumb)
      open (150, file = 'SUB_EL_CONFIG.temp')

     endif

!
!     Main cycle
!

    do CurrConf = 1, ElConfNumb

!
!     Electron configuration reading 
!
      write (140, *) ' '
      write (140, *) ' '
      read (150, *) (ConfFunc(i), i = 1, OrbNumb)
      write (140, 25) (ConfFunc(i), i = 1, OrbNumb)
      write (140, *) ' '

     
!     Compliance formation

!     Notations: (physical vacuum ---> configuration function)
!     "10"/"20" - no changes
!       "1" - one orbital to be annihilated
!       "2" - one orbital to be created 
!      "30" - one of two spin-orbitals (alfa or beta) to be annihilated
!  "3421"/"3422" - special case 1(alpha/beta) ---> 2
!   "341"/"342"  - special case 1(alpha/beta) ---> 1
!  "3411"/"3412" - special case 1(alpha/beta) ---> 0
!       "4" - one of two spin-orbitals (alpha or beta) to be created
       
      MaxSpin = 0
      ElemOpNumb = 0
      Gap = 0
      DestNumb = 0
      
       do i = 1, OrbNumb
         if (PhysVac(i) .eq. 2) then
            if (ConfFunc(i) .eq. 2) then
               Compliance(i) = 10
            elseif (ConfFunc(i) .eq. 0) then
               Compliance(i) = 1
               ElemOpNumb = ElemOpNumb + 2
            elseif (ConfFunc(i) .eq. 1) then
               Compliance(i) = 30
               MaxSpin = MaxSpin + 1
               ElemOpNumb = ElemOpNumb + 1
               DestNumb = DestNumb + 1
            endif
          elseif (PhysVac(i) .eq. 1) then
            if (ConfFunc(i) .eq. 2) then
               Compliance(i) = 3421
               ElemOpNumb = ElemOpNumb + 1
            elseif (ConfFunc(i) .eq. 1) then
               Compliance(i) = 341
               ElemOpNumb = ElemOpNumb + 2
               MaxSpin = MaxSpin + 1
               Gap = Gap + 1
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
               MaxSpin = MaxSpin + 1
               Gap = Gap + 1
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
               MaxSpin = MaxSpin + 1
               ElemOpNumb = ElemOpNumb + 1
            elseif (ConfFunc(i) .eq. 0) then
               Compliance(i) = 20
            endif
         endif
       enddo
    
!
!      Path generation
!

       MaxPathNumb = 2**MaxSpin

        if (MaxSpin .ne. 0) then

           allocate (Paths(MaxPathNumb, MaxSpin))
           allocate (CurrSValue(MaxPathNumb), NextSValue(MaxPathNumb))

           CurrSValue = 0
           NextSValue = 0
           PathNumb = 1
           counter = 1
             do i = 1, MaxSpin
              NextSValue = 0
              do j = 1, PathNumb
                if (CurrSValue(j) - 1 .ge. 0 .and. CurrSValue(j) + 1 .le. MaxSpin) then
                  NextSValue(j) = CurrSValue(j) + 1
                  Paths(j, i) = 1
                  counter = counter + 1
                  Paths(counter, i) = -1
                  NextSValue(counter) = CurrSValue(j) - 1
                   if (i > 1) then
                      do k = 1, i - 1
                        Paths(counter, k) = Paths(j, k)
                      enddo
                   endif
                elseif (CurrSValue(j) - 1 .lt. 0 .and. CurrSValue(j) + 1 .le. MaxSpin) then
                  NextSValue(j) = CurrSValue(j) + 1
                  Paths(j, i) = 1
                elseif (CurrSValue(j) - 1 .ge. 0 .and. CurrSValue(j) + 1 .gt. MaxSpin) then
                  NextSValue(j) = CurrSValue(j) - 1
                  Paths(j, i) = -1
                else
                  cycle
                endif
              enddo
              PathNumb = counter
              CurrSValue = NextSValue
             enddo

           deallocate (CurrSValue)

          if (SingVac) then

             SingPathNumb = 0
             TripPathNumb = 0
             QuinPathNumb = 0
             SepPathNumb = 0

              do i = 1, PathNumb
                if (NextSValue(i) .eq. 0 .and. target_S .eq. 0) then
                  SingPathNumb = SingPathNumb + 1
                 elseif (NextSValue(i) .eq. 2 .and. target_S .eq. 2) then
                  TripPathNumb = TripPathNumb + 1
                 elseif (NextSValue(i) .eq. 4 .and. target_S .eq. 4) then
                  QuinPathNumb = QuinPathNumb + 1
                 elseif (NextSValue(i) .eq. 6 .and. target_S .eq. 6) then
                  SepPathNumb = SepPathNumb + 1
                endif
              enddo
            
              if (SingPathNumb .ne. 0) then
                allocate (Singlet_Paths(SingPathNumb, MaxSpin))
                coun1 = 1
              endif
              if (TripPathNumb .ne. 0) then
                allocate (Triplet_Paths(TripPathNumb, MaxSpin))
                coun2 = 1
              endif
              if (QuinPathNumb .ne. 0) then
                allocate (Quintet_Paths(QuinPathNumb, MaxSpin))
                coun3 = 1
              endif

              if (SepPathNumb .ne. 0) then
                allocate (Septet_Paths(SepPathNumb, MaxSpin))
                coun4 = 1
              endif

                do i = 1, PathNumb
                  if (NextSValue(i) .eq. 0 .and. target_S .eq. 0) then
                    do j = 1, MaxSpin
                      Singlet_Paths(coun1, j) = Paths(i, j)
                    enddo
                   coun1 = coun1 + 1
                  elseif (NextSValue(i) .eq. 2 .and. target_S .eq. 2) then
                    do j = 1, MaxSpin
                      Triplet_Paths(coun2, j) = Paths(i, j)
                    enddo
                   coun2 = coun2 + 1
                  elseif (NextSValue(i) .eq. 4 .and. target_S .eq. 4) then
                    do j = 1, MaxSpin
                      Quintet_Paths(coun3, j) = Paths(i, j)
                    enddo
                   coun3 = coun3 + 1
                  elseif (NextSValue(i) .eq. 6 .and. target_S .eq. 6) then
                    do j = 1, MaxSpin
                      Septet_Paths(coun4, j) = Paths(i, j)
                    enddo
                   coun4 = coun4 + 1
                  endif
               enddo
                
            deallocate (Paths)
            deallocate (NextSValue)
 
          else

            DublPathNumb = 0
            QuarPathNumb = 0
            SextPathNumb = 0
            OctPathNumb = 0

              do i = 1, PathNumb
                if (NextSValue(i) .eq. 1 .and. target_S .eq. 1) then
                  DublPathNumb = DublPathNumb + 1
                 elseif (NextSValue(i) .eq. 3 .and. target_S .eq. 3) then
                  QuarPathNumb = QuarPathNumb + 1
                 elseif (NextSValue(i) .eq. 5 .and. target_S .eq. 5) then
                  SextPathNumb = SextPathNumb + 1
                 elseif (NextSValue(i) .eq. 7 .and. target_S .eq. 7) then
                  OctPathNumb = OctPathNumb + 1
                endif
              enddo

              if (DublPathNumb .ne. 0) then
                allocate (Dublet_Paths(DublPathNumb, MaxSpin))
                coun1 = 1
              endif
              if (QuarPathNumb .ne. 0) then
                allocate (Quartet_Paths(QuarPathNumb, MaxSpin))
                coun2 = 1
              endif
              if (SextPathNumb .ne. 0) then
                allocate (Sextet_Paths(SextPathNumb, MaxSpin))
                coun3 = 1
              endif
              if (OctPathNumb .ne. 0) then
                allocate (Octet_Paths(OctPathNumb, MaxSpin))
                coun4 = 1
              endif

                do i = 1, PathNumb
                  if (NextSValue(i) .eq. 1 .and. target_S .eq. 1) then
                    do j = 1, MaxSpin
                      Dublet_Paths(coun1, j) = Paths(i, j)
                    enddo
                   coun1 = coun1 + 1
                  elseif (NextSValue(i) .eq. 3 .and. target_S .eq. 3) then
                    do j = 1, MaxSpin
                      Quartet_Paths(coun2, j) = Paths(i, j)
                    enddo
                   coun2 = coun2 + 1
                  elseif (NextSValue(i) .eq. 5 .and. target_S .eq. 5) then
                    do j = 1, MaxSpin
                      Sextet_Paths(coun3, j) = Paths(i, j)
                    enddo
                   coun3 = coun3 + 1
                  elseif (NextSValue(i) .eq. 7 .and. target_S .eq. 7) then
                    do j = 1, MaxSpin
                      Octet_Paths(coun4, j) = Paths(i, j)
                    enddo
                   coun4 = coun4 + 1
                  endif
               enddo

            deallocate (Paths)
            deallocate (NextSValue)

          endif

        endif 

!
!       SPIN-TENSOR OPERATOR GENERATION
!

        length = 2**MaxSpin
         
         if (SingVac) then

           if (SingPathNumb .ne. 0 .or. MaxSpin .eq. 0) then

              if (MaxSpin .eq. 0) then
                SingPathNumb = 1
                allocate (Singlet_Paths(1, 0))
              endif

              allocate (q_sing_val(SingPathNumb * length, ElemOpNumb))
              allocate (k_q_sing_val(SingPathNumb * length))
              allocate (cf_q_sing_val(SingPathNumb))
              allocate (coef_q_sing_val(SingPathNumb * length))
              allocate (coef_q_sing_val_config(SingPathNumb * length))

              q_sing_val = 0
              k_q_sing_val = 0
              cf_q_sing_val = 0
              coef_q_sing_val = 1.00D0
              coef_q_sing_val_config = 1.00D0

              coun2 = 0

!              write (*, *) ' '
!              write (*, *) 'Total number of singlet paths: ', SingPathNumb

                 do i = 1, SingPathNumb
                    write (140, *) ' '
                    write (140, 50), ConFuncNumb + 1, (Singlet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, SingPathNumb, length, 0, &   
                         0, coun1, coun2, Gap, DestNumb, Compliance, Singlet_Paths(i, :), & 
                         q_sing_val, k_q_sing_val, coef_q_sing_val, coef_q_sing_val_config)
                    cf_q_sing_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_sing_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_sing_val(i)
!                       write (140, 40) (k_q_sing_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2
                             write (140, 70) coef_q_sing_val(j), (q_sing_val(j, k), k = 1, k_q_sing_val(j)/2)
                             write (140, 80) (q_sing_val(j, k), k = 1 + k_q_sing_val(j)/2, k_q_sing_val(j))
                             write (110, 10) k_q_sing_val(j)
                             write (130, 14) coef_q_sing_val(j)
                               if (k_q_sing_val(j) .gt. 0) then
                                 write (120, 12) (q_sing_val(j, k), k = 1, k_q_sing_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo

              deallocate (q_sing_val)
              deallocate (k_q_sing_val)
              deallocate (cf_q_sing_val)
              deallocate (coef_q_sing_val)
              deallocate (coef_q_sing_val_config)
              deallocate (Singlet_Paths)

           endif

           if (TripPathNumb .ne. 0) then

              allocate (q_trip_val(TripPathNumb * length, ElemOpNumb))
              allocate (k_q_trip_val(TripPathNumb * length))
              allocate (cf_q_trip_val(TripPathNumb))
              allocate (coef_q_trip_val(TripPathNumb * length))
              allocate (coef_q_trip_val_config(TripPathNumb * length))

              q_trip_val = 0
              k_q_trip_val = 0
              cf_q_trip_val = 0
              coef_q_trip_val = 1.00D0
              coef_q_trip_val_config = 1.00D0

              coun2 = 0

!              write (*, *) ' '
!              write (*, *) 'Total number of triplet paths: ', TripPathNumb

                 do i = 1, TripPathNumb
                    write (140, *) ' '
                    write (140, 52), ConFuncNumb + 1, (Triplet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, TripPathNumb, length, 2, &  
                         0, coun1, coun2, Gap, DestNumb, Compliance, Triplet_Paths(i, :), &  
                         q_trip_val, k_q_trip_val, coef_q_trip_val, coef_q_trip_val_config)
                    cf_q_trip_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_trip_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_trip_val(i)
!                       write (140, 40) (k_q_trip_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2
                             write (140, 70) coef_q_trip_val(j), (q_trip_val(j, k), k = 1, k_q_trip_val(j)/2)
                             write (140, 80) (q_trip_val(j, k), k = 1 + k_q_trip_val(j)/2, k_q_trip_val(j))
                             write (110, 10) k_q_trip_val(j)
                             write (130, 14) coef_q_trip_val(j)
                               if (k_q_trip_val(j) .gt. 0) then
                                 write (120, 12) (q_trip_val(j, k), k = 1, k_q_trip_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo

              deallocate (q_trip_val)
              deallocate (k_q_trip_val)
              deallocate (cf_q_trip_val)
              deallocate (coef_q_trip_val)
              deallocate (coef_q_trip_val_config)
              deallocate (Triplet_Paths)

           endif

           if (QuinPathNumb .ne. 0) then

              allocate (q_quin_val(QuinPathNumb * length, ElemOpNumb))
              allocate (k_q_quin_val(QuinPathNumb * length))
              allocate (cf_q_quin_val(QuinPathNumb))
              allocate (coef_q_quin_val(QuinPathNumb * length))
              allocate (coef_q_quin_val_config(QuinPathNumb * length))

              q_quin_val = 0
              k_q_quin_val = 0
              cf_q_quin_val = 0
              coef_q_quin_val = 1.00D0
              coef_q_quin_val_config = 1.00D0

              coun2 = 0
              
!              write (*, *) ' '
!              write (*, *) 'Total number of quintet paths: ', QuinPathNumb

                 do i = 1, QuinPathNumb
                    write (140, *) ' '
                    write (140, 54), ConFuncNumb + 1, (Quintet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, QuinPathNumb, length, 4, &  
                         0, coun1, coun2, Gap, DestNumb, Compliance, Quintet_Paths(i, :), &  
                         q_quin_val, k_q_quin_val, coef_q_quin_val, coef_q_quin_val_config)
                    cf_q_quin_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_quin_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_quin_val(i)
!                       write (140, 40) (k_q_quin_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2
                             write (140, 70) coef_q_quin_val(j), (q_quin_val(j, k), k = 1, k_q_quin_val(j)/2)
                             write (140, 80) (q_quin_val(j, k), k = 1 + k_q_quin_val(j)/2, k_q_quin_val(j))
                             write (110, 10) k_q_quin_val(j)
                             write (130, 14) coef_q_quin_val(j)
                               if (k_q_quin_val(j) .gt. 0) then
                                 write (120, 12) (q_quin_val(j, k), k = 1, k_q_quin_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo
             
              deallocate (q_quin_val)
              deallocate (k_q_quin_val)
              deallocate (cf_q_quin_val)
              deallocate (coef_q_quin_val)
              deallocate (coef_q_quin_val_config)
              deallocate (Quintet_Paths)

           endif

           if (SepPathNumb .ne. 0) then

              allocate (q_sep_val(SepPathNumb * length, ElemOpNumb))
              allocate (k_q_sep_val(SepPathNumb * length))
              allocate (cf_q_sep_val(SepPathNumb))
              allocate (coef_q_sep_val(SepPathNumb * length))
              allocate (coef_q_sep_val_config(SepPathNumb * length))

              q_sep_val = 0
              k_q_sep_val = 0
              cf_q_sep_val = 0
              coef_q_sep_val = 1.00D0
              coef_q_sep_val_config = 1.00D0

              coun2 = 0
              
!              write (*, *) ' '
!              write (*, *) 'Total number of septet paths: ', SepPathNumb

                 do i = 1, SepPathNumb
                    write (140, *) ' '
                    write (140, 56), ConFuncNumb + 1, (Septet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, SepPathNumb, length, 6, &  
                         0, coun1, coun2, Gap, DestNumb, Compliance, Septet_Paths(i, :), &  
                         q_sep_val, k_q_sep_val, coef_q_sep_val, coef_q_sep_val_config)
                    cf_q_sep_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_sep_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_sep_val(i)
!                       write (140, 40) (k_q_sep_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2
                             write (140, 70) coef_q_sep_val(j), (q_sep_val(j, k), k = 1, k_q_sep_val(j)/2)
                             write (140, 80) (q_sep_val(j, k), k = 1 + k_q_sep_val(j)/2, k_q_sep_val(j))
                             write (110, 10) k_q_sep_val(j)
                             write (130, 14) coef_q_sep_val(j)
                               if (k_q_sep_val(j) .gt. 0) then
                                 write (120, 12) (q_sep_val(j, k), k = 1, k_q_sep_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo
             
              deallocate (q_sep_val)
              deallocate (k_q_sep_val)
              deallocate (cf_q_sep_val)
              deallocate (coef_q_sep_val)
              deallocate (coef_q_sep_val_config)
              deallocate (Septet_Paths)

           endif

         else

           if (DublPathNumb .ne. 0) then

              allocate (q_dubl_val(DublPathNumb * length, ElemOpNumb))
              allocate (k_q_dubl_val(DublPathNumb * length))
              allocate (cf_q_dubl_val(DublPathNumb))
              allocate (coef_q_dubl_val(DublPathNumb * length))
              allocate (coef_q_dubl_val_config(DublPathNumb * length))

              q_dubl_val = 0
              k_q_dubl_val = 0
              cf_q_dubl_val = 0
              coef_q_dubl_val = 1.00D0
              coef_q_dubl_val_config = 1.00D0

              coun2 = 0

!              write (*, *) ' '
!              write (*, *) 'Total number of dublet paths: ', DublPathNumb

               do i = 1, DublPathNumb
                    write (140, *) ' '
                    write (140, 51), ConFuncNumb + 1, (Dublet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, DublPathNumb, length, 1, &  
                         1, coun1, coun2, Gap, DestNumb, Compliance, Dublet_Paths(i, :), &  
                         q_dubl_val, k_q_dubl_val, coef_q_dubl_val, coef_q_dubl_val_config)
                    cf_q_dubl_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_dubl_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_dubl_val(i)
!                       write (140, 40) (k_q_dubl_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2 
                             write (140, 70) coef_q_dubl_val(j), (q_dubl_val(j, k), k = 1, k_q_dubl_val(j)/2)
                             write (140, 80) (q_dubl_val(j, k), k = 1 + k_q_dubl_val(j)/2, k_q_dubl_val(j))
                             write (110, 10) k_q_dubl_val(j)
                             write (130, 14) coef_q_dubl_val(j)
                               if (k_q_dubl_val(j) .gt. 0) then
                                 write (120, 12) (q_dubl_val(j, k), k = 1, k_q_dubl_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo

              deallocate (q_dubl_val)
              deallocate (k_q_dubl_val)
              deallocate (cf_q_dubl_val)
              deallocate (coef_q_dubl_val)
              deallocate (coef_q_dubl_val_config)
              deallocate (Dublet_Paths)

           endif

           if (QuarPathNumb .ne. 0) then

              allocate (q_quar_val(QuarPathNumb * length, ElemOpNumb))
              allocate (k_q_quar_val(QuarPathNumb * length))
              allocate (cf_q_quar_val(QuarPathNumb))
              allocate (coef_q_quar_val(QuarPathNumb * length))
              allocate (coef_q_quar_val_config(QuarPathNumb * length))

              q_quar_val = 0
              k_q_quar_val = 0
              cf_q_quar_val = 0
              coef_q_quar_val = 1.00D0
              coef_q_quar_val_config = 1.00D0

              coun2 = 0

!              write (*, *) ' '
!              write (*, *) 'Total number of quartet paths: ', QuarPathNumb

               do i = 1, QuarPathNumb
                    write (140, *) ' '
                    write (140, 53), ConFuncNumb + 1, (Quartet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, QuarPathNumb, length, 3, &  
                         1, coun1, coun2, Gap, DestNumb, Compliance, Quartet_Paths(i, :), &  
                         q_quar_val, k_q_quar_val, coef_q_quar_val, coef_q_quar_val_config)
                    cf_q_quar_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_quar_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_quar_val(i)
!                       write (140, 40) (k_q_quar_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2
                             write (140, 70) coef_q_quar_val(j), (q_quar_val(j, k), k = 1, k_q_quar_val(j)/2)
                             write (140, 80) (q_quar_val(j, k), k = 1 + k_q_quar_val(j)/2, k_q_quar_val(j))
                             write (110, 10) k_q_quar_val(j)
                             write (130, 14) coef_q_quar_val(j)
                               if (k_q_quar_val(j) .gt. 0) then
                                 write (120, 12) (q_quar_val(j, k), k = 1, k_q_quar_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo

              deallocate (q_quar_val)
              deallocate (k_q_quar_val)
              deallocate (cf_q_quar_val)
              deallocate (coef_q_quar_val)
              deallocate (coef_q_quar_val_config)
              deallocate (Quartet_Paths)

           endif

           if (SextPathNumb .ne. 0) then

              allocate (q_sext_val(SextPathNumb * length, ElemOpNumb))
              allocate (k_q_sext_val(SextPathNumb * length))
              allocate (cf_q_sext_val(SextPathNumb))
              allocate (coef_q_sext_val(SextPathNumb * length))
              allocate (coef_q_sext_val_config(SextPathNumb * length))

              q_sext_val = 0
              k_q_sext_val = 0
              cf_q_sext_val = 0
              coef_q_sext_val = 1.00D0
              coef_q_sext_val_config = 1.00D0

              coun2 = 0

!              write (*, *) ' '
!              write (*, *) 'Total number of sextet paths: ', SextPathNumb

               do i = 1, SextPathNumb
                    write (140, *) ' '
                    write (140, 55), ConFuncNumb + 1, (Sextet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, SextPathNumb, length, 5, &  
                         1, coun1, coun2, Gap, DestNumb, Compliance, Sextet_Paths(i, :), &  
                         q_sext_val, k_q_sext_val, coef_q_sext_val, coef_q_sext_val_config)
                    cf_q_sext_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_sext_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_sext_val(i)
!                       write (140, 40) (k_q_sext_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2
                             write (140, 70) coef_q_sext_val(j), (q_sext_val(j, k), k = 1, k_q_sext_val(j)/2)
                             write (140, 80) (q_sext_val(j, k), k = 1 + k_q_sext_val(j)/2, k_q_sext_val(j))
                             write (110, 10) k_q_sext_val(j)
                             write (130, 14) coef_q_sext_val(j)
                               if (k_q_sext_val(j) .gt. 0) then
                                 write (120, 12) (q_sext_val(j, k), k = 1, k_q_sext_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo

              deallocate (q_sext_val)
              deallocate (k_q_sext_val)
              deallocate (cf_q_sext_val)
              deallocate (coef_q_sext_val)
              deallocate (coef_q_sext_val_config)
              deallocate (Sextet_Paths)

           endif

           if (OctPathNumb .ne. 0) then

              allocate (q_oct_val(OctPathNumb * length, ElemOpNumb))
              allocate (k_q_oct_val(OctPathNumb * length))
              allocate (cf_q_oct_val(OctPathNumb))
              allocate (coef_q_oct_val(OctPathNumb * length))
              allocate (coef_q_oct_val_config(OctPathNumb * length))

              q_oct_val = 0
              k_q_oct_val = 0
              cf_q_oct_val = 0
              coef_q_oct_val = 1.00D0
              coef_q_oct_val_config = 1.00D0

              coun2 = 0

!              write (*, *) ' '
!              write (*, *) 'Total number of octet paths: ', OctPathNumb

               do i = 1, OctPathNumb
                    write (140, *) ' '
                    write (140, 57), ConFuncNumb + 1, (Octet_Paths(i, j), j = 1, MaxSpin)
                    write (140, *) '--------------------------------------------------------------------------------------------'
                    write (140, 60)
                    call mediator (OrbNumb, ElemOpNumb, MaxSpin, OctPathNumb, length, 7, &  
                         1, coun1, coun2, Gap, DestNumb, Compliance, Octet_Paths(i, :), &  
                         q_oct_val, k_q_oct_val, coef_q_oct_val, coef_q_oct_val_config)
                    cf_q_oct_val(i) = coun1
                     if (coun1 .gt. 0) then
                       ConFuncNumb = ConFuncNumb + 1
                       write (100, 10) cf_q_oct_val(i)
                       write (140, *) ' '
!                       write (140, 30) cf_q_oct_val(i)
!                       write (140, 40) (k_q_oct_val(k), k = 1 + coun2, coun1 + coun2)
                       write (140, *) ' '
                          do j = 1 + coun2, coun1 + coun2
                             write (140, 70) coef_q_oct_val(j), (q_oct_val(j, k), k = 1, k_q_oct_val(j)/2)
                             write (140, 80) (q_oct_val(j, k), k = 1 + k_q_oct_val(j)/2, k_q_oct_val(j))
                             write (110, 10) k_q_oct_val(j)
                             write (130, 14) coef_q_oct_val(j)
                               if (k_q_oct_val(j) .gt. 0) then
                                 write (120, 12) (q_oct_val(j, k), k = 1, k_q_oct_val(j))
                               endif
                         enddo
                     endif
                    coun2 = coun2 + coun1
                 enddo

              deallocate (q_oct_val)
              deallocate (k_q_oct_val)
              deallocate (cf_q_oct_val)
              deallocate (coef_q_oct_val)
              deallocate (coef_q_oct_val_config)
              deallocate (Octet_Paths)

           endif

        endif

    enddo

       close (100)
       close (110)
       close (120)
       close (130)
       close (140)
       close (150)

    10 format (I4)
    12 format (90I4)
    14 format (F16.9)

    20 format ('  Physical Vacuum : ', 40I4)
    25 format ('  Target electron configuration : ', 40I4)
!    30 format ('  Number of complex operators :   ', I4)
!    40 format ('  Number of elementary operators in each complex operator : ', 50I4)
    50 format ('  Singlet configuration function corresponding to Path # ', I4, ':     ', 50I3)
    51 format ('  Dublet configuration function corresponding to Path #  ', I4, ':     ', 50I3)
    52 format ('  Triplet configuration function corresponding to Path # ', I4, ':     ', 50I3)
    53 format ('  Quartet configuration function corresponding to Path # ', I4, ':     ', 50I3)
    54 format ('  Quintet configuration function corresponding to Path # ', I4, ':     ', 50I3)
    55 format ('  Sextet configuration function corresponding to Path #  ', I4, ':     ', 50I3)
    56 format ('  Septet  configuration function corresponding to Path # ', I4, ':     ', 50I3)
    57 format ('  Octet  configuration function corresponding to Path #  ', I4, ':     ', 50I3)
    60 format ('  Coefficient                      Complex operator   ')
    70 format ('   ', F9.6,'               annihilation part   ', 40I4)
    80 format ('                        ', '       creation part   ', 40I4)
       end subroutine
         



