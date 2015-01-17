!
!          VARIABLES DESCRIPTION
!
      use integrals_amplitudes
      implicit none
      integer i, j, k           !dummy variables

      integer ElNumb            ! Number of electrons

      integer OrbNumb           ! Number of orbitals
      integer CoreNumb          ! Number of core orbitals
      integer ValNumb           ! Number of valent orbitals
      integer OccValNumb        ! Number of valent spin-orbitals occupied in a physical vacuum
      integer NonOccValNumb     ! Number of valent spin-orbitals not occupied in a physical vacuum
      integer VirtNumb          ! Number of virtual orbitals
      integer FzrNumb           ! Number of frozen core orbitals

      integer ValElConfNumb     ! Number of valent electron configurations
      integer ValElConfInitNumb ! Number of valent electron configurations for zero-order initial state eigenvalue problem generation
      integer ValConfFuncNumb   ! Number of valent configuration functions
      integer ValConfFuncInitNumb ! Number of valent configuration functions for zero-order initial state eigenvalue problem generation
      integer SubElConfNumb     ! Number of electron configurations for subsidiary space
      integer SubConfFuncNumb   ! Number of configuration functions for subsidiary space
      integer AllSubElConfNumb  ! Total number of electron configurations for subsidiary space
      integer SubExOpNumb       ! Total number of excitation operators in subsidiary space

      integer S_init_state      ! 2 x S for initial state: 0 -- singlet, 1 -- dublet, 2 -- triplet, etc.
      integer S_final_state     ! 2 x S for final state: 0 -- singlet, 1 -- dublet, 2 -- triplet, etc.

      logical S_init_final_same ! S_init_state = S_final_state --> true
                                ! S_init_state != S_final_state --> false
      logical S_curr_final      ! if S current is equal to S_final, then --> true
                                ! if S current is equal to S_initial, then --> false
      logical SingVac           ! Even total number of electrons --> true
                                ! Odd total number of electrons --> false
      logical FirstOrder_On_Off ! Calculations within first-order MRPT --> true (user defined)
                                ! Calculations only within zero-order MRPT --> false (user defined)
      logical FirstOrderIdent   ! Calculation within first-order MRPT --> true
                                ! Calculation within zero-order MRPT --> false
      logical SubAddIdent       ! Calculation of all electron configurations in subsidiary space is required --> true
                                ! Calculation of all electron configurations in subsidiary space is not required --> false
      integer InitStateNumb     ! Serial number of the state in zero-order eigenvalue problem output to be taken as initial state

      integer, allocatable, dimension(:) :: PhysVac   ! 1D -Array with Physical Vacuum Structure
      
      real time1, time2         ! CPU time points
      
!**********************************************************************************************
!********************************* INTRODUCTION ***********************************************
!**********************************************************************************************

      call introduction ()

!**********************************************************************************************
!********************************* INPUT READING **********************************************
!**********************************************************************************************

      write (*, *) '**********************************************************************************************'
      write (*, *) '********************************* INPUT READING **********************************************'
      write (*, *) '**********************************************************************************************'

      write (*, *) ' '
      write (*, *) ' ...Input data reading... '
      write (*, *) ' '
      open (100, file = 'INPUT')

      read (100, *) CoreNumb, ValNumb, VirtNumb, FzrNumb
      write (*, 11) CoreNumb
      write (*, 12) ValNumb
      write (*, 13) VirtNumb
      write (*, 14) FzrNumb
      write (*, *) ' '

      read (100, *) S_init_state
      write (*, 15) S_init_state + 1

      read (100, *) S_final_state
      write (*, 16) S_final_state + 1
      write (*, *) ' ' 

        if (S_init_state .eq. S_final_state) then
           S_init_final_same = .true.
        else
           S_init_final_same = .false.
        endif

      OrbNumb = CoreNumb + ValNumb + VirtNumb
      SpOrbNumb = 2 * OrbNumb

      allocate (PhysVac(OrbNumb))
      read (100, *) (PhysVac(i), i = 1, OrbNumb)
      write (*, 17) (PhysVac(i), i = 1, OrbNumb)
      write (*, *) ' '
    
      OccValNumb = 0
        do i = CoreNumb + 1, CoreNumb + ValNumb 
           OccValNumb = OccValNumb + abs(PhysVac(i))
        enddo
      NonOccValNumb = 2 * ValNumb - OccValNumb

        if (OccValNumb/2*2 .eq. OccValNumb) then
           SingVac = .true.
         else
           SingVac = .false.
        endif

      OnePartIntNumb = OrbNumb * (OrbNumb + 1) / 2                                       ! This integer type variable is defined in the integrals_amplitudes module
      TwoPartIntNumb = (OrbNumb**4 + 2 * OrbNumb**3 + 3 * OrbNumb**2 + 2 * OrbNumb) / 8  ! This integer type variable is defined in the integrals_amplitudes module
      write (*, 18) OnePartIntNumb
      write (*, 19) TwoPartIntNumb
      call integrals_amplitudes_read (OrbNumb)

      read (100, *) MeanEner                 ! This real*8 type variable is defined in the integtals_amplitudes module
      write (*, 20) MeanEner

      read (100, *) NuclRepEner              ! This real*8 type variable is defined in the integrals_amplitudes module
      write (*, 21) NuclRepEner

      read (100, *) FirstOrder_ON_OFF

      read (100, *) InitStateNumb
      write (*, 22) InitStateNumb

      call PhysVacEnerCalc (OrbNumb, CoreNumb + ValNumb, PhysVac)
      call ZeroOrderVacuumEnergy (OrbNumb, PhysVac)
      call integral_transformation (OrbNumb, CoreNumb + ValNumb, PhysVac)

      close (100)

      write (*, *) ' '
      write (*, *) ' ...DONE'
      write (*, *) ' '


!**********************************************************************************************
!********************************** ZERO ORDER ***********************************************
!**********************************************************************************************

!
!     GENERATION OF ZERO-ORDER INITIAL STATE WAVE FUNCTION (IF NECESSARY)
!

      if (.not. S_init_final_same) then

      write (*, *) '**********************************************************************************************'
      write (*, *) '************************** ZERO ORDER (INITIAL STATE) ****************************************'
      write (*, *) '**********************************************************************************************'

      S_curr_final = .false.

      write (*, *) ' '
!
!          ELECTRON CONFIGURATION GENERATION (MODEL SPACE) 
!
      write (*, *) ' '
      write (*, *) ' ...Valent electron configuration generation for zero-order initial state eigenvalue problem... '
      write (*, *) ' '
      call CPU_TIME(time1)

      call val_elconf_generator(OrbNumb, CoreNumb, ValNumb, ElNumb, OccValNumb, & 
                                S_init_state, S_curr_final, ValElConfInitNumb)

      call CPU_TIME(time2)
      write (*, *) ' ...DONE '
      write (*, *) ' '
      write (*, 10) time2 - time1
      write (*, *) ' '
      write (*, 30) ValElConfInitNumb
      write (*, *) ' '
      write (*, *) ' '

!
!     all the configurations generated have been saved in a temporary file "Val_EL_CONFIG_INIT.temp",
!     which will be removed after the end of the program
!
!     the total number of configurations amounts to "ValElConfInitNumb"
!

!
!          SPIN-TENSOR OPERATOR GENERATION (MODEL SPACE)
!

      FirstOrderIdent = .false.
      write (*, *) ' '
      write (*, *) ' ...Valent configuration function generation for zero-order initial state eigenvalue problem... '
      write (*, *) ' '
      call CPU_TIME(time1)

      call spin_tensor_up_level(OrbNumb, ValElConfInitNumb, S_init_state, PhysVac, & 
                                   FirstOrderIdent, S_curr_final, ValConfFuncInitNumb)

      call CPU_TIME(time2)
      write (*, *) ' ...DONE '
      write (*, *) ' '
      write (*, 10) time2 - time1
      write (*, *) ' '
      write (*, 31) ValConfFuncInitNumb
      write (*, *) ' '
      write (*, *) ' '

!
!    all the spin-tensor operators generated have been saved in temporary unformatted files "Q_VAL_INIT.temp",
!    "K_Q_VAL_INIT.temp", "CF_Q_VAL_INIT.temp", "COEF_Q_VAL_INIT.temp". For user conveniency, the formatted file "CONF_FUNC_STRUCTURE_INIT",
!    containing all the above mentioned data, was created. All .temp files will be removed after the end of the program
!

!
!           ZERO-ORDER EIGENVALUE PROBLEM (MODEL SPACE)
!

      write (*, *) ' '
      write (*, *) ' ...Zero-Order Initial State Eigenvalue Problem Generation and & 
                    Solution... '

      call zero_order_EVP (ValConfFuncInitNumb, InitStateNumb, S_init_final_same, S_curr_final)

      write (*, *) ' '
      write (*, *) ' '

!
!     the initial state wave function structure has been saved in temporary unformatted files "Q_VAL_INIT.temp",
!     "K_Q_VAL_INIT.temp", "COEF_Q_VAL_INIT.temp".
!
      endif

!********************************************************************************************************

      write (*, *) '**********************************************************************************************'
      write (*, *) '*************************** ZERO ORDER (FINAL STATE) *****************************************'
      write (*, *) '**********************************************************************************************'

!
!          ELECTRON CONFIGURATION GENERATION (MODEL SPACE)
!

      S_curr_final = .true.
      write (*, *) ' '
      write (*, *) ' ...Valent electron configuration generation... '
      write (*, *) ' '
      call CPU_TIME(time1)

      call val_elconf_generator(OrbNumb, CoreNumb, ValNumb, ElNumb, OccValNumb, & 
                                S_final_state, S_curr_final, ValElConfNumb)

      call CPU_TIME(time2)
      write (*, *) ' ...DONE '
      write (*, *) ' '
      write (*, 10) time2 - time1
      write (*, *) ' '
      write (*, 30) ValElConfNumb
      write (*, *) ' '
      write (*, *) ' '


!
!     all the configurations generated have been saved in a temporary file "Val_EL_CONFIG.temp",
!     which is deleted after execution of the program
!
!     the total number of configurations amounts to "ValElConfNumb"
!

!
!          SPIN-TENSOR OPERATOR GENERATION (MODEL SPACE)
!

      FirstOrderIdent = .false.
      write (*, *) ' '
      write (*, *) ' ...Valent configuration function generation... '
      write (*, *) ' '
      call CPU_TIME(time1)

      call spin_tensor_up_level(OrbNumb, ValElConfNumb, S_final_state, PhysVac, FirstOrderIdent, S_curr_final, ValConfFuncNumb)

      call CPU_TIME(time2)
      write (*, *) ' ...DONE '
      write (*, *) ' '
      write (*, 10) time2 - time1
      write (*, *) ' '
      write (*, 31) ValConfFuncNumb
      write (*, *) ' '
      write (*, *) ' '

!
!     all the spin-tensor operators generated have been saved in temporary unformatted files "Q_VAL.temp",
!     "K_Q_VAL.temp", "CF_Q_VAL.temp", "COEF_Q_VAL.temp". For user conveniency, the formatted file "CONF_FUNC_STRUCTURE",
!     containing all the above mentioned data, was created. All .temp files will be removed after the end of the program
!

!
!           ZERO-ORDER EIGENVALUE PROBLEM (MODEL SPACE)
!

      write (*, *) ' '
      write (*, *) ' ...Zero-Order Eigenvalue Problem Generation and & 
                    Solution... '

      call zero_order_EVP (ValConfFuncNumb, InitStateNumb, S_init_final_same, S_curr_final)

      write (*, *) ' '
      write (*, *) ' '

!
!     the initial state wave function structure has been saved in temporary unformatted files "Q_VAL_INIT.temp",
!     "K_Q_VAL_INIT.temp", "COEF_Q_VAL_INIT.temp".
!

!**********************************************************************************************
!********************************** FIRST ORDER ***********************************************
!**********************************************************************************************

      if (FirstOrder_ON_OFF) then

      write (*, *) '**********************************************************************************************'
      write (*, *) '********************************* FIRST ORDER ************************************************'
      write (*, *) '**********************************************************************************************'

!
!          ELECTRON CONFIGURATION GENERATION (SUBSIDIARY SPACE)
!

      write (*, *) ' '
      write (*, *) ' ...Subsidiary electron configuration generation... '
      write (*, *) ' '
      SubAddIdent = .false.
      call CPU_TIME(time1)

      call sub_elconf_generator(OrbNumb, CoreNumb, ValNumb, VirtNumb, FzrNumb, & 
                                      ElNumb, OccValNumb, S_final_state, SubAddIdent, SubElConfNumb)
      call CPU_TIME(time2)
      write (*, *) ' ...DONE '
      write (*, *) ' '
      write (*, 10) time2 - time1
      write (*, *) ' '
      write (*, 32) SubElConfNumb
      write (*, *) ' '
      write (*, *) ' '

!
!     all the configurations generated have been saved in a temporary file "SUB_EL_CONFIG.temp",
!     which is deleted after execution of the program
!
!     the total number of configurations amounts to "SubElConfNumb"
!

!
!          ETA COEFFICIENT GENERATION (SUBSIDIARY SPACE)
!

      write (*, *) ' '
      write (*, *) ' ...Eta coefficient generation... '
      write (*, *) ' '
      call CPU_TIME(time1)

         if (S_final_state .ne. 0 .and. S_final_state .ne. 1) then
            SubAddIdent = .true.
            if (SingVac) then
               call sub_elconf_generator(OrbNumb, CoreNumb, ValNumb, VirtNumb, FzrNumb, & 
                                          ElNumb, OccValNumb, 0, SubAddIdent, AllSubElConfNumb)
            else
               call sub_elconf_generator(OrbNumb, CoreNumb, ValNumb, VirtNumb, FzrNumb, & 
                                          ElNumb, OccValNumb, 1, SubAddIdent, AllSubElConfNumb)
            endif
         else
            AllSubElConfNumb = SubElConfNumb
         endif

      call sub_eta_coefficients (OrbNumb, AllSubElConfNumb, PhysVac, S_final_state, SubExOpNumb)

      call CPU_TIME(time2)
      write (*, *) ' ...DONE '
      write (*, *) ' '
      write (*, 10) time2 - time1
      write (*, *) ' '
      write (*, 33) SubExOpNumb
      write (*, *) ' '
      write (*, *) ' '

!
!     all the configurations generated have been saved in a temporary file "SUB_EL_CONFIG.temp",
!     which is deleted after execution of the program
!
!     the total number of configurations amounts to "SubElConfNumb"
!

!
!          SPIN-TENSOR OPERATOR GENERATION (SUBSIDIARY SPACE)
!

      FirstOrderIdent = .true.
      write (*, *) ' '
      write (*, *) ' ...Subsidiary space configuration function generation... '
      write (*, *) ' '
      call CPU_TIME(time1)

      call spin_tensor_up_level(OrbNumb, SubElConfNumb, S_final_state, PhysVac, FirstOrderIdent, S_curr_final, SubConfFuncNumb)

      call CPU_TIME(time2)
      write (*, *) ' ...DONE '
      write (*, *) ' '
      write (*, 10) time2 - time1
      write (*, *) ' '
      write (*, 34) SubConfFuncNumb
      write (*, *) ' '
      write (*, *) ' '

!
!     all the spin-tensor operators generated have been saved in temporary unformatted files "Q_SUB.temp",
!     "K_Q_SUB.temp", "CF_Q_SUB.temp", "COEF_Q_SUB.temp". For user conveniency, the formatted file "CONF_FUNC_STRUCTURE_SUB",
!     containing all the above mentioned data, was created.
!


      endif

   10 format (' EXECUTION TIME............................................................', F9.6)
   11 format (' NUMBER OF CORE ORBITALS.....................................................', I2)
   12 format (' NUMBER OF VALENT ORBITALS...................................................', I2)
   13 format (' NUMBER OF VIRTUAL ORBITALS..................................................', I2)
   14 format (' NUMBER OF FROZEN CORE ORBITALS..............................................', I2)
   15 format (' MULTIPLICITY OF INITIAL STATE...............................................', I2)
   16 format (' MULTIPLICITY OF FINAL STATES................................................', I2)
   17 format (' PHYSICAL VACUUM STRUCTURE............', 90I3)
   18 format (' TOTAL NUMBER OF ONE ELECTRON INTEGRALS......................................', I6)
   19 format (' TOTAL NUMBER OF TWO ELECTRON INTEGRALS......................................', I6)
   20 format (' MEAN ENERGY IN THE LOWEST ENERGY CLUSTER....................................', F12.6)
   21 format (' NUCLEAR REPULSION ENERGY....................................................', F12.6)
   22 format (' STATE WITH THE FOLLOWING NUMBER WILL BE TAKEN AS ZERO-ORDER INITIAL STATE.. ', I2)

   30 format (' TOTAL NUMBER OF VALENT ELECTRON CONFIGURATIONS .............................', I4)
   31 format (' TOTAL NUMBER OF VALENT CONFIGURATION FUNCTIONS .............................', I4)
   32 format (' TOTAL NUMBER OF SUBSIDIARY ELECTRON CONFIGURATIONS..........................', I4)
   33 format (' TOTAL NUMBER OF EXCITATION OPERATORS IN SUBSIDIARY SPACE....................', I4)
   34 format (' TOTAL NUMBER OF SUBSIDIARY CONFIGURATION FUNCTIONS..........................', I4)
      end 























