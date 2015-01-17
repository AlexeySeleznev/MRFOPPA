      subroutine val_elconf_generator(OrbNumb, CoreNumb, ValNumb, ElNumb, OccValNumb, S, S_curr_final, ElConfNumb)
      use elconfig
      implicit none
      integer i, j                                 !dummy variables
      integer counter
      integer S                                    ! 2 * S^2 eigenvalue
      integer ElNumb
      integer OrbNumb, CoreNumb, ValNumb           ! number of orbitals
      integer OccValNumb                           ! number of occupied valent spin-orbitals
      integer MinSpin, MaxSpin
      integer ElConfNumb
      integer OneOccNumb, DublOccNumb
      logical Doub_ON_OFF, One_ON_OFF
      logical S_curr_final

      integer, target :: ElConf(OrbNumb)
      integer, pointer :: ElConfVal(:)

      integer, allocatable, dimension(:) :: OnePosition
      integer, allocatable, dimension(:, :) :: DublMagOrder, OneMagOrder

      ElConfNumb = 0
!
!            BOUNDS
!
      MinSpin = S
        if (OccValNumb > ValNumb) then
           MaxSpin = 2 * ValNumb - OccValNumb
         else
           MaxSpin = OccValNumb
        endif
!
!            GENERATOR
!
      if (MaxSpin .ge. MinSpin) then

         if (.not. S_curr_final) then
            open (20, file = 'VAL_EL_CONFIG_INIT.temp')
         elseif (S_curr_final) then
            open (20, file = 'VAL_EL_CONFIG.temp')
         endif
         ElConf = 0

         do i = 1, CoreNumb
          ElConf(i) = 2
         enddo

         ElConfVal => ElConf(CoreNumb + 1 : CoreNumb + ValNumb)
!
!
!      
        do OneOccNumb = MinSpin, MaxSpin, 2

           DublOccNumb = OccValNumb - OneOccNumb

           if (OneOccNumb .eq. 0) then
              if (DublOccNumb .gt. 0) then
                 allocate (DublMagOrder(DublOccNumb/2, 3))
                 do j = 1, DublOccNumb/2
                    DublMagOrder(j, 1) = j
                    DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                    DublMagOrder(j, 3) = j
                 enddo
              endif
              if (DublOccNumb .gt. 0) then
                 Doub_ON_OFF = .true.
                   do while (Doub_ON_OFF)
                     ElConfVal = 0
                      do i = 1, DublOccNumb/2
                         ElConfVal(DublMagOrder(i, 3)) = 2
                      enddo
                     ElConfNumb = ElConfNumb + 1
                     write (20, 50) (ElConf(j), j = 1, OrbNumb)
                     call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                   enddo
               deallocate(DublMagOrder)  
              endif
 

           elseif (OneOccNumb .gt. 0) then
              if (DublOccNumb .gt. 0) then
                 allocate (DublMagOrder(DublOccNumb/2, 3))
                 do j = 1, DublOccNumb/2
                    DublMagOrder(j, 1) = j
                    DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                    DublMagOrder(j, 3) = j
                 enddo
              endif

              if (DublOccNumb .gt. 0) then
                 Doub_ON_OFF = .true.
                   do while (Doub_ON_OFF)
                      ElConfVal = 0
                      do i = 1, DublOccNumb/2
                         ElConfVal(DublMagOrder(i, 3)) = 2
                      enddo
                      allocate (OnePosition(ValNumb - DublOccNumb/2))
                      OnePosition = 0
                      counter = 0
                      do i = 1, ValNumb
                         if (ElConfVal(i) .eq. 0) then
                           counter = counter + 1
                           OnePosition(counter) = i
                         endif
                      enddo
                      allocate(OneMagOrder(OneOccNumb, 3))
                      do j = 1, OneOccNumb
                         OneMagOrder(j, 1) = j
                         OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - DublOccNumb/2 - j + 1
                         OneMagOrder(j, 3) = j
                      enddo
                         One_ON_OFF = .true.
                            do while (One_ON_OFF)
                               do i = 1, ValNumb - DublOccNumb/2
                                  ElConfVal(OnePosition(i)) = 0
                               enddo
                               do i = 1, OneOccNumb
                                  ElConfVal(OnePosition(OneMagOrder(i, 3))) = 1
                               enddo
                               ElConfNumb = ElConfNumb + 1
                               write (20, 50) (ElConf(j), j = 1, OrbNumb)
                               call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                            enddo
                    deallocate (OneMagOrder)
                    deallocate (OnePosition)
                    call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                   enddo
                 deallocate (DublMagOrder)

              elseif (DublOccNumb .eq. 0) then
                      allocate(OneMagOrder(OneOccNumb, 3))
                      do j = 1, OneOccNumb
                         OneMagOrder(j, 1) = j
                         OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - j + 1
                         OneMagOrder(j, 3) = j
                      enddo
                         One_ON_OFF = .true.
                            do while (One_ON_OFF)
                               ElConfVal = 0
                               do i = 1, OneOccNumb
                                  ElConfVal(OneMagOrder(i, 3)) = 1
                               enddo
                               ElConfNumb = ElConfNumb + 1
                               write (20, 50) (ElConf(j), j = 1, OrbNumb)
                               call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                            enddo
                    deallocate (OneMagOrder)

              endif

           endif

        enddo
       close (20)
      endif
   50 format(80I4)
      end subroutine

      subroutine sub_elconf_generator(OrbNumb, CoreNumb, ValNumb, VirtNumb, FzrNumb, & 
                                      ElNumb, OccValNumb, S, SubAddIdent, ElConfNumb)
      use elconfig
      implicit none
      integer i, j, virt_i, core_i                 !dummy variables
      integer counter
      integer S                                    ! 2 * S^2 eigenvalue
      integer ElNumb
      integer OrbNumb, CoreNumb, ValNumb           ! number of orbitals
      integer VirtNumb, FzrNumb                    ! number of orbitals
      integer OccValNumb                           ! number of occupied valent spin-orbitals
      integer MinSpin, MaxSpin
      integer ElConfNumb
      integer OneOccNumb, DublOccNumb
      logical Doub_ON_OFF, One_ON_OFF
      logical SubAddIdent

      integer, target :: ElConf(OrbNumb)
      integer, pointer :: ElConfSubInner(:)

      integer, allocatable, dimension(:) :: OnePosition
      integer, allocatable, dimension(:, :) :: DublMagOrder, OneMagOrder

      ElConf = 0
      ElConfNumb = 0

      if (CoreNumb .ne. 0) then
        do i = 1, CoreNumb
           ElConf(i) = 2
        enddo
      endif

      if (.not. SubAddIdent) then
         open (20, file = 'SUB_EL_CONFIG.temp')
      elseif (SubAddIdent) then
         open (20, file = 'SUB_EL_CONFIG_ALL.temp')
      endif

      ElConfSubInner => ElConf(CoreNumb + 1 : CoreNumb + ValNumb)

      if (VirtNumb .gt. 0) then

          if (S .eq. 0) then
            MinSpin = 1
          elseif (S .ne. 0) then
            MinSpin = S - 1
          endif

          if (OccValNumb - 1 > ValNumb) then
            MaxSpin = 2 * ValNumb - OccValNumb + 1
          else
            MaxSpin = OccValNumb - 1
          endif

          if (MaxSpin .ge. MinSpin) then

            do virt_i = CoreNumb + ValNumb + 1, OrbNumb
              ElConf(virt_i) = 1

                do OneOccNumb = MinSpin, MaxSpin, 2
                  DublOccNumb = OccValNumb - 1 - OneOccNumb

                  if (OneOccNumb .eq. 0) then
                      if (DublOccNumb .gt. 0) then
                         allocate (DublMagOrder(DublOccNumb/2, 3))
                           do j = 1, DublOccNumb/2
                             DublMagOrder(j, 1) = j
                             DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                             DublMagOrder(j, 3) = j
                           enddo
                      endif
                      if (DublOccNumb .gt. 0) then
                        Doub_ON_OFF = .true.
                          do while (Doub_ON_OFF)
                            ElConfSubInner = 0
                              do i = 1, DublOccNumb/2
                                 ElConfSubInner(DublMagOrder(i, 3)) = 2
                              enddo
                            ElConfNumb = ElConfNumb + 1
                            write (20, 50) (ElConf(j), j = 1, OrbNumb)
                            call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                          enddo
                        deallocate(DublMagOrder)
                      endif


                  elseif (OneOccNumb .gt. 0) then
                      if (DublOccNumb .gt. 0) then
                         allocate (DublMagOrder(DublOccNumb/2, 3))
                           do j = 1, DublOccNumb/2
                             DublMagOrder(j, 1) = j
                             DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                             DublMagOrder(j, 3) = j
                           enddo
                      endif
                      if (DublOccNumb .gt. 0) then
                          Doub_ON_OFF = .true.
                            do while (Doub_ON_OFF)
                               ElConfSubInner = 0
                                  do i = 1, DublOccNumb/2
                                     ElConfSubInner(DublMagOrder(i, 3)) = 2
                                  enddo
                               allocate (OnePosition(ValNumb - DublOccNumb/2))
                               OnePosition = 0
                               counter = 0
                                  do i = 1, ValNumb
                                     if (ElConfSubInner(i) .eq. 0) then
                                         counter = counter + 1
                                         OnePosition(counter) = i
                                     endif
                                  enddo
                               allocate(OneMagOrder(OneOccNumb, 3))
                                  do j = 1, OneOccNumb
                                     OneMagOrder(j, 1) = j
                                     OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - & 
                                                                DublOccNumb/2 - j + 1
                                     OneMagOrder(j, 3) = j
                                  enddo
                               One_ON_OFF = .true.
                                  do while (One_ON_OFF)
                                     do i = 1, ValNumb - DublOccNumb/2
                                        ElConfSubInner(OnePosition(i)) = 0
                                     enddo
                                     do i = 1, OneOccNumb
                                        ElConfSubInner(OnePosition(OneMagOrder(i, 3))) = 1
                                     enddo
                                     ElConfNumb = ElConfNumb + 1
                                     write (20, 50) (ElConf(j), j = 1, OrbNumb)
                                     call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                                  enddo
                               deallocate (OneMagOrder)
                               deallocate (OnePosition)
                               call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                            enddo
                         deallocate (DublMagOrder)

                      elseif (DublOccNumb .eq. 0) then
                         allocate(OneMagOrder(OneOccNumb, 3))
                           do j = 1, OneOccNumb
                              OneMagOrder(j, 1) = j
                              OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - j + 1
                              OneMagOrder(j, 3) = j
                           enddo
                              One_ON_OFF = .true.
                                 do while (One_ON_OFF)
                                    ElConfSubInner = 0
                                      do i = 1, OneOccNumb
                                        ElConfSubInner(OneMagOrder(i, 3)) = 1
                                      enddo
                                    ElConfNumb = ElConfNumb + 1
                                    write (20, 50) (ElConf(j), j = 1, OrbNumb)
                                    call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                                 enddo
                              deallocate (OneMagOrder)

                      endif
                  endif
                enddo
              ElConf(virt_i) = 0
            enddo
          endif
      endif

      if (CoreNumb - FzrNumb .gt. 0) then

          if (S .eq. 0) then
            MinSpin = 1
          elseif (S .ne. 0) then
            MinSpin = S + 1
          endif

          if (OccValNumb + 1 > ValNumb) then
            MaxSpin = 2 * ValNumb - OccValNumb - 1
          else
            MaxSpin = OccValNumb + 1
          endif

          if (MaxSpin .ge. MinSpin) then

            do core_i = CoreNumb, 1 + FzrNumb, -1
              ElConf(core_i) = 1

                do OneOccNumb = MinSpin, MaxSpin, 2
                  DublOccNumb = OccValNumb + 1 - OneOccNumb

                  if (OneOccNumb .eq. 0) then
                      if (DublOccNumb .gt. 0) then
                         allocate (DublMagOrder(DublOccNumb/2, 3))
                           do j = 1, DublOccNumb/2
                             DublMagOrder(j, 1) = j
                             DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                             DublMagOrder(j, 3) = j
                           enddo
                      endif
                      if (DublOccNumb .gt. 0) then
                        Doub_ON_OFF = .true.
                          do while (Doub_ON_OFF)
                            ElConfSubInner = 0
                              do i = 1, DublOccNumb/2
                                 ElConfSubInner(DublMagOrder(i, 3)) = 2
                              enddo
                            ElConfNumb = ElConfNumb + 1
                            write (20, 50) (ElConf(j), j = 1, OrbNumb)
                            call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                          enddo
                        deallocate(DublMagOrder)
                      endif

                  elseif (OneOccNumb .gt. 0) then
                      if (DublOccNumb .gt. 0) then
                         allocate (DublMagOrder(DublOccNumb/2, 3))
                           do j = 1, DublOccNumb/2
                             DublMagOrder(j, 1) = j
                             DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                             DublMagOrder(j, 3) = j
                           enddo
                      endif
                      if (DublOccNumb .gt. 0) then
                          Doub_ON_OFF = .true.
                            do while (Doub_ON_OFF)
                               ElConfSubInner = 0
                                  do i = 1, DublOccNumb/2
                                     ElConfSubInner(DublMagOrder(i, 3)) = 2
                                  enddo
                               allocate (OnePosition(ValNumb - DublOccNumb/2))
                               OnePosition = 0
                               counter = 0
                                  do i = 1, ValNumb
                                     if (ElConfSubInner(i) .eq. 0) then
                                         counter = counter + 1
                                         OnePosition(counter) = i
                                     endif
                                  enddo
                               allocate(OneMagOrder(OneOccNumb, 3))
                                  do j = 1, OneOccNumb
                                     OneMagOrder(j, 1) = j
                                     OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - & 
                                                                DublOccNumb/2 - j + 1
                                     OneMagOrder(j, 3) = j
                                  enddo
                               One_ON_OFF = .true.
                                  do while (One_ON_OFF)
                                     do i = 1, ValNumb - DublOccNumb/2
                                        ElConfSubInner(OnePosition(i)) = 0
                                     enddo
                                     do i = 1, OneOccNumb
                                        ElConfSubInner(OnePosition(OneMagOrder(i, 3))) = 1
                                     enddo
                                     ElConfNumb = ElConfNumb + 1
                                     write (20, 50) (ElConf(j), j = 1, OrbNumb)
                                     call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                                  enddo
                               deallocate (OneMagOrder)
                               deallocate (OnePosition)
                               call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                            enddo
                         deallocate (DublMagOrder)

                     elseif (DublOccNumb .eq. 0) then
                         allocate(OneMagOrder(OneOccNumb, 3))
                           do j = 1, OneOccNumb
                              OneMagOrder(j, 1) = j
                              OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - j + 1
                              OneMagOrder(j, 3) = j
                           enddo
                              One_ON_OFF = .true.
                                 do while (One_ON_OFF)
                                    ElConfSubInner = 0
                                      do i = 1, OneOccNumb
                                        ElConfSubInner(OneMagOrder(i, 3)) = 1
                                      enddo
                                    ElConfNumb = ElConfNumb + 1
                                    write (20, 50) (ElConf(j), j = 1, OrbNumb)
                                    call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                                 enddo
                              deallocate (OneMagOrder)
                      endif
                  endif
                enddo
              ElConf(core_i) = 2
            enddo
          endif
      endif

      if (VirtNumb .gt. 0 .and. CoreNumb - FzrNumb .gt. 0) then

            MinSpin = S
          if (OccValNumb > ValNumb) then
            MaxSpin = 2 * ValNumb - OccValNumb
          else
            MaxSpin = OccValNumb
          endif

          if (MaxSpin .ge. MinSpin) then

            do core_i = CoreNumb, 1 + FzrNumb, -1
              ElConf(core_i) = 1

              do virt_i = CoreNumb + ValNumb + 1, OrbNumb
                ElConf(virt_i) = 1

                  do OneOccNumb = MinSpin, MaxSpin, 2

                    DublOccNumb = OccValNumb - OneOccNumb

                    if (OneOccNumb .eq. 0) then
                        if (DublOccNumb .gt. 0) then
                           allocate (DublMagOrder(DublOccNumb/2, 3))
                             do j = 1, DublOccNumb/2
                               DublMagOrder(j, 1) = j
                               DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                               DublMagOrder(j, 3) = j
                             enddo
                        endif
                        if (DublOccNumb .gt. 0) then
                          Doub_ON_OFF = .true.
                            do while (Doub_ON_OFF)
                              ElConfSubInner = 0
                                do i = 1, DublOccNumb/2
                                   ElConfSubInner(DublMagOrder(i, 3)) = 2
                                enddo
                              ElConfNumb = ElConfNumb + 1
                              write (20, 50) (ElConf(j), j = 1, OrbNumb)
                              call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                            enddo
                          deallocate(DublMagOrder)
                        endif

                    elseif (OneOccNumb .gt. 0) then
                        if (DublOccNumb .gt. 0) then
                           allocate (DublMagOrder(DublOccNumb/2, 3))
                             do j = 1, DublOccNumb/2
                               DublMagOrder(j, 1) = j
                               DublMagOrder(DublOccNumb/2 - j + 1, 2) = ValNumb - j + 1
                               DublMagOrder(j, 3) = j
                             enddo
                        endif
                        if (DublOccNumb .gt. 0) then
                            Doub_ON_OFF = .true.
                              do while (Doub_ON_OFF)
                                 ElConfSubInner = 0
                                    do i = 1, DublOccNumb/2
                                       ElConfSubInner(DublMagOrder(i, 3)) = 2
                                    enddo
                                 allocate (OnePosition(ValNumb - DublOccNumb/2))
                                 OnePosition = 0
                                 counter = 0
                                    do i = 1, ValNumb
                                       if (ElConfSubInner(i) .eq. 0) then
                                           counter = counter + 1
                                           OnePosition(counter) = i
                                       endif
                                    enddo
                                 allocate(OneMagOrder(OneOccNumb, 3))
                                    do j = 1, OneOccNumb
                                       OneMagOrder(j, 1) = j
                                       OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - & 
                                                                DublOccNumb/2 - j + 1
                                       OneMagOrder(j, 3) = j
                                    enddo
                                 One_ON_OFF = .true.
                                    do while (One_ON_OFF)
                                       do i = 1, ValNumb - DublOccNumb/2
                                          ElConfSubInner(OnePosition(i)) = 0
                                       enddo
                                       do i = 1, OneOccNumb
                                          ElConfSubInner(OnePosition(OneMagOrder(i, 3))) = 1
                                       enddo
                                       ElConfNumb = ElConfNumb + 1
                                       write (20, 50) (ElConf(j), j = 1, OrbNumb)
                                       call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                                    enddo
                                 deallocate (OneMagOrder)
                                 deallocate (OnePosition)
                                 call switcher_ElConf (DublOccNumb/2, DublMagOrder, Doub_ON_OFF)
                              enddo
                           deallocate (DublMagOrder)

                       elseif (DublOccNumb .eq. 0) then
                           allocate(OneMagOrder(OneOccNumb, 3))
                             do j = 1, OneOccNumb
                                OneMagOrder(j, 1) = j
                                OneMagOrder(OneOccNumb - j + 1, 2) = ValNumb - j + 1
                                OneMagOrder(j, 3) = j
                             enddo
                                One_ON_OFF = .true.
                                   do while (One_ON_OFF)
                                      ElConfSubInner = 0
                                        do i = 1, OneOccNumb
                                          ElConfSubInner(OneMagOrder(i, 3)) = 1
                                        enddo
                                      ElConfNumb = ElConfNumb + 1
                                      write (20, 50) (ElConf(j), j = 1, OrbNumb)
                                      call switcher_ElConf (OneOccNumb, OneMagOrder, One_ON_OFF)
                                   enddo
                                deallocate (OneMagOrder)
                        endif
                    endif
                  enddo
                ElConf(virt_i) = 0
              enddo
              ElConf(core_i) = 2
            enddo
          endif
      endif
      nullify(ElConfSubInner) 
      close (20)
   50 format(80I4)
      end subroutine


