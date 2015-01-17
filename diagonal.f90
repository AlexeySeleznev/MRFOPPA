      subroutine diagonal (ki, q_cr_I, q_an_I, results)
      use integrals_amplitudes
      implicit none
      integer i, j
      integer ki
      integer q_cr_I(ki/2), q_an_I(ki/2)
      real*8 results

!       print *, 'q_cr_I = ', q_cr_I(:)
!       print *, 'q_an_I = ', q_an_I(:)
       if (ki .eq. 0) then
         results = PhysVacEner - MeanEner

       elseif (ki .ne. 0) then

         results = 0.00D0

         do i = 1, ki/2
            results = results + amplitudes((q_cr_I(i) + 1)/2) - amplitudes((q_an_I(i) + 1)/2)
         enddo

         if (ki .eq. 2) then
            results = results + PhysVacEner - MeanEner + oneparticle1(q_cr_I(1), q_cr_I(1)) - & 
            amplitudes((q_cr_I(1) + 1)/2) + amplitudes((q_an_I(1) + 1)/2) - oneparticle1(q_an_I(1), q_an_I(1)) - & 
            twoparticle1(q_cr_I(1), q_an_I(1), q_cr_I(1), q_an_I(1)) + twoparticle1(q_cr_I(1), q_an_I(1), q_an_I(1), q_cr_I(1))  

          elseif (ki .gt. 2) then 
            results = results + PhysVacEner - MeanEner

            do j = 1, ki/2
               results = results + amplitudes((q_an_I(j) + 1)/2) - oneparticle1(q_an_I(j), q_an_I(j))
            enddo

            do j = 1, ki/2
               results = results + oneparticle1(q_cr_I(j), q_cr_I(j)) - amplitudes((q_cr_I(j) + 1)/2)
            enddo

            do i = 1, ki/2
              do j= 1, ki/2
                 results = results + twoparticle1(q_cr_I(i), q_an_I(j), q_an_I(j), q_cr_I(i)) - & 
                                     twoparticle1(q_cr_I(i), q_an_I(j), q_cr_I(i), q_an_I(j))  
              enddo
            enddo

            do i = 1, ki/2 - 1
              do j = i + 1, ki/2
                 results = results + twoparticle1(q_cr_I(i), q_cr_I(j), q_cr_I(i), q_cr_I(j)) - & 
                                     twoparticle1(q_cr_I(i), q_cr_I(j), q_cr_I(j), q_cr_I(i))
              enddo
            enddo

            do i = 1, ki/2 - 1
              do j = i + 1, ki/2
                 results = results + twoparticle1(q_an_I(i), q_an_I(j), q_an_I(i), q_an_I(j)) - & 
                                     twoparticle1(q_an_I(i), q_an_I(j), q_an_I(j), q_an_I(i))
              enddo
            enddo

         endif
       endif
      end subroutine
