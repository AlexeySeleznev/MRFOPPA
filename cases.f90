       subroutine case2 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       integer q_cr_J(ki/2), q_an_J(ki/2)
       real*8 results

         if (ki .eq. 2) then
            results = - oneparticle1(q_an_J(ki/2), q_an_I(ki/2))
              do i = 1, ki/2
                 results = results + twoparticle1(q_cr_I(i), q_an_J(ki/2), q_an_I(ki/2), q_cr_I(i)) - & 
                                     twoparticle1(q_cr_I(i), q_an_J(ki/2), q_cr_I(i), q_an_I(ki/2))
              enddo
          elseif (ki .gt. 2) then
            results = - oneparticle1(q_an_J(ki/2), q_an_I(ki/2))
              do i = 1, ki/2
                 results = results + twoparticle1(q_cr_I(i), q_an_J(ki/2), q_an_I(ki/2), q_cr_I(i)) - & 
                                     twoparticle1(q_cr_I(i), q_an_J(ki/2), q_cr_I(i), q_an_I(ki/2))
              enddo
              do j = 1, ki/2 - 1
                 results = results + twoparticle1(q_an_J(ki/2), q_an_I(j), q_an_I(ki/2), q_an_I(j)) - & 
                                     twoparticle1(q_an_J(ki/2), q_an_I(j), q_an_I(j), q_an_I(ki/2))
              enddo
         endif
       end subroutine

       subroutine case3 (ki, q_cr_I, q_an_I, q_cr_J, q_an_j, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       integer q_cr_J(ki/2), q_an_J(ki/2)
       real*8 results

         if (ki .eq. 2) then
            results = oneparticle1(q_cr_I(ki/2), q_cr_J(ki/2))
              do i = 1, ki/2
                 results = results + twoparticle1(q_cr_I(ki/2), q_an_I(i), q_an_I(i), q_cr_J(ki/2)) - & 
                                     twoparticle1(q_cr_I(ki/2), q_an_I(i), q_cr_J(ki/2), q_an_I(i))
              enddo
          elseif (ki .gt. 2) then
            results = oneparticle1(q_cr_I(ki/2), q_cr_J(ki/2))
              do i = 1, ki/2
                 results = results + twoparticle1(q_cr_I(ki/2), q_an_I(i), q_an_I(i), q_cr_J(ki/2)) - & 
                                     twoparticle1(q_cr_I(ki/2), q_an_I(i), q_cr_J(ki/2), q_an_I(i))
              enddo
              do j = 1, ki/2 - 1
                 results = results + twoparticle1(q_cr_I(j), q_cr_I(ki/2), q_cr_I(j), q_cr_J(ki/2)) - & 
                                     twoparticle1(q_cr_I(j), q_cr_I(ki/2), q_cr_J(ki/2), q_cr_I(j))
              enddo
         endif
       end subroutine

       subroutine case4 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       integer q_cr_J(ki/2), q_an_J(ki/2)
       real*8 results

            results = twoparticle1(q_cr_I(ki/2), q_an_J(ki/2), q_an_I(ki/2), q_cr_J(ki/2)) - & 
                      twoparticle1(q_cr_I(ki/2), q_an_J(ki/2), q_cr_J(ki/2), q_an_I(ki/2))
       end subroutine

       subroutine case5 (ki, q_an_I, q_an_J, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_an_I(ki/2)
       integer q_an_J(ki/2)
       real*8 results

            results = twoparticle1(q_an_J(ki/2), q_an_J(ki/2 - 1), q_an_I(ki/2), q_an_I(ki/2 - 1)) - & 
                      twoparticle1(q_an_J(ki/2), q_an_J(ki/2 - 1), q_an_I(ki/2 - 1), q_an_I(ki/2))
       end subroutine

       subroutine case6 (ki, q_cr_I, q_cr_J, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2)
       integer q_cr_J(ki/2)
       real*8 results

            results = twoparticle1(q_cr_I(ki/2), q_cr_I(ki/2 - 1), q_cr_J(ki/2), q_cr_J(ki/2 - 1)) - & 
                      twoparticle1(q_cr_I(ki/2), q_cr_I(ki/2 - 1), q_cr_J(ki/2 - 1), q_cr_J(ki/2))
       end subroutine

       subroutine case7 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
! here ki = ki, kj = ki - 2
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       integer q_cr_J(ki/2 - 1), q_an_J(ki/2 - 1)
       real*8 results

         if (ki .eq. 2) then
            results = oneparticle1(q_cr_I(ki/2), q_an_I(ki/2))

          elseif (ki .gt. 2) then
            results = oneparticle1(q_cr_I(ki/2), q_an_I(ki/2))

              do i = 1, ki/2 - 1
                 results = results - twoparticle1(q_cr_I(ki/2), q_an_J(i), q_an_I(ki/2), q_an_J(i)) + & 
                                     twoparticle1(q_cr_I(ki/2), q_an_J(i), q_an_J(i), q_an_I(ki/2))
              enddo
              do j = 1, ki/2 - 1
                 results = results + twoparticle1(q_cr_I(ki/2), q_cr_J(j), q_an_I(ki/2), q_cr_J(j)) - & 
                                     twoparticle1(q_cr_I(ki/2), q_cr_J(j), q_cr_J(j), q_an_I(ki/2))
              enddo

         endif
       end subroutine

       subroutine case8 (ki, q_cr_I, q_an_I, q_an_J, results)
! here ki = ki, kj = ki - 2
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       integer q_an_J(ki/2 - 1)
       real*8 results

                 results = twoparticle1(q_cr_I(ki/2), q_an_J(ki/2 - 1), q_an_I(ki/2 - 1), q_an_I(ki/2)) - & 
                           twoparticle1(q_cr_I(ki/2), q_an_J(ki/2 - 1), q_an_I(ki/2), q_an_I(ki/2 - 1))
       end subroutine

       subroutine case9 (ki, q_cr_J, q_an_I, q_cr_I, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       integer q_cr_J(ki/2 - 1)
       real*8 results

                 results = twoparticle1(q_cr_I(ki/2), q_cr_I(ki/2 - 1), q_an_I(ki/2), q_cr_J(ki/2 - 1)) - & 
                           twoparticle1(q_cr_I(ki/2), q_cr_I(ki/2 - 1), q_cr_J(ki/2 - 1), q_an_I(ki/2))
       end subroutine

       subroutine case10 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
! ki = ki, kj = ki + 2
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       integer q_cr_J(ki/2 + 1), q_an_J(ki/2 + 1)
       real*8 results

         if (ki .eq. 0) then
            results = oneparticle1(q_an_J(ki/2 + 1), q_cr_J(ki/2 + 1))

          elseif (ki .gt. 0) then
            results = oneparticle1(q_an_J(ki/2 + 1), q_cr_J(ki/2 + 1))
              do i = 1, ki/2
                 results = results - twoparticle1(q_an_J(ki/2 + 1), q_an_I(i), q_cr_J(ki/2 + 1), q_an_I(i)) + & 
                                     twoparticle1(q_an_J(ki/2 + 1), q_an_I(i), q_an_I(i), q_cr_J(ki/2 + 1))
              enddo
              do j = 1, ki/2 
                 results = results + twoparticle1(q_an_J(ki/2 + 1), q_cr_I(j), q_cr_J(ki/2 + 1), q_cr_I(j)) - & 
                                     twoparticle1(q_an_J(ki/2 + 1), q_cr_I(j), q_cr_I(j), q_cr_J(ki/2 + 1))
              enddo
         endif
       end subroutine

       subroutine case11 (ki, q_an_I, q_cr_J, q_an_J, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_an_I(ki/2)
       integer q_cr_J(ki/2 + 1), q_an_J(ki/2 + 1)
       real*8 results

                 results = twoparticle1(q_an_J(ki/2 + 1), q_an_J(ki/2), q_an_I(ki/2), q_cr_J(ki/2 + 1)) - & 
                           twoparticle1(q_an_J(ki/2 + 1), q_an_J(ki/2), q_cr_J(ki/2 + 1), q_an_I(ki/2))
       end subroutine

       subroutine case12 (ki, q_cr_I, q_cr_J, q_an_J, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2)
       integer q_cr_J(ki/2 + 1), q_an_J(ki/2 + 1)
       real*8 results

                 results = twoparticle1(q_an_J(ki/2 + 1), q_cr_I(ki/2), q_cr_J(ki/2 + 1), q_cr_J(ki/2)) - & 
                           twoparticle1(q_an_J(ki/2 + 1), q_cr_I(ki/2), q_cr_J(ki/2), q_cr_J(ki/2 + 1))
       end subroutine

       subroutine case13 (ki, q_cr_I, q_an_I, results)
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_I(ki/2), q_an_I(ki/2)
       real*8 results

                 results = twoparticle1(q_cr_I(ki/2), q_cr_I(ki/2 - 1), q_an_I(ki/2), q_an_I(ki/2 - 1)) - & 
                           twoparticle1(q_cr_I(ki/2), q_cr_I(ki/2 - 1), q_an_I(ki/2 - 1), q_an_I(ki/2))
       end subroutine

       subroutine case14 (ki, q_cr_J, q_an_J, results)
! ki = ki, kj = ki + 4
       use integrals_amplitudes
       implicit none
       integer i, j
       integer ki
       integer q_cr_J(ki/2 + 2), q_an_J(ki/2 + 2)
       real*8 results

                 results = twoparticle1(q_an_J(ki/2 + 2), q_an_J(ki/2 + 1), q_cr_J(ki/2 + 2), q_cr_J(ki/2 + 1)) - & 
                           twoparticle1(q_an_J(ki/2 + 2), q_an_J(ki/2 + 1), q_cr_J(ki/2 + 1), q_cr_J(ki/2 + 2))
       end subroutine

