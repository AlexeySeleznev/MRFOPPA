      subroutine V_element (ki, kj, q_I, q_J, MaxOpNumb, results)
      implicit none
      integer i, j
      integer ki, kj, MaxOpNumb 
      integer q_I(MaxOpNumb), q_J(MaxOpNumb)
      integer numer1(2), numer2(2), numer3(2), numer4(2)
      integer q_cr_I(ki/2), q_an_I(ki/2), q_cr_J(kj/2), q_an_J(kj/2)
      integer teta1, teta2, teta3, teta4, teta
      integer multiplier
      real*8 phase
      real*8 results

      if (abs(ki - kj) .gt. 4) then
        results = 0.00D0

      else
        multiplier = (-1) ** ((ki/2 * (ki/2 - 1) + kj/2 * (kj/2 - 1))/2)
        
          do i = 1, ki/2
             q_an_I(i) = q_I(i)
          enddo
          do i = 1, ki/2
             q_cr_I(i) = q_I(ki/2 + i)
          enddo
          do j = 1, kj/2
             q_an_J(j) = q_J(j)
          enddo
          do j = 1, kj/2
             q_cr_J(j) = q_J(kj/2 + j)
          enddo

          call vergleich (ki/2, kj/2, q_cr_I, q_cr_J, teta1, numer1)
          call vergleich (ki/2, kj/2, q_an_I, q_an_J, teta3, numer3)
          call vergleich (kj/2, ki/2, q_cr_J, q_cr_I, teta2, numer2)
          call vergleich (kj/2, ki/2, q_an_J, q_an_I, teta4, numer4)

          teta = teta1 + teta2 + teta3 + teta4

            if (teta .gt. 4) then                                 
               results = 0.00D0

              else

                if (ki .eq. kj .and. teta .eq. 0) then
                   call diagonal (ki, q_cr_I, q_an_I, results)

                elseif (ki .eq. kj .and. teta3 .eq. 1 .and. teta4 .eq. 1 .and. teta .eq. 2) then
                   multiplier = multiplier * (-1)**(numer3(1) + numer4(1))
                   call case2 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj .and. teta2 .eq. 1 .and. teta1 .eq. 1 .and. teta .eq. 2) then
                   multiplier = multiplier * (-1)**(numer1(1) + numer2(1))
                   call case3 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj .and. teta .eq. 4 .and. teta1 .eq. 1 .and. teta2 .eq. 1 .and. & 
                        teta3 .eq. 1 .and. teta4 .eq. 1) then
                   multiplier = multiplier * (-1)**(numer1(1) + numer2(1) + numer3(1) + numer4(1))
                   call case4 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj .and. teta .eq. 4 .and. teta4 .eq. 2 .and. teta3 .eq. 2) then
                   multiplier = multiplier * (-1)**(numer3(1) + numer3(2) + numer4(1) + numer4(2))
                   call case5 (ki, q_an_I, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj .and. teta .eq. 4 .and. teta2 .eq. 2 .and. teta1 .eq. 2) then
                   multiplier = multiplier * (-1)**(numer2(1) + numer2(2) + numer1(1) + numer1(2))
                   call case6 (ki, q_cr_I, q_cr_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj + 2 .and. teta .eq. 2 .and. teta3 .eq. 1 .and. teta1 .eq. 1) then
                   multiplier = multiplier * (-1)**(numer3(1) + numer1(1))
                   call case7 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj + 2 .and. teta .eq. 4 .and. teta3 .eq. 2 .and. teta1 .eq. 1 .and. teta4 .eq. 1) then
                   multiplier = multiplier * (-1)**(numer3(1) + numer3(2) + numer1(1) + numer4(1))
                   call case8 (ki, q_cr_I, q_an_I, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj + 2 .and. teta1 .eq. 2 .and. teta3 .eq. 1 .and. teta2 .eq. 1 .and. teta .eq. 4) then
                   multiplier = multiplier * (-1)**(numer1(1) + numer1(2) + numer3(1) + numer2(1))
                   call case9 (ki, q_cr_J, q_an_I, q_cr_I, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki + 2 .eq. kj .and. teta .eq. 2 .and. teta4 .eq. 1 .and. teta2 .eq. 1) then
                   multiplier = multiplier * (-1)**(numer4(1) + numer2(1))
                   call case10 (ki, q_cr_I, q_an_I, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki + 2 .eq. kj .and. teta .eq. 4 .and. teta4 .eq. 2 .and. teta2 .eq. 1 .and. teta3 .eq. 1) then
                   multiplier = multiplier * (-1)**(numer4(1) + numer2(1) + numer4(2) + numer3(1))
                   call case11 (ki, q_an_I, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki + 2 == kj .and. teta == 4 .and. teta2 == 2 .and. teta4 == 1 .and. teta1 == 1) then
                   multiplier = multiplier * (-1)**(numer2(1) + numer2(2) + numer4(1) + numer1(1))
                   call case12 (ki, q_cr_I, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki .eq. kj + 4 .and. teta1 .eq. 2 .and. teta3 .eq. 2 .and. teta .eq. 4) then
                   multiplier = multiplier * (-1)**(numer1(1) + numer1(2) + numer3(1) + numer3(2))
                   call case13 (ki, q_cr_I, q_an_I, results)
                   phase = real(multiplier)
                   results = results * phase

                elseif (ki + 4 .eq. kj .and. teta .eq. 4 .and. teta2 .eq. 2 .and. teta4 .eq. 2) then
                   multiplier = multiplier * (-1)**(numer2(1) + numer2(2) + numer4(1) + numer4(2))
                   call case14(ki, q_cr_J, q_an_J, results)
                   phase = real(multiplier)
                   results = results * phase
                                  
                else 
                    results = 0.00D0

                endif
            endif 
      endif
      end subroutine  
