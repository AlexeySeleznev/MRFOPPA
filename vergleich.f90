       subroutine vergleich (a_dim, b_dim, a, b, DifNumb, Dif_Numer)
       implicit none
       integer i, j
       integer a_dim, b_dim
       integer a(a_dim), b(b_dim)
       integer DifNumb
       integer res, num1, temp
       integer Dif_Numer(2)

       DifNumb = 0
       Dif_Numer = 0
       num1 = 1

       do i = 1, a_dim
          res = 0
          do j = 1, b_dim
             if (a(i) .eq. b(j)) then
                res = res + 1
                exit
              endif
          enddo

          if (res .eq. 0) then 
            DifNumb = DifNumb + 1
              if (DifNumb .gt. 2) then
                 exit
               else
                 Dif_Numer(num1) = i
                 num1 = num1 + 1
              endif
          endif
       enddo
      
       if (DifNumb .eq. 1) then
          temp = a(Dif_Numer(1))
          a(Dif_Numer(1)) = a(a_dim)
          a(a_dim) = temp

       elseif (DifNumb .eq. 2) then
         if (Dif_Numer(2) .eq. a_dim) then
          temp = a(Dif_Numer(1))
          a(Dif_Numer(1)) = a(a_dim - 1)
          a(a_dim - 1) = temp
         elseif (Dif_Numer(1) .eq. a_dim) then
          temp = a(Dif_Numer(2))
          a(Dif_Numer(2)) = a(a_dim - 1)
          a(a_dim - 1) = temp
         else
          temp = a(Dif_Numer(2))
          a(Dif_Numer(2)) = a(a_dim)
          a(a_dim) = temp
          temp = a(Dif_Numer(1))
          a(Dif_Numer(1)) = a(a_dim - 1)
          a(a_dim - 1) = temp
         endif
       endif

       end subroutine
