      module coefficients
      contains
      real*8 function coeff(S, M, t, p)
      integer S, M
      integer t, p
      real*8 S_real, M_real

      S_real = real(S) * 0.50D0
      M_real = real(M) * 0.50D0

       if (t .eq. 1) then
         if (p .eq. 1) then
           if (S .gt. 0) then
             coeff = sqrt((S_real + M_real)/(2.00D0 * S_real))
            elseif (S .eq. 0) then
             coeff = sqrt(1.00D0/2.00D0)
           endif
          elseif (p .eq. -1) then
           if (S .gt. 0) then
             coeff = sqrt((S_real - M_real)/(2.00D0 * S_real))
            elseif (S .eq. 0) then
             coeff = sqrt(1.00D0/2.00D0)
           endif
         endif
        elseif (t .eq. -1) then
         if (p .eq. 1) then
           if (S .gt. 0) then
             coeff = (-1.00D0) * sqrt((S_real + 1.00D0 - M_real)/(2.00D0 * (S_real + 1.00D0)))
            elseif (S .eq. 0) then
             coeff = (-1.00D0) * sqrt(1.00D0/2.00D0)
           endif
          elseif (p .eq. -1) then
           if (S .gt. 0) then
            coeff = sqrt((S_real + 1.00D0 + M_real)/(2.00D0*(S_real + 1.00D0)))
           elseif (S .eq. 0) then
            coeff = sqrt(1.00D0/2.00D0)
           endif
         endif
       endif
      end function
      end module                                        
