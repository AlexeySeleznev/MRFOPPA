      module spin_tensor
      contains

      recursive subroutine spin_tensor_down (S, M, n_initial, n, a, path, q_val, coef_q_val, length)
      use coefficients
      implicit none
      integer n_initial, n
      integer length
      integer :: a(n_initial), path(n_initial)
      integer, target :: q_val(length, n_initial)
      real*8, target :: coef_q_val(length)
      integer M, S
      integer, pointer :: semi_q_val1(:, :), semi_q_val2(:, :)
      real*8, pointer :: semi_coef_q_val1(:), semi_coef_q_val2(:)
     
       if (n .eq. 1) then

         if (S .eq. 1 .and. M .eq. 1) then
            q_val(1, n) = 2 * a(n)
            coef_q_val(:) = -1.0D0 * coef_q_val(:)
         elseif (S .eq. 1 .and. M .eq. -1) then
            q_val(1, n) = 2 * a(n) - 1
            coef_q_val(:) = 1.0D0 * coef_q_val(:)
         endif
       
       elseif (n .ne. 1) then

         if (abs(M + 1) .le. abs(S - path(n)) .and. abs(M - 1) .gt. abs(S - path(n))) then
            q_val(:, n) = 2 * a(n) - 1
            coef_q_val(:) = coef_q_val(:) * coeff(S, M, path(n), -1)
            call spin_tensor_down (S - path(n), M + 1, n_initial, n - 1, a, path, q_val, coef_q_val, length)

          elseif (abs(M + 1) .gt. abs(S - path(n)) .and. abs(M - 1) .le. abs(S - path(n))) then
            q_val(:, n) = 2 * a(n)
            coef_q_val(:) = (-1) * coef_q_val(:) * coeff(S, M, path(n), 1)
            call spin_tensor_down (S - path(n), M - 1, n_initial, n - 1, a, path, q_val, coef_q_val, length)

          elseif (abs(M + 1) .le. abs(S - path(n)) .and. abs(M - 1) .le. abs(S - path(n))) then
            semi_q_val1 => q_val(1 : length/2, :)
            semi_q_val2 => q_val(length/2 + 1 : length, :)
            semi_coef_q_val1 => coef_q_val(1 : length/2)
            semi_coef_q_val2 => coef_q_val(length/2 + 1 : length)
            semi_q_val1(:, n) = 2 * a(n) - 1
            semi_coef_q_val1(:) = semi_coef_q_val1(:) * coeff(S, M, path(n), -1)
            call spin_tensor_down (S - path(n), M + 1, n_initial, n - 1, a, path, semi_q_val1, semi_coef_q_val1, length/2)
            semi_q_val2(:, n) = 2 * a(n)
            semi_coef_q_val2(:) = (-1) * semi_coef_q_val2(:) * coeff(S, M, path(n), 1)
            call spin_tensor_down (S - path(n), M - 1, n_initial, n - 1, a, path, semi_q_val2, semi_coef_q_val2, length/2)
         endif
    
       endif


      end subroutine

      recursive subroutine spin_tensor_up (S, M, n_initial, n, a, path, q_val, coef_q_val, length, up_to_down)
      use coefficients
      implicit none
      integer n_initial, n
      integer length
      integer up_to_down
      integer :: a(n_initial), path(n_initial)
      integer, target :: q_val(length, n_initial)
      real*8, target :: coef_q_val(length)
      integer M, S
      integer, pointer :: semi_q_val1(:, :), semi_q_val2(:, :)
      real*8, pointer :: semi_coef_q_val1(:), semi_coef_q_val2(:)
      integer i, j

       if (n .le. up_to_down) then

            call spin_tensor_down (S, M, n_initial, n, a, path, q_val, coef_q_val, length)

       elseif (n .gt. up_to_down) then

        if (n .eq. 1) then
          if (S .eq. 1 .and. M .eq. 1) then
             q_val(1, n) = 2 * a(n) - 1
             coef_q_val(:) = 1.0D0 * coef_q_val(:)
          elseif (S .eq. 1 .and. M .eq. -1) then
             q_val(1, n) = 2 * a(n)
             coef_q_val(:) = 1.0D0 * coef_q_val(:)
          endif
      
        elseif (n .ne. 1) then

         if (abs(M - 1) .le. abs(S - path(n)) .and. abs(M + 1) .gt. abs(S - path(n))) then
            q_val(:, n) = 2 * a(n) - 1
            coef_q_val(:) = coef_q_val(:) * coeff(S, M, path(n), 1)
            call spin_tensor_up (S - path(n), M - 1, n_initial, n - 1, a, path, q_val, coef_q_val, length, up_to_down)

          elseif (abs(M - 1) .gt. abs(S - path(n)) .and. abs(M + 1) .le. abs(S - path(n))) then
            q_val(:, n) = 2 * a(n)
            coef_q_val(:) = coef_q_val(:) * coeff(S, M, path(n), -1)
            call spin_tensor_up (S - path(n), M + 1, n_initial, n - 1, a, path, q_val, coef_q_val, length, up_to_down)

          elseif (abs(M - 1) .le. abs(S - path(n)) .and. abs(M + 1) .le. abs(S - path(n))) then
            semi_q_val1 => q_val(1 : length/2, :)
            semi_q_val2 => q_val(length/2 + 1 : length, :)
            semi_coef_q_val1 => coef_q_val(1 : length/2)
            semi_coef_q_val2 => coef_q_val(length/2 + 1 : length)
            semi_q_val1(:, n) = 2 * a(n) - 1
            semi_coef_q_val1(:) = semi_coef_q_val1(:) * coeff(S, M, path(n), 1)
            call spin_tensor_up (S - path(n), M - 1, n_initial, n - 1, a, path, semi_q_val1, semi_coef_q_val1, length/2, up_to_down)
            semi_q_val2(:, n) = 2 * a(n)
            semi_coef_q_val2(:) = semi_coef_q_val2(:) * coeff(S, M, path(n), -1)
            call spin_tensor_up (S - path(n), M + 1, n_initial, n - 1, a, path, semi_q_val2, semi_coef_q_val2, length/2, up_to_down)
         endif

        endif

       endif

      end subroutine

      end module

