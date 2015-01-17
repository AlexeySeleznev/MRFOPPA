      module elconfig
      contains

      recursive subroutine switcher_ElConf (L_dim, L, LOOP_ON_OFF)
      implicit none
      integer i, j
      integer L_dim
      logical LOOP_ON_OFF
      integer, target :: L(L_dim, 3)
      integer, pointer :: L_next(:, :)
         if (L_dim .eq. 1) then
            if (L(1, 3) .eq. L(1, 2)) then
              LOOP_ON_OFF = .false.
            elseif (L(1, 3) .lt. L(1, 2)) then
              L(1, 3) = L(1, 3) + 1
            endif

         elseif (L_dim .gt. 1) then
            if (L(L_dim, 3) .eq. L(L_dim, 2)) then
!             if (L(L_dim, 1) + 1 .le. L(L_dim, 2)) then
!                 L(L_dim, 1) = L(L_dim, 1) + 1
!             elseif (L(L_dim, 1) .eq. L(L_dim, 2) .and. L(L_dim - 1, 1) + 1
!             endif
!              L(L_dim, 3) = L(L_dim, 1)
              L_next => L(1 : L_dim - 1, :)
              call switcher_ElConf (L_dim - 1, L_next, LOOP_ON_OFF)
              nullify (L_next)
             if (L(L_dim, 1) + 1 .le. L(L_dim, 2)) then
                 L(L_dim, 1) = L(L_dim, 1) + 1
             elseif (L(L_dim, 1) .eq. L(L_dim, 2) .and. L(L_dim - 1, 1) + 1 .le. L(L_dim, 2)) then
                 L(L_dim, 1) = L(L_dim - 1, 1) + 1
             endif
              L(L_dim, 3) = L(L_dim, 1)
            elseif (L(L_dim, 3) .lt. L(L_dim, 2)) then
              L(L_dim, 3) = L(L_dim, 3) + 1
            endif
         endif
      end subroutine

      recursive subroutine operator_constructor(OrbNumb, Compliance, operators_cr, & 
                                              operators_an, length, length_cr, length_an, height, CurrNumb)
      implicit none
      integer OrbNumb
      integer CurrNumb
      integer Compliance(OrbNumb)
      integer length, length_cr, length_an, height
      integer, target :: operators_cr(height, length)
      integer, target :: operators_an(height, length)
      integer, pointer :: operators_cr_1(:, :), operators_cr_2(:, :)
      integer, pointer :: operators_an_1(:, :), operators_an_2(:, :)

      if (CurrNumb .gt. 0) then
         if (Compliance(CurrNumb) .eq. 1) then
            operators_an(:, length_an) = 2 * CurrNumb
            operators_an(:, length_an - 1) = 2 * CurrNumb - 1
            call operator_constructor(OrbNumb, Compliance, operators_cr, operators_an, length, & 
                                      length_cr, length_an - 2, height, CurrNumb - 1)
         elseif (Compliance(CurrNumb) .eq. 0) then
            call operator_constructor(OrbNumb, Compliance, operators_cr, operators_an, length, & 
                                      length_cr, length_an, height, CurrNumb - 1)
         elseif (Compliance(CurrNumb) .eq. 2) then
            operators_cr(:, length_cr) = 2 * CurrNumb
            operators_cr(:, length_cr - 1) = 2 * CurrNumb - 1
            call operator_constructor(OrbNumb, Compliance, operators_cr, operators_an, length, & 
                                      length_cr - 2, length_an, height, CurrNumb - 1)
         elseif (Compliance(CurrNumb) .eq. 30) then
            operators_an_1 => operators_an(1 : height/2, :)
            operators_an_2 => operators_an(height/2 + 1 : height, :)
            operators_cr_1 => operators_cr(1 : height/2, :)
            operators_cr_2 => operators_cr(height/2 + 1 : height, :)
            operators_an_1(:, length_an) = 2 * CurrNumb
            operators_an_2(:, length_an) = 2 * CurrNumb - 1
            call operator_constructor(OrbNumb, Compliance, operators_cr_1, operators_an_1, length, & 
                                      length_cr, length_an - 1, height/2, CurrNumb - 1)      
            call operator_constructor(OrbNumb, Compliance, operators_cr_2, operators_an_2, length, & 
                                      length_cr, length_an - 1, height/2, CurrNumb - 1)
            nullify(operators_an_1, operators_an_2, operators_cr_1, operators_cr_2)
         elseif (Compliance(CurrNumb) .eq. 4) then
            operators_an_1 => operators_an(1 : height/2, :)
            operators_an_2 => operators_an(height/2 + 1 : height, :)
            operators_cr_1 => operators_cr(1 : height/2, :)
            operators_cr_2 => operators_cr(height/2 + 1 : height, :)
            operators_cr_1(:, length_cr) = 2 * CurrNumb
            operators_cr_2(:, length_cr) = 2 * CurrNumb - 1
            call operator_constructor(OrbNumb, Compliance, operators_cr_1, operators_an_1, length, & 
                                      length_cr - 1, length_an, height/2, CurrNumb - 1)      
            call operator_constructor(OrbNumb, Compliance, operators_cr_2, operators_an_2, length, & 
                                      length_cr - 1, length_an, height/2, CurrNumb - 1)
            nullify(operators_an_1, operators_an_2, operators_cr_1, operators_cr_2)
         elseif (Compliance(CurrNumb) .eq. 3421) then
            operators_cr(:, length_cr) = 2 * CurrNumb
            call operator_constructor(OrbNumb, Compliance, operators_cr, operators_an, length, & 
                                      length_cr - 1, length_an, height, CurrNumb - 1)
         elseif (Compliance(CurrNumb) .eq. 3411) then
            operators_an(:, length_an) = 2 * CurrNumb - 1
            call operator_constructor(OrbNumb, Compliance, operators_cr, operators_an, length, & 
                                      length_cr, length_an - 1, height, CurrNumb - 1)
         elseif (Compliance(CurrNumb) .eq. 3422) then
            operators_cr(:, length_cr) = 2 * CurrNumb - 1
            call operator_constructor(OrbNumb, Compliance, operators_cr, operators_an, length, & 
                                       length_cr - 1, length_an, height, CurrNumb - 1)
         elseif (Compliance(CurrNumb) .eq. 3412) then
            operators_an(:, length_an) = 2 * CurrNumb
            call operator_constructor(OrbNumb, Compliance, operators_cr, operators_an, length, & 
                                      length_cr, length_an - 1, height, CurrNumb - 1)
         elseif (Compliance(CurrNumb) .eq. 341) then
            operators_an_1 => operators_an(1 : height/2, :)
            operators_an_2 => operators_an(height/2 + 1 : height, :)
            operators_cr_1 => operators_cr(1 : height/2, :)
            operators_cr_2 => operators_cr(height/2 + 1 : height, :)
            operators_an_1(:, length_an) = 0
            operators_cr_1(:, length_cr) = 0
            operators_an_2(:, length_an) = 2 * CurrNumb - 1
            operators_cr_2(:, length_cr) = 2 * CurrNumb
            call operator_constructor(OrbNumb, Compliance, operators_cr_1, operators_an_1, length, & 
                                      length_cr - 1, length_an - 1, height/2, CurrNumb - 1)      
            call operator_constructor(OrbNumb, Compliance, operators_cr_2, operators_an_2, length, & 
                                      length_cr - 1, length_an - 1, height/2, CurrNumb - 1)
            nullify(operators_an_1, operators_an_2, operators_cr_1, operators_cr_2)
         elseif (Compliance(CurrNumb) .eq. 342) then
            operators_an_1 => operators_an(1 : height/2, :)
            operators_an_2 => operators_an(height/2 + 1 : height, :)
            operators_cr_1 => operators_cr(1 : height/2, :)
            operators_cr_2 => operators_cr(height/2 + 1 : height, :)
            operators_an_1(:, length_an) = 2 * CurrNumb
            operators_cr_1(:, length_cr) = 2 * CurrNumb - 1
            operators_an_2(:, length_an) = 0
            operators_cr_2(:, length_cr) = 0
            call operator_constructor(OrbNumb, Compliance, operators_cr_1, operators_an_1, length, & 
                                      length_cr - 1, length_an - 1, height/2, CurrNumb - 1)      
            call operator_constructor(OrbNumb, Compliance, operators_cr_2, operators_an_2, length, & 
                                      length_cr - 1, length_an - 1, height/2, CurrNumb - 1)
            nullify(operators_an_1, operators_an_2, operators_cr_1, operators_cr_2)
         endif
      endif
      end subroutine

      end module

