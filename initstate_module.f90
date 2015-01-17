      module InitialStateData
      integer, save :: InitExOpNumb
      integer, save :: MaxInitDim
      logical, allocatable, save, dimension(:) :: cf_q_init_logical
      integer, allocatable, save, dimension(:) :: k_q_init
      real*8, allocatable, save, dimension(:) :: ksi
      integer, allocatable, save, dimension(:, :) :: q_init
      end module