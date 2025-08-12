program myprogram
  implicit none
  integer :: i
  real*8 :: xr, yr, zr, sigmawalkerwalker, sigmachromatinchromatin, sigmawalkerchromatin, totalfrac
  real*8, allocatable :: allx(:)
  real*8, allocatable :: ally(:)
  real*8, allocatable :: allz(:)
  real*8 :: boxlength,xlength,ylength,zlength,volume, xlowerbound, ylowerbound, zlowerbound
  integer :: npart
  integer :: j
  real*8 :: r
  integer :: startline
  integer :: k
  integer :: linecount
  integer :: walkers
  integer, allocatable :: seed(:)
  integer :: s
  real*8 :: randnum
  real :: pi
  real*8, allocatable :: xi(:)
  real*8, allocatable :: yi(:)
  real*8, allocatable :: zi(:)
  real*8 :: xt
  real*8 :: yt
  real*8 :: zt
  character(len=50) :: configfilename
  character(len=100) :: wccommand
  
  configfilename = 'config-relaxed-5.dump'
  pi = dacos(-1.0d0)
  xlength=200
  ylength=200
  zlength=200
  xlowerbound = -100
  ylowerbound = -100
  zlowerbound = -100
  !totalfrac = 0.12001
  sigmawalkerwalker = 1.0
  sigmachromatinchromatin = 1.0
  sigmawalkerchromatin = (sigmawalkerwalker+sigmachromatinchromatin)/2.0
  walkers = 10000

  CALL execute_command_line('rm walkerinput.in')
  wccommand = trim('wc -l < ' // trim(configfilename) // ' > wc.txt')
  print *, wccommand
  CALL execute_command_line(wccommand)
  OPEN(unit=8, file='wc.txt')
  READ(8, *) linecount
  print *, linecount
  rewind(8)
  CALL execute_command_line('rm wc.txt')
  close(8)
  open(unit = 4, file = trim(configfilename))
  npart = 0
  startline = 10
  do i = 1, startline-1
     if (i == 4) then
        read(4,*) npart
     else
        read(4,*)
     endif
  enddo

  volume=xlength*ylength*zlength
  !walkers = anint(((totalfrac*volume)*0.75*(1/pi)-npart*(0.5*sigmachromatinchromatin)**3)*(1/(0.5*sigmawalkerwalker)**3))

  print *, 'sigmawalkerwalker:', sigmawalkerwalker

  print *, 'npart:', npart
  print *, 'nwalkers:', walkers

  allocate (allx(npart), ally(npart), allz(npart))

  do i = 1, npart
     read(4,*) k,k,allx(i),ally(i),allz(i)
     !print *, allx(i), ally(i), allz(i)
     if (allx(i) .gt. (xlowerbound+xlength)) then
        print *, 'out of bounds'
     else if (allx(i) .lt. xlowerbound) then
        print *, 'out of bounds'
     else if (ally(i) .gt. (ylowerbound+ylength)) then
        print *, 'out of bounds'
     else if (ally(i) .lt. ylowerbound) then
        print *, 'out of bounds'
     else if (allz(i) .gt. (zlowerbound+zlength)) then
        print *, 'out of bounds'
     else if (allz(i) .lt. zlowerbound) then
        print *, 'out of bounds'
     endif
  enddo

  close(4)
  
!  call exit(0)

  allocate(xi(walkers))
  allocate(yi(walkers))
  allocate(zi(walkers))

  print *, 'Starting random seed call'

  CALL RANDOM_SEED (size=s) ! using RANDOM_SEED & RANDOM_NUMBER() to generate the random number since some sites recommend that over RAND: https://gcc.gnu.org/onlinedocs/gfortran/RAND.html
  allocate(seed(s)) ! RANDOM_SEED has thing where seed must be array of some minimum size I won't know until I call RANDOM_SEED
  seed = 1

  open(5, file = 'walkerinput.in', status = 'new')
  CALL RANDOM_SEED(PUT=seed)
  do i = 1, walkers
     if (MOD(i,100) .eq. 0) then
        print *, i
     endif
71   CALL RANDOM_NUMBER(randnum)  ! 71 is label for GOTO statements
     xt = randnum*xlength + xlowerbound
     CALL RANDOM_NUMBER(randnum)
     yt = randnum*ylength + ylowerbound
     CALL RANDOM_NUMBER(randnum)
     zt = randnum*zlength + zlowerbound
     do j = 1, npart ! check against SR-EV beads
        xr = xt - allx(j)
        yr = yt - ally(j)
        zr = zt - allz(j)
        xr = xr - xlength*ANINT(xr/xlength) ! pbc
        yr = yr - ylength*ANINT(yr/ylength)
        zr = zr - zlength*ANINT(zr/zlength)
        r = sqrt(xr**2 + yr**2 + zr**2)
        if (r < (0.9*sigmawalkerchromatin)) then
           go to 71
        endif
     enddo
     do k = 1, i-1 ! check against already initialized walkers
        xr = xt-xi(k)
        yr = yt-yi(k)
        zr = zt-zi(k)
        xr = xr - xlength*ANINT(xr/xlength) ! pbc
        yr = yr - ylength*ANINT(yr/ylength)
        zr = zr - zlength*ANINT(zr/zlength)
        r = sqrt(xr**2 + yr**2 + zr**2)
        if (r < (0.9*sigmawalkerwalker)) then
           go to 71
        endif
     enddo
     xi(i) = xt
     yi(i) = yt
     zi(i) = zt
     write(5,fmt="(a,1x,i7.1,1x,a,1x,f16.8,1x,f16.8,1x,f16.8,1x)") "create_atoms", 2, "single", xt,yt,zt
  enddo

  close(5)

end program myprogram
