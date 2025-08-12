program immobilconversion
  implicit none
  real*8 :: wx, wy, wz, rcutoff, sigmawalkerwalker, sigmachromatinchromatin, sigmawalkerchromatin, cutoffwalkerchromatin 
  real*8 :: xr, yr, zr, r, boxlength, savesreal, walkerreal, npartreal, corr, delta, voxellength, voxeldelta
  real*8 :: totalwalkdens, protradius, x, y, z, kr
  integer :: linecount, startline, saves, i, j, k, n, w, npart, walkercount, nucnumber, xc, yc, zc, deltaint
  integer :: voxeldeltaint, nt, binding, count, type, counter, slabcounter, extractionpoint,srevtype
  real*8 :: walkervol, srevvol, walkradius, srevradius, PI, totalwalkvol, totalsrevvol, quatw, quati, quatj, quatk
  integer, allocatable :: localcount(:), timesteps(:), counts(:)
  integer, allocatable :: bindarray(:)
  real*8, allocatable :: srevx(:), srevy(:), srevz(:), srevquatw(:), srevquati(:), srevquatj(:), srevquatk(:)
  real*8 :: shapex, shapey, shapez, walkshapex, walkshapey, walkshapez, walkquatw, walkquati, walkquatj, walkquatk
  real*8 :: srevshapex, srevshapey, srevshapez
  real*8, allocatable :: walkx(:), walky(:), walkz(:)
  real*8, allocatable :: slabx(:), slaby(:), slabz(:), slabquatw(:), slabquati(:), slabquatj(:), slabquatk(:)
  character(len=100) :: syncommand, wccommand
  character(len=50) :: configfilename, rcutoffstring
  character(len=15) :: tmp, tmp1
  

  boxlength = 130d0 ! in ru
  voxellength = 10d0 ! in ru
  voxeldelta = boxlength/voxellength
  voxeldeltaint = INT(voxeldelta)
  print*,voxeldelta,voxeldeltaint
  srevradius = 1.00d0/2.0d0
  protradius = 1.00d0/2.0d0
  PI = dacos(-1.0d0)
  srevvol = (4.0d0/3.0d0)*PI*(srevradius)**3
  print*,srevradius,PI,srevvol
  ! moving past the lines in the SREV file that don't contain the nucleosome positions
  saves = 21 !
  extractionpoint = 20000
  configfilename = 'edited-config-5.dump'

  walkquatw = 1
  walkquati = 0
  walkquatj = 0
  walkquatk = 0

  walkshapex = 1.0
  walkshapey = 1.0
  walkshapez = 1.0

  srevshapex = 1.0
  srevshapey = 1.0
  srevshapez = 0.5


  counter = 0
  slabcounter = 0

!  open(unit=1,file='dump.srev_and_protein')

  open(unit = 1, file = trim(configfilename))
  do i = 1, 9
     if (i .eq. 4) then
        read(1,*) npart
     else
        read(1,*)
     endif
  enddo
!  npart=534975
  walkercount = 10000
  allocate(srevx(npart), srevy(npart), srevz(npart), srevquatw(npart), srevquati(npart), srevquatj(npart), srevquatk(npart))
  allocate(slabx(npart), slaby(npart), slabz(npart), slabquatw(npart), slabquati(npart), slabquatj(npart), slabquatk(npart))

!  do i = 1, 9
!     read(1,*)
!  enddo

  do i = 1, (npart)
     read(1,*) k, srevtype, x, y, z,kr,kr,kr,kr,kr,kr,kr !, quatw, quati, quatj, quatk, shapex, shapey, shapez
     !if (type .eq. 1) then
        !counter = counter + 1
        srevx(i) = x
        srevy(i) = y
        srevz(i) = z
     !   srevquatw(counter) = quatw
     !   srevquati(counter) = quati
     !   srevquatj(counter) = quatj
     !   srevquatk(counter) = quatk
     !   if (abs(z) .le. 5) then
     !      slabcounter = slabcounter + 1
     !      slabx(slabcounter) = x
     !      slaby(slabcounter) = y
     !      slabz(slabcounter) = z
     !      slabquatw(slabcounter) = quatw
     !      slabquati(slabcounter) = quati
     !      slabquatj(slabcounter) = quatj
     !      slabquatk(slabcounter) = quatk
     !   endif
     !endif
  enddo

  CALL execute_command_line('rm fort.2')
  CALL execute_command_line('ln -s walkers.xyz fort.2') ! read in trajectory file of dye labels


  allocate(walkx(walkercount), walky(walkercount), walkz(walkercount))
  allocate(timesteps(saves), counts(saves))
  allocate(bindarray(walkercount)) ! tracks binding statuses of dye labels
  do j = 1, saves
     print *, "j", j
     read(2,*) walkercount
     read(2,*) tmp, tmp1, nt
     timesteps(j) = nt
     print *, "entering walker loop"
     if (nt .eq. extractionpoint) then
        count = 0
     do i = 1, walkercount
!        print *, "i", i
        read(2,*) k, wx, wy, wz
        if (k .eq. 3) then
           bindarray(i) = k
           walkx(i) = wx
           walky(i) = wy
           walkz(i) = wz
        else
     binding = 2
     do n = 1, npart ! check against SR-EV beads
        xr = wx - srevx(n)
        yr = wy - srevy(n)
        zr = wz - srevz(n)
        xr = xr - boxlength*ANINT(xr/boxlength) ! pbc
        yr = yr - boxlength*ANINT(yr/boxlength)
        zr = zr - boxlength*ANINT(zr/boxlength)
        r = sqrt(xr**2 + yr**2 + zr**2)
        if (r < (srevradius+protradius)) then
 !          print *, "in r if"
           binding = 3
        endif
     enddo

     if (binding .eq. 3) then
        count = count + 1
     endif

     walkx(i) = wx
     walky(i) = wy
     walkz(i) = wz
 !    print *, "binding", binding
     bindarray(i) = binding
     endif
     enddo
     else
        do i = 1, walkercount
           read(2,*) k, wx, wy, wz
        enddo
     endif
enddo

print *, 'finished iterating'

!!do j = 1, saves
!!   count = counts(j)
!!   nt = timesteps(j)
!!   print *, count
!!   print *, nt
!!enddo

print *, "count", count

open(unit = 9, file = "post-immobil-input.in")
do j = 1, saves !! TRADE OUT THE UPPERBOUND FOR SAVES
  if (j .eq. INT(extractionpoint/10002)) then
  do i = 1, walkercount
     write(9,*) "create_atoms", bindarray(i), "single", walkx(i), walky(i), walkz(i)
  enddo
  endif
enddo

print *, 'finished writing'

 deallocate(srevx, srevy, srevz, walkx, walky, walkz)
  
end program immobilconversion

  
