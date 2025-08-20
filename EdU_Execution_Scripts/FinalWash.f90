program finalwash
  implicit none
  real*8 :: wx, wy, wz, rcutoff, sigmawalkerwalker, sigmachromatinchromatin, sigmawalkerchromatin, cutoffwalkerchromatin 
  real*8 :: xr, yr, zr, r, boxlength, savesreal, walkerreal, npartreal
  real*8 :: totalwalkdens, protradius, x, y, z, kr
  integer :: linecount, startline, saves, i, j, k, n, w, npart, walkercount, xc, yc, zc
  integer :: voxeldeltaint, nt, binding, extractionpoint
  real*8 :: walkervol, srevvol, walkradius, srevradius, PI
  integer, allocatable :: localcount(:), timesteps(:), counts(:)
  integer, allocatable :: bindarray(:)
  real*8, allocatable :: srevx(:), srevy(:), srevz(:), srevquatw(:), srevquati(:), srevquatj(:), srevquatk(:)
  real*8 :: shapex, shapey, shapez, walkshapex, walkshapey, walkshapez
  real*8 :: srevshapex, srevshapey, srevshapez
  real*8, allocatable :: walkx(:), walky(:), walkz(:)
  character(len=50) :: configfilename
  character(len=15) :: tmp, tmp1
  

  boxlength = 130d0 ! in ru; CHANGE THIS TO MATCH BOXLENGTH PARAMS IN SREV CONFIG FILE
  srevradius = 1.00d0/2.0d0
  protradius = 0.20d0/2.0d0
  PI = dacos(-1.0d0)
  saves = 21
  extractionpoint = 20000
  configfilename = 'edited-config-5.dump'

  open(unit = 1, file = trim(configfilename))
  do i = 1, 9
     if (i .eq. 4) then
        read(1,*) npart
     else
        read(1,*)
     endif
  enddo
  walkercount = 10000
  allocate(srevx(npart), srevy(npart), srevz(npart))

  do i = 1, (npart)
     read(1,*) k, srevtype, x, y, z,kr,kr,kr,kr,kr,kr,kr
     srevx(i) = x
     srevy(i) = y
     srevz(i) = z
  enddo

  CALL execute_command_line('rm fort.2')
  CALL execute_command_line('ln -s walkers.xyz fort.2')


  allocate(walkx(walkercount), walky(walkercount), walkz(walkercount))
  allocate(timesteps(saves), counts(saves))
  allocate(bindarray(walkercount))
  do j = 1, saves
     print *, "j", j
     read(2,*) walkercount
     read(2,*) tmp, tmp1, nt
     print *, "entering walker loop"
     if (nt .eq. extractionpoint) then
     do i = 1, walkercount
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
        xr = xr - boxlength*ANINT(xr/boxlength) ! adjust computed distance for periodic boundary conditions in sim
        yr = yr - boxlength*ANINT(yr/boxlength)
        zr = zr - boxlength*ANINT(zr/boxlength)
        r = sqrt(xr**2 + yr**2 + zr**2)
        if (r < (srevradius+protradius)) then
           binding = 3
        endif
     enddo
     walkx(i) = wx
     walky(i) = wy
     walkz(i) = wz
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

open(unit = 9, file = "EdU-postwash-config-7.xyz")
do i = 1, walkercount
	if (bindarray(i) .eq. 3) then
     	write(9,*) bindarray(i), walkx(i), walky(i), walkz(i)
	endif
enddo

print *, 'finished writing'

deallocate(srevx, srevy, srevz, walkx, walky, walkz, bindarray)
  
end program finalwash

  
