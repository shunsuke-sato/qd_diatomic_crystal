module global_variables

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical parameter
  real(8),parameter :: fs=0.024189d0


! physical system
  integer :: nsite, nelec
  real(8) :: delta_gap, t_hop
  real(8) :: mass, lattice_a
  real(8),allocatable :: potential(:),ham(:,:)
  complex(8),allocatable :: zpsi(:,:)
  real(8),allocatable :: phi_gs(:,:),sp_energy(:)
  integer,allocatable :: nelem_site(:)


! time propagation
  integer :: nt
  real(8) :: Tprop, dt

! laser fields
  real(8) :: omega0, Epulse0, Tpulse0
  real(8) :: Edc0, Tdc0
  real(8),allocatable :: Efield_t(:), Afield_t(:)


end module global_variables
!-------------------------------------------------------------------------------
program main
  use global_variables
  implicit none

  call input
  call preparation

!  stop
  call time_propagation

end program main
!-------------------------------------------------------------------------------
subroutine input
  use global_variables
  implicit none
  integer :: i,j

! system parameters
! GaAs
  lattice_a = 5.65d0/0.529d0
  mass = 1d0/(1d0/0.067d0+1d0/0.08d0)
  delta_gap = 1.52d0/27.2114d0

! SiO2
!  lattice_a = 5.405d0/0.529d0
!  mass = 0.2d0
!  delta_gap = 9d0/27.2114d0

  write(*,*)"lattice_a=",lattice_a
  write(*,*)"mass     =",mass
  write(*,*)"delta_gap=",delta_gap

  t_hop     = 0.5d0/lattice_a*sqrt(delta_gap/mass)
  write(*,*)"t_hop    =",t_hop


!  open(20, file="final_spins.dat", status="old", action="read")
!  read(20, *) nsite
!  allocate(nelem_site(0:nsite-1))
!  do i = 0, nsite-1
!    read(20, *) nelem_site(i)
!  end do
!  close(20)
!  write(*,*) nelem_site


  nsite = 2048
  allocate(nelem_site(0:nsite-1))
  do j = 0, nsite-1
    nelem_site(j) = -(-1)**j 
  end do
  nelec = 0
  do j = 0, nsite-1
    if(nelem_site(j) == -1)nelec = nelec + 1
  end do

  write(*,*)"nsite=",nsite
  write(*,*)"nelec=",nelec

! laser fields
  omega0 = 1.55d0/27.2114d0 ! 3.5 mum
  Epulse0 = 0d4*(0.529d-8/27.2114d0) ! MV/cm

  Tpulse0 = 0.5d0*pi*(10d0/fs)/acos((0.5d0)**(1d0/8d0))
  write(*,"(A,2x,e26.16e3)")"Tpulse0 [fs]=",Tpulse0*fs

  Edc0 = 1d0*(0.529d-8/27.2114d0) ! MV/cm
  Tdc0 = 10d0/fs ! fs

! time-propagation
  Tprop = Tpulse0
  dt = 0.2d0
  nt = aint(Tprop/dt)+1
  dt = Tprop/nt
  write(*,"(A,2x,e26.16e3)")"refined dt=",dt
  write(*,"(A,2x,I7)")"nt        =",nt


end subroutine input
!-------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer,parameter :: nw = 100
  integer :: j,i, iw
  real(8) :: max_ene, min_ene, ww, ss, dw
  real(8) :: dos(0:nw), dos_s(0:nw)
  real(8) :: rr
! for LAPACK    =====
  integer :: lwork,info
  real(8),allocatable :: work(:)
  real(8),allocatable    :: w(:)

  lwork = min(64, nsite)*max(1,2*nsite-1)+1024
  allocate(w(0:nsite-1))
  allocate(work(lwork))
! for LAPACK    =====

  allocate(potential(0:nsite-1),ham(0:nsite-1,0:nsite-1))
  allocate(zpsi(0:nsite-1,nelec), phi_gs(0:nsite-1, 0:nsite-1),sp_energy(0:nsite-1))


  potential = 0.5d0*delta_gap*nelem_site

  do j = 0 , nsite-1
    call random_number(rr)
    rr = (rr-0.5d0)
    ham(j,j) = potential(j)
    i = mod(j + 1,nsite)
    ham(j,i) = t_hop*(1d0+rr)
    i = mod(j - 1 + nsite,nsite)
    ham(j,i) = t_hop*(1d0+rr)
  end do


  call dsyev('V','U', nsite, ham(0:nsite-1, 0:nsite-1),nsite, &
      sp_energy(0:nsite-1),work,lwork,info)

!  write(*,*)sp_energy
  max_ene = maxval(sp_energy)
  min_ene = minval(sp_energy)
  dw = (max_ene-min_ene)/nw
  max_ene = max_ene + dw*10
  min_ene = min_ene - dw*10
  dw = (max_ene-min_ene)/nw


  dos = 0d0
  do j = 0, nsite-1
!    w = min_ene + iw*dw
    iw = aint( (sp_energy(j)-min_ene)/dw )
    dos(iw) = dos(iw) + 1d0/dw
  end do
  ss = sum(dos)*dw
  dos = dos/ss

  dos_s = 0d0
  do i = 0, nw
    
    ss = 0d0
    do j = 0, nw
      ss = ss + dos(j)*exp(-0.5d0*(i-j)**2)
    end do
    dos_s(i) = ss
  end do
  ss = sum(dos_s)*dw
  dos_s = dos_s/ss

  open(30,file="dos.out")
  do iw = 0, nw
    ww = min_ene + iw*dw
    write(30,"(99e26.16e3)")ww, dos(iw), dos_s(iw)
  end do
  close(30)


end subroutine preparation
!-------------------------------------------------------------------------------
subroutine time_propagation
  use global_variables
  implicit none
  integer :: it
  real(8) :: current
  real(8) :: nex_bloch,nex_houston,nex_pol_houston

  call init_laser_field

  open(20,file='current.out')
  open(21,file='nex.out')
  do it = 0,nt

    call calc_current(current,it)
    write(20,"(999e26.16e3)")it*dt,Afield_t(it),Efield_t(it),current

!    call calc_nex(nex_bloch,nex_houston,nex_pol_houston,it)

    call dt_evolve(it)

  end do
  close(20)
  close(21)

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it


end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine calc_current(jt_t,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: jt_t


end subroutine calc_current
!-------------------------------------------------------------------------------
subroutine init_laser_field
  use global_variables
  implicit none
  integer :: it
  real(8) :: tt, ss

  allocate(Efield_t(-1:nt+1),Afield_t(-1:nt+2))
  Efield_t = 0d0
  Afield_t = 0d0


  if(Epulse0 /= 0d0)then
    do it = 0, nt+2
      tt = dt*it
      ss = (tt - 0.5d0*Tpulse0)
      if(abs(ss)<= 0.5d0*Tpulse0)then
        Afield_t(it) = Afield_t(it) &
            -(Epulse0/omega0)*sin(omega0*ss)*cos(pi*ss/Tpulse0)**4
      end if
    end do
  end if

  if(Edc0 /= 0d0)then
    do it = 0, nt+2
      tt = dt*it
      if(tt<= Tdc0)then
        ss = tt/Tdc0
        Afield_t(it) = Afield_t(it) &
            -Edc0*Tdc0*(ss**3 -0.5d0*ss**4)
      else
        Afield_t(it) = Afield_t(it) &
            -Edc0*(tt-Tdc0) -0.5d0*Edc0*Tdc0
      end if
    end do
  end if


  do it = 0, nt+1
    Efield_t(it) = -0.5d0*(Afield_t(it+1)-Afield_t(it-1))/dt
  end do

end subroutine init_laser_field
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
