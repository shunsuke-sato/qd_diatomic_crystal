module global_variables

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical parameter
  real(8),parameter :: fs=0.024189d0


! physical system
  integer :: nsite, nelec
  real(8) :: delta_gap, t_hop
  real(8) :: mass, lattice_a, sigma_t_hop
  real(8),allocatable :: potential(:),ham(:,:)
  real(8),allocatable :: t_hop_site(:)
  complex(8),allocatable :: zpsi(:,:)
  complex(8),allocatable :: zpsi_t(:,:),zhpsi_t(:,:)
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
  sigma_t_hop = 1.0d0

! SiO2
!  lattice_a = 5.405d0/0.529d0
!  mass = 0.2d0
!  delta_gap = 9d0/27.2114d0

  write(*,*)"lattice_a=",lattice_a
  write(*,*)"mass     =",mass
  write(*,*)"delta_gap=",delta_gap
  write(*,*)"delta_gap(ev) =",delta_gap*27.2114d0

  t_hop     = 0.5d0/lattice_a*sqrt(delta_gap/mass)
!  t_hop     = 0d0 ! debug
  write(*,*)"t_hop    =",t_hop


!  open(20, file="final_spins.dat", status="old", action="read")
!  read(20, *) nsite
!  allocate(nelem_site(0:nsite-1))
!  do i = 0, nsite-1
!    read(20, *) nelem_site(i)
!  end do
!  close(20)
!  write(*,*) nelem_site


!  nsite = 64
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
  Epulse0 = 1d4*(0.529d-8/27.2114d0) ! MV/cm

  Tpulse0 = 0.5d0*pi*(10d0/fs)/acos((0.5d0)**(1d0/8d0))
  write(*,"(A,2x,e26.16e3)")"Tpulse0 [fs]=",Tpulse0*fs

  Edc0 = 0d0*(0.529d-8/27.2114d0) ! MV/cm
  Tdc0 = 10d0/fs ! fs

! time-propagation
  Tprop = Tpulse0
  dt = 0.1d0
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
  real(8) :: rr,r1,r2,rg1,rg2
! for LAPACK    =====
  integer :: lwork,info
  real(8),allocatable :: work(:)
  real(8),allocatable    :: w(:)

  lwork = min(64, nsite)*max(1,2*nsite-1)+1024
  allocate(w(0:nsite-1))
  allocate(work(lwork))
! for LAPACK    =====

  allocate(potential(0:nsite-1),ham(0:nsite-1,0:nsite-1))
  allocate(t_hop_site(0:nsite-1))
  allocate(zpsi(0:nsite-1,nelec), phi_gs(0:nsite-1, 0:nsite-1),sp_energy(0:nsite-1))
  allocate(zpsi_t(0:nsite-1,1:nelec),zhpsi_t(0:nsite-1,1:nelec))

  do i = 0, nsite-1
    call random_number(r1)
    call random_number(r2)
    if(r1 == 0)then
      rg1 = 0d0
      rg2 = 0d0
    else
      rg1 = sqrt(-2d0*log(r1))*cos(2d0*pi*r2)
      rg2 = sqrt(-2d0*log(r1))*sin(2d0*pi*r2)
    end if
    rg1 = rg1 * sigma_t_hop
    rg2 = rg2 * sigma_t_hop

    t_hop_site(i) = t_hop*(1d0+rg1)
  end do

  potential = 0.5d0*delta_gap*nelem_site

  do j = 0 , nsite-1
    ham(j,j) = potential(j)
    i = mod(j + 1,nsite)
    ham(j,i) = -t_hop_site(j)
    ham(i,j) = -t_hop_site(j)
  end do


  call dsyev('V','U', nsite, ham(0:nsite-1, 0:nsite-1),nsite, &
      sp_energy(0:nsite-1),work,lwork,info)

  zpsi(0:nsite-1,1:nelec) = ham(0:nsite-1, 0:nelec-1)

  write(*,*)"gap=",sp_energy(nsite/2)-sp_energy(nsite/2-1)

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
  write(*,*)"start norm=",sum(abs(zpsi)**2)/nelec
  do it = 0,nt

    call calc_current(current,it)
    write(20,"(999e26.16e3)")it*dt,Afield_t(it),Efield_t(it),current


!    call calc_nex(nex_bloch,nex_houston,nex_pol_houston,it)

    call dt_evolve(it)

  end do
  write(*,*)"final norm=",sum(abs(zpsi)**2)/nelec
  close(20)
  close(21)

end subroutine time_propagation
!-------------------------------------------------------------------------------
subroutine dt_evolve(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer,parameter :: nexp = 4
  integer :: iexp
  real(8) :: dt_t
  complex(8) :: zfact

  dt_t = dt*0.5d0

  zpsi_t = zpsi
  zfact = 1d0
  do iexp = 1, nexp
    zfact = zfact*(-zi*dt_t)/iexp
    call zhpsi(it)
    zpsi = zpsi + zfact*zhpsi_t
    zpsi_t = zhpsi_t
  end do


  zpsi_t = zpsi
  zfact = 1d0
  do iexp = 1, nexp
    zfact = zfact*(-zi*dt_t)/iexp
    call zhpsi(it+1)
    zpsi = zpsi + zfact*zhpsi_t
    zpsi_t = zhpsi_t
  end do


end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine zhpsi(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ib, i, j, ip, im
  complex(8) :: zexp_ac_p,zexp_ac_m

  zexp_ac_m = -exp(-zi*0.5d0*lattice_a*Afield_t(it))
  zexp_ac_p = -exp( zi*0.5d0*lattice_a*Afield_t(it))

!$omp parallel do private(ib,i,ip,im) 
  do ib = 1, nelec

    zhpsi_t(0:nsite-1,ib) = potential(0:nsite-1)*zpsi_t(0:nsite-1,ib)

    i = 0
    ip = i + 1
    im = i - 1 + nsite
    zhpsi_t(i,ib) = zhpsi_t(i,ib) &
        + t_hop_site(i)*zexp_ac_p*zpsi_t(ip,ib) &
        + t_hop_site(im)*zexp_ac_m*zpsi_t(im,ib)


    do i = 1, nsite-1-1
      ip = i + 1
      im = i - 1 
      zhpsi_t(i,ib) = zhpsi_t(i,ib) &
          + t_hop_site(i)*zexp_ac_p*zpsi_t(ip,ib) &
          + t_hop_site(im)*zexp_ac_m*zpsi_t(im,ib)
    end do


    i = nsite-1
    ip = i + 1 - nsite
    im = i - 1
    zhpsi_t(i,ib) = zhpsi_t(i,ib) &
        + t_hop_site(i)*zexp_ac_p*zpsi_t(ip,ib) &
        + t_hop_site(im)*zexp_ac_m*zpsi_t(im,ib)



  end do

end subroutine zhpsi
!-------------------------------------------------------------------------------
subroutine calc_current(jt_t,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: jt_t
  complex(8),allocatable :: zp_psi_t(:,:)
  integer :: i,ip,im, ib
  complex(8) :: zexp_ac_p,zexp_ac_m

  allocate(zp_psi_t(0:nsite-1,1:nelec))

  zexp_ac_m =-exp(-zi*0.5d0*lattice_a*Afield_t(it))*(-zi*0.5d0*lattice_a)
  zexp_ac_p =-exp( zi*0.5d0*lattice_a*Afield_t(it))*( zi*0.5d0*lattice_a)

!$omp parallel do private(ib,i,ip,im) 
  do ib = 1, nelec

    i = 0
    ip = i + 1
    im = i - 1 + nsite

    zp_psi_t(i,ib) = t_hop_site(i)*zexp_ac_p*zpsi_t(ip,ib) &
                    +t_hop_site(im)*zexp_ac_m*zpsi_t(im,ib) 



    do i = 1, nsite-1-1
      ip = i + 1
      im = i - 1
      zp_psi_t(i,ib) = t_hop_site(i)*zexp_ac_p*zpsi_t(ip,ib) &
                      +t_hop_site(im)*zexp_ac_m*zpsi_t(im,ib) 


    end do


    i = nsite-1
    ip = i + 1 - nsite
    im = i - 1
    zp_psi_t(i,ib) = t_hop_site(i)*zexp_ac_p*zpsi_t(ip,ib) &
                    +t_hop_site(im)*zexp_ac_m*zpsi_t(im,ib) 



  end do


  jt_t = sum(conjg(zpsi)*zp_psi_t)/(nelec*nsite)

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
