module global_variables

! mathematical constants
  real(8),parameter :: pi = 4d0*atan(1d0)
  complex(8),parameter :: zi = (0d0, 1d0)

! physical parameter
  real(8),parameter :: fs=0.024189d0


! physical system
  integer :: nk
  real(8) :: delta_gap, t_hop
  real(8) :: mass, lattice_a
  real(8),allocatable :: ham_kt(:,:,:)
  complex(8),allocatable :: zpsi(:,:)
  real(8),allocatable :: phi_gs(:,:,:),sp_energy(:,:)
  real(8),allocatable :: kn(:)


! time propagation
  integer :: nt
  real(8) :: Tprop, dt

! laser fields
  real(8) :: omega0, Epulse0, Tpulse0
  real(8) :: Edc0, Tdc0
  real(8),allocatable :: Efield_t(:), Afield_t(:)

! quantum master equation
  logical :: if_q_master_eq 
  real(8) :: T1_qm, T2_qm
  complex(8),allocatable :: zrho_k(:,:,:)
  integer, parameter :: N_PROJ_BLOCH = 0, &
                        N_PROJ_HOUSTON = 1, &
                        N_PROJ_P_HOUSTON = 2
  integer,parameter :: N_PROJ_METHOD = N_PROJ_P_HOUSTON

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

! system parameters
! GaAs
  lattice_a = 5.65d0/0.529d0
  mass = 1d0/(1d0/0.067d0+1d0/0.08d0)
  delta_gap = 1.52d0/27.2114d0

  write(*,*)"lattice_a=",lattice_a
  write(*,*)"mass     =",mass
  write(*,*)"delta_gap=",delta_gap

  t_hop     = 0.5d0/lattice_a*sqrt(delta_gap/mass)
  write(*,*)"t_hop    =",t_hop

  nk = 2

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


! quantum master equation
  if_q_master_eq = .false.
  T1_qm = 10d0/fs
  T2_qm = 10d0/fs

end subroutine input
!-------------------------------------------------------------------------------
subroutine preparation
  use global_variables
  implicit none
  integer :: ik
  real(8) :: ham(2,2)

  allocate(ham_kt(2,2,0:nk-1))
  allocate(zpsi(2,0:nk-1), phi_gs(2,2,0:nk-1),sp_energy(2,0:nk-1))
  allocate(kn(0:nk-1))

  do ik = 0, nk-1
    kn(ik) = (2d0/lattice_a)*(pi/nk)*ik
  end do


  do ik = 0, nk-1
    ham(1,1) = -0.5d0*delta_gap
    ham(1,2) = -2d0*t_hop*cos(0.5d0*lattice_a*kn(ik))
    ham(2,1) = ham(1,2)
    ham(2,2) = 0.5d0*delta_gap

    call diag_2x2(ham, phi_gs(:,:,ik), sp_energy(:,ik))
    if(sp_energy(1,ik) > sp_energy(2,ik))then
      zpsi(:,ik) = phi_gs(:,2,ik)
    else
      zpsi(:,ik) = phi_gs(:,1,ik)
    end if
    write(*,*)phi_gs(:,2,ik) ! debug
  end do


  if(if_q_master_eq)then
    allocate(zrho_k(2,2,0:nk-1))
    do ik = 0, nk-1
      if(sp_energy(1,ik) > sp_energy(2,ik))then
        zrho_k(1,1,ik) = phi_gs(1,2,ik)*phi_gs(1,2,ik)
        zrho_k(2,1,ik) = phi_gs(2,2,ik)*phi_gs(1,2,ik)
        zrho_k(1,2,ik) = phi_gs(1,2,ik)*phi_gs(2,2,ik)
        zrho_k(2,2,ik) = phi_gs(2,2,ik)*phi_gs(2,2,ik)
      else
        zrho_k(1,1,ik) = phi_gs(1,1,ik)*phi_gs(1,1,ik)
        zrho_k(2,1,ik) = phi_gs(2,1,ik)*phi_gs(1,1,ik)
        zrho_k(1,2,ik) = phi_gs(1,1,ik)*phi_gs(2,1,ik)
        zrho_k(2,2,ik) = phi_gs(2,1,ik)*phi_gs(2,1,ik)
      end if
    end do
  end if

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

    if(.not. if_q_master_eq)then
      call calc_nex(nex_bloch,nex_houston,nex_pol_houston,it)
      write(21,"(999e26.16e3)")it*dt,nex_bloch,nex_houston,nex_pol_houston
    end if

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

  if(if_q_master_eq)then
    call dt_evolve_q_master(it)
  else
    call dt_evolve_tdse(it)
  end if

end subroutine dt_evolve
!-------------------------------------------------------------------------------
subroutine dt_evolve_q_master(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ik
  real(8) :: ham(2,2), vec(2,2), lambda(2)
  complex(8) :: zvec(2), zrho_tmp(2,2), zrho_tmp_delta(2,2)
  real(8) :: rho_eq(2,2)
  real(8) :: tt, kt, dkdt
  complex(8) :: zref_state(2,2)
  real(8) :: ref_energy(2)

  do ik = 0, nk-1

! propagation from dt*it to dt*(it+0.5)
    tt = dt*it

    kt = kn(ik) + Afield_t(it)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)

    zrho_tmp = matmul(transpose(vec),matmul(zrho_k(:,:,ik),vec(:,:)))

    zrho_tmp(1,2) = zrho_tmp(1,2)*exp(zi*(lambda(2)-lambda(1))*dt*0.5d0)
    zrho_tmp(2,1) = conjg(zrho_tmp(1,2))

    zrho_k(:,:,ik) = matmul(vec, matmul(zrho_tmp,transpose(vec)))


! relaxation from dt*it to dt*(it+0.5)
    tt = dt*it
    kt = kn(ik) + Afield_t(it)
    dkdt = 0.5d0*(Afield_t(it+1)-Afield_t(it-1))/dt

    call calc_reference_states(zref_state, ref_energy, kn(ik), kt, dkdt)

    zrho_tmp = matmul(transpose(conjg(zref_state)),matmul(zrho_k(:,:,ik),zref_state(:,:)))

    rho_eq = 0d0
    if(ref_energy(1) >= ref_energy(2))then
      rho_eq(2,2) = 1d0
    else
      rho_eq(1,1) = 1d0
    end if
    zrho_tmp_delta = zrho_tmp - rho_eq
    zrho_tmp_delta(1,1) = zrho_tmp_delta(1,1)*exp(-0.5d0*dt/T1_qm)
    zrho_tmp_delta(2,1) = zrho_tmp_delta(2,1)*exp(-0.5d0*dt/T2_qm)
    zrho_tmp_delta(1,2) = zrho_tmp_delta(1,2)*exp(-0.5d0*dt/T2_qm)
    zrho_tmp_delta(2,2) = zrho_tmp_delta(2,2)*exp(-0.5d0*dt/T1_qm)

    zrho_tmp = zrho_tmp_delta + rho_eq
    zrho_k(:,:,ik) = matmul(zref_state, matmul(zrho_tmp,transpose(conjg(zref_state))))


! relaxation from dt*(it+0.5) to dt*(it+1)
    tt = dt*(it+1)
    kt = kn(ik) + Afield_t(it+1)
    dkdt = 0.5d0*(Afield_t(it+1+1)-Afield_t(it-1+1))/dt

    call calc_reference_states(zref_state, ref_energy, kn(ik), kt, dkdt)

    zrho_tmp = matmul(transpose(conjg(zref_state)),matmul(zrho_k(:,:,ik),zref_state(:,:)))

    rho_eq = 0d0
    if(ref_energy(1) >= ref_energy(2))then
      rho_eq(2,2) = 1d0
    else
      rho_eq(1,1) = 1d0
    end if
    zrho_tmp_delta = zrho_tmp - rho_eq
    zrho_tmp_delta(1,1) = zrho_tmp_delta(1,1)*exp(-0.5d0*dt/T1_qm)
    zrho_tmp_delta(2,1) = zrho_tmp_delta(2,1)*exp(-0.5d0*dt/T2_qm)
    zrho_tmp_delta(1,2) = zrho_tmp_delta(1,2)*exp(-0.5d0*dt/T2_qm)
    zrho_tmp_delta(2,2) = zrho_tmp_delta(2,2)*exp(-0.5d0*dt/T1_qm)

    zrho_tmp = zrho_tmp_delta + rho_eq
    zrho_k(:,:,ik) = matmul(zref_state, matmul(zrho_tmp,transpose(conjg(zref_state))))

! propagation from dt*(it+0.5) to dt*(it+1)
    tt = dt*(it+1)

    kt = kn(ik) + Afield_t(it+1)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)

    zrho_tmp = matmul(transpose(vec),matmul(zrho_k(:,:,ik),vec(:,:)))

    zrho_tmp(1,2) = zrho_tmp(1,2)*exp(zi*(lambda(2)-lambda(1))*dt*0.5d0)
    zrho_tmp(2,1) = conjg(zrho_tmp(1,2))

    zrho_k(:,:,ik) = matmul(vec, matmul(zrho_tmp,transpose(vec)))


  end do

end subroutine dt_evolve_q_master
!-------------------------------------------------------------------------------
subroutine calc_reference_states(zref_state, ref_energy, k0, kt, dkdt)
  use global_variables
  implicit none
  complex(8),intent(out) :: zref_state(2,2)
  real(8),intent(out) :: ref_energy(2)
  real(8),intent(in) :: k0, kt, dkdt
  real(8) :: ref_state(2,2)

  select case(N_PROJ_METHOD)
  case(N_PROJ_BLOCH)
    call calc_bloch_states(ref_state, ref_energy, k0)
    zref_state = ref_state
  case(N_PROJ_HOUSTON)
    call calc_houston_states(ref_state, ref_energy, kt)
    zref_state = ref_state
  case(N_PROJ_P_HOUSTON)
    call calc_pol_houston_states(zref_state, ref_energy, kt, dkdt)
  case default
    stop 'Error in calc_reference_states'
  end select

end subroutine calc_reference_states
!-------------------------------------------------------------------------------
subroutine calc_bloch_states(ref_state, ref_energy, k0)
  use global_variables
  implicit none
  real(8),intent(out) :: ref_state(2,2), ref_energy(2)
  real(8),intent(in) :: k0
  real(8) :: ham(2,2)

  ham(1,1) = -0.5d0*delta_gap
  ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*k0)
  ham(1,2) = ham(2,1)
  ham(2,2) =  0.5d0*delta_gap
  
  call diag_2x2(ham, ref_state, ref_energy)
  

end subroutine calc_bloch_states
!-------------------------------------------------------------------------------
subroutine calc_houston_states(ref_state, ref_energy, kt)
  use global_variables
  implicit none
  real(8),intent(out) :: ref_state(2,2), ref_energy(2)
  real(8),intent(in) :: kt
  real(8) :: ham(2,2)

  ham(1,1) = -0.5d0*delta_gap
  ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
  ham(1,2) = ham(2,1)
  ham(2,2) =  0.5d0*delta_gap
  
  call diag_2x2(ham, ref_state, ref_energy)
  
  
end subroutine calc_houston_states
!-------------------------------------------------------------------------------
subroutine calc_pol_houston_states(zref_state, ref_energy, kt, dkdt)
  use global_variables
  implicit none
  complex(8),intent(out) :: zref_state(2,2)
  real(8),intent(out) :: ref_energy(2)
  real(8),intent(in) :: kt, dkdt
  real(8) :: phi, xx, yy, eps_c, eps_v, factor
  real(8) :: duc_dk(2), uv(2), uc(2)
  complex(8) :: zham(2,2), zvec(2,2)

  phi = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
  xx =  phi/(0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))
  yy = -phi/(0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))

  eps_c =  sqrt(delta_gap**2/4d0 + phi**2)
  eps_v = -sqrt(delta_gap**2/4d0 + phi**2)

  duc_dk(1)=-xx/(sqrt(1d0+xx**2))**3*xx + 1d0/sqrt(1d0+xx**2)**3
  duc_dk(2)=-xx/(sqrt(1d0+xx**2))**3 
  factor = 1d0/(0.5d0*delta_gap + sqrt(delta_gap**2/4d0 + phi**2))
  factor = factor - phi**2/( &
      (0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))**2 &
      *sqrt(delta_gap**2/4d0+phi**2) &
      )

  factor = factor *lattice_a*t_hop*sin(0.5d0*lattice_a*kt)

  duc_dk = duc_dk*factor

  uc(1) =  xx/sqrt(1d0+xx**2)
  uc(2) = 1d0/sqrt(1d0+xx**2)

  uv(1) = 1d0/sqrt(1d0+yy**2)
  uv(2) =  yy/sqrt(1d0+yy**2)
    
  zham(1,1) = eps_v
  zham(1,2) = -zi*sum(uv*duc_dk)*dkdt
  zham(2,1) = conjg(zham(1,2))
  zham(2,2) = eps_c

  call diag_2x2_complex(zham, zvec, ref_energy)

  zref_state(:,1) = uv*zvec(1,1) + uc*zvec(2,1)
  zref_state(:,2) = uv*zvec(1,2) + uc*zvec(2,2)


end subroutine calc_pol_houston_states
!-------------------------------------------------------------------------------
subroutine dt_evolve_tdse(it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  integer :: ik
  real(8) :: ham(2,2), vec(2,2), lambda(2)
  complex(8) :: zvec(2)
  real(8) :: tt, kt




  do ik = 0, nk-1

! propagation from dt*it to dt*(it+0.5)
    tt = dt*it

    kt = kn(ik) + Afield_t(it)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)
    
    zvec(1) = vec(1,1)*zpsi(1,ik)+ vec(2,1)*zpsi(2,ik)
    zvec(2) = vec(1,2)*zpsi(1,ik)+ vec(2,2)*zpsi(2,ik)
    zvec(1) = exp(-zi*0.5*dt*lambda(1))*zvec(1)
    zvec(2) = exp(-zi*0.5*dt*lambda(2))*zvec(2)

    zpsi(1,ik) = vec(1,1)*zvec(1) + vec(1,2)*zvec(2)
    zpsi(2,ik) = vec(2,1)*zvec(1) + vec(2,2)*zvec(2)


! propagation from dt*(it+0.5) to dt*(it+1)
    tt = dt*(it+1)

    kt = kn(ik) + Afield_t(it+1)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)
    
    zvec(1) = vec(1,1)*zpsi(1,ik)+ vec(2,1)*zpsi(2,ik)
    zvec(2) = vec(1,2)*zpsi(1,ik)+ vec(2,2)*zpsi(2,ik)
    zvec(1) = exp(-zi*0.5*dt*lambda(1))*zvec(1)
    zvec(2) = exp(-zi*0.5*dt*lambda(2))*zvec(2)

    zpsi(1,ik) = vec(1,1)*zvec(1) + vec(1,2)*zvec(2)
    zpsi(2,ik) = vec(2,1)*zvec(1) + vec(2,2)*zvec(2)

  end do


end subroutine dt_evolve_tdse
!-------------------------------------------------------------------------------
subroutine calc_current(jt_t,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: jt_t

  if(if_q_master_eq)then
    call calc_current_q_master(jt_t,it)
  else
    call calc_current_tdse(jt_t,it)
  end if

end subroutine calc_current
!-------------------------------------------------------------------------------
subroutine calc_current_q_master(jt_t,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: jt_t
  integer :: ik
  real(8) :: pmat
  real(8) :: tt, kt
  complex(8) :: zmat(2,2)

  jt_t = 0d0
  do ik = 0, nk-1

    kt = kn(ik) + Afield_t(it)
    pmat = 2d0*t_hop*sin(0.5d0*lattice_a*kt)*0.5d0*lattice_a
    zmat = 0d0
    zmat(1,2) = pmat
    zmat(2,1) = pmat
    zmat = matmul(zmat, zrho_k(:,:,ik))

    jt_t = jt_t + real(zmat(1,1)+zmat(2,2))
  end do

  jt_t = jt_t/nk


end subroutine calc_current_q_master
!-------------------------------------------------------------------------------
subroutine calc_current_tdse(jt_t,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: jt_t
  integer :: ik
  real(8) :: pmat
  real(8) :: tt, kt

  tt = dt*it


  jt_t = 0d0
  do ik = 0, nk-1

    kt = kn(ik) + Afield_t(it)
    pmat = 2d0*t_hop*sin(0.5d0*lattice_a*kt)*0.5d0*lattice_a

    jt_t = jt_t + conjg(zpsi(1,ik))*pmat*zpsi(2,ik) + conjg(zpsi(2,ik))*pmat*zpsi(1,ik)
  end do

  jt_t = jt_t/nk

end subroutine calc_current_tdse
!-------------------------------------------------------------------------------
subroutine calc_nex(nex_bloch,nex_houston,nex_pol_houston,it)
  use global_variables
  implicit none
  integer,intent(in) :: it
  real(8),intent(out) :: nex_bloch, nex_houston,nex_pol_houston
  integer :: ik
  real(8) :: kt, dkdt, phi, xx, yy, factor
  real(8) :: ham(2,2), vec(2,2), lambda(2)
  real(8) :: duc_dk(2), uv(2), uc(2)
  real(8) :: eps_v, eps_c
  complex(8) :: zham(2,2), zvec(2,2)
  complex(8) :: zp_houston_states(2,2)

! Bloch projection
  nex_bloch = 0d0
  do ik = 0, nk-1
    nex_bloch = nex_bloch &
        + abs(phi_gs(1,1,ik)*zpsi(1,ik)+phi_gs(2,1,ik)*zpsi(2,ik))**2
  end do
  nex_bloch = nex_bloch/nk

! Houston projection
  nex_houston = 0d0
  do ik = 0, nk-1

    kt = kn(ik) + Afield_t(it)
    ham(1,1) = -0.5d0*delta_gap
    ham(2,1) = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    ham(1,2) = ham(2,1)
    ham(2,2) =  0.5d0*delta_gap

    call diag_2x2(ham, vec, lambda)

    nex_houston = nex_houston + abs(vec(1,1)*zpsi(1,ik)+vec(2,1)*zpsi(2,ik))**2
  end do
  nex_houston = nex_houston/nk

! polarized Houston projection
  nex_pol_houston = 0d0
  do ik = 0, nk-1

    kt = kn(ik) + Afield_t(it)
    dkdt = 0.5d0*(Afield_t(it+1)-Afield_t(it-1))/dt
    phi = -2d0*t_hop*cos(0.5d0*lattice_a*kt)
    xx =  phi/(0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))
    yy = -phi/(0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))

    eps_c =  sqrt(delta_gap**2/4d0 + phi**2)
    eps_v = -sqrt(delta_gap**2/4d0 + phi**2)

    duc_dk(1)=-xx/(sqrt(1d0+xx**2))**3*xx + 1d0/sqrt(1d0+xx**2)**3
    duc_dk(2)=-xx/(sqrt(1d0+xx**2))**3 
    factor = 1d0/(0.5d0*delta_gap + sqrt(delta_gap**2/4d0 + phi**2))
    factor = factor - phi**2/( &
        (0.5d0*delta_gap+sqrt(delta_gap**2/4d0+phi**2))**2 &
        *sqrt(delta_gap**2/4d0+phi**2) &
        )

    factor = factor *lattice_a*t_hop*sin(0.5d0*lattice_a*kt)

    duc_dk = duc_dk*factor

    uc(1) =  xx/sqrt(1d0+xx**2)
    uc(2) = 1d0/sqrt(1d0+xx**2)

    uv(1) = 1d0/sqrt(1d0+yy**2)
    uv(2) =  yy/sqrt(1d0+yy**2)
    
    zham(1,1) = eps_v
    zham(1,2) = -zi*sum(uv*duc_dk)*dkdt
    zham(2,1) = conjg(zham(1,2))
    zham(2,2) = eps_c

    call diag_2x2_complex(zham, zvec, lambda)

    zp_houston_states(:,1) = uv*zvec(1,1) + uc*zvec(2,1)
    zp_houston_states(:,2) = uv*zvec(1,2) + uc*zvec(2,2)

    nex_pol_houston = nex_pol_houston + abs(&
        sum(conjg(zp_houston_states(:,1))*zpsi(:,ik)) &
        )**2

  end do
  nex_pol_houston = nex_pol_houston/nk


end subroutine calc_nex
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
subroutine diag_2x2(mat, vec, lambda)
  implicit none
  real(8),intent(in) :: mat(2,2)
  real(8),intent(out) :: vec(2,2)
  real(8),intent(out) :: lambda(2)
  real(8) :: a, b, c
  real(8) :: ss

  vec = 0d0
  lambda = 0d0

  a = mat(1,1)
  b = mat(1,2)
  c = mat(2,2)
  

  lambda(1) = 0.5d0*((a+c) + sqrt((a-c)**2 + 4d0*b**2)) 
  lambda(2) = 0.5d0*((a+c) - sqrt((a-c)**2 + 4d0*b**2)) 


  if( abs(lambda(1) - a) > abs(lambda(1) - c)  ) then
    vec(2,1) = 1d0
    vec(1,1) = b/(lambda(1)-a)

    vec(1,2) = 1d0
    vec(2,2) = b/(lambda(2)-c)
  else
    vec(1,1) = 1d0
    vec(2,1) = b/(lambda(1)-c)

    vec(2,2) = 1d0
    vec(1,2) = b/(lambda(2)-a)
  end if

  ss = sum(abs(vec(:,1))**2)
  vec(:,1) = vec(:,1)/sqrt(ss)

  ss = sum(abs(vec(:,2))**2)
  vec(:,2) = vec(:,2)/sqrt(ss)

end subroutine diag_2x2
!-------------------------------------------------------------------------------
subroutine diag_2x2_complex(zmat, zvec, lambda)
  implicit none
  complex(8),intent(in) :: zmat(2,2)
  complex(8),intent(out) :: zvec(2,2)
  real(8),intent(out) :: lambda(2)
  real(8) :: a, c
  complex(8) :: zb
  real(8) :: ss

  zvec = 0d0
  lambda = 0d0

  a  = zmat(1,1)
  c  = zmat(2,2)
  zb = zmat(1,2)

  lambda(1) = 0.5d0*((a+c) + sqrt((a-c)**2 + 4d0*abs(zb)**2)) 
  lambda(2) = 0.5d0*((a+c) - sqrt((a-c)**2 + 4d0*abs(zb)**2)) 


  if( abs(lambda(1) - a) > abs(lambda(1) - c)  ) then
    zvec(2,1) = 1d0
    zvec(1,1) = zb/(lambda(1)-a)

    zvec(1,2) = 1d0
    zvec(2,2) = conjg(zb)/(lambda(2)-c)
  else
    zvec(1,1) = 1d0
    zvec(2,1) = conjg(zb)/(lambda(1)-c)

    zvec(2,2) = 1d0
    zvec(1,2) = zb/(lambda(2)-a)
  end if

  
  ss = sum(abs(zvec(:,1))**2)
  zvec(:,1) = zvec(:,1)/sqrt(ss)

  ss = sum(abs(zvec(:,2))**2)
  zvec(:,2) = zvec(:,2)/sqrt(ss)

end subroutine diag_2x2_complex
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
!-------------------------------------------------------------------------------
