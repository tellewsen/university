module evolution_mod
  use healpix_types
  use params
  use time_mod
  use ode_solver
  use rec_mod
  use spline_2D_mod
  implicit none



  !Use j,k,l as global variable
  integer(i4b) :: j,k,l

  ! Accuracy parameters
  real(dp),     parameter, private :: k_min    = 0.1d0 * H_0 / c
  real(dp),     parameter, private :: k_max    = 1.d3  * H_0 / c
  integer(i4b), parameter          :: n_k      = 100

  integer(i4b), parameter, private :: lmax_int = 6
  real(dp),allocatable, dimension(:) :: dydx

  ! Perturbation quantities
  real(dp), allocatable, dimension(:,:,:) :: Theta
  real(dp), allocatable, dimension(:,:)   :: delta
  real(dp), allocatable, dimension(:,:)   :: delta_b
  real(dp), allocatable, dimension(:,:)   :: Phi
  real(dp), allocatable, dimension(:,:)   :: Psi
  real(dp), allocatable, dimension(:,:)   :: v
  real(dp), allocatable, dimension(:,:)   :: v_b
  real(dp), allocatable, dimension(:,:)   :: dPhi
  real(dp), allocatable, dimension(:,:)   :: dPsi
  real(dp), allocatable, dimension(:,:)   :: dv_b
  real(dp), allocatable, dimension(:,:,:) :: dTheta

  !Precomputing variables
  real(dp), allocatable, dimension(:) :: dtau
  real(dp), allocatable, dimension(:) :: ddtau
  real(dp), allocatable, dimension(:) :: H_p
  real(dp), allocatable, dimension(:) :: dH_p
  real(dp), allocatable, dimension(:) :: eta_precomp

  ! Fourier mode list
  real(dp), allocatable, dimension(:) :: ks

  ! Book-keeping variables
  real(dp),     private :: k_current,ck_current,ckH_p
  integer(i4b), private :: npar = 6+lmax_int

  !Milestone 4 variables
  real(dp),     allocatable,     dimension(:,:,:,:) :: S_coeff
  real(dp), allocatable, dimension(:,:) :: S_lores

  !real(dp), pointer, dimension(:) :: x_high,k_high
  integer(i4b), parameter             :: n_k_highres = 5000
  integer(i4b), parameter             :: n_x_highres = 5000
contains


  ! NB!!! New routine for 4th milestone only; disregard until then!!!
  subroutine get_hires_source_function(x_hires, k_hires, S)
    implicit none

    real(dp), allocatable, dimension(:),   intent(out) :: x_hires, k_hires
    real(dp), allocatable, dimension(:,:), intent(out) :: S

    integer(i4b) :: i,k
    real(dp)     :: g, dg, ddg, dt, tau, ddt, Pi, dPi, ddPi

    ! Task: Output a pre-computed 2D array (over k and x) for the 
    !       source function, S(k,x). Remember to set up (and allocate) output 
    !       k and x arrays too. 

    allocate(x_hires(n_x_highres),k_hires(n_k_highres))

    do i=1,n_x_highres    
        do k=1,n_k_highres
            x_hires(i) = x_init + (x_0-x_init)*(i-1.d0)/(n_x_highres-1.d0)
            k_hires(k) = k_min  + (k_max -k_min)*((k-1.d0)/(n_k_highres-1.d0))**2
        end do 
    end do

    !Test of x and k grid.
    !write(*,*) 'x_hires'
    !write(*,*) x_hires(1),x_hires(n_x_highres)
    !write(*,*) 'k_hires'
    !write(*,*) k_hires(1),k_hires(n_k_highres)
    !write(*,*) 'ks'
    !write(*,*) ks(1),ks(n_k)

    ! Substeps:
    !   1) First compute the source function over the existing k and x
    !      grids
    allocate(S_lores(1:n_t,1:n_k))
    allocate(S_coeff(4,4,n_t,n_k))

    allocate(S(n_x_highres,n_k_highres))

    do k=1,n_k
        k_current = ks(k)
        ck_current= c*k_current

        do i=1,n_t
            g     = get_g(x_t(i))
            dg    = get_dg(x_t(i))
            ddg   = get_ddg(x_t(i))
            tau   = get_tau(x_t(i))
            dt    = dtau(i)
            ddt   = ddtau(i)
            !H_p   = get_H_p(x_t(i))
            !dH_p  = get_dH_p(x_t(i))
            Pi    = Theta(i,2,k)
            dPi   = dTheta(i,2,k)

            ddPi  = 2.d0*ck_current/(5.d0*H_p(i))*(-dH_p(i)/H_p(i)*Theta(i,1,k) + dTheta(i,1,k)) &
                    +0.3d0*(ddt*Pi+dt*dPi) &
                    -3.d0*ck_current/(5.d0*H_p(i))*(-dH_p(i)/H_p(i)*Theta(i,3,k) + dTheta(i,3,k))

            S_lores(i,k) = g*(Theta(i,0,k) +Psi(i,k) + .25d0*Pi) &
                           +exp(-tau)*(dPsi(i,k)-dPhi(i,k)) &
                           -1.d0/ck_current*(H_p(i)*(g*dv_b(i,k) + v_b(i,k)*dg) + g*v_b(i,k)*dH_p(i)) &
                           +.75d0/ck_current**2*((H_0**2/2.d0*((Omega_m+Omega_b)/exp(x_t(i))+4.d0*Omega_r/exp(2.d0*x_t(i)) +4.d0*Omega_lambda*exp(2.d0*x_t(i))))*g*Pi &
                           +3.d0*H_p(i)*dH_p(i)*(dg*Pi+g*dPi)+H_p(i)**2*(ddg*Pi +2.d0*dg*dPi+g*ddPi))
        end do
    end do

    !2) Then spline this function with a 2D spline
    call splie2_full_precomp(x_t, ks, S_lores,S_coeff)

    !3) Finally, resample the source function on a high-resolution uniform
    !      5000 x 5000 grid and return this, together with corresponding
    !      high-resolution k and x arrays
    do k=1,n_k_highres
        do i=1,n_x_highres
            S(i,k) = splin2_full_precomp(x_t, ks, S_coeff, x_hires(i), k_hires(k))
        end do
    end do
  end subroutine get_hires_source_function


  ! Routine for initializing and solving the Boltzmann and Einstein equations
  subroutine initialize_perturbation_eqns
    implicit none
    integer(i4b) :: i

  !Initialize k-grid, ks; quadratic between k_min and k_max
    allocate(ks(n_k))
    do k=1,n_k
        ks(k) = k_min +(k_max -k_min)*((k-1.d0)/(n_k-1.d0))**2
    end do

    !Allocate arrays for perturbation quantities
    allocate(delta(1:n_t, n_k))
    allocate(delta_b(1:n_t, n_k))
    allocate(v(1:n_t, n_k))
    allocate(v_b(1:n_t, n_k))
    allocate(Phi(1:n_t, n_k))
    allocate(Theta(1:n_t, 0:lmax_int, n_k))
    allocate(Psi(1:n_t, n_k))

    allocate(dPhi(1:n_t, n_k))
    allocate(dPsi(1:n_t, n_k))
    allocate(dv_b(1:n_t, n_k))
    allocate(dTheta(1:n_t, 0:lmax_int, n_k))

    !Set these to zero so we only have to fill in those that are non zero later
    Theta(:,:,:) = 0.d0
    dTheta(:,:,:) = 0.d0
    dPhi(:,:) = 0.d0
    dPsi(:,:) = 0.d0

    !Allocate arrays for precomputed variables
    allocate(dtau(n_t),H_p(n_t),dH_p(n_t))
    allocate(ddtau(n_t),eta_precomp(n_t))

    !Precompute useful variables
    do i=1,n_t
       dtau(i)  = get_dtau(x_t(i))
       ddtau(i) = get_ddtau(x_t(i))
       H_p(i)   = get_H_p(x_t(i))
       dH_p(i)  = get_dH_p(x_t(i))
       eta_precomp(i)   = eta_t(i) !was calculated in time_mod
    end do
    !write(*,*) 'dtau(1)'
    !write(*,*) dtau(1)
    !stop


    !write(*,'(*(2X, ES14.6))') H_p(1),dH_p(1),ddtau(1),dtau(1),ks(1)
    !Set up initial conditions for the Boltzmann and Einstein equations
    Phi(1,:)     = 1.d0
    delta(1,:)   = 1.5d0*Phi(1,:)
    delta_b(1,:) = delta(1,:)
    Theta(1,0,:) = 0.5d0*Phi(1,:)
    do k = 1, n_k
        v(1,k)       = c*ks(k)/(2.d0*H_p(1))*Phi(1,k)
        v_b(1,k)     = v(1,k)
        Theta(1,1,k) = -c*ks(k)/(6.d0*H_p(1))*Phi(1,k)
        Theta(1,2,k) = -20.d0*c*ks(k)/(45.d0*H_p(1)*dtau(1))*Theta(1,1,k) !without polarization
        do l = 3, lmax_int
            Theta(1,l,k) = -l/(2.d0*l+1.d0)*c*ks(k)/(H_p(1)*dtau(1))*Theta(1,l-1,k)
        end do
        Psi(1,k)     = -Phi(1,k) - 12.d0*H_0**2/(ks(k)*c*a_t(1))**2*Omega_r*Theta(1,2,k)
    end do
    !write(*,*) 'Theta(1,1,50),Theta(1,2,50),H_p(1),dtau(1),Psi(1,1),Theta(1,3,50),v(1,50)'
    !write(*,*) Theta(1,1,50),Theta(1,2,50),H_p(1),dtau(1),Psi(1,1),Theta(1,3,50),v(1,50)

  end subroutine initialize_perturbation_eqns

  subroutine integrate_perturbation_eqns
    implicit none
    real(dp)     :: x1, x2, x_init
    real(dp)     :: eps, hmin, h1, x_tc, j_tc, dt, t1, t2
    real(dp)     :: R,d_v,d_v_b,q
    real(dp), allocatable, dimension(:) :: y, y_tight_coupling
    logical(lgt)                        :: exist

    eps    = 1.d-8
    hmin   = 0.d0
    h1     = 1.d-5
    allocate(y(npar))
    allocate(dydx(npar))
    allocate(y_tight_coupling(7))
    y_tight_coupling = 0.d0
    y = 0.d0
    dydx = 0.d0

    ! Check if data files exist. If not, compute data.
    inquire(file='precomp_perturb.unf', exist=exist)
    if (exist) then
        open(58,file='precomp_perturb.unf', form='unformatted')
            read(58) Theta
            read(58) delta
            read(58) delta_b
            read(58) Phi
            read(58) Psi
            read(58) v
            read(58) v_b
            read(58) dPhi
            read(58) dPsi
            read(58) dv_b
            read(58) dTheta
        close(58)
    else

        ! Propagate each k-mode independently
        do k = 1, n_k
           write(*,*) 'Current k', k
           k_current = ks(k)  ! Store k_current as a global module variable
           ck_current = c*ks(k) !store c*k

           ! Initialize equation set for tight coupling
           y_tight_coupling(1) = delta(1,k)
           y_tight_coupling(2) = delta_b(1,k)
           y_tight_coupling(3) = v(1,k)
           y_tight_coupling(4) = v_b(1,k)
           y_tight_coupling(5) = Phi(1,k)
           y_tight_coupling(6) = Theta(1,0,k)
           y_tight_coupling(7) = Theta(1,1,k)
       
           ! Find the time to which tight coupling is assumed, 
           ! and integrate equations to that time
           x_tc = get_tight_coupling_time(k_current)

           !write(*,*) 'x_tc=',x_tc
           !write(*,*) 'under x_tc'

           !Integrate from x_init until the end of tight coupling, using
           !the tight coupling equations

           !write(*,*) 'Start of tight coupling'
           !write (*,'(*(2X, ES14.6))') delta(1,k), delta_b(1,k), &
           !v(1,k), v_b(1,k), Phi(1,k), Theta(1,0,k), Theta(1,1,k),Psi(1,k)
           !write (*,'(*(2X, ES14.6))') x_t(1),dv_b(1,k),dPsi(1,k), &
           !dPhi(1,k),dTheta(1,0,k),dTheta(1,1,k),dTheta(1,2,k)

           do j=2,n_t
               if (x_t(j)< x_tc) then 
                   !precompute some variables
                   ckH_p = ck_current/H_p(j)

                   !Solve next step
                   call odeint(y_tight_coupling,x_t(j-1),x_t(j),eps,h1,hmin,derivs_tc, bsstep, output)

                   !Save variables
                   delta(j,k)   = y_tight_coupling(1)
                   delta_b(j,k) = y_tight_coupling(2)
                   v(j,k)       = y_tight_coupling(3)
                   v_b(j,k)     = y_tight_coupling(4)
                   Phi(j,k)     = y_tight_coupling(5)
                   Theta(j,0,k) = y_tight_coupling(6)
                   Theta(j,1,k) = y_tight_coupling(7)
                   Theta(j,2,k) = Theta(1,2,k)
                   do l = 3, lmax_int
                      Theta(j,l,k) = Theta(1,l,k)
                   end do	
                   Psi(j,k)      = -Phi(j,k) - 12.d0*H_0**2.d0/(ck_current*a_t(j))**2.d0*Omega_r*Theta(j,2,k)

                   !Store derivatives that are required for C_l estimation
                   call derivs_tc(x_t(j),y_tight_coupling,dydx)
                   dv_b(j,k)     = dydx(4)
                   dPhi(j,k)     = dydx(5)
                   dTheta(j,0,k) = dydx(6)
                   dTheta(j,1,k) = dydx(7)

                   if(k==40) then
                       write(*,*) 'dPhi(',j,k,') =', dPhi(j,k)
                   endif                   
                   dPsi(j,k)     = -dPhi(j,k) - 12.d0*H_0**2.d0/(ck_current*a_t(j))**2.d0&
                                   *Omega_r*(-2.d0*Theta(j,2,k)+dTheta(j,2,k))


                   !write (*,'(*(2X, ES14.6))') delta(j,k), delta_b(j,k), &
                   !v(j,k), v_b(j,k), Phi(j,k),  Theta(j,0,k), Theta(j,1,k),Psi(j,k)
                   !write (*,'(*(2X, ES14.6))') x_t(j),dPsi(j,k),dPhi(j,k),dv_b(j,k),&
                   !dTheta(j,0,k),dTheta(j,1,k),dTheta(j,2,k)

               else
                   j_tc = j
                   exit
               end if
           end do
           !write(*,*) 'End of tight coupling'

           !Set up variables for integration from the end of tight coupling 
           !until today
           y(1:7) = y_tight_coupling(1:7)
           y(8)   = Theta(j_tc-1,2,k)
           do l = 3, lmax_int
              y(6+l) = Theta(j_tc-1,l,k)
           end do
           !write(*,*) 'Theta(j_tc-1,2,50)'
           !write(*,*) Theta(j_tc-1,2,50)

           !Continue after tight coupling
           !write(*,*) 'start of rec'       

           do j = j_tc, n_t
              
              !Precompute some variables
              ckH_p = ck_current/H_p(j)

              !Integrate equations from tight coupling to today
              !write(*,*) 'running odeint with j =', j
              call odeint(y, x_t(j-1) ,x_t(j), eps, h1, hmin, derivs, bsstep, output)
    
              !Store variables at time step j in global variables
              delta(j,k)   = y(1)
              delta_b(j,k) = y(2)
              v(j,k)       = y(3)
              v_b(j,k)     = y(4)
              Phi(j,k)     = y(5)
              
              do l = 0, lmax_int
                 Theta(j,l,k) = y(6+l)
              end do
              Psi(j,k)     =  - Phi(j,k) - (12.d0*H_0**2.d0)/(ck_current*a_t(j))**2.d0*Omega_r*Theta(j,2,k)
    
              !Store derivatives that are required for C_l estimation
              call derivs(x_t(j),y,dydx)
              dv_b(j,k)     = dydx(4)
              dPhi(j,k)     = dydx(5)
              if(k==40) then
                  write(*,*) 'dPhi(',j,k,') =', dPhi(j,k)
              endif            
              dTheta(j,0,k) = dydx(6)
              dTheta(j,1,k) = dydx(7)                          
              dTheta(j,2,k) = dydx(8)
    
              do l=3,lmax_int-1
                  dTheta(j,l,k) = dydx(6+l)
              end do
              dTheta(j,lmax_int,k) = dydx(6+lmax_int)
    
              dPsi(j,k)     = -dPhi(j,k) - 12.d0*H_0**2.d0/(ck_current*a_t(j))**2.d0*&
                               Omega_r*(-2.d0*Theta(j,2,k)+dTheta(j,2,k)) 
           end do
        end do

        deallocate(y_tight_coupling)
        deallocate(y)
        deallocate(dydx)

        ! Write to file
        open(58,file='precomp_perturb.unf', form='unformatted')
            write(58) Theta
            write(58) delta
            write(58) delta_b
            write(58) Phi
            write(58) Psi
            write(58) v
            write(58) v_b
            write(58) dPhi
            write(58) dPsi
            write(58) dv_b
            write(58) dTheta
        close(58)
    end if
  end subroutine integrate_perturbation_eqns


  subroutine derivs_tc(x,y_tc, dydx)
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y_tc
      real(dp), dimension(:), intent(out) :: dydx

      real(dp) :: d_delta
      real(dp) :: d_delta_b
      real(dp) :: d_v
      real(dp) :: q,R

      real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2
      real(dp) :: Psi,dPhi,dTheta0,dv_b,dTheta1

      delta   = y_tc(1)
      delta_b = y_tc(2)
      v       = y_tc(3)
      v_b     = y_tc(4)
      Phi     = y_tc(5)
      Theta0  = y_tc(6)
      Theta1  = y_tc(7)

      Theta2    = -20.d0*ckH_p/(45.d0*dtau(j))*Theta1

      R         = (4.d0*Omega_r)/(3.d0*Omega_b*a_t(j))

      Psi       = -Phi - 12.d0*(H_0/ck_current/a_t(j))**2.d0*Omega_r*Theta2

      dPhi      = Psi -(ckH_p**2.d0)/3.d0*Phi + H_0**2.d0/(2.d0*H_p(j)**2.d0)*(Omega_m/a_t(j)*delta+Omega_b/a_t(j)*delta_b+4.d0*Omega_r*Theta0/a_t(j)**2.d0)

      dTheta0   = -ckH_p*Theta1 - dPhi

      d_delta   = ckH_p*v   - 3.d0*dPhi

      d_delta_b = ckH_p*v_b - 3.d0*dPhi

      d_v       = -v -ckH_p*Psi

      q         = ( -((1.d0-2.d0*R)*dtau(j) + (1.d0+R)*ddtau(j)) *&
                  (3.d0*Theta1+v_b) - ckH_p*Psi +(1.d0-dH_p(j)/H_p(j))*&
                  ckH_p*(-Theta0 + 2.d0*Theta2) - ckH_p*dTheta0) / &
                  ((1.d0+R)*dtau(j)+dH_p(j)/H_p(j) -1.d0)

      dv_b      = (1.d0/(1.d0+R)) *(-v_b - ckH_p*Psi + &
                  R*(q+ckH_p*(-Theta0 + 2.d0*Theta2)-ckH_p*Psi))

      dTheta1   = (1.d0/3.d0)*(q-dv_b)

      dydx(1) = d_delta
      dydx(2) = d_delta_b
      dydx(3) = d_v
      dydx(4) = dv_b
      dydx(5) = dPhi
      dydx(6) = dTheta0
      dydx(7) = dTheta1

  end subroutine derivs_tc

  subroutine derivs(x,y, dydx) 
      use healpix_types
      implicit none
      real(dp),               intent(in)  :: x
      real(dp), dimension(:), intent(in)  :: y
      real(dp), dimension(:), intent(out) :: dydx

      real(dp) :: d_delta
      real(dp) :: d_delta_b
      real(dp) :: d_v
      real(dp) :: q,R
      integer(i4b) :: l
      real(dp) :: delta,delta_b,v,v_b,Phi,Theta0,Theta1,Theta2,Theta3,Theta4,Theta5,Theta6
      real(dp) :: Psi,dPhi,dTheta0,dv_b,dTheta1,dTheta2

      delta   = y(1)
      delta_b = y(2)
      v       = y(3)
      v_b     = y(4)
      Phi     = y(5)
      Theta0  = y(6)
      Theta1  = y(7)
      Theta2  = y(8)
      Theta3  = y(9)
      Theta4  = y(10)
      Theta5  = y(11)
      Theta6  = y(12)

      R         = (4.d0*Omega_r)/(3.d0*Omega_b*a_t(j))
      Psi       = -Phi - 12.d0*(H_0/ck_current/a_t(j))**2.d0*Omega_r*Theta2

      dPhi      = Psi -(ckH_p**2.d0)/3.d0*Phi + H_0**2.d0/(2.d0*H_p(j)**2.d0)*(Omega_m/a_t(j)*delta+Omega_b/a_t(j)*delta_b+4.d0*Omega_r*Theta0/a_t(j)**2.d0)

      dTheta0   = -ckH_p*Theta1 - dPhi
      d_delta   = ckH_p*v   - 3.d0*dPhi
      d_delta_b = ckH_p*v_b - 3.d0*dPhi
      d_v       = -v -ckH_p*Psi

      dv_b      = -v_b -ckH_p*Psi +dtau(j)*R*(3.d0*Theta1+v_b)
      dTheta1   = ckH_p/3.d0*Theta0 -2.d0/3.d0*ckH_p*Theta2 + &
                  ckH_p/3.d0*Psi +dtau(j)*(Theta1+v_b/3.d0)

      dTheta2   = 2.d0/5.d0*ckH_p*Theta1 - 3.d0/5.d0*ckH_p*Theta3+dtau(j)*0.9d0*Theta2
      do l=3,lmax_int-1
          dydx(6+l) = l/(2.d0*l+1.d0)*ckH_p*y(5+l) - &
                      (l+1.d0)/(2.d0*l+1.d0)*ckH_p*y(7+l) +dtau(j)*y(6+l)
      end do

      dydx(6+lmax_int) = ckH_p*y(6+lmax_int-1) -c*(lmax_int+1.d0)/H_p(j)/eta_t(j)*y(6+lmax_int) +dtau(j)*y(lmax_int)

      dydx(1) = d_delta
      dydx(2) = d_delta_b
      dydx(3) = d_v
      dydx(4) = dv_b
      dydx(5) = dPhi
      dydx(6) = dTheta0
      dydx(7) = dTheta1
      dydx(8) = dTheta2

      !write(*,*)'x,k,dtau'
      !write(*,*) x,k_current,dtau(j)
!      write(*,*)'y'
 !     write(*,*) y
  !    write(*,*)'Psi,Theta2,H_0'
   !   write(*,*) Psi,Theta2,H_0
     ! write(*,*)'dydx'
    !  write(*,*) dydx
!      write(*,*) 'X_e(1),X_e(100),sum(X_e)'
 !     write(*,*) X_e(1),X_e(100),sum(X_e)
  !    write(*,*) 'ckH_p'
   !   write(*,*) ckH_p
    !  stop

  end subroutine derivs

  function get_tight_coupling_time(k)
    implicit none

    real(dp), intent(in)  :: k
    real(dp)              :: get_tight_coupling_time
    integer(i4b)          :: i,n
    real(dp)              :: x
    n =1d4
    do i=0,n
        x = x_init +i*(0.d0-x_init)/n
        if (x < x_start_rec .and. &
            abs(c*k/(get_H_p(x)*get_dtau(x))) <= 0.1d0 .and.& 
            abs(get_dtau(x)) > 10.d0) then 
            get_tight_coupling_time = x
        end if
    end do
  end function get_tight_coupling_time

end module evolution_mod
