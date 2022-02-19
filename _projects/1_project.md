---
layout: page
title: Crystallisation
description: Simulation of 'Model C' phase transition
img: assets/img/1.jpg
importance: 1
autosize: true
category: fun
---

The transition of a given solid material between regions of its phase diagram is a complicated phenomena to understand. Many attempts have been made over the years, with the first analytical description of transient behaviour of the system's order parameter $$\varphi$$ given by what is called 'Model A'. 


In this model, we consider non-conserving fields, where enthalpy gradients are the primary driving force in the transition across phase space. Consider a free energy functional with a $$\varphi^4$$ potential that allows for two equilibrium phases. We can introduce a symetry breaking term $$\phi(\mathbf{x})$$ that allows preferential occupation of a given state. We can write this free energy functional as:

$$
F[\varphi] = \int d\mathbf{x}\bigg(\frac{\varepsilon_\phi}{2}(\nabla\phi(\mathbf{x})^2)+H\{\underbrace{-\frac{1}{2}\phi^2(\mathbf{x})+\frac{1}{4}\phi^4(\mathbf{x})}_{g(\phi)}\}-\underbrace{B\phi(\mathbf{x})}_{f_b(\phi,B)}      \bigg)
$$

This simple model is essentially a sharp-interface model that can only describe macro- or meso-scale behaviour due to the interfacial coupling of the distinct faces across the interface. We can write the phase field equation in dimensionaless units as:

$$
\tau\frac{\partial\phi}{\partial t}=-\frac{\delta F}{\delta\phi} = -\frac{\varepsilon_\phi^2}{H}\nabla^2\phi - \frac{dg}{d\phi}-\frac{1}{H}\frac{\partial f_b}{\partial\phi}+\underbrace{\eta(\mathbf{x},t)}_\text{thermal noise}
$$

This non-conserved field takes only changes in chemical potential into account, but it does not include thermal affects, an obvious constraint of real crystallisation that must be considered.

We can begin to create more advanced theoretical frameworks. The code outlined here describes Model C, which posits that the changes to the order parameter are highly correlated to changes in the temperature field across various dendrites. It can furthermore include anisotropic tip growth to create fascinating patterns and beautiful patterns.

We can start to formulate a new construction of the phase parameter. In this case, we write the free energy as:

$$
F[\varphi] = \int d\mathbf{x}\bigg(\frac{\varepsilon^2_\phi}{2H}|\nabla\phi(\mathbf{x})^2|+g(\phi)-\lambda\frac{T-T_m}{L/c_p}P(\phi)      \bigg)
$$

where the final term is a similar symmetry breaking condition, related to the heat capacity and latent heat in the sample. We can write the coupled equations as:

$$
\tau\frac{\partial\phi}{\partial t}=-\frac{\delta F(\phi,T)}{\delta\phi}+\zeta(\mathbf{x},t)
$$

$$
\frac{\partial u}{\partial t} = \frac{k}{\rho c_p}\nabla^2u + \frac{1}{2}\frac{\partial P(\phi)}{\partial t}
$$

By seeding the order parameter, we propagate these equations to determine how the crystal structure will develop. By making the assumption that $$\varepsilon_\phi=\varepsilon_\phi(\mathbf{x})$$, we can introduce anisotropy in the rate of dendrite growth. The example presented below includes 4-fold symmetry, where the rate of growth is dependent on the relative angle of the direction normal to the surface to the positive $$x$$-axis.

<div class="row">
    <div class="col-sm mt-3 mt-md-0">
    </div>
    <div class="col-sm mt-3 mt-md-0">
        {% responsive_image path: assets/img/1.jpg title: "example image" class: "img-fluid rounded z-depth-1" %}
    </div>
    <div class="col-sm mt-3 mt-md-0">
    </div>
</div>
<div class="caption">
   The order parameter of a characteristic material, with anisotropic surface tension that creates preferential dendritic growth.
</div>

Owing to the complex computation of coupled partial derivatives given in the equations above, the simulation is computationally expensive. An array-optimised language is best to determine how the heat field will impact the order parameter. For efficiency, a `Fortran` script is included below that creates the data. It is written to include OpenMP optimisation, and it uses the [m_npy](https://raw.githubusercontent.com/MRedies/NPY-for-Fortran/master/src/npy.f90) module to write the Fortran arrays into [NumPy](https://numpy.org) arrays for visualisation.

<div align="center" class="embed-responsive embed-responsive-4by3">
    <video autoplay controls loop class="embed-responsive-item">
        <source src="../../assets/video/test_anim.mp4" type="video/mp4">
    </video>
</div>
<div class="caption">
   The order parameter $\varphi$ as a function of time time, starting with circular seed, including anisotropic tip growth.
</div>
---
Below is the Fortran script used to generate the data. Visualisation is done with `Python`, which created the `mp4` above.
{% raw %}
```fortran

PROGRAM Crystallisation
    USE mpi
    USE omp_lib

    USE m_npy
    IMPLICIT NONE

    INTEGER, PARAMETER :: dp = SELECTED_REAL_KIND(15)
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! DECLARATIONS OF SIM PARAMS
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    REAL(dp) :: h 
    LOGICAL :: onBoundary, file_exists
    INTEGER :: i, j, kk, num_elem, maxk, N ! = 126*5; !N=800;
    INTEGER :: Im1_Jp1, I_Jp1, Ip1_Jp1
    INTEGER :: Im1_J,   I_J,   Ip1_J
    INTEGER :: Im1_Jm1, I_Jm1, Ip1_Jm1
    INTEGER, ALLOCATABLE :: edge(:)
    REAL(dp), ALLOCATABLE :: x(:) !x=(0:h:(N+2)*h)';
    REAL(dp), ALLOCATABLE :: y(:) !y=0:h:(N+2)*h;
    REAL(dp), ALLOCATABLE :: phi(:,:), phi0(:,:), phi_next(:,:), U(:,:), tip_pos(:), velocity(:)
    REAL(dp) :: eps4, U0, lambda !degree of anisotropy, init temp, and slope of growth
    REAL(dp) :: as, eps, a12, tol, a1, a2, D, dt, r, T, noise_coeff
    REAL(dp) :: x0, y0, r0; !seed center and radius
    REAL(dp) :: DERX_r, DERX_l, DERX_t, DERX_b, DERY_t, DERY_r, DERY_l, DERY_b
    REAL(dp) :: MAG2_R, A_R, dA_R, MAG2_L, A_L, dA_L, MAG2_T, A_T, dA_T, MAG2_B, A_B, dA_B
    REAL(dp) :: DERX, DERY, MAG2, A, JR, JL, JT, JB
    REAL(dp) :: g_prime, h_prime, P_prime, noise

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! DECLARATIONS of MPI
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    INTEGER :: ierr, myrank, nprocs, count, remainder, idx, idy, nthread
    INTEGER, ALLOCATABLE :: istart(:), istop(:), numElems(:), recvcounts(:)
    REAL(dp) :: wtime
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! VARIABLES OF Model C
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL omp_set_num_threads(6)
    nthread= omp_get_thread_num()
    WRITE(6,*) "Running on "//trim(int_to_char(nthread))//" OpenMP threads."
    eps4 = 0.06
    lambda = 3.19
    h = 0.4
    N = 126*5
    as = 1-3*eps4; eps = 3*eps4/as; a12 = 4*as*eps; tol=1e-8; 
    a1=0.8836; a2 = 0.6267; D = a2*lambda;
    
    dt = 0.25*h*h/D
    r = D*dt/(h*h)
    T = 3000
    maxk = CEILING(T/dt)
    noise_coeff = 16/1500*0.1

    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Seed and initial params
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    CALL MPI_Init ( ierr )
    call MPI_Comm_rank ( MPI_COMM_WORLD, myrank, ierr )
    call MPI_Comm_size ( MPI_COMM_WORLD, nprocs, ierr )
    wtime = MPI_Wtime( )

    WRITE(6, '(5x, a21)', advance="no") "Allocating arrays ..."
    ALLOCATE(x(1:N+2), y(1:N+2))
    ALLOCATE(phi(N+2, N+2), phi0(N+2, N+2), phi_next(N+2, N+2), U(N+2, N+2))
    x0 = h; y0 = h; r0 = 10*h;
    U0=0.7
    U=-U0
!dir$ loop count min(256)
    DO i=1, N+2
        x(i) = i*h
        y(i) = i*h
        edge(i) = i
    ENDDO
    phi_next = 0.d0
!dir$ loop count min(256)
    DO i=1, N+2
        DO j=1, N+2
            phi0(i,j) = -tanh((dsqrt((x(i)-x0)*(x(i)-x0)+(y(j)-y0)*(y(j)-y0))-r0)/dsqrt(2.0_dp))
        ENDDO
    ENDDO
    WRITE(6, '(5x, a5)') "DONE"
    phi = phi0
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    ! Main time loop
    !%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
    count = (N+2)**2/nprocs
    remainder = MOD((N+2)**2, nprocs)
    
    open (unit=6, carriagecontrol='fortran')
    DO kk=1, maxk
        !main grid loop
        !Boundary conditions (flux absorbing U(N+1)=U(N)+q*dx etc
        !!mirror boundary conditions on phi
        U(1,:) = U(2,:); 
        U(N+2,:) = U(N+1,:); 
        U(:,1) = U(:,2); 
        U(:,N+2) = U(:,N+1);
        U(1,1)=U(2,2); U(1,N+2)=U(2,N+1); U(N+2,1)=U(N+1,2); U(N+2,N+2)=U(N+1,N+1);
        phi(1,:) = phi(2,:); phi(N+2,:) = phi(N+1,:); phi(:,1) = phi(:,2); phi(:,N+2) = phi(:,N+1);
        phi(1,1)=phi(2,2); phi(1,N+2)=phi(2,N+1); phi(N+2,1)=phi(N+1,2); phi(N+2,N+2)=phi(N+1,N+1);

!$omp parallel do schedule(guided) default(shared) private(i, j)
        DO i = 2, N+1
            DO j = 2, N+1
                !define derivative on edges of finite volume in x
                DERX_r = phi(i+1,j)-phi(i,j); DERX_l = phi(i,j)-phi(i-1,j)
                DERX_t = 0.25*(phi(i+1,j+1)-phi(i-1,j+1)+phi(i+1,j)-phi(i-1,j))
                DERX_b = 0.25*(phi(i+1,j)-phi(i-1,j)+phi(i+1,j-1)-phi(i-1,j-1))
                
                !same in y
                DERY_t = phi(i,j+1)-phi(i,j); DERY_b = phi(i,j)-phi(i,j-1)
                DERY_r = 0.25*(phi(i+1,j+1)-phi(i+1,j-1)+phi(i,j+1)-phi(i,j-1))
                DERY_l = 0.25*(phi(i,j+1)-phi(i,j-1)+phi(i-1,j+1)-phi(i-1,j-1))
                
                !anisotropy
                MAG2_R = (DERX_r**2+DERY_r**2)**2;
                IF(MAG2_R.le.tol) THEN
                    A_R = as; dA_R = 0;
                ELSE
                    A_R = as*(1+eps*(DERX_r**4+DERY_r**4)/MAG2_R); dA_R = -a12*DERX_r*DERY_r*(DERX_r**2-DERY_r**2)/MAG2_R;
                ENDIF
                
                MAG2_L = (DERX_l**2+DERY_l**2)**2;
                IF (MAG2_L.le.tol) THEN
                    A_L = as; dA_L = 0;
                ELSE
                    A_L = as*(1+eps*(DERX_l**4+DERY_l**4)/MAG2_L); dA_L = -a12*DERX_l*DERY_l*(DERX_l**2-DERY_l**2)/MAG2_L;
                ENDIF
                
                MAG2_T = (DERX_t**2+DERY_t**2)**2;
                IF (MAG2_T.le.tol) THEN
                    A_T = as; dA_T = 0;
                ELSE
                    A_T = as*(1+eps*(DERX_t**4+DERY_t**4)/MAG2_T); dA_T = -a12*DERX_t*DERY_t*(DERX_t**2-DERY_t**2)/MAG2_T;
                ENDIF     
                
                MAG2_B = (DERX_b**2+DERY_b**2)**2;
                IF (MAG2_B.le.tol) THEN
                    A_B = as; dA_B = 0;
                ELSE
                    A_B = as*(1+eps*(DERX_b**4+DERY_b**4)/MAG2_B); dA_B = -a12*DERX_b*DERY_b*(DERX_b**2-DERY_b**2)/MAG2_B;
                ENDIF
                
                DERX = 0.5*(phi(i+1,j)-phi(i-1,j)); DERY = 0.5*(phi(i,j+1)-phi(i,j-1)); 
                MAG2 = (DERX**2+DERY**2)**2;
                IF (MAG2.le.tol) THEN
                    A = as;
                ELSE
                    A = as*(1+eps*(DERX**4+DERY**4)/MAG2);
                ENDIF
                
                !define fluxes at each side of the finite volume
                JR = A_R*(A_R*DERX_r-dA_R*DERY_r); JL = A_L*(A_L*DERX_l-dA_L*DERY_l);
                JT = A_T*(A_T*DERY_t+dA_T*DERX_t); JB = A_B*(A_B*DERY_b+dA_B*DERX_b);
                !update phi
                g_prime =  -phi(i,j) + phi(i,j)**3
                h_prime = 0.5
                P_prime = (1-phi(i,j)*phi(i,j))*(1-phi(i,j)*phi(i,j))
                
                noise = rnorm()
                phi_next(i,j) = phi(i,j) + dt/(A*A)*( (JR - JL + JT - JB)/(h*h) - g_prime - lambda*U(i,j)*P_prime ) + dt*noise
                !ENDIF !if not on boudnary
            ENDDO !x
        ENDDO !y
!$omp end parallel do

        DO i=2, N+1
            DO j=2, N+1
                U(i,j) = (1-4*r)*U(i,j) +r*( U(i-1,j) + U(i+1,j) + U(i,j-1) + U(i,j+1) ) + h_prime*(phi_next(i,j)-phi(i,j))
            ENDDO
        ENDDO

        phi = phi_next
        IF ((mod(kk, 1000).eq.0).or.(kk.eq.1).or.(kk.eq.maxk)) THEN
            CALL save_npy('data/Nt='//trim(int_to_char(kk))//'.npy', phi)
            CALL save_npy('data/Nt='//trim(int_to_char(kk))//'_U.npy', U)
        ENDIF
        CALL progress(kk,maxk)
    ENDDO !kk
    wtime = MPI_Wtime() - wtime
    CALL MPI_Finalize ( ierr )

    DEALLOCATE(x, y)
    DEALLOCATE(phi, phi0, phi_next, U)

    CONTAINS
        SUBROUTINE rnorm_2D(arr)
        REAL(DP), INTENT(INOUT), ALLOCATABLE :: arr(:,:)
        INTEGER :: i,j
        
        DO i=LBOUND(arr,1), UBOUND(arr,1)
            DO j=LBOUND(arr,2), UBOUND(arr,2)
                arr(i,j) = rnorm()
            ENDDO
        ENDDO
        END SUBROUTINE rnorm_2D
        
        FUNCTION rnorm() RESULT( fn_val )

        !   Generate a random normal deviate using the polar method.
        !   Reference: Marsaglia,G. & Bray,T.A. 'A convenient method for generating
        !              normal variables', Siam Rev., vol.6, 260-264, 1964.

            IMPLICIT NONE
            REAL(DP)  :: fn_val

            ! Local variables

            REAL            :: u, sum
            REAL, SAVE      :: v, sln
            LOGICAL, SAVE   :: second = .FALSE.
            REAL, PARAMETER :: one = 1.0, vsmall = TINY( one )

            IF (second) THEN
            ! If second, use the second random number generated on last call

              second = .false.
              fn_val = v*sln

            ELSE
            ! First call; generate a pair of random normals

              second = .true.
              DO
                CALL RANDOM_NUMBER( u )
                CALL RANDOM_NUMBER( v )
                u = SCALE( u, 1 ) - one
                v = SCALE( v, 1 ) - one
                sum = u*u + v*v + vsmall         ! vsmall added to prevent LOG(zero) / zero
                IF(sum < one) EXIT
              END DO
              sln = SQRT(- SCALE( LOG(sum), 1 ) / sum)
              fn_val = u*sln
            END IF

            RETURN
        END FUNCTION rnorm
        
        FUNCTION IJ2N(i,j,nn) RESULT( fn_val )
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: i,j, nn
            INTEGER :: fn_val
            
            fn_val = i*nn+j
            
        END FUNCTION IJ2N
        
        FUNCTION N2IJ(n, nn) RESULT( fn_val )
            IMPLICIT NONE
            INTEGER, INTENT(IN) :: n, nn
            INTEGER :: fn_val(1:2)
            
            fn_val(1) = n / nn + 1
            fn_val(2) = MOD(n, nn) + 1
            
        END FUNCTION N2IJ
        
        
        SUBROUTINE progress(j, jmax)
          implicit none
          integer(kind=4)::j,k, jmax
          character(len=67)::bar="??????/?????? - ???% |                                                           |"
          write(unit=bar(1:6),fmt="(i6)") j
          write(unit=bar(7:13),fmt="(i6)") jmax
          write(unit=bar(17:20),fmt="(i3)") 10*j/jmax
          do k=1, j/jmax
            bar(6+k:6+k)="*"
          enddo
          ! print the progress bar.
          write(unit=6,fmt="(a1,a1,a17)") '+',char(13), bar
          return
        END SUBROUTINE progress
        
        FUNCTION int_to_char( i )
            !-----------------------------------------------------------------------
            !! Converts an integer number of up to 6 figures into a left-justifed
            !! character variable.
            !
            IMPLICIT NONE
            !
            INTEGER, INTENT(IN) :: i
            CHARACTER (LEN=6)   :: int_to_char
            CHARACTER :: c
            INTEGER   :: n, j, nc
            LOGICAL   :: neg
            !   
            nc = 6
            !
            IF( i < 0 ) then
               nc  = nc - 1
               n   = -i
               neg = .true.
            ELSE
               n   = i
               neg = .false.
            END IF
            !
            j = 1
            DO WHILE( j <= nc ) 
               int_to_char(j:j) = CHAR( MOD( n, 10 ) + ICHAR( '0' ) )
               n = n / 10
               IF( n == 0 ) EXIT
               j = j + 1
            END DO
            !
            IF( j <= nc ) THEN
               DO n = 1, j/2
                  c = int_to_char( n : n )
                  int_to_char( n : n ) = int_to_char( j-n+1 : j-n+1 )
                  int_to_char( j-n+1 : j-n+1 ) = c
               END DO
               IF( j < nc ) int_to_char(j+1:nc) = ' '
            ELSE
               int_to_char(:) = '*'
            END IF
            !
            IF( neg ) THEN
               DO n = nc+1, 2, -1
                  int_to_char(n:n) = int_to_char(n-1:n-1)
               END DO
               int_to_char(1:1) = '-'
            END IF
            !
            RETURN
            !
        END FUNCTION int_to_char

END PROGRAM Crystallisation
```
{% endraw %}
