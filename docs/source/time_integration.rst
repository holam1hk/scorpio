.. _ch:time_integration:

****************
Time Integration
****************



Introduction
============
``gridModule.f90`` calling ::   

    subroutine evolveGridRK2(this)
    subroutine evolveGridRK2sg(this)
    call rk2_1D(this,this%q,this%q1,this%q2)  from rk2.f03
    call rk2_2D(this,this%q,this%q1,this%q2)
    call rk2_3D(this,this%q,this%q1,this%q2)

End

The TR-BDF2 (trapezoidal rule and backward-difference formula of order two) scheme36 is employed to overcome the numerical stiffness of the ion-neutral collision 
source term. The Kernel-based method37-39 is adopted to solve the Poisson equation. The above procedures are coupled into the second-order Runge-Kutta (RK2) 
time-integration method40

36. Tilley, D. A., Balsara, D. S., & Meyer, C. (2012). A numerical scheme and benchmark tests for non-isothermal two-fluid ambipolar diffusion. New Astronomy, 17(3), 368-376.
37. Yen, Chien-Chang, et al. "Self-gravitational force calculation of infinitesimally thin gaseous disks." Journal of Computational Physics 231.24 (2012): 8246-8261.
38. Wang, Hsiang-Hsu, Ronald E. Taam, and David CC Yen. "Self-Gravitational Force Calculation of Infinitesimally Thin Gaseous Disks on Nested Grids." The Astrophysical Journal Supplement Series 224.2 (2016): 16.
39. Wang, Hsiang-Hsu, et al. "Self-gravitational Force Calculation of High-order Accuracy for Infinitesimally Thin Gaseous Disks." The Astrophysical Journal Supplement Series 242.2 (2019): 17.
40. Mignone, A., et al. "PLUTO: a numerical code for computational astrophysics." The Astrophysical Journal Supplement Series 170.1 (2007): 228.

CFL condition
