# Fourier-Laplace-plasma-response-in-driven-Gyrokinetics
Analytic/numeric solution of plasma response to antenna driving in gyro kinetic system. Comparison is made with simulation by AstroGK.

Directory:
python/
  - keeps jupyter notebooks. 
  - figures/
  - gk_solver/
    - gk_apar0 provides functions to obtain dispersion tensor,
      dispersion equation and its derivative.

slow/
  - Theory and derivations of plasma response to an antenna that
    drives slow mode (either by parallel magnetic field or by phase space
    density perturbation).
    
Alfven/
  - Theory and derivations of plasma response to an antenna that drives
    Alfven mode via parallel vector potential (parallel current/perpendicular
    magnetic field fluctuations).
