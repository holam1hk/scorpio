.. Scorpio Documentation master file, created by
   sphinx-quickstart on Sun Jul  3 02:18:39 2022.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to Scorpio's Documentation!
=================================================
**Scorpio** is developed by *SFG* in CUHK Physics to study star formation, planet formation and galaxy evolution. Scorpio uses finite volume method (FVM) for HD/MHD/ADMHD simulations. For the two-fluid (ion-neutral) model, 
The Harten-Lax-van Leer-Contact (HLLC) Riemann solver30 and Harten-Lax-van Leer-Discontinuities (HLLD) Riemann solver31 are combined with linear piecewise polynomial interpolation and Minimod limiter32 to calculate the fluxes of neutral and ion flows.

30. Toro, E. F., Spruce, M., & Speares, W. (1994). Restoration of the contact surface in the HLL-Riemann solver. Shock waves, 4(1), 25-34.
31. Miyoshi, T., & Kusano, K. (2005). A multi-state HLL approximate Riemann solver for ideal magnetohydrodynamics. Journal of Computational Physics, 208(1), 315-344.
32. Bryan, Greg L., et al. Enzo: An adaptive mesh refinement code for astrophysics. The Astrophysical Journal Supplement Series, 2014, 211.2: 19.

.. note::

	This is still under development.

.. toctree::
   :maxdepth: 1
   :caption: Contents:

   hydro
   turbulence
   selfgravity
   time_integration
   code_structure
   parallelization
   
Indices and tables
==================

* :ref:`genindex`
* :ref:`modindex`
* :ref:`search`
