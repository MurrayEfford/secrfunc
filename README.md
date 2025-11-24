# secrfunc

Helper Functions for Package 'secr'

Certain functions for area and transect search in secr use 
numerical integration code from RcppNumerical and RcppEigen that is slow 
to install. This package separates those functions, allowing them to be 
used in secr without slowing down installation.
  
The package is imported and used by secr >= 5.3.2. It will not be used independently.

Please report bugs as Issues on this GitHub page. 
