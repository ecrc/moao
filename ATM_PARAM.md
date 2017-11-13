Atmospheric parameters {#ATM_PARAM}
===================================
 
The parameter file is named "prof%d-atmos-night%d.txt", where the %d are integers representing:
* the index of the observation
* the index of the night

These files must contain the following parameters in this strict order.
   
Each parameter must be preceded by a blank line (used for comments).
  
* The number of layers
* r0 at wfs lambda
* profile strength (units as in  E-SPE-ESO-276-0206_atmosphericparameters)
* altitude of layers (meters)
* outer scale (meters)

See an example below:
~~~{.txt}
Nlayer
10
r0 @ wfs lambda
0.128916
cn2 ESO units
16.1481 3.93694 3.15623 2.16198 1.04095 0.800734 0.867462 1.72825 0.411711 0.144132  
h in meters
30 200 390 1880 4500 7500 10500 13500 16500 19500  
l0 in meters
25 25 25 25 25 25 25 25 25 25 
~~~
