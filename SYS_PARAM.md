System parameters {#SYS_PARAM}
==============================

 The parameter file "sys-param.txt" must contain the following parameters in this strict order.

 Each parameter must be preceded by a blank line (used for comments).

*  Telescope diameter ( meters)
*  Central obscuration (%)
*  Frame rate (seconds)
*  Number of wave front sensors (including the truth sensor)
*  Number of laser guided stars
*  Number of targets
*  Number of subapertures along the diameter (Nssp)
*  Array of the inverse of the guide star altitude (in 1/meters)
*  Type of the guide star (1:NGS 2:LGS)
*  Direction alphaX of the guide star (arcsec) 
*  Direction alphaY of the guide star (arcsec)
*  Pupil shift of the WFS, (meters)
*  Pupil shift of the WFS, (meters)
*  Rotation of microlenses
*  Rotation of camera
*  Sensitivity coefficient of the WFSs
*  telescope tracking error parameters (arcsec^2)
*  Precision of DPHI precomputation
*  Number of CPU used (for openMP)
*  Magnitude of NGS
*  Photon return at M1 (ph/m2/s)
*  Pixel size for ngs
*  Pixel size for lgs
*  Wave length for NGS (meters)
*  Wave length for LGS (meters)
*  Bandwidth for NGS (meters
*  Transmission for NGS
*  Transmission for LGS
*  Atmosphere transmission
*  Read Out Noise (nb of e-)
*  Lgs_cst
*  Width of the guide stars (arcsec)
*  Altitude of the guide stars (meters)
*  Thickness of the  (meters)
*  Pointing direction of the target(s) in X (arcsec)
*  Pointing direction of the target(s) in Y (arcsec)

See an example below:
~~~{.txt}
tel diam in meters
5
cent obs
0.2
Tframe (s)
0.004
N WFS
7
nlgs
4
nTarget
1
Nssp n subaps per wfs
10
GsAlt
  1.00000000e-05   1.00000000e-05   1.00000000e-05   1.00000000e-05   0.00000000e+00   0.00000000e+00   0.00000000e+00
type
 2.  2.  2.  2.  1.  1.  1.
alphaX in arcsec //the last value is for the truth sensor, overwrite by the function matcov_init_tomo_tiled
 2.87997509e+01   1.76347614e-15  -2.87997509e+01  -5.29042842e-15  -1.77683260e+01   7.03340438e+00   0.00000000e+00
alphaY in arcsec //the last value is for the truth sensor, overwrite by the function matcov_init_tomo_tiled
 0.00000000e+00   2.87997509e+01   3.52695228e-15  -2.87997509e+01   1.64365123e+01   1.61264671e+01   0.00000000e+00
XPup
 0.  0.  0.  0.  0.  0.  0.
YPup
 0.  0.  0.  0.  0.  0.  0.
thetaML
 0.  0.  0.  0.  0.  0.  0.
thetaCam
 0.  0.  0.  0.  0.  0.  0.
sensibilite
 1.  1.  1.  1.  1.  1.  1.
tracking
 1.  1.  1.
pasDphi
0.0001
ncpu
1
magnitude of NGS
 13.  13.  13.
Flux for LGS
7000000.0
pixel size for NGS
0.3
pixel size for LGS
1
wave length for NGS (meters)
6.5e-07
wave length for LGS (meters)
5.89e-07
bandwidth (meters)
3.3e-07
transmission for NGS
0.425
transmission for LGS
0.382
atmosphere transmission
0.84
Read Out Noise (nb of e-)
3
lgs_cst
0.1
spot_width (arcsec)
1
lgs_alt (meters)
100000
lgs_depth (meters)
5000
targetX  (arcsec)
 0.
targetY  (arcsec)
 0.
~~~
