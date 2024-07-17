# CPM-FEM

This is a version of the CPM-FEM code introduced in Van Oers, Rens et al. P Comput Biol 2014, [doi.org/10.1371/journal.pcbi.1003774](https://doi.org/10.1371/journal.pcbi.1003774) with Qt visualization. The Qt visualization replaces the Matlab visualization in the original code.

## Usage:
First compile:
`qmake`
`make`

Run:
`./cpmfem Network/`

## Further information
Biological, mathematical and implementation details are available in the accompanying paper (see above) and in its supplements.


# Parameters
CPM parameters:
SEED =  23020 / int
NVX = 100 / int 
NVY = 100 / int
MCS = NVX*NVY / int
TARGETVOLUME = 50 / int
CELLFORCE = 1 / double
VOXSIZE = 0.0000025 / double
NRINC = 501 / int
MAXNRITER = 10000 / int
ACCURACY = 0.00001 / double
POISSON = 0.45 / double

Deze parameters kan je aanpassen:
YOUNGS = 12000 / double
GLOBALSTRAIN = true / bool
LOADANGLE = 90 / double
LOAD = 0.3 / double

MOTILITY = 1 / double
INELASTICITY = 500.0 / double
NOSTICKJCM = 500000.0 / double
NOSTICKJCC = 1000000.0 / double
LAMBDADUR = 24 / double

THRESHOLDSTIFF = 15E3 / double
STIFFSENSITIVITY = 0.0005 / double
STIFFENINGSTIFF = 0.1 / double
COMPRESSINGSTIFF = false / bool
LAMBDADISS = 0 / double
COLLAGEN = 10 / double
PIXPERVOX = 5 / int
LINEWIDTH = 2 / int
CELLFORCES = true / bool
DUROTAXIS = true / bool
STRIDE = 10 / int
WSTRIDE = 100 / int
STRAINFIELD = true / bool
FORCEFIELD = true / bool
PRINCFIELD = true / bool

WHICHSIGMF = 1 / int
NODECONNECTION = true / bool
CELLCOLOUR = true / bool
ONECELL = false / bool ---------als je een cell wilt 
TWOCELL = false / bool--------- als je twee cellen wilt 
DISTWOCELLS = 7 / int
CELLCOL = false / bool ---------Deze parameter zet concentratie cellen al dan niet aan     
READCELLS = false / bool
CELLDENSITY = 0 / double
BOUNDARYDIS = 5 / int
NRcf = 0 / int
COLORBAR = true / bool
MAXCOLORBAR = 0 / double

WIDTHCOLORBAR = 100/ int
WRATIOPA = false / bool
WLENGTH = false / bool
WAREA = false / bool
WSQDIS = false / bool
WECC = false / bool
WANGLE = false / bool
WSIGMA = false / bool
WTWOCELLCONTACT = false / bool
WTWOCELLANGLECM = false / bool




