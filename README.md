# FIB
Fluctuating Immersed Boundary Method
by Steven Delong, Florencio Balboa
and Aleksandar Donev (donev@courant.nyu.edu)
Courant Institute of Mathematical Sciences

The IBBrownianBlobHierarchyIntegrator simulates Brownian Dynamics
using the IBAMR library (https://github.com/IBAMR)
and the temporal integrators from the paper (see doc/OverdampedIB.pdf):

"Brownian Dynamics without Green's Functions"
S. Delong, F. Balboa Usabiaga, R. Delgado-Buscalioni, B. E. Griffith and A. Donev
J. Chem. Phys., 140, 134110, 2014.
http://arxiv.org/abs/1401.4198

After compilation, (make main2d, make main3d  as appropriate), one
runs the code with ./main3d <input_file>.  Example input files
are included in subfolders.
