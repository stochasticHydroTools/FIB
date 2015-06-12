Example of using FIB and IBAMR to calculate the mobility of an immersed
particle.

First compile in this folder:
  make main2d

Then, run the mobility estimator tool with
  ./main2d input2d.mobility.

This will create a file mobility.dat with the mobility of the immersed
particle evaluated at the position given in points2d.vertex.
