 go to the build directory,

  cmake ..  -DCMAKE_INSTALL_PREFIX=../

if cmake gives an error that finding any packages fails, 

1) module show MYMODULE (zoltan, mpi, ...) 
2) open cmake gui via :  ccmake .. 
3) press t (toggle): 	set the corresponding directory (press  enter to edit,after editing press enter again)
4) press c (configure)
5) press q (quit ccmake)
6) cmake ..  -DCMAKE_INSTALL_PREFIX=../ 

  make -j Nprocs VERBOSE=1 (VERBOSE=1 will spit out the header and library paths which is helpful for debugging)
  make install

to run the app:
mpirun -np 4 bin/amrGem  Geom/bunny.stl 



