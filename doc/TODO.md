
1) 16 Jan 2017: 
	
	Cmake CMakeList.txt is currently developed for the linux systems
        need to make some adjustments to make it suitable for Mac and Windows

2) 17 Jan 2017:
	Add vesrion control such that if the version of the external package is less than the one recommended, generate a warning
 
3) 18 Jan 2017:
	FindZoltan.cmake has a hard time finding the zoltan package, I hard coded the ZOLTAN_DIR for now. all other packages are found easily
	I notices that the difference between the zoltan and MPI and Parmetis and .. is the lack of env. variable $PATH
	need to fix this before the final release
