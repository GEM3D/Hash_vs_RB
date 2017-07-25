#ifndef _TYPEDEFS_H_
#define _TYPEDEFS_H_

#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <vector>
#include <stack>
#include <iterator>
#include <algorithm> 
#include <assert.h>
#include <unistd.h>
#include <mpi.h>
#include <memory>
#include <time.h>
//#include "parmetis.h"
#include "zoltan.h"
#define Rma 1
#define DEBUG 0


typedef std::vector<int> Vector_Int;

typedef std::vector<unsigned int> Vector_Unint;

typedef std::vector<double>Vector_Dbl;



typedef struct{
	double x;
	double y;
	double z;	
} Center_Coords;

typedef std::vector<Center_Coords>Center_coords;

typedef std::stack<int> Stack_Int;

//typedef std::unique_ptr<int[]> Unique_Ptr_Int

// non local neighbor means non local to the process, it does not reside in the current node
/* this is the structure used to keep track of each 
 element that has a neighbor in another processor,
  the aim here is to not unnecessarily save a processoe tag for each element
*/
#define CPP 1


typedef struct {
  unsigned int face_tag;
  unsigned int proc_id;
  unsigned int elem_id;
} Nonlocal_Nbr;

typedef std::vector<Nonlocal_Nbr>Nonlocal_nbr;

typedef struct {
  unsigned int level;
  //double xyz[6];
  unsigned int centeroid_index;
  // this is the neighbor list local to each processor, set to some minus number if at the boundary
  Vector_Int nbr[6];
   double q;
  // this is the elem_id'th element that resides in proc: proc_id and is a neighbor to current face 
  Nonlocal_nbr nonlocal_nbr;	
  
} cube_data;

typedef std::vector<cube_data>Cube;

typedef struct {
  Vector_Int nbr[6];
}  Cart_proc_con;

typedef struct {
  unsigned int elem_id;
  
  double c1;
  double c2;
  /*
  float c1;
  float c2;
  */ 
} Message_Struct;

typedef std::vector<Message_Struct>Message_struct;

/*
typedef struct 
{

  double P;
  double U;
  double V;
  double W;
	
} unknowns;

typedef std::vector<unknowns>Qs;
*/

#define ANSI_COLOR_RED     "\033[01;31m"
#define ANSI_COLOR_GREEN   "\033[22;32m"
#define ANSI_COLOR_YELLOW  "\033[22;33m"
#define ANSI_COLOR_BLUE    "\033[22;34m"
#define ANSI_COLOR_MAGENTA "\033[22;35m"
#define ANSI_COLOR_CYAN    "\033[22;36m"
#define ANSI_COLOR_RESET   "\033[22;0m"
 
typedef struct { 
    int changes;
    int numGidEntries;
    int numLidEntries;
    int numImport;
    int numExport;
    unsigned int* importGlobalGids;
    unsigned int* importLocalGids;
    unsigned int* exportGlobalGids;
    unsigned int* exportLocalGids; 
    int *importProcs;
    int *importToPart;
    int *exportProcs;
    int *exportToPart;
    int *parts;  
} Zoltan_Out;

typedef struct vec
{
	double nx;
	double ny;
	double nz;
} Vec;

typedef struct normal
{
	Vec face_normal[6];
} Normal;	





void save(Cube& cube,Cube& cube_temp, int taged_cube);
void add_kids(Cube& cube,int id);
void AmrRefineParallel(const int my_rank,const MPI_Comm comm,Cube& cube,const Vector_Unint &ordered_refine_list,const int L,const int M, const int N,Center_coords &XYZ, double ancestor_length[3]);
void sphere(Cube& cube, Vector_Int& refine_list,double C);
void elem_conn_check(Cube& cube);
void circle(double **x, double **y,unsigned int n);
void identify(Cube& cube, Vector_Unint& refine_list,double *x,double *y,unsigned int n);
void check_coordinates(Cube&cube);
void WriteHdf5ParallelSpatialCollection(Cube& cube, int my_rank,unsigned int L,unsigned  int M,unsigned  int N, int npx,int npy,int npz,MPI_Comm comm, MPI_Comm com,MPI_Info info,int my_offset,int ncube_total,Center_coords &XYZ, double ancestor_length[3]);
void XdmfParallelSpatialCollection(Cube &cube,int offset,MPI_Comm comm,MPI_Info info,unsigned int L,unsigned int M,unsigned int N,int ncube_total);
void calculate_processor_offset_for_IO(int my_rank,int mycube,int RMA,int np,int *offset,int *ranks, int* off);
void find_boundary_face(Cube &cube,int cube_id,Vector_Int &bface_id);
void Cartesian_Proc_Conn(const int *Dim,const int *coords,const MPI_Comm Comm_cart, Cart_proc_con *cart_proc_con);
void AmrNewList(const int my_rank,const MPI_Comm comm,Cube& cube,Vector_Unint & refine_list,Vector_Unint &ordered_refine_list);
void WriteOut(int my_rank,const Cube &cube);
void ReadMesh(int my_rank,Cube &cube);
void Broadcast(int my_rank,int mycube,int RMA,int np,int offset,int *ranks,int *off, int* ncube_total);
void ran_list(Cube &cube,Vector_Unint &ranlist,unsigned int level);
//void CallParMetis(const Cube &cube,int np,int offset,int ncube_total ,int my_rank,MPI_Comm Comm,idx_t*part);
void ZoltanGeometricPartitioner(int argcs, char* pArgs[],int my_rank,int np,MPI_Comm Comm,Cube &cube,int ncube_total,int offset,int method,struct Zoltan_Struct *zz,Center_coords &XYZ, double ancestor_length[3],Zoltan_Out *zoltan_out);
void ZoltanGraphPartitioner(int argcs, char* pArgs[],int my_rank,int np,MPI_Comm Comm,const Cube &cube,int ncube_total,int offset,int method,struct Zoltan_Struct *zz,Center_coords &XYZ, double ancestor_length[3],Zoltan_Out *zoltan_out);
void ReadInGeom(double **xyz,int *geom_nn);
void IsInsideSolid(const Cube& cube,const Center_coords &XYZ, Vector_Unint& refine_list,const double *geom_xyz, double *ancestor_length,unsigned int n);
void ReadSTLGeom(int argc, char* argv[],double **triangle_center, int *nn,const double *xyz);
void CheckConnectivityConsistency(Cube &cube);
void CheckWaterTight(int my_rank,Cube &cube,const Center_coords  &XYZ,double *ancestor_length,Normal *normals);
void CalculateNormals(int my_rank,Cube &cube,const Center_coords &XYZ,double *ancestor_length,Normal *normals);
void elem_conn_check(Cube& cube);
void ConstructCRSFormat(int my_rank,Cube &cube,unsigned int **ia,unsigned int **ja, MPI_Comm Comm,int offset, int ncube_total,int np);
void WriteHdf5ParallelPolyvertex(const Cube& cube,const int my_rank,const int npx,const int npy,const int npz,const MPI_Comm comm,const MPI_Info info,const int my_offset,const int ncube_total,const Center_coords &XYZ);
void XdmfParallelPolyvertex(Cube &cube,int offset,MPI_Comm comm,MPI_Info info,int ncube_total);
void WriteHdf5ParallelUnstructured(const Cube& cube,const int my_rank,const int npx,const int npy,const int npz,const MPI_Comm comm,const MPI_Info info,const int my_offset,const int ncube_total,const Center_coords &XYZ,const double ancestor_length[3],Center_coords &XYZ2);


// bitwise operation for power of two
inline void TwoPowN(unsigned int b,double *result) {unsigned int two=1; two=two<<b;	*result=(double)two; };
//inline void TwoPowN(unsigned int b,double *result) {unsigned int two=1; pow(two,b);	*result=(double)two; };


/*
 * \033[22;30m - black
\033[22;31m - red
\033[22;32m - green
\033[22;33m - brown
\033[22;34m - blue
\033[22;35m - magenta
\033[22;36m - cyan
\033[22;37m - gray
\033[01;30m - dark gray
\033[01;31m - light red
\033[01;32m - light green
\033[01;33m - yellow
\033[01;34m - light blue
\033[01;35m - light magenta
\033[01;36m - light cyan
\033[01;37m - white.
*/


#endif
