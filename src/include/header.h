#ifndef _HEADER_H_
#define _HEADER_H_
#include <algorithm>
#include <bitset>
#include <cstdint>
#include <cstdio>
#include <cstdlib>
#include <functional>
#include <iostream>
#include <stack>
#include <unordered_map>
#include <vector>
#include <stdexcept>
#include <map>
#include <unordered_set>
#include <time.h>
#include <set>

using real = double;
using uint = unsigned int;
using integer = int;
using namespace std;

// setting hash=1 will use unordered_map (i.e. hash table) 
// setting hash=0 will use map (i.e. red-black tree)

#define hash 1

// Major data structure used in this research

template <size_t N>
using morton = std::bitset<N>; /*!\var typedef template of bitset class used for Morton encoding */

// trying to override the hash function of gcc
#if ( hash )
template <size_t N>
class morton_hash
{
    public:
    size_t operator()( const morton<N> key ) const // <-- don't forget const, means this function is not allowed to modify the object
    {
        morton<N> kt;
        /*
               for(uint i=0;i<N;i++)
               {
                kt[i]=key[N-i-1];
                }
        */
        //  kt=key;
        size_t hashval = 0;
        hashval = key.to_ulong();
        return hashval;
    }
};

template <size_t N, typename value>
using bitmap = std::unordered_map<morton<N>, value *, morton_hash<N>>; /*! \var typedef template unordered_map used in amrMesh Class*/

#else

template <size_t N>
class compare
{
    public:
    bool operator()( const morton<N> &a, const morton<N> &b ) const
    {
        return a.to_ullong() < b.to_ullong();
    }
};

template <size_t N, typename value>
using bitmap = std::map<morton<N>, value *, compare<N>>; /*! \var typedef template map (red-black tree) used in amrMesh Class*/

#endif
// using bitmapgeom=std::unordered_map<morton , uint >; /*! \var  typedef unordered_map used in amrMesh Class*/

template <size_t N>
using bitvector = std::vector<morton<N>>;

/*
template <size_t N>
using bitunorderedset = std::set<morton<N>,compare<N>>;
*/

// template <size_t N>
// using bitunorderedset = std::unordered_set<morton<N>>;

/**    \class Tree
 * \brief  This Class Generates a 4:1 balancerd AMR mesh
 *
 *
 *
 */

template <size_t N, typename value>
class Hdf5Xmf;

template <size_t N, typename value>
class Tree
{
    template <size_t N1, typename value2>
    friend class Hdf5Xmf; /*!< this is a friend class to write out in hdf5 format*/

    protected:
    real ancestorlength[3];    /*!< original length of the first generation (root) element*/
    real ancestorcoords[3];    /*!< centeroid of the  of the first generation (root) element*/
    morton<N> ancestorkey = 0; /*!< root value is always set as 0000000000000   */
    uint npx;                  /*!< discritization in x direction*/
    uint npy;                  /*!< discretization in y direction */
    uint npz;                  /*!< discretization in y direction */

    private:
    bitmap<N, value> mesh;                          /*!< base main container for hashmap */
    bitvector<N> refinelist;                        /*!<  list of elements tagged to be refined  */
                                                    //    bitunorderedset<N> derefinelist; /*!< list of elements tagged to be removed*/
    std::unordered_map<morton<N>, int> refinelist1; /*!< list of elements tagged to be removed*/
    bitvector<N> mortonSTL;
    std::unordered_map<morton<N>, int> derefinelist; /*!< list of elements tagged to be removed*/

    public:
    Tree( real *length, real *coords, uint nx, uint ny, uint nz ); /*!< constructor*/
    void level( morton<N> key, uint *level );                      /*!< obtains the level of the element from morton key */
    void centroid( morton<N> key, real *xyz );   /*!< calculates the centroid of the cube given the morton key of the element*/
    void enclosingBox( morton<N> key, real *X ); /*!< calculates the range that an element occupies in 3D space for a given Element*/
    typename bitmap<N, value>::iterator begin(); /*!< iterator returning the first object*/
    typename bitmap<N, value>::iterator end();   /*!< iterator returning the last object*/
    uint size();                                 /*!< returns the size of the mesh*/
    //    void reserve( uint *reservedsize );          /*!<this function reserves the memory given the reservedsize of the mesh */
    void siblings( morton<N> key, uint mylevel, morton<N> *sibkey ); /*!< extracts the siblings from morton code */
    void refine( morton<N> key );                                    /*!< perfomrs refinement for a tagged element given the Morton Key*/
    void derefine( morton<N> key );                                  /*! performs derefinement on a single element given a morton key*/
    void refineRefineList();                                         /*!<performs the refinement */
    void fourToOne(); /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
    void findFlipLevel( morton<N> key, uint *mylevel, uint *changedirectionlevel,
                        uint *direction ); /*!< detects the flip level, this info used in finding nonlocal neighbors*/
    void flipForNbr( morton<N> *key, uint *mylevel, uint *changedirectionlevel,
                     uint *direction );      /*!< perform the actual operation to identify the nonlocal nbr */
    uint IsInVectorList( morton<N> key );    /*!< checks to see if a given code is ialready in the list */
    void addToList( morton<N> key );         /*!<adds element to refinelist */
    uint count( morton<N> key );             /*!< counts the number of elements*/
    void addToDerefineList( morton<N> key ); /*!< adds the element to derefinelist */
    void derefineDerefineList();             /*!< Derefines the mesh*/
                                             // geometry
    uint isInsideSolid( const morton<N> key, const real *geom_xyz,
                        uint n ); /*!< tags the elements if any points of the gemoetry resides in the enclosing box*/
    typename bitmap<N, value>::iterator find( morton<N> key ); /*!< this function is to find a value given the key*/
    void convertStl2Morton( uint geom_size, real *geom_xyz );  /*!< converts all the stl points to morton code*/
    void convertCoordToMorton( real *xyz, morton<N> &key );    /*!<converts coordinates of a point to morton code */
    void pushToRefinelist( uint level );
    void printMesh();
    void refineRefineList2();
    void refineRefineList3();
    void fourToOneSet();
    void pushToRefinelistSet( uint nlevel );
    void refineRefineListSet();
    typename unordered_map<morton<N>,int>::iterator Dbegin();
    typename unordered_map<morton<N>,int>::iterator Dend();
    void retainFourToOne();
    void play( uint maxlevel, uint geom_nn, real *geom_xyz, double xx[3]);
    void pushToDerefinelist( uint nlevel );
    void removeFromDerefineList(typename  std::unordered_map<morton<N>,int>::iterator it  );

    // void derefineDerefineList();
    //    float load_factor();
    //
    ~Tree(); /*!< Destructor of the class, it frees the memeory pointed by pointer in the hashmap value  if allocated*/
};

/**    \class Hdf5Xmf
 *  * \brief  This Writes out Tree data in hdf5 format with *.xmf as metadata suitable for paraview and visit
 *
 */

template <size_t N, typename value>
class Hdf5Xmf
{
    public:
    // Hdf5Xmf;

    void xdmfPolyvertexSerial( Tree<N, value> &T, uint appx );
    void writeHdf5PolyvertexSerial( Tree<N, value> &T, uint appx );
    void writeP( Tree<N, value> &T, uint iterate ); /*!< write only the centroids of the cubes */

    void writeHdf5MultiBlockSerial( Tree<N, value> &T, uint appx ); /*!<writes out xdmf metadata required to read the hdf5 file */
    void xdmfMultiBlockSerial( Tree<N, value> &T, uint appx );      /*!<writes out the mesh in hdf5  */
    void write( Tree<N, value> &T, uint iterate );                  /*!< write all of the the elements */

    void writeHdf5MultiBlockHighLevel( Tree<N, value> &T, uint appx );
    void xdmfMultiBlockHighLevel( Tree<N, value> &T, uint appx );
    void writeH( Tree<N, value> &T, uint appx ); /*!< write only the elements with highest level*/
                                                 //~Hdf5Xmf
};
// the reason I am keeping these separate is because of the need for the solution vector in the base class
// template<size_t N>
// class Hdf5XmfV;
/**    \class Voxel
 * \brief  This Class Generates an unbalancerd Voxel to improve search by geometry partitioning
 *
*/
template <size_t N, typename value>
class Voxel : public Tree<N, value>
{
    private:
    uint maxlevel; /*!< maximum level of refinement*/
    uint numMax;   /*!< number of elements having the highest level*/
    bitmap<N, value> mesh;
    bitvector<N> lookup;

    public:
    Voxel( real *length, real *coords ) : Tree<N, value>( length, coords, 2, 2, 2 ) {}; /*!< constructor */
    void setLevel( uint *l );                                                           /*!< sets the maximum level for refinement*/
    void generateSearchTree( real *geom_xyz, uint n );                                  /*!< generates an initial tree*/
    void distributeGeomToLeaves( real *geom_xyz, uint n ); /*!< distributes geometry to different cells (leaves) */
    uint checkSiblingStatus( morton<N> key, morton<N> *sibkey );
    ;                        /*!<checks to see if siblings include any points and whether they have the same level */
    void derefineGeomTree(); /*!< derefines the tree based on geometry*/
    bool IsInsideSegment( morton<N> key, real *xyz );
    template <size_t N1, typename value1>
    friend class Hdf5XmfV; /*!< this is a friend class to write out in hdf5 format*/
    ~Voxel();              /*!< Destructor of this class*/
    // add erase because regular refine does not free the memory, though I reallocate in the next step, normally it should not leak
    // void  erase(morton key); /* rather specialized version of erase, since the memory has to be freed as well*/
};

/**    \class Hdf5XmfV
 *   *
 *  * \brief  This class writes out Voxel data in hdf5 format with *.xmf as metadata suitable for paraview and visit,
 *   the reason for keeping these separate is because Tree data will have a solution vector soon and
 *    that class will be much more different down the road
 *   */

template <size_t N, typename value>
class Hdf5XmfV
{
    public:
    // Hdf5Xmf;

    void xdmfPolyvertexSerial( Voxel<N, value> &V );
    void writeHdf5PolyvertexSerial( Voxel<N, value> &V );
    void writeP( Voxel<N, value> &V );                    /*!< write only the centroids of the cubes */
    void writeHdf5MultiBlockSerial( Voxel<N, value> &V ); /*!<writes out all the elements */
    void xdmfMultiBlockSerial( Voxel<N, value> &V );
    void write( Voxel<N, value> &V ); /*!<writes out only the elements with highest level in hdf5  */
};

void readSTLGeom( int argc, char *argv[], real **triangle_center, int *nn, const real *xyz );
inline void TwoPowN( uint b, real *result )
{
    uint two = 1;
    two = two << b;
    *result = (real)two;
};

#define RED "\033[01;31m"
#define GREEN "\033[22;32m"
#define YELLOW "\033[22;33m"
#define BLUE "\033[22;34m"
#define MAGENTA "\033[22;35m"
#define CYAN "\033[22;36m"
#define RESET "\033[22;0m"

#endif
