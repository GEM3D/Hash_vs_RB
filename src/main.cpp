#include "header.h"


/**
 * \main
 * @file main.cpp
 * @brief      Main program file for Serial AMR.
 *             "Serial implementation of Morton Encoding for AMR using unordered_map as wel as map classes of C++ STL library"              
 *             
 *             Part of NSF project
 *             
 *
 *             Required Libraries:
 *              
 *              * HDF5 
 *              * CMAKE
 *             Usage:
 *              progName input/geometry>.stl <run.txt
 *
 * @author     Jaber Hasbestan, Inanc Senocak
 * @date       25 July 2017
 *
 * Copyright (c) 2017
 * Mechanical and Bio-medical Engineering Department
 * Boise State University
 *
 * Permission is hereby granted, free of charge, to any person obtaining a copy
 * of this software and associated documentation files (the "Software"), to deal
 * in the Software without restriction, including without limitation the rights
 * to use, copy, modify, merge, publish, distribute, sublicense, and/or sell
 * copies of the Software, and to permit persons to whom the Software is
 * furnished to do so, subject to the following conditions:
 *
 * The above copyright notice and this permission notice shall be included in all
 * copies or substantial portions of the Software.
 *
 * THE SOFTWARE IS PROVIDED "AS IS", WITHOUT WARRANTY OF ANY KIND, EXPRESS OR
 * IMPLIED, INCLUDING BUT NOT LIMITED TO THE WARRANTIES OF MERCHANTABILITY,
 * FITNESS FOR A PARTICULAR PURPOSE AND NONINFRINGEMENT. IN NO EVENT SHALL THE
 * AUTHORS OR COPYRIGHT HOLDERS BE LIABLE FOR ANY CLAIM, DAMAGES OR OTHER
 * LIABILITY, WHETHER IN AN ACTION OF CONTRACT, TORT OR OTHERWISE, ARISING FROM,
 * OUT OF OR IN CONNECTION WITH THE SOFTWARE OR THE USE OR OTHER DEALINGS IN THE
 * SOFTWARE.
 */


int main(int argcs, char* pArgs[])
{


morton<64> A;

cout<<A<<endl;


double a=3;

double *c=&a;

bitmap<64,double> mesh;  /*!< \typedef Basic Object the Will be used for unordered_map */

bitmap<16,int>test;

mesh.insert({A,c});

int b=2;
int *d=&b;

morton<16>B;
B.flip(15);
cout<<B<<endl;
test.insert({B,d});

auto it=mesh.begin();
//cout<<it->first<<endl;
//cout<<it->second<<endl;


auto it2=test.begin();
//cout<it2->first<<endl;
//cout<<it2->second<<endl;

cout<<A.size()<<endl;
cout<< B.size()<<endl;

real ancestorlength[3] = {2.0f, 2.0f, 2.0f};

real ancestorcoords[3] = {0.0f, 0.0f, 0.0f};

real xyz[3];

uint npx = 2, npy = 2, npz = 2;

Tree<64, real> myTree( ancestorlength, ancestorcoords, npx, npy, npz );
Tree<64, uint> myTree1( ancestorlength, ancestorcoords, npx, npy, npz );
Tree<32, uint> myTree2( ancestorlength, ancestorcoords, npx, npy, npz );

cout << "size of my tree " << myTree2.size() << endl;

/*
for(auto it=myTree.begin();it!=myTree.end();it++)
{
cout<<"hello \n"<<endl;
}
*/

uint gestimate=100;

//myTree1.reserve(&gestimate);

// Hdf5Xmf<32,int> C();

bitvector<32> mylist;

uint co;

real *geom_xyz = NULL;
int geom_nn;

real xyz1[6] = {-1.0f, 1.0f, -1.0f, 1.0f, -1.0f, 1.0f};
// cout<<xyz1[0]<<endl;

readSTLGeom( argcs, pArgs, &geom_xyz, &geom_nn, xyz1 );

myTree1.refine( 0 );

uint mylevel;

Hdf5Xmf<64, uint> IO;

#if(0)

for ( uint j = 0; j < 7; j++ )
{
    co = 0;

    for ( auto local_it = myTree1.begin(); local_it != myTree1.end(); local_it++ )
    {
        myTree1.level( local_it->first, &mylevel );
        // cout<<local_it->first<<endl;
        if ( myTree1.isInsideSolid( local_it->first, geom_xyz, geom_nn ) == 1 )
        {
            // cout<<"yes" <<endl;
            myTree1.addToList( local_it->first );
        }
    }

    myTree1.refineRefineList();
    cout << "size=" << myTree1.size() << endl;

 //   IO.write( myTree1, j );
}
#endif



myTree1.convertStl2Morton( geom_nn, geom_xyz );

cout<<myTree1.size()<<endl;

uint meshlevel;

cout<<BLUE<<"enter the refinement level for tree" <<RESET<<endl;

cin>>meshlevel;

delete [] geom_xyz;

#if(1)

clock_t t0=clock();

for ( uint j = 0; j < meshlevel; j++ )
{
   
    myTree1.pushToRefinelistSet(j+1);
    myTree1.refineRefineListSet();
   cout<<j<<" : "<<myTree1.size()<<endl; 
 //   std::cout << "load_factor = " << myTree1.load_factor() << std::endl;
  // IO.write( myTree1, j );
}

clock_t t1=clock();

cout<<myTree1.size()<<" "<<float(t1-t0)/CLOCKS_PER_SEC<<endl;

#endif


#if(0)
for ( uint j = 0; j < meshlevel; j++ )
{

    myTree1.pushToRefinelistSet(j+1);
//    myTree1.refineRefineList2();

    //myTree1.refineRefineList3();
   myTree1.refineRefineListSet();
   cout<<j<<" : "<<myTree1.size()<<endl;
    //std::cout << "load_factor = " << myTree1.load_factor() << std::endl;
     // IO.write( myTree1, j );
     //myTree1.printMesh();
}   

IO.write( myTree1, 0 );
// move the geometry 
//

double xx[3]={0,0,0};
xx[0]=-0.5;

for(uint i=0;i<geom_nn;i++)
{
geom_xyz[3*i+0]+=xx[0];
geom_xyz[3*i+1]+=xx[1];
geom_xyz[3*i+2]+=xx[2];

}




myTree1.convertStl2Morton(geom_nn, geom_xyz);
myTree1.pushToDerefinelist(meshlevel+1);
myTree1.derefineDerefineList();
myTree1.pushToDerefinelist(meshlevel);
myTree1.derefineDerefineList();
cout<<"After derefinement"<<" : "<<myTree1.size()<<endl;



myTree1.convertStl2Morton(geom_nn, geom_xyz);
for ( uint j = 0; j < meshlevel; j++ )
{

    myTree1.pushToRefinelistSet(j+1);
//    myTree1.refineRefineList2();

    //myTree1.refineRefineList3();
   myTree1.refineRefineListSet();
   cout<<j<<" : "<<myTree1.size()<<endl;
    //std::cout << "load_factor = " << myTree1.load_factor() << std::endl;
     // IO.write( myTree1, j );
     //myTree1.printMesh();
}   

#endif   


double xx[3]={0,0,0};


xx[0]=xx[0]+0.1;

xx[1]=xx[0]+0.0;

xx[2]=xx[0]+0.0;

//myTree1.play(meshlevel,geom_nn, geom_xyz,xx);

//IO.write( myTree1,0);

//  IO.writeP( myTree1, 0 );
//cout<<myTree1.size()<<endl;


//myTree1.printMesh();

/*
morton<32>key;
real xx[3]={0.1,0,0};
myTree2.convertCoordToMorton( xx, key ); 
//cout<<key<<endl;

morton<5>mod;

for(uint i=0;i<5;i++)
{
mod[4-i]=key[32-3*i-1];
}




cout<<mod<<endl;
xx[0]=xx[0]+.1;

myTree2.convertCoordToMorton( xx, key ); 

//cout<<key<<endl;


for(uint i=0;i<5;i++)
{
mod[4-i]=key[32-3*i-1];
}
cout<<mod<<endl;
/*
while(1)
{

}
*/



//cout<<myTree1.size()<<endl;

// myTree1.IO;
//cout << sizeof( real ) << endl;
//cout << sizeof( double ) << endl;

//Voxel<32, real> myVoxel( ancestorlength, ancestorcoords );

/*
uint vlevel=5;

myVoxel.setLevel(&vlevel);

myVoxel.generateSearchTree(geom_xyz,geom_nn);

myVoxel.distributeGeomToLeaves(geom_xyz,geom_nn);

Hdf5XmfV<32,real> IOV;

IOV.write(myVoxel);
*/
}

