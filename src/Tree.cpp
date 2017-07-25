#include "header.h"

//====================================================================================
//
//                  Constructor
//
//====================================================================================

template <size_t N, typename value>
Tree<N, value>::Tree( real *length, real *coords, uint nx, uint ny, uint nz )
{
    for ( uint i = 0; i < 3; i++ )
    {
        ancestorlength[i] = length[i];
        ancestorcoords[i] = coords[i];
    }
    npx = nx;
    npy = ny;
    npz = nz;

    mesh.insert( {ancestorkey, nullptr} );
   // mesh.max_load_factor(0.8);
    // cout<<ancestorkey<<endl;
}

template <size_t N, typename value>
Tree<N, value>::~Tree()
{
    for ( auto it = begin(); it != end(); it++ )
    {
        if ( it->second != nullptr )
        {
            delete[] it->second;
        }
    }
}

//====================================================================================
//
//                  Calculates the Level of the ELement from Morton Code
//
//====================================================================================

template <size_t N, typename value>
void Tree<N, value>::level( morton<N> key, uint *level )
{
    *level = N / 3;
    uint rem = N % 3;
    // cout<<"rem="<<rem<<endl; /*!< This remaining will affect the number of bits left over, this offset will be used in the following
    // loop*/
    uint iend = N / 3;
    morton<N> keytemp;

    /*! to prevent unnecesary bit operation, the morton code is placed from starting from  left hand side*/

    for ( uint i = 0; i < iend; i++ )
    {
        if ( key[3 * i + rem] == false && key[3 * i + rem + 1] == false && key[3 * i + rem + 2] == false )
        {
            /*! now look and see if any siblings exist*/
            //
            // cout<<"i=:"<<i <<endl;
            // cout<<key[3*i+rem]<<<<key[3*i+rem]<<key[3*i+rem]<<endl;
            keytemp = key;
            if ( mesh.count( keytemp.flip( 3 * i + rem ) ) == 0 )
            {
                // cout<<"count="<<mesh.count(keytemp)<<endl;

                *level = *level - 1;

                // cout<<"here"<<endl;
            }
        }
        else
        {
            break;
        }
    }
}

//=====================================================================
//
//
// This routine calculated coordinates based on morton code
//  Tis is bsically a series of the form (-1)^boolean*2^-i
//
//=====================================================================

template <size_t N, typename value>
void Tree<N, value>::centroid( morton<N> key, real *xyz )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    // uint k;

    uint mylevel;
    level( key, &mylevel );

    // cout<<"key="<<key<<endl;
    // cout<<"level"<<level<<endl;
    for ( uint j = 0; j < 3; j++ )
    {
        dx = 0.0;
        for ( uint i = 0; i < mylevel; i++ )
        {
            if ( key[N - 3 * i - j - 1] == false )
            {
                sign = -1.0;
            }
            else
            {
                sign = 1.0;
            }

            i++;
            TwoPowN( i, &result );
            // cout<<"result:"<<result<<endl;
            dx = dx + sign * 1. / result;
            i--;
        }
        xyz[j] = ancestorcoords[j] + dx * ( ancestorlength[j] ) * 0.5;
    }
}

//====================================================================
//
//                             Enclosing Box
//
//====================================================================
template <size_t N, typename value>
void Tree<N, value>::enclosingBox( morton<N> key, real *X )
{
    real result = 0.0, sign = 0.0, dx = 0.0;
    real xyz[3];
    real dy, dz;

    uint mylevel;

    level( key, &mylevel );

    // cout<<"mylevel is "<<level<<endl;
    // get coords

    for ( uint j = 0; j < 3; j++ )
    {
        dx = 0.0;
        for ( uint i = 0; i < mylevel; i++ )
        {
            if ( key[N - 3 * i - j - 1] == false )
            {
                sign = -1.0;
            }
            else
            {
                sign = 1.0;
            }

            i++;
            TwoPowN( i, &result );
            // cout<<"result:"<<result<<endl;
            dx = dx + sign * 1. / result;
            i--;
        }
        xyz[j] = ancestorcoords[j] + dx * ( ancestorlength[j] ) * 0.5;
    }

    real idenum;
    TwoPowN( mylevel, &idenum );

    // cout<<"idenum"<<idenum<<endl;
    idenum = 1. / idenum;
    dx = ancestorlength[0] * idenum * 0.5;
    X[0] = xyz[0] - dx;
    X[1] = xyz[0] + dx;

    // denum=pow(2,level);
    dy = ancestorlength[1] * idenum * 0.5;
    X[2] = xyz[1] - dy;
    X[3] = xyz[1] + dy;

    // denum=pow(2,level);
    dz = ancestorlength[2] * idenum * 0.5;
    X[4] = xyz[2] - dz;
    X[5] = xyz[2] + dz;
}

//======================================================================
//
//                  these two functions return theiterator object for
//                      begining and end of the unordered_map
//
//======================================================================

template <size_t N, typename value>
typename bitmap<N, value>::iterator Tree<N, value>::begin()
{
    return ( mesh.begin() );
}

template <size_t N, typename value>
typename bitmap<N, value>::iterator Tree<N, value>::end()
{
    return ( mesh.end() );
}

//======================================================================
//
//    return the Size of mesh and reserve the amount you want,
//          reserve is important to prevent rehashing
//
//======================================================================

template <size_t N, typename value>
uint Tree<N, value>::size()
{
    return ( mesh.size() );
}

/* \brief this initial guess is important to prevent from rehashing*/
/*
template <size_t N, typename value>
void Tree<N, value>::reserve( uint *reservedsize )
{
    mesh.reserve( *reservedsize );
}
*/
//=========================================================
//
//                       Get Siblings
//
//==========================================================

template <size_t N, typename value>
void Tree<N, value>::siblings( morton<N> key, uint mylevel, morton<N> *sibkey )
{
    morton<N> kt = key;

    mylevel = mylevel - 1;

    sibkey[0] = kt.flip( N - 3 * ( mylevel ) - 1 );
    // cout<<"sibkeyi[0] "<<sibkey[0]<<endl;

    kt = key;
    sibkey[1] = kt.flip( N - 3 * ( mylevel ) - 2 );
    // cout<<"sibkeyi[1] "<<sibkey[1]<<endl;

    kt = key;
    sibkey[2] = kt.flip( N - 3 * ( mylevel ) - 3 );
    // cout<<"sibkeyi[2] "<<sibkey[2]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 1 );
    sibkey[3] = kt.flip( N - 3 * ( mylevel ) - 2 );
    // cout<<"sibkeyi[3] "<<sibkey[3]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 1 );
    sibkey[4] = kt.flip( N - 3 * ( mylevel ) - 3 );
    // cout<<"sibkeyi[4] "<<sibkey[4]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 3 );
    sibkey[5] = kt.flip( N - 3 * ( mylevel ) - 2 );
    // cout<<"sibkeyi[3] "<<sibkey[5]<<endl;

    kt = key;
    kt.flip( N - 3 * ( mylevel ) - 1 );
    kt.flip( N - 3 * ( mylevel ) - 2 );

    sibkey[6] = kt.flip( N - 3 * ( mylevel ) - 3 );
}

//==================================================================
//
//                  Refines a single elmenet
//
//==================================================================

template <size_t N, typename value>
void Tree<N, value>::refine( morton<N> key )
{
    uint mylevel;
    level( key, &mylevel );
    morton<N> temp;
    temp = key;
    /* \brief if the morton code does not exist in mesh, refinement is not permitted (refining a nonexsiting element not permitted)*/
    /* \brief for now I have assigned the values as NULL, assign it to what you like for solving your particular problem*/
    if ( mesh.count( key ) == 0 )
    {
        cout << RED "fatal error: key does not exist" RESET << endl;
        cout << key << endl;
        exit( 0 );
    }

    // cout<<"N= "<<N<<endl;
    // cout<<"mylevel "<<mylevel<<endl;
    // cout<<"key "<<key<<endl;
    // cout<<"index "<<N-3*mylevel-1<<endl;

    // too many type conversions but the asnwer will always be positive so should not be a problem

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 2 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    temp.flip( N - 3 * mylevel - 2 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 2 );
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );

    temp = key;
    temp.flip( N - 3 * mylevel - 1 );
    temp.flip( N - 3 * mylevel - 2 );
    temp.flip( N - 3 * mylevel - 3 );
    mesh.insert( {temp, nullptr} );
}

//===================================================================
// \brief
//@return perfomes 4:1 consistent refinement for the elements
// in the refinelist, this list is empty after return
//
//
//===================================================================

template <size_t N, typename value>
void Tree<N, value>::refineRefineList()
{
    morton<N> key;
    fourToOne();

    while ( refinelist.size() != 0 )
    {
        key = refinelist.at( refinelist.size() - 1 );
        refine( key );
        refinelist.pop_back();
    }
}

template <size_t N, typename value>
void Tree<N, value>::refineRefineListSet()
{
    morton<N> key;
    fourToOneSet();

for(auto it=refinelist1.begin();it!=refinelist1.end();it++)
{
        key = it->first;
        refine( key );
}
refinelist1.clear();

}
//=============================================================
//
//            Check and see if that is inside solid
//
//============================================================
template <size_t N, typename value>
uint Tree<N, value>::isInsideSolid( const morton<N> key, const real *geom_xyz, uint n )
{
    uint j;
    uint a, b, c;
    real xyz[6];
    uint bol;

    bol = 0;
    // get bounding box

    enclosingBox( key, xyz );

    for ( j = 0; j < n; j++ )
    {
        a = geom_xyz[3 * j] > xyz[0] && geom_xyz[3 * j] < xyz[1];
        b = geom_xyz[3 * j + 1] > xyz[2] && geom_xyz[3 * j + 1] < xyz[3];
        c = geom_xyz[3 * j + 2] > xyz[4] && geom_xyz[3 * j + 2] < xyz[5];

        if ( a && b && c )
        {
            bol = 1;
            break;
        }
    }

    return ( bol );
}

// ================================================================
//
//                       Four to One balance
//
// ================================================================

template <size_t N, typename value>
uint Tree<N, value>::IsInVectorList( morton<N> key )
{
    uint bol = 0;

    for ( uint i = 0; i < refinelist.size(); i++ )
    {
        if ( key == refinelist.at( i ) )
        {
            bol = 1;
            break;
        }
    }
    return ( bol );
}

//====================================================================
//
/*!< all we are interested is the nonlocal neighbors, i.e. the neighbors of the parents as siblings will have same level*/

//====================================================================

template <size_t N, typename value>
void Tree<N, value>::fourToOne() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey;
    uint istart = 0;
    uint iend = refinelist.size();

    a = 1;
    // the so-called ripple effect

    while ( a == 1 )
    {
        a = 0;

        //      printf("istart %d iend %d\n",istart,iend);

        for ( uint i = istart; i < iend; i++ )
        {
            mykey = refinelist.at( i );
            level( mykey, &mylevel );

            if ( mylevel > 1 )
            {
                for ( uint j = 0; j < 3; j++ )
                {
                    //      uint j=1;

                    kt = mykey;

                    //      cout<<"mykey\t"<<kt<<endl;
                    findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                    // if the change in signe does not happen that is a boundary cube
                    //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                    if ( changedirectionlevel != 0 )
                    {
                        flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                        mesh.count( kt );

                        // if this element exists, the level of nbr>=level of tagged element

                        //      cout<<"kt="<<mesh.count(kt)<<endl;
                        if ( mesh.count( kt ) == 0 )
                        {
                            kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                        }
                        //      cout<<"modified key\t"<<kt<<endl;
                        level( kt, &nbrlevel );

                        // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                       // if its is not in the list add to list
                     //if ( mylevel > nbrlevel && !IsInVectorList( kt ) )
                      // use std find is faster
                        if ( mylevel > nbrlevel && std::find(refinelist.begin(),refinelist.end(), kt)==refinelist.end()  )
                         {
                            refinelist.push_back( kt );
                            a = 1;
                        }
                    }
                }
            }
        }

        if ( a == 1 )
        {
            istart = iend;
            iend = refinelist.size();
            //  cout<<"istart=\t"<<istart<<"iend\t"<<iend<<endl;
        }
    }

    /*!< this approach eliminates search algorithm  as now we do not have the restrictions on cutting the cube that we had in the previous
     * approach*/
    // std::sort (refine_list.begin(), refine_list.end(), compare_level);
    // cout<<"============================================================"<<endl;
    // cout<<           "EXISTING 4:1 BALANCE CHECK" <<endl;
    // cout<<"============================================================"<<endl;
}

//====================================================================================

template <size_t N, typename value>
void Tree<N, value>::fourToOneSet() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey;

    a = 1;
    // the so-called ripple effect

    while ( a == 1 )
    {
        a = 0;

        //      printf("istart %d iend %d\n",istart,iend);

        for ( auto it =refinelist1.begin() ; it != refinelist1.end(); it++ )
        {
            mykey = it->first;
            level( mykey, &mylevel );

            if ( mylevel > 1 /*&& it->second!=1*/ )
            {
                for ( uint j = 0; j < 3; j++ )
                {
                    //      uint j=1;

                    kt = mykey;

                    //      cout<<"mykey\t"<<kt<<endl;
                    findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                    // if the change in signe does not happen that is a boundary cube
                    //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                    if ( changedirectionlevel != 0 )
                    {
                        flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                        mesh.count( kt );

                        // if this element exists, the level of nbr>=level of tagged element

                        //      cout<<"kt="<<mesh.count(kt)<<endl;
                        if ( mesh.count( kt ) == 0 )
                        {
                            kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                        }
                        //      cout<<"modified key\t"<<kt<<endl;
                        level( kt, &nbrlevel );

                        // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                       // if its is not in the list add to list
                     //if ( mylevel > nbrlevel && !IsInVectorList( kt ) )
                      // use std find is faster
                       
                       //edit 

                        if ( mylevel > nbrlevel  && it->second==0 )
                         {
                            refinelist1.insert( {kt,0} );
                            //  
                            a = 1;
                       }
                      
/*                 
      else
                       {
                          refine(mykey);
                          refinelist1.erase(mykey);                          
                        }
                       */  
                    }
                }
             it->second=1;
            // refine(mykey);
            // refinelist1.erase(mykey); 
             
            }
        }
/*
        if ( a == 1 )
        {
            istart = iend;
            iend = refinelist.size();
            //  cout<<"istart=\t"<<istart<<"iend\t"<<iend<<endl;
        }
*/
//cout<<"here"<<refinelist1.size()<<endl;
    }

    /*!< this approach eliminates search algorithm  as now we do not have the restrictions on cutting the cube that we had in the previous
     * approach*/
    // std::sort (refine_list.begin(), refine_list.end(), compare_level);
    // cout<<"============================================================"<<endl;
    // cout<<           "EXISTING 4:1 BALANCE CHECK" <<endl;
    // cout<<"============================================================"<<endl;
}




//=================================================================================================
//=================================================================================================
// gets the first change of direction
//
template <size_t N, typename value>
void Tree<N, value>::flipForNbr( morton<N> *key, uint *mylevel, uint *changedirectionlevel, uint *direction )
{
    // cout<<"mylevel, changelevel and direction = "<<*mylevel<<"\t"<<*changedirectionlevel<<"\t"<<*direction <<endl;

    for ( uint i = ( *changedirectionlevel ); i <= ( *mylevel ); i++ )
    {
        // cout<<(bit-3*(i-1)-(*direction)-1)<<endl;
        ( *key ).flip( N - 3 * ( i - 1 ) - ( *direction ) - 1 );
        // cout<<"in here"<<(*key)<<endl;
    }

    // cout<<"exiting flip NBR ..."<<endl;
}

template <size_t N, typename value>
void Tree<N, value>::addToList( morton<N> key )
{
    refinelist.push_back( key );

    // cout<<"refinesize"<<refinelist.size()<<endl;
}

//
//   Given the Key find the level to be traversed to find the nbr
//   \brief
//   a) there is not change in sign in a given direction for the entire code, this corresponds to boundary element
//   b) otherwise, find the level at which the change in sign appears
//   c) siblings normally exist
//   direction gets the values 1, 2,3 for x,y, azd z directions
//=============================================================================================

template <size_t N, typename value>
void Tree<N, value>::findFlipLevel( morton<N> key, uint *mylevel, uint *changedirectionlevel, uint *direction )
{
    bool bol;

    bol = key[N - 3 * ( ( *mylevel ) - 1 ) - ( *direction ) - 1];

    // cout<<"index is = "<<bit-3*((*mylevel)-1)-(*direction)-1<<endl;
    // cout<<"bol = "<<bol<<endl;

    *changedirectionlevel = 0;

    for ( uint i = ( *mylevel ) - 1; i > 0; i-- )
    {
        // cout<<"i="<<i <<key[bit-3*(i-1)-(*direction)-1]<<endl;

        if ( key[N - 3 * ( i - 1 ) - ( *direction ) - 1] != bol )
        {
            *changedirectionlevel = i;
            break;
        }
    }

    // cout<<*changedirectionlevel<<endl;

    if ( *changedirectionlevel == 0 )
    {
        // cout<<"element is a boundary element = "<<endl;
    }

    // cout<<"exiting flip level ..."<<endl;
}
//=================================================
template <size_t N, typename value>
uint Tree<N, value>::count( morton<N> key )
{
    return ( mesh.count( key ) );
}


#if(0)

//======================================================================================
/**
*  \brief If any of the siblings are listed in the dereffinement do not add to the list
*  as derefining one child means removing all the siblings
*
*
* */

template <size_t N, typename value>
void Tree<N, value>::addToDerefineList( morton<N> key )
{
    bool bol = true;

    uint mylevel, klevel;
    morton<N> sibkey[7], kt;
    level( key, &mylevel );
    siblings( key, mylevel, sibkey );

    for ( uint i = 0; i < 7; i++ )
    {
        for ( uint j = 0; j < derefinelist.size(); j++ )
        {
            kt = derefinelist.at( j );
            // level(kt,&klevel);
            if ( sibkey[i] == kt && mylevel == klevel )
            {
                bol = false;
                break;
            }
        }
    }

    if ( bol )
    {
        derefinelist.push_back( key );
    }
}





template <size_t N, typename value>
void Tree<N, value>::derefineDerefineList()
{
    morton<N> key;
    bool bol;
    bitvector<N> temp;
    morton<N> sibkey[8];

    // std::sort(derefinelist.begin(),derefinelist.end(), compare);

    uint mylevel;
    uint maxlevel = 1;

    // get max level

    for ( uint i = 0; i < derefinelist.size(); i++ )
    {
        key = derefinelist.at( i );
        level( key, &mylevel );
        if ( mylevel > maxlevel )
        {
            maxlevel = mylevel;
        }
    }

    for ( uint j = maxlevel; j > 0; j-- )
    {
        for ( uint i = 0; i < derefinelist.size(); i++ )
        {
            key = derefinelist.at( i );
            level( key, &mylevel );

            if ( mylevel == j )
            {
                temp.push_back( key );
            }
        }
    }
    /*
    for(uint i=0;i<temp.size();i++)
    {
    cout<<"temp"<<temp.at(i)<<endl;
    }
    */

    for ( uint i = 0; i < temp.size(); i++ )
    {
        key = temp.at( i );
        level( key, &mylevel );
        cout << "temp level= " << mylevel << endl;
    }

    if ( temp.size() != derefinelist.size() )
    {
        cout << RED "missing element after sorting" RESET << endl;
    }

    uint changedirectionlevel, nbrlevel;

    morton<N> kt;

    for ( uint i = 0; i < temp.size(); i++ )
    // while(temp.size()!=0)
    {
        key = temp.at( i );
        //
        //
        // key=temp.front();
        // cout<<"tagged elem"<<key<<endl;

        bol = true;

        sibkey[7] = key;

        level( key, &mylevel );

        siblings( key, mylevel, sibkey );

        /*
        for(uint l=0;l<8;l++)
        {
        cout<<sibkey[l]<<endl;
        }
        */
        for ( uint j = 0; j < 8; j++ )
        {
            /*
            cout<<"==========================================\n"<<endl;
            cout<<"key is " <<key<<endl;
            cout<<"level is " <<mylevel<<endl;
            */
            for ( uint k = 0; k < 3; k++ )
            {
                key = sibkey[j];
                findFlipLevel( key, &mylevel, &changedirectionlevel, &k );
                // printf("k= %d\n",k);
                // cout<<"chanege dir "<<changedirectionlevel<<endl;
                if ( changedirectionlevel != 0 )
                {
                    flipForNbr( &key, &mylevel, &changedirectionlevel, &k );
                    mesh.count( key );
                    // cout<<"nbr key "<<key<<endl;
                    // if this element exists, the level of nbr>=level of tagged element

                    //      cout<<"kt="<<mesh.count(kt)<<endl;
                    if ( mesh.count( key ) == 0 )
                    {
                        key[N - 3 * ( mylevel - 1 ) - 1] = 0;
                        key[N - 3 * ( mylevel - 1 ) - 2] = 0;
                        key[N - 3 * ( mylevel - 1 ) - 3] = 0;
                    }
                    //      cout<<"modified key\t"<<kt<<endl;
                    level( key, &nbrlevel );

                    // cout<<"nbr level =  "<<nbrlevel<<" my level = "<<mylevel<<endl;

                    // if its is not in the list add to list
                    if ( mylevel < nbrlevel )
                    {
                        bol = false;
                        break;
                    }
                }
            }
            // cout<<"==========================================\n"<<endl;
        }

        // cout<<"bol= "<<bol<<endl;

        if ( bol == true )
        {
            key = sibkey[7];
            derefine( key );
        }

        // cout<<"*************************************************************"<<endl;
        // printKey();
        // cout<<"*************************************************************"<<endl;
        // cout<<RED "Mesh Size" RESET<<mesh.size()<<endl;
    }
    derefinelist.clear();
    derefinelist.shrink_to_fit();
}
#endif


template <size_t N, typename value>
typename bitmap<N, value>::iterator Tree<N, value>::find( morton<N> key ) /**!< this function is to find a value given the key*/
{
    typename bitmap<N, value>::iterator temp;

    temp = mesh.find( key );
    /*
    if(temp==end())
    {
    cout<<"key not found"<<endl;
    exit(0);
    }
    else
    {
    return(temp);
    }
    */
    return ( temp );
}



template <size_t N, typename value>
void Tree<N, value>::convertStl2Morton( uint geom_size, real *geom_xyz ) /**!< this function is to find a value given the key*/
{
    real x, y, z, xmid, ymid, zmid;
    morton<N> key = 0;

       real xyz[3];

   mortonSTL.clear();

    // std::cout<<"xmax "<<xmax<<" ymax "<<ymax<<" zmax "<<zmax<<endl;

    // std::cout<<"xmin "<<xmin<<" ymin "<<ymin<<" zmin "<<zmin<<endl;
   
    for ( uint i = 0; i < geom_size; i++ )
    {

        xyz[0] = geom_xyz[3 * i];
        xyz[1] = geom_xyz[3 * i + 1];
        xyz[2] = geom_xyz[3 * i + 2];
        convertCoordToMorton( xyz, key );

       mortonSTL.push_back( key );
    }

// sort this in the future 
//std::sort(mortonSTL.begin(),mortonSTL.end());

#if(0)
    real X[6];

    bool bol1, bol2, bol3;

    for ( uint i = 0; i < mortonSTL.size(); i++ )
    {

        x = geom_xyz[3 * i];
        y = geom_xyz[3 * i + 1];
        z = geom_xyz[3 * i + 2];

        key = mortonSTL.at( i );
        enclosingBox( key, X );

        bol1 = X[0] <= x && x <= X[1];
        bol2 = X[2] <= y && y <= X[3];
        bol3 = X[4] <= z && z <= X[5];

       // cout<<bol1<< " " <<bol2 << " "<<bol3 <<endl;
  //      cout<< x << " " <<y<< " " <<z <<endl;
    //    cout<<"  "<<X[0]<<" "<<X[1]<<" " <<X[2]<< " "<<X[3]<<" "<<X[4]<<" "<<X[5]<<endl;
//	cout<< "key "<<key<<endl;
        if ( !( bol1 && bol2 && bol3 ) )
        {
            throw ::runtime_error( RED "error in generating morton code for STL geometry" RESET );
        }
    }

    std::cout << GREEN << "Morton Code Construction for STL Successful" << RESET << endl;

  // std::sort(mortonSTL.begin(), mortonSTL.end(), [](const morton<N> & lhs, const morton<N> & rhs) { return lhs.to_string() < rhs.to_string(); });

#endif
    /*
    std::cout<<"x "<<x<<"y "<<y<<"z "<<z<<endl;
    std::cout<<key<<endl;
    std::cout<<" xmin "<<X[0]<<" xmax "<<X[1]<<endl;

    std::cout<<" ymin "<<X[2]<<" zmax "<<X[3]<<endl;
    std::cout<<" zmin "<<X[4]<<" zmax "<<X[5]<<endl;
    */
}

template <size_t N, typename value>
void Tree<N, value>::convertCoordToMorton( real *xyz, morton<N> &key ) /**!< this function is to find a value given the key*/
{
    real x, y, z, xmid, ymid, zmid;

    static const real half = 0.5;

    real xmin = ancestorcoords[0] - half * ancestorlength[0];
    real xmax = ancestorcoords[0] + half * ancestorlength[0];

    real ymin = ancestorcoords[1] - half * ancestorlength[1];
    real ymax = ancestorcoords[1] + half * ancestorlength[1];

    real zmin = ancestorcoords[2] - half * ancestorlength[2];
    real zmax = ancestorcoords[2] + half * ancestorlength[2];



   

        x = xyz[0];
        y = xyz[1];
        z = xyz[ 2];
 key=0;


        for ( uint j = 0; j < N / 3; j++ )
        {
            xmid = half * ( xmin + xmax );

            ymid = half * ( ymin + ymax );

            zmid = half * ( zmin + zmax );

            if ( x > half * ( xmin + xmax ) )
            {
                key.flip( N - 1 - 3 * j );
                xmin = xmid;
            }
            else
            {
                xmax = xmid;
            }

            if ( y > half * ( ymin + ymax ) )
            {
                key.flip( N - 1 - 3 * j - 1 );
                ymin = ymid;
            }
            else
            {
                ymax = ymid;
            }

            if ( z > half * ( zmin + zmax ) )
            {
                key.flip( N - 1 - 3 * j - 2 );
                zmin = zmid;
            }
            else
            {
                zmax = zmid;
            }
        }
}






template <size_t N, typename value>
void             Tree<N, value>::printMesh()
{

//bitmap<N,typename value*>::hasher fn = mesh.hash_function();

    for ( auto it = mesh.begin(); it != mesh.end(); it++ )
    {
        cout <<it->first<<"   " <<it->first.to_ulong() << endl;
        // cout<<fn(it->first)<<endl; 
    }
/*
for (auto& x: mesh) {
    std::cout << "Element [" << x.first << ":" << x.second << "]";
    std::cout << " is in bucket #" <<(((double) mesh.bucket (x.first))/mesh.bucket_count())<<"  " <<(double)x.first.to_ulong() << std::endl;
  }
*/
}


/*
template <size_t N, typename value>
float             Tree<N, value>::load_factor()
{

return(mesh.load_factor());

}
*/




template <size_t N, typename value>
void Tree<N, value>::pushToRefinelistSet( uint nlevel ) /**!< this function is to find a value given the key*/
{

    morton<N> key, key1;

    key = 0;
    // cout<<mortonSTL.size()<<endl;

    for ( uint i = 0; i < mortonSTL.size(); i++ )
    {

        key = mortonSTL.at( i );

        for ( uint j = 0; j < 3 * nlevel; j++ )
        {
            key1[N - j - 1] = key[N - j - 1];
        }
        // cout<<key1<<endl;
        if ( mesh.count(key1)!=0 )
        {
            refinelist1.insert({key1,0});
        }
    }
    // cout<<refinelist.size()<<endl;
    /*
    for(uint i=0;i<refinelist.size();i++)
    {
    std::cout<<refinelist.at(i)<<endl;
    }
    */
}




template <size_t N, typename value>
void Tree<N, value>::pushToRefinelist( uint nlevel ) /**!< this function is to find a value given the key*/
{

    morton<N> key, key1;

    key = 0;
    // cout<<mortonSTL.size()<<endl;

    for ( uint i = 0; i < mortonSTL.size(); i++ )
    {

        key = mortonSTL.at( i );

        for ( uint j = 0; j < 3 * nlevel; j++ )
        {
            key1[N - j - 1] = key[N - j - 1];
        }
        // cout<<key1<<endl;
        if ( std::find( refinelist.begin(), refinelist.end(), key1 ) == refinelist.end() )
        {
            addToList( key1 );
        }
    }
    // cout<<refinelist.size()<<endl;
    /*
    for(uint i=0;i<refinelist.size();i++)
    {
    std::cout<<refinelist.at(i)<<endl;
    }
    */
}

template <size_t N, typename value>
void Tree<N, value>::refineRefineList2() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey;
    uint istart = 0;
    uint iend = refinelist.size();

    a = 1;
    // the so-called ripple effect

    while ( a == 1 )
    {
        a = 0;

        //      printf("istart %d iend %d\n",istart,iend);

        for ( uint i = istart; i < iend; i++ )
        {
            mykey = refinelist.at( i );
            level( mykey, &mylevel );

            if ( mylevel > 1 )
            {
                for ( uint j = 0; j < 3; j++ )
                {
                    //      uint j=1;

                    kt = mykey;
                    //      cout<<"mykey\t"<<kt<<endl;
                    findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                    // if the change in signe does not happen that is a boundary cube
                    //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                    if ( changedirectionlevel != 0 )
                    {
                        flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                        mesh.count( kt );

                        // if this element exists, the level of nbr>=level of tagged element

                        //      cout<<"kt="<<mesh.count(kt)<<endl;
                        if ( mesh.count( kt ) == 0 )
                        {
                            kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                        }
                        //      cout<<"modified key\t"<<kt<<endl;

		         level( kt, &nbrlevel );

                        // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                       // if its is not in the list add to list
                     //if ( mylevel > nbrlevel && !IsInVectorList( kt ) )
                      // use std find is faster 
                        if ( mylevel > nbrlevel && std::find(refinelist.begin()+i,refinelist.end(), kt)==refinelist.end()  )
                         {  
                            refinelist.push_back( kt );
                            a = 1;
                        }
                    }
                }
            }
           refine(mykey);
        }
       
       
        if ( a == 1 )
        {
            istart = iend;
            iend = refinelist.size();
            //  cout<<"istart=\t"<<istart<<"iend\t"<<iend<<endl;
        }
    }
    
   refinelist.clear();
}

template <size_t N, typename value>
void Tree<N, value>::refineRefineList3() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey;
    uint istart = 0;
    uint iend = refinelist.size();

    a = 1;
    // the so-called ripple effect

    while ( a == 1 )
    {
        a = 0;

        //      printf("istart %d iend %d\n",istart,iend);

        for ( uint i = istart; i < iend; i++ )
        {
            mykey = refinelist.at( i );
            level( mykey, &mylevel );

            if ( mylevel > 1 )
            {
                for ( uint j = 0; j < 3; j++ )
                {
                    //      uint j=1;

                    kt = mykey;
                    //      cout<<"mykey\t"<<kt<<endl;
                    findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                    // if the change in signe does not happen that is a boundary cube
                    //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                    if ( changedirectionlevel != 0 )
                    {
                        flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                        mesh.count( kt );

                        // if this element exists, the level of nbr>=level of tagged element

                        //      cout<<"kt="<<mesh.count(kt)<<endl;
                        if ( mesh.count( kt ) == 0 )
                        {
                            kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                            kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                        }
                        //      cout<<"modified key\t"<<kt<<endl;

		         level( kt, &nbrlevel );

                        // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                       // if its is not in the list add to list
                     //if ( mylevel > nbrlevel && !IsInVectorList( kt ) )
                      // use std find is faster 
                        if ( mylevel > nbrlevel && std::find(refinelist.begin()+i,refinelist.end(), kt)==refinelist.end()  )
                         {  
                            refinelist.push_back( kt );
                            a = 1;
                        }
                    }
                }
            }
           refine(mykey);
        }
       
        refinelist.erase(refinelist.begin(), refinelist.begin()+iend);              
       
        if ( a == 1 )
        {
//            istart = iend;
            istart=0;             
            iend = refinelist.size();
            //  cout<<"istart=\t"<<istart<<"iend\t"<<iend<<endl;
        }
    }
    
   refinelist.clear();
}

//*************************************************************************
//
//                            Derefinement Routines
//
//
//*************************************************************************
// derefinelist is private I need to access it by forest

template <size_t N, typename value>
typename unordered_map<morton<N>,int>::iterator Tree<N, value>::Dbegin()
{
    return ( derefinelist.begin() );
}

template <size_t N, typename value>
typename unordered_map<morton<N>,int>::iterator Tree<N, value>::Dend()
{
    return ( derefinelist.end() );
}

//======================================================================================
/**
*  \brief If any of the siblings are listed in the dereffinement do not add to the list
*  as derefining one child means removing all the siblings
*
*
* */
#if ( 1 )
template <size_t N, typename value>
void Tree<N, value>::addToDerefineList( morton<N> key )
{
    bool bol = true;

    uint mylevel, klevel;
    morton<N> sibkey[7], kt;
    level( key, &mylevel ); 
    siblings( key, mylevel, sibkey );

    // if any of the siblings has a higher level we can not remove that element
    //
    for ( uint i = 0; i < 7; i++ )
    {
        level( sibkey[i], &klevel );
        if ( klevel > mylevel )
        {
            bol = false;
            break;
        }
    }

    if ( bol && derefinelist.count( key ) == 0 )
    {
        derefinelist.insert( {key,0} );

        for ( uint i = 0; i < 7; i++ )
        {
            derefinelist.insert( {sibkey[i],0} );
        }
    }
}
#endif


#if ( 1 )
template <size_t N, typename value>
void Tree<N, value>::pushToDerefinelist( uint nlevel ) /**!< this function is to find a value given the key*/
{
    morton<N> key, key1, kt;

    key = 0;
    // cout<<mortonSTL.size()<<endl;
    uint mylevel;
    morton<N> sibkey[7];
    std::unordered_set<morton<N>> temp;

    // temp is to speed up the search, mortonSTL is the new geometry encoded based on the new location (moving body)
    /*
    save a temporary set of the mortonSTL (temp) with one level lower than the actual nlevel
    if that element is not inside mortonSTL, insert it to temp

    */
    // look at the parent element, thats why I do (nlevel-1)
    //
    for ( uint i = 0; i < mortonSTL.size(); i++ )
    {
        key = mortonSTL.at( i );

        for ( uint j = 0; j < 3 * ( nlevel - 1 ); j++ )
        {
            key1[N - j - 1] = key[N - j - 1];
        }
        if ( temp.count( key1 ) == 0 )
        {
            temp.insert( key1 );
        }
    }

    for ( auto it = mesh.begin(); it != mesh.end(); it++ )
    {

        key = it->first;
        level( key, &mylevel );
        kt = key;

        kt[N - 3 * ( mylevel - 2 ) - 1] = 0;
        kt[N - 3 * ( mylevel - 2 ) - 2] = 0;
        kt[N - 3 * ( mylevel - 2 ) - 3] = 0;

        if ( mylevel == nlevel && temp.count( kt ) == 0 )
        {
            // insert this element and all its siblings

            addToDerefineList( key );
        }
    }
// cout<<derefinelist.size()<<endl;

}
#endif

template <size_t N, typename value>
void Tree<N, value>::retainFourToOne() /*!< imposes 4:1 balance given the list of elments to be refined in the vector refine list*/
{
    morton<N> kt;
    uint mylevel, nbrlevel, a, changedirectionlevel;
    morton<N> mykey, sibkey[7];
    //    uint      istart = 0;
    //    uint      iend   = refinelist.size();
    //      printf("istart %d iend %d\n",istart,iend);

    //
    // check the balance for nonlocal neighbors
    //
    for ( auto it = derefinelist.begin(); it != derefinelist.end(); it++ )
    {
        mykey = ( it->first );
        level( mykey, &mylevel );

        if ( mylevel > 1 )
        {
            for ( uint j = 0; j < 3; j++ )
            {
                //      uint j=1;

                kt = mykey;

                //      cout<<"mykey\t"<<kt<<endl;
                findFlipLevel( kt, &mylevel, &changedirectionlevel, &j );

                // if the change in signe does not happen that is a boundary cube
                //       cout<<"changedirectionlevel"<<changedirectionlevel<<endl;
                if ( changedirectionlevel != 0 )
                {
                    flipForNbr( &kt, &mylevel, &changedirectionlevel, &j );
                    mesh.count( kt );

                    // if this element exists, the level of nbr>=level of tagged element

                    //      cout<<"kt="<<mesh.count(kt)<<endl;
                    if ( mesh.count( kt ) == 0 )
                    {
                        kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
                        kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
                        kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
                    }
                    //      cout<<"modified key\t"<<kt<<endl;
                    level( kt, &nbrlevel );

                    // cout<<"nbr level"<<nbrlevel<<"my level"<<mylevel<<endl;

                    // if its level is lower than the neighbor,remove it from the derefine list
                    if ( mylevel < nbrlevel )
                    {
                        // remove the element and all its siblings from the list
                        cout<<mylevel<<" "<<nbrlevel<<endl;
                        removeFromDerefineList( it );
                    }
                }
            }
        }
    }

 cout<<derefinelist.size()<<endl;
}

template <size_t N, typename value>
void Tree<N, value>::removeFromDerefineList(typename std::unordered_map<morton<N>,int>::iterator it ) /*!< E;iminate from derefinelist*/
{
    // this is to remove from the list so removes all 8 element
    uint mylevel;
    morton<N> key,sibkey[8];
    key=it->first; 
    level( key, &mylevel );
    siblings( key, mylevel, sibkey );
    sibkey[7]=key;
 // find the key and its siblings and put int to 1

    for(uint i=0;i<8;i++)
    {
    auto it=derefinelist.find (sibkey[i]);
    if(derefinelist.count(sibkey[i])!=0)
    {
    it->second=1;
    }
    }
   
        
}

//===========================================================
//
//         derefine: Eliminate An Element From Mesh
//
//===========================================================
template <size_t N, typename value>
void Tree<N, value>::derefine( morton<N> key )
{
    /*! \brief if the morton code does not exist in mesh, refinement is not permitted (derefining a nonexsiting element not permitted)
     *  Also, if any of the siblings have a higher level of refinement, derefinement is ignored*/

    // cout<<"key is :"<<key<<endl;

    /*!< \brief if the key does not exist simply igonre doing anything*/
    morton<N> sibkey[7], kt;
    uint mylevel;

     level( key, &mylevel );
      
    siblings( key, mylevel, sibkey );
 

    if ( mesh.count( key ) != 0 )
    { // if the siblings

        for ( uint i = 0; i < 7; i++ )
        {
            mesh.erase( sibkey[i] ); 
//            cout<<BLUE<<sibkey[i]<<endl;
        }

        mesh.erase( key );
        kt = key;
        kt[N - 3 * ( mylevel - 1 ) - 1] = 0;
        kt[N - 3 * ( mylevel - 1 ) - 2] = 0;
        kt[N - 3 * ( mylevel - 1 ) - 3] = 0;
        mesh.insert( {kt, nullptr} );
    }
    else
    {
        cout << RED "derefinement ignored, mesh not in the list" RESET << endl;
    }
}

template <size_t N, typename value>
void Tree<N, value>::derefineDerefineList(  )
{
    morton<N> key;
    bool bol;
    bitvector<N> temp;
    morton<N> sibkey[8];

    uint mylevel;
//    uint maxlevel = nlevel;
    uint changedirectionlevel, nbrlevel;
    morton<N> kt;
    retainFourToOne();


    for ( auto i = derefinelist.begin(); i != derefinelist.end(); i++ )
    {
        key = i->first;
//        cout<<key<<" "<<i->second<<endl;
        if(i->second==0)
        {
        derefine( key );

      // cout<<RED<<key<<RESET<<endl;
        removeFromDerefineList(i);
        }
    }
    derefinelist.clear();
}





template <size_t N, typename value>
void Tree<N, value>::play( uint maxlevel, uint geom_nn, real *geom_xyz, double xx[3] )
{
 
convertStl2Morton(geom_nn, geom_xyz);

for ( uint j = 0; j < maxlevel; j++ )
{
  pushToRefinelistSet(j+1);
  refineRefineListSet();
}


cout<<"initial size"<<size()<<endl;
for(uint i=0;i<geom_nn;i++)
{
geom_xyz[3*i+0]+=xx[0];
geom_xyz[3*i+1]+=xx[1];
geom_xyz[3*i+2]+=xx[2];

}

for(uint j=maxlevel;j>maxlevel-3;j--)
{
convertStl2Morton(geom_nn, geom_xyz);
pushToDerefinelist(j+1);
derefineDerefineList();
}


convertStl2Morton(geom_nn, geom_xyz);

for ( uint j = 0; j < maxlevel; j++ )
{

   pushToRefinelistSet(j+1);
   refineRefineListSet();
}   

cout<<"after move size"<<size()<<endl;





}





template class Tree<64, uint>;
template class Tree<64, real>;
template class Tree<32, real>;
template class Tree<32, uint>;
template class Tree<48, uint>;

template class Tree<128, uint>;
