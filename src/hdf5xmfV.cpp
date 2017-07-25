#include "hdf5.h"
#include "header.h"
#define H5FILE_NAME "soln/Pxdmf3d.h5"
#define XDMF_NAME "soln/Pxdmf3d.xmf"

/**
 *
 *              Multi-Block Presentation
 *
 */
//==================================================================
//
//
//
//==================================================================

template <size_t N, typename value>
void Hdf5XmfV<N, value>::writeHdf5MultiBlockSerial( Voxel<N, value> &V )
{
    uint i, j, k, l;
    i = 0;
    hid_t file_id;
    file_id     = H5Fcreate( H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );
    real *xtemp = NULL;
    real *ytemp = NULL;
    real *ztemp = NULL;
    uint  msize = V.npx * V.npy * V.npz;
    xtemp       = new real[msize];
    ytemp       = new real[msize];
    ztemp       = new real[msize];

    // for now put Q here

    real *q = NULL;
    q       = new real[msize];
    real         hx;
    real         hy;
    real         hz;
    real         Xh;
    real         Yh;
    real         Zh;
    unsigned int index;
    real         X[6];

    hid_t   dataset_id, dataspace_id;
    hsize_t dims[3];
    herr_t  status;

    char str0[50];
    char str1[50];
    char str2[50];
    char str3[50];

    uint mylevel;

    uint counter = 0;

    cout << V.maxlevel << endl;

    for ( auto it = V.begin(); it != V.end(); it++ )
    {
        //
        //
        // call get coordinates form
        //
        //

        V.enclosingBox( it->first, X );
        V.level( it->first, &mylevel );

        // cout<<"mylevel=\t"<<mylevel<<"maxlevel="<<maxlevel<<endl;

        if ( mylevel == V.maxlevel )
        {
            hx = V.npx - 1.0;
            hy = V.npy - 1.0;
            hz = V.npz - 1.0;
            Xh = ( X[1] - X[0] ) / ( hx );
            Yh = ( X[3] - X[2] ) / ( hy );
            Zh = ( X[5] - X[4] ) / ( hz );

            // cout<<"dx=\t"<<dx<<"dy\t"<<dy<<"dz"<<dz<<endl;

            // cout<<"i="<<i<<endl;
            index = 0;
            // printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);

            for ( j = 0; j < V.npx; j++ )
            {
                for ( k = 0; k < V.npy; k++ )
                {
                    for ( l = 0; l < V.npz; l++ )
                    {
                        xtemp[index] = X[0] + Xh * j;
                        // printf("%lf\n",xtemp[index]);
                        ytemp[index] = X[2] + Yh * k;
                        ztemp[index] = X[4] + Zh * l;
                        q[index]     = (real)mylevel;
                        index++;
                    }
                }
            }

            sprintf( str0, "/X%d", i );
            sprintf( str1, "/Y%d", i );
            sprintf( str2, "/Z%d", i );
            sprintf( str3, "/Q%d", i );
            i++;

            /*
            printf( "%s\n", str0);
            printf( "%s\n", str1);
            printf( "%s\n", str2);
            printf( "%s\n", str3);
            */

            /* Write separate coordinate arrays for the x and y and z coordinates. */

            dims[0] = V.npx;
            dims[1] = V.npy;
            dims[2] = V.npz;

            if ( sizeof( real ) == sizeof( double ) )
            {
                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xtemp );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );

                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str1, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ytemp );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );

                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str2, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ztemp );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );

                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str3, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );
            }
            else if ( sizeof( real ) == sizeof( float ) )
            {
                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str0, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xtemp );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );

                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str1, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ytemp );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );

                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str2, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ztemp );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );

                dataspace_id = H5Screate_simple( 3, dims, NULL );
                dataset_id   = H5Dcreate( file_id, str3, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
                status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, q );
                status       = H5Dclose( dataset_id );
                status       = H5Sclose( dataspace_id );
            }

            else
            {
                cout << RED "error in defining type in hdf5" RESET << endl;
                exit( 0 );
            }
        }
    }
    delete[] xtemp;
    delete[] ytemp;
    delete[] ztemp;
    delete[] q;
}

template <size_t N, typename value>
void Hdf5XmfV<N, value>::xdmfMultiBlockSerial( Voxel<N, value> &V )
{
    unsigned int i;
    char         str0[50];
    char         str1[50];
    char         str2[50];
    char         str3[50];
    char         str4[50];
    char         str5[50];

    FILE *xmf = NULL;
    xmf       = fopen( XDMF_NAME, "w" );
    fprintf( xmf, "<?xml version=\"1.0\" ?>\n" );
    fprintf( xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" );
    fprintf( xmf, "<Xdmf Version=\"2.0\">\n" );
    fprintf( xmf, " <Domain>\n" );
    fprintf( xmf, "   <Grid Name=\"AMR0\" GridType=\"Collection\" CollectionType=\"Spatial\">\n" );
    char fname[50] = "Pxdmf3d.h5";
    uint mylevel, counter = 0;

    for ( auto it = V.begin(); it != V.end(); it++ )
    {
        V.level( it->first, &mylevel );

        if ( mylevel == V.maxlevel )
        {
            sprintf( str0, "/X%d", counter );
            sprintf( str1, "/Y%d", counter );
            sprintf( str2, "/Z%d", counter );
            sprintf( str3, "mesh%d", counter );
            sprintf( str4, "/Q%d", counter );
            /* printf( "%s\n", str0);
             printf( "%s\n", str1);
             printf( "%s\n", str2);
             printf( "%s\n", str3);
             printf( "%s\n", str4);  */
            fprintf( xmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n", str3 );
            fprintf( xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", V.npx, V.npy, V.npz );
            fprintf( xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", V.npx, V.npy,
                     V.npz );
            fprintf( xmf, "        %s:%s\n", fname, str0 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", V.npx, V.npy,
                     V.npz );
            fprintf( xmf, "        %s:%s\n", fname, str1 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", V.npx, V.npy,
                     V.npz );
            fprintf( xmf, "        %s:%s\n", fname, str2 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "     </Geometry>\n" );
            // add your variables here
            fprintf( xmf, "     <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\">\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", V.npx, V.npy,
                     V.npz );
            fprintf( xmf, "        %s:%s\n", fname, str4 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "     </Attribute>\n" );
            fprintf( xmf, "   </Grid>\n" );
            counter++;
            // cout<<counter<<endl;
        }
    }
    fprintf( xmf, "   </Grid>\n" );
    fprintf( xmf, " </Domain>\n" );
    fprintf( xmf, "</Xdmf>\n" );
    fclose( xmf );
}

template <size_t N, typename value>
void Hdf5XmfV<N, value>::write( Voxel<N, value> &V )
{
    writeHdf5MultiBlockSerial( V );
    xdmfMultiBlockSerial( V );
}

//===================================
//
// polyvetrex
//
//=================================
template <size_t N, typename value>
void Hdf5XmfV<N, value>::writeHdf5PolyvertexSerial( Voxel<N, value> &V )
{
    hid_t file_id;
    file_id = H5Fcreate( H5FILE_NAME, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT );

    cout << "mesh size" << V.size() << endl;

    real *xtemp = NULL;
    real *ytemp = NULL;
    real *ztemp = NULL;
    xtemp       = new real[V.size()];
    ytemp       = new real[V.size()];
    ztemp       = new real[V.size()];

    real *q = NULL;
    q       = new real[V.size()];

    hsize_t index;

    hid_t   dataset_id, dataspace_id;
    hsize_t dims[1];
    herr_t  status;

    char str0[50] = "X";
    char str1[50] = "Y";
    char str2[50] = "Z";
    char str3[50] = "Q";
    real xyz[3];

    index = 0;
    uint mylevel;

    for ( auto it = V.begin(); it != V.end(); it++ )
    {
        //
        //
        // call get coordinates form
        //
        //

        V.level( it->first, &mylevel );

        // cout<<"mylevel=\t"<<mylevel<<"maxlevel="<<maxlevel<<endl;

        if ( mylevel == V.maxlevel )
        {
            //
            //
            // call get coordinates form
            //
            //

            V.centroid( it->first, xyz );
            // cout<<xyz[0]<<"\t"<<xyz[1]<<"\t"<<xyz[2]<<endl;

            // printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);

            xtemp[index] = xyz[0];
            ytemp[index] = xyz[1];
            ztemp[index] = xyz[2];
            q[index]     = (real)mylevel;
            index++;
        }

    } /*
  sprintf(str0, "/X%d", 0);
  sprintf(str1, "/Y%d", 0);
  sprintf(str2, "/Z%d", 0);
  sprintf(str3, "/Q%d", 0);
  /*
  printf( "%s\n", str0);
  printf( "%s\n", str1);
  printf( "%s\n", str2);
  printf( "%s\n", str3);
  */

    /* Write separate coordinate arrays for the x and y and z coordinates. */

    dims[0] = index;

    if ( sizeof( real ) == sizeof( double ) )
    {
        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str0, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, xtemp );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );

        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str1, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ytemp );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );

        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str2, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, ztemp );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );

        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str3, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );
    }

    else if ( sizeof( real ) == sizeof( float ) )
    {
        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str0, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, xtemp );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );

        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str1, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ytemp );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );

        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str2, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, ztemp );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );

        dataspace_id = H5Screate_simple( 1, dims, NULL );
        dataset_id   = H5Dcreate( file_id, str3, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT );
        status       = H5Dwrite( dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, q );
        status       = H5Dclose( dataset_id );
        status       = H5Sclose( dataspace_id );
    }

    else
    {
        cout << RED "error in defining type in hdf5" RESET << endl;
        exit( 0 );
    }

    V.numMax = index;

    delete[] xtemp;
    delete[] ytemp;
    delete[] ztemp;
    delete[] q;
}

template <size_t N, typename value>
void Hdf5XmfV<N, value>::xdmfPolyvertexSerial( Voxel<N, value> &V )
{
    FILE *fp = NULL;
    fp       = fopen( XDMF_NAME, "w" );

    const char *names[] = {"X", "Y", "Z"};

    char str[1000];
    uint ncube_total = V.numMax;

    //  if(my_rank==0)
    {
        sprintf( str, "<?xml version=\"1.0\" ?>\n" );
        fprintf( fp, "%s", str );

        sprintf( str, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n" );
        fprintf( fp, "%s", str );

        // MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
        sprintf( str, "<Xdmf Version=\"2.0\">\n" );
        fprintf( fp, "%s", str );

        sprintf( str, "<Domain>\n" );
        fprintf( fp, "%s", str );
        sprintf( str, "   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n" );
        fprintf( fp, "%s", str );

        sprintf( str, "        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n", ncube_total );
        fprintf( fp, "%s", str );
        sprintf( str, "          <Geometry GeometryType=\"X_Y_Z\">  \n" );
        fprintf( fp, "%s", str );
        sprintf( str, "          <DataItem Name=\"X\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                 ncube_total );
        fprintf( fp, "%s", str );
        sprintf( str, "          Pxdmf3d.h5:/%s\n", names[0] );
        fprintf( fp, "%s", str );
        sprintf( str, "         </DataItem>\n" );
        fprintf( fp, "%s", str );
        sprintf( str, "          <DataItem Name=\"Y\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                 ncube_total );
        fprintf( fp, "%s", str );
        sprintf( str, "          Pxdmf3d.h5:/%s\n", names[1] );
        fprintf( fp, "%s", str );
        sprintf( str, "         </DataItem>\n" );
        fprintf( fp, "%s", str );
        sprintf( str, "          <DataItem Name=\"Z\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                 ncube_total );
        fprintf( fp, "%s", str );
        sprintf( str, "          Pxdmf3d.h5:/%s\n", names[2] );
        fprintf( fp, "%s", str );
        sprintf( str, "         </DataItem>\n" );
        fprintf( fp, "%s", str );
        sprintf( str, "      </Geometry>   \n" );
        fprintf( fp, "%s", str );
        sprintf( str, "         <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\"> \n" );
        fprintf( fp, "%s", str );
        sprintf( str, "          <DataItem Name=\"Q\" Dimensions=\"%d \" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",
                 ncube_total );
        fprintf( fp, "%s", str );
        sprintf( str, "          Pxdmf3d.h5:/Q\n" );
        fprintf( fp, "%s", str );
        sprintf( str, "         </DataItem>\n" );
        fprintf( fp, "%s", str );
        sprintf( str, " </Attribute>\n" );
        fprintf( fp, "%s", str );
        sprintf( str, "  </Grid>\n" );
        fprintf( fp, "%s", str );
        sprintf( str, "      </Domain>   \n" );
        fprintf( fp, "%s", str );
        sprintf( str, "  </Xdmf>\n" );
        fprintf( fp, "%s", str );
    }

    fclose( fp );
}

/////////////
//
//
//
//

template <size_t N, typename value>
void Hdf5XmfV<N, value>::writeP( Voxel<N, value> &V )
{
    writeHdf5PolyvertexSerial( V );
    xdmfPolyvertexSerial( V );
}

template class Hdf5XmfV<32, real>;

template class Hdf5XmfV<64, real>;


template class Hdf5XmfV<64, uint>;
