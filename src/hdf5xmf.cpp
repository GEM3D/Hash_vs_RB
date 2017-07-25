#include "header.h"
#include "hdf5.h"
#define H5FILE_NAME  "soln/Pxdmf3d%u.h5"
#define XDMF_NAME    "soln/Pxdmf3d%u.xmf" 
#define H5FILE  "Pxdmf3d%u.h5"

/*
template<size_t N,typename value>
void Hdf5Xmf<N,value>::Msize(Tree<N,value>& C)
{
cout<<C.mesh.size()<<endl;
}
*/

#if(1)

template<size_t N, typename value>
void Hdf5Xmf<N,value>::xdmfPolyvertexSerial(Tree<N,value>& T ,uint appx)
{
  const char* names[] = {"X", "Y", "Z"};

  char str[1000],fname[1000];
  uint ncube_total=T.mesh.size();

  FILE *fp=NULL;
  sprintf(str,XDMF_NAME,appx);
  fp=fopen(str,"w" );
  //printf("%s\n",str);
  
  sprintf(fname,H5FILE,appx);
  

//  if(my_rank==0)
    {


      sprintf(str,"<?xml version=\"1.0\" ?>\n");
      fprintf(fp,"%s", str);

      sprintf(str,"<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []> \n");
      fprintf(fp,"%s", str);

      //MPI_File_write(fp,str,strlen(str), MPI_CHAR,&status);
      sprintf(str,"<Xdmf Version=\"2.0\">\n");
      fprintf(fp,"%s", str);

      sprintf(str,"<Domain>\n");
       fprintf(fp,"%s", str);
      sprintf(str,"   <Grid Name=\"mesh0\" GridType=\"Uniform\">\n");
      fprintf(fp,"%s", str);

          sprintf(str,"        <Topology TopologyType=\"Polyvertex\" NumberOfElements=\"%d\"/>\n",ncube_total);
          fprintf(fp,"%s", str);
          sprintf(str,"          <Geometry GeometryType=\"X_Y_Z\">  \n");
          fprintf(fp,"%s", str);
          sprintf(str,"          <DataItem Name=\"X\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
          fprintf(fp,"%s", str);
          sprintf(str,"          %s:/%s\n",fname,names[0]);
          fprintf(fp,"%s", str);
          sprintf(str,"         </DataItem>\n");
          fprintf(fp,"%s", str);
          sprintf(str,"          <DataItem Name=\"Y\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
          fprintf(fp,"%s", str);
          sprintf(str,"          %s:/%s\n",fname,names[1]);
          fprintf(fp,"%s", str);
          sprintf(str,"         </DataItem>\n");
          fprintf(fp,"%s", str);
          sprintf(str,"          <DataItem Name=\"Z\" Dimensions=\"%d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
          fprintf(fp,"%s", str);
          sprintf(str,"          %s:/%s\n",fname,names[2]);
          fprintf(fp,"%s", str);
          sprintf(str,"         </DataItem>\n");
          fprintf(fp,"%s", str);
          sprintf(str,"      </Geometry>   \n");
          fprintf(fp,"%s", str);
          sprintf(str,"         <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\"> \n");
          fprintf(fp,"%s", str);
          sprintf(str,"          <DataItem Name=\"Q\" Dimensions=\"%d \" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",ncube_total);
          fprintf(fp,"%s", str);
          sprintf(str,"          %s:/Q\n",fname);
          fprintf(fp,"%s", str);
          sprintf(str,"         </DataItem>\n");
          fprintf(fp,"%s", str);
          sprintf(str," </Attribute>\n");
          fprintf(fp,"%s", str);
          sprintf(str,"  </Grid>\n");
          fprintf(fp,"%s", str);
          sprintf(str,"      </Domain>   \n");
          fprintf(fp,"%s", str);
          sprintf(str,"  </Xdmf>\n");
          fprintf(fp,"%s", str);
    }

      fclose(fp);

}
#endif

//=================================================================
//
// polyvertex 
//================================================================

template<size_t N, typename value>
void Hdf5Xmf<N,value>::writeHdf5PolyvertexSerial(Tree<N,value>& T, uint appx)
{

    hid_t     file_id;
 
    char str[1000];
 
    sprintf(str,H5FILE_NAME,appx);
    //printf("%s\n",str);   
    file_id = H5Fcreate(str, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    real *xtemp =NULL;
    real *ytemp =NULL;
    real *ztemp =NULL;
    xtemp=new real[T.mesh.size()];
    ytemp=new real[T.mesh.size()];
    ztemp=new real[T.mesh.size()];

    real *q =NULL;
    q=new real[T.mesh.size()];

    hsize_t index;

    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[1];
    herr_t    status;

        char str0[50]="X";
        char str1[50]="Y";
        char str2[50]="Z";
        char str3[50]="Q";
        real xyz[3];

         index=0;
  uint mylevel;

        for(auto it=T.mesh.begin();it!=T.mesh.end();it++)
        {

    // 
    // call get coordinates form
    //

   T.level(it->first,&mylevel);
   T.centroid(it->first,xyz);
// cout<<xyz[0]<<"\t"<<xyz[1]<<"\t"<<xyz[2]<<endl;
//printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);

               xtemp[index]=xyz[0];
               ytemp[index]=xyz[1];
               ztemp[index]=xyz[2];
               q[index]=(real)mylevel;
              index++;
        }
    /* Write separate coordinate arrays for the x and y and z coordinates. */

      dims[0]=index;

if(sizeof(real)==sizeof(double))
{
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, xtemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate(file_id, str1, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ytemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate(file_id, str2, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ztemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

       dataspace_id = H5Screate_simple(1, dims, NULL);
       dataset_id = H5Dcreate(file_id, str3, H5T_NATIVE_DOUBLE,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
       status = H5Dclose(dataset_id);
       status = H5Sclose(dataspace_id);
}
else if(sizeof(real)==sizeof(float))
{
        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate(file_id, str0, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, xtemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate(file_id, str1, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, ytemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(1, dims, NULL);
        dataset_id = H5Dcreate(file_id, str2, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, ztemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

       dataspace_id = H5Screate_simple(1, dims, NULL);
       dataset_id = H5Dcreate(file_id, str3, H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
       status = H5Dclose(dataset_id);
       status = H5Sclose(dataspace_id);
}
else
{
 
//throw std::runtime_error("Error in defining type in hdf5");

cout<<RED "error in defining type in hdf5" RESET<<endl;
exit(0);
}


       // printf("%d\n",i);
        delete[] xtemp  ;
        delete[] ytemp  ;
        delete[] ztemp  ;
        delete[] q;
}



template<size_t N, typename value>
void Hdf5Xmf<N,value>::writeP(Tree<N,value> &T, uint appx)
{
writeHdf5PolyvertexSerial(T,appx);
xdmfPolyvertexSerial(T,appx);
}




template<size_t N, typename value>
void Hdf5Xmf<N,value>::writeHdf5MultiBlockSerial(Tree<N,value>&T,uint appx)
{

        uint i,j,k,l;
        i=0;

        hid_t     file_id;

 char str[1000];

 sprintf(str,H5FILE_NAME,appx);
   // printf("%s\n",str);

    file_id = H5Fcreate(str, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);


        real *xtemp =NULL;
    real *ytemp =NULL;
    real *ztemp =NULL;
    uint msize=T.npx*T.npy*T.npz;
    xtemp=new real[msize];
    ytemp=new real[msize];
    ztemp=new real[msize];

// for now put Q here

    real *q =NULL;
    q=new real[msize];
        real hx;
        real hy;
        real hz;
        real Xh;
        real Yh;
        real Zh;
    unsigned int index;
    real X[6];

    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[3];
    herr_t    status;

        char str0[50];
        char str1[50];
        char str2[50];
        char str3[50];

 uint mylevel;


  for(auto it=T.mesh.begin();it!=T.mesh.end();it++)
        {

       T.enclosingBox(it->first,X);
       T.level(it->first,&mylevel);
        hx=T.npx-1.0;
         hy=T.npy-1.0;
         hz=T.npz-1.0;
        Xh=(X[1]-X[0])/(hx);
         Yh=(X[3]-X[2])/(hy);
         Zh=(X[5]-X[4])/(hz);

//cout<<"dx=\t"<<dx<<"dy\t"<<dy<<"dz"<<dz<<endl;

//cout<<"i="<<i<<endl;  
    index=0;
   // printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);

        for(j=0;j<T.npx;j++)
        {
                for(k=0;k<T.npy;k++)
          {
                for(l=0;l<T.npz;l++)
            {
               xtemp[index]=X[0]+Xh*j;
               //printf("%lf\n",xtemp[index]);      
               ytemp[index]=X[2]+Yh*k;
               ztemp[index]=X[4]+Zh*l;
               q[index]=(real)mylevel;
            index++;
        }
      }
        }

    sprintf(str0, "/X%d", i);
    sprintf(str1, "/Y%d", i);
    sprintf(str2, "/Z%d", i);
    sprintf(str3, "/Q%d", i);
    i++;


    /*
    printf( "%s\n", str0);
    printf( "%s\n", str1);
    printf( "%s\n", str2);
    printf( "%s\n", str3);
    */

    /* Write separate coordinate arrays for the x and y and z coordinates. */

        dims[0] = T.npx;
        dims[1] = T.npy;
        dims[2] = T.npz;

// hdf5 does not use template for native datatyppe so I have to code  this section twice 
if(sizeof(real)==sizeof(double))
{
        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, xtemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str1, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ytemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str2, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ztemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);


       dataspace_id = H5Screate_simple(3, dims, NULL);
       dataset_id = H5Dcreate(file_id, str3, H5T_NATIVE_DOUBLE,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
       status = H5Dclose(dataset_id);
      status = H5Sclose(dataspace_id);
}
else if(sizeof(real)==sizeof(float))
   { 
       dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str0, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, xtemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str1, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, ytemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str2, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, ztemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);


       dataspace_id = H5Screate_simple(3, dims, NULL);
       dataset_id = H5Dcreate(file_id, str3, H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
       status = H5Dclose(dataset_id);
       status = H5Sclose(dataspace_id);
}
else
{
cout<<RED "error in defining type in hdf5" RESET<<endl;
exit(0);
}
       // printf("%d\n",i);

        }

//      XdmfMultiBlockSerial(mesh,L,M, N);

        //printf("%d\n",cube.size());
        delete[] xtemp  ;
        delete[] ytemp  ;
        delete[] ztemp  ;
        delete[] q;

}

template <size_t N, typename value>
void Hdf5Xmf<N,value>::xdmfMultiBlockSerial(Tree<N,value>&T,uint appx)
{
        unsigned int i;
        char str0[50];
        char str1[50];
        char str2[50];
        char str3[50];
        char str4[50];
        char str5[50];

    FILE *xmf=NULL;

  sprintf(str0,XDMF_NAME,appx);

char fname[1000];

  sprintf(fname,H5FILE,appx);

    xmf = fopen(str0, "w");
    

    fprintf(xmf, "<?xml version=\"1.0\" ?>\n");
    fprintf(xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n");
    fprintf(xmf, "<Xdmf Version=\"2.0\">\n");
    fprintf(xmf, " <Domain>\n");
    fprintf(xmf, "   <Grid Name=\"AMR0\" GridType=\"Collection\" CollectionType=\"Spatial\">\n");

    for(i=0;i<T.mesh.size();i++)
   //   for(i=0;i<2;i++)
    {
          sprintf(str0, "/X%d", i);
          sprintf(str1, "/Y%d", i);
      sprintf(str2, "/Z%d", i);
      sprintf(str3, "mesh%d", i);
      sprintf(str4, "/Q%d", i);
   /* printf( "%s\n", str0);
    printf( "%s\n", str1);
    printf( "%s\n", str2);
    printf( "%s\n", str3);  
    printf( "%s\n", str4);  */
        fprintf(xmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n",str3);
        fprintf(xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", T.npx, T.npy, T.npz);
    fprintf(xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n",  T.npx, T.npy, T.npz);
    fprintf(xmf, "        %s:%s\n",fname,str0);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", T.npx, T.npy, T.npz);
    fprintf(xmf, "        %s:%s\n",fname,str1);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", T.npx, T.npy, T.npz);
    fprintf(xmf, "        %s:%s\n",fname,str2);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Geometry>\n");
    // add your variables here
    fprintf(xmf, "     <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\">\n");
    fprintf(xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", T.npx, T.npy, T.npz);
    fprintf(xmf, "        %s:%s\n",fname,str4);
    fprintf(xmf, "       </DataItem>\n");
    fprintf(xmf, "     </Attribute>\n");
    fprintf(xmf, "   </Grid>\n");
        }
    fprintf(xmf, "   </Grid>\n");
    fprintf(xmf, " </Domain>\n");
    fprintf(xmf, "</Xdmf>\n");
    fclose(xmf);



}

template<size_t N,typename value>
void Hdf5Xmf<N,value>::write(Tree<N,value>&T,uint appx)
{

writeHdf5MultiBlockSerial(T,appx);
xdmfMultiBlockSerial(T,appx);

}


template<size_t N, typename value>
void Hdf5Xmf<N,value>::writeHdf5MultiBlockHighLevel(Tree<N,value>&T,uint appx)
{

        uint i,j,k,l;
        i=0;

        hid_t     file_id;
   
 char str[1000];

 sprintf(str,H5FILE_NAME,appx);
  //  printf("%s\n",str);

    file_id = H5Fcreate(str, H5F_ACC_TRUNC, H5P_DEFAULT, H5P_DEFAULT);

    real *xtemp =NULL;
    real *ytemp =NULL;
    real *ztemp =NULL;
    uint msize=T.npx*T.npy*T.npz;
    xtemp=new real[msize];
    ytemp=new real[msize];
    ztemp=new real[msize];

// for now put Q here

    real *q =NULL;
    q=new real[msize];
        real hx;
        real hy;
        real hz;
        real Xh;
        real Yh;
        real Zh;
    unsigned int index;
    real X[6];

    hid_t     dataset_id, dataspace_id;
    hsize_t   dims[3];
    herr_t    status;

        char str0[50];
        char str1[50];
        char str2[50];
        char str3[50];

 uint mylevel,maxlevel;
 i=0;

  maxlevel=0;
  for(auto it=T.mesh.begin();it!=T.mesh.end();it++)
{

    T.level(it->first,&mylevel);
   if(mylevel>maxlevel)
{
maxlevel=mylevel;
}

}

  for(auto it=T.mesh.begin();it!=T.mesh.end();it++)
{
     T.enclosingBox(it->first,X);

    T.level(it->first,&mylevel);

if(mylevel==maxlevel)
{
        hx=T.npx-1.0;
         hy=T.npy-1.0;
         hz=T.npz-1.0;
        Xh=(X[1]-X[0])/(hx);
         Yh=(X[3]-X[2])/(hy);
         Zh=(X[5]-X[4])/(hz);

//cout<<"dx=\t"<<dx<<"dy\t"<<dy<<"dz"<<dz<<endl;

//cout<<"i="<<i<<endl;  
    index=0;
   // printf("%lf %lf %lf %lf %lf %lf\n", Xa,Xb,Ya,Yb,Za,Zb);

        for(j=0;j<T.npx;j++)
        {
                for(k=0;k<T.npy;k++)
          {
                for(l=0;l<T.npz;l++)
            {
               xtemp[index]=X[0]+Xh*j;
               //printf("%lf\n",xtemp[index]);      
               ytemp[index]=X[2]+Yh*k;
               ztemp[index]=X[4]+Zh*l;
               q[index]=(real)mylevel;
            index++;
        }
      }
        }

    sprintf(str0, "/X%d", i);
    sprintf(str1, "/Y%d", i);
    sprintf(str2, "/Z%d", i);
    sprintf(str3, "/Q%d", i);
    i++;
   /* Write separate coordinate arrays for the x and y and z coordinates. */

        dims[0] = T.npx;
        dims[1] = T.npy;
        dims[2] = T.npz;

if(sizeof(real)==sizeof(double))
{
        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str0, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, xtemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str1, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ytemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

            dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str2, H5T_NATIVE_DOUBLE, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL,H5P_DEFAULT, ztemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);


       dataspace_id = H5Screate_simple(3, dims, NULL);
       dataset_id = H5Dcreate(file_id, str3, H5T_NATIVE_DOUBLE,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       status = H5Dwrite(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
       status = H5Dclose(dataset_id);
       status = H5Sclose(dataspace_id);
}
else if(sizeof(real)==sizeof(float))
{
        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str0, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, xtemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

        dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str1, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, ytemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

            dataspace_id = H5Screate_simple(3, dims, NULL);
        dataset_id = H5Dcreate(file_id, str2, H5T_NATIVE_FLOAT, dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
        status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL,H5P_DEFAULT, ztemp);
        status = H5Dclose(dataset_id);
        status = H5Sclose(dataspace_id);

       dataspace_id = H5Screate_simple(3, dims, NULL);
       dataset_id = H5Dcreate(file_id, str3, H5T_NATIVE_FLOAT,dataspace_id, H5P_DEFAULT, H5P_DEFAULT, H5P_DEFAULT);
       status = H5Dwrite(dataset_id, H5T_NATIVE_FLOAT, H5S_ALL, H5S_ALL, H5P_DEFAULT, q);
       status = H5Dclose(dataset_id);
       status = H5Sclose(dataspace_id);



}
else
{
cout<<RED "error in defining type in hdf5" RESET<<endl;
exit(0);

}
 }
}

//      XdmfMultiBlockSerial(mesh,L,M, N);

// printf("%d\n",cube.size());
delete[] xtemp;
delete[] ytemp;
delete[] ztemp;
delete[] q;
}

template <size_t N, typename value>
void Hdf5Xmf<N, value>::xdmfMultiBlockHighLevel( Tree<N, value> &T, uint appx )
{
    uint i = 0, mylevel;
    char str0[50];
    char str1[50];
    char str2[50];
    char str3[50];
    char str4[50];
    char str5[50];

    FILE *xmf = NULL;

    sprintf( str0, XDMF_NAME, appx );
    xmf = fopen( str0, "w" );

    char fname[1000];

    sprintf( fname, H5FILE, appx );

    fprintf( xmf, "<?xml version=\"1.0\" ?>\n" );
    fprintf( xmf, "<!DOCTYPE Xdmf SYSTEM \"Xdmf.dtd\" []>\n" );
    fprintf( xmf, "<Xdmf Version=\"2.0\">\n" );
    fprintf( xmf, " <Domain>\n" );
    fprintf( xmf, "   <Grid Name=\"AMR0\" GridType=\"Collection\" CollectionType=\"Spatial\">\n" );

    uint maxlevel = 0;
    for ( auto it = T.mesh.begin(); it != T.mesh.end(); it++ )
    {
        T.level( it->first, &mylevel );

        if ( mylevel > maxlevel )
        {
            maxlevel = mylevel;
        }
    }

    for ( auto it = T.begin(); it != T.end(); it++ )
    //   for(i=0;i<2;i++)
    {
        T.level( it->first, &mylevel );
        if ( mylevel == maxlevel )
        {
            sprintf( str0, "/X%d", i );
            sprintf( str1, "/Y%d", i );
            sprintf( str2, "/Z%d", i );
            sprintf( str3, "mesh%d", i );
            sprintf( str4, "/Q%d", i );
            /* printf( "%s\n", str0);
             printf( "%s\n", str1);
             printf( "%s\n", str2);
             printf( "%s\n", str3);
             printf( "%s\n", str4);  */
            fprintf( xmf, "   <Grid Name=\"%s\" GridType=\"Uniform\">\n", str3 );
            fprintf( xmf, "     <Topology TopologyType=\"3DSMesh\" NumberOfElements=\"%d %d %d\"/>\n", T.npx, T.npy, T.npz );
            fprintf( xmf, "     <Geometry GeometryType=\"X_Y_Z\">\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", T.npx, T.npy,
                     T.npz );
            fprintf( xmf, "        %s:%s\n", fname, str0 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", T.npx, T.npy,
                     T.npz );
            fprintf( xmf, "        %s:%s\n", fname, str1 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", T.npx, T.npy,
                     T.npz );
            fprintf( xmf, "        %s:%s\n", fname, str2 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "     </Geometry>\n" );
            // add your variables here
            fprintf( xmf, "     <Attribute Name=\"Q\" AttributeType=\"Scalar\" Center=\"Node\">\n" );
            fprintf( xmf, "       <DataItem Dimensions=\"%d %d %d\" NumberType=\"Float\" Precision=\"4\" Format=\"HDF\">\n", T.npx, T.npy,
                     T.npz );
            fprintf( xmf, "        %s:%s\n", fname, str4 );
            fprintf( xmf, "       </DataItem>\n" );
            fprintf( xmf, "     </Attribute>\n" );
            fprintf( xmf, "   </Grid>\n" );

            i++;
        }
    }
    fprintf( xmf, "   </Grid>\n" );
    fprintf( xmf, " </Domain>\n" );
    fprintf( xmf, "</Xdmf>\n" );
    fclose( xmf );
}

template <size_t N, typename value>
void Hdf5Xmf<N, value>::writeH( Tree<N, value> &T, uint appx )
{
    writeHdf5MultiBlockHighLevel( T, appx );
    xdmfMultiBlockHighLevel( T, appx );
}

template class Hdf5Xmf<32, uint>;

template class Hdf5Xmf<64, uint>;


template class Hdf5Xmf<48, uint>;
