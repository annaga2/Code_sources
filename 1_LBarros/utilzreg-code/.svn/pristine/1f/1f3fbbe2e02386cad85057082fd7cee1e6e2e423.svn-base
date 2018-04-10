/*=========================================================================
 
 Author: Laurent Risser, Francois-Xavier Vialard
 
 Disclaimer: This software has been developed for research purposes only, and hence should 
 not be used as a diagnostic tool. In no event shall the authors or distributors
 be liable to any direct, indirect, special, incidental, or consequential 
 damages arising of the use of this software, its documentation, or any 
 derivatives thereof, even if the authors have been advised of the possibility 
 of such damage. 
 
 
 =========================================================================*/

#include <SciCalcPack.h>



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           1:   FUNCTIONS FOR THE CLASS "ScalarField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///constructor
ScalarField::ScalarField(void){
  this->NX=0;
  this->NY=0;
  this->NZ=0;
  this->NT=0;
}

///destructor
ScalarField::~ScalarField(void){
  if ((this->ScalField!=NULL)&&(this->NX>0)) delete this->ScalField;
}


///put a value
//-> inline function in the .h file

/// add a value
//-> inline function in the .h file

///put a the same value at every points of the scalar field
void ScalarField::PutToAllVoxels(float cste,int t)
{
  int x,y,z;
  for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) { this->P(cste,x,y,z,t); }
}

///get a value
//-> inline function in the .h file

///get a value using linear interpolation
float ScalarField::G(float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
	//values out of the image
	if (x<0.) x=0.0001;
	if (x>=this->NX-1.) x=this->NX-1.0001;
	if (y<0.) y=0.0001;
	if (y>=this->NY-1.) y=this->NY-1.0001;
	if (z<0.) z=0.0001;
	if (z>=this->NZ-1.) z=this->NZ-1.0001;
	if (t<0) t=0;
	if (t>this->NT-1) t=this->NT-1;
	
	//closest entire value
	xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
	yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
	zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
	
	//interpolation
	if (this->NZ==1){ //2D IMAGE
		wmm=xwm*ywm;
		wmp=xwm*ywp;
		wpm=xwp*ywm;
		wpp=xwp*ywp;
		
		InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
	}
	else{//3D IMAGE
		wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
		wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
		
		InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
	}
	
	return InterpoGreyLevel;
}


///same as above
float ScalarField::G(double x,double y,double z,int t){
	return this->G((float)x,(float)y,(float)z,t);
}

///get a value using linear interpolation. No extrapolation if out of field coordinates.
float ScalarField::G_NoExtrapo(float x,float y,float z,int t){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
	//values out of the image
	if (x<0.) return 0;
	if (x>=this->NX-1.) return 0;
	if (y<0.)           return 0;
	if (y>=this->NY-1.) return 0;
	if (z<0.)           return 0;
	if (z>=this->NZ-1.) return 0;
	if (t<0)            return 0;
	if (t>this->NT-1)   return 0;
	
	//closest entire value
	xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
	yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
	zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
	
	//interpolation
	if (this->NZ==1){ //2D IMAGE
		wmm=xwm*ywm;
		wmp=xwm*ywp;
		wpm=xwp*ywm;
		wpp=xwp*ywp;
		
		InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
	}
	else{//3D IMAGE
		wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
		wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
		
		InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
	}
	
	return InterpoGreyLevel;
}


///same as above
float ScalarField::G_NoExtrapo(double x,double y,double z,int t){
	return this->G_NoExtrapo((float)x,(float)y,(float)z,t);
}
 

///get a value using linear interpolation
float ScalarField::G(float coordSpace2imageSpace[4][4],float x,float y,float z,int t,int NN){
  float InterpoGreyLevel;
  int xi,yi,zi;
  float xt,yt,zt;
  float xwm,ywm,zwm,xwp,ywp,zwp;
  float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
  float wmm,wmp,wpm,wpp;
  
  
  //transform x,y,z to the image space
  xt=x;
  yt=y;
  zt=z;
  
  x=xt*coordSpace2imageSpace[0][0]+yt*coordSpace2imageSpace[0][1]+zt*coordSpace2imageSpace[0][2]+coordSpace2imageSpace[0][3];
  y=xt*coordSpace2imageSpace[1][0]+yt*coordSpace2imageSpace[1][1]+zt*coordSpace2imageSpace[1][2]+coordSpace2imageSpace[1][3];
  z=xt*coordSpace2imageSpace[2][0]+yt*coordSpace2imageSpace[2][1]+zt*coordSpace2imageSpace[2][2]+coordSpace2imageSpace[2][3];
  
  if (NN==0){ //A) trilinear interpolation
    //values out of the image
    if (x<0.) x=0.0001;
    if (x>=this->NX-1.) x=this->NX-1.0001;
    if (y<0.) y=0.0001;
    if (y>=this->NY-1.) y=this->NY-1.0001;
    if (z<0.) z=0.0001;
    if (z>=this->NZ-1.) z=this->NZ-1.0001;
    if (t<0) t=0;
    if (t>this->NT-1) t=this->NT-1;
    
    //closest entire value
    xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
    yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
    zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
    
    //interpolation
    if (this->NZ==1){ //2D IMAGE
      wmm=xwm*ywm;
      wmp=xwm*ywp;
      wpm=xwp*ywm;
      wpp=xwp*ywp;
      
      InterpoGreyLevel= wmm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wpm*this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wpp*this->ScalField[ t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
    }
    else{//3D IMAGE
      wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
      wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
      
      InterpoGreyLevel= wmmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmpm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wmpp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
      InterpoGreyLevel+=wpmm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wpmp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wppm*this->ScalField[ t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
      InterpoGreyLevel+=wppp*this->ScalField[ t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
    }
  }
  else{//B) nearest neighbor interpolation
      //values out of the image
      xi=static_cast<int>(x+0.5);
      yi=static_cast<int>(y+0.5);
      zi=static_cast<int>(z+0.5);
      
      if (xi<0) xi=0;  if (xi>=this->NX) xi=this->NX-1;  
      if (yi<0) yi=0;  if (yi>=this->NY) yi=this->NY-1;  
      if (zi<0) zi=0;  if (zi>=this->NZ) zi=this->NZ-1;  
      
      InterpoGreyLevel=this->ScalField[ t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
  }
  
  return InterpoGreyLevel;
}


///same as above
float ScalarField::G(float coordSpace2imageSpace[4][4],double x,double y,double z,int t,int NN){
  return this->G(coordSpace2imageSpace,(float)x,(float)y,(float)z,t,NN);
}

///same as above
float ScalarField::G(float coordSpace2imageSpace[4][4],int x,int y,int z,int t,int NN){
  return this->G(coordSpace2imageSpace,(float)x,(float)y,(float)z,t,NN);
}



///get the maximum absolute values out of the scalar field
float ScalarField::GetMaxAbsVal(int t){
	float max=0.0;
	int x,y,z;
	for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
	{
		if(max<abs(this->G(x,y,z,t))){max = abs(this->G(x,y,z,t));}
	}
	return max;
}


///read a scalar field (in a nifti image)
void ScalarField::Read(char * ImageName){
  nifti_1_header hdr;
  FILE *fp;
  int ret,i;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z,t;
  float a,b,c,d,qfac;

  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  //to open some multichannel 3D images which have the channels in the 5th dimension instead of the 4th dimension
  if ((hdr.dim[4]==1)&&(hdr.dim[5]>1)){
    hdr.dim[4]=hdr.dim[5];
    hdr.dim[5]=1;
  }
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=hdr.dim[1])||(this->NY!=hdr.dim[2])||(this->NZ!=hdr.dim[3])||(this->NT!=hdr.dim[4]))
      cout << "WARNING: THE SIZE OF A NON-NULL SCALAR FIELD IS CHANGED\n";
  
  
  //fill the parameters of the class
  this->NX=static_cast<int>(hdr.dim[1]);
  this->NY=static_cast<int>(hdr.dim[2]);
  this->NZ=static_cast<int>(hdr.dim[3]);
  this->NT=static_cast<int>(hdr.dim[4]);
  if (this->NT<1) {this->NT=1; hdr.dim[4]=1;}
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;

  //cout << "Size of " << ImageName << ": " << this->NX << " " <<  this->NY << " " <<  this->NZ << endl;
  
  //Image to world matrix
  if (hdr.sform_code>0){ 
    //METHOD 3 of nifti1.h
    //cout << "Orientation of image " << ImageName << " opened using method 3 (alternative - normal)" << endl;
    this->Image2World[0][0]=hdr.srow_x[0];  this->Image2World[0][1]=hdr.srow_x[1];  this->Image2World[0][2]=hdr.srow_x[2];  this->Image2World[0][3]=hdr.srow_x[3];
    this->Image2World[1][0]=hdr.srow_y[0];  this->Image2World[1][1]=hdr.srow_y[1];  this->Image2World[1][2]=hdr.srow_y[2];  this->Image2World[1][3]=hdr.srow_y[3];
    this->Image2World[2][0]=hdr.srow_z[0];  this->Image2World[2][1]=hdr.srow_z[1];  this->Image2World[2][2]=hdr.srow_z[2];  this->Image2World[2][3]=hdr.srow_z[3];
    this->Image2World[3][0]=0;              this->Image2World[3][1]=0;              this->Image2World[3][2]=0;              this->Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      //cout << "Orientation of image " << ImageName << " opened using method 2 (normal)" << endl;
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      this->Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); this->Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       this->Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     this->Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   this->Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     this->Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       this->Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;                               this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                    this->Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      //set the dixel dimensions to 1 if they are absolutely unrealistic 
      //if ((hdr.pixdim[1]<0.000000001)&&(hdr.pixdim[1]>1000000000)) {hdr.pixdim[1]=1;  hdr.qoffset_x=0;}
      //if ((hdr.pixdim[2]<0.000000001)&&(hdr.pixdim[2]>1000000000)) {hdr.pixdim[2]=1;  hdr.qoffset_y=0;}
      //if ((hdr.pixdim[3]<0.000000001)&&(hdr.pixdim[3]>1000000000)) {hdr.pixdim[3]=1;  hdr.qoffset_z=0;}
      
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      this->Image2World[0][0]=hdr.pixdim[1];  this->Image2World[0][1]=0;             this->Image2World[0][2]=0;             this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=0;              this->Image2World[1][1]=hdr.pixdim[2]; this->Image2World[1][2]=0;             this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=0;              this->Image2World[2][1]=0;             this->Image2World[2][2]=hdr.pixdim[3]; this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;              this->Image2World[3][1]=0;             this->Image2World[3][2]=0;             this->Image2World[3][3]=1;
    }
  }

//  cout << "Image to world matrix of " <<  ImageName << ":" << endl;
//  for (i=0;i<4;i++){
//    for (j=0;j<4;j++){
//      cout << this->Image2World[i][j] << " ";
//    }
//    cout << endl;
//  }
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
//  cout << "World to image matrix of " <<  ImageName << ":" << endl;
//  for (i=0;i<4;i++){
//    for (j=0;j<4;j++){
//      cout << this->World2Image[i][j] << " ";
//    }
//    cout << endl;
//  }
  
  //print a little header information
  //fprintf(stderr, "\n%s header information:",ImageName);
  //fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr.dim[1],hdr.dim[2],hdr.dim[3],hdr.dim[4]);
  //fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr.datatype,hdr.bitpix);
  //fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr.scl_slope,hdr.scl_inter);
  //fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr.vox_offset));
  //fprintf(stderr, "\n");
  
  //message about the type of grey levels
  //cout << endl;
  //cout << ImageName << " contains";
  //if (hdr.datatype==2) cout << " unsigned char ";
  //if (hdr.datatype==4) cout << " signed short ";
  //if (hdr.datatype==8) cout << " signed int ";
  //if (hdr.datatype==16) cout << " float ";
  //if (hdr.datatype==64) cout << " double ";
  //if (hdr.datatype==256) cout << " signed char ";
  //if (hdr.datatype==512) cout << " unsigned short ";
  //if (hdr.datatype==768) cout << " unsigned int ";
  //cout << "pixels and has a resolution of " << hdr.dim[1] << "*"  << hdr.dim[2] << "*"  << hdr.dim[3] << "*"  << hdr.dim[4] << " voxels" << endl;
  //cout << endl;
  
  //2) read the image
  
  //allocate the memory for the image
  this->ScalField = new float [this->NXtYtZ*this->NT];

  // jump to data offset
  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);
  
  //test the multiplicatory value of hdr
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << ImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }

  //load the image
  if (hdr.datatype==2){ //2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data2, sizeof(unsigned char), 1, fp);
      this->P(static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==4){ //4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data4, sizeof(signed short), 1, fp);
      this->P(static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==8){  //8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data8, sizeof(signed int), 1, fp);
      this->P(static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==16){  //16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data16, sizeof(float), 1, fp);
      this->P(static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==64){  //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data64, sizeof(double), 1, fp);
      this->P(static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==256){  //256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data256, sizeof(signed char), 1, fp);
      this->P(static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==512){  //512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data512, sizeof(unsigned short), 1, fp);
      this->P(static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else  if (hdr.datatype==768){  //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
      ret = fread(&data768, sizeof(unsigned int), 1, fp);
      this->P(static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter),x,y,z,t);
    }
  }
  else{
    cout << "I can't open an image with grey levels of this type" << endl;
  }
  
  fclose(fp);

}


#define HEADER_N_VOX_OFFSET(x,y,z,t,s) (long)(hdr.vox_offset+((x)+(y*OrigNBX)+(z*OrigNBY*OrigNBX)+(t*OrigNBZ*OrigNBY*OrigNBX))*sizeof(s))

#define HEADER_N_VOX_OFFSET2(x,y,z,t,sizeof_s) (long)(hdr.vox_offset+((x)+(y*OrigNBX)+(z*OrigNBY*OrigNBX)+(t*OrigNBZ*OrigNBY*OrigNBX))*sizeof_s)

///read a scalar field in a ROI (in a nifti image)
/// -> advanced memory managment for large images: only allocate memory for the ROI
/// -> Inputs are different than in Read_ROI_Given_ImageToWorld_and_Size:
///     We give here the min and max {X,Y,Z,T} in the input image defining the outputed ROI
void ScalarField::Read_only_ROI(char * ImageName,int xMin,int xMax,int yMin,int yMax,int zMin,int zMax,int tMin,int tMax){
  nifti_1_header hdr;
  int tempInt;
  FILE *fp;
  int ret,i;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z,t;
  float a,b,c,d,qfac;
  int OrigNBX,OrigNBY,OrigNBZ,OrigNBT;
  float RealOriginX,RealOriginY,RealOriginZ;
  float VoxOriginX,VoxOriginY,VoxOriginZ;
  int sizeof_currentType;

  
  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) open the file and read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  
      //to open some multichannel 3D images which have the channels in the 5th dimension instead of the 4th dimension
  if ((hdr.dim[4]==1)&&(hdr.dim[5]>1)){
    hdr.dim[4]=hdr.dim[5];
    hdr.dim[5]=1;
  }

  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=hdr.dim[1])||(this->NY!=hdr.dim[2])||(this->NZ!=hdr.dim[3])||(this->NT!=hdr.dim[4]))
      cout << "WARNING: THE SIZE OF A NON-NULL SCALAR FIELD IS CHANGED\n";
  
  
  //fill the parameters of the class
  if (hdr.dim[4]<1) {hdr.dim[4]=1;}
  
  if (xMax<xMin) {tempInt=xMin; xMin=xMax; xMax=tempInt;}
  if (yMax<yMin) {tempInt=yMin; yMin=yMax; yMax=tempInt;}
  if (zMax<zMin) {tempInt=zMin; zMin=zMax; zMax=tempInt;}
  if (tMax<tMin) {tempInt=tMin; tMin=tMax; tMax=tempInt;}
  
  if (xMin<0) xMin=0;
  if (yMin<0) yMin=0;
  if (zMin<0) zMin=0;
  if (tMin<0) tMin=0;
   
  if (xMax>(static_cast<int>(hdr.dim[1]))) xMax=(static_cast<int>(hdr.dim[1]));
  if (yMax>(static_cast<int>(hdr.dim[2]))) yMax=(static_cast<int>(hdr.dim[2]));
  if (zMax>(static_cast<int>(hdr.dim[3]))) zMax=(static_cast<int>(hdr.dim[3]));
  if (tMax>(static_cast<int>(hdr.dim[4]))) tMax=(static_cast<int>(hdr.dim[4]));

  this->NX=xMax-xMin;
  this->NY=yMax-yMin;
  this->NZ=zMax-zMin;
  this->NT=tMax-tMin;
  
  OrigNBX=static_cast<int>(hdr.dim[1]);
  OrigNBY=static_cast<int>(hdr.dim[2]);
  OrigNBZ=static_cast<int>(hdr.dim[3]);
  OrigNBT=static_cast<int>(hdr.dim[4]);

  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;

  //Image to world matrix
  if (hdr.sform_code>0){ 
    this->Image2World[0][0]=hdr.srow_x[0];  this->Image2World[0][1]=hdr.srow_x[1];  this->Image2World[0][2]=hdr.srow_x[2];  this->Image2World[0][3]=hdr.srow_x[3];
    this->Image2World[1][0]=hdr.srow_y[0];  this->Image2World[1][1]=hdr.srow_y[1];  this->Image2World[1][2]=hdr.srow_y[2];  this->Image2World[1][3]=hdr.srow_y[3];
    this->Image2World[2][0]=hdr.srow_z[0];  this->Image2World[2][1]=hdr.srow_z[1];  this->Image2World[2][2]=hdr.srow_z[2];  this->Image2World[2][3]=hdr.srow_z[3];
    this->Image2World[3][0]=0;              this->Image2World[3][1]=0;              this->Image2World[3][2]=0;              this->Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      this->Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); this->Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       this->Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     this->Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   this->Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     this->Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       this->Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;                               this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                    this->Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      this->Image2World[0][0]=hdr.pixdim[1];  this->Image2World[0][1]=0;             this->Image2World[0][2]=0;             this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=0;              this->Image2World[1][1]=hdr.pixdim[2]; this->Image2World[1][2]=0;             this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=0;              this->Image2World[2][1]=0;             this->Image2World[2][2]=hdr.pixdim[3]; this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;              this->Image2World[3][1]=0;             this->Image2World[3][2]=0;             this->Image2World[3][3]=1;
    }
  }

  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //message about the type of grey levels
  cout << endl;
  cout << ImageName << " contains";
  if (hdr.datatype==2) {cout << " unsigned char ";         sizeof_currentType=static_cast<int>(sizeof(unsigned char));}
  else if (hdr.datatype==4) {cout << " signed short ";     sizeof_currentType=static_cast<int>(sizeof(signed short));}
  else if (hdr.datatype==8) {cout << " signed int ";       sizeof_currentType=static_cast<int>(sizeof(signed int));}
  else if (hdr.datatype==16) {cout << " float ";           sizeof_currentType=static_cast<int>(sizeof(float));}
  else if (hdr.datatype==64) {cout << " double ";          sizeof_currentType=static_cast<int>(sizeof(double));}
  else if (hdr.datatype==256) {cout << " signed char ";    sizeof_currentType=static_cast<int>(sizeof(signed char));}
  else if (hdr.datatype==512) {cout << " unsigned short "; sizeof_currentType=static_cast<int>(sizeof(unsigned short));}
  else if (hdr.datatype==768) {cout << " unsigned int ";   sizeof_currentType=static_cast<int>(sizeof(unsigned int));}
  else{
    cout << endl;
    cout << "I can't open an image with grey levels of this type" << endl;
    exit(0);
  }
  cout << "pixels and has a resolution of " << OrigNBX << "*"  << OrigNBY << "*"  << OrigNBZ << "*"  << OrigNBT << " voxels" << endl;
  cout << "The ROI has a resolution " << this->NX << "*" << this->NY << "*" << this->NZ << "*" << this->NT << endl;
  cout << endl;
  
  //test the multiplicatory value of hdr
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    //cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << ImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }

  
  //2) read the image -- OLD MEMORY EFFICIENT VERSION (WHICH SEEMS TO CONTAIN A BUG)
  /*
  //allocate the memory for the image
  this->ScalField = new float [this->NXtYtZ*this->NT];

  // jump to data offset
  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);

  
  //load the image  (MODIFIED FILE BUT STILL SEEMS TO CONTAIN A BUG)
  t=tMin;
  while (t<tMax){
    z=zMin;
    while (z<zMax){
      y=yMin;
      while (y<yMax){
        x=xMin;
        ret = fseek(fp, HEADER_N_VOX_OFFSET2(x,y,z,t,sizeof_currentType), SEEK_SET);
        while (x<xMax){
          if (hdr.datatype==16){//16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data16, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==2){//2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data2, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==8){//8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data8, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==4){//4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data4, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==64){ //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data64, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==256){//256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data256, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==512){//512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data512, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          else if (hdr.datatype==768){ //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
            ret = fread(&data768, sizeof_currentType, 1, fp);
            this->P(static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter),x-xMin,y-yMin,z-zMin,t-tMin);
          }
          x++;
        }
        y++;
      }
      z++;
    }
    t++;
  }
   
  fclose(fp);
  */
  
  //2) read the image -- HIGHLY MEMORY CONSUMING (BUT BUG FREE) VERSION
  
  cerr << "The original algorithm which only allocates memory for the ROI seems to contain a bug. " << endl;
  cerr << "An alternative version consuming more memory but bug free is used instead." << endl;

  //allocate the memory for the image
  this->ScalField = new float [this->NXtYtZ*this->NT];

  //read the orginal image
  ScalarField tmpImag;
  tmpImag.Read(ImageName);
  
  //copy the ROI
  for (t=0;t<this->NT;t++) for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++) 
    this->P(tmpImag.G(x+xMin,y+yMin,z+zMin,t+tMin),x,y,z,t);
  
   //3) update the image to world matrix
  VoxOriginX=static_cast<float>(xMin);
  VoxOriginY=static_cast<float>(yMin);
  VoxOriginZ=static_cast<float>(zMin);
  
  
  RealOriginX=this->Image2World[0][0]*VoxOriginX+this->Image2World[0][1]*VoxOriginY+this->Image2World[0][2]*VoxOriginZ+this->Image2World[0][3];
  RealOriginY=this->Image2World[1][0]*VoxOriginX+this->Image2World[1][1]*VoxOriginY+this->Image2World[1][2]*VoxOriginZ+this->Image2World[1][3];
  RealOriginZ=this->Image2World[2][0]*VoxOriginX+this->Image2World[2][1]*VoxOriginY+this->Image2World[2][2]*VoxOriginZ+this->Image2World[2][3];
  
  this->Image2World[0][3]=RealOriginX;
  this->Image2World[1][3]=RealOriginY;
  this->Image2World[2][3]=RealOriginZ;

  invert_4t4quaternion(this->Image2World,this->World2Image);
  
}




///read a scalar field in a ROI (in a nifti image)
/// -> advanced memory managment for large images: only allocate memory for the ROI
/// -> Inputs are different than in Read_only_ROI:
///     We give here the min and max {X,Y,Z,T} in the input image defining the outputed ROI
void ScalarField::Read_ROI_Given_ImageToWorld_and_Size(float ROI_Image2World[4][4],int NBX,int NBY,int NBZ,char * RefImageName){
  nifti_1_header hdr;
  FILE *fp;
  int ret;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z;
  float a,b,c,d,qfac;
  int OrigNBX,OrigNBY,OrigNBZ;
  int sizeof_currentType;
  float RefIma_Image2World[4][4];
  float RefIma_World2Image[4][4];
  float ImaCoord_ROI2RefIma[4][4];
  float Orig_coord[4];
  float ROI_coord[4];
  int x_ref,y_ref,z_ref;
  
  
  //1) generate a void image in this and set its image 2 world properties
  
  this->CreateVoidField(NBX,NBY,NBZ);
  
  this->Image2World[0][0]=ROI_Image2World[0][0]; this->Image2World[0][1]=ROI_Image2World[0][1]; this->Image2World[0][2]=ROI_Image2World[0][2]; this->Image2World[0][3]=ROI_Image2World[0][3];
  this->Image2World[1][0]=ROI_Image2World[1][0]; this->Image2World[1][1]=ROI_Image2World[1][1]; this->Image2World[1][2]=ROI_Image2World[1][2]; this->Image2World[1][3]=ROI_Image2World[1][3];
  this->Image2World[2][0]=ROI_Image2World[2][0]; this->Image2World[2][1]=ROI_Image2World[2][1]; this->Image2World[2][2]=ROI_Image2World[2][2]; this->Image2World[2][3]=ROI_Image2World[2][3];
  this->Image2World[3][0]=ROI_Image2World[3][0]; this->Image2World[3][1]=ROI_Image2World[3][1]; this->Image2World[3][2]=ROI_Image2World[3][2]; this->Image2World[3][3]=ROI_Image2World[3][3];
  
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  cout << this->Image2World[0][0] << " " << this->Image2World[0][1] <<  " " << this->Image2World[0][2] <<  " " << this->Image2World[0][3] << endl;
  cout << this->Image2World[1][0] << " " << this->Image2World[1][1] <<  " " << this->Image2World[1][2] <<  " " << this->Image2World[1][3] << endl;
  cout << this->Image2World[2][0] << " " << this->Image2World[2][1] <<  " " << this->Image2World[2][2] <<  " " << this->Image2World[2][3] << endl;
  cout << this->Image2World[3][0] << " " << this->Image2World[3][1] <<  " " << this->Image2World[3][2] <<  " " << this->Image2World[3][3] << endl;


  //2) open the reference file and read its properties

  //2.1) open the file, read the header and make some security treatments on the header
  fp = fopen(RefImageName,"rb");
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  if (hdr.dim[4]<1) {hdr.dim[4]=1;}
  
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    //cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << RefImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }

  
  //2.2) Image to world matrix of the reference image
  if (hdr.sform_code>0){ 
    RefIma_Image2World[0][0]=hdr.srow_x[0];  RefIma_Image2World[0][1]=hdr.srow_x[1];  RefIma_Image2World[0][2]=hdr.srow_x[2];  RefIma_Image2World[0][3]=hdr.srow_x[3];
    RefIma_Image2World[1][0]=hdr.srow_y[0];  RefIma_Image2World[1][1]=hdr.srow_y[1];  RefIma_Image2World[1][2]=hdr.srow_y[2];  RefIma_Image2World[1][3]=hdr.srow_y[3];
    RefIma_Image2World[2][0]=hdr.srow_z[0];  RefIma_Image2World[2][1]=hdr.srow_z[1];  RefIma_Image2World[2][2]=hdr.srow_z[2];  RefIma_Image2World[2][3]=hdr.srow_z[3];
    RefIma_Image2World[3][0]=0;              RefIma_Image2World[3][1]=0;              RefIma_Image2World[3][2]=0;              RefIma_Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      RefIma_Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); RefIma_Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       RefIma_Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     RefIma_Image2World[0][3]=hdr.qoffset_x;
      RefIma_Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     RefIma_Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   RefIma_Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     RefIma_Image2World[1][3]=hdr.qoffset_y;
      RefIma_Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     RefIma_Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       RefIma_Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); RefIma_Image2World[2][3]=hdr.qoffset_z;
      RefIma_Image2World[3][0]=0;                               RefIma_Image2World[3][1]=0;                                 RefIma_Image2World[3][2]=0;                                    RefIma_Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      
      //put the voxel dimensions in image to world
      cout << "Orientations of " << RefImageName << " were basically estimated..." << endl;
      RefIma_Image2World[0][0]=hdr.pixdim[1];  RefIma_Image2World[0][1]=0;             RefIma_Image2World[0][2]=0;             RefIma_Image2World[0][3]=hdr.qoffset_x;
      RefIma_Image2World[1][0]=0;              RefIma_Image2World[1][1]=hdr.pixdim[2]; RefIma_Image2World[1][2]=0;             RefIma_Image2World[1][3]=hdr.qoffset_y;
      RefIma_Image2World[2][0]=0;              RefIma_Image2World[2][1]=0;             RefIma_Image2World[2][2]=hdr.pixdim[3]; RefIma_Image2World[2][3]=hdr.qoffset_z;
      RefIma_Image2World[3][0]=0;              RefIma_Image2World[3][1]=0;             RefIma_Image2World[3][2]=0;             RefIma_Image2World[3][3]=1;
    }
  }

  invert_4t4quaternion(RefIma_Image2World,RefIma_World2Image);
  
  //2.3) Size of the reference image
  OrigNBX=static_cast<int>(hdr.dim[1]);
  OrigNBY=static_cast<int>(hdr.dim[2]);
  OrigNBZ=static_cast<int>(hdr.dim[3]);

  
  //2.4) compute the matrix to transform the voxel coordinates from the output ROI to the input reference image
  mult_quat4t4mat_quat4t4mat(RefIma_World2Image,this->Image2World,ImaCoord_ROI2RefIma);
  
  cout << "ROI to Ref. image in image coordinates:" << endl;
  cout << ImaCoord_ROI2RefIma[0][0] << " " << ImaCoord_ROI2RefIma[0][1] <<  " " << ImaCoord_ROI2RefIma[0][2] <<  " " << ImaCoord_ROI2RefIma[0][3] << endl;
  cout << ImaCoord_ROI2RefIma[1][0] << " " << ImaCoord_ROI2RefIma[1][1] <<  " " << ImaCoord_ROI2RefIma[1][2] <<  " " << ImaCoord_ROI2RefIma[1][3] << endl;
  cout << ImaCoord_ROI2RefIma[2][0] << " " << ImaCoord_ROI2RefIma[2][1] <<  " " << ImaCoord_ROI2RefIma[2][2] <<  " " << ImaCoord_ROI2RefIma[2][3] << endl;
  cout << ImaCoord_ROI2RefIma[3][0] << " " << ImaCoord_ROI2RefIma[3][1] <<  " " << ImaCoord_ROI2RefIma[3][2] <<  " " << ImaCoord_ROI2RefIma[3][3] << endl;

  
  //2.5) define sizeof(type in the ref image)
  if (hdr.datatype==2) {        sizeof_currentType=static_cast<int>(sizeof(unsigned char));}
  else if (hdr.datatype==4) {   sizeof_currentType=static_cast<int>(sizeof(signed short));}
  else if (hdr.datatype==8) {   sizeof_currentType=static_cast<int>(sizeof(signed int));}
  else if (hdr.datatype==16) {  sizeof_currentType=static_cast<int>(sizeof(float));}
  else if (hdr.datatype==64) {  sizeof_currentType=static_cast<int>(sizeof(double));}
  else if (hdr.datatype==256) { sizeof_currentType=static_cast<int>(sizeof(signed char));}
  else if (hdr.datatype==512) { sizeof_currentType=static_cast<int>(sizeof(unsigned short));}
  else if (hdr.datatype==768) { sizeof_currentType=static_cast<int>(sizeof(unsigned int));}
  else{
    cout << endl;
    cout << "I can't open an image with grey levels of this type" << endl;
    exit(0);
  }
  
  
  //3) read the image
  ROI_coord[3]=1;
  for (z=0;z<this->NZ;z++) for (y=0;y<this->NY;y++) for (x=0;x<this->NX;x++) {
    //3.1) go to where the point is in the reference image
    ROI_coord[0]=x;
    ROI_coord[1]=y;
    ROI_coord[2]=z;
    
    mult_4t4mat_4vec(ImaCoord_ROI2RefIma,ROI_coord,Orig_coord);
  
    x_ref=static_cast<int>(Orig_coord[0]+0.5);
    y_ref=static_cast<int>(Orig_coord[1]+0.5);
    z_ref=static_cast<int>(Orig_coord[2]+0.5);
    
    if (x_ref<0) x_ref=0;
    if (y_ref<0) y_ref=0;
    if (z_ref<0) z_ref=0;
    if (x_ref>=OrigNBX) x_ref=OrigNBX-1;
    if (y_ref>=OrigNBY) y_ref=OrigNBY-1;
    if (z_ref>=OrigNBZ) z_ref=OrigNBZ-1;
  
    ret = fseek(fp, HEADER_N_VOX_OFFSET2(x_ref,y_ref,z_ref,0,sizeof_currentType), SEEK_SET);
    
    //3.2) pick up the grey level
    if (hdr.datatype==16){//16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data16, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==2){//2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data2, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==8){//8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data8, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==4){//4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data4, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==64){ //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data64, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==256){//256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data256, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==512){//512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data512, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
    else if (hdr.datatype==768){ //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
      ret = fread(&data768, sizeof_currentType, 1, fp);
      this->P(static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter),x,y,z);
    }
  }
  


  fclose(fp);

}



///read a scalar field in a ROI (in a nifti image)
/// -> advanced memory managment for large images: undersample the image direcly by averaging blocks of voxels of size BlockS*BlockS*BlockS
///Use 'Read_and_Undersample' or 'Read_and_Interpolate' to have finer resamplings requiring more memory
void ScalarField::Read_directly_Undersampled(char * ImageName,int BlockS){
  nifti_1_header hdr;
  int tempInt;
  FILE *fp;
  int ret,i;
  unsigned char data2;      //probably not the
  signed short data4;      //nicest technique 
  signed int data8;       //to open 
  float data16;          // different kinds
  double data64;         // of images
  signed char data256;   // but
  unsigned short data512; // it
  unsigned int data768;    // works
  int x,y,z,t;
  float a,b,c,d,qfac;
  float BlockVolume;
  float TempValue;
  int LocTmpX,LocTmpY,LocTmpZ;
  
  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) open the file and read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  //message if there is already an image in InputImage with another size of the opened one
  if (this->NX!=0)
    if ((this->NX!=hdr.dim[1])||(this->NY!=hdr.dim[2])||(this->NZ!=hdr.dim[3])||(this->NT!=hdr.dim[4]))
      cout << "WARNING: THE SIZE OF A NON-NULL SCALAR FIELD IS CHANGED\n";
  
  
  //fill the parameters of the class
  if (hdr.dim[4]<1) {hdr.dim[4]=1;}
  
  this->NX=static_cast<int>(hdr.dim[1]/BlockS);
  this->NY=static_cast<int>(hdr.dim[2]/BlockS);
  this->NZ=static_cast<int>(hdr.dim[3]/BlockS);
  this->NT=hdr.dim[4];
  
  BlockVolume=static_cast<float>(BlockS*BlockS*BlockS);
  
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;

  //Image to world matrix
  if (hdr.sform_code>0){ 
    this->Image2World[0][0]=hdr.srow_x[0];  this->Image2World[0][1]=hdr.srow_x[1];  this->Image2World[0][2]=hdr.srow_x[2];  this->Image2World[0][3]=hdr.srow_x[3];
    this->Image2World[1][0]=hdr.srow_y[0];  this->Image2World[1][1]=hdr.srow_y[1];  this->Image2World[1][2]=hdr.srow_y[2];  this->Image2World[1][3]=hdr.srow_y[3];
    this->Image2World[2][0]=hdr.srow_z[0];  this->Image2World[2][1]=hdr.srow_z[1];  this->Image2World[2][2]=hdr.srow_z[2];  this->Image2World[2][3]=hdr.srow_z[3];
    this->Image2World[3][0]=0;              this->Image2World[3][1]=0;              this->Image2World[3][2]=0;              this->Image2World[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      this->Image2World[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); this->Image2World[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       this->Image2World[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     this->Image2World[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   this->Image2World[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     this->Image2World[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       this->Image2World[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;                               this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                    this->Image2World[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      this->Image2World[0][0]=hdr.pixdim[1];  this->Image2World[0][1]=0;             this->Image2World[0][2]=0;             this->Image2World[0][3]=hdr.qoffset_x;
      this->Image2World[1][0]=0;              this->Image2World[1][1]=hdr.pixdim[2]; this->Image2World[1][2]=0;             this->Image2World[1][3]=hdr.qoffset_y;
      this->Image2World[2][0]=0;              this->Image2World[2][1]=0;             this->Image2World[2][2]=hdr.pixdim[3]; this->Image2World[2][3]=hdr.qoffset_z;
      this->Image2World[3][0]=0;              this->Image2World[3][1]=0;             this->Image2World[3][2]=0;             this->Image2World[3][3]=1;
    }
  }

  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //message about the type of grey levels
  cout << endl;
  cout << ImageName << " contains";
  if (hdr.datatype==2) cout << " unsigned char ";
  if (hdr.datatype==4) cout << " signed short ";
  if (hdr.datatype==8) cout << " signed int ";
  if (hdr.datatype==16) cout << " float ";
  if (hdr.datatype==64) cout << " double ";
  if (hdr.datatype==256) cout << " signed char ";
  if (hdr.datatype==512) cout << " unsigned short ";
  if (hdr.datatype==768) cout << " unsigned int ";
  cout << "pixels and has a resolution of " << hdr.dim[1] << "*"  << hdr.dim[2] << "*"  << hdr.dim[3] << "*"  << hdr.dim[4] << " voxels" << endl;
  cout << "The undersampled image has a resolution " << this->NX << "*" << this->NY << "*" << this->NZ << "*" << this->NT << endl;
  cout << endl;
  
  //2) read the image
  
  //allocate the memory for the image and fill it with zeros
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(0,x,y,z,t);
  

  // jump to data offset
  ret = fseek(fp, (long)(hdr.vox_offset), SEEK_SET);

  //test the multiplicatory value of hdr
  if (hdr.scl_slope==0){
    hdr.scl_slope=1;
    cout << "Warning the multiplicatory factor of the grey levels (scl_slope) in " << ImageName << " is equal to 0. We set it to 1 in the opened image." << endl;
  }
  
  //load the image
  if (hdr.datatype==2){ //2  ->  unsigned char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data2, sizeof(unsigned char), 1, fp);
      
      TempValue=static_cast<float>((data2 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==4){ //4  ->  signed short +++++++++++++++++++++++++++++++++++++++++++++++++    
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data4, sizeof(signed short), 1, fp);
      
      TempValue=static_cast<float>((data4 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==8){  //8  ->  signed int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data8, sizeof(signed int), 1, fp);
      
      TempValue=static_cast<float>((data8 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==16){  //16  ->  float +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data16, sizeof(float), 1, fp);
      
      TempValue=static_cast<float>((data16 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==64){  //64  ->  double +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data64, sizeof(double), 1, fp);
      
      TempValue=static_cast<float>((data64 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==256){  //256  ->  signed char +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data256, sizeof(signed char), 1, fp);
      
      TempValue=static_cast<float>((data256 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==512){  //512  ->  unsigned short +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data512, sizeof(unsigned short), 1, fp);
      
      TempValue=static_cast<float>((data512 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else  if (hdr.datatype==768){  //768  ->  unsigned int +++++++++++++++++++++++++++++++++++++++++++++++++
    for(t=0;t<hdr.dim[4];t++) for(z=0;z<hdr.dim[3];z++) for(y=0;y<hdr.dim[2];y++) for(x=0;x<hdr.dim[1];x++){
      ret = fread(&data768, sizeof(unsigned int), 1, fp);
      
      TempValue=static_cast<float>((data768 * hdr.scl_slope) + hdr.scl_inter)/BlockVolume;
      LocTmpX=x/BlockS; LocTmpY=y/BlockS; LocTmpZ=z/BlockS;
      
      if ((LocTmpX<this->NX)&&(LocTmpY<this->NY)&&(LocTmpZ<this->NZ)) 
        this->Add(TempValue,LocTmpX,LocTmpY,LocTmpZ,t);
    }
  }
  else{
    cout << "I can't open an image with grey levels of this type" << endl;
  }
  
  fclose(fp);
  
  //3) update the image to world matrix
  this->Image2World[0][0]*=BlockS; this->Image2World[0][1]*=BlockS; this->Image2World[0][2]*=BlockS; 
  this->Image2World[1][0]*=BlockS; this->Image2World[1][1]*=BlockS; this->Image2World[1][2]*=BlockS; 
  this->Image2World[2][0]*=BlockS; this->Image2World[2][1]*=BlockS; this->Image2World[2][2]*=BlockS; 

  invert_4t4quaternion(this->Image2World,this->World2Image);
  
}


///read a scalar field and expend its domain
void ScalarField::ReadAndExpend(char * ImageName,int addX1,int addX2,int addY1,int addY2,int addZ1,int addZ2){
	ScalarField OrigSF;
	int x,y,z,t;
	float x2,y2,z2;
  float minImag;
  
  
	//read the scalar field at the original format
	OrigSF.Read(ImageName);
	
	//fill the parameters of the class and allocate the memory for the scalar field
	this->NX=OrigSF.NX+addX1+addX2;
	this->NY=OrigSF.NY+addY1+addY2;
	this->NZ=OrigSF.NZ+addZ1+addZ2;
	this->NT=OrigSF.NT;
	this->NXtY=this->NX*this->NY;
	this->NXtYtZ=this->NXtY*this->NZ;
	this->ScalField = new float [this->NXtYtZ*this->NT];
	
  
  this->Image2World[0][0]=OrigSF.Image2World[0][0];  this->Image2World[0][1]=OrigSF.Image2World[0][1];  this->Image2World[0][2]=OrigSF.Image2World[0][2];  this->Image2World[0][3]=OrigSF.Image2World[0][3]-addX1*OrigSF.Image2World[0][0];
  this->Image2World[1][0]=OrigSF.Image2World[1][0];  this->Image2World[1][1]=OrigSF.Image2World[1][1];  this->Image2World[1][2]=OrigSF.Image2World[1][2];  this->Image2World[1][3]=OrigSF.Image2World[1][3]-addY1*OrigSF.Image2World[1][1];
  this->Image2World[2][0]=OrigSF.Image2World[3][0];  this->Image2World[2][1]=OrigSF.Image2World[2][1];  this->Image2World[2][2]=OrigSF.Image2World[2][2];  this->Image2World[2][3]=OrigSF.Image2World[2][3]-addZ1*OrigSF.Image2World[2][2];
  this->Image2World[3][0]=0;                         this->Image2World[3][1]=0;                         this->Image2World[3][2]=0;                         this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
	//fill the image image...
  
  //...init
  minImag=OrigSF.G(0,0,0,0);
  
	for(t=0;t<OrigSF.NT;t++) for(z=0;z<OrigSF.NZ;z++) for(y=0;y<OrigSF.NY;y++) for(x=0;x<OrigSF.NX;x++) if (minImag>OrigSF.G(x,y,z,t)) minImag=OrigSF.G(x,y,z,t);

  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) this->P(minImag-1,x,y,z,t);
  
  //... original image
	for(t=0;t<OrigSF.NT;t++) for(z=0;z<OrigSF.NZ;z++) for(y=0;y<OrigSF.NY;y++) for(x=0;x<OrigSF.NX;x++)
    if ((x+addX1>=0)&&(x+addX1<this->NX)&& (y+addY1>=0)&&(y+addY1<this->NY)&& (z+addZ1>=0)&&(z+addZ1<this->NZ))
      this->P(OrigSF.G(x,y,z,t),x+addX1,y+addY1,z+addZ1,t);

  //... image extension
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) if (this->G(x,y,z,t)<minImag){
    x2=x;
    if (x<addX1) x2=addX1;
    if (x>=OrigSF.NX+addX1) x2=OrigSF.NX+addX1-1;
    y2=y;
    if (y<addY1) y2=addY1;
    if (y>=OrigSF.NY+addY1) y2=OrigSF.NY+addY1-1;
    z2=z;
    if (z<addZ1) z2=addZ1;
    if (z>=OrigSF.NZ+addZ1) z2=OrigSF.NZ+addZ1-1;
    
    this->P(this->G(x2,y2,z2,t),x,y,z,t);
  
  }
}



///read a scalar field and perform linear interpolation to give it a specific size
void ScalarField::Read_and_Interpolate(char * ImageName,int NBX,int NBY,int NBZ){
	ScalarField OrigSF;
	int x,y,z,t;
	float x2,y2,z2;
	float factorX,factorY,factorZ;
  
  
	//read the scalar field at the original format
	OrigSF.Read(ImageName);
	
	//fill the parameters of the class and allocate the memory for the scalar field
	this->NX=NBX;
	this->NY=NBY;
	this->NZ=NBZ;
	this->NT=OrigSF.NT;
	this->NXtY=this->NX*this->NY;
	this->NXtYtZ=this->NXtY*this->NZ;
	this->ScalField = new float [this->NXtYtZ*this->NT];
	
  factorX=(static_cast<float>(OrigSF.NX)/static_cast<float>(this->NX));
  factorY=(static_cast<float>(OrigSF.NY)/static_cast<float>(this->NY));
  factorZ=(static_cast<float>(OrigSF.NZ)/static_cast<float>(this->NZ));
  
  this->Image2World[0][0]=OrigSF.Image2World[0][0]*factorX;  this->Image2World[0][1]=OrigSF.Image2World[0][1]*factorX;  this->Image2World[0][2]=OrigSF.Image2World[0][2]*factorX;  this->Image2World[0][3]=OrigSF.Image2World[0][3];
  this->Image2World[1][0]=OrigSF.Image2World[1][0]*factorY;  this->Image2World[1][1]=OrigSF.Image2World[1][1]*factorY;  this->Image2World[1][2]=OrigSF.Image2World[1][2]*factorY;  this->Image2World[1][3]=OrigSF.Image2World[1][3];
  this->Image2World[2][0]=OrigSF.Image2World[3][0]*factorZ;  this->Image2World[2][1]=OrigSF.Image2World[2][1]*factorZ;  this->Image2World[2][2]=OrigSF.Image2World[2][2]*factorZ;  this->Image2World[2][3]=OrigSF.Image2World[2][3];
  this->Image2World[3][0]=0;                                 this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                 this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
	//interpolate the original image
	for(t=0;t<this->NT;t++){
		for(z=0;z<this->NZ;z++){ 
			z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
			for(y=0;y<this->NY;y++){ 
				y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
				for(x=0;x<this->NX;x++){
					x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
					this->P(OrigSF.G(x2,y2,z2,t),x,y,z,t);
				}
			}
		}
	}
}

///read a scalar field and perform undersampling by a factor 'factor'
void ScalarField::Read_and_Undersample(char * ImageName,float factor, float NN){
  ScalarField OrigSF;
  int x,y,z,t;
  float x2,y2,z2;
  float factorX,factorY,factorZ;
  
  
  //read the scalar field at the original format
  OrigSF.Read(ImageName);
  
  //fill the parameters of the class and allocate the memory for the scalar field
  this->NX=static_cast<int>(OrigSF.NX/factor);
  this->NY=static_cast<int>(OrigSF.NY/factor);
  if (OrigSF.NZ>1)
    this->NZ=static_cast<int>(OrigSF.NZ/factor);
  else
    this->NZ=static_cast<int>(OrigSF.NZ/1);
  this->NT=OrigSF.NT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  factorX=(static_cast<float>(OrigSF.NX)/static_cast<float>(this->NX));
  factorY=(static_cast<float>(OrigSF.NY)/static_cast<float>(this->NY));
  factorZ=(static_cast<float>(OrigSF.NZ)/static_cast<float>(this->NZ));
  
  this->Image2World[0][0]=OrigSF.Image2World[0][0]*factorX;  this->Image2World[0][1]=OrigSF.Image2World[0][1]*factorX;  this->Image2World[0][2]=OrigSF.Image2World[0][2]*factorX;  this->Image2World[0][3]=OrigSF.Image2World[0][3];
  this->Image2World[1][0]=OrigSF.Image2World[1][0]*factorY;  this->Image2World[1][1]=OrigSF.Image2World[1][1]*factorY;  this->Image2World[1][2]=OrigSF.Image2World[1][2]*factorY;  this->Image2World[1][3]=OrigSF.Image2World[1][3];
  this->Image2World[2][0]=OrigSF.Image2World[3][0]*factorZ;  this->Image2World[2][1]=OrigSF.Image2World[2][1]*factorZ;  this->Image2World[2][2]=OrigSF.Image2World[2][2]*factorZ;  this->Image2World[2][3]=OrigSF.Image2World[2][3];
  this->Image2World[3][0]=0;                                 this->Image2World[3][1]=0;                                 this->Image2World[3][2]=0;                                 this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //interpolate the original image
  if (NN==0){  //triliear interpolation
    for(t=0;t<this->NT;t++){
      for(z=0;z<this->NZ;z++){ 
        z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
        for(y=0;y<this->NY;y++){ 
          y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
          for(x=0;x<this->NX;x++){
            x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
            this->P(OrigSF.G(x2,y2,z2,t),x,y,z,t);
          }
        }
      }
    }
  }
  else{   //nearest neighbor
    for(t=0;t<this->NT;t++){
      for(z=0;z<this->NZ;z++){ 
        z2=floor(0.5+static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1));
        for(y=0;y<this->NY;y++){ 
          y2=floor(0.5+static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1));
          for(x=0;x<this->NX;x++){
            x2=floor(0.5+static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1));
            this->P(OrigSF.G(x2,y2,z2,t),x,y,z,t);
          }
        }
      }
    }
  }
}

///create a void scalar field. All values are initialized to 'cste' which is null by default. No message is printed if Verbose!=1.
void ScalarField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT,float cste,int Verbose){
  int x,y,z,t;
  
  //message if there is already an image in InputImage with another size of the opened one
  if (Verbose==1) if (this->NX!=0)
    if ((this->NX!=NBX)||(this->NY!=NBY)||(this->NZ!=NBZ)||(this->NT!=NBT))
      cout << "WARNING: THE SIZE OF A NON-NULL SCALAR FIELD IS CHANGED\n";
  
  //image size
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=NBT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //allocate memory to cast (and eventually transform) the original template and target images
  //    -->  ScalarField[ptSF(x,y,z)]= gray level at (x,y,z)
  this->ScalField = new float [this->NXtYtZ*this->NT];
  
  //set all entries of the field at 0.
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(cste,x,y,z,t);
}


///Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
void ScalarField::SlashFieldSize(int verbative){
	
	//image size
	this->NX=1;
	this->NY=1;
	this->NZ=1;
	this->NT=1;
	this->NXtY=this->NX*this->NY;
	this->NXtYtZ=this->NXtY*this->NZ;
	
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  if (verbative==1) cout << "Slash some memory" << endl;
  
	//re-allocate memory
  if (this->NX!=0) delete this->ScalField;
  this->ScalField = new float [1];
}



///write a scalar field in a nifti image
void ScalarField::Write(char * OutputImageName){
  int i;
  int x,y,z,t;
  nifti_1_header hdr;
  nifti1_extender pad={0,0,0,0};
  FILE *fp;
  int ret;
  float *data=NULL;
  
  
  //1) create the header
  memset((void *)&hdr,0, sizeof(hdr));
  hdr.sizeof_hdr = MIN_HEADER_SIZE;
  hdr.dim[0] = 4;
  hdr.dim[1] = this->NX;
  hdr.dim[2] = this->NY;
  hdr.dim[3] = this->NZ;
  hdr.dim[4] = this->NT;
  hdr.datatype = NIFTI_TYPE_FLOAT32;
  hdr.bitpix = 32; 
  hdr.qform_code=0; // should ideally be set to 1 but I don't set the values of 'quatern_b', 'quatern_c' and 'quatern_d'
  hdr.pixdim[1] = sqrt(this->Image2World[0][0]*this->Image2World[0][0]+this->Image2World[0][1]*this->Image2World[0][1]+this->Image2World[0][2]*this->Image2World[0][2]);
  hdr.pixdim[2] = sqrt(this->Image2World[1][0]*this->Image2World[1][0]+this->Image2World[1][1]*this->Image2World[1][1]+this->Image2World[1][2]*this->Image2World[1][2]);
  hdr.pixdim[3] = sqrt(this->Image2World[2][0]*this->Image2World[2][0]+this->Image2World[2][1]*this->Image2World[2][1]+this->Image2World[2][2]*this->Image2World[2][2]);
  hdr.qoffset_x=this->Image2World[0][3];
  hdr.qoffset_y=this->Image2World[1][3];
  hdr.qoffset_z=this->Image2World[2][3];
  hdr.pixdim[4] = 1.0;
  hdr.sform_code=1;
  hdr.srow_x[0]=this->Image2World[0][0];  hdr.srow_x[1]=this->Image2World[0][1];  hdr.srow_x[2]=this->Image2World[0][2];  hdr.srow_x[3]=this->Image2World[0][3];
  hdr.srow_y[0]=this->Image2World[1][0];  hdr.srow_y[1]=this->Image2World[1][1];  hdr.srow_y[2]=this->Image2World[1][2];  hdr.srow_y[3]=this->Image2World[1][3];
  hdr.srow_z[0]=this->Image2World[2][0];  hdr.srow_z[1]=this->Image2World[2][1];  hdr.srow_z[2]=this->Image2World[2][2];  hdr.srow_z[3]=this->Image2World[2][3];
  hdr.vox_offset = (float) NII_HEADER_SIZE;
  hdr.scl_inter = 0.0;
  hdr.scl_slope = 1.0;
  hdr.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
  strncpy(hdr.magic, "n+1\0", 4);

  //2) save the image OutputImageName
  
  /********** allocate and fill the buffer  */
  //data = (float *) malloc(sizeof(float) * hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]);
  data = new float [hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4]];
  
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(x,y,z,t);
    i++;
  }
  
  /********** write first 348 bytes of header   */
  fp = fopen(OutputImageName,"wb");
  
  ret = fwrite(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  /********** write extender pad and image data   */
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr.bitpix/8), hdr.dim[1]*hdr.dim[2]*hdr.dim[3]*hdr.dim[4], fp);
  
  
  fclose(fp);
}



///write a scalar field in a nifti image
//The 2nd file is an input file containing the headers 
void ScalarField::Write(char * OutputImageName, char * ImageForHeaderName){
  nifti_1_header hdr_ref;
  FILE *fp_header;
  int i;
  int x,y,z,t;
  nifti1_extender pad={0,0,0,0};
  FILE *fp;
  int ret;
  float *data=NULL;
  
  
  //1) read the header of ImageForHeaderName
  fp_header = fopen(ImageForHeaderName,"rb");
  fread(&hdr_ref, MIN_HEADER_SIZE, 1, fp_header);
  fclose(fp_header);
  

  
  //print a little header information
  //fprintf(stderr, "\n%s header information:",ImageForHeaderName);
  //fprintf(stderr, "\nXYZT dimensions: %d %d %d %d",hdr_ref.dim[1],hdr_ref.dim[2],hdr_ref.dim[3],hdr_ref.dim[4]);
  //fprintf(stderr, "\nDatatype code and bits/pixel: %d %d",hdr_ref.datatype,hdr_ref.bitpix);
  //fprintf(stderr, "\nScaling slope and intercept: %.6f %.6f",hdr_ref.scl_slope,hdr_ref.scl_inter);
  //fprintf(stderr, "\nByte offset to data in datafile: %ld",(long)(hdr_ref.vox_offset));
  //fprintf(stderr, "\n");
  

  //2) save the image OutputImageName
  hdr_ref.dim[0] = 4;
  hdr_ref.dim[1] = this->NX;
  hdr_ref.dim[2] = this->NY;
  hdr_ref.dim[3] = this->NZ;
  hdr_ref.dim[4] = this->NT;
  hdr_ref.datatype = NIFTI_TYPE_FLOAT32;
  hdr_ref.bitpix = 32; 
  hdr_ref.scl_inter = 0.0;
  hdr_ref.scl_slope = 1.0;
  hdr_ref.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;

  
  /********** allocate and fill the buffer  */
  //data = (float *) malloc(sizeof(float) * hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]);
  data = new  float [hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]];
  
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  /********** write first 348 bytes of header   */
  fp = fopen(OutputImageName,"wb");

  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  
  /********** write extender pad and image data   */
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  
  fclose(fp);
}


///grey levels alignment of the grey levels using optimal transport 
///-> optimal transportation minmizes the wassertein 1 distance between the linearly aligned histograms
void ScalarField::GreyLevAlignment(ScalarField * RefImage){
  int i,j,k,t;
  int NbBinsHisto;
  float * CumHisto_LocImage_x_axis;
  float * CumHisto_LocImage_y_axis;
  float * CumHisto_TrgImage_x_axis;
  float * CumHisto_TrgImage_y_axis;
  float * LUT_x_axis;
  float * LUT_y_axis;
  float tmpFl;
  int CorrespBin,CorrespBinPP;
  float alpha;
  
  NbBinsHisto=500;

  //1) allocate and compute the cumulative histograms
  CumHisto_LocImage_x_axis = new float [NbBinsHisto];
  CumHisto_LocImage_y_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_x_axis = new float [NbBinsHisto];
  CumHisto_TrgImage_y_axis = new float [NbBinsHisto];
  LUT_x_axis = new float [NbBinsHisto];
  LUT_y_axis = new float [NbBinsHisto];
  
  this->CptCumulativeHistogram(NbBinsHisto,CumHisto_LocImage_x_axis,CumHisto_LocImage_y_axis,1);   //remark: the log of the histogram is actually used
  RefImage->CptCumulativeHistogram(NbBinsHisto,CumHisto_TrgImage_x_axis,CumHisto_TrgImage_y_axis,1);
  
  //2) compute the look-up table (LUT) from the intensities of the Local image to those of the reference image
  for (i=0;i<NbBinsHisto;i++) LUT_x_axis[i]=CumHisto_LocImage_x_axis[i];
  
  for (i=0;i<NbBinsHisto;i++){
    //... tmpFl is between 0 and 1
    tmpFl=CumHisto_LocImage_y_axis[i]; 
    
    //... find the corresponding bin of tmpFl in the inverse of CumHisto_TrgImage
    CorrespBin=-2;
    for (j=0;j<NbBinsHisto;j++) if (CorrespBin==-2) if (CumHisto_TrgImage_y_axis[j]>tmpFl){
      CorrespBin=j-1;
      }
    
    //... fill the LUT
    if (CorrespBin==-1){
      LUT_y_axis[i]=CumHisto_TrgImage_x_axis[0];
      }
    else  if (CorrespBin==-2){
      LUT_y_axis[i]=CumHisto_TrgImage_x_axis[NbBinsHisto-1];
      }
    else LUT_y_axis[i]=CumHisto_TrgImage_x_axis[CorrespBin];   //nearest neighbor should be OK if there is a sufficient amount of bins
  }
  
  //... final refinement to exactely match the min and max grey levels  (x-axis of the histgrams is the center of the areas where the grey levels are considered)
  LUT_y_axis[0]=CumHisto_TrgImage_x_axis[0]-((CumHisto_TrgImage_x_axis[1]-CumHisto_TrgImage_x_axis[0])/2);
  LUT_y_axis[NbBinsHisto-1]=CumHisto_TrgImage_x_axis[NbBinsHisto-1]+((CumHisto_TrgImage_x_axis[1]-CumHisto_TrgImage_x_axis[0])/2);
  

  //3) resample the image grey levels
  for (t=0;t<this->NT;t++) for (i=0;i<this->NZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++){
    //current grey level
    tmpFl=this->G(k,j,i,t);
    
    //corresponding bin in the LUT
    if (tmpFl<LUT_x_axis[0]) CorrespBin=0;
    else if (tmpFl>LUT_x_axis[NbBinsHisto-1]) CorrespBin=NbBinsHisto-1;
    else CorrespBin=static_cast<int>(static_cast<float>(NbBinsHisto)*(tmpFl-LUT_x_axis[0])/(LUT_x_axis[NbBinsHisto-1]-LUT_x_axis[0]));
    
    if (CorrespBin>NbBinsHisto-1) CorrespBin=NbBinsHisto-1;
    
    //put the new grey level value
    this->P(LUT_y_axis[CorrespBin],k,j,i,t);
    }
}




///Compute the histogram. The histogram is normalized (sum of values equal to 1). 
/// -> Input_BinsNb is the number of bins in the histogram
/// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
void ScalarField::CptHistogram(int Input_BinsNb,float * Output_Histo_x_axis,float * Output_Histo_y_axis,int useLogHisto){
  int i,j,k,t;
  float GLmin,GLmax;
  int locNZ,tmpBin;
  float Input_BinsNbFl,halfDeltaBin,VoxelsNb,tmpFl;
  
  
  //1) initiate
  if (this->NZ<=1) locNZ=1;
  else locNZ=this->NZ;
  
  for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]=0;
  
  Input_BinsNbFl=static_cast<float>(Input_BinsNb);
  
  
  //2) get the min and max grey level
  GLmin=this->G(0,0,0);
  GLmax=this->G(0,0,0);
  
  for (t=0;t<this->NT;t++)  for (i=0;i<locNZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++){
    if (GLmin>this->G(k,j,i,t)) GLmin=this->G(k,j,i,t);
    if (GLmax<this->G(k,j,i,t)) GLmax=this->G(k,j,i,t);
    }
    
  //3) fill the histogram x-axis
  halfDeltaBin=((GLmax-GLmin)/Input_BinsNbFl)/2;
  for (i=0;i<Input_BinsNb;i++) Output_Histo_x_axis[i]=halfDeltaBin+GLmin+(GLmax-GLmin)*(static_cast<float>(i)/Input_BinsNbFl);  
  
  //4) fill the histogram y-axis
  for (t=0;t<this->NT;t++)   for (i=0;i<locNZ;i++) for (j=0;j<this->NY;j++) for (k=0;k<this->NX;k++){
    tmpBin=static_cast<int>(Input_BinsNbFl*(this->G(k,j,i,t)-GLmin)/(GLmax-GLmin));
    if (tmpBin<0) tmpBin=0; //in case of numerical problems
    if (tmpBin>=Input_BinsNb) tmpBin=Input_BinsNb-1;  //in case of numerical problems
    Output_Histo_y_axis[tmpBin]++;
    }
  
  //5) normalize the histogram y-axis and eventually cpt its log before normalizing
  //5.1) standard case
  
  if (useLogHisto==0) {
    VoxelsNb=static_cast<float>(locNZ*this->NY*this->NX);
    for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]/=VoxelsNb;
  }
  
  //5.2) case where the log of the histogram is considered
  if (useLogHisto!=0) {
    tmpFl=0;
    for (i=0;i<Input_BinsNb;i++){
      Output_Histo_y_axis[i]=log(Output_Histo_y_axis[i]+1);
      tmpFl+=Output_Histo_y_axis[i];
    }
    
    for (i=0;i<Input_BinsNb;i++) Output_Histo_y_axis[i]/=tmpFl;
  }
  
  
}



///Compute the cumulative histogram. The cumulative histogram is normalized (last value equal to 1). 
/// -> Input_BinsNb is the number of bins in the histogram
/// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
void ScalarField::CptCumulativeHistogram(int Input_BinsNb,float * Output_CumHisto_x_axis,float * Output_CumHisto_y_axis,int useLogHisto){
  int i;
  float cumulatedvalues;
  
  //compute the histogram
  this->CptHistogram(Input_BinsNb,Output_CumHisto_x_axis,Output_CumHisto_y_axis,useLogHisto);
  
  //compute the cumulative histogram
  cumulatedvalues=0;
  for (i=0;i<Input_BinsNb;i++){
    cumulatedvalues+=Output_CumHisto_y_axis[i];
    Output_CumHisto_y_axis[i]=cumulatedvalues;
  }
}


///get the 'Image 2 World matrix' of an image without loading the image  (load its header only)
void Get_Image2World(char * ImageName,float LocI2W[4][4]){
  nifti_1_header hdr;
  FILE *fp;
  int ret,i;
  int x,y,z,t;
  float a,b,c,d,qfac;

  //0 open the file
  fp = fopen(ImageName,"rb");

  //1) read the header
  ret = fread(&hdr, MIN_HEADER_SIZE, 1, fp);
  
  //2) Define the Image to world matrix
  if (hdr.sform_code>0){ //METHOD 3 of nifti1.h
    //cout << "Orientation of image " << ImageName << " opened using method 3 (alternative - normal)" << endl;
    LocI2W[0][0]=hdr.srow_x[0];  LocI2W[0][1]=hdr.srow_x[1];  LocI2W[0][2]=hdr.srow_x[2];  LocI2W[0][3]=hdr.srow_x[3];
    LocI2W[1][0]=hdr.srow_y[0];  LocI2W[1][1]=hdr.srow_y[1];  LocI2W[1][2]=hdr.srow_y[2];  LocI2W[1][3]=hdr.srow_y[3];
    LocI2W[2][0]=hdr.srow_z[0];  LocI2W[2][1]=hdr.srow_z[1];  LocI2W[2][2]=hdr.srow_z[2];  LocI2W[2][3]=hdr.srow_z[3];
    LocI2W[3][0]=0;              LocI2W[3][1]=0;              LocI2W[3][2]=0;              LocI2W[3][3]=1;
  }
  else {
    if (hdr.qform_code>0){ //METHOD 2 of nifti1.h
      //cout << "Orientation of image " << ImageName << " opened using method 2 (normal)" << endl;
      b=hdr.quatern_b;
      c=hdr.quatern_c;
      d=hdr.quatern_d;
      a=sqrt(1.0-(b*b+c*c+d*d));
      
      qfac=1;
      if (hdr.pixdim[0]==-1.0) qfac=-1;
             
      LocI2W[0][0]=hdr.pixdim[1]*(a*a+b*b-c*c-d*d); LocI2W[0][1]=hdr.pixdim[2]*(2*b*c-2*a*d);       LocI2W[0][2]=qfac*hdr.pixdim[3]*(2*b*d+2*a*c);     LocI2W[0][3]=hdr.qoffset_x;
      LocI2W[1][0]=hdr.pixdim[1]*(2*b*c+2*a*d);     LocI2W[1][1]=hdr.pixdim[2]*(a*a+c*c-b*b-d*d);   LocI2W[1][2]=qfac*hdr.pixdim[3]*(2*c*d-2*a*b);     LocI2W[1][3]=hdr.qoffset_y;
      LocI2W[2][0]=hdr.pixdim[1]*(2*b*d-2*a*c);     LocI2W[2][1]=hdr.pixdim[2]*(2*c*d+2*a*b);       LocI2W[2][2]=qfac*hdr.pixdim[3]*(a*a+d*d-c*c-b*b); LocI2W[2][3]=hdr.qoffset_z;
      LocI2W[3][0]=0;                               LocI2W[3][1]=0;                                 LocI2W[3][2]=0;                                    LocI2W[3][3]=1;
      }
    else{//METHOD 1 of nifti1.h
      //put the voxel dimensions in image to world
      cout << "Orientations of " << ImageName << " were basically estimated..." << endl;
      LocI2W[0][0]=hdr.pixdim[1];  LocI2W[0][1]=0;             LocI2W[0][2]=0;             LocI2W[0][3]=hdr.qoffset_x;
      LocI2W[1][0]=0;              LocI2W[1][1]=hdr.pixdim[2]; LocI2W[1][2]=0;             LocI2W[1][3]=hdr.qoffset_y;
      LocI2W[2][0]=0;              LocI2W[2][1]=0;             LocI2W[2][2]=hdr.pixdim[3]; LocI2W[2][3]=hdr.qoffset_z;
      LocI2W[3][0]=0;              LocI2W[3][1]=0;             LocI2W[3][2]=0;             LocI2W[3][3]=1;
    }
  }
}

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           2:   FUNCTIONS FOR THE CLASS "VectorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///constructor
VectorField::VectorField(void){
	this->NX=0;
	this->NY=0;
	this->NZ=0;
	this->NT=0;
}

///destructor
VectorField::~VectorField(void){
  if ((this->VecField!=NULL)&&(this->NX>0)) delete this->VecField;
}

///put a value
//-> inline function in the .h file

/// add a value
//-> inline function in the .h file

///put the same value at all entries of the vector field
void VectorField::PutToAllVoxels(float cste,int t){
	int i,x,y,z;
	for (i=0;i<3;i++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) { this->P(cste,i,x,y,z,t); }
}


///get a value
//-> inline function in the .h file

///get a value using linear interpolation
float VectorField::G(int IdDirec,float x,float y,float z,int t){
	float InterpoGreyLevel;
	int xi,yi,zi;
	float xwm,ywm,zwm,xwp,ywp,zwp;
	float wmmm,wmmp,wmpm,wmpp,wpmm,wpmp,wppm,wppp;
	float wmm,wmp,wpm,wpp;
	
	
	//values out of the image
	if (x<0.) x=0.0001;
	if (x>=this->NX-1.) x=this->NX-1.0001;
	if (y<0.) y=0.0001;
	if (y>=this->NY-1.) y=this->NY-1.0001;
	if (z<0.) z=0.0001;
	if (z>=this->NZ-1.) z=this->NZ-1.0001;
	if (t<0) t=0;
	if (t>this->NT-1) t=this->NT-1;
	
	//closest entire value
	xi=static_cast<int>(x);  xwm=1-(x-static_cast<float>(xi));  xwp=x-static_cast<float>(xi);
	yi=static_cast<int>(y);  ywm=1-(y-static_cast<float>(yi));  ywp=y-static_cast<float>(yi);
	zi=static_cast<int>(z);  zwm=1-(z-static_cast<float>(zi));  zwp=z-static_cast<float>(zi);
	
	//interpolation
	if (this->NZ==1){ //2D IMAGE
		wmm=xwm*ywm;
		wmp=xwm*ywp;
		wpm=xwp*ywm;
		wpp=xwp*ywp;
		
		InterpoGreyLevel= wmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wpm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wpp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (yi+1)*this->NX+ (xi+1) ];
	}
	else{//3D IMAGE
		wmmm=xwm*ywm*zwm; wmmp=xwm*ywm*zwp; wmpm=xwm*ywp*zwm; wmpp=xwm*ywp*zwp;
		wpmm=xwp*ywm*zwm; wpmp=xwp*ywm*zwp; wppm=xwp*ywp*zwm; wppp=xwp*ywp*zwp;
		
		InterpoGreyLevel= wmmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmpm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wmpp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+   (xi) ];
		InterpoGreyLevel+=wpmm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wpmp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+   (yi)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wppm*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+   (zi)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
		InterpoGreyLevel+=wppp*this->VecField[ IdDirec*this->NXtYtZtT + t*this->NXtYtZ+ (zi+1)*this->NXtY+ (yi+1)*this->NX+ (xi+1) ];
	}
	
	return InterpoGreyLevel;
}

///same as above
float VectorField::G(int IdDirec,double x,double y,double z,int t){
	return this->G(IdDirec,(float) x,(float) y,(float) z,t);
}




///get the maximum of the absolute values of the vector field
float VectorField::GetMaxAbsVal(int t)
{
	float max=0.0;
	int direc,x,y,z;
	for(direc=0;direc<3;direc++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
	{
		if(max<abs(this->G(direc,x,y,z,t))){max = abs(this->G(direc,x,y,z,t));}
	}
	return max;
}


///read a vector field (in 3 nifti images -> X, Y, Z)
void VectorField::Read(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z){
	ScalarField VecField_X;
	ScalarField VecField_Y;
	ScalarField VecField_Z;
	int x,y,z,t;

	//read the vector fields
	VecField_X.Read(NameVecField_X);
	VecField_Y.Read(NameVecField_Y);
	VecField_Z.Read(NameVecField_Z);
	
	//message if there is already an image in InputImage with another size of the opened one
	if (this->NX!=0)
		if ((this->NX!=VecField_X.NX)||(this->NY!=VecField_X.NY)||(this->NZ!=VecField_X.NZ)||(this->NT!=VecField_X.NT))
			cout << "WARNING: THE SIZE OF A NON-NULL VECTOR FIELD IS CHANGED\n";
	
	//fill the parameters of the class and allocate the memory for the image
	this->NX=VecField_X.NX;
	this->NY=VecField_X.NY;
	this->NZ=VecField_X.NZ;
	this->NT=VecField_X.NT;
	this->NXtY=this->NX*this->NY;
	this->NXtYtZ=this->NXtY*this->NZ;
	this->NXtYtZtT=this->NXtYtZ*this->NT;
	this->VecField = new float [this->NXtYtZtT*3];
  this->Image2World[0][0]=VecField_X.Image2World[0][0];  this->Image2World[0][1]=VecField_X.Image2World[0][1];  this->Image2World[0][2]=VecField_X.Image2World[0][2];  this->Image2World[0][3]=VecField_X.Image2World[0][3];
  this->Image2World[1][0]=VecField_X.Image2World[1][0];  this->Image2World[1][1]=VecField_X.Image2World[1][1];  this->Image2World[1][2]=VecField_X.Image2World[1][2];  this->Image2World[1][3]=VecField_X.Image2World[1][3];
  this->Image2World[2][0]=VecField_X.Image2World[2][0];  this->Image2World[2][1]=VecField_X.Image2World[2][1];  this->Image2World[2][2]=VecField_X.Image2World[2][2];  this->Image2World[2][3]=VecField_X.Image2World[2][3];
  this->Image2World[3][0]=VecField_X.Image2World[3][0];  this->Image2World[3][1]=VecField_X.Image2World[3][1];  this->Image2World[3][2]=VecField_X.Image2World[3][2];  this->Image2World[3][3]=VecField_X.Image2World[3][3];

  invert_4t4quaternion(this->Image2World,this->World2Image);
  
	//cast the image to the format used in the class
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
		this->P(VecField_X.G(x,y,z,t),0,x,y,z,t);
	
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
		this->P(VecField_Y.G(x,y,z,t),1,x,y,z,t);
	
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
		this->P(VecField_Z.G(x,y,z,t),2,x,y,z,t);
}


///read a scalar vector and perform linear interpolation to give it a specific size
void VectorField::Read_and_Interpolate(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z,int NBX,int NBY,int NBZ,int rescaleVF){
	ScalarField OrigSF;
	int x,y,z,t;
	float x2,y2,z2;
	float scaleFactor;
	
	//0) init
	scaleFactor=1.;
  
	//1) X DIRECTION
	//1.1) read the scalar field in direction X at the original format
	OrigSF.Read(NameVecField_X);
	
  //1.2) fill the parameters of the class and allocate the memory for the vector field
	//(the directions Y and Z are supposed in the same format)
	this->NX=NBX;
	this->NY=NBY;
	this->NZ=NBZ;
	this->NT=OrigSF.NT;
	this->NXtY=this->NX*this->NY;
	this->NXtYtZ=this->NXtY*this->NZ;
	this->NXtYtZtT=this->NXtYtZ*this->NT;
	this->VecField = new float [this->NXtYtZtT*3];
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;

  
	//1.3) manage the scale factor
	if (rescaleVF!=0) scaleFactor=((float)this->NX)/((float)OrigSF.NX);
  this->Image2World[0][0]=OrigSF.Image2World[0][0]/scaleFactor;  this->Image2World[0][1]=OrigSF.Image2World[0][1]/scaleFactor;  this->Image2World[0][2]=OrigSF.Image2World[0][2]/scaleFactor;  this->Image2World[0][3]=OrigSF.Image2World[0][3];
	
  //1.4) interpolate the original image to compute the vector field in direction X
	for(t=0;t<this->NT;t++){
		for(z=0;z<this->NZ;z++){ 
			z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
			for(y=0;y<this->NY;y++){ 
				y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
				for(x=0;x<this->NX;x++){
					x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
					this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,0,x,y,z,t);
				}
			}
		}
	}

	//2) Y DIRECTION
	//2.1) read the scalar field in direction Y
	OrigSF.Read(NameVecField_Y);
	
	//2.2) manage the scale factor
	if (rescaleVF!=0) scaleFactor=((float)this->NY)/((float)OrigSF.NY);
  this->Image2World[1][0]=OrigSF.Image2World[1][0]/scaleFactor;  this->Image2World[1][1]=OrigSF.Image2World[1][1]/scaleFactor;  this->Image2World[1][2]=OrigSF.Image2World[1][2]/scaleFactor;  this->Image2World[1][3]=OrigSF.Image2World[1][3];
	
	//2.3) interpolate the original image to compute the vector field in direction Y
	for(t=0;t<this->NT;t++){
		for(z=0;z<this->NZ;z++){ 
			z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
			for(y=0;y<this->NY;y++){ 
				y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
				for(x=0;x<this->NX;x++){
					x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
					this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,1,x,y,z,t);
				}
			}
		}
	}
	
	
	//3) Z DIRECTION
	//3.1) read the scalar field in direction Z
	OrigSF.Read(NameVecField_Z);
	
	//3.2) manage the scale factor
	if (rescaleVF!=0) scaleFactor=((float)this->NZ)/((float)OrigSF.NZ);
  this->Image2World[2][0]=OrigSF.Image2World[2][0]/scaleFactor;  this->Image2World[2][1]=OrigSF.Image2World[2][1]/scaleFactor;  this->Image2World[2][2]=OrigSF.Image2World[2][2]/scaleFactor;  this->Image2World[2][3]=OrigSF.Image2World[2][3];
	
	//3.3) interpolate the original image to compute the vector field in direction Z
	for(t=0;t<this->NT;t++){
		for(z=0;z<this->NZ;z++){ 
			z2=static_cast<float>(z)*static_cast<float>(OrigSF.NZ-1)/static_cast<float>(this->NZ-1);
			for(y=0;y<this->NY;y++){ 
				y2=static_cast<float>(y)*static_cast<float>(OrigSF.NY-1)/static_cast<float>(this->NY-1);
				for(x=0;x<this->NX;x++){
					x2=static_cast<float>(x)*static_cast<float>(OrigSF.NX-1)/static_cast<float>(this->NX-1);
					this->P(OrigSF.G(x2,y2,z2,t)*scaleFactor,2,x,y,z,t);
				}
			}
		}
	}
  
  //compute the invert World2Image matrix
  invert_4t4quaternion(this->Image2World,this->World2Image);

  
}



///create a void vector field. No message is printed if  Verbose!=1
void VectorField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT,int Verbose){
  int x,y,z,t,direc;
  
  //message if there is already an image in InputImage with another size of the opened one
  if (Verbose==1) if (this->NX!=0)
    if ((this->NX!=NBX)||(this->NY!=NBY)||(this->NZ!=NBZ)||(this->NT!=NBT))
      cout << "WARNING: THE SIZE OF A NON-NULL VECTOR FIELD IS CHANGED\n";
  
  //image size
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  this->NT=NBT;
  this->NXtY=this->NX*this->NY;
  this->NXtYtZ=this->NXtY*this->NZ;
  this->NXtYtZtT=this->NXtYtZ*this->NT;
  
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  //allocate memory to cast (and eventually transform) the original template and target images
  this->VecField = new float [this->NXtYtZtT*3];
  
  //fill the image with 0.
  for(direc=0;direc<3;direc++) for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
    this->P(0.,direc,x,y,z,t);
}


///Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
void VectorField::SlashFieldSize(int verbative){

	//image size
	this->NX=1;
	this->NY=1;
	this->NZ=1;
	this->NT=1;
	this->NXtY=this->NX*this->NY;
	this->NXtYtZ=this->NXtY*this->NZ;
	this->NXtYtZtT=this->NXtYtZ*this->NT;
	
  this->Image2World[0][0]=1;  this->Image2World[0][1]=0;  this->Image2World[0][2]=0;  this->Image2World[0][3]=0;
  this->Image2World[1][0]=0;  this->Image2World[1][1]=1;  this->Image2World[1][2]=0;  this->Image2World[1][3]=0;
  this->Image2World[2][0]=0;  this->Image2World[2][1]=0;  this->Image2World[2][2]=1;  this->Image2World[2][3]=0;
  this->Image2World[3][0]=0;  this->Image2World[3][1]=0;  this->Image2World[3][2]=0;  this->Image2World[3][3]=1;
  
  //World to image matrix  (using functions of nifti1_io)
  invert_4t4quaternion(this->Image2World,this->World2Image);
  
  if (verbative==1) cout << "Slash some memory" << endl;

	//re-allocate memory
  if (this->NX!=0) delete this->VecField;
	
  this->VecField = new float [3];
}



///write a vector field (from 3 nifti images -> X, Y, Z)
void VectorField::Write(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z){
	ScalarField OutputImage;
	int x,y,z,t;
	
	//image to cast
	OutputImage.CreateVoidField(this->NX, this->NY, this->NZ,this->NT);
	
	//write the vector field in X direction
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
		OutputImage.P(this->G(0,x,y,z,t),x,y,z,t);
	
	OutputImage.Write(NameVecField_X);
	
	//write the vector field in Y direction
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
		OutputImage.P(this->G(1,x,y,z,t),x,y,z,t);
	
	OutputImage.Write(NameVecField_Y);
	
	//write the vector field in Z direction
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++)
		OutputImage.P(this->G(2,x,y,z,t),x,y,z,t);
	
	OutputImage.Write(NameVecField_Z);
}



///write a vector field (in 3 nifti images -> X, Y, Z)
//The 4th file is an input file containing the headers 
void VectorField::Write(char * NameVecField_X,char * NameVecField_Y,char * NameVecField_Z, char * ImageForHeaderName){
  nifti_1_header hdr_ref;
  FILE *fp_header;
  int i;
  int x,y,z,t;
  nifti1_extender pad={0,0,0,0};
  FILE *fp;
  int ret;
  float *data=NULL;
  
  
  //1) read the header of ImageForHeaderName
  fp_header = fopen(ImageForHeaderName,"rb");
  fread(&hdr_ref, MIN_HEADER_SIZE, 1, fp_header);
  fclose(fp_header);
  
  //2)give image parameters to the header
  hdr_ref.dim[0] = 4;
  hdr_ref.dim[1] = this->NX;
  hdr_ref.dim[2] = this->NY;
  hdr_ref.dim[3] = this->NZ;
  hdr_ref.dim[4] = this->NT;
  hdr_ref.datatype = NIFTI_TYPE_FLOAT32;
  hdr_ref.bitpix = 32; 
  hdr_ref.scl_inter = 0.0;
  hdr_ref.scl_slope = 1.0;
  hdr_ref.xyzt_units = NIFTI_UNITS_MM | NIFTI_UNITS_SEC;
  
  
  //3) allocate memory for the buffer
  //data = (float *) malloc(sizeof(float) * hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]);
  data = new float  [hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4]];
  
  //4.1) write the vector field in X direction
  //... fill the buffer
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(0,x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  //... save the field
  fp = fopen(NameVecField_X,"wb");
  
  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  fclose(fp);
	

  //4.2) write the vector field in Y direction
  //... fill the buffer
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(1,x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  //... save the field
  fp = fopen(NameVecField_Y,"wb");
  
  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  fclose(fp);
	
  //4.3) write the vector field in Z direction
  //... fill the buffer
  i=0;
  for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++){
    data[i] = this->G(2,x,y,z,t)/hdr_ref.scl_slope;
    i++;
  }
  
  //... save the field
  fp = fopen(NameVecField_Z,"wb");
  
  ret = fwrite(&hdr_ref, MIN_HEADER_SIZE, 1, fp);
  ret = fwrite(&pad, 4, 1, fp);
  ret = fwrite(data, (size_t)(hdr_ref.bitpix/8), hdr_ref.dim[1]*hdr_ref.dim[2]*hdr_ref.dim[3]*hdr_ref.dim[4], fp);
  
  fclose(fp);
  
  //5) free the memory
  //free(data);
  delete data;

}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           3:   FUNCTIONS FOR THE CLASS "TensorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



///constructor
TensorField::TensorField(void){
	this->NX=0;
	this->NY=0;
	this->NZ=0;
	this->NT=0;
}

///destructor
TensorField::~TensorField(void){
  if ((this->TField!=NULL)&&(this->NX>0)) delete this->TField;
}


///constructor
void TensorField::CreateVoidField(int NBX,int NBY,int NBZ,int NBT){
	int x,y,z,t,direc1,direc2;
	
  //message if there is already an image in InputImage with another size of the opened one
	if (this->NX!=0)
		if ((this->NX!=NBX)||(this->NY!=NBY)||(this->NZ!=NBZ)||(this->NT!=NBT))
			cout << "WARNING: THE SIZE OF A NON-NULL TENSOR FIELD IS CHANGED\n";
	
	//image size
	this->NX=NBX;
	this->NY=NBY;
	this->NZ=NBZ;
	this->NT=NBT;
  this->NXt9=this->NX*9;
  this->NXtYt9=this->NX*NY*9;
  this->NXtYtZt9=this->NX*NY*NZ*9;
  
  
	//allocate memory to cast (and eventually transform) the original template and target images
	this->TField = new float [this->NXtYtZt9*NT];
	
	//fill the image with 0.
	for(t=0;t<this->NT;t++) for(z=0;z<this->NZ;z++) for(y=0;y<this->NY;y++) for(x=0;x<this->NX;x++) for(direc2=0;direc2<3;direc2++) for(direc1=0;direc1<3;direc1++) 
		this->P(0.,direc1,direc2,x,y,z,t);
}

///put a value
void TensorField::P(float value,int IdDirec1,int IdDirec2,int x,int y,int z,int t){
	this->TField[IdDirec1+3*IdDirec2+9*x+this->NXt9*y+this->NXtYt9*z+this->NXtYtZt9*t]=value;
}

///add a value
void TensorField::Add(float value,int IdDirec1,int IdDirec2,int x,int y,int z,int t){
	this->TField[IdDirec1+3*IdDirec2+9*x+this->NXt9*y+this->NXtYt9*z+this->NXtYtZt9*t]+=value;
}


///get a value
float TensorField::G(int IdDirec1,int IdDirec2,int x,int y,int z,int t){
	return this->TField[IdDirec1+3*IdDirec2+9*x+this->NXt9*y+this->NXtYt9*z+this->NXtYtZt9*t];
}

//Add a tensorised vector to the existing tensor
void TensorField::AddTensorisedVector(float vec[3],int x,int y,int z,int t){
  this->Add(vec[0]*vec[0],0,0,x,y,z,t);
  this->Add(vec[1]*vec[0],1,0,x,y,z,t);
  this->Add(vec[2]*vec[0],2,0,x,y,z,t);
  this->Add(vec[0]*vec[1],0,1,x,y,z,t);
  this->Add(vec[1]*vec[1],1,1,x,y,z,t);
  this->Add(vec[2]*vec[1],2,1,x,y,z,t);
  this->Add(vec[0]*vec[2],0,2,x,y,z,t);
  this->Add(vec[1]*vec[2],1,2,x,y,z,t);
  this->Add(vec[2]*vec[2],2,2,x,y,z,t);
}



//Perform a principal component analysis of the 3*3 tensor
//The outputs are:
// lambda1,lambda2,lambda3: the eigenvalues in decrasing order
// vec1: 1st eigenvector  (must be initialised as vec1[3])
// vec2: 2nd eigenvector  (must be initialised as vec2[3])
// vec3: 3rd eigenvector  (must be initialised as vec3[3])
void TensorField::PCA(float * vec1,float * vec2,float * vec3,float * lambda1, float * lambda2, float * lambda3, int x,int y,int z,int t){
  float ** a;
  float ** q;
  float *  d;
  int i;
  
  //alloc temporary variables
  a = new float * [3];
  for(i=0;i<3;i++) a[i]=new float [3];
  q = new float * [3];
  for(i=0;i<3;i++) q[i]=new float [3];
  d = new float [3];
  
  //copy input values
  a[0][0]=this->G(0,0,x,y,z,t);  a[0][1]=this->G(0,1,x,y,z,t);  a[0][2]=this->G(0,2,x,y,z,t);
  a[1][0]=this->G(1,0,x,y,z,t);  a[1][1]=this->G(1,1,x,y,z,t);  a[1][2]=this->G(1,2,x,y,z,t);
  a[2][0]=this->G(2,0,x,y,z,t);  a[2][1]=this->G(2,1,x,y,z,t);  a[2][2]=this->G(2,2,x,y,z,t);
  
  
  //eigenvalue decomposition
  jacobi3(a,d,q);
  
  //copy output values
  *lambda1=d[0];  *lambda2=d[1];  *lambda3=d[2];
  
  vec1[0]=q[0][0];  vec1[1]=q[1][0];  vec1[2]=q[2][0];
  vec2[0]=q[0][1];  vec2[1]=q[1][1];  vec2[2]=q[2][1];
  vec3[0]=q[0][2];  vec3[1]=q[1][2];  vec3[2]=q[2][2];
  
  //delete allocated variables
  for(i=0;i<3;i++) delete a[i];
  delete a;
  for(i=0;i<3;i++) delete q[i];
  delete q;
  delete d;
  
}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           4: CLASS TO PERFORM CONVOLUTION AND DECONVOLUTION USING FFT
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///Constructor
FFTconvolver3D::FFTconvolver3D(){
	this->NX=0;
	this->NY=0;
	this->NZ=0;
	
	this->NXfft=0;
	this->NYfft=0;
	this->NZfft=0;
}

///destructor
FFTconvolver3D::~FFTconvolver3D(void){}


///Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
///4 Gaussians (set some weights to 0 if less Gaussians are required)
///* NX, NY, NZ: is the size of the input image
///* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
///* w2,sX2,sY2,sZ2,: weight of the 2nd Gaussian kernel and std. dev. in direction X, Y, Z
///* w3,sX3,sY3,sZ3,: weight of the 3rd Gaussian kernel and std. dev. in direction X, Y, Z
///* w4,sX4,sY4,sZ4,: weight of the 4th Gaussian kernel and std. dev. in direction X, Y, Z
///* w5,sX5,sY5,sZ5,: weight of the 5th Gaussian kernel and std. dev. in direction X, Y, Z
///* w6,sX6,sY6,sZ6,: weight of the 6th Gaussian kernel and std. dev. in direction X, Y, Z
///* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void FFTconvolver3D::InitiateConvolver(int NBX,int NBY, int NBZ,float w1,float sX1,float sY1,float sZ1,float w2,float sX2,float sY2,float sZ2,float w3,float sX3,float sY3,float sZ3,float w4,float sX4,float sY4,float sZ4,float w5,float sX5,float sY5,float sZ5,float w6,float sX6,float sY6,float sZ6,float w7,float sX7,float sY7,float sZ7){
	
	//set the size of the original image
	this->NX=NBX;
	this->NY=NBY;
	this->NZ=NBZ;
	
	//set the size of images for the FFT
	this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
	this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
	this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
	

	//cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << "\n";
	
	//allocate memory for the images for the FFT
	this->RealSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - real part
	this->ImagSignalForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //image  - imaginary part
	this->RealFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - real part
	this->ImagFilterForFFT.CreateVoidField(this->NXfft, this->NYfft, this->NZfft); //filter - imaginary part
	
	//allocate memory for the temporary image
	this->ImageTemp.CreateVoidField(this->NXfft,this->NYfft,this->NZfft);
	
	//define the kernel and transform it in Fourier spaces
	this->MakeSumOf7AnisotropicGaussianFilters(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7);
}

///change the kernel of the convolver (same notations as the constructor)
///... here the new kernel is normalized (w1+...+w7 is normalized to 1)
void FFTconvolver3D::ChangeKernel(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,
                                  float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,
                                  float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,
                                  float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,
                                  float weight5,float sigmaX5,float sigmaY5,float sigmaZ5,
                                  float weight6,float sigmaX6,float sigmaY6,float sigmaZ6,
                                  float weight7,float sigmaX7,float sigmaY7,float sigmaZ7){
  
  //define the kernel and transform it in Fourier spaces
	this->MakeSumOf7AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,weight2,sigmaX2,sigmaY2,sigmaZ2,weight3,sigmaX3,sigmaY3,sigmaZ3,weight4,sigmaX4,sigmaY4,sigmaZ4,weight5,sigmaX5,sigmaY5,sigmaZ5,weight6,sigmaX6,sigmaY6,sigmaZ6,weight7,sigmaX7,sigmaY7,sigmaZ7);
}

///change the kernel of the convolver (same notations as the constructor)
///... here the new kernel is not normalized
void FFTconvolver3D::ChangeKernel_SingleScale(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1){
  
  //define the kernel and transform it in Fourier spaces
	this->MakeSumOf7AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0);
}



///design a kernel that is the sum of up to 7 Gaussians and transform it in Fourier spaces
//if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
void FFTconvolver3D::MakeSumOf7AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,float weight5,float sigmaX5,float sigmaY5,float sigmaZ5,float weight6,float sigmaX6,float sigmaY6,float sigmaZ6,float weight7,float sigmaX7,float sigmaY7,float sigmaZ7,int NormalizeWeights){
	int k,x,y,z;
	float SumLoc;
	float weight,sigmaX,sigmaY,sigmaZ,sumWeight;
	
  //set RealFilterForFFT and ImagFilterForFFT to 0 in case it contains something
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(0.,x,y,z);
	
  for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT.P(0.,x,y,z);

  sumWeight=0;
  
	//compute and save the 7 kernels
	for (k=0;k<7;k++){
		//parameters of the current kernel
		if (k==6){weight=weight7; sigmaX=sigmaX7; sigmaY=sigmaY7; sigmaZ=sigmaZ7;}
		if (k==5){weight=weight6; sigmaX=sigmaX6; sigmaY=sigmaY6; sigmaZ=sigmaZ6;}
		if (k==4){weight=weight5; sigmaX=sigmaX5; sigmaY=sigmaY5; sigmaZ=sigmaZ5;}
		if (k==3){weight=weight4; sigmaX=sigmaX4; sigmaY=sigmaY4; sigmaZ=sigmaZ4;}
		if (k==2){weight=weight3; sigmaX=sigmaX3; sigmaY=sigmaY3; sigmaZ=sigmaZ3;}
		if (k==1){weight=weight2; sigmaX=sigmaX2; sigmaY=sigmaY2; sigmaZ=sigmaZ2;}
		if (k==0){weight=weight1; sigmaX=sigmaX1; sigmaY=sigmaY1; sigmaZ=sigmaZ1;}
		
    sumWeight+=weight;
    
    for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImageTemp.P(0,x,y,z);
    
		//design the current kernel with no influence of the weight
    for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=0;z<this->NZfft/2;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp( -(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp( -(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=0;z<this->NZfft/2;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp( -(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)(z*z)/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=0;y<this->NYfft/2;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp(-(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)(y*y)/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=0;x<this->NXfft/2;x++)
      this->ImageTemp.P((float)(exp(-(float)(x*x)/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    for (z=this->NZfft/2;z<this->NZfft;z++) for (y=this->NYfft/2;y<this->NYfft;y++) for (x=this->NXfft/2;x<this->NXfft;x++)
      this->ImageTemp.P((float)(exp(-(float)((this->NXfft-x)*(this->NXfft-x))/(2.*sigmaX*sigmaX) -(float)((this->NYfft-y)*(this->NYfft-y))/(2.*sigmaY*sigmaY) -(float)((this->NZfft-z)*(this->NZfft-z))/(2.*sigmaZ*sigmaZ))),x,y,z);
    
    
		//normalization of the current filter
		SumLoc=0.;
		for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) SumLoc+=this->ImageTemp.G(x,y,z);
		
		//copy in RealFilterForFFT
    if (fabs(SumLoc)>0.01){
			for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.Add(weight*this->ImageTemp.G(x,y,z)/SumLoc,x,y,z);
    }
    else{
      cout << "One kernel appears to be null or almost null" << endl;
    }
  }
  
  //normalize the weights
  if (NormalizeWeights!=0) for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealFilterForFFT.P(this->RealFilterForFFT.G(x,y,z)/sumWeight,x,y,z);

	//Transform RealFilterForFFT and ImagFilterForFFT in Fourier spaces
	this->DirectFFT(&this->RealFilterForFFT,&this->ImagFilterForFFT);
}


///Fast Fourier Transform of numerical recipies (slighly modified)
void FFTconvolver3D::four1NR(float data[], unsigned long nn, int isign){
	unsigned long n,mmax,m,j,istep,i;
	double wtemp,wr,wpr,wpi,wi,theta;
	float tempr,tempi;
	
	n=nn << 1;
	j=1;
	for (i=1;i<n;i+=2){
		if (j>i){
			tempr=data[j]; data[j]=data[i]; data[i]=tempr;
			tempr=data[j+1]; data[j+1]=data[i+1]; data[i+1]=tempr;
		}
		m=n >> 1;
		while ((m>=2) && (j>m)){
			j -= m;
			m >>= 1;
		}
		j += m;
	}
	mmax=2;
	while (n > mmax) {
		istep=mmax << 1;
		theta=isign*(6.28318530717959/mmax);
		wtemp=sin(0.5*theta);
		wpr = -2.0*wtemp*wtemp;
		wpi=sin(theta);
		wr=1.0;
		wi=0.0;
		for (m=1;m<mmax;m+=2) {
			for (i=m;i<=n;i+=istep) {
				j=i+mmax;
				tempr=wr*data[j]-wi*data[j+1];
				tempi=wr*data[j+1]+wi*data[j];
				data[j]=data[i]-tempr;
				data[j+1]=data[i+1]-tempi;
				data[i] += tempr;
				data[i+1] += tempi;
			}
			wr=(wtemp=wr)*wpr-wi*wpi+wr;
			wi=wi*wpr+wtemp*wpi+wi;
		}
		mmax=istep;
	}
}


#ifdef COMPILE_WITH_OPENMP

///Fast Fourier Transform
void FFTconvolver3D::DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  int MaxSizeXSizeYSizeZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  MaxSizeXSizeYSizeZ=SizeX;
  if (SizeY>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeY;
  if (SizeZ>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,dataX) 
  {
    dataX = new float [MaxSizeXSizeYSizeZ*2+1];
  
    //2) perform the fft along x axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){
      for (z = 0; z < SizeZ; z++){
        for (x = 0; x < SizeX; x++){
          dataX[2*x+1]=RealSignal->G(x, y, z);
          dataX[2*x+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeX, 1);
        for (x = 0; x < SizeX; x++){
          RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
          ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX, x, y, z);
        }  
      }
    }

    //3) perform the fft along y axis
    #pragma omp for
    for (x = 0; x < SizeX; x++){ 
      for (z = 0; z < SizeZ; z++) {
        for (y = 0; y < SizeY; y++){
          dataX[2*y+1]=RealSignal->G(x, y, z);
          dataX[2*y+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeY, 1);
        for (y = 0; y < SizeY; y++){
          RealSignal->P(dataX[2*y+1]/SqrtSizeY,x, y, z);
          ImaginarySignal->P(dataX[2*y+2]/SqrtSizeY, x, y, z);
        }
      }
    }
  
    //4) perform the fft along z axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){
      for (x = 0; x < SizeX; x++){
        for (z = 0; z < SizeZ; z++){
          dataX[2*z+1]=RealSignal->G(x, y, z);
          dataX[2*z+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeZ, 1);
        for (z = 0; z < SizeZ; z++){
          RealSignal->P(dataX[2*z+1]/SqrtSizeZ,x, y, z);
          ImaginarySignal->P(dataX[2*z+2]/SqrtSizeZ, x, y, z);
        }
      }
    }

    delete dataX;
  //END FORK FOR THREADS
  }

}

#else

///Fast Fourier Transform
void FFTconvolver3D::DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  
  //2) perform the fft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->G(x, y, z);
      dataX[2*x+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataX, (unsigned long)SizeX, 1);
    for (x = 0; x < SizeX; x++){
      RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
      ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX, x, y, z);
    }
  }
  delete dataX;
  
  //3) perform the fft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->G(x, y, z);
      dataY[2*y+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataY, (unsigned long)SizeY, 1);
    for (y = 0; y < SizeY; y++){
      RealSignal->P(dataY[2*y+1]/SqrtSizeY,x, y, z);
      ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
    }
  }
  delete dataY;
  
  
  //4) perform the fft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->G(x, y, z);
      dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataZ, (unsigned long)SizeZ, 1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->P(dataZ[2*z+1]/SqrtSizeZ,x, y, z);
      ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ, x, y, z);
    }
  }
  delete dataZ;
}

#endif



#ifdef COMPILE_WITH_OPENMP


///Inverse Fast Fourier Transform
void FFTconvolver3D::InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  int MaxSizeXSizeYSizeZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  MaxSizeXSizeYSizeZ=SizeX;
  if (SizeY>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeY;
  if (SizeZ>MaxSizeXSizeYSizeZ) MaxSizeXSizeYSizeZ=SizeZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,dataX) 
  {
    dataX = new float [MaxSizeXSizeYSizeZ*2+1];
  
    //2) perform the ifft along z axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){ 
      for (x = 0; x < SizeX; x++){
        for (z = 0; z < SizeZ; z++){
          dataX[2*z+1]=RealSignal->G(x, y, z);
          dataX[2*z+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeZ, -1);
        for (z = 0; z < SizeZ; z++){
          RealSignal->P(dataX[2*z+1]/SqrtSizeZ, x, y, z);
          ImaginarySignal->P(dataX[2*z+2]/SqrtSizeZ,x, y, z);
        }
      }
    }
    
    //3) perform the ifft along y axis
    #pragma omp for
    for (x = 0; x < SizeX; x++){
      for (z = 0; z < SizeZ; z++){
        for (y = 0; y < SizeY; y++){
          dataX[2*y+1]=RealSignal->G(x, y, z);
          dataX[2*y+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeY, -1);
        for (y = 0; y < SizeY; y++){
          RealSignal->P(dataX[2*y+1]/SqrtSizeY, x, y, z);
          ImaginarySignal->P(dataX[2*y+2]/SqrtSizeY, x, y, z);
        }
      }
    }
    
    //4) perform the ifft along x axis
    #pragma omp for
    for (y = 0; y < SizeY; y++){
      for (z = 0; z < SizeZ; z++){
        for (x = 0; x < SizeX; x++){
          dataX[2*x+1]=RealSignal->G(x, y, z);
          dataX[2*x+2]=ImaginarySignal->G(x, y, z);
        }
        four1NR(dataX, (unsigned long)SizeX, -1);
        for (x = 0; x < SizeX; x++){
          RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
          ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX,x, y, z);
        }
      }
    }
    
    delete dataX;
  //END FORK FOR THREADS
  }
}


#else

///Inverse Fast Fourier Transform
void FFTconvolver3D::InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;
  
  SqrtSizeX=static_cast<float>(sqrt(static_cast<double>(SizeX)));
  SqrtSizeY=static_cast<float>(sqrt(static_cast<double>(SizeY)));
  SqrtSizeZ=static_cast<float>(sqrt(static_cast<double>(SizeZ)));
  
  
  //2) perform the ifft along z axis
  dataZ = new float [SizeZ*2+1];
  for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
    for (z = 0; z < SizeZ; z++){
      dataZ[2*z+1]=RealSignal->G(x, y, z);
      dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataZ, (unsigned long)SizeZ, -1);
    for (z = 0; z < SizeZ; z++){
      RealSignal->P(dataZ[2*z+1]/SqrtSizeZ, x, y, z);
      ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ,x, y, z);
    }
  }
  delete dataZ;
  
  //3) perform the ifft along y axis
  dataY = new float [SizeY*2+1];
  for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
    for (y = 0; y < SizeY; y++){
      dataY[2*y+1]=RealSignal->G(x, y, z);
      dataY[2*y+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataY, (unsigned long)SizeY, -1);
    for (y = 0; y < SizeY; y++){
      RealSignal->P(dataY[2*y+1]/SqrtSizeY, x, y, z);
      ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
    }
  }
  delete dataY;
  
  //4) perform the ifft along x axis
  dataX = new float [SizeX*2+1];
  for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
    for (x = 0; x < SizeX; x++){
      dataX[2*x+1]=RealSignal->G(x, y, z);
      dataX[2*x+2]=ImaginarySignal->G(x, y, z);
    }
    four1NR(dataX, (unsigned long)SizeX, -1);
    for (x = 0; x < SizeX; x++){
      RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
      ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX,x, y, z);
    }
  }
  delete dataX;
}

#endif


///convolution of a 3D scalar field using the predifined kernel
void FFTconvolver3D::Convolution(ScalarField * SF){
	int x,y,z;
	float a,b,c,d;
	float CoefMult;
	
	//1) Copy the orginal image in the image that will be transformed
	for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
	for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
	
	for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++) this->RealSignalForFFT.P(SF->G(x,y,z),x,y,z);


	//2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
	
	//3) filtering in Fourier spaces
	CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
	
	for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
		a=this->RealSignalForFFT.G(x, y, z);
		b=this->ImagSignalForFFT.G(x, y, z);
		c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
		d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
		
		this->RealSignalForFFT.P(a*c-b*d, x, y, z);
		this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
	}
	
	//4) IFFT
	this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
	
	//5) Copy the image that has been convolved in the orginal image
	for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
		SF->P(this->RealSignalForFFT.G(x,y,z),x,y,z);
	}
}

///convolution of a 3D vector field using the predifined kernel
void FFTconvolver3D::Convolution(VectorField * VF){
	int x,y,z,i;
	float a,b,c,d;
	float CoefMult;
	
  for (i=0;i<3;i++){
		//1) Copy the orginal image in the image that will be transformed
		for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
		for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
		
		for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++) this->RealSignalForFFT.P(VF->G(i,x,y,z),x,y,z);
	
	
		//2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
		this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
		
		//3) filtering in Fourier spaces
		CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
		
		for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
			a=this->RealSignalForFFT.G(x, y, z);
			b=this->ImagSignalForFFT.G(x, y, z);
			c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
			d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
			
			this->RealSignalForFFT.P(a*c-b*d, x, y, z);
			this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
		}
		
		//4) IFFT
		this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
		
		//5) Copy the image that has been convolved in the orginal image
		for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++){
			VF->P(this->RealSignalForFFT.G(x,y,z),i,x,y,z);
		}
  }
}



///convolution of the real scalar field defined inside of the class
void FFTconvolver3D::Convolution(){
	int x,y,z;
	float a,b,c,d;
	float CoefMult;
	
	//1) Set to 0. all values that cannot be accessed by outside of the class
	for (z=this->NZ;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++)        for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
	for (z=0;z<this->NZ;z++)           for (y=this->NY;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
	for (z=0;z<this->NZ;z++)           for (y=0;y<this->NY;y++)           for (x=this->NX;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
	
	for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) 
    this->ImagSignalForFFT.P(0.,x,y,z);
	
	//2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
	
	//3) filtering in Fourier spaces
	CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
	
	for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
		a=this->RealSignalForFFT.G(x, y, z);
		b=this->ImagSignalForFFT.G(x, y, z);
		c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
		d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
		
		this->RealSignalForFFT.P(a*c-b*d, x, y, z);
		this->ImagSignalForFFT.P(c*b+a*d,x, y, z);
	}
	
	//4) IFFT
	this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
}

///put a value in the real part of the field that is transformed by the class
void FFTconvolver3D::P(float value,int x,int y, int z){
	this->RealSignalForFFT.P(value,x,y,z);
}

///put a value in the real part of the field that is transformed by the class
float FFTconvolver3D::G(int x,int y, int z){
	return this->RealSignalForFFT.G(x,y,z);
}

///deconvolution of a 3D scalar field using the predifined kernel
/// !!! NOT VALIDATED !!!
void FFTconvolver3D::Deconvolution(ScalarField * SF){
	int x,y,z;
	float a,b,c,d;
	float CoefMult;
	
	cout << "DECONVOLUTION SHOULD BE USED CARREFULLY HERE - NOT VALIDATED!!!\n";
	
	//1) Copy the orginal image in the image that will be transformed
	for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->RealSignalForFFT.P(0.,x,y,z);
	for (z=0;z<this->NZfft;z++) for (y=0;y<this->NYfft;y++) for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT.P(0.,x,y,z);
	
	for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
		this->RealSignalForFFT.P(SF->G(x,y,z),x,y,z);
	}
	//2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
	
	//3) filtering in Fourier spaces
	CoefMult=(float)(sqrt((double)this->RealSignalForFFT.NX)*sqrt((double)this->RealSignalForFFT.NY)*sqrt((double)this->RealSignalForFFT.NZ));
	
	for (z = 0; z < this->RealSignalForFFT.NZ; z++) for (y = 0; y < this->RealSignalForFFT.NY; y++) for (x = 0; x < this->RealSignalForFFT.NX; x++){
		a=this->RealSignalForFFT.G(x, y, z);
		b=this->ImagSignalForFFT.G(x, y, z);
		c=this->RealFilterForFFT.G(x, y, z)*CoefMult;
		d=this->ImagFilterForFFT.G(x, y, z)*CoefMult;
		
		this->RealSignalForFFT.P((a*c+b*d)/(c*c+d*d), x, y, z);
		this->ImagSignalForFFT.P((c*b-a*d)/(c*c+d*d),x, y, z);
	}
	//4) IFFT
	this->InverseFFT(&this->RealSignalForFFT,&this->ImagSignalForFFT);
	
	//5) Copy the image that has been deconvolved in the orginal image
	for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
		SF->P(this->RealSignalForFFT.G(x,y,z),x,y,z);
	}
}




///compute the update vector field using FFT
void SmoothVFUsingFFT(VectorField * SmoothedField,FFTconvolver3D * FFTconvolver_loc){
  int x,y,z;
  int direction;
  
  //smooth the information in the 3 directions X, Y, Z  (3D IMAGE)
  if (SmoothedField->NZ>1) for (direction=0;direction<3;direction++){
    //copy the region in the convolver
    for (z = 2; z < SmoothedField->NZ-2; z++) for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++){
      FFTconvolver_loc->P(SmoothedField->G(direction,x,y,z),x,y,z);
    }
    
    //convolve the region
    FFTconvolver_loc->Convolution();
    
    //copy the smoothed region in SmoothedField
    for (z = 0; z < SmoothedField->NZ-0; z++) for (y = 0; y < SmoothedField->NY-0; y++)  for (x = 0; x < SmoothedField->NX-0; x++){
      SmoothedField->P(FFTconvolver_loc->G(x,y,z),direction,x,y,z);
    }
  }

  //smooth the information in the 2 directions X, Y  (2D IMAGE)
  if (SmoothedField->NZ==1)  for (direction=0;direction<2;direction++){
    //copy the region in the convolver
    for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++){
      FFTconvolver_loc->P(SmoothedField->G(direction,x,y),x,y);
    }
    
    //convolve the region
    FFTconvolver_loc->Convolution();
    
    //copy the smoothed region in SmoothedField
    for (y = 0; y < SmoothedField->NY-0; y++)  for (x = 0; x < SmoothedField->NX-0; x++){
      SmoothedField->P(FFTconvolver_loc->G(x,y),direction,x,y);
    }
  }
  
  
}










///4.3: light weight convolver

///Constructor
LightFFTconvolver3D::LightFFTconvolver3D(){
  this->NX=0;    this->NY=0;    this->NZ=0;
  this->NXfft=0; this->NYfft=0; this->NZfft=0;
}

///destructor
LightFFTconvolver3D::~LightFFTconvolver3D(void){
  this->NX=0;    this->NY=0;    this->NZ=0;
  this->NXfft=0; this->NYfft=0; this->NZfft=0;
}

///Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 4 Gaussians (set some weights to 0 if less Gaussians are required)
///* NX, NY, NZ: is the size of the input image
///* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
///  ...
///* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void LightFFTconvolver3D::InitiateConvolver(int NBX,int NBY, int NBZ,float w1,float sX1,float sY1,float sZ1,float w2,float sX2,float sY2,float sZ2,float w3,float sX3,float sY3,float sZ3,float w4,float sX4,float sY4,float sZ4,float w5,float sX5,float sY5,float sZ5,float w6,float sX6,float sY6,float sZ6,float w7,float sX7,float sY7,float sZ7,int NormalizeWeights){
  //set the size of the original image
  this->NX=NBX;
  this->NY=NBY;
  this->NZ=NBZ;
  
  //set the size of images for the FFT
  this->NXfft=(int)(pow(2.,floor((log((double)this->NX)/log(2.))+0.99999))+0.00001); //smaller size higher than 'this->NX' and being a power of 2
  this->NYfft=(int)(pow(2.,floor((log((double)this->NY)/log(2.))+0.99999))+0.00001); // ... 'this->NY' ...
  this->NZfft=(int)(pow(2.,floor((log((double)this->NZ)/log(2.))+0.99999))+0.00001); // ... 'this->NZ' ...
  

  //cout << "Images to perform FFTs: " << this->NXfft << " , " << this->NYfft  << " , " << this->NZfft  << "\n";
  
  //allocate memory for the images for the FFT
  this->RealSignalForFFT_X.CreateVoidField(this->NXfft, 1, 1); //image  - real part
  this->ImagSignalForFFT_X.CreateVoidField(this->NXfft, 1, 1); //image  - imaginary part
  this->RealSignalForFFT_Y.CreateVoidField(1, this->NYfft, 1); //image  - real part
  this->ImagSignalForFFT_Y.CreateVoidField(1, this->NYfft, 1); //image  - imaginary part
  this->RealSignalForFFT_Z.CreateVoidField(1, 1, this->NZfft); //image  - real part
  this->ImagSignalForFFT_Z.CreateVoidField(1, 1, this->NZfft); //image  - imaginary part

  this->RealFilterForFFT_X.CreateVoidField(this->NXfft, 1, 1); //filter - real part
  this->ImagFilterForFFT_X.CreateVoidField(this->NXfft, 1, 1); //filter - imaginary part
  this->RealFilterForFFT_Y.CreateVoidField(1, this->NYfft, 1); //filter - real part
  this->ImagFilterForFFT_Y.CreateVoidField(1, this->NYfft, 1); //filter - imaginary part
  this->RealFilterForFFT_Z.CreateVoidField(1, 1, this->NZfft); //filter - real part
  this->ImagFilterForFFT_Z.CreateVoidField(1, 1, this->NZfft); //filter - imaginary part
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf7AnisotropicGaussianFilters(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,NormalizeWeights);
}

///change the kernel of the convolver (same notations as the constructor)
void LightFFTconvolver3D::ChangeKernel(float w1,float sX1,float sY1,float sZ1,float w2,float sX2,float sY2,float sZ2,float w3,float sX3,float sY3,float sZ3,float w4,float sX4,float sY4,float sZ4,float w5,float sX5,float sY5,float sZ5,float w6,float sX6,float sY6,float sZ6,float w7,float sX7,float sY7,float sZ7,int NormalizeWeights){
  
  //define the kernel and transform it in Fourier spaces
  this->MakeSumOf7AnisotropicGaussianFilters(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,NormalizeWeights);
}

///change the kernel of the convolver (same notations as the constructor)
///... here the new kernel is not normalized
void LightFFTconvolver3D::ChangeKernel_SingleScale(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1){
  
  //define the kernel and transform it in Fourier spaces
	this->MakeSumOf7AnisotropicGaussianFilters(weight1,sigmaX1,sigmaY1,sigmaZ1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0,1,1,1,0);
}




#ifdef COMPILE_WITH_OPENMP

///convolution of a 3D scalar field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(ScalarField * SF,int TimeFrame){
  int x,y,z,i;
  float a,b,c,d;
  float CoefMult;
  int MinTimeFrame;
  int MaxTimeFrame;
  ScalarField loc_RealSignalForFFT;     //for openmp
  ScalarField loc_ImagSignalForFFT;     //for openmp
  
  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=SF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,a,b,c,d,i,CoefMult,loc_RealSignalForFFT,loc_ImagSignalForFFT) 
  {
    for (i=MinTimeFrame;i<=MaxTimeFrame;i++){ //loop on the time frames
      //1) convolution on X axis
      loc_RealSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);     //added
      loc_ImagSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);     //added
  
      #pragma omp for
      for (y = 0; y < SF->NY; y++){
        for (z = 0; z < SF->NZ; z++){
          //1.1) Copy the orginal image in the image that will be transformed
          for (x=0;x<this->NXfft;x++) loc_RealSignalForFFT.P(0.,x,0,0);
          for (x=0;x<this->NXfft;x++) loc_ImagSignalForFFT.P(0.,x,0,0);
          
          for (x = 0; x < SF->NX; x++) loc_RealSignalForFFT.P(SF->G(x,y,z,i),x,0,0);
    
    
          //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NX);
          
          for (x = 0; x < loc_RealSignalForFFT.NX; x++){
            a=loc_RealSignalForFFT.G(x, 0, 0);
            b=loc_ImagSignalForFFT.G(x, 0, 0);
            c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
            d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,x,0,0);
            loc_ImagSignalForFFT.P(c*b+a*d,x,0,0);
          }
          
          //1.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.5) Copy the image that has been convolved in the orginal image
          for (x = 0; x < SF->NX; x++){
            SF->P(loc_RealSignalForFFT.G(x,0,0),x,y,z,i);
          }
        }
      }
    
      //2) convolution on Y axis
      loc_RealSignalForFFT.CreateVoidField(1, this->NYfft, 1,1,0,0);     //added
      loc_ImagSignalForFFT.CreateVoidField(1, this->NYfft, 1,1,0,0);     //added
  
      #pragma omp for
      for (x = 0; x < SF->NX; x++){
        for (z = 0; z < SF->NZ; z++) {
          //2.1) Copy the orginal image in the image that will be transformed
          for (y=0;y<this->NYfft;y++) loc_RealSignalForFFT.P(0.,0,y,0);
          for (y=0;y<this->NYfft;y++) loc_ImagSignalForFFT.P(0.,0,y,0);
          
          for (y=0;y<SF->NY;y++) loc_RealSignalForFFT.P(SF->G(x,y,z,i),0,y,0);
    
    
          //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NY);
          
          for (y = 0; y < loc_RealSignalForFFT.NY; y++){
            a=loc_RealSignalForFFT.G(0, y, 0);
            b=loc_ImagSignalForFFT.G(0, y, 0);
            c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
            d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,y,0);
            loc_ImagSignalForFFT.P(c*b+a*d,0,y,0);
          }
          
          //2.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.5) Copy the image that has been convolved in the orginal image
          for (y = 0; y < SF->NY; y++){
            SF->P(loc_RealSignalForFFT.G(0,y,0),x,y,z,i);
          }
        }
      }


      //3) convolution on Z axis
      loc_RealSignalForFFT.CreateVoidField(1, 1, this->NZfft,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1, 1, this->NZfft,1,0,0);    //added
  
      #pragma omp for
      for (y = 0; y < SF->NY; y++){
        for (x = 0; x < SF->NX; x++){
          //3.1) Copy the orginal image in the image that will be transformed
          for (z=0;z<this->NZfft;z++) loc_RealSignalForFFT.P(0.,0,0,z);
          for (z=0;z<this->NZfft;z++) loc_ImagSignalForFFT.P(0.,0,0,z);
          
          for (z = 0; z < SF->NZ; z++) loc_RealSignalForFFT.P(SF->G(x,y,z,i),0,0,z);
    
    
          //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NZ);
          
          for (z = 0; z < loc_RealSignalForFFT.NZ; z++){
            a=loc_RealSignalForFFT.G(0, 0, z);
            b=loc_ImagSignalForFFT.G(0, 0, z);
            c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
            d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,0,z);
            loc_ImagSignalForFFT.P(c*b+a*d,0,0,z);
          }
          
          //3.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.5) Copy the image that has been convolved in the orginal image
          for (z = 0; z < SF->NZ; z++){
            SF->P(loc_RealSignalForFFT.G(0,0,z),x,y,z,i);
          }
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///convolution of a 3D scalar field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(ScalarField * SF,int TimeFrame){
  int x,y,z,i;
  float a,b,c,d;
  float CoefMult;
  int MinTimeFrame;
  int MaxTimeFrame;

  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=SF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  for (i=MinTimeFrame;i<=MaxTimeFrame;i++){ //loop on the time frames
    //1) convolution on X axis
    for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      for (x = 0; x < SF->NX; x++) this->RealSignalForFFT_X.P(SF->G(x,y,z,i),x,0,0);


      //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
      
      for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
        a=this->RealSignalForFFT_X.G(x, 0, 0);
        b=this->ImagSignalForFFT_X.G(x, 0, 0);
        c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
        d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
        
        this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
        this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
      }
      
      //1.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.5) Copy the image that has been convolved in the orginal image
      for (x = 0; x < SF->NX; x++){
        SF->P(this->RealSignalForFFT_X.G(x,0,0),x,y,z,i);
      }
    }
    
    //2) convolution on Y axis
    for (z = 0; z < SF->NZ; z++) for (x = 0; x < SF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      for (y=0;y<SF->NY;y++) this->RealSignalForFFT_Y.P(SF->G(x,y,z,i),0,y,0);


      //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
      
      for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
        a=this->RealSignalForFFT_Y.G(0, y, 0);
        b=this->ImagSignalForFFT_Y.G(0, y, 0);
        c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
        d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
        
        this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
        this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
      }
      
      //2.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.5) Copy the image that has been convolved in the orginal image
      for (y = 0; y < SF->NY; y++){
        SF->P(this->RealSignalForFFT_Y.G(0,y,0),x,y,z,i);
      }
    }

    //3) convolution on Z axis
    for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      for (z = 0; z < SF->NZ; z++) this->RealSignalForFFT_Z.P(SF->G(x,y,z,i),0,0,z);


      //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
      
      for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
        a=this->RealSignalForFFT_Z.G(0, 0, z);
        b=this->ImagSignalForFFT_Z.G(0, 0, z);
        c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
        d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
        
        this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
        this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
      }
      
      //3.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.5) Copy the image that has been convolved in the orginal image
      for (z = 0; z < SF->NZ; z++){
        SF->P(this->RealSignalForFFT_Z.G(0,0,z),x,y,z,i);
      }
    }
  }
}

#endif

#ifdef COMPILE_WITH_OPENMP

///convolution of a 3D vector field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(VectorField * VF,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;
  ScalarField loc_RealSignalForFFT;     //for openmp
  ScalarField loc_ImagSignalForFFT;     //for openmp

  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,a,b,c,d,i,j,CoefMult,loc_RealSignalForFFT,loc_ImagSignalForFFT) 
  {
    for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
      //1) convolution on X axis
      loc_RealSignalForFFT.CreateVoidField(this->NXfft, 1,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(this->NXfft, 1,1,1,0,0);    //added

      #pragma omp for
      for (y = 0; y < VF->NY; y++){
        for (z = 0; z < VF->NZ; z++){
          //1.1) Copy the orginal image in the image that will be transformed
          for (x=0;x<this->NXfft;x++) loc_RealSignalForFFT.P(0.,x,0,0);
          for (x=0;x<this->NXfft;x++) loc_ImagSignalForFFT.P(0.,x,0,0);
          
          for (x = 0; x < VF->NX; x++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),x,0,0);
    
    
          //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NX);
          
          for (x = 0; x < loc_RealSignalForFFT.NX; x++){
            a=loc_RealSignalForFFT.G(x, 0, 0);
            b=loc_ImagSignalForFFT.G(x, 0, 0);
            c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
            d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,x,0,0);
            loc_ImagSignalForFFT.P(c*b+a*d,x,0,0);
          }
          
          //1.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.5) Copy the image that has been convolved in the orginal image
          for (x = 0; x < VF->NX; x++){
            VF->P(loc_RealSignalForFFT.G(x,0,0),i,x,y,z,j);
          }
        }
      }
      
      //2) convolution on Y axis
      loc_RealSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      
      #pragma omp for
      for (x = 0; x < VF->NX; x++){ 
        for (z = 0; z < VF->NZ; z++){
          //2.1) Copy the orginal image in the image that will be transformed
          for (y=0;y<this->NYfft;y++) loc_RealSignalForFFT.P(0.,0,y,0);
          for (y=0;y<this->NYfft;y++) loc_ImagSignalForFFT.P(0.,0,y,0);
          
          for (y = 0; y < VF->NY; y++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,y,0);
    
    
          //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NY);
          
          for (y = 0; y < loc_RealSignalForFFT.NY; y++){
            a=loc_RealSignalForFFT.G(0, y, 0);
            b=loc_ImagSignalForFFT.G(0, y, 0);
            c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
            d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,y,0);
            loc_ImagSignalForFFT.P(c*b+a*d,0,y,0);
          }
          
          //2.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.5) Copy the image that has been convolved in the orginal image
          for (y = 0; y < VF->NY; y++){
            VF->P(loc_RealSignalForFFT.G(0,y,0),i,x,y,z,j);
          }
        }
      }
      
      //3) convolution on Z axis
      loc_RealSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      
      #pragma omp for
      for (y = 0; y < VF->NY; y++){ 
        for (x = 0; x < VF->NX; x++){
          //3.1) Copy the orginal image in the image that will be transformed
          for (z=0;z<this->NZfft;z++) loc_RealSignalForFFT.P(0.,0,0,z);
          for (z=0;z<this->NZfft;z++) loc_ImagSignalForFFT.P(0.,0,0,z);
          
          for (z = 0; z < VF->NZ; z++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,0,z);
    
    
          //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NZ);
          
          for (z = 0; z < loc_RealSignalForFFT.NZ; z++){
            a=loc_RealSignalForFFT.G(0, 0, z);
            b=loc_ImagSignalForFFT.G(0, 0, z);
            c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
            d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,0,z);
            loc_ImagSignalForFFT.P(c*b+a*d,0,0,z);
          }
          
          //3.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.5) Copy the image that has been convolved in the orginal image
          for (z = 0; z < VF->NZ; z++){
            VF->P(loc_RealSignalForFFT.G(0,0,z),i,x,y,z,j);
          }
        }
      }
    } 
  //END FORK FOR THREADS
  }
}

#else

///convolution of a 3D vector field using the predifined kernel
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::Convolution(VectorField * VF,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;

  //0) Define the time frames to smooth
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }

  
  for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
    //1) convolution on X axis
    for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      for (x = 0; x < VF->NX; x++) this->RealSignalForFFT_X.P(VF->G(i,x,y,z,j),x,0,0);


      //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
      
      for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
        a=this->RealSignalForFFT_X.G(x, 0, 0);
        b=this->ImagSignalForFFT_X.G(x, 0, 0);
        c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
        d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
        
        this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
        this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
      }
      
      //1.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.5) Copy the image that has been convolved in the orginal image
      for (x = 0; x < VF->NX; x++){
        VF->P(this->RealSignalForFFT_X.G(x,0,0),i,x,y,z,j);
      }
    }
    
    //2) convolution on Y axis
    for (z = 0; z < VF->NZ; z++) for (x = 0; x < VF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      for (y = 0; y < VF->NY; y++) this->RealSignalForFFT_Y.P(VF->G(i,x,y,z,j),0,y,0);


      //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
      
      for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
        a=this->RealSignalForFFT_Y.G(0, y, 0);
        b=this->ImagSignalForFFT_Y.G(0, y, 0);
        c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
        d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
        
        this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
        this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
      }
      
      //2.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.5) Copy the image that has been convolved in the orginal image
      for (y = 0; y < VF->NY; y++){
        VF->P(this->RealSignalForFFT_Y.G(0,y,0),i,x,y,z,j);
      }
    }
    
    //3) convolution on Z axis
    for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      for (z = 0; z < VF->NZ; z++) this->RealSignalForFFT_Z.P(VF->G(i,x,y,z,j),0,0,z);


      //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
      
      for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
        a=this->RealSignalForFFT_Z.G(0, 0, z);
        b=this->ImagSignalForFFT_Z.G(0, 0, z);
        c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
        d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
        
        this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
        this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
      }
      
      //3.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.5) Copy the image that has been convolved in the orginal image
      for (z = 0; z < VF->NZ; z++){
        VF->P(this->RealSignalForFFT_Z.G(0,0,z),i,x,y,z,j);
      }
    }
  } 
}

#endif

#ifdef COMPILE_WITH_OPENMP

///convolution of a 3D vector field using the predifined kernel. Convolution is performed in the ROI defined by (xmin, xmax, ymin, ymax, zmin, zmax) only.
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::ConvolutionInROI(VectorField * VF,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;
  ScalarField loc_RealSignalForFFT;     //for openmp
  ScalarField loc_ImagSignalForFFT;     //for openmp

  //0) Define the time frames to smooth and check the boundaries
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  if (xmin<0)      {xmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (ymin<0)      {ymin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (zmin<0)      {zmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (xmax>VF->NX) {xmax=VF->NX; cout << "Region boundaries are not well defined" << endl;}
  if (ymax>VF->NY) {ymax=VF->NY; cout << "Region boundaries are not well defined" << endl;}
  if (zmax>VF->NZ) {zmax=VF->NZ; cout << "Region boundaries are not well defined" << endl;}

  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,a,b,c,d,i,j,CoefMult,loc_RealSignalForFFT,loc_ImagSignalForFFT) 
  {
    for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
      //1) convolution on X axis
      loc_RealSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(this->NXfft,1,1,1,0,0);    //added
      
      #pragma omp for
      for (y = ymin; y < ymax; y++){ 
        for (z = zmin; z < zmax; z++){
          //1.1) Copy the orginal image in the image that will be transformed
          for (x=0;x<this->NXfft;x++) loc_RealSignalForFFT.P(0.,x,0,0);
          for (x=0;x<this->NXfft;x++) loc_ImagSignalForFFT.P(0.,x,0,0);
          
          for (x = 0; x < VF->NX; x++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),x,0,0);
  
          //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NX);
          
          for (x = 0; x < loc_RealSignalForFFT.NX; x++){
            a=loc_RealSignalForFFT.G(x, 0, 0);
            b=loc_ImagSignalForFFT.G(x, 0, 0);
            c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
            d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,x,0,0);
            loc_ImagSignalForFFT.P(c*b+a*d,x,0,0);
          }
          
          //1.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,0);
          
          //1.5) Copy the image that has been convolved in the orginal image
          for (x = xmin; x < xmax; x++){
            VF->P(loc_RealSignalForFFT.G(x,0,0),i,x,y,z,j);
          }
        }
      }
      
      
      //2) convolution on Y axis
      loc_RealSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,this->NYfft,1,1,0,0);    //added
      
      #pragma omp for
      for (x = xmin; x < xmax; x++){ 
        for (z = zmin; z < zmax; z++){
          //2.1) Copy the orginal image in the image that will be transformed
          for (y=0;y<this->NYfft;y++) loc_RealSignalForFFT.P(0.,0,y,0);
          for (y=0;y<this->NYfft;y++) loc_ImagSignalForFFT.P(0.,0,y,0);
          
          for (y = 0; y < VF->NY; y++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,y,0);
    
    
          //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NY);
          
          for (y = 0; y < loc_RealSignalForFFT.NY; y++){
            a=loc_RealSignalForFFT.G(0, y, 0);
            b=loc_ImagSignalForFFT.G(0, y, 0);
            c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
            d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,y,0);
            loc_ImagSignalForFFT.P(c*b+a*d,0,y,0);
          }
          
          //2.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,1);
          
          //2.5) Copy the image that has been convolved in the orginal image
          for (y = ymin; y < ymax; y++){
            VF->P(loc_RealSignalForFFT.G(0,y,0),i,x,y,z,j);
          }
        }
      }
      
      //3) convolution on Z axis
      loc_RealSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      loc_ImagSignalForFFT.CreateVoidField(1,1,this->NZfft,1,0,0);    //added
      
      #pragma omp for
      for (y = ymin; y < ymax; y++){
        for (x = xmin; x < xmax; x++){
          //3.1) Copy the orginal image in the image that will be transformed
          for (z=0;z<this->NZfft;z++) loc_RealSignalForFFT.P(0.,0,0,z);
          for (z=0;z<this->NZfft;z++) loc_ImagSignalForFFT.P(0.,0,0,z);
          
          for (z = 0; z < VF->NZ; z++) loc_RealSignalForFFT.P(VF->G(i,x,y,z,j),0,0,z);
    
    
          //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
          this->DirectFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.3) filtering in Fourier spaces
          CoefMult=sqrt((float)loc_RealSignalForFFT.NZ);
          
          for (z = 0; z < loc_RealSignalForFFT.NZ; z++){
            a=loc_RealSignalForFFT.G(0, 0, z);
            b=loc_ImagSignalForFFT.G(0, 0, z);
            c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
            d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
            
            loc_RealSignalForFFT.P(a*c-b*d,0,0,z);
            loc_ImagSignalForFFT.P(c*b+a*d,0,0,z);
          }
          
          //3.4) IFFT
          this->InverseFFT(&loc_RealSignalForFFT,&loc_ImagSignalForFFT,2);
          
          //3.5) Copy the image that has been convolved in the orginal image
          for (z = zmin; z < zmax; z++){
            VF->P(loc_RealSignalForFFT.G(0,0,z),i,x,y,z,j);
          }
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///convolution of a 3D vector field using the predifined kernel. Convolution is performed in the ROI defined by (xmin, xmax, ymin, ymax, zmin, zmax) only.
///If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
void LightFFTconvolver3D::ConvolutionInROI(VectorField * VF,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax,int TimeFrame){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i,j;
  int MinTimeFrame;
  int MaxTimeFrame;

  //0) Define the time frames to smooth and check the boundaries
  if (TimeFrame<0){
    MinTimeFrame=0;
    MaxTimeFrame=VF->NT-1;
  }
  else{
    MinTimeFrame=TimeFrame;
    MaxTimeFrame=TimeFrame;
  }
  
  if (xmin<0)      {xmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (ymin<0)      {ymin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (zmin<0)      {zmin=0;      cout << "Region boundaries are not well defined" << endl;}
  if (xmax>VF->NX) {xmax=VF->NX; cout << "Region boundaries are not well defined" << endl;}
  if (ymax>VF->NY) {ymax=VF->NY; cout << "Region boundaries are not well defined" << endl;}
  if (zmax>VF->NZ) {zmax=VF->NZ; cout << "Region boundaries are not well defined" << endl;}

  
  for (j=MinTimeFrame;j<=MaxTimeFrame;j++) for (i=0;i<3;i++){ // loop on the time frames and directions
    //1) convolution on X axis
    for (z = zmin; z < zmax; z++) for (y = ymin; y < ymax; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      for (x = 0; x < VF->NX; x++) this->RealSignalForFFT_X.P(VF->G(i,x,y,z,j),x,0,0);


      //1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
      
      for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
        a=this->RealSignalForFFT_X.G(x, 0, 0);
        b=this->ImagSignalForFFT_X.G(x, 0, 0);
        c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
        d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
        
        this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
        this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
      }
      
      //1.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
      
      //1.5) Copy the image that has been convolved in the orginal image
      for (x = xmin; x < xmax; x++){
        VF->P(this->RealSignalForFFT_X.G(x,0,0),i,x,y,z,j);
      }
    }
    
    //2) convolution on Y axis
    for (z = zmin; z < zmax; z++) for (x = xmin; x < xmax; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      for (y = 0; y < VF->NY; y++) this->RealSignalForFFT_Y.P(VF->G(i,x,y,z,j),0,y,0);


      //2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
      
      for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
        a=this->RealSignalForFFT_Y.G(0, y, 0);
        b=this->ImagSignalForFFT_Y.G(0, y, 0);
        c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
        d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
        
        this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
        this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
      }
      
      //2.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
      
      //2.5) Copy the image that has been convolved in the orginal image
      for (y = ymin; y < ymax; y++){
        VF->P(this->RealSignalForFFT_Y.G(0,y,0),i,x,y,z,j);
      }
    }
    
    //3) convolution on Z axis
    for (y = ymin; y < ymax; y++) for (x = xmin; x < xmax; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      for (z = 0; z < VF->NZ; z++) this->RealSignalForFFT_Z.P(VF->G(i,x,y,z,j),0,0,z);


      //3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
      this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.3) filtering in Fourier spaces
      CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
      
      for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
        a=this->RealSignalForFFT_Z.G(0, 0, z);
        b=this->ImagSignalForFFT_Z.G(0, 0, z);
        c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
        d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
        
        this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
        this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
      }
      
      //3.4) IFFT
      this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
      
      //3.5) Copy the image that has been convolved in the orginal image
      for (z = zmin; z < zmax; z++){
        VF->P(this->RealSignalForFFT_Z.G(0,0,z),i,x,y,z,j);
      }
    }
  } 
}
#endif


///Hack to perform convolution of the 3D vector field 'VF' in a masked region with mirror conditions.
///  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
///  -> sX1,sY1,sZ1: size of the smoothing kernel
void LightFFTconvolver3D::Convolution_Mask_Mirror(VectorField * VF,ScalarField * Mask, int MaskId){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i;
  int PreviousInMask,MirrorInMask;
  int Something;
  float MaskId_fl;
  
  MaskId_fl=static_cast<float>(MaskId);
  
  for (i=0;i<3;i++){ // loop on the time frames and directions
    //1) convolution on X axis
    for (z = 0; z < VF->NZ; z++) for (y = 0; y < VF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      //1.1.1) Init
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      Something=0;
      for (x = 0; x < VF->NX; x++)  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_X.P(VF->G(i,x,y,z),x,0,0);
	Something=1;
      }
      
      if (Something==1){
	//1.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (x = 0; x < VF->NX; x++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}  

	//1.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (x = VF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_X.G(x,0,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_X.G(x,0,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}
	
	//1.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (x = 0; x < VF->NX; x++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)<0)  this->RealSignalForFFT_X.P(VF->G(i,MirrorInMask,y,z),x,0,0);
	  }
	}

	//1.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (x = VF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>VF->NX-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)>0)  this->RealSignalForFFT_X.P(VF->G(i,MirrorInMask,y,z),x,0,0);
	  }
	}	

	//1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
	
	for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
	  a=this->RealSignalForFFT_X.G(x, 0, 0);
	  b=this->ImagSignalForFFT_X.G(x, 0, 0);
	  c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  
	  this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
	  this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
	}
	
	//1.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.5) Copy the image that has been convolved in the orginal image
	for (x = 0; x < VF->NX; x++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  VF->P(this->RealSignalForFFT_X.G(x,0,0),i,x,y,z);
	}
      }
    }
  
    
    //2) convolution on Y axis
    for (z = 0; z < VF->NZ; z++) for (x = 0; x < VF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      //2.1.1) Init
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      Something=0;
      for (y = 0; y < VF->NY; y++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_Y.P(VF->G(i,x,y,z),0,y,0);
	Something=1;
      }
      
      
      if (Something==1){
	//2.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (y = 0; y < VF->NY; y++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (y = VF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Y.G(0,y,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Y.G(0,y,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (y = 0; y < VF->NY; y++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)<0) this->RealSignalForFFT_Y.P(VF->G(i,x,MirrorInMask,z),0,y,0);
	  }
	}

	//2.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (y = VF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>VF->NY-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)>0)  this->RealSignalForFFT_Y.P(VF->G(i,x,MirrorInMask,z),0,y,0);
	  }
	}


	//2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
	
	for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
	  a=this->RealSignalForFFT_Y.G(0, y, 0);
	  b=this->ImagSignalForFFT_Y.G(0, y, 0);
	  c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  
	  this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
	  this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
	}
	
	//2.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.5) Copy the image that has been convolved in the orginal image
	for (y = 0; y < VF->NY; y++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  VF->P(this->RealSignalForFFT_Y.G(0,y,0),i,x,y,z);
	}
      }
    }
    
    //3) convolution on Z axis
    for (y = 0; y < VF->NY; y++) for (x = 0; x < VF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      
      //3.1.1) Init
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      Something=0;
      for (z = 0; z < VF->NZ; z++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){ 
	this->RealSignalForFFT_Z.P(VF->G(i,x,y,z),0,0,z);
	Something=1;
      }
      
      
            
      if (Something==1){
	//3.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (z = VF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Z.G(0,0,z)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Z.G(0,0,z))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)<0) this->RealSignalForFFT_Z.P(VF->G(i,x,y,MirrorInMask),0,0,z);
	  }
	}

	//3.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (z = VF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>VF->NZ-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)>0) this->RealSignalForFFT_Z.P(VF->G(i,x,y,MirrorInMask),0,0,z);
	  }
	}


	
	//3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
	
	for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
	  a=this->RealSignalForFFT_Z.G(0, 0, z);
	  b=this->ImagSignalForFFT_Z.G(0, 0, z);
	  c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  
	  this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
	  this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
	}
	
	//3.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.5) Copy the image that has been convolved in the orginal image
	for (z = 0; z < VF->NZ; z++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  VF->P(this->RealSignalForFFT_Z.G(0,0,z),i,x,y,z);
	}
      }
    }
  } 
}





///Hack to perform convolution of the 3D scalar field 'SF' in a masked region with mirror conditions.
///  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
void LightFFTconvolver3D::Convolution_Mask_Mirror(ScalarField * SF,ScalarField * Mask, int MaskId){
  int x,y,z;
  float a,b,c,d;
  float CoefMult;
  int i;
  int PreviousInMask,MirrorInMask;
  int Something;
  float MaskId_fl;
  
  MaskId_fl=static_cast<float>(MaskId);
  
    //1) convolution on X axis
    for (z = 0; z < SF->NZ; z++) for (y = 0; y < SF->NY; y++){
      //1.1) Copy the orginal image in the image that will be transformed
      //1.1.1) Init
      for (x=0;x<this->NXfft;x++) this->RealSignalForFFT_X.P(0.,x,0,0);
      for (x=0;x<this->NXfft;x++) this->ImagSignalForFFT_X.P(0.,x,0,0);
      
      Something=0;
      for (x = 0; x < SF->NX; x++)  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_X.P(SF->G(x,y,z),x,0,0);
	Something=1;
      }
      
      if (Something==1){
	//1.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (x = 0; x < SF->NX; x++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}  

	//1.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (x = SF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_X.G(x,0,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_X.G(x,0,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_X.P(static_cast<float>(MirrorInMask),x,0,0);
	  }
	}
	
	//1.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (x = 0; x < SF->NX; x++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)<0)  this->RealSignalForFFT_X.P(SF->G(MirrorInMask,y,z),x,0,0);
	  }
	}

	//1.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (x = SF->NX-1; x >=0 ; x--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=x;MirrorInMask=x-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>SF->NX-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(MirrorInMask,y,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_X.G(x,0,0)>0)  this->RealSignalForFFT_X.P(SF->G(MirrorInMask,y,z),x,0,0);
	  }
	}	

	//1.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_X.NX);
	
	for (x = 0; x < this->RealSignalForFFT_X.NX; x++){
	  a=this->RealSignalForFFT_X.G(x, 0, 0);
	  b=this->ImagSignalForFFT_X.G(x, 0, 0);
	  c=this->RealFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  d=this->ImagFilterForFFT_X.G(x, 0, 0)*CoefMult;
	  
	  this->RealSignalForFFT_X.P(a*c-b*d,x,0,0);
	  this->ImagSignalForFFT_X.P(c*b+a*d,x,0,0);
	}
	
	//1.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_X,&this->ImagSignalForFFT_X,0);
	
	//1.5) Copy the image that has been convolved in the orginal image
	for (x = 0; x < SF->NX; x++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  SF->P(this->RealSignalForFFT_X.G(x,0,0),x,y,z);
	}
      }
    }
  
    
    //2) convolution on Y axis
    for (z = 0; z < SF->NZ; z++) for (x = 0; x < SF->NX; x++){
      //2.1) Copy the orginal image in the image that will be transformed
      //2.1.1) Init
      for (y=0;y<this->NYfft;y++) this->RealSignalForFFT_Y.P(0.,0,y,0);
      for (y=0;y<this->NYfft;y++) this->ImagSignalForFFT_Y.P(0.,0,y,0);
      
      Something=0;
      for (y = 0; y < SF->NY; y++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	this->RealSignalForFFT_Y.P(SF->G(x,y,z),0,y,0);
	Something=1;
      }
      
      
      if (Something==1){
	//2.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (y = 0; y < SF->NY; y++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (y = SF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Y.G(0,y,0)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Y.G(0,y,0))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Y.P(static_cast<float>(MirrorInMask),0,y,0);
	  }
	}

	//2.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (y = 0; y < SF->NY; y++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)<0) this->RealSignalForFFT_Y.P(SF->G(x,MirrorInMask,z),0,y,0);
	  }
	}

	//2.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (y = SF->NY-1; y >=0 ; y--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=y;MirrorInMask=y-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>SF->NY-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,MirrorInMask,z)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask; //should be OK
	    if (this->RealSignalForFFT_Y.G(0,y,0)>0)  this->RealSignalForFFT_Y.P(SF->G(x,MirrorInMask,z),0,y,0);
	  }
	}


	//2.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Y.NY);
	
	for (y = 0; y < this->RealSignalForFFT_Y.NY; y++){
	  a=this->RealSignalForFFT_Y.G(0, y, 0);
	  b=this->ImagSignalForFFT_Y.G(0, y, 0);
	  c=this->RealFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  d=this->ImagFilterForFFT_Y.G(0, y, 0)*CoefMult;
	  
	  this->RealSignalForFFT_Y.P(a*c-b*d,0,y,0);
	  this->ImagSignalForFFT_Y.P(c*b+a*d,0,y,0);
	}
	
	//2.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Y,&this->ImagSignalForFFT_Y,1);
	
	//2.5) Copy the image that has been convolved in the orginal image
	for (y = 0; y < SF->NY; y++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  SF->P(this->RealSignalForFFT_Y.G(0,y,0),x,y,z);
	}
      }
    }
    
    //3) convolution on Z axis
    for (y = 0; y < SF->NY; y++) for (x = 0; x < SF->NX; x++){
      //3.1) Copy the orginal image in the image that will be transformed
      
      //3.1.1) Init
      for (z=0;z<this->NZfft;z++) this->RealSignalForFFT_Z.P(0.,0,0,z);
      for (z=0;z<this->NZfft;z++) this->ImagSignalForFFT_Z.P(0.,0,0,z);
      
      Something=0;
      for (z = 0; z < SF->NZ; z++)   if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){ 
	this->RealSignalForFFT_Z.P(SF->G(x,y,z),0,0,z);
	Something=1;
      }
      
      
            
      if (Something==1){
	//3.1.2) check closest region boundary - forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.3) check closest region boundary - backward side
	PreviousInMask=-1;
	for (z = SF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=0;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if ((-this->RealSignalForFFT_Z.G(0,0,z)>static_cast<float>(MirrorInMask))||(fabs(this->RealSignalForFFT_Z.G(0,0,z))<0.001)) //closer than the boundary in the other direction or no boundary in the other direction
	      this->RealSignalForFFT_Z.P(static_cast<float>(MirrorInMask),0,0,z);
	  }
	}
	
	//3.1.4) mirror conditions outside of the mask -> forward side
	PreviousInMask=-1;
	for (z=0;z<this->NZfft;z++){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z+1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask--;
	    if (MirrorInMask<0) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)<0) this->RealSignalForFFT_Z.P(SF->G(x,y,MirrorInMask),0,0,z);
	  }
	}

	//3.1.5) mirror conditions outside of the mask -> backward side
	PreviousInMask=-1;
	for (z = SF->NZ-1; z >=0 ; z--){  
	  //In the mask
	  if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01) {PreviousInMask=z;MirrorInMask=z-1;}
	  
	  //Has been in the mask
	  if ((fabs(Mask->G(x,y,z)-MaskId_fl)>0.01)&&(PreviousInMask>-1)){
	    MirrorInMask++;
	    if (MirrorInMask>SF->NZ-1) MirrorInMask=PreviousInMask;
	    if (fabs(Mask->G(x,y,MirrorInMask)-MaskId_fl)>0.01) MirrorInMask=PreviousInMask;
	    if (this->RealSignalForFFT_Z.G(0,0,z)>0) this->RealSignalForFFT_Z.P(SF->G(x,y,MirrorInMask),0,0,z);
	  }
	}


	
	//3.2) Transform RealSignalForFFT and ImagSignalForFFT in Fourier spaces
	this->DirectFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.3) filtering in Fourier spaces
	CoefMult=sqrt((float)this->RealSignalForFFT_Z.NZ);
	
	for (z = 0; z < this->RealSignalForFFT_Z.NZ; z++){
	  a=this->RealSignalForFFT_Z.G(0, 0, z);
	  b=this->ImagSignalForFFT_Z.G(0, 0, z);
	  c=this->RealFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  d=this->ImagFilterForFFT_Z.G(0, 0, z)*CoefMult;
	  
	  this->RealSignalForFFT_Z.P(a*c-b*d,0,0,z);
	  this->ImagSignalForFFT_Z.P(c*b+a*d,0,0,z);
	}
	
	//3.4) IFFT
	this->InverseFFT(&this->RealSignalForFFT_Z,&this->ImagSignalForFFT_Z,2);
	
	//3.5) Copy the image that has been convolved in the orginal image
	for (z = 0; z < SF->NZ; z++) if (fabs(Mask->G(x,y,z)-MaskId_fl)<0.01){
	  SF->P(this->RealSignalForFFT_Z.G(0,0,z),x,y,z);
	}
      }
    }
}









///design a kernel that is the sum of up to 7 Gaussians and transform it in Fourier spaces
//if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
void LightFFTconvolver3D::MakeSumOf7AnisotropicGaussianFilters(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1,float weight2,float sigmaX2,float sigmaY2,float sigmaZ2,float weight3,float sigmaX3,float sigmaY3,float sigmaZ3,float weight4,float sigmaX4,float sigmaY4,float sigmaZ4,float weight5,float sigmaX5,float sigmaY5,float sigmaZ5,float weight6,float sigmaX6,float sigmaY6,float sigmaZ6,float weight7,float sigmaX7,float sigmaY7,float sigmaZ7,int NormalizeWeights){
  int k,x,y,z;
  float SumLoc;
  float weight,sigmaX,sigmaY,sigmaZ,sumWeight,CubicRootOfWeight;
  
  //1) set RealFilterForFFT and ImagFilterForFFT to 0 in case it contains something
  for (x=0;x<this->NXfft;x++) this->RealFilterForFFT_X.P(0.,x,0,0);
  for (y=0;y<this->NYfft;y++) this->RealFilterForFFT_Y.P(0.,0,y,0);
  for (z=0;z<this->NZfft;z++) this->RealFilterForFFT_Z.P(0.,0,0,z);
  
  for (x=0;x<this->NXfft;x++) this->ImagFilterForFFT_X.P(0.,x,0,0);
  for (y=0;y<this->NYfft;y++) this->ImagFilterForFFT_Y.P(0.,0,y,0);
  for (z=0;z<this->NZfft;z++) this->ImagFilterForFFT_Z.P(0.,0,0,z);
  
  sumWeight=0;
  
  //2) compute and save the 7 kernels
  for (k=0;k<7;k++){
    //parameters of the current kernel
    if (k==6){weight=weight7; sigmaX=sigmaX7; sigmaY=sigmaY7; sigmaZ=sigmaZ7;}
    if (k==5){weight=weight6; sigmaX=sigmaX6; sigmaY=sigmaY6; sigmaZ=sigmaZ6;}
    if (k==4){weight=weight5; sigmaX=sigmaX5; sigmaY=sigmaY5; sigmaZ=sigmaZ5;}
    if (k==3){weight=weight4; sigmaX=sigmaX4; sigmaY=sigmaY4; sigmaZ=sigmaZ4;}
    if (k==2){weight=weight3; sigmaX=sigmaX3; sigmaY=sigmaY3; sigmaZ=sigmaZ3;}
    if (k==1){weight=weight2; sigmaX=sigmaX2; sigmaY=sigmaY2; sigmaZ=sigmaZ2;}
    if (k==0){weight=weight1; sigmaX=sigmaX1; sigmaY=sigmaY1; sigmaZ=sigmaZ1;}
    
    CubicRootOfWeight=static_cast<float>(pow(static_cast<double>(weight),0.3333));
    sumWeight+=weight;
    
    //2.1) design the current kernel in Z direction...
    //...normalisation factor
    SumLoc=0.;
    for (z=0;z<this->NZfft/2;z++)
      SumLoc+=exp( -((float)(z*z))/(2.*sigmaZ*sigmaZ));
    
    for (z=this->NZfft/2;z<this->NZfft;z++)
      SumLoc+=exp( -((float)((this->NZfft-z)*(this->NZfft-z)))/(2.*sigmaZ*sigmaZ));
    
    //...fill the normalised values
    if (SumLoc>=0.01){
      for (z=0;z<this->NZfft/2;z++)
        this->RealFilterForFFT_Z.Add(CubicRootOfWeight*exp( -((float)(z*z)/(2.*sigmaZ*sigmaZ)))/SumLoc,0,0,z);
    
      for (z=this->NZfft/2;z<this->NZfft;z++)
        this->RealFilterForFFT_Z.Add(CubicRootOfWeight*exp( -((float)((this->NZfft-z)*(this->NZfft-z)))/(2.*sigmaZ*sigmaZ))/SumLoc,0,0,z);
    }
    else{
        cout << "Kernel on z axis has a problem" << endl;
      }
    
    
    //2.2) design the current kernel in Y direction
    //...normalisation factor
    SumLoc=0.;
    for (y=0;y<this->NYfft/2;y++)
      SumLoc+=exp( -((float)(y*y))/(2.*sigmaY*sigmaY));
    
    for (y=this->NYfft/2;y<this->NYfft;y++)
      SumLoc+=exp( -((float)((this->NYfft-y)*(this->NYfft-y)))/(2.*sigmaY*sigmaY));
    
    //...fill the normalised values
    if (SumLoc>=0.01){
      for (y=0;y<this->NYfft/2;y++)
	      this->RealFilterForFFT_Y.Add(CubicRootOfWeight*exp( -((float)(y*y))/(2.*sigmaY*sigmaY))/SumLoc,0,y,0);
      
      for (y=this->NYfft/2;y<this->NYfft;y++)
	      this->RealFilterForFFT_Y.Add(CubicRootOfWeight*exp( -((float)((this->NYfft-y)*(this->NYfft-y)))/(2.*sigmaY*sigmaY))/SumLoc,0,y,0);
    }
    else{
        cout << "Kernel on y axis has a problem" << endl;
      }

    
    //2.3) design the current kernel in X direction
    //...normalisation factor
    SumLoc=0.;
    for (x=0;x<this->NXfft/2;x++)
      SumLoc+=exp( -((float)(x*x))/(2.*sigmaX*sigmaX));
    
    for (x=this->NXfft/2;x<this->NXfft;x++)
      SumLoc+=exp( -((float)((this->NXfft-x)*(this->NXfft-x)))/(2.*sigmaX*sigmaX));

    //...fill the normalised values
    if (SumLoc>=0.01){
      for (x=0;x<this->NXfft/2;x++)
	this->RealFilterForFFT_X.Add(CubicRootOfWeight*exp( -((float)(x*x))/(2.*sigmaX*sigmaX))/SumLoc,x,0,0);
      
      for (x=this->NXfft/2;x<this->NXfft;x++)
	this->RealFilterForFFT_X.Add(CubicRootOfWeight*exp( -((float)((this->NXfft-x)*(this->NXfft-x)))/(2.*sigmaX*sigmaX))/SumLoc,x,0,0);
        }
    else{
        cout << "Kernel on x axis has a problem" << endl;
      }

  }
  
  //3) normalize the weights
  if (NormalizeWeights!=0){
    for (x=0;x<this->NXfft;x++) this->RealFilterForFFT_X.P(this->RealFilterForFFT_X.G(x,0,0)/sumWeight,x,0,0);
    for (y=0;y<this->NYfft;y++) this->RealFilterForFFT_Y.P(this->RealFilterForFFT_Y.G(0,y,0)/sumWeight,0,y,0);
    for (z=0;z<this->NZfft;z++) this->RealFilterForFFT_Z.P(this->RealFilterForFFT_Z.G(0,0,z)/sumWeight,0,0,z);
    
    //cout << "Note that the weights of the kernel are normalized -> play with the -VFpenalizer option to put more or less weight on the regularisation energy" << endl;
  }

  //4) Transform the RealFilterForFFT_. and ImagFilterForFFT_. in Fourier spaces
  this->DirectFFT(&this->RealFilterForFFT_X,&this->ImagFilterForFFT_X,0);
  this->DirectFFT(&this->RealFilterForFFT_Y,&this->ImagFilterForFFT_Y,1);
  this->DirectFFT(&this->RealFilterForFFT_Z,&this->ImagFilterForFFT_Z,2);

}


///Fast Fourier Transform
/// -> if axis == 0 -> FFT on X axis
/// -> if axis == 1 -> FFT on Y axis
/// -> if axis == 2 -> FFT on Z axis
void LightFFTconvolver3D::DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SizeY=RealSignal->NY;
  SizeZ=RealSignal->NZ;

  SqrtSizeX=sqrt(SizeX);
  SqrtSizeY=sqrt(SizeY);
  SqrtSizeZ=sqrt(SizeZ);
  
  
  //2) perform the fft along x axis
  if (axis==0){
    dataX = new float [SizeX*2+1];
    for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
      for (x = 0; x < SizeX; x++){
        dataX[2*x+1]=RealSignal->G(x, y, z);
        dataX[2*x+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataX, (unsigned long)SizeX, 1);
      for (x = 0; x < SizeX; x++){
        RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
        ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX, x, y, z);
      }
    }
    delete dataX;
  }
  
  //3) perform the fft along y axis
  if (axis==1){
    dataY = new float [SizeY*2+1];
    for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
      for (y = 0; y < SizeY; y++){
        dataY[2*y+1]=RealSignal->G(x, y, z);
        dataY[2*y+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataY, (unsigned long)SizeY, 1);
      for (y = 0; y < SizeY; y++){
        RealSignal->P(dataY[2*y+1]/SqrtSizeY,x, y, z);
        ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
      }
    }
    delete dataY;
  }
  
  //4) perform the fft along z axis
  if (axis==2){
    dataZ = new float [SizeZ*2+1];
    for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
      for (z = 0; z < SizeZ; z++){
        dataZ[2*z+1]=RealSignal->G(x, y, z);
        dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataZ, (unsigned long)SizeZ, 1);
      for (z = 0; z < SizeZ; z++){
        RealSignal->P(dataZ[2*z+1]/SqrtSizeZ,x, y, z);
        ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ, x, y, z);
      }
    }
    delete dataZ;
  }
}


///Inverse Fast Fourier Transform
/// -> if axis == 0 -> IFFT on X axis
/// -> if axis == 1 -> IFFT on Y axis
/// -> if axis == 2 -> IFFT on Z axis
void LightFFTconvolver3D::InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis){
  int SizeX,SizeY,SizeZ;
  float SqrtSizeX,SqrtSizeY,SqrtSizeZ;
  int x,y,z;
  float * dataX;
  float * dataY;
  float * dataZ;
  
  //1) extract the size of the images
  SizeX=RealSignal->NX;
  SqrtSizeX=sqrt(SizeX);
  
  SizeY=RealSignal->NY;
  SqrtSizeY=sqrt(SizeY);
  
  SizeZ=RealSignal->NZ;
  SqrtSizeZ=sqrt(SizeZ);
  
  
  //2) perform the ifft along z axis
  if (axis==2){
    dataZ = new float [SizeZ*2+1];
    for (y = 0; y < SizeY; y++) for (x = 0; x < SizeX; x++){
      for (z = 0; z < SizeZ; z++){
        dataZ[2*z+1]=RealSignal->G(x, y, z);
        dataZ[2*z+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataZ, (unsigned long)SizeZ, -1);
      for (z = 0; z < SizeZ; z++){
        RealSignal->P(dataZ[2*z+1]/SqrtSizeZ, x, y, z);
        ImaginarySignal->P(dataZ[2*z+2]/SqrtSizeZ,x, y, z);
      }
    }
    delete dataZ;
  }
  
  //3) perform the ifft along y axis
  if (axis==1){
    dataY = new float [SizeY*2+1];
    for (z = 0; z < SizeZ; z++) for (x = 0; x < SizeX; x++){
      for (y = 0; y < SizeY; y++){
        dataY[2*y+1]=RealSignal->G(x, y, z);
        dataY[2*y+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataY, (unsigned long)SizeY, -1);
      for (y = 0; y < SizeY; y++){
        RealSignal->P(dataY[2*y+1]/SqrtSizeY, x, y, z);
        ImaginarySignal->P(dataY[2*y+2]/SqrtSizeY, x, y, z);
      }
    }
    delete dataY;
  }
  
  //4) perform the ifft along x axis
  if (axis==0){
    dataX = new float [SizeX*2+1];
    for (z = 0; z < SizeZ; z++) for (y = 0; y < SizeY; y++){
      for (x = 0; x < SizeX; x++){
        dataX[2*x+1]=RealSignal->G(x, y, z);
        dataX[2*x+2]=ImaginarySignal->G(x, y, z);
      }
      four1NR(dataX, (unsigned long)SizeX, -1);
      for (x = 0; x < SizeX; x++){
        RealSignal->P(dataX[2*x+1]/SqrtSizeX, x, y, z);
        ImaginarySignal->P(dataX[2*x+2]/SqrtSizeX,x, y, z);
      }
    }
    delete dataX;
  }
}  


///Fast Fourier Transform of numerical recipies (slighly modified)
void LightFFTconvolver3D::four1NR(float data[], unsigned long nn, int isign){
  unsigned long n,mmax,m,j,istep,i;
  double wtemp,wr,wpr,wpi,wi,theta;
  float tempr,tempi;
  
  n=nn << 1;
  j=1;
  for (i=1;i<n;i+=2){
    if (j>i){
      tempr=data[j]; data[j]=data[i]; data[i]=tempr;
      tempr=data[j+1]; data[j+1]=data[i+1]; data[i+1]=tempr;
    }
    m=n >> 1;
    while ((m>=2) && (j>m)){
      j -= m;
      m >>= 1;
    }
    j += m;
  }
  mmax=2;
  while (n > mmax) {
    istep=mmax << 1;
    theta=isign*(6.28318530717959/mmax);
    wtemp=sin(0.5*theta);
    wpr = -2.0*wtemp*wtemp;
    wpi=sin(theta);
    wr=1.0;
    wi=0.0;
    for (m=1;m<mmax;m+=2) {
      for (i=m;i<=n;i+=istep) {
        j=i+mmax;
        tempr=wr*data[j]-wi*data[j+1];
        tempi=wr*data[j+1]+wi*data[j];
        data[j]=data[i]-tempr;
        data[j+1]=data[i+1]-tempi;
        data[i] += tempr;
        data[i+1] += tempi;
      }
      wr=(wtemp=wr)*wpr-wi*wpi+wr;
      wi=wi*wpr+wtemp*wpi+wi;
    }
    mmax=istep;
  }
}



///4.4: class where we consider N regions to smooth


///Constructor
MultiRegionFFTConvolver2::MultiRegionFFTConvolver2(){
  NumberOfRegions=0;
}

///destructor
MultiRegionFFTConvolver2::~MultiRegionFFTConvolver2(){
}



//Initiate the convolver in all regions using same kernel
//-> 'Part_Of_Unity' is a 3D mask which define different ROIs. It only contains integer values each of them associated to a ROI. Partition of unity is first 
//   defined by splitting this mask into several channels, each of them having an intensity equal to 1 in the corresponding ROI and 0 otherwise. Each channel
//   is then smoothed with a Gaussian kernel of stddev 'sigmaPOI'.
//-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
//     * NX, NY, NZ: is the size of the input image
//     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
//     * ...
//     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void MultiRegionFFTConvolver2::InitiateConvolver(ScalarField * Part_Of_Unity, float sigmaPOI, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){
  int i,j;
  float tmpFl,maxSum,minSum,maxValue;
  int NBX,NBY,NBZ;
  int z,y,x;
  float * RegionsFound;
  int RegionsNb;
  float epsilon;
  int tmpInt;
  LightFFTconvolver3D TmpLightFFTconvolver;
  float x_mm,y_mm,z_mm;
  
  NBX=Part_Of_Unity->NX;
  NBY=Part_Of_Unity->NY;
  NBZ=Part_Of_Unity->NZ;
  
  //1) find the number of regions and check that all values are integers
  RegionsFound = new float [100];
  RegionsNb=0;
  tmpInt=1;
  epsilon=0.000001;
  
  //1.1) detect if the mask is an actual mask
  for(z=0;z<NBZ;z++) for(y=0;y<NBY;y++) for(x=0;x<NBX;x++) if (fabs(Part_Of_Unity->G(x,y,z))>epsilon)
      if (fabs(1-(floor(Part_Of_Unity->G(x,y,z))/Part_Of_Unity->G(x,y,z)))>epsilon)
        tmpInt=0;
  
  if ((Part_Of_Unity->NT>1)||(tmpInt==0)){
    cout << endl;
    cout << "WARNING: The mask used to generate the partition of unity has not only integer values and/or has several channels -> considered as the actual partition of unity" << endl;
    cout << endl;
    this->InitiateConvolverWithActualPOI(Part_Of_Unity,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7);
    return;
    }
  
  //1.2) find the number of regions and their identifier
  for(z=0;z<NBZ;z++) for(y=0;y<NBY;y++) for(x=0;x<NBX;x++) {
    tmpInt=0;
    for (i=0;i<RegionsNb;i++) if (fabs(Part_Of_Unity->G(x,y,z)-RegionsFound[i])<epsilon) tmpInt=1;
    
    if ((tmpInt==0)&&(RegionsNb<99)){
      RegionsFound[RegionsNb]=Part_Of_Unity->G(x,y,z);
      RegionsNb++;
      }
    }
    
    if (RegionsNb==99)
      cout << "Too much regions found (more than 99!). Further results may be false." << endl;
  
  this->NumberOfRegions=RegionsNb;
  
  
  for (i=0;i<RegionsNb-1;i++) for (j=i+1;j<RegionsNb;j++) if (RegionsFound[i]>RegionsFound[j]){
      tmpFl=RegionsFound[i];
      RegionsFound[i]=RegionsFound[j];
      RegionsFound[j]=tmpFl;
    }
  
  cout << endl;
  cout << RegionsNb << " regions found: "; for (i=0;i<RegionsNb;i++) cout << RegionsFound[i] << " "; cout << endl; 
  cout << endl;
  
  
  //2) initiate and treat the partition of unity
  
  this->PartitionOfUnity.CreateVoidField(NBX,NBY,NBZ,this->NumberOfRegions);
  
  for(z=0;z<NBZ;z++) for(y=0;y<NBY;y++) for(x=0;x<NBX;x++){
    //2.1) find the local ROI
    tmpInt=0;
    for (i=0;i<RegionsNb;i++) if (fabs(Part_Of_Unity->G(x,y,z)-RegionsFound[i])<epsilon) tmpInt=i;
    
    //2.2) set it in the partition of unity
    this->PartitionOfUnity.P(1,x,y,z,tmpInt);
  }
  
  
  //2.3) smooth the partition of unity
  x_mm=sqrt(this->PartitionOfUnity.Image2World[0][0]*this->PartitionOfUnity.Image2World[0][0]+this->PartitionOfUnity.Image2World[0][1]*this->PartitionOfUnity.Image2World[0][1]+this->PartitionOfUnity.Image2World[0][2]*this->PartitionOfUnity.Image2World[0][2]);
  y_mm=sqrt(this->PartitionOfUnity.Image2World[1][0]*this->PartitionOfUnity.Image2World[1][0]+this->PartitionOfUnity.Image2World[1][1]*this->PartitionOfUnity.Image2World[1][1]+this->PartitionOfUnity.Image2World[1][2]*this->PartitionOfUnity.Image2World[1][2]);
  z_mm=sqrt(this->PartitionOfUnity.Image2World[2][0]*this->PartitionOfUnity.Image2World[2][0]+this->PartitionOfUnity.Image2World[2][1]*this->PartitionOfUnity.Image2World[2][1]+this->PartitionOfUnity.Image2World[2][2]*this->PartitionOfUnity.Image2World[2][2]);
 
  TmpLightFFTconvolver.InitiateConvolver(NBX,NBY,NBZ,1,sigmaPOI/x_mm,sigmaPOI/y_mm,sigmaPOI/z_mm);
  TmpLightFFTconvolver.Convolution(&this->PartitionOfUnity);
  
  
  //3) initiate the region convolvers
  this->Region_convolver = new LightFFTconvolver3D [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++)
    this->Region_convolver[i].InitiateConvolver(NBX,NBY,NBZ,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,0);
  
  //4) initiate the temporary vector fields
  this->TempVF1.CreateVoidField(NBX,NBY,NBZ);
  this->TempVF2.CreateVoidField(NBX,NBY,NBZ);

  //5) initiate the ROI containing information in each region
  xmin= new int [this->NumberOfRegions];
  xmax= new int [this->NumberOfRegions];
  ymin= new int [this->NumberOfRegions];
  ymax= new int [this->NumberOfRegions];
  zmin= new int [this->NumberOfRegions];
  zmax= new int [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++) xmin[i]=NBX-1;
  for (i=0;i<this->NumberOfRegions;i++) ymin[i]=NBY-1;
  for (i=0;i<this->NumberOfRegions;i++) zmin[i]=NBZ-1;
  for (i=0;i<this->NumberOfRegions;i++) xmax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) ymax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) zmax[i]=1;
  
  for (i=0;i<this->NumberOfRegions;i++){
    maxValue=0;
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++)
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue) maxValue=this->PartitionOfUnity.G(x,y,z,i);
    
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue/20){
	      if (xmin[i]>x) xmin[i]=x;
	      if (ymin[i]>y) ymin[i]=y;
	      if (zmin[i]>z) zmin[i]=z;
      	
	      if (xmax[i]<x) xmax[i]=x+1;
	      if (ymax[i]<y) ymax[i]=y+1;
	      if (zmax[i]<z) zmax[i]=z+1;
      }
    }
    
    cout << "ROI containing region " << i << ": ";
    cout << "x=["<< xmin[i] << "," << xmax[i] << "], ";
    cout << "y=["<< ymin[i] << "," << ymax[i] << "], ";
    cout << "z=["<< zmin[i] << "," << zmax[i] << "]" << endl;
  }
}


///Initiate the convolver in all regions using same kernel
///-> Part_Of_Unity is a '3D + channels' scalar field which encodes the partition of unity in the different channels.
///   * Its size and number of channels (NBT actually) defines the size and the number of regions of the convolver 
///   * The maximum point-wise sum of the probabilities may be different to 1: normalisation will be automatically performed
///   * Point-wise sum of the probabilities may vary in space. If so, a background region will be automatically defined
///-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
///     * NX, NY, NZ: is the size of the input image
///     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
///     * ...
///     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
void MultiRegionFFTConvolver2::InitiateConvolverWithActualPOI(ScalarField * Part_Of_Unity, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){
  int i;
  float tmpFl,maxSum,minSum,maxValue;
  int NBX,NBY,NBZ;
  int z,y,x;
  
  //1) initiate the size of the convolver and the number of regions
  NBX=Part_Of_Unity->NX;
  NBY=Part_Of_Unity->NY;
  NBZ=Part_Of_Unity->NZ;
  this->NumberOfRegions=Part_Of_Unity->NT;
  
  //2) initiate and treat the partition of unity
  
  //2.1) check whether an additional 'background channel' is necessary
  for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
    
    tmpFl=0;
    
    for (i=0;i<this->NumberOfRegions;i++) tmpFl+=Part_Of_Unity->G(x,y,z,i);
    
    if ((x==0)&&(y==0)&&(z==0)){ maxSum=tmpFl; minSum=tmpFl; }
    if (maxSum<tmpFl) maxSum=tmpFl;
    if (minSum>tmpFl) minSum=tmpFl;
  }
  
  if ((maxSum-minSum)/maxSum>0.9){ //more than 10% difference between the point-wise partions of unity
    this->NumberOfRegions++;  // => a 'background channel' will be generated
    cout << "+++ A background channel is created +++" << endl;
  }
  
  //2.2) generate and fill the partition of unity
  this->PartitionOfUnity.CreateVoidField(NBX,NBY,NBZ,this->NumberOfRegions);
  
  for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
    tmpFl=0;
    for (i=0;i<Part_Of_Unity->NT;i++){
      this->PartitionOfUnity.P(Part_Of_Unity->G(x,y,z,i)/maxSum,x,y,z,i);
      tmpFl+=Part_Of_Unity->G(x,y,z,i)/maxSum;
    }
    
    if (Part_Of_Unity->NT<this->NumberOfRegions){
       this->PartitionOfUnity.P(1-tmpFl,x,y,z,this->NumberOfRegions-1);
    }
  }
  
  //2) initiate the region convolvers
  this->Region_convolver = new LightFFTconvolver3D [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++)
    this->Region_convolver[i].InitiateConvolver(NBX,NBY,NBZ,w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,0);
  
  //3) initiate the temporary vector fields
  this->TempVF1.CreateVoidField(NBX,NBY,NBZ);
  this->TempVF2.CreateVoidField(NBX,NBY,NBZ);

  //4) initiate the ROI containing information in each region
  xmin= new int [this->NumberOfRegions];
  xmax= new int [this->NumberOfRegions];
  ymin= new int [this->NumberOfRegions];
  ymax= new int [this->NumberOfRegions];
  zmin= new int [this->NumberOfRegions];
  zmax= new int [this->NumberOfRegions];
  
  for (i=0;i<this->NumberOfRegions;i++) xmin[i]=NBX-1;
  for (i=0;i<this->NumberOfRegions;i++) ymin[i]=NBY-1;
  for (i=0;i<this->NumberOfRegions;i++) zmin[i]=NBZ-1;
  for (i=0;i<this->NumberOfRegions;i++) xmax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) ymax[i]=1;
  for (i=0;i<this->NumberOfRegions;i++) zmax[i]=1;
  
  for (i=0;i<this->NumberOfRegions;i++){
    maxValue=0;
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++)
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue) maxValue=this->PartitionOfUnity.G(x,y,z,i);
    
    for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++){
      if (this->PartitionOfUnity.G(x,y,z,i)>maxValue/20){
	      if (xmin[i]>x) xmin[i]=x;
	      if (ymin[i]>y) ymin[i]=y;
	      if (zmin[i]>z) zmin[i]=z;
      	
	      if (xmax[i]<x) xmax[i]=x+1;
	      if (ymax[i]<y) ymax[i]=y+1;
	      if (zmax[i]<z) zmax[i]=z+1;
      }
    }
    
    cout << "ROI containing region " << i << ": ";
    cout << "x=["<< xmin[i] << "," << xmax[i] << "], ";
    cout << "y=["<< ymin[i] << "," << ymax[i] << "], ";
    cout << "z=["<< zmin[i] << "," << zmax[i] << "]" << endl;
  }
}

///save the actual partition of unity (after potential treatments in InitiateConvolver or undersampling)
///-> 1st char* is the name in which the image is saved
///-> 2nd char* is the name of the image that will be used in the header of the saved image
void MultiRegionFFTConvolver2::SaveActualParitionOfUnity(char * OutputImageName, char * ImageForHeaderName){
  this->PartitionOfUnity.Write(OutputImageName,ImageForHeaderName);
  }


///Update the parition of unity 
///-> it must have the same size and number of layers/times as in the current POU (tested)
///-> to make sense, the new POU must have a sum of intensities across layers/times equals to 1 at each voxel (not tested)
void MultiRegionFFTConvolver2::UpdatePartitionOfUnity(ScalarField * Part_Of_Unity){
  int z,y,x,t,i;
  
  //1) Check the size of the POU and nb of regions
  if (this->PartitionOfUnity.NX!=Part_Of_Unity->NX){cout << "Partition of unity not changed" << endl; return;}
  if (this->PartitionOfUnity.NY!=Part_Of_Unity->NY){cout << "Partition of unity not changed" << endl; return;}
  if (this->PartitionOfUnity.NZ!=Part_Of_Unity->NZ){cout << "Partition of unity not changed" << endl; return;}
  if (this->PartitionOfUnity.NT!=Part_Of_Unity->NT){cout << "Partition of unity not changed" << endl; return;}
  
  //2) fill the partition of unity
  for(t=0;t<Part_Of_Unity->NT;t++) for(z=0;z<Part_Of_Unity->NZ;z++) for(y=0;y<Part_Of_Unity->NY;y++) for(x=0;x<Part_Of_Unity->NX;x++)
      this->PartitionOfUnity.P(Part_Of_Unity->G(x,y,z,t),x,y,z,t);
  
  //3) change region-based ROIs  (maybe not optimal but it is not re-estimated, so some time is gained here for sure)
  for (i=0;i<this->NumberOfRegions;i++){
    this->xmin[i]=0;   this->ymin[i]=0;  this->zmin[i]=0;
    this->xmax[i]=Part_Of_Unity->NX;  this->ymax[i]=Part_Of_Unity->NY;  this->zmax[i]=Part_Of_Unity->NZ;
  }
}


///change the smoothing kernel in one region
void MultiRegionFFTConvolver2::ChangeKernelInOneRegion(int IdRegion, float w1,float sX1,float sY1,float sZ1, float w2,float sX2,float sY2,float sZ2, float w3,float sX3,float sY3,float sZ3, float w4,float sX4,float sY4,float sZ4, float w5,float sX5,float sY5,float sZ5, float w6,float sX6,float sY6,float sZ6, float w7,float sX7,float sY7,float sZ7){

  this->Region_convolver[IdRegion].ChangeKernel(w1,sX1,sY1,sZ1,w2,sX2,sY2,sZ2,w3,sX3,sY3,sZ3,w4,sX4,sY4,sZ4,w5,sX5,sY5,sZ5,w6,sX6,sY6,sZ6,w7,sX7,sY7,sZ7,0);

}



///convolution of a 3D vector field using the predifined kernel
void MultiRegionFFTConvolver2::Convolution(VectorField * VF){
  int z,y,x,l;
  int IdRegion;
  
  //1) init
  this->TempVF2.PutToAllVoxels(0);
  
  //2) multi-region convolution...
  for (IdRegion=0;IdRegion<this->NumberOfRegions;IdRegion++){
    //... init temp VF
    DeepCopy(VF,&this->TempVF1);
    
    //... multiply by the PartitionOfUnity
    for(z=0;z<this->TempVF1.NZ;z++) for(y=0;y<this->TempVF1.NY;y++) for(x=0;x<this->TempVF1.NX;x++) for (l=0;l<3;l++){
      this->TempVF1.P(this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,x,y,z);
    }
    
    //... convolve with local convolver
    this->Region_convolver[IdRegion].ConvolutionInROI(&this->TempVF1,this->xmin[IdRegion],this->xmax[IdRegion],this->ymin[IdRegion],this->ymax[IdRegion],this->zmin[IdRegion],this->zmax[IdRegion]);
    
    //... update the final velocity field 
    for(z=this->zmin[IdRegion];z<this->zmax[IdRegion];z++) for(y=this->ymin[IdRegion];y<this->ymax[IdRegion];y++) for(x=this->xmin[IdRegion];x<this->xmax[IdRegion];x++) for (l=0;l<3;l++){
      this->TempVF2.Add(this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,x,y,z);
      
      //+++ begin hack for symmetry in one region +++
      /*int SymmetricRegion;
      float FactorSymmetricRegion;
      
      SymmetricRegion=0;
      FactorSymmetricRegion=0.1;
      
      if (IdRegion==SymmetricRegion){
        //if (l==0) this->TempVF2.Add(0.5*FactorSymmetricRegion*this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,x,y,z);
        //if (l==1) this->TempVF2.Add(-0.5*FactorSymmetricRegion*this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,x,y,z);
        //if (l==2) this->TempVF2.Add(-0.5*FactorSymmetricRegion*this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,x,y,z);

        if (l==0) this->TempVF2.Add(-FactorSymmetricRegion*this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,this->TempVF2.NX-1-x,y,z);
        if (l==1) this->TempVF2.Add(FactorSymmetricRegion*this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,this->TempVF2.NX-1-x,y,z);
        if (l==2) this->TempVF2.Add(FactorSymmetricRegion*this->TempVF1.G(l,x,y,z)*this->PartitionOfUnity.G(x,y,z,IdRegion),l,this->TempVF2.NX-1-x,y,z);
        }*/
      //+++  end hack for symmetry in one region +++
    }
  }
  
  //3) copy the temporary smoothed VF
  DeepCopy(&this->TempVF2,VF);

}



///return the number of regions considered
int MultiRegionFFTConvolver2::GetRegionNb(){
  return this->NumberOfRegions;
}





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                  5: CLASS TO MANAGE THE MUTUAL INFORMATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///Constructor
MImanager::MImanager(){
	this->NX=0;
	this->NY=0;
	this->NZ=0;
  this->NumberOfBinsS=0;
  this->NumberOfBinsT=0;
  this->MinGreyLevelsS=0;
  this->MinGreyLevelsT=0;
  this->SizeStepsGreyLevelsS=0; 
  this->SizeStepsGreyLevelsT=0;
  this->MI=0;
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
  this->JointEntropy=0;
  this->MarginalEntropyS=0;
  this->MarginalEntropyT=0;
  this->indicatorMaskDefined=0;
  
}

///destructor
MImanager::~MImanager(void){}


//Initiate the MI manager without any mask
void MImanager::Initiate(ScalarField * SourceImage,ScalarField * TargetImage,int NbBinsSrc,int NbBinsTrg, int LocMargin){
	int i;
  int x,y,z;
  float minval,maxval;
  
  this->SrcImage=SourceImage;
  this->TrgImage=TargetImage;
  
  this->indicatorMaskDefined=0;
  
  //number of bins and margin
  this->NumberOfBinsS=NbBinsSrc+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->NumberOfBinsT=NbBinsTrg+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->Margin=LocMargin;
  
  //allocate the histograms
  this->MarginalHistogramS= new float [this->NumberOfBinsS];
  this->MarginalHistogramT= new float [this->NumberOfBinsT];
  this->JointHistogram= new float * [this->NumberOfBinsS];
  for (i=0; i<this->NumberOfBinsS; i++) 
    this->JointHistogram[i]= new float [this->NumberOfBinsT];
  
  //put to 0 the 'up-to-date indicators' (to allow the first estimations of the histograms and the MI)
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
  
  //save the size of the treated images
  this->NX=SrcImage->NX;
  this->NY=SrcImage->NY;
  this->NZ=SrcImage->NZ;
  
  //compute the minimal and maximal intensities of the source image ...
  minval=SrcImage->G(0,0,0);
  maxval=SrcImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>SrcImage->G(x,y,z)) minval=SrcImage->G(x,y,z);
      if (maxval<SrcImage->G(x,y,z)) maxval=SrcImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>SrcImage->G(x,y,0)) minval=SrcImage->G(x,y,0);
      if (maxval<SrcImage->G(x,y,0)) maxval=SrcImage->G(x,y,0);
    }
  }
  MinGreyLevelsS=minval;
  SizeStepsGreyLevelsS=(maxval*1-minval)/(this->NumberOfBinsS-4);
  
  //compute the minimal and maximal intensities of the target image
  minval=TrgImage->G(0,0,0);
  maxval=TrgImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>TrgImage->G(x,y,z)) minval=TrgImage->G(x,y,z);
      if (maxval<TrgImage->G(x,y,z)) maxval=TrgImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      if (minval>TrgImage->G(x,y,0)) minval=TrgImage->G(x,y,0);
      if (maxval<TrgImage->G(x,y,0)) maxval=TrgImage->G(x,y,0);
    }
  }
  MinGreyLevelsT=minval;
  SizeStepsGreyLevelsT=(maxval*1-minval)/(this->NumberOfBinsT-4); 
}


//Initiate the MI manager with a mask. The MI and MI gradients will only be computed where the mask equals 1
void MImanager::Initiate(ScalarField * SourceImage,ScalarField * TargetImage,ScalarField * ROI_Mask,int NbBinsSrc,int NbBinsTrg, int LocMargin){
	int i;
  int x,y,z;
  float minval,maxval;
  
  //load the images and the mask
  this->SrcImage=SourceImage;
  this->TrgImage=TargetImage;
  this->Mask=ROI_Mask;
  
  //number of bins and margin
  this->NumberOfBinsS=NbBinsSrc+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->NumberOfBinsT=NbBinsTrg+4; //add two bins at the low intensities and two bins at the high ones to manage the boundary conditions
  this->Margin=LocMargin;
  
  //save the size of the treated images
  this->NX=SrcImage->NX;
  this->NY=SrcImage->NY;
  this->NZ=SrcImage->NZ;
  
  //check that the mask has at least 30 values equals to 1. If not, initiate the MImanager without a mask
  this->indicatorMaskDefined=0;
  for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++)
    if (fabs(Mask->G(x,y,z)-1)<0.0001) 
      this->indicatorMaskDefined++;
  
  if (this->indicatorMaskDefined<30){
    cout << "The ROI in the mask (where intensities equal 1) is too small to be considered -> MImanager initiated without a mask" << endl;
    this->Initiate(SourceImage,TargetImage,NbBinsSrc,NbBinsTrg,LocMargin);
    this->indicatorMaskDefined=0;
    return;
  }
  else{
    this->indicatorMaskDefined=1;
  }
  
  //allocate the histograms
  this->MarginalHistogramS= new float [this->NumberOfBinsS];
  this->MarginalHistogramT= new float [this->NumberOfBinsT];
  this->JointHistogram= new float * [this->NumberOfBinsS];
  for (i=0; i<this->NumberOfBinsS; i++) 
    this->JointHistogram[i]= new float [this->NumberOfBinsT];
  
  //put to 0 the 'up-to-date indicators' (to allow the first estimations of the histograms and the MI)
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
  
  
  //compute the minimal and maximal intensities of the source image ...
  minval=SrcImage->G(0,0,0);
  maxval=SrcImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,z)-1)<0.0001){
      if (minval>SrcImage->G(x,y,z)) minval=SrcImage->G(x,y,z);
      if (maxval<SrcImage->G(x,y,z)) maxval=SrcImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,0)-1)<0.0001){
      if (minval>SrcImage->G(x,y,0)) minval=SrcImage->G(x,y,0);
      if (maxval<SrcImage->G(x,y,0)) maxval=SrcImage->G(x,y,0);
    }
  }
  MinGreyLevelsS=minval;
  SizeStepsGreyLevelsS=(maxval*1-minval)/(this->NumberOfBinsS-4);
  
  //compute the minimal and maximal intensities of the target image
  minval=TrgImage->G(0,0,0);
  maxval=TrgImage->G(0,0,0);
  if (this->NZ>1){ //... 3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,z)-1)<0.0001){
      if (minval>TrgImage->G(x,y,z)) minval=TrgImage->G(x,y,z);
      if (maxval<TrgImage->G(x,y,z)) maxval=TrgImage->G(x,y,z);
    }
  }
  else{ //... 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,0)-1)<0.0001){
      if (minval>TrgImage->G(x,y,0)) minval=TrgImage->G(x,y,0);
      if (maxval<TrgImage->G(x,y,0)) maxval=TrgImage->G(x,y,0);
    }
  }
  MinGreyLevelsT=minval;
  SizeStepsGreyLevelsT=(maxval*1-minval)/(this->NumberOfBinsT-4); 
}






//Indicate to the MI manager that the Source image has changed
void MImanager::IndicateSrcHasChanged(){
  this->indicatorUpdatedSrcHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
}

//Indicate to the MI manager that the Target image has changed
void MImanager::IndicateTrgHasChanged(){
  this->indicatorUpdatedTrgHisto=0;
  this->indicatorUpdatedJointHisto=0;
  this->indicatorUpdatedMI=0;
}


//normalized the intensities adapted to the bins of the histograms
float MImanager::GiveFloatBinT(float intensity){
  return ((intensity-MinGreyLevelsT)/SizeStepsGreyLevelsT)+2;
}

float MImanager::GiveFloatBinS(float intensity){
  return ((intensity-MinGreyLevelsS)/SizeStepsGreyLevelsS)+2;
}


//when computing the histograms give the contribution of 'intensity' in the bin 'IdBin' for the source image
float MImanager::GiveValueParzenWindowT(float intensity,int IdBin){
  float GFBi,flId,tmp,tmp2;
  
  GFBi=GiveFloatBinT(intensity);
  flId=static_cast<float>(IdBin);
  
  if ((flId-2<GFBi)&&(GFBi<=flId-1)){
    return pow(GFBi-flId+2,3)/6;
  }
  else if ((flId-1<GFBi)&&(GFBi<=flId)){
    tmp=GFBi-flId+1;
    tmp2=tmp*tmp;
    return (-3*tmp2*tmp+3*tmp2+3*tmp+1)/6;
  }
  else if ((flId<GFBi)&&(GFBi<=flId+1)){
    tmp=GFBi-flId;
    tmp2=tmp*tmp;
    return (3*tmp2*tmp-6*tmp2+4)/6;
  }
  else if ((flId+1<GFBi)&&(GFBi<flId+2)){
    return pow(flId-GFBi+2,3)/6;
  }
  else{
    return 0;
  }
  
  //TO KEEP (fast approximation which can be inline)
  //return (0.5*(2-fabs(GiveFloatBinT(intensity)-static_cast<float>(IdBin))))*static_cast<float>(fabs(GiveFloatBinT(intensity)-static_cast<float>(IdBin))<2);
}


//when computing the histograms give the contribution of 'intensity' in the bin 'IdBin' for the target image
float MImanager::GiveValueParzenWindowS(float intensity,int IdBin){
  float GFBi,flId,tmp,tmp2;
  
  GFBi=GiveFloatBinS(intensity);
  flId=static_cast<float>(IdBin);
  
  if ((flId-2<GFBi)&&(GFBi<=flId-1)){
    return pow(GFBi-flId+2,3)/6;
  }
  else if ((flId-1<GFBi)&&(GFBi<=flId)){
    tmp=GFBi-flId+1;
    tmp2=tmp*tmp;
    return (-3*tmp2*tmp+3*tmp2+3*tmp+1)/6;
  }
  else if ((flId<GFBi)&&(GFBi<=flId+1)){
    tmp=GFBi-flId;
    tmp2=tmp*tmp;
    return (3*tmp2*tmp-6*tmp2+4)/6;
  }
  else if ((flId+1<GFBi)&&(GFBi<flId+2)){
    return pow(flId-GFBi+2,3)/6;
  }
  else{
    return 0;
  }
  
  
  //TO KEEP (fast approximation which can be inline)
  //return (0.5*(2-fabs(GiveFloatBinS(intensity)-static_cast<float>(IdBin))))*static_cast<float>(fabs(GiveFloatBinS(intensity)-static_cast<float>(IdBin))<2);
}



//Compute Joint Histogram And Entropy
void MImanager::ComputeJointHistogramAndEntropy(){
  int i,j;
  int x,y,z;
  int Total;
  int MinBinS,MinBinT;
  
  //histo
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    JointHistogram[i][j]=0.01;  //0.01 to avoid problems when computing the 
  
  
  
  if (this->indicatorMaskDefined==1){ // a mask is defined
    if (this->NZ>1){// ... 3D image
      for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,z)-1)<0.0001){
	//evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,z)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,z)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,z),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,z),j);
      }
    }
    else{// ... 2D image
      for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++) if (fabs(Mask->G(x,y,0)-1)<0.0001){
        //evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,0)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,0)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,0),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,0),j);
      }
    }
  }
  else{// no mask is defined
    if (this->NZ>1){// ... 3D image
      for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
        //evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,z)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,z)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,z),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,z),j);
      }
    }
    else{// ... 2D image
      for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
        //evaluate the bins in which the values of the parzen window are not null
        MinBinT=static_cast<int>(this->GiveFloatBinT(this->TrgImage->G(x,y,0)))-1;
        MinBinS=static_cast<int>(this->GiveFloatBinS(this->SrcImage->G(x,y,0)))-1;
        
        if (MinBinT<0) MinBinT=0;
        if (MinBinS<0) MinBinS=0;
        if (MinBinT>=NumberOfBinsT-3) MinBinT=NumberOfBinsT-4;
        if (MinBinS>=NumberOfBinsS-3) MinBinS=NumberOfBinsS-4;
      
        for (i = MinBinS; i < MinBinS+4; i++) for (j = MinBinT; j < MinBinT+4; j++)
          JointHistogram[i][j]+=GiveValueParzenWindowS(this->SrcImage->G(x,y,0),i)*GiveValueParzenWindowT(this->TrgImage->G(x,y,0),j);
      }
    }
  }
  
  Total=0;
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    Total+=JointHistogram[i][j];
  
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    JointHistogram[i][j]/=Total;
  
  //entropy
  this->JointEntropy=0;
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    JointEntropy+=JointHistogram[i][j]*log(JointHistogram[i][j]);
  
  //indicator
  indicatorUpdatedJointHisto=1;
}



//Compute Marginal Histogram And Entropy S
void MImanager::ComputeMarginalHistogramAndEntropyS(){
  int i,j;
  
  if (indicatorUpdatedJointHisto==0)
    this->ComputeJointHistogramAndEntropy();
  
  //histo
  for (i = 0; i < this->NumberOfBinsS; i++)
    MarginalHistogramS[i]=0;  
  
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    MarginalHistogramS[i]+=JointHistogram[i][j];
  
  //entropy
  this->MarginalEntropyS=0;
  for (i = 0; i < this->NumberOfBinsS; i++)
    MarginalEntropyS+=MarginalHistogramS[i]*log(MarginalHistogramS[i]);
  
  
  //indicator
  indicatorUpdatedSrcHisto=1;
}


//Compute Marginal Histogram And Entropy T
void MImanager::ComputeMarginalHistogramAndEntropyT(){
  int i,j;
  
  if (indicatorUpdatedJointHisto==0)
    this->ComputeJointHistogramAndEntropy();
  
  //histo
  for (i = 0; i < this->NumberOfBinsT; i++)
    MarginalHistogramT[i]=0;  
  
  for (i = 0; i < this->NumberOfBinsS; i++) for (j = 0; j < this->NumberOfBinsT; j++)
    MarginalHistogramT[j]+=JointHistogram[i][j];
  
  //entropy
  this->MarginalEntropyT=0;
  for (i = 0; i < this->NumberOfBinsT; i++)
    MarginalEntropyT+=MarginalHistogramT[i]*log(MarginalHistogramT[i]);
  
  
  //indicator
  indicatorUpdatedTrgHisto=1;
}  



//returns the normalized mutual information (plus update all the histograms)
float MImanager::EvaluateMI(){

  //compute the joint histogram and joint entropy...
  if ((indicatorUpdatedMI!=1)||(indicatorUpdatedJointHisto!=1)||(indicatorUpdatedTrgHisto!=1)||(indicatorUpdatedSrcHisto!=1))
    this->ComputeJointHistogramAndEntropy();
  
  
  //compute the histgram related to the source image and related marginal entropy...
  if ((indicatorUpdatedMI!=1)||(indicatorUpdatedSrcHisto!=1))
    ComputeMarginalHistogramAndEntropyS();
  
  //compute the histgram related to the target image and related marginal entropy...
  if ((indicatorUpdatedMI!=1)||(indicatorUpdatedTrgHisto!=1))
    ComputeMarginalHistogramAndEntropyT();
  
  //compute the normalised mutual information
  if ((indicatorUpdatedMI!=1)){
    this->MI=JointEntropy-MarginalEntropyT-MarginalEntropyS;  //mutual information actually
    
    indicatorUpdatedMI=1;
  }
  
  return this->MI;
  
}



//returns the estimated gradient of normalized mutual information
//In practice: 
//  * local intensity gradients of the deformed source image are first evaluated
//  * At each point, the demons-style update is computed and return in 'Gradient'
void MImanager::EvaluateGradMI(VectorField * Gradient){
  int i, j, k;
  int x,y,z;
  float LocNorm,direcX,direcY,direcZ;
  float tmpFl,tmpX,tmpY,tmpZ;
  int SizeNgbh,SizeCheckForwardBackward;
  float S_loc,S_fwd,S_bwd,Sprime;
  float epsilon;
  float * ListIntensities_Target;
  float * ListIntensities_Centered;
  float * ListIntensities_Forward;
  float * ListIntensities_Backward;
  
  int * MinBin_Target;
  int * MinBin_Centered;
  int * MinBin_Forward;
  int * MinBin_Backward;
  
  float lambdaX;
  
  lambdaX=1;
  
  //1) Initiate and allocate
  epsilon=0.0001;
  SizeNgbh=3;
  SizeCheckForwardBackward=1;
  
  ListIntensities_Target= new float [SizeNgbh];
  ListIntensities_Centered = new float [SizeNgbh];
  ListIntensities_Forward = new float [SizeNgbh];
  ListIntensities_Backward = new float [SizeNgbh];
  
  MinBin_Target = new int [SizeNgbh];
  MinBin_Centered = new int [SizeNgbh];
  MinBin_Forward = new int [SizeNgbh];
  MinBin_Backward = new int [SizeNgbh];
  
  this->EvaluateMI();
  
  //2) Compute the gradient of the intensities in the (normaly deformed) source image   (d f(g(x,mu)) / d mu) for all x where mu is the direction in which the Gateaux derivative is the highest
  Cpt_Grad_ScalarField(this->SrcImage,Gradient);
  
  //3) Compute the amplitude of the gradient for each point in the image (minus the margins)
  if (this->NZ>1){ //3D image
    for (z = this->Margin; z < this->NZ-this->Margin; z++) for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      
      //The estimated direction of the gradients is the one of mu.
      
      //3.1) direction in which we look
      LocNorm=sqrt(Gradient->G(0,x,y,z)*Gradient->G(0,x,y,z)+Gradient->G(1,x,y,z)*Gradient->G(1,x,y,z)+Gradient->G(2,x,y,z)*Gradient->G(2,x,y,z));
      
      if (LocNorm<epsilon){ // no intensity gradient
	Gradient->P(0,0,x,y,z);
	Gradient->P(0,1,x,y,z);
	Gradient->P(0,2,x,y,z);
      }
      else{// intensity gradient OK
	direcX=Gradient->G(0,x,y,z)/LocNorm;
	direcY=Gradient->G(1,x,y,z)/LocNorm;
	direcZ=Gradient->G(2,x,y,z)/LocNorm;
	
	//3.2) intensities which will be considered to compute the gradient  ... an idea: the deformations may be weighted with the kernel size
	for (i=0;i<SizeNgbh;i++){
	  tmpX=static_cast<float>(x)+(i-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
	  tmpY=static_cast<float>(y)+(i-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
	  tmpZ=static_cast<float>(z)+(i-(SizeNgbh/2))*direcZ*SizeCheckForwardBackward;
	  ListIntensities_Target[i]=this->TrgImage->G(tmpX,tmpY,tmpZ);
	  ListIntensities_Centered[i]=this->SrcImage->G(tmpX,tmpY,tmpZ);
	}
	for (i=0;i<SizeNgbh;i++){
	  if (i>0) ListIntensities_Forward[i]=ListIntensities_Centered[i-1];
	  if (i<SizeNgbh-1) ListIntensities_Backward[i]=ListIntensities_Centered[i+1];
	}
	
	tmpX=static_cast<float>(x)+(-1-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
	tmpY=static_cast<float>(y)+(-1-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
	tmpZ=static_cast<float>(z)+(-1-(SizeNgbh/2))*direcZ*SizeCheckForwardBackward;
	ListIntensities_Forward[0]=this->SrcImage->G(tmpX,tmpY,tmpZ);
	
	tmpX=static_cast<float>(x)+(SizeNgbh-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
	tmpY=static_cast<float>(y)+(SizeNgbh-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
	tmpZ=static_cast<float>(z)+(SizeNgbh-(SizeNgbh/2))*direcZ*SizeCheckForwardBackward;
	ListIntensities_Backward[SizeNgbh-1]=this->SrcImage->G(tmpX,tmpY,tmpZ);
	
	
	//3.3) min bins corresponding to the intensities of interest (to only check the information where it is not null)
	for (i=0;i<SizeNgbh;i++){
	  MinBin_Target[i]=static_cast<int>(this->GiveFloatBinT(ListIntensities_Target[i]))-1;
	  if (MinBin_Target[i]<0) MinBin_Target[i]=0;
	  if (MinBin_Target[i]>=NumberOfBinsT-3) MinBin_Target[i]=NumberOfBinsT-4;
	  
	  MinBin_Centered[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Centered[i]))-1;
	  if (MinBin_Centered[i]<0) MinBin_Centered[i]=0;
	  if (MinBin_Centered[i]>=NumberOfBinsS-3) MinBin_Centered[i]=NumberOfBinsS-4;
	  
	  MinBin_Forward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Forward[i]))-1;
	  if (MinBin_Forward[i]<0) MinBin_Forward[i]=0;
	  if (MinBin_Forward[i]>=NumberOfBinsS-3) MinBin_Forward[i]=NumberOfBinsS-4;
	  
	  MinBin_Backward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Backward[i]))-1;
	  if (MinBin_Backward[i]<0) MinBin_Backward[i]=0;
	  if (MinBin_Backward[i]>=NumberOfBinsS-3) MinBin_Backward[i]=NumberOfBinsS-4;
	}
	
	
	//3.4) evaluate the derivative of the mutual information at (x,y,z) in the direction mu
	
	//3.4.1) local mutual information (in one direction and irrespective to the intensities gradient)
	S_loc=0;
	for (k=0;k<SizeNgbh;k++){
	  tmpFl=0;
	  for (i = MinBin_Centered[k]; i < MinBin_Centered[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
	    //log_2 of the conditional probability of the bin pair [i][j]
	    tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
	    
	    //weight of the local intensity of the target image in the parzen window
	    tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
	    
	    //weight of the local intensity of the deformed source image in the parzen window
	    tmpFl*=GiveValueParzenWindowS(ListIntensities_Centered[k],i);
	    
	    //update S_loc
	    S_loc+=tmpFl;
	  }
	}
	
	//3.4.2) local forward mutual information (in one direction and irrespective to the intensities gradient)
	S_fwd=0;
	for (k=0;k<SizeNgbh;k++){
	  tmpFl=0;
	  for (i = MinBin_Forward[k]; i < MinBin_Forward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
	    //log_2 of the conditional probability of the bin pair [i][j]
	    tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
	    
	    //weight of the local intensity of the target image in the parzen window
	    tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
	    
	    //weight of the local intensity of the deformed source image in the parzen window
	    tmpFl*=GiveValueParzenWindowS(ListIntensities_Forward[k],i);
	    
	    //update S_fwd
	    S_fwd+=tmpFl;
	  }
	}
	
	//3.4.3) local backward mutual information (in one direction and irrespective to the intensities gradient)
	S_bwd=0;
	for (k=0;k<SizeNgbh;k++){
	  tmpFl=0;
	  for (i = MinBin_Backward[k]; i < MinBin_Backward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
	    //log_2 of the conditional probability of the bin pair [i][j]
	    tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
	    
	    //weight of the local intensity of the target image in the parzen window
	    tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
	    
	    //weight of the local intensity of the deformed source image in the parzen window
	    tmpFl*=GiveValueParzenWindowS(ListIntensities_Backward[k],i);
	    
	    //update S_bwd
	    S_bwd+=tmpFl;
	  }
	}
	
	
	//3.5) define the local gradient of mutual information   (the S_... values are supposed of the same sign / normaly < 0)
	if ((S_fwd-S_loc)<(S_bwd-S_loc)){ //best score in the forward direction
	  Sprime=S_fwd-S_loc;
	}
	else{
	  Sprime=S_loc-S_bwd;
	}
	
	Gradient->P(Gradient->G(0,x,y,z)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),0,x,y,z);
	Gradient->P(Gradient->G(1,x,y,z)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),1,x,y,z);
	Gradient->P(Gradient->G(2,x,y,z)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),2,x,y,z);
      }
    }
  }
  else{ //3.bis: 2D image
    for (y = this->Margin; y < this->NY-this->Margin; y++) for (x = this->Margin; x < this->NX-this->Margin; x++){
      
      //The estimated direction of the gradients is the one of mu.
      
      //3.bis.1) direction in which we look
      LocNorm=sqrt(Gradient->G(0,x,y,0)*Gradient->G(0,x,y,0)+Gradient->G(1,x,y,0)*Gradient->G(1,x,y,0));
      
      if (LocNorm<epsilon){ // no intensity gradient
	Gradient->P(0,0,x,y,0);
	Gradient->P(0,1,x,y,0);
	Gradient->P(0,2,x,y,0);
      }
      else{// intensity gradient OK
	direcX=Gradient->G(0,x,y,0)/LocNorm;
	direcY=Gradient->G(1,x,y,0)/LocNorm;
	direcZ=0;
	
	//3.bis.2) intensities which will be considered to compute the gradient  ... an idea: the deformations may be weighted with the kernel size
	for (i=0;i<SizeNgbh;i++){
	  tmpX=static_cast<float>(x)+(i-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
	  tmpY=static_cast<float>(y)+(i-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
	  tmpZ=0;
	  ListIntensities_Target[i]=this->TrgImage->G(tmpX,tmpY,tmpZ);
	  ListIntensities_Centered[i]=this->SrcImage->G(tmpX,tmpY,tmpZ);
	}
	for (i=0;i<SizeNgbh;i++){
	  if (i>0) ListIntensities_Forward[i]=ListIntensities_Centered[i-1];
	  if (i<SizeNgbh-1) ListIntensities_Backward[i]=ListIntensities_Centered[i+1];
	}
	
	tmpX=static_cast<float>(x)+(-1-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
	tmpY=static_cast<float>(y)+(-1-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
	tmpZ=0;
	ListIntensities_Forward[0]=this->SrcImage->G(tmpX,tmpY,tmpZ);
	
	tmpX=static_cast<float>(x)+(SizeNgbh-(SizeNgbh/2))*direcX*SizeCheckForwardBackward;
	tmpY=static_cast<float>(y)+(SizeNgbh-(SizeNgbh/2))*direcY*SizeCheckForwardBackward;
	tmpZ=0;
	ListIntensities_Backward[SizeNgbh-1]=this->SrcImage->G(tmpX,tmpY,tmpZ);
	
	
	//3.bis.3) min bins corresponding to the intensities of interest (to only check the information where it is not null)
	for (i=0;i<SizeNgbh;i++){
	  MinBin_Target[i]=static_cast<int>(this->GiveFloatBinT(ListIntensities_Target[i]))-1;
	  if (MinBin_Target[i]<0) MinBin_Target[i]=0;
	  if (MinBin_Target[i]>=NumberOfBinsT-3) MinBin_Target[i]=NumberOfBinsT-4;
	  
	  MinBin_Centered[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Centered[i]))-1;
	  if (MinBin_Centered[i]<0) MinBin_Centered[i]=0;
	  if (MinBin_Centered[i]>=NumberOfBinsS-3) MinBin_Centered[i]=NumberOfBinsS-4;
	  
	  MinBin_Forward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Forward[i]))-1;
	  if (MinBin_Forward[i]<0) MinBin_Forward[i]=0;
	  if (MinBin_Forward[i]>=NumberOfBinsS-3) MinBin_Forward[i]=NumberOfBinsS-4;
	  
	  MinBin_Backward[i]=static_cast<int>(this->GiveFloatBinS(ListIntensities_Backward[i]))-1;
	  if (MinBin_Backward[i]<0) MinBin_Backward[i]=0;
	  if (MinBin_Backward[i]>=NumberOfBinsS-3) MinBin_Backward[i]=NumberOfBinsS-4;
	}
	
	
	//3.bis.4) evaluate the derivative of the mutual information at (x,y,0) in the direction mu
	
	//3.bis.4.1) local mutual information (in one direction and irrespective to the intensities gradient)
	S_loc=0;
	for (k=0;k<SizeNgbh;k++){
	  tmpFl=0;
	  for (i = MinBin_Centered[k]; i < MinBin_Centered[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
	    //log_2 of the conditional probability of the bin pair [i][j]
	    tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
	    
	    //weight of the local intensity of the target image in the parzen window
	    tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
	    
	    //weight of the local intensity of the deformed source image in the parzen window
	    tmpFl*=GiveValueParzenWindowS(ListIntensities_Centered[k],i);
	    
	    //update S_loc
	    S_loc+=tmpFl;
	  }
	}
	
	//3.bis.4.2) local forward mutual information (in one direction and irrespective to the intensities gradient)
	S_fwd=0;
	for (k=0;k<SizeNgbh;k++){
	  tmpFl=0;
	  for (i = MinBin_Forward[k]; i < MinBin_Forward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
	    //log_2 of the conditional probability of the bin pair [i][j]
	    tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
	    
	    //weight of the local intensity of the target image in the parzen window
	    tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
	    
	    //weight of the local intensity of the deformed source image in the parzen window
	    tmpFl*=GiveValueParzenWindowS(ListIntensities_Forward[k],i);
	    
	    //update S_fwd
	    S_fwd+=tmpFl;
	  }
	}
	
	//3.bis.4.3) local backward mutual information (in one direction and irrespective to the intensities gradient)
	S_bwd=0;
	for (k=0;k<SizeNgbh;k++){
	  tmpFl=0;
	  for (i = MinBin_Backward[k]; i < MinBin_Backward[k]+4; i++) for (j = MinBin_Target[k]; j < MinBin_Target[k]+4; j++){
	    //log_2 of the conditional probability of the bin pair [i][j]
	    tmpFl=log(JointHistogram[i][j]/MarginalHistogramS[i])/log(2.);
	    
	    //weight of the local intensity of the target image in the parzen window
	    tmpFl*=GiveValueParzenWindowT(ListIntensities_Target[k],j);
	    
	    //weight of the local intensity of the deformed source image in the parzen window
	    tmpFl*=GiveValueParzenWindowS(ListIntensities_Backward[k],i);
	    
	    //update S_bwd
	    S_bwd+=tmpFl;
	  }
	}
	
	
	//3.bis.5) define the local gradient of mutual information   (the S_... values are supposed of the same sign / normaly < 0)
	if ((S_fwd-S_loc)<(S_bwd-S_loc)){ //best score in the forward direction
	  Sprime=S_fwd-S_loc;
	}
	else{
	  Sprime=S_loc-S_bwd;
	}
	
	Gradient->P(Gradient->G(0,x,y,0)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),0,x,y,0);
	Gradient->P(Gradient->G(1,x,y,0)*Sprime/((LocNorm)+(fabs(Sprime)/lambdaX)),1,x,y,0);
      }
    }
  }
  
  
  //4) put to zero the gradients that are not taken into account...
  //4.1) ... in the margin
  if (this->NZ>1){//... 3D
    for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      if ((z<this->Margin)||(z>=this->NZ-this->Margin)||(y<this->Margin)||(y>=this->NY-this->Margin)||(x<this->Margin)||(x>=this->NX-this->Margin)){
        Gradient->P(0,0,x,y,z);
        Gradient->P(0,1,x,y,z);
        Gradient->P(0,2,x,y,z);
      }
    }
  }
  else{//... 2D
    for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      if ((y<this->Margin)||(y>=this->NY-this->Margin)||(x<this->Margin)||(x>=this->NX-this->Margin)){
        Gradient->P(0,0,x,y,0);
        Gradient->P(0,1,x,y,0);
        Gradient->P(0,2,x,y,0);
      }
    }
  }
  //4.2) ... outside of the ROI
  if (indicatorMaskDefined==1){
    if (this->NZ>1){//... 3D
      for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (fabs(Mask->G(x,y,z)-1)>=0.0001){
        Gradient->P(0,0,x,y,z);
        Gradient->P(0,1,x,y,z);
        Gradient->P(0,2,x,y,z);
      }
    }
    else{//... 2D
      for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (fabs(Mask->G(x,y,0)-1)>=0.0001){
        Gradient->P(0,0,x,y,0);
        Gradient->P(0,1,x,y,0);
        Gradient->P(0,2,x,y,0);
      }
    }
  }
  
  //5) delete allocated variables
  delete ListIntensities_Target;
  delete ListIntensities_Centered;
  delete ListIntensities_Forward;
  delete ListIntensities_Backward;
  
  delete MinBin_Target;
  delete MinBin_Centered;
  delete MinBin_Forward;
  delete MinBin_Backward;
}






///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           6: LOW LEVEL FUNCTIONS MAKING USE OF THE CLASSES ScalarField AND VectorField 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///Isotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
void Diffusion_3D(ScalarField * SField, float alpha, float dTau,int ITERATIONS_NB, float dx,float dy, float dz){
	int x, y, z;
	int NBX,NBY,NBZ;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	int iteration;
	float Var1x,Var2x,Var3x;
	float Var1y,Var2y,Var3y;
	float Var1z,Var2z,Var3z;
  
  
  //IF alpha == 1, it can be interesting to use the following strategy (will however diffuse the image boundaries)
  //FFTconvolver3D convolver;
  //if (alpha==1){
  //  convolver.InitiateConvolver(SField->NX,SField->NY,SField->NZ,1,(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB));
  //  convolver.Convolution(SField);
  //}
  //return;
  
  
	//1) INITIALISATION
	
	//variables definition
	NBX=SField->NX;
	NBY=SField->NY;
	NBZ=SField->NZ;
	
	//precomputed values
  Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
  Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
  Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
  
	//temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
	n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
	Va= new float [n];
	Vb= new float [n];
	Vc= new float [n];
	Vd= new float [n];
	Vx= new float [n];
  
	//2) ISOTROPIC DIFFUSION
  
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		//cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		
		//2.1) diffusion - x direction
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
			for (x = 0; x < NBX; x++){
				Va[x+1]=Var2x;
				Vb[x+1]=Var1x;
				Vc[x+1]=Var2x;
        Vd[x+1]=SField->G(x,y,z);
			}
			Va[0]=Va[1]; Va[NBX+1]=Va[NBX]; //to avoid boundary effects
			Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX]; //to avoid boundary effects
			Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX]; //to avoid boundary effects
			Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
			for (x = 0; x < NBX; x++) SField->P(Vx[x+1],x,y,z);
		}
    
		//2.2) diffusion - y direction
		for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
			for (y = 0; y < NBY; y++){
				Va[y+1]=Var2y;
				Vb[y+1]=Var1y;
				Vc[y+1]=Var2y;
        Vd[y+1]=SField->G(x,y,z);
			}
			Va[0]=Va[1]; Va[NBY+1]=Va[NBY]; //to avoid boundary effects
			Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY]; //to avoid boundary effects
			Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY]; //to avoid boundary effects
			Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
			for (y = 0; y < NBY; y++) SField->P(Vx[y+1],x,y,z);
		}
    
		//2.3) diffusion - z direction
		if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
			for (z = 0; z < NBZ; z++){
				Va[z+1]=Var2z;
				Vb[z+1]=Var1z;
				Vc[z+1]=Var2z;
        Vd[z+1]=SField->G(x,y,z);
			}
			Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ]; //to avoid boundary effects
			Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ]; //to avoid boundary effects
			Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ]; //to avoid boundary effects
			Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
			for (z = 0; z < NBZ; z++) SField->P(Vx[z+1],x,y,z);
		}
		
	}
  
  
}


///Isotropic diffusion of 'SField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
void Diffusion_3D(ScalarField * SField,ScalarField * Mask,int MaskId, float alpha, float dTau,int ITERATIONS_NB,int optNoMaskNoDef, float dx,float dy, float dz){
	int x, y, z;
	int NBX,NBY,NBZ;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	int iteration;
	float Var1x,Var2x,Var3x;
	float Var1y,Var2y,Var3y;
	float Var1z,Var2z,Var3z;
  float epsilon;
  float FlMaskId;
	
	//1) INITIALISATION
	
	//variables definition
	NBX=SField->NX;
	NBY=SField->NY;
	NBZ=SField->NZ;
	epsilon=0.00001;  //intensities below this value are considered as null in the mask
  FlMaskId=static_cast<float>(MaskId);
  
	//precomputed values
  Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
  Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
  Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
  
	//temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
	n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
	Va= new float [n];
	Vb= new float [n];
	Vc= new float [n];
	Vd= new float [n];
	Vx= new float [n];
  
	//2) ISOTROPIC DIFFUSION
  
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		//cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		
		//2.1) diffusion - x direction
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
      //main values of the diagonal matrix
			for (x = 0; x < NBX; x++){
        if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
          Va[x+1]=0; 
          Vb[x+1]=1; 
          Vc[x+1]=0;
        }
        else {
          Va[x+1]=Var2x; 
          Vb[x+1]=Var1x; 
          Vc[x+1]=Var2x;
        }
        Vd[x+1]=SField->G(x,y,z);
        
			}
      //boundaries
			Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
			Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
			Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
			Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
			//matrix inversion
      TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
			//store result
			for (x = 0; x < NBX; x++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) SField->P(Vx[x+1],x,y,z);
		}
    
		//2.2) diffusion - y direction
		for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
      //main values of the diagonal matrix
			for (y = 0; y < NBY; y++){
        if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
          Va[y+1]=0; 
          Vb[y+1]=1; 
          Vc[y+1]=0;
        }
        else {
          Va[y+1]=Var2y; 
          Vb[y+1]=Var1y; 
          Vc[y+1]=Var2y;
        }
        Vd[y+1]=SField->G(x,y,z);
        
			}
      //boundaries
			Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
			Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
			Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
			Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
			//matrix inversion
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
			//store result
			for (y = 0; y < NBY; y++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) SField->P(Vx[y+1],x,y,z);
		}
    
		//2.3) diffusion - z direction
		if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
      //main values of the diagonal matrix
			for (z = 0; z < NBZ; z++){
        if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
          Va[z+1]=0; 
          Vb[z+1]=1; 
          Vc[z+1]=0;
        }
        else {
          Va[z+1]=Var2z; 
          Vb[z+1]=Var1z; 
          Vc[z+1]=Var2z;
        }
        Vd[z+1]=SField->G(x,y,z);
			}
      //boundaries
			Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
			Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
			Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
			Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
			//matrix inversion
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
			//store result
			for (z = 0; z < NBZ; z++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) SField->P(Vx[z+1],x,y,z);
		}
	}
  
  
  //3) Manage the option 'optNoMaskNoDef'
  if (optNoMaskNoDef==1)
    for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
      if(Mask->G(x,y,z)<epsilon)
        SField->P(0,x,y,z);
}







///Isotropic diffusion of 'VField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
void Diffusion_3D(VectorField * VField,float alpha, float dTau,int ITERATIONS_NB,float dx,float dy, float dz){
	int x, y, z;
	int NBX,NBY,NBZ;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	int iteration;
	float Var1x,Var2x,Var3x;
	float Var1y,Var2y,Var3y;
	float Var1z,Var2z,Var3z;
  float epsilon;
	int IdDirec;
  
  
  
  
  //IF alpha == 1, it can be interesting to use the following strategy (will however diffuse the image boundaries)
  //  FFTconvolver3D convolver;
  //  if (alpha==1){
  //  cout << "fft" << endl;
  //    convolver.InitiateConvolver(VField->NX,VField->NY,VField->NZ,1,(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB),(float)sqrt(2*dTau*ITERATIONS_NB));
  //
  //    for (IdDirec=0;IdDirec<3;IdDirec++){
  //      
  //     for (z = 0; z < VField->NZ; z++) for (y = 0; y < VField->NY; y++) for (x = 0; x < VField->NX; x++) 
  //       convolver.P(VField->G(IdDirec,x,y,z), x, y, z);
  //       
  //     convolver.Convolution();
  //     
  //      for (z = 0; z < VField->NZ; z++) for (y = 0; y < VField->NY; y++) for (x = 0; x < VField->NX; x++) 
  ////        VField->P(convolver.G(x,y,z),IdDirec, x, y,  z);
  //    }
  //  }
  // return;
  
  
	//1) INITIALISATION
	
	//variables definition
	NBX=VField->NX;
	NBY=VField->NY;
	NBZ=VField->NZ;
	epsilon=0.00001;  //intensities below this value are considered as null in the mask
  
	//precomputed values
	Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
	Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
	Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
	
	//temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
	n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
	Va= new float [n];
	Vb= new float [n];
	Vc= new float [n];
	Vd= new float [n];
	Vx= new float [n];
  
	//2) ISOTROPIC DIFFUSION
  
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		//cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		for (IdDirec=0;IdDirec<3;IdDirec++){
			//2.1) diffusion - x direction
			for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
  	    //main values of the diagonal matrix
				for (x = 0; x < NBX; x++){
	        Va[x+1]=Var2x; 
	        Vb[x+1]=Var1x; 
 	        Vc[x+1]=Var2x;
          Vd[x+1]=VField->G(IdDirec,x,y,z);
				}
        //boundaries
				Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
				Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
				Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
				Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
				//matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
				//store result
				for (x = 0; x < NBX; x++) VField->P(Vx[x+1],IdDirec,x,y,z); //only store the result in the current ROI
			}
			//2.2) diffusion - y direction
			for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
        //main values of the diagonal matrix
				for (y = 0; y < NBY; y++){
          Va[y+1]=Var2y; 
 	        Vb[y+1]=Var1y; 
 	        Vc[y+1]=Var2y;
  	      Vd[y+1]=VField->G(IdDirec,x,y,z);
				}
  	    //boundaries
				Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
				Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
				Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
				Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
				//matrix inversion
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
				//store result
				for (y = 0; y < NBY; y++)  VField->P(Vx[y+1],IdDirec,x,y,z);
			}
	    
			//2.3) diffusion - z direction
			if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
	      //main values of the diagonal matrix
				for (z = 0; z < NBZ; z++){
 	        Va[z+1]=Var2z; 
          Vb[z+1]=Var1z; 
          Vc[z+1]=Var2z;
  	      Vd[z+1]=VField->G(IdDirec,x,y,z);
				}
  	    //boundaries
				Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
				Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
				Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
				Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
				//matrix inversion
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
	  		//store result
	  		for (z = 0; z < NBZ; z++) VField->P(Vx[z+1],IdDirec,x,y,z);
  		}
		}
	}
  
}





///Isotropic diffusion of 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
///If direction > 0 then the smoothing is only performed in the direction {1=X, 2=Y or 3=Z}
void Diffusion_3D(VectorField * VField,ScalarField * Mask,int MaskId, float alpha, float dTau,int ITERATIONS_NB,  int optNoMaskNoDef,int direction, float dx,float dy, float dz){
	int x, y, z;
	int NBX,NBY,NBZ;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	int iteration;
	float Var1x,Var2x,Var3x;
	float Var1y,Var2y,Var3y;
	float Var1z,Var2z,Var3z;
  float epsilon;
  float FlMaskId;
	int IdDirec;
  
  
  

  
	//1) INITIALISATION
	
	//variables definition
	NBX=VField->NX;
	NBY=VField->NY;
	NBZ=VField->NZ;
	epsilon=0.00001;  //intensities below this value are considered as null in the mask
  FlMaskId=static_cast<float>(MaskId);
  
	//precomputed values
	Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
	Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
	Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
	
	//temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
	n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
	Va= new float [n];
	Vb= new float [n];
	Vc= new float [n];
	Vd= new float [n];
	Vx= new float [n];
  
	//2) ISOTROPIC DIFFUSION
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		//cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		for (IdDirec=0;IdDirec<3;IdDirec++) if ((direction<0)||(direction==IdDirec)) {
			//2.1) diffusion - x direction
			for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) {
  	    //main values of the diagonal matrix
				for (x = 0; x < NBX; x++){
  	      if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
  	        Va[x+1]=0; 
  	        Vb[x+1]=1; 
  	        Vc[x+1]=0;
            Vd[x+1]=VField->G(0,x,y,z);
  	      }
  	      else {
  	        Va[x+1]=Var2x; 
  	        Vb[x+1]=Var1x; 
  	        Vc[x+1]=Var2x;
            Vd[x+1]=VField->G(IdDirec,x,y,z);
          }
          //region boundaries
          if ((fabs(Mask->G(x,y,z)-FlMaskId)>epsilon)&&(x<NBX-3)&&(x>2)){
            if      (fabs(Mask->G(x-1,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x-1,y,z);
            else if (fabs(Mask->G(x+1,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x+1,y,z);
            else if (fabs(Mask->G(x-2,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x-2,y,z);
            else if (fabs(Mask->G(x+2,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x+2,y,z);
            else if (fabs(Mask->G(x-3,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x-3,y,z);
            else if (fabs(Mask->G(x+3,y,z)-FlMaskId)<epsilon) Vd[x+1]=VField->G(IdDirec,x+3,y,z);
          }
        }
        //boundaries
				Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
				Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
				Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
				Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
				//matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
				//store result
				for (x = 0; x < NBX; x++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[x+1],IdDirec,x,y,z); //only store the result in the current ROI
			}
			//2.2) diffusion - y direction
			for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
        //main values of the diagonal matrix
				for (y = 0; y < NBY; y++){
  	      if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
  	        Va[y+1]=0; 
            Vb[y+1]=1; 
  	        Vc[y+1]=0;
            Vd[y+1]=VField->G(0,x,y,z);
  	      }
  	      else {
            Va[y+1]=Var2y; 
  	        Vb[y+1]=Var1y; 
  	        Vc[y+1]=Var2y;
            Vd[y+1]=VField->G(IdDirec,x,y,z);
  	      }
          //region boundaries
          if ((fabs(Mask->G(x,y,z)-FlMaskId)>epsilon)&&(y<NBY-3)&&(y>2)){
            if      (fabs(Mask->G(x,y-1,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y-1,z);
            else if (fabs(Mask->G(x,y+1,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y+1,z);
            else if (fabs(Mask->G(x,y-2,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y-2,z);
            else if (fabs(Mask->G(x,y+2,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y+2,z);
            else if (fabs(Mask->G(x,y-3,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y-3,z);
            else if (fabs(Mask->G(x,y+3,z)-FlMaskId)<epsilon) Vd[y+1]=VField->G(IdDirec,x,y+3,z);
          }
				}
  	    //boundaries
				Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
				Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
				Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
				Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
				//matrix inversion
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
				//store result
				for (y = 0; y < NBY; y++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon)  VField->P(Vx[y+1],IdDirec,x,y,z);
			}
	    
			//2.3) diffusion - z direction
			if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
	      //main values of the diagonal matrix
				for (z = 0; z < NBZ; z++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[z+1]=0; 
            Vb[z+1]=1; 
  	        Vc[z+1]=0;
            Vd[z+1]=VField->G(0,x,y,z);
  	      }
  	      else {
  	        Va[z+1]=Var2z; 
            Vb[z+1]=Var1z; 
            Vc[z+1]=Var2z;
            Vd[z+1]=VField->G(IdDirec,x,y,z);
  	      }
          //region boundaries
          if ((fabs(Mask->G(x,y,z)-FlMaskId)>epsilon)&&(z<NBZ-3)&&(z>2)){
            if      (fabs(Mask->G(x,y,z-1)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z-1);
            else if (fabs(Mask->G(x,y,z+1)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z+1);
            else if (fabs(Mask->G(x,y,z-2)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z-2);
            else if (fabs(Mask->G(x,y,z+2)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z+2);
            else if (fabs(Mask->G(x,y,z-3)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z-3);
            else if (fabs(Mask->G(x,y,z+3)-FlMaskId)<epsilon) Vd[z+1]=VField->G(IdDirec,x,y,z+3);
          }
				}

  	    //boundaries
				Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
				Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
				Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
				Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
				//matrix inversion
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
	  		//store result
	  		for (z = 0; z < NBZ; z++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon)  VField->P(Vx[z+1],IdDirec,x,y,z);
  		}
		}
	}

  //3) Manage the option 'optNoMaskNoDef'
  if (optNoMaskNoDef==1)
    for (IdDirec=0;IdDirec<3;IdDirec++) for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
      if(Mask->G(x,y,z)<epsilon)
        VField->P(0,IdDirec,x,y,z);
}









///Isotropic diffusion of 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///!!! -> Dirichlet conditions are considered at the domain boundaries !!!
void Diffusion_3Dbis(VectorField * VField,ScalarField * Mask,int MaskId,float alpha, float dTau,int ITERATIONS_NB, int optNoMaskNoDef, float dx,float dy, float dz){
	int x, y, z;
	int NBX,NBY,NBZ;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	int iteration;
	float Var1x,Var2x,Var3x;
	float Var1y,Var2y,Var3y;
	float Var1z,Var2z,Var3z;
  float epsilon;
  float FlMaskId;
	int IdDirec;
  
  
  
  
	//1) INITIALISATION
	
	//variables definition
	NBX=VField->NX;
	NBY=VField->NY;
	NBZ=VField->NZ;
	epsilon=0.00001;  //intensities below this value are considered as null in the mask
  FlMaskId=static_cast<float>(MaskId);
  
	//precomputed values
	Var1x=(1/dTau)+2*(alpha/(dx*dx));  Var2x=-alpha/(dx*dx);  Var3x=1/dTau;  Var1x=Var1x/Var3x;  Var2x=Var2x/Var3x;
	Var1y=(1/dTau)+2*(alpha/(dy*dy));  Var2y=-alpha/(dy*dy);  Var3y=1/dTau;  Var1y=Var1y/Var3y;  Var2y=Var2y/Var3y;
	Var1z=(1/dTau)+2*(alpha/(dz*dz));  Var2z=-alpha/(dz*dz);  Var3z=1/dTau;  Var1z=Var1z/Var3z;  Var2z=Var2z/Var3z;
	
	//temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
	n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
	Va= new float [n];
	Vb= new float [n];
	Vc= new float [n];
	Vd= new float [n];
	Vx= new float [n];
  
	//2) ISOTROPIC DIFFUSION
  
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		//cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		for (IdDirec=0;IdDirec<3;IdDirec++){
			//2.1) diffusion - x direction
			for (z = 2; z < NBZ; z++) for (y = 0; y < NBY; y++) {
  	    //main values of the diagonal matrix
				for (x = 0; x < NBX; x++){
  	      if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
  	        Va[x+1]=0; 
  	        Vb[x+1]=1; 
  	        Vc[x+1]=0;
            if (optNoMaskNoDef==1) Vd[x+1]=0;
            else Vd[x+1]=VField->G(IdDirec,x,y,z);
  	      }
  	      else {
  	        Va[x+1]=Var2x; 
  	        Vb[x+1]=Var1x; 
  	        Vc[x+1]=Var2x;
            Vd[x+1]=VField->G(IdDirec,x,y,z);
          }
          
				}
        //boundaries
				Va[0]=Va[1]; Va[NBX+1]=Va[NBX];
				Vb[0]=Vb[1]; Vb[NBX+1]=Vb[NBX];
				Vc[0]=Vc[1]; Vc[NBX+1]=Vc[NBX];
				Vd[0]=Vd[1]; Vd[NBX+1]=Vd[NBX];
				//matrix inversion
        TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
				//store result
				for (x = 0; x < NBX; x++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[x+1],IdDirec,x,y,z);
			}
			//2.2) diffusion - y direction
			for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++){
        //main values of the diagonal matrix
				for (y = 0; y < NBY; y++){
  	      if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
  	        Va[y+1]=0; 
            Vb[y+1]=1; 
  	        Vc[y+1]=0;
            if (optNoMaskNoDef==1) Vd[y+1]=0;
            else Vd[y+1]=VField->G(IdDirec,x,y,z);
  	      }
  	      else {
            Va[y+1]=Var2y; 
  	        Vb[y+1]=Var1y; 
  	        Vc[y+1]=Var2y;
            Vd[y+1]=VField->G(IdDirec,x,y,z);
  	      }
				}
  	    //boundaries
				Va[0]=Va[1]; Va[NBY+1]=Va[NBY];
				Vb[0]=Vb[1]; Vb[NBY+1]=Vb[NBY];
				Vc[0]=Vc[1]; Vc[NBY+1]=Vc[NBY];
				Vd[0]=Vd[1]; Vd[NBY+1]=Vd[NBY];
				//matrix inversion
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
				//store result
        for (y = 0; y < NBY; y++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[y+1],IdDirec,x,y,z);
			}
	    
			//2.3) diffusion - z direction
      if (NBZ>1) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) {
	      //main values of the diagonal matrix
				for (z = 0; z < NBZ; z++){
          if(fabs(Mask->G(x,y,z)-FlMaskId)>epsilon){
            Va[z+1]=0; 
            Vb[z+1]=1;
  	        Vc[z+1]=0;
            if (optNoMaskNoDef==1) Vd[z+1]=0;
            else Vd[z+1]=VField->G(IdDirec,x,y,z);
  	      }
  	      else {
  	        Va[z+1]=Var2z; 
            Vb[z+1]=Var1z; 
            Vc[z+1]=Var2z;
            Vd[z+1]=VField->G(IdDirec,x,y,z);
  	      }
				}
  	    //boundaries
				Va[0]=Va[1]; Va[NBZ+1]=Va[NBZ];
				Vb[0]=Vb[1]; Vb[NBZ+1]=Vb[NBZ];
				Vc[0]=Vc[1]; Vc[NBZ+1]=Vc[NBZ];
				Vd[0]=Vd[1]; Vd[NBZ+1]=Vd[NBZ];
				//matrix inversion
				TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
	  		//store result
        for (z = 0; z < NBZ; z++) if(fabs(Mask->G(x,y,z)-FlMaskId)<epsilon) VField->P(Vx[z+1],IdDirec,x,y,z);
  		}
		}
	}
  
  //3) Manage the option 'optNoMaskNoDef'
  if (optNoMaskNoDef==1)
    for (IdDirec=0;IdDirec<3;IdDirec++) for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) 
      if(Mask->G(x,y,z)<epsilon)
        VField->P(0,IdDirec,x,y,z);
  
}








///Anisotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Directions and intensities of diffusion are defined using ax, ay, az.
///The size of each voxel is defined by dx, dy, dz.
///Remark: very close to the classic paper of Perona-Malik
void anisoDiff_3D(ScalarField * SField,float ax, float ay, float az, float dTau,int ITERATIONS_NB, float dx,float dy, float dz){
	int x, y, z;
	ScalarField imageE;
	ScalarField imageO;
	int NBX,NBY,NBZ;
	double dIdx,dIdy,dIdz;
	float DivDgradI;
	float *Va;
	float *Vb; 
	float *Vc;
	float *Vd;
	float *Vx;
	int n;
	float Dxx_div_dxSq,Dyy_div_dySq,Dzz_div_dzSq;
	int iteration;
	float DivPowDxSqu,DivPowDySqu,DivPowDzSqu;
	
	//1) INITIALISATION
	
	//variables definition
	NBX=SField->NX+2;  //for boundary effects
	NBY=SField->NY+2;  //for boundary effects
	NBZ=SField->NZ+2;  //for boundary effects
	
	//temporary input and output images
  imageE.CreateVoidField(NBX,NBY,NBZ);
  imageO.CreateVoidField(NBX,NBY,NBZ);
	
	//precomputed values
	DivPowDxSqu=1./pow(dx,2);
	DivPowDySqu=1./pow(dy,2);
	DivPowDzSqu=1./pow(dz,2);
  
	
	//temporary variables dedicated to the semi implicit scheme (+4 is to avoid boundary effects)
	n=max(max(NBX,NBY),NBZ)+4; //for boundary effects
	Va= new float [n];
	Vb= new float [n];
	Vc= new float [n];
	Vd= new float [n];
	Vx= new float [n];
  
	//2) ANISOTROPIC DIFFUSION
  
  //2.1) copy the values of the input image in a temporary 3D image
	for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
		imageE.P(SField->G(x, y, z),x+1,y+1,z+1);
	
	for (z = 0; z < NBZ-2; z++)  for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
		imageO.P(SField->G(x, y, z),x+1,y+1,z+1);
  
	
	//image extension to avoid boundary effects
	for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE.P(imageE.G(1,y,z),0,y,z);
	for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageE.P(imageE.G(NBX-2,y,z),NBX-1,y,z);
	for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO.P(imageO.G(1,y,z),0,y,z);
	for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) imageO.P(imageO.G(NBX-2,y,z),NBX-1,y,z);
	
	for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,1,z),x,0,z);
	for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,NBY-2,z),x,NBY-1,z);
	for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,1,z),x,0,z);
	for (z = 0; z < NBZ; z++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,NBY-2,z),x,NBY-1,z);
	
	for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,y,1),x,y,0);
	for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageE.P(imageE.G(x,y,NBZ-2),x,y,NBZ-1);
	for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,y,1),x,y,0);
	for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++) imageO.P(imageO.G(x,y,NBZ-2),x,y,NBZ-1);
	
	//2.2) diffusion in the temporary 3D image - ADI semi implicit scheme
	for (iteration=0 ; iteration<ITERATIONS_NB; iteration++){
		cout << "Iteration " << iteration+1 << " / " << ITERATIONS_NB << "\n";
		
		//2.2.2) diffusion - x implicit / y,z explicit
		//2.2.2.1) explicit part
		for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdy=(imageE.G(x,y+1,z)-imageE.G(x,y-1,z))/(2*dy);
			Dyy_div_dySq=(1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu;
			dIdz=(imageE.G(x,y,z+1)-imageE.G(x,y,z-1))/(2*dz);
			Dzz_div_dzSq=(1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu;
			//new value of the voxel
			DivDgradI=(imageE.G(x,y+1,z)-2*imageE.G(x,y,z)+imageE.G(x,y-1,z))*Dyy_div_dySq+
      (imageE.G(x,y,z+1)-2*imageE.G(x,y,z)+imageE.G(x,y,z-1))*Dzz_div_dzSq;
			
			imageO.P(imageE.G(x,y,z)+(dTau/3.)*DivDgradI,x,y,z);
		}
    
		//2.2.2.2) implicit part
		for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) {
			for (x = 1; x < NBX-1; x++){
				dIdx=(imageE.G(x+1,y,z)-imageE.G(x-1,y,z))/(2*dx);
				Dxx_div_dxSq=(1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu;
				Va[x+1]=(dTau/3.)*Dxx_div_dxSq;
				Vb[x+1]=-1-2*(dTau/3.)*Dxx_div_dxSq;
				Vc[x+1]=(dTau/3.)*Dxx_div_dxSq;
        Vd[x+1]=imageE.G(x,y,z); //why not imageO ???
			}
			Va[1]=Va[3]; Va[0]=Va[4]; Va[NBX]=Va[NBX-2]; Va[NBX+1]=Va[NBX-3]; //to avoid boundary effects
			Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBX]=Vb[NBX-2]; Vb[NBX+1]=Vb[NBX-3]; //to avoid boundary effects
			Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBX]=Vc[NBX-2]; Vc[NBX+1]=Vc[NBX-3]; //to avoid boundary effects
			Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBX]=Vd[NBX-2]; Vd[NBX+1]=Vd[NBX-3]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBX+2);
			for (x = 1; x < NBX-1; x++) imageO.P(-Vx[x+1],x,y,z);
		}
		
		//2.2.3) diffusion - y implicit / x,z explicit
		//2.2.3.1) explicit part
		for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdx=(imageO.G(x+1,y,z)-imageO.G(x-1,y,z))/(2*dx);
			Dxx_div_dxSq=(1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu;
			dIdz=(imageO.G(x,y,z+1)-imageO.G(x,y,z-1))/(2*dz);
			Dzz_div_dzSq=(1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu;
			//new value of the voxel
			DivDgradI=(imageO.G(x+1,y,z)-2*imageO.G(x,y,z)+imageO.G(x-1,y,z))*Dxx_div_dxSq+
      (imageO.G(x,y,z+1)-2*imageO.G(x,y,z)+imageO.G(x,y,z-1))*Dzz_div_dzSq;
			
			imageE.P(imageO.G(x,y,z)+(dTau/3.)*DivDgradI,x,y,z);
		}
		
		//2.2.3.2) implicit part
		for (z = 1; z < NBZ-1; z++) for (x = 1; x < NBX-1; x++){
			for (y = 1; y < NBY-1; y++){
				dIdy=(imageO.G(x,y+1,z)-imageO.G(x,y-1,z))/(2*dy);
				Dyy_div_dySq=(1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu;
				Va[y+1]=(dTau/3.)*Dyy_div_dySq;
				Vb[y+1]=-1-2*(dTau/3.)*Dyy_div_dySq;
				Vc[y+1]=(dTau/3.)*Dyy_div_dySq;
				Vd[y+1]=imageO.G(x,y,z);
			}
			Va[1]=Va[3]; Va[0]=Va[4]; Va[NBY]=Va[NBY-2]; Va[NBY+1]=Va[NBY-3]; //to avoid boundary effects
			Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBY]=Vb[NBY-2]; Vb[NBY+1]=Vb[NBY-3]; //to avoid boundary effects
			Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBY]=Vc[NBY-2]; Vc[NBY+1]=Vc[NBY-3]; //to avoid boundary effects
			Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBY]=Vd[NBY-2]; Vd[NBY+1]=Vd[NBY-3]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBY+2);
			for (y = 1; y < NBY-1; y++) imageE.P(-Vx[y+1],x,y,z);
		}
    
		//2.2.4) diffusion - z implicit / x,y explicit
		//2.2.4.1) explicit part
		for (z = 1; z < NBZ-1; z++) for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			dIdx=(imageE.G(x+1,y,z)-imageE.G(x-1,y,z))/(2*dx);
			Dxx_div_dxSq=(1-exp(-3.314/pow((dIdx/ax),4)))*DivPowDxSqu;
			dIdy=(imageE.G(x,y+1,z)-imageE.G(x,y-1,z))/(2*dy);
			Dyy_div_dySq=(1-exp(-3.314/pow((dIdy/ay),4)))*DivPowDySqu;
			//new value of the voxel
			DivDgradI=(imageE.G(x+1,y,z)-2*imageE.G(x,y,z)+imageE.G(x-1,y,z))*Dxx_div_dxSq+
      (imageE.G(x,y+1,z)-2*imageE.G(x,y,z)+imageE.G(x,y-1,z))*Dyy_div_dySq;
      
			imageO.P(imageE.G(x,y,z)+(dTau/3.)*DivDgradI,x,y,z);
		}
		
		//2.2.4.2) implicit part
		for (y = 1; y < NBY-1; y++) for (x = 1; x < NBX-1; x++) {
			for (z = 1; z < NBZ-1; z++){
				dIdz=(imageE.G(x,y,z+1)-imageE.G(x,y,z-1))/(2*dz);
				Dzz_div_dzSq=(1-exp(-3.314/pow((dIdz/az),4)))*DivPowDzSqu;
				Va[z+1]=(dTau/3.)*Dzz_div_dzSq;
				Vb[z+1]=-1-2*(dTau/3.)*Dzz_div_dzSq;
				Vc[z+1]=(dTau/3.)*Dzz_div_dzSq;
				Vd[z+1]=imageE.G(x,y,z);
			}
			Va[1]=Va[3]; Va[0]=Va[4]; Va[NBZ]=Va[NBZ-2]; Va[NBZ+1]=Va[NBZ-3]; //to avoid boundary effects
			Vb[1]=Vb[3]; Vb[0]=Vb[4]; Vb[NBZ]=Vb[NBZ-2]; Vb[NBZ+1]=Vb[NBZ-3]; //to avoid boundary effects
			Vc[1]=Vc[3]; Vc[0]=Vc[4]; Vc[NBZ]=Vc[NBZ-2]; Vc[NBZ+1]=Vc[NBZ-3]; //to avoid boundary effects
			Vd[1]=Vd[3]; Vd[0]=Vd[4]; Vd[NBZ]=Vd[NBZ-2]; Vd[NBZ+1]=Vd[NBZ-3]; //to avoid boundary effects
			TridiagonalSolveFloat(Va,Vb,Vc,Vd,Vx,NBZ+2);
			for (z = 1; z < NBZ-1; z++) imageO.P(-Vx[z+1],x,y,z);
		}
		
		//2.2.5) temporary output image is reinjected in temporary input image
		for (z = 0; z < NBZ; z++) for (y = 0; y < NBY; y++) for (x = 0; x < NBX; x++)
			imageE.P(imageO.G(x,y,z),x,y,z);
    
	}
	//2.3) save the filtered temporary 3D image in VoxelType in the  output image at time t
	for (z = 0; z < NBZ-2; z++) for (y = 0; y < NBY-2; y++) for (x = 0; x < NBX-2; x++)
    SField->P(imageE.G(x+1,y+1,z+1),x,y,z);
  
}



//compute the distance map of the edges in a 3D Mask. Greedy algorithm updating the distances in a band.
//-> 'NeaBoun' contains a vector field pointing to the nearest boundaries
//-> 'Distance' contains the distance to the nearest boundary
void Cpt_NearestBoundary(ScalarField * Mask,VectorField * NeaBoun,ScalarField * Distance){
  int x,y,z,direc,i,j,k;
  float epsilon;
  int *ListBx;
  int *ListBy;
  int *ListBz;
  int *ListBx2;
  int *ListBy2;
  int *ListBz2;
  int NbB;
  int NbB2;
  float TmpDist;
  int ix,iy,iz;
  int changes;
  int NbMaxBoundaryPts;
  float mindxdydz;
  float MaxDist;
  float dx,dy,dz;
  float ThreshDistanceNotAnyMoreTreated;

  epsilon=0.00001;
  NbMaxBoundaryPts=30000000;
  
  ListBx = new int [NbMaxBoundaryPts];
  ListBy = new int [NbMaxBoundaryPts];
  ListBz = new int [NbMaxBoundaryPts];
  
  ListBx2 = new int [NbMaxBoundaryPts];
  ListBy2 = new int [NbMaxBoundaryPts];
  ListBz2 = new int [NbMaxBoundaryPts];
  
  //initate dx,dy,dz, MaxDist, NeaBoun and Distance
  dx=sqrt(Mask->Image2World[0][0]*Mask->Image2World[0][0]+Mask->Image2World[0][1]*Mask->Image2World[0][1]+Mask->Image2World[0][2]*Mask->Image2World[0][2]);
  dy=sqrt(Mask->Image2World[1][0]*Mask->Image2World[1][0]+Mask->Image2World[1][1]*Mask->Image2World[1][1]+Mask->Image2World[1][2]*Mask->Image2World[1][2]);
  dz=sqrt(Mask->Image2World[2][0]*Mask->Image2World[2][0]+Mask->Image2World[2][1]*Mask->Image2World[2][1]+Mask->Image2World[2][2]*Mask->Image2World[2][2]);
  
  MaxDist=sqrt((Mask->NZ*Mask->NZ*dz*dz)+(Mask->NY*Mask->NY*dy*dy)+(Mask->NX*Mask->NX*dx*dx));
  
  mindxdydz=dx;
  if (dy<mindxdydz) mindxdydz=dy;
  if (dz<mindxdydz) mindxdydz=dz;
  
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    Distance->P(MaxDist*MaxDist,x,y,z);
  
  for (direc = 0; direc < 3; direc++) for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++) 
    NeaBoun->P(0,direc,x,y,z);

  
  //1) FIRST BOUNDARIES
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX-1; x++){
    if(fabs(Mask->G(x,y,z)-Mask->G(x+1,y,z))>epsilon){
      Distance->P(1,x,y,z);      
      NeaBoun->Add(0.5*dx,0,x,y,z);   
      Distance->P(1,x+1,y,z);    
      NeaBoun->Add(-0.5*dx,0,x+1,y,z);
    }
  }
    
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY-1; y++)  for (x = 0; x < Mask->NX; x++){
    if(fabs(Mask->G(x,y,z)-Mask->G(x,y+1,z))>epsilon){
      Distance->P(1,x,y,z);       
      NeaBoun->Add(0.5*dy,1,x,y,z);   
      Distance->P(1,x,y+1,z);     
      NeaBoun->Add(-0.5*dy,1,x,y+1,z);
    }
  }
    
  for (z = 0; z < Mask->NZ-1; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
    if(fabs(Mask->G(x,y,z)-Mask->G(x,y,z+1))>epsilon){
      Distance->P(1,x,y,z);     
      NeaBoun->Add(0.5*dz,2,x,y,z);   
      Distance->P(1,x,y,z+1);   
      NeaBoun->Add(-0.5*dz,2,x,y,z+1);
    }
  }
  
  NbB=0;
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++){
    if (Distance->G(x,y,z)<2){
      TmpDist=sqrt(NeaBoun->G(0,x,y,z)*NeaBoun->G(0,x,y,z)+NeaBoun->G(1,x,y,z)*NeaBoun->G(1,x,y,z)+NeaBoun->G(2,x,y,z)*NeaBoun->G(2,x,y,z));
      Distance->P(TmpDist,x,y,z);
      ListBx[NbB]=x; 
      ListBy[NbB]=y; 
      ListBz[NbB]=z;
      NbB++;
    }
    else{
      NeaBoun->P(MaxDist*MaxDist,0,x,y,z);
      NeaBoun->P(MaxDist*MaxDist,1,x,y,z);
      NeaBoun->P(MaxDist*MaxDist,2,x,y,z);
    }
  }
  
  //2) DISTANCE PROPAGATION
  changes=1;
  
  ThreshDistanceNotAnyMoreTreated=-3*mindxdydz;
  while (changes==1){
    changes=0;
    NbB2=0;
    
    
    //one iteration of the propagation
    for (i=0; i<NbB; i++) {
      //init
      x=ListBx[i]; y=ListBy[i]; z=ListBz[i];
      
      //update the list with the current point for the next iteration
      if (Distance->G(x,y,z)>ThreshDistanceNotAnyMoreTreated){
	ListBx2[NbB2]=x;
	ListBy2[NbB2]=y;
	ListBz2[NbB2]=z;
	NbB2++;
      }
      
      //update the 6 ngbh of the points in the list
      ix=0;iy=0;iz=0;
      for (j=0;j<6;j++){
	if (j==0) {ix=1;iy=0;iz=0;}
	if (j==1) {ix=-1;}
	if (j==2) {ix=0;iy=1;}
	if (j==3) {iy=-1;}
	if (j==4) {iy=0;iz=1;}
	if (j==5) {iz=-1;};
	
	TmpDist=sqrt((NeaBoun->G(0,x,y,z)-ix*dx)*(NeaBoun->G(0,x,y,z)-ix*dx)+(NeaBoun->G(1,x,y,z)-iy*dy)*(NeaBoun->G(1,x,y,z)-iy*dy)+(NeaBoun->G(2,x,y,z)-iz*dz)*(NeaBoun->G(2,x,y,z)-iz*dz));
	
	if ((x+ix>=0)&&(x+ix<NeaBoun->NX)&&(y+iy>=0)&&(y+iy<NeaBoun->NY)&&(z+iz>=0)&&(z+iz<NeaBoun->NZ)) if ((TmpDist<Distance->G(x+ix,y+iy,z+iz))&&(TmpDist<MaxDist)){
	  if (Distance->G(x+ix,y+iy,z+iz)>MaxDist*MaxDist-2){
	    ListBx2[NbB2]=x+ix;
	    ListBy2[NbB2]=y+iy;
	    ListBz2[NbB2]=z+iz;
	    NbB2++;
	    if (NbB2>NbMaxBoundaryPts-10){
	      cout << "Too much points at the boundary -> limit extended" << endl;
	      delete ListBx;
	      delete ListBy;
	      delete ListBz;
	      ListBx = new int [Mask->NX*Mask->NY*Mask->NZ];
	      ListBy = new int [Mask->NX*Mask->NY*Mask->NZ];
	      ListBz = new int [Mask->NX*Mask->NY*Mask->NZ];
	      
	      for (k=0;k<NbB2;k++) ListBx[k]=ListBx2[k];
	      for (k=0;k<NbB2;k++) ListBy[k]=ListBy2[k];
	      for (k=0;k<NbB2;k++) ListBz[k]=ListBz2[k];
	      
	      delete ListBx2;
	      delete ListBy2;
	      delete ListBz2;
	      ListBx2 = new int [Mask->NX*Mask->NY*Mask->NZ];
	      ListBy2 = new int [Mask->NX*Mask->NY*Mask->NZ];
	      ListBz2 = new int [Mask->NX*Mask->NY*Mask->NZ];
	      
	      for (k=0;k<NbB2;k++) ListBx2[k]=ListBx[k];
	      for (k=0;k<NbB2;k++) ListBy2[k]=ListBy[k];
	      for (k=0;k<NbB2;k++) ListBz2[k]=ListBz[k];
	      
	      NbMaxBoundaryPts=Mask->NX*Mask->NY*Mask->NZ;
	    } 
	    changes=1;
	  }
	  NeaBoun->P(NeaBoun->G(0,x,y,z)-ix*dx,0,x+ix,y+iy,z+iz);
	  NeaBoun->P(NeaBoun->G(1,x,y,z)-iy*dy,1,x+ix,y+iy,z+iz);
	  NeaBoun->P(NeaBoun->G(2,x,y,z)-iz*dz,2,x+ix,y+iy,z+iz);
	  Distance->P(TmpDist,x+ix,y+iy,z+iz);
	}
	
      }
      
    }
    
    //prepare the next iteration of the propagation
    for(i=0;i<NbB2;i++){
      ListBx[i]=ListBx2[i];
      ListBy[i]=ListBy2[i];
      ListBz[i]=ListBz2[i];
    }
    NbB=NbB2;
    
    ThreshDistanceNotAnyMoreTreated+=mindxdydz;
  }
  
  
  //put 1 in the points of Distance in which the computations have been done and 0 elsewhere
  for (z = 0; z < Mask->NZ; z++) for (y = 0; y < Mask->NY; y++)  for (x = 0; x < Mask->NX; x++)if (Distance->G(x,y,z)>MaxDist*MaxDist-2){
      Distance->P(0,x,y,z);
      NeaBoun->P(0,0,x,y,z);
      NeaBoun->P(0,1,x,y,z);
      NeaBoun->P(0,2,x,y,z);
  }

  delete ListBx;
  delete ListBy;
  delete ListBz;
  delete ListBx2;
  delete ListBy2;
  delete ListBz2;
  
}



///Remove the normal contributions of a velocity field that are too close to a boundary
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which SmoothedField has reduced normal contirbutions
void RemoveNormalContributions(VectorField * SmoothedField,VectorField * NeaBoun,ScalarField * TempSF,float BoundaMargin, float dx,float dy, float dz,int SetBoundaryToZero){
  int x,y,z;
  int i,j;
  float epsilon;
  float normLoc;
  float normLoc3,ScalProd;
  float DeltaTime;
  float MaxDist_DistMap;
  int direction;
  int x2,y2,z2;
  float NCx,NCy,NCz;
  float NBlocX,NBlocY,NBlocZ;
  int dxyz;
  
  dxyz=dx;
  if (dy<dxyz) dxyz=dy;
  if (dz<dxyz) dxyz=dz;
  
  //cout << "Remove normal contributions" << endl;
  
  if (SmoothedField->NZ>1){ //3D image
    for (z = 2; z < SmoothedField->NZ-2; z++) for (y = 2; y < SmoothedField->NY-2; y++)  for (x = 2; x < SmoothedField->NX-2; x++) if (TempSF->G(x,y,z)>0.5){
      
      
      //1) INIT
      //put the vector representing the nearest boundary in NCx,NCy,NCz
      NCx=NeaBoun->G(0,x,y,z);
      NCy=NeaBoun->G(1,x,y,z);
      NCz=NeaBoun->G(2,x,y,z);
      
      normLoc=sqrt(NCx*NCx+NCy*NCy+NCz*NCz); //distance to the boundary in millimiters
      
      //2) REMOVE THE NORMAL CONTRIBUTION
      if ((normLoc<=0.5*dxyz)&&(SetBoundaryToZero==1)){
	NCx=0; NCy=0; NCz=0; normLoc=0;
	if (TempSF->G(x+1,y,z)>0.5) NCx=NeaBoun->G(0,x+1,y,z); NCy=NeaBoun->G(1,x+1,y,z); NCz=NeaBoun->G(2,x+1,y,z); normLoc++;
	if (TempSF->G(x-1,y,z)>0.5) NCx=NeaBoun->G(0,x-1,y,z); NCy=NeaBoun->G(1,x-1,y,z); NCz=NeaBoun->G(2,x-1,y,z); normLoc++;
	if (TempSF->G(x,y+1,z)>0.5) NCx=NeaBoun->G(0,x,y+1,z); NCy=NeaBoun->G(1,x,y+1,z); NCz=NeaBoun->G(2,x,y+1,z); normLoc++;
	if (TempSF->G(x,y-1,z)>0.5) NCx=NeaBoun->G(0,x,y-1,z); NCy=NeaBoun->G(1,x,y-1,z); NCz=NeaBoun->G(2,x,y-1,z); normLoc++;
	if (TempSF->G(x,y,z+1)>0.5) NCx=NeaBoun->G(0,x,y,z+1); NCy=NeaBoun->G(1,x,y,z+1); NCz=NeaBoun->G(2,x,y,z+1); normLoc++;
	if (TempSF->G(x,y,z-1)>0.5) NCx=NeaBoun->G(0,x,y,z-1); NCy=NeaBoun->G(1,x,y,z-1); NCz=NeaBoun->G(2,x,y,z-1); normLoc++;
	
	if (normLoc>0.5){
	  NCx/=normLoc;
	  NCy/=normLoc;
	  NCz/=normLoc;
	  normLoc=sqrt(NCx*NCx+NCy*NCy+NCz*NCz);
	}
      }
      
      if ((normLoc<BoundaMargin)&&(normLoc>0.5*dxyz)){
        //normalise NeaBoun   (in mm)
        NCx/=normLoc;
        NCy/=normLoc;
        NCz/=normLoc;
        
        
        //compute the scalar product with SmoothedField    (in mm)
        ScalProd=SmoothedField->G(0,x,y,z)*dx*NCx+SmoothedField->G(1,x,y,z)*dy*NCy+SmoothedField->G(2,x,y,z)*dz*NCz;
        
        //add a weight to ScalProd to take into account the distance to the bounadry
        ScalProd*=((normLoc/BoundaMargin)-1)*((normLoc/BoundaMargin)-1);
        
        //compute the normal contribution of SmoothedField * its weight and transform it in voxels
        NCx*=(ScalProd)/dx;
        NCy*=(ScalProd)/dy;
        NCz*=(ScalProd)/dz;
        
        //remove the normal contribution
        SmoothedField->P(SmoothedField->G(0,x,y,z)-NCx,0,x,y,z);
        SmoothedField->P(SmoothedField->G(1,x,y,z)-NCy,1,x,y,z);
        SmoothedField->P(SmoothedField->G(2,x,y,z)-NCz,2,x,y,z);
        
      }
      
      if ((normLoc<=0.5*dxyz)&&(SetBoundaryToZero==1)){
        SmoothedField->P(0,0,x,y,z);
        SmoothedField->P(0,1,x,y,z);
        SmoothedField->P(0,2,x,y,z);
      }

    }
  }
  else{ //2D image
    cout << "Treatment of normal contributions in the 2D case -> TO DO" << endl;
  }
}






//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
void Cpt_Grad_ScalarField(ScalarField * SField,VectorField * Gradient,int SpecificTimeFrame,float DeltaX){
	int NBX,NBY,NBZ,NBT;
	int x,y,z,t;
	NBX=SField->NX;
	NBY=SField->NY;
	NBZ=SField->NZ;
	NBT=SField->NT;
	
	//COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
	if (SpecificTimeFrame>=0){
		//1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
		
		//1.1) allocate memory in Gradient if not done
		if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)){
			Gradient->CreateVoidField(NBX,NBY,NBZ);
			cout << "Gradient added in Cpt_Grad_ScalarField\n";
		}
		//1.2) Calculations
		t=SpecificTimeFrame;
		
		//1.2.1) gradient in direction x, y, z
		for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
			Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z);
			Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z);
			Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z);
		}
		
		//1.2.2) boundaries at 0.
		z=0;
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		z=NBZ-1;
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		y=0;
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		y=NBY-1;
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		x=0;
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
		x=NBX-1;
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
		
		//1.2.3) 2D image case
		if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
			Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0);
			Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0);
		}
	}
	else{
		//2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
		//1.1) allocate memory in Gradient if not done
		if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)||(NBT!=Gradient->NT))
			Gradient->CreateVoidField(NBX,NBY,NBZ,NBT);
		
		//1.2) Calculations
		for (t=0;t<NBT;t++){
			//gradient in direction x, y, z
			for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
				Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z,t);
				Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z,t);
				Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z,t);
			}
			
			//boundaries at 0.
			z=0;
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			z=NBZ-1;
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			y=0;
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			y=NBY-1;
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			x=0;
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
			x=NBX-1;
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
			
			//2D image case
			if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
				Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0,t);
				Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0,t);
			}
		}
	}
}




//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
void Cpt_Grad_MaskedScalarField(ScalarField * SField,VectorField * Gradient,ScalarField * Mask,int SpecificTimeFrame,float DeltaX){
	int NBX,NBY,NBZ,NBT;
	int x,y,z,t;
	NBX=SField->NX;
	NBY=SField->NY;
	NBZ=SField->NZ;
	NBT=SField->NT;
	float epsilon;
  
  epsilon=0.0001;
  
	//COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
	if (SpecificTimeFrame>=0){
		//1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
		
		//1.1) allocate memory in Gradient if not done
		if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)){
			Gradient->CreateVoidField(NBX,NBY,NBZ);
			cout << "Gradient added in Cpt_Grad_ScalarField\n";
		}
		//1.2) Calculations
		t=SpecificTimeFrame;
		
		//1.2.1) gradient in direction x, y, z
		for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
			Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z);
			Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z);
			Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z);
      
      if ((fabs(Mask->G(x,y,z,t)-Mask->G(x+1,y,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x-1,y,z,t))>epsilon)) 
        Gradient->P(0,0,x,y,z);
      
      if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y+1,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y-1,z,t))>epsilon)) 
        Gradient->P(0,1,x,y,z);
      
      if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z+1,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z-1,t))>epsilon)) 
        Gradient->P(0,2,x,y,z);
      
      if (Mask->G(x,y,z,t)<epsilon){
        Gradient->P(0,0,x,y,z);
        Gradient->P(0,1,x,y,z);
        Gradient->P(0,2,x,y,z);
      } 
    }
		
		//1.2.2) boundaries at 0.
		z=0;
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		z=NBZ-1;
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		y=0;
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		y=NBY-1;
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z);
		x=0;
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
		x=NBX-1;
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z);
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z);
		
		//1.2.3) 2D image case
		if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
			Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0);
			Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0);
		}
	}
	else{
		//2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
		//1.1) allocate memory in Gradient if not done
		if ((NBX!=Gradient->NX)||(NBY!=Gradient->NY)||(NBZ!=Gradient->NZ)||(NBT!=Gradient->NT))
			Gradient->CreateVoidField(NBX,NBY,NBZ,NBT);
		
		//1.2) Calculations
		for (t=0;t<NBT;t++){
			//gradient in direction x, y, z
			for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
				Gradient->P((SField->G(x+1,y,z,t)-SField->G(x-1,y,z,t))/(2.*DeltaX),0,x,y,z,t);
				Gradient->P((SField->G(x,y+1,z,t)-SField->G(x,y-1,z,t))/(2.*DeltaX),1,x,y,z,t);
				Gradient->P((SField->G(x,y,z+1,t)-SField->G(x,y,z-1,t))/(2.*DeltaX),2,x,y,z,t);
        
        if ((fabs(Mask->G(x,y,z,t)-Mask->G(x+1,y,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x-1,y,z,t))>epsilon)) 
          Gradient->P(0,0,x,y,z);
        
        if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y+1,z,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y-1,z,t))>epsilon)) 
          Gradient->P(0,1,x,y,z);
        
        if ((fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z+1,t))>epsilon)||(fabs(Mask->G(x,y,z,t)-Mask->G(x,y,z-1,t))>epsilon)) 
          Gradient->P(0,2,x,y,z);
        
        if (Mask->G(x,y,z,t)<epsilon){
          Gradient->P(0,0,x,y,z);
          Gradient->P(0,1,x,y,z);
          Gradient->P(0,2,x,y,z);
        } 
      }
			
			//boundaries at 0.
			z=0;
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			z=NBZ-1;
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			y=0;
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			y=NBY-1;
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Gradient->P(0.,2,x,y,z,t);
			x=0;
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
			x=NBX-1;
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,0,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,1,x,y,z,t);
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Gradient->P(0.,2,x,y,z,t);
			
			//2D image case
			if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
				Gradient->P((SField->G(x+1,y,0,t)-SField->G(x-1,y,0,t))/(2.*DeltaX),0,x,y,0,t);
				Gradient->P((SField->G(x,y+1,0,t)-SField->G(x,y-1,0,t))/(2.*DeltaX),1,x,y,0,t);
			}
		}
	}
}






//Compute (d VField(X) / d x) + (d VField(Y) / d y) + (d VField(Z) / d z) and put the result in 'GradScalVF'
//where 'VField' is a vector field and 'GradScalVF' a scalar field
void Cpt_Grad_Scal_VectorField(VectorField * VField,ScalarField * GradScalVF,int SpecificTimeFrame,float DeltaX){
	int NBX,NBY,NBZ,NBT;
	int x,y,z,t;
	float GradX,GradY,GradZ;
	
	NBX=VField->NX;
	NBY=VField->NY;
	NBZ=VField->NZ;
	NBT=VField->NT;
	
	
	//COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
	if (SpecificTimeFrame>=0){
		//1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
		
		//1.1) allocate memory in Gradient if not done
		if ((NBX!=GradScalVF->NX)||(NBY!=GradScalVF->NY)||(NBZ!=GradScalVF->NZ))
			GradScalVF->CreateVoidField(NBX,NBY,NBZ);
		
		//1.2) Calculations
		t=SpecificTimeFrame;
		
		//1.2.1) sum of gradients in direction x, y, z
		for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
			GradX=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
			GradY=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
			GradZ=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
			GradScalVF->P(GradX+GradY+GradZ,x,y,z);
		}
		
		//1.2.2) boundaries at 0.
		z=0;
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
		z=NBZ-1;
		for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
		y=0;
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
		y=NBY-1;
		for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z);
		x=0;
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z);
		x=NBX-1;
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z);
		
		//1.2.3) 2D image case
		if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
			GradX=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
			GradY=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
			GradScalVF->P(GradX+GradY,x,y,0);
		}
	}
	else{
		//2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
		//1.1) allocate memory in Gradient if not done
		if ((NBX!=GradScalVF->NX)||(NBY!=GradScalVF->NY)||(NBZ!=GradScalVF->NZ)||(NBT!=GradScalVF->NT))
			GradScalVF->CreateVoidField(NBX,NBY,NBZ,NBT);
		
		//1.2) Calculations
		for (t=0;t<NBT;t++){
			//sum of gradients in direction x, y, z
			for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
				GradX=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
				GradY=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
				GradZ=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
				GradScalVF->P(GradX+GradY+GradZ,x,y,z,t);
			}
			
			//boundaries at 0.
			z=0;
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
			z=NBZ-1;
			for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
			y=0;
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
			y=NBY-1;
			for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) GradScalVF->P(0.,x,y,z,t);
			x=0;
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z,t);
			x=NBX-1;
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) GradScalVF->P(0.,x,y,z,t);
			
			//2D image case
			if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
				GradX=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
				GradY=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
				GradScalVF->P(GradX+GradY,x,y,0,t);
			}
		}
	}
}


#ifdef COMPILE_WITH_OPENMP

//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;
  
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  NBT=VField->NT;
  
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
    
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ))
      DetJ->CreateVoidField(NBX,NBY,NBZ);
    
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //BEGIN FORK FOR THREADS
    #pragma omp parallel default(shared) private(d11,d12,d13,d21,d22,d23,d31,d32,d33,x,y,z) 
    {
      //1.2.1) sum of gradients in direction x, y, z
      #pragma omp for
      for (y=1;y<NBY-1;y++){ 
         for (z=1;z<NBZ-1;z++) for (x=1;x<NBX-1;x++){
          d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
          d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
          d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
          d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
          d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
          d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
          d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
          d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
          d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
          DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z);
        }
      }
    //END FORK FOR THREADS
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
      DetJ->P(d11*d22-d21*d12,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ)||(NBT!=DetJ->NT))
      DetJ->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //BEGIN FORK FOR THREADS
    #pragma omp parallel default(shared) private(d11,d12,d13,d21,d22,d23,d31,d32,d33,x,y,z,t) 
    {
      //1.2) Calculations
      #pragma omp for
      for (t=0;t<NBT;t++){
        //sum of gradients in direction x, y, z
        for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
          d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
          d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
          d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
          d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
          d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
          d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
          d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
          d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
          d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
          DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z,t);
        }
        
        //boundaries at 0.
        z=0;
        for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        z=NBZ-1;
        for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        y=0;
        for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        y=NBY-1;
        for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
        x=0;
        for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
        x=NBX-1;
        for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
        
        //2D image case
        if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
          d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
          d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
          d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
          d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
          DetJ->P(d11*d22-d21*d12,x,y,0,t);
        }
      }
    //END FORK FOR THREADS
    }
  }
}

#else

//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame,float DeltaX){
  int NBX,NBY,NBZ,NBT;
  int x,y,z,t;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;
  
  NBX=VField->NX;
  NBY=VField->NY;
  NBZ=VField->NZ;
  NBT=VField->NT;
  
  
  //COMPUTATIONS IN ONE TIME FRAME OR ALL OF THEM?
  if (SpecificTimeFrame>=0){
    //1) COMPUTATIONS IN ONLY ONE TIME FRAME -> 3D image returned
    
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ))
      DetJ->CreateVoidField(NBX,NBY,NBZ);
    
    //1.2) Calculations
    t=SpecificTimeFrame;
    
    //1.2.1) sum of gradients in direction x, y, z
    for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
      d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
      d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
      d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
      d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
      d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
      DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z);
    }
    
    //1.2.2) boundaries at 0.
    z=0;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    z=NBZ-1;
    for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=0;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    y=NBY-1;
    for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z);
    x=0;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    x=NBX-1;
    for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z);
    
    //1.2.3) 2D image case
    if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
      d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
      d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
      d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
      d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
      DetJ->P(d11*d22-d21*d12,x,y,0);
    }
  }
  else{
    //2) COMPUTATIONS IN ALL TIME FRAMES -> 4D image returned
    //1.1) allocate memory in Gradient if not done
    if ((NBX!=DetJ->NX)||(NBY!=DetJ->NY)||(NBZ!=DetJ->NZ)||(NBT!=DetJ->NT))
      DetJ->CreateVoidField(NBX,NBY,NBZ,NBT);
    
    //1.2) Calculations
    for (t=0;t<NBT;t++){
      //sum of gradients in direction x, y, z
      for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        d11=(VField->G(0,x+1,y,z,t)-VField->G(0,x-1,y,z,t))/(2.*DeltaX);
        d12=(VField->G(0,x,y+1,z,t)-VField->G(0,x,y-1,z,t))/(2.*DeltaX);
        d13=(VField->G(0,x,y,z+1,t)-VField->G(0,x,y,z-1,t))/(2.*DeltaX);
        d21=(VField->G(1,x+1,y,z,t)-VField->G(1,x-1,y,z,t))/(2.*DeltaX);
        d22=(VField->G(1,x,y+1,z,t)-VField->G(1,x,y-1,z,t))/(2.*DeltaX);
        d23=(VField->G(1,x,y,z+1,t)-VField->G(1,x,y,z-1,t))/(2.*DeltaX);
        d31=(VField->G(2,x+1,y,z,t)-VField->G(2,x-1,y,z,t))/(2.*DeltaX);
        d32=(VField->G(2,x,y+1,z,t)-VField->G(2,x,y-1,z,t))/(2.*DeltaX);
        d33=(VField->G(2,x,y,z+1,t)-VField->G(2,x,y,z-1,t))/(2.*DeltaX);
        DetJ->P(d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13),x,y,z,t);
      }
      
      //boundaries at 0.
      z=0;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      z=NBZ-1;
      for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      y=0;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      y=NBY-1;
      for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) DetJ->P(1.,x,y,z,t);
      x=0;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
      x=NBX-1;
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) DetJ->P(1.,x,y,z,t);
      
      //2D image case
      if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++){
        d11=(VField->G(0,x+1,y,0,t)-VField->G(0,x-1,y,0,t))/(2.*DeltaX);
        d12=(VField->G(0,x,y+1,0,t)-VField->G(0,x,y-1,0,t))/(2.*DeltaX);
        d21=(VField->G(1,x+1,y,0,t)-VField->G(1,x-1,y,0,t))/(2.*DeltaX);
        d22=(VField->G(1,x,y+1,0,t)-VField->G(1,x,y-1,0,t))/(2.*DeltaX);
        DetJ->P(d11*d22-d21*d12,x,y,0,t);
      }
    }
  }
}

#endif






#ifdef COMPILE_WITH_OPENMP

///Integrate a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptDefFromSteadyVeloField(VectorField * VeloField,VectorField * DeformationField,int log2TimeStepNb){
  int x,y,z;
  int i;
  float VecTempX,VecTempY,VecTempZ;
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=1/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  cout << "totot" << endl;
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,VecTempX,VecTempY,VecTempZ,i) 
  {
    //initiate DeformationField
    #pragma omp for
    for (y=0;y<VeloField->NY;y++) {
      for (z=0;z<VeloField->NZ;z++) for (x=0;x<VeloField->NX;x++){
        DeformationField->P(VeloField->G(0,x,y,z)*coef,0,x,y,z);
        DeformationField->P(VeloField->G(1,x,y,z)*coef,1,x,y,z);
        DeformationField->P(VeloField->G(2,x,y,z)*coef,2,x,y,z);
      }
    }
    
    //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
    for (i=0;i<log2TimeStepNb;i++){
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          VecTempX=x+DeformationField->G(0,x,y,z);
          VecTempY=y+DeformationField->G(1,x,y,z);
          VecTempZ=z+DeformationField->G(2,x,y,z);
          
          TempVeloField1.P(DeformationField->G(0,x,y,z)+DeformationField->G(0,VecTempX,VecTempY,VecTempZ),0,x,y,z);
          TempVeloField1.P(DeformationField->G(1,x,y,z)+DeformationField->G(1,VecTempX,VecTempY,VecTempZ),1,x,y,z);
          TempVeloField1.P(DeformationField->G(2,x,y,z)+DeformationField->G(2,VecTempX,VecTempY,VecTempZ),2,x,y,z);
        }
      }
      
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          DeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
          DeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
          DeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///Integrate a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptDefFromSteadyVeloField(VectorField * VeloField,VectorField * DeformationField,int log2TimeStepNb){
  int x,y,z;
  int i;
  float VecTemp[3];
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=1/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  
  //initiate DeformationField
  for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
    DeformationField->P(VeloField->G(0,x,y,z)*coef,0,x,y,z);
    DeformationField->P(VeloField->G(1,x,y,z)*coef,1,x,y,z);
    DeformationField->P(VeloField->G(2,x,y,z)*coef,2,x,y,z);
  }
  
  //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
  for (i=0;i<log2TimeStepNb;i++){
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      VecTemp[0]=x+DeformationField->G(0,x,y,z);
      VecTemp[1]=y+DeformationField->G(1,x,y,z);
      VecTemp[2]=z+DeformationField->G(2,x,y,z);
      
      TempVeloField1.P(DeformationField->G(0,x,y,z)+DeformationField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
      TempVeloField1.P(DeformationField->G(1,x,y,z)+DeformationField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
      TempVeloField1.P(DeformationField->G(2,x,y,z)+DeformationField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
    }
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      DeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
      DeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
      DeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
    }
  }
}

#endif

#ifdef COMPILE_WITH_OPENMP

///Integrate the inverse of a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptInvDefFromSteadyVeloField(VectorField * VeloField,VectorField * InvDeformationField,int log2TimeStepNb){
  int x,y,z;
  int i;
  float VecTempX,VecTempY,VecTempZ;
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=1/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,VecTempX,VecTempY,VecTempZ,i) 
  {
    //initiate InvDeformationField
    #pragma omp for
    for (y=0;y<VeloField->NY;y++){
      for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
        InvDeformationField->P(-VeloField->G(0,x,y,z)*coef,0,x,y,z);
        InvDeformationField->P(-VeloField->G(1,x,y,z)*coef,1,x,y,z);
        InvDeformationField->P(-VeloField->G(2,x,y,z)*coef,2,x,y,z);
      }
    }
    
    //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
    for (i=0;i<log2TimeStepNb;i++){
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          VecTempX=x+InvDeformationField->G(0,x,y,z);
          VecTempY=y+InvDeformationField->G(1,x,y,z);
          VecTempZ=z+InvDeformationField->G(2,x,y,z);
          
          TempVeloField1.P(InvDeformationField->G(0,x,y,z)+InvDeformationField->G(0,VecTempX,VecTempY,VecTempZ),0,x,y,z);
          TempVeloField1.P(InvDeformationField->G(1,x,y,z)+InvDeformationField->G(1,VecTempX,VecTempY,VecTempZ),1,x,y,z);
          TempVeloField1.P(InvDeformationField->G(2,x,y,z)+InvDeformationField->G(2,VecTempX,VecTempY,VecTempZ),2,x,y,z);
        }
      }
      
      #pragma omp for
      for (y=0;y<VeloField->NY;y++){
        for (z=0;z<VeloField->NZ;z++)  for (x=0;x<VeloField->NX;x++){
          InvDeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
          InvDeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
          InvDeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else


///Integrate the inverse of a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptInvDefFromSteadyVeloField(VectorField * VeloField,VectorField * InvDeformationField,int log2TimeStepNb){
  int x,y,z;
  int i;
  float VecTemp[3];
  float coef;
  VectorField TempVeloField1;
  
  //init
  TempVeloField1.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  coef=1/(pow(static_cast<float>(2),static_cast<float>(log2TimeStepNb)));
  
  //initiate InvDeformationField
  for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
    InvDeformationField->P(-VeloField->G(0,x,y,z)*coef,0,x,y,z);
    InvDeformationField->P(-VeloField->G(1,x,y,z)*coef,1,x,y,z);
    InvDeformationField->P(-VeloField->G(2,x,y,z)*coef,2,x,y,z);
  }
  
  //integrate the velocity field in 2^log2TimeStepNb time steps using the technique of Arsigny MICCAI 2006
  for (i=0;i<log2TimeStepNb;i++){
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      VecTemp[0]=x+InvDeformationField->G(0,x,y,z);
      VecTemp[1]=y+InvDeformationField->G(1,x,y,z);
      VecTemp[2]=z+InvDeformationField->G(2,x,y,z);
      
      TempVeloField1.P(InvDeformationField->G(0,x,y,z)+InvDeformationField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
      TempVeloField1.P(InvDeformationField->G(1,x,y,z)+InvDeformationField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
      TempVeloField1.P(InvDeformationField->G(2,x,y,z)+InvDeformationField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
    }
    for (z=0;z<VeloField->NZ;z++) for (y=0;y<VeloField->NY;y++) for (x=0;x<VeloField->NX;x++){
      InvDeformationField->P(TempVeloField1.G(0,x,y,z),0,x,y,z);
      InvDeformationField->P(TempVeloField1.G(1,x,y,z),1,x,y,z);
      InvDeformationField->P(TempVeloField1.G(2,x,y,z),2,x,y,z);
    }
  }
}

#endif



///compute the Lie Braket of VF1 and VF2. Put the result in VF3.   (subfunction of ComposeTwoLogFieldsUsingBCH)
void LieBracket(VectorField * VF1,VectorField * VF2,VectorField * VF3){
  int x,y,z;
  //put to 0 the boundaries of VF3
  x=0;
  for (z=0;z<VF1->NZ;z++) for (y=0;y<VF1->NY;y++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  x=VF1->NX-1;
  for (z=0;z<VF1->NZ;z++) for (y=0;y<VF1->NY;y++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  y=0;
  for (z=0;z<VF1->NZ;z++)  for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  y=VF1->NY-1;
  for (z=0;z<VF1->NZ;z++)  for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  z=0;
  for (y=0;y<VF1->NY;y++) for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  z=VF1->NZ-1;
  for (y=0;y<VF1->NY;y++) for (x=0;x<VF1->NX;x++){
    VF3->P(0,0,x,y,z); VF3->P(0,1,x,y,z); VF3->P(0,2,x,y,z);
  }
  
  //treat the rest of the data
  for (z=1;z<VF1->NZ-1;z++) for (y=1;y<VF1->NY-1;y++) for (x=1;x<VF1->NX-1;x++){
    VF3->P(-((VF2->G(0,x+1,y,z)-VF2->G(0,x-1,y,z))*VF1->G(0,x,y,z))-(VF2->G(0,x,y+1,z)-(VF2->G(0,x,y-1,z))*VF1->G(1,x,y,z))-((VF2->G(0,x,y,z+1)-VF2->G(0,x,y,z-1))*VF1->G(2,x,y,z)),0,x,y,z);
    VF3->P(-((VF2->G(1,x+1,y,z)-VF2->G(1,x-1,y,z))*VF1->G(0,x,y,z))-(VF2->G(1,x,y+1,z)-(VF2->G(1,x,y-1,z))*VF1->G(1,x,y,z))-((VF2->G(1,x,y,z+1)-VF2->G(1,x,y,z-1))*VF1->G(2,x,y,z)),1,x,y,z);
    VF3->P(-((VF2->G(2,x+1,y,z)-VF2->G(2,x-1,y,z))*VF1->G(0,x,y,z))-(VF2->G(2,x,y+1,z)-(VF2->G(2,x,y-1,z))*VF1->G(1,x,y,z))-((VF2->G(2,x,y,z+1)-VF2->G(2,x,y,z-1))*VF1->G(2,x,y,z)),2,x,y,z);
    
    VF3->Add(((VF1->G(0,x+1,y,z)-VF1->G(0,x-1,y,z))*VF2->G(0,x,y,z))+((VF1->G(0,x,y+1,z)-VF1->G(0,x,y-1,z))*VF2->G(1,x,y,z))+((VF1->G(0,x,y,z+1)-VF1->G(0,x,y,z-1))*VF2->G(2,x,y,z)),0,x,y,z);
    VF3->Add(((VF1->G(1,x+1,y,z)-VF1->G(1,x-1,y,z))*VF2->G(0,x,y,z))+((VF1->G(1,x,y+1,z)-VF1->G(1,x,y-1,z))*VF2->G(1,x,y,z))+((VF1->G(1,x,y,z+1)-VF1->G(1,x,y,z-1))*VF2->G(2,x,y,z)),1,x,y,z);
    VF3->Add(((VF1->G(2,x+1,y,z)-VF1->G(2,x-1,y,z))*VF2->G(0,x,y,z))+((VF1->G(2,x,y+1,z)-VF1->G(2,x,y-1,z))*VF2->G(1,x,y,z))+((VF1->G(2,x,y,z+1)-VF1->G(2,x,y,z-1))*VF2->G(2,x,y,z)),2,x,y,z);
    
    VF3->P(VF3->G(0,x,y,z)/2,0,x,y,z);
    VF3->P(VF3->G(1,x,y,z)/2,1,x,y,z);
    VF3->P(VF3->G(2,x,y,z)/2,2,x,y,z);
  }
}






///RefVeloField and UpdateVeloField are two steady velocity fields. Their exponentials are respectively the current deformation
///and the update of the deformation. This function approximates the velocity field which is the log of the composition 
///between the current deformation and the update deformation (Baker-Campbell-Hausdorff formula) (cf Vercauteren MICCAI 2008)
///In output, RefVeloField is the updated velocity field. UpdateVeloField is also modified for algorithmic reasons but it
///represents nothing pertinent as an output.
void ComposeTwoLogFieldsUsingBCH(VectorField * RefVeloField,VectorField * UpdateVeloField){
  int x,y,z;
  VectorField TempVeloField1;
  VectorField TempVeloField2;
  
  //init
  TempVeloField1.CreateVoidField(RefVeloField->NX,RefVeloField->NY,RefVeloField->NZ);
  TempVeloField2.CreateVoidField(RefVeloField->NX,RefVeloField->NY,RefVeloField->NZ);
  
  
  //compute LieBracket(RefVeloField,UpdateVeloField)
  LieBracket(RefVeloField,UpdateVeloField,&TempVeloField1);
  
  //compute LieBracket(RefVeloField,LieBracket(RefVeloField,UpdateVeloField))
  LieBracket(RefVeloField,&TempVeloField1,&TempVeloField2);
  
  
  
  // estimate the new RefVeloField
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(0,x,y,z)+(TempVeloField1.G(0,x,y,z)/2)+(TempVeloField2.G(0,x,y,z)/12),0,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(1,x,y,z)+(TempVeloField1.G(1,x,y,z)/2)+(TempVeloField2.G(1,x,y,z)/12),1,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(2,x,y,z)+(TempVeloField1.G(2,x,y,z)/2)+(TempVeloField2.G(2,x,y,z)/12),2,x,y,z);
}


///RefVeloField and UpdateVeloField are two steady velocity fields. Their exponentials are respectively the current deformation
///and the update of the deformation. This function approximates the velocity field which is the log of the composition 
///between the current deformation and the update deformation
///In output, RefVeloField is the updated velocity field. UpdateVeloField is also modified for algorithmic reasons but it
///represents nothing pertinent as an output.
void ComposeTwoLogFieldsUsingSum(VectorField * RefVeloField,VectorField * UpdateVeloField){
  int x,y,z;
  
  
  // estimate the new RefVeloField
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(0,x,y,z),0,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(1,x,y,z),1,x,y,z);
  
  for (z=0;z<RefVeloField->NZ;z++) for (y=0;y<RefVeloField->NY;y++) for (x=0;x<RefVeloField->NX;x++)
    RefVeloField->Add(UpdateVeloField->G(2,x,y,z),2,x,y,z);
}



///Compose RefField with UpdateField. The result is saved in RefField
void DisplacementFieldCompose(VectorField * RefField,VectorField * UpdateField){
  int x,y,z;
  float VecTemp[3];
  
  for (z=0;z<RefField->NZ;z++) for (y=0;y<RefField->NY;y++) for (x=0;x<RefField->NX;x++){
    VecTemp[0]=x+RefField->G(0,x,y,z);
    VecTemp[1]=y+RefField->G(1,x,y,z);
    VecTemp[2]=z+RefField->G(2,x,y,z);
    
    RefField->Add(UpdateField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
    RefField->Add(UpdateField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
    RefField->Add(UpdateField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
  }
}


///Compose InvUpdateField with InvRefField. The result is saved in InvRefField (and InvRefField)
///(we consider that an inv. disp. f. points from the target source to the source)
void InvDisplacementFieldCompose(VectorField * InvUpdateField,VectorField * InvRefField){
  int x,y,z;
  float VecTemp[3];
  
  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++){
    VecTemp[0]=x+InvUpdateField->G(0,x,y,z);
    VecTemp[1]=y+InvUpdateField->G(1,x,y,z);
    VecTemp[2]=z+InvUpdateField->G(2,x,y,z);
    
    InvUpdateField->Add(InvRefField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]),0,x,y,z);
    InvUpdateField->Add(InvRefField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]),1,x,y,z);
    InvUpdateField->Add(InvRefField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]),2,x,y,z);
  }

  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++) InvRefField->P(InvUpdateField->G(0,x,y,z),0,x,y,z);
  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++) InvRefField->P(InvUpdateField->G(1,x,y,z),1,x,y,z);
  for (z=0;z<InvUpdateField->NZ;z++) for (y=0;y<InvUpdateField->NY;y++) for (x=0;x<InvUpdateField->NX;x++) InvRefField->P(InvUpdateField->G(2,x,y,z),2,x,y,z);

}




///We consider here that the vectors of DispField point from 'StaticImage' to 'DeformedImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
///remark: 'DefoField' has to be sufficiently smooth to allow an accurate inversion
void ProjectImageUsingDispField(VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage, int NearestNgbh){
  int x,y,z;
  float VecTemp[3];
  float VecTemp2[3];
  int NbIt;
  int i;
  
  NbIt=3;
  
  for (z=0;z<DispField->NZ;z++) for (y=0;y<DispField->NY;y++) for (x=0;x<DispField->NX;x++){
    VecTemp[0]=x-DispField->G(0,x,y,z);
    VecTemp[1]=y-DispField->G(1,x,y,z);
    VecTemp[2]=z-DispField->G(2,x,y,z);
    
    for (i=0; i<NbIt; i++) {
      VecTemp2[0]=x-DispField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[1]=y-DispField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[2]=z-DispField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]);
      
      VecTemp[0]=VecTemp2[0];
      VecTemp[1]=VecTemp2[1];
      VecTemp[2]=VecTemp2[2];
    }
    
    if (NearestNgbh!=1) DeformedImage->P(StaticImage->G(VecTemp[0],VecTemp[1],VecTemp[2]),x,y,z);
    else DeformedImage->P(StaticImage->G(static_cast<int>(VecTemp[0]+0.5),static_cast<int>(VecTemp[1]+0.5),static_cast<int>(VecTemp[2]+0.5)),x,y,z);
  }
}


///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndDispField(float ProjectCS_2_OriginCS[4][4],VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage){
  int x,y,z;
  float VecTemp[3];
  float VecTemp2[3];
  int NbIt;
  int i;
  
  NbIt=3;
  
  for (z=0;z<DispField->NZ;z++) for (y=0;y<DispField->NY;y++) for (x=0;x<DispField->NX;x++){
    VecTemp[0]=x-DispField->G(0,x,y,z);
    VecTemp[1]=y-DispField->G(1,x,y,z);
    VecTemp[2]=z-DispField->G(2,x,y,z);
    
    for (i=0; i<NbIt; i++) {
      VecTemp2[0]=x-DispField->G(0,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[1]=y-DispField->G(1,VecTemp[0],VecTemp[1],VecTemp[2]);
      VecTemp2[2]=z-DispField->G(2,VecTemp[0],VecTemp[1],VecTemp[2]);
      
      VecTemp[0]=VecTemp2[0];
      VecTemp[1]=VecTemp2[1];
      VecTemp[2]=VecTemp2[2];
    }
    
    DeformedImage->P(StaticImage->G(ProjectCS_2_OriginCS,VecTemp[0],VecTemp[1],VecTemp[2]),x,y,z);
  }
}




///We consider here that the vectors of InvDispField point from 'DeformedImage' to 'StaticImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
void ProjectImageUsingInvDispField(VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh){
  int x,y,z;
  float VecTemp[3];
  
  
  //project the image
  for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
    VecTemp[0]=x+InvDispField->G(0,x,y,z);
    VecTemp[1]=y+InvDispField->G(1,x,y,z);
    VecTemp[2]=z+InvDispField->G(2,x,y,z);
    
    
    if (NearestNgbh!=1) DeformedImage->P(StaticImage->G(VecTemp[0],VecTemp[1],VecTemp[2]),x,y,z);
    else DeformedImage->P(StaticImage->G(static_cast<int>(VecTemp[0]+0.5),static_cast<int>(VecTemp[1]+0.5),static_cast<int>(VecTemp[2]+0.5)),x,y,z);
  }
  
}



#ifdef COMPILE_WITH_OPENMP

///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndInvDispField(float ProjectCS_2_OriginCS[4][4],VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  int x,y,z;
  float x2,y2,z2;
  int x3,y3,z3;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,x2,y2,z2) 
  {
    #pragma omp for
    for (y=0;y<DeformedImage->NY;y++){
      for (z=0;z<DeformedImage->NZ;z++)  for (x=0;x<DeformedImage->NX;x++){
        x2=x+InvDispField->G(0,x,y,z);
        y2=y+InvDispField->G(1,x,y,z);
        z2=z+InvDispField->G(2,x,y,z);
        
        DeformedImage->P(StaticImage->G(ProjectCS_2_OriginCS,x2,y2,z2,0,NN),x,y,z);
      }
    }
  //END FORK FOR THREADS
  }
}

#else

///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndInvDispField(float ProjectCS_2_OriginCS[4][4],VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  int x,y,z;
  float x2,y2,z2;
  
  
  //project the image
  for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
    x2=x+InvDispField->G(0,x,y,z);
    y2=y+InvDispField->G(1,x,y,z);
    z2=z+InvDispField->G(2,x,y,z);
    
    DeformedImage->P(StaticImage->G(ProjectCS_2_OriginCS,x2,y2,z2,0,NN),x,y,z);
  }
}

#endif




///... explicit name ... here DispField represents the mapping from the target c.s to the source c.s. in voxels
void ProjectImageUsingDispFieldAndInvDispField(VectorField * DispField,VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  int x,y,z;
  float x2,y2,z2;
  int x3,y3,z3;
  float VecTemp[3];
  
  //project the image using ...
  if (NN==0){  //... trilinear interpolation
    for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
      x2=x+InvDispField->G(0,x,y,z);
      y2=y+InvDispField->G(1,x,y,z);
      z2=z+InvDispField->G(2,x,y,z);
      
      DeformedImage->P(StaticImage->G(DispField->G(0,x2,y2,z2),DispField->G(1,x2,y2,z2),DispField->G(2,x2,y2,z2)),x,y,z);
    }
  }
  else{ //... nearest neighbor interpolation
    for (z=0;z<DeformedImage->NZ;z++) for (y=0;y<DeformedImage->NY;y++) for (x=0;x<DeformedImage->NX;x++){
      x2=x+InvDispField->G(0,x,y,z);
      y2=y+InvDispField->G(1,x,y,z);
      z2=z+InvDispField->G(2,x,y,z);
      
      x3=static_cast<int>(DispField->G(0,x2,y2,z2)+0.5);
      y3=static_cast<int>(DispField->G(1,x2,y2,z2)+0.5);
      z3=static_cast<int>(DispField->G(2,x2,y2,z2)+0.5);
      
      if (x3<0) x3=0; if (x3>=StaticImage->NX) x3=StaticImage->NX-1; 
      if (y3<0) y3=0; if (y3>=StaticImage->NY) y3=StaticImage->NY-1; 
      if (z3<0) z3=0; if (z3>=StaticImage->NZ) z3=StaticImage->NZ-1; 
      
      DeformedImage->P(StaticImage->G(x3,y3,z3),x,y,z);
    }
  }
}







///... explicit name ....
void ProjectImageUsingSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage, int NearestNgbh){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptInvDefFromSteadyVeloField(VeloField,&TempField,5);
  
  //project the image
  ProjectImageUsingInvDispField(&TempField,StaticImage,DeformedImage,NearestNgbh);
  
}



///... explicit name ....
void ProjectImageUsingAffineTransfoAndSteadyVeloField(float ProjectCS_2_OriginCS[4][4],VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptInvDefFromSteadyVeloField(VeloField,&TempField,5);
  
  //project the image
  ProjectImageUsingAffineTransfoAndInvDispField(ProjectCS_2_OriginCS,&TempField,StaticImage,DeformedImage,NN);
}


///... explicit name .... here DispField represents the mapping from the target c.s to the source c.s. in voxels
void ProjectImageUsingDispFieldAndSteadyVeloField(VectorField * DispField,VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptInvDefFromSteadyVeloField(VeloField,&TempField,5);
  
  //project the image
  ProjectImageUsingDispFieldAndInvDispField(DispField,&TempField,StaticImage,DeformedImage,NN);
  
}







///... explicit name ....
void ProjectImageUsingInverseSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage, int NearestNgbh){
  VectorField TempField;
  
  //compute the inverse mapping
  TempField.CreateVoidField(VeloField->NX,VeloField->NY,VeloField->NZ);
  CptDefFromSteadyVeloField(VeloField,&TempField,5);
  
  //project the image
  ProjectImageUsingInvDispField(&TempField,StaticImage,DeformedImage,NearestNgbh);

}


///Image transformation using the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfo(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,ScalarField * TransformedImage){
	int x,y,z;
	
	for (z = 0; z < TransformedImage->NZ; z++) for (y = 0; y < TransformedImage->NY; y++) for (x = 0; x < TransformedImage->NX; x++){
		TransformedImage->P(ImagToPropag->G(ProjectCS_2_OriginCS,x,y,z),x,y,z);
	}
}


//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'
void CptMappingFromVeloField(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps,float DeltaX){
	float VecTemp[3];
	float VecTemp2[3];
	int NBX,NBY,NBZ,NBT;
	int i,x,y,z;
	int t;
	float DeltaT_div_DeltaX;
	
	//0) INITIALISATION
	
	//initialisation
	NBX=VeloField->NX;
	NBY=VeloField->NY;
	NBZ=VeloField->NZ;
	NBT=VeloField->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
	//allocate memory in GradScalVF if not done
	if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ)||(NBT!=Map->NT))
		Map->CreateVoidField(NBX,NBY,NBZ,NBT);
	
	
	//1) MAPPING AT THE REFERENCE TIME SUBDIVISION
	if ((refTimeStep<0)||(refTimeStep>NBT-1)) refTimeStep=0;
	
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
		Map->P(MappingAtRefTimeStep->G(0,x,y,z),0,x,y,z,refTimeStep);
		Map->P(MappingAtRefTimeStep->G(1,x,y,z),1,x,y,z,refTimeStep);
		Map->P(MappingAtRefTimeStep->G(2,x,y,z),2,x,y,z,refTimeStep);
	}
	
	
	
	//2) FORWARD MAPPING for the time subdivisions > refTimeStep
	for (t=refTimeStep+1;t<NBT;t++){
		if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				VecTemp[0]=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[1]=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[2]=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				
				//find the original coordinates
				Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
				Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
				Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
			}
		}
		else{ // leap frog scheme
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				//init
				VecTemp[0]=0.; 
				VecTemp[1]=0.;
				VecTemp[2]=0.;
				
				//convergence
				for (i=0;i<ConvergenceSteps;i++){
					VecTemp2[0]=VeloField->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					VecTemp2[1]=VeloField->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					VecTemp2[2]=VeloField->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					
					VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				}
				
				//find the original coordinates
				Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
				Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
				Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
			}
		}
	}
	
	//3) BACKWARD MAPPING for the time subdivisions < refTimeStep
	for (t=refTimeStep-1;t>=0;t--){
		if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				VecTemp[0]=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[1]=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[2]=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				
				Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
				Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
				Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
			}
		}
		else{ // leap frog scheme
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				//init
				VecTemp[0]=0.; 
				VecTemp[1]=0.;
				VecTemp[2]=0.;
				
				//convergence
				for (i=0;i<ConvergenceSteps;i++){
					VecTemp2[0]=VeloField->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
					VecTemp2[1]=VeloField->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
					VecTemp2[2]=VeloField->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
					
					VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				}
				
				//find the original coordinates
				Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
				Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
				Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
			}
		}
	}
}




#ifdef COMPILE_WITH_OPENMP

//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'
void CptMappingFromVeloField_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps,float DeltaX){
  float VecTempX,VecTempY,VecTempZ;
  int NBX,NBY,NBZ,NBT;
  int x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  
  //0) INITIALISATION
  
  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //allocate memory in GradScalVF if not done
  if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ)||(NBT!=Map->NT))
    Map->CreateVoidField(NBX,NBY,NBZ,NBT);
  

  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,VecTempX,VecTempY,VecTempZ,t) 
  {
    //1) MAPPING AT THE REFERENCE TIME SUBDIVISION
    if ((refTimeStep<0)||(refTimeStep>NBT-1)) refTimeStep=0;
    
    #pragma omp for
    for (y=0;y<NBY;y++){
      for (z=0;z<NBZ;z++)  for (x=0;x<NBX;x++){
        Map->P(static_cast<float>(x),0,x,y,z,refTimeStep);
        Map->P(static_cast<float>(y),1,x,y,z,refTimeStep);
        Map->P(static_cast<float>(z),2,x,y,z,refTimeStep);
      }
    }
    
    //2) FORWARD MAPPING for the time subdivisions > refTimeStep
    for (t=refTimeStep+1;t<NBT;t++){
      #pragma omp for
      for (y=0;y<NBY;y++){ 
        for (z=0;z<NBZ;z++)  for (x=0;x<NBX;x++){
          VecTempX=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempY=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempZ=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
          
          //find the original coordinates
          Map->P(Map->G(0,x-VecTempX,y-VecTempY,z-VecTempZ,t-1),0,x,y,z,t);
          Map->P(Map->G(1,x-VecTempX,y-VecTempY,z-VecTempZ,t-1),1,x,y,z,t);
          Map->P(Map->G(2,x-VecTempX,y-VecTempY,z-VecTempZ,t-1),2,x,y,z,t);
        }
      }
    }
    
    //3) BACKWARD MAPPING for the time subdivisions < refTimeStep
    for (t=refTimeStep-1;t>=0;t--){
      #pragma omp for
      for (y=0;y<NBY;y++){
        for (z=0;z<NBZ;z++)  for (x=0;x<NBX;x++){
          VecTempX=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempY=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTempZ=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
          
          Map->P(Map->G(0,x+VecTempX,y+VecTempY,z+VecTempZ,t+1),0,x,y,z,t);
          Map->P(Map->G(1,x+VecTempX,y+VecTempY,z+VecTempZ,t+1),1,x,y,z,t);
          Map->P(Map->G(2,x+VecTempX,y+VecTempY,z+VecTempZ,t+1),2,x,y,z,t);
        }
      }
    }
  //END FORK FOR THREADS
  }
}

#else

//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'
void CptMappingFromVeloField_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps,float DeltaX){
  float VecTemp[3];
  float VecTemp2[3];
  int NBX,NBY,NBZ,NBT;
  int i,x,y,z;
  int t;
  float DeltaT_div_DeltaX;
  
  //0) INITIALISATION
  
  //initialisation
  NBX=VeloField->NX;
  NBY=VeloField->NY;
  NBZ=VeloField->NZ;
  NBT=VeloField->NT;
  DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
  
  //allocate memory in GradScalVF if not done
  if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ)||(NBT!=Map->NT))
    Map->CreateVoidField(NBX,NBY,NBZ,NBT);
  
  
  //1) MAPPING AT THE REFERENCE TIME SUBDIVISION
  if ((refTimeStep<0)||(refTimeStep>NBT-1)) refTimeStep=0;
  
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    Map->P(static_cast<float>(x),0,x,y,z,refTimeStep);
    Map->P(static_cast<float>(y),1,x,y,z,refTimeStep);
    Map->P(static_cast<float>(z),2,x,y,z,refTimeStep);
  }
  
  
  
  //2) FORWARD MAPPING for the time subdivisions > refTimeStep
  for (t=refTimeStep+1;t<NBT;t++){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        //find the original coordinates
        Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[1]=VeloField->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          VecTemp2[2]=VeloField->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Map->P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),0,x,y,z,t);
        Map->P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),1,x,y,z,t);
        Map->P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1),2,x,y,z,t);
      }
    }
  }
  
  //3) BACKWARD MAPPING for the time subdivisions < refTimeStep
  for (t=refTimeStep-1;t>=0;t--){
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        VecTemp[0]=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[1]=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
        VecTemp[2]=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        
        Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
    else{ // leap frog scheme
      for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
        //init
        VecTemp[0]=0.; 
        VecTemp[1]=0.;
        VecTemp[2]=0.;
        
        //convergence
        for (i=0;i<ConvergenceSteps;i++){
          VecTemp2[0]=VeloField->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[1]=VeloField->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          VecTemp2[2]=VeloField->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
          
          VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
          VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
        }
        
        //find the original coordinates
        Map->P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),0,x,y,z,t);
        Map->P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),1,x,y,z,t);
        Map->P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1),2,x,y,z,t);
      }
    }
  }
}

#endif



//Do the same thing as the initialisation of CptMappingFromVeloField -> load the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
//The function 'CptMappingFromVeloField2_Increment' then allows to increment the field backward or forward according to 'VeloField'
void CptMappingFromVeloField2_Init(VectorField * MappingAtRefTimeStep,VectorField * Map){
	int NBX,NBY,NBZ;
	int x,y,z;
	
	//0) INITIALISATION
	
	//initialisation
	NBX=MappingAtRefTimeStep->NX;
	NBY=MappingAtRefTimeStep->NY;
	NBZ=MappingAtRefTimeStep->NZ;
	
	//allocate memory in GradScalVF if not done
	if ((NBX!=Map->NX)||(NBY!=Map->NY)||(NBZ!=Map->NZ))
		Map->CreateVoidField(NBX,NBY,NBZ);
	
	
	//1) MAPPING AT THE REFERENCE TIME SUBDIVISION
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
		Map->P(MappingAtRefTimeStep->G(0,x,y,z),0,x,y,z);
		Map->P(MappingAtRefTimeStep->G(1,x,y,z),1,x,y,z);
		Map->P(MappingAtRefTimeStep->G(2,x,y,z),2,x,y,z);
	}
	
}



//same as CptMappingFromVeloField2_Init with an identity mapping
void CptMappingFromVeloField2_Init_IniIdMap(VectorField * Map){
	int x,y,z;
	
	for (z=0;z<Map->NZ;z++) for (y=0;y<Map->NY;y++) for (x=0;x<Map->NX;x++){
		Map->P(static_cast<float>(x),0,x,y,z);
		Map->P(static_cast<float>(y),1,x,y,z);
		Map->P(static_cast<float>(z),2,x,y,z);
	}
}




//Do the same thing as an incrementation of CptMappingFromVeloField from 'CurrentTimeStep' and backward (BackwardOrForward==-1) or forward (BackwardOrForward==1) 
//The function 'CptMappingFromVeloField2_Init' is then supposed to have loaded the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
void CptMappingFromVeloField2_Increment(VectorField * VeloField,VectorField * Map,int CurrentTimeStep,int BackwardOrForward,int ConvergenceSteps,float DeltaX){
	float VecTemp[3];
	float VecTemp2[3];
	int NBX,NBY,NBZ,NBT;
	int i,x,y,z;
	int t;
	float DeltaT_div_DeltaX;
	VectorField NewMap;
  
	//1) INITIALISATION
	
	//initialisation
	NBX=VeloField->NX;
	NBY=VeloField->NY;
	NBZ=VeloField->NZ;
	NBT=VeloField->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
  NewMap.CreateVoidField(NBX,NBY,NBZ);
  
  t=CurrentTimeStep;
  
	
  if (BackwardOrForward==1){  //2.1) FORWARD MAPPING for the time subdivisions > refTimeStep
    if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				VecTemp[0]=(VeloField->G(0,x,y,z,t-1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[1]=(VeloField->G(1,x,y,z,t-1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[2]=(VeloField->G(2,x,y,z,t-1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				
				//find the original coordinates
				NewMap.P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),0,x,y,z);
				NewMap.P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),1,x,y,z);
				NewMap.P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),2,x,y,z);
			}
		}
		else{ // leap frog scheme
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				//init
				VecTemp[0]=0.; 
				VecTemp[1]=0.;
				VecTemp[2]=0.;
				
				//convergence
				for (i=0;i<ConvergenceSteps;i++){
					VecTemp2[0]=VeloField->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					VecTemp2[1]=VeloField->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					VecTemp2[2]=VeloField->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					
					VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				}
				
				//find the original coordinates
				NewMap.P(Map->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),0,x,y,z);
				NewMap.P(Map->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),1,x,y,z);
				NewMap.P(Map->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2]),2,x,y,z);
			}
		}
	}
	else{	//2.2) BACKWARD MAPPING for the time subdivisions < refTimeStep
		if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				VecTemp[0]=(VeloField->G(0,x,y,z,t+1)+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[1]=(VeloField->G(1,x,y,z,t+1)+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[2]=(VeloField->G(2,x,y,z,t+1)+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				
				NewMap.P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),0,x,y,z);
				NewMap.P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),1,x,y,z);
				NewMap.P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),2,x,y,z);
			}
		}
		else{ // leap frog scheme
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				//init
				VecTemp[0]=0.; 
				VecTemp[1]=0.;
				VecTemp[2]=0.;
				
				//convergence
				for (i=0;i<ConvergenceSteps;i++){
					VecTemp2[0]=VeloField->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
					VecTemp2[1]=VeloField->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
					VecTemp2[2]=VeloField->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
					
					VecTemp[0]=(VecTemp2[0]+VeloField->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[1]=(VecTemp2[1]+VeloField->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[2]=(VecTemp2[2]+VeloField->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				}
				
				//find the original coordinates
				NewMap.P(Map->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),0,x,y,z);
				NewMap.P(Map->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),1,x,y,z);
				NewMap.P(Map->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),2,x,y,z);
			}
		}
	}
  
  
  //3) Compute the new map
  for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
    Map->P(NewMap.G(0,x,y,z),0,x,y,z);
    Map->P(NewMap.G(1,x,y,z),1,x,y,z);
    Map->P(NewMap.G(2,x,y,z),2,x,y,z);
  }
  
}





//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialMapping' which is the partial mapping of 'MappingAtRefTimeStep' from the time 
//subdivision 'refTimeStep' due to the contribution of 'PartialVeloField'. Note, that an Identity mapping 'MappingId' is //also defined in the inputs (to avoid defining it each time the function is used)
void CptPartialMappingFromVeloFields(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * MappingId,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps,float DeltaX){
	int NBX,NBY,NBZ,NBT;
	int i,x,y,z;
	float x1,y1,z1;
	float x2,y2,z2;
	float x3,y3,z3;
	float x4,y4,z4;
	float xS,yS,zS;
	int t;
	float DeltaT_div_DeltaX;
	VectorField TotalBmap;  //total Backward mapping from the [new time_sub-1] to 0
	VectorField TotalBmap2;  //total Backward mapping from the [new time_sub] to 0
	VectorField TotalFmap;  //total forward mapping from the [new time_sub-1] to 0
	VectorField TotalFmap2;  //total forwrd mapping from the [new time_sub] to 0
	
	//1) INITIALISATION
	//1.1) constants
	NBX=VeloField->NX;
	NBY=VeloField->NY;
	NBZ=VeloField->NZ;
	NBT=VeloField->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
	//1.2) allocate memory in GradScalVF if not done
	if ((NBX!=PartialMapping->NX)||(NBY!=PartialMapping->NY)||(NBZ!=PartialMapping->NZ)||(NBT!=PartialMapping->NT))
		PartialMapping->CreateVoidField(NBX,NBY,NBZ,NBT);
	
	//1.3) allocate memory for TotalBmapand TotalFmap
	TotalBmap.CreateVoidField(NBX,NBY,NBZ,NBT);
	TotalBmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
	TotalFmap.CreateVoidField(NBX,NBY,NBZ,NBT);
	TotalFmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
	
	
	//1.4) PartialMapping at the reference time subdivision
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
		PartialMapping->P(MappingAtRefTimeStep->G(0,x,y,z),0,x,y,z,refTimeStep);
		PartialMapping->P(MappingAtRefTimeStep->G(1,x,y,z),1,x,y,z,refTimeStep);
		PartialMapping->P(MappingAtRefTimeStep->G(2,x,y,z),2,x,y,z,refTimeStep);
	}
	
	//2) PARTIAL FORWARD MAPPING FOR t>refTimeStep
	for (t=refTimeStep+1;t<NBT;t++){
		
		//2.1) compute the total backward mapping from t-1 to 0
		CptMappingFromVeloField(t-1,MappingId,VeloField,&TotalBmap,ConvergenceSteps);
		
		//2.2) compute the total backward mapping from t to 0
		CptMappingFromVeloField(t,MappingId,VeloField,&TotalBmap2,ConvergenceSteps);
		
		//2.3) compute partial forward map at t from the one at t-1
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
			//2.3.1) first estimation
			//2.3.1.a) first guess of where the information comes from at t-1
			x1=PartialMapping->G(0,x,y,z,t-1); y1=PartialMapping->G(1,x,y,z,t-1); z1=PartialMapping->G(2,x,y,z,t-1);
			x2=TotalBmap.G(0,x1,y1,z1,0);   y2=TotalBmap.G(1,x1,y1,z1,0);   z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
			
			xS=x-PartialVeloField->G(0,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
			yS=y-PartialVeloField->G(1,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
			zS=z-PartialVeloField->G(2,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
			
			//2.3.1.b) first transport of the information
			PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
			PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
			PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
			
			
			//2.3.1) leap frog style improvement of the estimation
			for (i=0;i<ConvergenceSteps*2;i++){
				//2.3.2.a) where the information comes from at t-1
				x1=PartialMapping->G(0,xS,yS,zS,t-1); y1=PartialMapping->G(1,xS,yS,zS,t-1); z1=PartialMapping->G(2,xS,yS,zS,t-1);
				x2=TotalBmap.G(0,x1,y1,z1,0);      y2=TotalBmap.G(1,x1,y1,z1,0);      z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
				
				x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
				x4=TotalBmap2.G(0,x3,y3,z3,0);  y4=TotalBmap2.G(1,x3,y3,z3,0);  z4=TotalBmap2.G(2,x3,y3,z3,0); //TotalBmap2 -> t
				
				xS=x-(PartialVeloField->G(0,x2,y2,z2,t-1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				yS=y-(PartialVeloField->G(1,x2,y2,z2,t-1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				zS=z-(PartialVeloField->G(2,x2,y2,z2,t-1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				
				//2.3.1.b) update the transport of the information
				PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
				PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
				PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
			}
		}
	}
	
	//3) PARTIAL BACWARD MAPPING FOR t<refTimeStep
	for (t=refTimeStep-1;t>=0;t--){      //not tested
		//3.1) compute the total forward mapping from 0 to t+1
		CptMappingFromVeloField(t+1,MappingId,VeloField,&TotalFmap,ConvergenceSteps);
		
		//3.2) compute the total forward mapping from 0 to t
		CptMappingFromVeloField(t,MappingId,VeloField,&TotalFmap2,ConvergenceSteps);
		
		//3.3) compute partial forward map at t from the one at t+1
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
			//3.3.1) first estimation
			//3.3.1.a) first guess of where the information comes from at t-1
			
			x1=PartialMapping->G(0,x,y,z,t+1); y1=PartialMapping->G(1,x,y,z,t+1); z1=PartialMapping->G(2,x,y,z,t+1);
			x2=TotalFmap.G(0,x1,y1,z1,NBT-1);   y2=TotalFmap.G(1,x1,y1,z1,NBT-1);   z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t+1
			
			xS=x+PartialVeloField->G(0,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
			yS=y+PartialVeloField->G(1,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
			zS=z+PartialVeloField->G(2,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
			
			//3.3.1.b) first transport of the information
			PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
			PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
			PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
			
			//3.3.1) leap frog style improvement of the estimation
			for (i=0;i<ConvergenceSteps*2;i++){
				//3.3.2.a) where the information comes from at t-1
				x1=PartialMapping->G(0,xS,yS,zS,t+1); y1=PartialMapping->G(1,xS,yS,zS,t+1); z1=PartialMapping->G(2,xS,yS,zS,t+1);
				x2=TotalFmap.G(0,x1,y1,z1,NBT-1);      y2=TotalFmap.G(1,x1,y1,z1,NBT-1);      z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t-1
				
				x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
				x4=TotalFmap2.G(0,x3,y3,z3,NBT-1);  y4=TotalFmap2.G(1,x3,y3,z3,NBT-1);  z4=TotalFmap2.G(2,x3,y3,z3,NBT-1); //TotalFmap2 -> t
				
				xS=x+(PartialVeloField->G(0,x2,y2,z2,t+1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				yS=y+(PartialVeloField->G(1,x2,y2,z2,t+1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				zS=z+(PartialVeloField->G(2,x2,y2,z2,t+1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				
				//3.3.1.b) update the transport of the information
				PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
				PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
				PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
			}
		}
	}
}



//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialMapping' which is the partial mapping of 'MappingAtRefTimeStep' from the time 
//subdivision 'refTimeStep' due to the contribution of 'PartialVeloField'.
void CptPartialMappingFromVeloFields_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps,float DeltaX){
	int NBX,NBY,NBZ,NBT;
	int i,x,y,z;
	float x1,y1,z1;
	float x2,y2,z2;
	float x3,y3,z3;
	float x4,y4,z4;
	float xS,yS,zS;
	int t;
	float DeltaT_div_DeltaX;
	VectorField TotalBmap;  //total Backward mapping from the [new time_sub-1] to 0
	VectorField TotalBmap2;  //total Backward mapping from the [new time_sub] to 0
	VectorField TotalFmap;  //total forward mapping from the [new time_sub-1] to 0
	VectorField TotalFmap2;  //total forwrd mapping from the [new time_sub] to 0
	
	//1) INITIALISATION
	//1.1) constants
	NBX=VeloField->NX;
	NBY=VeloField->NY;
	NBZ=VeloField->NZ;
	NBT=VeloField->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
	//1.2) allocate memory in GradScalVF if not done
	if ((NBX!=PartialMapping->NX)||(NBY!=PartialMapping->NY)||(NBZ!=PartialMapping->NZ)||(NBT!=PartialMapping->NT))
		PartialMapping->CreateVoidField(NBX,NBY,NBZ,NBT);
	
	//1.3) allocate memory for TotalBmapand TotalFmap
	TotalBmap.CreateVoidField(NBX,NBY,NBZ,NBT);
	TotalBmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
	TotalFmap.CreateVoidField(NBX,NBY,NBZ,NBT);
	TotalFmap2.CreateVoidField(NBX,NBY,NBZ,NBT);
	
	
	//1.4) PartialMapping at the reference time subdivision
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
		PartialMapping->P(static_cast<float>(x),0,x,y,z,refTimeStep);
		PartialMapping->P(static_cast<float>(y),1,x,y,z,refTimeStep);
		PartialMapping->P(static_cast<float>(z),2,x,y,z,refTimeStep);
	}
	
	//2) PARTIAL FORWARD MAPPING FOR t>refTimeStep
	for (t=refTimeStep+1;t<NBT;t++){
		
		//2.1) compute the total backward mapping from t-1 to 0
		CptMappingFromVeloField_IniIdMap(t-1,VeloField,&TotalBmap,ConvergenceSteps);
		
		//2.2) compute the total backward mapping from t to 0
		CptMappingFromVeloField_IniIdMap(t,VeloField,&TotalBmap2,ConvergenceSteps);
		
		//2.3) compute partial forward map at t from the one at t-1
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
			//2.3.1) first estimation
			//2.3.1.a) first guess of where the information comes from at t-1
			x1=PartialMapping->G(0,x,y,z,t-1); y1=PartialMapping->G(1,x,y,z,t-1); z1=PartialMapping->G(2,x,y,z,t-1);
			x2=TotalBmap.G(0,x1,y1,z1,0);   y2=TotalBmap.G(1,x1,y1,z1,0);   z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
			
			xS=x-PartialVeloField->G(0,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
			yS=y-PartialVeloField->G(1,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
			zS=z-PartialVeloField->G(2,x2,y2,z2,t-1)*DeltaT_div_DeltaX;
			
			//2.3.1.b) first transport of the information
			PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
			PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
			PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
			
			
			//2.3.1) leap frog style improvement of the estimation
			for (i=0;i<ConvergenceSteps*2;i++){
				//2.3.2.a) where the information comes from at t-1
				x1=PartialMapping->G(0,xS,yS,zS,t-1); y1=PartialMapping->G(1,xS,yS,zS,t-1); z1=PartialMapping->G(2,xS,yS,zS,t-1);
				x2=TotalBmap.G(0,x1,y1,z1,0);      y2=TotalBmap.G(1,x1,y1,z1,0);      z2=TotalBmap.G(2,x1,y1,z1,0); //TotalBmap -> t-1
				
				x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
				x4=TotalBmap2.G(0,x3,y3,z3,0);  y4=TotalBmap2.G(1,x3,y3,z3,0);  z4=TotalBmap2.G(2,x3,y3,z3,0); //TotalBmap2 -> t
				
				xS=x-(PartialVeloField->G(0,x2,y2,z2,t-1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				yS=y-(PartialVeloField->G(1,x2,y2,z2,t-1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				zS=z-(PartialVeloField->G(2,x2,y2,z2,t-1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				
				//2.3.1.b) update the transport of the information
				PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t-1),0,x,y,z,t);
				PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t-1),1,x,y,z,t);
				PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t-1),2,x,y,z,t);
			}
		}
	}
	
	//3) PARTIAL BACWARD MAPPING FOR t<refTimeStep
	for (t=refTimeStep-1;t>=0;t--){      //not tested
		//3.1) compute the total forward mapping from 0 to t+1
		CptMappingFromVeloField_IniIdMap(t+1,VeloField,&TotalFmap,ConvergenceSteps);
		
		//3.2) compute the total forward mapping from 0 to t
		CptMappingFromVeloField_IniIdMap(t,VeloField,&TotalFmap2,ConvergenceSteps);
		
		//3.3) compute partial forward map at t from the one at t+1
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
			//3.3.1) first estimation
			//3.3.1.a) first guess of where the information comes from at t-1
			
			x1=PartialMapping->G(0,x,y,z,t+1); y1=PartialMapping->G(1,x,y,z,t+1); z1=PartialMapping->G(2,x,y,z,t+1);
			x2=TotalFmap.G(0,x1,y1,z1,NBT-1);   y2=TotalFmap.G(1,x1,y1,z1,NBT-1);   z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t+1
			
			xS=x+PartialVeloField->G(0,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
			yS=y+PartialVeloField->G(1,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
			zS=z+PartialVeloField->G(2,x2,y2,z2,t+1)*DeltaT_div_DeltaX;
			
			//3.3.1.b) first transport of the information
			PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
			PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
			PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
			
			//3.3.1) leap frog style improvement of the estimation
			for (i=0;i<ConvergenceSteps*2;i++){
				//3.3.2.a) where the information comes from at t-1
				x1=PartialMapping->G(0,xS,yS,zS,t+1); y1=PartialMapping->G(1,xS,yS,zS,t+1); z1=PartialMapping->G(2,xS,yS,zS,t+1);
				x2=TotalFmap.G(0,x1,y1,z1,NBT-1);      y2=TotalFmap.G(1,x1,y1,z1,NBT-1);      z2=TotalFmap.G(2,x1,y1,z1,NBT-1); //TotalFmap -> t-1
				
				x3=PartialMapping->G(0,x,y,z,t);   y3=PartialMapping->G(1,x,y,z,t);   z3=PartialMapping->G(2,x,y,z,t);
				x4=TotalFmap2.G(0,x3,y3,z3,NBT-1);  y4=TotalFmap2.G(1,x3,y3,z3,NBT-1);  z4=TotalFmap2.G(2,x3,y3,z3,NBT-1); //TotalFmap2 -> t
				
				xS=x+(PartialVeloField->G(0,x2,y2,z2,t+1)+PartialVeloField->G(0,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				yS=y+(PartialVeloField->G(1,x2,y2,z2,t+1)+PartialVeloField->G(1,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				zS=z+(PartialVeloField->G(2,x2,y2,z2,t+1)+PartialVeloField->G(2,x4,y4,z4,t))*DeltaT_div_DeltaX/2.;
				
				//3.3.1.b) update the transport of the information
				PartialMapping->P(PartialMapping->G(0,xS,yS,zS,t+1),0,x,y,z,t);
				PartialMapping->P(PartialMapping->G(1,xS,yS,zS,t+1),1,x,y,z,t);
				PartialMapping->P(PartialMapping->G(2,xS,yS,zS,t+1),2,x,y,z,t);
			}
		}
	}
}




//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialLocMap' which is the partial mapping ONLY AT 'TargetSubdiv' FROM 'SourceSubdiv' due to the contribution of PartialVeloField.
//-> PartialLocMap therefore represents where are the coordinates of the points of time subdivision 'SourceSubdiv' when transported on time subdivision 'TargetSubdiv'. Note that we consider here an identity mapping at 'SourceSubdiv'.
void ComputeLagrangianPartialMapping(int SourceSubdiv,int TargetSubdiv,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialLocMap,float DeltaX){
	int NBX,NBY,NBZ,NBT;
	int x,y,z,t;
	float x1,y1,z1;
	float x2,y2,z2;
	float x3,y3,z3;
	float DeltaT_div_DeltaX;
	float DX,DY,DZ;
	float DX2,DY2,DZ2;
	
	//1) initialisation
	//1.1) constants
	NBX=VeloField->NX;
	NBY=VeloField->NY;
	NBZ=VeloField->NZ;
	NBT=VeloField->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
	//1.2) allocate memory in GradScalVF if not done
	if ((NBX!=PartialLocMap->NX)||(NBY!=PartialLocMap->NY)||(NBZ!=PartialLocMap->NZ)||(NBT!=1))
		PartialLocMap->CreateVoidField(NBX,NBY,NBZ,1);
	
	//2) Compute the transportation of the points of time SourceSubdiv to TargetSubdiv
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
		
		//initial coordinates
		x1=x*1.;  y1=y*1.;  z1=z*1.;  //for the transportation in the complete velocity field
		x2=x*1.;  y2=y*1.;  z2=z*1.;  //for the transportation in the partial velocity field
		
		//transportation...
		//...forward
		if (SourceSubdiv<TargetSubdiv) for (t=SourceSubdiv;t<TargetSubdiv;t++){
			x3=x1; y3=y1; z3=z1;
			
			DX=VeloField->G(0,x3,y3,z3,t)*DeltaT_div_DeltaX;
			DY=VeloField->G(1,x3,y3,z3,t)*DeltaT_div_DeltaX;
			DZ=VeloField->G(2,x3,y3,z3,t)*DeltaT_div_DeltaX;
			
			DX2=(VeloField->G(0,x3,y3,z3,t)+VeloField->G(0,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
			DY2=(VeloField->G(1,x3,y3,z3,t)+VeloField->G(1,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
			DZ2=(VeloField->G(2,x3,y3,z3,t)+VeloField->G(2,x3+DX,y3+DY,z3+DZ,t+1))*DeltaT_div_DeltaX/2.;
			DX=DX2;
			DY=DY2;
			DZ=DZ2;
			
			x1=x1+DX;
			y1=y1+DY;
			z1=z1+DZ;
			
			x2=x2+(PartialVeloField->G(0,x3,y3,z3,t)+PartialVeloField->G(0,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
			y2=y2+(PartialVeloField->G(1,x3,y3,z3,t)+PartialVeloField->G(1,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
			z2=z2+(PartialVeloField->G(2,x3,y3,z3,t)+PartialVeloField->G(2,x1,y1,z1,t+1))*DeltaT_div_DeltaX/2.;
		}
		
		//...backward
		if (SourceSubdiv<TargetSubdiv) for (t=SourceSubdiv;t>TargetSubdiv;t--){ //not tested
			x3=x1; y3=y1; z3=z1;
			
			DX=VeloField->G(0,x3,y3,z3,t)*DeltaT_div_DeltaX;
			DY=VeloField->G(1,x3,y3,z3,t)*DeltaT_div_DeltaX;
			DZ=VeloField->G(2,x3,y3,z3,t)*DeltaT_div_DeltaX;
			
			DX2=(VeloField->G(0,x3,y3,z3,t)+VeloField->G(0,x3-DX,y3-DY,z3-DZ,t-1))*DeltaT_div_DeltaX/2.;
			DY2=(VeloField->G(1,x3,y3,z3,t)+VeloField->G(1,x3-DX,y3-DY,z3-DZ,t-1))*DeltaT_div_DeltaX/2.;
			DZ2=(VeloField->G(2,x3,y3,z3,t)+VeloField->G(2,x3-DX,y3-DY,z3-DZ,t-1))*DeltaT_div_DeltaX/2.;
			DX=DX2;
			DY=DY2;
			DZ=DZ2;
			
			x1=x1-DX;
			y1=y1-DY;
			z1=z1-DZ;
			
			x2=x2-(PartialVeloField->G(0,x3,y3,z3,t)+PartialVeloField->G(0,x1,y1,z1,t-1))*DeltaT_div_DeltaX/2.;
			y2=y2-(PartialVeloField->G(1,x3,y3,z3,t)+PartialVeloField->G(1,x1,y1,z1,t-1))*DeltaT_div_DeltaX/2.;
			z2=z2-(PartialVeloField->G(2,x3,y3,z3,t)+PartialVeloField->G(2,x1,y1,z1,t-1))*DeltaT_div_DeltaX/2.;
		}
		
		
		//save where the point x,y,z is transported
		PartialLocMap->P(x2,0,x,y,z);
		PartialLocMap->P(y2,1,x,y,z);
		PartialLocMap->P(z2,2,x,y,z);
	}
}


#ifdef COMPILE_WITH_OPENMP

//Compute the projection of a 3D image 'ImagToPropag' using Mapping 'Map'.
//The image is projected at the time step 'TimeStepProj' of 'Map' and stored in 'ImageTimeT'.
//
//Importantly, the Mapping 'Map' should be an identity transformation at the time step 't' where 'ImagToPropag' is.
//It should also represent a forward mapping after 't' and a backward mapping before 't'.
void Project3Dimage(ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
  //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,t) 
  {
    #pragma omp for
    for (y = 0; y < ImageTimeT->NY; y++){
      for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++)  for (x = 0; x < ImageTimeT->NX; x++){
        ImageTimeT->P(ImagToPropag->G(Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
      }
    }
  //END FORK FOR THREADS
  }
}

#else

//Compute the projection of a 3D image 'ImagToPropag' using Mapping 'Map'.
//The image is projected at the time step 'TimeStepProj' of 'Map' and stored in 'ImageTimeT'.
//
//Importantly, the Mapping 'Map' should be an identity transformation at the time step 't' where 'ImagToPropag' is.
//It should also represent a forward mapping after 't' and a backward mapping before 't'.
void Project3Dimage(ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
  for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    ImageTimeT->P(ImagToPropag->G(Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
  }
}

#endif


#ifdef COMPILE_WITH_OPENMP

//same as above but the coordinates are transformed with the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfoAndTimeDepVF(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
 //BEGIN FORK FOR THREADS
  #pragma omp parallel default(shared) private(x,y,z,t) 
  {
    #pragma omp for
    for (y = 0; y < ImageTimeT->NY; y++){
      for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++)  for (x = 0; x < ImageTimeT->NX; x++){
        ImageTimeT->P(ImagToPropag->G(ProjectCS_2_OriginCS,Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
      }
    }
  //END FORK FOR THREADS
  }
}

#else

//same as above but the coordinates are transformed with the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfoAndTimeDepVF(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
  int x,y,z,t;
  
  for (t = 0; t < ImageTimeT->NT; t++) for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    ImageTimeT->P(ImagToPropag->G(ProjectCS_2_OriginCS,Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj),t),x,y,z,t);
  }
}

#endif


//same as above but the coordinates are transformed with a displacement field in voxels (from 'Map' c. s. to 'ImagToPropag' c. s.)
void Project3DImageUsingDispFieldAndTimeDepVF(VectorField * DispField,ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj){
	int x,y,z;
	float x2,y2,z2;
	
	for (z = 0; z < ImageTimeT->NZ; z++) for (y = 0; y < ImageTimeT->NY; y++) for (x = 0; x < ImageTimeT->NX; x++){
    x2=Map->G(0,x,y,z,TimeStepProj);
    y2=Map->G(1,x,y,z,TimeStepProj);
    z2=Map->G(2,x,y,z,TimeStepProj);
    
		ImageTimeT->P(ImagToPropag->G(DispField->G(0,x2,y2,z2),DispField->G(1,x2,y2,z2),DispField->G(2,x2,y2,z2)),x,y,z);
	}
}






//same as above but with vector fields (NOT REORIENTED!!!)
void Project3Dimage(VectorField * VFToPropag,VectorField * Map,VectorField * VFTimeT,int TimeStepProj){
	int x,y,z,i;
	
  for (i = 0; i < 3; i++) for (z = 0; z < VFTimeT->NZ; z++) for (y = 0; y < VFTimeT->NY; y++) for (x = 0; x < VFTimeT->NX; x++){
		VFTimeT->P(VFToPropag->G(i,Map->G(0,x,y,z,TimeStepProj), Map->G(1,x,y,z,TimeStepProj),Map->G(2,x,y,z,TimeStepProj)),i,x,y,z);
	}
}





///By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
///'VeloField4Measure' in the length of the flow from each point of the field. The length of flow
///is projected AT T=0 and returned in the 3D scalar field 'LengthOfFlow'
/// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
/// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
///   is computed.
void CptLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps,float DeltaX){
	float VecTemp[3];
	float VecTemp2[3];
	float NormVecTemp;
	int NBX,NBY,NBZ,NBT;
	int i,x,y,z;
	int t;
	float DeltaT_div_DeltaX;
	ScalarField PrevLengthOfFlow;
	
	//1) INITIALISATION
	//1.1) field size
	NBX=VeloField4Flow->NX;
	NBY=VeloField4Flow->NY;
	NBZ=VeloField4Flow->NZ;
	NBT=VeloField4Flow->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
	//1.2) allocate memory in LengthOfFlow at time t if not done
	if ((NBX!=LengthOfFlow->NX)||(NBY!=LengthOfFlow->NY)||(NBZ!=LengthOfFlow->NZ))
		LengthOfFlow->CreateVoidField(NBX,NBY,NBZ);
	
	//1.3) allocate memory of the PrevLengthOfFlow at time t+1
	PrevLengthOfFlow.CreateVoidField(NBX,NBY,NBZ);
	
	//1.4) PrevLengthOfFlow at the last time subdivision
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) PrevLengthOfFlow.P(0.,x,y,z);
	
	//JO at the other time subdivisions
	for (t=NBT-2;t>=0;t--){
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
			//init
			VecTemp[0]=0.; 
			VecTemp[1]=0.;
			VecTemp[2]=0.;
			
			//convergence
			for (i=0;i<ConvergenceSteps;i++){
				VecTemp2[0]=VeloField4Flow->G(0,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
				VecTemp2[1]=VeloField4Flow->G(1,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
				VecTemp2[2]=VeloField4Flow->G(2,x+VecTemp[0],y+VecTemp[1],z+VecTemp[2],t+1);
				
				VecTemp[0]=(VecTemp2[0]+VeloField4Flow->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[1]=(VecTemp2[1]+VeloField4Flow->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[2]=(VecTemp2[2]+VeloField4Flow->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
			}
			
			//compute the lenght
			NormVecTemp =(float)pow((double)VeloField4Measure->G(0,x,y,z,t)*DeltaT_div_DeltaX,2.);
			NormVecTemp+=(float)pow((double)VeloField4Measure->G(1,x,y,z,t)*DeltaT_div_DeltaX,2.);
			NormVecTemp+=(float)pow((double)VeloField4Measure->G(2,x,y,z,t)*DeltaT_div_DeltaX,2.);
			NormVecTemp=sqrt(NormVecTemp);
			
			LengthOfFlow->P(NormVecTemp+PrevLengthOfFlow.G(x+VecTemp[0],y+VecTemp[1],z+VecTemp[2]),x,y,z);
		}
		for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
			PrevLengthOfFlow.P(LengthOfFlow->G(x,y,z),x,y,z);
	}
}











///By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
///'VeloField4Measure' AT THE CURRENT TIME in the length of the flow from each point of the field. The length of flow
///is returned in the 3D+t scalar field 'LengthOfFlow'
/// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
/// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
///   is computed.
void CptEvoLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps,float DeltaX){
	float VecTemp[3];
	float VecTemp2[3];
	float VecTemp3[3];
	float TmpFl;
	int NBX,NBY,NBZ,NBT;
	int i,x,y,z;
	int t;
	float DeltaT_div_DeltaX;
	
	//1) INITIALISATION
	//1.1) field size
	NBX=VeloField4Flow->NX;
	NBY=VeloField4Flow->NY;
	NBZ=VeloField4Flow->NZ;
	NBT=VeloField4Flow->NT;
	DeltaT_div_DeltaX=1./((NBT-1.)*DeltaX);
	
	//1.2) allocate memory in LengthOfFlow at time t if not done
	if ((NBX!=LengthOfFlow->NX)||(NBY!=LengthOfFlow->NY)||(NBZ!=LengthOfFlow->NZ)||(NBT!=LengthOfFlow->NT))
		LengthOfFlow->CreateVoidField(NBX,NBY,NBZ,NBT);
	
	//1.3) JO at the first time subdivision
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
		LengthOfFlow->P(0.,x,y,z,0);
	
	//   ScalarField toto;
	//   toto.Read("AOD_Deformation1to2.nii");
	//   for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++)
	//         LengthOfFlow->P(toto.G(x,y,z,14),x,y,z,0);
	
	
	//2) Computation of the amplitude of the deformations
	for (t=1;t<NBT;t++){
		if (ConvergenceSteps<=1){ // simple integration scheme (centered in time)
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				VecTemp[0]=(VeloField4Flow->G(0,x,y,z,t-1)+VeloField4Flow->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[1]=(VeloField4Flow->G(1,x,y,z,t-1)+VeloField4Flow->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp[2]=(VeloField4Flow->G(2,x,y,z,t-1)+VeloField4Flow->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				
				//propagate the lenght of flow
				VecTemp3[0]=(VeloField4Measure->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp3[1]=(VeloField4Measure->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp3[2]=(VeloField4Measure->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				
				TmpFl =(float)pow(VecTemp3[0],2);
				TmpFl+=(float)pow(VecTemp3[1],2);
				TmpFl+=(float)pow(VecTemp3[2],2);
				TmpFl=sqrt(TmpFl);
				TmpFl+=LengthOfFlow->G(x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
				LengthOfFlow->P(TmpFl,x,y,z,t);
			}
		}
		else{ // leap frog scheme
			for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) for (x=0;x<NBX;x++){
				//init
				VecTemp[0]=0.; 
				VecTemp[1]=0.;
				VecTemp[2]=0.;
				
				//convergence
				for (i=0;i<ConvergenceSteps;i++){
					VecTemp2[0]=VeloField4Flow->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					VecTemp2[1]=VeloField4Flow->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					VecTemp2[2]=VeloField4Flow->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
					
					VecTemp[0]=(VecTemp2[0]+VeloField4Flow->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[1]=(VecTemp2[1]+VeloField4Flow->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
					VecTemp[2]=(VecTemp2[2]+VeloField4Flow->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				}
				
				//propagate the lenght of flow
				VecTemp3[0]=(VeloField4Measure->G(0,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(0,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp3[1]=(VeloField4Measure->G(1,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(1,x,y,z,t))*DeltaT_div_DeltaX/2;
				VecTemp3[2]=(VeloField4Measure->G(2,x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1)+VeloField4Measure->G(2,x,y,z,t))*DeltaT_div_DeltaX/2;
				
				
				TmpFl =(float)pow(VecTemp3[0],2);
				TmpFl+=(float)pow(VecTemp3[1],2);
				TmpFl+=(float)pow(VecTemp3[2],2);
				TmpFl=sqrt(TmpFl);
				TmpFl+=LengthOfFlow->G(x-VecTemp[0],y-VecTemp[1],z-VecTemp[2],t-1);
				LengthOfFlow->P(TmpFl,x,y,z,t);
			}
		}
	}
}


///compute the L_2 norm of the difference between two scalar fields
float CalcSqrtSumOfSquaredDif(ScalarField * I1,ScalarField * I2){
	int x,y,z;
	float L2_norm,tmp;
	
	L2_norm=0.;
	
	for (z=0;z<I1->NZ;z++) for (y=0;y<I1->NY;y++) for (x=0;x<I1->NX;x++){
		tmp=(I1->G(x,y,z)-I2->G(x,y,z));
		L2_norm+=tmp*tmp;
	}
	
	L2_norm=sqrt(L2_norm);
	
	return L2_norm;
}




///...
void TransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,int t)
{	int x,y,z;
	int NBX,NBY,NBZ;
	float d11,d12,d13,d21,d22,d23,d31,d32,d33;
	float temp;
	NBX=TempInvDiffeo->NX;
	NBY=TempInvDiffeo->NY;
	NBZ=TempInvDiffeo->NZ;
	for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
		d11=(TempInvDiffeo->G(0,x+1,y,z,t)-TempInvDiffeo->G(0,x-1,y,z,t))/(2.*DeltaX);
		d12=(TempInvDiffeo->G(0,x,y+1,z,t)-TempInvDiffeo->G(0,x,y-1,z,t))/(2.*DeltaX);
		d13=(TempInvDiffeo->G(0,x,y,z+1,t)-TempInvDiffeo->G(0,x,y,z-1,t))/(2.*DeltaX);
		d21=(TempInvDiffeo->G(1,x+1,y,z,t)-TempInvDiffeo->G(1,x-1,y,z,t))/(2.*DeltaX);
		d22=(TempInvDiffeo->G(1,x,y+1,z,t)-TempInvDiffeo->G(1,x,y-1,z,t))/(2.*DeltaX);
		d23=(TempInvDiffeo->G(1,x,y,z+1,t)-TempInvDiffeo->G(1,x,y,z-1,t))/(2.*DeltaX);
		d31=(TempInvDiffeo->G(2,x+1,y,z,t)-TempInvDiffeo->G(2,x-1,y,z,t))/(2.*DeltaX);
		d32=(TempInvDiffeo->G(2,x,y+1,z,t)-TempInvDiffeo->G(2,x,y-1,z,t))/(2.*DeltaX);
		d33=(TempInvDiffeo->G(2,x,y,z+1,t)-TempInvDiffeo->G(2,x,y,z-1,t))/(2.*DeltaX);
		temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0);
		Momentum->P( temp* (d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13)),x,y,z);
	}
	//1.2.2) boundaries at 0.
	z=0;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	z=NBZ-1;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	y=0;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	y=NBY-1;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	x=0;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	x=NBX-1;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->P(InitialMomentum->G(x,y,z,0),x,y,z);
	
	//1.2.3) 2D image case
	//float max=0.0;
	if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
		d11=(TempInvDiffeo->G(0,x+1,y,0,t)-TempInvDiffeo->G(0,x-1,y,0,t))/(2.*DeltaX);
		d12=(TempInvDiffeo->G(0,x,y+1,0,t)-TempInvDiffeo->G(0,x,y-1,0,t))/(2.*DeltaX);
		d21=(TempInvDiffeo->G(1,x+1,y,0,t)-TempInvDiffeo->G(1,x-1,y,0,t))/(2.*DeltaX);
		d22=(TempInvDiffeo->G(1,x,y+1,0,t)-TempInvDiffeo->G(1,x,y-1,0,t))/(2.*DeltaX);
		temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,0,t),TempInvDiffeo->G(1,x,y,0,t),TempInvDiffeo->G(2,x,y,0,t),0);
		//if (max<abs(temp*(d11*d22-d21*d12))){max=abs(temp*(d11*d22-d21*d12));}
		Momentum->P(temp*(d11*d22-d21*d12),x,y,0);
	}
	//cout << "c'est penible  "<< Momentum->GetMaxAbsVal() <<"\n";
}

/// Bug - Not VaLidated !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
void AddTransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,float cste, int t)
{	int x,y,z;
	int NBX,NBY,NBZ;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33;
	float temp;
	NBX=TempInvDiffeo->NX;
	NBY=TempInvDiffeo->NY;
	NBZ=TempInvDiffeo->NZ;
	for (z=1;z<NBZ-1;z++) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
		d11=(TempInvDiffeo->G(0,x+1,y,z,t)-TempInvDiffeo->G(0,x-1,y,z,t))/(2.*DeltaX);
		d12=(TempInvDiffeo->G(0,x,y+1,z,t)-TempInvDiffeo->G(0,x,y-1,z,t))/(2.*DeltaX);
		d13=(TempInvDiffeo->G(0,x,y,z+1,t)-TempInvDiffeo->G(0,x,y,z-1,t))/(2.*DeltaX);
		d21=(TempInvDiffeo->G(1,x+1,y,z,t)-TempInvDiffeo->G(1,x-1,y,z,t))/(2.*DeltaX);
		d22=(TempInvDiffeo->G(1,x,y+1,z,t)-TempInvDiffeo->G(1,x,y-1,z,t))/(2.*DeltaX);
		d23=(TempInvDiffeo->G(1,x,y,z+1,t)-TempInvDiffeo->G(1,x,y,z-1,t))/(2.*DeltaX);
		d31=(TempInvDiffeo->G(2,x+1,y,z,t)-TempInvDiffeo->G(2,x-1,y,z,t))/(2.*DeltaX);
		d32=(TempInvDiffeo->G(2,x,y+1,z,t)-TempInvDiffeo->G(2,x,y-1,z,t))/(2.*DeltaX);
		d33=(TempInvDiffeo->G(2,x,y,z+1,t)-TempInvDiffeo->G(2,x,y,z-1,t))/(2.*DeltaX);
		temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0);
		Momentum->Add(cste* temp* (d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13)),x,y,z);
	}
	//1.2.2) boundaries at 0.
	z=0;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	z=NBZ-1;
	for (y=0;y<NBY;y++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	y=0;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	y=NBY-1;
	for (z=0;z<NBZ;z++) for (x=0;x<NBX;x++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	x=0;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	x=NBX-1;
	for (z=0;z<NBZ;z++) for (y=0;y<NBY;y++) Momentum->Add(cste*InitialMomentum->G(x,y,z,0),x,y,z);
	
	//1.2.3) 2D image case
	if (NBZ==1) for (y=1;y<NBY-1;y++) for (x=1;x<NBX-1;x++)
	{
		d11=(TempInvDiffeo->G(0,x+1,y,0,t)-TempInvDiffeo->G(0,x-1,y,0,t))/(2.*DeltaX);
		d12=(TempInvDiffeo->G(0,x,y+1,0,t)-TempInvDiffeo->G(0,x,y-1,0,t))/(2.*DeltaX);
		d21=(TempInvDiffeo->G(1,x+1,y,0,t)-TempInvDiffeo->G(1,x-1,y,0,t))/(2.*DeltaX);
		d22=(TempInvDiffeo->G(1,x,y+1,0,t)-TempInvDiffeo->G(1,x,y-1,0,t))/(2.*DeltaX);
		temp = InitialMomentum->G(TempInvDiffeo->G(0,x,y,0,t),TempInvDiffeo->G(1,x,y,0,t));
		Momentum->Add(cste*temp*(d11*d22-d21*d12),x,y,0);
	}
}

///...
void TransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image, int t)
{
	int x,y,z;
  for (z = 0; z < InitialImage->NZ; z++) for (y = 0; y < InitialImage->NY; y++) for (x = 0; x < InitialImage->NX; x++)
	{
		Image->P(InitialImage->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0),x,y,z,0);
	}
}

///...
void AddTransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image, float cste, int t)
{
	int x,y,z;
  for (z = 0; z < InitialImage->NZ; z++) for (y = 0; y < InitialImage->NY; y++) for (x = 0; x < InitialImage->NX; x++)
	{
		Image->Add(cste * InitialImage->G(TempInvDiffeo->G(0,x,y,z,t),TempInvDiffeo->G(1,x,y,z,t),TempInvDiffeo->G(2,x,y,z,t),0),x,y,z,0);
	}
}

///...
void DeepCopy(VectorField *VectorField1,VectorField *VectorField2,int t)
{
	int i,x,y,z;
	for (i=0;i<3;i++)
	{
		for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
		{
			VectorField2->P(VectorField1->G(i,x,y,z),i,x,y,z,t);
		}
	}
}
void DeepCopy(VectorField *VectorField,ScalarField* ScalarField,int direc,int t)
{
	int i,x,y,z;
	for (i=0;i<3;i++)
	{
		for (z = 0; z < ScalarField->NZ; z++) for (y = 0; y < ScalarField->NY; y++) for (x = 0; x < ScalarField->NX; x++)
		{
			ScalarField->P(VectorField->G(direc,x,y,z,t),x,y,z,0);
		}
	}
}


///...
void DeepCopy(ScalarField *ScalarField1,ScalarField *ScalarField2,int t)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField2->P(ScalarField1->G(x,y,z),x,y,z,t);
	}	
}


/// Compute the scalar product between the vector fields and put it in ScalarField0 (for which NT=1)
void ScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t,float cste)
{
	int i,x,y,z;
	float temp;
	for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		temp=0.0;
		for (i=0;i<3;i++){temp+=VectorField1->G(i,x,y,z,t)*VectorField2->G(i,x,y,z,t);}
		ScalarField0->P(cste*temp,x,y,z);
	}
}


/// Compute the scalar product between the scalar fields and put it in ScalarField0 (for which NT=1)
void ScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t,float cste)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField0->P(cste * ScalarField1->G(x,y,z,t)*ScalarField2->G(x,y,z,t),x,y,z);
	}	
}


/// Add the scalar product between the vector fields to ScalarField0 (for which NT=1)
void AddScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t)
{
	int i,x,y,z;
	for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		for (i=0;i<3;i++){ScalarField0->Add(VectorField1->G(i,x,y,z,t)*VectorField2->G(i,x,y,z,t),x,y,z);}
	}	
}


/// Add the scalar product between the scalar fields to ScalarField0 (for which NT=1)
void AddScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField0->Add(ScalarField1->G(x,y,z,t)*ScalarField2->G(x,y,z,t),x,y,z);
	}	
}


/// Add  ScalarField1 at time t1 to ScalarField2 at time t2
void AddScalarField(ScalarField *ScalarField1, ScalarField *ScalarField2,float cste, int t1,int t2)
{
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		ScalarField2->Add(cste*ScalarField1->G(x,y,z,t1),x,y,z,t2);
	}	
}
/// Add  ScalarField1 at time t1 to ScalarField2 at time t2
void AddVectorField(VectorField *VectorField1, VectorField *VectorField2,float cste, int t1,int t2)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		VectorField2->Add(cste*VectorField1->G(i,x,y,z,t1),i,x,y,z,t2);
	}	
}
/// Multiply a vector field by the cste
void MultiplyVectorField(VectorField *VectorField1, float cste,int t)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		VectorField1->P(cste*VectorField1->G(x,y,z,t),i,x,y,z,t);
	}	
}

/// Sum two vector fields and put it in Output.
void SumVectorField(VectorField *VectorField1, VectorField *VectorField2, VectorField *Output, int t1, int t2, int t3, float cste1 ,float cste2)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < VectorField1->NZ; z++) for (y = 0; y < VectorField1->NY; y++) for (x = 0; x < VectorField1->NX; x++)
	{
		Output->P(cste1 * VectorField1->G(i,x,y,z,t1) + cste2 * VectorField2->G(i,x,y,z,t2),i,x,y,z,t3);
	}
}
/// compute the product element by element along each dimension
void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2)
{
	int i,x,y,z;
	for (i=0;i<3;i++) for (z = 0; z < ScalarField->NZ; z++) for (y = 0; y < ScalarField->NY; y++) for (x = 0; x < ScalarField->NX; x++)
	{
		VectorField2->P(VectorField1->G(i,x,y,z)*ScalarField->G(x,y,z),i,x,y,z);
	}
}


///...
float DotProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, int t1,int t2)
{
	float result=0.0;
	int x,y,z;
	for (z = 0; z < ScalarField1->NZ; z++) for (y = 0; y < ScalarField1->NY; y++) for (x = 0; x < ScalarField1->NX; x++)
	{
		result += ScalarField1->G(x,y,z,t1) * ScalarField2->G(x,y,z,t2);
	}
	return result;
}







///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                   7: OTHER FUNCTIONS OF SCIENTIFIC COMPUTATION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// Solve the problem: MX=D where D is a known vector, M a tridiagonal matrix and X the unknown vector.
/// Inputs are a,b,c,d,n where M(i,i)=b(i), M(i,i-1)=a(i), M(i,i+1)=c(i), D(i)=d(i), D in R^n and M in R^n*R^n.
/// Output is X where X in R^n.  Warning: will modify c and d! */
void TridiagonalSolveFloat(const float *a, const float *b, float *c, float *d, float *x, int n){
  int i;
  double id;
  
  /* Modify the coefficients. */
  c[0] /= b[0];                       /* Division by zero risk. */
  d[0] /= b[0];                       /* Division by zero would imply a singular matrix. */
  for(i = 1; i < n; i++){
    id = (b[i] - c[i-1] * a[i]);      /* Division by zero risk. */
    c[i] /= id;                       /* Last value calculated is redundant. */
    d[i] = (d[i] - d[i-1] * a[i])/id;
  }
  
  /* Now back substitute. */
  x[n - 1] = d[n - 1];
  for(i = n - 2; i >= 0; i--)
    x[i] = d[i] - c[i] * x[i + 1];
}


//Perform the eigenvalue decomposition of a 3*3 matrix
//Adapted from the algorithm having the same name in 'numerical recipes'.
//Input:  The 3*3 matrix 'MatIni' that has to be symetric.
//Ouputs: 'ValP' is a vector of size 3 which contains the eigen values (in decreasing order). 
//        'VecP' is a 3*3 matrix containg the eigen vectors in columns.
#define ROTATE(a,i,j,k,l) g=a[i][j];h=a[k][l];a[i][j]=g-s*(h+g*tau);a[k][l]=h+s*(g-h*tau);
void jacobi3(float **MatIni,float *ValP, float **VecP){
  int j,iq,ip,i; 
  float tresh,theta,tau,t,sm,s,h,g,c;
  float b[4];
  float z[4];
  float a[4][4];   //correspond a MatIni
  float d[4];    //correspond a ValP
  float v[4][4];   //correspond a VecP
  int vTri1,vTri2;
  float TempF;
  int n;
  
  
  n=3;
  for(i=0;i<n;i++) for(j=0;j<n;j++) a[i+1][j+1]=MatIni[i][j];
  
  
  //algo de numerical recipes
  for (ip=1;ip<=n;ip++) {
    for (iq=1;iq<=n;iq++) v[ip][iq]=0.0;
    v[ip][ip]=1.0;
  }
	
  for (ip=1;ip<=n;ip++) {
    b[ip]=d[ip]=a[ip][ip];
    z[ip]=0.0;
  }
	
  for (i=1;i<=50;i++) {
    sm=0.0;
    for (ip=1;ip<=n-1;ip++)
      for (iq=ip+1;iq<=n;iq++)
        sm += fabs(a[ip][iq]);
    
    if (sm == 0.0) {
      //adaptation des valeurs de l'algo de numerical recipes aux valeurs de sortie
      for(i=0;i<n;i++) ValP[i]=d[i+1];
      for(i=0;i<n;i++) for(j=0;j<n;j++) MatIni[i][j]=a[i+1][j+1];
      for(i=0;i<n;i++) for(j=0;j<n;j++) VecP[i][j]=v[i+1][j+1];
      
      //tri des donnees
      for(vTri1=0;vTri1<n-1;vTri1++) for(vTri2=vTri1+1;vTri2<n;vTri2++) if (ValP[vTri1]<ValP[vTri2]){
        TempF=ValP[vTri1]; ValP[vTri1]=ValP[vTri2]; ValP[vTri2]=TempF;
        for(i=0;i<n;i++) { TempF=VecP[i][vTri1]; VecP[i][vTri1]=VecP[i][vTri2]; VecP[i][vTri2]=TempF;}
      }
      
      return;
    }
    if (i < 4) tresh=0.2*sm/(n*n);
    else tresh=0.0;
    
    for (ip=1;ip<=n-1;ip++) {
      for (iq=ip+1;iq<=n;iq++) {
        g=100.0*fabs(a[ip][iq]);
        
        if (i > 4 && (float)(fabs(d[ip])+g) == (float)fabs(d[ip])&& (float)(fabs(d[iq])+g) == (float)fabs(d[iq]))
          a[ip][iq]=0.0;
        else if (fabs(a[ip][iq]) > tresh) {
          h=d[iq]-d[ip];
          if ((float)(fabs(h)+g) == (float)fabs(h))
            t=(a[ip][iq])/h;
          else {
            theta=0.5*h/(a[ip][iq]);
            t=1.0/(fabs(theta)+sqrt(1.0+theta*theta));
            if (theta < 0.0) t = -t;
          }
          c=1.0/sqrt(1+t*t); s=t*c; tau=s/(1.0+c); h=t*a[ip][iq];
          z[ip] -= h; z[iq] += h; d[ip] -= h; d[iq] += h; a[ip][iq]=0.0;
          for (j=1;j<=ip-1;j++) { ROTATE(a,j,ip,j,iq) }
          for (j=ip+1;j<=iq-1;j++) { ROTATE(a,ip,j,j,iq) }
          for (j=iq+1;j<=n;j++) { ROTATE(a,ip,j,iq,j) }
          for (j=1;j<=n;j++) { ROTATE(v,j,ip,j,iq)}
        }
      }
    }
    for (ip=1;ip<=n;ip++) { b[ip] += z[ip]; d[ip]=b[ip]; z[ip]=0.0; }
  }
  printf("Too many iterations in the routine jacobi\n");
  
  //adaptation des valeurs de l'algo de numerical recipes aux valeurs de sortie
  for(i=0;i<n;i++) ValP[i]=d[i+1];
  for(i=0;i<n;i++) for(j=0;j<n;j++) MatIni[i][j]=a[i+1][j+1];
  for(i=0;i<n;i++) for(j=0;j<n;j++) VecP[i][j]=v[i+1][j+1];
  
  //tri des donnees
  for(vTri1=0;vTri1<n-1;vTri1++) for(vTri2=vTri1+1;vTri2<n;vTri2++) if (ValP[vTri1]<ValP[vTri2]){
    TempF=ValP[vTri1]; ValP[vTri1]=ValP[vTri2]; ValP[vTri2]=TempF;
    for(i=0;i<n;i++) { TempF=VecP[i][vTri1]; VecP[i][vTri1]=VecP[i][vTri2]; VecP[i][vTri2]=TempF;}
  }
  
}



///compute two orthogonal vectors tvec1 and tvec2 in R^3 which are orthogonal to nvec
///the norm of tvec1 and tvec2 is defined as equal to the one of nvec
void CptVecsTangentPlane(float nvec[3],float tvec1[3],float tvec2[3]){
  float epsilon,dist;
  
  dist=sqrt(nvec[0]*nvec[0]+nvec[1]*nvec[1]+nvec[2]*nvec[2]);
  epsilon=dist/100;
  
  //define two orthogonal directions in the plan where transformations are allowed...
  //... vec1
  if (fabs(nvec[0])<epsilon){
    tvec1[0]=1; tvec1[1]=0; tvec1[2]=0;
  }
  else{
    tvec1[0]=0;
    
    if (fabs(nvec[1])<epsilon){
      tvec1[1]=1; tvec1[2]=0; 
    }
    else{
      tvec1[1]=-nvec[2]/nvec[1];
      tvec1[2]=1;
    }
  }
  
  VecNormalize(tvec1,dist);
  
  //... vec2
  tvec2[0]=(nvec[1]*tvec1[2])-(nvec[2]*tvec1[1]);
  tvec2[1]=-(nvec[0]*tvec1[2])+(nvec[2]*tvec1[0]);
  tvec2[2]=(nvec[0]*tvec1[1])-(nvec[1]*tvec1[0]);
  
  VecNormalize(tvec2,dist);
}



///normalize a vector
void VecNormalize(float vec[3],float norm){
  float tmpfl;
        
  tmpfl=sqrt(vec[0]*vec[0]+vec[1]*vec[1]+vec[2]*vec[2]);
  vec[0]*=norm/tmpfl;
  vec[1]*=norm/tmpfl;
  vec[2]*=norm/tmpfl;
}


///compute the determinant of a 3*3 matrix
float determinant_3t3matrix(float m[3][3]){
  float deter;
  deter=m[0][0]*(m[1][1]*m[2][2]-m[1][2]*m[2][1]);
  deter-=m[1][0]*(m[0][1]*m[2][2]-m[2][1]*m[0][2]);
  deter+=m[2][0]*(m[0][1]*m[1][2]-m[1][1]*m[0][2]);
  return deter;
}

///compute the comatrix of a 3*3 matrix
void comatrix_3t3matrix(float m1[3][3],float m2[3][3]){
  m2[0][0]=(m1[1][1]*m1[2][2])-(m1[2][1]*m1[1][2]);
  m2[1][0]=(m1[2][1]*m1[0][2])-(m1[0][1]*m1[2][2]);
  m2[2][0]=(m1[0][1]*m1[1][2])-(m1[1][1]*m1[0][2]);
  m2[0][1]=(m1[2][0]*m1[1][2])-(m1[1][0]*m1[2][2]);
  m2[1][1]=(m1[0][0]*m1[2][2])-(m1[2][0]*m1[0][2]);
  m2[2][1]=(m1[1][0]*m1[0][2])-(m1[0][0]*m1[1][2]);
  m2[0][2]=(m1[1][0]*m1[2][1])-(m1[2][0]*m1[1][1]);
  m2[1][2]=(m1[2][0]*m1[0][1])-(m1[0][0]*m1[2][1]);
  m2[2][2]=(m1[0][0]*m1[1][1])-(m1[1][0]*m1[0][1]);
  
}

///Estimate the exponential of a 3*3 matrix   (Checked with matlab and OK)
void Exponential_3t3matrix(float m1[3][3],float m2[3][3]){
  float tempMat[3][3];
  float tempMat2[3][3];
  int k,facto;
  
  //init  (k=0)
  m2[0][0]=1; m2[0][1]=0; m2[0][2]=0; 
  m2[1][0]=0; m2[1][1]=1; m2[1][2]=0; 
  m2[2][0]=0; m2[2][1]=0; m2[2][2]=1; 
  
  tempMat[0][0]=m1[0][0]; tempMat[0][1]=m1[0][1]; tempMat[0][2]=m1[0][2];
  tempMat[1][0]=m1[1][0]; tempMat[1][1]=m1[1][1]; tempMat[1][2]=m1[1][2];
  tempMat[2][0]=m1[2][0]; tempMat[2][1]=m1[2][1]; tempMat[2][2]=m1[2][2];
  
  //estimation...
  facto=1;
  for (k=1;k<15;k++){
    //...update m2
    m2[0][0]+=tempMat[0][0]/facto; m2[0][1]+=tempMat[0][1]/facto; m2[0][2]+=tempMat[0][2]/facto;
    m2[1][0]+=tempMat[1][0]/facto; m2[1][1]+=tempMat[1][1]/facto; m2[1][2]+=tempMat[1][2]/facto;
    m2[2][0]+=tempMat[2][0]/facto; m2[2][1]+=tempMat[2][1]/facto; m2[2][2]+=tempMat[2][2]/facto;
    
    //...update facto
    facto*=k+1;
    
    //...update tempMat
    mult_3t3mat_3t3mat(tempMat,m1,tempMat2);
    
    tempMat[0][0]=tempMat2[0][0]; tempMat[0][1]=tempMat2[0][1]; tempMat[0][2]=tempMat2[0][2];
    tempMat[1][0]=tempMat2[1][0]; tempMat[1][1]=tempMat2[1][1]; tempMat[1][2]=tempMat2[1][2];
    tempMat[2][0]=tempMat2[2][0]; tempMat[2][1]=tempMat2[2][1]; tempMat[2][2]=tempMat2[2][2];
  }
}


///transpose a 3*3 matrix
void transpose_3t3matrix(float m1[3][3],float m2[3][3]){
  int i,j;
  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      m2[i][j]=m1[j][i];
}


///multiply two 3*3 matrices
void  mult_3t3mat_3t3mat(float m1[3][3], float m2[3][3], float MatRes[3][3]){
  int i,j,k;
  
  for (i=0;i<3;i++)
    for(j=0;j<3;j++){
      MatRes[i][j]=0;
      for(k=0;k<3;k++) MatRes[i][j]+=m1[i][k]*m2[k][j];
    }
}



///inverse of a 3*3 matrix
void invert_3t3matrix(float m1[3][3],float m2[3][3]){
  float det;
  int i,j;
  
  det=determinant_3t3matrix(m1);
  comatrix_3t3matrix(m1,m2);
  transpose_3t3matrix(m2,m1);
  
  for (i=0;i<3;i++)
    for(j=0;j<3;j++)
      m2[i][j]=m1[i][j]/det;
  
}

///multiply a vector of size 3 and a 3*3 matrix
void mult_3t3mat_3vec(float mat[3][3], float vectIni[3], float vectRes[3]){
  int lig,col;
  
  for (lig=0;lig<3;lig++){
    vectRes[lig]=0;
    for(col=0;col<3;col++)
      vectRes[lig]+=mat[lig][col]*vectIni[col];
	}
}


///inverse of a quaternion
void invert_4t4quaternion(float q1[4][4],float q2[4][4]){
  float r11,r12,r13,r21,r22,r23,r31,r32,r33,v1,v2,v3 , deti ;
  
  //algorithm inspired from the one of nifti_io.c
  
  /*  INPUT MATRIX IS:  */
  r11 = q1[0][0]; r12 = q1[0][1]; r13 = q1[0][2];  /* [ r11 r12 r13 v1 ] */
  r21 = q1[1][0]; r22 = q1[1][1]; r23 = q1[1][2];  /* [ r21 r22 r23 v2 ] */
  r31 = q1[2][0]; r32 = q1[2][1]; r33 = q1[2][2];  /* [ r31 r32 r33 v3 ] */
  v1  = q1[0][3]; v2  = q1[1][3]; v3  = q1[2][3];  /* [  0   0   0   1 ] */
  
  deti = r11*r22*r33-r11*r32*r23-r21*r12*r33
  +r21*r32*r13+r31*r12*r23-r31*r22*r13 ;
  
  if( deti != 0.0 ) deti = 1.0 / deti ;
  
  q2[0][0] = deti*( r22*r33-r32*r23) ;
  q2[0][1] = deti*(-r12*r33+r32*r13) ;
  q2[0][2] = deti*( r12*r23-r22*r13) ;
  q2[0][3] = deti*(-r12*r23*v3+r12*v2*r33+r22*r13*v3
                    -r22*v1*r33-r32*r13*v2+r32*v1*r23) ;
  
  q2[1][0] = deti*(-r21*r33+r31*r23) ;
  q2[1][1] = deti*( r11*r33-r31*r13) ;
  q2[1][2] = deti*(-r11*r23+r21*r13) ;
  q2[1][3] = deti*( r11*r23*v3-r11*v2*r33-r21*r13*v3
                    +r21*v1*r33+r31*r13*v2-r31*v1*r23) ;
  
  q2[2][0] = deti*( r21*r32-r31*r22) ;
  q2[2][1] = deti*(-r11*r32+r31*r12) ;
  q2[2][2] = deti*( r11*r22-r21*r12) ;
  q2[2][3] = deti*(-r11*r22*v3+r11*r32*v2+r21*r12*v3
                    -r21*r32*v1-r31*r12*v2+r31*r22*v1) ;
  
  q2[3][0] = q2[3][1] = q2[3][2] = 0.0 ;
  q2[3][3] = (deti == 0.0) ? 0.0 : 1.0 ; /* failure flag if deti == 0 */
  
}


///multiply a vector of size 4 and a 4*4 matrix
void mult_4t4mat_4vec(float mat[4][4], float vectIni[4], float vectRes[4]){
  int lig,col;
  
  for (lig=0;lig<4;lig++){
    vectRes[lig]=0;
    for(col=0;col<4;col++)
      vectRes[lig]+=mat[lig][col]*vectIni[col];
	}
}

///multiply two 4*4 matrix: mat_i1 * mat_i2 -> mat_o
void mult_quat4t4mat_quat4t4mat(float mat_i1[4][4], float mat_i2[4][4], float mat_o[4][4]){
  int o1,o2,i;
  
  for (o1=0;o1<4;o1++) for (o2=0;o2<4;o2++){
    mat_o[o1][o2]=0;
    for (i=0;i<4;i++) mat_o[o1][o2]+=mat_i1[o1][i]*mat_i2[i][o2];
  }
}


///read a 4*4 matrix in a text file
void Read_quat4t4mat(char * FileName,float locmat[4][4]){
  int i,j;
  FILE *DataFile;
  
  //open file
  DataFile=fopen(FileName,"r");
  
  //go to the file beginning
  fseek(DataFile,0,SEEK_SET);
  
  //read the 4*4 matrix
  fscanf(DataFile,"%f	%f	%f	%f",&locmat[0][0],&locmat[0][1],&locmat[0][2],&locmat[0][3]);
  fscanf(DataFile,"%f	%f	%f	%f",&locmat[1][0],&locmat[1][1],&locmat[1][2],&locmat[1][3]);
  fscanf(DataFile,"%f	%f	%f	%f",&locmat[2][0],&locmat[2][1],&locmat[2][2],&locmat[2][3]);
  fscanf(DataFile,"%f	%f	%f	%f",&locmat[3][0],&locmat[3][1],&locmat[3][2],&locmat[3][3]);
  
  //close file
  fclose(DataFile);
}

///write a 4*4 matrix in a text file
void Write_quat4t4mat(char * FileName,float locmat[4][4]){
  int i,j;
  FILE *DataFile;
  
  //open file
  DataFile=fopen(FileName,"w");
  
  //read the 4*4 matrix
  fprintf(DataFile,"%f	%f	%f	%f\n",locmat[0][0],locmat[0][1],locmat[0][2],locmat[0][3]);
  fprintf(DataFile,"%f	%f	%f	%f\n",locmat[1][0],locmat[1][1],locmat[1][2],locmat[1][3]);
  fprintf(DataFile,"%f	%f	%f	%f\n",locmat[2][0],locmat[2][1],locmat[2][2],locmat[2][3]);
  fprintf(DataFile,"%f	%f	%f	%f\n",locmat[3][0],locmat[3][1],locmat[3][2],locmat[3][3]);
  
  //close file
  fclose(DataFile);
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                            8: LANDMARKS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// 8.1) ++++++++++++++++++ point landmarks ++++++++++++++++++ 

///constructor
LDMK_Points::LDMK_Points(void){
	this->LDMK_Points_Nb=0;
}

///destructor
LDMK_Points::~LDMK_Points(void){
  if ((this->Lx!=NULL)&&(this->LDMK_Points_Nb>0)) delete this->Lx;
  if ((this->Ly!=NULL)&&(this->LDMK_Points_Nb>0)) delete this->Ly;
  if ((this->Lz!=NULL)&&(this->LDMK_Points_Nb>0)) delete this->Lz;
  this->LDMK_Points_Nb=0;
}


///return the X coordinate of a Landmark
float LDMK_Points::GetX(int Id){
	
  if (Id<0) Id=0;
  if (Id>=this->LDMK_Points_Nb) Id=this->LDMK_Points_Nb-1;
  
  if (this->LDMK_Points_Nb==0) return 0;
  else return Lx[Id];
}  

///return the Y coordinate of a Landmark
float LDMK_Points::GetY(int Id){
	
  if (Id<0) Id=0;
  if (Id>=this->LDMK_Points_Nb) Id=this->LDMK_Points_Nb-1;
  
  if (this->LDMK_Points_Nb==0) return 0;
  else return Ly[Id];
}  

///return the Z coordinate of a Landmark
float LDMK_Points::GetZ(int Id){
	
  if (Id<0) Id=0;
  if (Id>=this->LDMK_Points_Nb) Id=this->LDMK_Points_Nb-1;
  
  if (this->LDMK_Points_Nb==0) return 0;
  else return Lz[Id];
}  

///Get the number of LDMK_Points
int LDMK_Points::Get_LDMK_PointsNumber(void){
  return this->LDMK_Points_Nb;
}  


///read LDMK_Points in a CSV file
///File should have this format:
///[Point 1: x]\t[Point 1: y]\t[Point 1:  z]
///[Point 2: x]\t[Point 2: y]\t[Point 2:  z]
///...
void LDMK_Points::Read(char * FileName){
  FILE *DataFile;
  float LocValue;
  int Nb_LDMK_Points;
  
  //allocate the memory for the LDMK_Points
  this->Lx = new float [10000];
  this->Ly = new float [10000];
  this->Lz = new float [10000];

  //open file
  DataFile=fopen(FileName,"r");
  
  //go to the file beginning
  fseek(DataFile,0,SEEK_SET);
  
  //read the file
  Nb_LDMK_Points=0;
  while(!feof(DataFile)){
    fscanf(DataFile,"%f	%f	%f",&this->Lx[Nb_LDMK_Points],&this->Ly[Nb_LDMK_Points],&this->Lz[Nb_LDMK_Points]);
    Nb_LDMK_Points++;
  }
  
  this->LDMK_Points_Nb=Nb_LDMK_Points-1;
  
  //check the read values
  cout << "Data read in " << FileName << ":" << endl; 
  for (Nb_LDMK_Points=0;Nb_LDMK_Points<this->LDMK_Points_Nb;Nb_LDMK_Points++)
    cout << this->Lx[Nb_LDMK_Points] << " "  << this->Ly[Nb_LDMK_Points] << " "  << this->Lz[Nb_LDMK_Points] << endl;

}


/// 8.2) ++++++++++++++++++ curve landmarks ++++++++++++++++++ 



///constructor
LDMK_Curves::LDMK_Curves(void){
	this->NbSeg=0;
}

///alternative constructor
LDMK_Curves::LDMK_Curves(int SegNb, int ElNb){
  int i,j;
  this->NbSeg=SegNb;
  
  //allocate the segments
  this->NbEl = new int [2*SegNb];
  this->x = new float* [2*SegNb];
  this->y = new float* [2*SegNb];
  this->z = new float* [2*SegNb];
  this->d = new float* [2*SegNb];
  
  //allocate the elements
  for (i=0;i<SegNb;i++){
    this->NbEl[i]=ElNb;
    this->x[i] = new float [ElNb];
    this->y[i] = new float [ElNb];
    this->z[i] = new float [ElNb];
    this->d[i] = new float [ElNb];
  }
  
  //give default values to the elements
  for (i=0;i<SegNb;i++) for (j=0;j<this->NbEl[i];j++){
    this->x[i][j] = 0;
    this->y[i][j] = 0;
    this->z[i][j] = 0;
    this->d[i][j] = 1;
  }
}

///destructor
LDMK_Curves::~LDMK_Curves(void){
  int i,j;
  
  for (i=0;i<this->NbSeg;i++) if (this->NbEl[i]>0){
      delete this->x[i];
      delete this->y[i];
      delete this->z[i];
      delete this->d[i];
  }
  
  delete this->NbEl;
  this->NbSeg=0;
}

///read a LDMK_Curves in a mv3d file
void LDMK_Curves::Read(char * Mv3dName){
  FILE *DataFile;
  char CT;
  char CharTestPrec;
  int NoLoc,i,j,k,m,n;
  double xLoc,yLoc,zLoc,dLoc;

  //open the file containing the curves (at the mv3d format)
  DataFile=fopen(Mv3dName,"rb");

  //1) count the number of segments
  cout << "Load the curves of "<< Mv3dName << endl;

  fseek(DataFile,0,SEEK_SET);  //start at the beginning of the file
  CT=fgetc(DataFile);
    
  while(!feof(DataFile)){
    //read the current character and store the previous one
    CharTestPrec=CT;
    CT=fgetc(DataFile);   //move to next character in DataFile
    
    //test
    if (CharTestPrec=='\n')
      if ((CT=='0')||(CT=='1')||(CT=='2')||(CT=='3')||(CT=='4')||(CT=='5')||(CT=='6')||(CT=='7')||(CT=='8')||(CT=='9')){
        fseek(DataFile,-1,SEEK_CUR);
        fscanf(DataFile,"%d	%lf	%lf	%lf	%lf",&NoLoc,&xLoc,&yLoc,&zLoc,&dLoc);  //move to next raw
        fseek(DataFile,-2,SEEK_CUR);   //move backward of 2 characters to reinitiate properly CT and CharTestPrec
        }
    }

  NoLoc++; //the last NoLoc is the number of segments minus one
  this->NbSeg=NoLoc;   //contains the number of segments

  //2) memory allocation (1/2)
  this->NbEl = new int [2*this->NbSeg];     //we allocate two times the amount of required memory, just in case...
  this->x = new float* [2*this->NbSeg];
  this->y = new float* [2*this->NbSeg];
  this->z = new float* [2*this->NbSeg];
  this->d = new float* [2*this->NbSeg];
  
  for (i=0;i<2*this->NbSeg;i++) this->NbEl[i]=0;
  
  //3) read the number of elements in each segment

  fseek(DataFile,0,SEEK_SET);  //go back to the beginning of the file
  CT=fgetc(DataFile);
    
  while(!feof(DataFile)){
    //read the current character and store the previous one
    CharTestPrec=CT;
    CT=fgetc(DataFile);   //move to next character in DataFile
    
    //test
    if (CharTestPrec=='\n')
      if ((CT=='0')||(CT=='1')||(CT=='2')||(CT=='3')||(CT=='4')||(CT=='5')||(CT=='6')||(CT=='7')||(CT=='8')||(CT=='9')){
        fseek(DataFile,-1,SEEK_CUR);
        fscanf(DataFile,"%d	%lf	%lf	%lf	%lf",&NoLoc,&xLoc,&yLoc,&zLoc,&dLoc);  //move to next raw
        fseek(DataFile,-2,SEEK_CUR);   //move backward of 2 characters to reinitiate properly CT and CharTestPrec
        this->NbEl[NoLoc]++;  //we add an element to the current segment
        }
    }


  //4) memory allocation (2/2)

  for(i=0;i<this->NbSeg;i++){
      this->x[i] = new float [this->NbEl[i]];
      this->y[i] = new float [this->NbEl[i]];
      this->z[i] = new float [this->NbEl[i]];
      this->d[i] = new float [this->NbEl[i]];
    }
    
  for (i=0;i<this->NbSeg;i++)
    this->NbEl[i]=0;
    
  //5) load the values of each element

  fseek(DataFile,0,SEEK_SET);  //go back to the beginning of the file
  CT=fgetc(DataFile);
  fseek(DataFile,1,SEEK_CUR);
    
  while(!feof(DataFile)){
    //read the current character and store the previous one
    CharTestPrec=CT;
    CT=fgetc(DataFile);   //move to next character in DataFile
    
    //test
    if (CharTestPrec=='\n')
      if ((CT=='0')||(CT=='1')||(CT=='2')||(CT=='3')||(CT=='4')||(CT=='5')||(CT=='6')||(CT=='7')||(CT=='8')||(CT=='9')){
        fseek(DataFile,-1,SEEK_CUR);
        fscanf(DataFile,"%d	%lf	%lf	%lf	%lf",&NoLoc,&xLoc,&yLoc,&zLoc,&dLoc);  //move to next raw
        fseek(DataFile,-2,SEEK_CUR);   //move backward of 2 characters to reinitiate properly CT and CharTestPrec
        this->x[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(xLoc);
        this->y[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(yLoc);
        this->z[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(zLoc);
        this->d[NoLoc][this->NbEl[NoLoc]]=static_cast<float>(dLoc);
        this->NbEl[NoLoc]++;
        }
    }


  //close DataFile
  fclose(DataFile);
  
}



///write a LDMK_Curves in a mv3d file.
///Set Preserve_IDs to 1 to preserve the original segment identifiers (default). They are optimaly resampled otherwise.
void LDMK_Curves::Write(char * Mv3dName, int Preserve_IDs){
  int i,iBis,j,tempInt;
  FILE *DataFile;
  int LocNbSeg,LocNbEl;
  int temp;
  int * NbSegEnds;
  int CaseSegEndI,CaseSegEndJ;
  int SegEndI,SegEndJ;
  int Nb1,Nb2,Nb3,Nb4,Nb5,Nb6,Nb7;
  int NbIntersections;

  
  //1) open the file in which the data will be saved
  DataFile=fopen(Mv3dName,"w");

  //2) header

  //2.1 - total number of segments and elements
  LocNbEl=0;
  LocNbSeg=0;
  for (i=0;i<2*this->NbSeg;i++){    //the "2 *" is in case there are additional segments
    if (this->NbEl[i]!=0){
      for (j=0;j<this->NbEl[i];j++)
	LocNbEl++;
      LocNbSeg++;
      }
    }

  //2.2 - count the number of intersections
  NbSegEnds = new int [4*(this->NbSeg)];  // we consider a maximum of 4*[the segments number]  for the number of segment-ends

  for (i=0;i<2*this->NbSeg;i++){
    NbSegEnds[2*i]=0;
    NbSegEnds[2*i+1]=0;
    }

  for (i=0;i<2*this->NbSeg-1;i++) if (this->NbEl[i]!=0) for(CaseSegEndI=0;CaseSegEndI<2;CaseSegEndI++){
    for (j=i+1;j<2*this->NbSeg;j++) if (this->NbEl[j]!=0) for(CaseSegEndJ=0;CaseSegEndJ<2;CaseSegEndJ++){
      //number of the current element
      SegEndI=(this->NbEl[i]-1)*CaseSegEndI;
      SegEndJ=(this->NbEl[j]-1)*CaseSegEndJ;
      
      //add the number of segment-ends if an intersection is found
      if (fabs(this->x[i][SegEndI]-this->x[j][SegEndJ])<0.1)
      if (fabs(this->y[i][SegEndI]-this->y[j][SegEndJ])<0.1)
      if (fabs(this->z[i][SegEndI]-this->z[j][SegEndJ])<0.1){
	  NbSegEnds[2*i+CaseSegEndI]++;
	NbSegEnds[2*j+CaseSegEndJ]++;
	}
      }
    }


  //for each segment-end, the number N of segment arriving on a node is: NbSegEnds[X]= sum_{i=1}^{i=N-1} i = N*(N-1)/2
  Nb2=0; Nb3=0; Nb4=0; Nb5=0; Nb6=0; Nb7=0;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for(CaseSegEndI=0;CaseSegEndI<2;CaseSegEndI++){
    if (NbSegEnds[2*i+CaseSegEndI]==0)  NbSegEnds[2*i+CaseSegEndI]=0;
    else if (NbSegEnds[2*i+CaseSegEndI]==1)   {NbSegEnds[2*i+CaseSegEndI]=2; Nb2++;}   //
    else if (NbSegEnds[2*i+CaseSegEndI]<=3)   {NbSegEnds[2*i+CaseSegEndI]=3; Nb3++;}   //
    else if (NbSegEnds[2*i+CaseSegEndI]<=6)   {NbSegEnds[2*i+CaseSegEndI]=4; Nb4++;}   // -> we use <= and not == for potential loops
    else if (NbSegEnds[2*i+CaseSegEndI]<=10)  {NbSegEnds[2*i+CaseSegEndI]=5; Nb5++;}  //
    else if (NbSegEnds[2*i+CaseSegEndI]<=15)  {NbSegEnds[2*i+CaseSegEndI]=6; Nb6++;}  //
    else if (NbSegEnds[2*i+CaseSegEndI]<=21)  {NbSegEnds[2*i+CaseSegEndI]=7; Nb7++;}  //
    }

  NbIntersections=Nb2/2+Nb3/3+Nb4/4+Nb5/5+Nb6/6+Nb7/7;
    

  //2.3) write header
  fprintf(DataFile,"# MicroVisu3D file\n");
  fprintf(DataFile,"# Number of lines   %d\n",LocNbSeg);
  fprintf(DataFile,"# Number of points  %d\n",LocNbEl);
  fprintf(DataFile,"# Number of inter.  %d\n",NbIntersections);
  fprintf(DataFile,"#\n");
  fprintf(DataFile,"# No		x		y		z		d\n");
  fprintf(DataFile,"#\n");

  //3) save data


  if (Preserve_IDs==1){
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
      for (j=0;j<this->NbEl[i];j++)
	fprintf(DataFile,"%d	%lf	%lf	%lf	%lf\n",i,static_cast<double>(this->x[i][j]),static_cast<double>(this->y[i][j]),static_cast<double>(this->z[i][j]),static_cast<double>(this->d[i][j]));
      fprintf(DataFile,"\n");
      }
    }
  else{
    tempInt=0;
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
      for (j=0;j<this->NbEl[i];j++)
	fprintf(DataFile,"%d	%lf	%lf	%lf	%lf\n",tempInt,static_cast<double>(this->x[i][j]),static_cast<double>(this->y[i][j]),static_cast<double>(this->z[i][j]),static_cast<double>(this->d[i][j]));
      tempInt++;
      fprintf(DataFile,"\n");
      }
    }

  fclose(DataFile);
  

}



///Delete a segment
void LDMK_Curves::DeleteSegment(int Seg){
  delete this->x[Seg];
  delete this->y[Seg];
  delete this->z[Seg];
  delete this->d[Seg];
  this->NbEl[Seg]=0;
}



///Merge two segments at their nearest extremity. The new segment is saved in Seg1.
void LDMK_Curves::MergeSegments(int Seg1, int Seg2){
  float * LocX;
  float * LocY;
  float * LocZ;
  float * LocD;
  int Side1,Side2,Side1_star,Side2_star;
  float tempFL,BestTempFL;
  int i;
  int NewElNb;
  
  if ((this->NbEl[Seg1]==0)||(this->NbEl[Seg2]==0)) return;

  
  //find the nearest connection
  for (Side1=0;Side1<2;Side1++) for (Side2=0;Side2<2;Side2++){
    tempFL=(this->x[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->x[Seg2][(this->NbEl[Seg2]-1)*Side2])*(this->x[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->x[Seg2][(this->NbEl[Seg2]-1)*Side2]);
    tempFL+=(this->y[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->y[Seg2][(this->NbEl[Seg2]-1)*Side2])*(this->y[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->y[Seg2][(this->NbEl[Seg2]-1)*Side2]);
    tempFL+=(this->z[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->z[Seg2][(this->NbEl[Seg2]-1)*Side2])*(this->z[Seg1][(this->NbEl[Seg1]-1)*Side1]-this->z[Seg2][(this->NbEl[Seg2]-1)*Side2]);
  
    if ((Side1==0)&&(Side2==0)){
      Side1_star=Side1;
      Side2_star=Side2;
      BestTempFL=tempFL;
    }
    
    if (BestTempFL>tempFL){
      Side1_star=Side1;
      Side2_star=Side2;
      BestTempFL=tempFL;
    }
  }
  
  
  //allocate memory for the temporary segment
  NewElNb=this->NbEl[Seg1]+this->NbEl[Seg2]-1;
  LocX=new float [NewElNb];
  LocY=new float [NewElNb];
  LocZ=new float [NewElNb];
  LocD=new float [NewElNb];
  
  //fill the temporary segment
  if (Side1_star==1){
    for (i=0;i<this->NbEl[Seg1];i++){
      LocX[i]=this->x[Seg1][i];
      LocY[i]=this->y[Seg1][i];
      LocZ[i]=this->z[Seg1][i];
      LocD[i]=this->d[Seg1][i];
    }
  }
  else{
    for (i=0;i<this->NbEl[Seg1];i++){
      LocX[i]=this->x[Seg1][this->NbEl[Seg1]-1-i];
      LocY[i]=this->y[Seg1][this->NbEl[Seg1]-1-i];
      LocZ[i]=this->z[Seg1][this->NbEl[Seg1]-1-i];
      LocD[i]=this->d[Seg1][this->NbEl[Seg1]-1-i];
    }
  }

  if (Side2_star==0){
    for (i=1;i<this->NbEl[Seg2];i++){
      LocX[this->NbEl[Seg1]+i-1]=this->x[Seg2][i];
      LocY[this->NbEl[Seg1]+i-1]=this->y[Seg2][i];
      LocZ[this->NbEl[Seg1]+i-1]=this->z[Seg2][i];
      LocD[this->NbEl[Seg1]+i-1]=this->d[Seg2][i];
    }
  }
  else{
    for (i=1;i<this->NbEl[Seg2];i++){
      LocX[this->NbEl[Seg1]+i-1]=this->x[Seg2][this->NbEl[Seg2]-1-i];
      LocY[this->NbEl[Seg1]+i-1]=this->y[Seg2][this->NbEl[Seg2]-1-i];
      LocZ[this->NbEl[Seg1]+i-1]=this->z[Seg2][this->NbEl[Seg2]-1-i];
      LocD[this->NbEl[Seg1]+i-1]=this->d[Seg2][this->NbEl[Seg2]-1-i];
    }
  }
  
  
  //delete Seg1 and Seg2
  this->DeleteSegment(Seg1);
  this->DeleteSegment(Seg2);
  
  //rebuild Seg1
  this->NbEl[Seg1]=NewElNb;
  this->x[Seg1] = new float [NewElNb];
  this->y[Seg1] = new float [NewElNb];
  this->z[Seg1] = new float [NewElNb];
  this->d[Seg1] = new float [NewElNb];
  
  for (i=0;i<NewElNb;i++){
    this->x[Seg1][i]=LocX[i];
    this->y[Seg1][i]=LocY[i];
    this->z[Seg1][i]=LocZ[i];
    this->d[Seg1][i]=LocD[i];
  }
  
  //dealloc temporary data
  delete LocX;
  delete LocY;
  delete LocZ;
  delete LocD;

}


///count the number of segments related to the segment end SegEnd (= 0 pr 1) of the segment Seg
int LDMK_Curves::CountLinkedSeg(int Seg,int SegEnd,float epsilon){
  int i,SideI;
  float LocX,LocY,LocZ;
  float LocX2,LocY2,LocZ2;
  int count;
  float tmpFl;
  
  if (SegEnd==0){ LocX=this->x[Seg][0];                 LocY=this->y[Seg][0];                 LocZ=this->z[Seg][0]; }
  if (SegEnd==1){ LocX=this->x[Seg][this->NbEl[Seg]-1]; LocY=this->y[Seg][this->NbEl[Seg]-1]; LocZ=this->z[Seg][this->NbEl[Seg]-1]; }
  
  count=0;
  for (i=0;i<2*this->NbSeg;i++) if ((this->NbEl[i]!=0)&&(i!=Seg)) for (SideI=0;SideI<2;SideI++){
    if (SideI==0){ LocX2=this->x[i][0];               LocY2=this->y[i][0];               LocZ2=this->z[i][0]; }
    if (SideI==1){ LocX2=this->x[i][this->NbEl[i]-1]; LocY2=this->y[i][this->NbEl[i]-1]; LocZ2=this->z[i][this->NbEl[i]-1]; }

    tmpFl=sqrt(((LocX-LocX2)*(LocX-LocX2))+((LocY-LocY2)*(LocY-LocY2))+((LocZ-LocZ2)*(LocZ-LocZ2)));
    
    if (tmpFl<epsilon)count++;
  }

  return count;

}

///Clean-up a moderately large LDMK_Curves structure
/// -> two segments ends are supposed linked if their distance is less than epsilon
/// -> merge the segments linked by a node with only two segments
/// -> remove the segments linked to only one node and for which the node has more than two segments related to other nodes
void LDMK_Curves::CleanUp(float epsilon){
  int changes;
  int i,j,j_star;
  int SideI, SideJ;
  float LocX,LocY,LocZ;
  float LocX2,LocY2,LocZ2;
  int count,count2;
  float tmpFl;
  
  //1) Remove isolated points
  //cout << "Remove isolated points" << endl;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]==1) this->DeleteSegment(i);
  
  //2) remove the segments linked to only one node and for which the node has more than two segments related to other nodes
  //cout << "Useless segments removal" << endl;
  int NgbhNbS0,NgbhNbS1,RefSide;
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
    NgbhNbS0=CountLinkedSeg(i,0,epsilon);
    NgbhNbS1=CountLinkedSeg(i,1,epsilon);
    
    //only one node related to something
    if ( ((NgbhNbS0==0)&&(NgbhNbS1!=0))  || ((NgbhNbS1==0)&&(NgbhNbS0!=0)) ){
      if (NgbhNbS0==0) RefSide=1;
      else RefSide=0;
      
      //find the ngbhs of 'i/RefSide' and treat them
      LocX=this->x[i][(this->NbEl[i]-1)*RefSide];
      LocY=this->y[i][(this->NbEl[i]-1)*RefSide];
      LocZ=this->z[i][(this->NbEl[i]-1)*RefSide];
      
      count=0;
      
      for (j=0;j<2*this->NbSeg;j++) if ((this->NbEl[j]!=0)&&(i!=j)) for (SideJ=0;SideJ<2;SideJ++){
	LocX2=this->x[j][(this->NbEl[j]-1)*SideJ];
	LocY2=this->y[j][(this->NbEl[j]-1)*SideJ];
	LocZ2=this->z[j][(this->NbEl[j]-1)*SideJ];
	
	tmpFl=sqrt(((LocX-LocX2)*(LocX-LocX2))+((LocY-LocY2)*(LocY-LocY2))+((LocZ-LocZ2)*(LocZ-LocZ2)));
	
	//'i/RefSide' is linked to 'j/SideJ'
	if (tmpFl<epsilon){
	  if (SideJ==0) count2=CountLinkedSeg(j,1,epsilon);
	  if (SideJ==1) count2=CountLinkedSeg(j,0,epsilon);
	  
	  //'j/abs(1-SideJ)' is related to something
	  if (count2>0) count++;
	}
	
      }
      
      if (count>=2){
	this->DeleteSegment(i);
	//cout << "Delete " << i << endl;
      }
    }
  }
  
  //3) merge the segments linked by a node with only two segments
  
  changes=1;
  while (changes==1){
    //cout << "Iteration of segments merging" << endl;
    changes=0;
    
    for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0){
      for (SideI=0;SideI<2;SideI++){
        
	//define the segment end treated
        if (SideI==0){ LocX=this->x[i][0];               LocY=this->y[i][0];               LocZ=this->z[i][0]; }
        if (SideI==1){ LocX=this->x[i][this->NbEl[i]-1]; LocY=this->y[i][this->NbEl[i]-1]; LocZ=this->z[i][this->NbEl[i]-1]; }
        
        //count the number of segment-ends linked to LocX,LocY,LocZ
	count=0;
	for (j=i+1;j<2*this->NbSeg;j++) if (this->NbEl[j]!=0) for (SideJ=0;SideJ<2;SideJ++){
	  //define the segment end treated
          if (SideJ==0){ LocX2=this->x[j][0];               LocY2=this->y[j][0];               LocZ2=this->z[j][0]; }
          if (SideJ==1){ LocX2=this->x[j][this->NbEl[j]-1]; LocY2=this->y[j][this->NbEl[j]-1]; LocZ2=this->z[j][this->NbEl[j]-1]; }
	  
	  //compute the distance and make the test
	  tmpFl=sqrt(((LocX-LocX2)*(LocX-LocX2))+((LocY-LocY2)*(LocY-LocY2))+((LocZ-LocZ2)*(LocZ-LocZ2)));
	  
	  if (tmpFl<epsilon){
	    count++;
	    j_star=j;
	  }
	}
	
	//merge the segments if only one segment is related to segment i
	if (count==1){
	  this->MergeSegments(i,j_star);
	  //cout << "Merge " << i << " " << j_star << endl;
	  changes=1;
	  SideI=2;
	  SideJ=2;
	}
      }
    }
  }
}

///transform the LDMK_Curves coordinates and diameters from voxels to mm according to the image 2 world properties of RefSF 
void LDMK_Curves::VoxelsToMillimeters(ScalarField * RefSF){
  int i,j;
  float x_new,y_new,z_new;
  float x_mm,y_mm,z_mm,vox_mm,three;
  
  three=3;
  
  x_mm=sqrt(RefSF->Image2World[0][0]*RefSF->Image2World[0][0]+RefSF->Image2World[0][1]*RefSF->Image2World[0][1]+RefSF->Image2World[0][2]*RefSF->Image2World[0][2]);
  y_mm=sqrt(RefSF->Image2World[1][0]*RefSF->Image2World[1][0]+RefSF->Image2World[1][1]*RefSF->Image2World[1][1]+RefSF->Image2World[1][2]*RefSF->Image2World[1][2]);
  z_mm=sqrt(RefSF->Image2World[2][0]*RefSF->Image2World[2][0]+RefSF->Image2World[2][1]*RefSF->Image2World[2][1]+RefSF->Image2World[2][2]*RefSF->Image2World[2][2]);
  vox_mm=(x_mm+y_mm+z_mm)/three;
  
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    x_new=this->x[i][j]*RefSF->Image2World[0][0]+this->y[i][j]*RefSF->Image2World[0][1]+this->z[i][j]*RefSF->Image2World[0][2]+RefSF->Image2World[0][3];
    y_new=this->x[i][j]*RefSF->Image2World[1][0]+this->y[i][j]*RefSF->Image2World[1][1]+this->z[i][j]*RefSF->Image2World[1][2]+RefSF->Image2World[1][3];
    z_new=this->x[i][j]*RefSF->Image2World[2][0]+this->y[i][j]*RefSF->Image2World[2][1]+this->z[i][j]*RefSF->Image2World[2][2]+RefSF->Image2World[2][3];

    this->x[i][j]=x_new;
    this->y[i][j]=y_new;
    this->z[i][j]=z_new;
    this->d[i][j]=this->d[i][j]*vox_mm;
  }

}

///transform the LDMK_Curves coordinates and diameters from mm to voxels according to the image 2 world properties of RefSF 
void LDMK_Curves::MillimetersToVoxels(ScalarField * RefSF){
  int i,j;
  float x_new,y_new,z_new;
  float x_mm,y_mm,z_mm,vox_mm,three;
  
  three=3;

  x_mm=sqrt(RefSF->Image2World[0][0]*RefSF->Image2World[0][0]+RefSF->Image2World[0][1]*RefSF->Image2World[0][1]+RefSF->Image2World[0][2]*RefSF->Image2World[0][2]);
  y_mm=sqrt(RefSF->Image2World[1][0]*RefSF->Image2World[1][0]+RefSF->Image2World[1][1]*RefSF->Image2World[1][1]+RefSF->Image2World[1][2]*RefSF->Image2World[1][2]);
  z_mm=sqrt(RefSF->Image2World[2][0]*RefSF->Image2World[2][0]+RefSF->Image2World[2][1]*RefSF->Image2World[2][1]+RefSF->Image2World[2][2]*RefSF->Image2World[2][2]);
  vox_mm=(x_mm+y_mm+z_mm)/three;
  
  for (i=0;i<2*this->NbSeg;i++) if (this->NbEl[i]!=0) for (j=0;j<this->NbEl[i];j++){
    x_new=this->x[i][j]*RefSF->World2Image[0][0]+this->y[i][j]*RefSF->World2Image[0][1]+this->z[i][j]*RefSF->World2Image[0][2]+RefSF->World2Image[0][3];
    y_new=this->x[i][j]*RefSF->World2Image[1][0]+this->y[i][j]*RefSF->World2Image[1][1]+this->z[i][j]*RefSF->World2Image[1][2]+RefSF->World2Image[1][3];
    z_new=this->x[i][j]*RefSF->World2Image[2][0]+this->y[i][j]*RefSF->World2Image[2][1]+this->z[i][j]*RefSF->World2Image[2][2]+RefSF->World2Image[2][3];

    this->x[i][j]=x_new;
    this->y[i][j]=y_new;
    this->z[i][j]=z_new;
    this->d[i][j]=this->d[i][j]/vox_mm;
  }
}


