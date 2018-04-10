/*=========================================================================
 
 Date      : $Date: 18.01.2011$
 Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$
 
 =========================================================================*/

#include <Demons_plain.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                   CONSTRUCTOR AND DESTRUCTOR
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
LargeDefDemons_plain::LargeDefDemons_plain(void){
	
  //default parameters
	iteration_nb=20;
  sigmaFluid=8;
  sigmaDiff=0;
  MaxUpdateAllowed=1;
  IndicatorMI=0;
  lambdaX=1;
	strcpy(PrefixInputs,"Null");
	strcpy(PrefixOutputs,"Outputs");
	strcpy(MovingImFile,"Null");
	strcpy(FixedImFile,"Null");
  World_FixedIm2MovingIm[0][0]=1; World_FixedIm2MovingIm[0][1]=0;   World_FixedIm2MovingIm[0][2]=0;    World_FixedIm2MovingIm[0][3]=0;   
  World_FixedIm2MovingIm[1][0]=0; World_FixedIm2MovingIm[1][1]=1;   World_FixedIm2MovingIm[1][2]=0;    World_FixedIm2MovingIm[1][3]=0;   
  World_FixedIm2MovingIm[2][0]=0; World_FixedIm2MovingIm[2][1]=0;   World_FixedIm2MovingIm[2][2]=1;    World_FixedIm2MovingIm[2][3]=0;   
  World_FixedIm2MovingIm[3][0]=0; World_FixedIm2MovingIm[3][1]=0;   World_FixedIm2MovingIm[3][2]=0;    World_FixedIm2MovingIm[3][3]=1;
  x_mm=1;
  y_mm=1;
  y_mm=1;
  UnderSampleFixedImFactor=1;
  strcpy(IniDispFieldX,"Null");
  strcpy(IniDispFieldY,"Null");
  strcpy(IniDispFieldZ,"Null");
  IniDispFieldDefined=0;
}

LargeDefDemons_plain::~LargeDefDemons_plain(void){}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                        SUB-FUNCTIONS TO PERFORM THE REGISTRATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///initiate the gradient descent
void LargeDefDemons_plain::ReadAndTreatInputImages(void){
	int x, y, z;
	int DistClosestEdge;
	int i,j;
	double mean1,mean2,std1,std2;
  float tempQuat[4][4];
  float TmpFl;
  ScalarField TempFixedIm;
	char FileTreatedFixedIm[256];
  char DotNii[256];
  
  
	//2) READ INPUTS
  //2.1.a) read moving image file
  this->MovingIm.Read(this->MovingImFile);
    
  //2.1.b) read fixed image file
  if (fabs(this->UnderSampleFixedImFactor-1)<=0.01){
    
    // the images do not required any downsampling or domain expansion
      this->FixedIm.Read(this->FixedImFile);
  }
  else{ // there is a special treatment to do with the FixedIm images
    //name of the treated FixedIm file
    strcpy(FileTreatedFixedIm,this->PrefixOutputs);
    strcpy(DotNii,"_tmp.nii");
    strcat(FileTreatedFixedIm,DotNii);
    
    cout << "After treatments (resampling and/or expansion): " <<  this->FixedImFile << " becomes " << FileTreatedFixedIm << endl;
      
    //undersample the image
    if (fabs(this->UnderSampleFixedImFactor-1)>0.01){
      TempFixedIm.Read_and_Undersample(this->FixedImFile,this->UnderSampleFixedImFactor);
      TempFixedIm.Write(FileTreatedFixedIm);
      
      strcpy(this->FixedImFile,FileTreatedFixedIm);  // a bit dirty but useful as the results are saved in the FixedIm coordinate space
      this->FixedIm.Read(this->FixedImFile);
    }
  }
  
  //2.1.a) allocate the memory for the deformed MovingIm
  this->DeformedMovingIm.Read(this->FixedImFile);
  
  //2.2) check whether  3D or 2D images are opened
  if (this->MovingIm.NT>1) cout << "Moving image  depends on time!!!";
  if (this->FixedIm.NT>1) cout << "Fixed image depends on time!!!";
  
	//2.3) variables containing the size of the Fixed image (and all the other images, except the moving image)
	this->NX=this->FixedIm.NX;
	this->NY=this->FixedIm.NY;
	this->NZ=this->FixedIm.NZ;
	this->NT=1;
  
  //2.4) variables containing the size of the moving image
	this->NXs=this->MovingIm.NX;
	this->NYs=this->MovingIm.NY;
	this->NZs=this->MovingIm.NZ;
  
	cout << "Image size: " << this->NX <<  "*"  <<  this->NY  <<  "*"  << this->NZ  << " (Moving image: " << this->NXs <<  "*"  <<  this->NYs  <<  "*"  << this->NZs  << ")\n";
	
  //2.5) compute the quaternion to convert fixed image coordinates into moving image coordinates
  mult_quat4t4mat_quat4t4mat(World_FixedIm2MovingIm,this->FixedIm.Image2World,tempQuat);
  mult_quat4t4mat_quat4t4mat(this->MovingIm.World2Image,tempQuat,FixedIm2MovingImCoord);
  
  cout << endl;
  
  if (IniDispFieldDefined==1){
    cout << "FixedIm to MovingIm encoded in world coordinate by " << this->IniDispFieldX << ", " <<  this->IniDispFieldY << " and " <<  this->IniDispFieldZ << endl;
  }
  else{
    cout << "FixedIm to MovingIm in voxel coordinates:" << endl;
    for (i=0;i<4;i++){
      for (j=0;j<4;j++){
        cout << FixedIm2MovingImCoord[i][j] << " ";
      }
      cout << endl;
    }
  }
  
  
  //2.6 compute the voxels size in mm
  this->x_mm=sqrt(this->FixedIm.Image2World[0][0]*this->FixedIm.Image2World[0][0]+this->FixedIm.Image2World[0][1]*this->FixedIm.Image2World[0][1]+this->FixedIm.Image2World[0][2]*this->FixedIm.Image2World[0][2]);
  this->y_mm=sqrt(this->FixedIm.Image2World[1][0]*this->FixedIm.Image2World[1][0]+this->FixedIm.Image2World[1][1]*this->FixedIm.Image2World[1][1]+this->FixedIm.Image2World[1][2]*this->FixedIm.Image2World[1][2]);
  this->z_mm=sqrt(this->FixedIm.Image2World[2][0]*this->FixedIm.Image2World[2][0]+this->FixedIm.Image2World[2][1]*this->FixedIm.Image2World[2][1]+this->FixedIm.Image2World[2][2]*this->FixedIm.Image2World[2][2]);
  
  cout << endl;
  cout << "Fixed image resolution: " << this->x_mm << " "  << this->y_mm << " "  << this->z_mm << endl;

  
  
	//4) COMPUTE THE MAPPING (in from the fixed image c.s. to the moving image c.s.)
  if (IniDispFieldDefined==1) this->ReadAndConvertMapping();
  
  

  //7) LOAD THE DEFORMATION FIELD OR DISPLACEMENT FIELD
  //    -->  VelocityField.G(0,x,y,z,i)= direction ex of the vector at (x,y,z)
	//    -->  VelocityField.G(1,x,y,z,i)= direction ey of the vector at (x,y,z)
	//    -->  VelocityField.G(2,x,y,z,i)= direction ez of the vector at (x,y,z)
	if (strcmp(PrefixInputs,"Null")!=0){
      this->LoadVelocityField(PrefixInputs);   //! the velocity field is put in VelocityField !
  }
	else{
      this->VelocityField.CreateVoidField(this->NX,this->NY,this->NZ);
  }
  
  
  
  
}

///allocate all variables used for the gradient descent
void LargeDefDemons_plain::AllocateAllVariables(void){
  int i;
	
  //contains the gradient of Energy. Is directely used to update the Displacement Field
  this->GradE.CreateVoidField(this->NX,this->NY,this->NZ);
  
  
  //manage the mutual information
  if (IndicatorMI==1){
    if (IniDispFieldDefined==1){
        ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);
    }
    else{
        ProjectImageUsingAffineTransfoAndSteadyVeloField(FixedIm2MovingImCoord,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);
    }
    
    this->NorMutInfMan.Initiate(&this->DeformedMovingIm,&this->FixedIm,50,50,3);
  }
  
  //initiate the FFT convolvers

    //allocate the memory
    if (this->sigmaFluid>0.05) this->FFTconvolver_fluid.InitiateConvolver(this->NX,this->NY,this->NZ,1,sigmaFluid/this->x_mm,sigmaFluid/this->y_mm,sigmaFluid/this->z_mm);
    
    if (this->sigmaDiff>0.05) this->FFTconvolver_Diff.InitiateConvolver(this->NX,this->NY,this->NZ,1,sigmaDiff/this->x_mm,sigmaDiff/this->y_mm,sigmaDiff/this->z_mm);
}



///save the result of the gradient descent
void LargeDefDemons_plain::ReadAndConvertMapping(void){
  float flX,flY,flZ;
  float WTX,WTY,WTZ;
  float WSX,WSY,WSZ;
  float ISX,ISY,ISZ;
  int x,y,z;
  
  //read the mapping in the world coordinate system   (IniDispField:  points from FixedIm to MovingIm in millimeters)
  this->IniDispField.Read_and_Interpolate(this->IniDispFieldX,this->IniDispFieldY,this->IniDispFieldZ,this->NX,this->NY,this->NZ,0);
  
  //convert the mapping in voxel coordonate system (points from FixedIm to moving image in voxels)
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    
    flX=static_cast<float>(x); flY=static_cast<float>(y); flZ=static_cast<float>(z);
    
    WTX=flX*this->FixedIm.Image2World[0][0]+flY*this->FixedIm.Image2World[0][1]+flZ*this->FixedIm.Image2World[0][2]+this->FixedIm.Image2World[0][3];
    WTY=flX*this->FixedIm.Image2World[1][0]+flY*this->FixedIm.Image2World[1][1]+flZ*this->FixedIm.Image2World[1][2]+this->FixedIm.Image2World[1][3];
    WTZ=flX*this->FixedIm.Image2World[2][0]+flY*this->FixedIm.Image2World[2][1]+flZ*this->FixedIm.Image2World[2][2]+this->FixedIm.Image2World[2][3];
    
    WSX=WTX+this->IniDispField.G(0,x,y,z);
    WSY=WTY+this->IniDispField.G(1,x,y,z);
    WSZ=WTZ+this->IniDispField.G(2,x,y,z);
    
    ISX=WSX*this->MovingIm.World2Image[0][0]+WSY*this->MovingIm.World2Image[0][1]+WSZ*this->MovingIm.World2Image[0][2]+this->MovingIm.World2Image[0][3];
    ISY=WSX*this->MovingIm.World2Image[1][0]+WSY*this->MovingIm.World2Image[1][1]+WSZ*this->MovingIm.World2Image[1][2]+this->MovingIm.World2Image[1][3];
    ISZ=WSX*this->MovingIm.World2Image[2][0]+WSY*this->MovingIm.World2Image[2][1]+WSZ*this->MovingIm.World2Image[2][2]+this->MovingIm.World2Image[2][3];
    
    this->IniDispField.P(ISX,0,x,y,z);
    this->IniDispField.P(ISY,1,x,y,z);
    this->IniDispField.P(ISZ,2,x,y,z);
  }
  
  //IniDispField.Write("IniDispFieldDemVoxX.nii","IniDispFieldDemVoxY.nii","IniDispFieldDemVoxZ.nii",this->FixedImFile);
}


///save the result of the gradient descent
void LargeDefDemons_plain::SaveResultGradientDescent(void){
  VectorField RealDisplacementField;
  int i,j;
  float TempQuat_T[4][4];
  float TempQuat_S[4][4];
  
  
  //2) save the results...
  
    //save the velocity field
  	this->SaveVelocityField(&this->VelocityField,this->PrefixOutputs);
    
    //save the final deformation
    this->SaveDeformations(this->PrefixOutputs);
    
    //save the inverse deformation field
    RealDisplacementField.CreateVoidField(this->NX,this->NY,this->NZ);
    CptInvDefFromSteadyVeloField(&this->VelocityField,&RealDisplacementField,7);
 
    //CptDefFromSteadyVeloField(&this->VelocityField,&RealDisplacementField,7);
    
    this->SaveInvTotalDisplacement(&RealDisplacementField,this->PrefixOutputs);
    

    
  
	

}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                          FUNCTIONS TO SAVE AND LOAD THE VARIOUS STRUCTURES
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




///load the velocity fields
void LargeDefDemons_plain::LoadVelocityField(char Prefix[256]){
	char FileNameX[256];
	char FileNameY[256];
	char FileNameZ[256];
	char DisplacementField_X[256];
	char DisplacementField_Y[256];
	char DisplacementField_Z[256];

	//1) intialisation
	strcpy(FileNameX,Prefix);
	strcpy(DisplacementField_X,"_VelocityField_X.nii");
	strcat(FileNameX,DisplacementField_X);
	strcpy(FileNameY,Prefix);
	strcpy(DisplacementField_Y,"_VelocityField_Y.nii");
	strcat(FileNameY,DisplacementField_Y);
	strcpy(FileNameZ,Prefix);
	strcpy(DisplacementField_Z,"_VelocityField_Z.nii");
	strcat(FileNameZ,DisplacementField_Z);
	this->VelocityField.Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
}




///save the inverse deformation field
void LargeDefDemons_plain::SaveInvTotalDisplacement(VectorField * LocDisplacementField,char Prefix[256]){
	char FileNameX[256];
	char FileNameY[256];
	char FileNameZ[256];
	char InvDisplacementField_X[256];
	char InvDisplacementField_Y[256];
	char InvDisplacementField_Z[256];
  float srcX,srcY,srcZ;
  float trgX,trgY,trgZ;
  float tmpX,tmpY,tmpZ;
  float tempQuat[4][4];
  int x,y,z;
  float flX,flY,flZ;
  float tmpX2,tmpY2,tmpZ2;
  
	//intialisation
	strcpy(FileNameX,Prefix);
	strcpy(InvDisplacementField_X,"_DispField_Fix2Mov_X.nii");
	strcat(FileNameX,InvDisplacementField_X);
	strcpy(FileNameY,Prefix);
	strcpy(InvDisplacementField_Y,"_DispField_Fix2Mov_Y.nii");
	strcat(FileNameY,InvDisplacementField_Y);
	strcpy(FileNameZ,Prefix);
	strcpy(InvDisplacementField_Z,"_DispField_Fix2Mov_Z.nii");
	strcat(FileNameZ,InvDisplacementField_Z);
	
  
  
  
  
  if (IniDispFieldDefined!=1){  //+++++++ Case A: Only an initial affine alignment is considered between the images +++++++
    
    //compute the inverse of the total displacement field in mm (world coordinate)
    mult_quat4t4mat_quat4t4mat(World_FixedIm2MovingIm,this->FixedIm.Image2World,tempQuat);
    
    //compute the mapping
    for (z = 0; z < LocDisplacementField->NZ; z++)  for (y = 0; y < LocDisplacementField->NY; y++) for (x = 0; x < LocDisplacementField->NX; x++){
      trgX=x*this->FixedIm.Image2World[0][0]+y*this->FixedIm.Image2World[0][1]+z*this->FixedIm.Image2World[0][2]+this->FixedIm.Image2World[0][3];
      trgY=x*this->FixedIm.Image2World[1][0]+y*this->FixedIm.Image2World[1][1]+z*this->FixedIm.Image2World[1][2]+this->FixedIm.Image2World[1][3];
      trgZ=x*this->FixedIm.Image2World[2][0]+y*this->FixedIm.Image2World[2][1]+z*this->FixedIm.Image2World[2][2]+this->FixedIm.Image2World[2][3];
      
      tmpX=x+LocDisplacementField->G(0,x,y,z);
      tmpY=y+LocDisplacementField->G(1,x,y,z);
      tmpZ=z+LocDisplacementField->G(2,x,y,z);
      
      srcX=tmpX*tempQuat[0][0]+tmpY*tempQuat[0][1]+tmpZ*tempQuat[0][2]+tempQuat[0][3];
      srcY=tmpX*tempQuat[1][0]+tmpY*tempQuat[1][1]+tmpZ*tempQuat[1][2]+tempQuat[1][3];
      srcZ=tmpX*tempQuat[2][0]+tmpY*tempQuat[2][1]+tmpZ*tempQuat[2][2]+tempQuat[2][3];
      
      LocDisplacementField->P(srcX-trgX,0,x,y,z);
      LocDisplacementField->P(srcY-trgY,1,x,y,z);
      LocDisplacementField->P(srcZ-trgZ,2,x,y,z);
    }
  
  }
  else{ //+++++++ Case B: The initial alignment between the images is encoded in a displacement field +++++++
   
    for (z = 0; z < LocDisplacementField->NZ; z++)  for (y = 0; y < LocDisplacementField->NY; y++) for (x = 0; x < LocDisplacementField->NX; x++){
      flX=static_cast<float>(x); flY=static_cast<float>(y); flZ=static_cast<float>(z);
      
      trgX=x*this->FixedIm.Image2World[0][0]+y*this->FixedIm.Image2World[0][1]+z*this->FixedIm.Image2World[0][2]+this->FixedIm.Image2World[0][3];
      trgY=x*this->FixedIm.Image2World[1][0]+y*this->FixedIm.Image2World[1][1]+z*this->FixedIm.Image2World[1][2]+this->FixedIm.Image2World[1][3];
      trgZ=x*this->FixedIm.Image2World[2][0]+y*this->FixedIm.Image2World[2][1]+z*this->FixedIm.Image2World[2][2]+this->FixedIm.Image2World[2][3];
      
      tmpX=flX+LocDisplacementField->G(0,x,y,z);
      tmpY=flY+LocDisplacementField->G(1,x,y,z);
      tmpZ=flZ+LocDisplacementField->G(2,x,y,z);
      
      tmpX2=this->IniDispField.G(0,tmpX,tmpY,tmpZ);
      tmpY2=this->IniDispField.G(1,tmpX,tmpY,tmpZ);
      tmpZ2=this->IniDispField.G(2,tmpX,tmpY,tmpZ);
      
      srcX=tmpX2*this->MovingIm.Image2World[0][0]+tmpY2*this->MovingIm.Image2World[0][1]+tmpZ2*this->MovingIm.Image2World[0][2]+this->MovingIm.Image2World[0][3];
      srcY=tmpX2*this->MovingIm.Image2World[1][0]+tmpY2*this->MovingIm.Image2World[1][1]+tmpZ2*this->MovingIm.Image2World[1][2]+this->MovingIm.Image2World[1][3];
      srcZ=tmpX2*this->MovingIm.Image2World[2][0]+tmpY2*this->MovingIm.Image2World[2][1]+tmpZ2*this->MovingIm.Image2World[2][2]+this->MovingIm.Image2World[2][3];
      
      LocDisplacementField->P(srcX-trgX,0,x,y,z);
      LocDisplacementField->P(srcY-trgY,1,x,y,z);
      LocDisplacementField->P(srcZ-trgZ,2,x,y,z);
    }
    
  }
	
  
  //save the deformation field
	LocDisplacementField->Write(FileNameX,FileNameY,FileNameZ,FixedImFile);
}

///save the Velocity field
void LargeDefDemons_plain::SaveVelocityField(VectorField * LocVelocityField,char Prefix[256]){
	char FileNameX[256];
	char FileNameY[256];
	char FileNameZ[256];
	char VelocityField_X[256];
	char VelocityField_Y[256];
	char VelocityField_Z[256];
  
	//intialisation
	strcpy(FileNameX,Prefix);
	strcpy(VelocityField_X,"_VelocityField_X.nii");
	strcat(FileNameX,VelocityField_X);
	strcpy(FileNameY,Prefix);
	strcpy(VelocityField_Y,"_VelocityField_Y.nii");
	strcat(FileNameY,VelocityField_Y);
	strcpy(FileNameZ,Prefix);
	strcpy(VelocityField_Z,"_VelocityField_Z.nii");
	strcat(FileNameZ,VelocityField_Z);

	//save the Velocity field
	LocVelocityField->Write(FileNameX,FileNameY,FileNameZ,FixedImFile);
}

///save the deformations in time subdivisions (not the convergence)
void LargeDefDemons_plain::SaveDeformations(char Prefix[256]){
  char FileName[256];
	char Deformation[256];

  
	//deformed image
  strcpy(Deformation,"_DeformedMovIma.nii");
  strcpy(FileName,Prefix);
  strcat(FileName,Deformation);
    
  if (IniDispFieldDefined==1){
      ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);
  }
  else{
      ProjectImageUsingAffineTransfoAndSteadyVeloField(FixedIm2MovingImCoord,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);
  }
  
  
  
  
  this->DeformedMovingIm.Write(FileName,this->FixedImFile);
}



///compute the update vector field
void LargeDefDemons_plain::ComputeUpdateFieldSSD(){
  float locCoef;
  float locCoef2;
  int x,y,z;
  float epsilon;
  
  
  //init variables and gradient
  epsilon=0.0001;
  
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->GradE.P(0,0,x,y,z);
    this->GradE.P(0,1,x,y,z);
    this->GradE.P(0,2,x,y,z);
  }
  
  
  //project the moving image in DeformedMovingIm 
  //compute the temporary deformed MovingIm
  
  if (IniDispFieldDefined==1)
      ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);
  else
      ProjectImageUsingAffineTransfoAndSteadyVeloField(FixedIm2MovingImCoord,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);

  
  //compute the temporary gradient of the deformed MovingIm
  Cpt_Grad_ScalarField(&this->DeformedMovingIm,&this->GradE);
  
  float Av1st_term;
  float Av2nd_term;
  int NbPts;
  
  Av1st_term=0;
  Av2nd_term=0;
  NbPts=0;
  
  //multiply the temporary gradient by the difference between the FixedIm and deformed images AND normalize the temporary gradient (3D IMAGE)
  if (this->NZ>1) for (z = 1; z < this->NZ-1; z++)  for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++){
    locCoef=-(this->FixedIm.G(x,y,z)-this->DeformedMovingIm.G(x,y,z));
    locCoef2=this->GradE.G(0,x,y,z)*this->GradE.G(0,x,y,z)+this->GradE.G(1,x,y,z)*this->GradE.G(1,x,y,z)+this->GradE.G(2,x,y,z)*this->GradE.G(2,x,y,z);
    
    
    if ((locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)))>epsilon){
      locCoef = locCoef/(locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)));
      Av1st_term+=locCoef2;
      Av2nd_term+=locCoef*locCoef;
      NbPts++;
    }
    else 
      locCoef=0;
      
    this->GradE.P(this->GradE.G(0,x,y,z)*locCoef,0,x,y,z);
    this->GradE.P(this->GradE.G(1,x,y,z)*locCoef,1,x,y,z);
    this->GradE.P(this->GradE.G(2,x,y,z)*locCoef,2,x,y,z);
  }
  
  //multiply the temporary gradient by the difference between the FixedIm and deformed images AND normalize the temporary gradient (2D IMAGE)
  if (this->NZ==1) for (y = 1; y < this->NY-1; y++) for (x = 1; x < this->NX-1; x++){
    locCoef=-(this->FixedIm.G(x,y)-this->DeformedMovingIm.G(x,y));
    locCoef2=this->GradE.G(0,x,y)*this->GradE.G(0,x,y)+this->GradE.G(1,x,y)*this->GradE.G(1,x,y)+this->GradE.G(2,x,y)*this->GradE.G(2,x,y);
    
    
    if ((locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)))>epsilon){
      locCoef = locCoef/(locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)));
      Av1st_term+=locCoef2;
      Av2nd_term+=locCoef*locCoef;
      NbPts++;
    }
    else 
      locCoef=0;
    
    this->GradE.P(this->GradE.G(0,x,y)*locCoef,0,x,y);
    this->GradE.P(this->GradE.G(1,x,y)*locCoef,1,x,y);
    this->GradE.P(this->GradE.G(2,x,y)*locCoef,2,x,y);
  }
  
}





///compute the update vector field
void LargeDefDemons_plain::ComputeUpdateFieldMI(){
  float evaMI;
  int x,y,z;
  
  //compute the current deformed moving image
  if (IniDispFieldDefined==1){
      ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);
  }
  else{
      ProjectImageUsingAffineTransfoAndSteadyVeloField(FixedIm2MovingImCoord,&this->VelocityField,&this->MovingIm,&this->DeformedMovingIm);
  }  
  
  this->NorMutInfMan.IndicateSrcHasChanged();
  
  //compute the mutual information and update the histograms
  this->NorMutInfMan.IndicateSrcHasChanged();
  evaMI=this->NorMutInfMan.EvaluateMI(); 
  cout << "MI: " << evaMI << endl;
  
  
  //evaluate the gradient of mutual information
  this->NorMutInfMan.EvaluateGradMI(&this->GradE); 
  
  
}



///control the maximum update. 
///* if the max update is larger than the initial one -> normalisation of all the update field to this->MaxUpdateAllowed
///* Othewise -> normalisation of all the update field to this->MaxUpdateAllowed*(max/InitMaxUpdate)
void LargeDefDemons_plain::ControlMaxUpdate(int IterationNb){
  float max;
  int direc,x,y,z;
  float factor;
  
  max=0;
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
    if (max<this->GradE.G(0,x,y,z)*this->GradE.G(0,x,y,z)+this->GradE.G(1,x,y,z)*this->GradE.G(1,x,y,z)+this->GradE.G(2,x,y,z)*this->GradE.G(2,x,y,z))
      max=this->GradE.G(0,x,y,z)*this->GradE.G(0,x,y,z)+this->GradE.G(1,x,y,z)*this->GradE.G(1,x,y,z)+this->GradE.G(2,x,y,z)*this->GradE.G(2,x,y,z);
  
  max=sqrt(max);
  
  if (IterationNb==0) InitMaxUpdate=max;
  
  if (max>InitMaxUpdate) factor=this->MaxUpdateAllowed/max;
  else factor=this->MaxUpdateAllowed/InitMaxUpdate;
  
  for (direc = 0; direc < 3; direc++) for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) 
    this->GradE.P(this->GradE.G(direc,x,y,z)*factor,direc,x,y,z);
}



///Estimate the reference update scale for a given sigma (in the smoothing kernel) 
float LargeDefDemons_plain::EstimRefUpdateScale(float sigmaloc,LightFFTconvolver3D * ConvolverLoc){
  float max,refWght;
  int d,x,y,z;
  
  //1 Compute the gradient
  if (IndicatorMI==1)
    this->ComputeUpdateFieldMI();
  else 
    this->ComputeUpdateFieldSSD();
  
  //2 smooth the gradient with the proposed sigma
  ConvolverLoc->ChangeKernel(1,sigmaloc/this->x_mm,sigmaloc/this->y_mm,sigmaloc/this->z_mm);
  
  ConvolverLoc->Convolution(&this->GradE);
  
  //3 compute the maximum gradient of energy
  max=0;
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
    if (max<this->GradE.G(0,x,y,z)*this->GradE.G(0,x,y,z)+this->GradE.G(1,x,y,z)*this->GradE.G(1,x,y,z)+this->GradE.G(2,x,y,z)*this->GradE.G(2,x,y,z))
      max=this->GradE.G(0,x,y,z)*this->GradE.G(0,x,y,z)+this->GradE.G(1,x,y,z)*this->GradE.G(1,x,y,z)+this->GradE.G(2,x,y,z)*this->GradE.G(2,x,y,z);
  
  max=sqrt(max);
  refWght=1/max;

  return refWght;
}




///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                      RUN FUNCTIONS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



///Function to solve the registration using the gradient descent algorithm of Beg 05
void LargeDefDemons_plain::Run(void){
	int IterationNb;
  
	
	//1) INITIALISATION
	
	//1.1) Pre-treatment of the inuput images 
	this->ReadAndTreatInputImages();
	
	//1.2) Allocations of the scalar and vector fields + definition of global parameters
	this->AllocateAllVariables();
	
	
	//2) GRADIENT DESCENT
	for (IterationNb=0;IterationNb<this->iteration_nb;IterationNb++){
		cout << "Iteration Number " << IterationNb+1 << " / " << this->iteration_nb << "\n";
		
    //2.1 Compute the gradient
    if (IndicatorMI==1)
      this->ComputeUpdateFieldMI();
    else{
      this->ComputeUpdateFieldSSD();
    }
    
    //2.2 smooth the gradient... 
    if (this->sigmaFluid>0.05) this->FFTconvolver_fluid.Convolution(&this->GradE);
    
    //2.3 control the maximum amplitude of update
    this->ControlMaxUpdate(IterationNb);
    
    //2.4 update the transformation
    //ComposeTwoLogFieldsUsingBCH(&this->VelocityField,&this->GradE);
    ComposeTwoLogFieldsUsingSum(&this->VelocityField,&this->GradE);
    
    //2.5 smooth the velocities
    if (this->sigmaDiff>0.05) this->FFTconvolver_Diff.Convolution(&this->VelocityField);
  }
	
	//3) SAVE THE RESULTS
	this->SaveResultGradientDescent();
}



