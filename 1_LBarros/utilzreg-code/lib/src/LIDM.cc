/*=========================================================================
 
 Authors: Laurent Risser, Francois-Xavier Vialard
 
 =========================================================================*/

#include <LIDM.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                   CONSTRUCTOR AND DESTRUCTOR
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
LIDM::LIDM(void){
  int i;
  
  //default parameters
  iteration_nb=10;
  NbTimeSubdiv=7;
  MaxVelocityUpdate=0.49;  //rem: Delta Voxels = 1
  Margin=0;
  VFpenalizer=0.999;
  RefMaxGrad=-1.;
  FinalDefInvVec=0;
  ShowSSD=0;
  PreserveWeights=0;
  MovingPOU=0;
  strcpy(PrefixInputs,"Null");
  strcpy(PrefixOutputs,"Outputs");
  strcpy(SourceFile,"Null");
  strcpy(TargetFile,"Null");
  strcpy(PartiOfUnityFile,"Null");
  World_Target2Template[0][0]=1; World_Target2Template[0][1]=0;   World_Target2Template[0][2]=0;    World_Target2Template[0][3]=0;   
  World_Target2Template[1][0]=0; World_Target2Template[1][1]=1;   World_Target2Template[1][2]=0;    World_Target2Template[1][3]=0;   
  World_Target2Template[2][0]=0; World_Target2Template[2][1]=0;   World_Target2Template[2][2]=1;    World_Target2Template[2][3]=0;   
  World_Target2Template[3][0]=0; World_Target2Template[3][1]=0;   World_Target2Template[3][2]=0;    World_Target2Template[3][3]=1; 
  x_mm=1;
  y_mm=1;
  y_mm=1;
  UnderSampleFactor=1;
  SigmaPOU=5;
}

LIDM::~LIDM(void){}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                        SUB-FUNCTIONS TO PERFORM THE REGISTRATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///initiate the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
void LIDM::ReadAndTreatInputImages(void){
  int x, y, z;
  int DistClosestEdge;
  int i,j;
  double mean1,mean2,std1,std2;
  ScalarField Mask;
  float tempQuat[4][4];
  float World_Template2Target[4][4];
  ScalarField TempImSrc;
  char FileTreatedSrc[256];
  char DotNii[256];
  

  //2) READ INPUTS
  //2.1) read files
  
  ImTarget.Read(this->TargetFile);
  
  if (fabs(this->UnderSampleFactor-1)<=0.01){ // no undersampling
    ImTemplate.Read(this->SourceFile);
    PartiOfUnity.Read(this->PartiOfUnityFile);
  }
  else{ // undersample the source/template image and the partition of unity ...
    
    // ... partition of unity
    if (SigmaPOU>0) PartiOfUnity.Read_and_Undersample(this->PartiOfUnityFile,this->UnderSampleFactor,1);  //nearest neighbor interpolation
    else PartiOfUnity.Read_and_Undersample(this->PartiOfUnityFile,this->UnderSampleFactor,0);    //trilinear interpolation
    
    // ... source image  (this trick is made as the results are saved in the source image coordinate system)
    //strcpy(FileTreatedSrc,"tmp_");
    //strcpy(DotNii,".nii");
    //strcat(FileTreatedSrc,this->PrefixOutputs);
    //strcat(FileTreatedSrc,DotNii);
    strcpy(FileTreatedSrc,this->PrefixOutputs);
    strcpy(DotNii,"_tmp.nii");
    strcat(FileTreatedSrc,DotNii);
    
    TempImSrc.Read_and_Undersample(this->SourceFile,this->UnderSampleFactor);
    TempImSrc.Write(FileTreatedSrc);
        
    strcpy(this->SourceFile,FileTreatedSrc);
    this->ImTemplate.Read(this->SourceFile);
    
    cout << "FYI: Resampled " <<  this->SourceFile << " is saved in " << FileTreatedSrc << endl;
  }
  
    
  //2.2) check whether  3D or 2D images are opened
  if (ImTemplate.NT>1) cout << "Source image " << " depends on time!!!";
  if (ImTarget.NT>1) cout << "Target image "  << " depends on time!!!";
    
  
  //2.5) variables containing the size of the image
  this->NX=ImTemplate.NX;
  this->NY=ImTemplate.NY;
  this->NZ=ImTemplate.NZ;
  this->NT=1;
  
  cout << "Image size: " << this->NX << "*" << this->NY << "*" << this->NZ << " voxels  (target: " << this->ImTarget.NX <<  "*"  <<  this->ImTarget.NY  <<  "*"  << this->ImTarget.NZ  << " voxels)\n";
  
  
  //2.6) compute the quaternion to convert target coordinates into template coordinates
  
  invert_4t4quaternion(this->World_Target2Template,World_Template2Target);
  
  mult_quat4t4mat_quat4t4mat(World_Template2Target,ImTemplate.Image2World,tempQuat);
  mult_quat4t4mat_quat4t4mat(ImTarget.World2Image,tempQuat,Template2TargetCoord);
  
  cout << endl;
  cout << "Template to target:" << endl;
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      cout << Template2TargetCoord[i][j] << " ";
    }
    cout << endl;
  }
  
  //2.7 compute the voxels size in mm
  this->x_mm=sqrt(ImTemplate.Image2World[0][0]*ImTemplate.Image2World[0][0]+ImTemplate.Image2World[0][1]*ImTemplate.Image2World[0][1]+ImTemplate.Image2World[0][2]*ImTemplate.Image2World[0][2]);
  this->y_mm=sqrt(ImTemplate.Image2World[1][0]*ImTemplate.Image2World[1][0]+ImTemplate.Image2World[1][1]*ImTemplate.Image2World[1][1]+ImTemplate.Image2World[1][2]*ImTemplate.Image2World[1][2]);
  this->z_mm=sqrt(ImTemplate.Image2World[2][0]*ImTemplate.Image2World[2][0]+ImTemplate.Image2World[2][1]*ImTemplate.Image2World[2][1]+ImTemplate.Image2World[2][2]*ImTemplate.Image2World[2][2]);
  
  cout << endl;
  cout << "Template image resolution: " << this->x_mm << " "  << this->y_mm << " "  << this->z_mm << " millimeters" << endl;
  
}


///allocate all variables used for the gradient descent (Beg 2005) of the current 3D image from the treated 4D time sequence.
///Compute also the dimension of the scalar and vector fields in use
void LIDM::AllocateAllVariables(void){
  int i,j,k;
  float s1_x,s1_y,s1_z,s2_x,s2_y,s2_z,s3_x,s3_y,s3_z,s4_x,s4_y,s4_z,s5_x,s5_y,s5_z,s6_x,s6_y,s6_z,s7_x,s7_y,s7_z,w1,w2,w3,w4,w5,w6,w7;
  LightFFTconvolver3D  TmpLightFFTconvolver;
  int NeedForInit;
  float normalFactor;
  char fileName[256];
  char APOU[256];
  
  //time step between two subdivision
  this->DeltaTimeSubdiv=1./(static_cast<float>(NbTimeSubdiv-1));
  
  //4) initiate the velocity field
  //... velocity field
  //    -->  VelocityField.G(0,x,y,z,i)= direction ex of the vector at (x,y,z)
  //    -->  VelocityField.G(1,x,y,z,i)= direction ey of the vector at (x,y,z)
  //    -->  VelocityField.G(2,x,y,z,i)= direction ez of the vector at (x,y,z)
  //    -->  where n is the id of the velocity field
  if (strcmp(PrefixInputs,"Null")!=0)
    this->LoadVelocityFields(PrefixInputs);  //NbTimeSubdiv should be checked
  else
    this->VelocityField.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
  
  //... forward mapping
  //    -->  ForwardMapping.G(0,x,y,z,i) = coordinate x at time i corresponding to (x,y,z) at time 0
  //    -->  ForwardMapping.G(1,x,y,z,i) = coordinate y at time i corresponding to (x,y,z) at time 0
  //    -->  ForwardMapping.G(2,x,y,z,i) = coordinate z at time i corresponding to (x,y,z) at time 0
  this->ForwardMapping.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
  
  //... backward mapping
  //    -->  BackwardMapping.G(0,x,y,z,i) = coordinate x at time i corresponding to (x,y,z) at time 1
  //    -->  BackwardMapping.G(1,x,y,z,i) = coordinate y at time i corresponding to (x,y,z) at time 1
  //    -->  BackwardMapping.G(2,x,y,z,i) = coordinate z at time i corresponding to (x,y,z) at time 1
  this->BackwardMapping.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);
  
  //... temporary image transformed using the forward mapping from time 0
  //    -->  J0.G(x,y,z) = gray level of the transformed image J0 at (x,y,z)
  this->J0.CreateVoidField(this->NX,this->NY,this->NZ);
  
  //... temporary image transformed using the backward mapping from time 1
  //    -->  J1.G(x,y,z) = gray level of the transformed image J1 at (x,y,z)
  this->J1.CreateVoidField(this->NX,this->NY,this->NZ);
  
  //... gradient of J
  //    -->  GradJ.G(0,x,y,z)= gradient of J0 in direction ex at (x,y,z)
  //    -->  GradJ.G(1,x,y,z)= gradient of J0 in direction ey at (x,y,z)
  //    -->  GradJ.G(2,x,y,z)= gradient of J0 in direction ez at (x,y,z)
  this->GradJ.CreateVoidField(this->NX,this->NY,this->NZ);
  
  //... determinent of the Jacobians  
  //    -->  Jacobians.G(x,y,z)= determinant of the jacobian at (x,y,z)
  this->DetJacobians.CreateVoidField(this->NX,this->NY,this->NZ);
  
  //... Energy Gradient
  //    -->  GradE.G(0,i,x,y,z) = Energy gradient at time i in direction ex at (x,y,z)
  //    -->  GradE.G(1,i,x,y,z) = Energy gradient at time i in direction ey at (x,y,z)
  //    -->  GradE.G(2,i,x,y,z) = Energy gradient at time i in direction ez at (x,y,z)
  this->GradE.CreateVoidField(this->NX,this->NY,this->NZ,this->NbTimeSubdiv);




  //7) Initiate the class to smooth the images
  
  //7.1) initiate the convolver with the partition of unity
  if (this->SigmaPOU<=0) 
    LIDMConvolver.InitiateConvolverWithActualPOI(&this->PartiOfUnity);
  else
    LIDMConvolver.InitiateConvolver(&this->PartiOfUnity,this->SigmaPOU);
  
  //7.1.2) Save the actual Partition of Unity
  strcpy(fileName,this->PrefixOutputs); strcpy(APOU,"_ActualPOU.nii"); strcat(fileName,APOU);
  LIDMConvolver.SaveActualParitionOfUnity(fileName,this->SourceFile);

  //7.1.3) Store in memory the actual Partition of Unity and allocate memory for its projection
  if (this->MovingPOU==1){
    this->ActualPOU.Read(fileName);
    this->ProjectionActualPOU.Read(fileName);
  }

  //7.2) tune the different kernels in the different regions
  TmpLightFFTconvolver.InitiateConvolver(this->PartiOfUnity.NX,this->PartiOfUnity.NY,this->PartiOfUnity.NZ);
  
  NeedForInit=1;
  
  for (i=0;i<LIDMConvolver.GetRegionNb();i++){ //loop on the regions
    cout << "Region " << i << ": ";
    
    //7.2.1) cast the scales from millimeters to voxels in region $i$
    s1_x=this->stdDev[i][0]/this->x_mm; s1_y=this->stdDev[i][0]/this->y_mm; s1_z=this->stdDev[i][0]/this->z_mm;
    s2_x=this->stdDev[i][1]/this->x_mm; s2_y=this->stdDev[i][1]/this->y_mm; s2_z=this->stdDev[i][1]/this->z_mm;
    s3_x=this->stdDev[i][2]/this->x_mm; s3_y=this->stdDev[i][2]/this->y_mm; s3_z=this->stdDev[i][2]/this->z_mm;
    s4_x=this->stdDev[i][3]/this->x_mm; s4_y=this->stdDev[i][3]/this->y_mm; s4_z=this->stdDev[i][3]/this->z_mm;
    s5_x=this->stdDev[i][4]/this->x_mm; s5_y=this->stdDev[i][4]/this->y_mm; s5_z=this->stdDev[i][4]/this->z_mm;
    s6_x=this->stdDev[i][5]/this->x_mm; s6_y=this->stdDev[i][5]/this->y_mm; s6_z=this->stdDev[i][5]/this->z_mm;
    s7_x=this->stdDev[i][6]/this->x_mm; s7_y=this->stdDev[i][6]/this->y_mm; s7_z=this->stdDev[i][6]/this->z_mm;
    
    //7.2.2) compute the real weights from the apparent weights...
    if (PreserveWeights==1){
		w1=this->weight[i][0]; w2=this->weight[i][1]; w3=this->weight[i][2]; w4=this->weight[i][3]; w5=this->weight[i][4]; w6=this->weight[i][5]; w7=this->weight[i][6];
		cout << "s=" << this->stdDev[i][0] << "mm,  aw=" << w1;
		if (this->weight[i][1]>0.00011) cout << " / s=" << this->stdDev[i][1] << "mm,  aw=" << w2;
		if (this->weight[i][2]>0.00011) cout << " / s=" << this->stdDev[i][2] << "mm,  aw=" << w3;
		if (this->weight[i][3]>0.00011) cout << " / s=" << this->stdDev[i][3] << "mm,  aw=" << w4;
		if (this->weight[i][4]>0.00011) cout << " / s=" << this->stdDev[i][4] << "mm,  aw=" << w5;
		if (this->weight[i][5]>0.00011) cout << " / s=" << this->stdDev[i][5] << "mm,  aw=" << w6;
		if (this->weight[i][6]>0.00011) cout << " / s=" << this->stdDev[i][6] << "mm,  aw=" << w7;
		cout << endl;
	  }
	else{
      //... region $i$ / kernel 1
      if (NeedForInit==1){
        normalFactor=this->EstimRefWeight(this->stdDev[i][0],&TmpLightFFTconvolver,NeedForInit);
        w1=this->weight[i][0];  // = this->weight[i][0]*EstimRefWeight(...)/normalFactor
        NeedForInit=0;
      }
      else{
        w1=this->weight[i][0]*this->EstimRefWeight(this->stdDev[i][0],&TmpLightFFTconvolver,NeedForInit)/normalFactor;
	  }
      cout << "s=" << this->stdDev[i][0] << "mm,  w=" << w1 << " (aw=" << this->weight[i][0] << ")";
      
      //... region $i$ / kernel 2
      if (this->weight[i][1]>0.00011){
        w2=this->weight[i][1]*this->EstimRefWeight(this->stdDev[i][1],&TmpLightFFTconvolver,NeedForInit)/normalFactor;
        cout << " / s=" << this->stdDev[i][1] << "mm,  w=" << w2 << " (aw=" << this->weight[i][1] << ")";
      }
      else w2=0;
      
      //... region $i$ / kernel 3
      if (this->weight[i][2]>0.00011){
        w3=this->weight[i][2]*this->EstimRefWeight(this->stdDev[i][2],&TmpLightFFTconvolver,NeedForInit)/normalFactor;
        cout << " / s=" << this->stdDev[i][2] << "mm,  w=" << w3 << " (aw=" << this->weight[i][2] << ")";
      }
      else w3=0;
      
      //... region $i$ / kernel 4
      if (this->weight[i][3]>0.00011){
        w4=this->weight[i][3]*this->EstimRefWeight(this->stdDev[i][3],&TmpLightFFTconvolver,NeedForInit)/normalFactor;
        cout << " / s=" << this->stdDev[i][3] << "mm,  w=" << w4 << " (aw=" << this->weight[i][3] << ")";
      }
      else w4=0;
      
      //... region $i$ / kernel 5
      if (this->weight[i][4]>0.00011){
        w5=this->weight[i][4]*this->EstimRefWeight(this->stdDev[i][4],&TmpLightFFTconvolver,NeedForInit)/normalFactor;
        cout << " / s=" << this->stdDev[i][4] << "mm,  w=" << w5 << " (aw=" << this->weight[i][4] << ")";
      }
      else w5=0;
      
      //... region $i$ / kernel 6
      if (this->weight[i][5]>0.00011){
        w6=this->weight[i][5]*this->EstimRefWeight(this->stdDev[i][5],&TmpLightFFTconvolver,NeedForInit)/normalFactor;
        cout << " / s=" << this->stdDev[i][5] << "mm,  w=" << w6 << " (aw=" << this->weight[i][5] << ")";
      }
      else w6=0;
      
      //... region $i$ / kernel 7
      if (this->weight[i][6]>0.00011){
        w7=this->weight[i][6]*this->EstimRefWeight(this->stdDev[i][6],&TmpLightFFTconvolver,NeedForInit)/normalFactor;
        cout << " / s=" << this->stdDev[i][6] << "mm,  w=" << w7 << " (aw=" << this->weight[i][6] << ")";
      }
      else w7=0;
      
      cout << endl;
    }
    
    //7.2.3) Tune the kernel in region $i$
    //cout << "Reg " << i << ": "  << w1 << " " << s1_x << " " << s1_y << " " << s1_z << " / ";
    //cout << w2 << " " << s2_x << " " << s2_y << " " << s2_z << " / ";
    //cout << w3 << " " << s3_x << " " << s3_y << " " << s3_z << endl;
    LIDMConvolver.ChangeKernelInOneRegion(i,w1,s1_x,s1_y,s1_z,w2,s2_x,s2_y,s2_z,w3,s3_x,s3_y,s3_z,w4,s4_x,s4_y,s4_z,w5,s5_x,s5_y,s5_z,w6,s6_x,s6_y,s6_z,w7,s7_x,s7_y,s7_z);
  }

}



///Compute the energy gradients
void LIDM::ComputeEnergyGradient(int timeSubdiv){
  int x,y,z,i,k;
  float temp;
    
  //loop on the vector directions (x,y,z)
  if (this->NZ>1){
	  for (i=0;i<3;i++) for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      if ((z<this->Margin)||(z>this->NZ-this->Margin-1)||(y<this->Margin)||(y>this->NY-this->Margin-1)||(x<this->Margin)||(x>this->NX-this->Margin-1))
        temp=0;
      else 
        temp=(this->J0.G(x,y,z) - this->J1.G(x,y,z)) * this->DetJacobians.G(x,y,z) * this->GradJ.G(i,x,y,z);
      this->GradJ.P(static_cast<float>(temp),i,x,y,z);
    }
  }
  else{
	  for (i=0;i<3;i++) for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
        if ((y<this->Margin)||(y>this->NY-this->Margin-1)||(x<this->Margin)||(x>this->NX-this->Margin-1))
          temp=0;
        else 
          temp=(this->J0.G(x,y,z) - this->J1.G(x,y,z)) * this->DetJacobians.G(x,y,z) * this->GradJ.G(i,x,y,z);
        this->GradJ.P(static_cast<float>(temp),i,x,y,z);
      }
    }
  
  this->LIDMConvolver.Convolution(&this->GradJ);
  
  for (i=0;i<3;i++) for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
    this->GradE.P(-2*this->GradJ.G(i,x,y,z),i,x,y,z,timeSubdiv);
}



///project the partition of unity
void LIDM::ProjectPOU(int timeSubdiv){
    
  Project3Dimage(&this->ActualPOU,&this->ForwardMapping,&this->ProjectionActualPOU,timeSubdiv);
  this->LIDMConvolver.UpdatePartitionOfUnity(&this->ProjectionActualPOU);
  
  //if (timeSubdiv==this->NbTimeSubdiv-1) this->ProjectionActualPOU.Write("DefPOU.nii");
}





///Update VelocityField -> part 1: multiply the current field with 'this->VFpenalizer'
void LIDM::UpdateVelocityField_part1(){
  int x, y, z, i, k;
  float MaxGrad,MultFactor;
  double LocGrad;
  
  //...3D images
  for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    this->VelocityField.P(this->VFpenalizer*this->VelocityField.G(0,x,y,z,i),0,x,y,z,i);
    this->VelocityField.P(this->VFpenalizer*this->VelocityField.G(1,x,y,z,i),1,x,y,z,i);
    this->VelocityField.P(this->VFpenalizer*this->VelocityField.G(2,x,y,z,i),2,x,y,z,i);
  }
  
  //...2D images
  if (this->NZ==1){
    for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
      this->VelocityField.P(this->VFpenalizer*this->VelocityField.G(0,x,y,0,i),0,x,y,0,i);
      this->VelocityField.P(this->VFpenalizer*this->VelocityField.G(1,x,y,0,i),1,x,y,0,i);
    }
  }
  
}


///Update VelocityField -> part 2: add the energy gradients.
///Return MaxGrad/this->RefMaxGrad
float LIDM::UpdateVelocityField_part2(int IterationNb){
  int x, y, z, i, k;
  float MaxGrad,MultFactor;
  double LocGrad;
  
  //1) Compute the maximum of gradient in all time frames...
  //...3D images
  MaxGrad=0;
  for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    LocGrad=sqrt(pow((double)this->GradE.G(0,x,y,z,i),2.)+pow((double)this->GradE.G(1,x,y,z,i),2.)+pow((double)this->GradE.G(2,x,y,z,i),2.));
    if (MaxGrad<(float)LocGrad) MaxGrad=(float)LocGrad;
  }
  
  //...2D images
  if (this->NZ==1){
    for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
      LocGrad=sqrt(pow((double)this->GradE.G(0,x,y,0,i),2.)+pow((double)this->GradE.G(1,x,y,0,i),2.));
      if (MaxGrad<(float)LocGrad) MaxGrad=(float)LocGrad;
    }
  }
  
  //2) maximum update control at the first iteration
  if ((IterationNb==0)&&(RefMaxGrad<0.)) {
    this->RefMaxGrad=MaxGrad;
    if (this->RefMaxGrad==0){
      cout << "It seems that the registered images are identical\n";
      this->RefMaxGrad=1;
    }
    cout << "\n\nRefMaxGrad is set to " << this->RefMaxGrad << ". Keep this value if you continue these\n";
    cout << "computations (using -PrefixInputs) to manage well the convergence.\n\n";
  }
  
  //3) compute the MultFactor
  if (MaxGrad>this->RefMaxGrad) MultFactor=this->MaxVelocityUpdate/MaxGrad;
  else MultFactor=this->MaxVelocityUpdate/(this->RefMaxGrad);
  
  //4) Messages
  cout << " -> MaxGrad/RefMaxGrad=" << MaxGrad/this->RefMaxGrad  << "  |  "; //endl; //"   (update forces considered only)\n";
  
  //5) update the vector field...
  //...3D images
  for (i=0;i<this->NbTimeSubdiv;i++) for (z = 1; z < this->NZ-2; z++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
    this->VelocityField.P(this->VelocityField.G(0,x,y,z,i)-this->GradE.G(0,x,y,z,i)*MultFactor,0,x,y,z,i);
    this->VelocityField.P(this->VelocityField.G(1,x,y,z,i)-this->GradE.G(1,x,y,z,i)*MultFactor,1,x,y,z,i);
    this->VelocityField.P(this->VelocityField.G(2,x,y,z,i)-this->GradE.G(2,x,y,z,i)*MultFactor,2,x,y,z,i);
  }
  
  //...2D images
  if (this->NZ==1){
    for (i=0;i<this->NbTimeSubdiv;i++) for (y = 1; y < this->NY-2; y++) for (x = 1; x < this->NX-2; x++){
      this->VelocityField.P(this->VelocityField.G(0,x,y,0,i)-this->GradE.G(0,x,y,0,i)*MultFactor,0,x,y,0,i);
      this->VelocityField.P(this->VelocityField.G(1,x,y,0,i)-this->GradE.G(1,x,y,0,i)*MultFactor,1,x,y,0,i);
    }
  }
  
  return MaxGrad/this->RefMaxGrad;
}


///save the result of the gradient descent (Beg 2005) for the current 3D image of the 4D time sequence
void LIDM::SaveResultGradientDescent(void){
  //init -> compute the forward mapping and import the original input template (non pre-treated)
  CptMappingFromVeloField_IniIdMap(0,&this->VelocityField,&this->ForwardMapping);
  
  this->SaveVecDeformation(this->PrefixOutputs);
  
  //whole transformations
  this->SaveVelocityFields(&this->VelocityField,this->PrefixOutputs);
  this->SaveFinalDeformation(this->PrefixOutputs);
  if (this->FinalDefInvVec==1) this->SaveInvVecDeformation(this->PrefixOutputs);
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                          FUNCTIONS TO SAVE AND LOAD THE VARIOUS STRUCTURES
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///load the velocity fields
void LIDM::LoadVelocityFields(char Prefix[256]){
  char FileNameX[256];
  char FileNameY[256];
  char FileNameZ[256];
  char VelocityField_X[256];
  char VelocityField_Y[256];
  char VelocityField_Z[256];
  
  //1) intialisation
  strcpy(FileNameX,Prefix);
  strcpy(VelocityField_X,"_VelocityField_X.nii");
  strcat(FileNameX,VelocityField_X);
  strcpy(FileNameY,Prefix);
  strcpy(VelocityField_Y,"_VelocityField_Y.nii");
  strcat(FileNameY,VelocityField_Y);
  strcpy(FileNameZ,Prefix);
  strcpy(VelocityField_Z,"_VelocityField_Z.nii");
  strcat(FileNameZ,VelocityField_Z);
  
  this->VelocityField.Read_and_Interpolate(FileNameX,FileNameY,FileNameZ,this->NX,this->NY,this->NZ,1);
}





///save the velocity fields
void LIDM::SaveVelocityFields(VectorField * VelocityFieldLoc,char Prefix[256]){
  char FileNameX[256];
  char FileNameY[256];
  char FileNameZ[256];
  char VelocityField_X[256];
  char VelocityField_Y[256];
  char VelocityField_Z[256];
  
  //intialisation
  strcpy(FileNameX,Prefix);//this->PrefixOutputs
  strcpy(VelocityField_X,"_VelocityField_X.nii");
  strcat(FileNameX,VelocityField_X);
  strcpy(FileNameY,Prefix);
  strcpy(VelocityField_Y,"_VelocityField_Y.nii");
  strcat(FileNameY,VelocityField_Y);
  strcpy(FileNameZ,Prefix);
  strcpy(VelocityField_Z,"_VelocityField_Z.nii");
  strcat(FileNameZ,VelocityField_Z);
  
  //save the velocity field
  VelocityFieldLoc->Write(FileNameX,FileNameY,FileNameZ,this->SourceFile);
}




///save the deformations in time subdivisions (not the convergence)
void LIDM::SaveFinalDeformation(char Prefix[256]){
  int x, y, z;
  char FileName[256];
  char FinalDef[256];
  ScalarField source_image;
  
  strcpy(FinalDef,"_DeformedSource.nii");
  
  source_image.Read(this->SourceFile);
  Project3Dimage(&source_image,&this->ForwardMapping,&this->J0,this->NbTimeSubdiv-1);
  
  strcpy(FileName,Prefix);
  strcat(FileName,FinalDef);
  this->J0.Write(FileName,this->SourceFile);
}





///save the vector field that transforms [source] into [target]
void LIDM::SaveVecDeformation(char Prefix[256]){
  int x, y, z;
  VectorField Temp3DField;
  char FileName_X[256];
  char FileName_Y[256];
  char FileName_Z[256];
  char VecDef_X[256];
  char VecDef_Y[256];
  char VecDef_Z[256];
  float flX,flY,flZ;
  float srcX,srcY,srcZ;
  float trgX,trgY,trgZ;
  float tmpX,tmpY,tmpZ;
  float tmpX2,tmpY2,tmpZ2;
  float World_Template2Target[4][4];
  
  //intialisation
  CptMappingFromVeloField_IniIdMap(0,&this->VelocityField,&this->BackwardMapping);
  
  Temp3DField.CreateVoidField(this->NX, this->NY, this->NZ);
  
  strcpy(VecDef_X,"_DispField_Src2Trg_X.nii");
  strcpy(VecDef_Y,"_DispField_Src2Trg_Y.nii");
  strcpy(VecDef_Z,"_DispField_Src2Trg_Z.nii");
  
  invert_4t4quaternion(this->World_Target2Template,World_Template2Target);
  
  //save the forward mapping
  for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    flX=static_cast<float>(x); flY=static_cast<float>(y); flZ=static_cast<float>(z);
    
    srcX=flX*ImTemplate.Image2World[0][0]+flY*ImTemplate.Image2World[0][1]+flZ*ImTemplate.Image2World[0][2]+ImTemplate.Image2World[0][3];
    srcY=flX*ImTemplate.Image2World[1][0]+flY*ImTemplate.Image2World[1][1]+flZ*ImTemplate.Image2World[1][2]+ImTemplate.Image2World[1][3];
    srcZ=flX*ImTemplate.Image2World[2][0]+flY*ImTemplate.Image2World[2][1]+flZ*ImTemplate.Image2World[2][2]+ImTemplate.Image2World[2][3];
    
    tmpX=this->BackwardMapping.G(0,x,y,z,this->BackwardMapping.NT-1);
    tmpY=this->BackwardMapping.G(1,x,y,z,this->BackwardMapping.NT-1);
    tmpZ=this->BackwardMapping.G(2,x,y,z,this->BackwardMapping.NT-1);

    tmpX2=tmpX*ImTemplate.Image2World[0][0]+tmpY*ImTemplate.Image2World[0][1]+tmpZ*ImTemplate.Image2World[0][2]+ImTemplate.Image2World[0][3];
    tmpY2=tmpX*ImTemplate.Image2World[1][0]+tmpY*ImTemplate.Image2World[1][1]+tmpZ*ImTemplate.Image2World[1][2]+ImTemplate.Image2World[1][3];
    tmpZ2=tmpX*ImTemplate.Image2World[2][0]+tmpY*ImTemplate.Image2World[2][1]+tmpZ*ImTemplate.Image2World[2][2]+ImTemplate.Image2World[2][3];
    
    trgX=tmpX2*World_Template2Target[0][0]+tmpY2*World_Template2Target[0][1]+tmpZ2*World_Template2Target[0][2]+World_Template2Target[0][3];
    trgY=tmpX2*World_Template2Target[1][0]+tmpY2*World_Template2Target[1][1]+tmpZ2*World_Template2Target[1][2]+World_Template2Target[1][3];
    trgZ=tmpX2*World_Template2Target[2][0]+tmpY2*World_Template2Target[2][1]+tmpZ2*World_Template2Target[2][2]+World_Template2Target[2][3];
    
    Temp3DField.P(trgX-srcX,0,x,y,z);
    Temp3DField.P(trgY-srcY,1,x,y,z);
    Temp3DField.P(trgZ-srcZ,2,x,y,z);
    
  }
  
  strcpy(FileName_X,Prefix);
  strcat(FileName_X,VecDef_X);
  strcpy(FileName_Y,Prefix);
  strcat(FileName_Y,VecDef_Y);
  strcpy(FileName_Z,Prefix);
  strcat(FileName_Z,VecDef_Z);
  
  Temp3DField.Write(FileName_X,FileName_Y,FileName_Z,this->SourceFile);
}





///save the vector field that transforms [target] into [source]
void LIDM::SaveInvVecDeformation(char Prefix[256]){
  int x, y, z;
  VectorField Temp3DField;
  char FileName_X[256];
  char FileName_Y[256];
  char FileName_Z[256];
  char VecDef_X[256];
  char VecDef_Y[256];
  char VecDef_Z[256];
  float flX,flY,flZ;
  float srcX,srcY,srcZ;
  float trgX,trgY,trgZ;
  float tmpX,tmpY,tmpZ;
  float tmpX2,tmpY2,tmpZ2;
  float World_Template2Target[4][4];
  
  //intialisation
  CptMappingFromVeloField_IniIdMap(this->VelocityField.NT-1,&this->VelocityField,&this->BackwardMapping);
  
  Temp3DField.CreateVoidField(this->ImTarget.NX, this->ImTarget.NY, this->ImTarget.NZ);
  
  invert_4t4quaternion(this->World_Target2Template,World_Template2Target);
  
  //save the forward mapping
  for (z = 0; z < this->ImTarget.NZ; z++) for (y = 0; y < this->ImTarget.NY; y++) for (x = 0; x < this->ImTarget.NX; x++){
    flX=static_cast<float>(x); flY=static_cast<float>(y); flZ=static_cast<float>(z);
    
    trgX=flX*ImTarget.Image2World[0][0]+flY*ImTarget.Image2World[0][1]+flZ*ImTarget.Image2World[0][2]+ImTarget.Image2World[0][3];
    trgY=flX*ImTarget.Image2World[1][0]+flY*ImTarget.Image2World[1][1]+flZ*ImTarget.Image2World[1][2]+ImTarget.Image2World[1][3];
    trgZ=flX*ImTarget.Image2World[2][0]+flY*ImTarget.Image2World[2][1]+flZ*ImTarget.Image2World[2][2]+ImTarget.Image2World[2][3];
    
    tmpX=trgX*this->World_Target2Template[0][0]+trgY*this->World_Target2Template[0][1]+trgZ*this->World_Target2Template[0][2]+this->World_Target2Template[0][3];
    tmpY=trgX*this->World_Target2Template[1][0]+trgY*this->World_Target2Template[1][1]+trgZ*this->World_Target2Template[1][2]+this->World_Target2Template[1][3];
    tmpZ=trgX*this->World_Target2Template[2][0]+trgY*this->World_Target2Template[2][1]+trgZ*this->World_Target2Template[2][2]+this->World_Target2Template[2][3];
    
    tmpX2=tmpX*ImTemplate.World2Image[0][0]+tmpY*ImTemplate.World2Image[0][1]+tmpZ*ImTemplate.World2Image[0][2]+ImTemplate.World2Image[0][3];
    tmpY2=tmpX*ImTemplate.World2Image[1][0]+tmpY*ImTemplate.World2Image[1][1]+tmpZ*ImTemplate.World2Image[1][2]+ImTemplate.World2Image[1][3];
    tmpZ2=tmpX*ImTemplate.World2Image[2][0]+tmpY*ImTemplate.World2Image[2][1]+tmpZ*ImTemplate.World2Image[2][2]+ImTemplate.World2Image[2][3];
    
    tmpX=this->BackwardMapping.G(0,tmpX2,tmpY2,tmpZ2,0);
    tmpY=this->BackwardMapping.G(1,tmpX2,tmpY2,tmpZ2,0);
    tmpZ=this->BackwardMapping.G(2,tmpX2,tmpY2,tmpZ2,0);
    
    srcX=tmpX*ImTemplate.Image2World[0][0]+tmpY*ImTemplate.Image2World[0][1]+tmpZ*ImTemplate.Image2World[0][2]+ImTemplate.Image2World[0][3];
    srcY=tmpX*ImTemplate.Image2World[1][0]+tmpY*ImTemplate.Image2World[1][1]+tmpZ*ImTemplate.Image2World[1][2]+ImTemplate.Image2World[1][3];
    srcZ=tmpX*ImTemplate.Image2World[2][0]+tmpY*ImTemplate.Image2World[2][1]+tmpZ*ImTemplate.Image2World[2][2]+ImTemplate.Image2World[2][3];
    
    
    
    Temp3DField.P(srcX-trgX,0,x,y,z);
    Temp3DField.P(srcY-trgY,1,x,y,z);
    Temp3DField.P(srcZ-trgZ,2,x,y,z);
  }
  
  strcpy(VecDef_X,"_DispField_Trg2Src_X.nii");
  strcpy(VecDef_Y,"_DispField_Trg2Src_Y.nii");
  strcpy(VecDef_Z,"_DispField_Trg2Src_Z.nii");
  strcpy(FileName_X,Prefix);
  strcat(FileName_X,VecDef_X);
  strcpy(FileName_Y,Prefix);
  strcat(FileName_Y,VecDef_Y);
  strcpy(FileName_Z,Prefix);
  strcat(FileName_Z,VecDef_Z);
  
  Temp3DField.Write(FileName_X,FileName_Y,FileName_Z,this->SourceFile);
}




///Estimate a reference weight for a given sigma in the smoothing kernel
float LIDM::EstimRefWeight(float sigmaloc,LightFFTconvolver3D * ConvolverLoc,int NeedForInit){
  float max,refWght,tmpFl,temp;
  int d,x,y,z;
  int i;
  
  //1) Initiate the mappings
  if (NeedForInit==1){
    CptMappingFromVeloField_IniIdMap(0,&this->VelocityField,&this->ForwardMapping);
    CptMappingFromVeloField_IniIdMap(this->VelocityField.NT-1,&this->VelocityField,&this->BackwardMapping);
  }
  
  //2) do the process equivalent to an iteration (without the smoothing)
  if (NeedForInit==1){
    Cpt_JacobianDeterminant(&BackwardMapping,&DetJacobians,0);
    Project3Dimage(&this->ImTemplate,&this->ForwardMapping,&this->J0,0);
    Project3DImageUsingAffineTransfoAndTimeDepVF(this->Template2TargetCoord,&this->ImTarget,&this->BackwardMapping,&this->J1,0);
  }
  Cpt_Grad_ScalarField(&this->J0,&this->GradJ);
  
  //3) compute typical energy gradient2d image
  for (i=0;i<3;i++) for (z = 0; z < this->NZ; z++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      if ((z<this->Margin)||(z>this->NZ-this->Margin-1)||(y<this->Margin)||(y>this->NY-this->Margin-1)||(x<this->Margin)||(x>this->NX-this->Margin-1))
        temp=0;
      else 
        temp=(this->J0.G(x,y,z) - this->J1.G(x,y,z)) * this->DetJacobians.G(x,y,z) * this->GradJ.G(i,x,y,z);
      this->GradJ.P(static_cast<float>(temp),i,x,y,z);
  }
  
  if (this->NZ==1){//...2d image
    for (i=0;i<3;i++) for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
        if ((y<this->Margin)||(y>this->NY-this->Margin-1)||(x<this->Margin)||(x>this->NX-this->Margin-1))
          temp=0;
        else 
          temp=(this->J0.G(x,y,0) - this->J1.G(x,y,0)) * this->DetJacobians.G(x,y,0) * this->GradJ.G(i,x,y,0);
        this->GradJ.P(static_cast<float>(temp),i,x,y,0);
    }
  }
  
  //4) change the convolver
  ConvolverLoc->ChangeKernel(1,sigmaloc/this->x_mm,sigmaloc/this->y_mm,sigmaloc/this->z_mm);
  
  //5) convolution
  ConvolverLoc->Convolution(&this->GradJ);
  
  //6 compute the maximum gradient of energy
  max=0;
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    tmpFl=this->GradJ.G(0,x,y,z)*this->GradJ.G(0,x,y,z)+this->GradJ.G(1,x,y,z)*this->GradJ.G(1,x,y,z)+this->GradJ.G(2,x,y,z)*this->GradJ.G(2,x,y,z);
    if (max<tmpFl) max=tmpFl;
  }
  
  //7) Return the update scale
  refWght=1/sqrt(max);
  
  return refWght;
}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                      RUN FUNCTIONS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




///Function to solve the registration using the gradient descent algorithm of Beg 05
void LIDM::Run_Default(void){
  int IterationStopper;
  int IterationNb;
  int TimeSubdiv;
  float SqrtSSD;
  int i;
  
  //1) INITIALISATION
  
  //1.1) Pre-treatment of the inuput images (grey level alignment + margins)
  this->ReadAndTreatInputImages();
  
  //1.2) Allocations of the scalar and vector fields + definition of global parameters
  this->AllocateAllVariables();
  
  
  //2) GRADIENT DESCENT
  if (this->iteration_nb==0) IterationStopper=1;
  else IterationStopper=0;
  IterationNb=0;
  
  
  while (IterationStopper==0){
    cout << "Iteration Number " << IterationNb+1 << " / " << this->iteration_nb << " ";// "\n";
    
    //2.1) update the velocity fields (part 1: multiply the VF by the VFpenalizer)
    this->UpdateVelocityField_part1();
    
    //2.2) compute the forward and backward mappings on space
    CptMappingFromVeloField_IniIdMap(0,&this->VelocityField,&this->ForwardMapping);
    CptMappingFromVeloField_IniIdMap(this->VelocityField.NT-1,&this->VelocityField,&this->BackwardMapping);
    
    //2.3) LOOP ON THE TIME SUBDIVISIONS AND THE CHANNELS
    for (TimeSubdiv=0;TimeSubdiv<this->NbTimeSubdiv;TimeSubdiv++){//LOOP ON THE TIME SUBDIVISIONS
      
      //2.3.1) compute the determinant of the jacobian of the transformation
      Cpt_JacobianDeterminant(&BackwardMapping,&DetJacobians,TimeSubdiv);
      
      //2.3.2) compute the temporary image transformed using the forward mapping from time 0 -> J0
      Project3Dimage(&this->ImTemplate,&this->ForwardMapping,&this->J0,TimeSubdiv);
      
      //2.3.2.bis) project also the partition of unity if asked
      if (this->MovingPOU==1) this->ProjectPOU(TimeSubdiv);
      
      //2.3.3) compute the temporary image transformed using the backward mapping from time 1 -> J1
      Project3DImageUsingAffineTransfoAndTimeDepVF(Template2TargetCoord,&this->ImTarget,&this->BackwardMapping,&this->J1,TimeSubdiv);
      
      //2.3.4) compute gradient of J
      Cpt_Grad_ScalarField(&this->J0,&this->GradJ);
      
      //2.3.5) compute the gradient of energy
      this->ComputeEnergyGradient(TimeSubdiv);
      
      //2.3.7) Measure the convergence of the similarity between the source and target image
      if (ShowSSD==1) if (TimeSubdiv==this->NbTimeSubdiv-1){
        SqrtSSD=CalcSqrtSumOfSquaredDif(&this->J0,&this->J1);
        cout << "Sqrt SSD = " << SqrtSSD  << "\n";
      }
    }
    
    //2.4) update the velocity fields (part 2: add the energy gradients)
    this->UpdateVelocityField_part2(IterationNb);
    
    //2.5) end of the convergence...
    //...controled by the number of iterations
    IterationNb++;
    if (IterationNb>=this->iteration_nb) IterationStopper=1;
  }
  
  //3) SAVE THE RESULTS
  this->SaveResultGradientDescent();
  
  cout << endl;
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                        MAIN RUN FUNCTION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///run function
void LIDM::Run(void)
{
    this->Run_Default();
}

