/*=========================================================================
 
 Author: Laurent Risser
 
 Disclaimer: This software has been developed for research purposes only, and hence should 
 not be used as a diagnostic tool. In no event shall the authors or distributors
 be liable to any direct, indirect, special, incidental, or consequential 
 damages arising of the use of this software, its documentation, or any 
 derivatives thereof, even if the authors have been advised of the possibility 
 of such damage. 
 
 =========================================================================*/

#include <Demons.h>

///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                   CONSTRUCTOR AND DESTRUCTOR
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
LargeDefDemons::LargeDefDemons(void){
  
  //default parameters
  iteration_nb=10;
  NbKernelScales=0;
  sigma1=2;
  sigma2=1;
  sigma3=1;
  sigma4=1;
  sigma5=1;
  sigma6=1;
  sigma7=1;
  w1=1;
  w2=0;
  w3=0;
  w4=0;
  w5=0;
  w6=0;
  w7=0;
  sigmaDiff=0;
  MaxUpdateAllowed=1;
  Margin=0;
  GreyLevAlign=0;
  IndicatorMI=0;
  GLA_Padding_Src=-1.;
  GLA_Padding_Trg=-1.;
  lambdaX=1;
  strcpy(PrefixInputs,"Null");
  strcpy(PrefixOutputs,"Outputs");
  strcpy(SourceFiles[0],"Null");
  strcpy(TargetFiles[0],"Null");
  NbChannels=0;
  strcpy(MaskFile,"Null");
  MaskDefined=0;
  MaskLabel=1;
  World_Target2Template[0][0]=1; World_Target2Template[0][1]=0;   World_Target2Template[0][2]=0;    World_Target2Template[0][3]=0;   
  World_Target2Template[1][0]=0; World_Target2Template[1][1]=1;   World_Target2Template[1][2]=0;    World_Target2Template[1][3]=0;   
  World_Target2Template[2][0]=0; World_Target2Template[2][1]=0;   World_Target2Template[2][2]=1;    World_Target2Template[2][3]=0;   
  World_Target2Template[3][0]=0; World_Target2Template[3][1]=0;   World_Target2Template[3][2]=0;    World_Target2Template[3][3]=1;
  BoundaMargin=30;
  x_mm=1;
  y_mm=1;
  y_mm=1;
  ExtendTrgImag_LowerX=0;
  ExtendTrgImag_UpperX=0;
  ExtendTrgImag_LowerY=0;
  ExtendTrgImag_UpperY=0;
  ExtendTrgImag_LowerZ=0;
  ExtendTrgImag_UpperZ=0;
  UnderSampleTrgFactor=1;
  strcpy(IniDispFieldX,"Null");
  strcpy(IniDispFieldY,"Null");
  strcpy(IniDispFieldZ,"Null");
  IniDispFieldDefined=0;
  IndicatorSaveMemory=0;
  FinalDefVec=0;
}

LargeDefDemons::~LargeDefDemons(void){}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                        SUB-FUNCTIONS TO PERFORM THE REGISTRATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///initiate the gradient descent
void LargeDefDemons::ReadAndTreatInputImages(void){
  int x, y, z;
  int DistClosestEdge;
  int i,j;
  double mean1,mean2,std1,std2;
  float tempQuat[4][4];
  float TmpFl;
  ScalarField TempImTarget;
  char FileTreatedTarget[256];
  char DotNii[256];
  
  //1) CREATE THE SOURCE AND TARGET IMAGES 3D FOR THE CALCULATIONS
  //    -->  ImTemplate[c].G(x,y,z) = gray level at (x,y,z)
  //    -->  ImTarget[c].G(x,y,z)  = gray level at (x,y,z)
  this->ImTemplate = new ScalarField [this->NbChannels];
  this->ImTarget = new ScalarField [this->NbChannels];
  
  //2) READ INPUTS
  //2.1.a) read source files
  for (i=0;i<this->NbChannels;i++){
    this->ImTemplate[i].Read(this->SourceFiles[i]);
  }
    
  //2.1.b) read target files
  if ((fabs(this->UnderSampleTrgFactor-1)<=0.01)&&(ExtendTrgImag_LowerX<=0)&&(ExtendTrgImag_UpperX<=0)&&(ExtendTrgImag_LowerY<=0)&&(ExtendTrgImag_UpperY<=0)&&(ExtendTrgImag_LowerZ<=0)&&(ExtendTrgImag_UpperZ<=0)){
    
    // the images do not required any downsampling or domain expansion
    for (i=0;i<this->NbChannels;i++){
      this->ImTarget[i].Read(this->TargetFiles[i]);
    }
  }
  else{ // there is at leat a special treatment to do with the target images
    //name of the treated target file
    strcpy(FileTreatedTarget,this->PrefixOutputs);
    strcpy(DotNii,"_tmp.nii");
    strcat(FileTreatedTarget,DotNii);
    
    cout << "After treatments (resampling and/or expansion): " <<  this->TargetFiles[0] << " becomes " << FileTreatedTarget << endl;
    
    //treat all the channels
    for (i=0;i<this->NbChannels;i++){
      
      //expend the image
      if (!((ExtendTrgImag_LowerX<=0)&&(ExtendTrgImag_UpperX<=0)&&(ExtendTrgImag_LowerY<=0)&&(ExtendTrgImag_UpperY<=0)&&(ExtendTrgImag_LowerZ<=0)&&(ExtendTrgImag_UpperZ<=0))){
        TempImTarget.ReadAndExpend(this->TargetFiles[i],ExtendTrgImag_LowerX,ExtendTrgImag_UpperX,ExtendTrgImag_LowerY,ExtendTrgImag_UpperY,ExtendTrgImag_LowerZ,ExtendTrgImag_UpperZ);
        TempImTarget.Write(FileTreatedTarget);
        
        strcpy(this->TargetFiles[i],FileTreatedTarget);  // a bit dirty but useful as the results are saved in the target coordinate space
        if (fabs(this->UnderSampleTrgFactor-1)<=0.01) this->ImTarget[i].Read(this->TargetFiles[i]);
      }
      
      //undersample the image
      if (fabs(this->UnderSampleTrgFactor-1)>0.01){
        TempImTarget.Read_and_Undersample(this->TargetFiles[i],this->UnderSampleTrgFactor);
        TempImTarget.Write(FileTreatedTarget);
        
        strcpy(this->TargetFiles[i],FileTreatedTarget);  // a bit dirty but useful as the results are saved in the target coordinate space
        this->ImTarget[i].Read(this->TargetFiles[i]);
      }
    }
  }
  
  //2.1.a) allocate the memory for the deformed template
  this->DeformedTemplate.Read(this->TargetFiles[0]);
  
  //2.2) check whether  3D or 2D images are opened
  if (this->ImTemplate[0].NT>1) cout << "Source image  depends on time!!!";
  if (this->ImTarget[0].NT>1) cout << "Target image depends on time!!!";
  
  //2.3) variables containing the size of the target image (and all the other images, except the source image)
  this->NX=this->ImTarget[0].NX;
  this->NY=this->ImTarget[0].NY;
  this->NZ=this->ImTarget[0].NZ;
  this->NT=1;
  
  //2.4) variables containing the size of the source image
  this->NXs=this->ImTemplate[0].NX;
  this->NYs=this->ImTemplate[0].NY;
  this->NZs=this->ImTemplate[0].NZ;
  
  cout << "Image size: " << this->NX <<  "*"  <<  this->NY  <<  "*"  << this->NZ  << " (source: " << this->NXs <<  "*"  <<  this->NYs  <<  "*"  << this->NZs  << ")\n";
  
  //2.5) compute the quaternion to convert target coordinates into template coordinates
  mult_quat4t4mat_quat4t4mat(World_Target2Template,this->ImTarget[0].Image2World,tempQuat);
  mult_quat4t4mat_quat4t4mat(this->ImTemplate[0].World2Image,tempQuat,Target2TemplateCoord);
  
  cout << endl;
  
  if (IniDispFieldDefined==1){
    cout << "Target to template encoded in world coordinate by " << this->IniDispFieldX << ", " <<  this->IniDispFieldY << " and " <<  this->IniDispFieldZ << endl;
  }
  else{
    cout << "Target to template in voxel coordinates:" << endl;
    for (i=0;i<4;i++){
      for (j=0;j<4;j++){
        cout << Target2TemplateCoord[i][j] << " ";
      }
      cout << endl;
    }
  }
  
  
  //2.6 compute the voxels size in mm
  this->x_mm=sqrt(this->ImTarget[0].Image2World[0][0]*this->ImTarget[0].Image2World[0][0]+this->ImTarget[0].Image2World[0][1]*this->ImTarget[0].Image2World[0][1]+this->ImTarget[0].Image2World[0][2]*this->ImTarget[0].Image2World[0][2]);
  this->y_mm=sqrt(this->ImTarget[0].Image2World[1][0]*this->ImTarget[0].Image2World[1][0]+this->ImTarget[0].Image2World[1][1]*this->ImTarget[0].Image2World[1][1]+this->ImTarget[0].Image2World[1][2]*this->ImTarget[0].Image2World[1][2]);
  this->z_mm=sqrt(this->ImTarget[0].Image2World[2][0]*this->ImTarget[0].Image2World[2][0]+this->ImTarget[0].Image2World[2][1]*this->ImTarget[0].Image2World[2][1]+this->ImTarget[0].Image2World[2][2]*this->ImTarget[0].Image2World[2][2]);
  
  cout << endl;
  cout << "Target image resolution: " << this->x_mm << " "  << this->y_mm << " "  << this->z_mm << endl;

  //3) TAKING INTO ACCOUNT THE MARGINS:  !!! this->Margin  !!!
  
  // to do ?
  
  
  //4) COMPUTE THE MAPPING (in from the target c.s. to the source c.s.)
  if (IniDispFieldDefined==1) this->ReadAndConvertMapping();
  
  
  
  //5) LINEAR ALIGNMENT OF THE GREY LEVELS OF ImTarget ON THOSE OF ImTemplate
  int NbVoxelsOK;
  
  if (GreyLevAlign!=0) for (i=0;i<this->NbChannels;i++){
    //compute mean and std dev of the source and target images
    mean1=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZs; z++)  for (y = 0; y < this->NYs; y++) for (x = 0; x < this->NXs; x++) if (this->ImTemplate[i].G(x,y,z)>GLA_Padding_Src){
      mean1+=(double)this->ImTemplate[i].G(x,y,z);
      NbVoxelsOK++;
    }
    mean1/=(double)(NbVoxelsOK);
    
    mean2=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTarget[i].G(x,y,z)>GLA_Padding_Trg){
      mean2+=(double)this->ImTarget[i].G(x,y,z);
      NbVoxelsOK++;
    }
    mean2/=(double)(NbVoxelsOK);
    
    std1=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZs; z++)  for (y = 0; y < this->NYs; y++) for (x = 0; x < this->NXs; x++) if (this->ImTemplate[i].G(x,y,z)>GLA_Padding_Src){
      std1+=pow((double)this->ImTemplate[i].G(x,y,z)-mean1,2.);
      NbVoxelsOK++;
    }
    std1/=(double)(NbVoxelsOK);
    std1=sqrt(std1);
    
    std2=0.;
    NbVoxelsOK=0;
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (this->ImTarget[i].G(x,y,z)>GLA_Padding_Trg){
      std2+=pow((double)this->ImTarget[i].G(x,y,z)-mean2,2.);
      NbVoxelsOK++;
    }
    std2/=(double)(NbVoxelsOK);
    std2=sqrt(std2);
    
    cout << "Template: mean=" << mean1 << ", stddev=" << std1 << ".    Target: mean=" << mean2 << ", stddev=" << std2 << "\n";
    
    
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
      this->ImTarget[i].P((this->ImTarget[i].G(x,y,z)-(float)mean2)*(((float)std1)/((float)std2))+(float)mean1,x,y,z);
    
    
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++)
      if ((this->ImTarget[i].G(x,y,z)<(GLA_Padding_Trg-(float)mean2)*(((float)std1)/((float)std2))+(float)mean1)||(this->ImTarget[i].G(x,y,z)<GLA_Padding_Src))
        this->ImTarget[i].P(0.,x,y,z);
  }
  

  
  //6) READ THE MASK AND PROJECT IT
  //6.1) All kinds of masks
  if (this->MaskDefined>0){ 
    this->Mask.Read(this->MaskFile);
    
    if ((ImTemplate[0].NX!=this->Mask.NX)) cout << "The source image and the mask do not have the same size!!!";
    if ((ImTemplate[0].NY!=this->Mask.NY)) cout << "The source image and the mask do not have the same size!!!";
    if ((ImTemplate[0].NZ!=this->Mask.NZ)) cout << "The source image and the mask do not have the same size!!!";
    
    
    this->ProjMask.CreateVoidField(this->NX,this->NY,this->NZ);
    this->VelocityField.CreateVoidField(this->NX,this->NY,this->NZ);
    
    if (IniDispFieldDefined==1) ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->Mask,&this->ProjMask);
    if (IniDispFieldDefined!=1) ProjectImageUsingAffineTransfoAndSteadyVeloField(Target2TemplateCoord,&this->VelocityField,&this->Mask,&this->ProjMask);
  }
  
  
  //6.2) mask to perform sliding motion registration
  if ((this->MaskDefined==2)&&(this->NZ==1)){
    cout << "Sliding motion encoded for 3D images only" << endl;
    this->MaskDefined=1;
  } 
  
  if (this->MaskDefined==2){
    //6.2.1) distance to the boundary if discontinuous deformations are considered
    
    this->NearestBoundary.CreateVoidField(this->NX,this->NY,this->NZ);
    this->TempSF.CreateVoidField(this->NX,this->NY,this->NZ);
    
    //6.2.2) compute the distance map
    cout << "Compute the distance map" << endl;
    
    Cpt_NearestBoundary(&this->ProjMask,&this->NearestBoundary,&this->TempSF);
   
   //6.2.3) smooth the distance map and change TempSF
   LightFFTconvolver3D LocConvolver;
   LocConvolver.InitiateConvolver(this->TempSF.NX,this->TempSF.NY,this->TempSF.NZ,1,2,2,2);
   
   //LocConvolver.Convolution(&this->TempSF);  // to improve: smoothing with no influence of the boundaries
   LocConvolver.Convolution(&this->NearestBoundary); //remark here that some problem may arise if BoundaMargin is too large: Imagine boundaries = "|       |" -> NearestBoundary = "|1 2 3 -3 -2 -1|"  -> smoothed NearestBoundary = "|1 2 1 -1 -2 -1|"
   
   //6.2.4) Define the mask in which the reorientations will performed
   for (z = 0; z < this->TempSF.NZ; z++)  for (y = 0; y < this->TempSF.NY; y++) for (x = 0; x < this->TempSF.NX; x++) {
      if ((fabs(this->ProjMask.G(x,y,z)-this->MaskLabel)>0.001)||(this->TempSF.G(x,y,z)>this->BoundaMargin)){
	this->TempSF.P(0,x,y,z);
	this->NearestBoundary.P(0,0,x,y,z);
	this->NearestBoundary.P(0,1,x,y,z);
	this->NearestBoundary.P(0,2,x,y,z);
      }
      else{
        this->TempSF.P(1,x,y,z);
      }
      if (fabs(this->ProjMask.G(x,y,z)-this->MaskLabel)>0.001){
	this->ProjMask.P(0,x,y,z);
      }
      else{
        this->ProjMask.P(1,x,y,z);
      }
    }
    
    //At this point: 
    // -> TempSF is equal to 1 where reorientation should be considered
    // -> NearestBoundary points to the nearest boundaries of ProjMask in the whole domain
    // -> ProjMask is equal to 1 in the ROI where the computations are made and 0 otherwise
    
    //this->ProjMask.Write("ProjMask.nii",this->TargetFiles[0]);
    //this->TempSF.Write("TempSF.nii",this->TargetFiles[0]);
    //this->NearestBoundary.Write("DMX.nii","DMY.nii","DMZ.nii",this->TargetFiles[0]);
    
    cout << "Distance map computed" << endl;
    
  }
  
  
  //7) LOAD THE DEFORMATION FIELD
  //    -->  VelocityField.G(0,x,y,z,i)= direction ex of the vector at (x,y,z)
  //    -->  VelocityField.G(1,x,y,z,i)= direction ey of the vector at (x,y,z)
  //    -->  VelocityField.G(2,x,y,z,i)= direction ez of the vector at (x,y,z)
  if (strcmp(PrefixInputs,"Null")!=0){
      this->LoadVelocityField(PrefixInputs);
  }
  else{
      this->VelocityField.CreateVoidField(this->NX,this->NY,this->NZ);
  }
  
  
  
  
}

///allocate all variables used for the gradient descent
void LargeDefDemons::AllocateAllVariables(void){
  int i;
  
  //contains the gradient of Energy. Is directely used to update the Displacement Field
  this->GradE.CreateVoidField(this->NX,this->NY,this->NZ);
  
  if (this->NbChannels>1) this->TempVF.CreateVoidField(this->NX,this->NY,this->NZ);  // usefull to compute the gradients when there are several channels
  
  
  //manage the mutual information
  if (IndicatorMI==1){
    if (IniDispFieldDefined==1){
        ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
    }
    else{
        ProjectImageUsingAffineTransfoAndSteadyVeloField(Target2TemplateCoord,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
    }
    
    //this->Mask.Read("Mask.nii");//IF I WANT TO TAKE INTO ACCOUNT A MASK
    //this->NorMutInfMan.Initiate(&this->DeformedTemplate,&this->ImTarget[0],&this->Mask,50,50,3);//IF I WANT TO TAKE INTO ACCOUNT A DEFORMING MASK
    
    this->NorMutInfMan.Initiate(&this->DeformedTemplate,&this->ImTarget[0],50,50,3);//TO COMMENT IF I WANT TO TAKE INTO ACCOUNT A DEFORMING MASK
  }
  
  //initiate the FFT convolvers
  if (sigma1>0.005){
    //allocate the memory
    this->FFTconvolver_fluid.InitiateConvolver(this->NX,this->NY,this->NZ);
    
    //esimate the weights
    if (fabs(this->w1)<0.0001){
      for (i=0;i<this->NbKernelScales;i++){
        if (i==0) { this->w1=this->EstimRefUpdateScale(sigma1,&this->FFTconvolver_fluid); cout << "Weight for scale 1 is " <<  100     << endl; }
        if (i==1) { this->w2=this->EstimRefUpdateScale(sigma2,&this->FFTconvolver_fluid); cout << "Weight for scale 2 is " <<  100*this->w2/this->w1 << endl; }
        if (i==2) { this->w3=this->EstimRefUpdateScale(sigma3,&this->FFTconvolver_fluid); cout << "Weight for scale 3 is " <<  100*this->w3/this->w1 << endl; } 
        if (i==3) { this->w4=this->EstimRefUpdateScale(sigma4,&this->FFTconvolver_fluid); cout << "Weight for scale 4 is " <<  100*this->w4/this->w1 << endl; }
        if (i==4) { this->w5=this->EstimRefUpdateScale(sigma5,&this->FFTconvolver_fluid); cout << "Weight for scale 5 is " <<  100*this->w5/this->w1 << endl; }
        if (i==5) { this->w6=this->EstimRefUpdateScale(sigma6,&this->FFTconvolver_fluid); cout << "Weight for scale 6 is " <<  100*this->w6/this->w1 << endl; }
        if (i==6) { this->w7=this->EstimRefUpdateScale(sigma7,&this->FFTconvolver_fluid); cout << "Weight for scale 7 is " <<  100*this->w7/this->w1 << endl; }
      }
    }
    
    w7/=w1;    w6/=w1;    w5/=w1;    w4/=w1;    w3/=w1;    w2/=w1;    w1=1;
    
    //set the kernels
    this->FFTconvolver_fluid.InitiateConvolver(this->NX,this->NY,this->NZ,
                                               w1,sigma1/this->x_mm,sigma1/this->y_mm,sigma1/this->z_mm,
                                               w2,sigma2/this->x_mm,sigma2/this->y_mm,sigma2/this->z_mm,
                                               w3,sigma3/this->x_mm,sigma3/this->y_mm,sigma3/this->z_mm,
                                               w4,sigma4/this->x_mm,sigma4/this->y_mm,sigma4/this->z_mm,
                                               w5,sigma5/this->x_mm,sigma5/this->y_mm,sigma5/this->z_mm,
                                               w6,sigma6/this->x_mm,sigma6/this->y_mm,sigma6/this->z_mm,
                                               w7,sigma7/this->x_mm,sigma7/this->y_mm,sigma7/this->z_mm);
  }
  
  if (sigmaDiff>0.005) this->FFTconvolver_Diff.InitiateConvolver(this->NX,this->NY,this->NZ,1,sigmaDiff/this->x_mm,sigmaDiff/this->y_mm,sigmaDiff/this->z_mm);
  
}



///save the result of the gradient descent
void LargeDefDemons::ReadAndConvertMapping(void){
  float flX,flY,flZ;
  float WTX,WTY,WTZ;
  float WSX,WSY,WSZ;
  float ISX,ISY,ISZ;
  int x,y,z;
  
  //read the mapping in the world coordinate system   (IniDispField:  points from target to source in millimeters)
  this->IniDispField.Read_and_Interpolate(this->IniDispFieldX,this->IniDispFieldY,this->IniDispFieldZ,this->NX,this->NY,this->NZ,0);
  
  //convert the mapping in voxel coordonate system (points from target to source in voxels)
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    
    flX=static_cast<float>(x); flY=static_cast<float>(y); flZ=static_cast<float>(z);
    
    WTX=flX*this->ImTarget[0].Image2World[0][0]+flY*this->ImTarget[0].Image2World[0][1]+flZ*this->ImTarget[0].Image2World[0][2]+this->ImTarget[0].Image2World[0][3];
    WTY=flX*this->ImTarget[0].Image2World[1][0]+flY*this->ImTarget[0].Image2World[1][1]+flZ*this->ImTarget[0].Image2World[1][2]+this->ImTarget[0].Image2World[1][3];
    WTZ=flX*this->ImTarget[0].Image2World[2][0]+flY*this->ImTarget[0].Image2World[2][1]+flZ*this->ImTarget[0].Image2World[2][2]+this->ImTarget[0].Image2World[2][3];
    
    WSX=WTX+this->IniDispField.G(0,x,y,z);
    WSY=WTY+this->IniDispField.G(1,x,y,z);
    WSZ=WTZ+this->IniDispField.G(2,x,y,z);
    
    ISX=WSX*this->ImTemplate[0].World2Image[0][0]+WSY*this->ImTemplate[0].World2Image[0][1]+WSZ*this->ImTemplate[0].World2Image[0][2]+this->ImTemplate[0].World2Image[0][3];
    ISY=WSX*this->ImTemplate[0].World2Image[1][0]+WSY*this->ImTemplate[0].World2Image[1][1]+WSZ*this->ImTemplate[0].World2Image[1][2]+this->ImTemplate[0].World2Image[1][3];
    ISZ=WSX*this->ImTemplate[0].World2Image[2][0]+WSY*this->ImTemplate[0].World2Image[2][1]+WSZ*this->ImTemplate[0].World2Image[2][2]+this->ImTemplate[0].World2Image[2][3];
    
    this->IniDispField.P(ISX,0,x,y,z);
    this->IniDispField.P(ISY,1,x,y,z);
    this->IniDispField.P(ISZ,2,x,y,z);
  }
  
  //IniDispField.Write("IniDispFieldDemVoxX.nii","IniDispFieldDemVoxY.nii","IniDispFieldDemVoxZ.nii",this->TargetFiles[0]);
}


///save the result of the gradient descent
void LargeDefDemons::SaveResultGradientDescent(void){
  VectorField RealDisplacementField;
  int i,j;
  float TempQuat_T[4][4];
  float TempQuat_S[4][4];
  
  //1) free some memory
  if (IndicatorSaveMemory==1){
    cout << "DEALOCATE SOME MEMORY" << endl;
    for (i=1;i<this->NbChannels;i++){
      this->ImTarget[i].SlashFieldSize();
      this->ImTemplate[i].SlashFieldSize();
    }
    
    if (this->MaskDefined==1) this->Mask.SlashFieldSize();
    
    this->GradE.SlashFieldSize();
    
    if (this->NbChannels>1) this->TempVF.SlashFieldSize();
    
    this->FFTconvolver_fluid.InitiateConvolver(1,1,1);
  }
  
  //2) save the results...

  //save the velocity field
  this->SaveVelocityField(&this->VelocityField,this->PrefixOutputs);
  
  //save the final deformation
  if (IndicatorSaveMemory!=1) this->SaveDeformations(this->PrefixOutputs);
  else cout << "Final deformation not saved" << endl;
  
  //dealocate more memory
  if (IndicatorSaveMemory==1){
    cout << "DEALOCATE MORE MEMORY" << endl;
    
    for (i=0;i<4;i++) for (j=0;j<4;j++) TempQuat_T[i][j]=this->ImTarget[0].Image2World[i][j];
    for (i=0;i<4;i++) for (j=0;j<4;j++) TempQuat_S[i][j]=this->ImTemplate[0].Image2World[i][j];
    
    this->ImTemplate[0].SlashFieldSize();
    this->ImTarget[0].SlashFieldSize();
    this->DeformedTemplate.SlashFieldSize();
    
    for (i=0;i<4;i++) for (j=0;j<4;j++) this->ImTarget[0].Image2World[i][j]=TempQuat_T[i][j];
    for (i=0;i<4;i++) for (j=0;j<4;j++) this->ImTemplate[0].Image2World[i][j]=TempQuat_S[i][j];
  }
  
  
  //save the inverse displacement field
  RealDisplacementField.CreateVoidField(this->NX,this->NY,this->NZ);
  CptInvDefFromSteadyVeloField(&this->VelocityField,&RealDisplacementField,7);
  
  if (IndicatorSaveMemory==1){
    cout << "DEALOCATE MORE MEMORY" << endl;
    this->VelocityField.SlashFieldSize();
  }
  
  this->SaveInvTotalDisplacement(&RealDisplacementField,this->PrefixOutputs);
    
  //save the displacement field
  if (this->FinalDefVec==1){
     if (this->IndicatorSaveMemory==1){
          cout << "Can't save the displacement field if the -SaveMemory option is turned on (inverse disp. field is however saved)" << endl;
        }
     else if (this->IniDispFieldDefined==1){
         cout << "Can't save the displacement field if -IniDispF option is used  (inverse disp. field is however saved)" << endl;
       }
     else{
          CptDefFromSteadyVeloField(&this->VelocityField,&RealDisplacementField,7);
          this->SaveTotalDisplacement(&RealDisplacementField,this->PrefixOutputs);
        }
    }

}


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                          FUNCTIONS TO SAVE AND LOAD THE VARIOUS STRUCTURES
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++




///load the velocity fields
void LargeDefDemons::LoadVelocityField(char Prefix[256]){
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


///save the inverse displacement field
void LargeDefDemons::SaveInvTotalDisplacement(VectorField * LocDisplacementField,char Prefix[256]){
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
  strcpy(InvDisplacementField_X,"_DispField_Trg2Src_X.nii");
  strcat(FileNameX,InvDisplacementField_X);
  strcpy(FileNameY,Prefix);
  strcpy(InvDisplacementField_Y,"_DispField_Trg2Src_Y.nii");
  strcat(FileNameY,InvDisplacementField_Y);
  strcpy(FileNameZ,Prefix);
  strcpy(InvDisplacementField_Z,"_DispField_Trg2Src_Z.nii");
  strcat(FileNameZ,InvDisplacementField_Z);
  
  if (IniDispFieldDefined!=1){  //+++++++ Case A: Only an initial affine alignment is considered between the images +++++++
    
    //compute the inverse of the total displacement field in mm (world coordinate)
    mult_quat4t4mat_quat4t4mat(World_Target2Template,this->ImTarget[0].Image2World,tempQuat);
    
    //compute the mapping
    for (z = 0; z < LocDisplacementField->NZ; z++)  for (y = 0; y < LocDisplacementField->NY; y++) for (x = 0; x < LocDisplacementField->NX; x++){
      trgX=x*this->ImTarget[0].Image2World[0][0]+y*this->ImTarget[0].Image2World[0][1]+z*this->ImTarget[0].Image2World[0][2]+this->ImTarget[0].Image2World[0][3];
      trgY=x*this->ImTarget[0].Image2World[1][0]+y*this->ImTarget[0].Image2World[1][1]+z*this->ImTarget[0].Image2World[1][2]+this->ImTarget[0].Image2World[1][3];
      trgZ=x*this->ImTarget[0].Image2World[2][0]+y*this->ImTarget[0].Image2World[2][1]+z*this->ImTarget[0].Image2World[2][2]+this->ImTarget[0].Image2World[2][3];
      
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
      
      trgX=x*this->ImTarget[0].Image2World[0][0]+y*this->ImTarget[0].Image2World[0][1]+z*this->ImTarget[0].Image2World[0][2]+this->ImTarget[0].Image2World[0][3];
      trgY=x*this->ImTarget[0].Image2World[1][0]+y*this->ImTarget[0].Image2World[1][1]+z*this->ImTarget[0].Image2World[1][2]+this->ImTarget[0].Image2World[1][3];
      trgZ=x*this->ImTarget[0].Image2World[2][0]+y*this->ImTarget[0].Image2World[2][1]+z*this->ImTarget[0].Image2World[2][2]+this->ImTarget[0].Image2World[2][3];
      
      tmpX=flX+LocDisplacementField->G(0,x,y,z);
      tmpY=flY+LocDisplacementField->G(1,x,y,z);
      tmpZ=flZ+LocDisplacementField->G(2,x,y,z);
      
      tmpX2=this->IniDispField.G(0,tmpX,tmpY,tmpZ);
      tmpY2=this->IniDispField.G(1,tmpX,tmpY,tmpZ);
      tmpZ2=this->IniDispField.G(2,tmpX,tmpY,tmpZ);
      
      srcX=tmpX2*this->ImTemplate[0].Image2World[0][0]+tmpY2*this->ImTemplate[0].Image2World[0][1]+tmpZ2*this->ImTemplate[0].Image2World[0][2]+this->ImTemplate[0].Image2World[0][3];
      srcY=tmpX2*this->ImTemplate[0].Image2World[1][0]+tmpY2*this->ImTemplate[0].Image2World[1][1]+tmpZ2*this->ImTemplate[0].Image2World[1][2]+this->ImTemplate[0].Image2World[1][3];
      srcZ=tmpX2*this->ImTemplate[0].Image2World[2][0]+tmpY2*this->ImTemplate[0].Image2World[2][1]+tmpZ2*this->ImTemplate[0].Image2World[2][2]+this->ImTemplate[0].Image2World[2][3];
      
      LocDisplacementField->P(srcX-trgX,0,x,y,z);
      LocDisplacementField->P(srcY-trgY,1,x,y,z);
      LocDisplacementField->P(srcZ-trgZ,2,x,y,z);
    }
    
  }
  
  
  //save the deformation field
  LocDisplacementField->Write(FileNameX,FileNameY,FileNameZ,TargetFiles[0]);
}




///save the inverse displacement field
void LargeDefDemons::SaveTotalDisplacement(VectorField * LocDisplacementField,char Prefix[256]){
  char FileNameX[256];
  char FileNameY[256];
  char FileNameZ[256];
  char DisplacementField_X[256];
  char DisplacementField_Y[256];
  char DisplacementField_Z[256];
  float srcX,srcY,srcZ;
  float trgX,trgY,trgZ;
  float tmpX,tmpY,tmpZ;
  int x,y,z;
  
  //intialisation
  strcpy(FileNameX,Prefix);
  strcpy(DisplacementField_X,"_DispField_Src2Trg_X.nii");
  strcat(FileNameX,DisplacementField_X);
  strcpy(FileNameY,Prefix);
  strcpy(DisplacementField_Y,"_DispField_Src2Trg_Y.nii");
  strcat(FileNameY,DisplacementField_Y);
  strcpy(FileNameZ,Prefix);
  
  strcpy(DisplacementField_Z,"_DispField_Src2Trg_Z.nii");
  strcat(FileNameZ,DisplacementField_Z);
  
  if (IniDispFieldDefined==1){
    cout << "Can't save the displacement field if -IniDispF option is used" << endl;
    return;
    }
  
  //compute the mapping
  for (z = 0; z < LocDisplacementField->NZ; z++)  for (y = 0; y < LocDisplacementField->NY; y++) for (x = 0; x < LocDisplacementField->NX; x++){
    tmpX=x+LocDisplacementField->G(0,x,y,z);
    tmpY=y+LocDisplacementField->G(1,x,y,z);
    tmpZ=z+LocDisplacementField->G(2,x,y,z);

    trgX=tmpX*this->ImTarget[0].Image2World[0][0]+tmpY*this->ImTarget[0].Image2World[0][1]+tmpZ*this->ImTarget[0].Image2World[0][2]+this->ImTarget[0].Image2World[0][3];
    trgY=tmpX*this->ImTarget[0].Image2World[1][0]+tmpY*this->ImTarget[0].Image2World[1][1]+tmpZ*this->ImTarget[0].Image2World[1][2]+this->ImTarget[0].Image2World[1][3];
    trgZ=tmpX*this->ImTarget[0].Image2World[2][0]+tmpY*this->ImTarget[0].Image2World[2][1]+tmpZ*this->ImTarget[0].Image2World[2][2]+this->ImTarget[0].Image2World[2][3];
    
    tmpX=x*this->ImTarget[0].Image2World[0][0]+y*this->ImTarget[0].Image2World[0][1]+z*this->ImTarget[0].Image2World[0][2]+this->ImTarget[0].Image2World[0][3];
    tmpY=x*this->ImTarget[0].Image2World[1][0]+y*this->ImTarget[0].Image2World[1][1]+z*this->ImTarget[0].Image2World[1][2]+this->ImTarget[0].Image2World[1][3];
    tmpZ=x*this->ImTarget[0].Image2World[2][0]+y*this->ImTarget[0].Image2World[2][1]+z*this->ImTarget[0].Image2World[2][2]+this->ImTarget[0].Image2World[2][3];

    srcX=tmpX*this->World_Target2Template[0][0]+tmpY*this->World_Target2Template[0][1]+tmpZ*this->World_Target2Template[0][2]+this->World_Target2Template[0][3];
    srcY=tmpX*this->World_Target2Template[1][0]+tmpY*this->World_Target2Template[1][1]+tmpZ*this->World_Target2Template[1][2]+this->World_Target2Template[1][3];
    srcZ=tmpX*this->World_Target2Template[2][0]+tmpY*this->World_Target2Template[2][1]+tmpZ*this->World_Target2Template[2][2]+this->World_Target2Template[2][3];
    
    LocDisplacementField->P(trgX-srcX,0,x,y,z);
    LocDisplacementField->P(trgY-srcY,1,x,y,z);
    LocDisplacementField->P(trgZ-srcZ,2,x,y,z);
  }
  
  //save the deformation field
  LocDisplacementField->Write(FileNameX,FileNameY,FileNameZ,TargetFiles[0]);
  cout << "Remark: The " << Prefix << "_DispField_Src2Trg* files are saved in the target image domain."  << endl;
}





///save the Velocity field
void LargeDefDemons::SaveVelocityField(VectorField * LocVelocityField,char Prefix[256]){
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
  LocVelocityField->Write(FileNameX,FileNameY,FileNameZ,TargetFiles[0]);
}

///save the deformations in time subdivisions (not the convergence)
void LargeDefDemons::SaveDeformations(char Prefix[256]){
  char FileName[256];
  char Deformation[256];

  
  //deformed image
  strcpy(Deformation,"_FinalDefSrc.nii");
  strcpy(FileName,Prefix);
  strcat(FileName,Deformation);
    
  if (IniDispFieldDefined==1){
      ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
  }
  else{
      ProjectImageUsingAffineTransfoAndSteadyVeloField(Target2TemplateCoord,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
  }
  
  
  
  
  this->DeformedTemplate.Write(FileName,this->TargetFiles[0]);
}



///compute the update vector field
void LargeDefDemons::ComputeUpdateFieldSSD(){
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
  
  
  //project the source image in DeformedTemplate 
  //compute the temporary deformed template
  
  if (IniDispFieldDefined==1){
      if (IndicatorSaveMemory==1) this->GradE.SlashFieldSize();
      ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
      if (IndicatorSaveMemory==1) this->GradE.CreateVoidField(this->NX,this->NY,this->NZ);
  }
  else{
      if (IndicatorSaveMemory==1) this->GradE.SlashFieldSize();
      ProjectImageUsingAffineTransfoAndSteadyVeloField(Target2TemplateCoord,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
      if (IndicatorSaveMemory==1) this->GradE.CreateVoidField(this->NX,this->NY,this->NZ);
  }
  
  //compute the temporary gradient of the deformed template
  Cpt_Grad_ScalarField(&this->DeformedTemplate,&this->GradE);
  
  float Av1st_term;
  float Av2nd_term;
  int NbPts;
  
  Av1st_term=0;
  Av2nd_term=0;
  NbPts=0;
  
  //multiply the temporary gradient by the difference between the target and deformed images AND normalize the temporary gradient
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    locCoef=-(this->ImTarget[0].G(x,y,z)-this->DeformedTemplate.G(x,y,z));
    locCoef2=this->GradE.G(0,x,y,z)*this->GradE.G(0,x,y,z)+this->GradE.G(1,x,y,z)*this->GradE.G(1,x,y,z)+this->GradE.G(2,x,y,z)*this->GradE.G(2,x,y,z);
    
    
    if ((locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)))>epsilon){
      locCoef = locCoef/(locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)));
      Av1st_term+=locCoef2;
      Av2nd_term+=locCoef*locCoef;
      NbPts++;
    }
    else 
      locCoef=0;
    
    //take the margin into account  (3D IMAGE)
    if (this->NZ>1) if ((z<1+this->Margin)||(z>=this->NZ-1-this->Margin)||(y<1+this->Margin)||(y>=this->NY-1-this->Margin)||(x<1+this->Margin)||(x>=this->NX-1-this->Margin)) locCoef=0;

    //take the margin into account  (2D IMAGE)
    if (this->NZ==1) if ((y<1+this->Margin)||(y>=this->NY-1-this->Margin)||(x<1+this->Margin)||(x>=this->NX-1-this->Margin)) locCoef=0;

    //compute the update field
    this->GradE.P(this->GradE.G(0,x,y,z)*locCoef,0,x,y,z);
    this->GradE.P(this->GradE.G(1,x,y,z)*locCoef,1,x,y,z);
    this->GradE.P(this->GradE.G(2,x,y,z)*locCoef,2,x,y,z);
  }
  
  //mask the image if required
  if (this->MaskDefined==1)
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) 
      if (fabs(ProjMask.G(x,y,z)-MaskID)>0.01){
        this->GradE.P(0,0,x,y);
        this->GradE.P(0,1,x,y);
        this->GradE.P(0,2,x,y);
    }
}

///compute the update vector field
void LargeDefDemons::ComputeUpdateFieldSSD_multiChannel(){
  float locCoef;
  float locCoef2;
  int x,y,z;
  float epsilon;
  float Av1st_term;
  float Av2nd_term;
  int NbPts;
  int LocChannel;
  
  //1) init variables and gradient
  epsilon=0.0001;
  
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->GradE.P(0,0,x,y,z);
    this->GradE.P(0,1,x,y,z);
    this->GradE.P(0,2,x,y,z);
  }

  //2) loop on the channels
  for (LocChannel=0;LocChannel<this->NbChannels;LocChannel++){
    
    //2.1) initialize the temporary velocity field containing the gradients related to the current channel
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      this->TempVF.P(0,0,x,y,z);
      this->TempVF.P(0,1,x,y,z);
      this->TempVF.P(0,2,x,y,z);
    }

    //2.2) project the source image (or current source channel) in DeformedTemplate 
    if (IniDispFieldDefined==1){
        ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->ImTemplate[LocChannel],&this->DeformedTemplate);
    }
    else{
        ProjectImageUsingAffineTransfoAndSteadyVeloField(Target2TemplateCoord,&this->VelocityField,&this->ImTemplate[LocChannel],&this->DeformedTemplate);
    }    
    //2.3) compute the temporary gradient of the deformed template (... the case of 3D and 2D images is distinguished because of the boundary conditions)
    
    //2.3.1) compute the gradients of the deformed template
    Cpt_Grad_ScalarField(&this->DeformedTemplate,&this->TempVF);
    
    Av1st_term=0;
    Av2nd_term=0;
    NbPts=0;

    //2.3.2) Multiply the temporary gradient by the difference between the target and deformed images AND normalize the temporary gradient
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      locCoef=-(this->ImTarget[LocChannel].G(x,y,z)-this->DeformedTemplate.G(x,y,z));
      locCoef2=this->TempVF.G(0,x,y,z)*this->TempVF.G(0,x,y,z)+this->TempVF.G(1,x,y,z)*this->TempVF.G(1,x,y,z)+this->TempVF.G(2,x,y,z)*this->TempVF.G(2,x,y,z);
      
      //compute the coef
      if ((locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)))>epsilon){
        locCoef = locCoef/(locCoef2+(locCoef*locCoef/(this->lambdaX*this->lambdaX)));
        Av1st_term+=locCoef2;
        Av2nd_term+=locCoef*locCoef;
        NbPts++;
      }
      else 
        locCoef=0;
      
      //take the margin into account  (3D IMAGE)
      if (this->NZ>1) if ((z<1+this->Margin)||(z>=this->NZ-1-this->Margin)||(y<1+this->Margin)||(y>=this->NY-1-this->Margin)||(x<1+this->Margin)||(x>=this->NX-1-this->Margin)) locCoef=0;
      
      //take the margin into account  (2D IMAGE)
      if (this->NZ==1) if ((y<1+this->Margin)||(y>=this->NY-1-this->Margin)||(x<1+this->Margin)||(x>=this->NX-1-this->Margin)) locCoef=0;
      
      //compute the update
      this->TempVF.P(this->TempVF.G(0,x,y,z)*locCoef,0,x,y,z);
      this->TempVF.P(this->TempVF.G(1,x,y,z)*locCoef,1,x,y,z);
      this->TempVF.P(this->TempVF.G(2,x,y,z)*locCoef,2,x,y,z);
    }
    
    //2.4) update the gradients with the contribution of the current channel
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
      this->GradE.Add(this->weightsChannels[LocChannel]*this->TempVF.G(0,x,y,z),0,x,y,z);
      this->GradE.Add(this->weightsChannels[LocChannel]*this->TempVF.G(1,x,y,z),1,x,y,z);
      this->GradE.Add(this->weightsChannels[LocChannel]*this->TempVF.G(2,x,y,z),2,x,y,z);
    }
  }

  //3) compute the final gradient
  for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++){
    this->GradE.P(this->GradE.G(0,x,y,z)/this->NbChannels,0,x,y,z);
    this->GradE.P(this->GradE.G(1,x,y,z)/this->NbChannels,1,x,y,z);
    this->GradE.P(this->GradE.G(2,x,y,z)/this->NbChannels,2,x,y,z);  
  }
  

  //mask the image if required
  if (this->MaskDefined==1)
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) 
      if (fabs(ProjMask.G(x,y,z)-MaskID)>0.01){
        this->GradE.P(0,0,x,y);
        this->GradE.P(0,1,x,y);
        this->GradE.P(0,2,x,y);
      }
  
}



///compute the update vector field
void LargeDefDemons::ComputeUpdateFieldMI(){
  float evaMI;
  int x,y,z;
  
  //compute the current deformed source image
  if (IniDispFieldDefined==1){
      ProjectImageUsingDispFieldAndSteadyVeloField(&this->IniDispField,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
  }
  else{
      ProjectImageUsingAffineTransfoAndSteadyVeloField(Target2TemplateCoord,&this->VelocityField,&this->ImTemplate[0],&this->DeformedTemplate);
  }  
  
  this->NorMutInfMan.IndicateSrcHasChanged();
  
  //compute the mutual information and update the histograms
  this->NorMutInfMan.IndicateSrcHasChanged();
  evaMI=this->NorMutInfMan.EvaluateMI(); //!!! TO CHANGE !!! -> adding a mask
  cout << "MI: " << evaMI << endl;
  
  
  //evaluate the gradient of mutual information
  this->NorMutInfMan.EvaluateGradMI(&this->GradE); //!!! TO CHANGE !!! -> adding a mask
  
  
  //mask the image if required
  if (this->MaskDefined==1)
    for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) 
      if (fabs(ProjMask.G(x,y,z)-MaskID)>0.01){
        this->GradE.P(0,0,x,y);
        this->GradE.P(0,1,x,y);
        this->GradE.P(0,2,x,y);
      }
  
}



///control the maximum update. 
///* if the max update is larger than the initial one -> normalisation of all the update field to this->MaxUpdateAllowed
///* Othewise -> normalisation of all the update field to this->MaxUpdateAllowed*(max/InitMaxUpdate)
void LargeDefDemons::ControlMaxUpdate(int IterationNb){
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
float LargeDefDemons::EstimRefUpdateScale(float sigmaloc,LightFFTconvolver3D * ConvolverLoc){
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
void LargeDefDemons::Run_Default(void){
  int IterationNb;
  int x,y,z;
  
  //1) INITIALISATION
  
  //1.1) Pre-treatment of the inuput images (grey level alignment + margins)
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
      if (this->NbChannels==1) this->ComputeUpdateFieldSSD();
      else  this->ComputeUpdateFieldSSD_multiChannel();
    }
    
    //2.2 smooth the gradient... 
    if (sigma1>0.005){
      if (this->MaskDefined!=2){ //standard smoothing
        this->FFTconvolver_fluid.Convolution(&this->GradE);
      }
      else{//region-wise smoothing with sliding motion
        RemoveNormalContributions(&this->GradE,&this->NearestBoundary,&this->TempSF,this->BoundaMargin,this->x_mm,this->y_mm,this->z_mm);
        this->FFTconvolver_fluid.Convolution_Mask_Mirror(&this->GradE,&this->ProjMask);   //remark: ProjMask is not exactly the projected mask but equals to 1 in the smoothed ROI and equals to 0 otherwise
	for (z = 0; z < this->NZ; z++)  for (y = 0; y < this->NY; y++) for (x = 0; x < this->NX; x++) if (fabs(this->ProjMask.G(x,y,z)-1)>0.01){this->GradE.P(0,0,x,y,z); this->GradE.P(0,1,x,y,z); this->GradE.P(0,2,x,y,z);}
        RemoveNormalContributions(&this->GradE,&this->NearestBoundary,&this->TempSF,this->BoundaMargin,this->x_mm,this->y_mm,this->z_mm);
      }
    }
    
    //2.3 control the maximum amplitude of update
    this->ControlMaxUpdate(IterationNb);
    
    //2.4 update the transformation
      //ComposeTwoLogFieldsUsingBCH(&this->VelocityField,&this->GradE);
      ComposeTwoLogFieldsUsingSum(&this->VelocityField,&this->GradE);


    
    //2.5 smooth the deformation
    int x,y,z;
    if (sigmaDiff>0.005){
      if (this->MaskDefined!=2){ //standard smoothing
        this->FFTconvolver_Diff.Convolution(&this->VelocityField);
      }
      else{//region-wise smoothing with sliding motion
        RemoveNormalContributions(&this->VelocityField,&this->NearestBoundary,&this->TempSF,this->BoundaMargin,this->x_mm,this->y_mm,this->z_mm);
        this->FFTconvolver_Diff.Convolution_Mask_Mirror(&this->VelocityField,&this->ProjMask); //remark: ProjMask is not exactly the projected mask but equals to 1 in the smoothed ROI and equals to 0 otherwise
        RemoveNormalContributions(&this->VelocityField,&this->NearestBoundary,&this->TempSF,this->BoundaMargin,this->x_mm,this->y_mm,this->z_mm);
      }
    }
    
    
  }
  
  if (this->MaskDefined==2)
    RemoveNormalContributions(&this->VelocityField,&this->NearestBoundary,&this->TempSF,this->BoundaMargin,this->x_mm,this->y_mm,this->z_mm,1);
  
  //3) SAVE THE RESULTS
  this->SaveResultGradientDescent();
}



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                        MAIN RUN FUNCTION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///run function
void LargeDefDemons::Run(void)
{
  this->Run_Default();
}

