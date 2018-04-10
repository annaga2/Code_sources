/*=========================================================================
 
 Author: Laurent Risser
 
 Disclaimer: This software has been developed for research purposes only, and hence should 
 not be used as a diagnostic tool. In no event shall the authors or distributors
 be liable to any direct, indirect, special, incidental, or consequential 
 damages arising of the use of this software, its documentation, or any 
 derivatives thereof, even if the authors have been advised of the possibility 
 of such damage. 
 
 =========================================================================*/


#include <SciCalcPack.h>





/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 THE TOOLS
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ image transport +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void transportImage(char * ImageInit,char * ImageFinal,int TimeInit,int TimeFinal,char * VFX,char * VFY,char * VFZ){
  ScalarField ImInit;
  ScalarField ImFinal;
  VectorField VelocityField;
  VectorField Map;
  VectorField IdMap;
  int x,y,z;
  
  //init
  ImInit.Read(ImageInit);
  ImFinal.Read(ImageInit);  //to allocate it
  VelocityField.Read(VFX,VFY,VFZ);
  IdMap.CreateVoidField(VelocityField.NX,VelocityField.NY,VelocityField.NZ);
  
  for (z = 0; z < IdMap.NZ; z++)  for (y = 0; y < IdMap.NY; y++) for (x = 0; x < IdMap.NX; x++){
    IdMap.P(static_cast<float>(x),0,x,y,z);
    IdMap.P(static_cast<float>(y),1,x,y,z);
    IdMap.P(static_cast<float>(z),2,x,y,z);
  }
  
  
  //compute the mapping
  CptMappingFromVeloField(TimeInit,&IdMap,&VelocityField,&Map);
  
  //project the image
  Project3Dimage(&ImInit,&Map,&ImFinal,TimeFinal);
  
  //save the projected image
  ImFinal.Write(ImageFinal,ImageInit);
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ transport several images  +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//largedeformationsTools -MultiTransport [VFX][VFY][VFZ] [Src TimeSub] [Nb Ima] [Src Ima 1][Trg TimeSub 1][Trg Ima 1] ...
//-> Transport [Nb Ima] images [Src Ima i] from the time subdivision [Src TimeSub] to the time subdivision [Trg TimeSub i]
//through the velocity field [VFX][VFY][VFZ]. The transported images are saved in [Trg Ima i].
void MultiTransport(int argc, char **argv){
  char ImageInit[256];
  char ImageFinal[256];
  ScalarField ImInit;
  ScalarField ImFinal;
  VectorField VelocityField;
  VectorField Map;
  VectorField IdMap;
  int x,y,z;
  char VFX[256];
  char VFY[256];
  char VFZ[256];
  int NbImages;
  int TimeFinal;
  int TimeInit;
  int i;
  
  //read global parameters
  argc--; argv++;
  strcpy(VFX,argv[1]); //velocity field X
  argc--; argv++;
  strcpy(VFY,argv[1]); //velocity field Y
  argc--; argv++;
  strcpy(VFZ,argv[1]); //velocity field Z
  argc--; argv++;
  TimeInit = atoi(argv[1]);
  argc--;  argv++;
  NbImages = atoi(argv[1]);
  argc--;  argv++;
  
  //compute the mapping
  VelocityField.Read(VFX,VFY,VFZ);
  IdMap.CreateVoidField(VelocityField.NX,VelocityField.NY,VelocityField.NZ);
  
  for (z = 0; z < IdMap.NZ; z++)  for (y = 0; y < IdMap.NY; y++) for (x = 0; x < IdMap.NX; x++){
    IdMap.P(static_cast<float>(x),0,x,y,z);
    IdMap.P(static_cast<float>(y),1,x,y,z);
    IdMap.P(static_cast<float>(z),2,x,y,z);
  }
  
  CptMappingFromVeloField(TimeInit,&IdMap,&VelocityField,&Map);
  
  //transport all the images
  for (i=0;i<NbImages;i++){
    //read the parameters
    strcpy(ImageInit,argv[1]);
    argc--; argv++;
    TimeFinal = atoi(argv[1]);
    argc--;  argv++;
    strcpy(ImageFinal,argv[1]);
    argc--; argv++;
    
    //read and allocate the images
    ImInit.Read(ImageInit);
    ImFinal.Read(ImageInit);
    
    //project the image
    Project3Dimage(&ImInit,&Map,&ImFinal,TimeFinal);
    
    //save the projected image
    ImFinal.Write(ImageFinal,ImageInit);
  }
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ weighted sum +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int weightedSum(int argc, char **argv)
{
  ScalarField in;
  ScalarField out;
  int NbImages;
  double weight;
  char *input_name = NULL, *output_name = NULL;
  int x,y,z,t,i;
  
  
  argc--;
  argv++;
  NbImages = atoi(argv[1]);
  argc--;
  argv++;
  
  for (i=0;i<NbImages;i++){
    if (i==0){
      // first image ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);
      
      out.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      
      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.P(in.G(x,y,z,t)*weight,x,y,z,t);
    }
    else{
      // other images ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);

      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.Add(in.G(x,y,z,t)*weight,x,y,z,t);
    }
  }

  output_name = argv[1];
  argc--;
  argv++;
  out.Write(output_name,input_name);
  
  return 0;
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++++ weighted sum +++++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int weightedMean(int argc, char **argv)
{
  ScalarField in;
  ScalarField out;
  int NbImages;
  double weight;
  double SumWeight;
  char *input_name = NULL, *output_name = NULL;
  int x,y,z,t,i;
  
  
  argc--;
  argv++;
  NbImages = atoi(argv[1]);
  argc--;
  argv++;
  
  SumWeight=0.;
  
  for (i=0;i<NbImages;i++){
    if (i==0){
      // first image ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);
      
      out.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      
      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.P(in.G(x,y,z,t)*weight,x,y,z,t);
      
      SumWeight+=weight;
    }
    else{
      // other images ...
      //... read parameters
      weight = atof(argv[1]);
      argc--;
      argv++;
      input_name  = argv[1];
      argc--;
      argv++;
      //... do the treatment
      in.Read(input_name);
      
      for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
        out.Add(in.G(x,y,z,t)*weight,x,y,z,t);
      
      SumWeight+=weight;
    }
  }
  
  //compute the mean
  for (t = 0; t < out.NT; t++) for (z = 0; z < out.NZ; z++)  for (y = 0; y < out.NY; y++) for (x = 0; x < out.NX; x++)
    out.P(out.G(x,y,z,t)/SumWeight,x,y,z,t);
  
  
  output_name = argv[1];
  argc--;
  argv++;
  out.Write(output_name,input_name);
  
  return 0;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++ weighted standard deviation ++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

int weightedStd(int argc, char **argv)
{
  ScalarField in;
  ScalarField std;
  ScalarField mean;
  int NbImages;
  double weight;
  double SumWeight;
  char *input_name = NULL, *output_name = NULL;
  int x,y,z,t,i;
  
  argc--;
  argv++;
  NbImages = atoi(argv[1]);
  argc--;
  argv++;
  
  
  //1) compute the mean
  SumWeight=0.;
  
  for (i=0;i<NbImages;i++){
    //... read parameters
    weight = atof(argv[2*i+1]);
    input_name  = argv[2*i+2];
    
    //... do the treatment
    in.Read(input_name);
    
    if (i==0){
      mean.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      for (t = 0; t < mean.NT; t++) for (z = 0; z < mean.NZ; z++)  for (y = 0; y < mean.NY; y++) for (x = 0; x < mean.NX; x++)
        mean.P(in.G(x,y,z,t)*weight,x,y,z,t);
    }
    else{
      for (t = 0; t < mean.NT; t++) for (z = 0; z < mean.NZ; z++)  for (y = 0; y < mean.NY; y++) for (x = 0; x < mean.NX; x++)
        mean.Add(in.G(x,y,z,t)*weight,x,y,z,t);
    }
    
    SumWeight+=weight;
  }
  
  //divide by the sum of the weights to obtain a mean
  for (t = 0; t < mean.NT; t++) for (z = 0; z < mean.NZ; z++)  for (y = 0; y < mean.NY; y++) for (x = 0; x < mean.NX; x++)
    mean.P(mean.G(x,y,z,t)/SumWeight,x,y,z,t);
  
  
  //2) compute the standard deviation
  
  SumWeight=0.;
  
  for (i=0;i<NbImages;i++){
    //... read parameters
    weight = atof(argv[2*i+1]);
    input_name  = argv[2*i+2];
    
    //... do the treatment
    in.Read(input_name);
    
    if (i==0){
      std.CreateVoidField(in.NX,in.NY,in.NZ,in.NT);
      for (t = 0; t < std.NT; t++) for (z = 0; z < std.NZ; z++)  for (y = 0; y < std.NY; y++) for (x = 0; x < std.NX; x++)
        std.P((in.G(x,y,z,t)-mean.G(x,y,z,t))*(in.G(x,y,z,t)-mean.G(x,y,z,t))*weight,x,y,z,t);
    }
    else{
      for (t = 0; t < std.NT; t++) for (z = 0; z < std.NZ; z++)  for (y = 0; y < std.NY; y++) for (x = 0; x < std.NX; x++)
        std.Add((in.G(x,y,z,t)-mean.G(x,y,z,t))*(in.G(x,y,z,t)-mean.G(x,y,z,t))*weight,x,y,z,t);
    }
    
    SumWeight+=weight;
  }
  
  //divide by the sum of the weights and compute the square root to obtain a std
  for (t = 0; t < std.NT; t++) for (z = 0; z < std.NZ; z++)  for (y = 0; y < std.NY; y++) for (x = 0; x < std.NX; x++)
    std.P(sqrt(std.G(x,y,z,t)/SumWeight),x,y,z,t);
  
  
  output_name = argv[2*NbImages+1];
  argc--;
  argv++;
  std.Write(output_name,input_name);
  
  return 0;
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++  basic information about an image +++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Info(int argc, char **argv)
{
  ScalarField image1;
  int x,y,z,t;
  int i,j;
  char File1[256];
  float minGL,maxGL;
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  
  
  //intialisation
  image1.Read(File1);
  
  //image size
  cout << "Image size: " << image1.NX << " " << image1.NY << " " << image1.NZ << " " << image1.NT << endl;
  
  //image to world matrix
  cout << endl;
  cout << "Image to world matrix:" << endl;
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      cout << image1.Image2World[i][j] << " ";
    }
    cout << endl;
  }
  
  //image to world matrix
  cout << endl;
  cout << "World to image matrix:" << endl;
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      cout << image1.World2Image[i][j] << " ";
    }
    cout << endl;
  }
  
  //min and max grey levels (for relatively small images only)
  if (image1.NX*image1.NY*image1.NZ*image1.NT<200*200*200){
    minGL=image1.G(0,0,0);
    maxGL=image1.G(0,0,0);
    for (t = 0; t < image1.NT; t++) for (z = 0; z < image1.NZ; z++)  for (y = 0; y < image1.NY; y++) for (x = 0; x < image1.NX; x++){
      if (image1.G(x,y,z,t)>maxGL) maxGL=image1.G(x,y,z,t);
      if (image1.G(x,y,z,t)<minGL) minGL=image1.G(x,y,z,t);
    }
    
    cout << endl;
    cout << "Min grey level: " <<  minGL << endl;
    cout << "Max grey level: " <<  maxGL << endl;
  }
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++ sum of squared differences between two images +++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double MeasureSSD(int argc, char **argv)
{
  ScalarField image1;
  ScalarField image2;
  double SSD;
  int x,y,z,t;
  int nbPts;
  char File1[256];
  char File2[256];
  int  margin,testMargin;
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  argc--; argv++;
  
  margin=0;
  if (argc>1){
    margin = atoi(argv[1]);
    argc--; argv++;
  }
  
  
  //intialisation
  image1.Read(File1);
  image2.Read(File2);
  
  cout << File1 << " / " << File2 << ":" << endl;
  
  
  SSD=0.;
  nbPts=0;
  
  //computation
  for (t = 0; t < image1.NT; t++) for (z = 0; z < image1.NZ; z++)  for (y = 0; y < image1.NY; y++) for (x = 0; x < image1.NX; x++){
    testMargin=1;
    
    if ((y<margin)||(y>image1.NY-1-margin)||(x<margin)||(x>image1.NX-1-margin)) testMargin=0;
    if (image1.NZ>1) if ((z<margin)||(z>image1.NZ-1-margin))  testMargin=0;
    
    if (testMargin==1){
      SSD+=(image1.G(x,y,z,t)-image2.G(x,y,z,t))*(image1.G(x,y,z,t)-image2.G(x,y,z,t));
      nbPts++;
    }
  }
  
  SSD=sqrt(SSD);
  SSD/=(double)nbPts;
  
  cout << "SSD=" << SSD << "\n";
  
  
  
  //BEGIN BONUS: GIVE ALSO THE MUTUAL INFORMATION
  MImanager MiMa;
  float LocMI;
  
  if (margin>0) cout << "MI computed in the whole images. The margin is only taken into account for the SSD." << endl;
  
  MiMa.Initiate(&image1,&image2,40,40,3);
  LocMI=MiMa.EvaluateMI();
  cout << "MI=" << LocMI << endl;
  //END BONUS: GIVE ALSO THE MUTUAL INFORMATION
  
  
  
  
  return SSD;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++        Dice metric between two images         +++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double MeasureOverlap(int argc, char **argv)
{
  ScalarField image1;
  ScalarField image2;
  int x,y,z,t;
  char File1[256];
  char File2[256];
  int  margin,testMargin;
  float thresh;
  int NbAndPoints;
  int NbOrPoints;
  double DiceM;
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  argc--; argv++;
  thresh = atof(argv[1]);
  argc--; argv++;
  
  margin=0;
  if (argc>1){
    margin = atoi(argv[1]);
    argc--; argv++;
  }
  
  
  //intialisation
  image1.Read(File1);
  image2.Read(File2);
  
  cout << File1 << " / " << File2 << ":" << endl;
  
  
  NbAndPoints=0.;
  NbOrPoints=0;
  
  //computation
  for (t = 0; t < image1.NT; t++) for (z = 0; z < image1.NZ; z++)  for (y = 0; y < image1.NY; y++) for (x = 0; x < image1.NX; x++){
    testMargin=1;
    
    if ((y<margin)||(y>image1.NY-1-margin)||(x<margin)||(x>image1.NX-1-margin)) testMargin=0;
    if (image1.NZ>1) if ((z<margin)||(z>image1.NZ-1-margin))  testMargin=0;
    
    if (testMargin==1){
      if ((image1.G(x,y,z,t)>thresh)&&(image2.G(x,y,z,t)>thresh)) NbAndPoints++;
      if ((image1.G(x,y,z,t)>thresh)||(image2.G(x,y,z,t)>thresh)) NbOrPoints++;
    }
  }
  
  DiceM=static_cast<double>(NbAndPoints)/static_cast<double>(NbOrPoints);
  
  
  cout << NbAndPoints << "  " << NbOrPoints << "\n";
  
  cout << "Overlap=" << DiceM << "\n";
  
  return DiceM;
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++        Dice metric between two images         +++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MeasureOverlap2(int argc, char **argv)
{
  ScalarField image1;
  ScalarField image2;
  int x,y,z;
  int i;
  char File1[256];
  char File2[256];
  int LocID,LocID_S,LocID_T,MaxIdRegion;
  int * Card_Id_Intersection;
  int * Card_Id_Source;
  int * Card_Id_Target;
  double TargetOverlap,MeanOverlap,tmp1,tmp2;
  char OutputFileValueOverlap[256];
  int MakeOutputfile;
  FILE *outFile;

  //1) read parameters / init
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  argc--; argv++;
  
  //read the name of the optional output file
  MakeOutputfile=0;
  if (argc>1) {
    MakeOutputfile=1;
    strcpy(OutputFileValueOverlap,argv[1]);
    argc--; argv++;
  }
  
  image1.Read(File1);  // Source
  image2.Read(File2);  // Target
  
  //2) Evaluate Regions...
  //... identify the max ID
  MaxIdRegion=0;
  for (z = 0; z < image1.NZ; z++)  for (y = 0; y < image1.NY; y++) for (x = 0; x < image1.NX; x++){
    if (image1.G(x,y,z)>0.5) if ( fabs(image1.G(x,y,z)-fabs(image1.G(x,y,z)))<0.01){ //strictly positive integer value
      LocID=static_cast<int>(image1.G(x,y,z)+0.1);
      if (MaxIdRegion<LocID) MaxIdRegion=LocID;
    }
    if (image2.G(x,y,z)>0.5) if ( fabs(image2.G(x,y,z)-fabs(image2.G(x,y,z)))<0.01){ //strictly positive integer value
      LocID=static_cast<int>(image2.G(x,y,z)+0.1);
      if (MaxIdRegion<LocID) MaxIdRegion=LocID;
    }
  }
  
  if (MaxIdRegion>10000){
    cout << "The maximum region identifier is too high. Please resample region identifiers" << endl;
    return;
  }
  
  //... compute the number of points in (unions, source, target) for each region
  Card_Id_Intersection= new int[MaxIdRegion];
  Card_Id_Source= new int[MaxIdRegion];
  Card_Id_Target= new int[MaxIdRegion];
  
  for (i=0;i<MaxIdRegion;i++) Card_Id_Intersection[i]=0;
  for (i=0;i<MaxIdRegion;i++) Card_Id_Source[i]=0;
  for (i=0;i<MaxIdRegion;i++) Card_Id_Target[i]=0;
  
  for (z = 0; z < image1.NZ; z++)  for (y = 0; y < image1.NY; y++) for (x = 0; x < image1.NX; x++){
    LocID_S=-1;
    LocID_T=-1;
    if (image1.G(x,y,z)>0.5) if ( fabs(image1.G(x,y,z)-fabs(image1.G(x,y,z)))<0.01) LocID_S=static_cast<int>(image1.G(x,y,z)+0.1);
    if (image2.G(x,y,z)>0.5) if ( fabs(image2.G(x,y,z)-fabs(image2.G(x,y,z)))<0.01) LocID_T=static_cast<int>(image2.G(x,y,z)+0.1);
    
    if (LocID_S>0) Card_Id_Source[LocID_S]++;
    if (LocID_T>0) Card_Id_Target[LocID_T]++;
    if ((LocID_S>0)&&(LocID_S==LocID_T)) Card_Id_Intersection[LocID_S]++;
  }
  
  //3) compute and return the target overlap and the mean overlap (dice coef.) ...
  cout << "Source: " << File1 << " /  Target: " << File2 << " :" << endl;
    
  //... target overlap
  tmp1=0; tmp2=0;
  for (i=0;i<MaxIdRegion;i++) tmp1+=static_cast<double>(Card_Id_Intersection[i]);
  for (i=0;i<MaxIdRegion;i++) tmp2+=static_cast<double>(Card_Id_Target[i]);
  
  TargetOverlap=tmp1/tmp2;
  cout << "Target overlap=" << TargetOverlap << "\n";
  
  //... mean overlap
  tmp1=0; tmp2=0;
  for (i=0;i<MaxIdRegion;i++) tmp1+=static_cast<double>(Card_Id_Intersection[i]);
  for (i=0;i<MaxIdRegion;i++) tmp2+=static_cast<double>(Card_Id_Target[i]+Card_Id_Source[i]);
  
  MeanOverlap=2*tmp1/tmp2;
  cout << "Mean overlap (Dice coef)=" << MeanOverlap << "\n";
  
  
  if (MakeOutputfile==1){
    outFile=fopen(OutputFileValueOverlap,"w");
    fprintf(outFile,"%f\n",MeanOverlap);
    //fclose(outFile);
  }
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++ surfacic average and std values ++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double Surf_Av_Std_values(int argc, char **argv)
{
  ScalarField RefIma;
  ScalarField SalarsIma;
  double Thresh;
  int x,y,z,t;
  int nbPts;
  char File1[256];
  char File2[256];
  double av,std;
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  Thresh = atof(argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  
  //intialisation
  RefIma.Read(File1);
  SalarsIma.Read(File2);
  
  av=0.;
  nbPts=0;
  
  //computation
  for (t = 0; t < RefIma.NT; t++) for (z = 0; z < RefIma.NZ; z++)  for (y = 0; y < RefIma.NY; y++) for (x = 0; x < RefIma.NX; x++)
    if (RefIma.G(x,y,z,t)>Thresh)
      if ((RefIma.G(x+1,y,z,t)<Thresh)||(RefIma.G(x-1,y,z,t)<Thresh)||(RefIma.G(x,y+1,z,t)<Thresh)||(RefIma.G(x,y-1,z,t)<Thresh)||(RefIma.G(x,y,z+1,t)<Thresh)||(RefIma.G(x,y,z-1,t)<Thresh)){
        av+=SalarsIma.G(x,y,z,t);
        nbPts++;
      }
  
  av/=(double)nbPts;
  
  std=0;
  
  for (t = 0; t < RefIma.NT; t++) for (z = 0; z < RefIma.NZ; z++)  for (y = 0; y < RefIma.NY; y++) for (x = 0; x < RefIma.NX; x++)
    if (RefIma.G(x,y,z,t)>Thresh)
      if ((RefIma.G(x+1,y,z,t)<Thresh)||(RefIma.G(x-1,y,z,t)<Thresh)||(RefIma.G(x,y+1,z,t)<Thresh)||(RefIma.G(x,y-1,z,t)<Thresh)||(RefIma.G(x,y,z+1,t)<Thresh)||(RefIma.G(x,y,z-1,t)<Thresh)){
        std+=(SalarsIma.G(x,y,z,t)-av)*(SalarsIma.G(x,y,z,t)-av);
      }
  
  std/=(double)nbPts;
  std=sqrt(std);
  
  cout << "av=" << av <<  "   std=" << std << "\n";
  return av;
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++ quantify a deformation in a ROI ++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double DispFieldStudy(int argc, char **argv)
{
  ScalarField DefX;
  ScalarField DefY;
  ScalarField DefZ;
  ScalarField Mask;
  ScalarField DispField;
  ScalarField DivergenceField;
  ScalarField DetJField;
  int indicatorMask;
  float IDMask,epsilon;
  int x,y,z,t;
  int nbPts;
  char File1[256];
  char File2[256];
  char File3[256];
  char File4[256];
  float av,std,minDJac,maxDJac,NbNegJacobians,tmp;
  float d11,d12,d13,d21,d22,d23,d31,d32,d33,zttt;
  char Output[256];
  float x_mm,y_mm,z_mm;
  char OutputFileMaxDetJ[256];
  FILE *outFile;
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  argc--; argv++;
  strcpy(File3,argv[1]);
  
  //load the mask if defined
  indicatorMask=0;
  if (argc==4) {
    indicatorMask=1;
    argc--; argv++;
    strcpy(File4,argv[1]);
    argc--; argv++;
    IDMask = atof(argv[1]);
  }
  
  //argc--; argv++;
  //strcpy(File4,argv[1]);
  //argc--; argv++;
  //IDMask = atof(argv[1]);
  
  
  
  //intialisation
  DefX.Read(File1);
  DefY.Read(File2);
  DefZ.Read(File3);
  if (indicatorMask==1){
    Mask.Read(File4);
  }
  else {
    Mask.Read(File3);
    for (z = 0; z < Mask.NZ; z++)  for (y = 0; y < Mask.NY; y++) for (x = 0; x < Mask.NX; x++) Mask.P(1,x,y,z);
    IDMask=1;
  }
  
  
  DispField.Read(File1);
  for (z = 0; z < DefX.NZ; z++)  for (y = 0; y < DefX.NY; y++) for (x = 0; x < DefX.NX; x++) DispField.P(0,x,y,z);
  DivergenceField.Read(File1);
  for (z = 0; z < DefX.NZ; z++)  for (y = 0; y < DefX.NY; y++) for (x = 0; x < DefX.NX; x++) DivergenceField.P(0,x,y,z);
  DetJField.Read(File1);
  for (z = 0; z < DefX.NZ; z++)  for (y = 0; y < DefX.NY; y++) for (x = 0; x < DefX.NX; x++) DetJField.P(0,x,y,z);
  
  epsilon=0.001;
  
  //1) norm displacement
  std=0.;
  av=0.;
  nbPts=0;
  
  if (DefX.NZ!=1){ //3D image
    for (z = 1; z < DefX.NZ-1; z++)  for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      tmp=sqrt((DefX.G(x,y,z)*DefX.G(x,y,z))+(DefY.G(x,y,z)*DefY.G(x,y,z))+(DefZ.G(x,y,z)*DefZ.G(x,y,z)));
      
      DispField.P(tmp,x,y,z);
      
      av+=tmp;
      nbPts++;
    }
    av/=nbPts;
    
    for (z = 1; z < DefX.NZ-1; z++)  for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      tmp=sqrt((DefX.G(x,y,z)*DefX.G(x,y,z))+(DefY.G(x,y,z)*DefY.G(x,y,z))+(DefZ.G(x,y,z)*DefZ.G(x,y,z)));
      
      std+=(tmp-av)*(tmp-av);
    }
    std/=nbPts;
    std=sqrt(std);
  }
  else{ //2D image
    z=0;
    for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      tmp=sqrt((DefX.G(x,y,z)*DefX.G(x,y,z))+(DefY.G(x,y,z)*DefY.G(x,y,z)));
      
      DispField.P(tmp,x,y,z);
      
      av+=tmp;
      nbPts++;
    }
    av/=nbPts;
    
    for (z = 1; z < DefX.NZ-1; z++)  for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      tmp=sqrt((DefX.G(x,y,z)*DefX.G(x,y,z))+(DefY.G(x,y,z)*DefY.G(x,y,z)));
      
      std+=(tmp-av)*(tmp-av);
    }
    std/=nbPts;
    std=sqrt(std);
  }
  
  
  cout << "Norm of the displacement: av=" << av <<  " /  std=" << std << "  (Nb pts=" <<  nbPts <<  ")\n";
  
  //2) divergence
  x_mm=sqrt(DefX.Image2World[0][0]*DefX.Image2World[0][0]+DefX.Image2World[0][1]*DefX.Image2World[0][1]+DefX.Image2World[0][2]*DefX.Image2World[0][2]);
  y_mm=sqrt(DefX.Image2World[1][0]*DefX.Image2World[1][0]+DefX.Image2World[1][1]*DefX.Image2World[1][1]+DefX.Image2World[1][2]*DefX.Image2World[1][2]);
  z_mm=sqrt(DefX.Image2World[2][0]*DefX.Image2World[2][0]+DefX.Image2World[2][1]*DefX.Image2World[2][1]+DefX.Image2World[2][2]*DefX.Image2World[2][2]);
  
  if (DefX.NZ!=1){ //3D image
    for (z = 1; z < DefX.NZ-1; z++)  for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon)
      DivergenceField.P(((DefX.G(x+1,y,z)-DefX.G(x-1,y,z))/x_mm)+((DefY.G(x,y+1,z)-DefY.G(x,y-1,z))/y_mm)+((DefZ.G(x,y,z+1)-DefZ.G(x,y,z-1))/z_mm),x,y,z);
  }
  else{ //2D image
    z=0;
    for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon)
      DivergenceField.P(((DefX.G(x+1,y,z)-DefX.G(x-1,y,z))/x_mm)+((DefY.G(x,y+1,z)-DefY.G(x,y-1,z))/y_mm),x,y,z);
  }
  
  
  //3) determinent of the Jacobian
  std=0.;
  av=0.;
  NbNegJacobians=0;
  cout << "toto" << endl;
  //cout << DefX.Image2World[0][0] << " " << DefY.Image2World[1][1] << " " << DefZ.Image2World[2][2]<< endl;
  if (DefX.NZ!=1){ //3D image
    for (z = 1; z < DefX.NZ-1; z++)  for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      d11=1+((DefX.G(x+1,y,z)-DefX.G(x-1,y,z))/(2.*DefX.Image2World[0][0]));
      d12=(DefX.G(x,y+1,z)-DefX.G(x,y-1,z))/(2.*DefX.Image2World[0][0]);
      d13=(DefX.G(x,y,z+1)-DefX.G(x,y,z-1))/(2.*DefX.Image2World[0][0]);
      
      d21=(DefY.G(x+1,y,z)-DefY.G(x-1,y,z))/(2.*DefY.Image2World[1][1]);
      d22=1+ ((DefY.G(x,y+1,z)-DefY.G(x,y-1,z))/(2.*DefY.Image2World[1][1]));
      d23=(DefY.G(x,y,z+1)-DefY.G(x,y,z-1))/(2.*DefY.Image2World[1][1]);
      
      d31=(DefZ.G(x+1,y,z)-DefZ.G(x-1,y,z))/(2.*DefZ.Image2World[2][2]);
      d32=(DefZ.G(x,y+1,z)-DefZ.G(x,y-1,z))/(2.*DefZ.Image2World[2][2]);
      d33=1+((DefZ.G(x,y,z+1)-DefZ.G(x,y,z-1))/(2.*DefZ.Image2World[2][2]));
      
      tmp=d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13);
      
      DetJField.P(tmp,x,y,z);
      
      av+=tmp;
      
      if (tmp<-0.01) NbNegJacobians++;
      
      if ((x==1)&&(y==1)&&(z==1)){
        minDJac=tmp;
        maxDJac=tmp;
      }
      else{
        if (minDJac>tmp) minDJac=tmp;
        if (maxDJac<tmp) maxDJac=tmp;
      }
    }
  }
  else{ //2D image
    z=0;
    for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      d11=(DefX.G(x+1,y,z)-DefX.G(x-1,y,z)+2.*DefX.Image2World[0][0])/(2.*DefX.Image2World[0][0]);
      d12=(DefX.G(x,y+1,z)-DefX.G(x,y-1,z))/(2.*DefX.Image2World[0][0]);
      d13=1;
      d21=(DefY.G(x+1,y,z)-DefY.G(x-1,y,z))/(2.*DefY.Image2World[1][1]);
      d22=(DefY.G(x,y+1,z)-DefY.G(x,y-1,z)+2.*DefX.Image2World[1][1])/(2.*DefY.Image2World[1][1]);
      d23=1;
      d31=(DefZ.G(x+1,y,z)-DefZ.G(x-1,y,z))/(2.*DefZ.Image2World[2][2]);
      d32=(DefZ.G(x,y+1,z)-DefZ.G(x,y-1,z))/(2.*DefZ.Image2World[2][2]);
      d33=1;
      tmp=d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13);
      
      DetJField.P(tmp,x,y,z);
      
      av+=tmp;
      
      if (tmp<-0.01) NbNegJacobians++;
      
      if ((x==1)&&(y==1)&&(z==1)){
        minDJac=tmp;
        maxDJac=tmp;
      }
      else{
        if (minDJac>tmp) minDJac=tmp;
        if (maxDJac<tmp) maxDJac=tmp;
      }
    }
  }
  av/=nbPts;
  
  if (DefX.NZ!=1){ //3D image
    for (z = 1; z < DefX.NZ-1; z++)  for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      d11=(DefX.G(x+1,y,z)-DefX.G(x-1,y,z)+2.*DefX.Image2World[0][0])/(2.*DefX.Image2World[0][0]);
      d12=(DefX.G(x,y+1,z)-DefX.G(x,y-1,z))/(2.*DefX.Image2World[0][0]);
      d13=(DefX.G(x,y,z+1)-DefX.G(x,y,z-1))/(2.*DefX.Image2World[0][0]);
      d21=(DefY.G(x+1,y,z)-DefY.G(x-1,y,z))/(2.*DefY.Image2World[1][1]);
      d22=(DefY.G(x,y+1,z)-DefY.G(x,y-1,z)+2.*DefX.Image2World[1][1])/(2.*DefY.Image2World[1][1]);
      d23=(DefY.G(x,y,z+1)-DefY.G(x,y,z-1))/(2.*DefY.Image2World[1][1]);
      d31=(DefZ.G(x+1,y,z)-DefZ.G(x-1,y,z))/(2.*DefZ.Image2World[2][2]);
      d32=(DefZ.G(x,y+1,z)-DefZ.G(x,y-1,z))/(2.*DefZ.Image2World[2][2]);
      d33=(DefZ.G(x,y,z+1)-DefZ.G(x,y,z-1)+2.*DefX.Image2World[2][2])/(2.*DefZ.Image2World[2][2]);
      tmp=d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13);
      
      std+=(tmp-av)*(tmp-av);
    }
  }
  else{ //2D image
    z=0;
    for (y = 1; y < DefX.NY-1; y++) for (x = 1; x < DefX.NX-1; x++) if (fabs(Mask.G(x,y,z)-IDMask)<epsilon){
      d11=(DefX.G(x+1,y,z)-DefX.G(x-1,y,z)+2.*DefX.Image2World[0][0])/(2.*DefX.Image2World[0][0]);
      d12=(DefX.G(x,y+1,z)-DefX.G(x,y-1,z))/(2.*DefX.Image2World[0][0]);
      d13=1;
      d21=(DefY.G(x+1,y,z)-DefY.G(x-1,y,z))/(2.*DefY.Image2World[1][1]);
      d22=(DefY.G(x,y+1,z)-DefY.G(x,y-1,z)+2.*DefX.Image2World[1][1])/(2.*DefY.Image2World[1][1]);
      d23=1;
      d31=(DefZ.G(x+1,y,z)-DefZ.G(x-1,y,z))/(2.*DefZ.Image2World[2][2]);
      d32=(DefZ.G(x,y+1,z)-DefZ.G(x,y-1,z))/(2.*DefZ.Image2World[2][2]);
      d33=1;
      tmp=d11*(d22*d33-d32*d23)-d21*(d12*d33-d32*d13)+d31*(d12*d23-d22*d13);
      
      std+=(tmp-av)*(tmp-av);
    }
  }
  
  std/=nbPts;
  std=sqrt(std);
  
  
  NbNegJacobians=100*NbNegJacobians/nbPts;
  
  
  //show and save the results
  
  cout << "Determinent of the Jacobians: av=" << av <<  " /  std=" << std <<   " /  min=" << minDJac <<   " /  max=" << maxDJac <<    " /  % neg. DetJ=" << NbNegJacobians << "\n";
  
  
  strcpy(Output,"Out_DetJ.nii");
  DetJField.Write(Output,File1);
  
  strcpy(Output,"Out_DefAmplitude.nii");
  DispField.Write(Output,File1);
  
  strcpy(Output,"Out_Divergence.nii");
  DivergenceField.Write(Output,File1);
  
  strcpy(OutputFileMaxDetJ,"Out_MaxDetJ.txt");
  outFile=fopen(OutputFileMaxDetJ,"w");
  fprintf(outFile,"%f\n",maxDJac);
  
  return av;
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++ ROI average and std values ++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

double ROI_Av_Std_values(int argc, char **argv)
{
  ScalarField SalarsIma;
  int x,y,z,t;
  int nbPts;
  char File1[256];
  double av,std;
  int Rx1,Rx2,Ry1,Ry2,Rz1,Rz2;
  
  //read parameters
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  Rx1 = atoi(argv[1]);
  argc--; argv++;
  Rx2 = atoi(argv[1]);
  argc--; argv++;
  Ry1 = atoi(argv[1]);
  argc--; argv++;
  Ry2 = atoi(argv[1]);
  argc--; argv++;
  Rz1 = atoi(argv[1]);
  argc--; argv++;
  Rz2 = atoi(argv[1]);
  
  //intialisation
  SalarsIma.Read(File1);
  
  av=0.;
  nbPts=0;
  
  //treat the negative values
  if (Rx1<0) Rx1=SalarsIma.NX-Rx1;
  if (Rx2<0) Rx2=SalarsIma.NX-Rx2;
  if (Ry1<0) Ry1=SalarsIma.NY-Ry1;
  if (Ry2<0) Ry2=SalarsIma.NY-Ry2;
  if (Rz1<0) Rz1=SalarsIma.NZ-Rz1;
  if (Rz2<0) Rz2=SalarsIma.NZ-Rz2;
  
  
  
  //computation
  for (t = 0; t < SalarsIma.NT; t++) for (z = 0; z < SalarsIma.NZ; z++)  for (y = 0; y < SalarsIma.NY; y++) for (x = 0; x < SalarsIma.NX; x++)
    if ((x>=Rx1)||(x<Rx2)||(y>=Ry1)||(y<Ry2)||(z>=Rz1)||(z<Rz2)){
      av+=SalarsIma.G(x,y,z,t);
      nbPts++;
    }
  
  av/=(double)nbPts;
  
  std=0;
  
  for (t = 0; t < SalarsIma.NT; t++) for (z = 0; z < SalarsIma.NZ; z++)  for (y = 0; y < SalarsIma.NY; y++) for (x = 0; x < SalarsIma.NX; x++)
    if ((x>=Rx1)||(x<Rx2)||(y>=Ry1)||(y<Ry2)||(z>=Rz1)||(z<Rz2)){
      std+=(SalarsIma.G(x,y,z,t)-av)*(SalarsIma.G(x,y,z,t)-av);
    }
  
  std/=(double)nbPts;
  std=sqrt(std);
  
  cout << "av=" << av <<  "   std=" << std << "\n";
  return av;
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++ Anisotropic diffusion +++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void AnisoDiff(int argc, char **argv)
{
  ScalarField SField;
  float ax,ay,az,dTau;
  int interationNb;
  char Input[256];
  char Output[256];
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  ax = atof(argv[1]);
  argc--; argv++;
  ay = atof(argv[1]);
  argc--; argv++;
  az = atof(argv[1]);
  argc--; argv++;
  interationNb = atoi(argv[1]);
  argc--; argv++;
  dTau = atof(argv[1]);
  
  
  //read, diffuse and write the image
  SField.Read(Input);
  anisoDiff_3D(&SField,ax,ay,az,dTau,interationNb);
  SField.Write(Output,Input);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++ Isotropic diffusion ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Diffusion(int argc, char **argv)
{
  ScalarField SField;
  ScalarField Mask;
  float alpha,dTau;
  int iterationNb;
  char Input[256];
  char Output[256];
  char MaskFile[256];
  int indicatorMask;
  int IdMask;
  float resX,resY,resZ;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  iterationNb = atoi(argv[1]);
  argc--; argv++;
  dTau = atof(argv[1]);
  
  alpha = 1;
  
  //load the mask if defined
  indicatorMask=0;
  if (argc==4) {
    indicatorMask=1;
    argc--; argv++;
    strcpy(MaskFile,argv[1]);
    Mask.Read(MaskFile);
    argc--; argv++;
    IdMask = atoi(argv[1]);
  }
  
  
  //read the image
  SField.Read(Input);
  
  
  resX=sqrt(SField.Image2World[0][0]*SField.Image2World[0][0]+SField.Image2World[0][1]*SField.Image2World[0][1]+SField.Image2World[0][2]*SField.Image2World[0][2]);
  resY=sqrt(SField.Image2World[1][0]*SField.Image2World[1][0]+SField.Image2World[1][1]*SField.Image2World[1][1]+SField.Image2World[1][2]*SField.Image2World[1][2]);
  resZ=sqrt(SField.Image2World[2][0]*SField.Image2World[2][0]+SField.Image2World[2][1]*SField.Image2World[2][1]+SField.Image2World[2][2]*SField.Image2World[2][2]);
  
  cout << "Voxel size: " << resX << " " <<  resY << " " <<  resZ << endl;
  
  //diffuse the image
  if (indicatorMask==0)
    Diffusion_3D(&SField,alpha,dTau,iterationNb,resX,resY,resZ);
  else {
    Diffusion_3D(&SField,&Mask,IdMask,alpha,dTau,iterationNb,0,resX,resY,resZ);
  }
  
  //write the image
  SField.Write(Output,Input);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++ Isotropic smoothing ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SmoothImage(int argc, char **argv)
{
  ScalarField SField;
  LightFFTconvolver3D FFTconvolver;
  float sigma;
  char Input[256];
  char Output[256];
  float resX,resY,resZ;
  float sigmaX,sigmaY,sigmaZ;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  sigma = atof(argv[1]);
  argc--; argv++;
  
  //read the image
  SField.Read(Input);
  
  resX=sqrt(SField.Image2World[0][0]*SField.Image2World[0][0]+SField.Image2World[0][1]*SField.Image2World[0][1]+SField.Image2World[0][2]*SField.Image2World[0][2]);
  resY=sqrt(SField.Image2World[1][0]*SField.Image2World[1][0]+SField.Image2World[1][1]*SField.Image2World[1][1]+SField.Image2World[1][2]*SField.Image2World[1][2]);
  resZ=sqrt(SField.Image2World[2][0]*SField.Image2World[2][0]+SField.Image2World[2][1]*SField.Image2World[2][1]+SField.Image2World[2][2]*SField.Image2World[2][2]);
  
  cout << "Voxel size: " << resX << " " <<  resY << " " <<  resZ << endl;
  
  sigmaX=sigma/resX;
  sigmaY=sigma/resY;
  sigmaZ=sigma/resZ;
  
  
  //diffuse the image
  FFTconvolver.InitiateConvolver(SField.NX,SField.NY,SField.NZ,1,sigmaX,sigmaY,sigmaZ);
  FFTconvolver.Convolution(&SField);
  
  //write the image
  SField.Write(Output,Input);
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++        Erosion      ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Erosion(int argc, char **argv)
{
  ScalarField SField;
  ScalarField TempField;
  ScalarField Mask;
  int iterationNb,it;
  char Input[256];
  char Output[256];
  char MaskFile[256];
  int indicatorMask;
  float IdMask;
  int x,y,z,t;
  float tmpfl;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  iterationNb = atoi(argv[1]);
  
  //load the mask if defined
  indicatorMask=0;
  if (argc==4) {
    indicatorMask=1;
    argc--; argv++;
    strcpy(MaskFile,argv[1]);
    Mask.Read(MaskFile);
    argc--; argv++;
    IdMask = atof(argv[1]);
  }
  
  
  //read the image
  SField.Read(Input);
  TempField.Read(Input);
  
  
  
  //treat the image
  if (indicatorMask==0){//no mask
    for (it=0;it<iterationNb;it++) for (t=0;t<SField.NT;t++){
      //x direction
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++){
        for (x=1;x<SField.NX-1;x++) {
          tmpfl=SField.G(x,y,z,t);
          if (tmpfl>SField.G(x+1,y,z,t)) tmpfl=SField.G(x+1,y,z,t);
          if (tmpfl>SField.G(x-1,y,z,t)) tmpfl=SField.G(x-1,y,z,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //y direction
      for (z=0;z<SField.NZ;z++) for (x=0;x<SField.NX;x++){
        for (y=1;y<SField.NY-1;y++){
          tmpfl=TempField.G(x,y,z,t);
          if (tmpfl>TempField.G(x,y+1,z,t)) tmpfl=TempField.G(x,y+1,z,t);
          if (tmpfl>TempField.G(x,y-1,z,t)) tmpfl=TempField.G(x,y-1,z,t);
          SField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //z direction
      for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++){
        for (z=1;z<SField.NZ-1;z++){
          tmpfl=SField.G(x,y,z,t);
          if (tmpfl>SField.G(x,y,z+1,t)) tmpfl=SField.G(x,y,z+1,t);
          if (tmpfl>SField.G(x,y,z-1,t)) tmpfl=SField.G(x,y,z-1,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      }
      
      //copy result
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) SField.P(TempField.G(x,y,z,t),x,y,z,t);
    }
    
  }
  else {//with mask
    for (it=0;it<iterationNb;it++) for (t=0;t<SField.NT;t++){
      //x direction
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) {
        for (x=1;x<SField.NX-1;x++) if (abs(Mask.G(x,y,z)-IdMask)<0.01){
          tmpfl=SField.G(x,y,z,t);
          if (tmpfl>SField.G(x+1,y,z,t))  if (abs(Mask.G(x+1,y,z)-IdMask)<0.01) tmpfl=SField.G(x+1,y,z,t);
          if (tmpfl>SField.G(x-1,y,z,t))  if (abs(Mask.G(x-1,y,z)-IdMask)<0.01)tmpfl=SField.G(x-1,y,z,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //y direction
      for (z=0;z<SField.NZ;z++) for (x=0;x<SField.NX;x++){
        for (y=1;y<SField.NY-1;y++) if (abs(Mask.G(x,y,z)-IdMask)<0.01){
          tmpfl=TempField.G(x,y,z,t);
          if (tmpfl>TempField.G(x,y+1,z,t)) if (abs(Mask.G(x,y+1,z)-IdMask)<0.01) tmpfl=TempField.G(x,y+1,z,t);
          if (tmpfl>TempField.G(x,y-1,z,t)) if (abs(Mask.G(x,y-1,z)-IdMask)<0.01) tmpfl=TempField.G(x,y-1,z,t);
          SField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //z direction
      for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++){
        for (z=1;z<SField.NZ-1;z++) if (abs(Mask.G(x,y,z)-IdMask)<0.01){
          tmpfl=SField.G(x,y,z,t);
          if (tmpfl>SField.G(x,y,z+1,t)) if (abs(Mask.G(x,y,z+1)-IdMask)<0.01) tmpfl=SField.G(x,y,z+1,t);
          if (tmpfl>SField.G(x,y,z-1,t)) if (abs(Mask.G(x,y,z-1)-IdMask)<0.01) tmpfl=SField.G(x,y,z-1,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      }
      
      //copy result
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) SField.P(TempField.G(x,y,z,t),x,y,z,t);
    }
  }
  
  //write the image
  SField.Write(Output,Input);
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++        Dilation      ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Dilation(int argc, char **argv)
{
  ScalarField SField;
  ScalarField TempField;
  ScalarField Mask;
  int iterationNb,it;
  char Input[256];
  char Output[256];
  char MaskFile[256];
  int indicatorMask;
  float IdMask;
  int x,y,z,t;
  float tmpfl;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  iterationNb = atoi(argv[1]);
  
  //load the mask if defined
  indicatorMask=0;
  if (argc==4) {
    indicatorMask=1;
    argc--; argv++;
    strcpy(MaskFile,argv[1]);
    Mask.Read(MaskFile);
    argc--; argv++;
    IdMask = atof(argv[1]);
  }
  
  
  //read the image
  SField.Read(Input);
  TempField.Read(Input);
  
  
  
  //treat the image
  if (indicatorMask==0){//no mask
    for (it=0;it<iterationNb;it++) for (t=0;t<SField.NT;t++){
      //x direction
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++){
        for (x=0;x<SField.NX;x++) {
          tmpfl=SField.G(x,y,z,t);
          if (x!=SField.NX-1) if (tmpfl<SField.G(x+1,y,z,t)) tmpfl=SField.G(x+1,y,z,t);
          if (x!=0)           if (tmpfl<SField.G(x-1,y,z,t)) tmpfl=SField.G(x-1,y,z,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //y direction
      for (z=0;z<SField.NZ;z++) for (x=0;x<SField.NX;x++){
        for (y=0;y<SField.NY;y++){
          tmpfl=TempField.G(x,y,z,t);
          if (y!=SField.NY-1) if (tmpfl<TempField.G(x,y+1,z,t)) tmpfl=TempField.G(x,y+1,z,t);
          if (y!=0)           if (tmpfl<TempField.G(x,y-1,z,t)) tmpfl=TempField.G(x,y-1,z,t);
          SField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //z direction
      for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++){
        for (z=0;z<SField.NZ;z++){
          tmpfl=SField.G(x,y,z,t);
          if (z!=SField.NZ-1) if (tmpfl<SField.G(x,y,z+1,t)) tmpfl=SField.G(x,y,z+1,t);
          if (z!=0)           if (tmpfl<SField.G(x,y,z-1,t)) tmpfl=SField.G(x,y,z-1,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      }
      
      //copy result
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) SField.P(TempField.G(x,y,z,t),x,y,z,t);
    }
    
  }
  else {//with mask
    for (it=0;it<iterationNb;it++) for (t=0;t<SField.NT;t++){
      //x direction
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) {
        for (x=1;x<SField.NX-1;x++) if (abs(Mask.G(x,y,z)-IdMask)<0.01){
          tmpfl=SField.G(x,y,z,t);
          if (x!=SField.NX-1) if (tmpfl<SField.G(x+1,y,z,t))  if (abs(Mask.G(x+1,y,z)-IdMask)<0.01) tmpfl=SField.G(x+1,y,z,t);
          if (x!=0)           if (tmpfl<SField.G(x-1,y,z,t))  if (abs(Mask.G(x-1,y,z)-IdMask)<0.01)tmpfl=SField.G(x-1,y,z,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //y direction
      for (z=0;z<SField.NZ;z++) for (x=0;x<SField.NX;x++){
        for (y=0;y<SField.NY;y++) if (abs(Mask.G(x,y,z)-IdMask)<0.01){
          tmpfl=TempField.G(x,y,z,t);
          if (y!=SField.NY-1) if (tmpfl<TempField.G(x,y+1,z,t)) if (abs(Mask.G(x,y+1,z)-IdMask)<0.01) tmpfl=TempField.G(x,y+1,z,t);
          if (y!=0)           if (tmpfl<TempField.G(x,y-1,z,t)) if (abs(Mask.G(x,y-1,z)-IdMask)<0.01) tmpfl=TempField.G(x,y-1,z,t);
          SField.P(tmpfl,x,y,z,t);
        }
      } 
      
      //z direction
      for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++){
        for (z=0;z<SField.NZ;z++) if (abs(Mask.G(x,y,z)-IdMask)<0.01){
          tmpfl=SField.G(x,y,z,t);
          if (z!=SField.NZ-1) if (tmpfl<SField.G(x,y,z+1,t)) if (abs(Mask.G(x,y,z+1)-IdMask)<0.01) tmpfl=SField.G(x,y,z+1,t);
          if (z!=0)           if (tmpfl<SField.G(x,y,z-1,t)) if (abs(Mask.G(x,y,z-1)-IdMask)<0.01) tmpfl=SField.G(x,y,z-1,t);
          TempField.P(tmpfl,x,y,z,t);
        }
      }
      
      //copy result
      for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) SField.P(TempField.G(x,y,z,t),x,y,z,t);
    }
  }
  
  
  //for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++){
  //   if (z<104) SField.P(1,x,y,z,t);
  //    else SField.P(0,x,y,z,t);
  //}
  
  
  //write the image
  SField.Write(Output,Input);
  
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++       Add noise     ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void AddNoise(int argc, char **argv)
{
  ScalarField SField;
  char Input[256];
  char Output[256];
  int x,y,z,t;
  float NoiseLevel,tmpfl;
  double u,v,r,c;
  int indicatorROI;
  int Xmin,Xmax,Ymin,Ymax,Zmin,Zmax,tmpInt;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  NoiseLevel = atof(argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  
  //load the ROI parameters if defined
  indicatorROI=0;
  cout << argc << endl;
  if (argc==8) {
    indicatorROI=1;
    argc--; argv++;
    Xmin = atoi(argv[1]);
    argc--; argv++;
    Xmax = atoi(argv[1]);
    argc--; argv++;
    Ymin = atoi(argv[1]);
    argc--; argv++;
    Ymax = atoi(argv[1]);
    argc--; argv++;
    Zmin = atoi(argv[1]);
    argc--; argv++;
    Zmax = atoi(argv[1]);
  }
  
  //read the image
  SField.Read(Input);

  //define the ROI
  if (indicatorROI==0){
    Xmin=0;
    Ymin=0;
    Zmin=0;
    Xmax=SField.NX;
    Ymax=SField.NY;
    Zmax=SField.NZ;
  }
  else{
    if (Xmin>Xmax){tmpInt=Xmax; Xmax=Xmin; Xmin=tmpInt;}
    if (Ymin>Ymax){tmpInt=Ymax; Ymax=Ymin; Ymin=tmpInt;}
    if (Zmin>Zmax){tmpInt=Zmax; Zmax=Zmin; Zmin=tmpInt;}
    if (Xmin<0) Xmin=0;
    if (Ymin<0) Ymin=0;
    if (Zmin<0) Zmin=0;
    if (Xmax>SField.NX) Xmax=SField.NX;
    if (Ymax>SField.NY) Ymax=SField.NY;
    if (Zmax>SField.NZ) Zmax=SField.NZ;
    }

  cout << "ROI:" << Xmin << " " << Xmax << " " << Ymin << " " << Ymax << " " << Zmin << " " << Zmax << endl;
  
  //treat the image
  srand (time(NULL));
  for (t=0;t<SField.NT;t++) for (z=Zmin;z<Zmax;z++) for (y=Ymin;y<Ymax;y++) for (x=Xmin;x<Xmax;x++) {
    //generate the noise
    u = (static_cast<double>(rand()) / (RAND_MAX)) * 2 - 1;
    v = (static_cast<double>(rand()) / (RAND_MAX)) * 2 - 1;
    r = u * u + v * v;
    if (r<0.1) r=0.1;
    if (r>0.9) r=0.9;
    c = sqrt(-2 * log(r) / r);
    tmpfl=NoiseLevel*(static_cast<float>(u * c));
    
    //noise the intensity
    SField.P(SField.G(x,y,z,t)+tmpfl,x,y,z,t);
  }
  
  //write the image
  SField.Write(Output,Input);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++     mask an image   ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MaskImage(int argc, char **argv)
{
  ScalarField SField;
  ScalarField Mask;
  char Input[256];
  char Output[256];
  char MaskFile[256];
  float IdMask;
  int x,y,z,t;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  strcpy(MaskFile,argv[1]);
  Mask.Read(MaskFile);
  argc--; argv++;
  IdMask = atof(argv[1]);
  
  
  //read the image
  SField.Read(Input);
  
  for (t=0;t<SField.NT;t++)  for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) if (abs(Mask.G(x,y,z)-IdMask)<0.01) SField.P(0,x,y,z,t);
  
  //write the image
  SField.Write(Output,Input);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++     max grey level   ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MaxGreyLev(int argc, char **argv)
{
  ScalarField Input1;
  ScalarField Output;
  char InputFile1[256];
  char InputFile2[256];
  char OutputFile[256];
  int x,y,z,t;
  
  //read parameters
  argc--; argv++;
  strcpy(InputFile1,argv[1]);
  argc--; argv++;
  strcpy(InputFile2,argv[1]);
  argc--; argv++;
  strcpy(OutputFile,argv[1]);
  
  //read the images
  Input1.Read(InputFile1);
  Output.Read(InputFile2);
  
  //do the job
  for (t=0;t<Output.NT;t++)  for (z=0;z<Output.NZ;z++) for (y=0;y<Output.NY;y++) for (x=0;x<Output.NX;x++) 
    if (Input1.G(x,y,z,t)>Output.G(x,y,z,t)) Output.P(Input1.G(x,y,z,t),x,y,z,t);
  
  //write the image
  Output.Write(OutputFile,InputFile1);
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++  cpt the norm of the gradient in an image +++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void NormGrad(int argc, char **argv)
{
  ScalarField SField;
  ScalarField SField2;
  char Input[256];
  char Output[256];
  int x,y,z,t;
  float temp;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  
  
  //read the image
  SField.Read(Input);
  SField2.Read(Input);
  
  for (t=0;t<SField.NT;t++)  for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) SField2.P(0,x,y,z,t);
  
  for (t=0;t<SField.NT;t++)  for (z=1;z<SField.NZ-1;z++) for (y=1;y<SField.NY-1;y++) for (x=1;x<SField.NX-1;x++){
    temp= (SField.G(x+1,y,z,t)-SField.G(x-1,y,z,t))*(SField.G(x+1,y,z,t)-SField.G(x-1,y,z,t));
    temp+=(SField.G(x,y+1,z,t)-SField.G(x,y-1,z,t))*(SField.G(x,y+1,z,t)-SField.G(x,y-1,z,t));
    temp+=(SField.G(x,y,z+1,t)-SField.G(x,y,z-1,t))*(SField.G(x,y,z+1,t)-SField.G(x,y,z-1,t));
    temp=sqrt(temp);
    SField2.P(temp,x,y,z,t);
  }
  
  //write the image
  SField2.Write(Output,Input);
  
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++  transform the image2world values of an image to be centered on another image +++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RegCenterImage(int argc, char **argv)
{
  ScalarField SFieldMov;
  ScalarField SFieldFix;
  char InputMov[256];
  char InputFix[256];
  char Output[256];
  float cxm,cym,czm,cxf,cyf,czf;
  
  
  //read parameters
  argc--; argv++;
  strcpy(InputMov,argv[1]);
  argc--; argv++;
  strcpy(InputFix,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  
  
  //read the image
  SFieldMov.Read(InputMov);
  SFieldFix.Read(InputFix);
  
  //center of the moving and fixed image
  cxm=(SFieldMov.NX/2)*SFieldMov.Image2World[0][0]+(SFieldMov.NY/2)*SFieldMov.Image2World[0][1]+(SFieldMov.NZ/2)*SFieldMov.Image2World[0][2]+SFieldMov.Image2World[0][3];
  cym=(SFieldMov.NX/2)*SFieldMov.Image2World[1][0]+(SFieldMov.NY/2)*SFieldMov.Image2World[1][1]+(SFieldMov.NZ/2)*SFieldMov.Image2World[1][2]+SFieldMov.Image2World[1][3];
  czm=(SFieldMov.NX/2)*SFieldMov.Image2World[2][0]+(SFieldMov.NY/2)*SFieldMov.Image2World[2][1]+(SFieldMov.NZ/2)*SFieldMov.Image2World[2][2]+SFieldMov.Image2World[2][3];
  
  cxf=(SFieldFix.NX/2)*SFieldFix.Image2World[0][0]+(SFieldFix.NY/2)*SFieldFix.Image2World[0][1]+(SFieldFix.NZ/2)*SFieldFix.Image2World[0][2]+SFieldFix.Image2World[0][3];
  cyf=(SFieldFix.NX/2)*SFieldFix.Image2World[1][0]+(SFieldFix.NY/2)*SFieldFix.Image2World[1][1]+(SFieldFix.NZ/2)*SFieldFix.Image2World[1][2]+SFieldFix.Image2World[1][3];
  czf=(SFieldFix.NX/2)*SFieldFix.Image2World[2][0]+(SFieldFix.NY/2)*SFieldFix.Image2World[2][1]+(SFieldFix.NZ/2)*SFieldFix.Image2World[2][2]+SFieldFix.Image2World[2][3];
  
  //image translation to have the same image center
  SFieldMov.Image2World[0][3]+=cxf-cxm;
  SFieldMov.Image2World[1][3]+=cyf-cym;
  SFieldMov.Image2World[2][3]+=czf-czm;
  
  cout << "The translation of " << InputMov << " is (" << cxf-cxm << ", " << cyf-cym << ", "  << czf-czm << ") mm." << endl; 
  
  //write the image
  SFieldMov.Write(Output);
  
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++                       make a sequence of 3D or 2D images                      +++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void MakeSequence(int argc, char **argv){
  int TimeFrameNb,i;
  ScalarField InputFrame;
  ScalarField OutputSequence;
  char InputFrameFile[256];
  char OutputSequenceFile[256];
  char Output[256];
  float cxm,cym,czm,cxf,cyf,czf;
  
  
  //1) read the number of time frames and the first frame of the sequence
  argc--; argv++;
  TimeFrameNb = atoi(argv[1]);
  argc--; argv++;
  strcpy(InputFrameFile,argv[1]);
  
  cout << "Import "  << InputFrameFile << endl;

  InputFrame.Read(InputFrameFile);
  
  OutputSequence.CreateVoidField(InputFrame.NX,InputFrame.NY,InputFrame.NZ,TimeFrameNb);
  
  OutputSequence.Image2World[0][0]=InputFrame.Image2World[0][0]; OutputSequence.Image2World[0][1]=InputFrame.Image2World[0][1]; OutputSequence.Image2World[0][2]=InputFrame.Image2World[0][2]; OutputSequence.Image2World[0][3]=InputFrame.Image2World[0][3];  
  OutputSequence.Image2World[1][0]=InputFrame.Image2World[1][0]; OutputSequence.Image2World[1][1]=InputFrame.Image2World[1][1]; OutputSequence.Image2World[1][2]=InputFrame.Image2World[1][2]; OutputSequence.Image2World[1][3]=InputFrame.Image2World[1][3];  
  OutputSequence.Image2World[2][0]=InputFrame.Image2World[2][0]; OutputSequence.Image2World[2][1]=InputFrame.Image2World[2][1]; OutputSequence.Image2World[2][2]=InputFrame.Image2World[2][2]; OutputSequence.Image2World[2][3]=InputFrame.Image2World[2][3];
  OutputSequence.Image2World[3][0]=InputFrame.Image2World[3][0]; OutputSequence.Image2World[3][1]=InputFrame.Image2World[3][1]; OutputSequence.Image2World[3][2]=InputFrame.Image2World[3][2]; OutputSequence.Image2World[3][3]=InputFrame.Image2World[3][3];

  OutputSequence.World2Image[0][0]=InputFrame.World2Image[0][0]; OutputSequence.World2Image[0][1]=InputFrame.World2Image[0][1]; OutputSequence.World2Image[0][2]=InputFrame.World2Image[0][2]; OutputSequence.World2Image[0][3]=InputFrame.World2Image[0][3];  
  OutputSequence.World2Image[1][0]=InputFrame.World2Image[1][0]; OutputSequence.World2Image[1][1]=InputFrame.World2Image[1][1]; OutputSequence.World2Image[1][2]=InputFrame.World2Image[1][2]; OutputSequence.World2Image[1][3]=InputFrame.World2Image[1][3];  
  OutputSequence.World2Image[2][0]=InputFrame.World2Image[2][0]; OutputSequence.World2Image[2][1]=InputFrame.World2Image[2][1]; OutputSequence.World2Image[2][2]=InputFrame.World2Image[2][2]; OutputSequence.World2Image[2][3]=InputFrame.World2Image[2][3];
  OutputSequence.World2Image[3][0]=InputFrame.World2Image[3][0]; OutputSequence.World2Image[3][1]=InputFrame.World2Image[3][1]; OutputSequence.World2Image[3][2]=InputFrame.World2Image[3][2]; OutputSequence.World2Image[3][3]=InputFrame.World2Image[3][3];
  
  DeepCopy(&InputFrame,&OutputSequence,0);
  
  //2) read and treat other frames
  for (i=1;i<TimeFrameNb;i++){
    argc--; argv++;
    strcpy(InputFrameFile,argv[1]);
    cout << "Import "  << InputFrameFile << endl;

    InputFrame.Read(InputFrameFile);
    DeepCopy(&InputFrame,&OutputSequence,i);
  }
  
  //3) save the sequence
  argc--; argv++;
  strcpy(OutputSequenceFile,argv[1]);
  
  cout << "Save "  << OutputSequenceFile << endl;
  OutputSequence.Write(OutputSequenceFile);
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++     mask an image   ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RemoveSmallReg(int argc, char **argv)
{
  ScalarField SField;
  ScalarField TestField;
  char Input[256];
  char Output[256];
  float IdMask;
  float IdBackG;
  int MaxSize;
  int x,y,z,i,j,x1,y1,z1;
  float tempfl;
  int * ListVoxX;
  int * ListVoxY;
  int * ListVoxZ;
  int * ListVoxX1;
  int * ListVoxY1;
  int * ListVoxZ1;
  int * ListVoxXtot;
  int * ListVoxYtot;
  int * ListVoxZtot;
  int NbVoxListVox,NbVoxListVox1,tempint;
  int Changes;
  int NbPtsRegion;
  int NbRegions;
  int NbRegionRemoved;
  int NbVoxelsRemoved;
  
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  IdMask = atof(argv[1]);
  argc--; argv++;
  IdBackG = atof(argv[1]);
  argc--; argv++;
  MaxSize = atoi(argv[1]);
  
  
  //read the image
  SField.Read(Input);
  
  //initiate the field containing the labels
  TestField.Read(Input);
  for (z=0;z<TestField.NZ;z++) for (y=0;y<TestField.NY;y++) for (x=0;x<TestField.NX;x++) TestField.P(0,x,y,z);
  
  NbRegions=0;
  NbRegionRemoved=0;
  NbVoxelsRemoved=0;
  
  
  //initiate the voxel lists
  ListVoxX = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxY = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxZ = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxX1 = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxY1 = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxZ1 = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxXtot = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxYtot = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxZtot = new int[SField.NZ*SField.NY*SField.NX];
  
  
  //big loop on the seed voxel
  tempfl=1;
  for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) if ((abs(SField.G(x,y,z)-IdMask)<0.01)&&(abs(TestField.G(x,y,z))<0.01)){
    
    //initiate the label propagation
    NbVoxListVox=0;
    NbPtsRegion=0;
    
    
    TestField.P(tempfl,x,y,z);
    ListVoxX[NbVoxListVox]=x;
    ListVoxY[NbVoxListVox]=y;
    ListVoxZ[NbVoxListVox]=z;
    ListVoxXtot[NbPtsRegion]=x;
    ListVoxYtot[NbPtsRegion]=y;
    ListVoxZtot[NbPtsRegion]=z;
    NbVoxListVox++;
    NbPtsRegion++;
    Changes=1;
    
    while (Changes==1){
      //propagate
      Changes=0;
      NbVoxListVox1=0;
      for (i=0;i<NbVoxListVox;i++) for (j=0;j<6;j++){
        if (j==0) {x1=ListVoxX[i]+1; y1=ListVoxY[i];   z1=ListVoxZ[i];   }
        if (j==1) {x1=ListVoxX[i]-1; y1=ListVoxY[i];   z1=ListVoxZ[i];   }
        if (j==2) {x1=ListVoxX[i];   y1=ListVoxY[i]+1; z1=ListVoxZ[i];   }
        if (j==3) {x1=ListVoxX[i];   y1=ListVoxY[i]-1; z1=ListVoxZ[i];   }
        if (j==4) {x1=ListVoxX[i];   y1=ListVoxY[i];   z1=ListVoxZ[i]+1; }
        if (j==5) {x1=ListVoxX[i];   y1=ListVoxY[i];   z1=ListVoxZ[i]-1; }
        
        if ((x1>0)&&(y1>0)&&(z1>0)&&(x1<TestField.NX)&&(y1<TestField.NY)&&(z1<TestField.NZ))
          if ((abs(SField.G(x1,y1,z1)-IdMask)<0.01)&&(abs(TestField.G(x1,y1,z1))<0.01)){
            TestField.P(tempfl,x1,y1,z1);
            ListVoxX1[NbVoxListVox1]=x1;
            ListVoxY1[NbVoxListVox1]=y1;
            ListVoxZ1[NbVoxListVox1]=z1;
            ListVoxXtot[NbPtsRegion]=x1;
            ListVoxYtot[NbPtsRegion]=y1;
            ListVoxZtot[NbPtsRegion]=z1;
            NbVoxListVox1++;
            NbPtsRegion++;
            Changes=1;
          }
      }
      //update the list of the points at the boundary of the propagation
      for (i=0;i<NbVoxListVox1;i++){
        ListVoxX[i]=ListVoxX1[i];
        ListVoxY[i]=ListVoxY1[i];
        ListVoxZ[i]=ListVoxZ1[i];
      }
      NbVoxListVox=NbVoxListVox1;
    }
    
    NbRegions++;
    
    //remove the region if to small
    if (NbPtsRegion<MaxSize){
      for (i=0;i<NbPtsRegion;i++) SField.P(IdBackG,ListVoxXtot[i],ListVoxYtot[i],ListVoxZtot[i]);
      
      NbRegionRemoved++;
      NbVoxelsRemoved+=NbPtsRegion;
    }
    
    
    
    
    tempfl++;
    
  }
  
  cout << NbRegionRemoved << " regions removed out of " << NbRegions << "  ("  << NbVoxelsRemoved << " voxels)" << endl;
  
  
  //write the image
  SField.Write(Output,Input);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++     mask an image   ++++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void KeepLargestReg(int argc, char **argv)
{
  ScalarField SField;
  ScalarField TestField;
  char Input[256];
  char Output[256];
  float IdMask;
  float IdBackG;
  int MaxSize;
  int x,y,z,i,j,x1,y1,z1;
  float tempfl;
  int * ListVoxX;
  int * ListVoxY;
  int * ListVoxZ;
  int * ListVoxX1;
  int * ListVoxY1;
  int * ListVoxZ1;
  int * ListVoxXtot;
  int * ListVoxYtot;
  int * ListVoxZtot;
  int NbVoxListVox,NbVoxListVox1,tempint;
  int Changes;
  int NbPtsRegion;
  int NbRegions;
  int NbRegionRemoved;
  int NbVoxelsRemoved;
  
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  IdMask = atof(argv[1]);
  argc--; argv++;
  IdBackG = atof(argv[1]);
  argc--; argv++;
  
  //1) Init
  
  //read the image
  SField.Read(Input);
  
  
  //2) Identify the largest connected set of voxels
  MaxSize=0;
  //initiate the field containing the labels
  TestField.Read(Input);
  for (z=0;z<TestField.NZ;z++) for (y=0;y<TestField.NY;y++) for (x=0;x<TestField.NX;x++) TestField.P(0,x,y,z);
  
  //initiate the voxel lists
  ListVoxX = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxY = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxZ = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxX1 = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxY1 = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxZ1 = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxXtot = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxYtot = new int[SField.NZ*SField.NY*SField.NX];
  ListVoxZtot = new int[SField.NZ*SField.NY*SField.NX];
  
  
  //big loop on the seed voxel
  tempfl=1;
  for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) if ((abs(SField.G(x,y,z)-IdMask)<0.01)&&(abs(TestField.G(x,y,z))<0.01)){
    
    //initiate the label propagation
    NbVoxListVox=0;
    NbPtsRegion=0;
    
    
    TestField.P(tempfl,x,y,z);
    ListVoxX[NbVoxListVox]=x;
    ListVoxY[NbVoxListVox]=y;
    ListVoxZ[NbVoxListVox]=z;
    ListVoxXtot[NbPtsRegion]=x;
    ListVoxYtot[NbPtsRegion]=y;
    ListVoxZtot[NbPtsRegion]=z;
    NbVoxListVox++;
    NbPtsRegion++;
    Changes=1;
    
    while (Changes==1){
      //propagate
      Changes=0;
      NbVoxListVox1=0;
      for (i=0;i<NbVoxListVox;i++) for (j=0;j<6;j++){
        if (j==0) {x1=ListVoxX[i]+1; y1=ListVoxY[i];   z1=ListVoxZ[i];   }
        if (j==1) {x1=ListVoxX[i]-1; y1=ListVoxY[i];   z1=ListVoxZ[i];   }
        if (j==2) {x1=ListVoxX[i];   y1=ListVoxY[i]+1; z1=ListVoxZ[i];   }
        if (j==3) {x1=ListVoxX[i];   y1=ListVoxY[i]-1; z1=ListVoxZ[i];   }
        if (j==4) {x1=ListVoxX[i];   y1=ListVoxY[i];   z1=ListVoxZ[i]+1; }
        if (j==5) {x1=ListVoxX[i];   y1=ListVoxY[i];   z1=ListVoxZ[i]-1; }
        
        if ((x1>0)&&(y1>0)&&(z1>0)&&(x1<TestField.NX)&&(y1<TestField.NY)&&(z1<TestField.NZ))
          if ((abs(SField.G(x1,y1,z1)-IdMask)<0.01)&&(abs(TestField.G(x1,y1,z1))<0.01)){
            TestField.P(tempfl,x1,y1,z1);
            ListVoxX1[NbVoxListVox1]=x1; ListVoxY1[NbVoxListVox1]=y1; ListVoxZ1[NbVoxListVox1]=z1;
            ListVoxXtot[NbPtsRegion]=x1; ListVoxYtot[NbPtsRegion]=y1; ListVoxZtot[NbPtsRegion]=z1;
            NbVoxListVox1++; NbPtsRegion++; Changes=1;
          }
      }
      //update the list of the points at the boundary of the propagation
      for (i=0;i<NbVoxListVox1;i++){
        ListVoxX[i]=ListVoxX1[i]; ListVoxY[i]=ListVoxY1[i]; ListVoxZ[i]=ListVoxZ1[i];
      }
      NbVoxListVox=NbVoxListVox1;
    }
    
    //update MaxSize if necessary
    if (NbPtsRegion>MaxSize)
      MaxSize=NbPtsRegion;
    
    tempfl++;
    
  }
  
  //3) remove all the connected set of voxels expect the largest one
  
  //initiate the field containing the labels
  
  for (z=0;z<TestField.NZ;z++) for (y=0;y<TestField.NY;y++) for (x=0;x<TestField.NX;x++) TestField.P(0,x,y,z);
  
  NbRegions=0;
  NbRegionRemoved=0;
  NbVoxelsRemoved=0;
  
  
  //big loop on the seed voxel
  tempfl=1;
  for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++) if ((abs(SField.G(x,y,z)-IdMask)<0.01)&&(abs(TestField.G(x,y,z))<0.01)){
    
    //initiate the label propagation
    NbVoxListVox=0;
    NbPtsRegion=0;
    
    
    TestField.P(tempfl,x,y,z);
    ListVoxX[NbVoxListVox]=x;
    ListVoxY[NbVoxListVox]=y;
    ListVoxZ[NbVoxListVox]=z;
    ListVoxXtot[NbPtsRegion]=x;
    ListVoxYtot[NbPtsRegion]=y;
    ListVoxZtot[NbPtsRegion]=z;
    NbVoxListVox++;
    NbPtsRegion++;
    Changes=1;
    
    while (Changes==1){
      //propagate
      Changes=0;
      NbVoxListVox1=0;
      for (i=0;i<NbVoxListVox;i++) for (j=0;j<6;j++){
        if (j==0) {x1=ListVoxX[i]+1; y1=ListVoxY[i];   z1=ListVoxZ[i];   }
        if (j==1) {x1=ListVoxX[i]-1; y1=ListVoxY[i];   z1=ListVoxZ[i];   }
        if (j==2) {x1=ListVoxX[i];   y1=ListVoxY[i]+1; z1=ListVoxZ[i];   }
        if (j==3) {x1=ListVoxX[i];   y1=ListVoxY[i]-1; z1=ListVoxZ[i];   }
        if (j==4) {x1=ListVoxX[i];   y1=ListVoxY[i];   z1=ListVoxZ[i]+1; }
        if (j==5) {x1=ListVoxX[i];   y1=ListVoxY[i];   z1=ListVoxZ[i]-1; }
        
        if ((x1>0)&&(y1>0)&&(z1>0)&&(x1<TestField.NX)&&(y1<TestField.NY)&&(z1<TestField.NZ))
          if ((abs(SField.G(x1,y1,z1)-IdMask)<0.01)&&(abs(TestField.G(x1,y1,z1))<0.01)){
            TestField.P(tempfl,x1,y1,z1);
            ListVoxX1[NbVoxListVox1]=x1; ListVoxY1[NbVoxListVox1]=y1; ListVoxZ1[NbVoxListVox1]=z1;
            ListVoxXtot[NbPtsRegion]=x1; ListVoxYtot[NbPtsRegion]=y1; ListVoxZtot[NbPtsRegion]=z1;
            NbVoxListVox1++; NbPtsRegion++; Changes=1;
          }
      }
      //update the list of the points at the boundary of the propagation
      for (i=0;i<NbVoxListVox1;i++){
        ListVoxX[i]=ListVoxX1[i]; ListVoxY[i]=ListVoxY1[i]; ListVoxZ[i]=ListVoxZ1[i];
      }
      NbVoxListVox=NbVoxListVox1;
    }
    
    NbRegions++;
    
    //remove the region if to small
    if (NbPtsRegion<MaxSize){
      for (i=0;i<NbPtsRegion;i++) SField.P(IdBackG,ListVoxXtot[i],ListVoxYtot[i],ListVoxZtot[i]);
      
      NbRegionRemoved++;
      NbVoxelsRemoved+=NbPtsRegion;
    }
    
    tempfl++;
    
  }
  
  cout << NbRegionRemoved << " regions removed out of " << NbRegions << "  ("  << NbVoxelsRemoved << " voxels)" << endl;
  
  
  //write the image
  SField.Write(Output,Input);
  
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++     Put a grey level in a ROI   ++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PutGLInROI(int argc, char **argv)
{
  ScalarField SField;
  char Input[256];
  char Output[256];
  char MaskFile[256];
  float GreyLev;
  int x,y,z,t;
  int Rx1,Rx2,Ry1,Ry2,Rz1,Rz2;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  Rx1 = atoi(argv[1]);
  argc--; argv++;
  Rx2 = atoi(argv[1]);
  argc--; argv++;
  Ry1 = atoi(argv[1]);
  argc--; argv++;
  Ry2 = atoi(argv[1]);
  argc--; argv++;
  Rz1 = atoi(argv[1]);
  argc--; argv++;
  Rz2 = atoi(argv[1]);
  argc--; argv++;
  GreyLev = atof(argv[1]);
  
  
  //read the image
  SField.Read(Input);
  
  
  //treat the negative values
  if (Rx1<0) {cout << Rx1 << " becomes " << SField.NX+Rx1 << endl; Rx1=SField.NX+Rx1;}
  if (Rx2<0) {cout << Rx2 << " becomes " << SField.NX+Rx2 << endl; Rx2=SField.NX+Rx2;}
  if (Ry1<0) {cout << Ry1 << " becomes " << SField.NY+Ry1 << endl; Ry1=SField.NY+Ry1;}
  if (Ry2<0) {cout << Ry2 << " becomes " << SField.NY+Ry2 << endl; Ry2=SField.NY+Ry2;}
  if (Rz1<0) {cout << Rz1 << " becomes " << SField.NZ+Rz1 << endl; Rz1=SField.NZ+Rz1;}
  if (Rz2<0) {cout << Rz2 << " becomes " << SField.NZ+Rz2 << endl; Rz2=SField.NZ+Rz2;}
  
  //do the job
  for (t=0;t<SField.NT;t++)  for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++)
    if ((z>=Rz1)&&(z<=Rz2)&&(y>=Ry1)&&(y<=Ry2)&&(x>=Rx1)&&(x<=Rx2))
      SField.P(GreyLev,x,y,z,t);
  
  //write the image
  SField.Write(Output,Input);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++     Put a grey level in a ROI   ++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void PutGLInROI2(int argc, char **argv)
{
  ScalarField SField;
  char Input[256];
  char Output[256];
  float GreyLevIn,GreyLevOut;
  int x,y,z,t;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  GreyLevIn = atof(argv[1]);
  argc--; argv++;
  GreyLevOut = atof(argv[1]);
  argc--; argv++;
  
  
  //read the image
  SField.Read(Input);
  
  
  //do the job
  for (t=0;t<SField.NT;t++)  for (z=0;z<SField.NZ;z++) for (y=0;y<SField.NY;y++) for (x=0;x<SField.NX;x++)
    if (fabs(SField.G(x,y,z,t)-GreyLevIn)<0.0001)
      SField.P(GreyLevOut,x,y,z,t);
  
  //write the image
  SField.Write(Output,Input);
  
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++++ Undersample an image +++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void UnderSample(int argc, char **argv){
  ScalarField SField;
  ScalarField SField2;
  float factor;
  int NX2,NY2,NZ2;
  char Input[256];
  char Output[256];
  int optionNN;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  factor = atof(argv[1]);
  
  optionNN=0;
  if (argc>1) optionNN=1;

  //read the image
  SField.Read(Input);
  
  //cpt the new resolution
  NX2=static_cast<int>(SField.NX/factor);
  NY2=static_cast<int>(SField.NY/factor);
  if (SField.NZ>1) NZ2=static_cast<int>(SField.NZ/factor);
  else  NZ2=static_cast<int>(SField.NZ);
  
  //read again the image but undersampled
  if (optionNN==0) SField2.Read_and_Interpolate(Input,NX2,NY2,NZ2);
  else SField2.Read_and_Undersample(Input,factor,1);
  
  //save the undersampled image
  SField2.Write(Output);
  
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++++ grey level alignment +++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GreyLevAlign(int argc, char **argv){
  ScalarField ImToAlign;
  ScalarField RefImage;
  char Input1[256];
  char Input2[256];
  char Output[256];
  int UsingTwoModes;
  
  //read parameters
  argc--; argv++;
  strcpy(Input1,argv[1]);
  argc--; argv++;
  strcpy(Input2,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  
  //read the images
  ImToAlign.Read(Input1);
  RefImage.Read(Input2);
  
  //grey level alignment
  ImToAlign.GreyLevAlignment(&RefImage);
  
  //save the undersampled image
  ImToAlign.Write(Output);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++++++ Undersample a very large image ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void LargeImageUndersample(int argc, char **argv){
  ScalarField SField;
  int factor;
  char Input[256];
  char Output[256];
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  factor = atoi(argv[1]);
  
  //read and undersample the image
  SField.Read_directly_Undersampled(Input,factor);
  
  //save the undersampled image
  SField.Write(Output);
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++ Cut a region of interest +++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

//cerr << "-ROI_Cut [ImageIn][ImageOut] [MinX][MaxX] [MinY][MaxY] [MinZ][MaxZ] [MinT][MaxT]\n";

void ROI_Cut(int argc, char **argv){
  ScalarField SField;
  int MinX,MaxX,MinY,MaxY,MinZ,MaxZ,MinT,MaxT;
  char Input[256];
  char Output[256];
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  MinX = atoi(argv[1]);
  argc--; argv++;
  MaxX = atoi(argv[1]);
  argc--; argv++;
  MinY = atoi(argv[1]);
  argc--; argv++;
  MaxY = atoi(argv[1]);
  argc--; argv++;
  MinZ = atoi(argv[1]);
  argc--; argv++;
  MaxZ = atoi(argv[1]);
  argc--; argv++;
  MinT = atoi(argv[1]);
  argc--; argv++;
  MaxT = atoi(argv[1]);
  
  
  //read the image
  SField.Read_only_ROI(Input,MinX,MaxX,MinY,MaxY,MinZ,MaxZ,MinT,MaxT);
  
  //save the ROI
  SField.Write(Output);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++ threshold of the grey levels +++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GreyLevThresh(int argc, char **argv){
  ScalarField SField;
  float minGL,maxGL,mi,ma;
  char Input[256];
  char Output[256];
  int i,j,k;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  minGL = atof(argv[1]);
  argc--; argv++;
  maxGL = atof(argv[1]);
  
  
  //read the image
  SField.Read(Input);
  
  
  
  //do the job
  if (abs(minGL-maxGL)<0.01){
    mi=SField.G(0,0,0);
    ma=SField.G(0,0,0);
    for (i=0;i<SField.NZ;i++) for (j=0;j<SField.NY;j++) for (k=0;k<SField.NX;k++) {
      if (SField.G(k,j,i)<mi) mi=SField.G(k,j,i);
      if (SField.G(k,j,i)>ma) ma=SField.G(k,j,i);
    }
    for (i=0;i<SField.NZ;i++) for (j=0;j<SField.NY;j++) for (k=0;k<SField.NX;k++) {
      if (SField.G(k,j,i)<minGL) SField.P(mi,k,j,i);
      else SField.P(ma,k,j,i);
    }
  }
  else{
    for (i=0;i<SField.NZ;i++) for (j=0;j<SField.NY;j++) for (k=0;k<SField.NX;k++) {
      if (SField.G(k,j,i)<minGL) SField.P(minGL,k,j,i);
      if (SField.G(k,j,i)>maxGL) SField.P(maxGL,k,j,i);
    }
  }
  
  //save the undersampled image
  SField.Write(Output);
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++ linearly resample the grey levels +++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GreyLevResample(int argc, char **argv){
  ScalarField SField;
  float minGL,maxGL;
  char Input[256];
  char Output[256];
  int i,j,k;
  float NewMinGL,NewMaxGL;
  float a,b,tmp;
  
  //read parameters
  argc--; argv++;
  strcpy(Input,argv[1]);
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  NewMinGL = atof(argv[1]);
  argc--; argv++;
  NewMaxGL = atof(argv[1]);
  
  
  //read the image
  SField.Read(Input);
  
  //do the job
  minGL=SField.G(0,0,0);
  maxGL=SField.G(0,0,0);
  for (i=0;i<SField.NZ;i++) for (j=0;j<SField.NY;j++) for (k=0;k<SField.NX;k++) {
    if (SField.G(k,j,i)<minGL) minGL=SField.G(k,j,i);
    if (SField.G(k,j,i)>maxGL) maxGL=SField.G(k,j,i);
  }
  
  a=(NewMinGL-NewMaxGL)/(minGL-maxGL);
  b=NewMaxGL-a*maxGL;
  
  if (minGL==maxGL){
    a=0;
    b=NewMinGL;
  }
  
  for (i=0;i<SField.NZ;i++) for (j=0;j<SField.NY;j++) for (k=0;k<SField.NX;k++){
    tmp=a*SField.G(k,j,i)+b;
    SField.P(tmp,k,j,i);
  }
  
  //post treatment in case there are minor errors when we want to binarise an image
  for (i=0;i<SField.NZ;i++) for (j=0;j<SField.NY;j++) for (k=0;k<SField.NX;k++){
    if (fabs(SField.G(k,j,i)-NewMinGL)<0.0001) SField.P(NewMinGL,k,j,i);
    if (fabs(SField.G(k,j,i)-NewMaxGL)<0.0001) SField.P(NewMaxGL,k,j,i);
  }
  
  //for (i=0;i<SField.NZ;i++) for (j=0;j<SField.NY;j++) for (k=0;k<SField.NX;k++)
  //  if ((i<30)||(i>SField.NZ-12)||(j<12)||(j>SField.NY-12)||(k<12)||(k>SField.NX-12))    SField.P(1,k,j,i);
  
  
  
  //save the undersampled image
  SField.Write(Output);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++++++++++ Create a new void image ++++++++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void CreateVoidImage(int argc, char **argv){
  ScalarField SField;
  char Output[256];
  int dx,dy,dz;
  float rx,ry,rz;
  int i,j,k;
  
  //read parameters
  argc--; argv++;
  strcpy(Output,argv[1]);
  argc--; argv++;
  dx = atoi(argv[1]);
  argc--; argv++;
  dy = atoi(argv[1]);
  argc--; argv++;
  dz = atoi(argv[1]);
  argc--; argv++;
  rx = atof(argv[1]);
  argc--; argv++;
  ry = atof(argv[1]);
  argc--; argv++;
  rz = atof(argv[1]);
  argc--; argv++;
  
  
  //generate the new image and set its values to zero
  SField.CreateVoidField(dx,dy,dz);
  
  //define the image resolution
  SField.Image2World[0][0]=rx; SField.Image2World[0][1]=0;  SField.Image2World[0][2]=0;  SField.Image2World[0][3]=0;
  SField.Image2World[1][0]=0;  SField.Image2World[1][1]=ry; SField.Image2World[1][2]=0;  SField.Image2World[1][3]=0;
  SField.Image2World[2][0]=0;  SField.Image2World[2][1]=0;  SField.Image2World[2][2]=rz; SField.Image2World[2][3]=0;
  SField.Image2World[3][0]=0;  SField.Image2World[3][1]=0;  SField.Image2World[3][2]=0;  SField.Image2World[3][3]=1;
  
  //save the undersampled image
  SField.Write(Output);
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++ evaluate the accuracy of a displacement field with 3D points +++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Compare_LDMK_Points(int argc, char **argv){
  LDMK_Points Target_LDMK_Points;
  LDMK_Points Source_LDMK_Points;
  char TargetLM[256];
  char SourceLM[256];
  VectorField DFi;
  char File1[256];
  char File2[256];
  char File3[256];
  int i,j,Nb_LDMK_Points;
  float flX,flY,flZ;
  float flX2,flY2,flZ2;
  float flX3,flY3,flZ3;
  float srX,srY,srZ;
  float Distance;
  
  //read parameters
  argc--; argv++;
  strcpy(TargetLM,argv[1]);
  argc--; argv++;
  strcpy(File1,argv[1]);
  argc--; argv++;
  strcpy(File2,argv[1]);
  argc--; argv++;
  strcpy(File3,argv[1]);
  argc--; argv++;
  strcpy(SourceLM,argv[1]);
  argc--; argv++;
  
  
  //read the LDMK_Points
  Target_LDMK_Points.Read(TargetLM);
  Source_LDMK_Points.Read(SourceLM);
  
  Nb_LDMK_Points=Target_LDMK_Points.Get_LDMK_PointsNumber();
  if (Nb_LDMK_Points!=Source_LDMK_Points.Get_LDMK_PointsNumber()) cout << "The two files have a different number of points";
  
  //read the displacement field
  DFi.Read(File1,File2,File3);
  
  
  //compare the LDMK_Points
  for (i=0;i<Nb_LDMK_Points;i++){
    flX=Target_LDMK_Points.GetX(i);
    flY=Target_LDMK_Points.GetY(i);
    flZ=Target_LDMK_Points.GetZ(i);
    
    flX2=flX*DFi.World2Image[0][0]+flY*DFi.World2Image[0][1]+flZ*DFi.World2Image[0][2]+DFi.World2Image[0][3];
    flY2=flX*DFi.World2Image[1][0]+flY*DFi.World2Image[1][1]+flZ*DFi.World2Image[1][2]+DFi.World2Image[1][3];
    flZ2=flX*DFi.World2Image[2][0]+flY*DFi.World2Image[2][1]+flZ*DFi.World2Image[2][2]+DFi.World2Image[2][3];
    
    
    flX3=flX+DFi.G(0,flX2,flY2,flZ2);
    flY3=flY+DFi.G(1,flX2,flY2,flZ2);
    flZ3=flZ+DFi.G(2,flX2,flY2,flZ2);
    
    srX=Source_LDMK_Points.GetX(i);
    srY=Source_LDMK_Points.GetY(i);
    srZ=Source_LDMK_Points.GetZ(i);
    
    Distance=sqrt((flX3-srX)*(flX3-srX) + (flY3-srY)*(flY3-srY) + (flZ3-srZ)*(flZ3-srZ));
    
    cout << flX3 << " " << flY3 << " " << flZ3 << "\t<-> "  << srX << " " << srY << " " << srZ << "\t Distance = "  << Distance << endl;
  }
  
  cout << "If unexpected results, check that the landmarks coordinates are nifti coordinates (and not vtk coordinates as often)" << endl;
  
  
}  


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++  compose two displacement fields  +++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void ComposeDispFields(int argc, char **argv){
  VectorField VField1;
  VectorField VFieldComp;
  char Input1_X[256];
  char Input2_X[256];
  char Output_X[256];
  char Input1_Y[256];
  char Input2_Y[256];
  char Output_Y[256];
  char Input1_Z[256];
  char Input2_Z[256];
  char Output_Z[256];
  int x,y,z;
  float flX,flY,flZ;
  float tmpX,tmpY,tmpZ;
  float tmpX2,tmpY2,tmpZ2;
  float srcX,srcY,srcZ;
  float trgX,trgY,trgZ;
  
  
  //read parameters
  argc--; argv++;
  strcpy(Input1_X,argv[1]);
  argc--; argv++;
  strcpy(Input1_Y,argv[1]);
  argc--; argv++;
  strcpy(Input1_Z,argv[1]);
  argc--; argv++;
  strcpy(Input2_X,argv[1]);
  argc--; argv++;
  strcpy(Input2_Y,argv[1]);
  argc--; argv++;
  strcpy(Input2_Z,argv[1]);
  argc--; argv++;
  strcpy(Output_X,argv[1]);
  argc--; argv++;
  strcpy(Output_Y,argv[1]);
  argc--; argv++;
  strcpy(Output_Z,argv[1]);
  argc--; argv++;
  
  
  //read the image
  VField1.Read(Input1_X,Input1_Y,Input1_Z);
  VFieldComp.Read(Input2_X,Input2_Y,Input2_Z);
  
  //do the job
  for (z=0;z<VFieldComp.NZ;z++) for (y=0;y<VFieldComp.NY;y++) for (x=0;x<VFieldComp.NX;x++) {
    flX=static_cast<float>(x);
    flY=static_cast<float>(y);
    flZ=static_cast<float>(z);
    
    srcX=flX*VFieldComp.Image2World[0][0]+flY*VFieldComp.Image2World[0][1]+flZ*VFieldComp.Image2World[0][2]+VFieldComp.Image2World[0][3];
    srcY=flX*VFieldComp.Image2World[1][0]+flY*VFieldComp.Image2World[1][1]+flZ*VFieldComp.Image2World[1][2]+VFieldComp.Image2World[1][3];
    srcZ=flX*VFieldComp.Image2World[2][0]+flY*VFieldComp.Image2World[2][1]+flZ*VFieldComp.Image2World[2][2]+VFieldComp.Image2World[2][3];
    
    tmpX=srcX+VFieldComp.G(0,x,y,z);
    tmpY=srcY+VFieldComp.G(1,x,y,z);
    tmpZ=srcZ+VFieldComp.G(2,x,y,z);
    
    tmpX2=tmpX*VField1.World2Image[0][0]+tmpY*VField1.World2Image[0][1]+tmpZ*VField1.World2Image[0][2]+VField1.World2Image[0][3];
    tmpY2=tmpX*VField1.World2Image[1][0]+tmpY*VField1.World2Image[1][1]+tmpZ*VField1.World2Image[1][2]+VField1.World2Image[1][3];
    tmpZ2=tmpX*VField1.World2Image[2][0]+tmpY*VField1.World2Image[2][1]+tmpZ*VField1.World2Image[2][2]+VField1.World2Image[2][3];
    
    trgX=tmpX+VField1.G(0,tmpX2,tmpY2,tmpZ2);
    trgY=tmpY+VField1.G(1,tmpX2,tmpY2,tmpZ2);
    trgZ=tmpZ+VField1.G(2,tmpX2,tmpY2,tmpZ2);
    
    VFieldComp.P(trgX-srcX,0,x,y,z);
    VFieldComp.P(trgY-srcY,1,x,y,z);
    VFieldComp.P(trgZ-srcZ,2,x,y,z);
  }
  
  
  //save the undersampled image
  VFieldComp.Write(Output_X,Output_Y,Output_Z,Input2_X);
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++ transform an image with a displacement field ++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TransfoImag(int argc, char **argv){
  ScalarField InputI;
  ScalarField OutputI;
  VectorField VField;
  char InputImage[256];
  char DF_X[256];
  char DF_Y[256];
  char DF_Z[256];
  char OutputImage[256];
  int x,y,z;
  int NN;
  float flX,flY,flZ;
  float tmpX,tmpY,tmpZ;
  float tmpX2,tmpY2,tmpZ2;
  float srcX,srcY,srcZ;
  float trgX,trgY,trgZ;
  
  
  //read parameters
  argc--; argv++;
  strcpy(InputImage,argv[1]);
  argc--; argv++;
  strcpy(DF_X,argv[1]);
  argc--; argv++;
  strcpy(DF_Y,argv[1]);
  argc--; argv++;
  strcpy(DF_Z,argv[1]);
  argc--; argv++;
  strcpy(OutputImage,argv[1]);
  argc--; argv++;
  
  NN=0;
  if (argc>1){
    NN=1;
    cout << "Nearest neighbor interpolation" << endl;
  }
  

  //read the image
  VField.Read(DF_X,DF_Y,DF_Z);
  InputI.Read(InputImage);
  OutputI.Read(DF_X);
  
  
  //do the job
  for (z=0;z<OutputI.NZ;z++) for (y=0;y<OutputI.NY;y++) for (x=0;x<OutputI.NX;x++) {
    flX=static_cast<float>(x);
    flY=static_cast<float>(y);
    flZ=static_cast<float>(z);
    
    srcX=flX*OutputI.Image2World[0][0]+flY*OutputI.Image2World[0][1]+flZ*OutputI.Image2World[0][2]+OutputI.Image2World[0][3];
    srcY=flX*OutputI.Image2World[1][0]+flY*OutputI.Image2World[1][1]+flZ*OutputI.Image2World[1][2]+OutputI.Image2World[1][3];
    srcZ=flX*OutputI.Image2World[2][0]+flY*OutputI.Image2World[2][1]+flZ*OutputI.Image2World[2][2]+OutputI.Image2World[2][3];
    
    tmpX=srcX+VField.G(0,x,y,z);
    tmpY=srcY+VField.G(1,x,y,z);
    tmpZ=srcZ+VField.G(2,x,y,z);
    
    tmpX2=tmpX*InputI.World2Image[0][0]+tmpY*InputI.World2Image[0][1]+tmpZ*InputI.World2Image[0][2]+InputI.World2Image[0][3];
    tmpY2=tmpX*InputI.World2Image[1][0]+tmpY*InputI.World2Image[1][1]+tmpZ*InputI.World2Image[1][2]+InputI.World2Image[1][3];
    tmpZ2=tmpX*InputI.World2Image[2][0]+tmpY*InputI.World2Image[2][1]+tmpZ*InputI.World2Image[2][2]+InputI.World2Image[2][3];
    
    if (NN!=1)
      OutputI.P(InputI.G(tmpX2,tmpY2,tmpZ2),x,y,z);
    else{
      if ((static_cast<int>(tmpX2+0.5)>=0)&&(static_cast<int>(tmpX2+0.5)<InputI.NX)&&(static_cast<int>(tmpY2+0.5)>=0)&&(static_cast<int>(tmpY2+0.5)<InputI.NY)&&(static_cast<int>(tmpZ2+0.5)>=0)&&(static_cast<int>(tmpZ2+0.5)<InputI.NZ))
        OutputI.P(InputI.G(static_cast<int>(tmpX2+0.5),static_cast<int>(tmpY2+0.5),static_cast<int>(tmpZ2+0.5)),x,y,z);
      else
        OutputI.P(0,x,y,z);
    }
  }
  
  //save the image after de-allocating sime memory
  VField.SlashFieldSize(0);
  
  OutputI.Write(OutputImage,DF_X);
  
}


//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++      Extract a 2D slice from a 3D image      ++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void Extract2Dslice(int argc, char **argv){
  char In3DImFile[256];
  char Out2DpTImFile[256];
  float centerCoord[3];
  float normalVector[3];
  float Tvec1[3];
  float Tvec2[3];
  ScalarField ImageIn;
  ScalarField ImageOut;
  float tmpFl;
  float VoxReso;
  int i,j,k;
  float tranfoMat[4][4];
  float Radius;
  int NbGreyLevDirec;
  
  //1) read parameters
  argc--; argv++;
  strcpy(In3DImFile,argv[1]);
  argc--; argv++;
  centerCoord[0] = atof(argv[1]);
  argc--; argv++;
  centerCoord[1]= atof(argv[1]);
  argc--; argv++;
  centerCoord[2]= atof(argv[1]);
  argc--; argv++;
  normalVector[0] = atof(argv[1]);
  argc--; argv++;
  normalVector[1]= atof(argv[1]);
  argc--; argv++;
  normalVector[2]= atof(argv[1]);
  argc--; argv++;
  Radius= atof(argv[1]);
  argc--; argv++;
  strcpy(Out2DpTImFile,argv[1]);
  argc--; argv++;
  
  
  //2) load input image
  ImageIn.Read(In3DImFile);

  //3) compute the voxel resolution of the output image (same as the finest one in the input image)
  VoxReso=sqrt(ImageIn.Image2World[0][0]*ImageIn.Image2World[0][0]+ImageIn.Image2World[0][1]*ImageIn.Image2World[0][1]+ImageIn.Image2World[0][2]*ImageIn.Image2World[0][2]);
  
  tmpFl=sqrt(ImageIn.Image2World[1][0]*ImageIn.Image2World[1][0]+ImageIn.Image2World[1][1]*ImageIn.Image2World[1][1]+ImageIn.Image2World[1][2]*ImageIn.Image2World[1][2]);
  if (VoxReso>tmpFl)
    VoxReso=tmpFl;
  
  tmpFl=sqrt(ImageIn.Image2World[2][0]*ImageIn.Image2World[2][0]+ImageIn.Image2World[2][1]*ImageIn.Image2World[2][1]+ImageIn.Image2World[2][2]*ImageIn.Image2World[2][2]);
  if (VoxReso>tmpFl)
    VoxReso=tmpFl;
  
  //4) compute the vectors defining the tangent plane in the 2D output image and normalize it with the output image voxel resolution
  CptVecsTangentPlane(normalVector,Tvec1,Tvec2);
  VecNormalize(normalVector,VoxReso);
  VecNormalize(Tvec1,VoxReso);
  VecNormalize(Tvec2,VoxReso);
  
  //5) allocate the output image and compute its properties
  NbGreyLevDirec=static_cast<int>(Radius/VoxReso);
  
  ImageOut.CreateVoidField(2*NbGreyLevDirec+1,2*NbGreyLevDirec+1);
  
  ImageOut.Image2World[0][0]=Tvec1[0]; ImageOut.Image2World[0][1]=Tvec2[0]; ImageOut.Image2World[0][2]=normalVector[0]; 
  ImageOut.Image2World[1][0]=Tvec1[1]; ImageOut.Image2World[1][1]=Tvec2[1]; ImageOut.Image2World[1][2]=normalVector[1]; 
  ImageOut.Image2World[2][0]=Tvec1[2]; ImageOut.Image2World[2][1]=Tvec2[2]; ImageOut.Image2World[2][2]=normalVector[2]; 
  ImageOut.Image2World[3][0]=0;        ImageOut.Image2World[3][1]=0;        ImageOut.Image2World[3][2]=0;        
  
  ImageOut.Image2World[0][3]=centerCoord[0]-NbGreyLevDirec*(Tvec1[0]+Tvec2[0]);
  ImageOut.Image2World[1][3]=centerCoord[1]-NbGreyLevDirec*(Tvec1[1]+Tvec2[1]);
  ImageOut.Image2World[2][3]=centerCoord[2]-NbGreyLevDirec*(Tvec1[2]+Tvec2[2]);
  ImageOut.Image2World[3][3]=1;
  
  invert_4t4quaternion(ImageOut.Image2World,ImageOut.World2Image);
  
  //6) fill the 2D image
  
  mult_quat4t4mat_quat4t4mat(ImageIn.World2Image,ImageOut.Image2World,tranfoMat);
  
  for (j=0;j<ImageOut.NY;j++) for (i=0;i<ImageOut.NX;i++) 
      ImageOut.P(ImageIn.G(tranfoMat,i,j,0),i,j);
  
  //7) save the 2D slice
  ImageOut.Write(Out2DpTImFile);
}



//++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++      Set the image to world matrix     ++++++++
//++++++++++++++++++++++++++++++++++++++++++++++++++++++

void SetImage2World(int argc, char **argv){
  ScalarField InputI;
  char InputImage[256];
  char OutputImage[256];
  float Loc_Image2World[4][4];
  
  //read parameters
  argc--; argv++;
  strcpy(InputImage,argv[1]);
  argc--; argv++;
  Loc_Image2World[0][0]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[0][1]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[0][2]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[0][3]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[1][0]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[1][1]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[1][2]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[1][3]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[2][0]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[2][1]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[2][2]=atof(argv[1]);
  argc--; argv++;
  Loc_Image2World[2][3]=atof(argv[1]);
  argc--; argv++;
  strcpy(OutputImage,argv[1]);
  argc--; argv++;
  
  
  //read the image
  InputI.Read(InputImage);
  
  InputI.Image2World[0][0]=Loc_Image2World[0][0]; InputI.Image2World[0][1]=Loc_Image2World[0][1]; InputI.Image2World[0][2]=Loc_Image2World[0][2]; InputI.Image2World[0][3]=Loc_Image2World[0][3];
  InputI.Image2World[1][0]=Loc_Image2World[1][0]; InputI.Image2World[1][1]=Loc_Image2World[1][1]; InputI.Image2World[1][2]=Loc_Image2World[1][2]; InputI.Image2World[1][3]=Loc_Image2World[1][3];
  InputI.Image2World[2][0]=Loc_Image2World[2][0]; InputI.Image2World[2][1]=Loc_Image2World[2][1]; InputI.Image2World[2][2]=Loc_Image2World[2][2]; InputI.Image2World[2][3]=Loc_Image2World[2][3];
  InputI.Image2World[3][0]=0;                     InputI.Image2World[3][1]=0;                     InputI.Image2World[3][2]=0;                     InputI.Image2World[3][3]=1;
  
  
  invert_4t4quaternion(InputI.Image2World,InputI.World2Image);
  
  //save the image
  InputI.Write(OutputImage);
  
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++      Transform the 'image 2 world' matrix of [SourceImage] so that voxel [i][j][k] is the center of the world coordinates    +++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void DefineCenterWorldCoord(int argc, char **argv){
  ScalarField InputI;
  char InputImage[256];
  char OutputImage[256];
  int i,j,k;
  
  //read parameters
  argc--; argv++;
  strcpy(InputImage,argv[1]);
  argc--; argv++;
  i=atoi(argv[1]);
  argc--; argv++;
  j=atoi(argv[1]);
  argc--; argv++;
  k=atoi(argv[1]);
  argc--; argv++;
  strcpy(OutputImage,argv[1]);
  argc--; argv++;
  
  
  //read the image
  InputI.Read(InputImage);
  
  InputI.Image2World[0][3]=-i*InputI.Image2World[0][0]-j*InputI.Image2World[0][1]-k*InputI.Image2World[0][2];
  InputI.Image2World[1][3]=-i*InputI.Image2World[1][0]-j*InputI.Image2World[1][1]-k*InputI.Image2World[1][2];
  InputI.Image2World[2][3]=-i*InputI.Image2World[2][0]-j*InputI.Image2World[2][1]-k*InputI.Image2World[2][2];
  
  
  invert_4t4quaternion(InputI.Image2World,InputI.World2Image);
  
  //save the image
  InputI.Write(OutputImage);
  
}




//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++   affine transformation of an image   ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TransfoImagAffine(int argc, char **argv){
  ScalarField InputI;
  ScalarField OutputI;
  char InputImage[256];
  char TargetImage[256];
  char OutputImage[256];
  int x,y,z;
  float flX,flY,flZ;
  float srcX,srcY,srcZ;
  float tmpX,tmpY,tmpZ;
  float tmpX2,tmpY2,tmpZ2;
  float World_Target2Template[4][4];
  int NoEtrapo;
  
  //read parameters
  argc--; argv++;
  strcpy(InputImage,argv[1]);
  argc--; argv++;
  strcpy(TargetImage,argv[1]);
  argc--; argv++;
  World_Target2Template[0][0]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[0][1]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[0][2]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[0][3]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[1][0]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[1][1]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[1][2]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[1][3]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[2][0]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[2][1]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[2][2]=atof(argv[1]);
  argc--; argv++;
  World_Target2Template[2][3]=atof(argv[1]);
  argc--; argv++;
  strcpy(OutputImage,argv[1]);
  argc--; argv++;
  
  NoEtrapo=0;
  if (argc>1){
    NoEtrapo=1;
  }
  
 
  
  World_Target2Template[3][0]=0;
  World_Target2Template[3][1]=0;
  World_Target2Template[3][2]=0;
  World_Target2Template[3][3]=1;
  
  
  //read the images
  InputI.Read(InputImage);
  OutputI.Read(TargetImage);
  
  
  //do the job
  for (z=0;z<OutputI.NZ;z++) for (y=0;y<OutputI.NY;y++) for (x=0;x<OutputI.NX;x++) {
    flX=static_cast<float>(x);
    flY=static_cast<float>(y);
    flZ=static_cast<float>(z);
    
    srcX=flX*OutputI.Image2World[0][0]+flY*OutputI.Image2World[0][1]+flZ*OutputI.Image2World[0][2]+OutputI.Image2World[0][3];
    srcY=flX*OutputI.Image2World[1][0]+flY*OutputI.Image2World[1][1]+flZ*OutputI.Image2World[1][2]+OutputI.Image2World[1][3];
    srcZ=flX*OutputI.Image2World[2][0]+flY*OutputI.Image2World[2][1]+flZ*OutputI.Image2World[2][2]+OutputI.Image2World[2][3];
    
    tmpX=srcX*World_Target2Template[0][0]+srcY*World_Target2Template[0][1]+srcZ*World_Target2Template[0][2]+World_Target2Template[0][3];
    tmpY=srcX*World_Target2Template[1][0]+srcY*World_Target2Template[1][1]+srcZ*World_Target2Template[1][2]+World_Target2Template[1][3];
    tmpZ=srcX*World_Target2Template[2][0]+srcY*World_Target2Template[2][1]+srcZ*World_Target2Template[2][2]+World_Target2Template[2][3];
    
    tmpX2=tmpX*InputI.World2Image[0][0]+tmpY*InputI.World2Image[0][1]+tmpZ*InputI.World2Image[0][2]+InputI.World2Image[0][3];
    tmpY2=tmpX*InputI.World2Image[1][0]+tmpY*InputI.World2Image[1][1]+tmpZ*InputI.World2Image[1][2]+InputI.World2Image[1][3];
    tmpZ2=tmpX*InputI.World2Image[2][0]+tmpY*InputI.World2Image[2][1]+tmpZ*InputI.World2Image[2][2]+InputI.World2Image[2][3];
    
    if (NoEtrapo==1) OutputI.P(InputI.G_NoExtrapo(tmpX2,tmpY2,tmpZ2),x,y,z);
    else OutputI.P(InputI.G(tmpX2,tmpY2,tmpZ2),x,y,z);
    
  }
  
  OutputI.Write(OutputImage,TargetImage);
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++   affine transformation of an image   ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TransfoImagAffine2(int argc, char **argv){
  ScalarField InputI;
  ScalarField OutputI;
  char InputImage[256];
  char TargetImage[256];
  char OutputImage[256];
  char TransfoMatrixFile[256];
  int x,y,z;
  float flX,flY,flZ;
  float srcX,srcY,srcZ;
  float tmpX,tmpY,tmpZ;
  float tmpX2,tmpY2,tmpZ2;
  float World_Target2Template[4][4];
  int NoEtrapo;
  
  //read parameters
  argc--; argv++;
  strcpy(InputImage,argv[1]);
  argc--; argv++;
  strcpy(TargetImage,argv[1]);
  argc--; argv++;
  strcpy(TransfoMatrixFile,argv[1]);
  argc--; argv++;
  strcpy(OutputImage,argv[1]);
  argc--; argv++;
  
  NoEtrapo=0;
  if (argc>1){
    NoEtrapo=1;
  }
  

  
  //read the images
  InputI.Read(InputImage);
  OutputI.Read(TargetImage);
  
  //read the transformation matrix
  Read_quat4t4mat(TransfoMatrixFile,World_Target2Template);
  
  
  //do the job
  for (z=0;z<OutputI.NZ;z++) for (y=0;y<OutputI.NY;y++) for (x=0;x<OutputI.NX;x++) {
    flX=static_cast<float>(x);
    flY=static_cast<float>(y);
    flZ=static_cast<float>(z);
    
    srcX=flX*OutputI.Image2World[0][0]+flY*OutputI.Image2World[0][1]+flZ*OutputI.Image2World[0][2]+OutputI.Image2World[0][3];
    srcY=flX*OutputI.Image2World[1][0]+flY*OutputI.Image2World[1][1]+flZ*OutputI.Image2World[1][2]+OutputI.Image2World[1][3];
    srcZ=flX*OutputI.Image2World[2][0]+flY*OutputI.Image2World[2][1]+flZ*OutputI.Image2World[2][2]+OutputI.Image2World[2][3];
    
    tmpX=srcX*World_Target2Template[0][0]+srcY*World_Target2Template[0][1]+srcZ*World_Target2Template[0][2]+World_Target2Template[0][3];
    tmpY=srcX*World_Target2Template[1][0]+srcY*World_Target2Template[1][1]+srcZ*World_Target2Template[1][2]+World_Target2Template[1][3];
    tmpZ=srcX*World_Target2Template[2][0]+srcY*World_Target2Template[2][1]+srcZ*World_Target2Template[2][2]+World_Target2Template[2][3];
    
    /*if ((x==50)&&(y==50)&&(z==50)){
     cout << srcX << " " << srcY << " " << srcZ << endl;
     cout << tmpX << " " << tmpY << " " << tmpZ << endl;
     }*/
    
    tmpX2=tmpX*InputI.World2Image[0][0]+tmpY*InputI.World2Image[0][1]+tmpZ*InputI.World2Image[0][2]+InputI.World2Image[0][3];
    tmpY2=tmpX*InputI.World2Image[1][0]+tmpY*InputI.World2Image[1][1]+tmpZ*InputI.World2Image[1][2]+InputI.World2Image[1][3];
    tmpZ2=tmpX*InputI.World2Image[2][0]+tmpY*InputI.World2Image[2][1]+tmpZ*InputI.World2Image[2][2]+InputI.World2Image[2][3];
    
    if (NoEtrapo==1) OutputI.P(InputI.G_NoExtrapo(tmpX2,tmpY2,tmpZ2),x,y,z);
    else OutputI.P(InputI.G(tmpX2,tmpY2,tmpZ2),x,y,z);
  }
  
  OutputI.Write(OutputImage,TargetImage);
  
}



//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//++++++++   affine transformation of an image   ++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void TransfoImagAffine3(int argc, char **argv){
  ScalarField OutputI;
  char InputImage[256];
  char OutputImage[256];
  char TransfoMatrixFile[256];
  float TransfoMatrix[4][4];
  float InvTransfoMatrix[4][4];
  float NewImage2World[4][4];
  int i,j;
  
  //read parameters
  argc--; argv++;
  strcpy(InputImage,argv[1]);
  argc--; argv++;
  strcpy(TransfoMatrixFile,argv[1]);
  argc--; argv++;
  strcpy(OutputImage,argv[1]);
  argc--; argv++;
  
  //read the images
  OutputI.Read(InputImage);
  
  //read the transformation matrix
  Read_quat4t4mat(TransfoMatrixFile,TransfoMatrix);
  
  //multiply the matrices
  invert_4t4quaternion(TransfoMatrix,InvTransfoMatrix);
  mult_quat4t4mat_quat4t4mat(InvTransfoMatrix, OutputI.Image2World,  NewImage2World);
  
  for (i=0;i<4;i++) for (j=0;j<4;j++)
    OutputI.Image2World[i][j]=NewImage2World[i][j];
  
  invert_4t4quaternion(OutputI.Image2World,OutputI.World2Image);
    
  OutputI.Write(OutputImage);
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++  remove the affine contrib of a displacement field  +++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void RemAffineContribDF(int argc, char **argv){
  VectorField DFi;
  VectorField DFo;
  char fileDFi_X[256];
  char fileDFi_Y[256];
  char fileDFi_Z[256];
  char fileDFo_X[256];
  char fileDFo_Y[256];
  char fileDFo_Z[256];
  int x,y,z;
  float flX,flY,flZ;
  float srcX,srcY,srcZ;
  float tmpX,tmpY,tmpZ;
  float trgX,trgY,trgZ;
  float LocalRigReg[4][4];
  float LocalInvRigReg[4][4];
  
  //read parameters
  argc--; argv++;
  strcpy(fileDFi_X,argv[1]);
  argc--; argv++;
  strcpy(fileDFi_Y,argv[1]);
  argc--; argv++;
  strcpy(fileDFi_Z,argv[1]);
  argc--; argv++;
  LocalRigReg[0][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][3]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][3]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][3]= atof(argv[1]);
  argc--; argv++;
  strcpy(fileDFo_X,argv[1]);
  argc--; argv++;
  strcpy(fileDFo_Y,argv[1]);
  argc--; argv++;
  strcpy(fileDFo_Z,argv[1]);
  argc--; argv++;
  
  LocalRigReg[3][0]=0;
  LocalRigReg[3][1]=0;
  LocalRigReg[3][2]=0;
  LocalRigReg[3][3]=1;
  
  //read the input field and allocate the output field
  DFi.Read(fileDFi_X,fileDFi_Y,fileDFi_Z);
  DFo.Read(fileDFi_X,fileDFi_Y,fileDFi_Z);
  
  //invert the rigid transformation
  invert_4t4quaternion(LocalRigReg,LocalInvRigReg);
  
  //do the job
  for (z=0;z<DFi.NZ;z++) for (y=0;y<DFi.NY;y++) for (x=0;x<DFi.NX;x++) {
    flX=static_cast<float>(x);
    flY=static_cast<float>(y);
    flZ=static_cast<float>(z);
    
    srcX=flX*DFi.Image2World[0][0]+flY*DFi.Image2World[0][1]+flZ*DFi.Image2World[0][2]+DFi.Image2World[0][3];
    srcY=flX*DFi.Image2World[1][0]+flY*DFi.Image2World[1][1]+flZ*DFi.Image2World[1][2]+DFi.Image2World[1][3];
    srcZ=flX*DFi.Image2World[2][0]+flY*DFi.Image2World[2][1]+flZ*DFi.Image2World[2][2]+DFi.Image2World[2][3];
    
    tmpX=srcX+DFi.G(0,x,y,z);
    tmpY=srcY+DFi.G(1,x,y,z);
    tmpZ=srcZ+DFi.G(2,x,y,z);
    
    trgX=tmpX*LocalInvRigReg[0][0]+tmpY*LocalInvRigReg[0][1]+tmpZ*LocalInvRigReg[0][2]+LocalInvRigReg[0][3];
    trgY=tmpX*LocalInvRigReg[1][0]+tmpY*LocalInvRigReg[1][1]+tmpZ*LocalInvRigReg[1][2]+LocalInvRigReg[1][3];
    trgZ=tmpX*LocalInvRigReg[2][0]+tmpY*LocalInvRigReg[2][1]+tmpZ*LocalInvRigReg[2][2]+LocalInvRigReg[2][3];
    
    DFo.P(trgX-srcX,0,x,y,z);
    DFo.P(trgY-srcY,1,x,y,z);
    DFo.P(trgZ-srcZ,2,x,y,z);
  }
  
  DFo.Write(fileDFo_X,fileDFo_Y,fileDFo_Z,fileDFi_X);
  
  cout <<  endl;
  cout << "To obtain the amplitude of the deformations and the Jacobians, type:   " <<  "'uTIlzReg_Tools -DispFieldStudy " << fileDFo_X << " "  << fileDFo_Y << " "  << fileDFo_Z << "'" << endl;
  cout <<  endl;
  
}

//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++  generate a displacement field using an affine deformation  +++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void GenerAffineDefDF(int argc, char **argv){
  VectorField DF;
  char fileDF_X[256];
  char fileDF_Y[256];
  char fileDF_Z[256];
  char TargetIma[256];
  int x,y,z;
  float flX,flY,flZ;
  float srcX,srcY,srcZ;
  float trgX,trgY,trgZ;
  float LocalRigReg[4][4];
  
  //read parameters
  argc--; argv++;
  strcpy(TargetIma,argv[1]);
  argc--; argv++;
  LocalRigReg[0][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][3]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][3]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][3]= atof(argv[1]);
  argc--; argv++;
  strcpy(fileDF_X,argv[1]);
  argc--; argv++;
  strcpy(fileDF_Y,argv[1]);
  argc--; argv++;
  strcpy(fileDF_Z,argv[1]);
  argc--; argv++;
  
  LocalRigReg[3][0]=0;
  LocalRigReg[3][1]=0;
  LocalRigReg[3][2]=0;
  LocalRigReg[3][3]=1;
  
  //allocate the memory for the output field and get its image to world properties
  DF.Read(TargetIma,TargetIma,TargetIma);
  
  //do the job
  for (z=0;z<DF.NZ;z++) for (y=0;y<DF.NY;y++) for (x=0;x<DF.NX;x++) {
    flX=static_cast<float>(x);
    flY=static_cast<float>(y);
    flZ=static_cast<float>(z);
    
    trgX=flX*DF.Image2World[0][0]+flY*DF.Image2World[0][1]+flZ*DF.Image2World[0][2]+DF.Image2World[0][3];
    trgY=flX*DF.Image2World[1][0]+flY*DF.Image2World[1][1]+flZ*DF.Image2World[1][2]+DF.Image2World[1][3];
    trgZ=flX*DF.Image2World[2][0]+flY*DF.Image2World[2][1]+flZ*DF.Image2World[2][2]+DF.Image2World[2][3];
    
    srcX=trgX*LocalRigReg[0][0]+trgY*LocalRigReg[0][1]+trgZ*LocalRigReg[0][2]+LocalRigReg[0][3];
    srcY=trgX*LocalRigReg[1][0]+trgY*LocalRigReg[1][1]+trgZ*LocalRigReg[1][2]+LocalRigReg[1][3];
    srcZ=trgX*LocalRigReg[2][0]+trgY*LocalRigReg[2][1]+trgZ*LocalRigReg[2][2]+LocalRigReg[2][3];
    
    DF.P(srcX-trgX,0,x,y,z);
    DF.P(srcY-trgY,1,x,y,z);
    DF.P(srcZ-trgZ,2,x,y,z);
  }
  
  DF.Write(fileDF_X,fileDF_Y,fileDF_Z,TargetIma);
  
}


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++                  Invert the affine deformation              +++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void InvertAffineDef(int argc, char **argv){
  float LocalRigReg[4][4];
  float LocalInvRigReg[4][4];
  int i,j;
  
  //read parameters
  argc--; argv++;
  LocalRigReg[0][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[0][3]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[1][3]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][0]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][1]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][2]= atof(argv[1]);
  argc--; argv++;
  LocalRigReg[2][3]= atof(argv[1]);
  argc--; argv++;
  
  LocalRigReg[3][0]=0;
  LocalRigReg[3][1]=0;
  LocalRigReg[3][2]=0;
  LocalRigReg[3][3]=1;
  
  
  //do the job
  invert_4t4quaternion(LocalRigReg,LocalInvRigReg);
  
  //ouputs
  
  cout << endl;
  cout << "Input matrix:" << endl;
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      cout << LocalRigReg[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
  
  cout << "Inverse of the input matrix:" << endl;
  for (i=0;i<4;i++){
    for (j=0;j<4;j++){
      cout << LocalInvRigReg[i][j] << " ";
    }
    cout << endl;
  }
  cout << endl;
  
  cout << "Inverse of the input matrix  ([r_xx r_xy r_xz t_x  r_yx ... t_z]):" << endl;
  for (i=0;i<3;i++) for (j=0;j<4;j++)
    cout << LocalInvRigReg[i][j] << " ";
  cout << endl;
  cout << endl;
  
  
}




/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 OTHER Stuff
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


void otherStuff(int argc, char **argv){
  
  
  //estimate the radii on a curve
  ScalarField ImageIn;
  LDMK_Curves PointsIn;
  char ImFile[256];
  char PtsFile[256];
  char OutPtsFile[256];
  int i,j,x,y,z;
  float TmpFl,TmpFl2;
  float VoxStep;
  float thresh;
  int RadiusLoc;
  int tstOK;
  int PtsIn,PtsAll;
  
  
  //set some parameters
  thresh=12500;
  
  //read parameters
  argc--; argv++;
  strcpy(ImFile,argv[1]);
  argc--; argv++;
  strcpy(PtsFile,argv[1]);
  argc--; argv++;
  strcpy(OutPtsFile,argv[1]);
  argc--; argv++;
  
  //read image and points of interest
  ImageIn.Read(ImFile);
  PointsIn.Read(PtsFile);
  
  cout << PointsIn.GetSegNumber() << endl;
  
  //we suppose an orthogonal and isotropic grid
  VoxStep=fabs(ImageIn.Image2World[0][0]);
  cout << "Voxel steps: " << VoxStep << endl;
  
  
  //read stuffs
  for (i=0;i<PointsIn.GetSegNumber();i++){
    for (j=0;j<PointsIn.GetElNumber(i);j++){
      
      RadiusLoc=3;
      tstOK=1;
      
      while (tstOK==1){
        //test current radius
        tstOK=0;
        PtsAll=0;
        PtsIn=0;
        for (z=-RadiusLoc;z<=RadiusLoc;z++) for (y=-RadiusLoc;y<=RadiusLoc;y++) for (x=-RadiusLoc;x<=RadiusLoc;x++){
          TmpFl=sqrt(static_cast<float>(x*x)+static_cast<float>(y*y)+static_cast<float>(z*z));
          if (TmpFl<static_cast<float>(RadiusLoc)){
            PtsAll++;
            TmpFl2=ImageIn.G(ImageIn.World2Image,PointsIn.GetX(i,j)+x*VoxStep,PointsIn.GetY(i,j)+y*VoxStep,PointsIn.GetZ(i,j)+z*VoxStep);
            if (TmpFl2<thresh) PtsIn++;
          }
        }
        
        //try a larger radius if still OK
        if (static_cast<float>(PtsIn)/static_cast<float>(PtsAll)>0.99){
          tstOK=1;
          RadiusLoc++;
        }
      }
      PointsIn.PutD(static_cast<float>(RadiusLoc)*VoxStep,i,j);
      
      cout <<  "Segment " << i <<  " / Element " << j << ": Radius = " << PointsIn.GetD(i,j) << endl;
    }
  }
  
  PointsIn.Write(OutPtsFile);
  
  
  
  /*
   ScalarField SFieldToPut;
   ScalarField SFieldOut;
   float minGL,maxGL;
   char Input1[256];
   char Input2[256];
   char Output[256];
   int x,y,z;
   float minX,minY,minZ;
   int minXi,minYi,minZi;
   int DecX,DecY,DecZ;
   
   //read parameters
   argc--; argv++;
   strcpy(Input1,argv[1]);
   argc--; argv++;
   strcpy(Input2,argv[1]);
   argc--; argv++;
   strcpy(Output,argv[1]);
   
   SFieldOut.Read(Input1);
   for (z=0;z<2;z++) for (y=0;y<SFieldOut.NY;y++) for (x=0;x<SFieldOut.NX;x++)
   SFieldOut.P(0,x,y,z);
   SFieldOut.Write(Input2,Input1);
   */
  
  /*argc--; argv++;
   minX = atof(argv[1]);
   argc--; argv++;
   minY = atof(argv[1]);
   argc--; argv++;
   minZ = atof(argv[1]);
   */
  /*
   SFieldOut.Read(Input1);
   SFieldToPut.Read(Input2);
   
   for (z=0;z<SFieldOut.NZ;z++) for (y=0;y<SFieldOut.NY;y++) for (x=0;x<SFieldOut.NX;x++)
   if ((z<126)&&(SFieldToPut.G(x,y,z)==1)&&(SFieldOut.G(x,y,z)<-100))
   SFieldOut.P(-22,x,y,z);
   
   SFieldOut.Write(Output);
   */
  
  /*
   SFieldOut.Read(Input1);
   for (z=0;z<SFieldOut.NZ;z++) for (y=0;y<SFieldOut.NY;y++) for (x=0;x<SFieldOut.NX;x++)
   if ((x==20)||(x==491)||(y==485)||(y==75)||(y==55)||(y==35)||(y==15)||(z==435)) 
   SFieldOut.P(2,x,y,z);
   
   SFieldOut.Write(Input2);
   */
  
  
  //read the image
  
  /* 
   SFieldOut.Read(Input1);
   SFieldToPut.Read(Input2);
   
   cout << SFieldOut.Image2World[0][0] << " " << SFieldOut.Image2World[1][1] << " " << SFieldOut.Image2World[2][2]<< endl;
   cout << SFieldToPut.Image2World[0][0] << " " << SFieldToPut.Image2World[1][1] << " " << SFieldToPut.Image2World[2][2]<< endl;
   
   
   DecX=static_cast<int>((SFieldOut.NX-SFieldToPut.NX)/2.); if (DecX<0) DecX=0;
   DecY=static_cast<int>((SFieldOut.NY-SFieldToPut.NY)/2.); if (DecY<0) DecY=0;
   DecZ=static_cast<int>((SFieldOut.NZ-SFieldToPut.NZ)/2.); if (DecZ<0) DecZ=0;
   
   
   
   //do the job
   for (z=0;z<SFieldOut.NZ;z++) for (y=0;y<SFieldOut.NY;y++) for (x=0;x<SFieldOut.NX;x++) SFieldOut.P(0,x,y,z);
   
   for (z=0;z<SFieldToPut.NZ;z++) for (y=0;y<SFieldToPut.NY;y++) for (x=0;x<SFieldToPut.NX;x++)
   if ((x+DecX<SFieldOut.NX)&&(y+DecY<SFieldOut.NY)&&(z+DecZ<SFieldOut.NZ))
   SFieldOut.P(SFieldToPut.G(x,y,z),x+DecX,y+DecY,z+DecZ);
   
   //save the undersampled image
   SFieldOut.Write(Output);
   */
  
}


/*++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
 INPUTS MANAGEMENT
 ++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++*/


//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
//+++++++++++ linearly resample the grey levels +++++++++++++
//+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

void AllOptions(int argc, char **argv){
  cerr << "Usage: uTIlzReg_Tools <options>\n";
  cerr << "\n";
  cerr << "Where <options> are the following:\n";
  cerr << "\n";
  cerr << "(1) PARALLEL TRANSPORT STUFFS:\n";
  cerr << "\n";
  cerr << "-transport [Time subdivision i][Image T_i] [Velo Field X][Velo Field Y][Velo Field Z] [Time subdivision j][Image T_j]\n";
  cerr << "    -> Transport the image [Image T_i] from the time subdivision i to j using the velocity field. The result is\n";
  cerr << "    saved in [Image T_j].\n";
  cerr << "\n";
  cerr << "-MultiTransport [VFX][VFY][VFZ] [Src TimeSub] [Nb Ima] [Src Ima 1][Trg TimeSub 1][Trg Ima 1] ...\n";
  cerr << "    -> Transport [Nb Ima] images [Src Ima i] from the time subdivision [Src TimeSub] to the time subdivision [Trg TimeSub i]\n";
  cerr << "    through the velocity field [VFX][VFY][VFZ]. The transported images are saved in [Trg Ima i].\n";
  cerr << "\n";
  cerr << "\n";
  cerr << "(2) ARITHMETIC ON IMAGES:\n";
  cerr << "\n";
  cerr << "-weightedSum [Nb Images] [Weight1][Image1] [Weight2][Image2]... [Output Image]\n";
  cerr << "    -> Compute the weighted sum of one or several images\n";
  cerr << "\n";
  cerr << "-weightedMean [Nb Images] [Weight1][Image1] [Weight2][Image2]... [Output Image]\n";
  cerr << "    -> Compute the weighted mean of several images\n";
  cerr << "\n";
  cerr << "-weightedStd [Nb Images] [Weight1][Image1] [Weight2][Image2]... [Output Image]\n";
  cerr << "    -> Compute the weighted standard deviation between several images\n";
  cerr << "\n";
  cerr << "-NormGrad [ImageIn][ImageOut]\n";
  cerr << "    -> Compute the norm of the gradients of [ImageIn] to detect the boundaries. Result is stored in [ImageOut].\n";
  cerr << "\n";
  cerr << "\n";
  cerr << "(3) IMAGE TREATMENT:\n";
  cerr << "\n";
  cerr << "-AnisoDiff [ImageIn][ImageOut] [ax][ay][az] [IterationNb][dTau]\n";
  cerr << "    -> Perform anisotropic diffusion of [ImageIn]. Spatial diffusion parameters are [ax][ay][az] and temporal (iteration\n";
  cerr << "    more precisely) parameters are [IterationNb][dTau]. Result is stored in [ImageOut].\n";
  cerr << "\n";
  cerr << "-Diffusion [ImageIn][ImageOut] [IterationNb][dTau] <[Mask][MaskId]>\n";
  cerr << "    -> Perform isotropic diffusion of [ImageIn]. Spatial diffusion parameter is 1 and temporal (iteration\n";
  cerr << "    more precisely) parameters are [IterationNb][dTau]. Result is stored in [ImageOut]. If the options <[Mask][MaskId]>\n";
  cerr << "    are defined, the diffusion is only performed where the intensities of <Mask> are equal to [MaskId].\n";
  cerr << "\n";
  cerr << "-SmoothImage [ImageIn][ImageOut] [sigma]\n";
  cerr << "    -> Smooth an image using Gaussian kernel convolution. Input image is [ImageIn] and output image is [ImageOut].\n";
  cerr << "    The standard deviation of the Gaussian kernel is [sigma] millimeters\n";
  cerr << "\n";
  cerr << "-Erosion [ImageIn][ImageOut] [IterationNb] <[Mask][MaskId]>\n";
  cerr << "    -> Erosion of [ImageIn] during [IterationNb] iterations and with a structuring element of 3*3*3 voxels. Result is \n";
  cerr << "    stored in [ImageOut]. If the options <[Mask][MaskId]> are defined, the erosion is only performed where the intensities\n";
  cerr << "    of <Mask> are equal to [MaskId].\n";
  cerr << "\n";
  cerr << "-Dilation [ImageIn][ImageOut] [IterationNb] <[Mask][MaskId]>\n";
  cerr << "    -> Dilation of [ImageIn] during [IterationNb] iterations and with a structuring element of 3*3*3 voxels. Result is \n";
  cerr << "    stored in [ImageOut]. If the options <[Mask][MaskId]> are defined, the dilation is only performed where the intensities\n";
  cerr << "    of <Mask> are equal to [MaskId].\n";
  cerr << "\n";
  cerr << "-AddNoise [ImageIn][NoiseLevel][ImageOut] <[Xmin][Xmax][Ymin][Ymax][Zmin][Zmax]>\n";
  cerr << "    -> Add white noise on [ImageIn], where [NoiseLevel] is the standard deviation of the noise\n";
  cerr << "\n";
  cerr << "-PutGLInROI [ImageIn][ImageOut] [Rx1][Rx2][Ry1][Ry2][Rz1][Rz2] [GreyLev]\n";
  cerr << "    -> Put the grey level [GreyLev] in the ROI defined by [Rx1][Rx2][Ry1][Ry2][Rz1][Rz2] (in voxels)\n";
  cerr << "\n";
  cerr << "-PutGLInROI2 [ImageIn][ImageOut] [GreyLevIn][GreyLevOut]\n";
  cerr << "    -> In [ImageIn] put the intensity [GreyLevOut] at the voxels where there is [GreyLevIn]. Result is saved in [ImageOut]\n";
  cerr << "\n";
  cerr << "-MaxGreyLev [InputImage1][InputImage2] [OutputImage]\n";
  cerr << "    -> At each voxel [OutputImage] has the largest intensity of [InputImage1] and [InputImage2] (which have the same size)\n";
  cerr << "\n";
  cerr << "-MaskImage [ImageIn][ImageOut] [Mask][MaskId]\n";
  cerr << "    -> Apply the mask [Mask] to image [ImageIn]. The values at [MaskId] in [Mask] are set to zero in [ImageOut]\n";
  cerr << "\n";
  cerr << "-RemoveSmallReg [ImageIn][ImageOut] [IdVoxShape][IdVoxBackG] [MaxSize]\n";
  cerr << "    -> Put to [IdVoxBackG] the connected sets of voxels of intensity [IdVoxShape] having a size of less than [MaxSize] voxels.\n";
  cerr << "\n";
  cerr << "-KeepLargestReg [ImageIn][ImageOut] [IdVoxShape][IdVoxBackG]\n";
  cerr << "    -> Put to [IdVoxBackG] all the connected sets of voxels of intensity [IdVoxShape] expect the largest one.\n";
  cerr << "\n";
  cerr << "-GreyLevThresh [ImageIn][ImageOut] [NewMinGL][NewMaxGL]\n";
  cerr << "    -> Set to [NewMinGL] all the grey levels under [NewMinGL] and to [NewMaxGL] all the grey levels above [NewMaxGL].\n";
  cerr << "    -> Remark: if [NewMinGL]==[NewMaxGL] a traditional threshold around [NewMinGL] is done instead.\n";
  cerr << "\n";
  cerr << "-GreyLevResample [ImageIn][ImageOut] [NewMinGL][NewMaxGL]\n";
  cerr << "    -> Linearly Resample the grey levels between [NewMinGL] and [NewMaxGL].\n";
  cerr << "\n";
  cerr << "-CreateVoidImage [ImName] [DX][DY][DZ] [RX][RY][RZ]\n";
  cerr << "    -> Create a void image of dimension [DX][DY][DZ] voxels having a resolution of [RX][RY][RZ]mm.\n";
  cerr << "\n";
  cerr << "-Extract2Dslice [ImageIn] [Cx][Cy][Cz] [Nx][Ny][Nz] [Rad] [ImageOut]\n";
  cerr << "    -> Cut [ImageIn] in the 2D ROI centered in [Cx][Cy][Cz] and normal to the vector [Nx][Ny][Nz] (all in world coordinates).\n";
  cerr << "       The output ROI [ImageOut] has the resolution of [ImageIn]. If [ImageIn] is anisotropic the finest resolution is considered.\n";
  cerr << "       Warning: The ROI may be not be openable under all viewers. It works using rview (irtk).\n";
  cerr << "\n";
  cerr << "-ROI_Cut [ImageIn][ImageOut] [MinX][MaxX] [MinY][MaxY] [MinZ][MaxZ] [MinT][MaxT]\n";
  cerr << "    -> Cut [ImageIn] in the region of interest [MinX][MaxX] [MinY][MaxY] [MinZ][MaxZ] for the space and [MinT][MaxT]\n";
  cerr << "       for the time (use MinT=0 and MaxT=1 for 3D images).\n";
  cerr << "\n";
  cerr << "\n";
  cerr << "(4) TREATMENT OF DEFORMATIONS AND IMAGE SPACES:\n";
  cerr << "\n";
  cerr << "-LargeImageUndersample [ImageIn][ImageOut] [factor]\n";
  cerr << "    -> Undersample the resolution of [ImageIn] with an integer factor [factor]. The result is saved in [ImageOut].\n";
  cerr << "       Remark: problem may occur on 2D images. For 2D images, use -UnderSample.\n";
  cerr << "\n";
  cerr << "-UnderSample [ImageIn][ImageOut] [factor] <[-nn]>\n";
  cerr << "    -> Undersample the resolution of [ImageIn] with a float factor [factor]. The result is saved in [ImageOut].\n";
  cerr << "       Nearest neighbor interpolation is used if the option -nn is entered (tri-linear interpolation by default).\n";
  cerr << "\n";
  cerr << "-GreyLevAlign [ImToAlign][ImRef] [AlignedIm]\n";
  cerr << "    -> Aligment of the Grey-levels of [ImToAlign] on those of [ImRef] using histogram matching. Transformed image is saved in [AlignedIm].\n";
  cerr << "       An optimal transport strategy based on the Wasserstein 1 distance between the histograms is used here.\n";
  cerr << "\n";
  cerr << "-MakeSequence [N] N*[Frame n] [OutputSequence]\n";
  cerr << "    -> Gather N 3D images/frames [Frame n] in a single 3D+time image [OutputSequence].\n";
  cerr << "\n";
  cerr << "-RegCenterImage [MovingIma][FixedIma][TransfoIma]\n";
  cerr << "    -> Transform the 'image 2 world' properties of [ImageIn] to be centered on [ImageOut]. Result is in [TransfoIma].\n";
  cerr << "\n";
  cerr << "-SetImage2World [SourceImage] [r_xx][r_xy][r_xz][t_x][r_yx][...][t_z] [TransfoImage]\n";
  cerr << "    -> Define manually the 'image 2 world' matrix. \n";
  cerr << "\n";
  cerr << "-DefineCenterWorldCoord [SourceImage] [i][j][k] [TransfoImage]\n";
  cerr << "    -> Transform the 'image 2 world' matrix of [SourceImage] so that voxel [i][j][k] is the center of the world coordinates.\n";
  cerr << "\n";
  cerr << "-TransfoImag [SourceImage] [DF_X][DF_Y][DF_Z] [TransfoImage] <-NN>\n";
  cerr << "    -> Transform the image [SourceImage] using the displacement field [DF_X][DF_Y][DF_Z] (in mm). \n";
  cerr << "    -> The transformed image [TransfoImage] is in the D.F. coord. system.\n";
  cerr << "    -> If the option -NN is entered, nearest neighbor is used instead of trilinear interpolation.\n";
  cerr << "\n";
  cerr << "-TransfoImagAffine [SourceImage][TargetImage] [r_xx][r_xy][r_xz][t_x][r_yx][...][t_z] [TransfoImage] <-NoExtrapo>\n";
  cerr << "    -> Transform the image [SourceImage] using the affine transformation [r_xx][...][t_z] (from trg to src). \n";
  cerr << "    -> The output transformed image [TransfoImage] is in the coord. system of [TargetImage].\n";
  cerr << "    -> No extrapolation of the grey levels is performed if -NoExtrapo is set.\n";
  cerr << "\n";
  cerr << "-TransfoImagAffine2 [SourceImage][TargetImage] [RotTransMat] [TransfoImage] <-NoExtrapo>\n";
  cerr << "    -> Same as TransfoImagAffine except that the 4*4 transformation matrix is in the file [RotTransMat].\n";
  cerr << "\n";
  cerr << "-TransfoImagAffine3 [SourceImage] [RotTransMat] [TransfoImage]\n";
  cerr << "    -> Same as TransfoImagAffine2 but only the 'image to world' matrix of [SourceImage] is modified.\n";
  cerr << "\n";
  cerr << "-ComposeDispFields [DF1_X][DF1_Y][DF1_Z] [DF2_X][DF2_Y][DF2_Z] [DFC_X][DC1_Y][DC1_Z]\n";
  cerr << "    -> Compute DispField1 o DispField2 (the displacement fields are in mm / DispFieldComp <- coord. sys. of DispField2)\n";
  cerr << "\n";
  cerr << "-RemAffineContribDF [DFi_X][DFi_Y][DFi_Z]  [r_xx r_xy r_xz t_x  r_yx ... t_z] [DFo_X][DFo_Y][DFo_Z]\n";
  cerr << "    -> Suppose [DFi] = (Affine def.) o (Non-affine def.) -> remove the affine contrib [r_xx...] and save the result in [DFo]\n";
  cerr << "\n";
  cerr << "-GenerAffineDefDF [TargetIma] [r_xx r_xy r_xz t_x  r_yx ... t_z] [DF_X][DF_Y][DF_Z]\n";
  cerr << "    -> Generate a disp. field [DF_X][DF_Y][DF_Z] in the coord. sys. of [TargetIma] using the affine def. [r_xx r_xy ...].\n";
  cerr << "\n";
  cerr << "-InvertAffineDef [r_xx r_xy r_xz t_x  r_yx ... t_z]\n";
  cerr << "    -> Invert the affine deformation [r_xx r_xy r_xz t_x  r_yx ... t_z].\n";
  cerr << "\n";
  cerr << "\n";
  cerr << "(5) IMAGE AND DEFORMATIONS ANALYSIS:\n";
  cerr << "\n";
  cerr << "-Info [Image]\n";
  cerr << "    -> Basic information about an image\n";
  cerr << "\n";
  cerr << "-SSD [Image1] [Image2] <[Margin]>\n";
  cerr << "    -> Measure the sum of squared differences and the mutual information between two images\n";
  cerr << "\n";
  cerr << "-Overlap [Image1] [Image2] [Thresh] <[Margin]>\n";
  cerr << "    -> Compute de overlap (Dice metric) between two images. In each image, the shape is the grey levels strictly above [Thresh].\n";
  cerr << "\n";
  cerr << "-Overlap2 [Image1] [Image2]  <[FileWithOverlap]>\n";
  cerr << "    -> Compute de overlap (Dice metric) between two images. Different regions are represented by strictly positive integer values.\n";
  cerr << "    -> If the file name [FileWithOverlap] is defined, the value of the average overlap is written in this file.\n";
  cerr << "\n";
  cerr << "-Surf_Av_Std_values [RefIma] [Thresh] [Values]\n";
  cerr << "    -> Mesure the average and std of the intensities in [Values] at the surface of [RefIma] at the threshold [Thresh]\n";
  cerr << "\n";
  cerr << "-ROI_Av_Std_values [Image] [Rx1][Rx2][Ry1][Ry2][Rz1][Rz2]\n";
  cerr << "    -> Mesure the average and std of the intensities in the ROI of [Image] defined by [Rx1][Rx2][Ry1][Ry2][Rz1][Rz2]\n";
  cerr << "\n";
  cerr << "-DispFieldStudy [DefX][DefY][DefZ] <[Mask][MaskID]>\n";
  cerr << "    -> Quantify the properties of the desiplacement field [DefX][DefY][DefZ]. If the options <[Mask][MaskId]> are defined,\n";
  cerr << "     the properties are only studied where the intensities of <Mask> are equal to [MaskId].\n";
  cerr << "\n";
  cerr << "-CompareLandmarkPoints [Target lm file] [DefX][DefY][DefZ] [Source lm file]\n";
  cerr << "    -> Measure how the displacement field [DefX][DefY][DefZ] matches the target points [Target lm file] to the source \n";
  cerr << "     ones of [Source lm file]\n";
  cerr << "\n";
}



void usage(){
  cerr << "Usage: uTIlzReg_Tools <options>\n";
  cerr << "\n";
  cerr << "Where <options> are the following:\n";
  cerr << "\n";
  cerr << "-SSD [Image1] [Image2]  <[Margin]>\n";
  cerr << "    -> Measure the sum of squared differences and the mutual information between two images\n";
  cerr << "\n";
  cerr << "-weightedSum [Nb Images] [Weight1][Image1] [Weight2][Image2]... [Output Image]\n";
  cerr << "    -> Compute the weighted sum of one or several images\n";
  cerr << "\n";
  cerr << "-SmoothImage [ImageIn][ImageOut] [sigma]\n";
  cerr << "    -> Smooth an image using Gaussian kernel convolution. Input image is [ImageIn] and output image is [ImageOut].\n";
  cerr << "    The standard deviation of the Gaussian kernel is [sigma] millimeters\n";
  cerr << "\n";
  cerr << "-UnderSample [ImageIn][ImageOut] [factor] <[-nn]>\n";
  cerr << "    -> Undersample the resolution of [ImageIn] with a factor [factor]. The result is saved in [ImageOut].\n";
  cerr << "       Nearest neighbor interpolation is used if the option -nn is entered (tri-linear interpolation by default).\n";
  cerr << "\n";
  cerr << "-GreyLevThresh [ImageIn][ImageOut] [NewMinGL][NewMaxGL]\n";
  cerr << "    -> Set to [NewMinGL] all the grey levels under [NewMinGL] and to [NewMaxGL] all the grey levels above [NewMaxGL].\n";
  cerr << "    -> Remark: if [NewMinGL]==[NewMaxGL] a traditional threshold around [NewMinGL] is done instead.\n";
  cerr << "\n";
  cerr << "-GreyLevResample [ImageIn][ImageOut] [NewMinGL][NewMaxGL]\n";
  cerr << "    -> Linearly Resample the grey levels between [NewMinGL] and [NewMaxGL].\n";
  cerr << "\n";
  cerr << "-GreyLevAlign [ImToAlign][ImRef] [AlignedIm]\n";
  cerr << "    -> Aligment of the Grey-levels of [ImToAlign] on those of [ImRef] using histogram matching. Transformed image is saved in [AlignedIm].\n";
  cerr << "\n";
  cerr << "-AllOptions     (or -A)\n";
  cerr << "    -> Show all options.\n";
  cerr << "\n";
  
  exit(1);
}

int main(int argc, char **argv){
  bool done;
  int int1,int2;
  char File1[256];
  char File2[256];
  char File3[256];
  char File4[256];
  char File5[256];
  
  
  done=false;
  
  // Check command line
  if (argc < 2) {
    usage();
  }
  else {
    //+++++++++++++++++++ image transport +++++++++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-transport") == 0) {
      argc--; argv++;
      int1 = atoi(argv[1]); //first time step
      argc--; argv++;
      strcpy(File1,argv[1]); //image to transport
      argc--; argv++;
      strcpy(File2,argv[1]); //velocity field X
      argc--; argv++;
      strcpy(File3,argv[1]); //velocity field Y
      argc--; argv++;
      strcpy(File4,argv[1]); //velocity field Z
      argc--; argv++;
      int2 = atoi(argv[1]); //final time step
      argc--; argv++;
      strcpy(File5,argv[1]); //transported image
      argc--; argv++;
      
      transportImage(File1,File5,int1,int2,File2,File3,File4);
      done = true;
    }
    //+++++++++++++   transport several images +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-MultiTransport") == 0) {
      MultiTransport(argc,argv);
      done = true;
    }
    //+++++++++++++   image weighted sum  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-weightedSum") == 0) {
      weightedSum(argc,argv);
      done = true;
    }
    //+++++++++++++   images weighted mean  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-weightedMean") == 0) {
      weightedMean(argc,argv);
      done = true;
    }
    //+++++++++++++   images weighted standard deviation  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-weightedStd") == 0) {
      weightedStd(argc,argv);
      done = true;
    }
    //+++++++++++++   basic information about an image  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Info") == 0) {
      Info(argc,argv);
      done = true;
    }
    //+++++++++++++   sum of squared differences between two images  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-SSD") == 0) {
      MeasureSSD(argc,argv);
      done = true;
    }
    //+++++++++++++   Dice metric between two images  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Overlap") == 0) {
      MeasureOverlap(argc,argv);
      done = true;
    }
    //+++++++++++++   Dice metric between two images  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Overlap2") == 0) {
      MeasureOverlap2(argc,argv);
      done = true;
    }
    //+++++++++++++   quantify deformations  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-DispFieldStudy") == 0) {
      DispFieldStudy(argc,argv);
      done = true;
    }
    
    //+++++++++++++   anisotropic diffusion  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-AnisoDiff") == 0) {
      AnisoDiff(argc,argv);
      done = true;
    }
    //++++++++++++++   isotropic diffusion  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Diffusion") == 0) {
      Diffusion(argc,argv);
      done = true;
    }
    //++++++++++++++   isotropic smoothing  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-SmoothImage") == 0) {
      SmoothImage(argc,argv);
      done = true;
    }
    //++++++++++++++   erosion  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Erosion") == 0) {
      Erosion(argc,argv);
      done = true;
    }
    //++++++++++++++   dilation  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Dilation") == 0) {
      Dilation(argc,argv);
      done = true;
    }
    //++++++++++++++   add noise  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-AddNoise") == 0) {
      AddNoise(argc,argv);
      done = true;
    }
    //++++++++++++++   NormGrad  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-NormGrad") == 0) {
      NormGrad(argc,argv);
      done = true;
    }
    //++++++++++++++   mask image  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-MaskImage") == 0) {
      MaskImage(argc,argv);
      done = true;
    }
    //++++++++++++++   max grey level  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-MaxGreyLev") == 0) {
      MaxGreyLev(argc,argv);
      done = true;
    }
    //++++++++++++++   mask image  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-PutGLInROI") == 0) {
      PutGLInROI(argc,argv);
      done = true;
    }
    //++++++++++++++   mask image  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-PutGLInROI2") == 0) {
      PutGLInROI2(argc,argv);
      done = true;
    }
    //++++++++++++++   mask image  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-RemoveSmallReg") == 0) {
      RemoveSmallReg(argc,argv);
      done = true;
    }
    //++++++++++++++   mask image  ++++++++++++++++
    if (done == false) if (strcmp(argv[1], "-KeepLargestReg") == 0) {
      KeepLargestReg(argc,argv);
      done = true;
    }
    
    //+++++++++++++  average and std values in a ROI  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-ROI_Av_Std_values") == 0) {
      ROI_Av_Std_values(argc,argv);
      done = true;
    }
    //+++++++++++++   surfacic average and std values  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Surf_Av_Std_values") == 0) {
      Surf_Av_Std_values(argc,argv);
      done = true;
    }
    //+++++++++++++   cut a region of interest  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-ROI_Cut") == 0) {
      ROI_Cut(argc,argv);
      done = true;
    }
    //+++++++++++++   undesample a very large image  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-LargeImageUndersample") == 0) {
      LargeImageUndersample(argc,argv);
      done = true;
    }
    //+++++++++++++   Grey level alignment  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-GreyLevAlign") == 0) {
      GreyLevAlign(argc,argv);
      done = true;
    }
    //+++++++++++++   undersample the image resolution  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-UnderSample") == 0) {
      UnderSample(argc,argv);
      done = true;
    }
    //+++++++++++++   threshold on the lowest and highest grey level  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-GreyLevThresh") == 0) {
      GreyLevThresh(argc,argv);
      done = true;
    }
    //+++++++++++++   linearly resample the grey levels  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-GreyLevResample") == 0) {
      GreyLevResample(argc,argv);
      done = true;
    }
    //+++++++++++++   create a new void image  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-CreateVoidImage") == 0) {
      CreateVoidImage(argc,argv);
      done = true;
    }
    //+++++++++++++   linearly resample the grey levels  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-ComposeDispFields") == 0) {
      ComposeDispFields(argc,argv);
      done = true;
    }
    //+++++++++++++      Set the image to world matrix      +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-SetImage2World") == 0) {
      SetImage2World(argc,argv);
      done = true;
    }
    //+++++++++++++      Set the center of world coordinate s+++++++++++++++
    if (done == false) if (strcmp(argv[1], "-DefineCenterWorldCoord") == 0) {
      DefineCenterWorldCoord(argc,argv);
      done = true;
    }
    //+++++++++++++   Extract a 2D slice from a 3D image  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-Extract2Dslice") == 0) {
      Extract2Dslice(argc,argv);
      done = true;
    }
    //+++++++++++++   transform an image with a displacement field  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-TransfoImag") == 0) {
      TransfoImag(argc,argv);
      done = true;
    }
    //+++++++++++++   affine transformation of an image   +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-TransfoImagAffine") == 0) {
      TransfoImagAffine(argc,argv);
      done = true;
    }
    //+++++++++++++   affine transformation of an image   +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-TransfoImagAffine2") == 0) {
      TransfoImagAffine2(argc,argv);
      done = true;
    }
    //+++++++++++++   affine transformation of an image   +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-TransfoImagAffine3") == 0) {
      TransfoImagAffine3(argc,argv);
      done = true;
    }
    //+++++++++++++   transform an image with a displacement field  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-RemAffineContribDF") == 0) {
      RemAffineContribDF(argc,argv);
      done = true;
    }
    //+++  generate a displacement field using an affine deformation  +++
    if (done == false) if (strcmp(argv[1], "-GenerAffineDefDF") == 0) {
      GenerAffineDefDF(argc,argv);
      done = true;
    }
    //+++  generate a displacement field using an affine deformation  +++
    if (done == false) if (strcmp(argv[1], "-InvertAffineDef") == 0) {
      InvertAffineDef(argc,argv);
      done = true;
    }
    //+++++++++++++   transform the image2world values of an image to be centered on another image   +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-RegCenterImage") == 0) {
      RegCenterImage(argc,argv);
      done = true;
    }
    //+++++++++++++   make a sequence of 3D or 2D images   +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-MakeSequence") == 0) {
      MakeSequence(argc,argv);
      done = true;
    }
    //+++++++++++++   evaluate the accuracy of a displacement field with landmarks   +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-CompareLandmarkPoints") == 0) {
      Compare_LDMK_Points(argc,argv);
      done = true;
    }
    //+++++++++++++   other stuffs  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-other") == 0) {
      otherStuff(argc,argv);
      done = true;
    }
    //+++++++++++++   all options  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-AllOptions") == 0) {
      AllOptions(argc,argv);
      done = true;
    }
    //+++++++++++++   all options (bis)  +++++++++++++++
    if (done == false) if (strcmp(argv[1], "-A") == 0) {
      AllOptions(argc,argv);
      done = true;
    }
    
  }
  
  
  return 0;
}
