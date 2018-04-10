/*=========================================================================

Authors: Francois-Xavier Vialard, Laurent Risser

Disclaimer: This software has been developed for research purposes only, and hence should 
not be used as a diagnostic tool. In no event shall the authors or distributors
be liable to any direct, indirect, special, incidental, or consequential 
damages arising of the use of this software, its documentation, or any 
derivatives thereof, even if the authors have been advised of the possibility 
of such damage. 

=========================================================================*/

#include <GeoShoot.h>

void usage(){
  cerr << "Usage: uTIlzReg_GeoShoot [Template] [Target] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "  Primary options:\n";
  cerr << "    <-iterations n>             Number of iterations (default=10)\n";
  cerr << "    <-subdivisions n>           Number of subdivisons between t=0 and t=1 (default=10)\n";
  cerr << "    <-UnderSampleTpl n>         Undersample the template image with a factor n (default n = 1)\n";
  cerr << "    <-alpha n>                  Weight of the norm in the cost function (Default=0.001)\n";
  cerr << "  Kernels (Default: -Gauss 3):\n";
  cerr << "    <-Gauss n>                  Gaussian kernel of std. dev. Sigma (in mm)\n";
  cerr << "    <-M_Gauss n>                Sum of Gaussian kernels (max 7)   --- n = k W1 S1 ... Wk Sk   (k=[#kernels], W.=weight, S.=Sigma in mm)\n" ;
  cerr << "    <-M_Gauss_easier n>         Sum of 7 linearly sampled Gaussian kernels with apparent weights = 1    --- n = Smax Smin  (S. in mm)\n" ;
  cerr << "  Inputs (default: nothing):\n";
  cerr << "    <-PrefixInputs n>           Prefix of the file containing the initial momentum (default=\"Null\")\n";
  cerr << "    <-InputInitMomentum n>      Initial Momentum to initiate the gradient descent (default=\"Null\")\n";
  cerr << "    <-InputIniMoTxt n>          Initial Momentum in an ascii file (instead of nifti) (default=\"Null\")\n";
  cerr << "    <-affineT n>                Affine transfo from Trg to Template in the world domain. The 4*3 parameters are: r_xx r_xy r_xz t_x  r_yx ... t_z\n";
  cerr << "    <-affineT_txt n>            Affine transfo from Trg to Template in the world domain. The 4*4 matrix is an ascii text file.\n";
  cerr << "  Outputs (default: initial momentum in a nifti file):\n";
  cerr << "    <-PrefixOutputs n>          Prefix of the output files (default=\"Outputs\")\n";
  cerr << "    <-OutFinalDef>              Outputs the deformed template (Nothing by default)\n";
  cerr << "    <-OutDispField>             Outputs the displacement field in mm (Nothing by default)\n";
  cerr << "    <-OutDispFieldEvo>          Outputs the evolution of the displacement field in mm along the diffeomorphism (Nothing by default)\n";
  cerr << "    <-OutIniMoTxt n>            Outputs the initial momentum in an ascii file (default=\"Null\")\n";
  cerr << "    <-OutVeloField>             Outputs the 3D+t velocity field in voxels (Nothing by default)\n";
  cerr << "    <-OutDistEnSim>             Outputs the distance, enrgy and similarity measure (Nothing by default)\n";
  cerr << "    <-OutDeformation>           Outputs the 3D+t deformation (Nothing by default)\n";
  cerr << "    <-OutDiff_Tpl_DefTrg>       Outputs the difference between the template and the deformed target image (Nothing by default)\n";
  cerr << "  Secondary options:\n";
  cerr << "    <-MaxVeloUpdt n>            Size of the maximum updates of the vector field (Default=0.5 voxels)\n";
  cerr << "    <-margins n>                Margin of the image where the calculations are reduced  (default=3 voxels)\n";
  cerr << "    <-GreyLevAlign n>           Grey level linear alignment of each channel -- n = [Padding Src] [Padding Trg]\n";
  cerr << "    <-Mask n>                   Mask in which the momenta are computed (values!=0  -- in the template image domain)\n";
  exit(1);
}

int main(int argc, char **argv){
  bool ok;
  EulerianShooting Shoot;
  int temp;
  char TempChars[256];

  // Check command line
  if (argc < 3){
    usage();
  }

  // read mandatory parameters
  strcpy(Shoot.SourceImageName,argv[1]);
  argc--;  argv++;
  strcpy(Shoot.TargetImageName,argv[1]);
  argc--;  argv++;
  
  // Parse remaining parameters
  while (argc > 1) {
    //1 - Primary options
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-subdivisions") == 0)) {
      argc--; argv++;
      Shoot.NbTimeSubdivisions = atoi(argv[1]);
      if (Shoot.NbTimeSubdivisions<2) Shoot.NbTimeSubdivisions=2;
      argc--; argv++;
      ok = true;
    }
  if ((ok == false) && (strcmp(argv[1], "-PrefixOutputs") == 0)) {
      argc--; argv++;
      strcpy(Shoot.PrefixOutputs,argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-PrefixInputs") == 0)) {
      argc--; argv++;
      strcpy(Shoot.PrefixInputs,argv[1]);
      Shoot.indicatorInitialMomentum = 1;
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--; argv++;
      Shoot.NbIterations = atoi(argv[1]);
      if (Shoot.NbIterations<0) Shoot.NbIterations=0;
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-UnderSampleTpl") == 0)) {
      argc--; argv++;
      Shoot.UnderSampleTplFactor = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-InputInitMomentum") == 0)) {
      argc--; argv++;
      strcpy(Shoot.InputInitialMomentumName,argv[1]);
      Shoot.indicatorInitialMomentum = 1;
      argc--; argv++;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-InputIniMoTxt") == 0)) {
      argc--; argv++;
      strcpy(Shoot.InputInitialMomentumName,argv[1]);
      Shoot.indicatorInitialMomentum = 2;
      argc--; argv++;
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-OutIniMoTxt") == 0)) {
      argc--; argv++;
      Shoot.OutIniMoTxt = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutVeloField") == 0)) {
      argc--; argv++;
      Shoot.OutVeloField = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutDistEnSim") == 0)) {
      argc--; argv++;
      Shoot.OutDistEnSim = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutDeformation") == 0)) {
      argc--; argv++;
      Shoot.OutDeformation = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutDiff_Tpl_DefTrg") == 0)) {
      argc--; argv++;
      Shoot.OutDiff_Tpl_DefTrg = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutFinalDef") == 0)) {
      argc--; argv++;
      Shoot.OutFinalDef = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutDispField") == 0)) {
      argc--; argv++;
      Shoot.OutDispField = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-OutDispFieldEvo") == 0)) {
      argc--; argv++;
      Shoot.OutDispFieldEvo = 1;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-affineT") == 0)) {
      argc--; argv++;
      Shoot.World_Target2Template[0][0]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[0][1]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[0][2]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[0][3]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[1][0]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[1][1]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[1][2]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[1][3]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[2][0]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[2][1]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[2][2]=atof(argv[1]);
      argc--; argv++;
      Shoot.World_Target2Template[2][3]=atof(argv[1]);
      argc--; argv++;

      Shoot.World_Target2Template[3][0]=0;
      Shoot.World_Target2Template[3][1]=0;
      Shoot.World_Target2Template[3][2]=0;
      Shoot.World_Target2Template[3][3]=1;
      ok = true;
    }

    if ((ok == false) && (strcmp(argv[1], "-affineT_txt") == 0)) {
      argc--; argv++;
      strcpy(TempChars,argv[1]);
      argc--; argv++;
      
      Read_quat4t4mat(TempChars,Shoot.World_Target2Template);
      Shoot.World_Target2Template[3][0]=0;
      Shoot.World_Target2Template[3][1]=0;
      Shoot.World_Target2Template[3][2]=0;
      Shoot.World_Target2Template[3][3]=1;
      
      ok = true;
    }
    
    
    
    if ((ok == false) && (strcmp(argv[1], "-Gauss") == 0)) {
      argc--; argv++;
      Shoot.sigmaX1 = atof(argv[1]);
      Shoot.sigmaY1 = atof(argv[1]);
      Shoot.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-M_Gauss") == 0)) {
      argc--; argv++;
      temp= atoi(argv[1]);
      argc--; argv++;
      if (temp>=1){
        Shoot.weight1 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX1 = atof(argv[1]); Shoot.sigmaY1 = atof(argv[1]); Shoot.sigmaZ1 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=2){
        Shoot.weight2 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX2 = atof(argv[1]); Shoot.sigmaY2 = atof(argv[1]); Shoot.sigmaZ2 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=3){
        Shoot.weight3 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX3 = atof(argv[1]); Shoot.sigmaY3 = atof(argv[1]); Shoot.sigmaZ3 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=4){
        Shoot.weight4 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX4 = atof(argv[1]); Shoot.sigmaY4 = atof(argv[1]); Shoot.sigmaZ4 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=5){
        Shoot.weight5 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX5 = atof(argv[1]); Shoot.sigmaY5 = atof(argv[1]); Shoot.sigmaZ5 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=6){
        Shoot.weight6 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX6 = atof(argv[1]); Shoot.sigmaY6 = atof(argv[1]); Shoot.sigmaZ6 = atof(argv[1]);
        argc--; argv++;
      }
      if (temp>=7){
        Shoot.weight7 = atof(argv[1]);
        argc--; argv++;
        Shoot.sigmaX7 = atof(argv[1]); Shoot.sigmaY7 = atof(argv[1]); Shoot.sigmaZ7 = atof(argv[1]);
        argc--; argv++;
      }
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-M_Gauss_easier") == 0)) {
      cout << "in" << endl;
      argc--; argv++;
      Shoot.sigmaX1 = atof(argv[1]);
      Shoot.sigmaY1 = atof(argv[1]);
      Shoot.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      Shoot.sigmaX7 = atof(argv[1]);
      Shoot.sigmaY7 = atof(argv[1]);
      Shoot.sigmaZ7 = atof(argv[1]);
      argc--; argv++;
      
      Shoot.weight1 = 0; Shoot.weight2 = 0; Shoot.weight3 = 0; Shoot.weight4 = 0;
      Shoot.weight5 = 0; Shoot.weight6 = 0; Shoot.weight7 = 0;

      float tmpFl;
      float a,b;
      int i;
      
      if (Shoot.sigmaX1<Shoot.sigmaX7){
        tmpFl=Shoot.sigmaX1; Shoot.sigmaX1=Shoot.sigmaX7; Shoot.sigmaX7=tmpFl;
        tmpFl=Shoot.sigmaY1; Shoot.sigmaY1=Shoot.sigmaY7; Shoot.sigmaY7=tmpFl;
        tmpFl=Shoot.sigmaZ1; Shoot.sigmaZ1=Shoot.sigmaZ7; Shoot.sigmaZ7=tmpFl;
      }
      
      a=(Shoot.sigmaY7-Shoot.sigmaY1)/6.;
      b=Shoot.sigmaY1-a;
      
      Shoot.sigmaX2 = a*2+b; Shoot.sigmaY2 = a*2+b; Shoot.sigmaZ2 = a*2+b;
      Shoot.sigmaX3 = a*3+b; Shoot.sigmaY3 = a*3+b; Shoot.sigmaZ3 = a*3+b;
      Shoot.sigmaX4 = a*4+b; Shoot.sigmaY4 = a*4+b; Shoot.sigmaZ4 = a*4+b;
      Shoot.sigmaX5 = a*5+b; Shoot.sigmaY5 = a*5+b; Shoot.sigmaZ5 = a*5+b;
      Shoot.sigmaX6 = a*6+b; Shoot.sigmaY6 = a*6+b; Shoot.sigmaZ6 = a*6+b;
      
      ok = true;
      cout << "out" << endl;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-alpha") == 0)) {
      argc--; argv++;
      Shoot.alpha = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-MaxVeloUpdt") == 0)) {
      argc--; argv++;
    Shoot.MaxUpdate = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
  if ((ok == false) && (strcmp(argv[1], "-margins") == 0)) {
      argc--; argv++;
      Shoot.Margin = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
  if ((ok == false) && (strcmp(argv[1], "-GreyLevAlign") == 0)) {
      argc--; argv++;
      Shoot.GreyLevAlign = 1;
      Shoot.GLA_Padding_Src = atof(argv[1]);
      argc--; argv++;
      Shoot.GLA_Padding_Trg = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
  if ((ok == false) && (strcmp(argv[1], "-Mask") == 0)) {
      argc--; argv++;
      strcpy(Shoot.InputMask,argv[1]);
      Shoot.indicatorMask = 1;
      argc--; argv++;
      ok = true;
    }

  if (ok == false) 
  {
    usage();
  }
  }
  
  //check the inputs
  printf("\nPARAMETERS:\n");
  printf("  Source image:           %s\n",Shoot.SourceImageName);
  printf("  Target image:           %s\n",Shoot.TargetImageName);
  printf("  Subdivisions number:    %d\n",Shoot.NbTimeSubdivisions);
  printf("  Gaussian Std. dev.:     %lf\n\n",Shoot.sigmaX1);
  
  
  //run process
  cout << "Large deformation registration using geodesic shooting... \n"; cout.flush();
  Shoot.Run();
  cout << "done" << endl;

  return 0;
}


