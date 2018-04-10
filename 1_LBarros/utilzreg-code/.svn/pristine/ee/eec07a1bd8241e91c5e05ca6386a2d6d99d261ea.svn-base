/*=========================================================================
 
Author: Laurent Risser

Disclaimer: This software has been developed for research purposes only, and hence should 
not be used as a diagnostic tool. In no event shall the authors or distributors
be liable to any direct, indirect, special, incidental, or consequential 
damages arising of the use of this software, its documentation, or any 
derivatives thereof, even if the authors have been advised of the possibility 
of such damage. 
 
 =========================================================================*/


#include <Demons_plain.h>

void usage(){
  cerr << "Usage: uTIlzReg_Demons_plain [Moving] [Fixed] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "  Primary options:\n";
  cerr << "    <-iterations n>          Number of iterations (default=20)\n";
  cerr << "    <-UnderSample n>         Undersample the fixed image with the factor n (default n = 1)\n";
  cerr << "  Inputs and Outputs:\n";
  cerr << "    <-PrefixInputs n>        Prefix of the files containing an initial velocity field (default=\"Null\")\n";
  cerr << "    <-PrefixOutputs n>       Prefix of the files containing the outputs (default=\"Outputs\")\n";
  cerr << "    <-affineT n>             Affine transfo from fixed image to moving image in the world domain. The 4*3 parameters are: r_xx r_xy r_xz t_x  r_yx ... t_z\n";
  cerr << "    <-affineT_txt n>         Affine transfo from fixed image to moving image in the world domain. The 4*4 matrix is an ascii text file.\n";
  cerr << "    <-IniDispF n>            Initial displacement field from fixed to moving in the world domain. (n = DX.nii DY.nii DZ.nii / affineT not considered)\n";
  cerr << "  Kernels    (Gaussian kernel of std dev 's'):\n";
  cerr << "    <-Gauss_fluid s>         Gaussian kernel for the fluid regularisation (default=[8])\n";
  cerr << "    <-Gauss_diffusion s>     Gaussian kernel for the diffusion reg (0 -> no dif.) (default=[0.0])\n";
  cerr << "  Secondary options:\n";
  cerr << "    <-MI>                    Minimize the mutual information instead of the sum of squared differences\n";
  cerr << "    <-MaxVeloUpdt n>         Maximum update at each iteration (default=0.5 voxels)\n";
  cerr << "    <-lambdaX n>             Value of lambdaX (default=1)\n";
  cerr << "  \n";
  
  exit(1);
}

int main(int argc, char **argv){
  LargeDefDemons_plain DefDemons;
  bool ok;
  int tmp,i;
  char TempChars[256];
  
  
  // Check command line
  if (argc < 3) {
    usage();
  }
  
  // Read the name of input and output images
  strcpy(DefDemons.MovingImFile,argv[1]);
  argc--;  argv++;
  strcpy(DefDemons.FixedImFile,argv[1]);
  argc--;  argv++;
  
  // Parse remaining parameters
  while (argc > 1) {
    //1 - Primary options
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--; argv++;
      DefDemons.iteration_nb = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Gauss_diffusion") == 0)){
      argc--; argv++;
      DefDemons.sigmaDiff = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Gauss_fluid") == 0)) {
      argc--; argv++;
      DefDemons.sigmaFluid = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-MaxVeloUpdt") == 0)) {
      argc--; argv++;
      DefDemons.MaxUpdateAllowed = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-lambdaX") == 0)) {
      argc--; argv++;
      DefDemons.lambdaX = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-MI") == 0)) {
      argc--; argv++;
      DefDemons.IndicatorMI = 1;
      ok = true;
    }
    //2 - Inputs and Outputs
    if ((ok == false) && (strcmp(argv[1], "-PrefixInputs") == 0)) {
      argc--; argv++;
      strcpy(DefDemons.PrefixInputs,argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-PrefixOutputs") == 0)) {
      argc--; argv++;
      strcpy(DefDemons.PrefixOutputs,argv[1]);
      argc--; argv++;
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-affineT") == 0)) {
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[0][0]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[0][1]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[0][2]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[0][3]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[1][0]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[1][1]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[1][2]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[1][3]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[2][0]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[2][1]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[2][2]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[2][3]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_FixedIm2MovingIm[3][0]=0;
      DefDemons.World_FixedIm2MovingIm[3][1]=0;
      DefDemons.World_FixedIm2MovingIm[3][2]=0;
      DefDemons.World_FixedIm2MovingIm[3][3]=1;
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-affineT_txt") == 0)) {
      argc--; argv++;
      strcpy(TempChars,argv[1]);
      argc--; argv++;
      
      Read_quat4t4mat(TempChars,DefDemons.World_FixedIm2MovingIm);
      DefDemons.World_FixedIm2MovingIm[3][0]=0;
      DefDemons.World_FixedIm2MovingIm[3][1]=0;
      DefDemons.World_FixedIm2MovingIm[3][2]=0;
      DefDemons.World_FixedIm2MovingIm[3][3]=1;
      
      ok = true;
    }

    
    //4 - Secondary options
    if ((ok == false) && (strcmp(argv[1], "-IniDispF") == 0)) {
      argc--; argv++;
      strcpy(DefDemons.IniDispFieldX,argv[1]);
      argc--; argv++;
      strcpy(DefDemons.IniDispFieldY,argv[1]);
      argc--; argv++;
      strcpy(DefDemons.IniDispFieldZ,argv[1]);
      argc--; argv++;
      DefDemons.IniDispFieldDefined=1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-UnderSample") == 0)) {
      argc--; argv++;
      DefDemons.UnderSampleFixedImFactor = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    
    if (ok == false) usage();
  }
  
  
  //run process
  cout << "Registration go"<< endl;
  DefDemons.Run();
  cout << "done" << endl;
  
  return 0;
}
