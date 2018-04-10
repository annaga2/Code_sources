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

void usage(){
  cerr << "Usage: uTIlzReg_Demons [Source] [Target] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "  Primary options:\n";
  cerr << "    <-iterations n>          Number of iterations (default=10)\n";
  cerr << "    <-UnderSampleTrg n>      Undersample the target image with the factor n (default n = 1)\n";
  cerr << "  Inputs and Outputs:\n";
  cerr << "    <-AddChannel w S T>      Add a channel, where w is its weight (w=1 for the 1st channel) and S,T are the src and trg channels \n";
  cerr << "    <-PrefixInputs n>        Prefix of the files containing an initial velocity field (default=\"Null\")\n";
  cerr << "    <-PrefixOutputs n>       Prefix of the files containing the outputs (default=\"Outputs\")\n";
  cerr << "    <-affineT n>             Affine transfo from Trg to Src in the world domain. The 4*3 parameters are: r_xx r_xy r_xz t_x  r_yx ... t_z\n";
  cerr << "    <-affineT_txt n>         Affine transfo from Trg to Src in the world domain. The 4*4 matrix is an ascii text file.\n";
  cerr << "    <-IniDispF n>            Initial displacement field from Trg to Src in the world domain. (n = DX.nii DY.nii DZ.nii / affineT not considered)\n";
  cerr << "    <-FinalDefVec>           Displacement field in mm from [Source] to [Target]   (only the inverse deformation is saved by default)\n";
  cerr << "  Kernels    (Gaussian kernel of std dev 's'):\n";
  cerr << "    <-Gauss_diffusion s>     Gaussian kernel for the diffusion reg (0 -> no dif.) (default=[0.0])\n";
  cerr << "    <-Gauss_fluid s>         Gaussian kernel for the fluid regularisation (default=[2])\n";
  cerr << "    <-MK_fluid N s1 ... sN>  Fluid regularisation using the sum of N gaussian kernels of std dev s1 ... sN  (as in [Risser et al, TMI 2011])\n";
  cerr << "  Secondary options:\n";
  cerr << "    <-MI>                    Minimize the mutual information instead of the sum of squared differences      (as in [Risser et al, ISBI 2012])\n";
  cerr << "    <-MaxVeloUpdt n>         Maximum update at each iteration (default=0.5 voxels)\n";
  cerr << "    <-Mask n>                Demons only considered where the mask is equal to ID  (n= [Mask.nii] ID / default=\"Null\")    \n";
  cerr << "    <-DiscMask mask l N>     Piecewise-diffeo. reg. where mask==l / sliding motion in a margin M (ex: M=30mm)   ([Risser et al, MedIA 2012])\n";
  cerr << "    <-ExpendDomain n>        Extend the domain of the target image, where the computations are done (n = x1 x2 y1 y2 z1 z2 / default = null)\n";
  cerr << "    <-SaveMemory>            Save as much memory as possible (should be only used on very large images to avoid swapping)\n";
  cerr << "    <-GreyLevAlign n>        Grey level linear alignment (Inputs: Padding Src - Padding Trg)\n";
  cerr << "    <-margins n>             Margin of the image where the calculations are reduced  (default=0 voxels)\n";
  cerr << "    <-lambdaX n>             Value of lambdaX (default=1)\n";
  cerr << "  \n";
  
  exit(1);
}

int main(int argc, char **argv){
  LargeDefDemons DefDemons;
  bool ok;
  int tmp,i;
  char TempChars[256];
  
  
  // Check command line
  if (argc < 3) {
    usage();
  }
  
  // Read the name of input and output images
  strcpy(DefDemons.SourceFiles[DefDemons.NbChannels],argv[1]);
  argc--;  argv++;
  strcpy(DefDemons.TargetFiles[DefDemons.NbChannels],argv[1]);
  argc--;  argv++;
  
  DefDemons.weightsChannels[DefDemons.NbChannels]=1;
  DefDemons.NbChannels++;
  
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
      DefDemons.sigma1 = atof(argv[1]);
      argc--; argv++;
      ok = true;
      DefDemons.w1=100;
    }
        if ((ok == false) && (strcmp(argv[1], "-MK_fluid") == 0)) {
      argc--; argv++;
      DefDemons.NbKernelScales= atoi(argv[1]);
      argc--; argv++;
      DefDemons.w1=0;  // will launch the automatic tuning of the weights
      if (DefDemons.NbKernelScales>=1){
        DefDemons.sigma1 = atof(argv[1]);
        argc--; argv++;
      }
      if (DefDemons.NbKernelScales>=2){
        DefDemons.sigma2 = atof(argv[1]);
        argc--; argv++;
      }
      if (DefDemons.NbKernelScales>=3){
        DefDemons.sigma3 = atof(argv[1]);
        argc--; argv++;
      }
      if (DefDemons.NbKernelScales>=4){
        DefDemons.sigma4 = atof(argv[1]);
        argc--; argv++;
      }
      if (DefDemons.NbKernelScales>=5){
        DefDemons.sigma5 = atof(argv[1]);
        argc--; argv++;
      }
      if (DefDemons.NbKernelScales>=6){
        DefDemons.sigma6 = atof(argv[1]);
        argc--; argv++;
      }
      if (DefDemons.NbKernelScales>=7){
        DefDemons.sigma7 = atof(argv[1]);
        argc--; argv++;
      }
      ok = true;
      DefDemons.IndicatorFFT_fluid = 1;
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
    if ((ok == false) && (strcmp(argv[1], "-SaveMemory") == 0)) {
      argc--; argv++;
      DefDemons.IndicatorSaveMemory = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-MI") == 0)) {
      argc--; argv++;
      DefDemons.IndicatorMI = 1;
      ok = true;
      if (DefDemons.NbChannels > 1) cout << "Gradients of mutual information only programmed for one channel!" << endl;
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
    if ((ok == false) && (strcmp(argv[1], "-AddChannel") == 0)) {
      argc--; argv++;
      DefDemons.weightsChannels[DefDemons.NbChannels]=atof(argv[1]);
      argc--;  argv++;
      strcpy(DefDemons.SourceFiles[DefDemons.NbChannels],argv[1]);
      argc--;  argv++;
      strcpy(DefDemons.TargetFiles[DefDemons.NbChannels],argv[1]);
      argc--;  argv++;
      DefDemons.NbChannels++;
      ok = true;
      if (DefDemons.IndicatorMI == 1) cout << "Gradients of mutual information only programmed for one channel!" << endl;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-affineT") == 0)) {
      argc--; argv++;
      DefDemons.World_Target2Template[0][0]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[0][1]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[0][2]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[0][3]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[1][0]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[1][1]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[1][2]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[1][3]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[2][0]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[2][1]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[2][2]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[2][3]=atof(argv[1]);
      argc--; argv++;
      DefDemons.World_Target2Template[3][0]=0;
      DefDemons.World_Target2Template[3][1]=0;
      DefDemons.World_Target2Template[3][2]=0;
      DefDemons.World_Target2Template[3][3]=1;
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-affineT_txt") == 0)) {
      argc--; argv++;
      strcpy(TempChars,argv[1]);
      argc--; argv++;
      
      Read_quat4t4mat(TempChars,DefDemons.World_Target2Template);
      DefDemons.World_Target2Template[3][0]=0;
      DefDemons.World_Target2Template[3][1]=0;
      DefDemons.World_Target2Template[3][2]=0;
      DefDemons.World_Target2Template[3][3]=1;
      
      ok = true;
    }

    
    //4 - Secondary options
    if ((ok == false) && (strcmp(argv[1], "-GreyLevAlign") == 0)) {
      argc--; argv++;
      DefDemons.GreyLevAlign = 1;
      DefDemons.GLA_Padding_Src = atof(argv[1]);
      argc--; argv++;
      DefDemons.GLA_Padding_Trg = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-margins") == 0)) {
      argc--; argv++;
      DefDemons.Margin = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-Mask") == 0)) {
      argc--; argv++;
      strcpy(DefDemons.MaskFile,argv[1]);
      DefDemons.MaskDefined = 1;
      argc--; argv++;
      DefDemons.MaskID = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-DiscMask") == 0)) {
      argc--; argv++;
      strcpy(DefDemons.MaskFile,argv[1]);
      argc--; argv++;
      DefDemons.MaskLabel = atof(argv[1]);
      argc--; argv++;
      DefDemons.BoundaMargin = atof(argv[1]);
      DefDemons.MaskDefined = 2;
      argc--; argv++;
      ok = true;
    }
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
    if ((ok == false) && (strcmp(argv[1], "-FinalDefVec") == 0)) {
      argc--; argv++;
      DefDemons.FinalDefVec = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-UnderSampleTrg") == 0)) {
      argc--; argv++;
      DefDemons.UnderSampleTrgFactor = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ExpendDomain") == 0)) {
      argc--; argv++;
      DefDemons.ExtendTrgImag_LowerX = atoi(argv[1]);
      argc--; argv++;
      DefDemons.ExtendTrgImag_UpperX = atoi(argv[1]);
      argc--; argv++;
      DefDemons.ExtendTrgImag_LowerY = atoi(argv[1]);
      argc--; argv++;
      DefDemons.ExtendTrgImag_UpperY = atoi(argv[1]);
      argc--; argv++;
      DefDemons.ExtendTrgImag_LowerZ = atoi(argv[1]);
      argc--; argv++;
      DefDemons.ExtendTrgImag_UpperZ = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
    
    if (ok == false) usage();
  }
  
  
  //run process
  cout << "Registration using ";
  if (DefDemons.IndicatorMI==1) cout << "mutual information-based ";
  cout << "LogDemons " << " ... \n"; cout.flush();
  DefDemons.Run();
  cout << "done" << endl;
  
  return 0;
}
