/*=========================================================================
 
 
Author: Laurent Risser, Francois-Xavier Vialard

Disclaimer: This software has been developed for research purposes only, and hence should 
not be used as a diagnostic tool. In no event shall the authors or distributors
be liable to any direct, indirect, special, incidental, or consequential 
damages arising of the use of this software, its documentation, or any 
derivatives thereof, even if the authors have been advised of the possibility 
of such damage. 


 =========================================================================*/

#include <LDDMM.h>

void usage(){
  cerr << "Usage: uTIlzReg_LDDMM_Beg [Source] [Target] <options>\n";
  cerr << "Where <options> are one or more of the following:\n";
  cerr << "  Primary options:\n";
  cerr << "    <-iterations n>      Number of iterations (default=10)\n";
  cerr << "    <-subdivisions n>    Number of subdivisons (default=10)\n";
  cerr << "    <-MaxVeloUpdt n>     Maximum velocity update at each iteration (default=0.4 voxels)\n";
  cerr << "  Inputs and Outputs:\n";
  cerr << "    <-PrefixInputs n>    Prefix of the files containing an initial velocity field (default=\"Null\")\n";
  cerr << "    <-PrefixOutputs n>   Prefix of the files containing the outputs (default=\"Outputs\")\n";
  cerr << "    <-AddChannel W S T>  Add a channel -> W=weight (wgt of ref channel is 1) S=Source T=Target\n";
  cerr << "    <-Mask n>            Definition of a mask (default=\"Null\")\n";
  cerr << "    <-affineT n>         Affine transfo from Trg to Src in the world domain. The 4*3 parameters are: r_xx r_xy r_xz t_x  r_yx ... t_z\n";
  cerr << "    <-affineT_txt n>     Affine transfo from Trg to Src in the world domain. The 4*4 matrix is an ascii text file.\n";
  cerr << "  Kernels (Default: -Gauss 1):\n";
  cerr << "    <-Gauss S>           Gaussian kernel (S = std. dev. in mm)\n";
  cerr << "    <-M_Gauss n>         Sum of Gaussian kernels (max 7) -- n = k W1 S1 ... Wk Sk   (k=[#kernels], W.=weight)\n" ;
  cerr << "    <-M_Gauss_easy n>    Sum of Gaussian kernels (max 7) with apparent weights = 1 -- n = k S1 ... Sk   (k=[#kernels])\n" ;
  cerr << "    <-M_Gauss_easier n>  Sum of 7 linearly sampled Gaussian kernels with apparent weights = 1 -- n = Smax Smin\n" ;
  cerr << "  Secondary options:\n";
  cerr << "    <-symmetric>         Symmetric registration\n";
  cerr << "    <-epsilon n>         Threshold on the normalized max update of the velicty field (default=0.2)\n";
  cerr << "    <-GreyLevAlign n>    Grey level linear alignment of each channel (Inputs: Padding Src - Padding Trg)\n";
  cerr << "    <-margins n>         Margin of the image where the calculations are reduced  (default=0 voxels)\n";
  cerr << "    <-WghtVeloField n>   Weight of the velocity field in the energy (default=0.001) \n";
  cerr << "    <-RefMaxGrad n>      Value to manage the convergence. Automatically configured if <0 (default=-1.)\n";
  cerr << "  Special Outputs:\n";
  cerr << "    <-FinalDefVec>       Displacement field in mm from [Source] to [Target]   (by default for the whole deformation)\n";
  cerr << "    <-FinalDefInvVec>    Displacement field in mm from [Target] to [Source]\n";
  cerr << "    <-SplitKernels>      Split the contribution of each kernel\n";
  cerr << "    <-AOD>               Amplitude of the deformations from each voxel of the source image\n";
  cerr << "    <-DetJacobian>       Determinant of the Jacobian at each voxel\n";
  cerr << "    <-InitMomentum>      Estimated initial momentum\n";
  cerr << "    <-ShowSSD>           Show the Sum of the Squared Differences at t=1 ieration after iteration\n";
  cerr << "  \n";
  cerr << "  Or: Typical inverse amplitude of the updates with a Gaussian kernel of stddev Sigma: <-TypicInvAmp Sigma>.  \n";
  cerr << "   \n";
  
  exit(1);
}

int main(int argc, char **argv){
  LargeDefGradLagrange LargeDef;
  bool ok;
  char TempChars[256];
  
  // Check command line
  if (argc < 3) {
    usage();
  }
  
  // Read the name of input and output images (= 1st channel)
  strcpy(LargeDef.SourceFiles[LargeDef.NbChannels],argv[1]);
  argc--;  argv++;
  strcpy(LargeDef.TargetFiles[LargeDef.NbChannels],argv[1]);
  argc--;  argv++;
  LargeDef.NbChannels++;
  
  // Parse remaining parameters
  while (argc > 1) {
    //1 - Primary options
    ok = false;
    if ((ok == false) && (strcmp(argv[1], "-iterations") == 0)) {
      argc--; argv++;
      LargeDef.iteration_nb = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-subdivisions") == 0)) {
      argc--; argv++;
      LargeDef.NbTimeSubdiv = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-MaxVeloUpdt") == 0)) {
      argc--; argv++;
      LargeDef.MaxVelocityUpdate = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    //2 - Inputs and Outputs
    if ((ok == false) && (strcmp(argv[1], "-PrefixInputs") == 0)) {
      argc--; argv++;
      strcpy(LargeDef.PrefixInputs,argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-PrefixOutputs") == 0)) {
      argc--; argv++;
      strcpy(LargeDef.PrefixOutputs,argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-AddChannel") == 0) && (LargeDef.NbChannels<100)) {
      argc--; argv++;
      LargeDef.weightChannel[LargeDef.NbChannels] = atof(argv[1]);
      argc--; argv++;
      strcpy(LargeDef.SourceFiles[LargeDef.NbChannels],argv[1]);
      argc--; argv++;
      strcpy(LargeDef.TargetFiles[LargeDef.NbChannels],argv[1]);
      argc--; argv++;
      LargeDef.NbChannels++;
      ok = true;
      if (LargeDef.NbChannels==100) cout << "\n \n MAXIMUM NUMBER OF 100 CHANNELS IS REACHED !!!\n \n ";
    }
    if ((ok == false) && (strcmp(argv[1], "-Mask") == 0)) {
      argc--; argv++;
      strcpy(LargeDef.MaskFile,argv[1]);
      argc--; argv++;
      ok = true;
    }
    //3 - Kernels
    if ((ok == false) && (strcmp(argv[1], "-Gauss") == 0)) {
      LargeDef.weight1 = 100.;
      argc--; argv++;
      LargeDef.sigmaX1 = atof(argv[1]);
      LargeDef.sigmaY1 = atof(argv[1]);
      LargeDef.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDef.NbKernels=1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-M_Gauss") == 0)) {
      argc--; argv++;
      LargeDef.NbKernels= atoi(argv[1]);
      argc--; argv++;
      if (LargeDef.NbKernels>=1){
        LargeDef.weight1 = atof(argv[1]);
        argc--; argv++;
        LargeDef.sigmaX1 = atof(argv[1]);
        LargeDef.sigmaY1 = atof(argv[1]);
        LargeDef.sigmaZ1 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=2){
        LargeDef.weight2 = atof(argv[1]);
        argc--; argv++;
        LargeDef.sigmaX2 = atof(argv[1]);
        LargeDef.sigmaY2 = atof(argv[1]);
        LargeDef.sigmaZ2 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=3){
        LargeDef.weight3 = atof(argv[1]);
        argc--; argv++;
        LargeDef.sigmaX3 = atof(argv[1]);
        LargeDef.sigmaY3 = atof(argv[1]);
        LargeDef.sigmaZ3 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=4){
        LargeDef.weight4 = atof(argv[1]);
        argc--; argv++;
        LargeDef.sigmaX4 = atof(argv[1]);
        LargeDef.sigmaY4 = atof(argv[1]);
        LargeDef.sigmaZ4 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=5){
        LargeDef.weight5 = atof(argv[1]);
        argc--; argv++;
        LargeDef.sigmaX5 = atof(argv[1]);
        LargeDef.sigmaY5 = atof(argv[1]);
        LargeDef.sigmaZ5 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=6){
        LargeDef.weight6 = atof(argv[1]);
        argc--; argv++;
        LargeDef.sigmaX6 = atof(argv[1]);
        LargeDef.sigmaY6 = atof(argv[1]);
        LargeDef.sigmaZ6 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=7){
        LargeDef.weight7 = atof(argv[1]);
        argc--; argv++;
        LargeDef.sigmaX7 = atof(argv[1]);
        LargeDef.sigmaY7 = atof(argv[1]);
        LargeDef.sigmaZ7 = atof(argv[1]);
        argc--; argv++;
      }
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-M_Gauss_easy") == 0)) {
      argc--; argv++;
      LargeDef.NbKernels= atoi(argv[1]);
      argc--; argv++;
      if (LargeDef.NbKernels>=1){
        LargeDef.weight1 = 0;
        LargeDef.sigmaX1 = atof(argv[1]);
        LargeDef.sigmaY1 = atof(argv[1]);
        LargeDef.sigmaZ1 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=2){
        LargeDef.weight2 = 0;
        LargeDef.sigmaX2 = atof(argv[1]);
        LargeDef.sigmaY2 = atof(argv[1]);
        LargeDef.sigmaZ2 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=3){
        LargeDef.weight3 = 0;
        LargeDef.sigmaX3 = atof(argv[1]);
        LargeDef.sigmaY3 = atof(argv[1]);
        LargeDef.sigmaZ3 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=4){
        LargeDef.weight4 = 0;
        LargeDef.sigmaX4 = atof(argv[1]);
        LargeDef.sigmaY4 = atof(argv[1]);
        LargeDef.sigmaZ4 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=5){
        LargeDef.weight5 = 0;
        LargeDef.sigmaX5 = atof(argv[1]);
        LargeDef.sigmaY5 = atof(argv[1]);
        LargeDef.sigmaZ5 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=6){
        LargeDef.weight6 = 0;
        LargeDef.sigmaX6 = atof(argv[1]);
        LargeDef.sigmaY6 = atof(argv[1]);
        LargeDef.sigmaZ6 = atof(argv[1]);
        argc--; argv++;
      }
      if (LargeDef.NbKernels>=7){
        LargeDef.weight7 = 0;
        LargeDef.sigmaX7 = atof(argv[1]);
        LargeDef.sigmaY7 = atof(argv[1]);
        LargeDef.sigmaZ7 = atof(argv[1]);
        argc--; argv++;
      }
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-M_Gauss_easier") == 0)) {
      argc--; argv++;
      LargeDef.sigmaX1 = atof(argv[1]);
      LargeDef.sigmaY1 = atof(argv[1]);
      LargeDef.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDef.sigmaX7 = atof(argv[1]);
      LargeDef.sigmaY7 = atof(argv[1]);
      LargeDef.sigmaZ7 = atof(argv[1]);
      argc--; argv++;
      
      LargeDef.NbKernels=7;
      
      LargeDef.weight1 = 0; LargeDef.weight2 = 0; LargeDef.weight3 = 0; LargeDef.weight4 = 0;
      LargeDef.weight5 = 0; LargeDef.weight6 = 0; LargeDef.weight7 = 0;

      float tmpFl;
      float a,b;
      int i;
      
      if (LargeDef.sigmaX1<LargeDef.sigmaX7){
        tmpFl=LargeDef.sigmaX1; LargeDef.sigmaX1=LargeDef.sigmaX7; LargeDef.sigmaX7=tmpFl;
        tmpFl=LargeDef.sigmaY1; LargeDef.sigmaY1=LargeDef.sigmaY7; LargeDef.sigmaY7=tmpFl;
        tmpFl=LargeDef.sigmaZ1; LargeDef.sigmaZ1=LargeDef.sigmaZ7; LargeDef.sigmaZ7=tmpFl;
      }
      
      a=(LargeDef.sigmaY7-LargeDef.sigmaY1)/6.;
      b=LargeDef.sigmaY1-a;
      
      LargeDef.sigmaX2 = a*2+b; LargeDef.sigmaY2 = a*2+b; LargeDef.sigmaZ2 = a*2+b;
      LargeDef.sigmaX3 = a*3+b; LargeDef.sigmaY3 = a*3+b; LargeDef.sigmaZ3 = a*3+b;
      LargeDef.sigmaX4 = a*4+b; LargeDef.sigmaY4 = a*4+b; LargeDef.sigmaZ4 = a*4+b;
      LargeDef.sigmaX5 = a*5+b; LargeDef.sigmaY5 = a*5+b; LargeDef.sigmaZ5 = a*5+b;
      LargeDef.sigmaX6 = a*6+b; LargeDef.sigmaY6 = a*6+b; LargeDef.sigmaZ6 = a*6+b;
      
      ok = true;
    }
    
    
    if ((ok == false) && (strcmp(argv[1], "-affineT") == 0)) {
      argc--; argv++;
      LargeDef.World_Target2Template[0][0]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[0][1]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[0][2]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[0][3]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[1][0]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[1][1]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[1][2]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[1][3]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[2][0]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[2][1]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[2][2]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[2][3]=atof(argv[1]);
      argc--; argv++;
      LargeDef.World_Target2Template[3][0]=0;
      LargeDef.World_Target2Template[3][1]=0;
      LargeDef.World_Target2Template[3][2]=0;
      LargeDef.World_Target2Template[3][3]=1;
      ok = true;
    }
    
    if ((ok == false) && (strcmp(argv[1], "-affineT_txt") == 0)) {
      argc--; argv++;
      strcpy(TempChars,argv[1]);
      argc--; argv++;
      
      Read_quat4t4mat(TempChars,LargeDef.World_Target2Template);
      LargeDef.World_Target2Template[3][0]=0;
      LargeDef.World_Target2Template[3][1]=0;
      LargeDef.World_Target2Template[3][2]=0;
      LargeDef.World_Target2Template[3][3]=1;
      
      ok = true;
    }

    
    //4 - Secondary options
    if ((ok == false) && (strcmp(argv[1], "-GreyLevAlign") == 0)) {
      argc--; argv++;
      LargeDef.GreyLevAlign = 1;
      LargeDef.GLA_Padding_Src = atof(argv[1]);
      argc--; argv++;
      LargeDef.GLA_Padding_Trg = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-SplitKernels") == 0)) {
      argc--; argv++;
      LargeDef.SplitKernels = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-symmetric") == 0)) {
      argc--; argv++;
      LargeDef.symmetric = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-margins") == 0)) {
      argc--; argv++;
      LargeDef.Margin = atoi(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-WghtVeloField") == 0)) {
      argc--; argv++;
      LargeDef.WghtVelField = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-RefMaxGrad") == 0)) {
      argc--; argv++;
      LargeDef.RefMaxGrad = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-epsilon") == 0)) {
      argc--; argv++;
      LargeDef.epsilon = atof(argv[1]);
      argc--; argv++;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-TypicInvAmp") == 0)) {
      argc--; argv++;
      LargeDef.sigmaX1 = atof(argv[1]);
      LargeDef.sigmaY1 = atof(argv[1]);
      LargeDef.sigmaZ1 = atof(argv[1]);
      argc--; argv++;
      LargeDef.NbKernels=1;
      LargeDef.MeasureTypicAmp=1;
      ok = true;
    }
    //5) Special outputs
    if ((ok == false) && (strcmp(argv[1], "-AOD") == 0)) {
      argc--; argv++;
      LargeDef.FlowLength = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-DetJacobian") == 0)) {
      argc--; argv++;
      LargeDef.DetJacobian = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-FinalDefVec") == 0)) {
      argc--; argv++;
      LargeDef.FinalDefVec = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-InitMomentum") == 0)) {
      argc--; argv++;
      LargeDef.CptInitMomentum = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-FinalDefInvVec") == 0)) {
      argc--; argv++;
      LargeDef.FinalDefInvVec = 1;
      ok = true;
    }
    if ((ok == false) && (strcmp(argv[1], "-ShowSSD") == 0)) {
      argc--; argv++;
      LargeDef.ShowSSD = 1;
      ok = true;
    }
    if (ok == false) usage();
  }
  
  //run process
  cout << "Large Deformation registration using Beg 05's technique ... \n"; cout.flush();
  LargeDef.Run();
  cout << "done" << endl;
  
  return 0;
}
