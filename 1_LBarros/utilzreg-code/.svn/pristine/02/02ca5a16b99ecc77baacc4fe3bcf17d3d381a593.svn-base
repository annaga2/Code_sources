/*=========================================================================
 
 Authors: Francois-Xavier Vialard, Laurent Risser
 
 Disclaimer: This software has been developed for research purposes only, and hence should 
 not be used as a diagnostic tool. In no event shall the authors or distributors
 be liable to any direct, indirect, special, incidental, or consequential 
 damages arising of the use of this software, its documentation, or any 
 derivatives thereof, even if the authors have been advised of the possibility 
 of such damage. 
 
 =========================================================================*/

#ifndef _IRTKLARGEDEFORMATIONSHOOTING_H
#define _IRTKLARGEDEFORMATIONSHOOTING_H

#include <SciCalcPack.h>

class EulerianShooting
{
  private:
  VectorField InvDiffeo, Diffeo, NablaI, NablaAdM, VelocityField,TempInvDiffeoLocal;
  VectorField AdjointVectorField,TempDiffeo,TempInvDiffeo,TempDiffeoLocal,TempVectorField;
  ScalarField ImTemplate, ImTarget, InitialMomentum, Image, Momentum, AdjointMomentum, AdjointImage,InputInitialMomentum;
  ScalarField TransfoImTarget;
  ScalarField GradientMomentum, TempAdMomentum, TempAdImage, OptimizedMomentum,TempScalarField,TempScalarField3;
  ScalarField Mask;
  int NX;
  int NY;
  int NZ;
  int NT;
  float DeltaX;
  // time step of one iteration
  float DeltaTimeSubdiv;
  
  static inline float Limiter(float a, float b){    //MinModLimiter
    if (a*b>0.0){if (a<b){return a;} else {return b;}}
    return 0.0;
  }

  
  
  
  // Cost to optimize
  float Cost;
  float Energy;
  float MaxVectorField;
  public:

  // Constructor and Destructor
  EulerianShooting();
  ~EulerianShooting();


  void AllocateVariablesShooting();
  void InitializeAdjointVariables();
  void InitializeVariables();
  void ReadAndTreatInputImages();
  void ComputeVelocityField(ScalarField *Momentum,VectorField * NablaI);
  void ComputeAdjointVectorField();
  void SchemeStep();
  void Shooting();
  void ShootingShow();
  void GradientDescent(float);
  void CptGradient();
  void ComputeDiffeoToSave();
  void SaveResult();
  void SaveDispField();
  void SaveTimeDepDispField();
  void SaveInitialMomentum();
  void SaveVelocityField();
  void SaveEnergyAndSSD();
  void Save_Diff_Tpl_DefTrg();
  void ReInitiateConvolver_HomoAppaWeights();

  float SimilarityMeasure();
  void Run();
  /// private parameters
  
  //fft convolver (to smooth the images)
  FFTconvolver3D FFTconvolver;
  
  /// public Parameters
  char SourceImageName[256];
  char TargetImageName[256];
  char InputInitialMomentumName[256];
  char InputMask[256];
  // Number of times from time 0 to time T
  int NbTimeSubdivisions;
  int NbIterations;

  // weight of the norm in the cost function
  float alpha;
  
  double weight1,sigmaX1,sigmaY1,sigmaZ1;
  double weight2,sigmaX2,sigmaY2,sigmaZ2;
  double weight3,sigmaX3,sigmaY3,sigmaZ3;
  double weight4,sigmaX4,sigmaY4,sigmaZ4;
  double weight5,sigmaX5,sigmaY5,sigmaZ5;
  double weight6,sigmaX6,sigmaY6,sigmaZ6;
  double weight7,sigmaX7,sigmaY7,sigmaZ7;
  int GreyLevAlign, Margin;
  float GLA_Padding_Src;         //if grey level alignment: padding value for the source image
  float GLA_Padding_Trg;         //if grey level alignment: padding value for the target image
  int indicatorInitialMomentum;
  int OutIniMoTxt;
  int OutVeloField;
  int OutDistEnSim;
  int OutDeformation;
  int OutDiff_Tpl_DefTrg;
  int OutDispField;
  int OutDispFieldEvo;
  int OutFinalDef;
  // size of the update in the GradientDescent
  float MaxUpdate;
  int indicatorMask;
  char PrefixOutputs[256];       //Prefix of the output files
  char PrefixInputs[256];       //Prefix of the input files
  float World_Target2Template[4][4]; // in the world domain, affine registration that matches the target image to the template (source) image
  float Template2TargetCoord[4][4];  //Quaternion to convert target coordinates into template (=source) coordinates
  float x_mm,y_mm,z_mm;   //voxel size in mm
  float UnderSampleTplFactor;
};

#endif
