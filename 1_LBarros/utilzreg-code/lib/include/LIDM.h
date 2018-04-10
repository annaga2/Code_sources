/*=========================================================================
 
 Authors: Laurent Risser, Francois-Xavier Vialard
 
 =========================================================================*/

#ifndef _UTILZREGLIDM_H
#define _UTILZREGLIDM_H

#include <SciCalcPack.h>


/**
 * Class for Large Deformation registration using Beg 05's technique
 */

class LIDM{
private:
  
protected:
  /// Function to launch the default calculations
  virtual void Run_Default();
  
  /// functions to perform the registration
  virtual void AllocateAllVariables();
  virtual void ReadAndTreatInputImages();
  virtual void ComputeEnergyGradient(int);
  virtual void UpdateVelocityField_part1();
  virtual float UpdateVelocityField_part2(int);
  virtual void ProjectPOU(int);
  virtual void SaveResultGradientDescent();
  virtual float EstimRefWeight(float sigmaloc,LightFFTconvolver3D * ConvolverLoc,int NeedForInit=1);
  
  ///functions to save and load various structures
  virtual void LoadVelocityFields(char *);
  virtual void SaveVelocityFields(VectorField *, char *);
  virtual void SaveFinalDeformation(char *);
  virtual void SaveVecDeformation(char *);
  virtual void SaveInvVecDeformation(char *);
  
  /// Protected parameters
  //scalar and vector fields
  ScalarField ImTemplate;   //3D image
  ScalarField ImTarget;     //3D image
  ScalarField J0;           //projection of 'ImTemplate' at a time between 0 and 1
  ScalarField J1;           //projection of 'ImTarget' at a time between 0 and 1
  ScalarField PartiOfUnity;
  ScalarField ActualPOU;
  ScalarField ProjectionActualPOU;
  VectorField GradJ;
  ScalarField DetJacobians;
  VectorField ForwardMapping;
  VectorField BackwardMapping;
  VectorField VelocityField;
  VectorField GradE;
  

  //Quaternion to convert target coordinates into template (=source) coordinates
  float Template2TargetCoord[4][4];

  
  //size of the fields
  int NX;
  int NY;
  int NZ;
  int NT;
  
  //fft convolver (to smooth the images)
  MultiRegionFFTConvolver2 LIDMConvolver;
  
public:
  
  /// Constructor
  LIDM();
  
  /// Destructor
  ~LIDM();
  
  /// Run  Large Deformation registration
  virtual void Run();
  
  
  /// public Parameters
  int iteration_nb;              //number of iterations 
  int NbTimeSubdiv;              //Number of subdivision of virtual time steps
  float MaxVelocityUpdate;       //Maximum allowed velocity update at each subdivision and iteration (in voxels)
  float RefMaxGrad;              //Reference update gradient (typically the value of MaxGrad at the 1st iteration)
  float DeltaTimeSubdiv;         //time step between two subdivision
  int Margin;                    //Margin of the image in voxels where the calculations are reduced
  float ** stdDev;               //std. dev. of the Gaussian kernels
  float ** weight;               //weights of the Gaussian kernels
  float VFpenalizer;             //Factor to multiply the velocity field at each iteration. - related to the weight between the 2 energy terms
  int DetJacobian;               //Save Determinant of the Jacobian at each voxel if .==1
  int FinalDefInvVec;            //Save Vector field of the estimated deformation from [Target] to [Source] if .==1
  int ShowSSD;                   //Show the evolution of the Sum of Squared Differences at t=1
  int PreserveWeights;           //Force the weights to be considered as they are and not as apparent weights
  int MovingPOU;                 //Make move the partition of unity with the diffeomorphism (registration algorithm is not LIDM anymore)
  char PrefixInputs[256];        //Prefix of the files containing an initial velocity field
  char PrefixOutputs[256];       //Prefix of the files containing the final velocity field
  char SourceFile[256];          //name of the file containing the source image
  char TargetFile[256];          //name of the file containing the target image
  char PartiOfUnityFile[256];        //partition of unity to distinguish the influence of the two kernels in LIDM
  float World_Target2Template[4][4]; // in the world domain, affine registration that matches the target image to the template (source) image
  float x_mm,y_mm,z_mm;              //voxel size in mm
  float UnderSampleFactor;        //undersampling factor
  float SigmaPOU;                 //std of the gaussian kernel to smooth the partition of unity
};


#endif
