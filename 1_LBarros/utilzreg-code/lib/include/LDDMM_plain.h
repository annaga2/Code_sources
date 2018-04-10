/*=========================================================================
 
 Authors: Laurent Risser, Francois-Xavier Vialard
 
 =========================================================================*/

#ifndef _UTILZREGLDDMMPLAIN_H
#define _UTILZREGLDDMMPLAIN_H

#include <SciCalcPack.h>


/**
 * Class for Large Deformation registration using Beg 05's technique
 */

class LDDMM_Plain{
private:
  
protected:
  /// Function to launch the default calculations
  virtual void Run_Default();
  
  /// functions to perform the registration
  virtual void AllocateAllVariables();
  virtual void ReadAndTreatInputImages();
  virtual void ComputeEnergyGradient(int);
  virtual float UpdateVelocityField(int);
  virtual void SaveResultGradientDescent();
  
  ///functions to save and load various structures
  virtual void LoadVelocityFields(char *);
  virtual void SaveVelocityFields(VectorField *, char *);
  virtual void SaveDeformations(char *);
  virtual void SaveFlowLength(VectorField *,VectorField *,char *,char *);
  virtual void SaveEvoFlowLength(VectorField *,VectorField *,char *,char *);
  virtual void SaveVecDeformation(char *);
  virtual void SaveInvVecDeformation(char *);
  virtual void SaveInitMomentum(char *);
  virtual void SaveGlobalFlowLength(char *);
  virtual void SaveDetJacobian(char *);
  
  /// Protected parameters
  //scalar and vector fields
  ScalarField ImTemplate;   //3D image
  ScalarField ImTarget;     //3D image
  ScalarField J0;           //projection of 'ImTemplate' at a time between 0 and 1
  ScalarField J1;           //projection of 'ImTarget' at a time between 0 and 1
  ScalarField InitialMomentum;
  ScalarField GradInitialMomentum;
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
  LightFFTconvolver3D FFTconvolver;
  
public:
  
  /// Constructor
  LDDMM_Plain();
  
  /// Destructor
  ~LDDMM_Plain();
  
  /// Run  Large Deformation registration
  virtual void Run();
  
  
  /// public Parameters
  int iteration_nb;              //number of iterations 
  float epsilon;                 //Threshold on the energy gradient convergence
  int NbTimeSubdiv;              //Number of subdivision of virtual time steps
  float MaxVelocityUpdate;       //Maximum allowed velocity update at each subdivision and iteration (in voxels)
  float RefMaxGrad;              //Reference update gradient (typically the value of MaxGrad at the 1st iteration)
  float DeltaTimeSubdiv;         //time step between two subdivision
  int Margin;                    //Margin of the image in voxels where the calculations are reduced
  float sigmaX1,sigmaY1,sigmaZ1; //std. dev. of the 1st Gaussian kernel in direction x,y,z
  float WghtVelField;            //Weight of the velocity field in the energy
  int FlowLength;                //Save Length of the deformation flow from each voxel if .==1
  int DetJacobian;               //Save Determinant of the Jacobian at each voxel if .==1
  int FinalDefInvVec;            //Save Vector field of the estimated deformation from [Target] to [Source] if .==1
  int CptInitMomentum;              //Save the initial momentum if .==1
  int ShowSSD;                   //Show the evolution of the Sum of Squared Differences at t=1
  char PrefixInputs[256];        //Prefix of the files containing an initial velocity field
  char PrefixOutputs[256];       //Prefix of the files containing the final velocity field
  char SourceFile[256];          //name of the file containing the source image
  char TargetFile[256];          //name of the file containing the target image
  char MaskFile[256];            //name of the file containing the mask on the images
  float World_Target2Template[4][4]; // in the world domain, affine registration that matches the target image to the template (source) image
  float x_mm,y_mm,z_mm;   //voxel size in mm
};


#endif
