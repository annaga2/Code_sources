/*=========================================================================
 
 Date      : $Date: 18.01.2011$
 Changes   : $Authors: Laurent Risser, Francois-Xavier Vialard$
 
 =========================================================================*/

#ifndef _IRTKLARGEDEFORMATIONGRADIENTLAGRANGE_H
#define _IRTKLARGEDEFORMATIONGRADIENTLAGRANGE_H

#include <SciCalcPack.h>


/**
 * Class for basic Demons registration
 */

class LargeDefDemons_plain{
private:
	
protected:
	
	
	/// functions to perform the registration
	virtual void AllocateAllVariables();
	virtual void ReadAndTreatInputImages();
	virtual void SaveResultGradientDescent();
	virtual void ComputeUpdateFieldSSD();
	virtual void ComputeUpdateFieldMI();
	virtual void ReadAndConvertMapping();
	virtual void ControlMaxUpdate(int IterationNb);
  
	///functions to save and load various structures
	virtual void LoadVelocityField(char *);
	virtual void SaveInvTotalDisplacement(VectorField *, char *);
	virtual void SaveVelocityField(VectorField *, char *);
	virtual void SaveDeformations(char *);
  virtual float EstimRefUpdateScale(float,LightFFTconvolver3D *);
	
	/// Protected parameters
	//scalar and vector fields
	ScalarField MovingIm;       //3D image
	ScalarField FixedIm;         //3D image
	ScalarField DeformedMovingIm;
	VectorField VelocityField;
	VectorField GradE;
  VectorField IniDispField;
  
  
  //FFT convolver
  LightFFTconvolver3D FFTconvolver_fluid;
  LightFFTconvolver3D FFTconvolver_Diff;
  
  //parameters related to the mutual information
  MImanager NorMutInfMan;
	
  //Quaternion to convert FixedIm coordinates into MovingIm (=source) coordinates
  float FixedIm2MovingImCoord[4][4];
  
	//size of the fields
	int NX;
	int NY;
	int NZ;
	int NT;

  //size of the moving image
  int NXs;
	int NYs;
	int NZs;
  
  
	
public:
	
	/// Constructor
	LargeDefDemons_plain();
	
	/// Destructor
	~LargeDefDemons_plain();
	
	/// Run  Large Deformation registration
	virtual void Run();
	
	
	/// public Parameters
	int iteration_nb;              //number of iterations 
	float MaxVelocityUpdate;       //Maximum allowed velocity update at each subdivision and iteration (in voxels)
	float RefMaxGrad;              //Reference update gradient (typically the value of MaxGrad at the 1st iteration)
  int IndicatorMI;
  float sigmaFluid;
  float sigmaDiff;
  char PrefixInputs[256];        //Prefix of the files containing an initial velocity field
	char PrefixOutputs[256];       //Prefix of the files containing the final velocity field
	char MovingImFile[256];         //name of the file containing the moving image 
	char FixedImFile[256];         //name of the file containing the fixed image 
  float lambdaX;            //controls the matching uncertainty
	float InitMaxUpdate;
  float MaxUpdateAllowed;
  float World_FixedIm2MovingIm[4][4]; // in the world domain, affine registration that matches the fixed image to the moving image
  float x_mm,y_mm,z_mm;   //voxel size in mm
  float UnderSampleFixedImFactor;
  int IniDispFieldDefined;
  char IniDispFieldX[256];
  char IniDispFieldY[256];
  char IniDispFieldZ[256];
  
};


#endif
