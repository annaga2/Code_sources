

/*=========================================================================
 
 
 Author: Laurent Risser, Francois-Xavier Vialard
 
 Disclaimer: This software has been developed for research purposes only, and hence should 
 not be used as a diagnostic tool. In no event shall the authors or distributors
 be liable to any direct, indirect, special, incidental, or consequential 
 damages arising of the use of this software, its documentation, or any 
 derivatives thereof, even if the authors have been advised of the possibility 
 of such damage. 
 
 
 =========================================================================*/


#ifndef _LD_SCICALCPACK_H
#define _LD_SCICALCPACK_H

//Parallelization with openmp
#ifndef __APPLE__
  #define COMPILE_WITH_OPENMP
  #include <omp.h>
#endif

//related to the i/o nifti and added instead of 'irtkImage.h' in irtk
#define MIN_HEADER_SIZE 348
#define NII_HEADER_SIZE 352


// C++ header files
#include <iostream>
using namespace std;
#include <iomanip>
#include <fstream>
#include <complex>
#include <algorithm>
#include <string>
#include <limits>

// C header files
#include <stdio.h>
#include <stdlib.h>
#include <limits.h>
#include <float.h>
#include <math.h>
#include <string.h>


#include "nifti1.h"
#include "nifti1_io.h"


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           1:   FUNCTIONS FOR THE CLASS "ScalarField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class ScalarField{
  /// ******************************************************************************
private:
  //size of the fields
  int NXtY;      //NX*NY
  int NXtYtZ;    //NX*NY*NZ
  
  //scalar field
  float * ScalField;
  
  /// ******************************************************************************
public:
  //size of the fields
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  int NT;        //image size on time
  
  float Image2World[4][4];  //quaternion to go from the image coordinates to the world coordinates
  float World2Image[4][4];  //quaternion to go from the world coordinates to the image coordinates
  
  /// Constructor and destructor
  ScalarField();
  ~ScalarField();
  
  /// functions associated to the field
  //put a value in the scalar field
  virtual inline void P(float value,int x,int y,int z=0,int t=0){
    this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x]=value;
  }
  
  //Add a value to the scalar field
  virtual inline void Add(float value,int x,int y, int z=0,int t=0){
    this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x]+=value;
  }
  
  // put a the same value at every points of the scalar field
  virtual void PutToAllVoxels(float cste,int t=0);
  
  //get a value from the scalar field
  virtual inline float G(int x,int y,int z=0,int t=0){
    return this->ScalField[t*this->NXtYtZ+z*this->NXtY+y*this->NX+x];
  }
  
  //get a value from the scalar field by linear interpolation
  virtual float G(double x,double y,double z=0.,int t=0);
  virtual float G(float x, float y, float z=0., int t=0);

  //same as G but zero is returned if the coordinate is out of the image
  virtual float G_NoExtrapo(double x,double y,double z=0.,int t=0);
  virtual float G_NoExtrapo(float x, float y, float z=0., int t=0);

  
  //get a value from the scalar field by linear interpolation. 
  //-> Here the coordinate not defined in the image space but in 'coordSpace' 
  //-> The 4*4 matrix 'coordSpace2imageSpace' is a quaternion that allows to transform the coordinates from 'coordSpace' to the image space
  virtual float G(float coordSpace2imageSpace[4][4],double x,double y,double z=0.,int t=0,int NN=0);
  virtual float G(float coordSpace2imageSpace[4][4],float x,float y,float z=0.,int t=0,int NN=0);
  virtual float G(float coordSpace2imageSpace[4][4],int x,int y,int z=0.,int t=0,int NN=0);
  
  // get the maximum absolute values out of the scalar field
  virtual float GetMaxAbsVal(int t=0);
  
  //read a scalar field (in a nifti image)
  virtual void Read(char *);
  
  //read a scalar field in a ROI (in a nifti image)
  // -> advanced memory managment for large images: only allocate memory for the ROI
  // -> Inputs are different than in Read_ROI_Given_ImageToWorld:
  //     We give here the min and max {X,Y,Z,T} in the input image defining the outputed ROI
  virtual void Read_only_ROI(char *,int,int,int,int,int,int,int,int);

  //read a scalar field in a ROI (in a nifti image)
  // -> advanced memory managment for large images: only allocate memory for the ROI
  // -> Inputs are different than in Read_only_ROI:
  //     We give here the ImageToWorld coordinates and the size (in voxels) of the outputed ROI
  // -> remark that nearest ngbh intensity is considered to sample the output image
  virtual void Read_ROI_Given_ImageToWorld_and_Size(float ROI_Image2World[4][4],int NBX,int NBY,int NBZ,char * RefImageName);
  
  //read a scalar field in a ROI (in a nifti image)
  // -> advanced memory managment for large images: undersample the image direcly by averaging blocks of voxels of size BlockS*BlockS*BlockS
  //Use 'Read_and_Undersample' or 'Read_and_Interpolate' to have finer resamplings requiring more memory
  virtual void Read_directly_Undersampled(char * ImageName,int BlockS);
  
  //read a scalar field and perform linear interpolation to give it a specific size
  virtual void Read_and_Interpolate(char *,int,int,int);
  
  ///read a scalar field and perform undersampling by a factor 'factor'
  virtual void Read_and_Undersample(char * ImageName,float factor,float NN=0);

  ///read a scalar field and expend its domain
  virtual void ReadAndExpend(char * ,int ,int ,int ,int ,int ,int);
    
  //create a void scalar field. All values are initialized to 'cste' which is null by default. No message is printed if Verbose==1.
  virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1,float cste=0.0,int Verbose=1);
  
  //Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
  virtual void SlashFieldSize(int verbative=1);

  //write a scalar field in a nifti image
  //The optional 2nd file is an input file containing the headers 
  virtual void Write(char *);
  virtual void Write(char *,char *);
    
  //return the number of voxels in a ScalarField
  virtual inline int GetNbVoxels(void){
    return this->NX *this->NY*this->NZ*this->NT;
  }
  
  ///grey levels alignment of the grey levels using optimal transport 
  ///-> optimal transportation minmizes the wassertein 1 distance between the linearly aligned histograms
  virtual void GreyLevAlignment(ScalarField * RefImage);


  ///Compute the histogram. The histogram is normalized (sum of values equal to 1). 
  /// -> Input_BinsNb is the number of bins in the histogram
  /// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
  void CptHistogram(int Input_BinsNb,float * Output_Histo_x_axis,float * Output_Histo_y_axis,int useLogHisto=0);

  ///Compute the cumulative histogram. The cumulative histogram is normalized (last value equal to 1). 
  /// -> Input_BinsNb is the number of bins in the histogram
  /// -> Output_Histo_x_axis and Output_Histo_y_axis represent the histogram and must be allocated before calling the function
  void CptCumulativeHistogram(int Input_BinsNb,float * Output_CumHisto_x_axis,float * Output_CumHisto_y_axis,int useLogHisto=0);
};

///get the 'Image 2 World matrix' of an image without loading the image  (only its header)
void Get_Image2World(char *,float LocI2W[4][4]);



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           2:   FUNCTIONS FOR THE CLASS "VectorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class VectorField{
  /// ******************************************************************************
private:
  //size of the fields
  int NXtY;      //NX*NY
  int NXtYtZ;    //NX*NY*NZ
  int NXtYtZtT;  //NX*NY*NZ*NT
  
  //scalar field
  float * VecField;
  
  /// ******************************************************************************
public:
  //size of the fields
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  int NT;        //image size on time
  
  float Image2World[4][4];  //quaternion to go from the image coordinates to the world coordinates
  float World2Image[4][4];  //quaternion to go from the world coordinates to the image coordinates

  
  /// Constructor
  VectorField();
  
  /// Destructor
  ~VectorField();
  
  //put a value in the scalar field
  virtual inline void P(float value,int IdDirec,int x,int y,int z=0,int t=0){
    this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x]=value;
  }
  
  //Add a value to the scalar field
  virtual inline void Add(float value,int IdDirec,int x,int y, int z=0,int t=0){
    this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x]+=value;
  }
  
  //put the same value at all entries of the vector field
  virtual void PutToAllVoxels(float cste,int t=0);
  
  //get a value from the scalar field
  virtual inline float G(int IdDirec,int x,int y,int z=0,int t=0){
    return this->VecField[IdDirec*this->NXtYtZtT + t*NXtYtZ + z*this->NXtY + y*this->NX + x];
  }
  
  //get a value from the scalar field by linear interpolation
  virtual float G(int IdDirec,double x,double y,double z=0.,int t=0);
  virtual float G(int IdDirec,float x,float y,float z=0.,int t=0);
  
  
  // get the maximum of the absolute values of the vector field
  virtual float GetMaxAbsVal(int t=0);
  
  //read a vector field (from 3 nifti images -> X, Y, Z)
  virtual void Read(char *,char *,char *);
  
  //read a vector field and perform linear interpolation to give it a specific size
  //  If rescaleVF!=0, the values of the vector field are rescaled proportionally to
  //  the re-sizing (usefull for velocity fields)
  virtual void Read_and_Interpolate(char *,char *,char *,int,int,int,int rescaleVF=0);
  
  //create a void vector field. No message is printed if  Verbose!=1
  virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1,int Verbose=1);

  //Do not destruct 'this' but strongly reduce its size. As a result, it cannot be used any more until 'CreateVoidField' realoc all the memory.
  virtual void SlashFieldSize(int verbative=1);

  //write a vector field (in 3 nifti images -> X, Y, Z)
  //The optional 4th file is an input file containing the headers 
  virtual void Write(char *,char *,char *);
  virtual void Write(char *,char *,char *,char *);
};


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                           3:   FUNCTIONS FOR THE CLASS "TensorField"
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


class TensorField{
  /// ******************************************************************************
private:
  //size of the fields
  //int NXtY;      //NX*NY
  //int NXtYtZ;    //NX*NY*NZ
  //int NXtYtZtT;  //NX*NY*NZ*NT
  //int NXtYtZtTt3;  //NX*NY*NZ*NT*3
  
  int NXt9;         //NX*9
  int NXtYt9;       //NX*NY*9
  int NXtYtZt9;     //NX*NY*NZ*9  
  
  
  //scalar field
  float * TField;
  
  /// ******************************************************************************
public:
  //size of the fields
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  int NT;        //image size on time
  
  /// Constructor
  TensorField();
  
  /// Destructor
  ~TensorField();
  
  //create a void tensor field
  virtual void CreateVoidField(int NBX,int NBY,int NBZ=1,int NBT=1);
  
  //put a value in the tensor field
  virtual void P(float value,int IdDirec1,int IdDirec2,int x,int y,int z=0,int t=0);
  
  //add a value in the tensor field
  virtual void Add(float value,int IdDirec1,int IdDirec2,int x,int y,int z=0,int t=0);
  
  //get a value from the tensor field
  virtual float G(int IdDirec1,int IdDirec2,int x,int y,int z=0,int t=0);
  
  //Add a tensorised vector to the existing tensor
  virtual void AddTensorisedVector(float vec[3],int x,int y,int z=0,int t=0);
  
  //Perform a principal component analysis of the 3*3 tensor
  //The outputs are:
  // lambda1,lambda2,lambda3: the eigenvalues in decrasing order
  // vec1: 1st eigenvector  (must be initialised as vec1[3])
  // vec2: 2nd eigenvector  (must be initialised as vec2[3])
  // vec3: 3rd eigenvector  (must be initialised as vec3[3])
  virtual void PCA(float vec1[3],float vec2[3],float vec3[3],float * lambda1, float * lambda2, float * lambda3, int x,int y,int z=0,int t=0);
  
};




///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           4: CLASSES TO PERFORM CONVOLUTION AND DECONVOLUTION USING FFT
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

///4.1: Main class (sufficiently general to have kernels which are not separable but memory consuming)
class FFTconvolver3D{
  /// ******************************************************************************
private:
  //size of the inputs (that will be copied in images having sizes = 2^{...})
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  
  //size of the fields transformed by fft
  int NXfft;
  int NYfft;
  int NZfft;
  
  //fields transformed by fft
  ScalarField RealSignalForFFT;
  ScalarField ImagSignalForFFT;
  ScalarField RealFilterForFFT;
  ScalarField ImagFilterForFFT;
  
  //temporary scalar field
  ScalarField ImageTemp;
  
  //design a kernel that is the sum of up to 7 Gaussians
  //... if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
  void MakeSumOf7AnisotropicGaussianFilters(float weight1=100,float sigmaX1=1,float sigmaY1=1,float sigmaZ1=1,
                                              float weight2=0,float sigmaX2=1,float sigmaY2=1,float sigmaZ2=1,
                                              float weight3=0,float sigmaX3=1,float sigmaY3=1,float sigmaZ3=1,
                                              float weight4=0,float sigmaX4=1,float sigmaY4=1,float sigmaZ4=1,
                                              float weight5=0,float sigmaX5=1,float sigmaY5=1,float sigmaZ5=1,
                                              float weight6=0,float sigmaX6=1,float sigmaY6=1,float sigmaZ6=1,
                                              float weight7=0,float sigmaX7=1,float sigmaY7=1,float sigmaZ7=1,
                                              int NormalizeWeights=1);
  
  //Fast Fourier Transform
  void DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal);
  
  //Inverse Fast Fourier Transform
  void InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal);
  
  //Fast Fourier Transform of numerical recipies (slighly modified)
  void four1NR(float data[], unsigned long nn, int isign);
  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  FFTconvolver3D();
  
  //Destructor
  ~FFTconvolver3D();
  
  //Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
  //7 Gaussians (set some weights to 0 if less Gaussians are required)
  //* NX, NY, NZ: is the size of the input image
  //* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //* w2,sX2,sY2,sZ2,: weight of the 2nd Gaussian kernel and std. dev. in direction X, Y, Z
  //* w3,sX3,sY3,sZ3,: weight of the 3rd Gaussian kernel and std. dev. in direction X, Y, Z
  //* w4,sX4,sY4,sZ4,: weight of the 4th Gaussian kernel and std. dev. in direction X, Y, Z
  //* w5,sX5,sY5,sZ5,: weight of the 5th Gaussian kernel and std. dev. in direction X, Y, Z
  //* w6,sX6,sY6,sZ6,: weight of the 6th Gaussian kernel and std. dev. in direction X, Y, Z
  //* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolver(int NBX,int NBY, int NBZ, 
                                 float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., 
                                 float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., 
                                 float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., 
                                 float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., 
                                 float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., 
                                 float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., 
                                 float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
  
  //change the kernel of the convolver (same notations as the constructor)
  //... here the new kernel is normalized (w1+...+w7 is normalized to 1)
  virtual void ChangeKernel(float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., 
                            float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., 
                            float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., 
                            float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1.,
                            float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1.,
                            float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1.,
                            float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
  
  //change the kernel of the convolver (same notations as the constructor)
  //... here the new kernel is not normalized
  virtual void ChangeKernel_SingleScale(float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1.);
  
  
  //convolution of a 3D scalar field using the predifined kernel
  virtual void Convolution(ScalarField *);
  
  //convolution of a 3D vector field using the predifined kernel
  virtual void Convolution(VectorField *);
  
  //convolution of the real scalar field defined inside of the class using the predifined kernel
  virtual void Convolution();
  
  //deconvolution of a 3D scalar field using the predifined kernel
  // !!! NOT VALIDATED !!!
  virtual void Deconvolution(ScalarField *);
  
  //put a value in the real part of the field that is transformed by the class
  virtual void P(float value,int x,int y, int z=0);
  
  //put a value in the real part of the field that is transformed by the class
  virtual float G(int x,int y, int z=0);
};


void SmoothVFUsingFFT(VectorField * SmoothedField,FFTconvolver3D * FFTconvolver_loc);



///4.3: light weight convolver (requires little memory but can only have separable kernels. Contains the most common functions of FFTconvolver3D, plus a direct smoothing of velocity fields)
class LightFFTconvolver3D{
  /// ******************************************************************************
private:
  //size of the inputs (that will be copied in images having sizes = 2^{...})
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  
  //size of the fields transformed by fft
  int NXfft;
  int NYfft;
  int NZfft;
  
  //fields transformed by fft
  ScalarField RealSignalForFFT_X;
  ScalarField ImagSignalForFFT_X;
  ScalarField RealSignalForFFT_Y;
  ScalarField ImagSignalForFFT_Y;
  ScalarField RealSignalForFFT_Z;
  ScalarField ImagSignalForFFT_Z;
  
  ScalarField RealFilterForFFT_X;
  ScalarField ImagFilterForFFT_X;
  ScalarField RealFilterForFFT_Y;
  ScalarField ImagFilterForFFT_Y;
  ScalarField RealFilterForFFT_Z;
  ScalarField ImagFilterForFFT_Z;
  
  //design a kernel that is the sum of up to 7 Gaussians
  //... if the option NormalizeWeights == 0 then the different weights (and then the whole filter) are not normalized
  void MakeSumOf7AnisotropicGaussianFilters(float weight1=100,float sigmaX1=1,float sigmaY1=1,float sigmaZ1=1,
                                              float weight2=0,float sigmaX2=1,float sigmaY2=1,float sigmaZ2=1,
                                              float weight3=0,float sigmaX3=1,float sigmaY3=1,float sigmaZ3=1,
                                              float weight4=0,float sigmaX4=1,float sigmaY4=1,float sigmaZ4=1,
                                              float weight5=0,float sigmaX5=1,float sigmaY5=1,float sigmaZ5=1,
                                              float weight6=0,float sigmaX6=1,float sigmaY6=1,float sigmaZ6=1,
                                              float weight7=0,float sigmaX7=1,float sigmaY7=1,float sigmaZ7=1,
                                              int NormalizeWeights=1);

  //Fast Fourier Transform
  // -> if axis == 0 -> FFT on X axis
  // -> if axis == 1 -> FFT on Y axis
  // -> if axis == 2 -> FFT on Z axis
  void DirectFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis);
  
  //Inverse Fast Fourier Transform
  // -> if axis == 0 -> IFFT on X axis
  // -> if axis == 1 -> IFFT on Y axis
  // -> if axis == 2 -> IFFT on Z axis
  void InverseFFT(ScalarField * RealSignal,ScalarField * ImaginarySignal,int axis);
  
  //Fast Fourier Transform of numerical recipies (slighly modified)
  void four1NR(float data[], unsigned long nn, int isign);
  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  LightFFTconvolver3D();
  
  //Destructor
  ~LightFFTconvolver3D();
  
  //Initiate the complex fields for the FFT and the smoothing kernels being the sum of up to 
  //7 Gaussians (set some weights to 0 if less Gaussians are required)
  //* NX, NY, NZ: is the size of the input image
  //* w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //...
  //* w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolver(int NBX,int NBY, int NBZ, 
                                 float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., 
                                 float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., 
                                 float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., 
                                 float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., 
                                 float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., 
                                 float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., 
                                 float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.,
                                 int NormalizeWeights=1);
  
  //change the kernel of the convolver (same notations as the constructor)
  //... here the new kernel is normalized (w1+...+w7 is normalized to 1)
  virtual void ChangeKernel(float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1.,
                            float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1.,
                            float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1.,
                            float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1.,
                            float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1.,
                            float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1.,
                            float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.,
                            int NormalizeWeights=1);

  //change the kernel of the convolver (same notations as the constructor)
  //... here the new kernel is not normalized
  virtual void ChangeKernel_SingleScale(float weight1,float sigmaX1,float sigmaY1,float sigmaZ1);

  //convolution of a 3D scalar field using the predifined kernel
  //If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
  virtual void Convolution(ScalarField *,int TimeFrame=-1);
  
  //convolution of a 3D vector field using the predifined kernel
  //If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
  virtual void Convolution(VectorField *,int TimeFrame=-1);

  //convolution of a 3D vector field using the predifined kernel. Convolution is performed in the ROI defined by (xmin, xmax, ymin, ymax, zmin, zmax) only.
  //If (TimeFrame==-1), all time frames are smoothed. / If (TimeFrame>=0), only time frame 'TimeFrame' is smoothed.
  virtual void ConvolutionInROI(VectorField *,int xmin,int xmax,int ymin,int ymax,int zmin,int zmax,int TimeFrame=-1);

  //Hack to perform convolution of the 3D vector field 'VF' in a masked region with mirror conditions.
  //  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
  //  -> sX1,sY1,sZ1: size of the smoothing kernel
  virtual void Convolution_Mask_Mirror(VectorField * VF,ScalarField * Mask, int MaskId=1);
 
  //Hack to perform convolution of the 3D scalar field 'SF' in a masked region with mirror conditions.
  //  -> Mask: convolution is performed where the mask equals 'MaskId' only. Mirror conditions are applied at the boundary of the domain.
  virtual void Convolution_Mask_Mirror(ScalarField * SF,ScalarField * Mask, int MaskId=1);
};



///4.4: class where we consider N regions to smooth

class MultiRegionFFTConvolver2{
  /// ******************************************************************************
private:
  LightFFTconvolver3D * Region_convolver;
  ScalarField PartitionOfUnity; //3D + channels
  int * xmin;
  int * xmax;
  int * ymin;
  int * ymax;
  int * zmin;
  int * zmax;
  int NumberOfRegions;
  VectorField TempVF1;
  VectorField TempVF2;

  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  MultiRegionFFTConvolver2();
  
  //Destructor
  ~MultiRegionFFTConvolver2();
  //Initiate the convolver in all regions using same kernel
  //-> 'Part_Of_Unity' is a 3D mask which define different ROIs. It only contains integer values each of them associated to a ROI. Partition of unity is first 
  //   defined by splitting this mask into several channels, each of them having an intensity equal to 1 in the corresponding ROI and 0 otherwise. Each channel
  //   is then smoothed with a Gaussian kernel of stddev 'sigmaPOI'.
  //-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
  //     * NX, NY, NZ: is the size of the input image
  //     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //     * ...
  //     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolver(ScalarField * Part_Of_Unity, float sigmaPOI=5., float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);

  //Initiate the convolver in all regions using same kernel
  //-> Part_Of_Unity is a '3D + channels' scalar field which encodes the partition of unity in the different channels.
  //   * Its size and number of channels (NBT actually) defines the size and the number of regions of the convolver 
  //   * The maximum point-wise sum of the probabilities may be different to 1: normalisation will be automatically performed
  //   * Point-wise sum of the probabilities may vary in space. If so, a background region will be automatically defined
  //-> 7 Gaussians (set some weights to 0 if less Gaussians are required)
  //     * NX, NY, NZ: is the size of the input image
  //     * w1,sX1,sY1,sZ1,: weight of the 1st Gaussian kernel and std. dev. in direction X, Y, Z
  //     * ...
  //     * w7,sX7,sY7,sZ7,: weight of the 7th Gaussian kernel and std. dev. in direction X, Y, Z
  virtual void InitiateConvolverWithActualPOI(ScalarField * Part_Of_Unity, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);

  //save the actual partition of unity (after potential treatments in InitiateConvolver or undersampling)
  //-> 1st char* is the name in which the image is saved
  //-> 2nd char* is the name of the image that will be used in the header of the saved image
  virtual void SaveActualParitionOfUnity(char *,char *);

  //change the smoothing kernel in one region
  virtual void ChangeKernelInOneRegion(int IdRegion, float w1=1.,float sX1=1.,float sY1=1.,float sZ1=1., float w2=0.,float sX2=1.,float sY2=1.,float sZ2=1., float w3=0.,float sX3=1.,float sY3=1.,float sZ3=1., float w4=0.,float sX4=1.,float sY4=1.,float sZ4=1., float w5=0.,float sX5=1.,float sY5=1.,float sZ5=1., float w6=0.,float sX6=1.,float sY6=1.,float sZ6=1., float w7=0.,float sX7=1.,float sY7=1.,float sZ7=1.);
  
  //Update the parition of unity 
  //-> it must have the same size and number of layers/times as in the current POU (tested)
  //-> to make sense, the new POU must have a sum of intensities across layers/times equals to 1 at each voxel (not tested)
  virtual void UpdatePartitionOfUnity(ScalarField * Part_Of_Unity);

  //convolution of a 3D vector field using the predifined kernel
  virtual void Convolution(VectorField * VF);
  
  //return the number of regions considered
  virtual int GetRegionNb();
};



  



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                  5: CLASS TO MANAGE THE MUTUAL INFORMATION
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

class MImanager{
  /// ******************************************************************************
private:
  //size of the treated images
  int NX;        //image size on X axis
  int NY;        //image size on Y axis
  int NZ;        //image size on Z axis
  
  //nbBins
  int NumberOfBinsS;
  int NumberOfBinsT;
  float MinGreyLevelsS;
  float MinGreyLevelsT;
  float SizeStepsGreyLevelsS; 
  float SizeStepsGreyLevelsT; 
  
  //margin on which the mutual information is not computed, plus mask indicator
  int Margin;
  int indicatorMaskDefined;
  
  //registered images and mask (ROI is where the mask intensities are 1)
  ScalarField * SrcImage;
  ScalarField * TrgImage;
  ScalarField * Mask;
  
  //joint Histo
  float ** JointHistogram;
  float * MarginalHistogramS;
  float * MarginalHistogramT;
  
  //joint and marginal entropies
  float JointEntropy;
  float MarginalEntropyS;
  float MarginalEntropyT;
  
  //normalised mutual informations
  float MI;
  
  //indicator saying that the saved computations are made on up-to-date data
  int indicatorUpdatedSrcHisto;
  int indicatorUpdatedTrgHisto;
  int indicatorUpdatedJointHisto;
  int indicatorUpdatedMI;
  
  
  //normalized the intensities adapted to the bins of the histograms
  float GiveFloatBinT(float intensity);
  float GiveFloatBinS(float intensity);
  
  
  //when computing the histograms give the contribution of 'intensity' in the bin 'IdBin' for the source and target image
  float GiveValueParzenWindowS(float intensity,int IdBin);
  float GiveValueParzenWindowT(float intensity,int IdBin);
  
  //compute the histograms
  void ComputeJointHistogramAndEntropy();
  void ComputeMarginalHistogramAndEntropyS();
  void ComputeMarginalHistogramAndEntropyT();
  
  
  
  /// +++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
public:
  //Constructor
  MImanager();
  
  //Destructor
  ~MImanager();
  
  //Initiate the MI manager without any mask
  //Important: the current MImanager will always point to the source and target images defined here (with intensities that can change!!!)
  virtual void Initiate(ScalarField * SourceImage,ScalarField * TargetImage,int NbBinsSrc=30,int NbBinsTrg=30, int LocMargin=3);

  //Initiate the MI manager with a mask. The MI and MI gradients will only be computed where the mask equals 1
  virtual void Initiate(ScalarField * SourceImage,ScalarField * TargetImage,ScalarField * ROI_Mask,int NbBinsSrc=30,int NbBinsTrg=30, int LocMargin=3);

  //returns the normalized mutual information (plus update all the histograms)
  virtual float EvaluateMI();
  
  //returns the estimated gradient of normalized mutual information
  virtual void EvaluateGradMI(VectorField * Gradient);
  
  //Indicate to the MI manager that the Source image has changed
  virtual void IndicateSrcHasChanged();
  
  //Indicate to the MI manager that the Target image has changed
  virtual void IndicateTrgHasChanged();
};





///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///           6: LOW LEVEL FUNCTIONS MAKING USE OF THE CLASSES ScalarField AND VectorField 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


///+++++++++++++++++++++++++++++++++      6.1: Diffusion stuffs       +++++++++++++++++++++++++++++++++

///Isotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
void Diffusion_3D(ScalarField * SField, float alpha, float dTau,int ITERATIONS_NB, float dx=1,float dy=1, float dz=1);
void Diffusion_3D(VectorField * VField, float alpha, float dTau,int ITERATIONS_NB, float dx=1,float dy=1, float dz=1);

///Isotropic diffusion of 'SField' or 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///Neumann conditions are considered at the domain boundaries 
///If optNoMaskNoDef==1 then the values for which the mask is 0 are set to 0
void Diffusion_3D(ScalarField * SField,ScalarField * Mask, int MaskId, float alpha, float dTau,int ITERATIONS_NB, int optNoMaskNoDef=0, float dx=1,float dy=1, float dz=1);
void Diffusion_3D(VectorField * VField,ScalarField * Mask, int MaskId, float alpha, float dTau,int ITERATIONS_NB,  int optNoMaskNoDef=0,int direction=-1, float dx=1,float dy=1, float dz=1);


///Isotropic diffusion of 'VField' in the mask 'Mask' (values='MaskId') during 'ITERATIONS_NB' where 'dTau' is the time step.
///Diffusion intensity is defined in alpha.
///The size of each voxel is defined by dx, dy, dz.
///!!! -> Dirichlet conditions are considered at the domain boundaries !!!
void Diffusion_3Dbis(VectorField * VField,ScalarField * Mask,int MaskId, float alpha, float dTau,int ITERATIONS_NB, int optNoMaskNoDef=0, float dx=1,float dy=1, float dz=1);



///Anisotropic diffusion of 'SField' during 'ITERATIONS_NB' where 'dTau' is the time step.
///Directions and intensities of diffusion are defined using ax, ay, az.
///The size of each voxel is defined by dx, dy, dz.
///Remark: very similar to the classic paper of Perona-Malik
void anisoDiff_3D(ScalarField * SField,float ax, float ay, float az, float dTau,int ITERATIONS_NB, float dx=1,float dy=1, float dz=1);


///compute the distance map of the edges in a 3D Mask. Greedy algorithm updating the distances in a band.
///-> 'NeaBoun' contains a vector field pointing to the nearest boundaries
///-> 'Distance' contains the distance to the nearest boundary
void Cpt_NearestBoundary(ScalarField * Mask,VectorField * NeaBoun,ScalarField * Distance);


///Remove the normal contributions of a velocity field that are too close to a boundary
///-> 'NeaBoun' contains the vector map that points to the closest voxel at a boundary
///-> 'TempSF'  = 0 where no vector is defined in NeaBoun and = 1 otherwise
///-> BoundaMargin is the margin around the bounadries (in voxels) in which SmoothedField has reduced normal contirbutions
void RemoveNormalContributions(VectorField * SmoothedField,VectorField * NeaBoun,ScalarField * TempSF,float BoundaMargin, float dx,float dy, float dz,int SetBoundaryToZero=0);


///+++++++++++++++++++++++++++++++++      6.2: derivatives and diamonds      +++++++++++++++++++++++++++++++++

///6.2.1: derivatives

//Compute the gradient of the scalar field "SField" and put the result in "Gradient"
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_Grad_ScalarField(ScalarField * SField,VectorField * Gradient,int SpecificTimeFrame=-1,float DeltaX=1);
void Cpt_Grad_MaskedScalarField(ScalarField * SField,VectorField * Gradient,ScalarField * Mask,int SpecificTimeFrame=-1,float DeltaX=1);


//Compute (d VField(X) / d x) + (d VField(Y) / d y) + (d VField(Z) / d z) and put the result in 'GradScalVF'
//where 'VField' is a vector field and 'GradScalVF' a scalar field
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_Grad_Scal_VectorField(VectorField * VField,ScalarField * GradScalVF,int SpecificTimeFrame=-1,float DeltaX=1);

//Compute the determinant of the Jacobian of the vector field 'VField' and put the result in the scalar field 'DetJ'
//* 'DeltaX' is the spatial step between two voxels.
//* If SpecificTimeFrame<0, the calculations are done in all time frames. They are only done
//  in the time frame 'SpecificTimeFrame' otherwise.
void Cpt_JacobianDeterminant(VectorField * VField,ScalarField * DetJ,int SpecificTimeFrame=-1,float DeltaX=1);




///+++++++++++++++++++++++++++++++++      6.3: image deformation and mappings      +++++++++++++++++++++++++++++++++


///6.3.0: affine transformation

///Image transformation using the 4*4 -- quaternion like -- matrix
void Project3DImageUsingAffineTransfo(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,ScalarField * TransformedImage);


///6.3.1: displacement fields

///Compose RefField with UpdateField. The result is saved in RefField  (we consider that a disp. f. points from the source to the target)
void DisplacementFieldCompose(VectorField * RefField,VectorField * UpdateField);

///Compose InvUpdateField with InvRefField. The result is saved in InvRefField  (we consider that an inv. disp. f. points from the target source to the source)
void InvDisplacementFieldCompose(VectorField * InvUpdateField,VectorField * InvRefField);

///We consider here that the vectors of DispField point from 'StaticImage' to 'DeformedImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
///remark: 'DefoField' has to be sufficiently smooth to allow an accurate inversion
void ProjectImageUsingDispField(VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0);


///We consider here that the vectors of InvDispField point from 'DeformedImage' to 'StaticImage'
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
void ProjectImageUsingInvDispField(VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0);



///Project 'StaticImage' into 'DeformedImage'
//-> 'ProjectCS_2_OriginCS' first projects 'StaticImage' from its own coordinate system to the one of 'DeformedImage' (eventually by integrating an affine mapping)
//       (It actually encodes the affine transformation from 'DeformedImage' to 'StaticImage')
//-> 'InvDispField' then projects 'StaticImage' to 'DeformedImage' 
//       (It actually encodes the inverse transformation from 'DeformedImage' to 'StaticImage')
void ProjectImageUsingAffineTransfoAndInvDispField(float ProjectCS_2_OriginCS[4][4],VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN=0);
void ProjectImageUsingDispFieldAndInvDispField(VectorField * DispField,VectorField * InvDispField,ScalarField * StaticImage,ScalarField * DeformedImage, int NN=0);


void ProjectImageUsingAffineTransfoAndDispField(float ProjectCS_2_OriginCS[4][4],VectorField * DispField,ScalarField * StaticImage,ScalarField * DeformedImage);


///6.3.2: stationary (steady) velocity fields


///Integrate a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptDefFromSteadyVeloField(VectorField * VeloField,VectorField * DeformationField,int log2TimeStepNb);

///Integrate the inverse of a steady velocity field using 2^[log2TimeStepNb] time steps (technique of Arsigny MICCAI 2006)
void CptInvDefFromSteadyVeloField(VectorField * VeloField,VectorField * InvDeformationField,int log2TimeStepNb);


///compute the Lie Braket of VF1 and VF2. Put the result in VF3.   (subfunction of ComposeTwoLogFieldsUsingBCH)
void LieBracket(VectorField * VF1,VectorField * VF2,VectorField * VF3);

///RefVeloField and UpdateVeloField are two steady velocity fields. Their exponentials are respectively the current deformation
///and the update of the deformation. This function approximates the velocity field which is the log of the composition 
///between the current deformation and the update deformation (Baker-Campbell-Hausdorff formula) (cf Vercauteren MICCAI 2008)
///In output, RefVeloField is the updated velocity field. UpdateVeloField is also modified for algorithmic reasons but it
///represents nothing pertinent as an output.
void ComposeTwoLogFieldsUsingBCH(VectorField * RefVeloField,VectorField * UpdateVeloField);
void ComposeTwoLogFieldsUsingSum(VectorField * RefVeloField,VectorField * UpdateVeloField);

///... explicit name ....
///If NearestNgbh==1, the nearest neighbor is projected / otherwise tri-linear interpolation is used
void ProjectImageUsingInverseSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0);

void ProjectImageUsingSteadyVeloField(VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NearestNgbh=0);
void ProjectImageUsingAffineTransfoAndSteadyVeloField(float ProjectCS_2_OriginCS[4][4],VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN=0);
void ProjectImageUsingDispFieldAndSteadyVeloField(VectorField * DispField,VectorField * VeloField,ScalarField * StaticImage,ScalarField * DeformedImage,int NN=0);

///6.3.3: time-dependent velocity fields (no momenta)


//Compute the 3D+t mapping 'Map' from the time step 'refTimeStep' by following the velocity field 'VeloField'.
//'MappingAtRefTimeStep' is the mapping at refTimeStep (possibly not the identity).
//* An iterative leap-frog like technique is performed to compute the backward mapping. The more
//  iterations (='ConvergenceSteps') the more accurate the mapping. For most of the deformations
//  1 iteration is far enough but more iterations are suitable if the deformations have large 
//  Jacobians.
//* 'DeltaX' is the spatial step between two voxels.
void CptMappingFromVeloField(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps=1,float DeltaX=1);
void CptMappingFromVeloField_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * Map,int ConvergenceSteps=1,float DeltaX=1);

//Do the same thing as the initialisation of CptMappingFromVeloField -> load the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
//The function 'CptMappingFromVeloField2_Increment' then allows to increment the field backward or forward according to 'VeloField'
void CptMappingFromVeloField2_Init(VectorField * MappingAtRefTimeStep,VectorField * Map);
void CptMappingFromVeloField2_Init_IniIdMap(VectorField * Map);

//Do the same thing as an incrementation of CptMappingFromVeloField from 'CurrentTimeStep' and backward (BackwardOrForward==-1) or forward (BackwardOrForward==1) 
//The function 'CptMappingFromVeloField2_Init' is then supposed to have loaded the mapping 'MappingAtRefTimeStep' in the 3D FIELD 'Map'
void CptMappingFromVeloField2_Increment(VectorField * VeloField,VectorField * Map,int CurrentTimeStep,int BackwardOrForward,int ConvergenceSteps=1,float DeltaX=1);

//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialMapping' which is the partial mapping of 'MappingAtRefTimeStep' from the time 
//subdivision 'refTimeStep' due to the contribution of 'PartialVeloField'. Note, that an Identity mapping 'MappingId' is //also defined in the inputs (to avoid defining it each time the function is used)
void CptPartialMappingFromVeloFields(int refTimeStep,VectorField * MappingAtRefTimeStep,VectorField * MappingId,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps=1,float DeltaX=1);
void CptPartialMappingFromVeloFields_IniIdMap(int refTimeStep,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialMapping,int ConvergenceSteps=1,float DeltaX=1);

//We consider here that 'PartialVeloField' contributes to 'VeloField'   (VeloField= [A velocity field] + PartialVeloField).
//This function then computes 'PartialLocMap' which is the partial mapping ONLY AT 'TargetSubdiv' FROM 'SourceSubdiv' due to the contribution of PartialVeloField.
//-> PartialLocMap therefore represents where are the coordinates of the points of time subdivision 'SourceSubdiv' when transported on time subdivision 'TargetSubdiv'
void ComputeLagrangianPartialMapping(int SourceSubdiv,int TargetSubdiv,VectorField * VeloField,VectorField * PartialVeloField,VectorField * PartialLocBmap,float DeltaX=1);


//Compute the projection of a 3D image 'ImagToPropag' using Mapping 'Map'.
//The image is projected at the time step 'TimeStepProj' of 'Map' and stored in 'ImageTimeT'.
//
//Importantly, we consider that 'Map' has an identity transformation at the time step (say 't') where 'ImagToPropag' is.
//'Project3Dimage' can therefore perform a forward projection (t<TimeStepProj) or a backward projection (TimeStepProj<t).
void Project3Dimage(ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj);

void Project3DImageUsingAffineTransfoAndTimeDepVF(float ProjectCS_2_OriginCS[4][4],ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj);

void Project3DImageUsingDispFieldAndTimeDepVF(VectorField * DispField,ScalarField * ImagToPropag,VectorField * Map,ScalarField * ImageTimeT,int TimeStepProj);

void Project3Dimage(VectorField * VFToPropag,VectorField * Map,VectorField * VFTimeT,int TimeStepProj);



//By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
//'VeloField4Measure' in the length of the flow from each point of the field. The length of flow
//is projected AT T=0 and returned in the 3D scalar field 'LengthOfFlow'
// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
//   is computed.
void CptLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps=3,float DeltaX=1);

//By following the flow defined by the velocity field 'VeloField4Flow' measure the contribution of
//'VeloField4Measure' AT THE CURRENT TIME in the length of the flow from each point of the field. The length of flow
//is returned in the 3D+t scalar field 'LengthOfFlow'
// * 'VeloField4Measure' is assumed to be part of a linear decomposition of 'VeloField4Flow'.
// * If 'VeloField4Measure'=='VeloField4Flow' then the length of the flow defined by 'VeloField4Flow'
//   is computed.
void CptEvoLengthOfFlow(VectorField * VeloField4Flow,VectorField * VeloField4Measure,ScalarField * LengthOfFlow,int ConvergenceSteps=3,float DeltaX=1);



///6.3.4: time-dependent velocity fields (with momenta)

// Computes the transport of the initial momentum by the diffeo and stores it in Momentum
void TransportMomentum(ScalarField *InitialMomentum, VectorField *InvDiffeo, ScalarField *Momentum,float DeltaX,int t=0);

// Computes cste * transport momentum from the initial momentum by the diffeo and add it in Image.
void AddTransportMomentum(ScalarField *InitialMomentum,VectorField *TempInvDiffeo, ScalarField *Momentum,float DeltaX,float cste=1.0, int t=0);

// Computes the transport image from the initial image by the diffeo and stores it in Image.
void TransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image,int t=0);

// Computes cste * transport image from the initial image by the diffeo and add it in Image.
void AddTransportImage(ScalarField *InitialImage, VectorField *TempInvDiffeo, ScalarField *Image,float cste=1.0,int t=0);



///+++++++++++++++++++++++++++++++++      6.4: linear algebra       +++++++++++++++++++++++++++++++++


// Copies the values of a VectorField1(t=0) in VectorField2(t)
void DeepCopy(VectorField *VectorField1,VectorField *VectorField2,int t=0);
void DeepCopy(VectorField *VectorField,ScalarField* ScalarField,int direc,int t=0);

// Copies the values of a ScalarField1(t=0) in ScalarField2(t)
void DeepCopy(ScalarField *ScalarField1,ScalarField *ScalarField2,int t=0);

// Compute the L^2 scalar product and store it in ScalarField0
void ScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t=0,float cste = 1.0);
void ScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t=0, float cste = 1.0);

// Compute the L^2 scalar product between two vectorfields at time t and add it to ScalarField0
void AddScalarProduct(VectorField *VectorField1, VectorField *VectorField2, ScalarField *ScalarField0, int t=0);
void AddScalarProduct(ScalarField *ScalarField1, ScalarField *ScalarField2, ScalarField *ScalarField0, int t=0);

// Add  cste * ScalarField1 at time t1 to ScalarField2 at time t2
void AddScalarField(ScalarField *ScalarField1, ScalarField *ScalarField2,float cste,int t1 = 0,int t2=0);
// Add  cste * VectorField1 at time t1 to VectorField2 at time t2
void AddVectorField(VectorField *VectorField1, VectorField *VectorField2,float cste,int t1 = 0,int t2=0);
// Multiply a vector field by the cste
void MultiplyVectorField(VectorField *VectorField1, float cste,int t=0);
// 
void SumVectorField(VectorField *VectorField1, VectorField *VectorField2, VectorField *Output, int t1=0,int t2=0,int t3=0, float cste1 = 1.0,float cste2 =1.0);
// Compute the product element by element of ScalarField and VectorField and store it in VectorField2
void Product(ScalarField *ScalarField, VectorField *VectorField1, VectorField *VectorField2);

// Compute the dot product
float DotProduct(ScalarField *ScalarField1, ScalarField *ScalarField2,int t1=0,int t2=0);



///+++++++++++++++++++++++++++++++++      6.5: similarity measures       +++++++++++++++++++++++++++++++++


//compute the L_2 norm of the difference between two scalar fields
float CalcSqrtSumOfSquaredDif(ScalarField * I1,ScalarField * I2);



///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                   7: OTHER FUNCTIONS OF SCIENTIFIC COMPUTATION 
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


//Solve the problem: MX=D where D is a known vector, M a tridiagonal matrix and X the unknown vector.
// Inputs are a,b,c,d,n where  M(i,i-1)=a(i), M(i,i)=b(i), M(i,i+1)=c(i), D(i)=d(i), D in R^n and M in R^n*R^n.
// Output is X where X in R^n.  Warning: will modify c and d!
void TridiagonalSolveFloat(const float *, const float *, float *, float *, float *, int);


//Perform the eigenvalue decomposition of a 3*3 matrix
//Adapted from the algorithm having the same name in 'numerical recipes'.
//Input:  The 3*3 matrix 'MatIni' that has to be symetric.
//Ouputs: 'ValP' is a vector of size 3 which contains the eigen values (in decreasing order). 
//        'VecP' is a 3*3 matrix containg the eigen vectors in columns.
void jacobi3(float **MatIni,float *ValP, float **VecP);


///compute two orthogonal vectors tvec1 and tvec2 in R^3 which are orthogonal to nvec
///the norm of tvec1 and tvec2 is defined as equal to the one of nvec
void CptVecsTangentPlane(float nvec[3],float tvec1[3],float tvec2[3]);

///normalize a vector
void VecNormalize(float vec[3],float norm);



///compute the determinant of a 3*3 matrix
float determinant_3t3matrix(float m[3][3]);

///compute the comatrix of a 3*3 matrix
void comatrix_3t3matrix(float m1[3][3],float m2[3][3]);

///Estimate the exponential of a 3*3 matrix
void Exponential_3t3matrix(float m1[3][3],float m2[3][3]);

///transpose a 3*3 matrix
void transpose_3t3matrix(float m1[3][3],float m2[3][3]);

///inverse of a 3*3 matrix
void invert_3t3matrix(float m1[3][3],float m2[3][3]);

///multiply two 3*3 matrices
void mult_3t3mat_3t3mat(float m1[3][3], float m2[3][3], float MatRes[3][3]);

///multiply a vector of size 3 and a 3*3 matrix
void mult_3t3mat_3vec(float mat[3][3], float vectIni[3], float vectRes[3]);

///inverse of a quaternion
void invert_4t4quaternion(float q1[4][4],float q2[4][4]);

///multiply a vector of size 4 and a 4*4 matrix
void mult_4t4mat_4vec(float mat[4][4], float vectIni[4], float vectRes[4]);

///multiply two 4*4 matrix representing quaternions: mat_i1 * mat_i2 -> mat_o
void mult_quat4t4mat_quat4t4mat(float mat_i1[4][4], float mat_i2[4][4], float mat_o[4][4]);

///read a 4*4 matrix in a text file
void Read_quat4t4mat(char *,float locmat[4][4]);

///write a 4*4 matrix in a text file
void Write_quat4t4mat(char *,float locmat[4][4]);


///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
///                                            8: LANDMARKS
///++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

/// 8.1) ++++++++++++++++++ point landmarks ++++++++++++++++++ 
class LDMK_Points{
  /// ******************************************************************************
private:
  //number of LDMK_Points
  int LDMK_Points_Nb;
  
  //landmark coordinates
  float * Lx;
  float * Ly;
  float * Lz;
  
  /// ******************************************************************************
public:
  
  /// Constructor and destructor
  LDMK_Points();
  ~LDMK_Points();
  
  ///read LDMK_Points in a CSV file
  virtual void Read(char *);
  
  ///Get the X, Y, Z coordinates of the LDMK_Points
  virtual float GetX(int Id=0);
  virtual float GetY(int Id=0);
  virtual float GetZ(int Id=0);
  
  ///Get the number of LDMK_Points
  virtual int Get_LDMK_PointsNumber(void);
};
  
  
  
  
  
/// 8.2) ++++++++++++++++++ curve landmarks ++++++++++++++++++ 
class LDMK_Curves{
  /// ******************************************************************************
private:
  int NbSeg;       //segments number (amount of curves)
  int *NbEl;       //elements number in each segment (amount of nodes in each segment). NbEl[i]=0 means no memory is allocated at segment i
  float **x;      //x[i][j] coordinate of segment j / element i on the x axis
  float **y;      //y[i][j] coordinate of segment j / element i on the y axis
  float **z;      //z[i][j] coordinate of segment j / element i on the z axis
  float **d;      //d[i][j] diametre of segment j / element i
  
  /// ******************************************************************************
public:
  
  /// Constructors and destructor
  LDMK_Curves();
  LDMK_Curves(int SegNb, int ElNb);
  ~LDMK_Curves();
  
  ///read LDMK_Curves in a mv3d file
  virtual void Read(char *);
  
  ///write LDMK_Curves in a mv3d file
  ///Set Preserve_IDs to 1 to preserve the original segment identifiers (default). They are optimaly resampled otherwise.
  virtual void Write(char *, int Preserve_IDs=1);

  ///Get the coordinate value of x, y, z or get the diameter of the element 'IdEl' of segment 'IdSeg'
  virtual inline float GetX(int IdSeg, int IdEl){return x[IdSeg][IdEl];}
  virtual inline float GetY(int IdSeg, int IdEl){return y[IdSeg][IdEl];}
  virtual inline float GetZ(int IdSeg, int IdEl){return z[IdSeg][IdEl];}
  virtual inline float GetD(int IdSeg, int IdEl){return d[IdSeg][IdEl];}
  
  ///Put a coordinate value for x, y, z or put a diameter to the element 'IdEl' of segment 'IdSeg'
  virtual inline void PutX(float value, int IdSeg, int IdEl){x[IdSeg][IdEl]=value;}
  virtual inline void PutY(float value, int IdSeg, int IdEl){y[IdSeg][IdEl]=value;}
  virtual inline void PutZ(float value, int IdSeg, int IdEl){z[IdSeg][IdEl]=value;}
  virtual inline void PutD(float value, int IdSeg, int IdEl){d[IdSeg][IdEl]=value;}
  
  ///get the number of segments
  virtual inline int GetSegNumber(){return NbSeg;}
  
  ///get the number of elements in segment 'IdSeg'
  virtual inline int GetElNumber(int IdSeg){return NbEl[IdSeg];}
  
  ///Merge two segments at their nearest extremity. The new segment is saved in Seg1.
  virtual void MergeSegments(int Seg1, int Seg2);

  ///Delete a segment
  virtual void DeleteSegment(int Seg);

  ///count the number of segments related to the segment end SegEnd (= 0 pr 1) of the segment Seg
  virtual int CountLinkedSeg(int Seg,int SegEnd,float epsilon=0.0001);

  ///Clean-up a moderately large LDMK_Curves structure
  /// -> two segments ends are supposed linked if their distance is less than epsilon
  /// -> merge the segments linked by a node with only two segments
  /// -> remove the segments linked to only one node and for which the node has more than two segments related to other nodes
  virtual void CleanUp(float epsilon=0.0001);

  ///transform the LDMK_Curves coordinates and diameters from voxels to mm according to the image 2 world properties of RefSF 
  virtual void VoxelsToMillimeters(ScalarField * RefSF);

  ///transform the LDMK_Curves coordinates and diameters from mm to voxels according to the image 2 world properties of RefSF 
  virtual void MillimetersToVoxels(ScalarField * RefSF);

  
};
  
  


#endif
