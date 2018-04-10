# -*- coding: utf-8 -*-

#BEFORE USING THIS FILE:
#-> We recommend itksnap to vizualize the results


import numpy as np
import os
from numpy.linalg import inv
import sys
from xml.etree import ElementTree
import os.path as op
import glob
import SimpleITK as sitk


#import tkinter
if sys.version_info[0] < 3:
    import Tkinter as Tk
    import tkFileDialog as TkFd
else:
    import tkinter as Tk
    import filedialog as TkFd

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                         FUNCTIONS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def ShowSlice_x_y_z(SitkImage3D,SliceId_x,SliceId_y,SliceId_z):
  tmpArray=sitk.GetArrayFromImage(SitkImage3D)
  plt.figure()
  plt.imshow(tmpArray[SliceId_x,:,:],cmap='Greys')
  plt.colorbar()
  plt.title('x='+str(SliceId_x))
  plt.figure()
  plt.imshow(tmpArray[:,SliceId_y,:],cmap='Greys')
  plt.colorbar()
  plt.title('y='+str(SliceId_y))
  plt.figure()
  plt.imshow(tmpArray[:,:,SliceId_z],cmap='Greys')
  plt.colorbar()
  plt.title('z='+str(SliceId_z))
  plt.show()
  
def ShowImDiff_x_y_z(SitkImage3D_a,SitkImage3D_b,SliceId_x,SliceId_y,SliceId_z):
  tmpArray_a=sitk.GetArrayFromImage(SitkImage3D_a)
  tmpArray_b=sitk.GetArrayFromImage(SitkImage3D_b)
  plt.figure()
  plt.imshow(tmpArray_a[SliceId_x,:,:]-tmpArray_b[SliceId_x,:,:],cmap='bwr')
  plt.title('x='+str(SliceId_x))
  plt.colorbar()
  plt.figure()
  plt.imshow(tmpArray_a[:,SliceId_y,:]-tmpArray_b[:,SliceId_y,:],cmap='bwr')
  plt.colorbar()
  plt.title('y='+str(SliceId_y))
  plt.figure()
  plt.imshow(tmpArray_a[:,:,SliceId_z]-tmpArray_b[:,:,SliceId_z],cmap='bwr')
  plt.colorbar()
  plt.title('z='+str(SliceId_z))
  plt.show()


def ShowVectorField_x_y_z(SitkVectorField3D,SliceId_x,SliceId_y,SliceId_z):
  tmpArray=sitk.GetArrayFromImage(SitkVectorField3D)
  #SliceId_x
  Y, Z = np.meshgrid(np.arange(0, tmpArray.shape[1], 1), np.arange(0, tmpArray.shape[2], 1))
  vY=tmpArray[SliceId_x,Y,Z,1]
  vZ=tmpArray[SliceId_x,Y,Z,2]
  tmpNorm=np.sqrt(vY*vY+vZ*vZ)
  plt.figure()
  plt.quiver(Y, Z, vY, vZ)
  plt.title('x='+str(SliceId_x)+'  -- norm in ['+str(tmpNorm.min())+','+str(tmpNorm.max())+')')
  #SliceId_y
  X, Z = np.meshgrid(np.arange(0, tmpArray.shape[0], 1), np.arange(0, tmpArray.shape[2], 1))
  vX=tmpArray[X,SliceId_y,Z,0]
  vZ=tmpArray[X,SliceId_y,Z,2]
  tmpNorm=np.sqrt(vX*vX+vZ*vZ)
  plt.figure()
  plt.quiver(X, Z, vX, vZ)
  plt.title('y='+str(SliceId_y)+'  -- norm in ['+str(tmpNorm.min())+','+str(tmpNorm.max())+')')
  #SliceId_z
  X, Y = np.meshgrid(np.arange(0, tmpArray.shape[0], 1), np.arange(0, tmpArray.shape[1], 1))
  vX=tmpArray[X,Y,SliceId_z,0]
  vY=tmpArray[X,Y,SliceId_z,1]
  tmpNorm=np.sqrt(vX*vX+vY*vY)
  plt.figure()
  plt.quiver(X, Y, vX, vY)
  plt.title('z='+str(SliceId_z)+'  -- norm in ['+str(tmpNorm.min())+','+str(tmpNorm.max())+')')
  plt.show()


def Convert_Voxel2Millimeter(VoxCoord,RefSitkImage):
  """
  Convert a voxel coordinate of RefSitkImage into millimeter coordinate in the physcal space
  Inuts:
  -> VoxCoord is a list or a numpy array of size 3 (x,y,z)
  -> RefSitkImage is a simple itk image
  Output:
    -> The mm coordinates as a numpy array of size 3 and dtype=np.float64
  """
  
  Orig=RefSitkImage.GetOrigin()
  Spac=RefSitkImage.GetSpacing()
  DirX=RefSitkImage.GetDirection()[0:3]
  DirY=RefSitkImage.GetDirection()[3:6]
  DirZ=RefSitkImage.GetDirection()[6:9]
  
  xmm=Orig[0]+Spac[0]*(DirX[0]*VoxCoord[0]+DirX[1]*VoxCoord[1]+DirX[2]*VoxCoord[2])
  ymm=Orig[1]+Spac[1]*(DirY[0]*VoxCoord[0]+DirY[1]*VoxCoord[1]+DirY[2]*VoxCoord[2])
  zmm=Orig[2]+Spac[2]*(DirZ[0]*VoxCoord[0]+DirZ[1]*VoxCoord[1]+DirZ[2]*VoxCoord[2])
  
  return np.array([xmm , ymm , zmm],dtype=np.float64)



def DefineMaskWithOnes(InputSitkImage):
  """
  Create a mask containing ones everywhere.
  Input:
   -> InputSitkImage: a sitk image
  Output:
   -> the mask as a sitk image in the same physical domain as InputSitkImage
  """
  tmpMask=sitk.GetArrayFromImage(InputSitkImage)
  tmpMask=tmpMask*0.+1.
  SitkMask=sitk.GetImageFromArray(tmpMask,isVector=False)
  SitkMask=sitk.Cast(SitkMask,sitk.sitkInt8)   #to alternatively copy the type of InputSitkImage -> InputSitkImage.GetPixelID()
  SitkMask.SetDirection(InputSitkImage.GetDirection())
  SitkMask.SetOrigin(InputSitkImage.GetOrigin())
  SitkMask.SetSpacing(InputSitkImage.GetSpacing())
  
  return SitkMask


def SmoothImage(InputSitkImage,KernelSigma):
  """
  Smooth a sitk image with a gaussian kernel of stdev sigma.
  Inputs:
   -> InputSitkImage: a sitk image
   -> KernelSigma: stdev of the gaussian kernel
  Output:
   -> the smoothed image as a sitk image in the same physical domain as InputSitkImage
  """
  
  return sitk.SmoothingRecursiveGaussian(InputSitkImage,sigma=KernelSigma)


def ComputeLandmarkBasedRigid3DTransform(FixSitkImage,FixListLandmarks,MovSitkImage,MovListLandmarks):
  """
  Return a transform of type sitk.VersorRigid3DTransform which best aligns the paired landmarks of FixListLandmarks and MovListLandmarks.
  The list of landmark must have the followinf structure: [[x1,y1,z1],[x2,y2,z2],...]
  Inputs:
   -> FixSitkImage: the fixed sitk image
   -> FixListLandmarks: the list of landmarks out of FixSitkImage in VOXEL coordinates
   -> MovSitkImage: the moving sitk image
   -> MovListLandmarks: the list of landmarks out of MovSitkImage in VOXEL coordinates
  Output:
   -> the optimal transform of type sitk.VersorRigid3DTransform 
  """
  
  #concatenate the landmarks of FixSitkImage
  Fix_LDMKs=Convert_Voxel2Millimeter(FixListLandmarks[0],FixSitkImage)
  for i in range(1,len(FixListLandmarks)):
    Fix_LDMKs=np.concatenate((Fix_LDMKs,Convert_Voxel2Millimeter(FixListLandmarks[i],FixSitkImage)))

  #concatenate the landmarks of MovSitkImage
  Mov_LDMKs=Convert_Voxel2Millimeter(MovListLandmarks[0],MovSitkImage)
  for i in range(1,len(MovListLandmarks)):
    Mov_LDMKs=np.concatenate((Mov_LDMKs,Convert_Voxel2Millimeter(MovListLandmarks[i],MovSitkImage)))
  
  #create a vector of weights
  LDMKs_wgt=Mov_LDMKs*0.+1.

  #compute the transform
  Coarse_transform = sitk.LandmarkBasedTransformInitializer(sitk.VersorRigid3DTransform(), Fix_LDMKs,Mov_LDMKs,LDMKs_wgt,FixSitkImage)
  
  return Coarse_transform

def ComputeAutomaticRigid3DTransform(FixSitkImage,MovSitkImage,MovMask,init_transform):
  """
  Return a transform of type sitk.VersorRigid3DTransform which best aligns FixSitkImage and MovSitkImage.
  Inputs:
   -> FixSitkImage: the fixed sitk image
   -> MovSitkImage: the moving sitk image
   -> MovMask: mask on MovSitkImage
   -> init_transform: the initial tranform of type sitk.VersorRigid3DTransform   (WARNING: VERY LIKELY TO BE MODIFIED)
  Output:
   -> the optimal transform of type sitk.VersorRigid3DTransform 
  """
  
  #step 1
  RigRegPerformer = sitk.ImageRegistrationMethod()
  #RigRegPerformer.SetMetricAsMeanSquares()
  RigRegPerformer.SetMetricAsANTSNeighborhoodCorrelation(radius=5)
  
  RigRegPerformer.SetMetricMovingMask(MovMask)
  
  
  RigRegPerformer.SetInterpolator(sitk.sitkLinear)
  #RigRegPerformer.SetOptimizerAsGradientDescent(learningRate= 1.0, numberOfIterations=100, convergenceWindowSize=10)
  RigRegPerformer.SetOptimizerAsGradientDescent(learningRate=1.0, numberOfIterations=30)
  RigRegPerformer.SetOptimizerScalesFromIndexShift()
  RigRegPerformer.SetInitialTransform(init_transform)   #warning: RigReg_transform1 will be transformed when executing the class
  
  RigRegPerformer.SetSmoothingSigmasAreSpecifiedInPhysicalUnits(True)
  RigRegPerformer.SetShrinkFactorsPerLevel(shrinkFactors = [4,2,2])
  RigRegPerformer.SetSmoothingSigmasPerLevel(smoothingSigmas=[2.,2.,1.])
  
  tmp_transform=RigRegPerformer.Execute(FixSitkImage,MovSitkImage)
  
  #step 2
  RigRegPerformer = sitk.ImageRegistrationMethod()
  RigRegPerformer.SetMetricAsANTSNeighborhoodCorrelation(radius=3)
  
  RigRegPerformer.SetMetricMovingMask(MovMask)
  
  RigRegPerformer.SetInterpolator(sitk.sitkLinear)
  RigRegPerformer.SetOptimizerAsGradientDescent(learningRate= 0.3, numberOfIterations=10)
  RigRegPerformer.SetOptimizerScalesFromIndexShift()
  RigRegPerformer.SetInitialTransform(tmp_transform) 
  RigRegPerformer.SetSmoothingSigmasAreSpecifiedInPhysicalUnits(True)
  RigRegPerformer.SetShrinkFactorsPerLevel(shrinkFactors = [1])
  RigRegPerformer.SetSmoothingSigmasPerLevel(smoothingSigmas=[0.1])
  
  return RigRegPerformer.Execute(FixSitkImage,MovSitkImage)


def ComputeLogDemonsTransform(FixSitkImage,MovSitkImage,TypicalBrainWidth):   #MovSitkMask
  """
  LogDemons registration of FixSitkImage and MovSitkImage. Images are supposed to be in the same domain and to have the same size.
  Inputs:
   -> FixSitkImage: the fixed sitk image
   -> MovSitkImage: the moving sitk image
   -> TypicalBrainWidth: typical width of the studied brain in millimeters
  Output:
   -> the optimal mapping as a sitk.DisplacementFieldTransform
  """

  
  #In Klein et al.: Gauss[2,0]  (here typical brain width and lenght is supposed 4 * smaller than in human brain)
  #Here 140 is the typical brain width in humans and is the 2 std val used in Klein et al, NeuroImage 51(1), 2010
  stdGauss=2*TypicalBrainWidth/140   
  PixelRes=np.array(FixSitkImage.GetSpacing())   #[res. x, res. y, res. z]
  
  #step 1
  DiffDem=sitk.DiffeomorphicDemonsRegistrationFilter()
  DiffDem.SetNumberOfIterations(10)
  DiffDem.SetSmoothDisplacementField(True)
  DiffDem.SetStandardDeviations([2,2,2])
  DiffDem.SetSmoothUpdateField(True)
  DiffDem.SetUpdateFieldStandardDeviations([5.*stdGauss/PixelRes[0],5.*stdGauss/PixelRes[1],5.*stdGauss/PixelRes[2]])
  Deformation=DiffDem.Execute(FixSitkImage,MovSitkImage)
    
  return sitk.DisplacementFieldTransform(Deformation)



def ComputeDemonsTransform(FixSitkImage,MovSitkImage,MovMask):
    """
    Demons registration of FixSitkImage and MovSitkImage. Images are supposed to be in the same domain and to have the same size.
    Inputs:
     -> FixSitkImage: the fixed sitk image
     -> MovSitkImage: the moving sitk image
     -> MovMask: mask on MovSitkImage
    Output:
     -> the optimal mapping as a sitk.DisplacementFieldTransform
    """
    
    PixelRes=np.array(FixSitkImage.GetSpacing())   #[res. x, res. y, res. z]
    
    
    FixSitkImage=sitk.Cast(sitk.RescaleIntensity(FixSitkImage), MovSitkImage.GetPixelIDValue())
    
    #Step 1
    
    transform_to_displacment_field_filter = sitk.TransformToDisplacementFieldFilter()
    transform_to_displacment_field_filter.SetReferenceImage(FixSitkImage)
    
    initial_transform = sitk.DisplacementFieldTransform(transform_to_displacment_field_filter.Execute(sitk.Transform()))
    initial_transform.SetSmoothingGaussianOnUpdate(varianceForUpdateField=9./PixelRes[0], varianceForTotalField=0.3/PixelRes[0]) 
    
    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetInitialTransform(initial_transform)
    registration_method.SetMetricAsANTSNeighborhoodCorrelation(radius=3)
    registration_method.SetMetricMovingMask(MovMask)
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=0.5, numberOfIterations=10, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()
    
    Deformation=registration_method.Execute(FixSitkImage, MovSitkImage)
    
    
    #Step 2
    
    initial_transform = sitk.DisplacementFieldTransform(Deformation)
    initial_transform.SetSmoothingGaussianOnUpdate(varianceForUpdateField=3./PixelRes[0], varianceForTotalField=0.3/PixelRes[0]) 
    
    registration_method = sitk.ImageRegistrationMethod()
    registration_method.SetInitialTransform(initial_transform)
    registration_method.SetMetricAsANTSNeighborhoodCorrelation(radius=3)
    registration_method.SetMetricMovingMask(MovMask)
    registration_method.SetInterpolator(sitk.sitkLinear)
    registration_method.SetOptimizerAsGradientDescent(learningRate=0.5, numberOfIterations=5, convergenceMinimumValue=1e-6, convergenceWindowSize=10)
    registration_method.SetOptimizerScalesFromPhysicalShift()
    
    Deformation=registration_method.Execute(FixSitkImage, MovSitkImage)
    
    return sitk.DisplacementFieldTransform(Deformation)



def ResampleMovingImage(FixSitkImage,MovSitkImage,transfo,NN=0):
  """
  Resample the moving image MovSitkImage in the domain of FixSitkImage using the transform transfo
  inputs:
   -> FixSitkImage
   -> MovSitkImage
   -> transfo
   -> NN (optional): use nearest neighbor interpolation if == 1
  """
  
  if (NN==1):
    return sitk.Resample(MovSitkImage,FixSitkImage,transfo, sitk.sitkNearestNeighbor, 0.0, MovSitkImage.GetPixelIDValue())
  else:
    return sitk.Resample(MovSitkImage,FixSitkImage,transfo, sitk.sitkLinear, 0.0, MovSitkImage.GetPixelIDValue())


def ResampleSetOfMovingImages(FixSitkImage,PrefixMovImageFiles,RigRegTransform,DispFieldTransform):
  """
  Resample a set moving images in the domain of FixSitkImage using the rigid transform RigRegTransform 
  and then the displacement field DispFieldTransform 
  inputs:
   -> PrefixMovImageFiles: prefix of the images to resample
   -> FixSitkImage: fixed image 
   -> RigRegTransform: rigid alignement
   -> DispFieldTransform: non-rigid alignement
  """  
  
  
  Lst_ImageFilesToResample=glob.glob(PrefixMovImageFiles+'*')
  
  for InputImageFile in Lst_ImageFilesToResample:
    #print ('Resample image '+InputImageFile)
    OutputImageFile=InputImageFile.replace('.nii','_RSP.nii')
    tmpIm=ReadSitkImage(InputImageFile)
    tmpIm=ResampleMovingImage(FixSitkImage,tmpIm,RigRegTransform,NN=1)
    tmpIm=ResampleMovingImage(FixSitkImage,tmpIm,DispFieldTransform,NN=1)
    WriteSitkImage(tmpIm,OutputImageFile)
  



def ReadSitkImage(ImFile):
  """
  ReadSitkImage(ImFile) -> read the image in ImFile and return it as a sitk image of type sitk.sitkFloat32
  """
  
  return sitk.Cast(sitk.ReadImage(ImFile),sitk.sitkFloat32)


def WriteSitkImage(SitkImage,ImFile):
  """
  WriteSitkImage(SitkImage,ImFile) -> write the sitk image SitkImage in ImFile
  """
  
  sitk.WriteImage(SitkImage,ImFile)

  

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                         MAIN TREATMENTS
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++

def SaveTypicalParameterFiles(OutParamFile):
  """
  All is in the function name
  """
  
  fo = open(OutParamFile, 'w+')
  fo.write('<ImageAnalysisParams title=\'Parameters For Image Analysis\'>\n')
  fo.write('    <ImageInfo title=\'Image to analyse\'>\n')
  fo.write('        <ImageFile ImageFile=\''+op.join('..','MISC','ExampleData','Example1.nii.gz')+'\' />\n')
  fo.write('        <T2insteadofT1 Bool=\'0\' />\n')
  fo.write('    </ImageInfo>\n')
  fo.write('    <Landmarks title=\'Image landmarks\'>\n')
  fo.write('        <LDMK1 title=\'Landmark 1 in voxel coordinates: pituitary gland\'>\n')
  fo.write('            <LDMK1_x value=\'48\' />\n')
  fo.write('            <LDMK1_y value=\'67\' />\n')
  fo.write('            <LDMK1_z value=\'25\' />\n')
  fo.write('        </LDMK1>\n')
  fo.write('        <LDMK2 title=\'Landmark 2 in voxel coordinates: anterior commissure\'>\n')
  fo.write('            <LDMK2_x value=\'48\' />\n')
  fo.write('            <LDMK2_y value=\'69\' />\n')
  fo.write('            <LDMK2_z value=\'39\' />\n')
  fo.write('        </LDMK2>\n')
  fo.write('        <LDMK3 title=\'Landmark 3 in voxel coordinates: posterior corpus callosum\'>\n')
  fo.write('            <LDMK3_x value=\'48\' />\n')
  fo.write('            <LDMK3_y value=\'45\' />\n')
  fo.write('            <LDMK3_z value=\'45\' />\n')
  fo.write('        </LDMK3>\n')
  fo.write('        <LDMK4 title=\'Landmark 4 in voxel coordinates: anterior caudate nucleus head (left)\'>\n')
  fo.write('            <LDMK4_x value=\'40\' />\n')
  fo.write('            <LDMK4_y value=\'81\' />\n')
  fo.write('            <LDMK4_z value=\'47\' />\n')
  fo.write('        </LDMK4>\n')
  fo.write('        <LDMK5 title=\'Landmark 5 in voxel coordinates: anterior caudate nucleus head (right)\'>\n')
  fo.write('            <LDMK5_x value=\'55\' />\n')
  fo.write('            <LDMK5_y value=\'81\' />\n')
  fo.write('            <LDMK5_z value=\'47\' />\n')
  fo.write('        </LDMK5>\n')
  fo.write('        <LDMK6 title=\'Landmark 6 in voxel coordinates: orbital canal extremity (left)\'>\n')
  fo.write('            <LDMK6_x value=\'44\' />\n')
  fo.write('            <LDMK6_y value=\'77\' />\n')
  fo.write('            <LDMK6_z value=\'30\' />\n')
  fo.write('        </LDMK6>\n')
  fo.write('        <LDMK7 title=\'Landmark 7 in voxel coordinates: orbital canal extremity (right)\'>\n')
  fo.write('            <LDMK7_x value=\'52\' />\n')
  fo.write('            <LDMK7_y value=\'77\' />\n')
  fo.write('            <LDMK7_z value=\'30\' />\n')
  fo.write('        </LDMK7>\n')
  fo.write('    </Landmarks>\n')
  fo.write('    <GeneralParameters title=\'General parameters\'>\n')
  fo.write('        <PrefixOutputs str=\''+op.join('.','out_')+'\' />\n')
  fo.write('        <SkipCorticalRegionsResampling Bool=\'0\' />\n')
  fo.write('        <SaveIntermediateRes Bool=\'0\' />\n')
  fo.write('    </GeneralParameters>\n')
  fo.write('</ImageAnalysisParams>\n')
  fo.close()
  

  
def main_ImageAnalysis(parametersFile,RunFile):
  """
  All is in the function name. 
  -> ParametersFile is a string indicating the xml file which contains the parameters for the analysis.
  -> RunFile is the name of the present file given when running the coputations (allows to locate where ImpecBrainSeg is installed)
  """

  
  #1) Read input parameters and define general parameters
  
  #1.1) Parse input xml file and check whether the inputs make sense
  
  #1.1.1) Default parameters in case the input file is not complete
  SaveIntermediateRes=0
  SkipCorticalRegionsResampling=0
  PrefixOutputs=op.join('.','out_')
  CerCoDir=RunFile.replace("ImpecBrainSeg_CLI.py","")
  CerCoTemplateDir=op.join(CerCoDir,'..','Template')
  T2insteadofT1=0
  
  #1.1.2) parse the input file
  
  with open(op.join(parametersFile), 'rt') as f:
      tree = ElementTree.parse(f)
      
  for node in tree.iter():
      if node.tag=='ImageFile':
        File_ImageToSegment=node.attrib['ImageFile']
      if node.tag=='SaveIntermediateRes':
        SaveIntermediateRes=int(node.attrib['Bool'])
      if node.tag=='LDMK1_x':
        LDMK1_x=float(node.attrib['value'])
      if node.tag=='LDMK1_y':
        LDMK1_y=float(node.attrib['value'])
      if node.tag=='LDMK1_z':
        LDMK1_z=float(node.attrib['value'])
      if node.tag=='LDMK2_x':
        LDMK2_x=float(node.attrib['value'])
      if node.tag=='LDMK2_y':
        LDMK2_y=float(node.attrib['value'])
      if node.tag=='LDMK2_z':
        LDMK2_z=float(node.attrib['value'])
      if node.tag=='LDMK3_x':
        LDMK3_x=float(node.attrib['value'])
      if node.tag=='LDMK3_y':
        LDMK3_y=float(node.attrib['value'])
      if node.tag=='LDMK3_z':
        LDMK3_z=float(node.attrib['value'])
      if node.tag=='LDMK4_x':
        LDMK4_x=float(node.attrib['value'])
      if node.tag=='LDMK4_y':
        LDMK4_y=float(node.attrib['value'])
      if node.tag=='LDMK4_z':
        LDMK4_z=float(node.attrib['value'])
      if node.tag=='LDMK5_x':
        LDMK5_x=float(node.attrib['value'])
      if node.tag=='LDMK5_y':
        LDMK5_y=float(node.attrib['value'])
      if node.tag=='LDMK5_z':
        LDMK5_z=float(node.attrib['value'])
      if node.tag=='LDMK6_x':
        LDMK6_x=float(node.attrib['value'])
      if node.tag=='LDMK6_y':
        LDMK6_y=float(node.attrib['value'])
      if node.tag=='LDMK6_z':
        LDMK6_z=float(node.attrib['value'])
      if node.tag=='LDMK7_x':
        LDMK7_x=float(node.attrib['value'])
      if node.tag=='LDMK7_y':
        LDMK7_y=float(node.attrib['value'])
      if node.tag=='LDMK7_z':
        LDMK7_z=float(node.attrib['value'])
      if node.tag=='PrefixOutputs':
        PrefixOutputs=node.attrib['str']
      if node.tag=='SkipCorticalRegionsResampling':
        SkipCorticalRegionsResampling=int(node.attrib['Bool'])
      if node.tag=='T2insteadofT1':
        T2insteadofT1=int(node.attrib['Bool'])     
  
  #1.1.3) correct obvious mistakes in the inputs
  
    
  if op.exists(CerCoTemplateDir)==0:
    print ('ERROR: CerCo template directory '+CerCoTemplateDir+' does not exist')
    sys.exit()
  
  if op.exists(File_ImageToSegment)==0:
    print ('ERROR: Entered image '+File_ImageToSegment+' does not exist')
    sys.exit()
  
  if (SaveIntermediateRes!=0):
    SaveIntermediateRes=1
  
  if (T2insteadofT1!=0):
    T2insteadofT1=1
  
  if (SkipCorticalRegionsResampling!=0):
    SkipCorticalRegionsResampling=1
    
    

  #1.1.4) Linux/mac vs windows commands for system calls
  if "win" in sys.platform.lower():
    loc_mv='move '
    loc_rm='del '
    loc_cp='copy '
  else:
    loc_mv='mv '
    loc_rm='rm '
    loc_cp='cp '

  #1.2) Show input parameters
  print ('Image to segment: '+File_ImageToSegment)
  if (SkipCorticalRegionsResampling==0):
    print ('Skip Cortical Regions Resampling: No')
  else:
    print ('Skip Cortical Regions Resampling: Yes')
  if (T2insteadofT1==0):
    print ('Reference image is a T1-weighted MR image')
  else:
    print ('Reference image is a T2-weighted MR image')
  print ('Voxel location of landmark 1: ('+str(LDMK1_x)+','+str(LDMK1_y)+','+str(LDMK1_z)+')')
  print ('Voxel location of landmark 2: ('+str(LDMK2_x)+','+str(LDMK2_y)+','+str(LDMK2_z)+')')
  print ('Voxel location of landmark 3: ('+str(LDMK3_x)+','+str(LDMK3_y)+','+str(LDMK3_z)+')')
  print ('Voxel location of landmark 4: ('+str(LDMK4_x)+','+str(LDMK4_y)+','+str(LDMK4_z)+')')
  print ('Voxel location of landmark 5: ('+str(LDMK5_x)+','+str(LDMK5_y)+','+str(LDMK5_z)+')')
  print ('Voxel location of landmark 6: ('+str(LDMK6_x)+','+str(LDMK6_y)+','+str(LDMK6_z)+')')
  print ('Voxel location of landmark 7: ('+str(LDMK7_x)+','+str(LDMK7_y)+','+str(LDMK7_z)+')')
  print ('Prefix for the outputs: '+PrefixOutputs)
  
  
  #2) Registration of the template to the image to segment
  
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  print ('+    BEGIN TEMPLATE ALIGNMENT WITH THE SEGMENTED IMAGE        +')
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
  #2.1) semi automatic rigid alignement
  
  
  print('Coarse alignment')
  
  #2.1.1) Load the target image (image to segment - fixed image) and corresponding landmarks for the coarse rigid alignment
  
  #... L1 -> pituitary gland  / L2 -> anterior commissure  / L3 -> posterior corpus callosum  / L4 -> anterior caudate nucleus (left)  / L5 -> anterior caudate nucleus (right)  / L6 -> extremity orbital canal (left)  / L7 -> extremity orbital canal (right)
  FixImage=ReadSitkImage(File_ImageToSegment)
  FixListLandmarks=[[LDMK1_x,LDMK1_y,LDMK1_z],[LDMK2_x,LDMK2_y,LDMK2_z],[LDMK3_x,LDMK3_y,LDMK3_z],[LDMK4_x,LDMK4_y,LDMK4_z],[LDMK5_x,LDMK5_y,LDMK5_z],[LDMK6_x,LDMK6_y,LDMK6_z],[LDMK7_x,LDMK7_y,LDMK7_z]]
  
  #2.1.2) Reference (moving) image, its mask and landmarks for rigid alignment
  
  if(T2insteadofT1==0):
    File_MovImage=op.join(CerCoTemplateDir,'Template_T1.nii.gz')
  else:
    File_MovImage=op.join(CerCoTemplateDir,'template_T2.nii.gz')
  
  MovImage=ReadSitkImage(File_MovImage)
  
  #... mask
  File_MovMask=op.join(CerCoTemplateDir,'BrainStructures','AllMask.nii.gz')
  MovMask=ReadSitkImage(File_MovMask)
  MovMask=sitk.DilateObjectMorphology(MovMask,4)
  
  #... L1 -> pituitary gland  / L2 -> anterior commissure  / L3 -> posterior corpus callosum  / L4 -> anterior caudate nucleus (left)  / L5 -> anterior caudate nucleus (right)  / L6 -> extremity orbital canal (left)  / L7 -> extremity orbital canal (right)
  MovListLandmarks=[[62 , 56 , 28],[62 , 59 , 44],[61 ,  38 , 52],[54 , 70 , 52],[70 , 70 , 52],[58 , 65 , 33],[68 , 65 , 33]]
  
  #2.1.3) Establish a coarse rigid transformation to match the Reference image on the Target image
  RigReg_Fix2Mov_Coarse=ComputeLandmarkBasedRigid3DTransform(FixImage,FixListLandmarks,MovImage,MovListLandmarks)

  if SaveIntermediateRes==1:
    if op.exists(PrefixOutputs+'IntermediateResults')==0:
      os.mkdir(PrefixOutputs+'IntermediateResults')
  

  if SaveIntermediateRes==1:
    DefImage=ResampleMovingImage(FixImage,MovImage,RigReg_Fix2Mov_Coarse)
    WriteSitkImage(DefImage,PrefixOutputs+'1_aligned_DefImage.nii')
    os.system(loc_mv+' '+PrefixOutputs+'1_aligned_DefImage.nii '+PrefixOutputs+'IntermediateResults')
  
  
  #2.2) automatic refinement of the coarse rigid alignement
  print('Rigid registration')
  
  RigReg_Fix2Mov_Fine=ComputeAutomaticRigid3DTransform(FixImage,MovImage,MovMask,RigReg_Fix2Mov_Coarse)

  DefImage=ResampleMovingImage(FixImage,MovImage,RigReg_Fix2Mov_Fine)
  DefMask=ResampleMovingImage(FixImage,MovMask,RigReg_Fix2Mov_Fine)

  if SaveIntermediateRes==1:
    WriteSitkImage(DefImage,PrefixOutputs+'2_RigReg_DefImage.nii')
    os.system(loc_mv+' '+PrefixOutputs+'2_RigReg_DefImage.nii '+PrefixOutputs+'IntermediateResults')
    WriteSitkImage(DefMask,PrefixOutputs+'2_RigReg_DefMask.nii')
    os.system(loc_mv+' '+PrefixOutputs+'2_RigReg_DefMask.nii '+PrefixOutputs+'IntermediateResults')
  
  #2.3) non-rigid registration
  print('Non-rigid registration')
  
  DispField_Fix2RigRegMov=ComputeDemonsTransform(FixImage,DefImage,DefMask)
  
  
  DefImage2=ResampleMovingImage(FixImage,DefImage,DispField_Fix2RigRegMov)
  
  if SaveIntermediateRes==1:
    WriteSitkImage(DefImage2,PrefixOutputs+'3_DiffeoReg_DefImage.nii')
    os.system(loc_mv+' '+PrefixOutputs+'3_DiffeoReg_DefImage.nii '+PrefixOutputs+'IntermediateResults')
  
  
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  print ('+     END TEMPLATE ALIGNMENT WITH THE SEGMENTED IMAGE         +')
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  
  #3)  segmentation resampling
  
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  print ('+             BEGIN SEGMENTATION RESAMPLING                   +')
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  
  useNN=''   #alternatively set to ' --use-NN' if you want to perform nearest neighbor interpolation
  
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','BrainMask.nii.gz') 
  LocRsp=PrefixOutputs+'BrainMask.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine,NN=1)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov,NN=1)  #NN=1 for nearest neighbor interpolation |  linear by default
  WriteSitkImage(tmpIm,LocRsp)

  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','gm_all.nii.gz')
  LocRsp=PrefixOutputs+'gm_all.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)  
  
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','gm_cortex.nii.gz')
  LocRsp=PrefixOutputs+'gm_cortex.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
  
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','gm_subco.nii.gz')
  LocRsp=PrefixOutputs+'gm_subco.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
  
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','wm.nii.gz')
  LocRsp=PrefixOutputs+'wm.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
    
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','csf.nii.gz')
  LocRsp=PrefixOutputs+'csf.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
    
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','Cerebelum.nii.gz')
  LocRsp=PrefixOutputs+'Cerebelum.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
    
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','hippo.nii.gz')
  LocRsp=PrefixOutputs+'hippo.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
    
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','BrainStem.nii.gz')
  LocRsp=PrefixOutputs+'BrainStem.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
    
  LocSrc=op.join(CerCoTemplateDir,'BrainStructures','OlfactoryBulb.nii.gz')
  LocRsp=PrefixOutputs+'OlfactoryBulb.nii'
  tmpIm=ReadSitkImage(LocSrc)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,RigReg_Fix2Mov_Fine)
  tmpIm=ResampleMovingImage(FixImage,tmpIm,DispField_Fix2RigRegMov)
  WriteSitkImage(tmpIm,LocRsp)
  
  if op.exists(PrefixOutputs+'BrainStructures')==0:
      os.mkdir(PrefixOutputs+'BrainStructures')
  
  os.system(loc_mv+' '+PrefixOutputs+'BrainMask.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'gm_all.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'gm_cortex.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'gm_subco.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'wm.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'csf.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'Cerebelum.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'hippo.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'BrainStem.nii '+PrefixOutputs+'BrainStructures')
  os.system(loc_mv+' '+PrefixOutputs+'OlfactoryBulb.nii '+PrefixOutputs+'BrainStructures')
  
  print ('Estimated brain structures are in the directory '+PrefixOutputs+'BrainStructures')
    
 
  
  
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  print ('+             END  SEGMENTATION RESAMPLING                    +')
  print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')

  
  #4) resample the cortical regions
  if (SkipCorticalRegionsResampling==0):
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('+          BEGIN DEFINITION OF THE CORTICAL AREAS             +')
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    
    PrefixImageFilesToResample=op.join(CerCoTemplateDir,'CorticalRegions','Reg_L')
    ResampleSetOfMovingImages(FixImage,PrefixImageFilesToResample,RigReg_Fix2Mov_Fine,DispField_Fix2RigRegMov)
    
    PrefixImageFilesToResample=op.join(CerCoTemplateDir,'CorticalRegions','Reg_R')
    ResampleSetOfMovingImages(FixImage,PrefixImageFilesToResample,RigReg_Fix2Mov_Fine,DispField_Fix2RigRegMov)
    
    PrefixImageFilesToResample=op.join(CerCoTemplateDir,'CorticalRegions','Reg_GroupL')
    ResampleSetOfMovingImages(FixImage,PrefixImageFilesToResample,RigReg_Fix2Mov_Fine,DispField_Fix2RigRegMov)
    
    PrefixImageFilesToResample=op.join(CerCoTemplateDir,'CorticalRegions','Reg_GroupR')
    ResampleSetOfMovingImages(FixImage,PrefixImageFilesToResample,RigReg_Fix2Mov_Fine,DispField_Fix2RigRegMov)
    
    PrefixImageFilesToResample=op.join(CerCoTemplateDir,'CorticalRegions','All_')
    ResampleSetOfMovingImages(FixImage,PrefixImageFilesToResample,RigReg_Fix2Mov_Fine,DispField_Fix2RigRegMov)
    
    if op.exists(PrefixOutputs+'CorticalRegions')==0:
      os.mkdir(PrefixOutputs+'CorticalRegions')
    
    os.system(loc_mv+' '+op.join(CerCoTemplateDir,'CorticalRegions')+'/*RSP.* '+PrefixOutputs+'CorticalRegions')
    os.system(loc_cp+' '+op.join(CerCoTemplateDir,'CorticalRegions','CortialRegions_Labels.txt')+' '+PrefixOutputs+'CorticalRegions')
    os.system(loc_cp+' '+op.join(CerCoTemplateDir,'CorticalRegions','ITKSnap_AreasLabels.txt')+' '+PrefixOutputs+'CorticalRegions')
    os.system(loc_cp+' '+op.join(CerCoTemplateDir,'CorticalRegions','CortialRegions_Labels.txt')+' '+PrefixOutputs+'CorticalRegions')
    
    print ('Deformed Cortical Regions are in the directory '+PrefixOutputs+'CorticalRegions')
    
    
    
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
    print ('+           END DEFINITION OF THE CORTICAL AREAS              +')
    print ('+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++')
  
  
  print ('')
  print ('All results were saved using the prefix '+PrefixOutputs)
  sys.exit()
  
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                            INTERFACE
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++


def ShowHelpMessage(OutParamFile='ImageAnalysis.xml'):
  """
  show error message
  """
  SaveTypicalParameterFiles(OutParamFile)
  
  print ('')
  print ('A parameter file with default values was generated: \''+OutParamFile+'\'')
  print ('')
  print ('Now, please:')
  print ('  (1) Fill '+OutParamFile+' with your parameters')
  print ('  (2) Run again ImpecBrainSeg_CLI: \'python ImpecBrainSeg_CLI.py\'')
  print ('  (3) Press the button \'Launch computations\'')
  print ('')
  print ('Remark: Instead of steps (2) and (3), you can execute \'python ImpecBrainSeg_CLI.py '+OutParamFile+'\'')
  print ('')
  print ('')



class Application(Tk.Frame):
    
    #function to extract the frames as png
    def LaunchComputations(self):
        self.ParamFileToExec = self.Entry_ParamFileToExec.get()
        
        if self.ParamFileToExec.endswith('.xml')==0:
          ShowHelpMessage()
          sys.exit()

        
        
        with open(self.ParamFileToExec, 'rt') as f:
            tree = ElementTree.parse(f)
        root = tree.getroot()

        if (root.tag=='ImageAnalysisParams'):
          print ('run the image analysis pipeline')
          main_ImageAnalysis(self.ParamFileToExec,self.argv_0)
        else:
          print ('The format of '+self.ParamFileToExec+' is not recognized')
          ShowHelpMessage()
          sys.exit()
          
    #function to extract the frames as csv
    def GenerateParameterFiles(self):
        SvgParamFile=self.Entry_ParamFileToExec.get()
        ShowHelpMessage(OutParamFile=SvgParamFile)
        sys.exit()
    
    #interface
    def createWidgets(self):

        #raw with nothing
        self.Entry_DirTablesText = Tk.Label(self, text="")
        self.Entry_DirTablesText.grid(row=1, column=0)
        
        #button to lauch computations
        self.run_cpt = Tk.Button(self)
        self.run_cpt["text"] = "Launch computations",
        self.run_cpt["command"] = self.LaunchComputations
        self.run_cpt.grid(row=2, column=0)
        
        #button to generate typical parameter files
        self.run_cpt = Tk.Button(self)
        self.run_cpt["text"] = "Generate parameter file",
        self.run_cpt["command"] = self.GenerateParameterFiles
        self.run_cpt.grid(row=2, column=1)
        
        #quit button
        self.QUIT = Tk.Button(self)
        self.QUIT["text"] = "QUIT"
        self.QUIT["fg"]   = "red"
        self.QUIT["command"] =  self.quit
        self.QUIT.grid(row=2, column=2)
        
        #raw with nothing
        self.Entry_DirTablesText = Tk.Label(self, text="")
        self.Entry_DirTablesText.grid(row=3, column=0)

        #Parameter file to exectute
        self.Entry_ParamFileToExecText = Tk.Label(self, text="Parameter file for computations:")
        self.Entry_ParamFileToExecText.grid(row=4, column=0)
        
        defaultFile=op.join('.','ImageAnalysis.xml')
        
        self.Entry_ParamFileToExec = Tk.Entry(self,width=25)
        self.Entry_ParamFileToExec.grid(row=4, column=1)
        self.Entry_ParamFileToExec.delete(0, Tk.END)
        self.Entry_ParamFileToExec.insert(0, defaultFile)
        self.ParamFileToExec = self.Entry_ParamFileToExec.get()
        
        #raw with nothing
        self.Entry_DirTablesText = Tk.Label(self, text="")
        self.Entry_DirTablesText.grid(row=5, column=0)

        
        
    #run the interface
    def __init__(self, master=None):
        Tk.Frame.__init__(self, master)
        self.pack()
        self.createWidgets()

#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++
#                                             MAIN
#++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++



#Check whether an input xml file is given and write typical ones if this is not the case


if len(sys.argv)==1:   #launch the interface
  TKroot = Tk.Tk()
  TKroot.title("ImpecBrainSeg laucher") 
  app = Application(master=TKroot)
  app.argv_0=sys.argv[0]
  app.mainloop()
  TKroot.destroy()
  
if len(sys.argv)>1:  #run directly the computations according to the parameters file
  if sys.argv[1].endswith('.xml')==0:
    ShowHelpMessage()
    sys.exit()
  
  parametersFile=sys.argv[1]

  with open(parametersFile, 'rt') as f:
    tree = ElementTree.parse(f)
  root = tree.getroot()

  if (root.tag=='ImageAnalysisParams'):
    print ('run the image analysis pipeline')
    main_ImageAnalysis(parametersFile,sys.argv[0])
  else:
    print ('The format of '+parametersFile+' is not recognized')
    ShowHelpMessage()
    sys.exit()
  
