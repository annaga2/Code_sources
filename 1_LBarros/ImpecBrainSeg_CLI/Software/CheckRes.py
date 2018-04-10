import numpy as np
import sys
import os
import os.path as op
import SimpleITK as sitk


#1) get intputs
VoxCoordX=int(sys.argv[1])
VoxCoordY=int(sys.argv[2])
VoxCoordZ=int(sys.argv[3])

PrefixOutputs='out_'

if len(sys.argv)>4:
  PrefixOutputs=sys.argv[4]



#2) get the tissu type

print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+           Probabilistic segmentation of brain structures             +")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
for tissuType in ['gm_cortex','wm','csf','gm_subco','Cerebelum','BrainStem','OlfactoryBulb','hippo']:
  SegmentationFile=op.join('.',PrefixOutputs+'BrainStructures',PrefixOutputs+tissuType+'.nii')
  SegArray=sitk.GetArrayFromImage(sitk.ReadImage(SegmentationFile))
  LocProba=SegArray[VoxCoordZ-1,VoxCoordY-1,VoxCoordX-1]
  if LocProba>0.01:
    print(tissuType+': '+str(LocProba))

#3) get the brain area 

print("")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+            Probabilistic segmentation of cortical areas              +")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")

#3.1) load region labels information

RegionLabelsFile=op.join('.',PrefixOutputs+'CorticalRegions','CortialRegions_Labels.txt')

List_RegionLabels=[]

f = open(RegionLabelsFile, 'r')
sep=0
for line in f:
    line2=line.replace('\n', '')
    if line2=="#+++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++":
      sep+=1
    words=line2.split(' ')
    if sep==2 and len(words)>1:
      List_RegionLabels.append(words[1])
f.close()
  

#3.2) get the probabilities

for i in range(1,116):
    if i<10:
      zeros='00'
    elif i<100:
      zeros='0'
    else:
      zeros=''
    
    #3.2.1) load proba map L+str(i)
    LocProbaFile=op.join('.',PrefixOutputs+'CorticalRegions','Reg_L'+zeros+str(i)+'_RSP.nii.gz')
    SegArray=sitk.GetArrayFromImage(sitk.ReadImage(LocProbaFile))
    LocProba=SegArray[VoxCoordZ-1,VoxCoordY-1,VoxCoordX-1]
    if LocProba>0.01:
      print(List_RegionLabels[i]+' L (L'+zeros+str(i)+'): '+str(LocProba))
    
    #3.2.2) load proba map R+str(i)
    LocProbaFile=op.join('.',PrefixOutputs+'CorticalRegions','Reg_R'+zeros+str(i)+'_RSP.nii.gz')
    SegArray=sitk.GetArrayFromImage(sitk.ReadImage(LocProbaFile))
    LocProba=SegArray[VoxCoordZ-1,VoxCoordY-1,VoxCoordX-1]
    if LocProba>0.01:
      print(List_RegionLabels[i]+' R (R'+zeros+str(i)+'): '+str(LocProba))

#4) get the brain regions 

print("")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")
print("+           Probabilistic segmentation of cortical regions             +")
print("++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++++")


for i in range(1,8):
    
      
    #4.1) group location
    if i==1:
      groupLoc="Frontal -- pre-fontal"
    if i==2:
      groupLoc="Frontal -- motor"
    if i==3:
      groupLoc="Limbic"
    if i==4:
      groupLoc="Temporal"
    if i==5:
      groupLoc="Parietal"
    if i==6:
      groupLoc="Insular"
    if i==7:
      groupLoc="Occipital"
    
    #4.2) load proba map L+str(i)
    LocProbaFile=op.join('.',PrefixOutputs+'CorticalRegions','Reg_GroupL'+str(i)+'_RSP.nii.gz')
    SegArray=sitk.GetArrayFromImage(sitk.ReadImage(LocProbaFile))
    LocProba=SegArray[VoxCoordZ-1,VoxCoordY-1,VoxCoordX-1]
    if LocProba>0.01:
      print(groupLoc+' L: '+str(LocProba))
    
    #4.3) load proba map R+str(i)
    LocProbaFile=op.join('.',PrefixOutputs+'CorticalRegions','Reg_GroupR'+str(i)+'_RSP.nii.gz')
    SegArray=sitk.GetArrayFromImage(sitk.ReadImage(LocProbaFile))
    LocProba=SegArray[VoxCoordZ-1,VoxCoordY-1,VoxCoordX-1]
    if LocProba>0.01:
      print(groupLoc+' R: '+str(LocProba))





