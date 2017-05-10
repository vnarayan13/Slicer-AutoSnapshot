from __future__ import print_function
from __main__ import vtk, qt, ctk, slicer, os

import vtkITK
import numpy
import SimpleITK as sitk
import sitkUtils as su
import vtk.util.numpy_support
import collections
import itertools
import math
import string
import operator

import fnmatch
import pdb
import radiomicsDatabase


def Execute(self, inputdir, outputdir, fileIdentifierSettings, applyButton):
  from datetime import datetime
  progress_file = str(datetime.now().strftime(('%Y-%m-%d_%H-%M-%S')) + '_ImageSnapshot.txt')
  progress_filename = os.path.join(outputdir, progress_file)
  
  with open(progress_filename,mode='w') as printfile: printfile.write(str(fileIdentifierSettings[0]) + '\n')
  
  completed = [item.split('_')[0] for item in os.listdir(outputdir)]
  PatientDirs, PatientNames = radiomicsDatabase.getFolderList(self,inputdir,excludeDirName=False)
  for x,patientName in enumerate(PatientNames):
    if patientName in completed: continue
    print(x,patientName,PatientDirs[x])
    
    #get button widget by reference and change text to x value
    applyButton.text = ("Processing: " + str(x+1) + " of " + str(len(PatientNames)) + ' - ' + patientName)
    applyButton.repaint()
    slicer.app.processEvents()
    
    # Set ID current patient
    IDcurrPatient = radiomicsDatabase.getIDcurrPatient(self,2,PatientDirs[x])
    
    if(os.path.isfile(os.path.join(outputdir,IDcurrPatient + '_axial_sliceZ.png')) or 
       os.path.exists(os.path.join(outputdir,IDcurrPatient+'_flagged'))):
      radiomicsDatabase.deleteDataSlicer(self)
      slicer.mrmlScene.Clear(0)
      continue
    
    #if not os.path.exists(str(os.path.join(outputdir, IDcurrPatient))):
    with open(progress_filename,mode='a') as printfile: 
      printfile.write("\n" + str(x) + "- Loading Data: " + str(IDcurrPatient) + " at location: " + str(PatientDirs[x]) + '\n')
    
    ExtractionDirs = []
    includeSubfolders = True
    subfoldersIndependent = True
    # Normal or Treat each subfolder as independent scans (option GUI)
    if subfoldersIndependent: ExtractionDirs = radiomicsDatabase.getFolderList(self,PatientDirs[x])[0].values()
    else: ExtractionDirs.append(PatientDirs[x])

    # Load and Extract features for each selected folder
    for index in xrange(radiomicsDatabase.lenghtList(ExtractionDirs)):
      for settings in fileIdentifierSettings:
        print(patientName, ExtractionDirs[index])
        ImageFilePath, LabelFilePath = radiomicsDatabase.getDataFiles(ExtractionDirs[index],includeSubfolders,settings["Mask"],settings["selImage"],settings["selLabel"])
        if ImageFilePath and LabelFilePath:
          # Load Image Nodes into Slicer
          with open(progress_filename,mode='a') as printfile: 
            imageNode, labelNode = radiomicsDatabase.loadDataSlicer(self,ImageFilePath,LabelFilePath,printfile)
          maxZ=ImageModelSnapshot(self, imageNode, labelNode, outputdir, IDcurrPatient, ImageFilePath, progress_filename)   
          with open(progress_filename,mode='a') as printfile: 
            printfile.write('  Z index: ' + str(maxZ) + '\n')
        else: 
          with open(progress_filename,mode='a') as printfile: 
            printfile.write("ERROR: Could not find image or label path: " + patientName + '\n')
        radiomicsDatabase.deleteDataSlicer(self)
        slicer.mrmlScene.Clear(0)
    radiomicsDatabase.deleteDataSlicer(self)
    slicer.mrmlScene.Clear(0)
    
def ImageModelSnapshot(self, dataNode, labelNode, outputDir, IDcurrPatient, ImageFilePath, progress_filename, selLabels=False, levels=False):   
  cast1 = vtk.vtkImageCast()
  cast1.SetOutputScalarTypeToInt()
  cast1.SetInputConnection(labelNode.GetImageDataConnection())
  cast1.Update()
  labelNode.SetImageDataConnection(cast1.GetOutputPort())
  
  cast2 = vtk.vtkImageCast()
  cast2.SetOutputScalarTypeToInt()
  cast2.SetInputConnection(dataNode.GetImageDataConnection())
  cast2.Update()
  dataNode.SetImageDataConnection(cast2.GetOutputPort())
  
  labelNode.SetOrigin(dataNode.GetOrigin()) 
  try:
    imageNodeSITK = su.PullFromSlicer(dataNode.GetName())
    labelNodeSITK = su.PullFromSlicer(labelNode.GetName())
    
    labelArraySITK = sitk.GetArrayFromImage(labelNodeSITK)
    lbound, hbound = squeeze(self, labelArraySITK, labelNode)

    self.cmaj = sitk.CropImageFilter()
    self.cmaj.SetLowerBoundaryCropSize(lbound)
    self.cmaj.SetUpperBoundaryCropSize(hbound)

    zoomLabelNodeSITK = self.cmaj.Execute(labelNodeSITK)
    zoomImageNodeSITK = self.cmaj.Execute(imageNodeSITK)

    zoomLabelNodeName = 'zoomLabelNode'+'_'+IDcurrPatient
    zoomImageNodeName = 'zoomImageNode'+'_'+IDcurrPatient
    su.PushToSlicer(zoomLabelNodeSITK, zoomLabelNodeName)
    su.PushToSlicer(zoomImageNodeSITK, zoomImageNodeName)  
    zoomLabelNode = slicer.util.getNode(zoomLabelNodeName)
    zoomImageNode = slicer.util.getNode(zoomImageNodeName)  
    
    zoomLabelNodeArraySITK = sitk.GetArrayFromImage(zoomLabelNodeSITK)
    zoomLabelNodeArraySITK = numpy.where(zoomLabelNodeArraySITK == 0, numpy.NAN, zoomLabelNodeArraySITK)
    maxZind = find_LargestSlices(self, zoomLabelNodeArraySITK)
    
    imgOutputDir = str(outputDir)
    modelHandler(ImageFilePath, zoomImageNode, zoomLabelNode, imgOutputDir, IDcurrPatient, maxZind)

    return(maxZind)
    
  except Exception, e:
    with open(progress_filename,mode='a') as printfile: printfile.write('ERROR: ' + str(e) + '\n')
    os.makedirs(os.path.join(outputDir,IDcurrPatient+'_flagged'))
    return(-1,-1,-1)

def find_LargestSlices(self, matrix):
  maxZ = matrix[0]
  sizeMaxZ = maxZ[~numpy.isnan(maxZ)].size
  maxZind = 0
  for zind, z in enumerate(matrix):
    if z[~numpy.isnan(z)].size > sizeMaxZ: 
      maxZind = zind
      maxZ = z
      sizeMaxZ = maxZ[~numpy.isnan(maxZ)].size
  return (maxZind)
  
def squeeze(self, matrix, volume):
  zmin = 0
  zmax = matrix.shape[0]
  zmin = minfinder(matrix)
  zmax = maxfinder(matrix)

  Xmat = numpy.rollaxis(matrix,2)
  xmin = 0
  xmax = Xmat.shape[0]
  xmin = minfinder( Xmat )
  xmax = maxfinder( Xmat )
  
  Ymat = numpy.rollaxis(matrix,1)
  ymin = 0
  ymax = Ymat.shape[0]
  ymin = minfinder( Ymat)
  ymax = maxfinder( Ymat )
  
  cube = (200.00,200.00,200.00) # lung
  # use (100.00,100.00,100.00) pad cube for brain tumors
  
  dims = tuple(map( lambda x: x-1, volume.GetImageData().GetDimensions() ))
  spacing = volume.GetSpacing()
  
  minCoordinates = (xmin, ymin, zmin)
  maxCoordinates = (xmax, ymax, zmax)

  minCoordinates, maxCoordinates = padXYZ(dims, spacing, minCoordinates, maxCoordinates, cube=cube)
  lbound = minCoordinates
  hbound = tuple(map(lambda (x,y): x-y, zip(dims,maxCoordinates)))

  return lbound, hbound
  
def minfinder(array):
  ind_slice = ( enumerate(iter(array[:-1:])) )
  
  index,slice = next(ind_slice)
  while True:
    try:
      if numpy.all(slice==0):
        index, slice = next(ind_slice)
        if numpy.any(slice!=0):
          return (index)
    except:
      return 0

def maxfinder(array):
  ind_slice = ( ((int(array.shape[0]) - index - 1),slice) for index,slice in enumerate(iter(array[-1:0:-1])) )
  
  index,slice = next(ind_slice)
  while True:
    try:
      if numpy.all(slice==0):
        index, slice = next(ind_slice)
        if numpy.any(slice!=0):
          return (index)
    except:
      return array.shape[0]
      
def padXYZ(dims, spacing, minCoordinates, maxCoordinates, cube=(200.0,200.0,200.0)):
  padFx = lambda (minCoordinates,maxCoordinates,cube,spacing): int(numpy.floor( ((cube - spacing*(maxCoordinates - minCoordinates)) / spacing) / 2.0 ))
  pad = tuple(map(padFx, zip(minCoordinates,maxCoordinates,cube,spacing)))
  
  minCoordinatesPad = tuple(map(lambda (min,pad): min-pad if min-pad > 0 else 0, zip(minCoordinates, pad)))
  maxCoordinatesPad = tuple(map(lambda (max,pad,dims): max+pad if max+pad < dims else dims, zip(maxCoordinates, pad, dims)))
  
  return minCoordinatesPad, maxCoordinatesPad
  

def modelHandler(ImageFilePath, imageNode, labelNode, imgOutputDir, IDcurrPatient, Zind):
  """
  imageBuffer = vtk.vtkImageData()
  label = growCut(imageNode.GetImageData(),labelNode.GetImageData(),imageBuffer)
  q = labelNode.GetImageData()
  q.DeepCopy(label)
  if labelNode.GetImageDataConnection():
    labelNode.GetImageDataConnection().GetProducer().Update()
  if labelNode.GetImageData().GetPointData().GetScalars() != None:
    labelNode.GetImageData().GetPointData().GetScalars().Modified()
  labelNode.GetImageData().Modified()
  labelNode.Modified()
  pdb.set_trace()
  """
  imageSeriesDescription = os.path.basename(ImageFilePath).replace('.nrrd','')
  imagePatientID = os.path.basename(os.path.dirname(os.path.dirname(os.path.dirname((ImageFilePath)))))
  imageStudyDate = os.path.basename(os.path.dirname(os.path.dirname((ImageFilePath))))
  imagePatientID_StudyDate = imagePatientID + '_' + imageStudyDate
  
  inputImage = su.PullFromSlicer(imageNode.GetName())
  inputLabel = su.PullFromSlicer(labelNode.GetName())
  pixelID = inputImage.GetPixelIDValue()
  FlipYAxis = sitk.FlipImageFilter()
  FlipYAxis.SetFlipAxes([False, True, False])  
  sizeZ = list(inputImage.GetSize())
  sizeZ[2] = 0
  indexZ = [0, 0, Zind]
  ExtractorZ = sitk.ExtractImageFilter()
  ExtractorZ.SetSize(sizeZ)
  ExtractorZ.SetIndex(indexZ)
  sliceZImage = ExtractorZ.Execute(inputImage)
  sliceZLabel = ExtractorZ.Execute(inputLabel)
  imageName = 'Axial_Slice_Image' + IDcurrPatient
  labelName = 'Axial_Slice_Label' + IDcurrPatient
  su.PushToSlicer(sliceZImage, imageName)
  su.PushToSlicer(sliceZLabel, labelName)
  
  sliceZImageNode = slicer.util.getNode(imageName)
  sliceZLabelNode = slicer.util.getNode(labelName)
  volumesLogic = slicer.vtkSlicerVolumesLogic()
  volumesLogic.SetVolumeAsLabelMap(sliceZLabelNode, True)
  volumesLogic.SetVolumeAsLabelMap(labelNode, True)
   
  #sliceZLabelNode = binarizeLabelMapToValue(sliceZLabelNode, labelValue=296)
  #labelNode = binarizeLabelMapToValue(labelNode, labelValue=296)
  
  sliceZLabelNodeDisplay = sliceZLabelNode.GetDisplayNode()
  sliceZLabelNodeDisplay.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileGenericColors.txt')
  sliceZLabelNode.SetAndObserveDisplayNodeID(sliceZLabelNodeDisplay.GetID())
  labelNodeDisplay = labelNode.GetDisplayNode()
  labelNodeDisplay.SetAndObserveColorNodeID('vtkMRMLColorTableNodeFileGenericColors.txt')
  labelNode.SetAndObserveDisplayNodeID(labelNodeDisplay.GetID())
  
  ##
  appLogic = slicer.app.applicationLogic()
  selectionNode = appLogic.GetSelectionNode()
  selectionNode.SetReferenceActiveVolumeID(sliceZImageNode.GetID())
  selectionNode.SetReferenceActiveLabelVolumeID(sliceZLabelNode.GetID())  
  appLogic.PropagateVolumeSelection()
  
  ##
  lm = slicer.app.layoutManager()
  lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)
  redWidget = lm.sliceWidget('Red')
  redLogic = redWidget.sliceLogic()
  redView = redWidget.sliceView()
  #redLogic.SetBackgroundWindowLevel(*windowLevel)
   
  sln = slicer.util.getNode('vtkMRMLSliceNodeRed')
  dims = list(sliceZImageNode.GetImageData().GetDimensions())
  # dims[0] is x, dims[1] is y, dims [2] is Z 
  redWidget.setFixedSize(720,660)
  slncw = redWidget.sliceController()
  slncw.showLabelOutline(1)
  slncw.fitSliceToBackground()
  #sln.SetFieldOfView(dims[0],dims[1],1)
  #sln.SetDimensions(dims[0],dims[1],1)
  
  sliceannotations = slicer.modules.DataProbeInstance.infoWidget.sliceAnnoations
  if sliceannotations.sliceViewAnnotationsCheckBox.checkState() == 2:
    sliceannotations.sliceViewAnnotationsCheckBox.click()
  
  slicer.app.processEvents()
  wti=vtk.vtkWindowToImageFilter()  
  wti.SetInput(redView.renderWindow())
  wti.Update()
  imgDataRedSlice = wti.GetOutput()
 
  modelMakerCLI(labelNode) 
  imgData3D = GetModelSnapshot()
  
  append = vtk.vtkImageAppend() 
  append.SetAppendAxis(0)
  append.AddInputData(imgDataRedSlice)
  append.AddInputData(imgData3D)
  append.Update()
  finalImage = append.GetOutput()
  
  #finalImagePath = str(os.path.join(imgOutputDir, IDcurrPatient + '_axial_sliceZ.png'))
  finalImagePath = str(os.path.join(imgOutputDir, imagePatientID_StudyDate + '_' + imageSeriesDescription +'.png'))
  writer=vtk.vtkPNGWriter()
  writer.SetFileName(finalImagePath)
  writer.SetInputData(finalImage)
  writer.Write()
  
def GetModelSnapshot():
  lm = slicer.app.layoutManager()
  lm.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutFourUpView)
  view = lm.threeDWidget(0).threeDView()
  viewNode = view.mrmlViewNode()
  
  view.resetFocalPoint()
  viewNode.SetBackgroundColor((0,0,0))
  viewNode.SetBackgroundColor2((0,0,0))
  viewNode.SetAxisLabelsVisible(0)
  viewNode.SetBoxVisible(0)
  view.lookFromViewAxis(4)
  view.setPitchRollYawIncrement(90)
  view.pitch()
  view.setPitchRollYawIncrement(30)
  view.yaw()
  view.pitchDirection = 1
  view.pitch()
  view.setZoomFactor(10)
  view.zoomIn()
  #view.setZoomFactor(5)
  #view.zoomIn()
  #view.zoomIn()
  view.forceRender()
   
  rw=view.renderWindow()
  wti=vtk.vtkWindowToImageFilter()
  wti.SetInput(rw)
  slicer.app.processEvents()
  wti.Update() 
  imgData3D = wti.GetOutput()
  return(imgData3D)
  
def modelMakerCLI(inputvol):
  modelHNode = slicer.mrmlScene.CreateNodeByClass('vtkMRMLModelHierarchyNode')
  modelHNode.SetName('TumorModel')
  modelHNode = slicer.mrmlScene.AddNode(modelHNode)
  
  parameters = {}
  parameters["InputVolume"] = inputvol
  parameters['ModelSceneFile'] = modelHNode.GetID()
  parameters['Smooth'] = 8
  modelmaker = slicer.modules.modelmaker 
  return (slicer.cli.run(modelmaker, None, parameters, wait_for_completion = True))  

def binarizeLabelMapToValue(labelNode, labelValue=1):
  labelNodeImageData = labelNode.GetImageData()
  change = slicer.vtkImageLabelChange()
  change.SetInputData(labelNodeImageData)
  change.SetOutput(labelNodeImageData)
  change.SetOutputLabel(labelValue)
  
  for i in xrange(1,int(labelNodeImageData.GetScalarRange()[1])+1):
    change.SetInputLabel(i)
    change.Update()
    
  labelNode.SetAndObserveImageData(labelNodeImageData)
  return labelNode


def growCut(background,gestureInput,growCutOutput):
  growCutFilter = vtkITK.vtkITKGrowCutSegmentationImageFilter()

  # set the make a zero-valued volume for the output
  # TODO: maybe this should be done in numpy as a one-liner
  thresh = vtk.vtkImageThreshold()
  thresh.ReplaceInOn()
  thresh.ReplaceOutOn()
  thresh.SetInValue(0)
  thresh.SetOutValue(0)
  thresh.SetOutputScalarType( vtk.VTK_SHORT )

  thresh.SetInputData( gestureInput )
  
  thresh.SetOutput( growCutOutput )
  thresh.Update()
  growCutOutput.DeepCopy( gestureInput )

  growCutFilter.SetInputData( 0, background )
  growCutFilter.SetInputData( 1, gestureInput )
  growCutFilter.SetInputConnection( 2, thresh.GetOutputPort() )

  objectSize = 10. # TODO: this is a magic number
  contrastNoiseRatio = 0.8 # TODO: this is a magic number
  priorStrength = 0.003 # TODO: this is a magic number
  segmented = 2 # TODO: this is a magic number
  conversion = 1000 # TODO: this is a magic number

  spacing = gestureInput.GetSpacing()
  voxelVolume = reduce(lambda x,y: x*y, spacing)
  voxelAmount = objectSize / voxelVolume
  voxelNumber = round(voxelAmount) * conversion

  cubeRoot = 1./3.
  oSize = int(round(pow(voxelNumber,cubeRoot)))

  growCutFilter.SetObjectSize( oSize )
  growCutFilter.SetContrastNoiseRatio( contrastNoiseRatio )
  growCutFilter.SetPriorSegmentConfidence( priorStrength )
  growCutFilter.Update()

  growCutOutput.DeepCopy( growCutFilter.GetOutput() )

  return growCutOutput
  
def castVolumeCLI(self, InputVolume, OutputVolume, Type='Short'):
  parameters = {}
  parameters["InputVolume"] = InputVolume
  parameters["OutputVolume"] = OutputVolume
  parameters["Type"] = Type
  castVolume = slicer.modules.castscalarvolume
  return (slicer.cli.run(castVolume, None, parameters, wait_for_completion = True))