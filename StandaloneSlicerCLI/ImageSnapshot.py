from __main__ import vtk, qt, ctk, slicer

import os
import collections
import itertools
import glob
import operator
import fnmatch
import pdb

from datetime import datetime

#import vtkITK
import numpy as np
import SimpleITK as sitk
import sitkUtils as su
#import vtk.util.numpy_support

### USAGE ###
# "C:\Program Files\Slicer 4.4.0\Slicer.exe" --python-script "ImageSnapshot.py"


def main(): 
    inputDirectory = r"C:\Users\Vivek Narayan\Desktop\batch_script\02_SlicerData"
    outputDirectory = r"C:\Users\Vivek Narayan\Desktop\batch_script\test_output"
    filetype = ".nrrd"
    
    keywordSettings = {} 
    keywordSettings['image'] = ""
    keywordSettings['imageExclusion'] = "label"
    keywordSettings['label'] = "label"
    keywordSettings['labelExclusion'] = ""
    
    if not os.path.exists(outputDirectory): os.mkdir(outputDirectory)
    progress_file = str(datetime.now().strftime(('%Y-%m-%d_%H-%M-%S')) + '_ImageSnapshotLog.txt')
    progress_filename = os.path.join(outputDirectory, progress_file)
    with open(progress_filename,mode='w') as printfile: printfile.write(str(keywordSettings.items()) + '\n')
           
    dataHierarchyDict = ReadDatasetHierarchy(inputDirectory, imagetype=filetype)    
    #completed = [item.split('_')[0] for item in os.listdir(outputDirectory)]
    
    for patientIndex, patientDirectory in enumerate(dataHierarchyDict):
        patientID = os.path.basename(patientDirectory)
        
        # skip completed
        #if patientID in completed: continue
        
        for studyDirectory in dataHierarchyDict[patientDirectory]:
            studyDate = os.path.basename(studyDirectory)

            with open(progress_filename,mode='a') as printfile: 
                printfile.write("(%s/%s) Processing Patient: %s, Study: %s" %(str(patientIndex+1), str(len(dataHierarchyDict.keys())), patientID, studyDate) + '\n')
                print "(%s/%s) Processing Patient: %s, Study: %s" %(str(patientIndex+1), str(len(dataHierarchyDict.keys())), patientID, studyDate)
                
            imageFilepaths = dataHierarchyDict[patientDirectory][studyDirectory]["reconstructions"]
            labelFilepaths = dataHierarchyDict[patientDirectory][studyDirectory]["segmentations"]
            resourceFilepaths = dataHierarchyDict[patientDirectory][studyDirectory]["resources"]
            
            imageFilepath, labelFilepath = getImageAndLabelPair(imageFilepaths, labelFilepaths, keywordSettings)

            if (imageFilepath is not None) and (labelFilepath is not None):
                with open(progress_filename,mode='a') as printfile: 
                    imageNode, labelNode = loadDataSlicer(imageFilepath, labelFilepath, printfile)
                
                outputImagePath = str(os.path.join(outputDirectory, patientID + '_' + studyDate + '_' + os.path.basename(labelFilepath).replace(filetype,'') + '.png'))
                maxZ = ImageModelSnapshot(imageNode, labelNode, outputImagePath, progress_filename)   
                
                with open(progress_filename,mode='a') as printfile: 
                    printfile.write('\tMax Surface Area Z-index: ' + str(maxZ) + '\n')
            
            slicer.mrmlScene.Clear(0)

def ImageModelSnapshot(imageNode, labelNode, outputImagePath, progress_filename):   
    castScalarToInt(imageNode)
    castScalarToInt(labelNode)
     
    # center label to origin of image, if needed 
    #labelNode.SetOrigin(imageNode.GetOrigin()) 
    
    try:
        minBounds, maxBounds = getCropBoundaries(labelNode)
        
        zoomImageNode = cropToZoom(imageNode, minBounds, maxBounds, label=False)
        zoomLabelNode = cropToZoom(labelNode, minBounds, maxBounds, label=True)
        
        sitkZoomLabelNode = su.PullFromSlicer(zoomLabelNode.GetName())
        zoomLabelNodeArray = sitk.GetArrayFromImage(sitkZoomLabelNode)
        zoomLabelNodeArray = np.where(zoomLabelNodeArray == 0, np.NAN, zoomLabelNodeArray)
        
        # Find index of Z axis slice with largest surface area in the label node
        maxZind = findLargestSurfaceAreaSlice(zoomLabelNodeArray)
        
        outputImage = modelHandler(zoomImageNode, zoomLabelNode, maxZind)
        
        writer=vtk.vtkPNGWriter()
        writer.SetFileName(outputImagePath)
        writer.SetInputData(outputImage)
        writer.Write()
        
        return(maxZind)
      
    except Exception, e:
        with open(progress_filename,mode='a') as printfile: printfile.write('ERROR: ' + str(e) + '\n\n')
        print "ERROR: Error Creating", outputImagePath, str(e) 
        return(-1,-1,-1)
     
def castScalarToInt(volumeNode):
    castScalarToInt = vtk.vtkImageCast()
    castScalarToInt.SetOutputScalarTypeToInt()
    
    castScalarToInt.SetInputConnection(volumeNode.GetImageDataConnection())
    castScalarToInt.Update()
    volumeNode.SetImageDataConnection(castScalarToInt.GetOutputPort())
    
def cropToZoom(volumeNode, minBounds, maxBounds, label=False):     
    # bring image and label node to sitk, crop both to tumor region (for zooming in)
    # send cropped image and label back to slicer
    cmaj = sitk.CropImageFilter()
    cmaj.SetLowerBoundaryCropSize(minBounds)
    cmaj.SetUpperBoundaryCropSize(maxBounds)
    
    sitkVolumeNode = su.PullFromSlicer(volumeNode.GetName())
    sitkZoomVolumeNode = cmaj.Execute(sitkVolumeNode)
  
    zoomVolumeNodeName = volumeNode.GetName() + '_' + 'zoomed'   
    if label: su.PushLabel(sitkZoomVolumeNode, zoomVolumeNodeName)
    else: su.PushBackground(sitkZoomVolumeNode, zoomVolumeNodeName)    
    zoomVolumeNode = slicer.util.getNode(zoomVolumeNodeName)
    
    return zoomVolumeNode
    
def getCropBoundaries(labelNode):
    sitkLabelNode = su.PullFromSlicer(labelNode.GetName())
    labelArray = sitk.GetArrayFromImage(sitkLabelNode)
    
    zmin = 0
    zmax = labelArray.shape[0]
    zmin = minfinder(labelArray)
    zmax = maxfinder(labelArray)

    Xmat = np.rollaxis(labelArray,2)
    xmin = 0
    xmax = Xmat.shape[0]
    xmin = minfinder( Xmat )
    xmax = maxfinder( Xmat )
    
    Ymat = np.rollaxis(labelArray,1)
    ymin = 0
    ymax = Ymat.shape[0]
    ymin = minfinder( Ymat)
    ymax = maxfinder( Ymat )
    
    cube = (200.00,200.00,200.00) # lung
    # use (100.00,100.00,100.00) pad cube for brain tumors
    
    dims = tuple(map( lambda x: x-1, labelNode.GetImageData().GetDimensions() ))
    spacing = labelNode.GetSpacing()
    
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
            if np.all(slice==0):
                index, slice = next(ind_slice)
                if np.any(slice!=0):
                    return (index)
            else:
                return(index)          
        except:
            return 0

def maxfinder(array):
    ind_slice = ( ((int(array.shape[0]) - index - 1),slice) for index,slice in enumerate(iter(array[-1:0:-1])) )
    
    index,slice = next(ind_slice)
    while True:
        try:
            if np.all(slice==0):
                index, slice = next(ind_slice)
                if np.any(slice!=0):
                    return (index)
            else:
                return(index)      
        except:
            return array.shape[0]
      
def padXYZ(dims, spacing, minCoordinates, maxCoordinates, cube=(200.0,200.0,200.0)):
    padFx = lambda (minCoordinates,maxCoordinates,cube,spacing): int(np.floor( ((cube - spacing*(maxCoordinates - minCoordinates)) / spacing) / 2.0 ))
    pad = tuple(map(padFx, zip(minCoordinates,maxCoordinates,cube,spacing)))
    
    minCoordinatesPad = tuple(map(lambda (min,pad): min-pad if min-pad > 0 else 0, zip(minCoordinates, pad)))
    maxCoordinatesPad = tuple(map(lambda (max,pad,dims): max+pad if max+pad < dims else dims, zip(maxCoordinates, pad, dims)))
    
    return minCoordinatesPad, maxCoordinatesPad

def findLargestSurfaceAreaSlice(array):
    maxZ = array[0]
    sizeMaxZ = maxZ[~np.isnan(maxZ)].size
    maxZind = 0
    for zind, z in enumerate(array):
        if z[~np.isnan(z)].size > sizeMaxZ: 
            maxZind = zind
            maxZ = z
            sizeMaxZ = maxZ[~np.isnan(maxZ)].size
    return (maxZind)
  
def modelHandler(imageNode, labelNode, Zind): 
    sliceZImageNode = extractMaxZSlice(imageNode, Zind)
    sliceZLabelNode = extractMaxZSlice(labelNode, Zind) 

    volumesLogic = slicer.vtkSlicerVolumesLogic()
    volumesLogic.SetVolumeAsLabelMap(sliceZLabelNode, True)
    volumesLogic.SetVolumeAsLabelMap(labelNode, True)
    
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
    
    return finalImage

def extractMaxZSlice(volumeNode, Zind):
    sitkVolumeNode = su.PullFromSlicer(volumeNode.GetName())
    
    pixelID = sitkVolumeNode.GetPixelIDValue()
    FlipYAxis = sitk.FlipImageFilter()
    FlipYAxis.SetFlipAxes([False, True, False])  
    sizeZ = list(sitkVolumeNode.GetSize())
    sizeZ[2] = 0
    indexZ = [0, 0, Zind]
    ExtractorZ = sitk.ExtractImageFilter()
    ExtractorZ.SetSize(sizeZ)
    ExtractorZ.SetIndex(indexZ)
    
    sitkZSlice = ExtractorZ.Execute(sitkVolumeNode)
    ZSliceName = 'maxz_' + volumeNode.GetName()
    su.PushToSlicer(sitkZSlice, ZSliceName) 
    ZSliceNode = slicer.util.getNode(ZSliceName)
    return ZSliceNode
    
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
    
def loadDataSlicer(imageFilepath, labelFilepath, printfile):
    # loads Image and Label into Slicer
    
    image = list(slicer.util.loadVolume(imageFilepath, returnNode=True))
    if not (image[0]==True): 
        printfile.write("Error loading Image data\n")
        image = None    
    else:
        printfile.write("Image file loaded " + str(imageFilepath) + '\n')
        image = image[1]
        
    label = list(slicer.util.loadVolume(labelFilepath, properties={'labelmap':True}, returnNode=True))
    if not (label[0]==True): 
        printfile.write("Error loading Label data\n")
        label = None
    else:
        printfile.write("Label file loaded " + str(labelFilepath) + '\n')
        label = label[1]
        
    return image, label    
    
def getImageAndLabelPair(imageFilepaths, labelFilepaths, keywordSettings):
    keywordSettings = {k:[str(keyword.strip()) for keyword in v.split(',')] for (k,v) in keywordSettings.iteritems()}
    
    matchedImages = []
    for imageFilepath in imageFilepaths:
        imageFilename = str(os.path.basename(imageFilepath))
        if testString(imageFilename, keywordSettings['image'], keywordSettings['imageExclusion']): 
            matchedImages.append(imageFilepath)
     
    matchedLabels = []
    for labelFilepath in labelFilepaths:
        labelFilename = str(os.path.basename(labelFilepath))
        if testString(labelFilename, keywordSettings['label'], keywordSettings['labelExclusion']):
            matchedLabels.append(labelFilepath)
            
    if len(matchedImages) < 1: print "ERROR: No Images Matched"
    elif len(matchedImages) > 1: print "ERROR: Multiple Images Matched"
    
    if len(matchedLabels) < 1: print "ERROR: No Labels Matched"
    elif len(matchedLabels) > 1: print "ERROR: Multiple Labels Matched"
    
    if (len(matchedImages) == 1) and (len(matchedLabels) == 1):
        return matchedImages[0], matchedLabels[0]
    else:
        return None, None
    
def testString(fileName, inclusionKeywords, exclusionKeywords):
    fileName = fileName.upper()
    inclusionKeywords = [keyword.upper() for keyword in inclusionKeywords if (keyword != '')]
    exclusionKeywords = [keyword.upper() for keyword in exclusionKeywords if (keyword != '')]
    
    result = False
    if (len(inclusionKeywords) == 0) and (len(exclusionKeywords) > 0):
        if (not any(keyword in fileName for keyword in exclusionKeywords)):
            result = True
    elif (len(inclusionKeywords) > 0) and (len(exclusionKeywords) == 0):
        if (all(keyword in fileName for keyword in inclusionKeywords)):
            result = True        
    elif (len(inclusionKeywords) > 0) and (len(exclusionKeywords) > 0):
        if (all(keyword in fileName for keyword in inclusionKeywords)) and \
          (not any(keyword in fileName for keyword in exclusionKeywords)):
            result = True 
    elif (len(inclusionKeywords) == 0) and (len(exclusionKeywords) == 0):
        result = True 
    
    return result    
      
def ReadReconstructionsDir(studyDir, subfolders, imagetype, create=False):
    images = []
    recDir = "NONE"
    try:
        recDir = [item for item in subfolders if 'reconstructions' in os.path.basename(item).lower()][0]
        images = [item for item in glob.glob(os.path.join(recDir,"*")) if imagetype in os.path.basename(item)]
    except IndexError:
        if create:
            recDir = os.path.join(studyDir,"Reconstructions")
            if not os.path.exists(recDir):
                os.mkdir(recDir)
                print "\tCreated:", recDir
                
    return recDir, images
            
def ReadSegmentationsDir(studyDir, subfolders, imagetype, create=False):
    labels = []
    segDir = "NONE"
    try:
        segDir = [item for item in subfolders if 'segmentations' in os.path.basename(item).lower()][0]
        labels = [item for item in glob.glob(os.path.join(segDir,"*")) if imagetype in os.path.basename(item)]
    except IndexError:
        if create:
            segDir = os.path.join(studyDir,"Segmentations")
            if not os.path.exists(segDir):
                os.mkdir(segDir)
                print "\tCreated:", segDir

    return segDir, labels
    
def ReadResourcesDir(studyDir, subfolders, create=False):
    resources = []
    resDir = "NONE"
    try:
        resDir = [item for item in subfolders if 'resources' in os.path.basename(item).lower()][0]
        resources = [item for item in glob.glob(os.path.join(resDir,"*"))]
    except IndexError:
        if create:
            resDir = os.path.join(studyDir,"Resources")
            if not os.path.exists(resDir):
                os.mkdir(resDir)
                print "\tCreated:", resDir

    return resDir, resources
        
def ReadDatasetHierarchy(mainPatientDir, imagetype='.nrrd'):
    DatabaseHierarchyDict = collections.OrderedDict()
    
    patientDirs = glob.glob(os.path.join(mainPatientDir,'*'))
    
    for patientDir in patientDirs:
        patientID = os.path.basename(patientDir)
        studydateDirs = glob.glob(os.path.join(patientDir,'*'))
        DatabaseHierarchyDict[patientDir] = {}
        
        for studydateDir in studydateDirs:
            DatabaseHierarchyDict[patientDir][studydateDir] = {}
            studydatename = os.path.basename(studydateDir)
            
            subfolders = [dirpath for dirpath in glob.glob(os.path.join(studydateDir,'*')) if os.path.isdir(dirpath)]
            
            reconstructionsDir, images = ReadReconstructionsDir(studydateDir, subfolders, imagetype, create=False)
            DatabaseHierarchyDict[patientDir][studydateDir][reconstructionsDir] = images
            DatabaseHierarchyDict[patientDir][studydateDir]["reconstructions"] = images
            
            resourcesDir, resources = ReadResourcesDir(studydateDir, subfolders, create=False)
            DatabaseHierarchyDict[patientDir][studydateDir][resourcesDir] = resources
            DatabaseHierarchyDict[patientDir][studydateDir]["resources"] = resources
            
            segmentationsDir, labels = ReadSegmentationsDir(studydateDir, subfolders, imagetype, create=False)
            DatabaseHierarchyDict[patientDir][studydateDir][segmentationsDir] = labels
            DatabaseHierarchyDict[patientDir][studydateDir]["segmentations"] = labels
            
    return DatabaseHierarchyDict
       
main()    