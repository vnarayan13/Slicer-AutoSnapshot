from __main__ import vtk, qt, ctk, slicer, os, datetime

import string
import pdb

import IMGRadiomicsToolsLib

class radiomicsImageSnapshot:
  def __init__(self, parent):  
    parent.title = "RadiomicsImageSnapshot"
    parent.categories = ["Radiomics"]
    parent.contributors = ["Vivek Narayan / Hugo Aerts"]
    parent.helpText = """
    This module takes a screenshot of a label map overlaid on an image at the largest slice.
    A model of the label map is also generated and appended to the screenshot.
    """
    parent.acknowledgementText = """
    """ 
    self.parent = parent

#
# Widget
#

class radiomicsImageSnapshotWidget:
  def __init__(self, parent=None):
    if not parent:
      self.parent = slicer.qMRMLWidget()
      self.parent.setLayout(qt.QVBoxLayout())
      self.parent.setMRMLScene(slicer.mrmlScene)
    else:
      self.parent = parent
    self.layout = self.parent.layout()
    if not parent:
      self.setup()
      self.parent.show()

  def setup(self):
    # intialize variables
    # Data of current loaded patient
    self.IDcurrPatient = {}
    self.imageNodes = {}      # List of currently loaded image nodes
    self.labelNodes = {}      # List of currently loaded label nodes
    self.ModelNodes = {}
    # Data of database
    self.mainPatientdir = self.outputDir = self.datafile = {}

    #---------------------------------------------------------
    # 2D Slice Snapshots
    self.IMGSliceExtractCollapsibleButton = ctk.ctkCollapsibleButton()
    self.IMGSliceExtractCollapsibleButton.text = "2D Slice and Model Image Capture"
    self.layout.addWidget(self.IMGSliceExtractCollapsibleButton)
    self.IMGSliceExtractCollapsibleButton.collapsed = False
    self.IMGSliceExtractFormLayout = qt.QFormLayout(self.IMGSliceExtractCollapsibleButton)

    # Input 1: Input Directory selector
    self.input6Frame = qt.QFrame(self.IMGSliceExtractCollapsibleButton)
    self.input6Frame.setLayout(qt.QHBoxLayout())
    self.IMGSliceExtractFormLayout.addWidget(self.input6Frame)
    
    self.input6Selector = qt.QLabel("Input Directory (DICOM):  ", self.input6Frame)
    self.input6Frame.layout().addWidget(self.input6Selector)
    self.input6Button = qt.QPushButton("Select main directory DICOM files")
    self.input6Button.toolTip = "Select main directory with DICOM files (folder names are patient names)"
    self.input6Button.enabled = True
    self.input6Frame.layout().addWidget(self.input6Button)

    # Input 2: Output Directory selector
    self.input7Frame = qt.QFrame(self.IMGSliceExtractCollapsibleButton)
    self.input7Frame.setLayout(qt.QHBoxLayout())
    self.IMGSliceExtractFormLayout.addWidget(self.input7Frame)
    
    self.input7Selector = qt.QLabel("Output Directory:  ", self.input7Frame)
    self.input7Frame.layout().addWidget(self.input7Selector)
    self.input7Button = qt.QPushButton("Select main directory output NRRD or NFITI files")
    self.input7Button.toolTip = "Select main directory for output NRRD or NIFTI files (folder names are patient names)"
    self.input7Button.enabled = True
    self.input7Frame.layout().addWidget(self.input7Button)
	
    # Keywords Collapsible Button
    self.KeywordsCollapsibleButton = ctk.ctkCollapsibleButton(self.IMGSliceExtractCollapsibleButton)
    self.KeywordsCollapsibleButton.text = "Keyword Matching"
    self.KeywordsCollapsibleButtonLayout = qt.QHBoxLayout()
    self.KeywordsCollapsibleButton.setLayout(self.KeywordsCollapsibleButtonLayout)    
    self.IMGSliceExtractFormLayout.addWidget(self.KeywordsCollapsibleButton)
    self.KeywordsCollapsibleButton.collapsed = False
    
    self.keywordsFrame = qt.QFrame(self.KeywordsCollapsibleButton)
    self.keywordsFrame.setLayout(qt.QFormLayout())
    self.KeywordsCollapsibleButtonLayout.addWidget(self.keywordsFrame)
    
    # Radio Buttons Frame
    self.radioButtonFrame = qt.QFrame(self.keywordsFrame)
    self.radioButtonFrame.setLayout(qt.QFormLayout())
    self.fileFormatGroup = qt.QButtonGroup(self.radioButtonFrame)
    self.nrrdButton = qt.QRadioButton("NRRD")
    self.nrrdButton.checked = True
    self.niftiButton = qt.QRadioButton("NIFTI")
    self.fileFormatGroup.addButton(self.nrrdButton)
    self.fileFormatGroup.addButton(self.niftiButton)
    self.radioButtonFrame.layout().addRow(self.nrrdButton, self.niftiButton)
    self.keywordsHeader = qt.QLabel("Keyword Matching:", self.keywordsFrame)
    self.keywordsFrame.layout().addRow(self.keywordsHeader)
    
    # Keywords Frame
    self.inputMaskHeader = qt.QLabel("Input Image File Type:", self.keywordsFrame)
    self.keywordsFrame.layout().addRow(self.inputMaskHeader, self.radioButtonFrame)
    
    self.inputImageKeywords = qt.QLabel("Input Image Keywords:", self.keywordsFrame)
    self.inputImageKeywordsField = qt.QLineEdit("",self.keywordsFrame)
    self.keywordsFrame.layout().addRow(self.inputImageKeywords, self.inputImageKeywordsField )
    
    self.inputImageExclusionKeywords = qt.QLabel("Input Image Exclusion Keywords:", self.keywordsFrame)
    self.inputImageExclusionKeywordsField = qt.QLineEdit("",self.keywordsFrame)
    self.keywordsFrame.layout().addRow(self.inputImageExclusionKeywords, self.inputImageExclusionKeywordsField )
    
    self.inputLabelKeywords = qt.QLabel("Input Label Keywords:", self.keywordsFrame)
    self.inputLabelKeywordsField = qt.QLineEdit("",self.keywordsFrame)
    self.keywordsFrame.layout().addRow(self.inputLabelKeywords, self.inputLabelKeywordsField )
    
    self.inputLabelExclusionKeywords = qt.QLabel("Input Label Exclusion Keywords:", self.keywordsFrame)
    self.inputLabelExclusionKeywordsField = qt.QLineEdit("",self.keywordsFrame)
    self.keywordsFrame.layout().addRow(self.inputLabelExclusionKeywords, self.inputLabelExclusionKeywordsField )
     
    # IMG Slice Extract Button
    self.IMGSliceExtractButton = qt.QPushButton("Extract Images")
    self.IMGSliceExtractFormLayout.addWidget(self.IMGSliceExtractButton)
    #---------------------------------------------------------
    
    # Connections
    self.input6Button.connect('clicked(bool)', self.onInput6Button)
    self.input7Button.connect('clicked(bool)', self.onInput7Button)   
    self.IMGSliceExtractButton.connect('clicked(bool)', self.onIMGSliceExtract)
    # End Connections
    
  def onInput6Button(self):
    self.inputPatientdirIMG = qt.QFileDialog.getExistingDirectory()
    self.input6Button.text = self.inputPatientdirIMG

  def onInput7Button(self):
    self.outputPatientdirIMG = qt.QFileDialog.getExistingDirectory()
    self.input7Button.text = self.outputPatientdirIMG
    
  def onIMGSliceExtract(self):
    fileIdentifierSettings = []  
    setting = {}  
    if self.niftiButton.checked: setting["Mask"] = '*.nii'
    else: setting["Mask"] = '*.nrrd'
    
    inputImageKeywords = [str(keyword.strip()) for keyword in self.inputImageKeywordsField.text.split(',')]
    inputImageExcludeKeywords = [str(keyword.strip()) for keyword in self.inputImageExclusionKeywordsField.text.split(',')]

    inputLabelKeywords = [str(keyword.strip()) for keyword in self.inputLabelKeywordsField.text.split(',')]
    inputLabelExcludeKeywords = [str(keyword.strip()) for keyword in self.inputLabelExclusionKeywordsField.text.split(',')]
    
    setting["selImage"] = [ inputImageKeywords , inputImageExcludeKeywords ]
    setting["selLabel"] = [ inputLabelKeywords , inputLabelExcludeKeywords ]
    setting["levels"] = [1]
    fileIdentifierSettings.append(setting)

    IMGRadiomicsToolsLib.ImageSnapshot.Execute(self, self.inputPatientdirIMG, self.outputPatientdirIMG, fileIdentifierSettings, self.IMGSliceExtractButton)
   
  def printStatus(caller, event):
    print("Got a %s from a %s" % (event, caller.GetClassName()))
    if caller.IsA('vtkMRMLCommandLineModuleNode'):
      print("Status is %s" % caller.GetStatusString())