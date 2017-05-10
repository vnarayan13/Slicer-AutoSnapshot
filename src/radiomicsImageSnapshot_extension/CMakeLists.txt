#-----------------------------------------------------------------------------
set(MODULE_NAME radiomicsimagesnapshot)

#-----------------------------------------------------------------------------
set(MODULE_PYTHON_SCRIPTS
  ${MODULE_NAME}.py
  IMGRadiomicsToolsLib/__init__
  IMGRadiomicsToolsLib/radiomicsDatabase
  IMGRadiomicsToolsLib/path
  IMGRadiomicsToolsLib/ImageSnapshot
  )

set(MODULE_PYTHON_RESOURCES
  )

#-----------------------------------------------------------------------------
slicerMacroBuildScriptedModule(
  NAME ${MODULE_NAME}
  SCRIPTS ${MODULE_PYTHON_SCRIPTS}
  RESOURCES ${MODULE_PYTHON_RESOURCES}
  WITH_GENERIC_TESTS
  )
