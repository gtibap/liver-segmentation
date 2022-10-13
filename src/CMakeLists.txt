PROJECT(EVALPARAM)

cmake_minimum_required(VERSION 2.6)

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

# Found ITK.
FIND_PACKAGE(ITK REQUIRED)
INCLUDE(${ITK_USE_FILE})

FIND_PACKAGE(vtkINRIA3D REQUIRED)
INCLUDE(${vtkINRIA3D_USE_FILE})


#------------------------------
# Metodo propuesto de segmentación del hígado en imágenes de TAC

# ADD_EXECUTABLE( evalParam preDeformacion.cxx InfoEstadisticaCorteAxial.cxx InicializacionDeformacion.cxx FuncionesDeVisualizacion.cxx MetodoDeDeformacion.cxx FuncionesDeDeformacion.cxx FuncionesBasicas.cxx RemalladoNuevo.cxx StopEvaluation.cxx WriterSegmentation.cxx )
#  AjusteParametros.cxx InfoEstadisticaCorteAxial.cxx SeleccionRegionHigado.cxx FuncionesDeVisualizacion.cxx ConstruccionSuperficie.cxx AjusteManual.cxx  
# TARGET_LINK_LIBRARIES( evalParam
# vtkRenderingAddOn vtkVisuManagement vtkRendering vtkVolumeRendering vtkIO vtkImaging vtkCommon vtkWidgets vtkGraphics vtkHybrid vtkFiltering ITKBasicFilters ITKCommon ITKStatistics ITKIO )


#------------------------------
# Metodo propuesto de segmentación del hígado en imágenes de TAC

ADD_EXECUTABLE( segHigado preDeformacion.cxx ajusteSuperficie.cxx FuncionesDeVisualizacion.cxx AjusteManual.cxx InfoEstadisticaCorteAxial.cxx InicializacionDeformacion.cxx MetodoDeDeformacion.cxx FuncionesDeDeformacion.cxx FuncionesBasicas.cxx RemalladoNuevo.cxx StopEvaluation.cxx WriterSegmentation.cxx)

TARGET_LINK_LIBRARIES( segHigado
vtkRenderingAddOn vtkVisuManagement vtkRendering vtkVolumeRendering vtkIO vtkImaging vtkCommon vtkWidgets vtkGraphics vtkHybrid vtkFiltering ITKBasicFilters ITKCommon ITKStatistics ITKIO )

