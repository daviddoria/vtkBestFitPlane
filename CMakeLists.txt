cmake_minimum_required(VERSION 2.6)
set (CMAKE_CXX_STANDARD 11)

PROJECT(BestFitPlane)
ENABLE_TESTING()

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(vtkBestFitPlaneExample vtkBestFitPlaneExample.cxx
vtkBestFitPlane.cxx vtkCenterOfMass.cxx)
TARGET_LINK_LIBRARIES(vtkBestFitPlaneExample ${VTK_LIBRARIES})

ADD_EXECUTABLE(vtkCenterOfMassExample vtkCenterOfMassExample.cxx
vtkCenterOfMass.cxx)
TARGET_LINK_LIBRARIES(vtkCenterOfMassExample ${VTK_LIBRARIES})

############

ADD_EXECUTABLE(ComputeBestFitPlane ComputeBestFitPlane.cxx
vtkBestFitPlane.cxx vtkCenterOfMass.cxx)
TARGET_LINK_LIBRARIES(ComputeBestFitPlane ${VTK_LIBRARIES})
