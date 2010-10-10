cmake_minimum_required(VERSION 2.6)

PROJECT(BestFitPlane)
ENABLE_TESTING()

FIND_PACKAGE(VTK REQUIRED)
INCLUDE(${VTK_USE_FILE})

ADD_EXECUTABLE(vtkBestFitPlaneExample vtkBestFitPlaneExample.cxx
vtkBestFitPlane.cxx vtkCenterOfMass.cxx)
TARGET_LINK_LIBRARIES(vtkBestFitPlaneExample vtkHybrid )

ADD_EXECUTABLE(vtkCenterOfMassExample vtkCenterOfMassExample.cxx
vtkCenterOfMass.cxx)
TARGET_LINK_LIBRARIES(vtkCenterOfMassExample vtkHybrid )

