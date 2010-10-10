/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkBestFitPlane.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkBestFitPlane - Find the least squares plane through a set of points
// .SECTION Description
// vtkBestFitPlane finds the least squares plane through a set of points.

#ifndef __vtkBestFitPlane_h
#define __vtkBestFitPlane_h

#include "vtkPolyDataAlgorithm.h"

class vtkBestFitPlane : public vtkPolyDataAlgorithm
{
public:
  static vtkBestFitPlane *New();
  vtkTypeMacro(vtkBestFitPlane,vtkPolyDataAlgorithm);

protected:
  vtkBestFitPlane(){}

  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector );
  int FillInputPortInformation(int port, vtkInformation* info);

private:

  vtkBestFitPlane(const vtkBestFitPlane&);  // Not implemented.
  void operator=(const vtkBestFitPlane&);  // Not implemented.
};

#endif
