/*=========================================================================

  Program:   Visualization Toolkit
  Module:    vtkCenterOfMass.h

  Copyright (c) Ken Martin, Will Schroeder, Bill Lorensen
  All rights reserved.
  See Copyright.txt or http://www.kitware.com/Copyright.htm for details.

     This software is distributed WITHOUT ANY WARRANTY; without even
     the implied warranty of MERCHANTABILITY or FITNESS FOR A PARTICULAR
     PURPOSE.  See the above copyright notice for more information.

=========================================================================*/
// .NAME vtkCenterOfMass - Find the center of mass of a set of points.
// .SECTION Description
// vtkCenterOfMass finds the center of mass of a set of points.

#ifndef __vtkCenterOfMass_h
#define __vtkCenterOfMass_h

#include "vtkPolyDataAlgorithm.h"

class vtkCenterOfMass : public vtkPolyDataAlgorithm
{
public:
  static vtkCenterOfMass *New();
  vtkTypeMacro(vtkCenterOfMass,vtkPolyDataAlgorithm);

protected:
  vtkCenterOfMass(){}

  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector );

private:

  vtkCenterOfMass(const vtkCenterOfMass&);  // Not implemented.
  void operator=(const vtkCenterOfMass&);  // Not implemented.
};

#endif
