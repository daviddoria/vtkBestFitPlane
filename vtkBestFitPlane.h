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
  
  vtkBestFitPlane();

  int RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector ) override;
  int FillInputPortInformation(int port, vtkInformation* info) override;
  
  enum WeightEnum {UNIFORMWEIGHT, GAUSSIANWEIGHT};

  // Description:
  // Compute a weight based on the WeightMode.
  double WeightFunction(double distance);

  // Description:
  // Set the WeightMode to UNIFORMWEIGHT
  void SetWeightModeToUniform(){this->WeightMode = UNIFORMWEIGHT;}

  // Description:
  // Set the WeightMode to GAUSSIANWEIGHT
  void SetWeightModeToGaussian(){this->WeightMode = GAUSSIANWEIGHT;}
  
  // Description:
  // A flag specifying how to weight points.
  int WeightMode;

  // Description:
  // Accessor/mutator for GaussianVariance
  vtkSetMacro(GaussianVariance, double);
  vtkGetMacro(GaussianVariance, double);
  
  // Description:
  // The variance of the Gaussian function to use if
  // WeightMode == GAUSSIANWEIGHT
  double GaussianVariance;
private:

  vtkBestFitPlane(const vtkBestFitPlane&);  // Not implemented.
  void operator=(const vtkBestFitPlane&);  // Not implemented.
};

#endif
