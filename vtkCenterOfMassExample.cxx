#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include <limits>
#include <cmath>

#include "vtkCenterOfMass.h"

template<class A>
bool fuzzyCompare1D(A a, A b)
{
  return std::abs(a - b) < std::numeric_limits<A>::epsilon();
}

template<class A>
bool fuzzyCompare2D(A a[2], A b[2])
{
  return fuzzyCompare1D(a[0], b[0]) &&
         fuzzyCompare1D(a[1], b[1]);
}

template<class A>
bool fuzzyCompare3D(A a[3], A b[3])
{
  return fuzzyCompare1D(a[0], b[0]) &&
         fuzzyCompare1D(a[1], b[1]) &&
         fuzzyCompare1D(a[2], b[2]);
}

int main (int, char *[])
{
  vtkSmartPointer<vtkPoints> points =
    vtkSmartPointer<vtkPoints>::New();
  points->InsertNextPoint(5,0,0);
  points->InsertNextPoint(6,0,0);
  points->InsertNextPoint(5,1,0);
  points->InsertNextPoint(6,1,0);

  vtkSmartPointer<vtkPolyData> polydata =
    vtkSmartPointer<vtkPolyData>::New();
  polydata->SetPoints(points);
  
  // Compute the center of mass
  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter =
    vtkSmartPointer<vtkCenterOfMass>::New();
  centerOfMassFilter->SetInputData(polydata);
  centerOfMassFilter->Update();
  
  vtkDoubleArray* centerOfMass = vtkDoubleArray::SafeDownCast(centerOfMassFilter->GetOutput()->GetFieldData()->GetArray("CenterOfMass"));
  double center[3];
  centerOfMass->GetTypedTuple(0,center);

  std::cout << "Center: " << center[0] << " " << center[1] << " " << center[2] << std::endl;

  double correct[3] = {5.5, 0.5, 0};
  if(!fuzzyCompare3D(center, correct))
    {
    std::cout << "ERROR!" << std::endl;
    }
    
  return EXIT_SUCCESS;
}
