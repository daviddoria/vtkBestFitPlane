#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>

#include "vtkBestFitPlane.h"

#include <limits>
#include <cmath>

template<class A>
bool fuzzyCompare1D(A a, A b)
{
  return std::abs(a - b) < std::numeric_limits<A>::epsilon();
}

template<class A>
bool fuzzyCompare3D(A a[3], A b[3])
{
  return fuzzyCompare1D(a[0], b[0]) &&
         fuzzyCompare1D(a[1], b[1]) &&
         fuzzyCompare1D(a[2], b[2]);
}

int main(int, char *[])
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
  
  // Compute the best fit plane.
  vtkSmartPointer<vtkBestFitPlane> bestFitFilter =
    vtkSmartPointer<vtkBestFitPlane>::New();
  bestFitFilter->SetInputData(polydata);
  bestFitFilter->Update();

  double normal[3];
  vtkDoubleArray* normalArray = vtkDoubleArray::SafeDownCast(bestFitFilter->
                      GetOutput()->GetFieldData()->GetArray("BestFitPlaneNormal"));
  normalArray->GetTypedTuple(0,normal);
  std::cout << "Normal = (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")" << std::endl;

  double origin[3];
  vtkDoubleArray* originArray = vtkDoubleArray::SafeDownCast(bestFitFilter->
                      GetOutput()->GetFieldData()->GetArray("BestFitPlaneOrigin"));
  originArray->GetTypedTuple(0,origin);
  std::cout << "Origin = (" << origin[0] << ", " << origin[1] << ", " << origin[2] << ")" << std::endl;

  double correctOrigin[3] = {5.5, 0.5, 0};
  if(!fuzzyCompare3D(origin, correctOrigin))
    {
    std::cout << "ERROR!" << std::endl;
    }

  double correctNormal[3] = {0,0,1};
  if(!fuzzyCompare3D(normal, correctNormal))
    {
    std::cout << "ERROR!" << std::endl;
    }
    
  return EXIT_SUCCESS;
}
