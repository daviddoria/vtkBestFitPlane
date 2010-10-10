#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkPointSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>

#include "vtkCenterOfMass.h"

int main (int, char *[])
{
  vtkSmartPointer<vtkPointSource> pointSource =
    vtkSmartPointer<vtkPointSource>::New();
  pointSource->SetCenter(5,0,0);
  pointSource->Update();
  
  // Compute the center of mass
  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter =
    vtkSmartPointer<vtkCenterOfMass>::New();
  centerOfMassFilter->SetInputConnection(pointSource->GetOutputPort());
  centerOfMassFilter->Update();
  
  vtkDoubleArray* centerOfMass = vtkDoubleArray::SafeDownCast(centerOfMassFilter->GetOutput()->GetFieldData()->GetArray("CenterOfMass"));
  double center[3];
  centerOfMass->GetTupleValue(0,center);

  std::cout << "Center: " << center[0] << " " << center[1] << " " << center[2] << std::endl;

  return EXIT_SUCCESS;
}
