#include "vtkCenterOfMass.h"

#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkSmartPointer.h"

vtkStandardNewMacro(vtkCenterOfMass);

int vtkCenterOfMass::RequestData(vtkInformation* request,
                  vtkInformationVector** inputVector,
                  vtkInformationVector* outputVector )
{
  // Get the input and ouptut
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet* input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData* output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  double center[3];
  center[0] = 0.0;
  center[1] = 0.0;
  center[2] = 0.0;

  for(vtkIdType i = 0; i < input->GetNumberOfPoints(); i++)
    {
    double point[3];
    input->GetPoint(i, point);
    vtkMath::Add(center, point, center);
    }

  double numberOfPoints = static_cast<double>(input->GetNumberOfPoints());
  vtkMath::MultiplyScalar(center, 1./numberOfPoints);

  output->ShallowCopy(input);

  vtkSmartPointer<vtkDoubleArray> centerOfMass =
    vtkSmartPointer<vtkDoubleArray>::New();
  centerOfMass->SetNumberOfComponents(3);
  centerOfMass->SetName("CenterOfMass");
  centerOfMass->InsertNextTuple(center);

  output->GetFieldData()->AddArray(centerOfMass);

  return 1;
}
