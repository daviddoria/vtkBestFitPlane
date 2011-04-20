#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>
#include <vtkXMLPolyDataWriter.h>

#include "vtkBestFitPlane.h"

int main(int argc, char *argv[])
{
  // Verify arguments
  if(argc != 3)
    {
    std::cerr << "Required arguments: InputPoints(vtp) OutputPlane(vtp)" << std::endl;
    return EXIT_FAILURE;
    }

  // Parse arguments
  std::string inputPointsFileName = argv[1];
  std::string outputPlaneFileName = argv[2];

  vtkSmartPointer<vtkXMLPolyDataReader> reader =
    vtkSmartPointer<vtkXMLPolyDataReader>::New();
  reader->SetFileName(inputPointsFileName.c_str());
  reader->Update();
 
  // Compute the best fit plane.
  vtkSmartPointer<vtkBestFitPlane> bestFitFilter =
    vtkSmartPointer<vtkBestFitPlane>::New();
  bestFitFilter->SetInputConnection(reader->GetOutput()->GetProducerPort());
  bestFitFilter->Update();

  double normal[3];
  vtkDoubleArray* normalArray = vtkDoubleArray::SafeDownCast(bestFitFilter->
                      GetOutput()->GetFieldData()->GetArray("BestFitPlaneNormal"));
  normalArray->GetTupleValue(0,normal);
  std::cout << "Normal = (" << normal[0] << ", " << normal[1] << ", " << normal[2] << ")" << std::endl;

  double origin[3];
  vtkDoubleArray* originArray = vtkDoubleArray::SafeDownCast(bestFitFilter->
                      GetOutput()->GetFieldData()->GetArray("BestFitPlaneOrigin"));
  originArray->GetTupleValue(0,origin);
  std::cout << "Origin = (" << origin[0] << ", " << origin[1] << ", " << origin[2] << ")" << std::endl;

  vtkSmartPointer<vtkPlaneSource> planeSource =
    vtkSmartPointer<vtkPlaneSource>::New();
  planeSource->SetCenter(origin);
  planeSource->SetNormal(normal);
  planeSource->Update();
  
  vtkSmartPointer<vtkXMLPolyDataWriter> writer =  
    vtkSmartPointer<vtkXMLPolyDataWriter>::New();
  writer->SetFileName(outputPlaneFileName.c_str());
  writer->SetInputConnection(planeSource->GetOutputPort());
  writer->Write();
  
  return EXIT_SUCCESS;
}
