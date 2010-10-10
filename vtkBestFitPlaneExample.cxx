#include <vtkDoubleArray.h>
#include <vtkFieldData.h>
#include <vtkPlaneSource.h>
#include <vtkPolyData.h>
#include <vtkSmartPointer.h>
#include <vtkXMLPolyDataReader.h>

#include "vtkBestFitPlane.h"

int main (int argc, char *argv[])
{
  vtkSmartPointer<vtkPolyData> input;
  std::string inputFilename;
  std::string outputFilename;

  if(argc > 1) //If a file name is specified, open and use the file.
    {

    //verify command line arguements
    if(argc != 3)
      {
      std::cout << "Required arguments: InputFilename(vtp) OutputFilename(vtp)" << std::endl;
      return EXIT_FAILURE;
      }
    inputFilename = argv[1];
    outputFilename = argv[2];

    vtkSmartPointer<vtkXMLPolyDataReader> reader =
      vtkSmartPointer<vtkXMLPolyDataReader>::New();
    reader->SetFileName(argv[1]);
    reader->Update();
    input = reader->GetOutput();
    }
  else //If a file name is not specified, create a sphere
    {
    vtkSmartPointer<vtkPlaneSource> planeSource =
      vtkSmartPointer<vtkPlaneSource>::New();
    planeSource->SetXResolution(10);
    planeSource->SetYResolution(10);
    planeSource->Update();

    input = vtkSmartPointer<vtkPolyData>::New();
    input->SetPoints(planeSource->GetOutput()->GetPoints());
    }

  // Compute the best fit plane.
  vtkSmartPointer<vtkBestFitPlane> bestFitFilter =
    vtkSmartPointer<vtkBestFitPlane>::New();
  bestFitFilter->SetInputConnection(input->GetProducerPort());
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

  return EXIT_SUCCESS;
}
