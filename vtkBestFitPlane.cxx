#include "vtkBestFitPlane.h"
#include "vtkCenterOfMass.h"

#include "vtkDataObject.h"
#include "vtkDoubleArray.h"
#include "vtkFieldData.h"
#include "vtkInformation.h"
#include "vtkInformationVector.h"
#include "vtkMath.h"
#include "vtkObjectFactory.h"
#include "vtkPlane.h"
#include "vtkPlaneSource.h"
#include "vtkSmartPointer.h"
#include "vtkStreamingDemandDrivenPipeline.h"

vtkStandardNewMacro(vtkBestFitPlane);

vtkBestFitPlane::vtkBestFitPlane()
{
  this->WeightMode = UNIFORMWEIGHT;
  this->GaussianVariance = 1.0;
}

double vtkBestFitPlane::WeightFunction(double distance)
{
  if(this->WeightMode == UNIFORMWEIGHT)
    {
    return 1;
    }
  else if(this->WeightMode == GAUSSIANWEIGHT)
    {
    return vtkMath::GaussianWeight(this->GaussianVariance, distance);
    }
  else
    {
    vtkErrorMacro(<< "Invalid WeightMode specified!");
    }
}

namespace //anonymous
{

/* allocate memory for an nrow x ncol matrix */
template<class TReal>
    TReal **create_matrix ( long nrow, long ncol )
{
  typedef TReal* TRealPointer;
  TReal **m = new TRealPointer[nrow];

  TReal* block = ( TReal* ) calloc ( nrow*ncol, sizeof ( TReal ) );
  m[0] = block;
  for ( int row = 1; row < nrow; ++row )
  {
    m[ row ] = &block[ row * ncol ];
  }
  return m;
}

/* free a TReal matrix allocated with matrix() */
template<class TReal>
    void free_matrix ( TReal **m )
{
  free ( m[0] );
  delete[] m;
}

}//end anonymous namespace

int vtkBestFitPlane::RequestData(
                                          vtkInformation *vtkNotUsed(request),
    vtkInformationVector **inputVector,
    vtkInformationVector *outputVector)
{

  // Get the input and ouptut
  vtkInformation *inInfo = inputVector[0]->GetInformationObject(0);
  vtkDataSet* input = vtkDataSet::SafeDownCast(
    inInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkInformation *outInfo = outputVector->GetInformationObject(0);

  vtkPolyData* output = vtkPolyData::SafeDownCast(
    outInfo->Get(vtkDataObject::DATA_OBJECT()));

  vtkIdType numPoints = input->GetNumberOfPoints();
  double dNumPoints = static_cast<double>(numPoints);

  // Find the center of mass of the points
  vtkSmartPointer<vtkCenterOfMass> centerOfMassFilter =
    vtkSmartPointer<vtkCenterOfMass>::New();
  centerOfMassFilter->SetInputConnection(input->GetProducerPort());
  centerOfMassFilter->Update();

  double center[3];
  vtkDoubleArray::SafeDownCast(centerOfMassFilter->GetOutput()->GetFieldData()->GetArray("CenterOfMass"))->GetTupleValue(0, center);

  //std::cout << "Center of mass: " << Center[0] << " " << Center[1] << " " << Center[2] << vtkstd::endl;

  // Compute sample covariance matrix
  double **a = create_matrix<double> ( 3,3 );
  a[0][0] = 0; a[0][1] = 0;  a[0][2] = 0;
  a[1][0] = 0; a[1][1] = 0;  a[1][2] = 0;
  a[2][0] = 0; a[2][1] = 0;  a[2][2] = 0;

  unsigned int i, j;

  double weightSum = 0;
  
  for(vtkIdType pointId = 0; pointId < numPoints; pointId++ )
    {
    double x[3];
    double xp[3];
    double distanceFromCenter = sqrt(vtkMath::Distance2BetweenPoints(x, center));
    double WeightFunction(distanceFromCenter);
    input->GetPoint(pointId, x);

    for(j = 0; j < 3; j++ )
      {
      xp[j] = x[j] - center[j];

      double w = this->WeightFunction(xp[j]);
      // Take advantage of the fact that the scatter materix is symmetric
      for (i = j; i < 3; i++)
        {
        double value = w * xp[j] * xp[i];
        // Store the computed value in the current position (j,i) and its
        // reflection (i,j).
        a[j][i] += value;
        a[i][j] += value;
        }
      weightSum += w;
      }
    }

  // Divide by N-1 for an unbiased estimate
  for( j = 0; j < 3; j++ )
    {
    for(i = 0; i < 3; i++)
      {
      a[j][i] /= weightSum;
      }
    }

  // Extract eigenvectors from covariance matrix
  double **eigvec = create_matrix<double> ( 3,3 );

  double eigval[3];
  vtkMath::Jacobi(a,eigval,eigvec);

  // Set the plane normal to the eigenvector corresponding
  // to the smallest eigenvalue.
  double normal[3];
  normal[0] = eigvec[0][2];
  normal[1] = eigvec[1][2];
  normal[2] = eigvec[2][2];

  output->ShallowCopy(input);

  // Store the result in the FieldData
  vtkSmartPointer<vtkDoubleArray> normalFieldData =
    vtkSmartPointer<vtkDoubleArray>::New();
  normalFieldData->SetNumberOfComponents(3);
  normalFieldData->SetName("BestFitPlaneNormal");
  normalFieldData->InsertNextTuple(normal);

  output->GetFieldData()->AddArray(normalFieldData);

  vtkSmartPointer<vtkDoubleArray> originFieldData =
    vtkSmartPointer<vtkDoubleArray>::New();
  originFieldData->SetNumberOfComponents(3);
  originFieldData->SetName("BestFitPlaneOrigin");
  originFieldData->InsertNextTuple(center);

  output->GetFieldData()->AddArray(originFieldData);

  // Cleanup
  free_matrix(eigvec);
  free_matrix(a);

  return 1;
}


int vtkBestFitPlane::FillInputPortInformation( int port, vtkInformation* info )
{
  if (port == 0)
    {
    info->Set( vtkAlgorithm::INPUT_REQUIRED_DATA_TYPE(), "vtkDataSet" );
    return 1;
    }

  return 0;
}
