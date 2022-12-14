// #include "vtkSphereSource.h"
#include "vtkActor.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkCleanPolyData.h"
#include "vtkConeSource.h"
#include "vtkFloatArray.h"
#include "vtkGlyph3D.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkImageAccumulate.h"
#include "vtkImageCast.h"
#include "vtkImageData.h"
#include "vtkImageGradient.h"
#include "vtkImageMagnitude.h"
#include "vtkImageMedian3D.h"
#include "vtkImageShiftScale.h"
#include "vtkImageStencil.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkMaskPoints.h"
#include "vtkMetaImageWriter.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSimpleImageFilterExample.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTriangleFilter.h"
#include "vtkUnsignedCharArray.h"
#include "vtkFeatureEdges.h"
#include "vtkDataSetSurfaceFilter.h"

// #include "vtkPolyDataReader.h"
// #include "vtkProgrammableFilter.h"

#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
// #include "vtkPolyDataWriter.h"
// #include "vtkSTLReader.h"
#include <stdio.h>
#include <time.h>
#include <iostream>
using namespace std;

void statisticRegionImagesClosedMesh( vtkImageData *images, vtkPolyData *liverMesh, double &meanRegion, double &sigmaRegion );

//Funciones para el suavizado de la malla
void findNeighbors1( vtkPolyData *input, vtkIdList **nbs );
void findCentroids( vtkPolyData *input, vtkIdList **lnbs, vtkDoubleArray *lcentroids );
void tensionForce( vtkPolyData *input, vtkDoubleArray *lcentroids, vtkDoubleArray *lTForce );

void rigidityForce1( vtkPolyData *input, vtkDoubleArray *lTForce, vtkIdList **lnbs, vtkDoubleArray *lRForce );
void rigidityForce2( vtkPolyData *input, vtkDoubleArray *lTForce, vtkDataArray *lNormalVectors, vtkDoubleArray *lRForce );

void relationIdMeshIdImages( vtkPoints *lMeshPoints, vtkImageData *binaryIntensityImages, vtkIdTypeArray *idMeshIdImages );

void intensityForce( vtkDataArray *lIntensityImages, vtkDataArray *lNormalVectors, vtkIdTypeArray *lIdMeshIdImages, double &mean, double &sigma, vtkDoubleArray *lIntensityForce, int selFP );
// void gradientEdgeImages( vtkImageData *magnitudeSobelImages, double c, vtkDataArray *lGradientEdgeImages );
void edgeForce( vtkDataArray *lGradientEdgeImages, vtkDataArray *lNormalVectors, vtkIdTypeArray *lIdMeshIdImages, vtkDoubleArray *lEdgeForce );

void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, vtkPoints *lNewMeshPoints );
void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, vtkPoints *lNewMeshPoints );
void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, double &p3, vtkDoubleArray *lForces3, vtkPoints *lNewMeshPoints );
void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, double &p3, vtkDoubleArray *lForces3, double &p4, vtkDoubleArray *lForces4, vtkPoints *lNewMeshPoints );

void movePointsMesh( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, double &p3, vtkDoubleArray *lForces3, double limitImages[6], vtkPoints *lNewMeshPoints );
// void deformationMesh( vtkPoints *lMeshPoints, vtkDoubleArray *lBalloonForce, vtkPoints *lNewMeshPoints );
// int fusionDeVertices( vtkPolyData *input, float minimumDistance );
int initialRemesh( vtkPolyData *liverMesh, float minimumDistance );
int remesh( vtkPolyData *liverMesh, float minimumDistance, vtkPoints *initialPoints, vtkDoubleArray *errorDistance, vtkPolyData *outputMesh );

int normalsCalculation( vtkPolyData *liverMesh, vtkDataArray *lNormalVectors );

int stopEvaluation( vtkDoubleArray *errorDistance, FILE *fpt, int lCounts[3][3], int medianlCounts[3] );

void MeshVisualization( vtkPolyData *mesh );
void NormalsMeshVisualization( vtkPolyData *normalsMesh );
void SynchronizedViews( vtkImageData *imagesCT, vtkPolyData *meshLiverModelTransformed, int range[2]);
void ViewImages( vtkImageData *invertSobelImages, int range[2] );
void VisualizationMeshComparation( vtkPolyData *initialLiver, vtkPolyData *liverMesh );

int WriterSegmentation( vtkImageData *segImages, vtkPolyData *meshLiver, int *regionBB, int sizeOrig[3], char *outSegFileName );

// void times(double s, double A[3], double C[3]);
// double norm (double A[3]);
// void divide (double A[3], double s, double C[3]);
// void vplus (double A[3], double B[3], double C[3]);
// void vminus (double A[3], double B[3], double C[3]);

/**
Inicia el programa para la deformaci??n de la malla. 
Esta funci??n es llamada por la funci??n deformationInit()
 */
// int deformation (vtkPolyData *initialLiver, vtkPolyData *initialLiver1, vtkImageData *binaryIntensityImages, vtkImageData *magnitudeSobelImages, int iter ) {
int metodoDeDeformacion (vtkPolyData *inputLiverMesh, vtkImageData *intensityImages, double *infoIntensity, vtkImageData *magnitudeSobelImages, double *infoMagGradient, int iter, int *OrigDims, int *bb, char outSegFileName[], double cs, double cp, double cb, int selFP ) {
  
//   cout << "Deformaci??n de la malla por la ecuaci??n de movimiento... " << endl;
  vtkCellArray *polys;
  vtkDataArray *lNormalVectors;
  vtkDoubleArray *errorDistance;
  vtkIdTypeArray *lIdMeshIdImages;
  vtkPoints *lInitialMeshPoints;
  vtkPoints *lMeshPoints;
  vtkPoints *lNewMeshPoints;
  vtkPoints *points;
  vtkPolyData *liverMesh;
  vtkPolyData *outputMesh;
    
  vtkIdType numberOfPointsMesh;
  
  float minimumDistance;
  int range[2] = {-1000, 1000};
//   double p = 1.0, q = 1.0; // p y q son los factores de ponderaci??n de las fuerzas de balloon y edge, respectivamente.
  // -----------------------------------------  
  
  /** Remallado de la malla que representa el h??gado en funci??n del tama??o del voxel. En este caso, la distancia m??nima de los bordes de la malla ser?? la magnitud de la diagonal de cada voxel en las im??genes */
  double spacingVoxel[3];
  
  intensityImages->GetSpacing( spacingVoxel );
  cout << "spacing: " << spacingVoxel[0] << ", " << spacingVoxel[1] << ", " << spacingVoxel[2] << endl;
  
//   minimumDistance = sqrt( spacingVoxel[0]*spacingVoxel[0] + spacingVoxel[1]*spacingVoxel[1] + spacingVoxel[2]*spacingVoxel[2] );
  minimumDistance = 2.0;  //  2.0mm
  cout << "minDistance= " << minimumDistance << endl;
  
  // copy the input (only polys) to our working mesh
  liverMesh = vtkPolyData::New();
  points = vtkPoints::New();
  polys = vtkCellArray::New();
  
  points->DeepCopy(inputLiverMesh->GetPoints());
  liverMesh->SetPoints(points);
  points->Delete();
  polys->DeepCopy(inputLiverMesh->GetPolys());
  liverMesh->SetPolys(polys);
  polys->Delete();
  liverMesh->GetPointData()->DeepCopy(inputLiverMesh->GetPointData());
  liverMesh->GetFieldData()->PassData(inputLiverMesh->GetFieldData());
  liverMesh->BuildCells();
  liverMesh->BuildLinks();
  

  /**  Visualizaci??n de la malla. */
  cout << "Visualizaci??n de la malla copiada..." << endl;
//   MeshVisualization( liverMesh );
  cout << "Visualizaci??n realizada." << endl;
  /**  Visualizaci??n de la malla.*/
 
  /**  Remallado inicial */
  float longitud = 1.0;
  cout << "Remallado inicial..." << endl;
  while( longitud < minimumDistance ){
    initialRemesh( liverMesh, longitud );
    cout << "pointsMesh = " << liverMesh->GetNumberOfPoints() << endl;
    longitud+=0.2;
//     cout << "Distancia min: " << longitud << endl;
  }
  initialRemesh( liverMesh, minimumDistance );
  cout << "pointsMesh = " << liverMesh->GetNumberOfPoints() << endl;
  /**  Remallado inicial */
  
  /**--------------------------
  Estas listas de puntos servir??n para calcular la distancia de desplazamiento de los puntos de la malla en cada iteraci??n de la deformaci??n de la malla  */
  
//   vtkPoints *initialPoints = vtkPoints::New();
//   vtkIdType numberOfPointsMesh;
  
//   errorDistance->SetNumberOfComponents(1);
  /**------------------------------------
  Estas listas de puntos servir??n para calcular la distancia de desplazamiento de los puntos de la malla en cada iteraci??n de la deformaci??n de la malla  */

/**---------------------------
  Visualizaci??n de la malla... */
//   MeshVisualization( liverMesh );
  /**  Visualizaci??n de la malla.
  ---------------------------*/
  // -----------------------------------------
  cout << "Deformaci??n de la malla por la ecuaci??n de movimiento, la cu??l se compone b??sicamente de una fuerza de rigidez de la malla y una fuerza de acople a las im??genes ... " << endl;
  /** La ecuaci??n de movimiento de los puntos de la malla, es una versi??n discretizada de la ecuaci??n diferencial de primer orden definada por:
   dX
  ---- = FuerzaDeRigidez + FuerzaDeAcople
   dt
  La fuerza de rigidez busca acercar cada punto al centroide de sus vecinos mientras preserva la forma inicial.
  La fuerza de acople busca llevar cada punto hacia los bordes de las im??genes tomando en cuenta los valores de intensidad de los pixeles m??s cercanos a la posici??n de los puntos de la malla.
  */
  
  //-----------------------------  
  cout << "C??lculo de la fuerza de acople a las im??genes... " << endl;
  /**Se comienza calculando la fuerza de acople (que estar?? en direcci??n del vector normal asociado a cada v??rtice), la cual se compone de dos t??rminos:
  FuerzaDeAcople = FuerzaDeBordes + FuerzaDeIntensidad */
  
  /**C??lculo de la FuerzaDeBordes:
                                         1
  FuerzaDeBordes = [ Grad(  ------------------------------) ?? normalVector ] normalVector
                              1 + (MagnitudGradiente)^p
  
  Donde la MagnitudGradiente es la magnitud del resultado de la aplicaci??n del operador Sobel (magnitudeSobelImages), lo cual entrega la informaci??n de bordes de las im??genes suavizadas; p es un factor de pondraci??n (valdr?? 1 o 2), Grad es el operador gradiente, (??) representa el producto punto, y, normalVector es el vector normal asociado con cada v??rtice de la malla.
  */
  cout << "1. C??lculo de la fuerza de bordes... " << endl;
  //dado que el grad (1 / (1 + y*magnitudeSobelImages)) es siempre constante, se calcula primero y, por supuesto, una sola vez 
  cout << "1.1 C??lculo del gradiente del inverso de uno m??s la magnitud del gradiente elevado a p... " << endl;  
  //El primer paso para generar las im??genes de la FuerzaDeBordes es localizar el espacio de memoria del mismo tama??o de las im??genes de entrada "magnitudeSobelImages" (es necesario generar las im??genes para la aplicaci??n del operador gradiente).

  /**-----------------------------------*/
  
  // en esta parte se visualizan las im??genes de sobel
  cout << "Visualizaci??n del resultado del operador Sobel sobre las Im??genes... " << endl;
//   range[0] = int(infoMagGradient[0] + infoMagGradient[1]);
/*  range[0] = 0;
  range[1] = int(infoMagGradient[0] + 3*infoMagGradient[1]);
  ViewImages( magnitudeSobelImages, range );*/
  cout << "Visualizaci??n realizada." << endl;
  
  /**-----------------------------------*/
  
  double meanSobel, sigmaSobel, medianSobel;
  
  /**----------------------------------------
  Informaci??n estad??stica (media y desviaci??n estandar) de la magnitud del gradiente (calculada por la aplicaci??n del operador Sobel) de un solo corte segmentado manualmente */
  
  cout << "Informaci??n de media y desviaci??n estandar de la magnitud del resultado del operador sobel sobre la regi??n de las im??genes encerrada por malla registrada... " << endl; 
//   statisticRegionImagesClosedMesh( magnitudeSobelImages, liverMesh, meanSobel, sigmaSobel );
  meanSobel = infoMagGradient[0];
  sigmaSobel = infoMagGradient[1];
  medianSobel = infoMagGradient[2];
  cout << "meanSobel = " << meanSobel << endl;
  cout << "sigmaSobel = " << sigmaSobel << endl;
  cout << "medianSobel = " << medianSobel << endl;
  //Esta informaci??n estad??stica permite hacer un filtrado en g(I) de los datos de magnitudeSobelImages para la estimaci??n del gradiente
  
  /**Esta informaci??n servir?? para eliminar valores de gradiente peque??os en las im??genes.
  ----------------------------------------*/

  /**C??lculo de la funci??n de velocidad g(I) = 1 / (1 + magSobel^p), donde la magSobel depende de la informaci??n estad??stica estimada, de tal forma que ser??n cero para los valores de magnitudeSobelImages menores a la media, multiplicados por valores entre 0 y 1 de forma lineal para los valores entre la media y la media + la desviaci??n estandar, y se mantiene el valor en otro caso.*/
  cout << "C??lculo de la funci??n de velocidad g(I) la cual esta en funci??n del gradiente (operador sobel)..." << endl;
  
  int dims[3];
  double spacing[3], origin[3];
  
  
  vtkIdType nov; //nov number of voxels
  vtkImageData *velocityFunctionImages = vtkImageData::New();
  
  velocityFunctionImages->SetScalarTypeToDouble();
  velocityFunctionImages->SetNumberOfScalarComponents( 1 );
  
  magnitudeSobelImages->GetDimensions( dims );
  magnitudeSobelImages->GetSpacing( spacing );
  magnitudeSobelImages->GetOrigin( origin );
  
  velocityFunctionImages->SetDimensions( dims );
  velocityFunctionImages->SetSpacing( spacing );
  velocityFunctionImages->SetOrigin( origin );
  
  cout << "dims: " << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;
  cout << "spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << endl;
  cout << "origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << endl;

  /** ----------------------
  Definici??n coeficientes de ponderaci??n fuerzas de deformaci??n (cp y cb) **/
  double mdv;
  mdv = spacing[0];
  if (mdv > spacing[1])
    mdv = spacing[1];
  if(mdv > spacing[2])
    mdv = spacing[2];
  cout << "minima distancia entre voxeles: " << mdv << endl;
  /** Definici??n coeficientes de ponderaci??n fuerzas de deformaci??n (cp y cb) 
  -------------------------**/
  velocityFunctionImages->AllocateScalars();
  //en este punto ya se tiene el espacio en memoria para el volumen de las im??genes, para escribir los datos de la funci??n de velocidad g(I)
  
  // se almacenan los valores de magnitudeSobelImages en la lista lMagnitudeSobelImages
  vtkDataArray *lMagnitudeSobelImages = magnitudeSobelImages->GetPointData()->GetScalars();
  
  nov = lMagnitudeSobelImages->GetNumberOfTuples();
  cout << "nov: " << nov << endl;
  
  //se define la lista lVelocityFunctionImages para almacenar el resultado del inverso del gradiente
  vtkDoubleArray *lVelocityFunctionImages = vtkDoubleArray::New();
  lVelocityFunctionImages->SetNumberOfTuples( nov );
  lVelocityFunctionImages->SetNumberOfComponents( 1 );
  
  double valueSobel[1], valueVelocityFunction[1];
  double y;
//   int cuenta=0, cuenta1=0;
  //aqu?? se define una funci??n de ponderaci??n para eliminar valores peque??os de la magnitud del gradiente obtenidos por el operador sobel, la cual depende de la media y la desviaci??n estandar de estos valores, tal que: Son cero los valores menores que la media, se multiplican por un valor entre cero y 1 (de forma lineal) si los valores est??n entre la media y la media+sigma, y, son los mismos (multiplicados por 1) para los valores mayores a la media+sigma.
  for( int idVoxel=0; idVoxel < nov; idVoxel++ )
  {
    lMagnitudeSobelImages->GetTuple( idVoxel, valueSobel );
    y = valueSobel[0] / (medianSobel);
    valueVelocityFunction[0] = 1.0 / ( 1.0 + y*y );
    lVelocityFunctionImages->SetTuple( idVoxel, valueVelocityFunction );
  }
  
  velocityFunctionImages->GetPointData()->SetScalars( lVelocityFunctionImages );
  velocityFunctionImages->Modified();
  velocityFunctionImages->Update();
  
/*  cout << "cuenta= " << cuenta << endl;
  cout << "cuenta1= " << cuenta1 << endl;*/
  cout << "Funci??n de velocidad g(I) calculada." << endl;
  
  /**--------------
  Guardar las im??genes de la velocidad g(I)*/
/*  cout << "Guardando los resultados de las im??genes de velocidad..." << endl;
  vtkImageShiftScale *scaleVImagesFilter = vtkImageShiftScale::New();
  scaleVImagesFilter->SetScale(255.0);
  scaleVImagesFilter->SetOutputScalarTypeToUnsignedChar();
  
  vtkMetaImageWriter *writerVImages = vtkMetaImageWriter::New();
  
  scaleVImagesFilter->SetInput( velocityFunctionImages );
  writerVImages->SetInputConnection( scaleVImagesFilter->GetOutputPort() );
  writerVImages->SetFileName( "imagesRegions/velTest03.mhd" );
  writerVImages->SetRAWFileName( "imagesRegions/velTest03.raw" );
  try{
    writerVImages->Write();
  }
  catch( exception& ){
    cerr << "Error en la escritura." << endl;
  }
  cout << "Realizada." << endl;
  writerVImages->Delete();
  scaleVImagesFilter->Delete();*/
  /**Guardar las im??genes de la velocidad g(I)
  -------------*/
  /**-----------
  Visualizaci??n funci??n de velocidad */
  cout << "Visualizaci??n de la funci??n de velocidad g(I)... ";
/*  range[0] = 0;
  range[1] = 1;
  ViewImages( velocityFunctionImages, range );*/
  cout << "realizada." << endl;
  /** Visualizaci??n funci??n de velocidad */
  
  cout << "lVelocityFunctionImages components: " << lVelocityFunctionImages->GetNumberOfComponents() << endl;
  cout << "lVelocityFunctionImages tuples: " << lVelocityFunctionImages->GetNumberOfTuples() << endl;
  
  //luego de calcular g(I) = 1 / (1 + y*magnitudSobelImages), se calcular?? el gradiente de este resultado
  cout << "C??lculo del gradiente de la funci??n de velocidad g(I)... " << endl;
  vtkImageGradient *gradientImages = vtkImageGradient::New();
  gradientImages->SetDimensionality( 3 );
  gradientImages->SetInput( velocityFunctionImages );
  
  gradientImages->Update();
  cout << "gradiente calculado." << endl;
  
  //-------
  lVelocityFunctionImages->Delete();
  velocityFunctionImages->Delete();
  //-------
  /**--------------
  Guardar las im??genes de la velocidad g(I)*/
/*  vtkImageMagnitude *magGradientImages = vtkImageMagnitude::New();
  magGradientImages->SetInputConnection( gradientImages->GetOutputPort() );

  cout << "Guardando los resultados de las im??genes del gradiente de las im??genes de velocidad..." << endl;
  vtkImageShiftScale *scaleVImagesFilter2 = vtkImageShiftScale::New();
  scaleVImagesFilter2->SetScale(255.0);
  scaleVImagesFilter2->SetOutputScalarTypeToUnsignedChar();
  
  vtkMetaImageWriter *writerVImages2 = vtkMetaImageWriter::New();
  
  scaleVImagesFilter2->SetInputConnection( magGradientImages->GetOutputPort() );
  writerVImages2->SetInputConnection( scaleVImagesFilter2->GetOutputPort() );
  writerVImages2->SetFileName( "imagesRegions/gradVelTest03.mhd" );
  writerVImages2->SetRAWFileName( "imagesRegions/gradVelTest03.raw" );
  try{
  writerVImages2->Write();
}
  catch( exception& ){
  cerr << "Error en la escritura." << endl;
}
  cout << "Realizada." << endl;
  writerVImages2->Delete();
  scaleVImagesFilter2->Delete();
  magGradientImages->Delete();*/
  /**Guardar las im??genes de la velocidad g(I)
  -------------*/
  
  /** Visualizaci??n de la magnitud de los vectores gradiente de las im??genes de la funci??n de velocidad */
  cout << "Visualizaci??n magnitud del gradiente de la funci??n de velocidad... ";
/*  vtkImageMagnitude *magGradientImages = vtkImageMagnitude::New();
  magGradientImages->SetInputConnection( gradientImages->GetOutputPort() );
  magGradientImages->Update();
  range[0] = 0;
  range[1] = 1;
  ViewImages( magGradientImages->GetOutput(), range );
  magGradientImages->Delete();*/
  cout << "realizada." << endl;
  /** Visualizaci??n de la magnitud de los vectores gradiente de las im??genes de la funci??n de velocidad */
  
  //-----------------------------
  /**Los valores de los vectores del gradiente de las im??genes se pasan a la lista
  !!!   lGradientEdgeImages   ??????
  */  
  //se pasan los valores escalares del resultado del gradiente a la lista lGradientEdgeImages
  vtkDataArray *lGradientEdgeImages = gradientImages->GetOutput()->GetPointData()->GetScalars();
  cout << "tuples gradient = " << lGradientEdgeImages->GetNumberOfTuples() << endl;
  cout << "components gradient = " << lGradientEdgeImages->GetNumberOfComponents() << endl;
  //-----------------------------
  
  /**
  C??lculo de la FuerzaDeIntensidad
  La Fuerza de Intensidad estar?? en direcci??n normal, y ponderada por la estimaci??n de la funci??n gausiana que representa la distribuci??n de intensidad del volumen del h??gado. Esta ponderaci??n estar?? dada por
  
  p = exp( -(x-b)^2 / (2*c^2) ),
  
  donde b es la media (mean) y c es la desviaci??n estandar (sigma), de la distribuci??n gausiana de probabilidad de la intensidad.
  
  Se cuenta con medidas "promedio" de la media y de la desviaci??n estandar estimadas a partir de m??s de 20 h??gados segmentados manualmente; y adem??s, puede extraerse informaci??n de la media y de la desviaci??n estandar directamente sobre el volumen de las im??genes que encierra la malla promedio que ha sido previamente registrada manualmente.
  */
  
  /** La fuerza basada en las intensidades se calcular?? iterativamente para cada v??rtice de la malla. La ponderaci??n de esta fuerza depender?? de la media y la desviaci??n estandar de la distribuci??n de probabilidad de la intensidad en las im??genes
  */

  /**----------------------------------------
  Informaci??n estad??stica (media y desviaci??n estandar) de las intensidades de las im??genes (luego del filtro de mediana) de un solo corte segmentado manualmente */

  cout << "Informaci??n de la distribuci??n de probabilidad de las intensidades del h??gado, a partir de la segmentaci??n manual de ??ste en un corte axial..." << endl;
  double meanRegion, sigmaRegion, medianRegion;
  //  aqu?? se usan las clases vtkPolyDataToImageStencil y vtkImageAcumulate, y tal vez vtkImageStencil
//   statisticRegionImagesClosedMesh( intensityImages, liverMesh, meanRegion, sigmaRegion );
  
  meanRegion = infoIntensity[0];
  sigmaRegion = infoIntensity[1];
  medianRegion = infoIntensity[2];
  
  cout << "meanRegion = " << meanRegion << endl;
  cout << "sigmaRegion  = " << sigmaRegion << endl;
  cout << "medianRegion  = " << medianRegion << endl;
  /** Esta informaci??n es utilizada para definir la funci??n de fuerza basada en la intensidad (del modelo deformable)
  ----------------------------------------*/
  
  /** referencia a la lista de los valores de intensidad de las im??genes por el puntero lIntensityImages
  */
  vtkDataArray *lIntensityImages = intensityImages->GetPointData()->GetScalars();
  cout << "intensity tuples: " << lIntensityImages->GetNumberOfTuples() << endl;
  cout << "intensity components: " << lIntensityImages->GetNumberOfComponents() << endl;
  
  // -----------------------------------------  
  
  /**Ahora se encontrar?? la relaci??n entre los ids de los v??rtices de la malla con los ids m??s cercanos a estos puntos sobre las im??genes. Este proceso tambi??n debe hacerse recursivamente (es decir, cada vez que son desplazados los v??rtices de la malla). La relaci??n se referencia a la lista
  lIdMeshIdImages
   */
    //  Paso de las coordenadas de los puntos de la malla a la lista referenciada por lMeshPoints  
//   lMeshPoints = vtkPoints::New();
//   lMeshPoints->DeepCopy( liverMesh->GetPoints() );
  //Referencia realizada
//   cout << "points Mesh= " << lMeshPoints->GetNumberOfPoints() << endl;

//   cout << "Relaci??n entre los ids de la malla con los ids de las im??genes... ";
  
//   lIdMeshIdImages = vtkIdTypeArray::New();
//   lIdMeshIdImages->SetNumberOfComponents( 1 );

//   relationIdMeshIdImages( lMeshPoints, intensityImages, lIdMeshIdImages );
//   cout << "realizada." << endl;
  /**Realizada la relaci??n entre los ids de la malla con los de las im??genes*/

  // -----------------------------------------  
  /**Definici??n de las variables, c??lculo iterativo de las fuerzas de deformaci??n y de las nuevas posiciones de cada v??rtice de la malla*/ 
  cout << "Definici??n de las variables necesarias adicionales... " << endl;
 
  //   C??lculo de la fuerza de Balloon, la c??al es la encargada de mover los v??rtices de la malla en funci??n del vector normal unitario (asociado a cada v??rtice) y al rango de intensidades del h??gado, para que la malla se mueva hacia adentro o hacia afuera, acercandose a los bordes de este ??rgano.
  
  // Se inicializa una lista para almacenar los valores que representan los vectores de la fuerza de Balloon
  vtkDoubleArray *lIntensityForce;
  vtkDoubleArray *lEdgeForce;
  //inicializaci??n de variables
//   vtkIdType nop = liverMesh->GetNumberOfPoints();
//   cout << "nop= " << nop << endl;
  
//   vtkIdList **lnbs = new vtkIdList*[nop];
  vtkIdList **lnbs;
  vtkDoubleArray *lcentroids;
  vtkDoubleArray *lTensionForce;
  
/*  vtkDoubleArray *lRigidityForce;
  lRigidityForce = vtkDoubleArray::New();
  lRigidityForce->SetNumberOfComponents(3);*/
  
  cout << "n??mero de iteraciones: " << iter << endl;
  
  int lCounts[3][3] = {}; //  Cuentas se inicializa en ceros
  int medianCounts[3];

  
  double w1, w3, w4; //valores de ponderaci??n de cada una de las fuerzas de deformaci??n
  w1= cs;// 0.7 ponderaci??n fuerza de tensi??n
//   double  w2=-1.0;//ponderaci??n fuerza de rigidez
  w3=-cb*mdv;  //  0.9mm ponderaci??n fuerza de bordes
  w4= cp*mdv;  //  0.5mm ponderaci??n fuerza de intensidad (basada en la fuerza de balloon)
  minimumDistance = 2.0;  //  2.0mm
  
  cout << "coeficientes cs, cp, cb: " << w1 << ", " << w4 << ", " << w3 << endl; 
  /**------------------------------------
  escritura de los datos obtenidos en un archivo de texto  */
  
  FILE *fpt;
//   char infoFileName[] = "/media/disk/Tesis/ImagesTraining/MICCAI2007/meshRegions/InfoRegistroManual.dat";
  char infoFileName[] = "RespuestaDelSistema.dat";
  
  fpt = fopen( infoFileName , "a" );
  
  fprintf( fpt, "\n%s\n", outSegFileName );
  fprintf( fpt, "cs=%f, cp=%f, cb=%f\n", w1, w4, w3);
  fprintf( fpt, "iter\t numVertices\t N??mero de vertices de acuerdo con sus desplazamientos entre [0, 1] con pasos de 0.1 (primera columna, vertices que se desplazan entre 0 y 0.1 mm, en la segunda, los que se desplazan entre 0.1 y 0.2, y as?? hasta la columna 10, de los v??rtices que se desplazan entre 0.9 a 1.0 )\n");
  /** realizada la escritura de los datos obtenidos en un archivo de texto 
  ------------------------------------------------------*/
  // Se obtiene la informaci??n de los l??mites de las im??genes 
  double originImages[3], spacingImages[3], limitImages[6];
  int dimensionImages[3];
  
  intensityImages->GetOrigin(originImages);
  cout << "originImages: " << originImages[0] <<", " << originImages[1] << ", " << originImages[2] << endl;
  
  intensityImages->GetSpacing(spacingImages);
  cout << "spacingImages: " << spacingImages[0] <<", " << spacingImages[1] << ", " << spacingImages[2] << endl;

  intensityImages->GetDimensions(dimensionImages);
  cout << "dimensionImages: " << dimensionImages[0] <<", " << dimensionImages[1] << ", " << dimensionImages[2] << endl;
  
  // se definen los l??mites de las im??genes
  limitImages[0] = originImages[0];
  limitImages[1] = originImages[0] + spacingImages[0]*(dimensionImages[0] - 1.5);
  limitImages[2] = originImages[1];
  limitImages[3] = originImages[1] + spacingImages[1]*(dimensionImages[1] - 1.5);
  limitImages[4] = originImages[2];
  limitImages[5] = originImages[2] + spacingImages[2]*(dimensionImages[2] - 1.5);
  
  // si la malla tiene puntos ubicados fuera de las im??genes, entonces se llevan a los l??mites
  lMeshPoints = vtkPoints::New();
  lNewMeshPoints = vtkPoints::New();
  double coordPointMesh[3];
  
  lMeshPoints->DeepCopy( liverMesh->GetPoints() );
  numberOfPointsMesh = lMeshPoints->GetNumberOfPoints();
  cout << "PointsMesh: " << numberOfPointsMesh << endl;
    
  lNewMeshPoints->SetNumberOfPoints( numberOfPointsMesh );
    
  for(int i=0; i < numberOfPointsMesh; i++){
    lMeshPoints->GetPoint(i, coordPointMesh);
    if( coordPointMesh[0] < limitImages[0] || coordPointMesh[0] > limitImages[1] ||         coordPointMesh[1] < limitImages[2] || coordPointMesh[1] > limitImages[3] ||         coordPointMesh[2] < limitImages[4] || coordPointMesh[2] > limitImages[5] ){
      
      if(coordPointMesh[0] < limitImages[0]){
        coordPointMesh[0] = limitImages[0];
      }else{
        if(coordPointMesh[0] > limitImages[1]){
          coordPointMesh[0] = limitImages[1];
        }
      };
      if(coordPointMesh[1] < limitImages[2]){
        coordPointMesh[1] = limitImages[2];
      }else{
        if(coordPointMesh[1] > limitImages[3]){
          coordPointMesh[1] = limitImages[3];
        }
      };
      if(coordPointMesh[2] < limitImages[4]){
        coordPointMesh[2] = limitImages[4];
      }else{
        if(coordPointMesh[2] > limitImages[5]){
          coordPointMesh[2] = limitImages[5];
        }
      };
    };
    lNewMeshPoints->SetPoint(i, coordPointMesh);
  }
  
  liverMesh->SetPoints( lNewMeshPoints );
  
  lNewMeshPoints->Delete();
  lMeshPoints->Delete();
/**---------------------------
  Visualizaci??n de la malla... */
//   MeshVisualization( liverMesh );
  /**  Visualizaci??n de la malla.
  ---------------------------*/
  
//   Antes del proceso de deformaci??n es necesario guardar el volumen inicial encerrado por la superficie manualmente ajustada. Para esto se utiliza el procedimiento para escribir las im??genes segmentadas a partir de la superficie incial.
  /**-----------------------
  Escritura de la segmentaci??n resultante de la superficie del h??gado "promedio" en im??genes binarias  */
/*  int flagIMWriter=0;
  cout << "Escribir la segmentaci??n de la superficie inicial?(0/1): ";
  cin >> flagIMWriter;
  if(flagIMWriter==1)
  {
    cout << " Escritura de la segmentaci??n en im??genes binarias... ";
    flagIMWriter = WriterSegmentation( intensityImages, liverMesh, bb, OrigDims, outIMSegFileName );
    if(flagIMWriter==0) 
      cout << "Realizada." << endl;
    else
      cout << "Error en la escritura." << endl;
  }
  int continuar;
  cout << "Continuar ? (0/1): ";
  cin >> continuar;
  if(continuar == 0)
    return 0;*/
  /**  Escritura de la segmentaci??n resultante en im??genes binarias  
  -----------------------*/
//   inicializaci??n de variables para el registro del tiempo
  time_t comienzo, final;
 /**-------------------------------------
    Inicia el proceso iterativo de deformaci??n
  -------------------------------------*/
  int i=0;
  comienzo = time( NULL );// inicio de la cuenta      
  
  do{
    cout << "iter = " << i << endl; 
    
    //  inicializaci??n de las variables para calcular la mag. de las distancias del desplazamiento de los puntos 
    lInitialMeshPoints = vtkPoints::New();
    lInitialMeshPoints->DeepCopy( liverMesh->GetPoints() );

    lMeshPoints = vtkPoints::New();
    lMeshPoints->DeepCopy( liverMesh->GetPoints() );

    numberOfPointsMesh = lMeshPoints->GetNumberOfPoints();
    cout << "PointsMesh: " << numberOfPointsMesh << endl;
    
//     fprintf( fpt, "%d, %d, ", i, numberOfPointsMesh );

    //actualizaci??n de la lista de los vectores normales lNormalVectors
    lNormalVectors = vtkDoubleArray::New();
//     lNormalVectors->SetNumberOfComponents( 3 );
//     lNormalVectors->SetNumberOfTuples( numberOfPointsMesh );
    // En este punto ya se conoce el desplazamiento de cada punto en errorDistance.
//     cout << "c??lculo de las normales..." << endl;
    normalsCalculation( liverMesh, lNormalVectors );
//     cout << "normales calculadas." << endl;
    
    //actualizaci??n de la relaci??n entre los ids de la malla con los ids de las im??genes
    lIdMeshIdImages = vtkIdTypeArray::New();
    lIdMeshIdImages->SetNumberOfValues( numberOfPointsMesh );
    
    relationIdMeshIdImages( lMeshPoints, intensityImages, lIdMeshIdImages );

    //C??lculo de la fuerza basada en la intensidad y en direcci??n de la normal (basada en la fuerza de balloon)
//     cout << "Fuerza de intensidad..." << endl;
    lIntensityForce = vtkDoubleArray::New();
    lIntensityForce->SetNumberOfComponents( 3 );
    lIntensityForce->SetNumberOfTuples( numberOfPointsMesh );
    
    intensityForce( lIntensityImages, lNormalVectors, lIdMeshIdImages, meanRegion, sigmaRegion, lIntensityForce, selFP );
//     cout << "Fuerza de intensidad calculada." << endl;
    //------------------------------
    
    /**C??lculo de la fuerza de bordes*/
//     cout << "Fuerza de bordes..." << endl;
    lEdgeForce = vtkDoubleArray::New();
    lEdgeForce->SetNumberOfComponents( 3 );
    lEdgeForce->SetNumberOfTuples( numberOfPointsMesh );
    
    edgeForce( lGradientEdgeImages, lNormalVectors, lIdMeshIdImages, lEdgeForce );
//     cout << "Fuerza de bordes calculada." << endl;
    //--------------------------------
    lNormalVectors->Delete();
    lIdMeshIdImages->Delete();

    /**C??lculo de la fuerza de tensi??n como parte de la fuerza de suavizado*/
    // c??lculo de la fuerza de suavizado en funci??n de la segunda derivada, conocida como fuerza de Tensi??n. (esta fuerza actua acercando cada punto Xi hacia el centroide de sus vecinos)
//     cout << "Fuerza de tensi??n..." << endl;
    lnbs = new vtkIdList*[numberOfPointsMesh];
    
    lcentroids = vtkDoubleArray::New();
    lcentroids->SetNumberOfComponents(3);
    lcentroids->SetNumberOfTuples( numberOfPointsMesh );
    
    lTensionForce = vtkDoubleArray::New();
    lTensionForce->SetNumberOfComponents(3);
    lTensionForce->SetNumberOfTuples( numberOfPointsMesh );
    
    findNeighbors1( liverMesh, lnbs );
    findCentroids( liverMesh, lnbs, lcentroids );
    
    delete[] lnbs;
    
    tensionForce( liverMesh, lcentroids, lTensionForce );
    
    lcentroids->Delete();
//     cout << "Fuerza de tensi??n calculada." << endl;
    //--------------------------------
    /**Movimento de los puntos de la malla por el efecto de las fuerzas de deformaci??n */  
    //extracci??n de las posiciones actuales de los v??rtices de la malla
//     cout << "Movimiento de los puntos..." << endl;
    lNewMeshPoints = vtkPoints::New();
    lNewMeshPoints->SetNumberOfPoints( numberOfPointsMesh );
    
    movePointsMesh( lMeshPoints, w1, lTensionForce, w3, lEdgeForce, w4, lIntensityForce, limitImages, lNewMeshPoints );
    liverMesh->SetPoints( lNewMeshPoints );
//     cout << "Movimiento de los puntos realizado." << endl;
    
    lNewMeshPoints->Delete();
    lMeshPoints->Delete();
    lIntensityForce->Delete();
    lEdgeForce->Delete();
    lTensionForce->Delete();
//     cout << "paso de los puntos a liverMesh..." << endl;
//     liverMesh->Reset();
//     points1->DeepCopy( tempLiverMesh->GetPoints() );
//     liverMesh->SetPoints( points1 );
//     polys1->DeepCopy( tempLiverMesh->GetPolys() );
//     liverMesh->SetPolys( polys1 );
//     cout << "paso realizado." << endl;
    
//     liverMesh->BuildCells();
//     liverMesh->BuildLinks();
//     liverMesh->Update();
    
    
    //lNewMeshPoints contiene las nuevas coordenadas de los puntos de la malla, que son asignados con SetPoints()
//     liverMesh->SetPoints( lNewMeshPoints );
//     liverMesh->BuildLinks();
//     liverMesh->Squeeze(); //  reclama memoria no usada
//     cout << "Movimiento de los puntos realizado." << endl;
    /** se evaluan las long. de los bordes de los tri??ngulos; si ent??n muy cerca se fusionan los puntos y se eliminan las celdas que comparten este borde, realizando antes de esto y, posteriormente una etapa de suavizado  */
    outputMesh = vtkPolyData::New();
    
    errorDistance = vtkDoubleArray::New();
    errorDistance->SetNumberOfValues( numberOfPointsMesh );

//     cout << "Remallado..." << endl;
    remesh( liverMesh, minimumDistance, lInitialMeshPoints, errorDistance, outputMesh );
//     cout << "Remallado realizado." << endl;
    
    lInitialMeshPoints->Delete();
    
    liverMesh->DeleteLinks();
    liverMesh->Delete();
    liverMesh = vtkPolyData::New();
    
    points = vtkPoints::New();
    polys = vtkCellArray::New();
  
  // copia de la malla de entrada (solo los polys) a la malla de trabajo  
    points->DeepCopy( outputMesh->GetPoints() );
    liverMesh->SetPoints( points );
    points->Delete();
    polys->DeepCopy( outputMesh->GetPolys() );
    liverMesh->SetPolys( polys );
    polys->Delete();
    liverMesh->GetPointData()->DeepCopy(outputMesh->GetPointData());
    liverMesh->GetFieldData()->PassData(outputMesh->GetFieldData());
    liverMesh->BuildCells();
    liverMesh->BuildLinks();
  
    outputMesh->Delete();

    
    // Funci??n que evalua el criterio de parada
    fprintf( fpt, "%d\t", i+1 );
    // medimos el tiempo de cada iteraci??n y almacenarlo en un archivo, en segundos
//     cout << "Registro desplazamiento v??rtices." << endl;
    stopEvaluation( errorDistance, fpt, lCounts, medianCounts );
    
    errorDistance->Delete();

    i++;

    
  }while( i < iter );
//   (medianCounts[0] < 0.9*numberOfPointsMesh || medianCounts[0]+medianCounts[1] < 0.98*numberOfPointsMesh) &&

  final = time( NULL );
  fprintf( fpt, "tiempo %f\t", difftime(final, comienzo) );

  fclose(fpt);
  cout << "Actualizaci??n realizada." << endl;
  /**Fin del proceso iterativo de deformaci??n
  -------------------------- */

  //visualizaci??n de la superficie del h??gado
  cout <<"\tVisualizaci??n de la malla... ";
//   MeshVisualization( liverMesh );
  cout << "realizada." << endl;
  //--------------------------
  //  visualizaci??n de los vectores normales sobre la superficie del h??gado
  cout <<"\tVisualizaci??n de las normales... ";
/*  vtkPolyDataNormals *normalsMesh = vtkPolyDataNormals::New();
  normalsMesh->SplittingOff();
  normalsMesh->SetInput( liverMesh );
  
  normalsMesh->Update();

  NormalsMeshVisualization( normalsMesh->GetOutput() );
  
  normalsMesh->Delete();*/
  cout << "realizada." << endl;
  //----------------------------


  // copy the input (only polys) to our working mesh
//   vtkPolyData *finalLiverMesh = vtkPolyData::New();
//   
//   polys->DeepCopy(liverMesh->GetPolys());
//   points->DeepCopy(liverMesh->GetPoints());
//   
//   finalLiverMesh->SetPolys(polys);
//   finalLiverMesh->SetPoints(points);
//   
//   finalLiverMesh->BuildCells();
//   finalLiverMesh->BuildLinks();
//   finalLiverMesh->Update();
//   


  /**   visualizaci??n de las curvas sobre las im??genes  */
//   cout << "Visualizaci??n de la segmentaci??n sobre las im??genes..." << endl;
//   range[0] = int( meanRegion - 3*sigmaRegion);
//   range[1] = int( meanRegion + 3*sigmaRegion);
// //    liverMesh
//   SynchronizedViews( intensityImages, liverMesh, range);
//   cout << "visualizaci??n realizada" << endl;

/**--------------
Dado que pueden generarse auto-intersecciones en la deformaci??n de la superficie, este bloque de c??digo busca obtener solo la envovente. */

//   vtkFeatureEdges *nuevaMalla = vtkFeatureEdges::New();
//   nuevaMalla->SetInput(liverMesh);
//   nuevaMalla->ColoringOff();
//   nuevaMalla->BoundaryEdgesOn();
//   nuevaMalla->ManifoldEdgesOn();
//   nuevaMalla->NonManifoldEdgesOn();	
//   nuevaMalla->Update();
/**-------------- */

//   vtkDataSetSurfaceFilter *nuevaMalla = vtkDataSetSurfaceFilter::New();
//   nuevaMalla->SetInput(liverMesh);
//   nuevaMalla->Update();

  /**   visualizaci??n de las curvas sobre las im??genes  */
  cout << "Visualizaci??n de la segmentaci??n sobre las im??genes..." << endl;
  range[0] = int( meanRegion - 3*sigmaRegion);
  range[1] = int( meanRegion + 3*sigmaRegion);
//     nuevaMalla->GetOutput()
  SynchronizedViews( intensityImages, liverMesh, range);
  cout << "visualizaci??n realizada" << endl;


 
    /**-----------------------
  Escritura de la segmentaci??n resultante en im??genes binarias  */
  
  int flagWriter=0;
  cout << "Escribir el resultado de la segmentaci??n?(0/1): ";
  cin >> flagWriter;
  if(flagWriter==1)
  {
    cout << " Escritura de la segmentaci??n resultante en im??genes binarias... ";
//     nuevaMalla->GetOutput()
    flagWriter = WriterSegmentation( intensityImages, liverMesh, bb, OrigDims, outSegFileName );
    if(flagWriter==0) 
      cout << "Realizada." << endl;
    else
      cout << "Error en la escritura." << endl;
  }
  
  /**  Escritura de la segmentaci??n resultante en im??genes binarias  
  -----------------------*/

  

  /** Visualizaci??n de la magnitud de los vectores gradiente de las im??genes de la funci??n de velocidad */
//   cout << "Visualizaci??n magnitud del gradiente de la funci??n de velocidad... ";
/*  vtkImageMagnitude *magGradientImages = vtkImageMagnitude::New();
  magGradientImages->SetInputConnection( gradientImages->GetOutputPort() );
  magGradientImages->Update();
  range[0] = 0;
  range[1] = 1;
  SynchronizedViews( magGradientImages->GetOutput(), finalLiverMesh, range);
//   ViewImages( magGradientImages->GetOutput(), range );
  magGradientImages->Delete();
  cout << "realizada." << endl;*/
  /** Visualizaci??n de la magnitud de los vectores gradiente de las im??genes de la funci??n de velocidad */
  //--------------------------
  
//   lEdgeForce->Delete();
//   lIntensityForce->Delete();
//   lIdMeshIdImages->Delete();
  gradientImages->Delete();
//   errorDistance->Delete();

//   nuevaMalla->Delete();
  liverMesh->DeleteLinks();
  liverMesh->Delete();
  
  return 0;
}

int normalsCalculation( vtkPolyData *liverMesh, vtkDataArray *lNormalVectors )
{
//   vtkPoints *points = vtkPoints::New();
//   vtkCellArray *polys = vtkCellArray::New();
//   vtkPolyData *temporalMesh = vtkPolyData::New();
//   
//   points->DeepCopy( liverMesh->GetPoints() );
//   temporalMesh->SetPoints( points );
//   points->Delete();
//   
//   polys->DeepCopy( liverMesh->GetPolys() );
//   temporalMesh->SetPolys( polys );
//   polysMesh->Delete();
//   
//   temporalMesh->BuildCells();
//   temporalMesh->BuildLinks();

  vtkPolyDataNormals *normalsMesh = vtkPolyDataNormals::New();
  normalsMesh->SplittingOff();
  normalsMesh->SetInput( liverMesh );
  normalsMesh->Update();

//   cout << "normalspts: " << normalsMesh->GetOutput()->GetNumberOfPoints() << endl;
  lNormalVectors->DeepCopy( normalsMesh->GetOutput()->GetPointData()->GetNormals() );

  normalsMesh->Delete();
//   temporalMesh->DeleteLinks();
//   temporalMesh->Delete();

  return 0;
}
