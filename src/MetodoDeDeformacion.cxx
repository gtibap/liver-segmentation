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
Inicia el programa para la deformación de la malla. 
Esta función es llamada por la función deformationInit()
 */
// int deformation (vtkPolyData *initialLiver, vtkPolyData *initialLiver1, vtkImageData *binaryIntensityImages, vtkImageData *magnitudeSobelImages, int iter ) {
int metodoDeDeformacion (vtkPolyData *inputLiverMesh, vtkImageData *intensityImages, double *infoIntensity, vtkImageData *magnitudeSobelImages, double *infoMagGradient, int iter, int *OrigDims, int *bb, char outSegFileName[], double cs, double cp, double cb, int selFP ) {
  
//   cout << "Deformación de la malla por la ecuación de movimiento... " << endl;
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
//   double p = 1.0, q = 1.0; // p y q son los factores de ponderación de las fuerzas de balloon y edge, respectivamente.
  // -----------------------------------------  
  
  /** Remallado de la malla que representa el hígado en función del tamaño del voxel. En este caso, la distancia mínima de los bordes de la malla será la magnitud de la diagonal de cada voxel en las imágenes */
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
  

  /**  Visualización de la malla. */
  cout << "Visualización de la malla copiada..." << endl;
//   MeshVisualization( liverMesh );
  cout << "Visualización realizada." << endl;
  /**  Visualización de la malla.*/
 
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
  Estas listas de puntos servirán para calcular la distancia de desplazamiento de los puntos de la malla en cada iteración de la deformación de la malla  */
  
//   vtkPoints *initialPoints = vtkPoints::New();
//   vtkIdType numberOfPointsMesh;
  
//   errorDistance->SetNumberOfComponents(1);
  /**------------------------------------
  Estas listas de puntos servirán para calcular la distancia de desplazamiento de los puntos de la malla en cada iteración de la deformación de la malla  */

/**---------------------------
  Visualización de la malla... */
//   MeshVisualization( liverMesh );
  /**  Visualización de la malla.
  ---------------------------*/
  // -----------------------------------------
  cout << "Deformación de la malla por la ecuación de movimiento, la cuál se compone básicamente de una fuerza de rigidez de la malla y una fuerza de acople a las imágenes ... " << endl;
  /** La ecuación de movimiento de los puntos de la malla, es una versión discretizada de la ecuación diferencial de primer orden definada por:
   dX
  ---- = FuerzaDeRigidez + FuerzaDeAcople
   dt
  La fuerza de rigidez busca acercar cada punto al centroide de sus vecinos mientras preserva la forma inicial.
  La fuerza de acople busca llevar cada punto hacia los bordes de las imágenes tomando en cuenta los valores de intensidad de los pixeles más cercanos a la posición de los puntos de la malla.
  */
  
  //-----------------------------  
  cout << "Cálculo de la fuerza de acople a las imágenes... " << endl;
  /**Se comienza calculando la fuerza de acople (que estará en dirección del vector normal asociado a cada vértice), la cual se compone de dos términos:
  FuerzaDeAcople = FuerzaDeBordes + FuerzaDeIntensidad */
  
  /**Cálculo de la FuerzaDeBordes:
                                         1
  FuerzaDeBordes = [ Grad(  ------------------------------) ° normalVector ] normalVector
                              1 + (MagnitudGradiente)^p
  
  Donde la MagnitudGradiente es la magnitud del resultado de la aplicación del operador Sobel (magnitudeSobelImages), lo cual entrega la información de bordes de las imágenes suavizadas; p es un factor de pondración (valdrá 1 o 2), Grad es el operador gradiente, (°) representa el producto punto, y, normalVector es el vector normal asociado con cada vértice de la malla.
  */
  cout << "1. Cálculo de la fuerza de bordes... " << endl;
  //dado que el grad (1 / (1 + y*magnitudeSobelImages)) es siempre constante, se calcula primero y, por supuesto, una sola vez 
  cout << "1.1 Cálculo del gradiente del inverso de uno más la magnitud del gradiente elevado a p... " << endl;  
  //El primer paso para generar las imágenes de la FuerzaDeBordes es localizar el espacio de memoria del mismo tamaño de las imágenes de entrada "magnitudeSobelImages" (es necesario generar las imágenes para la aplicación del operador gradiente).

  /**-----------------------------------*/
  
  // en esta parte se visualizan las imágenes de sobel
  cout << "Visualización del resultado del operador Sobel sobre las Imágenes... " << endl;
//   range[0] = int(infoMagGradient[0] + infoMagGradient[1]);
/*  range[0] = 0;
  range[1] = int(infoMagGradient[0] + 3*infoMagGradient[1]);
  ViewImages( magnitudeSobelImages, range );*/
  cout << "Visualización realizada." << endl;
  
  /**-----------------------------------*/
  
  double meanSobel, sigmaSobel, medianSobel;
  
  /**----------------------------------------
  Información estadística (media y desviación estandar) de la magnitud del gradiente (calculada por la aplicación del operador Sobel) de un solo corte segmentado manualmente */
  
  cout << "Información de media y desviación estandar de la magnitud del resultado del operador sobel sobre la región de las imágenes encerrada por malla registrada... " << endl; 
//   statisticRegionImagesClosedMesh( magnitudeSobelImages, liverMesh, meanSobel, sigmaSobel );
  meanSobel = infoMagGradient[0];
  sigmaSobel = infoMagGradient[1];
  medianSobel = infoMagGradient[2];
  cout << "meanSobel = " << meanSobel << endl;
  cout << "sigmaSobel = " << sigmaSobel << endl;
  cout << "medianSobel = " << medianSobel << endl;
  //Esta información estadística permite hacer un filtrado en g(I) de los datos de magnitudeSobelImages para la estimación del gradiente
  
  /**Esta información servirá para eliminar valores de gradiente pequeños en las imágenes.
  ----------------------------------------*/

  /**Cálculo de la función de velocidad g(I) = 1 / (1 + magSobel^p), donde la magSobel depende de la información estadística estimada, de tal forma que serán cero para los valores de magnitudeSobelImages menores a la media, multiplicados por valores entre 0 y 1 de forma lineal para los valores entre la media y la media + la desviación estandar, y se mantiene el valor en otro caso.*/
  cout << "Cálculo de la función de velocidad g(I) la cual esta en función del gradiente (operador sobel)..." << endl;
  
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
  Definición coeficientes de ponderación fuerzas de deformación (cp y cb) **/
  double mdv;
  mdv = spacing[0];
  if (mdv > spacing[1])
    mdv = spacing[1];
  if(mdv > spacing[2])
    mdv = spacing[2];
  cout << "minima distancia entre voxeles: " << mdv << endl;
  /** Definición coeficientes de ponderación fuerzas de deformación (cp y cb) 
  -------------------------**/
  velocityFunctionImages->AllocateScalars();
  //en este punto ya se tiene el espacio en memoria para el volumen de las imágenes, para escribir los datos de la función de velocidad g(I)
  
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
  //aquí se define una función de ponderación para eliminar valores pequeños de la magnitud del gradiente obtenidos por el operador sobel, la cual depende de la media y la desviación estandar de estos valores, tal que: Son cero los valores menores que la media, se multiplican por un valor entre cero y 1 (de forma lineal) si los valores están entre la media y la media+sigma, y, son los mismos (multiplicados por 1) para los valores mayores a la media+sigma.
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
  cout << "Función de velocidad g(I) calculada." << endl;
  
  /**--------------
  Guardar las imágenes de la velocidad g(I)*/
/*  cout << "Guardando los resultados de las imágenes de velocidad..." << endl;
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
  /**Guardar las imágenes de la velocidad g(I)
  -------------*/
  /**-----------
  Visualización función de velocidad */
  cout << "Visualización de la función de velocidad g(I)... ";
/*  range[0] = 0;
  range[1] = 1;
  ViewImages( velocityFunctionImages, range );*/
  cout << "realizada." << endl;
  /** Visualización función de velocidad */
  
  cout << "lVelocityFunctionImages components: " << lVelocityFunctionImages->GetNumberOfComponents() << endl;
  cout << "lVelocityFunctionImages tuples: " << lVelocityFunctionImages->GetNumberOfTuples() << endl;
  
  //luego de calcular g(I) = 1 / (1 + y*magnitudSobelImages), se calculará el gradiente de este resultado
  cout << "Cálculo del gradiente de la función de velocidad g(I)... " << endl;
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
  Guardar las imágenes de la velocidad g(I)*/
/*  vtkImageMagnitude *magGradientImages = vtkImageMagnitude::New();
  magGradientImages->SetInputConnection( gradientImages->GetOutputPort() );

  cout << "Guardando los resultados de las imágenes del gradiente de las imágenes de velocidad..." << endl;
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
  /**Guardar las imágenes de la velocidad g(I)
  -------------*/
  
  /** Visualización de la magnitud de los vectores gradiente de las imágenes de la función de velocidad */
  cout << "Visualización magnitud del gradiente de la función de velocidad... ";
/*  vtkImageMagnitude *magGradientImages = vtkImageMagnitude::New();
  magGradientImages->SetInputConnection( gradientImages->GetOutputPort() );
  magGradientImages->Update();
  range[0] = 0;
  range[1] = 1;
  ViewImages( magGradientImages->GetOutput(), range );
  magGradientImages->Delete();*/
  cout << "realizada." << endl;
  /** Visualización de la magnitud de los vectores gradiente de las imágenes de la función de velocidad */
  
  //-----------------------------
  /**Los valores de los vectores del gradiente de las imágenes se pasan a la lista
  !!!   lGradientEdgeImages   ¡¡¡
  */  
  //se pasan los valores escalares del resultado del gradiente a la lista lGradientEdgeImages
  vtkDataArray *lGradientEdgeImages = gradientImages->GetOutput()->GetPointData()->GetScalars();
  cout << "tuples gradient = " << lGradientEdgeImages->GetNumberOfTuples() << endl;
  cout << "components gradient = " << lGradientEdgeImages->GetNumberOfComponents() << endl;
  //-----------------------------
  
  /**
  Cálculo de la FuerzaDeIntensidad
  La Fuerza de Intensidad estará en dirección normal, y ponderada por la estimación de la función gausiana que representa la distribución de intensidad del volumen del hígado. Esta ponderación estará dada por
  
  p = exp( -(x-b)^2 / (2*c^2) ),
  
  donde b es la media (mean) y c es la desviación estandar (sigma), de la distribución gausiana de probabilidad de la intensidad.
  
  Se cuenta con medidas "promedio" de la media y de la desviación estandar estimadas a partir de más de 20 hígados segmentados manualmente; y además, puede extraerse información de la media y de la desviación estandar directamente sobre el volumen de las imágenes que encierra la malla promedio que ha sido previamente registrada manualmente.
  */
  
  /** La fuerza basada en las intensidades se calculará iterativamente para cada vértice de la malla. La ponderación de esta fuerza dependerá de la media y la desviación estandar de la distribución de probabilidad de la intensidad en las imágenes
  */

  /**----------------------------------------
  Información estadística (media y desviación estandar) de las intensidades de las imágenes (luego del filtro de mediana) de un solo corte segmentado manualmente */

  cout << "Información de la distribución de probabilidad de las intensidades del hígado, a partir de la segmentación manual de éste en un corte axial..." << endl;
  double meanRegion, sigmaRegion, medianRegion;
  //  aquí se usan las clases vtkPolyDataToImageStencil y vtkImageAcumulate, y tal vez vtkImageStencil
//   statisticRegionImagesClosedMesh( intensityImages, liverMesh, meanRegion, sigmaRegion );
  
  meanRegion = infoIntensity[0];
  sigmaRegion = infoIntensity[1];
  medianRegion = infoIntensity[2];
  
  cout << "meanRegion = " << meanRegion << endl;
  cout << "sigmaRegion  = " << sigmaRegion << endl;
  cout << "medianRegion  = " << medianRegion << endl;
  /** Esta información es utilizada para definir la función de fuerza basada en la intensidad (del modelo deformable)
  ----------------------------------------*/
  
  /** referencia a la lista de los valores de intensidad de las imágenes por el puntero lIntensityImages
  */
  vtkDataArray *lIntensityImages = intensityImages->GetPointData()->GetScalars();
  cout << "intensity tuples: " << lIntensityImages->GetNumberOfTuples() << endl;
  cout << "intensity components: " << lIntensityImages->GetNumberOfComponents() << endl;
  
  // -----------------------------------------  
  
  /**Ahora se encontrará la relación entre los ids de los vértices de la malla con los ids más cercanos a estos puntos sobre las imágenes. Este proceso también debe hacerse recursivamente (es decir, cada vez que son desplazados los vértices de la malla). La relación se referencia a la lista
  lIdMeshIdImages
   */
    //  Paso de las coordenadas de los puntos de la malla a la lista referenciada por lMeshPoints  
//   lMeshPoints = vtkPoints::New();
//   lMeshPoints->DeepCopy( liverMesh->GetPoints() );
  //Referencia realizada
//   cout << "points Mesh= " << lMeshPoints->GetNumberOfPoints() << endl;

//   cout << "Relación entre los ids de la malla con los ids de las imágenes... ";
  
//   lIdMeshIdImages = vtkIdTypeArray::New();
//   lIdMeshIdImages->SetNumberOfComponents( 1 );

//   relationIdMeshIdImages( lMeshPoints, intensityImages, lIdMeshIdImages );
//   cout << "realizada." << endl;
  /**Realizada la relación entre los ids de la malla con los de las imágenes*/

  // -----------------------------------------  
  /**Definición de las variables, cálculo iterativo de las fuerzas de deformación y de las nuevas posiciones de cada vértice de la malla*/ 
  cout << "Definición de las variables necesarias adicionales... " << endl;
 
  //   Cálculo de la fuerza de Balloon, la cúal es la encargada de mover los vértices de la malla en función del vector normal unitario (asociado a cada vértice) y al rango de intensidades del hígado, para que la malla se mueva hacia adentro o hacia afuera, acercandose a los bordes de este órgano.
  
  // Se inicializa una lista para almacenar los valores que representan los vectores de la fuerza de Balloon
  vtkDoubleArray *lIntensityForce;
  vtkDoubleArray *lEdgeForce;
  //inicialización de variables
//   vtkIdType nop = liverMesh->GetNumberOfPoints();
//   cout << "nop= " << nop << endl;
  
//   vtkIdList **lnbs = new vtkIdList*[nop];
  vtkIdList **lnbs;
  vtkDoubleArray *lcentroids;
  vtkDoubleArray *lTensionForce;
  
/*  vtkDoubleArray *lRigidityForce;
  lRigidityForce = vtkDoubleArray::New();
  lRigidityForce->SetNumberOfComponents(3);*/
  
  cout << "número de iteraciones: " << iter << endl;
  
  int lCounts[3][3] = {}; //  Cuentas se inicializa en ceros
  int medianCounts[3];

  
  double w1, w3, w4; //valores de ponderación de cada una de las fuerzas de deformación
  w1= cs;// 0.7 ponderación fuerza de tensión
//   double  w2=-1.0;//ponderación fuerza de rigidez
  w3=-cb*mdv;  //  0.9mm ponderación fuerza de bordes
  w4= cp*mdv;  //  0.5mm ponderación fuerza de intensidad (basada en la fuerza de balloon)
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
  fprintf( fpt, "iter\t numVertices\t Número de vertices de acuerdo con sus desplazamientos entre [0, 1] con pasos de 0.1 (primera columna, vertices que se desplazan entre 0 y 0.1 mm, en la segunda, los que se desplazan entre 0.1 y 0.2, y así hasta la columna 10, de los vértices que se desplazan entre 0.9 a 1.0 )\n");
  /** realizada la escritura de los datos obtenidos en un archivo de texto 
  ------------------------------------------------------*/
  // Se obtiene la información de los límites de las imágenes 
  double originImages[3], spacingImages[3], limitImages[6];
  int dimensionImages[3];
  
  intensityImages->GetOrigin(originImages);
  cout << "originImages: " << originImages[0] <<", " << originImages[1] << ", " << originImages[2] << endl;
  
  intensityImages->GetSpacing(spacingImages);
  cout << "spacingImages: " << spacingImages[0] <<", " << spacingImages[1] << ", " << spacingImages[2] << endl;

  intensityImages->GetDimensions(dimensionImages);
  cout << "dimensionImages: " << dimensionImages[0] <<", " << dimensionImages[1] << ", " << dimensionImages[2] << endl;
  
  // se definen los límites de las imágenes
  limitImages[0] = originImages[0];
  limitImages[1] = originImages[0] + spacingImages[0]*(dimensionImages[0] - 1.5);
  limitImages[2] = originImages[1];
  limitImages[3] = originImages[1] + spacingImages[1]*(dimensionImages[1] - 1.5);
  limitImages[4] = originImages[2];
  limitImages[5] = originImages[2] + spacingImages[2]*(dimensionImages[2] - 1.5);
  
  // si la malla tiene puntos ubicados fuera de las imágenes, entonces se llevan a los límites
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
  Visualización de la malla... */
//   MeshVisualization( liverMesh );
  /**  Visualización de la malla.
  ---------------------------*/
  
//   Antes del proceso de deformación es necesario guardar el volumen inicial encerrado por la superficie manualmente ajustada. Para esto se utiliza el procedimiento para escribir las imágenes segmentadas a partir de la superficie incial.
  /**-----------------------
  Escritura de la segmentación resultante de la superficie del hígado "promedio" en imágenes binarias  */
/*  int flagIMWriter=0;
  cout << "Escribir la segmentación de la superficie inicial?(0/1): ";
  cin >> flagIMWriter;
  if(flagIMWriter==1)
  {
    cout << " Escritura de la segmentación en imágenes binarias... ";
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
  /**  Escritura de la segmentación resultante en imágenes binarias  
  -----------------------*/
//   inicialización de variables para el registro del tiempo
  time_t comienzo, final;
 /**-------------------------------------
    Inicia el proceso iterativo de deformación
  -------------------------------------*/
  int i=0;
  comienzo = time( NULL );// inicio de la cuenta      
  
  do{
    cout << "iter = " << i << endl; 
    
    //  inicialización de las variables para calcular la mag. de las distancias del desplazamiento de los puntos 
    lInitialMeshPoints = vtkPoints::New();
    lInitialMeshPoints->DeepCopy( liverMesh->GetPoints() );

    lMeshPoints = vtkPoints::New();
    lMeshPoints->DeepCopy( liverMesh->GetPoints() );

    numberOfPointsMesh = lMeshPoints->GetNumberOfPoints();
    cout << "PointsMesh: " << numberOfPointsMesh << endl;
    
//     fprintf( fpt, "%d, %d, ", i, numberOfPointsMesh );

    //actualización de la lista de los vectores normales lNormalVectors
    lNormalVectors = vtkDoubleArray::New();
//     lNormalVectors->SetNumberOfComponents( 3 );
//     lNormalVectors->SetNumberOfTuples( numberOfPointsMesh );
    // En este punto ya se conoce el desplazamiento de cada punto en errorDistance.
//     cout << "cálculo de las normales..." << endl;
    normalsCalculation( liverMesh, lNormalVectors );
//     cout << "normales calculadas." << endl;
    
    //actualización de la relación entre los ids de la malla con los ids de las imágenes
    lIdMeshIdImages = vtkIdTypeArray::New();
    lIdMeshIdImages->SetNumberOfValues( numberOfPointsMesh );
    
    relationIdMeshIdImages( lMeshPoints, intensityImages, lIdMeshIdImages );

    //Cálculo de la fuerza basada en la intensidad y en dirección de la normal (basada en la fuerza de balloon)
//     cout << "Fuerza de intensidad..." << endl;
    lIntensityForce = vtkDoubleArray::New();
    lIntensityForce->SetNumberOfComponents( 3 );
    lIntensityForce->SetNumberOfTuples( numberOfPointsMesh );
    
    intensityForce( lIntensityImages, lNormalVectors, lIdMeshIdImages, meanRegion, sigmaRegion, lIntensityForce, selFP );
//     cout << "Fuerza de intensidad calculada." << endl;
    //------------------------------
    
    /**Cálculo de la fuerza de bordes*/
//     cout << "Fuerza de bordes..." << endl;
    lEdgeForce = vtkDoubleArray::New();
    lEdgeForce->SetNumberOfComponents( 3 );
    lEdgeForce->SetNumberOfTuples( numberOfPointsMesh );
    
    edgeForce( lGradientEdgeImages, lNormalVectors, lIdMeshIdImages, lEdgeForce );
//     cout << "Fuerza de bordes calculada." << endl;
    //--------------------------------
    lNormalVectors->Delete();
    lIdMeshIdImages->Delete();

    /**Cálculo de la fuerza de tensión como parte de la fuerza de suavizado*/
    // cálculo de la fuerza de suavizado en función de la segunda derivada, conocida como fuerza de Tensión. (esta fuerza actua acercando cada punto Xi hacia el centroide de sus vecinos)
//     cout << "Fuerza de tensión..." << endl;
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
//     cout << "Fuerza de tensión calculada." << endl;
    //--------------------------------
    /**Movimento de los puntos de la malla por el efecto de las fuerzas de deformación */  
    //extracción de las posiciones actuales de los vértices de la malla
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
    /** se evaluan las long. de los bordes de los triángulos; si entán muy cerca se fusionan los puntos y se eliminan las celdas que comparten este borde, realizando antes de esto y, posteriormente una etapa de suavizado  */
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

    
    // Función que evalua el criterio de parada
    fprintf( fpt, "%d\t", i+1 );
    // medimos el tiempo de cada iteración y almacenarlo en un archivo, en segundos
//     cout << "Registro desplazamiento vértices." << endl;
    stopEvaluation( errorDistance, fpt, lCounts, medianCounts );
    
    errorDistance->Delete();

    i++;

    
  }while( i < iter );
//   (medianCounts[0] < 0.9*numberOfPointsMesh || medianCounts[0]+medianCounts[1] < 0.98*numberOfPointsMesh) &&

  final = time( NULL );
  fprintf( fpt, "tiempo %f\t", difftime(final, comienzo) );

  fclose(fpt);
  cout << "Actualización realizada." << endl;
  /**Fin del proceso iterativo de deformación
  -------------------------- */

  //visualización de la superficie del hígado
  cout <<"\tVisualización de la malla... ";
//   MeshVisualization( liverMesh );
  cout << "realizada." << endl;
  //--------------------------
  //  visualización de los vectores normales sobre la superficie del hígado
  cout <<"\tVisualización de las normales... ";
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


  /**   visualización de las curvas sobre las imágenes  */
//   cout << "Visualización de la segmentación sobre las imágenes..." << endl;
//   range[0] = int( meanRegion - 3*sigmaRegion);
//   range[1] = int( meanRegion + 3*sigmaRegion);
// //    liverMesh
//   SynchronizedViews( intensityImages, liverMesh, range);
//   cout << "visualización realizada" << endl;

/**--------------
Dado que pueden generarse auto-intersecciones en la deformación de la superficie, este bloque de código busca obtener solo la envovente. */

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

  /**   visualización de las curvas sobre las imágenes  */
  cout << "Visualización de la segmentación sobre las imágenes..." << endl;
  range[0] = int( meanRegion - 3*sigmaRegion);
  range[1] = int( meanRegion + 3*sigmaRegion);
//     nuevaMalla->GetOutput()
  SynchronizedViews( intensityImages, liverMesh, range);
  cout << "visualización realizada" << endl;


 
    /**-----------------------
  Escritura de la segmentación resultante en imágenes binarias  */
  
  int flagWriter=0;
  cout << "Escribir el resultado de la segmentación?(0/1): ";
  cin >> flagWriter;
  if(flagWriter==1)
  {
    cout << " Escritura de la segmentación resultante en imágenes binarias... ";
//     nuevaMalla->GetOutput()
    flagWriter = WriterSegmentation( intensityImages, liverMesh, bb, OrigDims, outSegFileName );
    if(flagWriter==0) 
      cout << "Realizada." << endl;
    else
      cout << "Error en la escritura." << endl;
  }
  
  /**  Escritura de la segmentación resultante en imágenes binarias  
  -----------------------*/

  

  /** Visualización de la magnitud de los vectores gradiente de las imágenes de la función de velocidad */
//   cout << "Visualización magnitud del gradiente de la función de velocidad... ";
/*  vtkImageMagnitude *magGradientImages = vtkImageMagnitude::New();
  magGradientImages->SetInputConnection( gradientImages->GetOutputPort() );
  magGradientImages->Update();
  range[0] = 0;
  range[1] = 1;
  SynchronizedViews( magGradientImages->GetOutput(), finalLiverMesh, range);
//   ViewImages( magGradientImages->GetOutput(), range );
  magGradientImages->Delete();
  cout << "realizada." << endl;*/
  /** Visualización de la magnitud de los vectores gradiente de las imágenes de la función de velocidad */
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
