//En esta parte del programa se abre meshAtlasLiver.stl, se transforma con la información translate[] y matrixTransform[], se visualizan las curvas sobre las imágenes, y se deforma la susperficie. 
 
// Autor: Gerardo Tibamoso Pedraza. 

// librerias ITK
#include "itkCastImageFilter.h"
#include "itkImage.h"
#include "itkImageFileReader.h"
#include "itkImageFileWriter.h"
#include "itkMedianImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"
#include "itkSobelEdgeDetectionImageFilter.h"

// librerias VTK
#include "vtkImageData.h"
#include "vtkSTLReader.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"

// librerias vtkINRIA3D
#include "itkImageToVTKImageFilter.h"

// librerias básicas de C++
#include <string.h>
#include <iostream>
using namespace std;

typedef itk::Image< unsigned char, 3 >   UCharImagesType;
typedef itk::Image< short, 3 >   ShortImagesType;
typedef itk::Image< double, 3 >  DoubleImagesType;


int SynchronizedViews( vtkImageData *imagesCT, vtkPolyData *meshLiverModelTransformed, int range[2] );

int metodoDeDeformacion (vtkPolyData *inputLiverMesh, vtkImageData *intensityImages, double *infoIntensity, vtkImageData *magnitudeSobelImages, double *infoMagGradient, int iter, int *OrigDims, int *bb, char outSegFileName[], double cs, double cp, double cb, int selFP);
//char outIMSegFileName[]

/**--------------------------------------
Esta función es llamada por MetodoCompleto.cxx  */ 
int InicializacionDeformacion( vtkImageData *medianImages, vtkImageData *magGradImages, vtkPolyData *meshLiver, double translate1[3], double translate2[3], double matrixTransform[16], int *dims, int *bb, double *infoIntensity, double *infoMagGradient, int iter, char outSegFileName[], double cs, double cp, double cb, int selFP )
{
  
  cout << "iter = " << iter << endl;
  cout << outSegFileName << endl;
  
 // -----------------------------------------  
 
  /**---------------------------
  Transformación de la malla... */
  
  cout << "Transformación de la malla de acuerdo al registro manual realizado previamente sobre las imágenes... " << endl;
  
  vtkTransform *transform1 = vtkTransform::New();
  transform1->Translate( translate1 );
  
  vtkTransform *transform2 = vtkTransform::New();
  transform2->SetMatrix( matrixTransform );
  
  vtkTransform *transform3 = vtkTransform::New();
  transform3->Translate( -translate2[0], -translate2[1], -translate2[2] );

  vtkTransformPolyDataFilter *transformPDFilter1 = vtkTransformPolyDataFilter::New();
  transformPDFilter1->SetTransform( transform1 );

  vtkTransformPolyDataFilter *transformPDFilter2 = vtkTransformPolyDataFilter::New();
  transformPDFilter2->SetTransform( transform2 );

  vtkTransformPolyDataFilter *transformPDFilter3 = vtkTransformPolyDataFilter::New();
  transformPDFilter3->SetTransform( transform3 );

  transformPDFilter1->SetInput( meshLiver );
  transformPDFilter2->SetInputConnection( transformPDFilter1->GetOutputPort() );
  transformPDFilter3->SetInputConnection( transformPDFilter2->GetOutputPort() );
  
  transformPDFilter3->Update();
  
  transformPDFilter1->Delete();
  transform1->Delete();
  transformPDFilter2->Delete();
  transform2->Delete();
  
  cout << "realizada." << endl;
  
  /**Transformación de la malla
  ----------------------------------------*/
  /**-------------------------
  Visualización de las curvas de la malla sobre las imágenes */

  int range[2];
  
  range[0] = int(infoIntensity[0] - 5.0*infoIntensity[1]);
  range[1] = int(infoIntensity[0] + 5.0*infoIntensity[1]);
  SynchronizedViews( medianImages, transformPDFilter3->GetOutput(), range);

  /**Visualización realizada.
  ----------------------------*/  
  
  /**-------------------------
  Visualización de las curvas de la malla sobre las imágenes */

//   range[0] = 0;
//   range[1] = int(infoMagGradient[0] + 3*infoMagGradient[1]);
//   SynchronizedViews( connectorSobelFilter->GetOutput(), transformPDFilter2->GetOutput(), range);

  /**Visualización realizada.
  ----------------------------*/  

  cout << "infoIntensity (media y desviacion estandar): " << infoIntensity[0] << ", " << infoIntensity[1] << endl;
  cout << "infoMagGradient (mediana): " << infoMagGradient[2] << ", " << infoMagGradient[1] << endl; 

  /**----------------------------
  Proceso de la deformación  */ 
  
  //   visualización de las curvas sobre las imágenes  

  int flagDef=0;
  
  cout << "Realizar deformación?(0/1): ";
  cin >> flagDef;
  if(flagDef == 1)
  {
    metodoDeDeformacion( transformPDFilter3->GetOutput(), medianImages, infoIntensity, magGradImages, infoMagGradient, iter, dims, bb, outSegFileName, cs, cp, cb, selFP );
  }
  
  /**Proceso de deformación
  ------------------------------*/
  
  transformPDFilter3->Delete();
  transform3->Delete();

  cout << "Salida Normal." << endl;
  return EXIT_SUCCESS;

}
