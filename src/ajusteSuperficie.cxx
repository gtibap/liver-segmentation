/**-----------------------
En este programa se ajusta manualmente la superficie del hígado promedio, sobre las imágenes de TAC, por medio de otra superficie.
--------------------------**/
  
// Autor: Gerardo Tibamoso Pedraza. 

// librerias VTK
// #include "vtkPolyData.h"
#include "vtkSTLReader.h"

// librerias vtkINRIA3D
#include "itkImageToVTKImageFilter.h"


// librerias básicas de C++
#include <string.h>
#include <stdio.h>
#include <time.h>
#include <iostream>
using namespace std;


// int SeleccionRegionHigado( ShortImagesType::Pointer  , double *infoIntensity, char apertureFileName[], char meshRegionFileName[], vtkPolyData *meshLiver, double translate1[3], double translate2[3], double matrixTransform[16]  );

// int initialRemesh( vtkPolyData *liverMesh, float minimumDistance );

void MeshVisualization( vtkPolyData *mesh );

int AjusteManual( vtkPolyData *meshRegion, vtkPolyData *meshLiver, double translate1[3], double translate2[3], double matrixTransform[16] );

// int InicializacionDeformacion( vtkImageData *medianImages, vtkImageData *magGradImages, vtkPolyData *meshLiver, double translate1[3], double translate2[3], double matrixTransform[16], int *dims, int *bb, double *infoIntensity, double *infoMagGradient, int iter, char outSegFileName[] );
/**---------------------------------
Función principal  */ 
int ajusteSuperficie( char * meshLiverFileName, char * meshRegionFileName, double * translate1, double * matrixTransform, double * translate2 )
{
  
  /**-- lectura del modelo promedio del hígado ----*/
  cout << "Lectura de la superficie del hígado promedio... " << endl;

  vtkSTLReader *readerMeshLiver = vtkSTLReader::New();
  readerMeshLiver->SetFileName( meshLiverFileName );
  try{
    readerMeshLiver->Update();
  }
  catch( exception& ){
    cerr << "Error en la lectura de la superficie." << endl;
    return EXIT_FAILURE;
  }
  cout << "Lectura realizada." << endl;
  /**-- lectura del modelo promedio del hígado ----*/
  
  /**-- lectura de la superficie de la región del hígado de interés ----*/
  cout << "Lectura de la superficie de la región del hígado de interés... " << endl;

  vtkSTLReader *readerMeshRegion = vtkSTLReader::New();
  readerMeshRegion->SetFileName( meshRegionFileName );
  try{
    readerMeshRegion->Update();
  }
  catch( exception& ){
    cerr << "Error en la lectura de la superficie." << endl;
    return EXIT_FAILURE;
  }
  cout << "Lectura realizada." << endl;
  /**-- lectura del modelo promedio del hígado ----*/

  /** ********* **/
  /**-- Ubicación manual de la superficie promedio del hígado sobre las imágenes--*/

  AjusteManual( readerMeshRegion->GetOutput(), readerMeshLiver->GetOutput(), translate1, translate2,  matrixTransform );
  
  /**-- Ubicación manual de la superficie promedio del hígado sobre las imágenes--*/  
  /** ********* **/

  
  //   int flagRegion = SeleccionRegionHigado( medianImages, infoIntensity, apertureFileName, meshRegionFileName, readerMeshLiver->GetOutput(), translate1, translate2, matrixTransform );

  
  readerMeshRegion->Delete();
  readerMeshLiver->Delete();
  
  cout << "Salida Normal." << endl;
  return EXIT_SUCCESS;

}
