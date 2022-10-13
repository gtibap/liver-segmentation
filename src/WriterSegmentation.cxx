// En este programa se escriben las imágenes del resultado de la segmentación realizada por la malla que representa el hígado

//  librerías ITK
#include "itkAffineTransform.h"
#include "itkIdentityTransform.h"
#include "itkImage.h"
#include "itkImageFileWriter.h"
#include "itkLinearInterpolateImageFunction.h"
#include "itkNearestNeighborInterpolateImageFunction.h"
#include "itkResampleImageFilter.h"

  
// Librerías VTK
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkImageConstantPad.h"
#include "vtkImageData.h"
#include "vtkImagePadFilter.h"
#include "vtkImageReslice.h"
#include "vtkImageStencil.h"
#include <vtkImageTranslateExtent.h>
#include "vtkMetaImageWriter.h"
#include "vtkPointData.h"
#include "vtkPoints.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkPolyDataToImageStencil.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkUnsignedCharArray.h"

//  Librerías vtkINRIA3D
#include "itkVTKImageToImageFilter.h"


#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
// #include "vtkPolyDataWriter.h"
// #include "vtkSTLReader.h"

#include <iostream>
using namespace std;

typedef itk::Image< unsigned char, 3 > UCharImagesType;

/**
Inicia el programa para la escritura de las imágenes que representan las segmentaciones realizada por el método propuesto.
 */
int WriterSegmentation( vtkImageData *segImages, vtkPolyData *meshLiver, int *regionBB, int sizeOrig[3], char *outSegFileName ){

  int dims[3];
  double spacing[3], origin[3];
  
  
  vtkIdType nov; //nov number of voxels
  vtkImageData *newSegImages = vtkImageData::New();
  
  newSegImages->SetScalarTypeToUnsignedChar();
  newSegImages->SetNumberOfScalarComponents( 1 );
  
  segImages->GetDimensions( dims );
  segImages->GetSpacing( spacing );
  segImages->GetOrigin( origin );
  
  newSegImages->SetDimensions( dims );
  newSegImages->SetSpacing( spacing );
  newSegImages->SetOrigin( origin );
  
//   cout << "dims: " << dims[0] << ", " << dims[1] << ", " << dims[2] << endl;
//   cout << "spacing: " << spacing[0] << ", " << spacing[1] << ", " << spacing[2] << endl;
//   cout << "origin: " << origin[0] << ", " << origin[1] << ", " << origin[2] << endl;

  newSegImages->AllocateScalars();
  //en este punto ya se tiene el espacio en memoria para el volumen de las imágenes, para escribir los datos de la segmentación
  
  // se almacenan los valores de magnitudeSobelImages en la lista lMagnitudeSobelImages
  nov = segImages->GetPointData()->GetScalars()->GetNumberOfTuples();
//   cout << "nov: " << nov << endl;
  
  //se define la lista lVelocityFunctionImages para almacenar el resultado del inverso del gradiente
  vtkUnsignedCharArray *lNewSegImages = vtkUnsignedCharArray::New();
  lNewSegImages->SetNumberOfValues( nov );

  //  aquí se escriben las imágenes con 1's
  for( int idVoxel=0; idVoxel < nov; ++idVoxel ){
    lNewSegImages->SetValue( idVoxel, 1 );
  }

  newSegImages->GetPointData()->SetScalars( lNewSegImages );
  newSegImages->Modified();
  newSegImages->Update();
  
  /** Ya preparada las imágenes en blanco, ahora se deja unicamente la región que encierra la malla */
  

  vtkPolyDataNormals *normalsMesh = vtkPolyDataNormals::New();
  normalsMesh->SetInput( meshLiver );
  
  vtkPolyDataToImageStencil *dataToStencil = vtkPolyDataToImageStencil::New();
  dataToStencil->SetInformationInput( newSegImages );
  dataToStencil->SetInput( normalsMesh->GetOutput() );
  
  //---------------------------
  //extracción de la región que encierra la malla
  
  vtkImageStencil *stencil = vtkImageStencil::New();
  stencil->SetInput( newSegImages );
  stencil->SetStencil( dataToStencil->GetOutput() );
  stencil->ReverseStencilOff();
  stencil->SetBackgroundValue( 0 );
  
  stencil->Update();
  //  En este punto ya se tienen las imágenes de la segmentación generada por la malla, aunque tienen el tamaño de la región recortada de las imágenes originales
  double originStencil[3];
//   int extension[6];
  
  stencil->GetOutput()->GetOrigin( originStencil );
//   stencil->GetOutput()->GetExtent( extension );
  
//   cout << "origen stencil: " << originStencil[0] << ", " << originStencil[1] << ", " << originStencil[2] << endl;
  
  vtkTransform *transformImages = vtkTransform::New();
  vtkImageReslice *regionExtent = vtkImageReslice::New();
  
  transformImages->Translate( originStencil[0] - spacing[0]*(regionBB[0])  , originStencil[1] - spacing[1]*(regionBB[2]), originStencil[2] - spacing[2]*(regionBB[4]) );
//   
  regionExtent->SetInputConnection( stencil->GetOutputPort() );
  regionExtent->SetResliceTransform( transformImages );
  regionExtent->SetOutputExtent( 0, sizeOrig[0]-1, 0, sizeOrig[1]-1, 0, sizeOrig[2]-1 );
  regionExtent->SetOutputSpacing(spacing[0], spacing[1], spacing[2]);
  regionExtent->SetOutputOrigin( 0, 0, 0 );
  
//   cout << "Escritura de las imágenes de las segmentaciones... " << endl;
  
  char MHDFileName[ strlen(outSegFileName) + 5 ];
  char RAWFileName[ strlen(outSegFileName) + 5 ];
  char extMHD[] = ".mhd";
  char extRAW[] = ".raw";
  
  strcpy( MHDFileName, outSegFileName );
  strcpy( RAWFileName, outSegFileName );
  
  strcat( MHDFileName, extMHD );
  strcat( RAWFileName, extRAW );
  
  vtkMetaImageWriter *writerImages = vtkMetaImageWriter::New();
  writerImages->SetInputConnection( regionExtent->GetOutputPort() );
  writerImages->SetCompression(false);
  writerImages->SetFileName( MHDFileName );
  writerImages->SetRAWFileName( RAWFileName );
    
  try{
  writerImages->Write();
}
  catch( exception& ){
//   cerr << "Error en la escritura." << endl;
  return 1;
}
//   cout << "Realizada." << endl;
  
  writerImages->Delete();
  regionExtent->Delete();
  transformImages->Delete();
  stencil->Delete();
  dataToStencil->Delete();
  normalsMesh->Delete();
  newSegImages->Delete();
  lNewSegImages->Delete();


  return 0;
}

