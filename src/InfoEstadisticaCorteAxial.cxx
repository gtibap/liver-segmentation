// Programa Principal para la segmentación del hígado. 
//En este programa se extrae la media y la desviación estandar de la distribución de intensidad de un corte axial en las imágenes de TAC, segmentado manualmente. Para esto es necesario abrir las imágenes de intensidad y las imágenes de las segmentaciones manuales, y es necesario seleccionar el corte axial, el cual, deberá ser, el que aproximadamente comprenda la mayor área.

// Autor: Gerardo Tibamoso Pedraza. 

// librerias ITK
#include "itkConnectedComponentImageFilter.h"
#include "itkExtractImageFilter.h"
#include "itkImage.h"
#include "itkLabelStatisticsImageFilter.h"
#include "itkRelabelComponentImageFilter.h"
#include "itkRegionOfInterestImageFilter.h"


// librerias básicas de C++
#include <iostream>
#include <stdio.h>
using namespace std;

typedef itk::Image< short, 3 >   ShortImagesType;
typedef itk::Image< unsigned short, 3 >   UShortImagesType;
typedef itk::Image< unsigned char, 3 >   UCharImagesType;
typedef itk::Image< double, 3 >   DoubleImagesType;

// inicio del función principal
void InfoIntensidadCorteAxial( ShortImagesType::Pointer readerOrgImages, UCharImagesType::Pointer readerSegImages, int sliceNumber, double &mean2DRegion, double &sigma2DRegion, double &median2DRegion, int labelNumber )
{
  
  typedef char            CharPixelType;
  typedef unsigned char   UCharPixelType;
  typedef short           ShortPixelType;
  typedef double          DoublePixelType;

  const int Dimension = 3;

  typedef itk::Image< CharPixelType, Dimension >    CharImagesType;
  typedef itk::Image< ShortPixelType, Dimension >   ShortImagesType;
  typedef itk::Image< UCharPixelType, Dimension >   UCharImagesType;
  
  typedef itk::Image< ShortPixelType, 2 > Short2DImagesType;
  typedef itk::Image< UCharPixelType, 2 > UChar2DImagesType;
  
  
  typedef itk::ExtractImageFilter< ShortImagesType, Short2DImagesType > ExtractSliceShortType;
  typedef itk::ExtractImageFilter< UCharImagesType, UChar2DImagesType > ExtractSliceUCharType;
  
  typedef itk::ConnectedComponentImageFilter< UChar2DImagesType, UChar2DImagesType > ConnectedComponent2DFilterType;
  typedef itk::RelabelComponentImageFilter< UChar2DImagesType, UChar2DImagesType > RelabelComponent2DFilterType;
  typedef itk::LabelStatisticsImageFilter< Short2DImagesType, UChar2DImagesType > LabelStatistics2DImagesType;
  
// ----------------------------------------------------------------------  
// Extracción del corte axial seleccionado tanto de las imágenes originales como de las imágenes de las segmentaciones manuales
  ExtractSliceShortType::Pointer  extractOrgFilter = ExtractSliceShortType::New();
  ExtractSliceUCharType::Pointer  extractSegFilter = ExtractSliceUCharType::New();
  
  ShortImagesType::RegionType inputOrgRegion;
  UCharImagesType::RegionType inputSegRegion;
  
  inputOrgRegion = readerOrgImages->GetLargestPossibleRegion();
  inputSegRegion = readerSegImages->GetLargestPossibleRegion();
  
  ShortImagesType::SizeType sizeOrg = inputOrgRegion.GetSize();
  UCharImagesType::SizeType sizeSeg = inputSegRegion.GetSize();
  
  sizeOrg[2] = 0;
  sizeSeg[2] = 0;
  cout << "sizeOrg: " <<   sizeOrg[0] << ", " <<  sizeOrg[1] << ", " <<   sizeOrg[2] << endl;
  cout << "startSeg: " <<  sizeSeg[0] << ", " <<   sizeSeg[1] << ", " <<   sizeSeg[2] << endl;
  
  ShortImagesType::IndexType startOrg = inputOrgRegion.GetIndex();
  UCharImagesType::IndexType startSeg = inputSegRegion.GetIndex();
  
  startOrg[2] = sliceNumber;
  startSeg[2] = sliceNumber;
  cout << "startOrg: " <<   startOrg[0] << ", " <<   startOrg[1] << ", " <<   startOrg[2] << endl;
  cout << "startSeg: " <<   startSeg[0] << ", " <<   startSeg[1] << ", " <<   startSeg[2] << endl;
  cout << "sliceNumber: " << sliceNumber << endl;
  
  ShortImagesType::RegionType desiredOrgRegion;
  UCharImagesType::RegionType desiredSegRegion;
  
  desiredOrgRegion.SetSize( sizeOrg );
  desiredSegRegion.SetSize( sizeSeg );
  
  desiredOrgRegion.SetIndex( startOrg );
  desiredSegRegion.SetIndex( startSeg );
  
  extractOrgFilter->SetExtractionRegion( desiredOrgRegion );
  extractSegFilter->SetExtractionRegion( desiredSegRegion );
  
  extractOrgFilter->SetInput( readerOrgImages );
  extractSegFilter->SetInput( readerSegImages );
  
  //-----------------------
  
  cout << "Estimación de la media y la desviación estandar del corte seleccionado..." << endl;

//   ConnectedComponent2DFilterType::Pointer connectedComponent2DFilter = ConnectedComponent2DFilterType::New();
//   RelabelComponent2DFilterType::Pointer relabelComponent2DFilter = RelabelComponent2DFilterType::New();
  LabelStatistics2DImagesType::Pointer statistics2DFilter = LabelStatistics2DImagesType::New();
  
//   connectedComponent2DFilter->SetInput( extractSegFilter->GetOutput() );
//   relabelComponent2DFilter->SetInput( connectedComponent2DFilter->GetOutput() );
  
  statistics2DFilter->SetInput( extractOrgFilter->GetOutput() );
//   statistics2DFilter->SetLabelInput( relabelComponent2DFilter->GetOutput() );
  statistics2DFilter->SetLabelInput( extractSegFilter->GetOutput() );
  
  statistics2DFilter->Update();
  
  mean2DRegion = statistics2DFilter->GetMean( labelNumber );
  sigma2DRegion = statistics2DFilter->GetSigma( labelNumber );
  
  
  double max2DRegion, min2DRegion;
  
  max2DRegion = statistics2DFilter->GetMaximum( labelNumber );
  min2DRegion = statistics2DFilter->GetMinimum( labelNumber );
  
//   cout << "Máx: " << max2DRegion << "  Mín: " << min2DRegion << endl;
  
//   LabelStatistics2DImagesType::HistogramPointer  histogramImages;
//   LabelStatistics2DImagesType::HistogramType  tHistogram;
  
  int numBins = int(max2DRegion - min2DRegion + 1);
//   int nbin;
  
//   cout << "numBins: " << numBins << endl;
  
  statistics2DFilter->UseHistogramsOn();
  statistics2DFilter->SetHistogramParameters ( numBins, min2DRegion, max2DRegion );
  statistics2DFilter->Update();
  
  median2DRegion = statistics2DFilter->GetMedian( labelNumber );
//   cout << "median= " << median2DRegion << endl;
  
//   histogramImages = statistics2DFilter->GetHistogram( 1 );
//   LabelStatistics2DImagesType::HistogramType::SizeType size;
  
//   size = histogramImages->GetSize();
//   cout << "size: " << size[0] << endl;
  
  /**------------------------------------
  escritura de los datos del histograma en un archivo de texto  */
  
  FILE *fpt;
  
  fpt = fopen( "histogramaImagenes.dat", "a" );

  fprintf(fpt, "media, sigma, mediana, slice:%d\n", sliceNumber);
  fprintf( fpt, "%f, %f, %f\n", mean2DRegion, sigma2DRegion, median2DRegion );
//   fprintf( fpt, "sigma = %f\n", sigma2DRegion );
//   fprintf( fpt, "mediana = %f\n", median2DRegion );
//   fprintf( fpt, "rango = [%d, %d]\n", int(min2DRegion), int( max2DRegion) );
//   for( nbin = 0 ; nbin < int(size[0]) ; nbin++ )
//   {
//     fprintf( fpt, "%d, ", int( histogramImages->GetFrequency( nbin, 0 ) ) );
//   }
//   
  fprintf( fpt, "\n" );
  fclose(fpt);

  /** realizada la escritura de los datos obtenidos en un archivo de texto */
  
  return;

}
