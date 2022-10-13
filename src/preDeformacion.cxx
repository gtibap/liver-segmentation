/**-----------------------
En este programa se preparan las imágenes y la superficie para su deformación
--------------------------**/
  
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
#include "vtkCellArray.h"
#include "vtkImageData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataWriter.h"
#include "vtkSTLReader.h"
#include "vtkSTLWriter.h"

// librerias vtkINRIA3D
#include "itkImageToVTKImageFilter.h"


// librerias básicas de C++
#include <string.h>
#include <iostream>
#include <stdio.h>
using namespace std;

typedef itk::Image< char, 3 >   CharImagesType;
typedef itk::Image< unsigned char, 3 >   UCharImagesType;
typedef itk::Image< short, 3 >   ShortImagesType;
typedef itk::Image< unsigned short, 3 >   UShortImagesType;
typedef itk::Image< float, 3 >  FloatImagesType;
typedef itk::Image< double, 3 >  DoubleImagesType;

int ajusteSuperficie( char * meshLiverFileName, char * meshRegionFileName, double * translate1, double * matrixTransform, double * translate2 );

void InfoIntensidadCorteAxial( ShortImagesType::Pointer readerOrgImages, UCharImagesType::Pointer readerSegImages, int sliceNumber, double &mean2DRegion, double &sigma2DRegion, double &median2DRegion, int labelNumber );

// int SeleccionRegionHigado( ShortImagesType::Pointer  , double *infoIntensity, char apertureFileName[], char meshRegionFileName[], vtkPolyData *meshLiver, double translate1[3], double translate2[3], double matrixTransform[16]  );

// int initialRemesh( vtkPolyData *liverMesh, float minimumDistance );

// void MeshVisualization( vtkPolyData *mesh );


int InicializacionDeformacion( vtkImageData *medianImages, vtkImageData *magGradImages, vtkPolyData *meshLiver, double translate1[3], double translate2[3], double matrixTransform[16], int *dims, int *bb, double *infoIntensity, double *infoMagGradient, int iter, char outSegFileName[], double cs, double cp, double cb, int selFP );
/**---------------------------------
Función principal  */ 
int main( int argc , char **argv )
{
  //---------------------------------
  cout << "Función principal: AjusteParametros.cxx" << endl;
  //---------------------------------

  
  if(argc < 9){
    cerr << "used: " << argv[0] << " DirectorySourceImagesFiles/  numberFile  iter name-coef_Val cs_value cp_value cb_value selFP " << endl;
//     InputBinaryRegionCTImagesFile  InputLiverModelPolyDataFile  OutputImagesFile
    return EXIT_FAILURE;
  }
  
  int iter = atoi(argv[3]);
  int selFP = atoi(argv[8]);
  double cs, cp, cb;
  
  cs = atof(argv[5]);
  cp = atof(argv[6]);
  cb = atof(argv[7]);
  
  cout << "cs, cp, cb: " << cs << cp << cb << endl;
  cout << "cuadrada (1), triangular (2), trapezoidal (3). función: " << selFP << endl;
  
  char origFileName[ strlen(argv[1]) + 20 ];
  char segFileName[ strlen(argv[1]) + 20 ];
  char medianFileName[ strlen(argv[1]) + 40 ];
  char magGradFileName[ strlen(argv[1]) + 40 ];
  char apertureFileName[ strlen(argv[1]) + 40 ];
  char meshRegionFileName[ strlen(argv[1]) + 40 ];
  char meshLiverFileName[ strlen(argv[1]) + 40 ];
//   char infoFileName[strlen(argv[1]) + 50 ];
  char outSegFileName[ strlen(argv[1]) + 70 ];
  int *dimImages, *bb;
  double *infoIntensity, *infoMagGradient;
  int sliceNumber;


  char
			seg00[] = "pat0-label-ucRAI.nii.gz",

      seg01[] = "liver-segSlice001.mhd", 
      seg02[] = "liver-segSlice002.mhd", 
      seg03[] = "liver-segSlice003.mhd", 
      seg04[] = "liver-segSlice004.mhd",
      seg05[] = "liver-segSlice005.mhd",
      seg06[] = "liver-segSlice006.mhd",
      seg07[] = "liver-segSlice007.mhd",
      seg08[] = "liver-segSlice008.mhd",
      seg09[] = "liver-segSlice009.mhd",
      seg10[] = "liver-segSlice010.mhd",
      seg11[] = "liver-segSlice011.mhd",
      seg12[] = "liver-segSlice012.mhd",
      seg13[] = "liver-segSlice013.mhd",
      seg14[] = "liver-segSlice014.mhd",
      seg15[] = "liver-segSlice015.mhd",
      seg16[] = "liver-segSlice016.mhd",
      seg17[] = "liver-segSlice017.mhd",
      seg18[] = "liver-segSlice018.mhd",
      seg19[] = "liver-segSlice019.mhd",
      seg20[] = "liver-segSlice020.mhd";

  char
      median00[] = "imagesRegions/imagesMedian00.mhd", 	
		
      median01[] = "imagesRegions/imagesMedian01.mhd", 
      median02[] = "imagesRegions/imagesMedian02.mhd", 
      median03[] = "imagesRegions/imagesMedian03.mhd", 
      median04[] = "imagesRegions/imagesMedian04.mhd",
      median05[] = "imagesRegions/imagesMedian05.mhd",
      median06[] = "imagesRegions/imagesMedian06.mhd",
      median07[] = "imagesRegions/imagesMedian07.mhd",
      median08[] = "imagesRegions/imagesMedian08.mhd",
      median09[] = "imagesRegions/imagesMedian09.mhd",
      median10[] = "imagesRegions/imagesMedian10.mhd",
      median11[] = "imagesRegions/imagesMedian11.mhd",
      median12[] = "imagesRegions/imagesMedian12.mhd",
      median13[] = "imagesRegions/imagesMedian13.mhd",
      median14[] = "imagesRegions/imagesMedian14.mhd",
      median15[] = "imagesRegions/imagesMedian15.mhd",
      median16[] = "imagesRegions/imagesMedian16.mhd",
      median17[] = "imagesRegions/imagesMedian17.mhd",
      median18[] = "imagesRegions/imagesMedian18.mhd",
      median19[] = "imagesRegions/imagesMedian19.mhd",
      median20[] = "imagesRegions/imagesMedian20.mhd";
      
  char
      magGrad00[] = "imagesRegions/imagesMagGrad00.mhd", 

      magGrad01[] = "imagesRegions/imagesMagGrad01.mhd", 
      magGrad02[] = "imagesRegions/imagesMagGrad02.mhd", 
      magGrad03[] = "imagesRegions/imagesMagGrad03.mhd", 
      magGrad04[] = "imagesRegions/imagesMagGrad04.mhd",
      magGrad05[] = "imagesRegions/imagesMagGrad05.mhd",
      magGrad06[] = "imagesRegions/imagesMagGrad06.mhd",
      magGrad07[] = "imagesRegions/imagesMagGrad07.mhd",
      magGrad08[] = "imagesRegions/imagesMagGrad08.mhd",
      magGrad09[] = "imagesRegions/imagesMagGrad09.mhd",
      magGrad10[] = "imagesRegions/imagesMagGrad10.mhd",
      magGrad11[] = "imagesRegions/imagesMagGrad11.mhd",
      magGrad12[] = "imagesRegions/imagesMagGrad12.mhd",
      magGrad13[] = "imagesRegions/imagesMagGrad13.mhd",
      magGrad14[] = "imagesRegions/imagesMagGrad14.mhd",
      magGrad15[] = "imagesRegions/imagesMagGrad15.mhd",
      magGrad16[] = "imagesRegions/imagesMagGrad16.mhd",
      magGrad17[] = "imagesRegions/imagesMagGrad17.mhd",
      magGrad18[] = "imagesRegions/imagesMagGrad18.mhd",
      magGrad19[] = "imagesRegions/imagesMagGrad19.mhd",
      magGrad20[] = "imagesRegions/imagesMagGrad20.mhd";
      
  char
			meshHeart[] = "meshRegions/meshHeart2mm.stl",
      meshLiver[] = "meshRegions/meshLiver2mm.stl";

  char
      meshRegion00[] = "meshRegions/meshRegion00.stl",

      meshRegion01[] = "meshRegions/meshRegion01.stl",
      meshRegion02[] = "meshRegions/meshRegion02.stl",
      meshRegion03[] = "meshRegions/meshRegion03.stl",
      meshRegion04[] = "meshRegions/meshRegion04.stl",
      meshRegion05[] = "meshRegions/meshRegion05.stl",
      meshRegion06[] = "meshRegions/meshRegion06.stl",
      meshRegion07[] = "meshRegions/meshRegion07.stl",
      meshRegion08[] = "meshRegions/meshRegion08.stl",
      meshRegion09[] = "meshRegions/meshRegion09.stl",
      meshRegion10[] = "meshRegions/meshRegion10.stl",
      meshRegion11[] = "meshRegions/meshRegion11.stl",
      meshRegion12[] = "meshRegions/meshRegion12.stl",
      meshRegion13[] = "meshRegions/meshRegion13.stl",
      meshRegion14[] = "meshRegions/meshRegion14.stl",
      meshRegion15[] = "meshRegions/meshRegion15.stl",
      meshRegion16[] = "meshRegions/meshRegion16.stl",
      meshRegion17[] = "meshRegions/meshRegion17.stl",
      meshRegion18[] = "meshRegions/meshRegion18.stl",
      meshRegion19[] = "meshRegions/meshRegion19.stl",
      meshRegion20[] = "meshRegions/meshRegion20.stl";

  char
      outSeg00[] = "imagesRegions/EvalCoef/heart-outSeg01",

      outSeg01[] = "imagesRegions/EvalCoef/liver-outSeg01",
      outSeg02[] = "imagesRegions/EvalCoef/liver-outSeg02",
      outSeg03[] = "imagesRegions/EvalCoef/liver-outSeg03",
      outSeg04[] = "imagesRegions/EvalCoef/liver-outSeg04",
      outSeg05[] = "imagesRegions/EvalCoef/liver-outSeg05",
      outSeg06[] = "imagesRegions/EvalCoef/liver-outSeg06",
      outSeg07[] = "imagesRegions/EvalCoef/liver-outSeg07",
      outSeg08[] = "imagesRegions/EvalCoef/liver-outSeg08",
      outSeg09[] = "imagesRegions/EvalCoef/liver-outSeg09",
      outSeg10[] = "imagesRegions/EvalCoef/liver-outSeg10",
      outSeg11[] = "imagesRegions/EvalCoef/liver-outSeg11",
      outSeg12[] = "imagesRegions/EvalCoef/liver-outSeg12",
      outSeg13[] = "imagesRegions/EvalCoef/liver-outSeg13",
      outSeg14[] = "imagesRegions/EvalCoef/liver-outSeg14",
      outSeg15[] = "imagesRegions/EvalCoef/liver-outSeg15",
      outSeg16[] = "imagesRegions/EvalCoef/liver-outSeg16",
      outSeg17[] = "imagesRegions/EvalCoef/liver-outSeg17",
      outSeg18[] = "imagesRegions/EvalCoef/liver-outSeg18",
      outSeg19[] = "imagesRegions/EvalCoef/liver-outSeg19",
      outSeg20[] = "imagesRegions/EvalCoef/liver-outSeg20";

  int
      dims00[3] = {384, 384, 150},

      dims01[3] = {512, 512, 183},
      dims02[3] = {512, 512, 64},
      dims03[3] = {512, 512, 79},
      dims04[3] = {512, 512, 212},
      dims05[3] = {512, 512, 319},
      dims06[3] = {512, 512, 111},
      dims07[3] = {512, 512, 251},
      dims08[3] = {512, 512, 228},
      dims09[3] = {512, 512, 210},
      dims10[3] = {512, 512, 191},
      dims11[3] = {512, 512, 388},
      dims12[3] = {512, 512, 220},
      dims13[3] = {512, 512, 145},
      dims14[3] = {512, 512, 129},
      dims15[3] = {512, 512, 394},
      dims16[3] = {512, 512, 151},
      dims17[3] = {512, 512, 121},
      dims18[3] = {512, 512, 245},
      dims19[3] = {512, 512, 335},
      dims20[3] = {512, 512, 183};

  int 
      bb00[] = {120, 268, 30, 262, 0, 149},

      bb01[] = {35, 359, 55, 376, 50, 178},
      bb02[] = {14, 417, 55, 383, 0, 63},
      bb03[] = {37, 352, 59, 354, 7, 78},
      bb04[] = {41, 362, 92, 413, 41, 211},
      bb05[] = {0, 280, 73, 461, 95, 299},
      bb06[] = {39, 418, 98, 350, 10, 103},
      bb07[] = {36, 347, 96, 369, 40, 245},
      bb08[] = {17, 371, 72, 412, 10, 227},
      bb09[] = {19, 347, 117, 393, 4, 209},
      bb10[] = {29, 346, 57, 393, 0, 190},
      bb11[] = {20, 440, 111, 422, 187, 387},
      bb12[] = {36, 352, 100, 408, 27, 219},
      bb13[] = {37, 338, 84, 381, 12, 144},
      bb14[] = {28, 279, 91, 403, 41, 85},
      bb15[] = {17, 366, 19, 342, 229, 393},
      bb16[] = {28, 378, 88, 418, 0, 150},
      bb17[] = {37, 344, 86, 359, 2, 114},
      bb18[] = {41, 340, 120, 392, 17, 244},
      bb19[] = {7, 437, 57, 392, 64, 327},
      bb20[] = {15, 333, 116, 399, 101, 182};
      
  int
      sliceNumber00 = 57,

      sliceNumber01 = 137,
      sliceNumber02 = 43,
      sliceNumber03 = 54,
      sliceNumber04 = 145,
      sliceNumber05 = 216,
      sliceNumber06 = 68,
      sliceNumber07 = 167,
      sliceNumber08 = 143,
      sliceNumber09 = 131,
      sliceNumber10 = 122,
      sliceNumber11 = 329,
      sliceNumber12 = 161,
      sliceNumber13 = 110,
      sliceNumber14 = 68,
      sliceNumber15 = 320,
      sliceNumber16 = 105,
      sliceNumber17 = 73,
      sliceNumber18 = 162,
      sliceNumber19 = 217,
      sliceNumber20 = 150;
      
  double
      infoIntensity00[3] = {1553.17, 216.178, 1550.0},

      infoIntensity01[3] = {88.519852, 15.242842, 88.117241},
      infoIntensity07[3] = {97.917932, 23.643111, 98.030189},
      infoIntensity13[3] ={105.742819, 14.105102, 105.922581},
	
      infoIntensity02[3] ={95.415603, 33.112783, 92.883657},
      infoIntensity03[3] ={111.374480, 33.468011, 104.070000},
      infoIntensity04[3] ={188.198324, 48.604293, 183.941385},
      infoIntensity05[3] ={143.900469, 32.982319, 137.142349},
      infoIntensity06[3] ={146.291433, 16.580405, 143.975877},
    
      infoIntensity08[3] ={110.540526, 40.092671, 107.081871},
      infoIntensity09[3] ={175.252312, 31.236308, 174.064171},
      infoIntensity10[3] ={181.657906, 39.764360, 179.063529},
      infoIntensity11[3] ={150.423523, 26.727512, 147.975439},
      infoIntensity12[3] ={135.363025, 31.074490, 134.946250},
    
      infoIntensity14[3] ={90.575735, 15.539066, 91.792627},
      infoIntensity15[3] ={153.887665, 41.606842, 152.931947},
      infoIntensity16[3] ={107.866434, 19.464445, 106.947566},
      infoIntensity17[3] ={135.424498, 18.422857, 133.883162},
      infoIntensity18[3] ={141.031342, 37.706809, 137.974444},
      infoIntensity19[3] ={140.602963, 32.273074, 138.968153},
      infoIntensity20[3] ={122.656800, 13.347638, 122.865285};	

  double
      infoMagGradient00[3] = {407.604089, 243.489172, 371.446587},

      infoMagGradient01[3] = {407.604089, 243.489172, 371.446587},
      infoMagGradient07[3] = {634.938347, 341.964958, 579.355674},
      infoMagGradient13[3] = {377.876979, 211.933321, 337.377471},

	
      infoMagGradient02[3] = {645.290999, 740.230871, 421.454214},
      infoMagGradient03[3] = {579.721773, 606.171152, 374.443791},
      infoMagGradient04[3] = {1009.249678, 676.588576, 848.377547},
      infoMagGradient05[3] = {679.293500, 453.608857, 564.384656},
      infoMagGradient06[3] = {357.059133, 335.880933, 248.427237},
      
      infoMagGradient08[3] = {769.831243, 507.368015, 659.410573},
      infoMagGradient09[3] = {738.407742, 476.217676, 651.419952},
      infoMagGradient10[3] = {715.690829, 542.692702, 583.417937},
      infoMagGradient11[3] = {446.599508, 352.673281, 350.414211},
      infoMagGradient12[3] = {790.698610, 451.763074, 704.346163},
      
      infoMagGradient14[3] = {355.291101, 383.357734, 245.439012},
      infoMagGradient15[3] = {781.049510, 704.287696, 587.409791},
      infoMagGradient16[3] = {571.262821, 464.537324, 491.457700},
      infoMagGradient17[3] = {433.581018, 414.518326, 324.460183},
      infoMagGradient18[3] = {775.936164, 522.716179, 659.395073},
      infoMagGradient19[3] = {661.545690, 469.863723, 554.389872},
      infoMagGradient20[3] = {425.964205, 498.296433, 313.468616};
	
  strcpy( origFileName, argv[1] );
  strcpy( segFileName,  argv[1] );
  strcpy( medianFileName,  argv[1] );
  strcpy( magGradFileName,  argv[1] );
  strcpy( apertureFileName,  argv[1] );
  strcpy( meshRegionFileName,  argv[1] );
  strcpy( outSegFileName,  argv[1] );
  
  strcpy( meshLiverFileName,  argv[1] );
//  strcat( meshLiverFileName,  meshLiver );
  strcat( meshLiverFileName,  meshHeart );
  
  switch( atoi(argv[2]) ){

    case 0: 
      strcat( segFileName,  seg00 );
      strcat( medianFileName,  median00 );
      strcat(meshRegionFileName, meshRegion00 );
      strcat( magGradFileName,  magGrad00 );
      strcat( outSegFileName,  outSeg00 );
      strcat( outSegFileName,  argv[4] );
      bb = bb00;
      sliceNumber = sliceNumber00;
      infoIntensity = infoIntensity00;
      infoMagGradient = infoMagGradient00;
      dimImages = dims00;
      break;

    case 1: 
      strcat( segFileName,  seg01 );
      strcat( medianFileName,  median01 );
      strcat(meshRegionFileName, meshRegion01 );
      strcat( magGradFileName,  magGrad01 );
      strcat( outSegFileName,  outSeg01 );
      strcat( outSegFileName,  argv[4] );
      bb = bb01;
      sliceNumber = sliceNumber01;
      infoIntensity = infoIntensity01;
      infoMagGradient = infoMagGradient01;
      dimImages = dims01;
      break;

    case 2: 
      strcat( segFileName,  seg02 );
      strcat( medianFileName,  median02 );
      strcat(meshRegionFileName, meshRegion02 );
      strcat( magGradFileName,  magGrad02 );
      strcat( outSegFileName,  outSeg02 );
      strcat( outSegFileName,  argv[4] );
      bb = bb02;
      sliceNumber = sliceNumber02;
      infoIntensity = infoIntensity02;
      infoMagGradient = infoMagGradient02;
      dimImages = dims02;
      break;
      
    case 3: 
      strcat( segFileName,  seg03 );
      strcat( medianFileName,  median03 );
      strcat(meshRegionFileName, meshRegion03 );
      strcat( magGradFileName,  magGrad03 );
      strcat( outSegFileName,  outSeg03 );
      strcat( outSegFileName,  argv[4] );
      bb = bb03;
      sliceNumber = sliceNumber03;
      infoIntensity = infoIntensity03;
      infoMagGradient = infoMagGradient03;
      dimImages = dims03;
      break;
      
    case 4: 
      strcat( segFileName,  seg04 );
      strcat( medianFileName,  median04 );
      strcat(meshRegionFileName, meshRegion04 );
      strcat( magGradFileName,  magGrad04 );
      strcat( outSegFileName,  outSeg04 );
      strcat( outSegFileName,  argv[4] );
      bb = bb04;
      sliceNumber = sliceNumber04;
      infoIntensity = infoIntensity04;
      infoMagGradient = infoMagGradient04;
      dimImages = dims04;
      break;      
      
    case 5: 
      strcat( segFileName,  seg05 );
      strcat( medianFileName,  median05 );
      strcat(meshRegionFileName, meshRegion05 );
      strcat( magGradFileName,  magGrad05 );
      strcat( outSegFileName,  outSeg05 );
      strcat( outSegFileName,  argv[4] );
      bb = bb05;
      sliceNumber = sliceNumber05;
      infoIntensity = infoIntensity05;
      infoMagGradient = infoMagGradient05;
      dimImages = dims05;
      break;
      
    case 6: 
      strcat( segFileName,  seg06 );
      strcat( medianFileName,  median06 );
      strcat(meshRegionFileName, meshRegion06 );
      strcat( magGradFileName,  magGrad06 );
      strcat( outSegFileName,  outSeg06 );
      strcat( outSegFileName,  argv[4] );
      bb = bb06;
      sliceNumber = sliceNumber06;
      infoIntensity = infoIntensity06;
      infoMagGradient = infoMagGradient06;
      dimImages = dims06;
      break;
      
    case 7: 
      strcat( segFileName,  seg07 );
      strcat( medianFileName,  median07 );
      strcat( magGradFileName,  magGrad07 );
      strcat(meshRegionFileName, meshRegion07 );
      strcat( outSegFileName,  outSeg07 );
      strcat( outSegFileName,  argv[4] );
      bb = bb07;
      sliceNumber = sliceNumber07;
      infoIntensity = infoIntensity07;
      infoMagGradient = infoMagGradient07;
      dimImages = dims07;
      break;

    case 8: 
      strcat( segFileName,  seg08 );
      strcat( medianFileName,  median08 );
      strcat( magGradFileName,  magGrad08 );
      strcat(meshRegionFileName, meshRegion08 );
      strcat( outSegFileName,  outSeg08 );
      strcat( outSegFileName,  argv[4] );
      bb = bb08;
      sliceNumber = sliceNumber08;
      infoIntensity = infoIntensity08;
      infoMagGradient = infoMagGradient08;
      dimImages = dims08;
      break;
      
    case 9: 
      strcat( segFileName,  seg09 );
      strcat( medianFileName,  median09 );
      strcat(meshRegionFileName, meshRegion09 );
      strcat( magGradFileName,  magGrad09 );
      strcat( outSegFileName,  outSeg09 );
      strcat( outSegFileName,  argv[4] );
      bb = bb09;
      sliceNumber = sliceNumber09;
      infoIntensity = infoIntensity09;
      infoMagGradient = infoMagGradient09;
      dimImages = dims09;
      break;
      
    case 10: 
      strcat( segFileName,  seg10 );
      strcat( medianFileName,  median10 );
      strcat(meshRegionFileName, meshRegion10 );
      strcat( magGradFileName,  magGrad10 );
      strcat( outSegFileName,  outSeg10 );
      strcat( outSegFileName,  argv[4] );
      bb = bb10;
      sliceNumber = sliceNumber10;
      infoIntensity = infoIntensity10;
      infoMagGradient = infoMagGradient10;
      dimImages = dims10;
      break;
      
    case 11: 
      strcat( segFileName,  seg11 );
      strcat( medianFileName,  median11 );
      strcat(meshRegionFileName, meshRegion11 );
      strcat( magGradFileName,  magGrad11 );
      strcat( outSegFileName,  outSeg11 );
      strcat( outSegFileName,  argv[4] );
      bb = bb11;
      sliceNumber = sliceNumber11;
      infoIntensity = infoIntensity11;
      infoMagGradient = infoMagGradient11;
      dimImages = dims11;
      break;
      
    case 12: 
      strcat( segFileName,  seg12 );
      strcat( medianFileName,  median12 );
      strcat(meshRegionFileName, meshRegion12 );
      strcat( magGradFileName,  magGrad12 );
      strcat( outSegFileName,  outSeg12 );
      strcat( outSegFileName,  argv[4] );
      bb = bb12;
      sliceNumber = sliceNumber12;
      infoIntensity = infoIntensity12;
      infoMagGradient = infoMagGradient12;
      dimImages = dims12;
      break;
      
    case 13: 
      strcat( segFileName,  seg13 );
      strcat( medianFileName,  median13 );
      strcat(meshRegionFileName, meshRegion13 );
      strcat( magGradFileName,  magGrad13 );
      strcat( outSegFileName,  outSeg13 );
      strcat( outSegFileName,  argv[4] );
      bb = bb13;
      sliceNumber = sliceNumber13;
      infoIntensity = infoIntensity13;
      infoMagGradient = infoMagGradient13;
      dimImages = dims13;
      break;

    case 14: 
      strcat( segFileName,  seg14 );
      strcat( medianFileName,  median14 );
      strcat(meshRegionFileName, meshRegion14 );
      strcat( magGradFileName,  magGrad14 );
      strcat( outSegFileName,  outSeg14 );
      strcat( outSegFileName,  argv[4] );
      bb = bb14;
      sliceNumber = sliceNumber14;
      infoIntensity = infoIntensity14;
      infoMagGradient = infoMagGradient14;
      dimImages = dims14;
      break;
      
    case 15: 
      strcat( segFileName,  seg15 );
      strcat( medianFileName,  median15 );
      strcat(meshRegionFileName, meshRegion15 );
      strcat( magGradFileName,  magGrad15 );
      strcat( outSegFileName,  outSeg15 );
      strcat( outSegFileName,  argv[4] );
      bb = bb15;
      sliceNumber = sliceNumber15;
      infoIntensity = infoIntensity15;
      infoMagGradient = infoMagGradient15;
      dimImages = dims15;
      break;
      
    case 16: 
      strcat( segFileName,  seg16 );
      strcat( medianFileName,  median16 );
      strcat(meshRegionFileName, meshRegion16 );
      strcat( magGradFileName,  magGrad16 );
      strcat( outSegFileName,  outSeg16 );
      strcat( outSegFileName,  argv[4] );
      bb = bb16;
      sliceNumber = sliceNumber16;
      infoIntensity = infoIntensity16;
      infoMagGradient = infoMagGradient16;
      dimImages = dims16;
      break;
      
    case 17: 
      strcat( segFileName,  seg17 );
      strcat( medianFileName,  median17 );
      strcat(meshRegionFileName, meshRegion17 );
      strcat( magGradFileName,  magGrad17 );
      strcat( outSegFileName,  outSeg17 );
      strcat( outSegFileName,  argv[4] );
      bb = bb17;
      sliceNumber = sliceNumber17;
      infoIntensity = infoIntensity17;
      infoMagGradient = infoMagGradient17;
      dimImages = dims17;
      break;
      
    case 18: 
      strcat( segFileName,  seg18 );
      strcat( medianFileName,  median18 );
      strcat(meshRegionFileName, meshRegion18 );
      strcat( magGradFileName,  magGrad18 );
      strcat( outSegFileName,  outSeg18 );
      strcat( outSegFileName,  argv[4] );
      bb = bb18;
      sliceNumber = sliceNumber18;
      infoIntensity = infoIntensity18;
      infoMagGradient = infoMagGradient18;
      dimImages = dims18;
      break;
	
    case 19: 
      strcat( segFileName,  seg19 );
      strcat( medianFileName,  median19 );
      strcat(meshRegionFileName, meshRegion19 );
      strcat( magGradFileName,  magGrad19 );
      strcat( outSegFileName,  outSeg19 );
      strcat( outSegFileName,  argv[4] );
      bb = bb19;
      sliceNumber = sliceNumber19;
      infoIntensity = infoIntensity19;
      infoMagGradient = infoMagGradient19;
      dimImages = dims19;
      break;
      
    case 20: 
      strcat( segFileName,  seg20 );
      strcat( medianFileName,  median20 );
      strcat(meshRegionFileName, meshRegion20 );
      strcat( magGradFileName,  magGrad20 );
      strcat( outSegFileName,  outSeg20 );
      strcat( outSegFileName,  argv[4] );
      bb = bb20;
      sliceNumber = sliceNumber20;
      infoIntensity = infoIntensity20;
      infoMagGradient = infoMagGradient20;
      dimImages = dims20;
      break;
      
    default:
      cout << "Registro no encontrado" << endl;
      return EXIT_FAILURE;
  }

  cout << "segFileName: " << segFileName << endl;
  cout << "outSegFileName: " << outSegFileName << endl;
  
    /**--------------------------------------
  Definición de tipos para la lectura y/o escritura de las imágenes de TAC, de las segmentaciones manuales del hígado, y de la magnitud del gradiente */
  
  typedef itk::ImageFileReader< ShortImagesType > ReaderShortImagesType;
  typedef itk::ImageFileReader< UCharImagesType > ReaderUCharImagesType;
  typedef itk::ImageFileReader< FloatImagesType > ReaderFloatImagesType;
  typedef itk::ImageFileReader< DoubleImagesType > ReaderDoubleImagesType;
  
  typedef itk::ImageFileWriter< ShortImagesType >  WriterShortImagesType;
  typedef itk::ImageFileWriter< UShortImagesType > WriterUShortImagesType;
  
  typedef itk::ImageToVTKImageFilter< ShortImagesType > ConnectorShortFilterType;
  typedef itk::ImageToVTKImageFilter< FloatImagesType > ConnectorFloatFilterType;
  typedef itk::ImageToVTKImageFilter< DoubleImagesType > ConnectorDoubleFilterType;
  
//   typedef itk::RegionOfInterestImageFilter< ShortImagesType, ShortImagesType > ShortRegionFilter;
  typedef itk::RegionOfInterestImageFilter< UCharImagesType, UCharImagesType > UCharRegionFilter;

  cout << "\nLectura de las imágenes de la seg. manual... " << endl;
  
//   ReaderShortImagesType::Pointer readerOrgImages = ReaderShortImagesType::New();
  ReaderUCharImagesType::Pointer readerSegImages = ReaderUCharImagesType::New();
  
//   readerOrgImages->SetFileName( origFileName );
  readerSegImages->SetFileName( segFileName ); 
  //   lectura de las imágenes de entrada
//   try{
//     readerOrgImages->Update();
//   }
//   catch( itk::ExceptionObject & err ){
//     cerr << "Error en la lectura de las imágenes de intensidad." << endl;
//     return EXIT_FAILURE;
//   }
//   
  try{
    readerSegImages->Update();
  }
  catch( itk::ExceptionObject & err ){
    cerr << "Error en la lectura de las imágenes de las segmentaciones." << endl;
    return EXIT_FAILURE;
  }
  
  cout << "realizada." << endl;
  
  /** Lectura del tamaño original de las imágenes de entrada */
  //  dimensiones de las imágenes leidas
  UCharImagesType::RegionType selectRegion;
  UCharImagesType::SizeType sizeRegion;
  int sizeOrig[3];
  
  selectRegion = readerSegImages->GetOutput()->GetLargestPossibleRegion();
  sizeRegion = selectRegion.GetSize();
  
  sizeOrig[0] = sizeRegion.GetElement(0);
  sizeOrig[1] = sizeRegion.GetElement(1);
  sizeOrig[2] = sizeRegion.GetElement(2);

  cout << "Tamaño Orig: " << sizeOrig[0] << ", " << sizeOrig[1] << ", " << sizeOrig[2] << endl;
  
  /**Lectura realizada.
  ------------------------------------*/
  /**----------------------------------
  extracción de la región donde está el hígado, definida por el BoundingBox, tanto para las imágenes de TAC como para las de las segmentaciones manuales*/
  
  cout << "Extracción de la region seleccionada... " << endl;
//   ShortRegionFilter::Pointer orgRegionFilter = ShortRegionFilter::New();
  UCharRegionFilter::Pointer segRegionFilter = UCharRegionFilter::New();
//   
//   ShortImagesType::IndexType startOrgRegion;
  UCharImagesType::IndexType startSegRegion;
//   
//   startOrgRegion[0] = bb[0];
  startSegRegion[0] = bb[0];
//   
//   startOrgRegion[1] = bb[2];
  startSegRegion[1] = bb[2];
//   
//   startOrgRegion[2] = bb[4];
  startSegRegion[2] = bb[4];
//   //----------------------
//   
//   ShortImagesType::SizeType sizeOrgRegion;
  UCharImagesType::SizeType sizeSegRegion;
//   
//   sizeOrgRegion[0] = bb[1] - bb[0] + 1;
  sizeSegRegion[0] = bb[1] - bb[0] + 1;
//   
//   sizeOrgRegion[1] = bb[3] - bb[2] + 1;
  sizeSegRegion[1] = bb[3] - bb[2] + 1;
//   
//   sizeOrgRegion[2] = bb[5] - bb[4] + 1;
  sizeSegRegion[2] = bb[5] - bb[4] + 1;
//   //----------------------
//   ShortImagesType::RegionType desiredOrgRegion;
  UCharImagesType::RegionType desiredSegRegion;
//   
//   desiredOrgRegion.SetSize( sizeOrgRegion );
  desiredSegRegion.SetSize( sizeSegRegion );
//   
//   desiredOrgRegion.SetIndex( startOrgRegion );
  desiredSegRegion.SetIndex( startSegRegion );
//   //-----------------------
//   
//   orgRegionFilter->SetRegionOfInterest( desiredOrgRegion );
  segRegionFilter->SetRegionOfInterest( desiredSegRegion );
//   
//   orgRegionFilter->SetInput( readerOrgImages->GetOutput() );
  segRegionFilter->SetInput( readerSegImages->GetOutput() );
//   //------------------------
//   
//   orgRegionFilter->Update();
  segRegionFilter->Update();
  cout << "Extracción realizada." << endl;
//   
  /**Extracción realizada
  -------------------------*/
  /**-- Filtro de mediana aplicado a las imágenes----*/
  /**--------------------------------------
  Antes de leer las imágenes originales, hay que confirmar que la mediana y el operador sobel de las imágenes no se hayan calculado antes.Si ya lo fueron, entonces se leen estas imágenes y, por supuesto, ya no se recalcularan.   */
//   cout << "Filtro de Mediana... " << endl;
//   
//   typedef itk::MedianImageFilter< ShortImagesType, ShortImagesType > MedianImagesFilterType;  
//   MedianImagesFilterType::Pointer medianImagesFilter = MedianImagesFilterType::New();

  cout << "Lectura de las imágenes de mediana ... " << endl;
  int readerFlag = 1;
  ShortImagesType *medianImages;
  
  ReaderShortImagesType::Pointer readerMedianImages = ReaderShortImagesType::New();

  readerMedianImages->SetFileName( medianFileName );
  try{
    readerMedianImages->Update();
  }
  catch( itk::ExceptionObject & err ){
    cerr << "Error en la lectura de las imágenes de mediana." << endl;
    return EXIT_FAILURE;
  }
  medianImages = readerMedianImages->GetOutput();

  
  /**-- Informacion de media y desviación estandar del hígado en un corte axial segmentado manualmente--*/
    
  cout << "Información estadística de Intensidad del corte seleccionado..." << endl;
  double meanIntensity, sigmaIntensity, medianIntensity;
    
  /** ********* **/
  InfoIntensidadCorteAxial( medianImages, segRegionFilter->GetOutput(), sliceNumber-bb[4], meanIntensity, sigmaIntensity, medianIntensity, 1 );
  /** ********* **/
  
  infoIntensity[0] = meanIntensity;
  infoIntensity[1] = sigmaIntensity;
  infoIntensity[2] = medianIntensity;
    
  cout << "media Intensidad= " << infoIntensity[0] << endl;
  cout << "Desviación Estandar Intensidad= " << infoIntensity[1] << endl;
  cout << "mediana Intensidad= " << infoIntensity[2] << endl;
    
  /**-- Informacion de media y desviación estandar del hígado en un corte axial segmentado manualmente--*/
  
  /**--------------------------
  AJUSTE SUPERFICIE -----------*/
  double translate1[3], matrixTransform[16], translate2[3];
  
  int val = ajusteSuperficie( meshLiverFileName, meshRegionFileName, translate1, matrixTransform, translate2 );
  
  if (val ==0)
    cout<< "ajuste manual realizado."<< endl;
  else
    cout << "problemas con el ajuste manual." << endl;
  
  cout << "Translate1: ";
  for( int i=0; i<3; i++ ){
    cout << translate1[i] << ", ";
  }
  cout << endl;
  
  int cuenta;
  cout << "La matriz de transformacion final es: " << endl;
  cuenta=0;
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      cout << matrixTransform[cuenta] << " ";
      ++cuenta;
    }
    cout << endl;
  }
  cout << "Translate2: ";
  for( int i=0; i<3; i++ ){
    cout << translate2[i] << ", ";
  }
  cout << endl;

  /**-- lectura del modelo promedio del hígado ----*/
  cout << "Lectura de la malla modelo de forma del hígado... " << endl;

  vtkSTLReader *readerMeshLiver = vtkSTLReader::New();
  readerMeshLiver->SetFileName( meshLiverFileName );
  try{
    readerMeshLiver->Update();
  }
  catch( exception& ){
    cerr << "Error en la lectura de la superficie del hígado promedio." << endl;
    return EXIT_FAILURE;
  }
  cout << "Lectura realizada." << endl;
  /**-- lectura del modelo promedio del hígado ----*/
  
  /**-- Visualización de la malla. */
  //     MeshVisualization( meshLiver );
  /**-- Visualización de la malla. */
  
  
  /** Cálculo de la mag. del gradiente  */
  // Ahora se calculará el modulo del gradiente estimado con el operador Sobel. Este filtro requiere que los datos sean de tipo real, por lo tanto debe hacerse un casting (conversión) de entero corto(short) a real(double)
  cout << "Imágenes del Filtrado de Sobel... " << endl;
  
  FloatImagesType *magGradImages;
  readerFlag = 1;

  typedef itk::CastImageFilter< ShortImagesType, FloatImagesType > CastShortToFloatImagesType;
  typedef itk::SobelEdgeDetectionImageFilter< FloatImagesType, FloatImagesType > SobelEdgeDetectionImagesFilterType;  
  typedef itk::CastImageFilter< FloatImagesType, ShortImagesType > CastFloatToShortImagesType;
  
  CastShortToFloatImagesType::Pointer castShortToFloatImagesFilter = CastShortToFloatImagesType::New();
  
//   SobelEdgeDetectionImagesFilterType::Pointer sobelEdgeDetectionImagesFilter = SobelEdgeDetectionImagesFilterType::New();
    
//   CastFloatToShortImagesType::Pointer castRealToShortImagesFilter = CastFloatToShortImagesType::New();
    
  ReaderShortImagesType::Pointer readerSobelImages = ReaderShortImagesType::New();
  readerSobelImages->SetFileName( magGradFileName );
  try{
    readerSobelImages->Update();
  }
  catch( itk::ExceptionObject & err ){
    cerr << "Error en la lectura de las imágenes de magGrad." << endl;
    return EXIT_FAILURE;
  }
  /**-------Casting de MagGrad de Entero a Real ---------------*/
  castShortToFloatImagesFilter->SetInput( readerSobelImages->GetOutput() );
  castShortToFloatImagesFilter->Update();
  /**-------Casting de MagGrad de Entero a Real---------------*/

  /**-------Casting de MagGrad de Real a Entero---------------*/
//   castRealToShortImagesFilter->SetInput( readerSobelImages->GetOutput() );
//   castRealToShortImagesFilter->Update();
  /**-------Casting de MagGrad de Real a Entero---------------*/
  
  magGradImages = castShortToFloatImagesFilter->GetOutput();
  
  /**---Estimación de la información Estadística MagGrad------*/
  cout << "Información estadística de la magnitud del gradiente del corte seleccionado..." << endl;
  double meanSobel, sigmaSobel, medianSobel;
  
  /** ************ **/
  InfoIntensidadCorteAxial(readerSobelImages->GetOutput(), segRegionFilter->GetOutput(), sliceNumber-bb[4], meanSobel, sigmaSobel, medianSobel, 1);
  /** ************ **/

  infoMagGradient[0] = meanSobel;
  infoMagGradient[1] = sigmaSobel;
  infoMagGradient[2] = medianSobel;
      
  cout << "media Mag.Gradiente= " << infoMagGradient[0] << endl;
  cout << "Desviación Estandar Mag.Gradiente= " << infoMagGradient[1] << endl;
  cout << "median Mag.Gradiente= " << infoMagGradient[2] << endl;
  /**---Estimación de la información Estadística MagGrad------*/
  
  
  /**--- Exportación de las imágenes de itk a vtk... */
  cout << "Exportación de las imágenes del filtro de mediana ... " << endl;
  ConnectorShortFilterType::Pointer  medianImagesConnector = ConnectorShortFilterType::New();
  medianImagesConnector->SetInput( medianImages );
  medianImagesConnector->Update();
  cout << "realizada." << endl;

  cout << "Exportación de las imágenes de la mag. del gradiente... " << endl;
  ConnectorFloatFilterType::Pointer  sobelImagesConnector = ConnectorFloatFilterType::New();
  sobelImagesConnector->SetInput( magGradImages );
  sobelImagesConnector->Update();
  cout << "realizada." << endl;
  /**--- Exportación de las imágenes de itk a vtk... */

  
  
  
  
  
  
  
  
  
  
  
  
  
  /**------ lectura datos transformación superficie  -------**/
/*  FILE *lecturaDatos;
  char estudio[10];
  int numReg;
  float vtemp;
  double translate1[3], translate2[3], matrixTransform[16];
  lecturaDatos = fopen("transTemp.dat","r");
    
  fscanf(lecturaDatos, "%s %d",estudio, &numReg);
  cout <<"nombre "<< estudio << numReg << endl;
  cout << "Translate1: ";
  for( int i=0; i<3; i++ ){
    fscanf(lecturaDatos, "%f",&vtemp);
    translate1[i] = double(vtemp);
    cout << translate1[i] << ", ";
  }
  cout << endl;
  
  int cuenta;
  cout << "La matriz de transformacion final es: " << endl;
  cuenta=0;
  for(int i=0; i<4; i++){
    for(int j=0; j<4; j++){
      fscanf(lecturaDatos, "%f", &vtemp);
      matrixTransform[cuenta] = double (vtemp);
      cout << matrixTransform[cuenta] << " ";
      ++cuenta;
    }
    cout << endl;
  }
  
  cout << "Translate2: ";
  for( int i=0; i<3; i++ ){
    fscanf(lecturaDatos, "%f", &vtemp);
    translate2[i] = double (vtemp);
    cout << translate2[i] << ", ";
  }
  cout << endl;
  fclose(lecturaDatos);*/
  /**------ lectura datos transformación superficie  -------**/
  
//   if(atoi(argv[2]) == numReg)
//   {
    int flagRegion;
    cout << "Iniciar deformación ? (1/0): ";
    cin >> flagRegion;
    if( flagRegion == 1  )
    {
      /** ************* **/
      InicializacionDeformacion( medianImagesConnector->GetOutput(), sobelImagesConnector->GetOutput(), readerMeshLiver->GetOutput(), translate1, translate2, matrixTransform, dimImages, bb, infoIntensity, infoMagGradient, iter, outSegFileName, cs, cp, cb, selFP );
      /** ************* **/
    }
//   }
//   else{
//     cout << "Los datos de tranformación no son para el estudio especificado" << endl;
//   }
  /**  Exportación realizada.
      --------------------------------*/

  readerMeshLiver->Delete();
  
  cout << "Salida Normal." << endl;
  return EXIT_SUCCESS;

}
