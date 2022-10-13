// #include "vtkSphereSource.h"
#include "vtkActor.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkCamera.h"
#include "vtkCellArray.h"
#include "vtkConeSource.h"
#include "vtkFloatArray.h"
#include "vtkGlyph3D.h"
#include "vtkIdList.h"
#include "vtkIdTypeArray.h"
#include "vtkImageAccumulate.h"
#include "vtkImageData.h"
#include "vtkImageGradient.h"
#include "vtkImageMagnitude.h"
#include "vtkImageMedian3D.h"
#include "vtkImageShiftScale.h"
#include "vtkImageStencil.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkMaskPoints.h"
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

// #include "vtkPolyDataReader.h"
// #include "vtkProgrammableFilter.h"


#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
// #include "vtkPolyDataWriter.h"
// #include "vtkSTLReader.h"
#include <iostream>
using namespace std;

void times(double s, double A[3], double C[3]);
double vtimes (double A[3], double B[3]);
double norm (double A[3]);
void divide (double A[3], double s, double C[3]);
void vplus (double A[3], double B[3], double C[3]);
void vminus (double A[3], double B[3], double C[3]);

/**-------------------------------
Funciones
----------------------------------*/

/**Función para el cálculo de la media y la varianza de la región de las imágenes encerrada por la malla registrada manualmente */
void statisticRegionImagesClosedMesh( vtkImageData *images, vtkPolyData *liverMesh, double &mean, double &sigma ){
  
  vtkPolyDataNormals *normalsMesh = vtkPolyDataNormals::New();
  normalsMesh->SetInput( liverMesh );
  
  vtkPolyDataToImageStencil *dataToStencil = vtkPolyDataToImageStencil::New();
  dataToStencil->SetInformationInput( images );
  dataToStencil->SetInput( normalsMesh->GetOutput() );
  
  //---------------------------
  //extracción de la región que encierra la malla
/*
  vtkImageStencil *stencil = vtkImageStencil::New();
  stencil->SetInput( intensityImages );
  stencil->SetStencil( dataToStencil->GetOutput() );
  stencil->ReverseStencilOff();
  stencil->SetBackgroundValue( 0 );
  
  stencil->Update();
  
  range[0]=-100;
  range[1]= 400;
  ViewImages( stencil->GetOutput(), range );
  cout << "realizada." << endl;
  //-------  
  stencil->Delete();
*/
  //extracción de la región que encierra la malla
  //---------------------------
  //extracción de la información de la distribución estadística de la región de las imágenes que encierra la malla
  int numberOfVoxels;
  double meanRegion[3], sigmaRegion[3], minValueRegion[3], maxValueRegion[3];

  vtkImageAccumulate *statisticInformation = vtkImageAccumulate::New();
  statisticInformation->SetInput( images );
  statisticInformation->SetStencil( dataToStencil->GetOutput() );
  statisticInformation->ReverseStencilOff();
  statisticInformation->IgnoreZeroOff();
  
  statisticInformation->Update();
  
  statisticInformation->GetMean( meanRegion );
  statisticInformation->GetStandardDeviation( sigmaRegion );
  numberOfVoxels = statisticInformation->GetVoxelCount();
  statisticInformation->GetMin( minValueRegion );
  statisticInformation->GetMax( maxValueRegion );
  
  mean = meanRegion[0];
  sigma = sigmaRegion[0];
  
  statisticInformation->Delete();
  dataToStencil->Delete();
  normalsMesh->Delete();

  return;
}


//------------------------------------
/** función para acutalizar los puntos de la malla aplicando una sola fuerza */
void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, vtkPoints *lNewMeshPoints ){
  
  vtkIdType nop;
  double coordPoint[3], forceVector1[3];
  
  nop = lMeshPoints->GetNumberOfPoints();
  lNewMeshPoints->SetNumberOfPoints( nop );
  
  for(int idMesh=0; idMesh < nop; ++idMesh ){
    
    lForces1->GetTuple( idMesh, forceVector1 );
    
    times(p1, forceVector1, forceVector1);
    
    lMeshPoints->GetPoint( idMesh, coordPoint );
    
    vplus( coordPoint, forceVector1, coordPoint );
    
    lNewMeshPoints->SetPoint( idMesh, coordPoint );
  }

  return;
}

/** función para acutalizar los puntos de la malla aplicando la fuerza de tensión y la fuerza de bordes*/
void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, vtkPoints *lNewMeshPoints ){
  
  vtkIdType nop;
  double coordPoint[3], forceVector1[3], forceVector2[3], sumForces[3];
  
  nop = lMeshPoints->GetNumberOfPoints();
  lNewMeshPoints->SetNumberOfPoints( nop );
  
  for(int idMesh=0; idMesh < nop; ++idMesh ){
    
    lForces1->GetTuple( idMesh, forceVector1 );
    lForces2->GetTuple( idMesh, forceVector2 );
    
    times(p1, forceVector1, forceVector1);
    times(p2, forceVector2, forceVector2);
    
    sumForces[0] = sumForces[1] = sumForces[2] = 0.0;
    
    vplus( sumForces, forceVector1, sumForces );
    vplus( sumForces, forceVector2, sumForces );
    
    lMeshPoints->GetPoint( idMesh, coordPoint );
    
    vplus( coordPoint, sumForces, coordPoint );
    
    lNewMeshPoints->SetPoint( idMesh, coordPoint );
  }

  return;
}

/** función para acutalizar los puntos de la malla aplicando la fuerza de tensión, la fuerza de bordes y la fuerza de intensidad (basada en la fuerza de balloon)*/
void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, double &p3, vtkDoubleArray *lForces3, vtkPoints *lNewMeshPoints ){
  
  vtkIdType nop;
  double coordPoint[3], forceVector1[3], forceVector2[3], forceVector3[3], sumForces[3];
  
  nop = lMeshPoints->GetNumberOfPoints();
  lNewMeshPoints->SetNumberOfPoints( nop );
  
  for(int idMesh=0; idMesh < nop; ++idMesh ){
    
    lForces1->GetTuple( idMesh, forceVector1 );
    lForces2->GetTuple( idMesh, forceVector2 );
    lForces3->GetTuple( idMesh, forceVector3 );
    
    times(p1, forceVector1, forceVector1);
    times(p2, forceVector2, forceVector2);
    times(p3, forceVector3, forceVector3);
    
    sumForces[0] = sumForces[1] = sumForces[2] = 0.0;
    
    vplus( sumForces, forceVector1, sumForces );
    vplus( sumForces, forceVector2, sumForces );
    vplus( sumForces, forceVector3, sumForces );
    
    lMeshPoints->GetPoint( idMesh, coordPoint );
    
    vplus( coordPoint, sumForces, coordPoint );
    
    lNewMeshPoints->SetPoint( idMesh, coordPoint );
  }

  return;
}

/** función para acutalizar los puntos de la malla aplicando la fuerza de tensión, la fuerza de rigidez, la fuerza de bordes y la fuerza de intensidad (basada en la fuerza de balloon)*/
void moveMeshPoints( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, double &p3, vtkDoubleArray *lForces3, double &p4, vtkDoubleArray *lForces4, vtkPoints *lNewMeshPoints ){
  
  vtkIdType nop;
  double coordPoint[3], forceVector1[3], forceVector2[3], forceVector3[3], forceVector4[3], sumForces[3];
  
  nop = lMeshPoints->GetNumberOfPoints();
  lNewMeshPoints->SetNumberOfPoints( nop );
  
  for(int idMesh=0; idMesh < nop; ++idMesh ){
    
    lForces1->GetTuple( idMesh, forceVector1 );
    lForces2->GetTuple( idMesh, forceVector2 );
    lForces3->GetTuple( idMesh, forceVector3 );
    lForces4->GetTuple( idMesh, forceVector4 );
    
    times(p1, forceVector1, forceVector1);
    times(p2, forceVector2, forceVector2);
    times(p3, forceVector3, forceVector3);
    times(p4, forceVector4, forceVector4);
    
    sumForces[0] = sumForces[1] = sumForces[2] = 0.0;
    
    vplus( sumForces, forceVector1, sumForces );
    vplus( sumForces, forceVector2, sumForces );
    vplus( sumForces, forceVector3, sumForces );
    vplus( sumForces, forceVector4, sumForces );
    
    lMeshPoints->GetPoint( idMesh, coordPoint );
    
    vplus( coordPoint, sumForces, coordPoint );
    
    lNewMeshPoints->SetPoint( idMesh, coordPoint );
  }

  return;
}


/** función para acutalizar los puntos de la malla aplicando la fuerza de tensión, la fuerza de bordes y la fuerza de intensidad (basada en la fuerza de balloon)*/
void movePointsMesh( vtkPoints *lMeshPoints, double &p1, vtkDoubleArray *lForces1, double &p2, vtkDoubleArray *lForces2, double &p3, vtkDoubleArray *lForces3, double limitImages[6], vtkPoints *lNewMeshPoints ){
  
  vtkIdType nop;
  double coordPoint[3], forceVector1[3], forceVector2[3], forceVector3[3], sumForces[3];
  
  // copy the input (only polys) to our working mesh
//   vtkPoints *points = vtkPoints::New();
//   vtkPoints *newPoints = vtkPoints::New();
//   vtkCellArray *polys = vtkCellArray::New();
  
  nop = lMeshPoints->GetNumberOfPoints();
//   polys->DeepCopy(inputMesh->GetPolys());

//   nop = points->GetNumberOfPoints();
//   newPoints->SetNumberOfPoints( nop );
//   cout << "Reemplazo de los puntos de la malla..." << endl;
  for(int idMesh=0; idMesh < nop; ++idMesh ){
    
    lForces1->GetTuple( idMesh, forceVector1 );
    lForces2->GetTuple( idMesh, forceVector2 );
    lForces3->GetTuple( idMesh, forceVector3 );
    
    times(p1, forceVector1, forceVector1);
    times(p2, forceVector2, forceVector2);
    times(p3, forceVector3, forceVector3);
    
    sumForces[0] = sumForces[1] = sumForces[2] = 0.0;
    
    vplus( sumForces, forceVector1, sumForces );
    vplus( sumForces, forceVector2, sumForces );
    vplus( sumForces, forceVector3, sumForces );
    
    lMeshPoints->GetPoint(idMesh, coordPoint);
    vplus( coordPoint, sumForces, coordPoint );
    //Es necesario evaluar las coordenadas de los puntos, y si estan por fuera de las imágenes debe ajustarse a los límites de estas.
    if( coordPoint[0] < limitImages[0] || coordPoint[0] > limitImages[1] ||                 coordPoint[1] < limitImages[2] || coordPoint[1] > limitImages[3] ||                 coordPoint[2] < limitImages[4] || coordPoint[2] > limitImages[5] ){
      
      if(coordPoint[0] < limitImages[0]){
        coordPoint[0] = limitImages[0];
      }else{
        if(coordPoint[0] > limitImages[1]){
          coordPoint[0] = limitImages[1];
        }
      };
      if(coordPoint[1] < limitImages[2]){
        coordPoint[1] = limitImages[2];
      }else{
        if(coordPoint[1] > limitImages[3]){
          coordPoint[1] = limitImages[3];
        }
      };
      if(coordPoint[2] < limitImages[4]){
        coordPoint[2] = limitImages[4];
      }else{
        if(coordPoint[2] > limitImages[5]){
          coordPoint[2] = limitImages[5];
        }
      };
    };
    
    lNewMeshPoints->SetPoint(idMesh, coordPoint);
/*    inputMesh->GetPoints()->GetPoint( idMesh, coordPoint );
    vplus( coordPoint, sumForces, coordPoint );
    inputMesh->GetPoints()->SetPoint( idMesh, coordPoint );*/
  
    
  }
  
  
  
//   inputMesh->Reset();
//   inputMesh->Initialize();
//   inputMesh->SetPoints(points);
//   inputMesh->SetPolys(polys);
  
//   inputMesh->BuildCells();
//   inputMesh->BuildLinks();
//   inputMesh->Modified();
//   inputMesh->Update();

  //   outputMesh->Initialize();
/*  outputMesh->SetPoints(newPoints);
  outputMesh->SetPolys(polys);
  
  outputMesh->BuildCells();
  outputMesh->BuildLinks();
  outputMesh->Update();*/
//   cout << "Reemplazo realizado." << endl;
//   polys->Delete();
//   points->Delete();
//   newPoints->Delete();
  
  return;
}


//------------------------------------
/**Función para el cálculo de la fuerza de bordes en función de la función de velocidad*/
void edgeForce( vtkDataArray *lGradientEdgeImages, vtkDataArray *lNormalVectors, vtkIdTypeArray *lIdMeshIdImages, vtkDoubleArray *lEdgeForce ){
  
  vtkIdType idImages;
  double scalarValue;
  double gradientVector[3], normalVector[3], vectorEdgeForce[3];
    
  vtkIdType nop = lIdMeshIdImages->GetNumberOfTuples();
//   cout << "id nop: " << nop << endl;
  
//   lEdgeForce->SetNumberOfTuples( nop );
//   cout << "not: " << lEdgeForce->GetNumberOfTuples() << endl;
  
  for(int idMesh=0; idMesh < nop; idMesh++ ){
    //extracción del vector normal asociado al vértice idMesh de la malla
    lNormalVectors->GetTuple( idMesh, normalVector );
    //para relacionar el idMesh de la malla con el idImages de las imágenes
    idImages = lIdMeshIdImages->GetValue( idMesh );
    //con idImages se extrae el correspondiente vector gradiente asociado al vértice idMesh
    lGradientEdgeImages->GetTuple( idImages, gradientVector );
    //ahora se calculará el producto punto entre el vector normal y el vector gradiente (de la función de velocidad)
    scalarValue = vtimes( gradientVector, normalVector );
    //el scalarValue calculado ponderará el vector normal unitario para establecer la fuerza de bordes
    times( scalarValue, normalVector, vectorEdgeForce );
    //se almacena el vector resultante para cada vértice en lEdgeForce
    lEdgeForce->SetTuple( idMesh, vectorEdgeForce );
//     lEdgeForce->SetTuple( idMesh, gradientVector );
  }
    
  return;
}


//------------------------------------
/**Función para obtener la relación entre los puntos de la malla con los puntos de las imágenes*/
void relationIdMeshIdImages( vtkPoints *lMeshPoints, vtkImageData *intensityImages, vtkIdTypeArray *lIdMeshIdImages ){
    
  //Se obtiene la relación entre los id de los puntos de la malla con los puntos de las imágenes (tanto de intensidad como del gradiente)
  double coordPoint[3];
  vtkIdType nop;
  vtkIdType idImages;
  
  nop = lMeshPoints->GetNumberOfPoints();
  
//   lIdMeshIdImages->SetNumberOfValues( nop );

  for( int idMesh=0; idMesh < nop; idMesh++ ){
      //cout << "\tidMesh=" << idMesh; 
      //GetPoint obtiene las coordenadas globales del punto (coordenadas del espacio) de la Malla indentificado con idMesh
    lMeshPoints->GetPoint(idMesh, coordPoint); 
      //cout << " ("<< coordPoint[0] << ", " << coordPoint[1] << ", " << coordPoint[2] << "), ";
      //FindPoint encuentra el id asociado a las imágenes más cercano a las coodenadas almacenadas en coordPoint
    idImages = intensityImages->FindPoint( coordPoint );
      //idSobel = magnitudeSobelImages->FindPoint( coordPoint );
      //cout << " ( " << idImages << ", " << idSobel << " )\t";
    lIdMeshIdImages->SetValue( idMesh, idImages );
  }

  return;
}


//------------------------------------
/**Función para el cálculo de la fuerza de intensidad basada en la fuerza de balloon propuesta incialmente por Cohen (On Active Contour Models and Balloons)*/
void intensityForce( vtkDataArray *lIntensityImages, vtkDataArray *lNormalVectors, vtkIdTypeArray *lIdMeshIdImages, double &mean, double &sigma, vtkDoubleArray *lIntensityForce, int selFP ){
//Los valores de intensidad de las imágenes serán evaluados en cada vértice de la malla para establecer la magnitud y sentido de la fuerza en dirección del vector normal. Este cálculo depende de la distribución estadística de las intensidades de la región encerrada por la malla registrada manualmente (media y desviación estandar).
  
  double normalVector[3], intensityValue[1], vectorIntensityForce[3];
  double scalarValue;
  vtkIdType nop;
  vtkIdType idImages;
  
  nop = lIdMeshIdImages->GetNumberOfTuples();
//   cout << "\nNúmero de vértices de la malla: " << nop << endl;
//   lIntensityForce->SetNumberOfTuples( nop );
  
  for( int idMesh=0; idMesh < nop; idMesh++ ){
    //se encuentra la relación entre el id de la malla y el id de las imágenes
    idImages = lIdMeshIdImages->GetValue( idMesh );
    //se obtiene el valor de intensidad de las imágenes asociado a cada vértice idMesh
    lIntensityImages->GetTuple( idImages, intensityValue );
//     cout << "\t  intensityValue= " << intensityValue[0];
    //será evaluado el valor de intensidad obtenido con la función  triangular con máximo de 1 en la media, y cruces por cero en +/- 3*sigma, siguiendo la linealidad hasta +/- 6*sigma donde toma el valor de -1.
    
    switch(selFP)
    {
      /**------- Función rectangular -----------**/
      case 1:
	if( intensityValue[0] < mean-1.5*sigma  ||  intensityValue[0] > mean+1.5*sigma ){
	  scalarValue = -1.0;
	}
	else{
	    scalarValue = 1.0;
	}
	break;
      /**------- Función rectangular -----------**/
      /**------- Función triangular -----------**/
      case 2:
	if( intensityValue[0] < mean-3.0*sigma  ||  intensityValue[0] > mean+3.0*sigma ){
	  scalarValue = -1.0;
	}
	else{
	  if( intensityValue[0] >= mean-3.0*sigma  &&  intensityValue[0] < mean ){
	    scalarValue = ( intensityValue[0] - ( mean - 1.5*sigma ) ) / ( 1.5*sigma );
	  }
	  else{
	    scalarValue = -( intensityValue[0] - ( mean + 1.5*sigma ) ) / ( 1.5*sigma );
	  }
	}
	break;
	/**------- Función triangular -----------**/
        /**------- Función trapecio -----------**/
      case 3: // 2.25 = 9.0/4.0
// 	cout << "trapecio" << endl;
	if( intensityValue[0] < mean-2.25*sigma  ||  intensityValue[0] > mean+2.25*sigma ){
	  scalarValue = -1.0;
	}
	else{ // 0.75 = 3.0/4.0
	  if( intensityValue[0] >= mean-2.25*sigma  &&  intensityValue[0] < mean-0.75*sigma ){
	    scalarValue = ( intensityValue[0] - ( mean - 1.5*sigma ) ) / ( 0.75*sigma );
	  }
	  else{
	    if ( intensityValue[0] >= mean-0.75*sigma && intensityValue[0] < mean+0.75*sigma ){
	      scalarValue = 1.0;
	    }
	    else{
	      scalarValue = -( intensityValue[0] - ( mean + 1.5*sigma ) ) / (0.75*sigma);
	    }
	  }
	}
	break;
	/**------- Función trapecio -----------**/
      default:
	cout << "Ninguna función presión elegida." << endl;
	return;
    }

    //se obtiene el vector normal asociado al vértice idMesh de la malla
    lNormalVectors->GetTuple(idMesh, normalVector );
//     cout << "  normalVector: "<< normalVector[0] << ", " << normalVector[1] << ", " << normalVector[2] << "), ";
    //se multiplica el vector normal unitario por el factor de escala scaleValue
    times( scalarValue , normalVector, vectorIntensityForce );
    //se almacena el vector resultante para cada vértice de la malla
    lIntensityForce->SetTuple( idMesh, vectorIntensityForce );
//     cout << "normalVector1: "<< normalVector[0] << ", " << normalVector[1] << ", " << normalVector[2] << "), ";
  }

  return;
}

//------------------------------------
/**Función para establecer los vértices vecinos de cada vértice de la malla*/
void findNeighbors1 (vtkPolyData *input, vtkIdList **nbs) {
  
  vtkIdList *cells = vtkIdList::New();  //list to store cells
  vtkIdList *points = vtkIdList::New(); //list to store the cells points
  vtkIdType nop, nopc, noc, cell, point;
  
  nop = input->GetNumberOfPoints();
  //find neighbors for every vertex 
  for(int i=0; i<nop; i++){
    nbs[i] = vtkIdList::New();
    //obtain the cell list
    input->GetPointCells(i, cells);
    noc = cells->GetNumberOfIds();
    // for each cell, obtain the point list
    for(int j=0; j<noc; j++){
      cell = cells->GetId(j);
      input->GetCellPoints(cell, points);
      nopc = points->GetNumberOfIds();
      // for each point, try to insert it in the neighbors list
      for(int k=0; k<nopc; k++){
        point = points->GetId(k);
        if ( point != i  &&  nbs[i]->IsId(point) == -1 ) { 
          nbs[i]->InsertNextId(point);
        }
      }
    }
  }
  
  cells->Delete();
  points->Delete();

  return;
}


//------------------------------------
/**Función para estimar los centroides de cada vértice en función de las coordenadas de sus vecinos*/
void findCentroids( vtkPolyData *input, vtkIdList **lnbs, vtkDoubleArray *lcentroids ) {
  
  vtkIdType nop, non, nb;
  double x[3], centroid[3];
  double suma[3];
  
  nop = input->GetNumberOfPoints();
//   lcentroids->SetNumberOfTuples( nop );
  
  for(int i=0; i<nop; i++){
    non = lnbs[i]->GetNumberOfIds();
//     cout << "non= " << non << endl;
    suma[0] = suma[1] = suma[2] = 0;
    for(int j=0; j<non; j++){
      nb = lnbs[i]->GetId(j);
//       cout << "nb= " << nb << endl;
      input->GetPoint(nb, x);
//       cout << "x= " << x[0] <<", "<< x[1] <<", "<< x[2] << endl;
      vplus(suma, x, suma);
//       cout << "suma= " << suma[0] <<", "<< suma[1] <<", "<< suma[2] << endl;
    }
    divide(suma, non, centroid);
//     cout << i << " centroid= " << centroid[0] <<", "<< centroid[1] <<", "<< centroid[2] << endl;
//     lcentroids->InsertTupleValue(i, centroid);
    lcentroids->SetTuple(i, centroid);
  }
  
  return;
}


//------------------------------------
/**Función para el cálculo de la fuerza de tensión*/
void tensionForce( vtkPolyData *input, vtkDoubleArray *lcentroids, vtkDoubleArray *lTForce ){
  
  vtkIdType nop;
  double centroid[3], x[3], force[3];
  
  nop = input->GetNumberOfPoints();
//   lTForce->SetNumberOfTuples( nop );
  
  for(int i=0; i<nop; i++){
    lcentroids->GetTuple(i, centroid);
    input->GetPoint(i, x);
    vminus (centroid, x, force);
  //  lTForce->InsertTupleValue(i, force);
    lTForce->SetTuple(i, force);
  }

  return;
}


void rigidityForce1( vtkPolyData *input, vtkDoubleArray *lTForce, vtkIdList **lnbs, vtkDoubleArray *lRForce ){
  
  vtkIdType nop, non, nb;
  double ntforce[3], suma[3], vtforce[3], force[3];
  
  nop = input->GetNumberOfPoints();
  lRForce->SetNumberOfTuples( nop );
  
  for(int i=0; i<nop; i++){
    non = lnbs[i]->GetNumberOfIds();
    suma[0] = suma[1] = suma[2] = 0.0;
    for(int j=0; j<non; j++){
      nb = lnbs[i]->GetId(j);
      lTForce->GetTuple(nb, ntforce);
      vplus(suma, ntforce, suma);
    }
    divide(suma, non, ntforce);
    lTForce->GetTuple(i, vtforce);
    vminus(vtforce, ntforce, force);
    lRForce->SetTuple(i, force);
  }

  return;
}

/** Esta fuerza se propone con base en la tesis de Andreita  */
void rigidityForce2( vtkPolyData *input, vtkDoubleArray *lTForce, vtkDataArray *lNormalVectors, vtkDoubleArray *lRForce ){
  
  vtkIdType nop;
  double normalVector[3], vectorTensionForce[3], vectorRigidityForce[3];
  double scalarValue;
  
  nop = input->GetNumberOfPoints();
  lRForce->SetNumberOfTuples( nop );
  
  for(int idMesh=0; idMesh<nop; idMesh++){
    //  vector normal asociado a vertice idMesh
    lNormalVectors->GetTuple( idMesh, normalVector );
    //  vector de la fuerza de tensión
    lTForce->GetTuple(idMesh, vectorTensionForce);
    //ahora se calculará el producto punto entre el vector normal y el vector de la fuerza de tensión
    scalarValue = vtimes( vectorTensionForce, normalVector );
    //el scalarValue calculado ponderará el vector normal unitario para establecer la fuerza de bordes
    times( scalarValue, normalVector, vectorRigidityForce );
    //se almacena el vector resultante para cada vértice en lEdgeForce
    lRForce->SetTuple(idMesh, vectorRigidityForce);
  }

  return;
}
