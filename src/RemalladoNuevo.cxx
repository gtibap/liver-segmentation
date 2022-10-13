//La fusión de vértices se realiza para aquellos cuyo borde se menor a un valor establecido.
// La estrategia consiste en encontrar las celdas a borrar, cambiar los puntos de las celdas que comparten los puntos de las celdas a borrar, y solo cuando se haya recorrido toda la malla, las celdas y los puntos seleccionados serán borrados.
 
// Autor: Gerardo Tibamoso Pedraza. 

// librerias VTK
#include "vtkCell.h"
#include "vtkCellArray.h"
#include "vtkCellData.h"
#include "vtkCleanPolyData.h"
#include "vtkDoubleArray.h"
#include "vtkIdList.h"
#include "vtkMath.h"
#include "vtkPointData.h"
#include "vtkPolyData.h"
#include "vtkPolyDataNormals.h"
#include "vtkSmoothPolyDataFilter.h"
#include "vtkSTLReader.h"
#include "vtkSTLWriter.h"
#include "vtkTriangleFilter.h"
#include "vtkWindowedSincPolyDataFilter.h"


// librerias básicas de C++
#include <string.h>
#include <iostream>
using namespace std;

double distance (double A[3], double B[3]);
void vplus (double A[3], double B[3], double C[3]);
void divide (double A[3], double s, double C[3]);

void MeshVisualization( vtkPolyData *mesh );
void NormalsMeshVisualization( vtkPolyData *mesh );

int calculateErrorDistance( vtkPoints *initialPoints, vtkPoints *pointsSmoothMesh, vtkDoubleArray *errorDistance);

int reductionPointsMesh(vtkPolyData *meshLiver, double minimumDistance, vtkPolyData *outputMesh);
int CollapseEdge(vtkPolyData *Mesh, vtkIdType pt0Id, vtkIdType pt1Id);

int insertPointsMesh( vtkPolyData *inputMesh, double maximumDistance, vtkPolyData *outputMesh );
int insertNewCells(vtkPolyData *meshInput, vtkIdType idPoint1, vtkIdType idPoint2, vtkIdType idNewPoint);

/** Inicio función del remallado de la malla que representa el hígado para ajustarla a las imágenes particulares */
/** Inicio función del remallado de la malla que representa el hígado para ajustarla a las imágenes particulares */

/** Redistribuir los vértices sobre la malla con una operación de suavizado */
int initialRemesh( vtkPolyData *input, float minimumDistance )
{
  vtkPoints *pointsMesh = vtkPoints::New();
  vtkCellArray *polysMesh = vtkCellArray::New();
  vtkPolyData *temporalMesh = vtkPolyData::New();
  
  pointsMesh->DeepCopy( input->GetPoints() );
  polysMesh->DeepCopy( input->GetPolys() );
  
  temporalMesh->SetPoints( pointsMesh );
  pointsMesh->Delete();
  temporalMesh->SetPolys( polysMesh );
  polysMesh->Delete();
  temporalMesh->GetPointData()->DeepCopy(input->GetPointData());
  temporalMesh->GetFieldData()->PassData(input->GetFieldData());
  temporalMesh->BuildCells();
  temporalMesh->BuildLinks();
  
  vtkWindowedSincPolyDataFilter *initialSmoothMesh = vtkWindowedSincPolyDataFilter::New();
  initialSmoothMesh->NormalizeCoordinatesOn();
  initialSmoothMesh->SetNumberOfIterations( 10 );
  initialSmoothMesh->SetInput( temporalMesh );
  
  initialSmoothMesh->Update();
  
  // copy the input (only polys) to our working mesh
  vtkPolyData *temporalMesh1 = vtkPolyData::New();
  vtkPoints *pointsMesh1 = vtkPoints::New();
  vtkCellArray *polysMesh1 = vtkCellArray::New();

  pointsMesh1->DeepCopy(initialSmoothMesh->GetOutput()->GetPoints());
  polysMesh1->DeepCopy(initialSmoothMesh->GetOutput()->GetPolys());
  
  temporalMesh1->SetPoints(pointsMesh1);
  pointsMesh1->Delete();
  temporalMesh1->SetPolys(polysMesh1);
  polysMesh1->Delete();
  temporalMesh1->GetPointData()->DeepCopy(initialSmoothMesh->GetOutput()->GetPointData());
  temporalMesh1->GetFieldData()->PassData(initialSmoothMesh->GetOutput()->GetFieldData());
  temporalMesh1->BuildCells();
  temporalMesh1->BuildLinks();
  
  initialSmoothMesh->Delete();
  temporalMesh->DeleteLinks();
  temporalMesh->Delete();

  //  Lectura de cada uno de los triángulos de la malla y determinación de la longitud de cada uno de sus lados
  
  vtkIdType noc, nop;
  double point0[3], point1[3], point2[3], centerPoint[3];
  double edge0, edge1, edge2;
  double suma[3];
  int cuenta0, cuenta1, cuenta2, cuenta3, cuenta4;

  
  cuenta0=cuenta1=cuenta2=cuenta3=cuenta4=0;
  
  vtkIdList *idPoints = vtkIdList::New();

  nop = temporalMesh1->GetNumberOfPoints();
//   cout << "pointsMesh = " << nop << endl;
  noc = temporalMesh1->GetNumberOfCells();
//   cout << "cellsMesh = " << noc << endl;

  for( int idCell=0; idCell < noc; idCell++ )
  {
    
    if( temporalMesh1->GetCell(idCell)->GetCellType() != VTK_EMPTY_CELL )
    {
      temporalMesh1->GetCellPoints( idCell, idPoints );
      
  //     cout << "numberOfIds cell" << idCell << "=" << idPoints->GetNumberOfIds() << endl;
      
      temporalMesh1->GetPoint( idPoints->GetId(0), point0 );
  //     cout << "point0: " << point0[0] << ", " << point0[1] << ", " << point0[2] << endl;
      temporalMesh1->GetPoint( idPoints->GetId(1), point1 );
  //     cout << "point1: " << point1[0] << ", " << point1[1] << ", " << point1[2] << endl;
      temporalMesh1->GetPoint( idPoints->GetId(2), point2 );
  //     cout << "point2: " << point2[0] << ", " << point2[1] << ", " << point2[2] << endl;
      
      //cálculo de las long. de los lados del tríangulo
      edge0 = distance( point0, point1);
      edge1 = distance( point1, point2);
      edge2 = distance( point2, point0);
//       cout << "distances: " << endl;
//       cout << "edge0 = " << edge0 << endl;
//       cout << "edge1 = " << edge1 << endl;
//       cout << "edge2 = " << edge2 << endl;
      
      suma[0] = suma[1] = suma[2] = 0.0;
      if( edge0 < minimumDistance || edge1 < minimumDistance || edge2 < minimumDistance )
      {
        if( edge0 < minimumDistance )
        {
          vplus(point0, point1, suma);
          divide(suma, 2, centerPoint);

          temporalMesh1->GetPoints()->SetPoint( idPoints->GetId(0), centerPoint );

          CollapseEdge(temporalMesh1, idPoints->GetId(0), idPoints->GetId(1)); 

          cuenta1++;
        }
        else
        {
          if( edge1 < minimumDistance )
          {
            vplus(point1, point2, suma);
            divide(suma, 2, centerPoint);

            temporalMesh1->GetPoints()->SetPoint( idPoints->GetId(1), centerPoint );

            CollapseEdge(temporalMesh1, idPoints->GetId(1), idPoints->GetId(2)); 

            cuenta2++;
          }
          else
          {
            vplus(point2, point0, suma);
            divide(suma, 2, centerPoint);

            temporalMesh1->GetPoints()->SetPoint( idPoints->GetId(2), centerPoint );

            CollapseEdge(temporalMesh1, idPoints->GetId(2), idPoints->GetId(0)); 

            cuenta3++;
          }
        }
      }
    }
  }
  
//   cout << "Remallado realizado." << endl;
  
  idPoints->Delete();
  /** Aquí se pasan las celdas de la malla modificada a una nueva malla */
    
  vtkIdList *outputCellList = vtkIdList::New();
  vtkPolyData *output = vtkPolyData::New();
  
  // copy the simplified mesh from the working mesh to the output mesh
  for (int i = 0; i < temporalMesh1->GetNumberOfCells(); i++) 
  {
    if (temporalMesh1->GetCell(i)->GetCellType() != VTK_EMPTY_CELL) 
    {
      outputCellList->InsertNextId(i);
    }
  } 

  
  output->Reset();
  output->Allocate(temporalMesh1, outputCellList->GetNumberOfIds());
  output->GetPointData()->CopyAllocate(temporalMesh1->GetPointData(), 1);
  output->CopyCells(temporalMesh1, outputCellList);
  
  
//   nop = output->GetNumberOfPoints();
//   cout << "out-pointsMesh = " << nop << endl;
  
//   noc = output->GetNumberOfCells();
//   cout << "out-cellsMesh = " << noc << endl;

  /**-------------------------------------*/
  /** Escritura de la malla resultante */
/*  int flagw=0;
  cout << "Escribir malla?(1/0): ";
  cin >> flagw;
  if(flagw==1)
  {
    
    cout << "Escritura de la malla en la función initialRemesh..." << endl;
    vtkPolyDataWriter *writerMesh = vtkPolyDataWriter::New();
    writerMesh->SetFileTypeToASCII();
    writerMesh->SetFileName("../data/mesh02.vtk");
    writerMesh->SetInput(output);
    try{
      writerMesh->Write();
    }
    catch( exception& ){
      cerr << "error." << endl;
      return EXIT_FAILURE;
    }
    cout << "realizada." << endl;
    writerMesh->Delete();
  }*/
  /**-------------------------------------*/
  /** Redistribuir los vértices sobre la malla con una operación de suavizado */
  vtkWindowedSincPolyDataFilter *smoothMesh = vtkWindowedSincPolyDataFilter::New();
  smoothMesh->NormalizeCoordinatesOn();
  smoothMesh->SetNumberOfIterations( 10 );
  smoothMesh->SetInput( output );
  
  smoothMesh->Update();
  
  /** pasar los puntos de nuevo a la entrada  */
  vtkPoints *pointsMesh2 = vtkPoints::New();
  vtkCellArray *polysMesh2 = vtkCellArray::New();

  pointsMesh2->DeepCopy(smoothMesh->GetOutput()->GetPoints());
  polysMesh2->DeepCopy(smoothMesh->GetOutput()->GetPolys());
  
  input->Reset();
  input->SetPoints(pointsMesh2);
  pointsMesh2->Delete();
  input->SetPolys(polysMesh2);
  polysMesh2->Delete();
  input->GetPointData()->DeepCopy(smoothMesh->GetOutput()->GetPointData());
  input->GetFieldData()->PassData(smoothMesh->GetOutput()->GetFieldData());
  input->BuildCells();
  input->BuildLinks();
  
  smoothMesh->Delete();
  output->Delete();
  temporalMesh1->DeleteLinks();
  temporalMesh1->Delete();
  outputCellList->Delete();
  
  return 0;

}


/**-------------------------------------
Inicio función del remallado de la malla que representa el hígado para ajustarla a las imágenes particulares 
*/
int remesh( vtkPolyData *inputMesh, float minimumDistance, vtkPoints *initialPoints, vtkDoubleArray *errorDistance, vtkPolyData *newOutputMesh )
{
  vtkPoints *points;
  vtkPolyData *workMesh;
  vtkPolyData *outputMesh;
  vtkCellArray *polys;
  vtkWindowedSincPolyDataFilter *smoothMesh;
  vtkWindowedSincPolyDataFilter *smoothMesh1;

  // Redistribuir los vértices sobre la malla con una operación de suavizado 
//   cout << "smooth1...";
  smoothMesh = vtkWindowedSincPolyDataFilter::New();
  smoothMesh->NormalizeCoordinatesOn();
  smoothMesh->SetNumberOfIterations( 10 );
  smoothMesh->SetInput( inputMesh );
  
  smoothMesh->Update();
//   cout << " realizado." << endl;

  points = vtkPoints::New();
  polys = vtkCellArray::New();
  
  // copia de la malla de entrada (solo los polys) a la malla de trabajo  
  workMesh = vtkPolyData::New();
  points->DeepCopy( smoothMesh->GetOutput()->GetPoints() );
  workMesh->SetPoints( points );
  
  // Cálculo de la mag. de las distancias de desplazamiento de los puntos de la malla 
  calculateErrorDistance( initialPoints, points, errorDistance );

  points->Delete();
  
  polys->DeepCopy( smoothMesh->GetOutput()->GetPolys() );
  workMesh->SetPolys( polys );
  polys->Delete();
  workMesh->GetPointData()->DeepCopy(smoothMesh->GetOutput()->GetPointData());
  workMesh->GetFieldData()->PassData(smoothMesh->GetOutput()->GetFieldData());
  workMesh->BuildCells();
  workMesh->BuildLinks();

  smoothMesh->Delete();
  
  // Fusión de puntos con distancias menores a la distancia mínima 2mm distancia mín. entre vértices
  int npd, npi;
//   cout << "Fusión de puntos...";

  outputMesh = vtkPolyData::New();
  npd = reductionPointsMesh( workMesh, minimumDistance, outputMesh );
//   cout << " realizada.\n Puntos borrados = " << npd << endl;

  workMesh->Delete();
  workMesh = vtkPolyData::New();
  points = vtkPoints::New();
  polys = vtkCellArray::New();
  
  // copia de la malla de entrada (solo los polys) a la malla de trabajo  
  points->DeepCopy( outputMesh->GetPoints() );
  workMesh->SetPoints( points );
  points->Delete();
  polys->DeepCopy( outputMesh->GetPolys() );
  workMesh->SetPolys( polys );
  polys->Delete();
  workMesh->GetPointData()->DeepCopy(outputMesh->GetPointData());
  workMesh->GetFieldData()->PassData(outputMesh->GetFieldData());
  workMesh->BuildCells();
  workMesh->BuildLinks();
  
  outputMesh->Delete();

  //Inclusión de puntos en bordes con distancias mayores a la distancia máxima, y por tanto, inclusión de nuevas celdas
  
  double maximumDistance = 5.0; // 5mm de distancia máx. entre vértices
  outputMesh = vtkPolyData::New();
  
 // función para insertar puntos y celdas de acuerdo con la long. de los bordes
//   cout << "Inclusión de puntos y celdas...";
  npi = insertPointsMesh( workMesh, maximumDistance, outputMesh );
//   cout << "Realizada. \nPuntos insertados = " << npi << endl;
  
  workMesh->Delete();
  // copia de la malla de salida a la de trabajo (solo los polys)
/*  workMesh = vtkPolyData::New();
  points = vtkPoints::New();
  polys = vtkCellArray::New();

  points->DeepCopy( outputMesh->GetPoints() );
  workMesh->SetPoints( points );
  points->Delete();
  polys->DeepCopy( outputMesh->GetPolys() );
  workMesh->SetPolys( polys );
  polys->Delete();
  workMesh->GetFieldData()->PassData(outputMesh->GetFieldData());
  workMesh->BuildCells();
  workMesh->BuildLinks();

  outputMesh->Delete();  */
  
  // Redistribuir los vértices sobre la malla con una operación de suavizado
//   cout << "smooth2...";
  smoothMesh1 = vtkWindowedSincPolyDataFilter::New();
  smoothMesh1->NormalizeCoordinatesOn();
  smoothMesh1->SetNumberOfIterations( 10 );
  smoothMesh1->SetInput( outputMesh );
  
  smoothMesh1->Update();
//   cout << " realizado." << endl;
  outputMesh->Delete();
//   workMesh->DeleteLinks();
//   workMesh->Delete();
//   cout << " paso1 ";
    
  points = vtkPoints::New();
  polys = vtkCellArray::New();

  newOutputMesh->Reset();
//   newOutputMesh->Allocate(smoothMesh1->GetOutput());
  // copia de la malla de salida a la de entrada (solo los polys)
  points->DeepCopy( smoothMesh1->GetOutput()->GetPoints() );
  newOutputMesh->SetPoints( points );
  points->Delete();
  polys->DeepCopy( smoothMesh1->GetOutput()->GetPolys() );
  newOutputMesh->SetPolys( polys );
  polys->Delete();
  newOutputMesh->GetPointData()->DeepCopy(smoothMesh1->GetOutput()->GetPointData());
  newOutputMesh->GetFieldData()->PassData(smoothMesh1->GetOutput()->GetFieldData());
  newOutputMesh->BuildCells();
  newOutputMesh->BuildLinks();
  
//   cout << "paso2 ";
  smoothMesh1->Delete();

  return 0;
}

/**----------------------------------------
Función para calcular el desplazamiento de cada uno de los puntos dado por la aplicación de las fuerzas de deformación y suavizado. 
*/
int calculateErrorDistance( vtkPoints *initialPoints, vtkPoints *pointsSmoothMesh, vtkDoubleArray *errorDistance)
{
  vtkIdType numberOfPointsMesh;
  double initialCoordPoint[3], finalCoordPoint[3];
  double dist;
    
  numberOfPointsMesh = initialPoints->GetNumberOfPoints();
  //   cout << "numberPoints1= " << numberOfPointsMesh << endl;
  //   cout << "numberPoints2= " << points->GetNumberOfPoints() << endl;
    
  for(int idPoint=0; idPoint < numberOfPointsMesh; idPoint++)
  {
    initialPoints->GetPoint(idPoint, initialCoordPoint);
    pointsSmoothMesh->GetPoint(idPoint, finalCoordPoint);
    dist = distance(initialCoordPoint, finalCoordPoint);
    errorDistance->SetValue(idPoint, dist);
  }
    
  return 0;
}

/**-------------------------------------------
Función para la fusión de puntos con distancia menor a minimumDistance  
*/
int reductionPointsMesh(vtkPolyData *inputMesh, double minimumDistance, vtkPolyData *outputMesh)
{  
  //  Lectura de cada uno de los triángulos de la malla y determinación de la longitud de cada uno de sus lados
  vtkIdType noc, nop;
  double point0[3], point1[3], point2[3], centerPoint[3];
  double edge0, edge1, edge2;
  double suma[3];
  int cuenta0, cuenta1, cuenta2, cuenta3, cuenta4;
  
  // copy the input (only polys) to our working mesh
  vtkPolyData *temporalMesh = vtkPolyData::New();
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *polys = vtkCellArray::New();
  
  points->DeepCopy( inputMesh->GetPoints() );
  polys->DeepCopy( inputMesh->GetPolys() );
  
  temporalMesh->SetPoints( points );
  points->Delete();
  temporalMesh->SetPolys( polys );
  polys->Delete();
  temporalMesh->GetPointData()->DeepCopy(inputMesh->GetPointData());
  temporalMesh->GetFieldData()->PassData(inputMesh->GetFieldData());
  temporalMesh->BuildCells();
  temporalMesh->BuildLinks();
    
  cuenta0=cuenta1=cuenta2=cuenta3=cuenta4=0;
    
  vtkIdList *idPoints = vtkIdList::New();
    
  nop = temporalMesh->GetNumberOfPoints();
  //   cout << "pointsMesh = " << nop << endl;
  noc = temporalMesh->GetNumberOfCells();
  //   cout << "cellsMesh = " << noc << endl;
  
  for( int idCell=0; idCell < noc; idCell++ )
  {
      
    if( temporalMesh->GetCell(idCell)->GetCellType() != VTK_EMPTY_CELL )
    {
      temporalMesh->GetCellPoints( idCell, idPoints );
        
    //     cout << "numberOfIds cell" << idCell << "=" << idPoints->GetNumberOfIds() << endl;
        
      temporalMesh->GetPoint( idPoints->GetId(0), point0 );
    //     cout << "point0: " << point0[0] << ", " << point0[1] << ", " << point0[2] << endl;
      temporalMesh->GetPoint( idPoints->GetId(1), point1 );
    //     cout << "point1: " << point1[0] << ", " << point1[1] << ", " << point1[2] << endl;
      temporalMesh->GetPoint( idPoints->GetId(2), point2 );
    //     cout << "point2: " << point2[0] << ", " << point2[1] << ", " << point2[2] << endl;
        
        //cálculo de las long. de los lados del tríangulo
      edge0 = distance( point0, point1);
      edge1 = distance( point1, point2);
      edge2 = distance( point2, point0);
  //       cout << "distances: " << endl;
  //       cout << "edge0 = " << edge0 << endl;
  //       cout << "edge1 = " << edge1 << endl;
  //       cout << "edge2 = " << edge2 << endl;
        
      suma[0] = suma[1] = suma[2] = 0.0;
      if( edge0 < minimumDistance || edge1 < minimumDistance || edge2 < minimumDistance )
      {
        if( edge0 < minimumDistance )
        {
          vplus(point0, point1, suma);
          divide(suma, 2, centerPoint);
  
          temporalMesh->GetPoints()->SetPoint( idPoints->GetId(0), centerPoint );
  
          CollapseEdge(temporalMesh, idPoints->GetId(0), idPoints->GetId(1)); 
  
  //           cuenta1++;
        }
        else
        {
          if( edge1 < minimumDistance )
          {
            vplus(point1, point2, suma);
            divide(suma, 2, centerPoint);
  
            temporalMesh->GetPoints()->SetPoint( idPoints->GetId(1), centerPoint );
  
            CollapseEdge(temporalMesh, idPoints->GetId(1), idPoints->GetId(2)); 
  
  //             cuenta2++;
          }
          else
          {
            vplus(point2, point0, suma);
            divide(suma, 2, centerPoint);
  
            temporalMesh->GetPoints()->SetPoint( idPoints->GetId(2), centerPoint );
  
            CollapseEdge(temporalMesh, idPoints->GetId(2), idPoints->GetId(0)); 
  
  //             cuenta3++;
          }
        }
        cuenta1++;
      }
    }
  }
//   cout << cuenta1 << " puntos eliminados." << endl;
  //   cout << "Remallado realizado." << endl;
  idPoints->Delete();
  
  // Aquí se pasan las celdas de la malla modificada a una nueva malla 
      
  vtkIdList *outputCellList = vtkIdList::New();
    
    // copy the simplified mesh from the working mesh to the output mesh
  for (int i = 0; i < temporalMesh->GetNumberOfCells(); i++) 
  {
    if (temporalMesh->GetCell(i)->GetCellType() != VTK_EMPTY_CELL) 
    {
      outputCellList->InsertNextId(i);
    }
  } 
  
//   vtkPolyData *outputMesh = vtkPolyData::New();
  
  outputMesh->Reset();
  outputMesh->Allocate(temporalMesh, outputCellList->GetNumberOfIds());
  outputMesh->GetPointData()->CopyAllocate(temporalMesh->GetPointData(), 1);
  outputMesh->CopyCells(temporalMesh, outputCellList);
  
/*  vtkPoints *points1 = vtkPoints::New();
  vtkCellArray *polys1 = vtkCellArray::New();
  
  points1->DeepCopy( outputMesh->GetPoints() );
  polys1->DeepCopy( outputMesh->GetPolys() );
  
  inputMesh->Reset();
  inputMesh->SetPoints( points1 );
  points1->Delete();
  inputMesh->SetPolys( polys1 );
  polys1->Delete();
  inputMesh->GetPointData()->DeepCopy(outputMesh->GetPointData());
  inputMesh->GetFieldData()->PassData(outputMesh->GetFieldData());
  inputMesh->BuildCells();
  inputMesh->BuildLinks();
  
  outputMesh->Delete();*/
  temporalMesh->DeleteLinks();
  temporalMesh->Delete();
  outputCellList->Delete();
    
  return cuenta1;
}

/**-------------------------------------------
Función CollapseEdge de vtkQuadricDecimation*/

int CollapseEdge(vtkPolyData *Mesh, vtkIdType pt0Id, vtkIdType pt1Id) 
{
  int j, numDeleted=0;
  vtkIdType i, npts, *pts, cellId;
  vtkIdList *CollapseCellIds = vtkIdList::New();
//   int cont;
  
  Mesh->GetPointCells(pt0Id, CollapseCellIds);
  i = 0;
  while ( i < CollapseCellIds->GetNumberOfIds() && numDeleted < 2 ) 
  {
    cellId = CollapseCellIds->GetId(i);
    Mesh->GetCellPoints(cellId, npts, pts);
    for (j = 0; j < 3; j++) 
    {
      if (pts[j] == pt1Id) 
      {
        Mesh->RemoveCellReference(cellId);
        Mesh->DeleteCell(cellId);
        numDeleted++;
      }
    }
    i++;
  }

  Mesh->GetPointCells(pt1Id, CollapseCellIds);
  Mesh->ResizeCellList(pt0Id, CollapseCellIds->GetNumberOfIds());
  for (i=0; i < CollapseCellIds->GetNumberOfIds(); i++)
  {
    cellId = CollapseCellIds->GetId(i);
    Mesh->GetCellPoints(cellId, npts, pts);
    // making sure we don't already have the triangle we're about to
    // change this one to
    if ((pts[0] == pt1Id && Mesh->IsTriangle(pt0Id, pts[1], pts[2])) ||
         (pts[1] == pt1Id && Mesh->IsTriangle(pts[0], pt0Id, pts[2])) ||
         (pts[2] == pt1Id && Mesh->IsTriangle(pts[0], pts[1], pt0Id)))
    {
      Mesh->RemoveCellReference(cellId);
      Mesh->DeleteCell(cellId);
      numDeleted++;
//       cout << "entré!" << endl; cin >> cont; 
    }
    else
    {
      Mesh->AddReferenceToCell(pt0Id, cellId);
      Mesh->ReplaceCellPoint(cellId, pt1Id, pt0Id);
    }
  }
  Mesh->DeletePoint(pt1Id);
  
  CollapseCellIds->Delete();
  return numDeleted;
}

/**---------------------------------------------
Función para incertar puntos en la malla*/

int insertPointsMesh( vtkPolyData *inputMesh, double maximumDistance, vtkPolyData *outputMesh )
{
  //  Algunas cosas de esta parte están basadas en vtkInterpolatingSubdivisionFilter.cxx
  //  Lectura de cada uno de los triángulos de la malla y determinación de la longitud de cada uno de sus lados
  vtkIdType noc, nop;
  double point0[3], point1[3], point2[3], centerPoint[3];
  double edge0, edge1, edge2;
  double suma[3];
  int cuenta0, cuenta1, cuenta2, cuenta3, cuenta4;
  vtkIdType idNewPoint;
  vtkPolyData *temporalMesh;
  vtkPoints *points;
  vtkCellArray *polys;
  vtkIdType numCells;
  
  numCells = inputMesh->GetNumberOfCells();
  
  points = vtkPoints::New();
  points->GetData()->DeepCopy(inputMesh->GetPoints()->GetData());

  // se reserva un espacio de memoria para inclusión de puntos y celdas 
  temporalMesh = vtkPolyData::New();
  temporalMesh->GetPointData()->CopyAllocate(inputMesh->GetPointData(), 2*inputMesh->GetNumberOfPoints());
  temporalMesh->GetCellData()->CopyAllocate(inputMesh->GetCellData(),4*numCells);
  
  polys = vtkCellArray::New();
  temporalMesh->GetPolys()->Allocate(polys->EstimateSize(4*numCells,3));
  
  //  copy the input (only polys) to our working mesh  
  points->DeepCopy( inputMesh->GetPoints() );
  temporalMesh->SetPoints( points );
  points->Delete();
  polys->DeepCopy( inputMesh->GetPolys() );
  temporalMesh->SetPolys( polys );
  polys->Delete();
  temporalMesh->GetPointData()->DeepCopy(inputMesh->GetPointData());
  temporalMesh->GetFieldData()->PassData(inputMesh->GetFieldData());
  temporalMesh->BuildCells();
  temporalMesh->BuildLinks();
//   cout << "copia de la malla exitosa." << endl;
  
  cuenta0=cuenta1=cuenta2=cuenta3=cuenta4=0;
  
  vtkIdList *idPoints = vtkIdList::New();
  
  nop = temporalMesh->GetNumberOfPoints();
//   cout << "pointsMesh = " << nop << endl;
  
  noc = temporalMesh->GetNumberOfCells();
//   cout << "cellsMesh = " << noc << endl;
  
  for( int idCell=0; idCell < noc; idCell++ )
  {
    temporalMesh->GetCellPoints( idCell, idPoints );
      
  //     cout << "numberOfIds cell" << idCell << "=" << idPoints->GetNumberOfIds() << endl;
      
    temporalMesh->GetPoint( idPoints->GetId(0), point0 );
  //     cout << "point0: " << point0[0] << ", " << point0[1] << ", " << point0[2] << endl;
    temporalMesh->GetPoint( idPoints->GetId(1), point1 );
  //     cout << "point1: " << point1[0] << ", " << point1[1] << ", " << point1[2] << endl;
    temporalMesh->GetPoint( idPoints->GetId(2), point2 );
  //     cout << "point2: " << point2[0] << ", " << point2[1] << ", " << point2[2] << endl;
      
      //cálculo de las long. de los lados del tríangulo
    edge0 = distance( point0, point1);
    edge1 = distance( point1, point2);
    edge2 = distance( point2, point0);
//       cout << "distances: " << endl;
//       cout << "edge0 = " << edge0 << endl;
//       cout << "edge1 = " << edge1 << endl;
//       cout << "edge2 = " << edge2 << endl;
      
    suma[0] = suma[1] = suma[2] = 0.0;
    if( edge0 > maximumDistance || edge1 > maximumDistance || edge2 > maximumDistance )
    {
      if( edge0 > maximumDistance )
      {
        //  cálculo del punto medio
        vplus(point0, point1, suma);
        divide(suma, 2, centerPoint);
//         cout << "centerPoint: " << centerPoint[0] << centerPoint[1] << centerPoint[2] << endl;
//         cout << "insertando un punto a la malla..." << endl; 
        idNewPoint = temporalMesh->GetPoints()->InsertNextPoint( centerPoint );

//         temporalMesh->GetPoint(idNewPoint, pointNew);
//         cout << "pointNew: " << pointNew[0] << pointNew[1] << pointNew[2] << endl;
        
//         cout << "insertando celdas a la malla..." << endl;
        insertNewCells(temporalMesh, idPoints->GetId(0), idPoints->GetId(1), idNewPoint);
//         cout << "celdas insertadas." << endl;
      }
      else
      {
        if( edge1 > maximumDistance )
        {
          //  cálculo del punto medio
          vplus(point1, point2, suma);
          divide(suma, 2, centerPoint);

//           cout << "insertando un punto a la malla..." << endl;
          
          idNewPoint = temporalMesh->GetPoints()->InsertNextPoint( centerPoint );

//           cout << "insertando celdas a la malla..." << endl;
          insertNewCells(temporalMesh, idPoints->GetId(1), idPoints->GetId(2), idNewPoint);
//           cout << "celdas insertadas." << endl;
        }
        else
        {
          //  cálculo del punto medio
          vplus(point2, point0, suma);
          divide(suma, 2, centerPoint);

//           cout << "insertando un punto a la malla..." << endl;
          idNewPoint = temporalMesh->GetPoints()->InsertNextPoint( centerPoint );

//           cout << "insertando celdas a la malla..." << endl;
          insertNewCells(temporalMesh, idPoints->GetId(2), idPoints->GetId(0), idNewPoint);
//           cout << "celdas insertadas." << endl;
        }
      }
      cuenta1++;
    }
  }
  
  temporalMesh->Squeeze();
  
  outputMesh->Reset();
  
  points = vtkPoints::New();
  points->DeepCopy( temporalMesh->GetPoints() );
  outputMesh->SetPoints( points );
  points->Delete();
  
  polys = vtkCellArray::New();
  polys->DeepCopy( temporalMesh->GetPolys() );
  outputMesh->SetPolys( polys );
  polys->Delete();
  outputMesh->GetPointData()->DeepCopy(temporalMesh->GetPointData());
  outputMesh->GetFieldData()->PassData(temporalMesh->GetFieldData());
  outputMesh->BuildCells();
  outputMesh->BuildLinks();
  
  temporalMesh->DeleteLinks();
  temporalMesh->Delete();
  idPoints->Delete();

  return cuenta1;
}

/**--------------------------------------
Función que redefine las dos celdas que comparten el borde con mag. mayor a maximumDistance, e incluye dos celdas */

int insertNewCells(vtkPolyData *inputMesh, vtkIdType idPoint1, vtkIdType idPoint2, vtkIdType idNewPoint)
{
//   cout << "insertar celdas..." << endl;
  //  se determina cuales son las celdas que comparten el borde edge
  vtkIdList *idListCells1 = vtkIdList::New();
  vtkIdList *idListCells2 = vtkIdList::New();
  vtkIdList *idListCellsIntersection = vtkIdList::New();
  idListCellsIntersection->SetNumberOfIds(2);

//   // copy the input (only polys) to our working mesh
// //   vtkPolyData *temporalMesh = vtkPolyData::New();
//   vtkPoints *points = vtkPoints::New();
//   vtkCellArray *polys = vtkCellArray::New();
//   
//   points->DeepCopy( inputMesh->GetPoints() );
// //   temporalMesh->SetPoints( points );
//   polys->DeepCopy( inputMesh->GetPolys() );
// //   temporalMesh->SetPolys( polys );
// // 
// //   temporalMesh->BuildCells();
// //   temporalMesh->BuildLinks();
// //   temporalMesh->Update();

  inputMesh->GetPointCells( idPoint1, idListCells1);
  inputMesh->GetPointCells( idPoint2, idListCells2);

  vtkIdType mCells = idListCells1->GetNumberOfIds();
  vtkIdType nCells = idListCells2->GetNumberOfIds();

//   cout << "mCells= "<< mCells << " \tnCells= "<< nCells << endl;
  int j;
  int count=0;
  
//   cout << "Intersección: ";
  for(int i=0; i < mCells; i++)
  {
    j=0;
    while( idListCells1->GetId(i) != idListCells2->GetId(j) && j < nCells )
    {
      j++;
    }
    if( j < nCells )
    {
      idListCellsIntersection->SetId( count, idListCells1->GetId(i) );
//       cout << idListCells1->GetId(i) << ", ";
      count++;
      if( count >= 2 )
      {
        break;
      }
    }
  }
//   cout << endl;
  //  Después de identificar las dos celdas que comparten el borde, se copian estas como base para la inserción de las dos nuevas celdas
  vtkIdList *copyIdCellPoints1 = vtkIdList::New();
  vtkIdList *copyIdCellPoints2 = vtkIdList::New();

  inputMesh->GetCellPoints( idListCellsIntersection->GetId(0), copyIdCellPoints1 );
  inputMesh->GetCellPoints( idListCellsIntersection->GetId(1), copyIdCellPoints2 );

/*  cout << "lista de puntos de las celdas que coparten el borde:" << endl;
  cout << "cell1: " << copyIdCellPoints1->GetId(0) << ", " << copyIdCellPoints1->GetId(1) << ", " << copyIdCellPoints1->GetId(2) << endl;
  cout << "cell2: " << copyIdCellPoints2->GetId(0) << ", " << copyIdCellPoints2->GetId(1) << ", " << copyIdCellPoints2->GetId(2) << endl;
  cout << "idPoint1, idPoint2, idNewPoint: " << idPoint1 <<", "<< idPoint2 <<", "<< idNewPoint << endl;*/
  
  //  ahora se alistan los puntos que definen las dos celdas que deberán ser incluidas, con base en los puntos de las celdas que comparten el borde, reeplazando el punto idPoint1 por idNewPoint en cada celda
  vtkIdType idNewCell1, idNewCell2;
  vtkIdType npts = 3;
  vtkIdType ptsCell1[3], ptsCell2[3];

  vtkIdType tempIdPoint;
//   vtkIdType otherPoints[2];
  
  for(int i=0; i<3; i++)
  {
    tempIdPoint = copyIdCellPoints1->GetId(i);
   
    if( tempIdPoint == idPoint1 )
      ptsCell1[i] = idNewPoint;
    else
      ptsCell1[i] = tempIdPoint;
    
/*    if( tempIdPoint != idPoint1 && tempIdPoint != idPoint2 )
      otherPoints[0] = tempIdPoint;*/
    
    tempIdPoint = copyIdCellPoints2->GetId(i);
    
    if( tempIdPoint == idPoint1 )
      ptsCell2[i] = idNewPoint;
    else
      ptsCell2[i] = tempIdPoint;
    
/*    if( tempIdPoint != idPoint1 && tempIdPoint != idPoint2 )
      otherPoints[1] = tempIdPoint;*/
  }
  
  //  redefinición de las celdas que comparten el borde edge0
//   cout << "reemplazo de punto..."; 
  for(int i=0; i<2; i++){
    inputMesh->ReplaceCellPoint( idListCellsIntersection->GetId(i), idPoint2, idNewPoint );
  }
//   cout << " realizado." << endl;
  
/*  vtkIdList *listCellsOtherPoint = vtkIdList::New();
  
  for(int i=0; i<2; i++){
    inputMesh->GetPointCells(otherPoints[i], listCellsOtherPoint);
    inputMesh->ResizeCellList(otherPoints[i], listCellsOtherPoint->GetNumberOfIds() + 1);
  }
  listCellsOtherPoint->Delete();
  cout << "redimensión realizada." << endl;*/
  
  idNewCell1 = inputMesh->GetPolys()->InsertNextCell( npts, ptsCell1 );
  idNewCell2 = inputMesh->GetPolys()->InsertNextCell( npts, ptsCell2 );
//   cout << "celdas insertadas" << endl;
/*  inputMesh->Modified();
  inputMesh->Update();
  // copy the input (only polys) to our working mesh
  vtkPoints *points = vtkPoints::New();
  vtkCellArray *polys = vtkCellArray::New();
  
  points->DeepCopy( inputMesh->GetPoints() );
  polys->DeepCopy( inputMesh->GetPolys() );

  idNewCell1 = polys->InsertNextCell( npts, ptsCell1 );
  idNewCell2 = polys->InsertNextCell( npts, ptsCell2 );
  
  vtkIdType nocMesh = inputMesh->GetPolys()->GetNumberOfCells();
  inputMesh->Reset();
//   cout << "nocMesh: " << nocMesh << endl;
//   inputMesh->Allocate(nocMesh+2,nocMesh+2);

  inputMesh->SetPoints(points);
  inputMesh->SetPolys(polys);*/
  
  inputMesh->BuildCells();
  inputMesh->BuildLinks();
//   inputMesh->Update();

//   points->Delete();
//   polys->Delete();

  idListCells1->Delete();
  idListCells2->Delete();
  idListCellsIntersection->Delete();
  copyIdCellPoints1->Delete();
  copyIdCellPoints2->Delete();

  return 0;
}

