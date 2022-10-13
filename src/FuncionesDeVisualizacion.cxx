#include "vtkActor.h"
#include "vtkCellArray.h"
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"
#include "vtkCamera.h"
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

// Librerias de vtkINRIA3D
#include <vtkIsosurfaceManager.h>
#include <vtkViewImage2D.h>
#include <vtkViewImage3D.h>


void NormalsMeshVisualization( vtkPolyData *normalsMesh )
{
  vtkMaskPoints *ptMask = vtkMaskPoints::New();
// ptMask->SetOnRatio( 5 );
  ptMask->GenerateVerticesOn();
  ptMask->SingleVertexPerCellOn();
// ptMask->RandomModeOn();
  
  ptMask->SetInput( normalsMesh );
  
  vtkConeSource *cone = vtkConeSource::New();
  cone->SetResolution( 6 );
  
  vtkTransform *transform = vtkTransform::New();
  transform->Translate( 0.5, 0.0, 0.0 );
  
  vtkTransformPolyDataFilter *transformF = vtkTransformPolyDataFilter::New();
  transformF->SetTransform( transform );
  
  transformF->SetInput( cone->GetOutput() );
  
  vtkGlyph3D *glyph = vtkGlyph3D::New();
  glyph->SetInput( ptMask->GetOutput() );
  glyph->SetSourceConnection( transformF->GetOutputPort() );
  glyph->SetVectorModeToUseNormal();
  glyph->SetScaleModeToScaleByVector();
  glyph->SetScaleFactor( 3.0 );

//Pintado superficie desde aquí
  vtkPolyDataMapper *mapNormals = vtkPolyDataMapper::New();
  mapNormals->SetInput( normalsMesh );
    
  vtkPolyDataMapper *mapGlyph = vtkPolyDataMapper::New();
  mapGlyph->SetInputConnection(glyph->GetOutputPort());
// actor coordinates geometry, properties, transformation
  vtkActor *aNormals = vtkActor::New();
  aNormals->SetMapper(mapNormals);
  aNormals->GetProperty()->SetColor( 1.0, 0.49, 0.25 );
  aNormals->GetProperty()->SetOpacity( 1.0 );
  
  vtkActor *aGlyph = vtkActor::New();
  aGlyph->SetMapper(mapGlyph);
  aGlyph->GetProperty()->SetColor( 0.0, 0.0, 1.0 );
  
  vtkRenderer *ren = vtkRenderer::New();
  ren->AddActor(aNormals);
  ren->AddActor(aGlyph);
  ren->SetBackground( 0.0, 0.0, 0.0 );
  
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);
  
  renWin->SetSize(600,600);
  renWin->SetWindowName("Normales de la Superficie");

// an interactor
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  
  vtkInteractorStyleTrackballCamera *iCamera = vtkInteractorStyleTrackballCamera::New();
  iren->SetInteractorStyle( iCamera );
  
  renWin->Render();

// begin mouse interaction
  iren->Initialize();
  iren->Start();

  
  iren->Delete();
  renWin->Delete();
  iCamera->Delete();
  ren->Delete();
  aGlyph->Delete();
  aNormals->Delete();
  mapGlyph->Delete();
  mapNormals->Delete();
  glyph->Delete();
  transform->Delete();
  transformF->Delete();
  cone->Delete();
  ptMask->Delete();

  return;
}


void MeshVisualization( vtkPolyData *inputMesh ){

  vtkPolyData *mesh = vtkPolyData::New();
  vtkPoints *pointsMesh1 = vtkPoints::New();
  vtkCellArray *polysMesh1 = vtkCellArray::New();

  pointsMesh1->DeepCopy( inputMesh->GetPoints() );
  polysMesh1->DeepCopy( inputMesh->GetPolys() );
  
  mesh->Reset();
  mesh->SetPoints( pointsMesh1 );
  mesh->SetPolys( polysMesh1 );
  mesh->BuildCells();
  mesh->BuildLinks();
  
  pointsMesh1->Delete();
  polysMesh1->Delete();

//Pintado superficie desde aquí
  vtkPolyDataMapper *map = vtkPolyDataMapper::New();
  map->SetInput( mesh );
  
// actor coordinates geometry, properties, transformation
  vtkActor *actor = vtkActor::New();
  actor->SetMapper( map );
  actor->GetProperty()->SetColor( 1.0, 0.0, 0.0 );
  actor->GetProperty()->SetOpacity( 1.0 );
  
  
  vtkRenderer *ren = vtkRenderer::New();
  ren->AddActor( actor );
  ren->SetBackground( 0.0, 0.0, 0.0 );
  
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);

  renWin->SetSize(600,600);
  renWin->SetWindowName("Visualización de la Malla");

// an interactor
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  
  vtkInteractorStyleTrackballCamera *iCamera = vtkInteractorStyleTrackballCamera::New();
  iren->SetInteractorStyle( iCamera );
  
  renWin->Render();

// begin mouse interaction
  iren->Initialize();
  iren->Start();

  
  iren->Delete();
  renWin->Delete();
  iCamera->Delete();
  ren->Delete();
  actor->Delete();
  map->Delete();
  mesh->DeleteLinks();
  mesh->Delete();
  
  return;
}


void VisualizationMeshComparation( vtkPolyData *initialLiver, vtkPolyData *liverMesh ){

//Pintado superficie desde aquí
  vtkPolyDataMapper *map1 = vtkPolyDataMapper::New();
  map1->SetInput( initialLiver );
  
  vtkPolyDataMapper *map2 = vtkPolyDataMapper::New();
  map2->SetInput( liverMesh );
// actor coordinates geometry, properties, transformation
  vtkActor *actor1 = vtkActor::New();
  actor1->SetMapper( map1 );
  actor1->GetProperty()->SetColor( 1.0, 0.0, 0.0 );
  actor1->GetProperty()->SetOpacity( 0.8 );
  
  vtkActor *actor2 = vtkActor::New();
  actor2->SetMapper( map2 );
  actor2->GetProperty()->SetColor( 0.0, 0.0, 1.0 );
  actor2->GetProperty()->SetOpacity( 1.0 );
  
  vtkRenderer *ren = vtkRenderer::New();
  ren->AddActor( actor1 );
  ren->AddActor( actor2 );
  ren->SetBackground( 0.0, 0.0, 0.0 );
  
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  renWin->AddRenderer(ren);

  renWin->SetSize(600,600);
  renWin->SetWindowName("Comparación de las Superficies");

// an interactor
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  iren->SetRenderWindow(renWin);
  
  vtkInteractorStyleTrackballCamera *iCamera = vtkInteractorStyleTrackballCamera::New();
  iren->SetInteractorStyle( iCamera );
  
  renWin->Render();

// begin mouse interaction
  iren->Initialize();
  iren->Start();

  
  iren->Delete();
  renWin->Delete();
  iCamera->Delete();
  ren->Delete();
  actor1->Delete();
  actor2->Delete();
  map1->Delete();
  map2->Delete();

  return;
}


void ViewImages( vtkImageData *imagesCT, int range[2] )
{
//   cout << "\nFunción: SynchronizedViews()" << endl;
  /**
  Create 3 views, each of them will have a different orientation, .i.e.
  axial, sagittal and coronal.
  Create one 3D view.
   */
  vtkViewImage2D* view1 = vtkViewImage2D::New();
  vtkViewImage2D* view2 = vtkViewImage2D::New();
  vtkViewImage2D* view3 = vtkViewImage2D::New();
  vtkViewImage3D* view4 = vtkViewImage3D::New();

  vtkRenderWindowInteractor* iren1 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren2 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren3 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren4 = vtkRenderWindowInteractor::New();

  vtkRenderWindow* rwin1 = vtkRenderWindow::New();
  vtkRenderWindow* rwin2 = vtkRenderWindow::New();
  vtkRenderWindow* rwin3 = vtkRenderWindow::New();
  vtkRenderWindow* rwin4 = vtkRenderWindow::New();

  vtkRenderer* renderer1 = vtkRenderer::New();
  vtkRenderer* renderer2 = vtkRenderer::New();
  vtkRenderer* renderer3 = vtkRenderer::New();
  vtkRenderer* renderer4 = vtkRenderer::New();

  iren1->SetRenderWindow (rwin1);
  iren2->SetRenderWindow (rwin2);
  iren3->SetRenderWindow (rwin3);
  iren4->SetRenderWindow (rwin4);

  rwin1->AddRenderer (renderer1);
  rwin2->AddRenderer (renderer2);
  rwin3->AddRenderer (renderer3);
  rwin4->AddRenderer (renderer4);

  view1->SetRenderWindow ( rwin1 );
  view2->SetRenderWindow ( rwin2 );
  view3->SetRenderWindow ( rwin3 );
  view4->SetRenderWindow ( rwin4 );

  view1->SetRenderer ( renderer1 );
  view2->SetRenderer ( renderer2 );
  view3->SetRenderer ( renderer3 );
  view4->SetRenderer ( renderer4 );


  /**
  Set some properties to the views, like the interaction style, orientation and
  background color.
   */

  /*
  view1->SetInteractionStyle (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);
  view2->SetInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetInteractionStyle (vtkViewImage2D::ZOOM_INTERACTION);
  */

  // One can also associate to each button (left, middle, right and even wheel)
  // a specific interaction like this:
  
  view1->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view1->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view1->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view1->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);

  view2->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view2->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view2->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view2->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);

  view3->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view3->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);
  
  
  view1->SetLinkZoom (true);
  view2->SetLinkZoom (true);
  view3->SetLinkZoom (true);
  view4->SetLinkZoom (true);

  view1->SetOrientation (vtkViewImage2D::AXIAL_ID);
  view2->SetOrientation (vtkViewImage2D::CORONAL_ID);
  view3->SetOrientation (vtkViewImage2D::SAGITTAL_ID);

  view1->SetBackgroundColor (0.0,0.0,0.0);
  view2->SetBackgroundColor (0.0,0.0,0.0);
  view3->SetBackgroundColor (0.0,0.0,0.0);
  
  double color[3] = {0.0,0.0,0.0};
  view4->SetTextColor (color);
  view4->SetRenderingModeToPlanar();
  view4->SetCubeVisibility(1);
  

  view1->SetAboutData ("Powered by vtkINRIA3D");
  view2->SetAboutData ("Powered by vtkINRIA3D");
  view3->SetAboutData ("Powered by vtkINRIA3D");
  view4->SetAboutData ("Powered by vtkINRIA3D");

  /**
  Link the views together for synchronization.
   */
  view1->AddChild (view2);
  view2->AddChild (view3);
  view3->AddChild (view4);
  view4->AddChild (view1);


  view1->SetImage ( imagesCT );
  view2->SetImage ( imagesCT );
  view3->SetImage ( imagesCT );
  view4->SetImage ( imagesCT );


  /**
  Reset the window/level and the current position.
   */
  view1->SyncResetCurrentPoint();
  view1->SyncResetWindowLevel();


  int window, level;
    
  window = range[1] - range[0];
  level  = range[0] + window/2;

    
  view1->SetLevel( level );
  view1->SetWindow( window );
  view2->SetLevel( level );
  view2->SetWindow( window );
  view3->SetLevel( level );
  view3->SetWindow( window );
  view4->SetLevel( level );
  view4->SetWindow( window );
  
 
/*  vtkProperty *propertyPD = vtkProperty::New();
  propertyPD->SetColor( 0, 0, 0 );
  propertyPD->SetLineWidth( 3.0 );
  
  
  view1->AddDataSet( meshLiver, propertyPD);
  view2->AddDataSet( meshLiver, propertyPD);
  view3->AddDataSet( meshLiver, propertyPD);
  view4->AddDataSet( meshLiver, propertyPD);*/
  
  view1->Render();
  view2->Render();
  view3->Render();
  view4->Render();
  
  iren1->Start();
  iren2->Start();
  iren3->Start();
  iren4->Start();

  
  view1->Detach();
  view2->Detach();
  view3->Detach();
  view4->Detach();

  
  view1->Delete();
  view2->Delete();
  view3->Delete();
  view4->Delete();

  rwin1->Delete();
  rwin2->Delete();
  rwin3->Delete();
  rwin4->Delete();
  iren1->Delete();
  iren2->Delete();
  iren3->Delete();
  iren4->Delete();
  renderer1->Delete();
  renderer2->Delete();
  renderer3->Delete();
  renderer4->Delete();
  
  return;
}


void SynchronizedViews( vtkImageData *imagesCT, vtkPolyData *inputMesh, int range[2] )
{
  vtkPolyData *meshLiver = vtkPolyData::New();
  vtkPoints *pointsMesh1 = vtkPoints::New();
  vtkCellArray *polysMesh1 = vtkCellArray::New();

  pointsMesh1->DeepCopy( inputMesh->GetPoints() );
  polysMesh1->DeepCopy( inputMesh->GetPolys() );
  
  meshLiver->Reset();
  meshLiver->SetPoints( pointsMesh1 );
  meshLiver->SetPolys( polysMesh1 );
  meshLiver->BuildCells();
  meshLiver->BuildLinks();
  
  pointsMesh1->Delete();
  polysMesh1->Delete();

//   cout << "\nFunción: SynchronizedViews()" << endl;
  /**
  Create 3 views, each of them will have a different orientation, .i.e.
  axial, sagittal and coronal.
  Create one 3D view.
   */
  vtkViewImage2D* view1 = vtkViewImage2D::New();
  vtkViewImage2D* view2 = vtkViewImage2D::New();
  vtkViewImage2D* view3 = vtkViewImage2D::New();
  vtkViewImage3D* view4 = vtkViewImage3D::New();

  vtkRenderWindowInteractor* iren1 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren2 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren3 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren4 = vtkRenderWindowInteractor::New();

  vtkRenderWindow* rwin1 = vtkRenderWindow::New();
  vtkRenderWindow* rwin2 = vtkRenderWindow::New();
  vtkRenderWindow* rwin3 = vtkRenderWindow::New();
  vtkRenderWindow* rwin4 = vtkRenderWindow::New();

  vtkRenderer* renderer1 = vtkRenderer::New();
  vtkRenderer* renderer2 = vtkRenderer::New();
  vtkRenderer* renderer3 = vtkRenderer::New();
  vtkRenderer* renderer4 = vtkRenderer::New();

  iren1->SetRenderWindow (rwin1);
  iren2->SetRenderWindow (rwin2);
  iren3->SetRenderWindow (rwin3);
  iren4->SetRenderWindow (rwin4);

  rwin1->AddRenderer (renderer1);
  rwin2->AddRenderer (renderer2);
  rwin3->AddRenderer (renderer3);
  rwin4->AddRenderer (renderer4);

  view1->SetRenderWindow ( rwin1 );
  view2->SetRenderWindow ( rwin2 );
  view3->SetRenderWindow ( rwin3 );
  view4->SetRenderWindow ( rwin4 );

  view1->SetRenderer ( renderer1 );
  view2->SetRenderer ( renderer2 );
  view3->SetRenderer ( renderer3 );
  view4->SetRenderer ( renderer4 );


  /**
  Set some properties to the views, like the interaction style, orientation and
  background color.
   */

  /*
  view1->SetInteractionStyle (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);
  view2->SetInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetInteractionStyle (vtkViewImage2D::ZOOM_INTERACTION);
  */

  // One can also associate to each button (left, middle, right and even wheel)
  // a specific interaction like this:
  
  view1->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view1->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view1->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view1->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);

  view2->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view2->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view2->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view2->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);

  view3->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view3->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);
  
  
  view1->SetLinkZoom (true);
  view2->SetLinkZoom (true);
  view3->SetLinkZoom (true);
  view4->SetLinkZoom (true);

  view1->SetOrientation (vtkViewImage2D::AXIAL_ID);
  view2->SetOrientation (vtkViewImage2D::CORONAL_ID);
  view3->SetOrientation (vtkViewImage2D::SAGITTAL_ID);

  view1->SetBackgroundColor (0.0,0.0,0.0);
  view2->SetBackgroundColor (0.0,0.0,0.0);
  view3->SetBackgroundColor (0.0,0.0,0.0);
  
  double color[3] = {0.0,0.0,0.0};
  view4->SetTextColor (color);
  view4->SetRenderingModeToPlanar();
  view4->SetCubeVisibility(1);
  

  view1->SetAboutData ("Powered by vtkINRIA3D");
  view2->SetAboutData ("Powered by vtkINRIA3D");
  view3->SetAboutData ("Powered by vtkINRIA3D");
  view4->SetAboutData ("Powered by vtkINRIA3D");

  /**
  Link the views together for synchronization.
   */
  view1->AddChild (view2);
  view2->AddChild (view3);
  view3->AddChild (view4);
  view4->AddChild (view1);


  view1->SetImage ( imagesCT );
  view2->SetImage ( imagesCT );
  view3->SetImage ( imagesCT );
  view4->SetImage ( imagesCT );


  /**
  Reset the window/level and the current position.
   */
  view1->SyncResetCurrentPoint();
  view1->SyncResetWindowLevel();

  view1->Show2DAxis(false);
  view2->Show2DAxis(false);
  view3->Show2DAxis(false);

  view1->SetShowDirections(false);
  view2->SetShowDirections(false);
  view3->SetShowDirections(false);

  int window, level;
    
  window = range[1] - range[0];
  level  = range[0] + window/2;
//  window = 400;
//  level  = 70;

    
  view1->SetLevel( level );
  view1->SetWindow( window );
  view2->SetLevel( level );
  view2->SetWindow( window );
  view3->SetLevel( level );
  view3->SetWindow( window );
  view4->SetLevel( level );
  view4->SetWindow( window );
  
 
  vtkProperty *propertyPD = vtkProperty::New();
  propertyPD->SetColor( 1, 0, 0 );
  propertyPD->SetLineWidth( 2.0 );
  
  
  view1->AddDataSet( meshLiver, propertyPD);
  view2->AddDataSet( meshLiver, propertyPD);
  view3->AddDataSet( meshLiver, propertyPD);
  view4->AddDataSet( meshLiver, propertyPD);
  
  view1->Render();
  view2->Render();
  view3->Render();
  view4->Render();
  
  iren1->Start();
  iren2->Start();
  iren3->Start();

  propertyPD->SetOpacity( 0.7 );
  iren4->Start();
  


  view1->Detach();
  view2->Detach();
  view3->Detach();
  view4->Detach();

  propertyPD->Delete();

  view1->Delete();
  view2->Delete();
  view3->Delete();
  view4->Delete();

  rwin1->Delete();
  rwin2->Delete();
  rwin3->Delete();
  rwin4->Delete();
  iren1->Delete();
  iren2->Delete();
  iren3->Delete();
  iren4->Delete();
  renderer1->Delete();
  renderer2->Delete();
  renderer3->Delete();
  renderer4->Delete();
  
  meshLiver->DeleteLinks();
  meshLiver->Delete();
  
  return;
}

/**Aquí se visualizan las curvas de las segmentaciones de referencias y las obtenidas por el método propuesto sobre las imágenes de TAC*/
void SynchronizedViews( vtkImageData *imagesCT, vtkPolyData *inputRefMesh, vtkPolyData *inputSegMesh )
{
//   cout << "\nFunción: SynchronizedViews()" << endl;
  /**
  Create 3 views, each of them will have a different orientation, .i.e.
  axial, sagittal and coronal.
  Create one 3D view.
   */
  vtkViewImage2D* view1 = vtkViewImage2D::New();
  vtkViewImage2D* view2 = vtkViewImage2D::New();
  vtkViewImage2D* view3 = vtkViewImage2D::New();
  vtkViewImage3D* view4 = vtkViewImage3D::New();

  vtkRenderWindowInteractor* iren1 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren2 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren3 = vtkRenderWindowInteractor::New();
  vtkRenderWindowInteractor* iren4 = vtkRenderWindowInteractor::New();

  vtkRenderWindow* rwin1 = vtkRenderWindow::New();
  vtkRenderWindow* rwin2 = vtkRenderWindow::New();
  vtkRenderWindow* rwin3 = vtkRenderWindow::New();
  vtkRenderWindow* rwin4 = vtkRenderWindow::New();

  vtkRenderer* renderer1 = vtkRenderer::New();
  vtkRenderer* renderer2 = vtkRenderer::New();
  vtkRenderer* renderer3 = vtkRenderer::New();
  vtkRenderer* renderer4 = vtkRenderer::New();

  iren1->SetRenderWindow (rwin1);
  iren2->SetRenderWindow (rwin2);
  iren3->SetRenderWindow (rwin3);
  iren4->SetRenderWindow (rwin4);

  rwin1->AddRenderer (renderer1);
  rwin2->AddRenderer (renderer2);
  rwin3->AddRenderer (renderer3);
  rwin4->AddRenderer (renderer4);

  view1->SetRenderWindow ( rwin1 );
  view2->SetRenderWindow ( rwin2 );
  view3->SetRenderWindow ( rwin3 );
  view4->SetRenderWindow ( rwin4 );

  view1->SetRenderer ( renderer1 );
  view2->SetRenderer ( renderer2 );
  view3->SetRenderer ( renderer3 );
  view4->SetRenderer ( renderer4 );


  /**
  Set some properties to the views, like the interaction style, orientation and
  background color.
   */

  /*
  view1->SetInteractionStyle (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);
  view2->SetInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetInteractionStyle (vtkViewImage2D::ZOOM_INTERACTION);
  */

  // One can also associate to each button (left, middle, right and even wheel)
  // a specific interaction like this:
  
  view1->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view1->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view1->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view1->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);

  view2->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view2->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view2->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view2->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);

  view3->SetLeftButtonInteractionStyle   (vtkViewImage2D::ZOOM_INTERACTION);
  view3->SetMiddleButtonInteractionStyle (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetWheelInteractionStyle        (vtkViewImage2D::SELECT_INTERACTION);
  view3->SetRightButtonInteractionStyle  (vtkViewImage2D::WINDOW_LEVEL_INTERACTION);
  
  
  view1->SetLinkZoom (true);
  view2->SetLinkZoom (true);
  view3->SetLinkZoom (true);
  view4->SetLinkZoom (true);

  view1->SetOrientation (vtkViewImage2D::AXIAL_ID);
  view2->SetOrientation (vtkViewImage2D::CORONAL_ID);
  view3->SetOrientation (vtkViewImage2D::SAGITTAL_ID);

  view1->SetBackgroundColor (0.0,0.0,0.0);
  view2->SetBackgroundColor (0.0,0.0,0.0);
  view3->SetBackgroundColor (0.0,0.0,0.0);
  
  double color[3] = {0.0,0.0,0.0};
  view4->SetTextColor (color);
  view4->SetRenderingModeToPlanar();
  view4->SetCubeVisibility(1);
  

  view1->SetAboutData ("Powered by vtkINRIA3D");
  view2->SetAboutData ("Powered by vtkINRIA3D");
  view3->SetAboutData ("Powered by vtkINRIA3D");
  view4->SetAboutData ("Powered by vtkINRIA3D");

  /**
  Link the views together for synchronization.
   */
  view1->AddChild (view2);
  view2->AddChild (view3);
  view3->AddChild (view4);
  view4->AddChild (view1);


  view1->SetImage ( imagesCT );
  view2->SetImage ( imagesCT );
  view3->SetImage ( imagesCT );
  view4->SetImage ( imagesCT );


  /**
  Reset the window/level and the current position.
   */
  view1->SyncResetCurrentPoint();
  view1->SyncResetWindowLevel();


  int window, level;
    
//   window = range[1] - range[0];
//   level  = range[0] + window/2;
  window = 400;
  level  = 70;
  
    
  view1->SetLevel( level );
  view1->SetWindow( window );
  view2->SetLevel( level );
  view2->SetWindow( window );
  view3->SetLevel( level );
  view3->SetWindow( window );
  view4->SetLevel( level );
  view4->SetWindow( window );
  
 
  vtkProperty *propertyRef = vtkProperty::New();
  propertyRef->SetColor( 1, 0, 0 );
  propertyRef->SetLineWidth( 3.0 );
  
  view1->AddDataSet( inputRefMesh, propertyRef);
  view2->AddDataSet( inputRefMesh, propertyRef);
  view3->AddDataSet( inputRefMesh, propertyRef);
  view4->AddDataSet( inputRefMesh, propertyRef);
  
  vtkProperty *propertySeg = vtkProperty::New();
  propertySeg->SetColor( 0, 0, 1 );
  propertySeg->SetLineWidth( 3.0 );
  
  view1->AddDataSet( inputSegMesh, propertySeg);
  view2->AddDataSet( inputSegMesh, propertySeg);
  view3->AddDataSet( inputSegMesh, propertySeg);
  view4->AddDataSet( inputSegMesh, propertySeg);

  view1->Render();
  view2->Render();
  view3->Render();
  view4->Render();
  
  iren1->Start();
  iren2->Start();
  iren3->Start();

//   propertyRef->SetOpacity( 0.7 );
//   iren4->Start();

  view1->Detach();
  view2->Detach();
  view3->Detach();
  view4->Detach();

  propertyRef->Delete();
  propertySeg->Delete();
  
  view1->Delete();
  view2->Delete();
  view3->Delete();
  view4->Delete();

  rwin1->Delete();
  rwin2->Delete();
  rwin3->Delete();
  rwin4->Delete();
  iren1->Delete();
  iren2->Delete();
  iren3->Delete();
  iren4->Delete();
  renderer1->Delete();
  renderer2->Delete();
  renderer3->Delete();
  renderer4->Delete();
  
  return;
}
