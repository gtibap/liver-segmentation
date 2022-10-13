// Aquí se registrará una región (obtenida por una segmentacion inicial realizada por análisis del histograma y operaciones morfológicas) con el hígado promedio.
// Autor: Gerardo Tibamoso Pedraza. 
//     Claro con las ayudas y tutoriales correspondientes.

#include "vtkActor.h"
#include "vtkAxes.h"
#include "vtkBoxWidget.h"
#include "vtkCamera.h"
#include "vtkCommand.h"
#include "vtkCoordinate.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkMatrix4x4.h"
#include "vtkPolyData.h"
#include "vtkPolyDataMapper.h"
#include "vtkProperty.h"
#include "vtkRenderer.h"
#include "vtkRenderWindow.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkSTLReader.h"
#include "vtkSTLWriter.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkTubeFilter.h"

#include <stdio.h>
#include <iostream>
using namespace std;

// Callback for the interaction
class vtkMyCallback : public vtkCommand
{
  public:
    static vtkMyCallback *New() 
    { 
      return new vtkMyCallback; 
    }
    virtual void Execute(vtkObject *caller, unsigned long, void*)
    {
      vtkTransform *t = vtkTransform::New();
      vtkBoxWidget *widget = reinterpret_cast<vtkBoxWidget*>(caller);
      widget->GetTransform(t);
      widget->GetProp3D()->SetUserTransform(t);
      
      vtkMatrix4x4 *matriz = vtkMatrix4x4::New();
      t->GetMatrix( matriz );
      cuenta=0;
      for(int i=0; i<4; i++){
        for(int j=0; j<4; j++){
          paso[cuenta] = matriz->GetElement(i,j);
          ++cuenta;
        }
      }
      t->Delete();
      matriz->Delete();
    }
    
    void GetTransformMatrix( double tMatrix[16] )
    {
        for(int i=0; i<16; i++){
          tMatrix[i] = paso[i];
        }
    }
  
  private:
    int cuenta, test;
    double scale[3], orientation[3], position[3];
    static double paso[16]; 
};
double vtkMyCallback::paso[]={1,0,0,0, 0,1,0,0, 0,0,1,0, 0,0,0,1};

/**---------------------------------------
------------------------------------------
Esta función es llamada por ConstruccionSuperficie.cxx */

int AjusteManual( vtkPolyData *meshRegion, vtkPolyData *meshLiver, double translate1[3], double translate2[3], double matrixTransform[16] )
{
  
  vtkPolyDataMapper *mapper1 = vtkPolyDataMapper::New();
  vtkPolyDataMapper *mapper2 = vtkPolyDataMapper::New();
  vtkActor *actor1 = vtkActor::New();
  vtkActor *actor2 = vtkActor::New();
  vtkProperty *property1 = vtkProperty::New();
  vtkProperty *property2 = vtkProperty::New();
  vtkRenderer *ren1 = vtkRenderer::New();
  vtkRenderWindow *renWin = vtkRenderWindow::New();
  vtkRenderWindowInteractor *iren = vtkRenderWindowInteractor::New();
  vtkInteractorStyleTrackballCamera *iCamera = vtkInteractorStyleTrackballCamera::New();
  
  property1->SetColor( 0.85, 0.33, 0.15 );
  property1->SetOpacity( 1.0 );
  property2->SetColor( 0.15, 0.33, 0.85 );
  property2->SetOpacity( 1.0 );
//  property2->SetOpacity( 0.75 );
	//property2->SetDiffuseColor( 0.15, 0.33, 0.85 );
	//property2->SetAmbientColor( 0.15, 0.33, 0.85 );
	//property2->SetLighting(false);
	//property2->LightingOff();
	//property2->ShadingOn();
	//property2->FrontfaceCullingOff();
  
  actor1->SetProperty( property1 );
  actor2->SetProperty( property2 );
	

  ren1->SetBackground( 0.0 , 0.0 , 0.0 );

  renWin->SetSize( 600 , 600 );

  mapper1->SetInput( meshLiver ); 
  mapper2->SetInput( meshRegion );
  actor1->SetMapper( mapper1 );
  actor2->SetMapper( mapper2 );
  ren1->AddActor( actor1 );
  ren1->AddActor( actor2 );
  renWin->AddRenderer( ren1 );
  iren->SetRenderWindow( renWin );
  iren->SetInteractorStyle( iCamera );
  
  double bounds[6], positionActor[3], originActor[3];
  
  actor1->GetBounds( bounds );
  actor1->GetOrigin( originActor );
  actor1->GetPosition( positionActor );
  
/*  cout << "bounds actor 1: " << bounds[0] <<", " << bounds[1] <<", " << bounds[2] <<", " << bounds[3] <<", " << bounds[4] <<", " << bounds[5] << endl;
  cout << "originActor 1: " << originActor[0] << ", " << originActor[1] << ", " << originActor[2] << endl;
  cout << "positionActor 1: " << positionActor[0] << ", " << positionActor[1] << ", " << positionActor[2] << endl;*/
  
  actor2->GetBounds( bounds );
  actor2->GetOrigin( originActor );
  actor2->GetPosition( positionActor );

/*  cout << "bounds actor 2: " << bounds[0] <<", " << bounds[1] <<", " << bounds[2] <<", " << bounds[3] <<", " << bounds[4] <<", " << bounds[5] << endl;
  cout << "originActor 2: " << originActor[0] << ", " << originActor[1] << ", " << originActor[2] << endl;
  cout << "positionActor 2: " << positionActor[0] << ", " << positionActor[1] << ", " << positionActor[2] << endl;*/
  
// ajustar el centro del actor2 (en este caso el hígado promedio) sobre el centro del actor1 (que en este caso es la región donde está el hígado)  
//   double translate1[3], translate2[3];
  
  for( int i=0; i<3; i++ ){
//     translateAtlas[i] = -actor1->GetCenter()[i] + actor2->GetCenter()[i];
    translate1[i] = -actor1->GetCenter()[i];
    translate2[i] = -actor2->GetCenter()[i];
  }
//   actor1->SetPosition( translateAtlas );
//   double origen[3] = {0,0,0};
  actor1->SetPosition(translate1);
  actor2->SetPosition(translate2);
  
  // Here we use a vtkBoxWidget to transform the underlying coneActor (by
  // manipulating its transformation matrix). Many other types of widgets
  // are available for use, see the documentation for more details.
  //
  // The SetInteractor method is how 3D widgets are associated with the render
  // window interactor. Internally, SetInteractor sets up a bunch of callbacks
  // using the Command/Observer mechanism (AddObserver()). The place factor 
  // controls the initial size of the widget with respect to the bounding box
  // of the input to the widget.
  vtkBoxWidget *boxWidget = vtkBoxWidget::New();
  boxWidget->SetInteractor( iren );
  boxWidget->SetPlaceFactor( 1.5 );

  //
  // Place the interactor initially. The input to a 3D widget is used to 
  // initially position and scale the widget. The EndInteractionEvent is
  // observed which invokes the SelectPolygons callback.
  //
  boxWidget->SetProp3D( actor1 );
  boxWidget->PlaceWidget();
  vtkMyCallback *callback = vtkMyCallback::New();
  boxWidget->AddObserver( vtkCommand::InteractionEvent, callback );

  //
  // Normally the user presses the "i" key to bring a 3D widget to life. Here
  // we will manually enable it so it appears with the cone. 
  //
  boxWidget->On();
  
  
/*  vtkAxes *axes = vtkAxes::New();
  axes->SetOrigin( actor1->GetOrigin() );
  axes->SetScaleFactor( 50.0 );
  vtkTubeFilter *axesTubes = vtkTubeFilter::New();
  axesTubes->SetInputConnection(axes->GetOutputPort());
  axesTubes->SetRadius( 5.0 );
  axesTubes->SetNumberOfSides(6);
  vtkPolyDataMapper *axesMapper = vtkPolyDataMapper::New();
  axesMapper->SetInputConnection(axesTubes->GetOutputPort());
  vtkActor *axesActor = vtkActor::New();*/
  
//   axesActor->SetMapper(axesMapper);
//   ren1->AddActor( axesActor );
  
  
  iren->Initialize();
  iren->Start();
  
  //  se obtiene la matriz de tranformación
  callback->GetTransformMatrix( matrixTransform );
  
 
  boxWidget->Delete();
  callback->Delete();
  property1->Delete();
  property2->Delete();
//   axes->Delete();
//   axesTubes->Delete();
//   axesMapper->Delete();
//   axesActor->Delete();
  mapper1->Delete();
  mapper2->Delete();
  ren1->Delete();
  actor1->Delete();
  actor2->Delete();
  renWin->Delete();
  iCamera->Delete();
  iren->Delete();
  
  return 0;
}

