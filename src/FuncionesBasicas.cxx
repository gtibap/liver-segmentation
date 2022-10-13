/*#include "vtkSphereSource.h"
#include "vtkPolyDataMapper.h"
#include "vtkActor.h"
#include "vtkRenderWindow.h"
#include "vtkRenderer.h"
#include "vtkRenderWindowInteractor.h"
#include "vtkInteractorStyleTrackballCamera.h"
#include "vtkPointSet.h"
#include "vtkPolyData.h"
#include "vtkPointData.h"
#include "vtkIdList.h"
#include "vtkDataArray.h"
#include "vtkFloatArray.h"
#include "vtkDoubleArray.h"
#include "vtkPolyDataReader.h"
#include "vtkPolyDataNormals.h"
#include "vtkProgrammableFilter.h"

#include "vtkProperty.h"
#include "vtkTransform.h"
#include "vtkTransformPolyDataFilter.h"
#include "vtkMaskPoints.h"
#include "vtkConeSource.h"
#include "vtkGlyph3D.h"

#include "vtkTriangleFilter.h"
#include "vtkDecimatePro.h"
#include "vtkCamera.h"*/


#include <stdlib.h>
#include <float.h>
#include <math.h>
#include <string.h>
// #include "vtkPolyDataWriter.h"
// #include "vtkSTLReader.h"
#include <iostream>
using namespace std;


// double minimal (double a, double b);
// void vplus (double A[3], double B[3], double C[3]);
// void vminus (double A[3], double B[3], double C[3]);
// void divide (double A[3], double s, double C[3]);
// double vtimes (double A[3], double B[3]);
// void vcopy (double A[3], double B[3]);
// double norm (double A[3]);
// void times(double s, double A[3], double C[3]);
// double distance (double A[3], double B[3]);
// double oppositeAngle (double edgeAB, double edgeBC, double edgeCA);
// double cot (double angle);
// int trianglesArea (vtkPolyData *input, double* &areas);
// double patchArea (vtkPolyData *input, vtkIdType i);
// void findNeighbors (vtkPolyData *input, vtkIdType center, vtkIdList *&nbs);

// void findNeighbors1 (vtkPolyData *input, vtkIdList **nbs);
// void findCentroids(vtkPolyData *input, vtkIdList **lnbs, vtkDoubleArray *lcentroids );
// void tensionForce( vtkPolyData *input, vtkDoubleArray *lcentroids, vtkDoubleArray *lTForce );
/*void rigidityForce( vtkPolyData *input, vtkDoubleArray lTForce, vtkIdList **lnbs, vtkDoubleArray *lRForce );*/
// void ballonForce( vtkPolyData *input, vtkDoubleArray *lBForce );

// void localForce (vtkPolyData *input, vtkIdType center, vtkDoubleArray *lsF, vtkIdList *nbs);
// void globalForce (vtkPolyData *input, int k, vtkFloatArray *norm, vtkIdList **lnbs, double newpoint[3]);
// double radialForce (vtkPolyData *input, vtkIdType xi, double gcen[3], double radius, double newRpoint[3]);
// void expansionForce (vtkPolyData *input, vtkIdType i, double gcen[3], double radius, double Rext, double vH[3], double n[3], vtkDoubleArray *reF);
// void curvature (vtkPolyData *input, vtkIdType i, double vH[3], vtkDoubleArray *vK);
// void histogram (double *factorJ, int sizeJ, char *filename);
// 
// 
// vtkProgrammableFilter *pointf = vtkProgrammableFilter::New();
// int iter;
// 



double minimal (double a, double b) {
    if (a < b) 
        return a;
    else
        return b;
}

// vector operations
void vplus (double A[3], double B[3], double C[3]) {
    C[0] = A[0] + B[0];
    C[1] = A[1] + B[1];
    C[2] = A[2] + B[2];
    return;
}

void vminus (double A[3], double B[3], double C[3]) {
    C[0] = A[0] - B[0];
    C[1] = A[1] - B[1];
    C[2] = A[2] - B[2];
    return;
}

void divide (double A[3], double s, double C[3]) {
  
  if(s != 0.0){
    C[0] = A[0] / s;
    C[1] = A[1] / s;
    C[2] = A[2] / s;
  }
  else{
    C[0] = 0.0;
    C[1] = 0.0;
    C[2] = 0.0;
  }
    
    return;
}

double vtimes (double A[3], double B[3]) {
    return (A[0]*B[0] + A[1]*B[1] + A[2]*B[2]);
}

void vcopy (double A[3], double B[3]) {
    B[0] = A[0];
    B[1] = A[1];
    B[2] = B[2];
    return;
}

double norm (double A[3]) {
    return sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
}

void times(double s, double A[3], double C[3]) {
    C[0] = s * A[0];
    C[1] = s * A[1];
    C[2] = s * A[2];
    return;
}

double distance (double A[3], double B[3]) {
    return sqrt((B[0]-A[0])*(B[0]-A[0])+(B[1]-A[1])*(B[1]-A[1])+(B[2]-A[2])*(B[2]-A[2]));
}

double oppositeAngle (double edgeAB, double edgeBC, double edgeCA) {
    if (edgeBC != 0.0 && edgeCA != 0.0)
        return acos((edgeCA*edgeCA+edgeBC*edgeBC-edgeAB*edgeAB)/(2.0*edgeBC*edgeCA));
    else
        return 0.0;
}

double cot (double angle) {
    if (angle != 0.0)
        return 1.0/tan(angle);
    else
        return DBL_MAX;
}
// -----------------

// function used to find the triangles areas in the surface
// takes the list of cells in the polydata and calculates the area for each cell
// params: input polydata (surface), array of doubles (store each triangle area)
// returns: an integer with the number of cells in the polydata
// int trianglesArea (vtkPolyData *input, double* &areas) {
//     vtkIdList *points = vtkIdList::New(); //list to store the cell's points
// 	
// 	// get number of cells in the polydata
//     vtkIdType noc = input->GetNumberOfCells();
//     vtkIdType nop, point;
//     vtkIdType maxp = 0;
//     int i, j;
//     double A[3], B[3], C[3], plist[3][3], temp;
// 	
// 	// initialize array of triangles areas
//     areas = new double[noc];
// 	
// 	// for each cell in the polydata
//     for (i=0; i<noc; i++) {
// 		// get the points that define the cell
//         input->GetCellPoints(i, points);
//         nop = points->GetNumberOfIds();
// 		
// /*        if(nop > maxp)
//             maxp=nop;*/
//         
// 		// for each point in the cell
//         for (j=0; j<nop; j++) {
// 			// store point's coordinates
//             point = points->GetId(j);
//             input->GetPoint(point,plist[j]);
//         }
// 		
// 		// find the distances between each pair of points (lenght of each edge)
//         vminus(plist[2],plist[0],A);
//         vminus(plist[2],plist[1],B);
//         vminus(plist[1],plist[0],C);
// 		// calculates: |C|^2 * |A|^2 - (C·A)^2
//         temp = ((norm(C)*norm(C))*(norm(A)*norm(A)))-(vtimes(C,A)*vtimes(C,A));
// 		
// 		// if (temp == 0), then the triangle has area = 0
//         if (temp > 0.0) {
//             //comprobé que con esta formula se obtiene el área pero hay que dividirla en 2
//             areas[i] = sqrt(temp);
//         } else {
//             areas[i] = 0.0;
//         }
//     }
// //     cout << "número máx. de puntos de las celdas = " << maxp << endl;	
//     return (int)noc;
// }
// 
// // given a point ID, this function calculates the patch area around it
// // this patch area is the sum of the areas of each triangle that contains the given point
// // params: input polydata (surface), integer point ID, boolean
// // returns: value of patch area
// double patchArea (vtkPolyData *input, vtkIdType i) {
//     vtkIdList *cells = vtkIdList::New();  //list to store cells
//     vtkIdList *points = vtkIdList::New(); //list to store the cell's points
//     double pArea = 0.0;
//     double x[3], A[3], B[3], C[3], plist[3][3], temp;
//     int j, k;
// 	
// 	// obtain the point's coordinates
//     input->GetPoint(i, x);
// 	
// 	//obtain the cell list
//     input->GetPointCells(i, cells);
// 	
//     vtkIdType noc = cells->GetNumberOfIds();
//     vtkIdType nop, point;
// 	
//     if (noc < 3) {
//         cerr<<"Point "<<i<<" ("<<x[0]<<","<<x[1]<<","<<x[2]<<") has less than 3 triangles adjacents"<<endl;
//     }
// 	
// 	// for each cell in the patch
//     for (j=0; j<noc; j++) {
// 		// get the point that define a cell
//         input->GetCellPoints(j, points);
//         nop = points->GetNumberOfIds();
// 		
// 		// for each point in the cell
//         for (k=0; k<nop; k++) {
// 			// get the point's coordinates
//             point = points->GetId(k);
//             input->GetPoint(point,plist[k]);
//         }
// 		
// 		// find the distances between each pair of points (lenght of each edge)
//         vminus(plist[2],plist[0],A);
//         vminus(plist[2],plist[1],B);
//         vminus(plist[1],plist[0],C);
// 		// calculates: |C^2| * |A^2| - (C·A)^2
//         temp = ((norm(C)*norm(C))*(norm(A)*norm(A)))-(vtimes(C,A)*vtimes(C,A));
// 		
// 		// if (temp == 0), then the triangle has area = 0
//         if (temp > 0.0) {
//             pArea += sqrt(temp);
//         } else {
//             pArea += 0.0;
//         }
//     }
// 	
//     cells->Delete();
//     points->Delete();
// 	
//     return pArea;
// }
// 
// // procedure to obtain the list of neighbors of a point, called 'center'.
// // use the cells that share the point 'center', and the points that defines each cell
// void findNeighbors (vtkPolyData *input, vtkIdType center, vtkIdList *&nbs) {
//     vtkIdList *cells = vtkIdList::New();  //list to store cells
//     vtkIdList *points = vtkIdList::New(); //list to store the cells points
// 	
// 	//obtain the cell list
//     input->GetPointCells(center, cells);
// 	
//     vtkIdType noc = cells->GetNumberOfIds();
// 	
//     vtkIdType cell, nop, point, dummy;
//     int i, j;
// 	
//     for (i=0; i<noc; i++) {
//         cell = cells->GetId(i); 
// 		// for each cell, obtain the point list
//         input->GetCellPoints(cell, points);
// 		
//         nop = points->GetNumberOfIds();
// 		
// 		// for each point, try to insert it in the neighbors list
//         for (j=0; j<nop; j++) {
//             point = points->GetId(j); 
//             if (point != center && nbs->IsId(point)==(-1)) { 
//                 dummy = nbs->InsertNextId(point);
//             }
//         }
//     }
// 	
//     cells->Delete();
//     points->Delete();
//     return;
// }

// void findNeighbors1 (vtkPolyData *input, vtkIdList **nbs) {
//   
//   vtkIdList *cells = vtkIdList::New();  //list to store cells
//   vtkIdList *points = vtkIdList::New(); //list to store the cells points
//   vtkIdType nop, nopc, noc, cell, point;
//   	
//   nop = input->GetNumberOfPoints();
//   //find neighbors for every vertex 
//   for(int i=0; i<nop; i++){
//     nbs[i] = vtkIdList::New();
//     //obtain the cell list
//     input->GetPointCells(i, cells);
//     vtkIdType noc = cells->GetNumberOfIds();
//     // for each cell, obtain the point list
//     for(int j=0; j<noc; j++){
//       cell = cells->GetId(j);
//       input->GetCellPoints(cell, points);
//       nopc = points->GetNumberOfIds();
//       // for each point, try to insert it in the neighbors list
//       for(int k=0; k<nopc; k++){
//         point = points->GetId(k);
//         if (point != i && nbs[i]->IsId(point)==(-1)) { 
//           nbs[i]->InsertNextId(point);
//         }
//       }
//     }
//   }
// 	
//   cells->Delete();
//   points->Delete();
// 
//   return;
// }
// 
// void findCentroids( vtkPolyData *input, vtkIdList **lnbs, vtkDoubleArray *lcentroids ) {
//   
//   vtkIdType nop, non, nb;
//   double x[3], centroid[3];
//   double suma[3];
//   
//   nop = input->GetNumberOfPoints();
//   for(int i=0; i<nop; i++){
//     non = lnbs[i]->GetNumberOfIds();
// //     cout << "non= " << non << endl;
//     suma[0] = suma[1] = suma[2] = 0;
//     for(int j=0; j<non; j++){
//       nb = lnbs[i]->GetId(j);
// //       cout << "nb= " << nb << endl;
//       input->GetPoint(nb, x);
// //       cout << "x= " << x[0] <<", "<< x[1] <<", "<< x[2] << endl;
//       vplus(suma, x, suma);
// //       cout << "suma= " << suma[0] <<", "<< suma[1] <<", "<< suma[2] << endl;
//     }
//     divide(suma, non, centroid);
// //     cout << i << " centroid= " << centroid[0] <<", "<< centroid[1] <<", "<< centroid[2] << endl;
//     lcentroids->InsertTupleValue(i, centroid);
//   }
//   
//   return;
// }
// 
// void tensionForce( vtkPolyData *input, vtkDoubleArray *lcentroids, vtkDoubleArray *lTForce ){
//   
//   vtkIdType nop;
//   double centroid[3], x[3], force[3];
//   
//   nop = input->GetNumberOfPoints();
//   for(int i=0; i<nop; i++){
//     lcentroids->GetTupleValue(i, centroid);
//     input->GetPoint(i, x);
//     vminus (centroid, x, force);
//     lTForce->InsertTupleValue(i, force);
//   }
// 
// return;
// }

// void rigidityForce( vtkPolyData *input, vtkDoubleArray *lTForce, vtkIdList **lnbs, vtkDoubleArray *lRForce ){
//   
//   vtkIdType nop, non, nb;
//   double ntforce[3], suma[3], vtforce[3], force[3];
//   
//   nop = input->GetNumberOfPoints();
//   for(int i=0; i<nop; i++){
//     non = lnbs[i]->GetNumberOfIds();
//     suma[0] = suma[1] = suma[2] = 0.0;
//     for(int j=0; j<non; j++){
//       nb = lnbs[i]->GetId(j);
//       lTForce->GetTupleValue(nb, ntforce);
//       vplus(suma, ntforce, suma);
//     }
//     divide(suma, non, ntforce);
//     lTForce->GetTupleValue(i, vtforce);
//     vminus(vtforce, ntforce, force);
//     lRForce->InsertTupleValue(i, force);
//   }
// 
//   return;
// }

// void ballonForce( vtkPolyData *input, vtkDoubleArray *lBForce ){
//   
//   vtkPolyDataNormals *normalsPDFilter = vtkPolyDataNormals::New();
//   
//   normalsPDFilter->SetInput( input );
//   normalsPDFilter->Update();
//   
//   lBForce = (vtkDoubleArray *)( normalsPDFilter->GetOutput()->GetPointData()->GetNormals() );
//   
// //   normalsPDFilter->Delete();
// 
//   return;
// }

//procedure to calculate:  \frac{1}{N_k} \sum_{j \in N_k} (x_j - x_k)
//the local smoothing force that moves each point x_k in the direction of the centroid of its neighbors
// void localForce (vtkPolyData *input, vtkIdType center, vtkDoubleArray *lsF, vtkIdList *nbs) {
//     vtkIdType non;
// 	
//     double x[3], cen[3], T1[3], newpoint[3];
//     newpoint[0] = 0.0;
//     newpoint[1] = 0.0;
//     newpoint[2] = 0.0;
//     input->GetPoint(center, cen);
// 	
//     non = nbs->GetNumberOfIds();
// 	
// 	// sum all distances between center and neighbors
//     for (int i=0; i<non; i++) {
//         input->GetPoint(nbs->GetId(i), x);
//         vminus(x, cen, T1);
//         vplus(newpoint, T1, newpoint);
//     }
// 	//determinate the average distance
//     divide(newpoint, non, newpoint);
//     lsF->InsertTuple(center, newpoint);
//     return;
// }

// procedure to calculate: \frac{1}{N_p} \sum_i^V \sum_{j \in N_i} (n_i \cdot(x_j - x_i)) * n_i
// the global smoothing force
// void globalForce (vtkPolyData *input, int k, vtkFloatArray *norm, vtkIdList **lnbs, double newpoint[3]) {
//     int i, j;
//     double T1[3], T2[3], x[3], y[3];
//     double *n = new double[3];
//     vtkIdType nonb;
//     vtkIdType nop = input->GetNumberOfPoints();
//     newpoint[0] = 0.0;
//     newpoint[1] = 0.0;
//     newpoint[2] = 0.0;
// 	
//     for (i=0; i<nop; i++) {
// 		//get the point with id i, the "center", and store it in x
//         input->GetPoint(i, x);
// 			
//         if (k==0) {
//             lnbs[i] = vtkIdList::New();
//         }
// 		//find the list of neighbors of x
//         findNeighbors(input, i, lnbs[i]);
//         nonb = lnbs[i]->GetNumberOfIds();
// 		//get the normal in x
//         n = norm->GetTuple(i);
// 		//cout<<"n = ["<<n[0]<<","<<n[1]<<","<<n[2]<<"]"<<endl;
// 		
//         for (j=0; j<nonb; j++) {
// 			//get one of the neigbors, and store it in y
//             input->GetPoint(lnbs[i]->GetId(j), y);
// 			//calculate:  \sum_i^V \sum_{j \in N_i} (n_i \cdot(x_j - x_i)) * n_i
//             vminus(y, x, T1); // x_j - x_i
// 			//vminus(x, y, T1);
//             times(vtimes(n, T1), n, T2); // n_i  (x_j - x_i) * n_i
//             vplus(newpoint, T2, newpoint); // \sum (n_i  (x_j - x_i)) * n_i
// 			//cout<<"newpoint = ["<<newpoint[0]<<","<<newpoint[1]<<","<<newpoint[2]<<"]"<<endl;
//         }
//     }
// 	//divide the sum by the number of points in the surface
//     divide(newpoint, nop, newpoint);
// 	
// 	//cout<<"delete n 1"<<endl;
// 	//delete [] n;
//     return;
// }

//procedure to calculate: R_k - x_k
//the radial force that drives each point x_k toward the surface of a sphere
//R_k is the radial projection of x_k onto the sphere
// return the norm of the radial force vector
// double radialForce (vtkPolyData *input, vtkIdType xi, double gcen[3], double radius, double newRpoint[3]) {
//     double t, x[3], r[3];
// 	
//     newRpoint[0] = 0.0;
//     newRpoint[1] = 0.0;
//     newRpoint[2] = 0.0;
//     input->GetPoint(xi, x);
// 	
// 	// get distance between the point and the geometric center
//     vminus(x, gcen, r);
// 	
// 	// find the relation between the sphere radius and the distance r
//     t = ((radius / norm(r)) - 1.0);
// 	
// 	// multiply this factor to the distance
//     times(t, r, newRpoint);
// 	
// 	// return the norm of the radial force vector
//     return norm(newRpoint);
// }
// 
// void expansionForce (vtkPolyData *input, vtkIdType i, double gcen[3], double radius, double Rext, double vH[3], double n[3], vtkDoubleArray *reF) {
//     double x[3], y[3], coef, newEpoint[3];
//     vtkIdType nop = input->GetNumberOfPoints();
// 	
//     newEpoint[0] = 0.0;
//     newEpoint[1] = 0.0;
//     newEpoint[2] = 0.0;
// 	
// 	//get the point with id i and store it in x 
//     input->GetPoint(i, x);
// 	
// // 	if (i == 549) {
// // 		cout<<"radialForce = "<<radialForce(input, i, gcen, radius, y)<<endl;
// // 		cout<<"norm(vH) = "<<norm(vH)<<endl;
// // 		cout<<"vtimes(x, n) = "<<vtimes(x, n)<<endl;
// // 		cout<<"Rext - vtimes(x, n) = "<<(Rext - vtimes(x, n))<<endl;
// // 	}
//     coef = radialForce(input, i, gcen, radius, y) - norm(vH) + (Rext - vtimes(x, n));
// 	//coef = (Rext - vtimes(x, n)) - norm(vH);
// // 	if (i == 549) {
// // 		cout<<"coef = "<<coef<<endl<<endl;
// // 		cout<<"n = ["<<n[0]<<","<<n[1]<<","<<n[2]<<"]"<<endl;
// // 	}
//     times(coef, n, newEpoint);
// // 	if (i == 549) cout<<"eF_"<<i<<" = ["<<newEpoint[0]<<","<<newEpoint[1]<<","<<newEpoint[2]<<"]"<<endl;
//     reF->InsertTuple(i, newEpoint);
//     return;
// }
// 
// void curvature (vtkPolyData *input, vtkIdType i, double vH[3], vtkDoubleArray *vK) {
//     vtkIdList *cells = vtkIdList::New();  //list to store cells
//     vtkIdList *points = vtkIdList::New(); //list to store the cells points
// 	
//     double x[3], y[2][3], vKt[3], vN[3], vM[3], A[3], B[3], C[3], T1[3], T2[3], T3[3], T4[3], neig[3], cosAlpha;
//     int j,k,l;
// 	
//     vKt[0] = vKt[1] = vKt[2] = 0.0;
//     vH[0] = vH[1] = vH[2] = 0.0;
// 	
// 	// get the point's coordinates
//     input->GetPoint(i, x);
// 	
// 	//obtain the cell list
//     input->GetPointCells(i, cells);
// 	
//     vtkIdType noc = cells->GetNumberOfIds();
//     vtkIdType cell, nop, point;
// 	
// 	// for each cell in the patch
//     for (j=0; j<noc; j++) {
//         cell = cells->GetId(j); 
// 		// obtain the point list
//         input->GetCellPoints(cell, points);
// 		
//         nop = points->GetNumberOfIds();
// 		
//         l=0;
// 		// for each point in the cell
//         for (k=0; k<nop; k++) {
// 			// obtain the point's coordinates
//             point = points->GetId(k);
//             input->GetPoint(point,neig);
// 			
//             if (point != i && (x[0]!=neig[0] || x[1]!=neig[1] || x[2]!=neig[2])) { 
//                 input->GetPoint(point,y[l]);
//                 l++;
//             }
//         }
// 		
//         vminus(x, y[0], A); // x_i - x_j
//         vminus(y[1], y[0], B); // x_(j+1) - x_j
//         vminus(y[0], y[1], C); // x_j - x_(j+1)
//         divide(A, norm(A), vN); // n = (x_i - x_j) / |x_i - x_j|
//         divide(B, norm(B), vM); // m = (x_(j+1) - x_j) / |x_(j+1) - x_j|
//         cosAlpha = vtimes(vN, vM); // cos(alpha) = n · m
// 		
//         times(cosAlpha, vM, T1); // cos(alpha) * m
//         vminus(vN, T1, T2); // n - cos(alpha) * m
//         divide(T2, (2.0 * sqrt(1 - (cosAlpha * cosAlpha))), T3); // (n - cos(alpha) * m) / (2 * sin(alpha))
//         times(norm(C), T3, T4); // |x_j - x_(j+1)| * {(n - cos(alpha) * m) / (2 * sin(alpha))}
//         vplus(vKt, T4, vKt); // k = k + [|x_j - x_(j+1)| * {(n - cos(alpha) * m) / (2 * sin(alpha))}]
// 		//cout<<"T4 = ["<<T4[0]<<","<<T4[1]<<","<<T4[2]<<"]"<<endl;
//     }
// 	
// 	//cout<<"vK_"<<i<<" = ["<<vKt[0]<<","<<vKt[1]<<","<<vKt[2]<<"]"<<endl;
//     times(-1.0, vKt, vH);
// 	//cout<<"vH_"<<i<<" = ["<<vH[0]<<","<<vH[1]<<","<<vH[2]<<"]"<<endl;
//     vK->InsertTuple(i, vKt);
// 	
//     cells->Delete();
//     points->Delete();
//     return;
// }
// 
// // procedure to generate an histogram for factor j (variation between initial and actual area)
// // writes an file with the histogram
// // params: array of factors j, size of array
// void histogram (double *factorJ, int sizeJ, char *filename) {
//     int i, *hist, n;
//     double minJ, maxJ, x;
//     minJ = DBL_MAX;
//     maxJ = DBL_MIN;
//     int size = 50;
//     char str[25];
// 	
//     ofstream out;
// 	
//     strcpy(str, "histogram_");
//     strcat(str, filename);
//     strcat(str, ".dat");
//     out.open(str, ios::out);
// 	
// // 	for (i=0; i<sizeJ; i++) {
// // 		if (factorJ[i] < minJ) minJ = factorJ[i];
// // 		if (factorJ[i] > maxJ) maxJ = factorJ[i];
// // 	}
// 	
//     minJ = 0.0;
//     maxJ = 2.0;
// 	
//     hist = new int[size];
//     for (i=0; i<size; i++) {
//         hist[i] = 0;
//     }
// 	
//     for (i=0; i<sizeJ; i++) {
//         n = int( size * (factorJ[i] - minJ) / (maxJ - minJ) );
//         n=(n<size)?n:size;
//         hist[n]++;
//     }
// 	
//     for (i=0; i<size; i++) {
// 		//x=minJ+(i+0.5)*(maxJ-minJ)/size;
//         x=minJ+(i)*(maxJ-minJ)/size;
//         out<<i<<"\t"<<x<<"\t"<<hist[i]<<endl;
//     }
// 	
//     delete [] hist;
//     return;
// }


// void mapping (void * arg) {

//   vtkPolyData *input = pointf->GetPolyDataInput();
// 
//   vtkIdType nop = input->GetNumberOfPoints();
//   vtkIdType noc = input->GetNumberOfCells();
//   cout << "# Points = " << nop << endl;    
//   cout << "# Cells = " << noc << endl;
// 
//   vtkPoints *newpts = vtkPoints::New();
//   vtkPoints *oldpts = vtkPoints::New();
//   double x[3], y[3], min[3], max[3], r[3], gf[3], gcen[3], radius, T1[3], T2[3], T3[3], T4[3], Rext, lambda, vH[3], lI[1], xOld[3], patchAreasI[nop], patchAreasF[nop], prom[3], vKt[3], reFt[3], lambdaI_avg, predX[3], dxt[3], dx_1t[3], dx_2t[3], dpredX[3];
//   double *n = new double[3];
//   double sum1, sum2, patchA, dS;
//   double areaTotalI = 0.0;
//   double totalPatchAreasI = 0.0;
//   double totalPatchAreasF = 0.0;
//   double *areasI, *areasT;
//   int i, j, k, l, asize;
//   double p[3], q[3], s[3];
//   vtkIdType pId, qId, rId, sId;
//   double angleR, angleS, edgePQ, edgeQR, edgeRP, edgeQS, edgeSP, Dpp, areaPQR, areaPQS, temp;
//   double Hvector[nop][3], Hvalue[nop], meanHvalue, deltaMI, deltaMF;
//   vtkIdType non;
//   int NumNeighbs[nop];
//   vtkIdType nocellsnb, nopointsnb, pqrT;
//     
// 	
//   input->BuildLinks();
//   input->Update();
// 
//   vtkIdList **lnbs = new vtkIdList*[nop];
//   vtkDoubleArray *lcentroids = vtkDoubleArray::New();
//   vtkDoubleArray *lTForce = vtkDoubleArray::New();
//   vtkDoubleArray *lRForce = vtkDoubleArray::New();
//   vtkDoubleArray *lNormals = vtkDoubleArray::New();
//   vtkDoubleArray *lBForce = vtkDoubleArray::New();
//     
//   lcentroids->SetNumberOfComponents(3);
//   lTForce->SetNumberOfComponents(3);
//   lRForce->SetNumberOfComponents(3);
//   lNormals->SetNumberOfComponents(3);
//   lBForce->SetNumberOfComponents(3);
// 
//   findNeighbors1( input, lnbs );
//   findCentroids( input, lnbs, lcentroids );
//     
//   tensionForce( input, lcentroids, lTForce );
//   rigidityForce( input, lTForce, lnbs, lRForce );
// //     ballonForce( input, lBForce );
//     
//   vtkPolyDataNormals *normalsPDFilter = vtkPolyDataNormals::New();
//   
//   normalsPDFilter->SetInput( input );
//   normalsPDFilter->SplittingOff();
// //     normalsPDFilter->AutoOrientNormalsOn();
//     
// /*    cout << "computePointNormals= " << normalsPDFilter->GetComputePointNormals() << endl;
//   cout << "computeCellNormals= " << normalsPDFilter->GetComputeCellNormals() << endl;
//   cout << "splitting= " << normalsPDFilter->GetSplitting() << endl;
//   cout << "featureAngle= " << normalsPDFilter->GetFeatureAngle() << endl;
//   cout << "consistency= " << normalsPDFilter->GetConsistency() << endl;
//   cout << "autoOrientNormals= " << normalsPDFilter->GetAutoOrientNormals() << endl;
//   cout << "flipNormals= " << normalsPDFilter->GetFlipNormals() << endl;
//   cout << "NonManifoldTraversal= " << normalsPDFilter->GetNonManifoldTraversal() << endl;*/
//         
//   normalsPDFilter->Update();
//   
//     
//     
// //     lBForce = (vtkDoubleArray *)( normalsPDFilter->GetOutput()->GetPointData()->GetNormals() );
//   lNormals = (vtkDoubleArray *)( normalsPDFilter->GetOutput()->GetPointData()->GetNormals() );
//     
// //     double *scalar = new double[3];
//   double tuple[3];
//   int numberOfComponents;
//   vtkIdType numberOfTuples;
//     
//   numberOfComponents = lNormals->GetNumberOfComponents();
//   numberOfTuples = lNormals->GetNumberOfTuples();
//   cout << "components= " << numberOfComponents << "  tuples= " << numberOfTuples << endl;    
// //     se calculan los vectores normales unitarios numberOfTuples
//     
//   int inicio= iter, final= inicio + 10;
//   for( i=inicio; i <final ; i++ ){
//     cout << "paso " << i << endl;
//     lNormals->GetTupleValue( i, tuple );
// //       cout << "paso 2" << endl;
// //       cout << "tuple= " << tuple[0] << ", " << tuple[1] << ", " << tuple[2] << endl;
//     divide( tuple, norm( tuple ), tuple );
// //       cout << "tuple= " << tuple[0] << ", " << tuple[1] << ", " << tuple[2] << endl;
// //       cout << "paso 3" << endl;
//     lNormals->SetTupleValue( i, tuple );
//   }
// 
//   numberOfComponents = lNormals->GetNumberOfComponents();
//   numberOfTuples = lNormals->GetNumberOfTuples();
// 
//   for( i=inicio; i<final; i++ ){
//     lNormals->GetTupleValue( i, tuple );
//     cout << "tuple= " << tuple[0] << ", " << tuple[1] << ", " << tuple[2] << endl;
//     cout << "mag= " << sqrt( tuple[0]*tuple[0] + tuple[1]*tuple[1] + tuple[2]*tuple[2] ) << endl;
//   }
// // 

//     delete[] scalar;
    
//     vtkIdType nonb;    
//     double rf[3], tf[3], point[3], centroid[3];
//     for(i=0; i<10; i++){
//     lRForce->GetTupleValue(i, rf);
//     cout <<"rf= " << rf[0] <<", "<< rf[1] <<", "<< rf[2] << endl;
// }

    /*    double tf[3], point[3], centroid[3];
  for(i=100; i<110; i++){
  input->GetPoint(i, point);
  lcentroids->GetTupleValue(i, centroid);
  lTForce->GetTupleValue(i, tf);
  cout <<"point= " << point[0] <<", "<< point[1] <<", "<< point[2] << endl;
  cout <<"centroid= " << centroid[0] <<", "<< centroid[1] <<", "<< centroid[2] << endl;
  cout <<"tf= " << tf[0] <<", "<< tf[1] <<", "<< tf[2] << endl;
  cout << endl;
}*/
    
    
//     double *p1 = new double[3];
/*    double p1[3];
  for(i=0; i<10; i++){
  lcentroids->GetTupleValue(i, p1);
  cout <<"p1= " << p1[0] <<", "<< p1[1] <<", "<< p1[2] << endl;
}*/
    
/*    for(i=0; i<10; i++){
  cout << "nonb= " << lnbs[i]->GetNumberOfIds() << endl;
}*/
/*    for(i=0; i<10; i++){
  lnbs[i] = vtkIdList::New();
  findNeighbors1(input, i, lnbs[i]);
  nonb = lnbs[i]->GetNumberOfIds();
  cout << "nonb= " << nonb << endl;
}*/
    
    
//     lnbs->Delete();
//     lBForce->Delete();
/*  lRForce->Delete();
  lTForce->Delete();
  lcentroids->Delete();
  delete[] lnbs;*/
    
    /*
  double lambdaS = 1.2;
  double dt = 0.001;
  double beta = 1.0;
	
    
  vtkDoubleArray *dk_2 = vtkDoubleArray::New(); // derivative of x in time k-2
  dk_2->SetNumberOfComponents(3);
  vtkDoubleArray *dk_1 = vtkDoubleArray::New(); // derivative of x in time k-1
  dk_1->SetNumberOfComponents(3);
  vtkDoubleArray *dk = vtkDoubleArray::New(); // derivative of x in time k
  dk->SetNumberOfComponents(3);
  vtkDoubleArray *dvK = vtkDoubleArray::New(); // value of vK in time - 1
  dvK->SetNumberOfComponents(3);
	
  double zero[3];
  zero[0] = zero[1] = zero[2] = 0.0;
  for (i=0; i<nop; i++) {
  dk_2->InsertTuple(i, zero);
  dk_1->InsertTuple(i, zero);
  dk->InsertTuple(i, zero);
  dvK->InsertTuple(i, zero);
}
	
	//initialize auxiliary vectors
	
  for (i=0; i<nop; i++) {
  for (j=0; j<3; j++) {
  Hvector[i][j] = 0.0;
}
  NumNeighbs[i] = 0;
  Hvalue[i] = 0.0;
}
	
  vtkDoubleArray *vK = vtkDoubleArray::New(); // K vector
  vK->SetNumberOfComponents(3);
	
	//-----------------------------------------------------------------------------
	// Cálculo de \beta estimado (parte inicial)
  meanHvalue = 0.0;
  for (i=0; i<nop; i++) {
		
		// calculate: curvature
  curvature(input, i, vH, vK);
  vK->GetTuple(i, vKt);
		
			// calculate mean curvature vector in P
  Hvector[i][0] = - vKt[0] / (2.0 * patchAreasI[i]);
  Hvector[i][1] = - vKt[1] / (2.0 * patchAreasI[i]);
  Hvector[i][2] = - vKt[2] / (2.0 * patchAreasI[i]);
			// calculate mean curvature value in P
  Hvalue[i] = norm(Hvector[i]);
  cerr<<"Hvalue["<<i<<"]: "<<Hvalue[i]<<endl;
		
  meanHvalue += (Hvalue[i] * patchAreasI[i]);
}
  meanHvalue /= totalPatchAreasI;
  vK->Delete();
	
  deltaMI = 0.0;
  for (i=0; i<nop; i++) {
  deltaMI += ((Hvalue[i] - meanHvalue) * (Hvalue[i] - meanHvalue) * patchAreasI[i]);
}
  deltaMI /= totalPatchAreasI;
	
	//-----------------------------------------------------------------------------
	
  for (i=0; i<nop; i++) {
  for (j=0; j<3; j++) {
  Hvector[i][j] = 0.0;
}
  NumNeighbs[i] = 0;
  Hvalue[i] = 0.0;
}
	
  for (k=0; k<iter; k++) {
		// determine R_ext (extension reference value)
  Rext = 0.0;
  for (i=0; i<nop; i++) {
  input->GetPoint(i, x);
		//cout<<"distance(x_"<<i<<", cen) = "<<distance(x, gcen)<<endl;
  if (distance(x, gcen) > Rext) Rext = distance(x, gcen);
}
		//determine the sphere radius
  radius = Rext;
		//cout<<"R_ext = "<<Rext<<endl;
		
		//------------- Calculate the global force -------------
		//obtain the surface normals
  vtkPolyDataNormals *pdnorm = vtkPolyDataNormals::New();
  pdnorm->SetInput(input);
  pdnorm->Update();
		
  vtkFloatArray *norm = vtkFloatArray::New();
  norm->SetNumberOfComponents(3);
  norm = (vtkFloatArray *)(pdnorm->GetOutput()->GetPointData()->GetNormals());
		
  vtkDoubleArray *lsF = vtkDoubleArray::New(); // local smoothing force
  lsF->SetNumberOfComponents(3);
  vtkDoubleArray *vK = vtkDoubleArray::New(); // K vector
  vK->SetNumberOfComponents(3);
  vtkDoubleArray *reF = vtkDoubleArray::New(); // radial expansion Force
  reF->SetNumberOfComponents(3);
  vtkDoubleArray *lambdaI = vtkDoubleArray::New(); // local preservation value
  lambdaI->SetNumberOfComponents(1);
		
		//-- First step: smoothing forces
  gf[0] = gf[1] = gf[2] = 0.0;
  globalForce(input, k, norm, lnbs, gf);
		//cout<<"global force: ["<<gf[0]<<","<<gf[1]<<","<<gf[2]<<"]"<<endl;
		
  for (i=0; i<nop; i++) {
			//get the point with id i, and store it in x
  input->GetPoint(i, x);
  if (k==0) {
  oldpts->InsertPoint(i, x[0], x[1], x[2]);
} else {
  oldpts->SetPoint(i, x[0], x[1], x[2]);
}
			//cout<<"x_"<<i<<"(t) = ["<<x[0]<<","<<x[1]<<","<<x[2]<<"]"<<endl;
			//get the normal in x
  n = norm->GetTuple(i);
			//calculate local smoothing force:  \frac{1}{N_k} \sum_{j \in N_k} (x_j - x_k), and store it in lsF
  localForce(input, i, lsF, lnbs[i]);
			
			//update coordinates x^~(t) = x(t) + [F_s + F_sg] * dt
  vplus(lsF->GetTuple(i), gf, T1);
  times(lambdaS, T1, T1);
  times(dt, T1, T2);
			//cout<<"after first step, smoothing force = ["<<T2[0]<<","<<T2[1]<<","<<T2[2]<<"]"<<endl;
  vplus(x, T2, x);
			
  newpts->SetPoint(i, x[0], x[1], x[2]);
// 			oldpts->GetPoint(i, xOld);
// 			vminus(x, xOld, T1);
// 			//if (T1[0]>0.1 || T1[0]<-0.1 || T1[1]>0.1 || T1[1]<-0.1 || T1[2]>0.1 || T1[2]<-0.1)
// 				cout<<"(x - xOld)_"<<i<<" = ["<<T1[0]<<","<<T1[1]<<","<<T1[2]<<"]"<<endl;
}
		
  input->SetPoints(newpts);
		
		//-- Second step: radial expansion
		
  sum1 = sum2 = 0.0;
// 		cout<<endl<<endl<<"Iteración "<<k+1<<endl;
  for (i=0; i<nop; i++) {
			//get the normal in x
  n = norm->GetTuple(i);
			// calculate: curvature
  curvature(input, i, vH, vK);
			// calculate: local expansion force
  expansionForce(input, i, gcen, radius, Rext, vH, n, reF);
			// calculate: local preservation parameter
  vplus(lsF->GetTuple(i), gf, T1);
  vK->GetTuple(i, vKt);
  reF->GetTuple(i, reFt);
  sum1 += vtimes(vKt, T1);
  sum2 += vtimes(vKt, reFt);
// 			if (i==549) {
// 				cout<<"vtimes(vK, T1) = "<<vtimes(vKt, T1)<<endl;
// 				cout<<"vK_"<<i<<" = ["<<vKt[0]<<","<<vKt[1]<<","<<vKt[2]<<"]"<<endl;
// 				cout<<"reF_"<<i<<" = ["<<reFt[0]<<","<<reFt[1]<<","<<reFt[2]<<"]"<<endl;
// 				cout<<"vtimes(vK, reF) = "<<vtimes(vKt, reFt)<<endl;
// 			}
  if (vtimes(vKt, reFt)<0.001 && vtimes(vKt, reFt)>-0.001) {
  lI[0] = vtimes(dvK->GetTuple(i), T1) / vtimes(dvK->GetTuple(i), reFt);
} else {
  lI[0] = vtimes(vKt, T1) / vtimes(vKt, reFt);
}
  dvK->SetTuple(i, vKt);
  lambdaI->InsertTuple(i, lI);
}
		// calculate: global preservation parameter
  lambda = sum1 / sum2;
		
  prom[0] = prom[1] = prom[2] = 0.0;
  for (i=0; i<nop; i++) {
			// get the point with id i, and store it in x
  input->GetPoint(i, x);
			//get the normal in x
  n = norm->GetTuple(i);
			// calculate lambdaI_avg
// 			non = lnbs[i]->GetNumberOfIds();
// 			lambdaI_avg = 0.0;
// 			for (j=0; j<non; j++) {
// 				l = lnbs[i]->GetId(j);
// 				lambdaI->GetTuple(l, lI);
// 				lambdaI_avg += lI[0];
// 			}
// 			lambdaI_avg /= non;
  lambdaI->GetTuple(i, lI);
			// update coordinates
  times((((1.0 + beta) * lI[0]) - (beta * lambda)), reF->GetTuple(i), T1);
  times(dt, T1, T2);
			//if (T2[0]>0.5 || T2[0]<-0.5 || T2[1]>0.5 || T2[1]<-0.5 || T2[2]>0.5 || T2[2]<-0.5)
			//	cout<<"mov_"<<i<<"_"<<k<<" = ["<<T2[0]<<","<<T2[1]<<","<<T2[2]<<"]"<<endl;
  vplus(prom, T2, prom);
			
			//---- predictor-corrector
  dk_2->SetTuple(i, dk_1->GetTuple(i));
  dk_2->GetTuple(i, dx_2t);

  dk_1->SetTuple(i, dk->GetTuple(i));
  dk_1->GetTuple(i, dx_1t);
  dk->SetTuple(i, T1);
  dk->GetTuple(i, dxt);
  times(23, dxt, dxt); // 23 x'_n
  times(16, dx_1t, dx_1t); // 16 x'_{n-1}
  times(5, dx_2t, dx_2t); // 5 x'_{n-2}
  vminus(dxt, dx_1t, T3); // 23x'_n - 16 x'_{n-1}
  vplus(T3, dx_2t, T3); // 23x'_n - 16 x'_{n-1} + 5 x'_{n-2}
  times((dt / 12), T3, T3); // (h / 12) (23x'_n - 16 x'_{n-1} + 5 x'_{n-2})
  vplus(x, T3, predX); // x_{n+1} = x_n + (h / 12) (23x'_n - 16 x'_{n-1} + 5 x'_{n-2})
			
  j = 0;
  do {
  vminus(predX, x, T4);
  divide(T4, dt, dpredX);
				
  dk_1->GetTuple(i, dx_1t);
  dk->GetTuple(i, dxt);
  times(5, dpredX, dpredX); // 5 x'_{n+1}
  times(8, dxt, dxt); // 8 x'_n
  vplus(dpredX, dxt, T3); // 5 x'_{n+1} + 8 x'_n
  vminus(T3, dx_1t, T3); // 5 x'_{n+1} + 8 x'_n - x'_{n-1}
  times((dt / 12), T3, T3); // (h / 12) (5 x'_{n+1} + 8 x'_n - x'_{n-1})
  vplus(x, T3, predX);
				
  j++;
} while (j<1);
			
  vcopy(predX, x);
			
			//vplus(x, T2, x);
			//cout<<"x_"<<i<<"(t+1) = ["<<x[0]<<","<<x[1]<<","<<x[2]<<"]"<<endl;
			
			//i==1359 || i==549 || i==45 || i==118 || i==1862 || i==38 || i==122 || i==1999
// 			if (i==549) {
// 				cout<<endl<<"x_"<<i<<" = ["<<x[0]<<","<<x[1]<<","<<x[2]<<"]"<<endl;
// 				cout<<"n_"<<i<<" = ["<<n[0]<<","<<n[1]<<","<<n[2]<<"]"<<endl;
  // 				
// 				cout<<"lambdaI_avg_"<<i<<" = "<<lambdaI_avg<<endl;
// 				cout<<"lambda ="<<lambda<<endl;
// 				cout<<"((1.0 + beta) * lI[0]) - (beta * lambda) = "<<((1.0 + beta) * lambdaI_avg) - (beta * lambda)<<endl;
  // 				
// 				non = lnbs[i]->GetNumberOfIds();
// 				for (j=0; j<non; j++) {
// 					l = lnbs[i]->GetId(j);
// 					input->GetPoint(l, y);
// 					//get the normal in y
// 					n = norm->GetTuple(l);
// 					cout<<"x_"<<l<<" = ["<<y[0]<<","<<y[1]<<","<<y[2]<<"]"<<endl;
// 					cout<<"n_"<<l<<" = ["<<n[0]<<","<<n[1]<<","<<n[2]<<"]"<<endl;
// 				}
// 			}
			
  newpts->SetPoint(i, x[0], x[1], x[2]);
// 			oldpts->GetPoint(i, xOld);
// 			vminus(x, xOld, T1);
// 			cout<<"(x - xOld)_"<<i<<" = ["<<T1[0]<<","<<T1[1]<<","<<T1[2]<<"]"<<endl;
}
  input->SetPoints(newpts);
		
  divide(prom, nop, prom);
		//cout<<"prom-vel_"<<k<<" = ["<<prom[0]<<","<<prom[1]<<","<<prom[2]<<"]"<<endl;
		
  dS = 0.0;
  for (i=0; i<nop; i++) {
  patchA = patchArea (input, i);
  dS += ((patchA - patchAreasI[i]) * (patchA - patchAreasI[i])) / (patchAreasI[i] * patchAreasI[i]);
}
  dS = sqrt(dS) / nop;
		
  norm->Delete();
  lsF->Delete();
  vK->Delete();
  reF->Delete();
  lambdaI->Delete();
		//pdnorm->Delete();
}
// 	cout<<"dS = "<<dS<<endl;
// 	cout<<"numero de puntos: "<<nop<<endl;
	
  asize = trianglesArea (input, areasT);
	
  double areaTotal = 0.0;
  double factorJt_avg = 0.0;
  double factorJp_avg = 0.0;
  double factorJt[asize], factorJp[nop];
	
  for (i=0; i<asize; i++) {
  factorJt[i] = areasI[i] / areasT[i];
  factorJt_avg += factorJt[i];
  areaTotal += areasT[i];
		//cout<<"factorJ_"<<i<<" = "<<factorJ[i]<<endl;
}
  factorJt_avg = factorJt_avg / asize;
  for (i=0; i<asize; i++) {
  factorJt[i] = factorJt[i] / factorJt_avg;
}
  histogram (factorJt, asize, "triangles");
	
  cout<<"Area total inicial: "<<areaTotalI<<" unidades de area"<<endl;
  cout<<"Area total despues de "<<iter<<" iteraciones: "<<areaTotal<<" unidades de area"<<endl;
  cout<<"Cambio en el area total :"<<(areaTotalI - areaTotal)<<" unidades de area"<<endl;
// 	cout<<"Numero de triangulos en la superficie: "<<asize<<endl;
// 	cout<<"Factor de area promedio para los triangulos: "<<factorJt_avg<<endl;
	
  for (i=0; i<nop; i++) {
  patchAreasF[i] = patchArea (input, i);
  factorJp[i] = patchAreasI[i] / patchAreasF[i];
  factorJp_avg += factorJp[i];
  totalPatchAreasF += patchAreasF[i];
}
  factorJp_avg = factorJp_avg / nop;
  for (i=0; i<asize; i++) {
  factorJp[i] = factorJp[i] / factorJp_avg;
}
  histogram (factorJp, nop, "patches");
  cout<<"Numero de puntos en la superficie: "<<nop<<endl;
  cout<<"Factor de area promedio para los parches: "<<factorJp_avg<<endl;
	
	//-----------------------------------------------------------------------------
	// Cálculo de \beta estimado (parte final)
  for (i=0; i<nop; i++) {
  for (j=0; j<3; j++) {
  Hvector[i][j] = 0.0;
}
  NumNeighbs[i] = 0;
  Hvalue[i] = 0.0;
}
	
  vtkDoubleArray *vK2 = vtkDoubleArray::New(); // K vector
  vK2->SetNumberOfComponents(3);
	
  meanHvalue = 0.0;
  for (i=0; i<nop; i++) {
		
		// calculate: curvature
  curvature(input, i, vH, vK2);
  vK2->GetTuple(i, vKt);
		
			// calculate mean curvature vector in P
  Hvector[i][0] = - vKt[0] / (2.0 * patchAreasF[i]);
  Hvector[i][1] = - vKt[1] / (2.0 * patchAreasF[i]);
  Hvector[i][2] = - vKt[2] / (2.0 * patchAreasF[i]);
			// calculate mean curvature value in P
  Hvalue[i] = norm(Hvector[i]);
		
  meanHvalue += (Hvalue[i] * patchAreasF[i]);
}
  meanHvalue /= totalPatchAreasF;
  vK2->Delete();
	
  deltaMF = 0.0;
  for (i=0; i<nop; i++) {
  deltaMF += ((Hvalue[i] - meanHvalue) * (Hvalue[i] - meanHvalue) * patchAreasF[i]);
}
  deltaMF /= totalPatchAreasF;
	
  double betaEst;
  betaEst = sqrt(2 * (deltaMI - deltaMF));
  cout<<"beta estimado: "<<betaEst<<endl;
	
	//-----------------------------------------------------------------------------
	
  pointf->GetPolyDataOutput()->CopyStructure(input);
  pointf->GetPolyDataOutput()->SetPoints(newpts);
  newpts->Delete();
  oldpts->Delete();
  dk_2->Delete();
  dk_1->Delete();
  dk->Delete();
  dvK->Delete();
	//delete[] n;
  //lnbs->Delete();*/
//Código comentado hasta aquí    
//   return;
// } 
