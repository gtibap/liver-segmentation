//La fusión de vértices se realiza para aquellos cuyo borde se menor a un valor establecido.
// La estrategia consiste en encontrar las celdas a borrar, cambiar los puntos de las celdas que comparten los puntos de las celdas a borrar, y solo cuando se haya recorrido toda la malla, las celdas y los puntos seleccionados serán borrados.
 
// Autor: Gerardo Tibamoso Pedraza. 

// librerias VTK
#include "vtkDataArray.h"
#include "vtkDoubleArray.h"


// librerias básicas de C++
#include <stdio.h>
#include <string.h>
#include <iostream>
using namespace std;


int stopEvaluation( vtkDoubleArray *errorDistance, FILE *fpt, int lCounts[3][3], int medianlCounts[3] )
{
  
  int cuenta, temp;
  float factor, error;
  int countPoints[10] = {}; //  Se inicializa el arreglo countPoints[] en cero
  int copylCounts[3][3];
  
  vtkIdType nop = errorDistance->GetNumberOfTuples();
//   cout << "nop = " << nop << endl;
  fprintf( fpt, "%d\t", nop );
  
  for( int iPoint=0; iPoint < nop; iPoint++)
  {
    error = errorDistance->GetValue(iPoint);
    cuenta = 0;
    factor = 0.1;
    while( error > factor && factor < 1.0)
    {
      factor+=0.1;
      cuenta++;
    }
    countPoints[cuenta]++;
  }
    
  for(cuenta=0; cuenta<10; cuenta++)
  {
    fprintf( fpt, "%d\t", countPoints[cuenta] );
  }
  fprintf( fpt, "\n" );
  
//   for(int i=0; i<3; i++)
//   {
//     for(int j=0; j<2; j++)
//     {
//       lCounts[i][j] = lCounts[i][j+1];
//     }
//   }
//   lCounts[0][2] = countPoints[0];
//   lCounts[1][2] = countPoints[1];
//   lCounts[2][2] = countPoints[2];

  //  cálculo de la mediana. Organización de los datos acendentemente
  //copia de la matriz lCounts
/*  for(int i=0; i<3; i++)
  {
    for(int j=0; j<3; j++)
    {
      copylCounts[i][j] = lCounts[i][j];
    }
  }*/
  //  organización de los datos
//   for(int k=0; k<3; k++)
//   {
//     for(int i=0; i<2; i++)
//     {
//       for(int j=i+1;j<3; j++)
//       {
//         if( copylCounts[k][i] > copylCounts[k][j] )
//         {
//           temp = copylCounts[k][j];
//           copylCounts[k][j] = copylCounts[k][i];
//           copylCounts[k][i] = temp;
//         }
//       }
//     }
//   }
//   
//   medianlCounts[0] = copylCounts[0][1];
//   medianlCounts[1] = copylCounts[1][1];
//   medianlCounts[2] = copylCounts[2][1];
//   
//   cout << "PorcentajePuntos: " << 
//       (medianlCounts[0])*100/nop << "%, " << 
//       (medianlCounts[0]+medianlCounts[1])*100/nop << "%, " << (medianlCounts[0]+medianlCounts[1]+medianlCounts[2])*100/nop <<"%" << endl; 
// 
//   fprintf( fpt, "%d, %d, %d;\n", medianlCounts[0], medianlCounts[1], medianlCounts[2] );

  return 0;
}
