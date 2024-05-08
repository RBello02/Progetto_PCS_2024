#include "src/GeometryLibrary.hpp"
#include "src/Utils.hpp"

#include <iomanip>
#include <iostream>
#include <string>

using namespace GeometryLibrary;

int main()
{
    cout << setprecision(16);
    string path = "DFN";
    string filename = path + "/FR3_data.txt";
    Fractures frc;

    if (!importData(filename, frc)){return 1;}
    else{
        //stampo un poo' di roba per verificare che sia tutto giusto
        cout << "Ho un numero di fratture pari a " << frc.num_fractures << endl;

        for (unsigned int i = 0; i < frc.num_fractures ; i++){
            cout << "FRATTURA NUM " << i <<  " che ha " << frc.dim_fractures[i] << " vertici." << endl;

            for (unsigned int k = 0; k < frc.dim_fractures[i]; k++){

                 unsigned int id_vertice = frc.vertices_fractures[i][k];
                cout << "Il vertice " << k << ", che ha id = " << id_vertice << " , ha coord (" << frc.coordinates[id_vertice][0] << " , "  <<frc.coordinates[id_vertice][1] <<  " , " <<frc.coordinates[id_vertice][2] << " )" << endl;
             }
            cout << endl;
        }


    }

  return 0;
}
