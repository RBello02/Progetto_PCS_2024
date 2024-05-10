#include "src/GeometryLibrary.hpp"
#include "src/Utils.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include "Eigen/Eigen"
#include <cmath>

using namespace GeometryLibrary;
using namespace Eigen;
using namespace std;

inline Vector4d Equazione_piano(const Fractures& frc, unsigned int id_fract1)
{
    // volevo risolvere un sistema lineare del tipo:
    // |x1 y1 z1 1 | |a| |0|
    // |x2 y2 z2 1 | |b| |0|
    // |x3 y3 z3 1 | |c|=|0|
    //               |d|
    // prendendo le coordinate lei primi 3 punti di una frattura, ma non funziona :(
    MatrixXd A;
    A.resize(3,4);
    Vector3d b=Vector3d::Zero();
    for (unsigned int i=0; i<3; i++)
    {
        unsigned int id_vertice=frc.vertices_fractures[id_fract1][i];
        A.row(i) << frc.coordinates[id_vertice][0],frc.coordinates[id_vertice][1], frc.coordinates[id_vertice][2], 1 ;

    }
    // for (int i = 0; i < 3; ++i) { // la matrice è giusta
    //     for (int j = 0; j < 4; ++j) {
    //         std::cout << A(i,j) << " ";
    //     }
    //     std::cout << std::endl;}

    Vector4d coeff_piano = A.colPivHouseholderQr().solve(b);
    for(unsigned int i =0; i<4;i++)
    {
        cout << coeff_piano[i] << " ";
    }

    return coeff_piano;


}

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

    // ciclo sulle coppie di poligoni e vedo se sono vicini oppure no
    for(unsigned int i=0; i<frc.num_fractures; i++)
    {
        for (unsigned int j=i+1; j< frc.num_fractures; j++)
        {
        unsigned int id_fract1=i;
        unsigned int id_fract2=j;
        if( !NearFractures(frc, id_fract1,id_fract2))
        {cout << "non dobbiamo fare niente"<< endl;}
        else
        {
            cout<<"fare i piani e le rette" << endl;
        }

        }
    }

    // test equzione piano
    unsigned int id_fract=2;
    Vector4d a = Equazione_piano(frc,id_fract);


    return 0;
}
