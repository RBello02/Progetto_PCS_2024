#ifndef __FRACTURESLIBRARY_H // Header guards
#define __FRACTURESLIBRARY_H

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


//qui andremo a definire la struct della mesh che importiamo
namespace FracturesLibrary{

struct Trace
{
    unsigned int id; //id della traccia
    double len; //lunghezza della traccia
    MatrixXd coordinates_extremes; //matrice 3x2 con le coordinate degli estremi
    unsigned int id_frc1; //id della prima frattura a cui appartiene la traccia
    unsigned int id_frc2; //id della seconda frattura a cui appartiene la traccia

    //metodo che mi calcola la lunghezza di una traccia
    inline double calcolo_lunghezza(){
        Vector3d origin = coordinates_extremes.col(0);
        Vector3d end = coordinates_extremes.col(1);
        double val = sqrt( (origin[0]-end[0])*(origin[0]-end[0]) + (origin[1]-end[1])*(origin[1]-end[1]) + (origin[2]-end[2])*(origin[2]-end[2]));
        return val;
    }
};


struct Fracture
{
    unsigned int id; //id della frattura
    unsigned int num_vertici; //numero di vertici della frattura
    vector<unsigned int> vertices; //vector con gli id dei vertici

    //metodo per il calcolo del piano contenente la frattura
    inline Matrix3d calcolo_piano(const vector<Vector3d>& coord)
    {
        // salvo il piano in forma parametrica X = P0 + a1(P2-P0) + a2(P1-P0), dove a1 e a2 sono le parametrizzazioni
        // il piano è quindi identificato da 3 punti P0,P1,P2


        Matrix3d A;
        Vector3d P0 = coord[vertices[0]];
        Vector3d P1 = coord[vertices[1]];
        Vector3d P2 = coord[vertices[2]];

        A.row(0) = P0;
        A.row(1) = P2-P0;
        A.row(2) = P1-P0;


        return A; // restituisce una matrice con la seguente struttura   (P0;P1;P2)  dove Pi sono VETTORI RIGA
    }

};


}// namespace


#endif
