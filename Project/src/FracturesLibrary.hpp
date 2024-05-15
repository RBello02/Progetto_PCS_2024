#ifndef __FRACTURESLIBRARY_H // Header guards
#define __FRACTURESLIBRARY_H

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

//qui andremo a definire la struct della mesh che importiamo
namespace FracturesLibrary {



struct Trace
{
    unsigned int id; //id della traccia
    double len; //lunghezza della traccia
    MatrixXd coordinates_extremes; //matrice 3x2 con le coordinate degli estremi
    unsigned int id_frc1;
    unsigned int id_frc2;

    inline double calcolo_lunghezza(); //metodo che mi calcola la lunghezza di una traccia
};


struct Fracture
{
    unsigned int id;
    unsigned int num_vertici;
    vector<unsigned int> vertices;

    inline Matrix3d calcolo_piano(const vector<Vector3d>& coord);
};



bool importData(const string& path, vector<Fracture>& lista, vector<Vector3d>& coord);

bool NearFractures(const Fracture& frc1, const Fracture& frc2, const vector<Vector3d>& coord);

void IntersectionFractures(Fracture &frc1, Fracture &frc2, const vector<Vector3d>& coord, list<Trace>& list_traces, map<unsigned int, list<Trace>>& P_traces, map<unsigned int, list<Trace>>& NP_traces);



} //namespace

#endif
