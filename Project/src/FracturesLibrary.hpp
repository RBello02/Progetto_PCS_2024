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


struct Fractures
{
    vector<Vector3d> coordinates; //memorizzo tutti i vertici dei poligoni
    unsigned int num_fractures; //numero di fratture, aka poligoni
    vector<unsigned int> id_fractures; //vettore che contiene per ogni frattura il suo id
    vector<unsigned int> dim_fractures; //vettore che contiene per ogni frattura il numero di vertici della frattura
    vector<vector<unsigned int>> vertices_fractures; //per ogni frattura mi salvo la matrice con gli identificativi dei vertici
    map<unsigned int, list<Trace>> P_traces_of_fractures; //per ogni frattura memorizziamo una lista che contiene
        //degli array che contengono le tracce appartenenti alle fratture e se sono passanti o no
    map<unsigned int, list<Trace>> NP_traces_of_fractures; //analogo a sopra ma per tracce non passanti

};




bool importData(const string& path, Fractures& fract);

bool NearFractures(const Fractures& frc, unsigned int id_fract1, unsigned int id_fract2);

void IntersectionFractures(Fractures& frc, unsigned int id_fract1, unsigned int id_fract2, list<Trace>& list_traces);



} //namespace

#endif
