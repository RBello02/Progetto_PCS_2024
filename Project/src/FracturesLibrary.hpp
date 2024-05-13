#ifndef __FRACTURESLIBRARY_H // Header guards
#define __FRACTURESLIBRARY_H

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

//qui andremo a definire la struct della mesh che importiamo
namespace FracturesLibrary {


struct Fractures
{
    vector<Vector3d> coordinates; //memorizzo tutti i vertici dei poligoni
    unsigned int num_fractures; //numero di fratture, aka poligoni
    vector<unsigned int> id_fractures; //vettore che contiene per ogni frattura il suo id
    vector<unsigned int> dim_fractures; //vettore che contiene per ogni frattura il numero di vertici della frattura
    vector<vector<unsigned int>> vertices_fractures; //per ogni frattura mi salvo la matrice con gli identificativi dei vertici
    map<unsigned int, list<array<double, 2>>> P_traces_of_fractures; //per ogni frattura memorizziamo una lista che contiene
        //degli array che contengono le tracce appartenenti alle fratture e se sono passanti o no
    map<unsigned int, list<array<double, 2>>> NP_traces_of_fractures; //analogo a sopra ma per tracce non passanti

};

struct Traces
{
    unsigned int num_traces; //numero di tracce
    vector<map<unsigned int, bool>> fractures_of_traces; //per ogni traccia associo una mappa dove le chiavi sono le due fratture coinvolte e associo il bool per idre se passante o no
    vector<MatrixXd> extreme_traces; //vettore che per ogni traccia si salva la matrice con i due punti estremi

};

bool importData(const string& path, Fractures& fract);

bool NearFractures(const Fractures& frc, unsigned int id_fract1, unsigned int id_fract2);

bool IntersectionFractures(Fractures& frc, unsigned int id_fract1, unsigned int id_fract2);



} //namespace

#endif
