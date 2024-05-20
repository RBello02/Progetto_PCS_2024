#ifndef __FRACTURESLIBRARY_H // Header guards
#define __FRACTURESLIBRARY_H

#include "Eigen/Eigen"
#include "PolygonalMesh.hpp"
#include <limits>

using namespace std;
using namespace Eigen;
using namespace PolygonalLibrary;



//qui andremo a definire la struct della mesh che importiamo
namespace FracturesLibrary{

struct Trace
{
    unsigned int id; //id della traccia
    double len; //lunghezza della traccia
    MatrixXd coordinates_extremes; //matrice 3x2 con le coordinate degli estremi
    unsigned int id_frc1; //id della prima frattura a cui appartiene la traccia
    unsigned int id_frc2; //id della seconda frattura a cui appartiene la traccia

    inline double calcolo_lunghezza(); //metodo che mi calcola la lunghezza di una traccia
};


struct Fracture
{
    unsigned int id; //id della frattura
    unsigned int num_vertici; //numero di vertici della frattura
    vector<unsigned int> vertices; //vector con gli id dei vertici

    inline Matrix3d calcolo_piano(const vector<Vector3d>& coord); //metodo per il calcolo del piano contenente la frattura
};



struct FracturesFunctions {

    double eps_macchina  = numeric_limits<double>::epsilon();
    const double tolleranza1D = max(pow(10,-16), eps_macchina);
    const double tolleranza2D = max(eps_macchina, pow(tolleranza1D, 2) * 0.75);
    const double tolleranza3D = max(eps_macchina, sqrt(2)/12 * pow(tolleranza1D,3));

    bool importData(const string& path, vector<Fracture>& lista, vector<Vector3d>& coord);
    bool NearFractures(const Fracture& frc1, const Fracture& frc2, const vector<Vector3d>& coord);
    void IntersectionFractures(Fracture &frc1, Fracture &frc2, const vector<Vector3d>& coord, list<Trace>& list_traces, map<unsigned int, list<Trace>>& P_traces, map<unsigned int, list<Trace>>& NP_traces);
    inline bool Parallelismo(const Matrix3d& piano_1, const Matrix3d& piano_2);
    PolygonalMesh SottoPoligonazione(const Fracture& frattura, const list<Trace>& P_traces, const list<Trace>& NP_traces, const vector<Vector3d>& coord);

    //funzioni di supporto
    inline MatrixXd Retta_tra_piani(const Matrix3d& piano_1, const Matrix3d& piano_2);
    inline MatrixXd Retta_per_due_vertici_della_frattura(unsigned int id_vertice1, unsigned int id_vertice2, const vector<Vector3d>& coord);
    inline Vector2d alpha_di_intersezione(MatrixXd r_frattura, MatrixXd retta_intersez);
    inline MatrixXd Retta_per_due_punti(Vector3d& pt1, Vector3d& pt2);
    inline Vector2d intersezione_rette(MatrixXd& r_frattura, MatrixXd& r_traccia);
    inline bool appartiene_a_segmento(Vector3d& origin, Vector3d& end, Vector3d& pto);
    inline bool pto_unico(Vector3d& pto, vector<Vector3d>& punti, unsigned int& id);

};


}// namespace


#endif
