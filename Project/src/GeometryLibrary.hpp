#ifndef __GEOMETRYLIBRARY_H // Header guards
#define __GEOMETRYLIBRARY_H

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

//qui andremo a definire la struct della mesh che importiamo
namespace GeometryLibrary {

struct Fractures
{
    vector<array<double,3>> coordinates; //memorizzo tutti i vertici dei poligoni
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

struct PolygonalMesh{

    unsigned int NumberCell0D = 0; ///< number of Cell0D
    vector<unsigned int> Cell0DId = {}; ///< Cell0D id, size 1 x NumberCell0D
    vector<Vector3d> Cell0DCoordinates = {}; ///< Cell0D coordinates, size 3 x NumberCell0D (x,y,z)

    unsigned int NumberCell1D = 0; ///< number of Cell1D
    vector<unsigned int> Cell1DId = {}; ///< Cell1D id, size 1 x NumberCell1D
    vector<Vector2i> Cell1DVertices = {}; ///< Cell1D vertices indices, size 2 x NumberCell1D (fromId,toId)

    unsigned int NumberCell2D = 0; ///< number of Cell2D
    vector<unsigned int> Cell2DId = {}; ///< Cell2D id, size 1 x NumberCell2D
    vector<vector<unsigned int>> Cell2DVertices = {}; ///< Cell2D Vertices indices, size 1 x NumberCell2DVertices[NumberCell2D]
    vector<vector<unsigned int>> Cell2DEdges = {}; ///< Cell2D Cell1D indices, size 1 x NumberCell2DEdges[NumberCell2D]
};


} //namespace

#endif
