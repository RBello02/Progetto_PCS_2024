#ifndef __POLYGONALMESH_H // Header guards
#define __POLYGONALMESH_H

#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;


namespace PolygonalLibrary {

struct PolygonalMesh{

    unsigned int NumberCell0D = 0; ///< number of Cell0D
    vector<unsigned int> Cell0DId = {}; ///< Cell0D id, size 1 x NumberCell0D: contiene gli id dei punti
    vector<Vector3d> Cell0DCoordinates = {}; ///< Cell0D coordinates, size 3 x NumberCell0D (x,y,z): coordinate dei punti

    unsigned int NumberCell1D = 0; ///< number of Cell1D
    vector<unsigned int> Cell1DId = {}; ///< Cell1D id, size 1 x NumberCell1D: contiene gli id dei lati
    vector<Vector2i> Cell1DVertices = {}; ///< Cell1D vertices indices, size 2 x NumberCell1D (fromId,toId): contiene l'id della testa e della coda del segmento
    vector<bool> Cell1DStatus = {};

    unsigned int NumberCell2D = 0; ///< number of Cell2D
    vector<unsigned int> Cell2DId = {}; ///< Cell2D id, size 1 x NumberCell2D: id celle 2D
    vector<vector<unsigned int>> Cell2DVertices = {}; ///< Cell2D Vertices indices, size 1 x NumberCell2DVertices[NumberCell2D]
            ///< contine gli id dei vertici del poligono
    vector<vector<unsigned int>> Cell2DEdges = {}; ///< Cell2D Cell1D indices, size 1 x NumberCell2DEdges[NumberCell2D]
            ///< contine gli id dei segmenti del poligono
};


} //namespace

#endif
