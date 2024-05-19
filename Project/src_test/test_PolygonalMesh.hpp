#ifndef __TEST_POLYGONALMESH_H // Header guards
#define __TEST_POLYGONALMESH_H

#include "Eigen/Eigen"
#include <gtest/gtest.h>
#include <math.h>
#include "FracturesLibrary.hpp"
#include "PolygonalMesh.hpp"


using namespace std;
using namespace Eigen;
using namespace FracturesLibrary;




//************************* ho lasciato questi test per fare copia incolla, sono da eliminare


//****************+parallelismo funziona***************+
TEST(Parallelismo_test, generale){
    Matrix3d A1, A2;
    A1 << 0,-4,0,
        4,-4,0,
        4,4,0;
    A2 << 0,0,-2,
        0,0,4,
        2,0,2;
    bool b;
    FracturesFunctions f;
    b=f.Parallelismo(A1,A2);

    ASSERT_FALSE(b);
}
TEST(Parallelismo_test, prova){
    Matrix3d A1, A2;
    A1 << 0,-4,0,
        4,-4,0,
        4,4,0;
    A2 << 0,-2,0,
        2,-2,0,
        2,2,0;
    bool b;
    FracturesFunctions f;
    b=f.Parallelismo(A1,A2);

    ASSERT_TRUE( b);
}

// ******************++calcolo_lunghezza funziona********************
TEST(clacolo_lung_test, generale){

    Trace t;
    t.coordinates_extremes.resize(3,2);
    t.coordinates_extremes << 2, -2,
        0, 0,
        0, 0;
    double b;
    b=t.calcolo_lunghezza();

    ASSERT_DOUBLE_EQ(b,4);
}
TEST(clacolo_lung_test, prova){

    Trace t;
    t.coordinates_extremes.resize(3,2);
    t.coordinates_extremes << 1, 0,
                              0, 1,
                              0, 0;
    double b;
    b=t.calcolo_lunghezza();

    ASSERT_DOUBLE_EQ(b,sqrt(2));
}

//************************intersezioni: funziona**************
TEST(NearFracture_test, generale){
    //bool NearFractures(const Fracture& frc1, const Fracture& frc2, const vector<Vector3d>& coord);
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    vector<Vector3d> coordinate;
    coordinate={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0},{0,0,-2},{2,0,0},{0,0,2},{-2,0,0}};
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    bool Risultato_funzione;
    Risultato_funzione=g.NearFractures(f1,f2,coordinate);
    ASSERT_TRUE(Risultato_funzione);

}
// nel mio cè un problema
TEST(NearFracture_test, prova){
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    vector<Vector3d> coordinate;
    coordinate={{0,1,0},{0,0,1},{1,0,0},{30,30,-30},{28,30,-28},{28,28,-28},{30,28,-30}};
    f1.num_vertici=3;
    f2.num_vertici=4;
    f1.vertices={0,1,2};
    f2.vertices={3,4,5,6};
    bool Risultato_funzione;
    Risultato_funzione=g.NearFractures(f1,f2,coordinate);
    ASSERT_FALSE(Risultato_funzione);

}

#endif
