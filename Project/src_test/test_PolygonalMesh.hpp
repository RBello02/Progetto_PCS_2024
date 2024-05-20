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

inline MatrixXd Retta_per_due_punti(Vector3d& pt1, Vector3d& pt2)
{
    // l'equazione parametrica è X = at+P

    // passo 1: salvo le coordinate dei due punti
    double x1 = pt1[0];
    double y1 = pt1[1];
    double z1 = pt1[2];

    double x2 = pt2[0];
    double y2 = pt2[1];
    double z2 = pt2[2];

    // passo 2: trovo direttrice e punto di partenza della retta

    Vector3d t = {x2-x1,y2-y1,z2-z1};
    Vector3d P = {x1,y1,z1};

    // salvo in un formato particolare
    MatrixXd X;
    X.resize(2,3);
    X.row(0) = t.transpose();
    X.row(1) = P.transpose();  // come una matrice 2x3

    return X;

}

inline Vector2d intersezione_rette(MatrixXd& r_frattura, MatrixXd& r_traccia)
{

    //imposto un sistema lineare per la ricerca dei parametri alpha e beta
    //primo parametro è la matrice della retta del poligono --> retta in funzione di alpha
    Vector3d t1 = r_frattura.row(0).transpose();
    Vector3d P1 = r_frattura.row(1).transpose();

    //secondo parametro è la matrice della retta della traccia --> retta in funzione di beta
    Vector3d t2 = r_traccia.row(0).transpose();
    Vector3d P2 = r_traccia.row(1).transpose();

    MatrixXd A = MatrixXd::Zero(3,2);
    Vector3d b = Vector3d::Zero();

    //imposto i coefficienti della matrice e del termine noto
    A.col(0) = t1;
    A.col(1) = -t2;

    for (unsigned int i = 0; i<3; i++){b[i] = P2[i]-P1[i];}

    Vector2d x = A.householderQr().solve(b); //x =[alpha; beta]

    Vector2d alpha_beta = x;
    return alpha_beta;

}

inline double appartiene_a_segmento(Vector3d& origin, Vector3d& end, Vector3d& pto, double toll){
    bool appartiene = false;

    //calcolo la distanza di pto dai due estremi del segmento, se la somma di queste due distanza è maggiore della distanza dei
    // due estremi allora il pto è esterno
    double lung_segmento = (origin-end).norm();
    double pto_o = (origin-pto).norm();
    double pto_e = (end-pto).norm();

    if(abs(pto_o + pto_e - lung_segmento) < toll){appartiene = true;}

    return appartiene;

}

inline double pto_unico(Vector3d& pto, vector<Vector3d>& punti, double toll, unsigned int& id){
    //l'id mi serve quando cerco un punto tra i vertici

    bool unico = true;

    if(punti.size() == 0){return unico;}

    for (unsigned int i = 0; i < punti.size(); i++){
        Vector3d elem = punti[i];
        bool uguagl_x = (abs(pto[0] - elem [0]) < toll);
        bool uguagl_y = (abs(pto[1] - elem [1]) < toll);
        bool uguagl_z = (abs(pto[2] - elem [2]) < toll);

        if (uguagl_x && uguagl_y && uguagl_z){unico = false; id = i; return unico;}
    }
    return unico;
}
    TEST(Retta_per_due_punti_test, generale){
    Vector3d pt1,pt2;
    pt1={0,-4,0};
    pt2={0,0,-2};
    MatrixXd risultato_funzione;
    risultato_funzione.resize(2,3);
    risultato_funzione=Retta_per_due_punti(pt1,pt2);
    MatrixXd Risultato_atteso;
    Risultato_atteso.resize(2,3);
    Risultato_atteso << 0,4,-2,
        0,-4,0;
    ASSERT_EQ(Risultato_atteso,risultato_funzione);


}
TEST(intersezione_rette_test,prova){
    FracturesFunctions g;
    MatrixXd A1, A2;
    A1.resize(2,3);
    A2.resize(2,3);
    A1 << 2,0,1,
        3,1,2;
    A2 << 5,4,0,
        1,1,1;
    Vector2d Risultato_funzione;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={-1,0};
    ASSERT_EQ(risultato_atteso,Risultato_funzione);

}
TEST(appartiene_a_segmento_test,generale){
    bool risultato_funzione;
    Vector3d origine,fine,punto;
    double toll=1e-15;
    origine={1,2,3};
    fine={9,8,7};
    punto={2,4,5};
    risultato_funzione=appartiene_a_segmento(origine,fine,punto,toll);
    ASSERT_FALSE(risultato_funzione);
}
TEST(punto_unico_test,generale){
    bool risultato_funzione;
    Vector3d pto;
    vector<Vector3d> punti;
    double toll=1e-15;
    pto={2,4,5};
    punti={{1,4,5},{5,7,8},{9,20,4},{2,4,5}};
    unsigned int id=0;
    risultato_funzione=pto_unico(pto,punti,toll, id);
    ASSERT_FALSE(risultato_funzione);
}


#endif
