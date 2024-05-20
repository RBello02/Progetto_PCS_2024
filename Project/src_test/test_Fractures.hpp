#ifndef __TEST_FRACTURES_H // Header guards
#define __TEST_FRACTURES_H

#include "Eigen/Eigen"
#include <gtest/gtest.h>
#include <math.h>
#include "FracturesLibrary.hpp"
#include <limits>


using namespace std;
using namespace Eigen;
using namespace FracturesLibrary;


inline MatrixXd Retta_tra_piani(const Matrix3d& piano_1, const Matrix3d& piano_2)
{
    // calcolo la retta di intersezione tra i piani in forma parametrica
    // X = at+P


    // NOMENCLATURA
    // u1 = piano_1.row(1) v1 = piano_1.row(2) sono i vettori direttori del primo piano
    // u2 = piano_2.row(1) v2 = piano_2.row(2) sono i vettori direttori del secondo piano
    // P0_1 = piano_1.row(0) punto iniziale per il piano 1
    // P0_2 = piano_2.row(0) punto iniziale per il piano 2

    // passo 1: ricavare la normale ai piani
    Vector3d n1 = piano_1.row(1).cross(piano_1.row(2))/(piano_1.row(1).norm()*piano_1.row(2).norm());
    Vector3d n2 = piano_2.row(1).cross(piano_2.row(2))/(piano_2.row(1).norm()*piano_2.row(2).norm());

    // passo 2 : ricavare la direttrice della retta
    Vector3d t = n1.cross(n2);

    // passo 3: ricavare lo scalare d per i due piani
    double d1 = n1.dot(piano_1.row(0));
    double d2 = n2.dot(piano_2.row(0));
    Vector3d d = {d1,d2,0};  // salvo queste quantita in un vettore

    // passo 4: ricavare la matrice A
    Matrix3d A;
    A.row(0) = n1;
    A.row(1) = n2;
    A.row(2) = t;

    // passo 5: risolvere il sistema per ricavare P
    Vector3d P = A.partialPivLu().solve(d);

    // salvo in un formato particolare
    MatrixXd X;
    X.resize(2,3);
    X.row(0) = t.transpose();
    X.row(1) = P.transpose();  // come una matrice 2x3

    return X;

}

inline MatrixXd Retta_per_due_vertici_della_frattura(unsigned int id_vertice1, unsigned int id_vertice2, const vector<Vector3d>& coord)
{
    // l'equazione parametrica è X = at+P

    // passo 1: salvo le coordinate dei due punti
    double x1 = coord[id_vertice1][0];
    double y1 = coord[id_vertice1][1];
    double z1 = coord[id_vertice1][2];

    double x2 = coord[id_vertice2][0];
    double y2 = coord[id_vertice2][1];
    double z2 = coord[id_vertice2][2];

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

inline Vector2d alpha_di_intersezione(MatrixXd r_frattura, MatrixXd retta_intersez)
{

    //imposto un sistema lineare per la ricerca dei parametri alpha e beta
    //primo parametro è la matrice della retta del poligono
    Vector3d t1 = r_frattura.row(0).transpose();
    Vector3d P1 = r_frattura.row(1).transpose();

    //secondo parametro è la matrice della retta di intersezione tra i piani
    Vector3d t2 = retta_intersez.row(0).transpose();
    Vector3d P2 = retta_intersez.row(1).transpose();

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

inline Matrix3d Fracture::calcolo_piano(const vector<Vector3d>& coord)
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
//*************************
// ***************retta_per_vertici_fuziona**************************+
TEST(Retta_per_due_vertici, poligono1){
    vector<Vector3d> coord={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0}};
    Fracture f;
    f.vertices={0,1,2,3};
    MatrixXd risultato_funzione;
    risultato_funzione.resize(2,3);
    risultato_funzione=Retta_per_due_vertici_della_frattura(0,1,coord);
    MatrixXd Risultato_atteso;
    Risultato_atteso.resize(2,3);
    Risultato_atteso << 4,4,0,
        0,-4,0;
    ASSERT_EQ(Risultato_atteso,risultato_funzione); //retta 1 BA
    risultato_funzione=Retta_per_due_vertici_della_frattura(1,2,coord);
    Risultato_atteso << -4,4,0,
        4,0,0;
    ASSERT_EQ(Risultato_atteso,risultato_funzione); // retta CB
    risultato_funzione=Retta_per_due_vertici_della_frattura(2,3,coord);
    Risultato_atteso << -4,-4,0,
        0,4,0;
    ASSERT_EQ(Risultato_atteso,risultato_funzione); // retta CD
    risultato_funzione=Retta_per_due_vertici_della_frattura(3,0,coord);
    Risultato_atteso << 4,-4,0,
        -4,0,0;
    ASSERT_EQ(Risultato_atteso,risultato_funzione); // retta DA

}
//*************calcolo_piano non funziona******************
TEST(Calcolo_piano_test, generale_poligono1){
    vector<Vector3d> coord={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0}};
    Fracture f;
    f.vertices={0,1,2};

    Matrix3d risultato;
    risultato=f.calcolo_piano(coord);
    Matrix3d A;
    A<<0,-4,0,
       0,8,0,
        4,4,0;
    ASSERT_EQ(A,risultato);}
TEST(Calcolo_piano_test, generale_poligono2){

    Fracture f;
    f.vertices={0,1,2};
    vector<Vector3d> coord1={{0,0,-2},{2,0,0},{0,0,2},{-2,0,0}};
    Matrix3d risultato;
    risultato=f.calcolo_piano(coord1);
    Matrix3d B;
    B<<0,0,-2,
        0,0,4,
        2,0,2;
    ASSERT_EQ(B,risultato);

}

//***********retta_tra_piani: falliscono entrambi**************
TEST(Retta_tra_piani_test, generale){
    Matrix3d A1, A2;
    A1 << 0,-4,0,
        4,-4,0,
        4,4,0;
    A2 <<0,0,-2,
        0,0,4,
        2,0,2;
    MatrixXd Risultato_funzione;
    Risultato_funzione.resize(2,3);
    Risultato_funzione=Retta_tra_piani(A1,A2);
    MatrixXd risultato_atteso;
    risultato_atteso.resize(2,3);
    risultato_atteso << -256,0,0,
        0,0,0;
    ASSERT_EQ(risultato_atteso, Risultato_funzione);

}

//**********+alpha_intersezione:falliscono quelli generali, da rivedere i risultati*************
TEST(alpha_intersez_test,generale_poligono1){
    FracturesFunctions g;
    MatrixXd A1, A2;
    A1.resize(2,3);
    A2.resize(2,3);
    // A2 è la retta di intersezione tra i piani
    A2 << -256,0,0,
        0,0,0;
    // A1 sono le rette dei lati del poligono
    //***************** primo lato:
    A1 << 4,4,0,
        0,-4,0;
    Vector2d Risultato_funzione;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={1,-0.015625};

    //3.33067e-16 vs 2.22045e-16 mi restituisce questi valori
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);
    //***************** secondo lato:
    A1 << -4,4,0,
        4,0,0;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    risultato_atteso={0,-0.015625};
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);
    //************************* terzo lato:
    A1 << -4,-4,0,
        0,4,0;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    risultato_atteso={1,0.015625}; // sono scambiati rispetto a quelli di rena
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);
    //************************** quarto lato:
    A1 << 4,-4,0,
        -4,0,0;
        risultato_atteso={0,0.015625};
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);

}
TEST(alpha_intersez_test,generale_poligono2){
    FracturesFunctions g;
    MatrixXd A1, A2;
    A1.resize(2,3);
    A2.resize(2,3);
    A2 << -256,0,0,
        0,0,0;
    //***************** secondo lato:
    A1 << 2,0,2,
        0,0,-2;
    Vector2d Risultato_funzione;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={1,-0.0078125}; // sono scambiati rispetto a quelli di renato
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);
    //***************** secondo lato:
    A1 << -2,0,2,
        2,0,0;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    risultato_atteso={0,-0.0078125};
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);
    //************************* terzo lato:
    A1 << -2,0,-2,
                 0,0,2;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    risultato_atteso={1, 0.0078125};
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);
    //************************** quarto lato:
    A1 << 2,0,-2,
        -2,0,0;
    risultato_atteso={0, 0.0078125};
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),1e-15);

}
//** il mio funziona ma forse perchè sono int(?)
TEST(alpha_intersez_test,prova){
    FracturesFunctions g;
    MatrixXd A1, A2;
    A1.resize(2,3);
    A2.resize(2,3);
    A1 << 2,0,1,
        3,1,2;
    //***************** secondo lato:
    A2 << 5,4,0,
        1,1,1;
    Vector2d Risultato_funzione;
    Risultato_funzione=alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={-1,0};
    ASSERT_LE((risultato_atteso-Risultato_funzione).norm(),g.tolleranza1D);

}

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
