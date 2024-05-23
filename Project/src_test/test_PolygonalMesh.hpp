#ifndef __TEST_POLYGONALMESH_H // Header guards
#define __TEST_POLYGONALMESH_H

#include "Eigen/Eigen"
#include <gtest/gtest.h>
#include "test_Fractures.hpp"


using namespace UtilsFunction;

TEST(Retta_per_due_punti_test, generale_reny){
    FracturesFunctions fx;
    Vector3d pt1,pt2;
    pt1={0,-4,0};
    pt2={0,0,-2};
    MatrixXd risultato_funzione;
    risultato_funzione.resize(2,3);
    risultato_funzione=fx.Retta_per_due_punti(pt1,pt2);
    MatrixXd Risultato_atteso;
    Risultato_atteso.resize(2,3);
    Risultato_atteso << 0,4,-2,
        0,-4,0;

    for (unsigned int i = 0; i < Risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < Risultato_atteso.rows(); j++){
            ASSERT_NEAR(risultato_funzione(j,i), Risultato_atteso(j,i), fx.tolleranza1D);
        }
    }
    ASSERT_EQ(Risultato_atteso,risultato_funzione);


}

TEST(Retta_per_due_punti_test, punti_uguali){
    FracturesFunctions fx;
    Vector3d pt1,pt2;
    pt1={0,-4,0};
    pt2={0,-4,0};
    MatrixXd risultato_funzione;
    risultato_funzione.resize(2,3);
    risultato_funzione=fx.Retta_per_due_punti(pt1,pt2);
    MatrixXd Risultato_atteso;
    Risultato_atteso.resize(2,3);
    Risultato_atteso << 0,0,0,
        0,-4,0;

    for (unsigned int i = 0; i < Risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < Risultato_atteso.rows(); j++){
            ASSERT_NEAR(risultato_funzione(j,i), Risultato_atteso(j,i), fx.tolleranza1D);
        }
    }

    // dato che la direzione è nulla non ha senso la retta, capire cosa fa il programma in questo caso, se lo gestisce prima di usare la funzione


}


TEST(intersezione_rette_test, generale_sofi){
    FracturesFunctions fx;
    MatrixXd A1, A2;
    A1.resize(2,3);
    A2.resize(2,3);
    A1 << 2,0,1,
        3,1,2;
    A2 << 5,4,0,
        1,1,1;
    Vector2d Risultato_funzione;
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={-1,0};
    for (unsigned int i = 0; i < risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < risultato_atteso.rows(); j++){
            ASSERT_NEAR(Risultato_funzione(j,i), risultato_atteso(j,i), fx.tolleranza1D);
        }
    }

}

TEST(appartiene_a_segmento_test, genarale_sofi){
    FracturesFunctions fx;
    bool risultato_funzione;
    Vector3d origine,fine,punto;
    origine={1,2,3};
    fine={9,8,7};
    punto={2,4,5};
    risultato_funzione=fx.appartiene_a_segmento(origine,fine,punto);
    ASSERT_FALSE(risultato_funzione);
}

//ho aggiunto questo test per controllare che il segmento comprenda sia orgine che fine
TEST(appartiene_a_segmento_test, origine_fine){
    FracturesFunctions fx;
    double perturbazione = fx.tolleranza1D*1./10;
    bool risultato_funzione;
    Vector3d origine,fine,punto;
    origine={1,2,3};
    fine={9,8,7};
    punto={1+perturbazione,2,3};
    risultato_funzione=fx.appartiene_a_segmento(origine,fine,punto);
    ASSERT_TRUE(risultato_funzione);
    punto={9,8,7};
    risultato_funzione=fx.appartiene_a_segmento(origine,fine,punto);
    ASSERT_TRUE(risultato_funzione);

}


TEST(punto_unico_test, generale_sofi){
    FracturesFunctions fx;
    bool risultato_funzione;
    Vector3d pto;
    vector<Vector3d> punti;
    double perturbazione = fx.tolleranza1D*1./10;
    pto={2 + perturbazione,4,5};
    punti={{1,4,5},{5,7,8},{9,20,4},{2,4,5}};
    unsigned int id=0;
    risultato_funzione=fx.pto_unico(pto,punti, id);
    ASSERT_FALSE(risultato_funzione);
    EXPECT_EQ(id, 3);
}

TEST(punto_unico_test, generale_marti){
    FracturesFunctions fx;
    bool risultato_funzione;
    Vector3d pto;
    vector<Vector3d> punti;
    pto={2,4,6};
    punti={{1,4,5},{5,7,8},{9,20,4},{2,4,5}};
    unsigned int id=0;
    risultato_funzione=fx.pto_unico(pto,punti, id);
    ASSERT_TRUE(risultato_funzione);
    EXPECT_EQ(id, 0);
}


//test sulla sottopoloigonazione
TEST(sottopoligonazione, generale_marti){
    FracturesFunctions fx;
    PolygonalMesh mesh;

    vector<Vector3d> coordinates;
    Fracture frattura;
    frattura.id = 0;
    frattura.num_vertici = 6;
    frattura.vertices.reserve(6);

    for(unsigned int i = 0 ; i<6; i++){frattura.vertices.push_back(i);}

    coordinates.reserve(6);
    coordinates.push_back({2,9,3});
    coordinates.push_back({5,9,2});
    coordinates.push_back({5.5,9,4.5});
    coordinates.push_back({4,9,6});
    coordinates.push_back({2.5,9,6});
    coordinates.push_back({1,9,5});


    Trace trc1;
    trc1.id = 0;
    trc1.id_frc1 = 0;
    trc1.id_frc2 = 4;

    MatrixXd extr1;
    extr1.resize(3,2);
    Vector3d o1 = {3.49999, 9,2.97};
    Vector3d e1 =  {4.67782,9,4.80000};
    extr1.col(0) = o1;
    extr1.col(1) = e1;

    trc1.coordinates_extremes = extr1;

    Trace trc2;
    trc2.id = 3;
    trc2.id_frc1 = 8;
    trc2.id_frc2 = 0;

    MatrixXd extr2;
    extr2.resize(3,2);
    Vector3d o2 = {2, 9,4.80};
    Vector3d e2 =  {3,9,4.80000};
    extr2.col(0) = o2;
    extr2.col(1) = e2;

    trc2.coordinates_extremes = extr2;

    list<Trace> P_traces;
    list<Trace> NP_traces;
    NP_traces.push_back(trc1);
    NP_traces.push_back(trc2);

    mesh = fx.FracturesFunctions::SottoPoligonazione(frattura, P_traces, NP_traces, coordinates);

    //faccio solo una verifica del numero di celle 0D e 2D
    EXPECT_EQ(mesh.NumberCell0D, 10);
    EXPECT_EQ(mesh.NumberCell2D, 3);
}

TEST(sottopoligonazione, estr_traccia_coincidente_con_vertice){
    FracturesFunctions fx;
    PolygonalMesh mesh;

    vector<Vector3d> coordinates;
    Fracture frattura;
    frattura.id = 0;
    frattura.num_vertici = 5;
    frattura.vertices.reserve(5);

    for(unsigned int i = 0 ; i<5; i++){frattura.vertices.push_back(i);}

    coordinates.reserve(5);
    coordinates.push_back({1,3,0});
    coordinates.push_back({2,1,0});
    coordinates.push_back({4,1,0});
    coordinates.push_back({5,3,0});
    coordinates.push_back({3,4.5,0});


    Trace trc1;
    trc1.id = 0;
    trc1.id_frc1 = 0;
    trc1.id_frc2 = 4;

    MatrixXd extr1;
    extr1.resize(3,2);
    Vector3d o1 = {3,4.5,0};
    Vector3d e1 =  {3,1,0};
    extr1.col(0) = o1;
    extr1.col(1) = e1;

    trc1.coordinates_extremes = extr1;

    Trace trc2;
    trc2.id = 3;
    trc2.id_frc1 = 8;
    trc2.id_frc2 = 0;

    MatrixXd extr2;
    extr2.resize(3,2);
    Vector3d o2 = {1.75, 2,0};
    Vector3d e2 =  {2.6,2.5,0};
    extr2.col(0) = o2;
    extr2.col(1) = e2;

    trc2.coordinates_extremes = extr2;

    list<Trace> P_traces;
    list<Trace> NP_traces;
    P_traces.push_back(trc1);
    NP_traces.push_back(trc2);

    mesh = fx.FracturesFunctions::SottoPoligonazione(frattura, P_traces, NP_traces, coordinates);

    //faccio solo una verifica del numero di celle 0D e 2D
    EXPECT_EQ(mesh.NumberCell0D, 8);
    EXPECT_EQ(mesh.NumberCell2D, 3);
}

TEST(sottopoligonazione, traccia_comune_a_due_sottofratt){
    FracturesFunctions fx;
    PolygonalMesh mesh;

    vector<Vector3d> coordinates;
    Fracture frattura;
    frattura.id = 0;
    frattura.num_vertici = 5;
    frattura.vertices.reserve(5);

    for(unsigned int i = 0 ; i<5; i++){frattura.vertices.push_back(i);}

    coordinates.reserve(5);
    coordinates.push_back({1,3,0});
    coordinates.push_back({2,1,0});
    coordinates.push_back({4,1,0});
    coordinates.push_back({5,3,0});
    coordinates.push_back({3,4.5,0});


    Trace trc1;
    trc1.id = 0;
    trc1.id_frc1 = 0;
    trc1.id_frc2 = 4;

    MatrixXd extr1;
    extr1.resize(3,2);
    Vector3d o1 = {3,4.5,0};
    Vector3d e1 =  {3,1,0};
    extr1.col(0) = o1;
    extr1.col(1) = e1;

    trc1.coordinates_extremes = extr1;

    Trace trc2;
    trc2.id = 3;
    trc2.id_frc1 = 8;
    trc2.id_frc2 = 0;

    MatrixXd extr2;
    extr2.resize(3,2);
    Vector3d o2 = {2, 3.25,0};
    Vector3d e2 =  {4,2.96,0};
    extr2.col(0) = o2;
    extr2.col(1) = e2;

    trc2.coordinates_extremes = extr2;

    list<Trace> P_traces;
    list<Trace> NP_traces;
    P_traces.push_back(trc1);
    NP_traces.push_back(trc2);

    mesh = fx.FracturesFunctions::SottoPoligonazione(frattura, P_traces, NP_traces, coordinates);

    //faccio solo una verifica del numero di celle 0D e 2D
    EXPECT_EQ(mesh.NumberCell0D, 9);
    EXPECT_EQ(mesh.NumberCell2D, 4);
}



#endif
