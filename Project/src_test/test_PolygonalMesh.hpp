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


#endif
