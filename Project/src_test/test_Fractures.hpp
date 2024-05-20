#ifndef __TEST_FRACTURES_H // Header guards
#define __TEST_FRACTURES_H

#include "Eigen/Eigen"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <math.h>
#include "FracturesLibrary.hpp"
#include <limits>


using namespace std;
using namespace Eigen;
using namespace FracturesLibrary;



// ***************retta_per_vertici**************************+
TEST(Retta_per_due_vertici, poligono1){
    FracturesFunctions fx;
    vector<Vector3d> coord={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0}};
    Fracture f;
    f.vertices={0,1,2,3};
    MatrixXd risultato_funzione;
    risultato_funzione.resize(2,3);
    risultato_funzione=fx.Retta_per_due_vertici_della_frattura(0,1,coord);
    MatrixXd Risultato_atteso;
    Risultato_atteso.resize(2,3);
    Risultato_atteso << 4,4,0,
        0,-4,0;

    //retta 1 BA
    for (unsigned int i = 0; i < Risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < Risultato_atteso.rows(); j++){
            ASSERT_NEAR(risultato_funzione(j,i), Risultato_atteso(j,i), fx.tolleranza1D);
            }
    }

    risultato_funzione=fx.Retta_per_due_vertici_della_frattura(1,2,coord);
    Risultato_atteso << -4,4,0,
        4,0,0;
    // retta CB
    for (unsigned int i = 0; i < Risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < Risultato_atteso.rows(); j++){
            ASSERT_NEAR(risultato_funzione(j,i), Risultato_atteso(j,i), fx.tolleranza1D);
        }
    }

    risultato_funzione=fx.Retta_per_due_vertici_della_frattura(2,3,coord);
    Risultato_atteso << -4,-4,0,
        0,4,0;
    // retta CD
    for (unsigned int i = 0; i < Risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < Risultato_atteso.rows(); j++){
            ASSERT_NEAR(risultato_funzione(j,i), Risultato_atteso(j,i), fx.tolleranza1D);
        }
    }

    risultato_funzione=fx.Retta_per_due_vertici_della_frattura(3,0,coord);
    Risultato_atteso << 4,-4,0,
        -4,0,0;
     // retta DA
    for (unsigned int i = 0; i < Risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < Risultato_atteso.rows(); j++){
            ASSERT_NEAR(risultato_funzione(j,i), Risultato_atteso(j,i), fx.tolleranza1D);
        }
    }

}

//*************calcolo_piano******************
TEST(Calcolo_piano_test, generale_poligono1){
    vector<Vector3d> coord={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0}};
    Fracture f;
    f.vertices={0,1,2};
    FracturesFunctions fx;

    Matrix3d risultato;
    risultato=f.calcolo_piano(coord);
    Matrix3d A;
    A<<0,-4,0,
       0,8,0,
        4,4,0;
    for (unsigned int i = 0; i < risultato.cols(); i++){
        for (unsigned int j = 0; j < risultato.rows(); j++){
            ASSERT_NEAR(risultato(j,i), A(j,i), fx.tolleranza1D);
        }
    }
}

TEST(Calcolo_piano_test, generale_poligono2){
    FracturesFunctions fx;
    Fracture f;
    f.vertices={0,1,2};
    vector<Vector3d> coord1={{0,0,-2},{2,0,0},{0,0,2},{-2,0,0}};
    Matrix3d risultato;
    risultato=f.calcolo_piano(coord1);
    Matrix3d B;
    B<<0,0,-2,
        0,0,4,
        2,0,2;
    for (unsigned int i = 0; i < risultato.cols(); i++){
        for (unsigned int j = 0; j < risultato.rows(); j++){
            ASSERT_NEAR(risultato(j,i), B(j,i), fx.tolleranza1D);
        }
    }

}

//***********retta_tra_piani**************
TEST(Retta_tra_piani_test, generale_reny){
    FracturesFunctions fx;
    Matrix3d A1, A2;
    A1 <<0,-4,0,
        0,8,0,
        4,-4,0;
    A2 <<0,0,-2,
        0,0,4,
        2,0,2;
    MatrixXd Risultato_funzione;
    Risultato_funzione.resize(2,3);
    Risultato_funzione=fx.Retta_tra_piani(A1,A2);
    MatrixXd risultato_atteso;
    risultato_atteso.resize(2,3);
    risultato_atteso << -1/sqrt(2),0,0,
        0,0,0;

    for (unsigned int i = 0; i < risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < risultato_atteso.rows(); j++){
            ASSERT_NEAR(Risultato_funzione(j,i), risultato_atteso(j,i), fx.tolleranza1D);
        }
    }

}

//**********+alpha_intersezione*************
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
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={1,-0.015625};


    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }

    //EXPECT_EQ(risultato_atteso, Risultato_funzione);

    //***************** secondo lato:
    A1 << -4,4,0,
        4,0,0;
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    risultato_atteso={0,-0.015625};

    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }

    //************************* terzo lato:
    A1 << -4,-4,0,
        0,4,0;
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    risultato_atteso={1,0.015625}; // sono scambiati rispetto a quelli di rena

    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }

    //************************** quarto lato:
    A1 << 4,-4,0,
        -4,0,0;
        risultato_atteso={0,0.015625};

    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }

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
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={1,-0.0078125}; // sono scambiati rispetto a quelli di renato

    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }


    //***************** secondo lato:
    A1 << -2,0,2,
        2,0,0;
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    risultato_atteso={0,-0.0078125};

    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }

    //************************* terzo lato:
    A1 << -2,0,-2,
                 0,0,2;
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    risultato_atteso={1, 0.0078125};

    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }

    //************************** quarto lato:
    A1 << 2,0,-2,
        -2,0,0;
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    risultato_atteso={0, 0.0078125};

    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }

}

TEST(alpha_intersez_test,generale_sofi){
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
    Risultato_funzione=fx.alpha_di_intersezione(A1,A2);
    Vector2d risultato_atteso={-1,0};

    for(unsigned int i = 0; i < Risultato_funzione.size(); i++){
        ASSERT_DOUBLE_EQ(risultato_atteso[i], Risultato_funzione[i]);
    }


}

//****************+parallelismo***************+
TEST(Parallelismo_test, generale_reny){
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

TEST(Parallelismo_test, generale_sofi){
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

// ******************++calcolo_lunghezza********************
TEST(clacolo_lung_test, generale_reny){

    Trace t;
    t.coordinates_extremes.resize(3,2);
    t.coordinates_extremes << 2, -2,
        0, 0,
        0, 0;
    double b;
    b=t.calcolo_lunghezza();

    ASSERT_DOUBLE_EQ(b,4);
}

TEST(clacolo_lung_test, generale_sofi){

    Trace t;
    t.coordinates_extremes.resize(3,2);
    t.coordinates_extremes << 1, 0,
                              0, 1,
                              0, 0;
    double b;
    b=t.calcolo_lunghezza();

    ASSERT_DOUBLE_EQ(b,sqrt(2));
}

//************************intersezioni**************
TEST(NearFracture_test, generale_reny){
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

TEST(NearFracture_test, genarale_sofi){
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
    EXPECT_FALSE(Risultato_funzione);

}

TEST(NearFracture_test, sfere_tangenti){
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    vector<Vector3d> coordinate;
    coordinate={{0,1,1}, {0, -1,1}, {0,1,-1},{0,-1,-1},{2*sqrt(2),1,1}, {2*sqrt(2), -1,1}, {2*sqrt(2),1,-1},{2*sqrt(2),-1,-1} };
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    bool Risultato_funzione;
    Risultato_funzione=g.NearFractures(f1,f2,coordinate);
    EXPECT_FALSE(Risultato_funzione);

}
#endif
