#ifndef __TEST_FRACTURES_H // Header guards
#define __TEST_FRACTURES_H

#include "Eigen/Eigen"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <math.h>
#include "FracturesLibrary.hpp"
#include "Utils.hpp"



using namespace std;
using namespace Eigen;
using namespace FracturesLibrary;
using namespace UtilsFunction;


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
        4,4,0;
    A2 <<0,0,-2,
        0,0,4,
        2,0,2;
    MatrixXd Risultato_funzione;
    Risultato_funzione.resize(2,3);
    Risultato_funzione=fx.Retta_tra_piani(A1,A2);
    MatrixXd risultato_atteso;
    risultato_atteso.resize(2,3);
    risultato_atteso << 1./2.,0,0,
        0,0,0;

    for (unsigned int i = 0; i < risultato_atteso.cols(); i++){
        for (unsigned int j = 0; j < risultato_atteso.rows(); j++){
            ASSERT_NEAR(Risultato_funzione(j,i), risultato_atteso(j,i), fx.tolleranza1D);
        }
    }

}



//**********+alpha_intersezione*************
TEST(alpha_intersez_test,generale_poligono1){
    FracturesFunctions fx;
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
    FracturesFunctions fx;
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
    FracturesFunctions fx;
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

// ******************calcolo_lunghezza********************
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

// ho aggiunto questo test per vedere che succede se do due punti uguali, quindi lunghezza della traccia nulla
TEST(clacolo_lung_test, lunghezza_nulla){

    Trace t;
    t.coordinates_extremes.resize(3,2);
    t.coordinates_extremes << 1, 1,
        0, 0,
        0, 0;
    double b;
    b=t.calcolo_lunghezza();

    ASSERT_DOUBLE_EQ(b,0);
}

//************************intersezioni sfere**************
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

//ho aggiunto questo test per vedere che succede usando due poligoni uguali
TEST(NearFracture_test, poligoni_uguali){
    //bool NearFractures(const Fracture& frc1, const Fracture& frc2, const vector<Vector3d>& coord);
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    vector<Vector3d> coordinate;
    coordinate={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0},{0,-4,0},{4,0,0},{0,4,0},{-4,0,0}};
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    bool Risultato_funzione;
    Risultato_funzione=g.NearFractures(f1,f2,coordinate);
    ASSERT_TRUE(Risultato_funzione);

}


// FracturesFunctions::IntersectionFractures

TEST(IntersectionFractures_test, generale_sofi){
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    Trace t;
    vector<Vector3d> coordinate;
    coordinate.resize(8);
    list<Trace> list_traces; //lista delle tracce
    map<unsigned int, list<Trace>> P_traces_of_fractures; //per ogni frattura memorizziamo una lista contenente gli id elle tracce passanti
    map<unsigned int, list<Trace>> NP_traces_of_fractures; //analogo a sopra ma per tracce non passanti

    coordinate={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0},{0,0,-2},{2,0,0},{0,0,2},{-2,0,0}};
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    f1.id=0;
    f2.id=1;
    g.IntersectionFractures(f1, f2, coordinate, list_traces, P_traces_of_fractures, NP_traces_of_fractures);
    MatrixXd estremi;
    estremi.resize(3,2);
    estremi << 2, -2,
        0,0,
        0,0;

    t = (*list_traces.begin());

    MatrixXd risultato_funzione;
    risultato_funzione.resize(3,2);
    risultato_funzione=t.coordinates_extremes;

    for (unsigned int i = 0; i < risultato_funzione.cols(); i++){
        for (unsigned int j = 0; j < risultato_funzione.rows(); j++){
            EXPECT_NEAR(risultato_funzione(j,i), estremi(j,i), g.tolleranza1D);
        }
    }
    EXPECT_EQ(t.id_frc1,f1.id);
    EXPECT_EQ(t.id_frc2,f2.id);
    ASSERT_DOUBLE_EQ(t.len,4);

    //verifico che la traccia sia passante per la frattura 1 e non passante per la frattura 0
    list<unsigned int> trc_per_0_P = {};
    list<unsigned int> trc_per_0_NP = {};
    list<unsigned int> trc_per_1_P = {};
    list<unsigned int> trc_per_1_NP = {};

    for (auto it = P_traces_of_fractures.begin(); it != P_traces_of_fractures.end(); it++){
        list<Trace> lista= it->second;
        unsigned int id_frc = it->first;
        for (const Trace &trc : lista){
            if (id_frc == 0){trc_per_0_P.push_back(trc.id);}
            if (id_frc == 1){trc_per_1_P.push_back(trc.id);}
        }
    }

    for (auto it = NP_traces_of_fractures.begin(); it != NP_traces_of_fractures.end(); it++){
        list<Trace> lista= it->second;
        unsigned int id_frc = it->first;
        for (const Trace &trc : lista){
            if (id_frc == 0){trc_per_0_NP.push_back(trc.id);}
            if (id_frc == 1){trc_per_1_NP.push_back(trc.id);}
        }
    }

    EXPECT_EQ(trc_per_0_P.size(), 0); //la fratt 0 non ha tracce passanti
    EXPECT_EQ(trc_per_1_NP.size(), 0); //la fratt 1 non ha tracce non passanti

    auto it_0 = find(trc_per_0_NP.begin(), trc_per_0_NP.end(), t.id);
    auto it_1 = find(trc_per_1_P.begin(), trc_per_1_P.end(), t.id);

    bool pr0 = (it_0 == trc_per_0_NP.end()); //ritorna vero se non ho trovato la traccia
    bool pr1 = (it_1 == trc_per_1_P.end()); //ritorna vero se non ho trovato la traccia

    EXPECT_FALSE(pr0);
    EXPECT_FALSE(pr1);
}

TEST(IntersectionFractures_test, traccia_di_lunghezza_nulla){
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    double perturbazione = g.tolleranza1D*1./10;
    vector<Vector3d> coordinate;
    list<Trace> list_traces; //lista delle tracce
    map<unsigned int, list<Trace>> P_traces; //per ogni frattura memorizziamo una lista contenente gli id elle tracce passanti
    map<unsigned int, list<Trace>> NP_traces; //analogo a sopra ma per tracce non passanti

    coordinate={{0,-0.43,0},{1,-0.78,0},{1,1,0},{0,1,0},{1,0,-1},{1,0,2},{3,0,2},{3,0,-1}};
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    f1.id=0;
    f2.id=1;
    g.IntersectionFractures(f1, f2, coordinate, list_traces, P_traces, NP_traces);

    //in questo test mi aspetto di trovare una traccia di lunghezza < toll qyuindi non la dovrei considerare come tale
    unsigned int num_tracce_per_0 = P_traces[0].size() + NP_traces[0].size();
    unsigned int num_tracce_per_1 = P_traces[1].size() + NP_traces[1].size();

    EXPECT_EQ(num_tracce_per_0, 0);
    EXPECT_EQ(num_tracce_per_1, 0);
}

TEST(IntersectionFractures_test, generale_marti){
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    Trace t;
    vector<Vector3d> coordinate;
    list<Trace> list_traces; //lista delle tracce
    map<unsigned int, list<Trace>> P_traces; //per ogni frattura memorizziamo una lista contenente gli id elle tracce passanti
    map<unsigned int, list<Trace>> NP_traces; //analogo a sopra ma per tracce non passanti

    coordinate={{0,-0.43,0},{1,-0.78,0},{3.75,1,0},{0,1,0},{1,0,-1},{1,0,2},{3,0,2},{3,0,-1}};
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    f1.id=0;
    f2.id=1;
    g.IntersectionFractures(f1, f2, coordinate, list_traces, P_traces, NP_traces);

    //in questo test mi aspetto di trovare una traccia non passante per entrambe le fratture
    MatrixXd estremi;
    estremi.resize(3,2);
    estremi << 1, 2.205056179775281,
        0,0,
        0,0;

    t = (*list_traces.begin());

    MatrixXd risultato_funzione;
    risultato_funzione.resize(3,2);
    risultato_funzione=t.coordinates_extremes;

    for (unsigned int i = 0; i < risultato_funzione.cols(); i++){
        for (unsigned int j = 0; j < risultato_funzione.rows(); j++){
            EXPECT_NEAR(risultato_funzione(j,i), estremi(j,i), g.tolleranza1D);
        }
    }
    EXPECT_EQ(t.id_frc1,f1.id);
    EXPECT_EQ(t.id_frc2,f2.id);
    ASSERT_DOUBLE_EQ(t.len,1.205056179775281);

    //verifico che la traccia sia passante per la frattura 1 e non passante per la frattura 0
    list<unsigned int> trc_per_0_P = {};
    list<unsigned int> trc_per_0_NP = {};
    list<unsigned int> trc_per_1_P = {};
    list<unsigned int> trc_per_1_NP = {};

    for (auto it = P_traces.begin(); it != P_traces.end(); it++){
        list<Trace> lista= it->second;
        unsigned int id_frc = it->first;
        for (const Trace &trc : lista){
            if (id_frc == 0){trc_per_0_P.push_back(trc.id);}
            if (id_frc == 1){trc_per_1_P.push_back(trc.id);}
        }
    }

    for (auto it = NP_traces.begin(); it != NP_traces.end(); it++){
        list<Trace> lista= it->second;
        unsigned int id_frc = it->first;
        for (const Trace &trc : lista){
            if (id_frc == 0){trc_per_0_NP.push_back(trc.id);}
            if (id_frc == 1){trc_per_1_NP.push_back(trc.id);}
        }
    }

    EXPECT_EQ(trc_per_0_P.size(), 0); //la fratt 0 non ha tracce passanti
    EXPECT_EQ(trc_per_1_P.size(), 0); //la fratt 1 non ha tracce passanti

    auto it_0 = find(trc_per_0_NP.begin(), trc_per_0_NP.end(), t.id);
    auto it_1 = find(trc_per_1_NP.begin(), trc_per_1_NP.end(), t.id);

    bool pr0 = (it_0 == trc_per_0_NP.end()); //ritorna vero se non ho trovato la traccia
    bool pr1 = (it_1 == trc_per_1_NP.end()); //ritorna vero se non ho trovato la traccia

    EXPECT_FALSE(pr0);
    EXPECT_FALSE(pr1);
}

TEST(IntersectionFractures_test, traccia_non_trovata_reny)
{
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    Trace t;
    vector<Vector3d> coordinate;
    coordinate.resize(8);
    list<Trace> list_traces; //lista delle tracce
    map<unsigned int, list<Trace>> P_traces_of_fractures; //per ogni frattura memorizziamo una lista contenente gli id elle tracce passanti
    map<unsigned int, list<Trace>> NP_traces_of_fractures; //analogo a sopra ma per tracce non passanti

    coordinate={{6.2146082912229295e-01,3.6933505910269426e-01 ,6.5491065861714337e-01},{9.3634260291962368e-02,1.3967317445425522e+00,6.5491065861714337e-01},{1.1880550392475797e-01,1.4096635082402353e+00,1.4244250759175627e+00},{6.4663207275508849e-01,3.8226682280037716e-01,1.4244250759175630e+00},
                {1.0134507097882688e+00,-4.2358095677780244e-01,1.3440567176761081e-01},{-2.9565547130751324e-02,1.3939531644777750e-01,1.3440567176761081e-01},{2.3802571986642296e-01,6.3515695856546772e-01,6.8846529129663980e-01},{1.2810419767854433e+00,7.2180685339887751e-02,6.8846529129663991e-01}};
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    f1.id=0;
    f2.id=1;
    g.IntersectionFractures(f1, f2, coordinate, list_traces, P_traces_of_fractures, NP_traces_of_fractures);
    MatrixXd estremi;
    estremi.resize(3,2);
    estremi << 0.581538, 0.60718,
        0.449743,0.397132,
        0.688465,0.654911;
    double lunghezza_vera = (estremi.col(0)- estremi.col(1)).norm();   // il valore vero della lunghezza della traccia

    t = (*list_traces.begin());

    MatrixXd risultato_funzione;
    risultato_funzione.resize(3,2);
    risultato_funzione=t.coordinates_extremes;

    for (unsigned int i = 0; i < risultato_funzione.cols(); i++){
        for (unsigned int j = 0; j < risultato_funzione.rows(); j++){
            EXPECT_NEAR(risultato_funzione(j,i), estremi(j,i), g.tolleranza1D);
        }
    }
    EXPECT_EQ(t.id_frc1,f1.id);
    EXPECT_EQ(t.id_frc2,f2.id);
    ASSERT_DOUBLE_EQ(t.len,lunghezza_vera);

    //verifico che la traccia sia passante per la frattura 1 e non passante per la frattura 0
    list<unsigned int> trc_per_0_P = {};
    list<unsigned int> trc_per_0_NP = {};
    list<unsigned int> trc_per_1_P = {};
    list<unsigned int> trc_per_1_NP = {};

    for (auto it = P_traces_of_fractures.begin(); it != P_traces_of_fractures.end(); it++){
        list<Trace> lista= it->second;
        unsigned int id_frc = it->first;
        for (const Trace &trc : lista){
            if (id_frc == 0){trc_per_0_P.push_back(trc.id);}
            if (id_frc == 1){trc_per_1_P.push_back(trc.id);}
        }
    }

    for (auto it = NP_traces_of_fractures.begin(); it != NP_traces_of_fractures.end(); it++){
        list<Trace> lista= it->second;
        unsigned int id_frc = it->first;
        for (const Trace &trc : lista){
            if (id_frc == 0){trc_per_0_NP.push_back(trc.id);}
            if (id_frc == 1){trc_per_1_NP.push_back(trc.id);}
        }
    }

    EXPECT_EQ(trc_per_0_P.size(), 0); //la fratt 0 non ha tracce passanti
    EXPECT_EQ(trc_per_1_P.size(), 0); //la fratt 1 non ha tracce passanti

    auto it_0 = find(trc_per_0_NP.begin(), trc_per_0_NP.end(), t.id);
    auto it_1 = find(trc_per_1_P.begin(), trc_per_1_P.end(), t.id);

    bool pr0 = (it_0 == trc_per_0_NP.end()); //ritorna vero se non ho trovato la traccia
    bool pr1 = (it_1 == trc_per_1_P.end()); //ritorna vero se non ho trovato la traccia

    EXPECT_FALSE(pr0);
    EXPECT_FALSE(pr1);
}
#endif
