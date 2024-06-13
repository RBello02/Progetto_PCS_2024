#ifndef __TEST_FRACTURES_H // Header guards
#define __TEST_FRACTURES_H

#include "Eigen/Eigen"
#include <gtest/gtest.h>
#include <gmock/gmock.h>
#include <math.h>
#include <iostream>
#include <fstream>
#include <string>
#include "FracturesLibrary.hpp"
#include "Utils.hpp"



using namespace std;
using namespace Eigen;
using namespace FracturesLibrary;
using namespace UtilsFunction;

inline bool compare_tracce(const Trace& trc1, const Trace& trc2){
    return trc1.len > trc2.len;
}


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
    Risultato_funzione=g.vicinanza_sfere(f1,f2,coordinate);
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
    Risultato_funzione=g.vicinanza_sfere(f1,f2,coordinate);
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
    Risultato_funzione=g.vicinanza_sfere(f1,f2,coordinate);
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
    Risultato_funzione=g.vicinanza_sfere(f1,f2,coordinate);
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
    g.ricerca_tracce(f1, f2, coordinate, list_traces, P_traces_of_fractures, NP_traces_of_fractures);
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
    //tecnicamente tra le due fratture c'è una traccia di lunghezza l tc eps_macch < l < tolleranza1D, verifico che non venga etta come traccia
    Fracture f1;
    Fracture f2;
    FracturesFunctions g;
    double perturbazione = g.tolleranza1D - 0.5*(g.tolleranza1D - g.eps_macchina);
    vector<Vector3d> coordinate;
    list<Trace> list_traces; //lista delle tracce
    map<unsigned int, list<Trace>> P_traces; //per ogni frattura memorizziamo una lista contenente gli id elle tracce passanti
    map<unsigned int, list<Trace>> NP_traces; //analogo a sopra ma per tracce non passanti

    coordinate={{0,-0.43,0},{1+perturbazione,-0.78,0},{1+perturbazione,1,0},{0,1,0},{1,0,-1},{1,0,2},{3,0,2},{3,0,-1}};
    f1.num_vertici=4;
    f2.num_vertici=4;
    f1.vertices={0,1,2,3};
    f2.vertices={4,5,6,7};
    f1.id=0;
    f2.id=1;
    g.ricerca_tracce(f1, f2, coordinate, list_traces, P_traces, NP_traces);

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
    g.ricerca_tracce(f1, f2, coordinate, list_traces, P_traces, NP_traces);

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


TEST(Tuttoilprogramma_test, Generale){

    cout << setprecision(16);
    FracturesFunctions g;
    Trace t;
    string path="TEST_data.txt";
    string file1="tracceTEST.txt";
    string file2="tracce_per_frattura_TEST.txt";
    vector<Fracture> list_fractures;
    list<Trace> list_traces;
    map<unsigned int, list<Trace>> P_traces_of_fractures; //per ogni frattura memorizziamo una lista contenente gli id elle tracce passanti
    map<unsigned int, list<Trace>> NP_traces_of_fractures; //analogo a sopra ma per tracce non passanti
    unsigned int num_fratt;
    vector<Vector3d> coordinates;
    bool l;
    l=g.lettura_da_file(path, list_fractures, coordinates);
    EXPECT_TRUE(l);
    num_fratt = list_fractures.size();
    EXPECT_EQ(num_fratt,2);
    EXPECT_EQ(list_fractures[0].id,0);
    EXPECT_EQ(list_fractures[1].id,1);

    EXPECT_EQ(list_fractures[0].num_vertici,4);
    EXPECT_EQ(list_fractures[1].num_vertici,4);

    vector<Vector3d> CoordinateAttese;
    CoordinateAttese={{0,-4,0},{4,0,0},{0,4,0},{-4,0,0},{0,0,-2},{2,0,0},{0,0,2},{-2,0,0}};
    EXPECT_EQ(coordinates,CoordinateAttese);


    for(unsigned int i=0; i<num_fratt; i++)
    {
        for (unsigned int j=i+1; j< num_fratt; j++)
        {
            Fracture frc1 = list_fractures[i];
            Fracture frc2 = list_fractures[j];

            if( g.vicinanza_sfere(frc1, frc2, coordinates)){
                g.ricerca_tracce(frc1, frc2, coordinates, list_traces, P_traces_of_fractures, NP_traces_of_fractures);
            }

        }
    }

    //primo file di output
    string filenameO_tracce = "tracce.txt";
    ofstream ofs;
    ofs.open(filenameO_tracce);



    ofs << "# Number of Traces" << endl << list_traces.size() << endl;
    ofs << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    const string del = " ; ";

    for (Trace& traccia : list_traces){
        ofs << traccia.id << del << traccia.id_frc1 << del << traccia.id_frc2 << del << traccia.coordinates_extremes(0,0) << del << traccia.coordinates_extremes(1,0) << del << traccia.coordinates_extremes(2,0) << del << traccia.coordinates_extremes(0,1) << del << traccia.coordinates_extremes(1,1) << del << traccia.coordinates_extremes(2,1) << endl;
    }

    ofs.close();


    //secondo file di output
    string filenameO_frc = "tracce_per_frattura.txt";
    ofstream ofs1;
    ofs1.open(filenameO_frc);


    string header_frc = "# FractureId; NumTraces";
    string header_trc = "# TraceId; Tips; Lenght";

    for (Fracture& fratt : list_fractures){

        ofs1 << header_frc << endl;
        ofs1 << fratt.id << del << P_traces_of_fractures[fratt.id].size() + NP_traces_of_fractures[fratt.id].size() << endl;

        if (P_traces_of_fractures[fratt.id].size() + NP_traces_of_fractures[fratt.id].size() != 0){
            ofs1 << header_trc << endl;

            //prima ordiniamo le liste
            (P_traces_of_fractures[fratt.id]).sort(compare_tracce);
            (NP_traces_of_fractures[fratt.id]).sort(compare_tracce);

            //scorriamo la lista delle tracce passanti
            for (Trace& traccia : P_traces_of_fractures[fratt.id]){
                ofs1 << traccia.id << del << "false" << del << traccia.len << endl;
            }

            //scorriamo la lista delle tracce non passanti
            for (Trace& traccia : NP_traces_of_fractures[fratt.id]){
                ofs1 << traccia.id << del << "true"  << del << traccia.len << endl;
            }
        }
        ofs1 << endl;
    }

    ofs1.close();

    ifstream ifs_test(file1);
    ifstream ifs_prova(filenameO_tracce);

    //verifico da averli aperti correttamente
    if(ifs_test.fail()){cerr << "Errore nell'apertura del file di test" << endl; return;}
    if(ifs_prova.fail()){cerr << "Errore nell'apertura del file di prova" << endl; return;}

    string linetest;
    string lineatt;
    //ciclo finché entrambi i file di output sono ancora aperti
    while ( (!ifs_prova.eof())|| (!ifs_test.eof())) {
        getline(ifs_test, linetest);
        getline(ifs_prova, lineatt);
        EXPECT_EQ(linetest,lineatt);
    }

    ifs_test.close();
    ifs_prova.close();


    ifstream ifs_test2(file2);
    ifstream ifs_prova2(filenameO_frc);
    //verifico da averli aperti correttamente
    if(ifs_test2.fail()){cerr << "Errore nell'apertura del file di test2" << endl; return;}
    if(ifs_prova2.fail()){cerr << "Errore nell'apertura del file di prova2" << endl; return;}

    //ciclo finché entrambi i file di output sono ancora aperti
    while ( (!ifs_prova2.eof())|| (!ifs_test2.eof())) {
        getline(ifs_test2, linetest);
        getline(ifs_prova2, lineatt);
        EXPECT_EQ(linetest,lineatt);

    }

    ifs_test2.close();
    ifs_prova2.close();


    //verifica mesh 1 (sofi scrivi sopra please)
    PolygonalMesh mesh1 = g.creazione_mesh(list_fractures[0], P_traces_of_fractures[0], NP_traces_of_fractures[0], coordinates);

    EXPECT_EQ(mesh1.NumberCell0D,4);
    EXPECT_EQ(mesh1.NumberCell1D, 5);
    EXPECT_EQ(mesh1.NumberCell2D, 2);

    for (unsigned int i =0; i < mesh1.NumberCell0D; i++){
        //verifica delle celle 0D
        for (unsigned int k = 0; k < 3; k++){
            EXPECT_NEAR(mesh1.Cell0DCoordinates[i][k], coordinates[i][k], g.tolleranza1D);
        }
    }


    //verifica delle celle 1D
    unsigned int cont = 0;
    Vector2i expect = {0,1};
    EXPECT_EQ(mesh1.Cell1DVertices[cont], expect);
    expect = {1,3};
    cont ++;
    EXPECT_EQ(mesh1.Cell1DVertices[cont], expect);
    expect = {3,0};
    cont ++;
    EXPECT_EQ(mesh1.Cell1DVertices[cont], expect);
    expect = {1,2};
    cont++;
    EXPECT_EQ(mesh1.Cell1DVertices[cont], expect);
    expect = {2,3};
    cont ++;
    EXPECT_EQ(mesh1.Cell1DVertices[cont], expect);


    vector vert_frc1 = mesh1.Cell2DVertices[0]; //verifica dei vertici della sotto fratt 1 della fratt 0
    EXPECT_EQ(vert_frc1[0], 0);
    EXPECT_EQ(vert_frc1[1], 1);
    EXPECT_EQ(vert_frc1[2], 3);

    vector vert_frc2 = mesh1.Cell2DVertices[1]; //verifica dei vertici della sotto fratt 2 della fratt 0
    EXPECT_EQ(vert_frc2[0], 1);
    EXPECT_EQ(vert_frc2[1], 2);
    EXPECT_EQ(vert_frc2[2], 3);

    vector edges_frc1 =mesh1.Cell2DEdges[0];
    EXPECT_EQ(edges_frc1[0], 0);
    EXPECT_EQ(edges_frc1[1], 1);
    EXPECT_EQ(edges_frc1[2], 2);

    vector edges_frc2 =mesh1.Cell2DEdges[1];
    EXPECT_EQ(edges_frc2[0], 3);
    EXPECT_EQ(edges_frc2[1], 4);
    EXPECT_EQ(edges_frc2[2], 1);



    //verifica mesh 1
    PolygonalMesh mesh2 = g.creazione_mesh(list_fractures[1], P_traces_of_fractures[1], NP_traces_of_fractures[1], coordinates);

    EXPECT_EQ(mesh2.NumberCell0D,4);
    EXPECT_EQ(mesh2.NumberCell1D, 5);
    EXPECT_EQ(mesh2.NumberCell2D, 2);

    for (unsigned int i =0; i < mesh2.NumberCell0D; i++){
        //verifica delle celle 0D
        for (unsigned int k = 0; k < 3; k++){
            EXPECT_NEAR(mesh2.Cell0DCoordinates[i][k], coordinates[4+i][k], g.tolleranza1D);
        }
    }


    //verifica delle celle 1D
    cont = 0;
    expect = {0,1};
    EXPECT_EQ(mesh2.Cell1DVertices[cont], expect);
    expect = {1,3};
    cont ++;
    EXPECT_EQ(mesh2.Cell1DVertices[cont], expect);
    expect = {3,0};
    cont ++;
    EXPECT_EQ(mesh2.Cell1DVertices[cont], expect);
    expect = {1,2};
    cont++;
    EXPECT_EQ(mesh2.Cell1DVertices[cont], expect);
    expect = {2,3};
    cont ++;
    EXPECT_EQ(mesh2.Cell1DVertices[cont], expect);



    vert_frc1 = mesh2.Cell2DVertices[0]; //verifica dei vertici della sotto fratt 1 della fratt 0
    EXPECT_EQ(vert_frc1[0], 0);
    EXPECT_EQ(vert_frc1[1], 1);
    EXPECT_EQ(vert_frc1[2], 3);

    vert_frc2 = mesh2.Cell2DVertices[1]; //verifica dei vertici della sotto fratt 2 della fratt 0
    EXPECT_EQ(vert_frc2[0], 1);
    EXPECT_EQ(vert_frc2[1], 2);
    EXPECT_EQ(vert_frc2[2], 3);

    edges_frc1 =mesh2.Cell2DEdges[0];
    EXPECT_EQ(edges_frc1[0], 0);
    EXPECT_EQ(edges_frc1[1], 1);
    EXPECT_EQ(edges_frc1[2], 2);

    edges_frc2 =mesh2.Cell2DEdges[1];
    EXPECT_EQ(edges_frc2[0], 3);
    EXPECT_EQ(edges_frc2[1], 4);
    EXPECT_EQ(edges_frc2[2], 1);



}

TEST(ascissa_curvilinea, Generale_Renato_1)
{
    FracturesFunctions g;
    double toll = g.tolleranza1D;
    MatrixXd retta;
    retta.resize(2,3);
    retta << 1,1,1,
             0,0,0;
    Vector3d x = {5,5,5};
    EXPECT_NEAR(g.ascissa_curvilinea(retta,x), 5, toll);
}

TEST(ascissa_curvilinea, Generale_Renato_2)
{
    FracturesFunctions g;
    double toll = g.tolleranza1D;
    MatrixXd retta;
    retta.resize(2,3);
    retta << 1,1,1,
        0,0,0;
    Vector3d x = {5,0,5};
    EXPECT_TRUE(isnan(g.ascissa_curvilinea(retta,x)));
}
#endif
