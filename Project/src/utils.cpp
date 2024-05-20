#include "FracturesLibrary.hpp"
#include "PolygonalMesh.hpp"

using namespace std;
using namespace  Eigen;
using namespace PolygonalLibrary;


namespace FracturesLibrary{

inline MatrixXd FracturesFunctions::Retta_tra_piani(const Matrix3d& piano_1, const Matrix3d& piano_2)
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

inline MatrixXd FracturesFunctions::Retta_per_due_vertici_della_frattura(unsigned int id_vertice1, unsigned int id_vertice2, const vector<Vector3d>& coord)
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

inline Vector2d FracturesFunctions::alpha_di_intersezione(MatrixXd r_frattura, MatrixXd retta_intersez)
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

inline MatrixXd FracturesFunctions::Retta_per_due_punti(Vector3d& pt1, Vector3d& pt2)
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

inline Vector2d FracturesFunctions::intersezione_rette(MatrixXd& r_frattura, MatrixXd& r_traccia)
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

inline bool FracturesFunctions::appartiene_a_segmento(Vector3d& origin, Vector3d& end, Vector3d& pto){
    bool appartiene = false;

    //calcolo la distanza di pto dai due estremi del segmento, se la somma di queste due distanza è maggiore della distanza dei
    // due estremi allora il pto è esterno
    double lung_segmento = (origin-end).norm();
    double pto_o = (origin-pto).norm();
    double pto_e = (end-pto).norm();

    if(abs(pto_o + pto_e - lung_segmento) < tolleranza1D){appartiene = true;}

    return appartiene;

}

inline bool FracturesFunctions::pto_unico(Vector3d& pto, vector<Vector3d>& punti, unsigned int& id){
    //l'id mi serve quando cerco un punto tra i vertici

    bool unico = true;

    if(punti.size() == 0){return unico;}

    for (unsigned int i = 0; i < punti.size(); i++){
        Vector3d elem = punti[i];
        bool uguagl_x = (abs(pto[0] - elem [0]) < tolleranza1D);
        bool uguagl_y = (abs(pto[1] - elem [1]) < tolleranza1D);
        bool uguagl_z = (abs(pto[2] - elem [2]) < tolleranza1D);

        if (uguagl_x && uguagl_y && uguagl_z){unico = false; id = i; return unico;}
    }
    return unico;
}

}//namespace