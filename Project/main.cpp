#include "src/GeometryLibrary.hpp"
#include "src/Utils.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include "Eigen/Eigen"
#include <cmath>

using namespace GeometryLibrary;
using namespace Eigen;
using namespace std;


inline Matrix3d Coef_piano(const Fractures& frc, unsigned int id_fract1)
{
    // salvo il piano in forma parametrica X = P0 + a1(P2-P0) + a2(P1-P0), dove a1 e a2 sono le parametrizzazioni
    // il piano è quindi identificato da 3 punti P0,P1,P2
    Matrix3d A;
    for (unsigned int i=0; i<3; i++) // ciclo su 3 punti del poligono
    {
        unsigned int id_vertice = frc.vertices_fractures[id_fract1][i];
        for (unsigned int j = 0; j< 3; j++)
        {
            A(i,j) = frc.coordinates[id_vertice][j];
            //cout << A(i,j)<< " ";
        }
        //cout << endl;
    }
    return A; // restituisce una matrice con la seguente struttura   (P0;P1;P2)  dove Pi sono VETTORI RIGA
}

inline bool Parallelismo(const Matrix3d& piano_1, const Matrix3d& piano_2)
{
    // questa funzione verifica se i piani sono paralleli

    bool risultato = false;

    // NOMENCLATURA
    // u1 = piano_1.row(1) v1 = piano_1.row(2) sono i vettori direttori del primo piano
    // u2 = piano_2.row(1) v2 = piano_2.row(2) sono i vettori direttori del secondo piano
    // P0_1 = piano_1.row(0) punto iniziale per il piano 1
    // P0_2 = piano_2.row(0) punto iniziale per il piano 2

    // passo 1: ricavare la normale ai piani
    Vector3d n1 = piano_1.row(1).cross(piano_1.row(2));
    Vector3d n2 = piano_2.row(1).cross(piano_2.row(2));

    // passo 2 : vedo se il prodotto vettoriale tra le normali è 0
    Vector3d t = n1.cross(n2);
    if (t.norm() == 0)
    {
        risultato = !risultato;
    }
    return risultato;
}

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
    A << n1,n2,t;

    // passo 5: risolvere il sistema per ricavare P
    Vector3d P = A.partialPivLu().solve(d);

    // salvo in un formato particolare
    MatrixXd X;
    X.resize(2,3);
    X.row(0) = t.transpose();
    X.row(1) = P.transpose();  // come una matrice 2x3

    return X;

}

inline MatrixXd Retta_per_due_vertici_della_frattura(const Fractures& frc, unsigned int id_vertice1, unsigned int id_vertice2)
{
    // l'equazione parametrica è X = at+P

    // passo 1: salvo le coordinate dei due punti
    double x1 = frc.coordinates[id_vertice1][0];
    double y1 = frc.coordinates[id_vertice1][1];
    double z1 = frc.coordinates[id_vertice1][2];

    double x2 = frc.coordinates[id_vertice2][0];
    double y2 = frc.coordinates[id_vertice2][1];
    double z2 = frc.coordinates[id_vertice2][2];

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
inline double alpha_di_intersezione(MatrixXd A, MatrixXd B)
{
    // PRIMA DI APPLICARE QUESTA FUNZIONE DEVO FARE UN CONTROLLO PER VEDERE CHE t' != t
    // dati i due vettori contenenti i dati delle rette vogliamo trovarne l'intersezione
    // dati X = at + P e X = at'+P', imponiamo il sistema,
    // ricavo il valore di a per cui si intersecano, basta farlo rispetto ad una componente dei vettori e lo faccio a mano
    // tra le componenti disponibili cerco quella tale per cui t-t' non è nulla

    Vector3d t1 = A.row(0);
    Vector3d P1 = A.row(1);

    Vector3d t2 = B.row(0);
    Vector3d P2 = B.row(1);

    double alpha = 0;   // inizializzo il valore di alpha

    // questi sono tutti i casi al fine di non avere alpha subito zero

    if (t1(0) == t2(0))
    {
        if (t1(1) == t2(1))
        {
            alpha = (P2(2) - P1(2))/(t2(2)-t1(2));
        }
        else
        {
            if (P2(1) == P2(1))
            {
                if (t1(2) == t2(2))
                {
                    alpha = (P2(1) - P1(1))/(t2(1)-t1(1));
                }
                else
                {
                    alpha = (P2(2) - P1(2))/(t2(2)-t1(2));
                }
            }
            else
            {
                alpha = (P2(1) - P1(1))/(t2(1)-t1(1));
            }

        }
    }
    else
    {
        if (P1(0) == P2(0))
        {
            if (t1(1) == t2(1))
            {
                if (t1(2) == t2(2))
                {
                    alpha = (P2(0) - P1(0))/(t2(0)-t1(0));
                }
                else
                {
                alpha = (P2(2) - P1(2))/(t2(2)-t1(2));
                }
            }
            else
            {
                if(P2(1) == P1(2))
                {
                    if (t1(2) == t2(2))
                    {
                        alpha = (P2(1) - P1(1))/(t2(1)-t1(1));
                    }
                    else
                    {
                        alpha = (P2(2) - P1(2))/(t2(2)-t1(2));
                    }
                }
                else
                {
                    alpha = (P2(1) - P1(1))/(t2(1)-t1(1));
                }

            }
        }
        else
        {
            alpha = (P2(0) - P1(0))/(t2(0)-t1(0));
        }
    }

    return alpha;

}


int main()
{
    cout << setprecision(16);
    string path = "DFN";
    string filename = path + "/FR3_data.txt";
    Fractures frc;

    if (!importData(filename, frc)){return 1;}
    else{
        //stampo un poo' di roba per verificare che sia tutto giusto
        cout << "Ho un numero di fratture pari a " << frc.num_fractures << endl;

        for (unsigned int i = 0; i < frc.num_fractures ; i++){
            cout << "FRATTURA NUM " << i <<  " che ha " << frc.dim_fractures[i] << " vertici." << endl;

            for (unsigned int k = 0; k < frc.dim_fractures[i]; k++){

                 unsigned int id_vertice = frc.vertices_fractures[i][k];
                cout << "Il vertice " << k << ", che ha id = " << id_vertice << " , ha coord (" << frc.coordinates[id_vertice][0] << " , "  <<frc.coordinates[id_vertice][1] <<  " , " <<frc.coordinates[id_vertice][2] << " )" << endl;
             }
            cout << endl;
        }


    }

    // ciclo sulle coppie di poligoni e vedo se sono vicini oppure no
    for(unsigned int i=0; i<frc.num_fractures; i++)
    {
        for (unsigned int j=i+1; j< frc.num_fractures; j++)
        {
            unsigned int id_fract1=i;
            unsigned int id_fract2=j;
            if( !NearFractures(frc, id_fract1,id_fract2))
            {
                cout << "non dobbiamo fare niente"<< endl;
            }
            else
            {
                cout<<"fare i piani e le rette" << endl;

            }

        }
    }

    // esempi per vedere se funziona tutto

    unsigned int id_fract1 = 0;
    unsigned int id_fract2 = 2;
    Matrix3d a = Coef_piano(frc,id_fract1);
    Matrix3d b = Coef_piano(frc,id_fract2);
    //cout << Parallelismo(a,b)<<endl;
    MatrixXd X;
    X.resize(2,3);
    X = Retta_tra_piani(a,b);
    unsigned int id_vertice1 = frc.vertices_fractures[id_fract1][0];
    unsigned int id_vertice2 = frc.vertices_fractures[id_fract1][1];
    MatrixXd X0;
    X0.resize(2,3);
    X0 = Retta_per_due_vertici_della_frattura(frc,id_vertice1,id_vertice2);
    //cout << X << endl;
    double A = alpha_di_intersezione(X,X0);
    cout << X << endl;
    cout << X0 << endl;
    cout << A << endl;



    return 0;
}
