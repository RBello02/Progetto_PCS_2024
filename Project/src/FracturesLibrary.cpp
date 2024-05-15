#include "FracturesLibrary.hpp"
#include "reshaping_array.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>
#include <algorithm>

#include <cmath>

#include "Eigen/Eigen"

using namespace Eigen;

namespace FracturesLibrary{
double tolleranza = pow(10,-10);

//funzioni inline di supporto
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

inline bool FracturesFunctions::Parallelismo(const Matrix3d& piano_1, const Matrix3d& piano_2)
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
    if (t.norm() < tolleranza1D)
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
    //primo parametro è la matrice della retta del poligono --> retta in funzione di alpha
    Vector3d t1 = r_frattura.row(0).transpose();
    Vector3d P1 = r_frattura.row(1).transpose();

    //secondo parametro è la matrice della retta di intersezione tra i piani --> retta in funzione di beta
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

inline bool compare_beta(const array<double,2>& arr1, const array<double,2>& arr2){
    return arr1[1] > arr2[1];
}

inline double Trace::calcolo_lunghezza(){
    Vector3d origin = this->coordinates_extremes.col(0);
    Vector3d end = this->coordinates_extremes.col(1);
    double val = (origin-end).norm();
    return val;
}
/****************************************************************************************************************/

bool FracturesFunctions::importData(const string& path, vector<Fracture>& lista, vector<Vector3d>& coord){

    //apro il file
    ifstream file;
    file.open(path);

    if(file.fail()){
        cerr << "errore nell'apertura del file" << endl;
        return false;}

    //la prima riga è di header
    string header1;
    getline(file, header1);

    // mi salvo il numero di fratture
    unsigned int num_fratture;
    string line;
    getline(file, line);
    istringstream convert(line);
    convert >> num_fratture;
    lista.reserve(num_fratture);

    //riservo al vettore di coordinate uno spazio superiore a quello che credo che mi servirà
    coord.reserve(num_fratture*7);

    unsigned int cont_vertici = 0;

    //sfrutto la correttezza del file per leggermi le righe a gruppi di 6
    for(unsigned int i = 0; i < num_fratture; i++){
        Fracture frc;

        //la prima riga è di header
        string header1;
        getline(file, header1);

        //la seconda contiene id e numero vertici
        string id_num;
        getline(file, id_num);
        istringstream convert1(id_num);
        char del;
        convert1 >> frc.id >> del >> frc.num_vertici;
        frc.vertices.reserve(frc.num_vertici);

        //riga di intestazione
        string header2;
        getline(file,header2);

        //adesso ho 3 righe che contengono i vertici
        MatrixXd vert = MatrixXd::Zero(3, frc.num_vertici);
        string line1;
        for (unsigned int k = 0; k<3; k++){
            getline(file,line1);
            istringstream convert2(line1);
            for (unsigned int l = 0; l < frc.num_vertici -1; l++){
                char del;
                convert2 >> vert(k,l) >> del;
            }
            convert2 >> vert(k,frc.num_vertici-1);
        }

        //salvo ora le coordinate nel vettore coord e i rispettivi id nella struct Fracture
        for (unsigned int k = 0; k < frc.num_vertici; k++){
            Vector3d v = vert.col(k);
            coord.push_back(v);
            ReshapingArray::VerificaRaddoppio(coord);
            frc.vertices.push_back(cont_vertici);
            cont_vertici += 1;
        }

        lista.push_back(frc);

    }

    coord.shrink_to_fit(); //elimino la capacità in eccesso
    file.close();
    return true;
}


bool FracturesFunctions::NearFractures(const Fracture& frc1, const Fracture& frc2, const vector<Vector3d>& coord){

    // vettori per le coordinate dei due baricentri (approssimativamente):
    array<double, 3> bar1;
    array<double,3> bar2;

    //calcolo le coordinate facendo somma/numvertici per ogni coordinata
    double sommax=0,sommay=0,sommaz=0;

    for (unsigned int k = 0; k < frc1.num_vertici; k++){    // ciclo sui vertici della frattura e sommo tutte le coordinate

        unsigned int id_vertice = frc1.vertices[k];
        sommax += coord[id_vertice][0];
        sommay += coord[id_vertice][1];
        sommaz += coord[id_vertice][2];

    };
    bar1[0]=sommax/frc1.num_vertici;                // calcolo la coordinata del baricentro dividendo la somma delle coordinate per il numero di vertici
    bar1[1]=sommay/frc1.num_vertici;
    bar1[2]=sommaz/frc1.num_vertici;


    for (unsigned int k = 0; k < frc2.num_vertici; k++){    // ripeto il calcolo anche per la seconda frattura

        unsigned int id_vertice = frc2.vertices[k];
        sommax += coord[id_vertice][0];
        sommay += coord[id_vertice][1];
        sommaz += coord[id_vertice][2] ;
    };
    bar2[0]=sommax/frc2.num_vertici;
    bar2[1]=sommay/frc2.num_vertici;
    bar2[2]=sommaz/frc2.num_vertici;

    //calcolo i raggi delle sfere(al quadrato) e la distanza tra i due baricentri
    double raggio1=0,raggio2=0;

    // cerco il vertice con distanza massima dal baricentro, in realtà mi interessa solo la distanza

    double raggio_da_confrontare_1 = 0;

    for (unsigned int j = 0; j<frc1.num_vertici; j++)   // ciclo su tutti i poligoni per calcolare la distanza massima dei vertici dal baricentro
    {
        unsigned int id_vertice = frc1.vertices[j];  // inizializzo l'id del vertice
        for (unsigned int i =0; i<3;i++)        // ciclo sulle 3 coordinate per calcolare il raggio del poligono
        {
            raggio1 += pow((bar1[i]-coord[id_vertice][i]),2);
        }
        if (raggio1 >= raggio_da_confrontare_1 -tolleranza1D)
        {
            raggio_da_confrontare_1 = raggio1;
        }
    }

    double raggio_da_confrontare_2 = 0;

    for (unsigned int j = 0; j<frc2.num_vertici; j++)   // ciclo su tutti i poligoni per calcolare la distanza massima dei vertici dal baricentro
    {
        unsigned int id_vertice = frc2.vertices[j];  // inizializzo l'id del vertice
        for (unsigned int i =0; i<3;i++)        // ciclo sulle 3 coordinate per calcolare il raggio del poligono
        {
            raggio2 += pow((bar1[i]-coord[id_vertice][i]),2);
        }
        if (raggio2 >= raggio_da_confrontare_2 - tolleranza1D)
        {
            raggio_da_confrontare_2 = raggio2;
        }
    }

    double distbb = 0;

    for (unsigned int j = 0; j<3; j++)
    {
        distbb += pow(bar2[j]-bar1[j],2);
    }


    if (raggio_da_confrontare_1+raggio_da_confrontare_2 < distbb + tolleranza1D)
        return false;

    else
        return true;

}

void FracturesFunctions::IntersectionFractures(Fracture &frc1, Fracture &frc2, const vector<Vector3d>& coord, list<Trace>& list_traces, map<unsigned int, list<Trace> > &P_traces, map<unsigned int, list<Trace> > &NP_traces){
    Matrix3d piano_frc1 = frc1.calcolo_piano(coord);
    Matrix3d piano_frc2 = frc2.calcolo_piano(coord);

    //procedo nella funzione solo nel caso in cui i due piani non siano paralleli
    if(Parallelismo(piano_frc1, piano_frc2)){
        cout << "La frattura " << frc1.id << " e la fratt " << frc2.id << " non si intersecano" << endl;}
    else{

        //calcolo la retta di intersezione tra i piani in forma parametrica
        MatrixXd retta_intersez_piani = Retta_tra_piani(piano_frc1,piano_frc2);
        Vector3d dir_retta_intersez_piani;
        dir_retta_intersez_piani = retta_intersez_piani.row(0);

        //Devo ora ciclare sulle coppie di vertici delle due fratture per trovare gli alpha di intersezione tra la retta di
        //intersezione tra i due piani determinati dalle fratture e le rette determinate dai vertici
        list<array<double,2>> beta_inters; //--> nel primo array avrò gli alpha della prima frc e nel secondo della seconda

        //Introduco un contatore che incremento ogni volta che la retta di intersezione tra i due piani interseca un
        // segnemnto della frattura. se alla fine avrò cont == 4, vuol dire che ho una traccia
        unsigned int cont = 0;

        for (unsigned int i = 0; i < frc1.num_vertici; i++){

            unsigned int j;
            if (i == frc1.num_vertici-1){j = 0;}
            else{j = i+1;}

            //mi trovo la retta per ogni coppia di vertici
            MatrixXd retta_tra_vertici = Retta_per_due_vertici_della_frattura(frc1.vertices[i],frc1.vertices[j], coord);
            Vector3d dir_retta_tra_vertici;
            dir_retta_tra_vertici = retta_tra_vertici.row(0);

            //cerco l'eventuale intersezione tra questa retta e la retta di intrsezione tra i piani, queste esiste se le due rette non
            //sono parallele e complanari
            bool non_parallele = !((dir_retta_tra_vertici.cross(dir_retta_intersez_piani)).norm() < tolleranza1D);

            if (non_parallele){
                Vector2d a_b = alpha_di_intersezione(retta_tra_vertici, retta_intersez_piani);

                if (a_b[0] >= -tolleranza1D && a_b[0] <= 1+tolleranza1D){
                    cont += 1;
                    beta_inters.push_back({static_cast<double>(frc1.id),a_b[1]});
                }

            }
        }

        for (unsigned int i = 0; i < frc2.num_vertici; i++){

            unsigned int j;
            if (i == frc2.num_vertici-1){j = 0;}
            else{j = i+1;}

            //mi trovo la retta per ogni coppia di vertici
            MatrixXd retta_tra_vertici = Retta_per_due_vertici_della_frattura(frc2.vertices[i],frc2.vertices[j], coord);
            Vector3d dir_retta_tra_vertici;
            dir_retta_tra_vertici = retta_tra_vertici.row(0);

            //cerco l'eventuale intersezione tra questa retta e la retta di intrsezione tra i piani, queste esiste se le due rette non
            //sono parallele e complanari
            bool non_parallele = !((dir_retta_tra_vertici.cross(dir_retta_intersez_piani)).norm() < tolleranza1D);

            if (non_parallele){
                Vector2d a_b = alpha_di_intersezione(retta_tra_vertici, retta_intersez_piani);
                if (a_b[0] >= - tolleranza1D && a_b[0] <= 1+ tolleranza1D){
                    cont += 1;
                    beta_inters.push_back({static_cast<double>(frc2.id),a_b[1]});
                }
                }
            }


        //se trovo una traccia la salvo
        if (cont == 4){
            cout << "Ho una traccia tra la frattura " << frc1.id << " e la fratt " << frc2.id << endl;
            //posso ordinare gli array tra loro in base al valore di beta
            beta_inters.sort(compare_beta);


            cout << endl;

            Trace traccia;
            traccia.id = list_traces.size();
            MatrixXd coord_estremi_traccia = MatrixXd::Zero(3,2);
            auto it_beta = beta_inters.begin();
            double beta0 = (*it_beta)[1];
            unsigned int id_frc_beta0 = (*it_beta)[0];
            ++ it_beta;
            Vector3d origin = retta_intersez_piani.row(1).transpose() + (*it_beta)[1] * dir_retta_intersez_piani;
            double beta1 = (*it_beta)[1];
            unsigned int id_frc_beta1 = (*it_beta)[0];
            ++ it_beta;
            Vector3d end = retta_intersez_piani.row(1).transpose() + (*it_beta)[1] * dir_retta_intersez_piani;
            double beta2 = (*it_beta)[1];
            unsigned int id_frc_beta2 = (*it_beta)[0];
            ++ it_beta;
            double beta3 = (*it_beta)[1];
            unsigned int id_frc_beta3 = (*it_beta)[0];

            coord_estremi_traccia.col(0) = origin;
            coord_estremi_traccia.col(1) = end;
            traccia.coordinates_extremes = coord_estremi_traccia;
            traccia.id_frc1 = frc1.id;
            traccia.id_frc2 = frc2.id;

            traccia.len = traccia.calcolo_lunghezza();
            list_traces.push_back(traccia);


            //determiniamo se la traccia è passante o no
            if (abs(beta0 -beta1) < tolleranza1D && abs(beta2-beta3) < tolleranza1D){
                //la traccia è passante per entrambe
                auto ret1 = P_traces.insert({frc1.id, {traccia}});
                if(!ret1.second)
                    (ret1.first)->second.push_back(traccia);

                auto ret2 = P_traces.insert({frc2.id, {traccia}});
                if(!ret2.second)
                    (ret2.first)->second.push_back(traccia);
            }
            else if(abs(beta0 -beta1) < tolleranza1D){
                auto ret1 = P_traces.insert({id_frc_beta2, {traccia}});
                if(!ret1.second)
                    (ret1.first)->second.push_back(traccia);

                auto ret2 = NP_traces.insert({id_frc_beta3, {traccia}});
                if(!ret2.second)
                    (ret2.first)->second.push_back(traccia);
            }
            else if(abs(beta2-beta3) < tolleranza1D){
                auto ret1 = P_traces.insert({id_frc_beta1, {traccia}});
                if(!ret1.second)
                    (ret1.first)->second.push_back(traccia);

                auto ret2 = NP_traces.insert({id_frc_beta0, {traccia}});
                if(!ret2.second)
                    (ret2.first)->second.push_back(traccia);
            }
            else if ( id_frc_beta1 == id_frc_beta2){
                auto ret1 = P_traces.insert({id_frc_beta1, {traccia}});
                if(!ret1.second)
                    (ret1.first)->second.push_back(traccia);

                auto ret2 = NP_traces.insert({id_frc_beta0, {traccia}});
                if(!ret2.second)
                    (ret2.first)->second.push_back(traccia);}

            else{
                auto ret1 = NP_traces.insert({frc1.id, {traccia}});
                if(!ret1.second)
                    (ret1.first)->second.push_back(traccia);

                auto ret2 = NP_traces.insert({frc2.id, {traccia}});
                if(!ret2.second)
                    (ret2.first)->second.push_back(traccia);}

        }



    } //chiude l'else

}


} //namespace
