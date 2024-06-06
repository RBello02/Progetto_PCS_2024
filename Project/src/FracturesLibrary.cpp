#include "FracturesLibrary.hpp"
#include "reshaping_array.hpp"
#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <array>

#include <cmath>

#include "Eigen/Eigen"

using namespace Eigen;
using namespace UtilsFunction;

namespace UtilsFunction{

inline bool compare_beta(const array<double,2>& arr1, const array<double,2>& arr2){
    return arr1[1] > arr2[1];
}

/****************************************************************************************************************/

bool FracturesFunctions::importData(const string& path, vector<Fracture>& lista, vector<Vector3d>& coord){

    //apro il file
    ifstream file;
    file.open(path);

    if(file.fail()){
        cerr << "errore nell'apertura del file di input" << endl;
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
            double z2;
            convert2 >> z2;
            vert(k,frc.num_vertici-1) = z2;
        }

        //salvo ora le coordinate nel vettore coord e i rispettivi id nella struct Fracture
        for (unsigned int k = 0; k < frc.num_vertici; k++)
        {
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

    bool flag = true;
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


    sommax=0;
    sommay=0;
    sommaz=0;
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


    // cerco il vertice con distanza massima dal baricentro, in realtà mi interessa solo la distanza

    double raggio_da_confrontare_1 = 0;

    for (unsigned int j = 0; j<frc1.num_vertici; j++)   // ciclo su tutti i poligoni per calcolare la distanza massima dei vertici dal baricentro
    {
        double raggio1=0;
        unsigned int id_vertice = frc1.vertices[j];  // inizializzo l'id del vertice
        for (unsigned int i =0; i<3;i++)        // ciclo sulle 3 coordinate per calcolare il raggio del poligono
        {
            raggio1 += (bar1[i]-coord[id_vertice][i])*(bar1[i]-coord[id_vertice][i]);
        }
        if (raggio1 >= raggio_da_confrontare_1 +tolleranza1D)
        {
            raggio_da_confrontare_1 = raggio1;
        }
    }

    double raggio_da_confrontare_2 = 0;

    for (unsigned int j = 0; j<frc2.num_vertici; j++)   // ciclo su tutti i poligoni per calcolare la distanza massima dei vertici dal baricentro
    {
        double raggio2=0;
        unsigned int id_vertice = frc2.vertices[j];  // inizializzo l'id del vertice
        for (unsigned int i =0; i<3;i++)        // ciclo sulle 3 coordinate per calcolare il raggio del poligono
        {
            raggio2 += (bar2[i]-coord[id_vertice][i])*(bar2[i]-coord[id_vertice][i]);
        }
        if (raggio2 >= raggio_da_confrontare_2 + tolleranza1D)
        {
            raggio_da_confrontare_2 = raggio2;
        }
    }

    double distbb = 0;

    for (unsigned int j = 0; j<3; j++)
    {
        distbb += (bar2[j]-bar1[j])*(bar2[j]-bar1[j]);
    }


    if (sqrt(raggio_da_confrontare_1)+sqrt(raggio_da_confrontare_2) < sqrt(distbb) + tolleranza1D){flag = false;}

    return flag;

}

void FracturesFunctions::IntersectionFractures(Fracture &frc1, Fracture &frc2, const vector<Vector3d>& coord, list<Trace>& list_traces, map<unsigned int, list<Trace> > &P_traces, map<unsigned int, list<Trace> > &NP_traces){
    Matrix3d piano_frc1 = frc1.calcolo_piano(coord);
    Matrix3d piano_frc2 = frc2.calcolo_piano(coord);

    FracturesFunctions fx;

    //procedo nella funzione solo nel caso in cui i due piani non siano paralleli
    if(Parallelismo(piano_frc1, piano_frc2)){
        cout << "La frattura " << frc1.id << " e la fratt " << frc2.id << " non si intersecano" << endl;}
    else{

        //calcolo la retta di intersezione tra i piani in forma parametrica
        MatrixXd retta_intersez_piani = fx.Retta_tra_piani(piano_frc1,piano_frc2);
        Vector3d dir_retta_intersez_piani;
        dir_retta_intersez_piani = retta_intersez_piani.row(0);

        //Devo ora ciclare sulle coppie di vertici delle due fratture per trovare gli alpha di intersezione tra la retta di
        //intersezione tra i due piani determinati dalle fratture e le rette determinate dai vertici
        list<array<double,2>> beta_inters; //--> nel primo array avrò gli alpha della prima frc e nel secondo della seconda

        for (unsigned int i = 0; i < frc1.num_vertici; i++){

            unsigned int j;
            if (i == frc1.num_vertici-1){j = 0;}
            else{j = i+1;}

            //mi trovo la retta per ogni coppia di vertici
            MatrixXd retta_tra_vertici = fx.Retta_per_due_vertici_della_frattura(frc1.vertices[i],frc1.vertices[j], coord);
            Vector3d dir_retta_tra_vertici;
            dir_retta_tra_vertici = retta_tra_vertici.row(0);

            //cerco l'eventuale intersezione tra questa retta e la retta di intrsezione tra i piani, queste esiste se le due rette non
            //sono parallele e complanari
            bool non_parallele = !((dir_retta_tra_vertici.cross(dir_retta_intersez_piani)).norm() < tolleranza1D);

            if (non_parallele){
                Vector2d a_b = fx.alpha_di_intersezione(retta_tra_vertici, retta_intersez_piani);

                if (a_b[0] >= -tolleranza1D && a_b[0] <= 1+tolleranza1D){
                    beta_inters.push_back({static_cast<double>(frc1.id),a_b[1]});
                }

            }
        }

        for (unsigned int i = 0; i < frc2.num_vertici; i++){

            unsigned int j;
            if (i == frc2.num_vertici-1){j = 0;}
            else{j = i+1;}

            //mi trovo la retta per ogni coppia di vertici
            MatrixXd retta_tra_vertici = fx.Retta_per_due_vertici_della_frattura(frc2.vertices[i],frc2.vertices[j], coord);
            Vector3d dir_retta_tra_vertici;
            dir_retta_tra_vertici = retta_tra_vertici.row(0);

            //cerco l'eventuale intersezione tra questa retta e la retta di intrsezione tra i piani, queste esiste se le due rette non
            //sono parallele e complanari
            bool non_parallele = !((dir_retta_tra_vertici.cross(dir_retta_intersez_piani)).norm() < tolleranza1D);

            if (non_parallele){
                Vector2d a_b = fx.alpha_di_intersezione(retta_tra_vertici, retta_intersez_piani);
                if (a_b[0] >= - tolleranza1D && a_b[0] <= 1+ tolleranza1D){
                    beta_inters.push_back({static_cast<double>(frc2.id),a_b[1]});
                }
                }
            }


        //se trovo una traccia la salvo
        beta_inters.sort(compare_beta);
        beta_inters.unique();  // tolgo i duplicati al fine di eliminare i casi degeneri (vedi nel test)


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

        bool non_traccia = (id_frc_beta0 == id_frc_beta1) && (id_frc_beta2 == id_frc_beta3);

        if (beta_inters.size()==4 && !non_traccia){
            cout << "Ho una traccia tra la frattura " << frc1.id << " e la fratt " << frc2.id << endl;
            //posso ordinare gli array tra loro in base al valore di beta
            beta_inters.sort(compare_beta);

            cout << endl;

            Trace traccia;
            traccia.id = P_traces.size() + NP_traces.size();
            MatrixXd coord_estremi_traccia = MatrixXd::Zero(3,2);


            coord_estremi_traccia.col(0) = origin;
            coord_estremi_traccia.col(1) = end;
            traccia.coordinates_extremes = coord_estremi_traccia;
            traccia.id_frc1 = frc1.id;
            traccia.id_frc2 = frc2.id;
            traccia.len = traccia.calcolo_lunghezza();

            list_traces.push_back(traccia);


            //determiniamo se la traccia è passante o no
            if(abs(beta0 -beta1) < tolleranza1D && abs(beta2-beta3) < tolleranza1D){
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
