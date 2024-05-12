#include "Utils.hpp"
#include "reshaping_array.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>


#include <cmath>

#include "Eigen/Eigen"
using namespace Eigen;

namespace GeometryLibrary{

bool importData(const string& path, Fractures& fract){

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
    file >> fract.num_fractures;

    //riservo lo spazio necessario ai vettori
    fract.id_fractures.reserve(fract.num_fractures);
    fract.dim_fractures.reserve(fract.num_fractures);
    fract.vertices_fractures.resize(fract.num_fractures);

    //non so a priori quanti vertici ho considerando tutte le fratture, considero una stima pessimistica basata sul fatto che ogni frattura ha 7 vertici
    fract.coordinates.reserve(fract.num_fractures*7);

    //mi creo delle variabili ausiliari che mi aiutano a salvarmi i dati
    unsigned int cont = 0; //contatore che mi aiuta a capire per ogni frattura che dati sto leggendo
    unsigned int id;
    unsigned int num_vertici = 0;
    MatrixXd vertices;
    string line;
    getline(file, line);

    unsigned int contatore_vertici = 0; //questo contatore mi serve ad identificare in maniera univoca tra tutte le fratture i punti

    while (!file.eof()){
        getline(file, line);
        istringstream convert(line);

        //faccio una casistica basata su che riga sto leggendo, le righe che coincidono con cont == 0 o cont == 2 sono di intestazione quindi le saltiamo
        if (cont == 1){
                        char del; //mi serve per fermare la conversione ogni volta che incontro un ;
                        convert >> id >> del >> num_vertici;
                        fract.id_fractures.push_back(id);
                        fract.dim_fractures.push_back(num_vertici);
                        fract.vertices_fractures[id].reserve(num_vertici);
                        vertices = MatrixXd::Zero(3, num_vertici);}

        else if (cont < 6 && cont > 2){
                for (unsigned int i = 0; i < num_vertici -1; i++){
                    char del; //mi serve per fermare la conversione ogni volta che incontro un ;
                    convert >> vertices(cont-3,i) >> del;
                }
            convert >> vertices(cont-3,num_vertici-1);


        }

        //incremento il contatore o lo reimposto a 0 quando finisco le info relative ad una frattura
        if (cont == 5){
            cont = 0;

            //arrivata a questo punto mi sono salvata tutte le info riguardanti una frattura quindi vado ad aggiungerle alla strutture dati
            for (unsigned int i = 0; i < num_vertici; i++){
                fract.coordinates.push_back({vertices(0,i), vertices(1,i), vertices(2,i)});
                fract.coordinates = ReshapingArray::VerificaRaddoppio(fract.coordinates);

                fract.vertices_fractures[id].push_back(contatore_vertici);

                contatore_vertici += 1;

            }
        }
        else {cont += 1;}

    }

    fract.coordinates.shrink_to_fit(); //elimino la capacità in eccesso
    file.close();
    return true;
}


bool NearFractures(const Fractures& frc, unsigned int id_fract1, unsigned int id_fract2){

    // vettori per le coordinate dei due baricentri (approssimativamente):
    array<double, 3> bar1;
    array<double,3> bar2;
    //calcolo le coordinate facendo somma/numvertici per ogni coordinata
    double sommax=0,sommay=0,sommaz=0;

    for (unsigned int k = 0; k < frc.dim_fractures[id_fract1]; k++){    // ciclo sui vertici della frattura e sommo tutte le coordinate

        unsigned int id_vertice = frc.vertices_fractures[id_fract1][k];
        sommax += frc.coordinates[id_vertice][0];
        sommay += frc.coordinates[id_vertice][1];
        sommaz += frc.coordinates[id_vertice][2];

    };
    bar1[0]=sommax/frc.dim_fractures[id_fract1];                // calcolo la coordinata del baricentro dividendo la somma delle coordinate per il numero di vertici
    bar1[1]=sommay/frc.dim_fractures[id_fract1];
    bar1[2]=sommaz/frc.dim_fractures[id_fract1];


    for (unsigned int k = 0; k < frc.dim_fractures[id_fract2]; k++){    // ripeto il calcolo anche per la seconda frattura

        unsigned int id_vertice = frc.vertices_fractures[id_fract2][k];
        sommax += frc.coordinates[id_vertice][0];
        sommay += frc.coordinates[id_vertice][1];
        sommaz += frc.coordinates[id_vertice][2] ;
    };
    bar2[0]=sommax/frc.dim_fractures[id_fract2];
    bar2[1]=sommay/frc.dim_fractures[id_fract2];
    bar2[2]=sommaz/frc.dim_fractures[id_fract2];

    //calcolo i raggi delle sfere(al quadrato) e la distanza tra i due baricentri
    double raggio1=0,raggio2=0;

    // cerco il vertice con distanza massima dal baricentro, in realtà mi interessa solo la distanza

    double raggio_da_confrontare_1 = 0;

    for (unsigned int j = 0; j<frc.dim_fractures[id_fract1]; j++)   // ciclo su tutti i poligoni per calcolare la distanza massima dei vertici dal baricentro
    {
        unsigned int id_vertice = frc.vertices_fractures[id_fract1][j];  // inizializzo l'id del vertice
        for (unsigned int i =0; i<3;i++)        // ciclo sulle 3 coordinate per calcolare il raggio del poligono
        {
            raggio1 += pow((bar1[i]-frc.coordinates[id_vertice][i]),2);
        }
        if (raggio1 >= raggio_da_confrontare_1)
        {
            raggio_da_confrontare_1 = raggio1;
        }
    }

    double raggio_da_confrontare_2 = 0;

    for (unsigned int j = 0; j<frc.dim_fractures[id_fract2]; j++)   // ciclo su tutti i poligoni per calcolare la distanza massima dei vertici dal baricentro
    {
        unsigned int id_vertice = frc.vertices_fractures[id_fract2][j];  // inizializzo l'id del vertice
        for (unsigned int i =0; i<3;i++)        // ciclo sulle 3 coordinate per calcolare il raggio del poligono
        {
            raggio2 += pow((bar1[i]-frc.coordinates[id_vertice][i]),2);
        }
        if (raggio2 >= raggio_da_confrontare_2)
        {
            raggio_da_confrontare_2 = raggio2;
        }
    }

    double distbb = 0;

    for (unsigned int j = 0; j<3; j++)
    {
        distbb += pow(bar2[j]-bar1[j],2);
    }


    if (raggio_da_confrontare_1+raggio_da_confrontare_2<distbb)
        return false;

    else
        return true;

}


} //namespace
