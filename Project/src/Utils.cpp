#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

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
    fract.coordinates.resize(fract.num_fractures*7);

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
            for (unsigned int i = 0; i < num_vertici; i++)
            {
                fract.coordinates[contatore_vertici] = {vertices(0,i), vertices(1,i), vertices(2,i)};
                fract.vertices_fractures[id].push_back(contatore_vertici);

                contatore_vertici += 1;

            }
        }
        else {cont += 1;}

    }


    file.close();
    return true;
}


} //namespace
