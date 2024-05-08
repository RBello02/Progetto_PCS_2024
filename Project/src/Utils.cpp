#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>

namespace GeometryLibrary{

bool importData(const string& path, Fractures& fract){

    ifstream file;
    file.open(path);

    if(file.fail()){
        return false;}

    //la prima riga è di header
    string header1;
    file >> header1;

    // mi salvo il numero di fratture
    file >> fract.num_fractures;

    //riservo lo spazio necessario ai vettori
    fract.id_fractures.reserve(fract.num_fractures);
    fract.dim_fractures.reserve(fract.num_fractures);
    fract.vertices_fractures.reserve(fract.num_fractures);

    fract.coordinates.reserve(fract.num_fractures);

    unsigned int cont = 0; //contatore che mi aiuta a salvarmi i dati
    unsigned int id;
    unsigned int num_vertici;
    char del; //mi serve per fermare la conversione ogni volta che incontro un ;
    MatrixXd vertices;

    while (file.is_open()){
        cout << "Sono enyrato nel file" << endl;
        if (cont == 0 || cont == 2){string header;file >> header;}
        else if (cont == 1){file >> id >> del >> num_vertici;
                            fract.id_fractures.push_back(id);
                            fract.dim_fractures.push_back(num_vertici);
                            vertices = MatrixXd::Zero(3, num_vertici);}
        else if (cont == 3){
            //leggo le x
            for (unsigned int i = 0; i < num_vertici -1; i++){
                cout << "RENA HO RAGIONE IO" << endl;
                cout << "size of mat" << vertices.size() << endl;
                file >> vertices(0,i) >> del;
                cout << vertices(0,0);
            }
            file >> vertices(0,num_vertici-1);

            //leggo le y
            for (unsigned int i = 0; i < num_vertici -1; i++){
                file >> vertices(1,i) >> del;
            }
            file >> vertices(1,num_vertici-1);

            //leggo le z
            for (unsigned int i = 0; i < num_vertici -1; i++){
                file >> vertices(2,i) >> del;
            }
            file >> vertices(2,num_vertici-1);

        }

        if (cont == 3){cont = 0;}
        else {cont += 1;}

    }


    file.close();
    return true;
}


} //namespace
