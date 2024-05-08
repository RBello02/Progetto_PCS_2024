#include "Utils.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>

namespace GeometryLibrary{

bool importData(const string& path, Fractures& fract){

    ifstream file;
    file.open(path);

    if(file.fail()){
        return false;}

    //la prima riga è di header
    string header1;
    getline(file, header1);

    // mi salvo il numero di fratture
    file >> fract.num_fractures;

    //riservo lo spazio necessario ai vettori
    fract.id_fractures.reserve(fract.num_fractures);
    fract.dim_fractures.reserve(fract.num_fractures);
    fract.vertices_fractures.reserve(fract.num_fractures);

    fract.coordinates.reserve(fract.num_fractures);

    unsigned int cont = 0; //contatore che mi aiuta a salvarmi i dati
    unsigned int id;
    unsigned int num_vertici = 0;
    MatrixXd vertices;
    string line;

    while (!file.eof()){
        getline(file, line);
        cout << "contatore  " << cont << endl;
        cout << "line " << line << endl;
        istringstream convert(line);

        if (cont == 2){
                        char del; //mi serve per fermare la conversione ogni volta che incontro un ;
                        convert >> id >> del >> num_vertici;
                        fract.id_fractures.push_back(id);
                        fract.dim_fractures.push_back(num_vertici);
                        cout << "num vertici " << num_vertici << endl;
                        vertices = MatrixXd::Zero(3, num_vertici);}

        else if (cont < 6 && cont > 3){
                for (unsigned int i = 0; i < num_vertici -1; i++){
                char del; //mi serve per fermare la conversione ogni volta che incontro un ;
                convert >> vertices(cont-4,i) >> del;
            }
            convert >> vertices(cont-4,num_vertici-1);


        }

        if (cont == 5){cont = 0;}
        else {cont += 1;}

    }


    file.close();
    return true;
}


} //namespace
