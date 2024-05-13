#include "src/FracturesLibrary.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include "Eigen/Eigen"
#include <fstream>

using namespace FracturesLibrary;
using namespace Eigen;
using namespace std;

inline bool compare_tracce(const Trace& trc1, const Trace& trc2){
    return trc1.len > trc2.len;
}


int main()
{
    cout << setprecision(16);
    string path = "DFN";
    string filenameI = path + "/FR82_data.txt";
    Fractures frc;
    list<Trace> list_traces; //lista delle tracce

    if (!importData(filenameI, frc)){return 1;}
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
            unsigned int id_fract1 = i;
            unsigned int id_fract2 = j;
            if( NearFractures(frc, id_fract1,id_fract2)){
                IntersectionFractures(frc,i,j, list_traces);
            }

        }
    }


    //primo file di output
    string filenameO_tracce = "tracce.txt";
    ofstream ofs;
    ofs.open(filenameO_tracce);

    if(ofs.fail()){cerr << "file opened fail." << endl; return 1;}

    ofs << "# Number of Traces" << endl << list_traces.size() << endl;
    ofs << "# TraceId; FractureId1; FractureId2; X1; Y1; Z1; X2; Y2; Z2" << endl;
    const string del = " ; ";

    for (auto traccia : list_traces){
        ofs << traccia.id << del << traccia.id_frc1 << del << traccia.id_frc2 << del << traccia.coordinates_extremes(0,0) << del << traccia.coordinates_extremes(1,0) << del << traccia.coordinates_extremes(2,0) << del << traccia.coordinates_extremes(0,1) << del << traccia.coordinates_extremes(1,1) << del << traccia.coordinates_extremes(2,1) << endl;
    }

    ofs.close();


    //secondo file di output
    string filenameO_frc = "tracce_per_frattura.txt";
    ofstream ofs1;
    ofs1.open(filenameO_frc);

    if(ofs1.fail()){cerr << "file opened fail." << endl; return 1;}

    string header_frc = "# FractureId; NumTraces";
    string header_trc = "# TraceId; Tips; Lenght";

    for (unsigned int i = 0; i < frc.num_fractures; i++){
        ofs1 << header_frc << endl;
        ofs1 << i << del << frc.P_traces_of_fractures[i].size() + frc.NP_traces_of_fractures[i].size() << endl;

        if (frc.P_traces_of_fractures[i].size() + frc.NP_traces_of_fractures[i].size() != 0){
            ofs1 << header_trc << endl;

            //prima ordiniamo le liste
            (frc.P_traces_of_fractures[i]).sort(compare_tracce);
            (frc.NP_traces_of_fractures[i]).sort(compare_tracce);

            //scorriamo la lista delle tracce passanti
            for (auto traccia : frc.P_traces_of_fractures[i]){
                ofs1 << traccia.id << del << "false" << del << traccia.len << endl;
            }

            //scorriamo la lista delle tracce non passanti
            for (auto traccia : frc.NP_traces_of_fractures[i]){
                ofs1 << traccia.id << del << "true"  << del << traccia.len << endl;
            }
        }
        ofs1 << endl;
    }

    ofs1.close();


    return 0;
}
