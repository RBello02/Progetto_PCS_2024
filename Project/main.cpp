#include "src/FracturesLibrary.hpp"
#include "src/PolygonalMesh.hpp"
#include "src/Utils.hpp"
#include "UCDUtilities.hpp"

#include <iomanip>
#include <iostream>
#include <string>
#include "Eigen/Eigen"
#include <fstream>

using namespace FracturesLibrary;
using namespace PolygonalLibrary;
using namespace UtilsFunction;
using namespace Eigen;
using namespace std;

inline bool compare_tracce(const Trace& trc1, const Trace& trc2){
    return trc1.len > trc2.len;
}


int main()
{
    cout << setprecision(16);

    string path = "DFN";
    string filenameI = path + "/FR3_data.txt";

    //definisco le liste che conterranno le tracce e le fratture
    vector<Fracture> list_fractures; //lista di fratture (è un vettore)
    list<Trace> list_traces; //lista delle tracce
    map<unsigned int, list<Trace>> P_traces_of_fractures; //per ogni frattura memorizziamo una lista contenente gli id elle tracce passanti
    map<unsigned int, list<Trace>> NP_traces_of_fractures; //analogo a sopra ma per tracce non passanti

    //mi creo un vettore di array dove mi salvo tutte le coordinate
    vector<Vector3d> coordinates;
    unsigned int num_fratt;

    FracturesFunctions f;


    if (!f.importData(filenameI, list_fractures, coordinates)){return 1;}   // importo i dati
    else{
        //stampo un poo' di roba per verificare che sia tutto giusto
        num_fratt = list_fractures.size();
        cout << "Ho un numero di fratture pari a " << num_fratt << endl;

        for (Fracture frc : list_fractures)
        {

            cout << "FRATTURA NUM " << frc.id <<  " che ha " << frc.num_vertici << " vertici." << endl;

            for (unsigned int k = 0; k < frc.num_vertici; k++){

                unsigned int id_vertice = frc.vertices[k];
                cout << "Il vertice " << k << ", che ha id = " << id_vertice << " , ha coord (" << coordinates[id_vertice][0] << " , "  << coordinates[id_vertice][1] <<  " , " << coordinates[id_vertice][2] << " )" << endl;
             }
            cout << endl;
        }


    }


    // ciclo sulle coppie di poligoni e determino le tracce
    for(unsigned int i=0; i<num_fratt; i++)
    {
        for (unsigned int j=i+1; j< num_fratt; j++)
        {
            Fracture frc1 = list_fractures[i];
            Fracture frc2 = list_fractures[j];

            if( f.NearFractures(frc1, frc2, coordinates)){
                f.IntersectionFractures(frc1, frc2, coordinates, list_traces, P_traces_of_fractures, NP_traces_of_fractures);
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

    for (Trace traccia : list_traces){
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

    for (Fracture& fratt : list_fractures){

        ofs1 << header_frc << endl;
        ofs1 << fratt.id << del << P_traces_of_fractures[fratt.id].size() + NP_traces_of_fractures[fratt.id].size() << endl;

        if (P_traces_of_fractures[fratt.id].size() + NP_traces_of_fractures[fratt.id].size() != 0){
            ofs1 << header_trc << endl;

            //prima ordiniamo le liste
            (P_traces_of_fractures[fratt.id]).sort(compare_tracce);
            (NP_traces_of_fractures[fratt.id]).sort(compare_tracce);

            //scorriamo la lista delle tracce passanti
            for (Trace traccia : P_traces_of_fractures[fratt.id]){
                ofs1 << traccia.id << del << "false" << del << traccia.len << endl;
            }

            //scorriamo la lista delle tracce non passanti
            for (Trace traccia : NP_traces_of_fractures[fratt.id]){
                ofs1 << traccia.id << del << "true"  << del << traccia.len << endl;
            }
        }
        ofs1 << endl;
    }

    ofs1.close();

    //PARTE 2 PROGETTO
    vector<PolygonalMesh> sottoPoligonazione_per_frattura;  // creo un vettore di oggetti poligonalmesh
    sottoPoligonazione_per_frattura.reserve(num_fratt);     // riservo la memoria corrispondente al numero di fratture
    //calcolo ora la sottopoligonazione per ogni frattura
    for(const Fracture &frattura : list_fractures)   // ciclo sulle fratture
    {
        PolygonalMesh mesh;     // inizializzo la mesh
        mesh = f.FracturesFunctions::SottoPoligonazione(frattura, P_traces_of_fractures[frattura.id], NP_traces_of_fractures[frattura.id], coordinates);

        //aggiungo la mesh creata al vettore
        sottoPoligonazione_per_frattura.push_back(mesh);
    }


    // ESPORTAZIONE IN PARAVIEWx
    // string name="Mesh_";
    //     for(unsigned int i=0; i< sottoPoligonazione_per_frattura.size();i++){
    //        PolygonalMesh mesh=sottoPoligonazione_per_frattura[i];
    //         name=name+to_string(i);
    //         MatrixXd punti;
    //         punti.resize(3,mesh.NumberCell0D);
    //         for(unsigned int j=0; j<mesh.NumberCell0D;j++){
    //             punti.col(j)=mesh.Cell0DCoordinates[j];
    //         }
    //         vector<Gedim::UCDCell> cells;
    //         VectorXi materials;
    //         Gedim::UCDUtilities::CreatePointCells(punti,materials);
    //     }
        

    return 0;
}
