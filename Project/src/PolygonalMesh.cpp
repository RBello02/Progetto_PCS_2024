#include "PolygonalMesh.hpp"
#include "FracturesLibrary.hpp"
#include "reshaping_array.hpp"

using namespace PolygonalLibrary;
using namespace FracturesLibrary;

namespace FracturesLibrary{

//funzioni di supporto
inline MatrixXd Retta_per_due_punti(Vector3d& pt1, Vector3d& pt2)
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

inline Vector2d intersezione_rette(MatrixXd& r_frattura, MatrixXd& r_traccia)
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

inline double appartiene_a_segmento(Vector3d& origin, Vector3d& end, Vector3d& pto, double toll){
    bool appartiene = false;

    //calcolo la distanza di pto dai due estremi del segmento, se la somma di queste due distanza è maggiore della distanza dei
    // due estremi allora il pto è esterno
    double lung_segmento = (origin-end).norm();
    double pto_o = (origin-pto).norm();
    double pto_e = (end-pto).norm();

    if(abs(pto_o + pto_e - lung_segmento) < toll){appartiene = true;}

    return appartiene;

}

/*****************************************************************************************************************************/

void divisione_sottopol(const Fracture& frattura, list<Trace> P_traces,  list<Trace> NP_traces, PolygonalMesh& mesh, const double& toll, list<unsigned int> lista_vert){
    /*questa procedura sarà definita ricorsivamente:
        * prima divido la frattura in due in base alla prima traccia che dobbiamo considerare
        * vado a dividere le tracce restanti in base a quale poligono appartengono
        * if(traccia i poligono  == 0){salvo il poligono nella mesh}
        * else{richiamo questa funzione}

    *se sono entrata in questa procedura è perché c'è almeno una traccia associata alla frattura
    */

    //ricordiamo: prima le tracce passanti, poi quelle non passanti in ordine decrescente (le liste sono già ordinate)
    Trace traccia_tagliante;
    if (P_traces.size() == 0){
        auto it = P_traces.begin();
        traccia_tagliante = *it;
        P_traces.pop_front();
    }
    else{
        auto it = P_traces.begin();
        traccia_tagliante = *it;
        NP_traces.pop_front();
    }

    Vector3d pt1 = traccia_tagliante.coordinates_extremes.col(0);
    Vector3d pt2 = traccia_tagliante.coordinates_extremes.col(1);
    MatrixXd retta_traccia = Retta_per_due_punti(pt1, pt2);
    Vector3d dir_t = retta_traccia.row(0);

    //ciclo sui lati della frattura per cercare i due punti di intersezione
    array<unsigned int, 2> nuovi_punti;
    unsigned int contatore = 0;

    for (unsigned int i = 0; i< frattura.num_vertici; i++){
        unsigned int id_o = frattura.vertices[i];
        unsigned int id_e;
        if(i==frattura.num_vertici-1){id_e = frattura.vertices[0];}
        else{id_e = frattura.vertices[i+1];}

        Vector3d origin = mesh.Cell0DCoordinates[id_o];
        Vector3d end = mesh.Cell0DCoordinates[id_e];
        MatrixXd retta_fratt = Retta_per_due_punti(origin, end);

        //interseco le rette se queste non sono parallele
        Vector3d dir_f = retta_fratt.row(0);
        bool non_parallele = !((dir_f.cross(dir_t)).norm() < toll);

        if(non_parallele){
            Vector2d a_b = intersezione_rette(retta_fratt, retta_traccia);
            //il punto deve stare sul segmento della frattura, ossia alpha appartiene a [0,1]
            if (a_b[0] >= -toll && a_b[0] <= 1+ toll){
                Vector3d pto = retta_traccia.row(1).transpose() + a_b[1] * dir_t;

                //salvo i due punti nella lista e le loro coordinate nella mesh
                unsigned int id = mesh.NumberCell0D;
                nuovi_punti[contatore] = id;
                auto pos = find(lista_vert.begin(), lista_vert.end(), id_e);
                lista_vert.insert(pos, id);

                //aggiorno le celle 0D della mesh
                mesh.NumberCell0D += 1;
                ReshapingArray::VerificaRaddoppio(mesh.Cell0DCoordinates);
                ReshapingArray::VerificaRaddoppio(mesh.Cell0DId);
                mesh.Cell0DCoordinates.push_back(pto);
                mesh.Cell0DId.push_back(id);

                contatore+=1;
            }
        }
    } //end for

    //creo le due sottofratture e mi ricalcolo le loro info
    Fracture frc1;
    frc1.vertices.reserve(10);
    Fracture frc2;
    frc2.vertices.reserve(10);

    //per suddividere i vertici nelle due fratt scorro la lista e avrò un booleano che switcho appena trovo uno dei due nuovi vertici
    bool flag = true;
    for(unsigned int elem : lista_vert){
        if (elem == nuovi_punti[0] || elem == nuovi_punti[1]){
            //assegno l'elemento ad entrambe le fratture
            frc1.vertices.push_back(elem);
            ReshapingArray::VerificaRaddoppio(frc1.vertices);
            frc2.vertices.push_back(elem);
            ReshapingArray::VerificaRaddoppio(frc2.vertices);

            //switcho la flag
            flag = !flag;
        }
        else{
            if (flag){
                frc1.vertices.push_back(elem);
                ReshapingArray::VerificaRaddoppio(frc1.vertices);
            }
            else{
                frc2.vertices.push_back(elem);
                ReshapingArray::VerificaRaddoppio(frc2.vertices);
            }

        }

    }

    frc1.vertices.shrink_to_fit();
    frc2.vertices.shrink_to_fit();
    frc1.num_vertici = frc1.vertices.size();
    frc2.num_vertici = frc1.vertices.size();

    //creo una lista di vertici per ogni sottofrattura
    list<unsigned int> vert_fract1;
    for (unsigned int i = 0; i< frc1.num_vertici; i++){
        vert_fract1.push_back(frc1.vertices[i]);
    }

    list<unsigned int> vert_fract2;
    for (unsigned int i = 0; i< frc2.num_vertici; i++){
        vert_fract1.push_back(frc2.vertices[i]);
    }

    //determino le tracce passanti e quelle non per ogni sotto-frattura
    list<Trace> P_traces1;
    list<Trace> NP_traces1;
    list<Trace> P_traces2;
    list<Trace> NP_traces2;


    //MI MANCA DA VERIFICRE QUALI TRACCE APPARTENGANO A QUALI LISTE


    list<unsigned int> list_vert1;
    for (unsigned int i = 0; i< frc1.num_vertici; i++){
        list_vert1.push_back(frc1.vertices[i]);
    }

    list<unsigned int> list_vert2;
    for (unsigned int i = 0; i< frc2.num_vertici; i++){
        list_vert2.push_back(frc2.vertices[i]);
    }

    //se non ho tracce interne interrompo la ricorsione e salvo la mesh
    if (P_traces1.size() + NP_traces1.size() == 0){
        //le celle 0D le ho salvate man mano

        //salvo le celle 1D
        mesh.NumberCell1D += frc1.num_vertici;
        for (unsigned int v = 0; v < frc1.num_vertici; v++){
            unsigned int n = mesh.Cell1DId.size();
            unsigned int id_origin = frc1.vertices[v];
            unsigned int id_end;
            if(v == frc1.num_vertici){id_end = frc1.vertices[0];}
            else{id_end = frc1.vertices[v+1];}

            ReshapingArray::VerificaRaddoppio(mesh.Cell1DId);
            mesh.Cell1DId.push_back(n+v);
            mesh.Cell1DVertices.push_back({id_origin, id_end});
        }

        //salvo le celle 2D
        mesh.NumberCell2D +=1;
        ReshapingArray::VerificaRaddoppio(mesh.Cell2DId);
        ReshapingArray::VerificaRaddoppio(mesh.Cell2DEdges);
        ReshapingArray::VerificaRaddoppio(mesh.Cell2DVertices);
        unsigned int id_pol = mesh.NumberCell2D -1;
        mesh.Cell2DId.push_back(id_pol);


    }
    else{divisione_sottopol(frc1, P_traces1, NP_traces1, mesh, toll, list_vert1);}

}

PolygonalMesh FracturesFunctions::SottoPoligonazione(const Fracture& frattura, const list<Trace>& P_traces, const list<Trace>& NP_traces, const vector<Vector3d>& coord){
    //ridefinisco la frattura in modo da iniziare l'indicizzazione dei vertici da 0
    Fracture frc;
    vector<Vector3d> coord_frc; //conterrà solo i vertici di questa frattura
    frc.num_vertici = frattura.num_vertici;
    frc.vertices.reserve(frc.num_vertici);
    coord_frc.reserve(frattura.num_vertici);

    PolygonalMesh mesh;
    mesh.NumberCell0D += frc.num_vertici;
    mesh.Cell0DId.reserve(50);
    mesh.Cell0DCoordinates.reserve(50);

    for (unsigned int i = 0; i< frattura.num_vertici; i++){
        Vector3d pto = coord[frattura.vertices[i]];
        coord_frc.push_back(pto);
        frc.vertices.push_back(i);
        mesh.Cell0DCoordinates.push_back(pto);
        mesh.Cell0DId.push_back(i);
    }

    //mi creo una lista di vertici per poterli ordinare in maniera più easy
    list<unsigned int> vert_fract;
    for (unsigned int i = 0; i< frc.num_vertici; i++){
        vert_fract.push_back(frc.vertices[i]);
    }

    divisione_sottopol(frc, P_traces, NP_traces, mesh, tolleranza1D, vert_fract);

    return mesh;

}


} //namespace
