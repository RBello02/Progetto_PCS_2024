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

inline double pto_unico(Vector3d& pto, vector<Vector3d>& punti, double toll, unsigned int& id){
    //l'id mi serve quando cerco un punto tra i vertici

    bool unico = true;

    if(punti.size() == 0){return unico;}

    for (unsigned int i = 0; i < punti.size(); i++){
        Vector3d elem = punti[i];
        bool uguagl_x = (abs(pto[0] - elem [0]) < toll);
        bool uguagl_y = (abs(pto[1] - elem [1]) < toll);
        bool uguagl_z = (abs(pto[2] - elem [2]) < toll);

        if (uguagl_x && uguagl_y && uguagl_z){unico = false; id = i; return unico;}
    }
    return unico;
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
    if (P_traces.size() != 0){
        auto it = P_traces.begin();
        traccia_tagliante = *it;
        P_traces.pop_front();
    }
    else{
        auto it = NP_traces.begin();
        traccia_tagliante = *it;
        NP_traces.pop_front();
    }

    Vector3d pt1 = traccia_tagliante.coordinates_extremes.col(0);
    Vector3d pt2 = traccia_tagliante.coordinates_extremes.col(1);
    MatrixXd retta_traccia = Retta_per_due_punti(pt1, pt2);
    Vector3d dir_t = retta_traccia.row(0);

    //ciclo sui lati della frattura per cercare i due punti di intersezione (devo gestire il caso di capotare su un vertice --> conta doppio)
    vector<Vector3d> nuovi_punti;
    vector<unsigned int> id_nuoviPunti;
    nuovi_punti.reserve(2); //male che vada ho 4 beta
    id_nuoviPunti.reserve(4);


    for (unsigned int i = 0; i< frattura.num_vertici; i++){
        unsigned int id_o = frattura.vertices[i];
        unsigned int id_e;
        if(i==frattura.num_vertici-1){id_e = frattura.vertices[0];}
        else{id_e = frattura.vertices[i+1];}

        Vector3d origin = mesh.Cell0DCoordinates[id_o];
        Vector3d end = mesh.Cell0DCoordinates[id_e];
        MatrixXd retta_fratt = Retta_per_due_punti(origin, end);
        Vector3d dir_f = retta_fratt.row(0);

        //interseco le rette se queste non sono parallele
        bool non_parallele = !((dir_f.cross(dir_t)).norm() < toll);

        if(non_parallele){
            Vector2d a_b = intersezione_rette(retta_fratt, retta_traccia);
            //il punto deve stare sul segmento della frattura, ossia alpha appartiene a [0,1]
            if (a_b[0] > -toll && a_b[0] < 1+ toll){
                Vector3d pto = retta_traccia.row(1).transpose() + a_b[1] * dir_t;
                //se pto è già in nuovi_punti lo ignoro sennò faccio quello che devo
                unsigned int a = 0;
                if(pto_unico(pto, nuovi_punti, toll, a)){

                    //devo gesture anche il caso in cui il punto coincide con un vertice della frattura
                    vector<Vector3d> coord_frc;
                    coord_frc.reserve(frattura.num_vertici);
                    for(unsigned int i = 0; i<frattura.num_vertici; i++){
                        coord_frc.push_back(mesh.Cell0DCoordinates[frattura.vertices[i]]);
                    }
                    unsigned int id_new_pto = 0;
                    if (!pto_unico(pto, coord_frc, toll, id_new_pto)){
                        //in questo caso non devo aggiungere nulla alla mesh
                        nuovi_punti.push_back(pto);
                        id_nuoviPunti.push_back(id_new_pto);
                    }
                    else{
                    nuovi_punti.push_back(pto);
                    id_new_pto = mesh.NumberCell0D;
                    id_nuoviPunti.push_back(id_new_pto);
                    ReshapingArray::VerificaRaddoppio(mesh.Cell0DId);
                    mesh.Cell0DId.push_back(id_new_pto);
                    ReshapingArray::VerificaRaddoppio(mesh.Cell0DCoordinates);
                    mesh.Cell0DCoordinates.push_back(pto);
                    mesh.NumberCell0D += 1;

                    //inserisco l'id del punto nella lista prima del vertice di end del segmento
                    auto pos = find(lista_vert.begin(), lista_vert.end(), id_e);
                    lista_vert.insert(pos, id_new_pto);
                    }
                }

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
        if (elem == id_nuoviPunti[0] || elem == id_nuoviPunti[1]){
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
    list<Trace> P_traces1 = {};
    list<Trace> NP_traces1 = {};
    list<Trace> P_traces2 = {};
    list<Trace> NP_traces2  = {};


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

        //salvo le celle 1D e 2D
        unsigned int num_cell1d = mesh.NumberCell1D;
        mesh.NumberCell1D += frc1.num_vertici;
        unsigned int num_celle2d = mesh.NumberCell2D;
        mesh.NumberCell2D +=1;

        //mi creo due vettori ausiliari
        vector<unsigned int> vertices;
        vertices.reserve(frc1.num_vertici);
        vector<unsigned int> edges;
        edges.reserve(frc1.num_vertici);

        for (unsigned int v = 0; v < frc1.num_vertici; v++){
            unsigned int id_origin = frc1.vertices[v];
            unsigned int id_end;
            if(v == frc1.num_vertici){id_end = frc1.vertices[0];}
            else{id_end = frc1.vertices[v+1];}

            //celle 1D
            unsigned int id_segm = num_cell1d+v;
            ReshapingArray::VerificaRaddoppio(mesh.Cell1DId);
            mesh.Cell1DId.push_back(id_segm);
            ReshapingArray::VerificaRaddoppio(mesh.Cell1DVertices);
            mesh.Cell1DVertices.push_back({id_origin, id_end});

            //celle2D
            vertices.push_back(id_origin);
            edges.push_back(id_segm);
        }

        ReshapingArray::VerificaRaddoppio(mesh.Cell2DId);
        mesh.Cell2DId.push_back(num_celle2d);
        ReshapingArray::VerificaRaddoppio(mesh.Cell2DEdges);
        mesh.Cell2DEdges.push_back(edges);
        ReshapingArray::VerificaRaddoppio(mesh.Cell2DVertices);
        mesh.Cell2DVertices.push_back(vertices);



    }
    else{divisione_sottopol(frc1, P_traces1, NP_traces1, mesh, toll, list_vert1);}

    if (P_traces2.size() + NP_traces2.size() == 0){
        //le celle 0D le ho salvate man mano

        //salvo le celle 1D e 2D
        unsigned int num_cell1d = mesh.NumberCell1D;
        mesh.NumberCell1D += frc2.num_vertici;
        unsigned int num_celle2d = mesh.NumberCell2D;
        mesh.NumberCell2D +=1;

        //mi creo due vettori ausiliari
        vector<unsigned int> vertices;
        vertices.reserve(frc2.num_vertici);
        vector<unsigned int> edges;
        edges.reserve(frc2.num_vertici);

        for (unsigned int v = 0; v < frc2.num_vertici; v++){
            unsigned int id_origin = frc2.vertices[v];
            unsigned int id_end;
            if(v == frc2.num_vertici){id_end = frc2.vertices[0];}
            else{id_end = frc2.vertices[v+1];}

            //celle 1D
            unsigned int id_segm = num_cell1d+v;
            ReshapingArray::VerificaRaddoppio(mesh.Cell1DId);
            mesh.Cell1DId.push_back(id_segm);
            ReshapingArray::VerificaRaddoppio(mesh.Cell1DVertices);
            mesh.Cell1DVertices.push_back({id_origin, id_end});

            //celle2D
            vertices.push_back(id_origin);
            edges.push_back(id_segm);
        }

        ReshapingArray::VerificaRaddoppio(mesh.Cell2DId);
        mesh.Cell2DId.push_back(num_celle2d);
        ReshapingArray::VerificaRaddoppio(mesh.Cell2DEdges);
        mesh.Cell2DEdges.push_back(edges);
        ReshapingArray::VerificaRaddoppio(mesh.Cell2DVertices);
        mesh.Cell2DVertices.push_back(vertices);



    }
    else{divisione_sottopol(frc2, P_traces2, NP_traces2, mesh, toll, list_vert2);}

}

PolygonalMesh FracturesFunctions::SottoPoligonazione(const Fracture& frattura, const list<Trace>& P_traces, const list<Trace>& NP_traces, const vector<Vector3d>& coord){
    //ridefinisco la frattura in modo da iniziare l'indicizzazione delle celle 0D da 0
    Fracture frc;
    vector<Vector3d> coord_frc; //conterrà solo i vertici di questa frattura
    frc.num_vertici = frattura.num_vertici;
    frc.vertices.reserve(frc.num_vertici);
    coord_frc.reserve(frattura.num_vertici);

    PolygonalMesh mesh;
    mesh.NumberCell0D += frc.num_vertici;
    mesh.Cell0DId.reserve(50);
    mesh.Cell0DCoordinates.reserve(50);

    //inizio già a salvare le coordinate 0D
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
