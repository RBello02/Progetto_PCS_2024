#include "PolygonalMesh.hpp"
#include "FracturesLibrary.hpp"
#include "reshaping_array.hpp"
#include "Utils.hpp"
#include "Eigen/Eigen"
#include <iostream>

using namespace PolygonalLibrary;
using namespace FracturesLibrary;
using namespace UtilsFunction;
using namespace Eigen;

namespace UtilsFunction
{

inline vector<unsigned int> verifica_celle_1D(PolygonalMesh& mesh, Fracture& frattura){

    vector<unsigned int> edges;
    edges.reserve(frattura.num_vertici+1);

    vector<unsigned int> vertices = frattura.vertices;

    for(auto itor = frattura.vertices.begin(); itor != frattura.vertices.end(); itor++){
        unsigned int id_origin = (*itor);
        unsigned int id_end;
        if(itor == frattura.vertices.end()-1){id_end = frattura.vertices[0];}
        else{id_end = *(itor+1);}

        //CASO O: evito di salvare la cella 1D {a,b} se esiste già la cella {b,a}
        bool caso0 = false;
        unsigned int lato_dupl;

        for (unsigned int i = 0; i < mesh.Cell1DVertices.size(); i++){
            Vector2i cella1 = mesh.Cell1DVertices[i];


            //CASO 0
            if ( (unsigned(cella1[0]) == id_origin && unsigned(cella1[1]) == id_end) ||   (unsigned(cella1[0]) == id_end && unsigned(cella1[1]) == id_origin)){
                caso0 = true;
                lato_dupl = i;
                break;
            }
        }

        //salvo opportunamente la cella 1D
        if(caso0){edges.push_back(lato_dupl);}
        else{
            //devo creare una nuova cella 1D
            unsigned int new_id = mesh.Cell1DId.size();
            ReshapingArray::VerificaRaddoppio(mesh.Cell1DId);
            mesh.Cell1DId.push_back(new_id);
            ReshapingArray::VerificaRaddoppio(mesh.Cell1DVertices);
            mesh.Cell1DVertices.push_back({id_origin, id_end});
            edges.push_back(new_id);
        }
    }

    frattura.vertices = vertices;
    return edges;

}




/*****************************************************************************************************************************/

void sottopoligonazione_ricorsiva(const Fracture& frattura, list<Trace>& P_traces,  list<Trace>& NP_traces, PolygonalMesh& mesh, FracturesFunctions& fx, list<unsigned int>& lista_vert)
    {

        /*questa procedura sarà definita ricorsivamente:
            * prima divido la frattura in due in base alla prima traccia che dobbiamo considerare
            * vado a dividere le tracce restanti in base a quale poligono appartengono
            * if(traccia in poligono  == 0){salvo il poligono nella mesh}
            * else{richiamo questa funzione}

        *se sono entrata in questa procedura è perché c'è almeno una traccia associata alla frattura
        */

        //ricordiamo: prima le tracce passanti, poi quelle non passanti in ordine decrescente (le liste sono già ordinate)
        Trace traccia_tagliante;  // definisco l'oggetto traccia
        if (P_traces.size() != 0)  // verifico che il numero di tracce passanti sia diverso da 0
        {
            auto it = P_traces.begin();  // inizializzo la prima traccia, che è la più lunga passante
            traccia_tagliante = *it;
            P_traces.pop_front();  // la tolgo dalla lista di tracce passanti
        }
        else  // se non ho tracce passanti ( il numero di tracce non passanti è != 0 perché ho fatto un controllo per verificare la correttezza)
        {
            auto it = NP_traces.begin();   // prendo il primo elemento ( il più lungo)
            traccia_tagliante = *it;
            NP_traces.pop_front();         // elimino il primo elemento dalla lista
        }

        Vector3d estr1_tr_tagliante = traccia_tagliante.coordinates_extremes.col(0);   // dopo aver salvato l'oggetto traccia estraggo le coordinate degli estremi della traccia
        Vector3d estr2_tr_tagliante = traccia_tagliante.coordinates_extremes.col(1);
        MatrixXd retta_traccia_tagliante = fx.Retta_per_due_punti(estr1_tr_tagliante, estr2_tr_tagliante);         // calcolo la retta passante per gli estremi della traccia
        Vector3d dir_tr_tagliante = retta_traccia_tagliante.row(0);              // ricavo la direttrice
        Vector3d pt0_tr_tagliante = retta_traccia_tagliante.row(1);

        //ciclo sui lati della frattura per cercare i due punti di intersezione (devo gestire il caso di capotare su un vertice --> conta doppio)
        vector<Vector3d> nuovi_punti;
        vector<unsigned int> id_nuoviPunti;
        nuovi_punti.reserve(4); //male che vada ho 4 beta
        id_nuoviPunti.reserve(4);


        for (unsigned int i = 0; i< frattura.num_vertici; i++)  // ciclo sui lati della frattura per intersecarla con i lati
        {
            // definisco i segmenti del poligono
            unsigned int id_o = frattura.vertices[i];  // salvo l'ID del vertice di inizio del segmento
            unsigned int id_e;  // inizializzo l'ID del vertice di fine del segmento
            if(i==frattura.num_vertici-1){id_e = frattura.vertices[0];}  // effettuo le casistiche per capire qual è il vertice di fine
            else{id_e = frattura.vertices[i+1];}

            Vector3d origin = mesh.Cell0DCoordinates[id_o];   // determino le coordinate di origine segmento
            Vector3d end = mesh.Cell0DCoordinates[id_e];        // determino le coordinate di fine segmento
            MatrixXd retta_fratt = fx.Retta_per_due_punti(origin, end);  // costruisco la retta per questi due punti
            Vector3d dir_f = retta_fratt.row(0);        // calcolo la direttrice

            //interseco le rette se queste non sono parallele
            bool non_parallele = !((dir_f.cross(dir_tr_tagliante)).norm() < fx.tolleranza1D); // verifico che non sono parallele

            if(non_parallele) // se le rette non sono parallele si intersecano
            {
                Vector2d a_b = fx.intersezione_rette(retta_fratt, retta_traccia_tagliante);   // vettore contenente alpha e beta di intersezione tra le due rette
                //il punto deve stare sul segmento della frattura, ossia alpha appartiene a [0,1]
                if (a_b[0] > -fx.tolleranza1D && a_b[0] < 1+ fx.tolleranza1D){  // la posizione 0 corrisponde agli alpha, faccio un controllo su questi
                    Vector3d pto = pt0_tr_tagliante + a_b[1] * dir_tr_tagliante;  // determino il punto di intersezione prendendo il beta di intersezione e sostituendolo nell'equazione della retta
                    //se pto è già in nuovi_punti lo ignoro sennò faccio quello che devo
                    unsigned int a = 0;  // questa è l'inizializzazione del nuovo id da inserire
                    if(fx.pto_unico(pto, nuovi_punti, a))  // vedo il nuovo pto coincide con i punti di intersezione gia trovati, al passo 0 nuovi punti è vuoto quindi restituisco true ed entro nell'if
                    {

                        // struttura di supporto in cui ci salviamo i vertici della frattura
                        //devo gestire anche il caso in cui il punto coincide con un vertice della frattura
                        vector<Vector3d> coord_frc;
                        coord_frc.reserve(frattura.num_vertici);
                        for(unsigned int k = 0; k<frattura.num_vertici; k++)
                        {
                            coord_frc.push_back(mesh.Cell0DCoordinates[frattura.vertices[k]]);  // mi salvo le coordinate dei vertici del poligono
                        }
                        unsigned int id_new_pto = 0;

                        if (!fx.pto_unico(pto, coord_frc, id_new_pto)){
                            //caso in cui il prolungamento della traccia tagliante cade su un vertice della frattura che consideriamo
                            nuovi_punti.push_back(pto);  // a nuovo punto appendo il punto di intersezione
                            id_nuoviPunti.push_back(id_new_pto); // appendo anche l'ID nuovo che coincide con l'ID del vertice
                        }
                        else if(!fx.pto_unico(pto, mesh.Cell0DCoordinates, id_new_pto)){
                            //caso in cui il prolungamento della traccia tagliante tocca una cella 0D già creata
                            nuovi_punti.push_back(pto);  // appendo il punto di intersezione
                            id_nuoviPunti.push_back(id_new_pto);  // pusho dentro l'id dei nuovi punti

                            //inserisco l'id del punto nella lista prima del vertice di end del segmento
                            auto pos = find(lista_vert.begin(), lista_vert.end(), id_e);
                            lista_vert.insert(pos, id_new_pto);
                        }
                        else
                        {
                        nuovi_punti.push_back(pto);  // appendo il punto di intersezione
                        id_new_pto = mesh.NumberCell0D;  // Il su puo ID diventa equivalente al numerodi vertici nel poligono
                        id_nuoviPunti.push_back(id_new_pto);  // pusho dentro l'id dei nuovi punti
                        ReshapingArray::VerificaRaddoppio(mesh.Cell0DId);  // raddoppio la dimensione del vettore contenente le celle 0D se necessario
                        mesh.Cell0DId.push_back(id_new_pto);  // pusho i nuovi punti nel vettorre di celle 0d
                        ReshapingArray::VerificaRaddoppio(mesh.Cell0DCoordinates);  // verifico il raddoppio della dimensione del vettore
                        mesh.Cell0DCoordinates.push_back(pto);  // aggiundo le coordinate
                        mesh.NumberCell0D += 1;  // aggiungo 1 al numero di vertici

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
        frc1.id = frattura.id;
        frc1.vertices.reserve(10);  // riservo i vertici
        Fracture frc2;
        frc2.id = frattura.id;
        frc2.vertices.reserve(10);  // riserco i vertici

        //per suddividere i vertici nelle due fratture scorro la lista e avrò un booleano che switcho appena trovo uno dei due nuovi vertici
        bool flag = true; // true è per la frattura 1
        for(unsigned int elem : lista_vert)  // ciclo sulla lista dei vertici della frattura
        {
            if (elem == id_nuoviPunti[0] || elem == id_nuoviPunti[1])  // questo è il caso in cui vado a finire su uno dei nuovi punti della frattura, quindi lo assegno a entrambe le sotto fratture
            {
                //assegno l'elemento ad entrambe le fratture
                ReshapingArray::VerificaRaddoppio(frc1.vertices);
                frc1.vertices.push_back(elem);
                ReshapingArray::VerificaRaddoppio(frc1.vertices);
                frc2.vertices.push_back(elem);

                //switcho la flag
                flag = !flag;  // perché da qua in poi incontro gli elementi della seconda frattura
            }
            else
            {
                if (flag)  // se il frag è True allora assegno i punti alla prima frattura
                {
                    ReshapingArray::VerificaRaddoppio(frc1.vertices);
                    frc1.vertices.push_back(elem);

                }
                else // se il flag è false assegno gli elementi alla seconda frattura
                {
                    ReshapingArray::VerificaRaddoppio(frc2.vertices);
                    frc2.vertices.push_back(elem);

                }

            }

        }

        frc1.vertices.shrink_to_fit();  // riduco la capacità
        frc2.vertices.shrink_to_fit();
        frc1.num_vertici = frc1.vertices.size();   // ovviamente il numero di vertici diventa equivalente alla dimensione attuale del vettore
        frc2.num_vertici = frc2.vertices.size();

        list<unsigned int> list_vert1;
        for (unsigned int i = 0; i< frc1.num_vertici; i++){
            list_vert1.push_back(frc1.vertices[i]);
        }

        list<unsigned int> list_vert2;
        for (unsigned int i = 0; i< frc2.num_vertici; i++){
            list_vert2.push_back(frc2.vertices[i]);
        }

        //determino le tracce passanti e quelle non per ogni sotto-frattura
        list<Trace> P_traces1 = {};   // inizializzo la nuova lista di tracce per le mie fratture
        list<Trace> NP_traces1 = {};
        list<Trace> P_traces2 = {};
        list<Trace> NP_traces2  = {};



        // LAVORO PRIMA SULLE TRACCE PASSANTI
        if(P_traces.size() != 0){
            for (auto itor = P_traces.begin(); itor != P_traces.end(); itor++)  // ciclo sulle tracce
            {
                Trace traccia = *itor;

                // adesso devo vedere se le tracce dove sto ciclando sono dentro frac1 o frac2, considero anche il caso in cui appartengono a entrambi

               // passo 1: ricavo le coodinate degli estremi
                Vector3d origine = traccia.coordinates_extremes.col(0);
                Vector3d fine = traccia.coordinates_extremes.col(1);   // definisco le coordinate degli estremi della traccia

                // passo 2: devo riuscire a capire se la traccia appartiene alla prima o la seconda frattura (o entrambe)

                    // tutti i punti che hanno la stessa direzione e verso nel prodotto vettoriale con il p_0 della traccia tagliante appartengono ad una particolare frattura,
                    // basta quindi confrontare i segni dei prodotti vettoriali
                    // per ricavare un vettore direzzione da confrontare basta vedere sui vertici della della frattura 1. Ne scelgo uno, possibilmente con distanza da p0 > tolleranza, altrimenti faccio casino

                // vado a definire una segnatura (ovvero un vettore direzione)
                bool rimango_dentro_while = true;
                unsigned int i = 0;
                Vector3d segnatura (0,0,0);
                while (rimango_dentro_while && i<frc1.num_vertici)
                {
                    unsigned int id = frc1.vertices[i];
                    Vector3d vettore_provvisorio = mesh.Cell0DCoordinates[id];
                    if ((dir_tr_tagliante.cross(vettore_provvisorio-pt0_tr_tagliante)).norm()>fx.tolleranza1D)
                    {
                        segnatura = dir_tr_tagliante.cross(vettore_provvisorio-pt0_tr_tagliante);   // segnatura relativa alla frattura 1, tutti i punti sulla frattura 1 seguono questa segnatura
                        rimango_dentro_while = !rimango_dentro_while; // esco da while
                    }
                    i++;
                }
                // vediamo come sono segnati gli estremi delle tracce
                double segnatura_origine = ((dir_tr_tagliante).cross(origine-pt0_tr_tagliante)).dot(segnatura);
                double segnatura_fine = (dir_tr_tagliante).cross(fine-pt0_tr_tagliante).dot(segnatura);

                if ( (segnatura_origine > fx.tolleranza1D && segnatura_fine > fx.tolleranza1D) || ( abs(segnatura_origine) < fx.tolleranza1D  && segnatura_fine > fx.tolleranza1D) || (segnatura_origine > fx.tolleranza1D && abs(segnatura_origine) < fx.tolleranza1D))  // entrambe stessa direzione della segnatura della frattura 1
                {
                    P_traces1.push_back(traccia);
                }
                else if ( (segnatura_origine < - fx.tolleranza1D && segnatura_fine < - fx.tolleranza1D) || (abs(segnatura_origine) < fx.tolleranza1D && segnatura_fine < - fx.tolleranza1D) || (segnatura_origine < -fx.tolleranza1D && abs(segnatura_origine) < fx.tolleranza1D)) // tutte queste condizioni servono per vedere i casi delle traccie
                {
                    P_traces2.push_back(traccia);
                }
                else
                {
                    // devo creare due tracce da riassegnare ai sotto poligoni

                    // passo 1: cerco il punto di intersezione tra le due tracce
                    MatrixXd retta_traccia_da_smistare;
                    retta_traccia_da_smistare.resize(2,3);
                    retta_traccia_da_smistare = fx.Retta_per_due_punti(origine,fine);

                    Vector2d gamma = fx.intersezione_rette(retta_traccia_da_smistare,retta_traccia_tagliante);
                    Vector3d intersezione_tra_tracce = gamma(0)*retta_traccia_da_smistare.row(0)+retta_traccia_da_smistare.row(1);
                    // creo la matrice con le coordinate dei nuovi punti da salvare poi nella traccia
                    MatrixXd coordinate_traccia_per_frattura_1;
                    coordinate_traccia_per_frattura_1.resize(3,2);
                    MatrixXd coordinate_traccia_per_frattura_2;
                    coordinate_traccia_per_frattura_2.resize(3,2);
                    if (segnatura_origine > fx.tolleranza1D)
                    {
                        coordinate_traccia_per_frattura_1.col(0) = origine;
                        coordinate_traccia_per_frattura_1.col(1) = intersezione_tra_tracce;

                        coordinate_traccia_per_frattura_2.col(0) = intersezione_tra_tracce;
                        coordinate_traccia_per_frattura_2.col(1) = fine;
                    }
                    else
                    {
                        coordinate_traccia_per_frattura_1.col(0) = fine;
                        coordinate_traccia_per_frattura_1.col(1) = intersezione_tra_tracce;

                        coordinate_traccia_per_frattura_2.col(0) = intersezione_tra_tracce;
                        coordinate_traccia_per_frattura_2.col(1) = origine;
                    }
                    // passo 2: creo due oggetti traccia con lo stesso ID della traccia originale e le riassegno alle due sotto fratture
                    Trace traccia_per_frattura_1; // inizializzo la prima traccia
                    traccia_per_frattura_1.coordinates_extremes = coordinate_traccia_per_frattura_1;
                    traccia_per_frattura_1.id_frc1 = traccia.id_frc1;
                    traccia_per_frattura_1.id_frc2 = traccia.id_frc2;
                    traccia_per_frattura_1.len = traccia.len;
                    traccia_per_frattura_1.id = traccia.id;

                    Trace traccia_per_frattura_2; // inizializzo la prima traccia
                    traccia_per_frattura_2.coordinates_extremes = coordinate_traccia_per_frattura_2;
                    traccia_per_frattura_2.id_frc1 = traccia.id_frc1;
                    traccia_per_frattura_2.id_frc2 = traccia.id_frc2;
                    traccia_per_frattura_2.len = traccia.len;
                    traccia_per_frattura_2.id = traccia.id;


                    P_traces2.push_back(traccia_per_frattura_2);
                    P_traces1.push_back(traccia_per_frattura_1);
                }
            }
        }
        // LAVORO ORA SULLE TRACCE NON PASSANTI
        if(NP_traces.size() != 0){
            for (auto itor = NP_traces.begin(); itor != NP_traces.end(); itor++)  // ciclo sulle tracce
            {
                Trace traccia = *itor;

                // adesso devo vedere se le tracce dove sto ciclando sono dentro frac1 o frac2, considero anche il caso in cui appartengono a entrambi
                // passo 1: ricavo le coodinate degli estremi
                Vector3d origine = traccia.coordinates_extremes.col(0);
                Vector3d fine = traccia.coordinates_extremes.col(1);   // definisco le coordinate degli estremi della traccia

                // passo 2: devo riuscire a capire se la traccia appartiene alla prima o la seconda frattura (o entrambe)

                // tutti i punti che hanno la stessa direzione e verso nel prodotto vettoriale con il p_0 della traccia tagliante appartengono ad una particolare frattura,
                // basta quindi confrontare i segni dei prodotti vettoriali
                // per ricavare un vettore direzzione da confrontare basta vedere sui vertici della della frattura 1. Ne scelgo uno, possibilmente con distanza da p0 > tolleranza, altrimenti faccio casino

                // vado a definire una segnatura (ovvero un vettore direzione)
                bool rimango_dentro = true;
                unsigned int i = 0;
                Vector3d segnatura (0,0,0);
                while (rimango_dentro && i<frc1.num_vertici)
                {
                    unsigned int id = frc1.vertices[i];
                    Vector3d vettore_provvisorio = mesh.Cell0DCoordinates[id];
                    if ((dir_tr_tagliante.cross(vettore_provvisorio-pt0_tr_tagliante)).norm()>fx.tolleranza1D)
                    {
                        segnatura = dir_tr_tagliante.cross(vettore_provvisorio-pt0_tr_tagliante);   // segnatura relativa alla frattura 1, tutti i punti sulla frattura 1 seguono questa segnatura
                        rimango_dentro = !rimango_dentro; // esco da while
                    }
                    i++;
                }
                // vediamo come sono segnati gli estremi delle tracce
                double segnatura_origine = ((dir_tr_tagliante).cross(origine-pt0_tr_tagliante)).dot(segnatura);
                double segnatura_fine = ((dir_tr_tagliante).cross(fine-pt0_tr_tagliante)).dot(segnatura);

                if ( (segnatura_origine > fx.tolleranza1D && segnatura_fine > fx.tolleranza1D) || (segnatura_origine >= 0 && segnatura_origine < fx.tolleranza1D && segnatura_fine > fx.tolleranza1D) || (segnatura_origine > fx.tolleranza1D && segnatura_fine >= 0 && segnatura_fine < fx.tolleranza1D))  // entrambe stessa direzione della segnatura della frattura 1
                {
                    NP_traces1.push_back(traccia);
                }
                else if ( (segnatura_origine < - fx.tolleranza1D && segnatura_fine < - fx.tolleranza1D) || (segnatura_origine <= 0 && segnatura_origine < -fx.tolleranza1D && segnatura_fine < - fx.tolleranza1D) || (segnatura_origine < -fx.tolleranza1D && segnatura_fine <= 0 && segnatura_fine < -fx.tolleranza1D))   // tutte queste condizioni servono per vedere i casi delle traccie
                {
                    NP_traces2.push_back(traccia);
                }
                else
                {
                    MatrixXd retta_traccia_da_smistare;
                    retta_traccia_da_smistare.resize(2,3);
                    retta_traccia_da_smistare = fx.Retta_per_due_punti(origine,fine);

                    Vector2d gamma = fx.intersezione_rette(retta_traccia_da_smistare,retta_traccia_tagliante);
                    Vector3d intersezione_tra_tracce = gamma(0)*retta_traccia_da_smistare.row(0)+retta_traccia_da_smistare.row(1);
                    // creo la matrice con le coordinate dei nuovi punti da salvare poi nella traccia
                    MatrixXd coordinate_traccia_per_frattura_1;
                    coordinate_traccia_per_frattura_1.resize(3,2);
                    MatrixXd coordinate_traccia_per_frattura_2;
                    coordinate_traccia_per_frattura_2.resize(3,2);
                    if (segnatura_origine > fx.tolleranza1D)
                    {
                        coordinate_traccia_per_frattura_1.col(0) = origine;
                        coordinate_traccia_per_frattura_1.col(1) = intersezione_tra_tracce;

                        coordinate_traccia_per_frattura_2.col(0) = intersezione_tra_tracce;
                        coordinate_traccia_per_frattura_2.col(1) = fine;
                    }
                    else
                    {
                        coordinate_traccia_per_frattura_1.col(0) = fine;
                        coordinate_traccia_per_frattura_1.col(1) = intersezione_tra_tracce;

                        coordinate_traccia_per_frattura_2.col(0) = intersezione_tra_tracce;
                        coordinate_traccia_per_frattura_2.col(1) = origine;
                    }
                    // passo 2: creo due oggetti traccia con lo stesso ID della traccia originale e le riassegno alle due sotto fratture
                    Trace traccia_per_frattura_1; // inizializzo la prima traccia
                    traccia_per_frattura_1.coordinates_extremes = coordinate_traccia_per_frattura_1;
                    traccia_per_frattura_1.id_frc1 = traccia.id_frc1;
                    traccia_per_frattura_1.id_frc2 = traccia.id_frc2;
                    traccia_per_frattura_1.len = traccia.len;
                    traccia_per_frattura_1.id = traccia.id;

                    Trace traccia_per_frattura_2; // inizializzo la prima traccia
                    traccia_per_frattura_2.coordinates_extremes = coordinate_traccia_per_frattura_2;
                    traccia_per_frattura_2.id_frc1 = traccia.id_frc1;
                    traccia_per_frattura_2.id_frc2 = traccia.id_frc2;
                    traccia_per_frattura_2.len = traccia.len;
                    traccia_per_frattura_2.id = traccia.id;


                    NP_traces2.push_back(traccia_per_frattura_2);
                    NP_traces1.push_back(traccia_per_frattura_1);
                }
            }
        }

        //se non ho tracce interne interrompo la ricorsione e salvo la mesh
        if (P_traces1.size() + NP_traces1.size() == 0){
            //le celle 0D le ho salvate man mano

            //salvo le celle 1D e 2D
            unsigned int num_celle2d = mesh.NumberCell2D;
            mesh.NumberCell2D +=1;

            //mi creo due vettori ausiliari
            vector<unsigned int> edges;
            edges.reserve(frc1.num_vertici+1);
            edges = verifica_celle_1D(mesh, frc1);

            vector<unsigned int> vertices = frc1.vertices;
            mesh.NumberCell1D = mesh.Cell1DId.size();


            ReshapingArray::VerificaRaddoppio(mesh.Cell2DId);
            mesh.Cell2DId.push_back(num_celle2d);
            ReshapingArray::VerificaRaddoppio(mesh.Cell2DEdges);
            mesh.Cell2DEdges.push_back(edges);
            ReshapingArray::VerificaRaddoppio(mesh.Cell2DVertices);
            mesh.Cell2DVertices.push_back(vertices);



        }
        else{sottopoligonazione_ricorsiva(frc1, P_traces1, NP_traces1, mesh, fx, list_vert1);}

        if (P_traces2.size() + NP_traces2.size() == 0){
            //le celle 0D le ho salvate man mano

            //salvo le celle 1D e 2D
            unsigned int num_celle2d = mesh.NumberCell2D;
            mesh.NumberCell2D +=1;

            //mi creo due vettori ausiliari
            vector<unsigned int> edges;
            edges.reserve(frc2.num_vertici+1);
            edges = verifica_celle_1D(mesh, frc2);

            vector<unsigned int> vertices = frc2.vertices;
            mesh.NumberCell1D = mesh.Cell1DId.size();


            ReshapingArray::VerificaRaddoppio(mesh.Cell2DId);
            mesh.Cell2DId.push_back(num_celle2d);
            ReshapingArray::VerificaRaddoppio(mesh.Cell2DEdges);
            mesh.Cell2DEdges.push_back(edges);
            ReshapingArray::VerificaRaddoppio(mesh.Cell2DVertices);
            mesh.Cell2DVertices.push_back(vertices);

        }
        else{sottopoligonazione_ricorsiva(frc2, P_traces2, NP_traces2, mesh, fx, list_vert2);}

    }

    PolygonalMesh FracturesFunctions::creazione_mesh(const Fracture& frattura,  list<Trace>& P_traces,  list<Trace>& NP_traces, const vector<Vector3d>& coord)
    {
        PolygonalMesh mesh;    // inizializzo un oggetto mesh
        //ridefinisco la frattura in modo da iniziare l'indicizzazione delle celle 0D da 0
        Fracture frc;  // definisco un nuovo oggetto frattura (è ancora vuoto, lo devo riempire)
        vector<Vector3d> coord_frc; //conterrà solo i vertici di questa frattura
        frc.num_vertici = frattura.num_vertici; // questo nuovo oggetto avrà lo stesso numero di vertici della frattura che gli passo in input
        frc.vertices.reserve(frc.num_vertici);  // riservo la memoria per salvare gli ID di questa frattura
        coord_frc.reserve(frattura.num_vertici);  // riservo la memoria per salvare le coordinate di questa frattura
        frc.id = frattura.id;

        mesh.NumberCell0D += frc.num_vertici;  // il numero di vertici è lo stesso della frattura frc
        mesh.Cell0DId.reserve(50);      // riservo 50 (numero molto grande) perché non so quanti ID mi generano i tagli
        mesh.Cell0DCoordinates.reserve(50);   // non so nemmeno quante coordinate nuove mi generano i tagli per questo riservo 50

        //inizio già a salvare le coordinate 0D
        for (unsigned int i = 0; i< frattura.num_vertici; i++)  // ciclo sul numero di vertici della frattura che passo in input
        {
            Vector3d pto = coord[frattura.vertices[i]];  // di questa frattura mi prendo le coordinate dei punti
            coord_frc.push_back(pto);    // le salvo nella matrice coord_frc
            frc.vertices.push_back(i);  // appendo i vertici all'oggetto frattura
            mesh.Cell0DCoordinates.push_back(pto);  // aggiungo le coordinate del punto alla mesh
            mesh.Cell0DId.push_back(i);    // aggiungo l'ID del punto alla mesh
        }

        if (P_traces.size() + NP_traces.size() == 0)
        {
            //salvo le celle 1D
            mesh.NumberCell1D += frattura.num_vertici;   // il numero di vertici è equivalente al numero di lati
            mesh.Cell1DId.reserve(frattura.num_vertici);  // riservo la memoria per salvare i lati
            for (unsigned int v = 0; v < frattura.num_vertici; v++)  // ciclo sul numero di lati
            {
                unsigned int id_origin = frattura.vertices[v];          // salvo gli ID dei vertici di origine e fine
                unsigned int id_end;
                if(v == frattura.num_vertici -1){id_end = frattura.vertices[0];} // il poligono è chiuso => l'ultimo vertice è collegato tramite un segmento al primo
                else{id_end = frattura.vertices[v+1];}
                mesh.Cell1DId.push_back(v);         // aggiungo v agli ID dei segmenti
                mesh.Cell1DVertices.push_back({id_origin, id_end}); // aggiungo origine e fine
            }

            //salvo le celle 2D --> in questo caso ne ho solo una, la frattura stessa
            mesh.NumberCell2D = 1;   // perché ne ho una
            mesh.Cell2DId.reserve(1);  // riservo una cella di memoria
            mesh.Cell2DId.push_back(0); // inizializzo (Renato chiedi a Martina)
            mesh.Cell2DVertices.reserve(1);  // riservo una cella di memoria
            mesh.Cell2DVertices.push_back(frattura.vertices);  // inserisco gli ID dei vertici
            mesh.Cell2DEdges.reserve(1);  // riservo una cella di memoria
            mesh.Cell2DEdges.push_back(mesh.Cell1DId);  // aggiungo gli ID dei segmenti
        }
        else{

            //mi creo una lista di vertici per poterli ordinare in maniera più easy
            list<unsigned int> vert_fract;
            for (unsigned int i = 0; i< frc.num_vertici; i++)
            {
                vert_fract.push_back(frc.vertices[i]);   // in questa lista pusho tutti i vertici della frattura frc
            }

            sottopoligonazione_ricorsiva(frc, P_traces, NP_traces, mesh, *this, vert_fract);  // do tutto in pasto alla funzione che fa la sotto poligonazione
            // usiamo this per avere le funzioni di FracturesFunctions legato alla tolleranza
        }


        for(unsigned int id_cella_2d=0; id_cella_2d<mesh.NumberCell2D; id_cella_2d ++){
            vector<unsigned int> vertici=mesh.Cell2DVertices[id_cella_2d];
            vector<unsigned int> lati=mesh.Cell2DEdges[id_cella_2d];

            for(auto itor = mesh.Cell2DEdges[id_cella_2d].begin(); itor != mesh.Cell2DEdges[id_cella_2d].end(); itor++ ){
                unsigned int id_lato= (*itor);
                Vector2i id_estremi_lato = mesh.Cell1DVertices(id_lato);
                (*this).Retta_per_due_punti(mesh.Cell0DCoordinates[id_estremi_lato[0]], );

                vector<array<double,2>> punti_interni;
                punti_interni.reserve(50);




                //cicliamo sulle celle 0d e cerchiamo i punti interni al segmento id_lato
                for (unsigned int k = 0 ; k < mesh.NumberCell0D; k++){

                }



            }
        }
        return mesh;

    }


} //namespace
