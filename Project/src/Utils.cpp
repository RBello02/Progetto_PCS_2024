#include "Utils.hpp"
#include "reshaping_array.hpp"
#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>


#include <cmath>

#include "Eigen/Eigen"

using namespace Eigen;

namespace GeometryLibrary{

//funzioni inline di supporto
inline Matrix3d Coef_piano(const Fractures& frc, unsigned int id_fract1)
{
    // salvo il piano in forma parametrica X = P0 + a1(P2-P0) + a2(P1-P0), dove a1 e a2 sono le parametrizzazioni
    // il piano è quindi identificato da 3 punti P0,P1,P2

    Matrix3d A;
    Vector3d P0 = frc.coordinates[frc.vertices_fractures[id_fract1][0]];
    Vector3d P1 = frc.coordinates[frc.vertices_fractures[id_fract1][1]];
    Vector3d P2 = frc.coordinates[frc.vertices_fractures[id_fract1][2]];

    A.col(0) = P0;
    A.col(1) = P2-P0;
    A.col(2) = P1-P0;


    return A.transpose(); // restituisce una matrice con la seguente struttura   (P0;P1;P2)  dove Pi sono VETTORI RIGA
}

inline bool Parallelismo(const Matrix3d& piano_1, const Matrix3d& piano_2)
{
    // questa funzione verifica se i piani sono paralleli

    bool risultato = false;

    // NOMENCLATURA
    // u1 = piano_1.row(1) v1 = piano_1.row(2) sono i vettori direttori del primo piano
    // u2 = piano_2.row(1) v2 = piano_2.row(2) sono i vettori direttori del secondo piano
    // P0_1 = piano_1.row(0) punto iniziale per il piano 1
    // P0_2 = piano_2.row(0) punto iniziale per il piano 2

    // passo 1: ricavare la normale ai piani
    Vector3d n1 = piano_1.row(1).cross(piano_1.row(2));
    Vector3d n2 = piano_2.row(1).cross(piano_2.row(2));

    // passo 2 : vedo se il prodotto vettoriale tra le normali è 0
    Vector3d t = n1.cross(n2);
    if (t.norm() == 0)
    {
        risultato = !risultato;
    }
    return risultato;
}

inline MatrixXd Retta_tra_piani(const Matrix3d& piano_1, const Matrix3d& piano_2)
{
    // calcolo la retta di intersezione tra i piani in forma parametrica
    // X = at+P


    // NOMENCLATURA
    // u1 = piano_1.row(1) v1 = piano_1.row(2) sono i vettori direttori del primo piano
    // u2 = piano_2.row(1) v2 = piano_2.row(2) sono i vettori direttori del secondo piano
    // P0_1 = piano_1.row(0) punto iniziale per il piano 1
    // P0_2 = piano_2.row(0) punto iniziale per il piano 2

    // passo 1: ricavare la normale ai piani
    Vector3d n1 = piano_1.row(1).cross(piano_1.row(2))/(piano_1.row(1).norm()*piano_1.row(2).norm());
    Vector3d n2 = piano_2.row(1).cross(piano_2.row(2))/(piano_2.row(1).norm()*piano_2.row(2).norm());

    // passo 2 : ricavare la direttrice della retta
    Vector3d t = n1.cross(n2);

    // passo 3: ricavare lo scalare d per i due piani
    double d1 = n1.dot(piano_1.row(0));
    double d2 = n2.dot(piano_2.row(0));
    Vector3d d = {d1,d2,0};  // salvo queste quantita in un vettore

    // passo 4: ricavare la matrice A
    Matrix3d A;
    A << n1,n2,t;

    // passo 5: risolvere il sistema per ricavare P
    Vector3d P = A.partialPivLu().solve(d);

    // salvo in un formato particolare
    MatrixXd X;
    X.resize(2,3);
    X.row(0) = t.transpose();
    X.row(1) = P.transpose();  // come una matrice 2x3

    return X;

}

inline MatrixXd Retta_per_due_vertici_della_frattura(const Fractures& frc, unsigned int id_vertice1, unsigned int id_vertice2)
{
    // l'equazione parametrica è X = at+P

    // passo 1: salvo le coordinate dei due punti
    double x1 = frc.coordinates[id_vertice1][0];
    double y1 = frc.coordinates[id_vertice1][1];
    double z1 = frc.coordinates[id_vertice1][2];

    double x2 = frc.coordinates[id_vertice2][0];
    double y2 = frc.coordinates[id_vertice2][1];
    double z2 = frc.coordinates[id_vertice2][2];

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

inline Vector2d alpha_di_intersezione(MatrixXd r_frattura, MatrixXd retta_intersez)
{

    //imposto un sistema lineare per la ricerca dei parametri alpha e beta
    //primo parametro è la matrice della retta del poligono
    Vector3d t1 = r_frattura.row(0).transpose();
    Vector3d P1 = r_frattura.row(1).transpose();

    //secondo parametro è la matrice della retta di intersezione tra i piani
    Vector3d t2 = retta_intersez.row(0).transpose();
    Vector3d P2 = retta_intersez.row(1).transpose();

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


/****************************************************************************************************************/

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

bool IntersectionFractures(Fractures& frc, unsigned int id_fract1, unsigned int id_fract2){
    Matrix3d piano_frc1 = Coef_piano(frc,id_fract1);
    Matrix3d piano_frc2 = Coef_piano(frc,id_fract2);

    //procedo nella funzione solo nel caso in cui i due piani non siano paralleli
    if(Parallelismo(piano_frc1, piano_frc2)){
        cout << "La frattura " << id_fract1 << " e la fratt " << id_fract2 << " non si intersecano" << endl;
        return false;}
    else{
        //calcolo la retta di intersezione tra i piani in forma parametrica
        MatrixXd retta_intersez_piani = Retta_tra_piani(piano_frc1,piano_frc2);
        Vector3d dir_retta_intersez_piani;
        dir_retta_intersez_piani = retta_intersez_piani.row(0);

        //Devo ora ciclare sulle coppie di vertici delle due fratture per trovare gli alpha di intersezione tra la retta di
        //intersezione tra i due piani determinati dalle fratture e le rette determinate dai vertici
        array<vector<double>,2> beta_inters; //--> nel primo array avrò gli alpha della prima frc e nel secondo della seconda
        (beta_inters[0]).reserve(2);
        (beta_inters[1]).reserve(2);

        unsigned int num_vertici_frc1 = frc.dim_fractures[id_fract1];
        unsigned int num_vertici_frc2 = frc.dim_fractures[id_fract2];

        //Introduco un contatore che incremento ogni volta che la retta di intersezione tra i due piani interseca un
        // segnemnto della frattura. se alla fine avrò cont == 4, vuol dire che ho una traccia
        unsigned int cont = 0;

        for (unsigned int i = 0; i < num_vertici_frc1; i++){

            unsigned int j;
            if (i == num_vertici_frc1-1){j = 0;}
            else{j = i+1;}

            //mi trovo la retta per ogni coppia di vertici
            MatrixXd retta_tra_vertici = Retta_per_due_vertici_della_frattura(frc, frc.vertices_fractures[id_fract1][i],frc.vertices_fractures[id_fract1][j]);
            Vector3d dir_retta_tra_vertici;
            dir_retta_tra_vertici = retta_tra_vertici.row(0);

            //cerco l'eventuale intersezione tra questa retta e la retta di intrsezione tra i piani, queste esiste se le due rette non
            //sono parallele
            if (!((dir_retta_tra_vertici.cross(dir_retta_intersez_piani)).norm() ==0)){
                Vector2d a_b = alpha_di_intersezione(retta_tra_vertici, retta_intersez_piani);

                if (a_b[0] >= -pow(10,-10) && a_b[0] <= 1+pow(10, -10)){
                    cont += 1;
                    (beta_inters[0]).push_back(a_b[1]);
                }

            }
        }

        for (unsigned int i = 0; i < num_vertici_frc2; i++){

            unsigned int j;
            if (i == num_vertici_frc2-1){j = 0;}
            else{j = i+1;}

            //mi trovo la retta per ogni coppia di vertici
            MatrixXd retta_tra_vertici = Retta_per_due_vertici_della_frattura(frc, frc.vertices_fractures[id_fract2][i],frc.vertices_fractures[id_fract2][j]);
            Vector3d dir_retta_tra_vertici;
            dir_retta_tra_vertici = retta_tra_vertici.row(0);

            //cerco l'eventuale intersezione tra questa retta e la retta di intrsezione tra i piani, queste esiste se le due rette non
            //sono parallele
            if (!((dir_retta_tra_vertici.cross(dir_retta_intersez_piani)).norm() ==0)){
                Vector2d a_b = alpha_di_intersezione(retta_tra_vertici, retta_intersez_piani);
                if (a_b[0] >= -pow(10,-10) && a_b[0] <= 1+pow(10, -10)){
                    cont += 1;
                    (beta_inters[1]).push_back(a_b[1]);
                }
                }
            }

        cout << "id frc: " << id_fract1 << ", alpha: ";
        for (const auto& elemento : beta_inters[0]) {
            cout << elemento << " ";
        }
        cout << endl;

        cout << "id frc: " << id_fract2 << ", alpha: ";
        for (const auto& elemento : beta_inters[1]) {
            cout << elemento << " ";
        }
        cout << endl;


        if (cont == 4){cout << "Ho una traccia tra la frattura " << id_fract1 << " e la fratt " << id_fract2 << endl; }


    } //chiude l'else
    cout << endl;
    return true;
}




} //namespace
