#ifndef __RESHAPING_H // Header guards
#define __RESHAPING_H

#include <vector>
#include "Eigen/Eigen"

using namespace std;
using namespace Eigen;

namespace ReshapingArray{

//questo header conterrà le funzioni necessarie alla gestione della dimensione di array dinamici
template<typename T>
void VerificaRaddoppio( vector<T>& vec){

    //vec è un vettore di lunghezza d con n elementi
    unsigned int d = vec.capacity();
    unsigned int n = vec.size();

    if( n == d) {
        vector<T> new_vec;
        new_vec.reserve(2*d);
        for (unsigned int i = 0; i < n ; i++){new_vec.push_back(vec[i]);}
        vec = new_vec;
    }
}

template<typename T>
void VerificaDimezzamento( vector<T>& vec){

    //vec è un vettore di lunghezza d con n elementi
    unsigned int d = vec.capacity();
    unsigned int n = vec.size();

    if ( d > 1  && n == d/4){
        vector<T> new_vec;
        new_vec.reserve(0.5*d);
        for (unsigned int i = 0; i < n; i++){new_vec.push_back(vec[i]);}
        vec = new_vec;
    }
}


} //namespace
#endif
