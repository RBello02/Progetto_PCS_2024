#include "src/GeometryLibrary.hpp"
#include "src/Utils.hpp"

#include <string>

using namespace GeometryLibrary;

int main()
{
    string path = "DFN";
    string finame = path + "/FR3_data.txt";
    Fractures frc;

    if (!importData(finame, frc)){return 1;}

  return 0;
}
