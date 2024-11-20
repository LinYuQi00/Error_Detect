#include <iostream>
#include "ErrorDetection.h"

using namespace std;


int main(){
    string ioe = "../Data/PRECAI.DAT";
    string gcp = "../Data/PRECKI.DAT";
    string pixel = "../Data/PREPHIy10.DAT";
    string output = "../Data/result.txt";

    ErrorDetection test;
    test.process(ioe, gcp, pixel, output);

    return 0;
}

