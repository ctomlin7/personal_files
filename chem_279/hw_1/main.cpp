#include "hw_1.h"

int main(void) {
    string filename = "/home/ctomlin/chem_279/hw_1/chem_x79_hw1/sample_input/Force/1.txt";
    Atoms result(filename);
    vector<Vector3d> vec = result.read_xyz(filename);

    LennardJones lj(filename);

    //lj.run_LJ(filename);

    //cout << vec << endl;

    FiniteDifference fd(filename);

    fd.run_FD(filename);

    return 0;
}