#include "hw_1.h"

int main(void) {
    string filename = "/home/ctomlin/chem_279/hw_1/chem_x79_hw1/sample_input/SD_with_line_search/2.txt";
    Atoms result(filename);
    vector<vector<double>> vec = result.read_xyz(filename);

    LennardJones lj(filename);
    lj.run_LJ(filename);
    cout << vec << endl;

    FiniteDifference fd(filename);
    fd.run_FD(filename);

    fd.steepest_descent(vec);

    return 0;
}