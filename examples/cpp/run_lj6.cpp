#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <random>
#include <sstream>


#include "lj.h"
#include "array.h"
#include "distance.h"
#include "neb.h"

using pele::Array;
using std::cout;

std::vector<double> vector_from_file(std::string fname)
{
    std::vector<double> x;
    std::ifstream fin;
    fin.open(fname.c_str());
    // Prepare a pair of iterators to read the data from cin
    std::istream_iterator<double> eos;
    std::istream_iterator<double> input(fin);
    // No loop is necessary, because you can use copy()
    std::copy(input, eos, std::back_inserter(x));

    fin.close();
    std::cout << "size " << x.size() << std::endl;
    return x;
}

std::vector< Array<double> > linear_interpoloation(Array<double> xi, Array<double> xf, size_t nimages)
{
    std::vector< Array<double> > path;
    path.push_back(xi.copy());
    Array<double> x(xi.size());

    for (size_t i = 1; i < nimages-1; ++i) {
        double f = double(i) / (nimages-1);
        for (size_t k = 0; k < xi.size(); ++k) {
            x[k] = xi[k] + f * (xf[k] - xi[k]);
        }
        path.push_back(x.copy());
    }
    path.push_back(xf.copy());
    return path;
}

void print_EofS(Array<double> energies,
        Array<double> distances,
        size_t count)
{
    std::stringstream ss;
    ss << "neb.EofS." << count;

    std::ofstream fout(ss.str());
    double cum_dist = 0;
    for (size_t i = 0; i < energies.size(); ++i) {
        fout << cum_dist << " " << energies[i] << "\n";
        cum_dist += distances[i];
    }
    fout.close();
}

int main()
{
    auto v = vector_from_file("../lj13_m1");
    Array<double> xi = Array<double>(v).copy();
    v = vector_from_file("../lj13_m2");
    Array<double> xf = Array<double>(v).copy();

    cout << "hi\n";
    cout << xf.size() << std::endl;


    pele::LJ lj(4., 4.);
    pele::CartesianNEBDistance neb_dist;

    cout << "energy1 " << lj.get_energy(xi) << std::endl;
    cout << "energy2 " << lj.get_energy(xf) << std::endl;

    pele::NEB neb(&lj, &neb_dist);

    auto path = linear_interpoloation(xi, xf, 30);
    neb.set_path(path);
    neb.set_k(100.);
    neb.set_double_nudging(true);
    neb.set_verbosity(10);

//    for (auto x : path) {
//        cout << "size " << x.size() << std::endl;
//        cout << "energy_path " << lj.get_energy(x) << std::endl;
//    }

    neb.start_with_lbfgs(1e-3, 10, 1., .01);

    for (size_t i = 0; i < 500; ++i) {
        bool success = neb.step();


        if (i % 10 == 0.) {
            print_EofS(neb.get_true_energies(), neb.get_distancese(), i);
        }

        if (i % 5 == 0.) {
            neb.adjust_k();
        }


        if (success) break;
    }


}
