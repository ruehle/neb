#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <algorithm>
#include <iterator>
#include <ctime>
#include <random>


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
        for (size_t k = 0; k < xi.size(); ++k) {
            double f = double(i) / (nimages-1);
            x[k] = xi[k] + f * (xi[k] - xf[k]);
            path.push_back(x.copy());
        }
    }
    path.push_back(xf.copy());
    return path;
}

int main()
{
    auto v = vector_from_file("../lj6_m1");
    Array<double> xi = Array<double>(v).copy();
    v = vector_from_file("../lj6_m2");
    Array<double> xf = Array<double>(v).copy();

    cout << "hi\n";
    cout << xf << std::endl;

    pele::LJ lj(4., 4.);
    pele::CartesianNEBDistance neb_dist;

    pele::NEB neb(&lj, &neb_dist);

    auto path = linear_interpoloation(xi, xf, 5);
    neb.set_path(path);

    neb.start_with_lbfgs(.1, 10, 1., .01);
    neb.step();

}
