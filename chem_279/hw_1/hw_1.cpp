#include "hw_1.h"

Atoms::Atoms(const string &filename) 
        :filename(filename)
        { };

LennardJones::LennardJones(const string &filename) 
        :Atoms(filename)
        { };

FiniteDifference::FiniteDifference(const string &filename) 
        :LennardJones(filename)
        { };

vector<vector<double>> Atoms::read_xyz(const string &filename) {
    string firstline;
    string line;
    int num_atoms;
    ifstream file(filename);
    if (!file.is_open()) {
        throw runtime_error("File not found");
    }

    if (file.is_open()) {
        if (getline(file, firstline)) {
            istringstream iss(firstline);
            iss >> num_atoms;
        }
        vector<vector<double>> result;
        for (int i = 0; i < num_atoms; i++) {
            if (!getline(file, line)) { throw runtime_error("Unexpected end of file."); }

            istringstream iss(line);

            int atom;
            double x, y, z;
            iss >> atom >> x >> y >> z;

            if (!(atom == 79)) { throw runtime_error("This is not a gold atom"); }

            vector<double> v = {x, y, z};
            result.push_back(v);
        }
        file.close();
        return result;
    }
    throw runtime_error("File not found");
}

double LennardJones::calculate_LJ(double r_ij) {
    double inv_r_ij = sigma / r_ij;
    double r6_term = pow(inv_r_ij, 6.0);
    double r12_term = r6_term * r6_term;
    double pairwise_energy = epsilon * (r12_term - (2 * r6_term));
    return pairwise_energy;
}

double LennardJones::calculate_distance(AtomCoord coord1, AtomCoord coord2) {
    double distance = 0;

    for (int i = 0; i < 3; i++) {
        double dim_dist = (coord1[i] - coord2[i]);
        dim_dist = dim_dist * dim_dist;
        distance += dim_dist;
    }
    distance = sqrt(distance);
    return distance;
}

double LennardJones::calculate_total_energy(Coordinates coordinates) {
    double total_energy = 0.0;

    for (int i = 0; i < coordinates.size(); i++) {
        for (int j = i + 1; j < coordinates.size(); j++) {
            double dist_ij = calculate_distance(coordinates[i], coordinates[j]);
            total_energy += calculate_LJ(dist_ij);
        }
    }
    return total_energy;
}

double LennardJones::calculate_pair_energy(Coordinates coordinates, int i_particle) {
    double e_total = 0.0;
    int num_atoms = coordinates.size();

    for (int j_particle = 0; j_particle < num_atoms; j_particle++) {
        if (j_particle != i_particle) {
            double dist_ij = calculate_distance(coordinates[i_particle], coordinates[j_particle]);
            double interaction_energy = calculate_LJ(dist_ij);
            e_total += interaction_energy;
        }
    }
    return e_total;
}

void LennardJones::run_LJ(string &filename) {
    Atoms result(filename);
    vector<vector<double>> vec = result.read_xyz(filename);
    LennardJones lj(filename);
    
    double total_energy = lj.calculate_total_energy(vec);

    cout << total_energy << endl;
}

vector<vector<double>> FiniteDifference::analytical_force(Coordinates coordinates) {
    int num_atoms = coordinates.size();
    Coordinates force_vec(num_atoms, vector<double>(3, 0.0));

    for (int i = 0; i < num_atoms; ++i) {
        for (int j = i + 1; j < num_atoms; ++j) {
            vector<double> r_vec(3, 0.0);
            double r_ij = 0.0;

            for (int k = 0; k < 3; ++k) {
                double dim_dist = (coordinates[i][k] - coordinates[j][k]);
                r_vec[k] = dim_dist;
                r_ij += dim_dist * dim_dist;
            }
            r_ij = sqrt(r_ij);

            double r7_term = (pow(sigma, 6.0) / pow(r_ij, 7.0)); 
            double r13_term = (pow(sigma, 12.0) / pow(r_ij, 13.0));
            double force_mag = epsilon * ((12 * r13_term) - (12 * r7_term));

            for (int k = 0; k < 3; ++k) {
                double force_component = -1 * force_mag * r_vec[k] / r_ij;
                force_vec[i][k] += force_component;
                force_vec[j][k] -= force_component;  
            }
        }
    }
    return force_vec;
}

double FiniteDifference::forward_difference(Coordinates coordinates, int i_particle, int dim, double h) {
    double e_total = calculate_pair_energy(coordinates, i_particle);
    double e_plus = 0.0;

    Coordinates coordinates_plus = coordinates;

    coordinates_plus[i_particle][dim] += h;

    e_plus = calculate_pair_energy(coordinates_plus, i_particle);

    double derivative = -1 * (e_plus - e_total) / h;
    return derivative;
}

double FiniteDifference::central_difference(Coordinates coordinates, int i_particle, int dim, double h) {
    double e_total = calculate_pair_energy(coordinates, i_particle);
    double e_plus = 0.0;
    double e_minus = 0.0;

    Coordinates coordinates_plus = coordinates;
    Coordinates coordinates_minus = coordinates;

    coordinates_plus[i_particle][dim] += h;
    coordinates_minus[i_particle][dim] -= h;

    e_plus = calculate_pair_energy(coordinates_plus, i_particle);
    e_minus = calculate_pair_energy(coordinates_minus, i_particle);

    double derivative = -1 * (e_plus - e_minus) / (2 * h);
    return derivative;
}

void FiniteDifference::run_FD(string &filename) {
    Atoms result(filename);
    Coordinates vec = result.read_xyz(filename);
    LennardJones lj(filename);
    FiniteDifference fd(filename);

    cout << "E_LJ: " << lj.calculate_total_energy(vec) << endl;
    
    printFormattedData("F_LJ_Analytical", analytical_force(vec), 0.0);

    for (int i = 0; i < h.size(); i++) {
        Coordinates F_LJ_forward_difference;
        Coordinates F_LJ_central_difference;
        for (int k = 0; k < vec.size(); k++) {
            vector<double> for_row;
            vector<double> cen_row;
            for (int j = 0; j < 3; j++) {
                double forward_derivative = fd.forward_difference(vec, k, j, h[i]);
                double central_derivative = fd.central_difference(vec, k, j, h[i]);

                for_row.push_back(forward_derivative);
                cen_row.push_back(central_derivative);
            }
            F_LJ_forward_difference.push_back(for_row);
            F_LJ_central_difference.push_back(cen_row);
        }
        printFormattedData("F_LJ_forward_difference", F_LJ_forward_difference, h[i]);
        printFormattedData("F_LJ_central_difference", F_LJ_central_difference, h[i]);
    }
}

void FiniteDifference::printFormattedData(const string & label, Coordinates data, double stepsize) {
    cout << label << "   Stepsize: " << setprecision(4) << defaultfloat << stepsize << endl;
    
    for (const auto & row : data) {
        cout << "(";
        for (size_t i = 0; i < row.size(); i++) {
            cout << setw(12) << setprecision(4) << scientific << row[i];
            
            if ((i + 1) % 3 == 0) {
                cout << " )";
            }
        }
        cout << endl;
    }
    cout.flush();
}

void FiniteDifference::steepest_descent(Coordinates coordinates) {
    Golden golden_search;
    double displacement;
    const int max_iteration = 100;
    int iteration = 0;

    cout << "Initial Energy: " << calculate_total_energy(coordinates) << endl;

    Coordinates new_coords = coordinates;
    for (int i = 0; i < coordinates.size(); ++i) {
        for (int j = 0; j < 3; ++j) {
            new_coords[i][j] -= displacement * central_difference(coordinates, i, j, 1e-4);
        }
    }
    printFormattedData("Central Difference Force", new_coords, 0.0);

    while (iteration < max_iteration) {
        Coordinates gradient = analytical_force(coordinates);
        auto line_search_function = [this, &coordinates, &gradient](double displacement) {
            Coordinates new_coords = coordinates;
            for (int i = 0; i < coordinates.size(); ++i) {
                for (int j = 0; j < 3; ++j) {
                    new_coords[i][j] -= displacement * gradient[i][j];
                }
            }
            return calculate_total_energy(new_coords);
        };

        golden_search.bracket(0.0, 1.0, line_search_function);
        displacement = golden_search.minimize(line_search_function);

        for (int i = 0; i < coordinates.size(); ++i) {
            for (int j = 0; j < 3; ++j) {
                coordinates[i][j] -= displacement * gradient[i][j];
            }
        }

        double energy = calculate_total_energy(coordinates);
        cout << "Iteration " << iteration << ": Energy = " << energy << endl;

        if (displacement < 1e-6) {
            break;
        }

        ++iteration;
    }

    cout << "Final Energy: " << calculate_total_energy(coordinates) << endl;
}

ostream & operator<<(ostream & os, const vector<vector<double>> & vec) {
    for (const auto& v : vec) {
        os << "79 (" << v[0] << ", " << v[1] << ", " << v[2] << ")" << std::endl;
    }
    return os;
}