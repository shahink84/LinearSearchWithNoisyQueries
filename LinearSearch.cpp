/******************************************************************************************
*******************************************************************************************
********************* Linear Search with Noisy Queries  ***********************************
*******************************************************************************************
*******************************************************************************************/


// LinearSearch.cpp : This file contains the 'main' function. Program execution begins and ends there.
//

using namespace std;

#include <iostream>
#include <fstream>
#include "ErroneousBinarySearch.h"
#include <string>
#include <ctime>
#include <cstdlib>
#include <algorithm>
#include <chrono> 


/* the cost of an algorithm which explores distances x, x.y, x.y^2, x, y^3, ... until it finds the 
target at position dist on the right of the origin*/
double get_gemetric_alg_cost(int dist, double x, double y)
{
    bool found = false;
    double cost = 0;
    double jump = x;
    while (!found)
    {
        // plan to move jump units to the right
        if (jump >= dist)
        {
            cost += dist;
            found = true;
        }
        else
        {
            cost += 2 * jump;
            jump = jump * y; // move new jump units to the left and back
            cost += 2 * jump;
            jump = jump * y; // the updated jump for the next effort
        }
    }
    return cost;
}

// get the value of b (exponent), as given by Equation 6 in the paper
double get_b (int k, double rho, double zeta2)
{
    double b = -1;
    double l = (double)(1 << k); // 2^k
    if (rho <= (1 + l) * (1 + l) / l)
        b = pow(zeta2, 1.0 / l);
    else
        b = pow(1 + l, 1 / l);
    return b;
}


/* Given a target at distance dist on the right of origin, return the index of the best strategy
among 2^k strategies X0,X1,...x_{2^k-1} for this location of the target. 
Here X_i = b^i, b^{i+l}, b^{i+2l}, b^{i+3l}, ...*/

int find_best_alg_from_family(int dist, int k, double b)
{
    int best_indx = -1;
    int best_cost = INT_MAX;
    int l = 1 << k; // 2^k
    for (int i = 0; i < l; i++)
    {
        int this_alg_cost = get_gemetric_alg_cost(dist, pow(b, i), pow(b, l));
        if (this_alg_cost < best_cost)
        {
            best_cost = this_alg_cost;
            best_indx = i;
        }
    }
//    cout << "best_indx, best_ratio: " << best_indx << "\t" << (double)best_cost/dist << endl;
    return best_indx;
}

/* Let eta be a uniform random number in tau/2, tau. 
Return a boolean vector of size k with eta 0s at random positions (and 1 in other indices). */
bool* create_response_vector(int k, int tau)
{
    // generate a random permutation to define the flipped bits
    int* random_perm = new int[k];
    for (int i = 0; i < k; i++)
        random_perm[i] = i;
    random_shuffle(random_perm, random_perm + k - 1);

    // array of answers
    bool* ans = new bool[k];
    for (int i = 0; i < k; i++)
        ans[i] = true;

    // eta is a random number in {tau/2, tau/2+1, ..., tau}.
    int eta = -1;
    
    eta = tau; // (rand() % ((tau + 1) / 2 + 1)) + tau / 2;
    // cout << "eta is: " << eta << endl;
    for (int i = 0; i < eta; i++)
        ans[random_perm[i]] = false;
  /*  for (int i = 0; i < k; i++)
        cout << ans[i];
    cout << endl;//*/

    delete random_perm;
    return ans;
}

// showing progress while the test is running
void show_progress(double progress)
{
    int barWidth = 70;
    std::cout << "[";
    int pos = barWidth * progress;
    for (int i = 0; i < barWidth; ++i) {
        if (i < pos) std::cout << "=";
        else if (i == pos) std::cout << ">";
        else std::cout << " ";
    }
    std::cout << "] " << int(progress * 100.0) << " %\r";
    std::cout.flush();
}

/* 
The main test is performed in the following method:
Input:
k: the number of queries (e.g., k = 16)
r: the desired robustness (a real value, e.g., r = 15.00)
tau: the upper bound for the number of incorrect answers
max_tau: maximum value of tau; tau in [0..max-tau] will be checked (max_tau <= k/4, e.g., tau = 4)\n
max_distance: distance range from the origin (e.g., max-distance = 100000)
n_test: number of tests when taking average (e.g., n_test = 100)
out_file: the file where the results are recorded.

OUTPUT:
The experimental performance ratio (competitive ratio) when running 
Robust Binary Interval Search (RBIS) 
using k queries on r-robust algorithms when tau queries are incorrectly answered, 
for tau in [0..maxtau]. 
For each value of tau, n-tests are made. 
The average/max/min competitive ratio is reported over n-tests.
*/

void linearSearchTest(int k, double r, int max_tau, int max_distance, int n_tests, ofstream& out_file)
{
    int l = 1 << k; // l = 2^k
    // erroneous binary search oracle (RBIS oracle). 
    ErroneousBinarySearch* ebs = new ErroneousBinarySearch(k, l);


    // setting some variables which are a function of k and r (see the paper).
    double rho = (r - 1) / 2;
    double zeta1 = (rho - sqrt(rho * rho - 4 * rho))/2;
    double zeta2 = (rho + sqrt(rho * rho - 4 * rho))/2;

    // cout << "rho, zeta1, and zeta2: " << rho << '\t' << zeta1 << '\t' << zeta2 << endl;

    // setting array costs for different values of tau
    double* avg_ratio = new double[max_tau + 1];
    double* max_ratio = new double[max_tau + 1];
    double* min_ratio = new double[max_tau + 1];

    out_file << "distance,\tbaseline1,\tbaseline2,\t";
    for (int tau = 0; tau <= max_tau; tau++)
        out_file << "avg. tau =" << tau << ",\t";
    for (int tau = 0; tau <= max_tau; tau++)
        out_file << "min. tau =" << tau << ",\t";
    for (int tau = 0; tau <= max_tau; tau++)
        out_file << "max. tau =" << tau << ",\t";
    out_file << endl;

    for (int dist = 1; dist <= max_distance; dist = ceil(dist * 1.01))
    {
        show_progress(log10((double)(dist)) / log10(max_distance));
        // reporting distance from the origin (x-axis)
        out_file << dist << ",\t";
        // reporting baseline costs
        double baseline1 = get_gemetric_alg_cost(dist, 1, zeta1) / dist;
        double baseline2 = get_gemetric_alg_cost(dist, 1, zeta2) / dist;
        out_file << baseline1 << ",\t" << baseline2 << ",\t";

        // first, find the index of the best algorithm.
        double b = get_b(k, rho, zeta2);
        int i_star = find_best_alg_from_family(dist, k, b);

        for (int tau = 0; tau <= max_tau; tau++)
        {
            int min_test_cost = INT_MAX;
            int max_test_cost = -1;

            int sum = 0;
            for (int test = 0; test < n_tests; test++)
            {
                bool* answers = create_response_vector(k, tau);

                int sel_alg = ebs->bs_limited_query_rightmost_child(2 * i_star + 1, answers, tau, false);
                sel_alg = (sel_alg - 1) / 2;

                int sel_alg_cost = get_gemetric_alg_cost(dist, pow(b, sel_alg), pow(b, l));
                sum += sel_alg_cost;
                if (sel_alg_cost < min_test_cost)
                    min_test_cost = sel_alg_cost;
                if (sel_alg_cost > max_test_cost)
                    max_test_cost = sel_alg_cost;
            }
            avg_ratio[tau] = (double)(sum) / (dist * n_tests);
            min_ratio[tau] = (double)(min_test_cost)/dist;
            max_ratio[tau] = (double)(max_test_cost)/dist;
        }

        for (int tau = 0; tau <= max_tau; tau++)
            out_file << avg_ratio[tau] << ",\t";
        for (int tau = 0; tau <= max_tau; tau++)
            out_file << min_ratio[tau] << ",\t";
        for (int tau = 0; tau <= max_tau; tau++)
            out_file << max_ratio[tau] << ",\t";
        out_file << endl;
    }
    delete[] avg_ratio;
    delete[] max_ratio;
    delete[] min_ratio;
    delete ebs;
}

/* reading the following values from an input file :
k(no.queries), r(robustness), max_tau(range of tau), 
max_distance(range of target), and n_tests(no.random tests for a given target)*/
bool read_input_file(string input_file_name, int& k, double& r, int& max_tau, int& max_distance, int& n_tests)
{
    string line;
    ifstream myfile(input_file_name);
    if (myfile.is_open())
    {
        getline(myfile, line);
        k = stoi(line);
        getline(myfile, line);
        r = stod(line);
        getline(myfile, line);
        max_tau = stoi(line);
        getline(myfile, line);
        max_distance = stoi(line);
        getline(myfile, line);
        n_tests = stoi(line);
        myfile.close();
        return true;
    }
    else cout << "Unable to open file " << input_file_name << endl;
    return false;
}

// the main method; reading the parameters from a file and passing them to the linearSearchTest(..) method
int main()
{
    cout << "Experiments with Linear Search with Noisy Queries !\n\n";
    cout << "INPUT:\n";
    cout << "k: the number of queries (e.g., k = 16)\n";
    cout << "r: the desired robustness (a real value, e.g., r = 15.00)\n";
    cout << "tau: the upper bound for the number of incorrect answers\n";
    cout << "max_tau: maximum value of tau; tau in [0..max-tau] will be checked (max_tau <= k/4, e.g., tau = 4)\n";
    cout << "max_distance: distance range from the origin (e.g., max-distance = 100000)\n";
    cout << "n_test: number of tests when taking average (e.g., n_test = 100)\n\n";
    cout << "OUTPUT:\n";
    cout << "The experimental performance ratio (competitive ratio) when running Robust Binary Interval Search (RBIS) using k queries on r-robust algorithms when tau queries are incorrectly answered, for tau in [0..maxtau]. For each value of tau, n-tests are made and the average/max/min competitive ratio is reported over n-tests.\n\n\n";
   
    cout << "please enter a file name with the values of k, r, max_tau, max_distance, n_test written in the same order, one in each line (e.g., input1.txt)\n";
    string input_file_name;
    cin >> input_file_name;

    int k, max_tau, max_distance, n_tests;
    double r;

    bool success = read_input_file(input_file_name, k, r, max_tau, max_distance, n_tests);
    if (success)
    {
        cout << "\nreading from " << input_file_name << ":\n";
        cout << "k is : " << k << "\nr is : " << r << "\nmax_tau is : " << max_tau << "\nmax_distance is : " << max_distance << "\nn_tests is : " << n_tests << "\n\n";
        srand(time(0));  // Initialize random number generator.

        string output_file_name;
        cout << "please enter the output file name (e.g., output1.csv)\n";
        cin >> output_file_name;
        cout << "\nwriting on " << output_file_name << ". please wait ...\n";

        ofstream out_file;
        out_file.open (output_file_name);

        chrono::steady_clock::time_point begin = std::chrono::steady_clock::now(); // to measure elapsd time
        linearSearchTest(k, r, max_tau, max_distance, n_tests, out_file);
        chrono::steady_clock::time_point end = std::chrono::steady_clock::now();

        cout << "\n\nthe results are written on " << output_file_name << ".\n" << endl;
        std::cout << "Time elapsed = " << chrono::duration_cast<std::chrono::microseconds>(end - begin).count() << "[micro-s]" << std::endl;
        out_file.close();
    }
} 

