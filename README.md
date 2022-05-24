# LinearSearchWithNoisyQueries
Implementation of robust algorithms that use noisy queries to search for a hidden target on an infinite line.

# Input
the following inputs should be written on a file and be given as input to the program.
k: the number of queries (e.g., k = 16)
r: the desired robustness (a real value, e.g., r = 15.00)
tau: the upper bound for the number of incorrect answers
max_tau: maximum value of tau; tau in [0..max-tau] will be checked (max_tau <= k/4, e.g., tau = 4)
max_distance: distance range from the origin (e.g., max-distance = 100000)
n_test: number of tests when taking average (e.g., n_test = 100)

# Output  
The experimental performance ratio (competitive ratio) when running Robust Binary Interval Search (RBIS) algorithm using k queries on r-robust algorithms when tau queries are incorrectly answered, for tau in [0..maxtau]. We use a variation of RBIS that returns the rightmost descendant of the node at which the search stops.
For each value of tau, n-tests are made and the average/max/min competitive ratio is reported (on a CSV file) over n-tests.

# Prerequisites
a c++ compiler, e.g., g++ 8.3.0 or higher

# Compiling the Project
type the followings: 
g++ HeapBst.h g++ ErroneousBinarySearch.h g++ LinearSearch.cpp -o LinSearch

# Running the Program 
type the following: LinSearch.exe

Follow the instructions to enter an input file name (e.g., "input.txt"). 
The input file should be located in the same folder.



