#include <algorithm>
#include <cstring>
#include <iostream>
#include <vector>
#include "genetic_algorithm.hpp"
#include "utimer.hpp"
#include "json.hpp"
#include <fstream>

using namespace std;

#define CANVAS_WIDTH 800        // each x-coordinate of a city is comprised in [0,CANVAS_WIDTH]
#define CANVAS_HEIGHT 600       // each y-coordinate of a city is comprised in [0,CANVAS_HEIGHT]
#define CROSSOVER_PROB 0.8      // probability that cycle crossover is applied
#define MUTATION_PROB 0.1       // probability that mutation occurs


int main(int argc, char *argv[]){

    // verify that the number of arguments from command line is correct
    if(argc==2 && (strcmp(argv[1],"-help")==0 || strcmp(argv[1],"-h")==0)){
        cout << "Usage is: " << argv[0] << " cities pop_size generations seed" << endl;
        return 0;
    }
    if(argc!=5){
        cout << "The program expects 4 arguments from command line." << endl  
             << "Usage is: " << argv[0] << " cities pop_size generations seed" << endl;
        return 0;
    }

    // read the input from command line
    int n = atoi(argv[1]);               // number of cities
    int m = atoi(argv[2]);               // population size
    int generations = atoi(argv[3]);     // number of iterations/generations
    unsigned int seed = atoi(argv[4]);   // seed used for random initialization

    if(m%2==1) m++;   // population size should be even

    // initialize cities positions at random with x in [0, CANVAS_WIDTH] and y in [0, CANVAS_HEIGHT]
    vector<pair<double,double>> points(n);
    for(int i=0; i<n; i++){
        points[i] = { random_double(seed, 0, CANVAS_WIDTH), random_double(seed, 0, CANVAS_HEIGHT) };
    }
    // precompute cities distances
    vector<vector<double>> dist(n, vector<double>(n));
    for(int i=0; i<n; i++){ 
        for(int j=0; j<n; j++){
            if(i==j) dist[i][j] = 0;
            else if(i<j){
                dist[i][j] = points_distance(points[i], points[j]);
            }
            else dist[i][j] = dist[j][i]; // the graph is undirected, so dist[i][j] = dist[j][i]
        }
    }

    // allocate memory for the current population and a generic next population
    vector<individual> population(m), population_next(m);
    for(int i=0; i<population.size(); i++){
        population[i].permutation.resize(n);
        population_next[i].permutation.resize(n);
    }

    // allocating the vectors needed for crossover
    vector<int> map_p1(n), cycle_number(n);

    // allocate memory for program statistics
    vector<gen_stats> statistics(generations);
    for(int i=0; i<statistics.size(); i++) statistics[i].best_path.resize(n);

    long completion_time;
    {
        utimer t("", &completion_time);

        // STEP 1: create an initial random population P
        for(int i=0; i<population.size(); i++){
            for(int j=0; j<n; j++) population[i].permutation.at(j)=j;
            random_shuffle(population[i].permutation.begin(), population[i].permutation.end());
        }

        // STEP 2: compute the fitness of each individual and fitness_sum
        double fitness_sum = 0;
        for(int i=0; i<population.size(); i++){
            population[i].path_distance = path_distance(population[i], dist);
            population[i].fitness = 1.0/(1.0+population[i].path_distance);
            fitness_sum += population[i].fitness;
        }

        // creating the indicated number of generations
        for(int curr_gen = 0; curr_gen<generations; curr_gen++){

            // program stats
            int best_path_index = 0;
            statistics[curr_gen].gen_number = curr_gen;
            statistics[curr_gen].avg_path_distance = 0;
            statistics[curr_gen].best_path_distance = __LONG_LONG_MAX__;

            // STEP 3: use crossover to create a new generation
            for(int i=0; i<population.size(); i+=2){
                int p1 = weighted_selection(population, seed, fitness_sum);
                int p2 = weighted_selection(population, seed, fitness_sum);
                crossover(population[p1], population[p2], population_next[i], population_next[i+1], CROSSOVER_PROB, seed, map_p1, cycle_number);
            }

            // STEP 4: apply mutation
            for(int i=0; i<population_next.size(); i++){
                mutate(population_next[i], MUTATION_PROB, seed);
            }

            // STEP 5: replace old population with new one
            swap(population, population_next);

            // STEP 6: compute the fitness of each individual and fitness_sum
            fitness_sum = 0;
            for(int i=0; i<population.size(); i++){
                population[i].path_distance = path_distance(population[i], dist);
                population[i].fitness = 1.0/(1.0+population[i].path_distance);
                fitness_sum += population[i].fitness; 

                // program stats
                statistics[curr_gen].avg_path_distance += population[i].path_distance;
                if(population[i].path_distance < statistics[curr_gen].best_path_distance){
                    statistics[curr_gen].best_path_distance = population[i].path_distance;
                    best_path_index = i;
                }
            }

            // program stats
            statistics[curr_gen].avg_path_distance /= population.size();
            for(int i=0; i<n; i++) statistics[curr_gen].best_path[i] = population[best_path_index].permutation[i];
        }
    }

    // save program stats in JSON format
    nlohmann::json jsonfile;
    jsonfile["width"] = CANVAS_WIDTH;
    jsonfile["height"] = CANVAS_HEIGHT;
    jsonfile["n_cities"] = n;
    jsonfile["points"] = points;
    jsonfile["population_size"] = m;
    jsonfile["n_generations"] = generations;

    vector<nlohmann::json_abi_v3_11_2::json> gen_stats;
    for(auto stat : statistics){
        gen_stats.push_back(
            {
                { "gen_number", stat.gen_number },
                { "avg_path_distance", stat.avg_path_distance },
                { "best_path", stat.best_path },
                { "best_path_distance", stat.best_path_distance }
            }
        );
    }
    jsonfile["generations"] = gen_stats;
    std::ofstream file("../results.json");
    file << jsonfile;

    cout << "Completion Time: " << completion_time << " usec" << endl;

    return 0;
}
