#include <algorithm>
#include <cstring>
#include <iostream>
#include <queue>
#include <thread> 
#include <vector>
#include "genetic_algorithm.hpp"
#include "utimer.hpp"
#include <sched.h>
#include "json.hpp"
#include <fstream>

using namespace std;
using nlohmann::json;

#define CANVAS_WIDTH 800        // each x-coordinate of a city is comprised in [0,CANVAS_WIDTH]
#define CANVAS_HEIGHT 600       // each y-coordinate of a city is comprised in [0,CANVAS_HEIGHT]
#define CROSSOVER_PROB 0.8      // probability that cycle crossover is applied
#define MUTATION_PROB 0.1       // probability that mutation occurs

// parallel worker's task type
#define EOS -1
#define WAIT 0 
#define GENERATE_CHILDREN 1


int main(int argc, char *argv[]){

    // verify that the number of arguments from command line is correct
    if(argc==2 && (strcmp(argv[1],"-help")==0 || strcmp(argv[1],"-h")==0)){
        cout << "Usage is: " << argv[0] << " cities pop_size generations seed n_w" << endl;
        return 0;
    }
    if(argc!=6){
        cout << "The program expects 5 arguments from command line." << endl 
             << "Usage is: " << argv[0] << " cities pop_size generations seed n_w" << endl;
        return 0;
    }       

    // read the input from command line
    int n = atoi(argv[1]);               // number of cities
    int m = atoi(argv[2]);               // population size
    int generations = atoi(argv[3]);     // number of iterations/generations
    unsigned int seed = atoi(argv[4]);   // seed used for random initialization
    int nw = atoi(argv[5]);  

    if(m%2==1) m++;   // population size should be even
    // given that the crossover generates 2 individuals, the number of workers can be at most m/2 
    nw = min(nw, m/2);

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


        // the master assigns a task to worker i on tasks[i] and waits over ACK[i] 
        // worker i waits on tasks[i] and adds an acknowledgement on ACK[i]
        vector<int> tasks(nw); 
        vector<local_stats*> results(nw, NULL);
            
        // load balancing: each worker is assigned a static interval 
        vector<pair<int,int>> lb(nw);
        int interval_size = m/nw - (m/nw%2==1?1:0);
        int surplus = (m - nw*interval_size)/2;
        int start = 0;
        for(int i=0; i<nw; i++){
            int end = start + interval_size + (i<surplus?2:0);
            lb[i] = {start,end};
            start = end; 
        }
        

        // depending on its index, a worker generates children in population_next[start,end]
        auto generate_children = [&](int index, unsigned int &seed, vector<int> &crossover_map, vector<int> &crossover_cycle_number, local_stats &result){

            int start = lb[index].first;
            int end = lb[index].second;

            // use the crossover operator to generate children and then apply mutation
            for(int i=start; i<end; i+=2){
                
                int p1 = weighted_selection(population, seed, fitness_sum);
                int p2 = weighted_selection(population, seed, fitness_sum);
                crossover(population[p1], population[p2], population_next[i], population_next[i+1], CROSSOVER_PROB, seed, crossover_map, crossover_cycle_number);

                // apply mutation
                mutate(population_next[i], MUTATION_PROB, seed);
                mutate(population_next[i+1], MUTATION_PROB, seed);
            }

            // compute the local stats
            result.fitness_sum = 0;
            result.path_distance_sum = 0;
            result.best_path_index = start;
            result.best_path_distance = __LONG_LONG_MAX__;

            for(int i=start; i<end; i++){
                population_next[i].path_distance = path_distance(population_next[i], dist);
                population_next[i].fitness = 1.0/(1.0+population_next[i].path_distance);
                result.fitness_sum += population_next[i].fitness;
                result.path_distance_sum += population_next[i].path_distance;
                if(population_next[i].path_distance < result.best_path_distance){
                    result.best_path_index = i;
                    result.best_path_distance = population_next[i].path_distance;
                }
            }
            
            return;
        };

        // a generic worker loop
        auto worker_body = [&](int id){
            
            unsigned int worker_seed = rand_r(&seed);
            // allocating the vectors needed for crossover beforehand
            vector<int> map_p1(n), cycle_number(n);
            // local stats for each task
            local_stats result;

            while(true){

                while(tasks[id] == WAIT) this_thread::yield(); // busy wait

                if(tasks[id] == EOS) break;
                else if(tasks[id] == GENERATE_CHILDREN) generate_children(id, worker_seed, map_p1, cycle_number, result);

                tasks[id] = WAIT;
                results[id] = &result; // notify master
            }
        };

        // initialize workers
        vector<thread> tids(nw);
        for(int i=0; i<nw; i++){
            tids[i] = thread(worker_body, i);
        }

        // creating the indicated number of generations
        for(int curr_gen = 0; curr_gen<generations; curr_gen++){

            // program stats
            int best_path_index = 0;
            statistics[curr_gen].gen_number = curr_gen;
            statistics[curr_gen].avg_path_distance = 0;
            statistics[curr_gen].best_path_distance = __LONG_LONG_MAX__;
            double acc_fitness_sum = 0;

            // setting up tasks to generate children
            for(int i=0; i<nw; i++) tasks[i] = GENERATE_CHILDREN;

            // waiting for all workers to end
            for(int i=0; i<nw; i++){
                
                while(results[i]==NULL) this_thread::yield();

                // reducing the results
                acc_fitness_sum += results[i]->fitness_sum;
                statistics[curr_gen].avg_path_distance += results[i]->path_distance_sum;
                if(results[i]->best_path_distance < statistics[curr_gen].best_path_distance){
                    statistics[curr_gen].best_path_distance = results[i]->best_path_distance;
                    best_path_index = results[i]->best_path_index;
                }
                results[i] = NULL;
            }
            
            // replace old population with new one
            swap(population, population_next);

            // program stats
            fitness_sum = acc_fitness_sum;
            statistics[curr_gen].avg_path_distance /= population.size();
            for(int i=0; i<n; i++) statistics[curr_gen].best_path[i] = population[best_path_index].permutation[i];
        }

        // stopping workers
        for(int i=0; i<nw; i++) tasks[i] = EOS;
        for(int i=0; i<nw; i++) tids[i].join();
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
