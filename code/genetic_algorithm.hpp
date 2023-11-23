#include <algorithm>
#include <functional>
#include <vector>
#include <cmath>


struct individual{
    std::vector<int> permutation;
    double path_distance;
    double fitness;
};

// statistics for a single generation
struct gen_stats{
    int gen_number;
    double avg_path_distance;
    std::vector<int> best_path;
    double best_path_distance;
};

// statistics computed by a parallel worker over its assigned individuals.
// the master worker will then reduce the local stats to obtain the generation stats.
struct local_stats{
    double path_distance_sum;
    double best_path_distance;
    int best_path_index;
    double fitness_sum;
};

// generates a random floating-point number comprised in [lower_bound, upper_bound]
// pre-consition: 'lower_bound' <= 'upper_bound'
double random_double(unsigned int &seed, double lower_bound, double upper_bound){
    double range = (upper_bound - lower_bound); 
    double div = RAND_MAX / range;
    return lower_bound + (rand_r(&seed) / div);
}

// computes the distance between two points 'a' and 'b'
double points_distance(std::pair<double, double> a, std::pair<double, double> b){
    return sqrt(pow(b.first - a.first, 2) + pow(b.second - a.second, 2) * 1.0); 
}

// compute the path distance of an individual
double path_distance(individual &p, std::vector<std::vector<double>> &dist){
    double result = 0;
    for(int i=0; i<p.permutation.size(); i++){
        int city1 = p.permutation[i];
        int city2 = p.permutation[(i+1) % p.permutation.size()];
        result += dist[city1][city2];
    }
    return result;
}

// weighted selection of an individual from the population
int weighted_selection(std::vector<individual> &population, unsigned int &seed, double fitness_sum){

  int index = 0;
  // generate a random double in [0,1]
  double r = random_double(seed, 0, 1);
  // we want to find the first p_k such that:
  //        r <= \sum_{i=0}^k (fitness(p_i)/fitness_sum)
  // this is equivalent to 
  //        fitness_sum * r <= \sum_{i=0}^k fitness(p_i)
  r *= fitness_sum; 
  while(r>0){
      r -= population[index].fitness;
      // because of floating-point arithmetics, population probabilities might not add to 1
      // for this reason we also check that index is always less than population.size()
      if(r>0 && index+1<population.size()) index++;
  }
  return index;
}

void crossover(individual &p1, individual &p2, individual &c1, individual &c2, 
               double crossover_probability, unsigned int &seed, 
               std::vector<int> &map_p1, std::vector<int> &cycle_number){
    
    double r = random_double(seed, 0, 1);
    int n = p1.permutation.size();
    
    // CASE 1: children are exact copies of parents
    if(r > crossover_probability){
        for(int i=0; i<n; i++){
            c1.permutation[i] = p1.permutation[i];
            c2.permutation[i] = p2.permutation[i];
        }
        return;
    }

    // CASE 2: apply cycle crossover
    // preprocessing: map_p1[i]=k means that p1.permutation[k]=i 
    for(int i=0; i<n; i++){
        map_p1[p1.permutation[i]] = i;  
    }
    // divide the permutation into cycles
    int curr_index = 0, curr_cycle = 1;

    for(int i=0; i<n; i++) cycle_number[i]=0;

    while(curr_index < n){
        int initial = curr_index;
        cycle_number[curr_index] = curr_cycle;
        curr_index = map_p1[p2.permutation[curr_index]];
        // while the cycle is not closed      
        while(curr_index!=initial){
            cycle_number[curr_index] = curr_cycle;
            curr_index = map_p1[p2.permutation[curr_index]];
        }
        curr_cycle++;
        // find the next index that has no cycle number assigned
        while(curr_index<n && cycle_number[curr_index]!=0) curr_index++; 
    }

    // generate offspring
    for(int i=0; i<n; i++){
        if(cycle_number[i]%2==1){
            c1.permutation[i] = p1.permutation[i];
            c2.permutation[i] = p2.permutation[i];
        }
        else{
            c1.permutation[i] = p2.permutation[i];
            c2.permutation[i] = p1.permutation[i];
        }  
    }
}

void mutate(individual &p, double mutation_prob, unsigned int &seed){
    
    double r = random_double(seed, 0, 1);
    int n = p.permutation.size();

    if(r <= mutation_prob){
        int start = rand_r(&seed) % n;
        int end = rand_r(&seed) % n;
        if(start>end) std::swap(start,end);
        reverse(p.permutation.begin()+start, p.permutation.begin()+end+1);
    }
}
