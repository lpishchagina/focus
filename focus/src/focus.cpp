
#include <Rcpp.h>
#include <iostream>
#include <vector>
#include <libqhullcpp/Qhull.h>
#include <libqhullcpp/QhullFacetList.h>
#include <libqhullcpp/QhullVertexSet.h>
#include <libqhullcpp/PointCoordinates.h>
using namespace Rcpp;
using namespace orgQhull;

/*Function returns integer number:
 1  - gaussian cost, FOCuS0;
 11 - gaussian cost, FOCuS;
 2  - poisson cost, FOCuS0;
 22 - poisson cost, FOCuS;
 */
unsigned int get_costType(std::string method, 
                          std::string cost) {
  unsigned int costType = 0;
  if (cost == "gauss") costType = 1;
  if (cost == "poisson") costType = 2;
  if (method == "FOCuS") costType = costType * 11;
  return costType;
}

/* Function returns the value of cost function for FOCuS0 or FOCuS */
double get_cost(unsigned int typeCost, 
                unsigned int p, 
                unsigned int change, 
                unsigned int t, 
                double** &cumsumMatrix) {
  //Focus0 gauss
  double cost =  0;
  //memory
  double* right_cumsum = new double[p+1];
  //[t-change, S(t)-S(change)]
  for (size_t k = 0; k <= p; k++)
    right_cumsum[k] = cumsumMatrix[t][k] - cumsumMatrix[change][k];
  //Gaussian model
  if (typeCost == 1 || typeCost == 11) {
    //FOCuS0
    for (size_t k = 0; k < p; k++)
      cost = cost - right_cumsum[k+1] * right_cumsum[k+1] / right_cumsum[0];
    //FOCuS
    if (typeCost == 11)
      for (size_t k = 0; k < p; k++)
        cost = cost - cumsumMatrix[change][k+1] * cumsumMatrix[change][k+1] / (cumsumMatrix[change][0]+1);
  }
  //Poisson model
  if (typeCost == 2 || typeCost == 22) {
    //FOCuS0
    for (size_t k = 0; k < p; k++)
      if(right_cumsum[k+1] != 0)
        cost = cost + right_cumsum[k+1] * ( 1 - log(right_cumsum[k+1] / right_cumsum[0]));
    //FOCuS
    if (typeCost == 22)
      for (size_t k = 0; k < p; k++)
        if(cumsumMatrix[change][k+1] != 0)
          cost = cost + cumsumMatrix[change][k+1] * (1- log(cumsumMatrix[change][k+1] / (cumsumMatrix[change][0]+1)));
    cost = 2*cost; 
  }
  //clean memory
  delete []right_cumsum;
  return cost;
}

/*Function calculates matrix cumsumMatrix. (n+1)x(p+1)
 * row i :(i, cumsum x^k), k=1,..,p
 * */
void get_cumsumMatrix(unsigned int typeCost, 
               unsigned int p, 
               unsigned int length, 
               double** &cumsumMatrix,  
               Rcpp::NumericMatrix data) {// 0,0...0//1,x11,...,x1p//...
    for (size_t i = 0; i < length; i++) cumsumMatrix[i][0] = i; //time t
    for (size_t k = 0; k < p; k++) cumsumMatrix[0][k+1] = data(0, k);//row 0
    //cumsum
    for (size_t i = 1; i < length; i++)
      for (size_t k = 0; k < p; k++) 
        cumsumMatrix[i][k+1] = cumsumMatrix[i-1][k+1] + data(i, k);
}

/* Function obtains candidate's indices, 
 * builds the convex hull 
 * and removes candidate's indices that are not vertices of convex hull*/
void doPruning(unsigned int p, 
               std::list<int> &candidates, 
               bool* & isPruning, 
               double** & cumsumMatrix) {
  std::vector<coordT> points;       // coordinate vector of all candidates
  for ( auto it = candidates.begin(); it != candidates.end(); it++) {
    isPruning[(*it)] = true;        // flag: do pruning for all candidates at time t
    //we generate coordinate vector of all candidates from cumsumMatrix for Qhull
    for (unsigned int k = 0; k <= p; k++)
      points.push_back(cumsumMatrix[(*it)][k]);
  }
  //do Qhull
  orgQhull::Qhull qhull;
  qhull.runQhull("", p+1, candidates.size(), points.data(), "s");
  //get all vertices
  const orgQhull::QhullVertexList& vertices = qhull.vertexList();
  //Attention :  vertex.point().coordinates()[0] is always candidate's index 
  int index;
  for (const orgQhull::QhullVertex& vertex : vertices) {
    const double* coords = vertex.point().coordinates();
    index = int(vertex.point().coordinates()[0]);
    if (isPruning[index] == true) isPruning[index] = false;// change pruning's flag : this candidate exists at time t
  }
  //do pruning
  auto pruning_it = candidates.begin();
  while (pruning_it != candidates.end()) {
    if ( isPruning[(*pruning_it)] == true) {
      pruning_it = candidates.erase(pruning_it);
      --pruning_it;
    }
    ++pruning_it;
  }
}



// [[Rcpp::export]]
List getChangePoints(Rcpp::NumericMatrix data, 
                     std::string method = "FOCuS0", 
                     std::string cost = "gauss",
                     int common_difference_step = 1,
                     int common_ratio_step = 1, 
                     int first_step_qhull = 5,
                     bool cand_nb = false, 
                     bool opt_changes = false, 
                     bool opt_costs = false,
                     bool cands = false) {
  //stop
  if (method != "FOCuS0" && method !="FOCuS") {throw std::range_error("Parameter method should be FOCuS0 or FOCuS !");}
  if (cost != "gauss" && cost != "poisson") {throw std::range_error("Parameter cost should be gauss or poisson !");}
  
  //data parameters
  unsigned int p = (unsigned int)data.ncol();
  unsigned int length = (unsigned int) data.nrow();
  // cost and method
  unsigned int typeCost = get_costType(method, cost);
  
  //memory help tables and vectors
  bool* isPruning = new bool[length];   //for pruning's flag
  double** cumsumMatrix = new double*[length];
  for (unsigned int i = 0; i < length; i++) cumsumMatrix[i] = new double[p + 1];
  
  //pre-processing
  for (int i  = 0;  i < length; i++) isPruning[i] = false;
  get_cumsumMatrix(typeCost, p, length, cumsumMatrix, data);
  
  //get first step for qhull
  int next_step_qhull = first_step_qhull + p;         // step > p
  if (cost == "poisson")
    next_step_qhull = next_step_qhull  + 10; //ATTENTION: do the constant for poisson model, now it is 10!!! 
  
  //initialization
  std::list<int> candidates;     // change point candidates at time t
  std::vector<unsigned int> nb_at_time;   // vector of candidate number at time t
  std::vector<double>  opt_cost;          // optimal cost at time t
  std::vector<int>  opt_change;            // optimal change at time t
  std::vector<unsigned int> candidate_set; // candidates at time t=length-1
  for (unsigned int i = 0; i < length; i++) {
   nb_at_time.push_back(0);
    opt_cost.push_back(INFINITY);
    opt_change.push_back(0);
  }
  opt_cost[0] = 0;
  double candidate_cost;
  //loop
  for (unsigned int i = 1; i < length; i++) {
    candidates.push_back(i-1);
    nb_at_time[i] = candidates.size();
    // optimal cost and change
     for(auto it = candidates.begin(); it != candidates.end(); it++){
       candidate_cost = get_cost(typeCost, p, (*it), i, cumsumMatrix);
       
       if (candidate_cost < opt_cost[i]) {                      //find min
         opt_cost[i] = candidate_cost;
         opt_change[i] = (*it) + 1; //for R
       }
     }
     //pruning
     if(nb_at_time[i] >= next_step_qhull ||  i == (length - 1)){
         doPruning(p, candidates, isPruning, cumsumMatrix);
        // update-next_step_qhull
         next_step_qhull = common_ratio_step * candidates.size() + common_difference_step;
     }
  }
  //get res
  unsigned int change  = opt_change[length-1];
  std::vector<double> pre_change_mean;
  std::vector<double> post_change_mean;
  for(unsigned int k = 1; k <= p; k++) {
    pre_change_mean.push_back(cumsumMatrix[change-1][k]/change); 
    post_change_mean.push_back((cumsumMatrix[length-1][k] - cumsumMatrix[change-1][k])/(length - change)); 
  }
  //add the last nb_cand
  nb_at_time.push_back(candidates.size());
  //memory
  for (unsigned int i = 0; i < length; i++) delete(cumsumMatrix[i]);
  delete [] isPruning;
  delete [] cumsumMatrix ;
  isPruning = NULL;
  cumsumMatrix = NULL;
  //res
  List res;
  res["change"] = change;
  res["cost"] = opt_cost[length-1];
  res["pre_change_mean"] = pre_change_mean;
  res["post_change_mean"] = post_change_mean;
  
  res["nb_candidates"] = nb_at_time[length];
  
  if (cands){
    std::vector<int> cands_n;
    auto pruning_it = candidates.begin();
    while (pruning_it != candidates.end()) {
      cands_n.push_back((*pruning_it)+1);//for R
      ++pruning_it;
    }
    res["candidates"] = cands_n;
  } 
  if (cand_nb) res["nb_candidates_at_time"] = nb_at_time;
  if (opt_changes) res["opt_changes_at_time"] = opt_change;
  if (opt_costs) res["opt_costs_at_time"] = opt_cost;
  return res;
}
    