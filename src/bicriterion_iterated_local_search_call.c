#include <stdlib.h> 
#include <stdio.h>
#include <math.h>
#include <stdbool.h>
#include <string.h>

//datastructure to store a linked-list of partitions
struct Pareto_element { 
  double diversity; 
  double dispersion;
  int *partition;
  struct Pareto_element *next;
};


void printList(size_t N, struct Pareto_element* n) { 
  
  int pos = 1;
  while (n != NULL) {
    int *array = n->partition;
    printf("%d. partition: ", pos);
    for(int i = 0; i < N; i++){
      printf("%d" , array[i]);
    }
    printf("  dis: %.17le  div:%.17le \n", n->dispersion, n->diversity);
    n = n->next;
    pos++;
  } 
  printf("\n \n");
}   


// Declare Functions
void bicriterion_iterated_local_search_call(double *distances, int *N, int *G, int *R,int *upper_bound, int *WL, double *W, double *Xi, double *result);
struct Pareto_element* multistart_bicriterion_pairwise_interchange(size_t N, double matrix[N][N], int G, int R, int WL, double weights[WL]);
struct Pareto_element* bicriterion_iterated_local_search(struct Pareto_element* head, size_t N, double matrix[N][N],int G, int R, int WL, double weights[WL], double neighbor_percent[2]);
double sample(size_t array_size,double array[array_size]);
void random_partition(size_t N, int G, int random[N]);
double get_diversity(size_t N, int partition[N], double matrix[N][N]);
double get_dispersion(size_t N, int partition[N], double matrix[N][N]);
void cluster_swap(size_t N, int i, int j, int partition[N]);
void update_pareto(struct Pareto_element** head_ref, size_t N, int partition[N], double diversity, double dispersion);
void compress(size_t N, int partition[N]);
bool paretodominated(struct Pareto_element* head, double diversity, double dispersion);
bool paretoincluded(struct Pareto_element* head, size_t N, int partition[N]);
void push(struct Pareto_element** head_ref, double diversity, double dispersion, size_t N, int partition[N]);
void delete_outdated(struct Pareto_element** head_ref, double diversity, double dispersion);
void linked_list_sample(struct Pareto_element* head, size_t N, int* partition);
int linked_list_length(struct Pareto_element* head);
double random_in_range(double min, double max);
double get_diversity_fast(double diversity, int x,int y, size_t N, int partition[N], double matrix[N][N]);
double get_dispersion_fast(double dispersion, int x,int y, size_t N, int partition[N], double matrix[N][N], int *counter1, int *counter2, int *counter3,int *counter4, int *counter);


//receive data from r, call bils algorithm, save results for r
void bicriterion_iterated_local_search_call(double *distances, int *N, int *G, int *R, int *upper_bound, int *WL, double *W, double *Xi, double *result) {
  
  const int n = *N; // number of elements
  const int g = *G; // number of clusters
  const int r = *R; // number of restarts
  const int u = *upper_bound; //max. length of result-list
  const int wl = *WL; // length of possible weights
  
  int distance_ptr[n]; // indexing columns
  double distance_pts[n][n];
  
  // Column offsets (to convert one-dimensional array to Row/Col major)
  for(int i = 0; i < n; i++) {
     distance_ptr[i] = i * n;
  }
  
  // Reconstruct the data points as N x N matrix
  for(int i = 0; i < n; i++) {
    for(int j = 0; j < n; j++) {
      distance_pts[i][j] = distances[distance_ptr[j]++];
    }
  }
  
  double weights[wl];
  for(int i = 0; i < wl; i++){
    weights[i] = W[i];
  }
  
  double neighbor_percent[2];
  neighbor_percent[0] = Xi[0];
  neighbor_percent[1] = Xi[1];
  
  //divide restarts for both parts of the algorithm equally
  int half_restarts = r/2 + (r%2);
  
  struct Pareto_element* head = multistart_bicriterion_pairwise_interchange(n, distance_pts, g, half_restarts, wl, weights);
  head = bicriterion_iterated_local_search(head, n, distance_pts, g, half_restarts, wl, weights, neighbor_percent);
  
  //fill result with -1 to identify the last partition in the matrix
  for(int i = 0; i < n*u; i++){
    result[i] = -1;
  }
  
  //update result with the solutions from the linked list
  int position = 0;
  while (head != NULL && position != n*u) {
    int *array = head->partition;
    for(int i = 0; i < n; i++){
      result[position] = array[i];
      position++;
    }
    head = head->next;
  } 
}
  

struct Pareto_element* multistart_bicriterion_pairwise_interchange(size_t N, double matrix[N][N], int G, int R, int WL, double weights[WL]) {
  
  int counter = 0;
  int counter1 = 0;
  int counter2 = 0;
  int counter3 = 0;
  int counter4 = 0;
  
  struct Pareto_element* head = NULL; // head pointing on the later linked list(paretoset)
  
  for(int a = 0; a < R; a++){
    double div_weight = sample(WL, weights); 
    double dis_weight = 1 - div_weight;
    int partition[N];
    random_partition(N,G,partition);
    double diversity = get_diversity(N, partition, matrix);
    double save_diversity = diversity;
    double dispersion = get_dispersion(N, partition, matrix);
    double save_dispersion = dispersion;
    double max_bicriterion = div_weight*diversity + dis_weight*dispersion;
    bool Flag = false;
    while(!Flag){
      Flag = true;
      for(int i = 0; i < N-1; i++){
        for(int j = i+1; j < N; j++){
          int g = partition[i];
          int h = partition[j];
          if(g != h){
            cluster_swap(N, i, j, partition);
            double current_diversity = get_diversity_fast(save_diversity, i, j, N, partition, matrix);
            //double current_diversity = get_diversity(N, partition, matrix);
            double current_dispersion = get_dispersion_fast(save_dispersion , i, j, N, partition, matrix, &counter1,&counter2,&counter3,&counter4,&counter);
            //double current_dispersion = get_dispersion(N, partition, matrix);
            update_pareto(&head, N, partition,current_diversity, current_dispersion);
            double current_bicriterion = div_weight*current_diversity + dis_weight*current_dispersion;
            if(current_bicriterion > max_bicriterion){
              save_diversity = current_diversity;
              save_dispersion = current_dispersion;
              max_bicriterion = current_bicriterion;
              Flag = false;
            }else{
              cluster_swap(N, i, j, partition);
            }
          }
        }
      }
    }
  }
  printf("%d  %d  %d  %d  %d \n", counter1,counter2, counter3, counter4, counter);
  return (head);
}  


struct Pareto_element* bicriterion_iterated_local_search(struct Pareto_element* head, size_t N, double matrix[N][N],int G, int R, int WL, double weights[WL], double neighbor_percent[2]){
  
  int bcounter = 0;
  int bcounter1 = 0;
  int bcounter2 = 0;
  int bcounter3 = 0;
  int bcounter4 = 0;
  
  for(int a = 0; a < R; a++){
    double div_weight = sample(10, weights); 
    double dis_weight = 1 - div_weight;
    double neighborhood_size = random_in_range(neighbor_percent[0], neighbor_percent[1]);
    int partition[N];
    linked_list_sample(head, N, partition);
    for(int i = 0; i < N-1; i++){
      for(int j = i + 1; j < N; j++){
        int g = partition[i];
        int h = partition[j];
        if(g != h){
          double random = random_in_range(0,1);
          if(random < neighborhood_size){
            cluster_swap(N, i, j, partition);
          }
        }
      }
    }
    double diversity = get_diversity(N, partition, matrix);
    double save_diversity = diversity;
    double dispersion = get_dispersion(N, partition, matrix);
    double save_dispersion = dispersion;
    double max_bicriterion = div_weight*diversity + dis_weight*dispersion;
    bool Flag = false;
    while(!Flag){
      Flag = true;
      for(int i = 0; i < N-1; i++){
        for(int j = i+1; j < N; j++){
          int g = partition[i];
          int h = partition[j];
          if(g != h){
            cluster_swap(N, i, j, partition);
            double current_diversity = get_diversity_fast(save_diversity, i, j, N, partition, matrix);
            //double current_diversity = get_diversity(N, partition, matrix);
            double current_dispersion = get_dispersion_fast(save_dispersion , i, j, N, partition, matrix, &bcounter1,&bcounter2,&bcounter3,&bcounter4,&bcounter);
            //double current_dispersion = get_dispersion(N, partition, matrix);
            update_pareto(&head, N, partition, current_diversity, current_dispersion);
            double current_bicriterion = div_weight*current_diversity + dis_weight*current_dispersion;
            if(current_bicriterion > max_bicriterion){
              save_diversity = current_diversity;
              save_dispersion = current_dispersion;
              max_bicriterion = current_bicriterion;
              Flag = false;
            }else{
              cluster_swap(N, i, j, partition);
            }
          }
        }
      }
    }
  }
  printf("%d  %d  %d  %d  %d \n", bcounter1,bcounter2, bcounter3, bcounter4, bcounter);
  //printList(N, head);
  return (head);
}     

//sample functin from R. Get one random element from an array
double sample(size_t array_size,double array[array_size]){
  
  int r = rand() % array_size;
  return(array[r]);
}


void random_partition(size_t N, int G, int partition[N]){
  
  int clustersize = N/G;
  int freeclusters[G];
  
  for(int i = 0; i < G; i++){
    freeclusters[i] = clustersize; 
  }
  
  for(int i = 0; i < N; i++){
    int r = rand() % G;
    if(freeclusters[r] > 0){
      partition[i] = r;
      freeclusters[r]--;
    }else{
      i--;
    }
  }
}


double get_diversity(size_t N, int partition[N], double matrix[N][N]){
  
  double sum = 0;
  
  for(size_t i = 0; i < N-1; i++){
    for(size_t j = i+1; j < N; j++){
      if(partition[i] == partition[j]){
        sum = sum + matrix[i][j];
      }
    }
  }
  
  return(sum);  
}


double get_dispersion(size_t N, int partition[N], double matrix[N][N]){
  
  double min = INFINITY;
  double distance;
  
  for(size_t i = 0; i < N-1; i++){
    for(size_t j = i+1; j < N; j++){
      if(partition[i] == partition[j]){
        distance = matrix[i][j];
        if(distance < min){
          min = distance;
        }
      }
    }
  }
  
  return(min);
}


void cluster_swap(size_t N, int i, int j, int partition[N]){
  
  int save = partition[i];
  partition[i] = partition[j];
  partition[j] = save;
}


void update_pareto(struct Pareto_element** head_ref, size_t N, int partition[N], double diversity, double dispersion){
  
  if(*head_ref == NULL){
    compress(N, partition);
    push(head_ref, diversity, dispersion, N, partition);
    
  }else{
    if(!(paretodominated(*head_ref, diversity, dispersion))){
      compress(N, partition);
      if(!(paretoincluded(*head_ref, N, partition))){
        delete_outdated(head_ref, diversity, dispersion);
        push(head_ref, diversity, dispersion, N, partition);
      }
    }
  }
}


//compress partition to uniform ascending order
void compress(size_t N, int partition[N]){
  
  int cluster = 0;
  int change;
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster){
      cluster++;
    }else if(partition[i] > cluster){
      change = partition[i];
      for(int j = i; j < N; j++){
        if(change == partition[j]){
          partition[j] = cluster;
        }else if(cluster == partition[j]){
          partition[j] = change;          
        }
      }
      cluster++;
    }
  } 
}


bool paretodominated(struct Pareto_element* head, double diversity, double dispersion){
  
  while(head != NULL){
    if((head->diversity >= diversity && head->dispersion > dispersion) || (head->diversity > diversity && head->dispersion >= dispersion)){
      return (true);
    }
    head = head->next;
  }
  
  return (false);
}


bool paretoincluded(struct Pareto_element* head, size_t N, int partition[N]){
  
  int counter;
  
  while(head != NULL){
    int * array = head->partition;
    counter = 0;
    for(int i = 0; i < N; i++){
      if(partition[i] == array[i]){
        counter++;
      }
    }
    if(N == counter){
      return(true);
    }
    head = head->next;
  }
  
  return(false);
}


void push(struct Pareto_element** head_ref, double diversity, double dispersion, size_t N, int partition[N]) { 
  
  /* 1. allocate node */
    struct Pareto_element* new_node = (struct Pareto_element*) malloc(sizeof(struct Pareto_element)); 
    
    /* 2. put in the data  */
      new_node->diversity  = diversity;
      new_node->dispersion = dispersion;
      new_node->partition =(int *)malloc(sizeof(int)*N ) ;  // vllt memcpy anstatt schleife?
        for(int b = 0;b<N;b++){ 
          new_node->partition[b]  = partition[b];
        }
      
      /* 3. Make next of new node as head */
        new_node->next = (*head_ref); 
        
        /* 4. move the head to point to the new node */
          (*head_ref)    = new_node; 
}  


void delete_outdated(struct Pareto_element** head_ref, double diversity, double dispersion){
  
  struct Pareto_element* temp = *head_ref;
  struct Pareto_element* prev = *head_ref;
  struct Pareto_element* safe = NULL;
  
  while(temp != NULL && ((diversity >= temp->diversity && dispersion > temp->dispersion) || (diversity > temp->diversity && dispersion >= temp->dispersion))){
    *head_ref = temp->next; 
    prev = temp->next;
    free(temp); 
    temp = prev;
  }
  if(temp != NULL){
    temp = temp->next;
  }
  
  while(temp != NULL){
    if((diversity >= temp->diversity && dispersion > temp->dispersion) || (diversity > temp->diversity && dispersion >= temp->dispersion)){
      prev->next = temp->next;
      safe = temp->next;  
      free(temp);
      temp = safe;
    }else{
      prev = prev->next;
      temp = temp->next;
    }
  }
} 


void linked_list_sample(struct Pareto_element* head, size_t N, int* partition){
  
  int count = linked_list_length(head);
  int r = rand() % count;
  struct Pareto_element* current = head; 
  
  for(int i = 0; i < r; i++){
    current = current->next;
  }
  
  for(int i = 0; i < N; i++){
    partition[i] = current->partition[i];
  }
  
}


double random_in_range(double min, double max) {
  
  double range = (max - min); 
  double div = RAND_MAX / range;
  return min + (rand() / div);
}


int linked_list_length(struct Pareto_element* head) {
  
  int count = 0;   
  struct Pareto_element* current = head;  
  while (current != NULL) 
  { 
    count++; 
    current = current->next; 
  } 
  return count; 
} 


double get_diversity_fast(double diversity, int x,int y, size_t N, int partition[N], double matrix[N][N]){
  
  int cluster_x = partition[x];
  int cluster_y = partition[y];
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster_y && i != x && i != y){
      diversity -= matrix[i][x];
    }
  }
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster_x && i != x){
      diversity += matrix[i][x];
    }
  }
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster_x && i != x && i != y){
      diversity -= matrix[i][y];
    }
  }
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster_y && i != y){
      diversity += matrix[i][y];
    }
  }
  
  return(diversity);  
}


double get_dispersion_fast(double dispersion, int x,int y, size_t N, int partition[N], double matrix[N][N], int *counter1, int *counter2, int *counter3,int *counter4, int *counter){
  
  *counter = *counter + 1;

  int cluster_x = partition[x];
  int cluster_y = partition[y];
  bool before = false;
  bool after = false;
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster_y  && i != x && i != y){
      if(dispersion == matrix[i][x]){
        before = true;
        break;
      }
    }
  }
  
  for(int i = 0; i < N && (!before); i++){
    if(partition[i] == cluster_x && i != x && i != y){
      if(dispersion == matrix[i][y]){
        before = true;
        break;
      }
    }
  }
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster_x && i != x){
      if(dispersion >= matrix[i][x]){
        dispersion = matrix[i][x];
        after = true;
      }
    }
  }
  
  for(int i = 0; i < N; i++){
    if(partition[i] == cluster_y && i != y){
      if(dispersion >= matrix[i][y]){
        dispersion = matrix[i][y];
        after = true;
      }
    }
  }
  
  
  if(!before && !after){
    *counter1 = *counter1 + 1;
    return(dispersion);
  }
  
  if(!before && after){
    *counter2 = *counter2 + 1;
    return(dispersion);
  }
  
  if(before && after){
    *counter3 = *counter3 + 1;
    return(dispersion);
  }
  
  if(before && !after){
    *counter4 = *counter4 + 1;
    dispersion = get_dispersion(N, partition, matrix);
    return(dispersion);
  }
  
  return(dispersion);
}