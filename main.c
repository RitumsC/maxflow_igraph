/*
* Programm implementing MCMF agorithm for sorbtion
*
* TODO
*/
#include <stdio.h>
#include <igraph.h>
#include <time.h>
#include <string.h>

#define NEWLINE printf("\n");
#define EPS_RO 1e-6
#define EPS_MU 1e-5
#define MAX(x, y) (((x) > (y)) ? (x) : (y))

/* Capacity of source-sink edge  */
double c_st = 0;

/* Array of matrix-fluid energies for each cell*/
double *mf;

/* Sum capacity for each cell */
double *sum_cap;

int Network_size, coordination, N_inner_edges;


/* Function prototypes */
igraph_t make_graph(const char *fname, const double y, igraph_vector_t *weights);

void add_st( const double mu, igraph_vector_t *weights );

int get_ro( const igraph_t *graph, igraph_vector_t *weights, const double mu );

void dichotomy_right(const igraph_t *graph, igraph_vector_t *weights,
		     const double mu_i, const double mu_f,
		     const int ro_i, const int ro_f);

void dichotomy_left(const igraph_t *graph, igraph_vector_t *weights,
		    const double mu_i, const double mu_f,
		    const int ro_i, const int ro_f);

void dichotomy(const igraph_t *graph, igraph_vector_t *weights,
	       const double mu_i, const double mu_f);


/* Main start */
int
main(int argc, char **argv){

#if DEBUG
  printf("Writing in debug mode! \n");
#endif

  clock_t start_t, end_t;
  start_t = clock();

  /*  Argument list   */
  char *fname = argv[1];
  char *oname = argv[2];
  double y    = atof(argv[3]);
  double mu_i  = atof(argv[4]);
  double mu_f  = atof(argv[5]);

  /* Main code */
  igraph_vector_t weights;
  igraph_t graph = make_graph(fname, y, &weights);
  (void)dichotomy(&graph, &weights, mu_i, mu_f);

  /* Free up used space */
  free(mf);
  free(sum_cap);
  igraph_destroy(&graph);
  igraph_vector_destroy(&weights);

  end_t = clock();
  printf("Total CPU time: %fs\n", (double)(end_t-start_t)/ CLOCKS_PER_SEC);

  return 0;
}



/*************************************************
 * Create graph from fname given wettability y.  *
 * initializes only internal edges               *
 *************************************************/
igraph_t
make_graph(const char *fname, const double y, igraph_vector_t *weights){

  clock_t start_graph_t, end_graph_t;
  start_graph_t = clock();

  /* Read the file and obtain size and coordination information */
  FILE *read;
  if ( (read = fopen(fname, "r")) == NULL){
    printf("Error reading network file.");
    exit(1);
  }
  fscanf(read, "%d", &Network_size);
  fscanf(read, "%d", &coordination);
  
  /* Given known size allocate space for mf and sum_cap */
  mf = (double*) malloc(Network_size * sizeof(double));
  sum_cap = (double*) malloc(Network_size * sizeof(double));
  
  /* Start construction the graph */
  igraph_vector_init(weights, 0);
  igraph_t graph;
  igraph_empty(&graph, Network_size + 2, 1);
  
  int cell_ID, Nedges, nn;
  igraph_vector_t node_to_node;
  /* Edge list of type (from, to), ...  */
  /* initialize as max possible size to reduce time */
  /* note that the length is 1/2 of max edges as graph is directed */
  igraph_vector_init(&node_to_node, Network_size*coordination);
  igraph_vector_init(weights, 0.5*Network_size*coordination);
  int half_actual_size = 0;
  for(int i = 0; i < Network_size; ++i){
    fscanf(read, "%d", &cell_ID);
    fscanf(read, "%d", &Nedges);
    mf[i] = (coordination - Nedges)*y;
    double sum_cap_node = 0;
#if DEBUG
    printf("Added mf to node %d with value of%lf\n", i, mf[i]);
#endif
    for(int j = 0; j < Nedges; ++j){
      fscanf(read, "%d", &nn);
      if ( nn > i ){                      /* only add edge if nn > i */
	sum_cap_node += 1;
	VECTOR(node_to_node)[2*half_actual_size]=i+1;
	VECTOR(node_to_node)[2*half_actual_size+1]=nn+1;
	VECTOR(*weights)[half_actual_size++] = 1;
#if DEBUG
	printf("Add edge from %d to %d with weight 1\n", i+1, nn+1);
#endif
      }
    }
    sum_cap[i] = sum_cap_node;
#if DEBUG
    printf("Added %lf for vertex %d to sum_cap\n", sum_cap_node, i);
#endif
  }
  
  /* shrink down the vector to real size */
  igraph_vector_resize(&node_to_node, half_actual_size*2);
  /* prepare weight vector to hold also source-node node-sink weigts */
  igraph_vector_resize(weights, half_actual_size+2*Network_size);
  N_inner_edges = half_actual_size;
  
  igraph_add_edges(&graph, &node_to_node, 0);
  igraph_vector_destroy(&node_to_node);
  /* Add Source to Node and Node to Source edges */
  /* We do this at the end so we can modify capacity vector easier */
  igraph_vector_t source_node_sink;
  igraph_vector_init(&source_node_sink, 4*Network_size);
  for( int i = 0; i < Network_size; ++i){
    /* Source - Node */
    VECTOR(source_node_sink)[i*4] = 0;
    VECTOR(source_node_sink)[i*4+1] = i+1;
    /* Node - Sink */
    VECTOR(source_node_sink)[i*4+2] = i+1;
    VECTOR(source_node_sink)[i*4+3] = Network_size+1;
  }
  igraph_add_edges(&graph, &source_node_sink, 0);
  igraph_vector_destroy(&source_node_sink);
  /* Done forming inner edges of graph */
  fclose(read);
  
  end_graph_t = clock();
  printf("Read in graph time: %fs\n", (double)(end_graph_t-start_graph_t)/
	 CLOCKS_PER_SEC);
  
  return graph;
}


/****************************************************
 * Adds source-node and node-sink edges as required *
 ****************************************************/
void
add_st( const double mu, igraph_vector_t *weights ){

#if DEBUG
  printf("Starting new add_st with mu=%.2lf \n", mu);
#endif

  c_st = 0; /* reminder of declaration */
  for( int i = 0; i < Network_size; ++i ){
    double w = -1.0*mu - mf[i] - sum_cap[i];
    if (w > 0){
      #if DEBUG
      printf("Added edge from %d to sink with weigth of %lf\n", i, w);
      #endif
      VECTOR(*weights)[N_inner_edges + 2*i] = 0;
      VECTOR(*weights)[N_inner_edges + 2*i + 1] = w;
    }
    else{
      #if DEBUG
      printf("Added edge from source to %d with weigth of %lf\n", i, -w);
      #endif
      VECTOR(*weights)[N_inner_edges + 2*i] = -w;
      VECTOR(*weights)[N_inner_edges + 2*i + 1] = 0;
      c_st += w;
    }
  }
  
#if DEBUG
  printf("# Capacity vector: \n", mu);
  for (int i = 0; i < igraph_vector_size(weights); ++i){
    printf("%lf ", (double)VECTOR(*weights)[i]);
  }
  NEWLINE
#endif

}


/*****************************************
 * Returns fluid density given mu and y  *
 *****************************************/
int
get_ro( const igraph_t *graph, igraph_vector_t *weights, const double mu ){

#if DEBUG_REBUILD_TIME
  clock_t start_graph_t, end_graph_t;
  start_graph_t = clock();
#endif

  add_st( mu, weights );

#if DEBUG_REBUILD_TIME
  printf("Recalculate outer edges time: %fs\n", (double)(clock()-start_graph_t)/
	 CLOCKS_PER_SEC);
#endif

#if DEBUG
  igraph_integer_t result;
  igraph_es_t all;
  igraph_es_all(&all, IGRAPH_EDGEORDER_ID);
  igraph_es_size(graph, &all, &result);
  printf("Capacity vector size is %d and number of edges is %d \n",
	 igraph_vector_size(weights), result);
  fflush(stdout);
  igraph_es_destroy(&all);
#endif

  igraph_real_t value;
  static igraph_vector_t source_part;
  static int first_time = 1;
  if (first_time-- == 1){
    igraph_vector_init(&source_part, Network_size);
  }
  /* int igraph_st_mincut(const igraph_t *graph, igraph_real_t *value,
                          igraph_vector_t *cut, igraph_vector_t *partition,
                          igraph_vector_t *partition2,
                          igraph_integer_t source, igraph_integer_t target,
                          const igraph_vector_t *capacity); */
  (void)igraph_st_mincut( graph, &value,
			  NULL, &source_part,
			  NULL,
			  0, Network_size + 1,
			  weights );

#if DEBUG
  printf("Source node is %d and sink node is %d", Network_size, Network_size+1);
  printf("Cut value is %lf", value);
  NEWLINE
    printf("Source size is %d\n",
	   igraph_vector_size(&source_part));
  for (int i = 0; i < igraph_vector_size(&source_part); ++i){
    printf("%d ", (int)VECTOR(source_part)[i]);
  }
#endif

 
  int Cells_filled = igraph_vector_size(&source_part)-1;

  return (int)Cells_filled;
}


/********************************************
 * Samples isotherm using dichotomy method  *
 ********************************************/
void
dichotomy(const igraph_t *graph, igraph_vector_t *weights,
	  const double mu_i, const double mu_f){
  int ro_i = get_ro(graph, weights, mu_i);
  int ro_f = get_ro(graph, weights, mu_f);

#if DEBUG2
  printf("Rho difference is %d, given mu_i is %f"
	 "and mu_f is %f as ro_i is %d and ro_f is %d",
	 (ro_f-ro_i), mu_i, mu_f, ro_i, ro_f);
  fflush(stdout);
  exit(1);
  NEWLINE
#endif
    if( (ro_f-ro_i) > 1){
      int ro_mid = get_ro(graph, weights, (mu_f+mu_i)/2.0);
      dichotomy_left(graph, weights, mu_i, (mu_f+mu_i)/2, ro_i, ro_mid);
      dichotomy_right(graph, weights, (mu_f+mu_i)/2, mu_f, ro_mid, ro_f);
    }
    else{
      printf("%.12f, %d\n", mu_i, ro_i);
    }
}

void
dichotomy_left(const igraph_t *graph, igraph_vector_t *weights,
	       const double mu_i, const double mu_f,
	       const int ro_i, const int ro_f){

  if( (mu_f - mu_i) < EPS_MU){
    printf("%.12f, %d\n", mu_f, MAX(ro_i, ro_f));
    return;
  }

  if ( (ro_f-ro_i) > 0){
    int ro_mid = get_ro(graph, weights, (mu_f+mu_i)/2.0);
    dichotomy_left(graph, weights, mu_i, (mu_f+mu_i)/2, ro_i, ro_mid);
    dichotomy_right(graph, weights, (mu_f+mu_i)/2, mu_f, ro_mid, ro_f);
  }
  else{
    printf("%.12f, %d\n", mu_i, ro_i);
  }

}


void
dichotomy_right(const igraph_t *graph, igraph_vector_t *weights,
		const double mu_i, const double mu_f,
		const int ro_i, const int ro_f){
  if( (mu_f - mu_i) < EPS_MU){
    printf("%.12f, %d\n", mu_f, MAX(ro_i, ro_f));
    return;
  }

  if ( (ro_f-ro_i) > 0){
    int ro_mid = get_ro(graph, weights, (mu_f+mu_i)/2.0);
    dichotomy_left(graph, weights, mu_i, (mu_f+mu_i)/2, ro_i, ro_mid);
    dichotomy_right(graph, weights, (mu_f+mu_i)/2, mu_f, ro_mid, ro_f);
  }
  else{
    printf("%.12f, %d\n", mu_i, ro_i);
  }
}


