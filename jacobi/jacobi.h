#ifndef jacobi_header
#define jacobi_header


#include <pthread.h>
#include <semaphore.h>
#include <stdbool.h>

typedef
struct worker_st {
  int startrow;
  int endrow;
  int thread_id;
  double maxdiff;
  bool* finished;
  int* arrived;
  int* finarr;
  double** grid_ptr;
  double** newgrid_ptr;
  sem_t* barrier;
  pthread_mutex_t* mutex; 
} worker_t;

#endif
