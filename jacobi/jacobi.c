/* CS322 - Jacobi Project
 *
 * Connor Freitas
 * Theo Floor
 *
 * This is an implementation of 
 * Jacobi Iteration as a mean of
 * solving Laplace's Equation.
 *
 */

#define NUM_THREADS 8
#define EPSILON 0.0001
#define N 2048

#include <stdio.h>
#include <stdlib.h>
#include <pthread.h>
#include <semaphore.h>
#include <string.h>
#include <stdbool.h>
#include "jacobi.h"

void worker_init(worker_t* worker, int thread_number, int num_threads,
                 void* grid_init_ptr, void* newgrid_init_ptr, 
                 pthread_mutex_t* mutex_ptr, bool* finished_ptr,
                 void* barrier_ptr, int* arrived_ptr, int* finarr_ptr);
void * work (void* thread_args); 
void barrier_wait (worker_t* worker_args);
void read_floats (char* file_name, void* init_matrix);


void main(int argc, char* argv[]) {
  float** grid;
  float** newgrid;
  void* retval;
  bool finished = false;
  int finarr[NUM_THREADS];
  int arrived = 0;
  pthread_mutex_t mutex;
  pthread_mutex_init(&mutex, NULL);
  sem_t barrier[NUM_THREADS];
  
  grid = malloc(N * sizeof(float*));
  newgrid = malloc(N * sizeof(float*));
  for (int i = 0; i <= N; ++i) {
    grid[i] = malloc(N * sizeof(float*));
    newgrid[i] = malloc(N * sizeof(float*));
  }

  read_floats(argv[1], grid);
  
  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      newgrid[i][j] = grid[i][j];
    }
  }

  for (int i = 0; i < NUM_THREADS; ++i) {
    sem_init(&barrier[i], 0, 0);
  } 

  worker_t worker[NUM_THREADS];
  pthread_t thds[NUM_THREADS];
  for (int i = 0; i < NUM_THREADS; ++i) {
    worker_init(&worker[i], i, NUM_THREADS, grid, newgrid, 
                &mutex, &finished, barrier, &arrived, finarr);
    
    if (pthread_create(&thds[i], NULL, work, &worker[i])) {
      printf("error\n");
    }
  }
  for (int i = 0; i < NUM_THREADS; ++i) {
    pthread_join(thds[i], &retval);
  }
  free(grid);
  free(newgrid);
  printf("Done!\n");
}

void worker_init(worker_t* worker, int thread_number, int num_threads,
                 void* grid_init_ptr, void* newgrid_init_ptr, 
                 pthread_mutex_t* mutex_ptr, bool* finished_ptr,
                 void* barrier_ptr, int* arrived_ptr, int* finarr_ptr) {
  worker->startrow = (N * thread_number) / num_threads;
  worker->endrow = (N * (thread_number + 1)) / num_threads - 1;
  worker->thread_id = thread_number;
  worker->grid_ptr = (float **) grid_init_ptr;
  worker->newgrid_ptr = (float **) newgrid_init_ptr;
  worker->mutex = mutex_ptr;
  worker->finished = finished_ptr;  
  worker->barrier = (sem_t*) barrier_ptr;
  worker->arrived = arrived_ptr; 
  worker->finarr = finarr_ptr;
 
  if (worker->startrow == 0) {
    worker->startrow = 1;
  }

  if (worker->endrow == N) {
    worker->endrow = N - 1;
  }
} 

void * work (void* thread_args) {
  worker_t* worker_args = (worker_t*)thread_args; 
  while (!*(worker_args->finished)) {
    worker_args->maxdiff = 0.0;
    
    for (int i = worker_args->startrow; i < worker_args->endrow; ++i) {
      for (int j = 1; j < N - 1; ++j) {
        worker_args->newgrid_ptr[i][j] = (worker_args->grid_ptr[i - 1][j] + 
        worker_args->grid_ptr[i][j - 1] +
        worker_args->grid_ptr[i + 1][j] +
        worker_args->grid_ptr[i][j + 1]) * 0.25;
        
        if ((worker_args->newgrid_ptr[i][j] - worker_args->grid_ptr[i][j]) > 
          worker_args->maxdiff) {
          worker_args->maxdiff = worker_args->newgrid_ptr[i][j] - worker_args->grid_ptr[i][j];
        }
      }
    }

    barrier_wait(worker_args);
    float** tempgrid = worker_args->newgrid_ptr;
    worker_args->newgrid_ptr = worker_args->grid_ptr;
    worker_args->grid_ptr = tempgrid;
    if (worker_args->maxdiff < EPSILON) {
      worker_args->finarr[worker_args->thread_id] = 1;
    }      
    else {
      worker_args->finarr[worker_args->thread_id] = 0;
    }
    barrier_wait(worker_args);
  }
  pthread_exit(NULL);
}


void barrier_wait (worker_t* worker_args) {
  pthread_mutex_lock(worker_args->mutex); 
  if (*(worker_args->arrived) == NUM_THREADS - 1) {
    *(worker_args->finished) = true;
    for (int k = 0; k < NUM_THREADS; ++k) {
      if (!worker_args->finarr[k]) {
        *(worker_args->finished) = false;
      }
    }  
    for (int i = 0; i < NUM_THREADS; ++i) {
      sem_post(&(worker_args->barrier[i]));
    }
    *(worker_args->arrived) = 0;
    pthread_mutex_unlock(worker_args->mutex);
    sem_wait(&(worker_args->barrier[worker_args->thread_id]));
  }
  else {
    *(worker_args->arrived) += 1; 
    pthread_mutex_unlock(worker_args->mutex);
    sem_wait(&(worker_args->barrier[worker_args->thread_id]));
  }
}     
 

void read_floats (char* file_name, void* init_matrix) {
  FILE* file = fopen(file_name, "r");
  float f = 0.0;
  float** matrix = (float **) init_matrix;

  for (int i = 0; i < N; ++i) {
    for (int j = 0; j < N; ++j) {
      fscanf(file, "%f", &f);
      matrix[i][j] = f;

    }
  }
  fclose(file);
}
