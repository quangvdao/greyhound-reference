#include <pthread.h>
#include <stdatomic.h>
#include <stddef.h>
#include <stdlib.h>
#include <unistd.h>
#include "parallel.h"

typedef struct {
  _Atomic size_t next;
  size_t count;
  parallel_for_fn fn;
  void *context;
} parallel_job;

static _Thread_local int inside_parallel;

static void *parallel_worker(void *arg) {
  parallel_job *job = arg;
  size_t index;
  int previous = inside_parallel;

  inside_parallel = 1;
  while((index = atomic_fetch_add_explicit(&job->next,1,memory_order_relaxed)) < job->count)
    job->fn(index,job->context);
  inside_parallel = previous;
  return NULL;
}

static size_t parallel_thread_count(void) {
  const char *value = getenv("LATTICE_DOGS_THREADS");
  char *end;
  long count;

  if(value && *value) {
    count = strtol(value,&end,10);
    if(!*end && count > 0)
      return (size_t)count;
  }
  count = sysconf(_SC_NPROCESSORS_ONLN);
  return count > 0 ? (size_t)count : 1;
}

void parallel_for(size_t count, size_t work_per_task,
                  parallel_for_fn fn, void *context) {
  size_t i, created = 0;
  size_t threads = parallel_thread_count();

  if(threads > count) threads = count;
  if(inside_parallel || threads < 2 || work_per_task < 4096) {
    for(i=0;i<count;i++) fn(i,context);
    return;
  }

  parallel_job job = { .next = 0, .count = count, .fn = fn, .context = context };
  pthread_t workers[threads-1];
  for(i=0;i<threads-1;i++) {
    if(pthread_create(&workers[i],NULL,parallel_worker,&job))
      break;
    created++;
  }
  parallel_worker(&job);
  for(i=0;i<created;i++)
    pthread_join(workers[i],NULL);
}
