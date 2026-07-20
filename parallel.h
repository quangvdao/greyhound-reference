#ifndef LABRADOR_PARALLEL_H
#define LABRADOR_PARALLEL_H

#include <stddef.h>

typedef void (*parallel_for_fn)(size_t index, void *context);

void parallel_for(size_t count, size_t work_per_task,
                  parallel_for_fn fn, void *context);

#endif
