#line 1 "isotropic-cpp.c"
#line 1 "<built-in>"
#line 1 "<command-line>"
#line 1 "/usr/include/stdc-predef.h"
#line 1 "<command-line>"
#line 1 "isotropic-cpp.c"
#if _XOPEN_SOURCE < 700
#undef _XOPEN_SOURCE
#define _XOPEN_SOURCE 700
#endif
#if _GNU_SOURCE
#include <stdint.h>
#include <string.h>
#include <fenv.h>
#endif


#line 1 "/home/davidtang/basilisk/src/common.h"
#include <stdlib.h>
#include <stdio.h>
#include <stddef.h>
#include <stdbool.h>
#include <stdarg.h>
#include <string.h>
#include <float.h>
#include <limits.h>
#ifndef assert
# include <assert.h>
#endif
#include <math.h>
#include <time.h>
#include <sys/time.h>
#include <sys/resource.h>

#if _OPENMP
# include <omp.h>
# define OMP(x) _Pragma(#x)
#elif 1

# define OMP(x)

# include <mpi.h>
static int mpi_rank, mpi_npe;
# define tid() mpi_rank
# define pid() mpi_rank
# define npe() mpi_npe

#else

# define OMP(x)

#endif

#if _CADNA
# include <cadna.h>
#endif

#if __cplusplus
# define delete delete_qcc
# define right right_qcc
# define left left_qcc
# define norm norm_qcc
# define new new_qcc
#endif

#define pi 3.14159265358979
#undef HUGE
#define HUGE ((double)1e30)
#define nodata HUGE
#define _NVARMAX 65536
#define is_constant(v) ((v).i >= _NVARMAX)
#define constant(v) (is_constant(v) ? _constant[(v).i - _NVARMAX] : nodata)

#define max(a,b) ((a) > (b) ? (a) : (b))
#define min(a,b) ((a) < (b) ? (a) : (b))
#define sq(x) ((x)*(x))
#define cube(x) ((x)*(x)*(x))
#define sign(x) ((x) > 0 ? 1 : -1)
#define noise() (1. - 2.*rand()/(double)RAND_MAX)
#define clamp(x,a,b) ((x) < (a) ? (a) : (x) > (b) ? (b) : (x))
#define swap(type,a,b) { type tmp = a; a = b; b = tmp; }
#define unmap(x,y)

#define trash(x)


#define systderr stderr
#define systdout stdout

#if 1
FILE * qstderr (void);
FILE * qstdout (void);
FILE * ferr = NULL, * fout = NULL;
#define not_mpi_compatible()\
do {\
  if (npe() > 1) {\
    fprintf (ferr, "%s() is not compatible with MPI (yet)\n", __func__);\
    exit (1);\
  }\
} while(0)\

#line 83

# define system(command) (pid() == 0 ? system(command) : 0)
#else
# define qstderr() stderr
# define qstdout() stdout
# define ferr stderr
# define fout stdout
# define not_mpi_compatible()
#endif



#define sysmalloc malloc
#define syscalloc calloc
#define sysrealloc realloc
#define sysfree free
#define systrdup strdup

#if MTRACE

struct {
  FILE * fp;
  size_t total, max;
  size_t overhead, maxoverhead;
  size_t nr;
  size_t startrss, maxrss;
  char * fname;
} pmtrace;

typedef struct {
  char * func, * file;
  size_t max, total;
  int line, id;
} pmfunc;

typedef struct {
  size_t id, size;
} pmdata;

static pmfunc * pmfuncs = NULL;
static int pmfuncn = 0;

static int pmfunc_index (const char * func, const char * file, int line)
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++)
    if (p->line == line && !strcmp(func, p->func) && !strcmp(file, p->file))
      return p->id;
  pmfuncn++;
  pmfuncs = (pmfunc *) sysrealloc (pmfuncs, pmfuncn*sizeof(pmfunc));
  p = &pmfuncs[pmfuncn - 1];
  memset (p, 0, sizeof(pmfunc));
  p->func = systrdup(func);
  p->file = systrdup(file);
  p->line = line;
  p->id = pmfuncn;
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "@ %d %s %s %d\n", pmfuncn, func, file, line);
  return pmfuncn;
}

static void pmfunc_trace (pmfunc * f, char c)
{
  if (pmtrace.fp)
    fprintf (pmtrace.fp, "%c %d %ld %ld %ld",
      c, f->id, pmtrace.nr, pmtrace.total, f->total);
#if _GNU_SOURCE
  if (pmtrace.nr % 1 == 0) {
    struct rusage usage;
    getrusage (RUSAGE_SELF, &usage);
    if (pmtrace.fp)
      fprintf (pmtrace.fp, " %ld", usage.ru_maxrss*1024);
    if (!pmtrace.nr)
      pmtrace.startrss = usage.ru_maxrss;
    if (usage.ru_maxrss > pmtrace.maxrss)
      pmtrace.maxrss = usage.ru_maxrss;
  }
#endif
  if (pmtrace.fp)
    fputc ('\n', pmtrace.fp);
  pmtrace.nr++;
}

static void * pmfunc_alloc (pmdata * d, size_t size,
       const char * func, const char * file, int line,
       char c)
{
  assert (d != NULL);
  OMP (omp critical)
  {
    d->id = pmfunc_index(func, file, line);
    d->size = size;
    pmfunc * f = &pmfuncs[d->id - 1];
    f->total += size;
    if (f->total > f->max)
      f->max = f->total;
    pmtrace.total += size;
    pmtrace.overhead += sizeof(pmdata);
    if (pmtrace.total > pmtrace.max) {
      pmtrace.max = pmtrace.total;
      pmtrace.maxoverhead = pmtrace.overhead;
    }
    pmfunc_trace (f, c);
  }
  return ((char *)d) + sizeof(pmdata);
}

static void * pmfunc_free (void * ptr, char c)
{
  if (!ptr)
    return ptr;
  pmdata * d = (pmdata *) (((char *)ptr) - sizeof(pmdata));
  if (d->id < 1 || d->id > pmfuncn) {
    fputs ("*** MTRACE: ERROR!: corrupted free()", ferr);
    if (d->size == 0)
      fputs (", possible double free()", ferr);
    else
      fputs (", not traced?", ferr);
    fputs (", aborting...\n", ferr);
    abort();
    return ptr;
  }
  else
  OMP (omp critical)
  {
    pmfunc * f = &pmfuncs[d->id - 1];
    if (f->total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        f->total, d->size);
      abort();
    }
    else
      f->total -= d->size;
    if (pmtrace.total < d->size) {
      fprintf (ferr, "*** MTRACE: ERROR!: %ld < %ld: corrupted free()?\n",
        pmtrace.total, d->size);
      abort();
    }
    else {
      pmtrace.total -= d->size;
      pmtrace.overhead -= sizeof(pmdata);
    }
    d->id = 0;
    d->size = 0;
    pmfunc_trace (f, c);
  }
  return d;
}

static void * pmalloc (size_t size,
         const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysmalloc (sizeof(pmdata) + size),
         size, func, file, line, '+');
}

static void * pcalloc (size_t nmemb, size_t size,
         const char * func, const char * file, int line)
{
  void * p = pmalloc (nmemb*size, func, file, line);
  return memset (p, 0, nmemb*size);
}

static void * prealloc (void * ptr, size_t size,
   const char * func, const char * file, int line)
{
  return pmfunc_alloc ((pmdata *) sysrealloc (pmfunc_free(ptr, '<'),
           sizeof(pmdata) + size),
         size, func, file, line, '>');
}

static void pfree (void * ptr,
     const char * func, const char * file, int line)
{
  sysfree (pmfunc_free (ptr, '-'));
}

static char * pstrdup (const char * s,
         const char * func, const char * file, int line)
{
  char * d = (char *) pmalloc (strlen(s) + 1, func, file, line);
  return strcpy (d, s);
}

#if MTRACE < 3
static int pmaxsort (const void * a, const void * b) {
  const pmfunc * p1 = a, * p2 = b;
  return p1->max < p2->max;
}
#endif

static int ptotalsort (const void * a, const void * b) {
  const pmfunc * p1 = (const pmfunc *) a, * p2 = (const pmfunc *) b;
  return p1->total < p2->total;
}

static void pmfuncs_free()
{
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn; i++, p++) {
    sysfree (p->func);
    sysfree (p->file);
  }
  sysfree (pmfuncs);
}

void pmuntrace (void)
{
#if MTRACE < 3
  fprintf (ferr,
    "*** MTRACE: max resident  set size: %10ld bytes\n"
    "*** MTRACE: max traced memory size: %10ld bytes"
    " (tracing overhead %.1g%%)\n"
    "%10s    %20s   %s\n",
    pmtrace.maxrss*1024,
    pmtrace.max, pmtrace.maxoverhead*100./pmtrace.max,
    "max bytes", "function", "file");
  qsort (pmfuncs, pmfuncn, sizeof(pmfunc), pmaxsort);
  pmfunc * p = pmfuncs;
  for (int i = 0; i < pmfuncn && p->max > 0; i++, p++)
    fprintf (ferr, "%10ld    %20s   %s:%d\n",
      p->max, p->func, p->file, p->line);

  if (pmtrace.fp) {
    char * fname = pmtrace.fname, * s;
    while ((s = strchr(fname,'/')))
      fname = s + 1;

    fputs ("load(\"`echo $BASILISK`/mtrace.plot\")\n", pmtrace.fp);
    fprintf (pmtrace.fp,
      "plot '%s' u 3:($6-%g) w l t 'ru_maxrss - %.3g',"
      "total(\"%s\") w l t 'total'",
      fname,
      pmtrace.startrss*1024.,
      pmtrace.startrss*1024.,
      fname);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->max > 0.01*pmtrace.max; i++, p++)
      fprintf (pmtrace.fp,
        ",func(\"%s\",%d) w l t '%s'",
        fname, p->id, p->func);
    fputc ('\n', pmtrace.fp);
    fprintf (ferr,
      "*** MTRACE: To get a graph use: tail -n 2 %s | gnuplot -persist\n",
      fname);
    fclose (pmtrace.fp);
    pmtrace.fp = NULL;
    sysfree (pmtrace.fname);
  }
#endif

  if (pmtrace.total > 0) {
    qsort (pmfuncs, pmfuncn, sizeof(pmfunc), ptotalsort);
    pmfunc * p = pmfuncs;
    for (int i = 0; i < pmfuncn && p->total > 0; i++, p++)
      fprintf (ferr, "%s:%d: error: %ld bytes leaked here\n",
        p->file, p->line, p->total);
    pmfuncs_free();
    exit(1);
  }
  else {
#if MTRACE < 3
    fputs ("*** MTRACE: No memory leaks\n", ferr);
#endif
    pmfuncs_free();
  }
}

#else
# define pmalloc(s,func,file,line) malloc(s)
# define pcalloc(n,s,func,file,line) calloc(n,s)
# define prealloc(p,s,func,file,line) realloc(p,s)
# define pfree(p,func,file,line) free(p)
# define pstrdup(s,func,file,line) strdup(s)
#endif







typedef struct {
  void * p;
  long max, len;
} Array;

Array * array_new()
{
  Array * a = ((Array *) pmalloc ((1)*sizeof(Array),__func__,__FILE__,__LINE__));
  a->p = NULL;
  a->max = a->len = 0;
  return a;
}

void array_free (Array * a)
{
  pfree (a->p,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
}

void array_append (Array * a, void * elem, size_t size)
{
  if (a->len + size >= a->max) {
    a->max += max (size, 4096);
    a->p = prealloc (a->p, a->max,__func__,__FILE__,__LINE__);
  }
  memcpy (((char *)a->p) + a->len, elem, size);
  a->len += size;
}

void * array_shrink (Array * a)
{
  void * p = prealloc (a->p, a->len,__func__,__FILE__,__LINE__);
  pfree (a,__func__,__FILE__,__LINE__);
  return p;
}



#if TRACE == 1
#include <extrae_user_events.h>

typedef struct {
  Array index, stack;
  extrae_type_t type;
} Trace;

Trace trace_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000010,
};

Trace trace_mpi_func = {
  {NULL, 0, 0}, {NULL, 0, 0},
  60000011,
};

static int lookup_func (Array * a, const char * func)
{
  for (int i = 0; i < a->len/sizeof(char *); i++) {
    char * s = ((char **)a->p)[i];
    if (!strcmp (func, s))
      return i + 1;
  }
  char * s = pstrdup (func,__func__,__FILE__,__LINE__);
  array_append (a, &s, sizeof(char *));
  return a->len;
}

static void trace_push (Trace * t, const char * func)
{
  int value = lookup_func (&t->index, func);
  Extrae_eventandcounters (t->type, value);
  array_append (&t->stack, &value, sizeof(int));
}

static void trace_pop (Trace * t, const char * func)
{
  assert (t->stack.len > 0);
  t->stack.len -= sizeof(int);
  int value = t->stack.len > 0 ?
    ((int *)t->stack.p)[t->stack.len/sizeof(int) - 1] : 0;
  Extrae_eventandcounters (t->type, value);
}

static void trace_define (Trace * t, char * description)
{
  if (t->index.len > 0) {
    extrae_value_t values[t->index.len/sizeof(char *) + 1];
    char * names[t->index.len/sizeof(char *) + 1],
      ** func = (char **) t->index.p;
    names[0] = "OTHER";
    values[0] = 0;
    unsigned len = 1;
    for (int i = 0; i < t->index.len/sizeof(char *); i++, func++) {
      names[len] = *func;
      values[len++] = i + 1;
    }
    Extrae_define_event_type (&t->type, description, &len, values, names);
  }
}

static void trace_free (Trace * t)
{
  char ** func = (char **) t->index.p;
  for (int i = 0; i < t->index.len/sizeof(char *); i++, func++)
    pfree (*func,__func__,__FILE__,__LINE__);
  pfree (t->index.p,__func__,__FILE__,__LINE__);
  pfree (t->stack.p,__func__,__FILE__,__LINE__);
}

static void trace_off()
{
  trace_define (&trace_func, "Basilisk functions");
  trace_define (&trace_mpi_func, "Basilisk functions (MPI-related)");
  trace_free (&trace_func);
  trace_free (&trace_mpi_func);
}






# define trace(func, file, line) trace_push (&trace_func, func)
# define end_trace(func, file, line) trace_pop (&trace_func, func)

#elif TRACE

typedef struct {
  char * func, * file;
  int line, calls;
  double total, self;
#if 1
  double min, max;
#endif
} TraceIndex;

struct {
  Array stack, index;
  double t0;
} Trace = {
  {NULL, 0, 0}, {NULL, 0, 0},
  -1
};

static void trace_add (const char * func, const char * file, int line,
         double total, double self)
{
  TraceIndex * t = (TraceIndex *) Trace.index.p;
  int i, len = Trace.index.len/sizeof(TraceIndex);
  for (i = 0; i < len; i++, t++)
    if (t->line == line && !strcmp (func, t->func) && !strcmp (file, t->file))
      break;
  if (i == len) {
    TraceIndex t = {pstrdup(func,__func__,__FILE__,__LINE__), pstrdup(file,__func__,__FILE__,__LINE__), line, 1, total, self};
    array_append (&Trace.index, &t, sizeof(TraceIndex));
  }
  else
    t->calls++, t->total += total, t->self += self;
}

static void trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  if (Trace.t0 < 0)
    Trace.t0 = tv.tv_sec + tv.tv_usec/1e6;
  double t[2] = {(tv.tv_sec - Trace.t0) + tv.tv_usec/1e6, 0.};
  array_append (&Trace.stack, t, 2*sizeof(double));




}

static void end_trace (const char * func, const char * file, int line)
{
  struct timeval tv;
  gettimeofday (&tv, NULL);
  double te = (tv.tv_sec - Trace.t0) + tv.tv_usec/1e6;
  double * t = (double *) Trace.stack.p;
  assert (Trace.stack.len >= 2*sizeof(double));
  t += Trace.stack.len/sizeof(double) - 2;
  Trace.stack.len -= 2*sizeof(double);
  double dt = te - t[0];




  trace_add (func, file, line, dt, dt - t[1]);
  if (Trace.stack.len >= 2*sizeof(double)) {
    t -= 2;
    t[1] += dt;
  }
}

static int compar_self (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  return t1->self < t2->self;
}

#if 1
static int compar_func (const void * p1, const void * p2)
{
  const TraceIndex * t1 = p1, * t2 = p2;
  if (t1->line != t2->line)
    return t1->line < t2->line;
  return strcmp (t1->file, t2->file);
}
#endif

void trace_print (FILE * fp, double threshold)
{
  int i, len = Trace.index.len/sizeof(TraceIndex);
  double total = 0.;
  TraceIndex * t;
  Array * index = array_new();
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    array_append (index, t, sizeof(TraceIndex)), total += t->self;
#if 1
  qsort (index->p, len, sizeof(TraceIndex), compar_func);
  double tot[len], self[len], min[len], max[len];
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    tot[i] = t->total, self[i] = t->self;
  MPI_Reduce (self, min, len, MPI_DOUBLE, MPI_MIN, 0, MPI_COMM_WORLD);
  MPI_Reduce (self, max, len, MPI_DOUBLE, MPI_MAX, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? self : MPI_IN_PLACE,
       self, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  MPI_Reduce (pid() ? tot : MPI_IN_PLACE,
       tot, len, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
  total = 0.;
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    t->total = tot[i]/npe(), t->self = self[i]/npe(),
      t->max = max[i], t->min = min[i], total += t->self;
#endif
  qsort (index->p, len, sizeof(TraceIndex), compar_self);
  fprintf (fp, "   calls    total     self   %% total   function\n");
  for (i = 0, t = (TraceIndex *) index->p; i < len; i++, t++)
    if (t->self*100./total > threshold) {
      fprintf (fp, "%8d   %6.2f   %6.2f     %4.1f%%",
        t->calls, t->total, t->self, t->self*100./total);
#if 1
      fprintf (fp, " (%4.1f%% - %4.1f%%)", t->min*100./total, t->max*100./total);
#endif
      fprintf (fp, "   %s():%s:%d\n", t->func, t->file, t->line);
    }
  fflush (fp);
  array_free (index);
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    t->calls = t->total = t->self = 0.;
}

static void trace_off()
{
  trace_print (fout, 0.);

  int i, len = Trace.index.len/sizeof(TraceIndex);
  TraceIndex * t;
  for (i = 0, t = (TraceIndex *) Trace.index.p; i < len; i++, t++)
    pfree (t->func,__func__,__FILE__,__LINE__), pfree (t->file,__func__,__FILE__,__LINE__);

  pfree (Trace.index.p,__func__,__FILE__,__LINE__);
  Trace.index.p = NULL;
  Trace.index.len = Trace.index.max = 0;

  pfree (Trace.stack.p,__func__,__FILE__,__LINE__);
  Trace.stack.p = NULL;
  Trace.stack.len = Trace.stack.max = 0;
}

#else
# define trace(...)
# define end_trace(...)
#endif



#if _OPENMP

#define tid() omp_get_thread_num()
#define pid() 0
#define npe() omp_get_num_threads()
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#elif 1

static bool in_prof = false;
static double prof_start, _prof;
#define prof_start(name)\
  assert (!in_prof); in_prof = true;\
  prof_start = MPI_Wtime();\

#line 658

#define prof_stop()\
  assert (in_prof); in_prof = false;\
  _prof = MPI_Wtime();\
  mpi_time += _prof - prof_start;\

#line 663


#if FAKE_MPI
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)
#else

int mpi_all_reduce0 (void *sendbuf, void *recvbuf, int count,
       MPI_Datatype datatype, MPI_Op op, MPI_Comm comm)
{ trace ("mpi_all_reduce0", "/home/davidtang/basilisk/src/common.h", 672);
  { int _ret =  MPI_Allreduce (sendbuf, recvbuf, count, datatype, op, comm); end_trace("mpi_all_reduce0", "/home/davidtang/basilisk/src/common.h", 673);  return _ret; }
 end_trace("mpi_all_reduce0", "/home/davidtang/basilisk/src/common.h", 674); }
#define mpi_all_reduce(v,type,op) {\
  prof_start ("mpi_all_reduce");\
  union { int a; float b; double c;} global;\
  mpi_all_reduce0 (&(v), &global, 1, type, op, MPI_COMM_WORLD);\
  memcpy (&(v), &global, sizeof (v));\
  prof_stop();\
}\

#line 682

#define mpi_all_reduce_double(v,op) {\
  prof_start ("mpi_all_reduce");\
  double global, tmp = v;\
  mpi_all_reduce0 (&tmp, &global, 1, MPI_DOUBLE, op, MPI_COMM_WORLD);\
  v = global;\
  prof_stop();\
}\

#line 690


#endif

#define QFILE FILE

FILE * qstderr (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "log-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systderr;
  }
  return fp;
}

FILE * qstdout (void)
{
  static QFILE * fp = NULL;
  if (!fp) {
    if (mpi_rank > 0) {
      char name[80];
      sprintf (name, "out-%d", mpi_rank);
      fp = fopen (name, "w");
    }
    else
      fp = systdout;
  }
  return fp;
}

static void finalize (void)
{
  MPI_Finalize();
}

void mpi_init()
{
  int initialized;
  MPI_Initialized (&initialized);
  if (!initialized) {
    MPI_Init (NULL, NULL);
    MPI_Comm_set_errhandler (MPI_COMM_WORLD, MPI_ERRORS_ARE_FATAL);
    atexit (finalize);
  }
  MPI_Comm_rank (MPI_COMM_WORLD, &mpi_rank);
  MPI_Comm_size (MPI_COMM_WORLD, &mpi_npe);
  srand (mpi_rank + 1);
  if (ferr == NULL) {
    if (mpi_rank > 0) {
      ferr = fopen ("/dev/null", "w");
      fout = fopen ("/dev/null", "w");
    }
    else {
      ferr = systderr;
      fout = systdout;
    }
    char * etrace = getenv ("MALLOC_TRACE"), name[80];
    if (etrace && mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      setenv ("MALLOC_TRACE", name, 1);
    }
#if MTRACE == 1
    etrace = getenv ("MTRACE");
    if (!etrace)
      etrace = "mtrace";
    if (mpi_rank > 0) {
      sprintf (name, "%s-%d", etrace, mpi_rank);
      pmtrace.fp = fopen (name, "w");
      pmtrace.fname = systrdup(name);
    }
    else {
      pmtrace.fp = fopen (etrace, "w");
      pmtrace.fname = systrdup(etrace);
    }
#endif
  }
}

#else

#define tid() 0
#define pid() 0
#define npe() 1
#define mpi_all_reduce(v,type,op)
#define mpi_all_reduce_double(v,op)

#endif

void init_solver()
{
#if _CADNA
  cadna_init (-1);
#endif
#if 1
  mpi_init();
#elif MTRACE == 1
  char * etrace = getenv ("MTRACE");
  pmtrace.fp = fopen (etrace ? etrace : "mtrace", "w");
  pmtrace.fname = systrdup (etrace ? etrace : "mtrace");
#endif
}

#define OMP_PARALLEL() OMP(omp parallel)

#define NOT_UNUSED(x) (void)(x)

#define VARIABLES ;
#define _index(a,m) (a.i)
#define val(a,k,l,m) data(k,l,m)[_index(a,m)]

double _val_higher_dimension = 0.;
#define _val_higher_dimension(x,a,b,c) _val_higher_dimension
#line 816 "/home/davidtang/basilisk/src/common.h"
#if (_GNU_SOURCE || __APPLE__) && !_OPENMP && !_CADNA
double undefined;
# if __APPLE__
# include <stdint.h>
# include "fp_osx.h"
# endif
# define enable_fpe(flags) feenableexcept (flags)
# define disable_fpe(flags) fedisableexcept (flags)
static void set_fpe (void) {
  int64_t lnan = 0x7ff0000000000001;
  assert (sizeof (int64_t) == sizeof (double));
  memcpy (&undefined, &lnan, sizeof (double));
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#else
# define undefined ((double) DBL_MAX)
# define enable_fpe(flags)
# define disable_fpe(flags)
static void set_fpe (void) {}
#endif


typedef struct {
  long n;
  long tn;
  int depth;
  int maxdepth;
} Grid;
Grid * grid = NULL;

double X0 = 0., Y0 = 0., Z0 = 0.;

double L0 = 1.;




int N = 16;


typedef struct { int i; } scalar;

typedef struct {
  scalar x;

  scalar y;


  scalar z;

} vector;

typedef struct {
  vector x;

  vector y;


  vector z;

} tensor;

struct { int x, y, z; } Period = {false, false, false};

typedef struct {
  double x, y, z;
} coord;
#line 895 "/home/davidtang/basilisk/src/common.h"
void normalize (coord * n)
{
  double norm = 0.;
  {
#line 898

    norm += sq(n->x);
#line 898

    norm += sq(n->y);
#line 898

    norm += sq(n->z);}
  norm = sqrt(norm);
  {
#line 901

    n->x /= norm;
#line 901

    n->y /= norm;
#line 901

    n->z /= norm;}
}

struct _origin { double x, y, z; };

void origin (struct _origin p) {
  X0 = p.x; Y0 = p.y; Z0 = p.z;
}

void size (double L) {
  L0 = L;
}

double zero (double s0, double s1, double s2) { return 0.; }
#line 924 "/home/davidtang/basilisk/src/common.h"
  enum { right, left, top, bottom, front, back };

int nboundary = 2*3;



#define dirichlet(expr) (2.*(expr) - val(_s,0,0,0))
#define dirichlet_homogeneous() (- val(_s,0,0,0))
#define neumann(expr) (Delta*(expr) + val(_s,0,0,0))
#define neumann_homogeneous() (val(_s,0,0,0))

double * _constant = NULL;
extern size_t datasize;
typedef struct _Point Point;

#line 1 "/home/davidtang/basilisk/src/grid/boundaries.h"


typedef struct _Boundary Boundary;

struct _Boundary {
  void (* destroy) (Boundary * b);
  void (* level) (const Boundary * b, scalar * list, int l);

  void (* restriction) (const Boundary * b, scalar * list, int l);
};

static Boundary ** boundaries = NULL;

void add_boundary (Boundary * b) {
  int len = 0;
  if (boundaries) {
    Boundary ** i = boundaries;
    while (*i++) len++;
  }
  boundaries = (Boundary * *) prealloc (boundaries, (len + 2)*sizeof(Boundary *),__func__,__FILE__,__LINE__);
  boundaries[len] = b;
  boundaries[len+1] = NULL;
}

void free_boundaries() {
  if (!boundaries)
    return;
  Boundary ** i = boundaries, * b;
  while ((b = *i++))
    if (b->destroy)
      b->destroy (b);
    else
      pfree (b,__func__,__FILE__,__LINE__);
  pfree (boundaries,__func__,__FILE__,__LINE__);
  boundaries = NULL;
}
#line 47 "/home/davidtang/basilisk/src/grid/boundaries.h"
typedef struct {
  Boundary parent;
  int d;
} BoxBoundary;
#line 940 "/home/davidtang/basilisk/src/common.h"



typedef struct {

#line 945 "/home/davidtang/basilisk/src/common.h"

  double (** boundary) (Point, Point, scalar, void *);
  double (** boundary_homogeneous) (Point, Point, scalar, void *);
  double (* gradient) (double, double, double);
  void (* delete) (scalar);
  char * name;
  struct {
    int x;

    int y;


    int z;

  } d;
  vector v;
  bool face, nodump, freed;
  int block;

#line 17 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

  void (* prolongation) (Point, scalar);
  void (* restriction) (Point, scalar);

} _Attributes;
_Attributes * _attribute;
#line 943 "/home/davidtang/basilisk/src/common.h"



























int list_len (scalar * list)
{
  if (!list) return 0;
  int ns = 0;
  if (list) for (scalar s = *list, *_i0 = list; ((scalar *)&s)->i >= 0; s = *++_i0) ns++;
  return ns;
}

scalar * list_append (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  list[len] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_prepend (scalar * list, scalar s)
{
  int len = list_len (list);
  list = (scalar *) prealloc (list, (len + 2)*sizeof(scalar),__func__,__FILE__,__LINE__);
  for (int i = len; i >= 1; i--)
    list[i] = list[i-1];
  list[0] = s;
  list[len + 1].i = -1;
  return list;
}

scalar * list_add (scalar * list, scalar s)
{
  if (list) for (scalar t = *list, *_i1 = list; ((scalar *)&t)->i >= 0; t = *++_i1)
    if (t.i == s.i)
      return list;
  return list_append (list, s);
}

int list_lookup (scalar * l, scalar s)
{
  if (l != NULL)
    if (l) for (scalar s1 = *l, *_i2 = l; ((scalar *)&s1)->i >= 0; s1 = *++_i2)
      if (s1.i == s.i)
 return true;
  return false;
}

scalar * list_copy (scalar * l)
{
  scalar * list = NULL;
  if (l != NULL)
    if (l) for (scalar s = *l, *_i3 = l; ((scalar *)&s)->i >= 0; s = *++_i3)
      list = list_append (list, s);
  return list;
}

scalar * list_concat (scalar * l1, scalar * l2)
{
  scalar * l3 = list_copy (l1);
  if (l2) for (scalar s = *l2, *_i4 = l2; ((scalar *)&s)->i >= 0; s = *++_i4)
    l3 = list_append (l3, s);
  return l3;
}

void list_print (scalar * l, FILE * fp)
{
  int i = 0;
  if (l) for (scalar s = *l, *_i5 = l; ((scalar *)&s)->i >= 0; s = *++_i5)
    fprintf (fp, "%s%s", i++ == 0 ? "{" : ",", _attribute[s.i].name);
  fputs (i > 0 ? "}\n" : "{}\n", fp);
}

int vectors_len (vector * list)
{
  if (!list) return 0;
  int nv = 0;
  if (list) for (vector v = *list, *_i6 = list; ((scalar *)&v)->i >= 0; v = *++_i6) nv++;
  return nv;
}

vector * vectors_append (vector * list, vector v)
{
  int len = vectors_len (list);
  list = (vector *) prealloc (list, (len + 2)*sizeof(vector),__func__,__FILE__,__LINE__);
  list[len] = v;
  list[len + 1] = (vector){{-1}};
  return list;
}

vector * vectors_add (vector * list, vector v)
{
  if (list) for (vector w = *list, *_i7 = list; ((scalar *)&w)->i >= 0; w = *++_i7) {
    bool id = true;
    {
#line 1061

      if (w.x.i != v.x.i)
 id = false;
#line 1061

      if (w.y.i != v.y.i)
 id = false;
#line 1061

      if (w.z.i != v.z.i)
 id = false;}
    if (id)
      return list;
  }
  return vectors_append (list, v);
}

vector * vectors_copy (vector * l)
{
  vector * list = NULL;
  if (l != NULL)
    if (l) for (vector v = *l, *_i8 = l; ((scalar *)&v)->i >= 0; v = *++_i8)
      list = vectors_append (list, v);
  return list;
}

vector * vectors_from_scalars (scalar * s)
{
  vector * list = NULL;
  while (s->i >= 0) {
    vector v;
    {
#line 1084
 {
      assert (s->i >= 0);
      v.x = *s++;
    }
#line 1084
 {
      assert (s->i >= 0);
      v.y = *s++;
    }
#line 1084
 {
      assert (s->i >= 0);
      v.z = *s++;
    }}
    list = vectors_append (list, v);
  }
  return list;
}

int tensors_len (tensor * list)
{
  if (!list) return 0;
  int nt = 0;
  if (list) for (tensor t = *list, *_i9 = list; ((scalar *)&t)->i >= 0; t = *++_i9) nt++;
  return nt;
}

tensor * tensors_append (tensor * list, tensor t)
{
  int len = tensors_len (list);
  list = (tensor *) prealloc (list, (len + 2)*sizeof(tensor),__func__,__FILE__,__LINE__);
  list[len] = t;
  list[len + 1] = (tensor){{{-1}}};
  return list;
}

tensor * tensors_from_vectors (vector * v)
{
  tensor * list = NULL;
  while (v->x.i >= 0) {
    tensor t;
    {
#line 1115
 {
      assert (v->x.i >= 0);
      t.x = *v++;
    }
#line 1115
 {
      assert (v->y.i >= 0);
      t.y = *v++;
    }
#line 1115
 {
      assert (v->z.i >= 0);
      t.z = *v++;
    }}
    list = tensors_append (list, t);
  }
  return list;
}

scalar * all = NULL;



scalar (* init_scalar) (scalar, const char *);
scalar (* init_vertex_scalar) (scalar, const char *);
vector (* init_vector) (vector, const char *);
tensor (* init_tensor) (tensor, const char *);
vector (* init_face_vector) (vector, const char *);





typedef struct _Event Event;
typedef int (* Expr) (int *, double *, Event *);

struct _Event {
  int last, nexpr;
  int (* action) (const int, const double, Event *);
  Expr expr[3];
  int * arrayi;
  double * arrayt;
  char * file;
  int line;
  char * name;
  double t;
  int i, a;
  void * data;
  Event * next;
};

static Event * Events = NULL;

int iter = 0, inext = 0;
double t = 0, tnext = 0;
void init_events (void);
void event_register (Event event);
void _init_solver (void);



#if 1
static double mpi_time = 0.;
#endif

typedef struct {
  clock_t c;
  struct timeval tv;
  double tm;
} timer;

timer timer_start (void)
{
  timer t;
  t.c = clock();
  gettimeofday (&t.tv, NULL);
#if 1
  t.tm = mpi_time;
#endif
  return t;
}

double timer_elapsed (timer t)
{
  struct timeval tvend;
  gettimeofday (&tvend, NULL);
  return ((tvend.tv_sec - t.tv.tv_sec) +
   (tvend.tv_usec - t.tv.tv_usec)/1e6);
}



vector zerof= {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}};
vector unityf= {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
scalar unity= {_NVARMAX + 6};
scalar zeroc= {_NVARMAX + 7};



 vector fm = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar cm = {(_NVARMAX + 6)};
#line 1219 "/home/davidtang/basilisk/src/common.h"
static FILE ** qpopen_pipes = NULL;

FILE * qpopen (const char * command, const char * type)
{
  if (pid() > 0)
    return fopen ("/dev/null", type);
  FILE * fp = popen (command, type);
  if (fp) {
    FILE ** i = qpopen_pipes;
    int n = 0;
    while (i && *i) { n++; i++; }
    qpopen_pipes = (FILE * *) prealloc (qpopen_pipes, (n + 2)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    qpopen_pipes[n] = fp;
    qpopen_pipes[n+1] = NULL;
  }
  return fp;
}

int qpclose (FILE * fp)
{
  if (pid() > 0)
    return fclose (fp);
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i == fp)
      *i = (FILE *) 1;
    i++;
  }
  return pclose (fp);
}

static void qpclose_all()
{
  FILE ** i = qpopen_pipes;
  while (i && *i) {
    if (*i != (FILE *) 1)
      pclose (*i);
    i++;
  }
  pfree (qpopen_pipes,__func__,__FILE__,__LINE__);
  qpopen_pipes = NULL;
}






FILE * lfopen (const char * name, const char * mode)
{
  char fname[80];
  sprintf (fname, "%s-%d", name, pid());
  return fopen (fname, mode);
}



void * matrix_new (int n, int p, size_t size)
{
  void ** m = ((void * *) pmalloc ((n)*sizeof(void *),__func__,__FILE__,__LINE__));
  char * a = ((char *) pmalloc ((n*p*size)*sizeof(char),__func__,__FILE__,__LINE__));
  for (int i = 0; i < n; i++)
    m[i] = a + i*p*size;
  return m;
}

double matrix_inverse (double ** m, int n, double pivmin)
{
  int indxc[n], indxr[n], ipiv[n];
  int i, icol = 0, irow = 0, j, k, l, ll;
  double big, dum, pivinv, minpiv = HUGE;

  for (j = 0; j < n; j++)
    ipiv[j] = -1;

  for (i = 0; i < n; i++) {
    big = 0.0;
    for (j = 0; j < n; j++)
      if (ipiv[j] != 0)
 for (k = 0; k < n; k++) {
   if (ipiv[k] == -1) {
     if (fabs (m[j][k]) >= big) {
       big = fabs (m[j][k]);
       irow = j;
       icol = k;
     }
   }
 }
    ipiv[icol]++;
    if (irow != icol)
      for (l = 0; l < n; l++)
 swap (double, m[irow][l], m[icol][l]);
    indxr[i] = irow;
    indxc[i] = icol;
    if (fabs (m[icol][icol]) <= pivmin)
      return 0.;
    if (fabs (m[icol][icol]) < minpiv)
      minpiv = fabs (m[icol][icol]);
    pivinv = 1.0/m[icol][icol];
    m[icol][icol] = 1.0;
    for (l = 0; l < n; l++) m[icol][l] *= pivinv;
    for (ll = 0; ll < n; ll++)
      if (ll != icol) {
 dum = m[ll][icol];
 m[ll][icol] = 0.0;
 for (l = 0; l < n; l++)
   m[ll][l] -= m[icol][l]*dum;
      }
  }
  for (l = n - 1; l >= 0; l--) {
    if (indxr[l] != indxc[l])
      for (k = 0; k < n; k++)
 swap (double, m[k][indxr[l]], m[k][indxc[l]]);
  }
  return minpiv;
}

void matrix_free (void * m)
{
  pfree (((void **) m)[0],__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
}



typedef void (* free_solver_func) (void);

static Array * free_solver_funcs = NULL;

void free_solver_func_add (free_solver_func func)
{
  if (!free_solver_funcs)
    free_solver_funcs = array_new();
  array_append (free_solver_funcs, &func, sizeof(free_solver_func));
}



static char * display_defaults = NULL;

struct _display {
  const char * commands;
  bool overwrite;
};

static void free_display_defaults() {
  pfree (display_defaults,__func__,__FILE__,__LINE__);
}

void display (struct _display p)
{
  if (display_defaults == NULL)
    free_solver_func_add (free_display_defaults);
  if (p.overwrite) {
    pfree (display_defaults,__func__,__FILE__,__LINE__);
    display_defaults = pmalloc (strlen(p.commands) + 2,__func__,__FILE__,__LINE__);
    strcpy (display_defaults, "@");
    strcat (display_defaults, p.commands);
  }
  else {
    if (!display_defaults)
      display_defaults = pstrdup ("@",__func__,__FILE__,__LINE__);
    display_defaults =
      prealloc (display_defaults,
        strlen(display_defaults) + strlen(p.commands) + 1,__func__,__FILE__,__LINE__);
    strcat (display_defaults, p.commands);
  }
}
#line 13 "isotropic-cpp.c"
static double _boundary0 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary1 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary2 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary3 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary4 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s, void * data);
static double _boundary5 (Point point, Point neighbor, scalar _s, void * data);
static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s, void * data);
#line 1 "isotropic.c"
#line 13 "isotropic.c"
#line 1 "grid/multigrid3D.h"
#line 1 "/home/davidtang/basilisk/src/grid/multigrid3D.h"

#line 1 "grid/multigrid.h"
#line 1 "/home/davidtang/basilisk/src/grid/multigrid.h"







# define BGHOSTS 1







typedef struct {
  Grid g;
  char ** d;
} Multigrid;

struct _Point {
  int i;

  int j;


  int k;

  int level, n;
#ifdef foreach_block
  int l;
#define _BLOCK_INDEX , point.l
#else
#define _BLOCK_INDEX
#endif
};
static Point last_point;
#line 49 "/home/davidtang/basilisk/src/grid/multigrid.h"
static size_t _size (size_t l)
{
  size_t n = (1 << l) + 2*2;
  return cube(n);
}
#line 67 "/home/davidtang/basilisk/src/grid/multigrid.h"
#define data(l,m,o)\
  ((double *)&((Multigrid *)grid)->d[point.level][((point.i + l)*sq((1 << point.level) +\
         2*2) +\
      (point.j + m)*((1 << point.level) +\
       2*2) +\
      (point.k + o))*datasize]) 
#line 71

#line 101 "/home/davidtang/basilisk/src/grid/multigrid.h"
#define allocated(a,l,m) (point.i+a >= 0 &&\
         point.i+a < (1 << point.level) + 2*2 &&\
         point.j+l >= 0 &&\
         point.j+l < (1 << point.level) + 2*2 &&\
         point.k+m >= 0 &&\
         point.k+m < (1 << point.level) + 2*2)\

#line 107


#define allocated_child(a,l,m) (level < depth() &&\
         point.i > 0 && point.i <= (1 << point.level) + 2 &&\
         point.j > 0 && point.j <= (1 << point.level) + 2 &&\
         point.k > 0 && point.k <= (1 << point.level) + 2)\

#line 113




#define depth() (grid->depth)
#line 159 "/home/davidtang/basilisk/src/grid/multigrid.h"
#define fine(a,l,m,o)\
((double *)\
 &((Multigrid *)grid)->d[point.level+1][((2*point.i-2 +l)*sq(2*((1 << point.level) +\
         2)) +\
          (2*point.j-2 +m)*2*((1 << point.level) +\
             2) +\
          (2*point.k-2 +o))*datasize])[_index(a,m)]\

#line 166

#define coarse(a,l,m,o)\
((double *)\
 &((Multigrid *)grid)->d[point.level-1][(((point.i+2)/2+l)*sq((1 << point.level)/2 +\
        2*2) +\
          ((point.j+2)/2+m)*((1 << point.level)/2 +\
             2*2) +\
          (point.k+2)/2+o)*datasize])[_index(a,m)]\

#line 174

#define POINT_VARIABLES\
  VARIABLES\
  int level = point.level; NOT_UNUSED(level);\
  struct { int x, y, z; } child = {\
    2*((point.i + 2)%2) - 1,\
    2*((point.j + 2)%2) - 1,\
    2*((point.k + 2)%2) - 1\
  }; NOT_UNUSED(child);\
  Point parent = point; NOT_UNUSED(parent);\
  parent.level--;\
  parent.i = (point.i + 2)/2;\
  parent.j = (point.j + 2)/2;\
  parent.k = (point.k + 2)/2;\

#line 188



#define foreach_level(l)\
OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
      for (point.k = 2; point.k < point.n + 2; point.k++)\
\
 {\
\
          POINT_VARIABLES\

#line 208

#define end_foreach_level()\
\
 }\
\
  }\
}\

#line 215


#define foreach()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k < point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j < point.n + 2; point.j++)\
\
      for (point.k = 2; point.k < point.n + 2; point.k++)\
\
 {\
\
          POINT_VARIABLES\

#line 234

#define end_foreach()\
\
 }\
\
  }\
}\

#line 241


#define is_active(cell) (true)
#define is_leaf(cell) (level == depth())
#define is_local(cell) (true)
#define leaf 2
#define refine_cell(...) do {\
  fprintf (ferr, "grid depths do not match. Aborting.\n");\
  assert (0);\
} while (0)\

#line 251

#define tree ((Multigrid *)grid)
#line 1 "grid/foreach_cell.h"
#line 1 "/home/davidtang/basilisk/src/grid/foreach_cell.h"
#line 66 "/home/davidtang/basilisk/src/grid/foreach_cell.h"
#define foreach_cell_root(root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
 POINT_VARIABLES;\
\

#line 89

#define end_foreach_cell_root()\
        if (point.level < grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
          { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
        }\
        break;\
      }\
\
      case 1: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 2: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 3: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 4: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 5: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
      case 6: { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
       { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; }; break;\
      case 7: { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; }; break;\
\
      }\
    }\
  }\

#line 123


#define foreach_cell() {\
\
\
\
\
\
  Point root = {2,2,2,0};\
\
  foreach_cell_root (root)\

#line 134

#define end_foreach_cell() end_foreach_cell_root() }

#define foreach_cell_all() {\
  Point root = {0};\
  for (root.i = 2*Period.x; root.i <= 2*(2 - Period.x); root.i++)\
\
    for (root.j = 2*Period.y; root.j <= 2*(2 - Period.y); root.j++)\
\
\
      for (root.k = 2*Period.z; root.k <= 2*(2 - Period.z); root.k++)\
\
 foreach_cell_root (root)\

#line 147

#define end_foreach_cell_all() end_foreach_cell_root() }

#define foreach_cell_post_root(condition, root)\
  {\
    int ig = 0, jg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg);\
    Point point = {0};\
\
\
\
\
\
    int kg = 0; NOT_UNUSED(kg);\
    struct { int l, i, j, k, stage; } stack[20];\
\
    int _s = -1;\
    { _s++; stack[_s].l = 0; stack[_s].i = root.i; stack[_s].j = root.j; stack[_s].k = root.k; stack[_s].stage = 0; };\
    while (_s >= 0) {\
      int stage;\
      { point.level = stack[_s].l; point.i = stack[_s].i; point.j = stack[_s].j; point.k = stack[_s].k; stage = stack[_s].stage; _s--; };\
      if (!allocated (0,0,0))\
 continue;\
      switch (stage) {\
      case 0: {\
        POINT_VARIABLES;\
 if (point.level == grid->depth) {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 }\
 else {\
   { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 1; };\
   if (condition)\
     { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 }\
 break;\
      }\
\
      case 1:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 2; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 2:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 3; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 3:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 4; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = (2*point.i - 2); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 4:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 5; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 5:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 6; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = (2*point.j - 2); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
      case 6:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 7; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = (2*point.k - 2); stack[_s].stage = 0; };\
 break;\
      case 7:\
 { _s++; stack[_s].l = point.level; stack[_s].i = point.i; stack[_s].j = point.j; stack[_s].k = point.k; stack[_s].stage = 8; };\
 if (condition)\
   { _s++; stack[_s].l = point.level + 1; stack[_s].i = ((2*point.i - 2) + 1); stack[_s].j = ((2*point.j - 2) + 1); stack[_s].k = ((2*point.k - 2) + 1); stack[_s].stage = 0; };\
 break;\
\
      default: {\
        POINT_VARIABLES;\
\

#line 244

#define end_foreach_cell_post_root()\
      }\
      }\
    }\
  }\

#line 250


#define foreach_cell_post(condition)\
  {\
\
\
\
\
\
    Point root = {2,2,2,0};\
\
    foreach_cell_post_root(condition, root)\

#line 262

#define end_foreach_cell_post() end_foreach_cell_post_root() }

#define foreach_cell_post_all(condition) {\
  Point root = {0};\
  for (root.i = 0; root.i <= 2*2; root.i++)\
\
    for (root.j = 0; root.j <= 2*2; root.j++)\
\
\
      for (root.k = 0; root.k <= 2*2; root.k++)\
\
 foreach_cell_post_root (condition, root)\

#line 275

#define end_foreach_cell_post_all() end_foreach_cell_post_root() }

#define foreach_leaf() foreach_cell()\
  if (is_leaf (cell)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 281

#define end_foreach_leaf() } continue; } end_foreach_cell()
#line 254 "/home/davidtang/basilisk/src/grid/multigrid.h"

#define foreach_face_generic()\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = depth(); point.n = 1 << point.level;\
  int _k;\
  OMP(omp for schedule(static))\
  for (_k = 2; _k <= point.n + 2; _k++) {\
    point.i = _k;\
\
    for (point.j = 2; point.j <= point.n + 2; point.j++)\
\
      for (point.k = 2; point.k <= point.n + 2; point.k++)\
\
        {\
\
   POINT_VARIABLES\

#line 272

#define end_foreach_face_generic()\
\
 }\
\
  }\
}\

#line 279


#define foreach_vertex()\
foreach_face_generic() {\
  x -= Delta/2.;\
\
  y -= Delta/2.;\
\
\
  z -= Delta/2.;\
\

#line 290

#define end_foreach_vertex() } end_foreach_face_generic()

#define is_coarse() (point.level < depth())
#line 342 "/home/davidtang/basilisk/src/grid/multigrid.h"
#define foreach_vertex_aux()\
foreach_vertex() {\
  struct { int x, y, z; } _a = {point.i, point.j, point.k};\

#line 345

#define end_foreach_vertex_aux() } end_foreach_vertex()






#define is_face_x() (point.j < point.n + 2 && point.k < point.n + 2)
#define is_face_y() (point.i < point.n + 2 && point.k < point.n + 2)
#define is_face_z() (point.i < point.n + 2 && point.j < point.n + 2)

#define foreach_child() {\
  int _i = 2*point.i - 2;\
  int _j = 2*point.j - 2;\
  int _k = 2*point.k - 2;\
  point.level++;\
  point.n *= 2;\
  for (int _l = 0; _l < 2; _l++)\
    for (int _m = 0; _m < 2; _m++)\
      for (int _n = 0; _n < 2; _n++) {\
 point.i = _i + _l; point.j = _j + _m; point.k = _k + _n;\
 POINT_VARIABLES;\

#line 368

#define end_foreach_child()\
  }\
  point.i = (_i + 2)/2;\
  point.j = (_j + 2)/2;\
  point.k = (_k + 2)/2;\
  point.level--;\
  point.n /= 2;\
}\

#line 377

#define foreach_child_break() _l = _m = _n = 2


#if TRASH
# undef trash
# define trash(list) reset(list, undefined)
#endif

#line 1 "grid/neighbors.h"
#line 1 "/home/davidtang/basilisk/src/grid/neighbors.h"
#line 35 "/home/davidtang/basilisk/src/grid/neighbors.h"
#define foreach_neighbor(_s) {\
  int _nn = _s + 0 ? _s + 0 : 2;\
  int _i = point.i, _j = point.j, _k = point.k;\
  for (int _l = - _nn; _l <= _nn; _l++) {\
    point.i = _i + _l;\
    for (int _m = - _nn; _m <= _nn; _m++) {\
      point.j = _j + _m;\
      for (int _n = - _nn; _n <= _nn; _n++) {\
 point.k = _k + _n;\
 POINT_VARIABLES;\

#line 45

#define end_foreach_neighbor()\
      }\
    }\
  }\
  point.i = _i; point.j = _j; point.k = _k;\
}\

#line 52

#define foreach_neighbor_break() _l = _m = _n = _nn + 1
#line 387 "/home/davidtang/basilisk/src/grid/multigrid.h"

void reset (void * alist, double val)
{
  scalar * list = (scalar *) alist;
  Point p;
  p.level = depth(); p.n = 1 << p.level;
  for (; p.level >= 0; p.n /= 2, p.level--)
    for (int i = 0; i < cube(p.n + 2*2); i++)
      if (list) for (scalar s = *list, *_i10 = list; ((scalar *)&s)->i >= 0; s = *++_i10) {
 if (!is_constant(s))
   for (int b = 0; b < _attribute[s.i].block; b++)
     ((double *)(&((Multigrid *)grid)->d[p.level][i*datasize]))[s.i + b] = val;
      }
}
#line 473 "/home/davidtang/basilisk/src/grid/multigrid.h"
#define foreach_boundary_dir(l,d)\
  OMP_PARALLEL() {\
  int ig = 0, jg = 0, kg = 0; NOT_UNUSED(ig); NOT_UNUSED(jg); NOT_UNUSED(kg);\
  Point point = {0};\
  point.level = l < 0 ? depth() : l;\
  point.n = 1 << point.level;\
  int * _i = &point.j, * _j = &point.k;\
  if (d == left) {\
    point.i = 2;\
    ig = -1;\
  }\
  else if (d == right) {\
    point.i = point.n + 2 - 1;\
    ig = 1;\
  }\
  else if (d == bottom) {\
    point.j = 2;\
    _i = &point.i;\
    jg = -1;\
  }\
  else if (d == top) {\
    point.j = point.n + 2 - 1;\
    _i = &point.i;\
    jg = 1;\
  }\
  else if (d == back) {\
    point.k = 2;\
    _i = &point.i; _j = &point.j;\
    kg = -1;\
  }\
  else if (d == front) {\
    point.k = point.n + 2 - 1;\
    _i = &point.i; _j = &point.j;\
    kg = 1;\
  }\
  int _l;\
  OMP(omp for schedule(static))\
  for (_l = 0; _l < point.n + 2*2; _l++) {\
    *_i = _l;\
    for (int _m = 0; _m < point.n + 2*2; _m++) {\
      *_j = _m;\
      POINT_VARIABLES\

#line 515

#define end_foreach_boundary_dir()\
    }\
  }\
}\

#line 520


#define neighbor(o,p,q)\
  ((Point){point.i+o, point.j+p, point.k+q, point.level, point.n _BLOCK_INDEX})\

#line 524

#define is_boundary(point) (point.i < 2 || point.i >= point.n + 2 ||\
    point.j < 2 || point.j >= point.n + 2 ||\
    point.k < 2 || point.k >= point.n + 2)\

#line 528




#define foreach_boundary(b)\
  if (default_scalar_bc[b] != periodic_bc)\
    foreach_boundary_dir (depth(), b)\
      if (!is_boundary(point)) {\

#line 536

#define end_foreach_boundary() } end_foreach_boundary_dir()

#define neighborp(k,l,o) neighbor(k,l,o)

static double periodic_bc (Point point, Point neighbor, scalar s, void * data);

static void box_boundary_level (const Boundary * b, scalar * scalars, int l)
{
  disable_fpe (FE_DIVBYZERO|FE_INVALID);
  for (int bghost = 1; bghost <= BGHOSTS; bghost++)
    for (int d = 0; d < 2*3; d++) {

      scalar * list = NULL, * listb = NULL;
      if (scalars) for (scalar s = *scalars, *_i11 = scalars; ((scalar *)&s)->i >= 0; s = *++_i11)
 if (!is_constant(s) && _attribute[s.i].block > 0) {
   scalar sb = s;

   if (_attribute[s.i].v.x.i >= 0) {

     int j = 0;
     while ((&_attribute[s.i].v.x)[j].i != s.i) j++;
     sb = (&_attribute[s.i].v.x)[(j - d/2 + 3) % 3];
   }

   if (_attribute[sb.i].boundary[d] && _attribute[sb.i].boundary[d] != periodic_bc) {
     list = list_append (list, s);
     listb = list_append (listb, sb);
   }
 }

      if (list) {
  { foreach_boundary_dir (l, d){

#line 568 "/home/davidtang/basilisk/src/grid/multigrid.h"
 {
   scalar s, sb;
   scalar * _i0 = list; scalar * _i1 = listb; if (list) for (s = *list, sb = *listb; ((scalar *)&s)->i >= 0; s = *++_i0, sb = *++_i1) {
     if (_attribute[s.i].face && sb.i == _attribute[s.i].v.x.i) {

       if (bghost == 1)
 
    val(s,(ig + 1)/2,(jg + 1)/2,(kg + 1)/2) =
    _attribute[sb.i].boundary[d] (point, neighborp(ig,jg,kg), s, NULL);
     }
     else

      
  val(s,bghost*ig,bghost*jg,bghost*kg) =
  _attribute[sb.i].boundary[d] (neighborp((1 - bghost)*ig,
       (1 - bghost)*jg,
       (1 - bghost)*kg),
    neighborp(bghost*ig,bghost*jg,bghost*kg),
    s, NULL);
   }
 } } end_foreach_boundary_dir(); }
 pfree (list,__func__,__FILE__,__LINE__);
 pfree (listb,__func__,__FILE__,__LINE__);
      }
    }
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}
#line 703 "/home/davidtang/basilisk/src/grid/multigrid.h"
void free_grid (void)
{
  if (!grid)
    return;
  free_boundaries();
  Multigrid * m = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++)
    pfree (m->d[l],__func__,__FILE__,__LINE__);
  pfree (m->d,__func__,__FILE__,__LINE__);
  pfree (m,__func__,__FILE__,__LINE__);
  grid = NULL;
}

int log_base2 (int n) {
  int m = n, r = 0;
  while (m > 1)
    m /= 2, r++;
  return (1 << r) < n ? r + 1 : r;
}

void init_grid (int n)
{
  free_grid();
  Multigrid * m = ((Multigrid *) pmalloc ((1)*sizeof(Multigrid),__func__,__FILE__,__LINE__));
  grid = (Grid *) m;
  grid->depth = grid->maxdepth = log_base2(n);
  N = 1 << depth();

  grid->n = grid->tn = 1 << 3*depth();

  Boundary * b = ((Boundary *) pcalloc (1, sizeof(Boundary),__func__,__FILE__,__LINE__));
  b->level = box_boundary_level;
  add_boundary (b);

  Boundary * mpi_boundary_new();
  mpi_boundary_new();
#line 748 "/home/davidtang/basilisk/src/grid/multigrid.h"
  m->d = (char **) pmalloc(sizeof(Point *)*(depth() + 1),__func__,__FILE__,__LINE__);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l)*datasize;
    m->d[l] = (char *) pmalloc (len,__func__,__FILE__,__LINE__);


    double * v = (double *) m->d[l];
    for (int i = 0; i < len/sizeof(double); i++)
      v[i] = undefined;
  }
  reset (all, 0.);
}

void realloc_scalar (int size)
{
  Multigrid * p = ((Multigrid *)grid);
  for (int l = 0; l <= depth(); l++) {
    size_t len = _size(l);
    p->d[l] = (char *) prealloc (p->d[l], (len*(datasize + size))*sizeof(char),__func__,__FILE__,__LINE__);
    char * data = p->d[l] + (len - 1)*datasize;
    for (int i = len - 1; i > 0; i--, data -= datasize)
      memmove (data + i*size, data, datasize);
  }
  datasize += size;
}


int mpi_dims[3], mpi_coords[3];
#line 786 "/home/davidtang/basilisk/src/grid/multigrid.h"
struct _locate { double x, y, z; };

Point locate (struct _locate p)
{
  Point point = {0};  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 790 "/home/davidtang/basilisk/src/grid/multigrid.h"

  point.level = -1, point.n = 1 << depth();

  point.i = (p.x - X0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[0]*point.n;
  if (point.i < 2 || point.i >= point.n + 2)
    return point;

  point.j = (p.y - Y0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[1]*point.n;
  if (point.j < 2 || point.j >= point.n + 2)
    return point;


  point.k = (p.z - Z0)/L0*point.n*mpi_dims[0] + 2 - mpi_coords[2]*point.n;
  if (point.k < 2 || point.k >= point.n + 2)
    return point;
#line 821 "/home/davidtang/basilisk/src/grid/multigrid.h"
  point.level = depth();
  return point;
}

#line 1 "grid/multigrid-common.h"
#line 1 "/home/davidtang/basilisk/src/grid/multigrid-common.h"


#line 1 "grid/cartesian-common.h"
#line 1 "/home/davidtang/basilisk/src/grid/cartesian-common.h"
#line 1 "grid/events.h"
#line 1 "/home/davidtang/basilisk/src/grid/events.h"







static void event_error (Event * ev, const char * s)
{
  fprintf (ferr, "%s:%d: error: %s\n", ev->file, ev->line, s);
  exit (1);
}

static void init_event (Event * ev)
{
  if (ev->arrayi || ev->arrayt) {
    ev->i = ev->t = -1;
    if (ev->arrayi)
      ev->i = ev->arrayi[0];
    else
      ev->t = ev->arrayt[0];
    ev->a = 1;
    ev->expr[1] = NULL;
  }
  else {
    if (ev->nexpr > 0) {
      Expr init = NULL, cond = NULL, inc = NULL;
      for (int j = 0; j < ev->nexpr; j++) {
 int i = -123456; double t = i;
 (* ev->expr[j]) (&i, &t, ev);
 if (i == -123456 && t == -123456) {

   if (cond)
     event_error (ev, "events can only use a single condition");
   cond = ev->expr[j];
 }
 else {

   int i1 = i; double t1 = t;
   (* ev->expr[j]) (&i1, &t1, ev);
   if (i1 == i && t1 == t) {


     if (init)
       event_error (ev, "events can only use a single initialisation");
     init = ev->expr[j];
   }
   else {

     if (inc)
       event_error (ev, "events can only use a single increment");
     inc = ev->expr[j];
   }
 }
      }
      ev->expr[0] = init;
      ev->expr[1] = cond;
      ev->expr[2] = inc;
      ev->nexpr = 0;
    }
    ev->i = ev->t = -1;
    if (ev->expr[0]) {
      (* ev->expr[0]) (&ev->i, &ev->t, ev);
      if (ev->i == 1234567890 || ev->t == 1234567890) {
 ev->i = 1234567890; ev->t = -1;
      }
    }
    else if (ev->expr[2]) {
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (ev->i != -1)
 ev->i = 0;
      if (ev->t != -1)
 ev->t = 0;
    }
  }
}

enum { event_done, event_alive, event_stop };

static int event_finished (Event * ev)
{
  ev->t = ev->i = -1;
  return event_done;
}

void event_register (Event event) {
  assert (Events);
  assert (!event.last);
  int n = 0, parent = -1;
  for (Event * ev = Events; !ev->last; ev++) {
    if (!strcmp (event.name, ev->name)) {
      assert (parent < 0);
      parent = n;
    }
    n++;
  }
  if (parent < 0) {
    Events = (Event *) prealloc (Events, (n + 2)*sizeof(Event),__func__,__FILE__,__LINE__);
    Events[n] = event;
    Events[n].next = NULL;
    Events[n + 1].last = true;
    init_event (&Events[n]);
  }
  else {
    Event * ev = ((Event *) pcalloc (1, sizeof(Event),__func__,__FILE__,__LINE__));
    *ev = Events[parent];
    Events[parent] = event;
    Events[parent].next = ev;
    init_event (&Events[parent]);
  }
}

static int event_cond (Event * ev, int i, double t)
{
  if (!ev->expr[1])
    return true;
  return (* ev->expr[1]) (&i, &t, ev);
}
#line 131 "/home/davidtang/basilisk/src/grid/events.h"
static int event_do (Event * ev, bool action)
{
  if ((iter > ev->i && t > ev->t) || !event_cond (ev, iter, t))
    return event_finished (ev);
  if (iter == ev->i || fabs (t - ev->t) <= 1e-9) {
    if (action) {
      bool finished = false;
      for (Event * e = ev; e; e = e->next) {



 if ((* e->action) (iter, t, e))
   finished = true;
      }
      if (finished) {
 event_finished (ev);
 return event_stop;
      }
    }
    if (ev->arrayi) {
      ev->i = ev->arrayi[ev->a++];
      if (ev->i < 0)
 return event_finished (ev);
    }
    if (ev->arrayt) {
      ev->t = ev->arrayt[ev->a++];
      if (ev->t < 0)
 return event_finished (ev);
    }
    else if (ev->expr[2]) {
      int i0 = ev->i;
      (* ev->expr[2]) (&ev->i, &ev->t, ev);
      if (i0 == -1 && ev->i != i0)
 ev->i += iter + 1;
      if (!event_cond (ev, iter + 1, ev->t))
 return event_finished (ev);
    }
    else if (ev->expr[0] && !ev->expr[1])
      return event_finished (ev);
  }
  return event_alive;
}

static void end_event_do (bool action)
{




  for (Event * ev = Events; !ev->last; ev++)
    if (ev->i == 1234567890 && action)
      for (Event * e = ev; e; e = e->next) {



 e->action (iter, t, e);
      }
}

int events (bool action)
{





  if (iter == 0)
    for (Event * ev = Events; !ev->last; ev++)
      init_event (ev);

  int cond = 0, cond1 = 0;
  inext = 1234567890; tnext = HUGE;
  for (Event * ev = Events; !ev->last && !cond; ev++)
    if (ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond = 1;
  for (Event * ev = Events; !ev->last; ev++) {
    int status = event_do (ev, action);
    if (status == event_stop) {
      end_event_do (action);
      return 0;
    }
    if (status == event_alive && ev->i != 1234567890 &&
 (ev->expr[1] || (ev->expr[0] && !ev->expr[1] && !ev->expr[2]) || ev->arrayi || ev->arrayt))
      cond1 = 1;
    if (ev->t > t && ev->t < tnext)
      tnext = ev->t;
    if (ev->i > iter && ev->i < inext)
      inext = ev->i;
  }
  if ((!cond || cond1) && (tnext != HUGE || inext != 1234567890)) {
    inext = iter + 1;
    return 1;
  }
  end_event_do (action);
  return 0;
}

void event (const char * name)
{
  for (Event * ev = Events; !ev->last; ev++)
    if (!strcmp (ev->name, name))
      for (Event * e = ev; e; e = e->next) {



 (* e->action) (0, 0, e);
      }
}

double dtnext (double dt)
{
  if (tnext != HUGE && tnext > t) {
    unsigned int n = (tnext - t)/dt;
    assert (n < INT_MAX);
    if (n == 0)
      dt = tnext - t;
    else {
      double dt1 = (tnext - t)/n;
      if (dt1 > dt + 1e-9)
 dt = (tnext - t)/(n + 1);
      else if (dt1 < dt)
 dt = dt1;
      tnext = t + dt;
    }
  }
  else
    tnext = t + dt;
  return dt;
}
#line 2 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

void (* debug) (Point);

#define _val_constant(a,k,l,m) ((const double) _constant[a.i -_NVARMAX])

#undef VARIABLES
#define VARIABLES\
  double Delta = L0*(1./(1 << point.level)/mpi_dims[0]);\
  \
    double Delta_x = Delta;\
    double Delta_y = Delta;\
    double Delta_z = Delta;\
\
  double x = (ig/2. + (point.i - 2 + mpi_coords[0]*(1 << point.level)) + 0.5)*Delta + X0; NOT_UNUSED(x);\
\
  double y = (jg/2. + (point.j - 2 + mpi_coords[1]*(1 << point.level)) + 0.5)*Delta + Y0;\
\
\
\
 NOT_UNUSED(y);\
\
  double z = (kg/2. + (point.k - 2 + mpi_coords[2]*(1 << point.level)) + 0.5)*Delta + Z0;\
\
\
\
  NOT_UNUSED(z);\
\
  NOT_UNUSED(Delta);\
  \
    NOT_UNUSED(Delta_x);\
    NOT_UNUSED(Delta_y);\
    NOT_UNUSED(Delta_z);\
\
  ;\

#line 32


#line 1 "grid/fpe.h"
#line 1 "/home/davidtang/basilisk/src/grid/fpe.h"


#include <signal.h>
#include <unistd.h>

static int gdb()
{
  if (last_point.level >= 0) {
    debug (last_point);
    fputc ('\n', ferr);
    fflush (ferr);
  }
  char command[80];
  sprintf (command, "exec xterm -e 'gdb -p %d' & xterm -e 'gnuplot plot -'",
    getpid());
  return system (command);
}

static void caught_abort (int sig)
{
  fprintf (ferr, "Caught signal %d (Aborted)\n", sig);
  gdb();
}

static void caught_fpe (int sig)
{
  fprintf (ferr, "Caught signal %d (Floating Point Exception)\n", sig);
  gdb();
  exit (1);
}

static void caught_segfault (int sig)
{
  fprintf (ferr, "Caught signal %d (Segmentation Fault)\n", sig);
  gdb();
  exit (2);
}

void catch_fpe (void)
{
  struct sigaction act;
  act.sa_handler = caught_fpe;
  sigemptyset (&act.sa_mask);
  act.sa_flags = 0;
  last_point.level = -1;
  sigaction (8, &act, NULL);
  act.sa_handler = caught_segfault;
  sigaction (11, &act, NULL);
  act.sa_handler = caught_abort;
  act.sa_flags = SA_RESETHAND;
  sigaction (6, &act, NULL);
}
#line 35 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

#define end_foreach_face()

static void init_block_scalar (scalar sb, const char * name, const char * ext,
          int n, int block)
{
  char bname[strlen(name) + strlen(ext) + 10];
  if (n == 0) {
    sprintf (bname, "%s%s", name, ext);
    init_scalar (sb, bname);
    _attribute[sb.i].block = block;
  }
  else {
    sprintf (bname, "%s%d%s", name, n, ext);
    init_scalar (sb, bname);
    _attribute[sb.i].block = - n;
  }
  all = list_append (all, sb);
}

scalar new_block_scalar (const char * name, const char * ext, int block)
{
  int nvar = datasize/sizeof(double);

  scalar s = {0};
  while (s.i < nvar) {
    int n = 0;
    scalar sb = s;
    while (sb.i < nvar && n < block && _attribute[sb.i].freed)
      n++, sb.i++;
    if (n >= block) {
      for (sb.i = s.i, n = 0; n < block; n++, sb.i++)
 init_block_scalar (sb, name, ext, n, block);
      trash (((scalar []){s, {-1}}));
      return s;
    }
    s.i = sb.i + 1;
  }


  s = (scalar){nvar};
  assert (nvar + block <= _NVARMAX);
  _attribute = (_Attributes *) prealloc (_attribute, (nvar + block)*sizeof(_Attributes),__func__,__FILE__,__LINE__);
  memset (&_attribute[nvar], 0, block*sizeof (_Attributes));
  for (int n = 0; n < block; n++, nvar++) {
    scalar sb = (scalar){nvar};
    init_block_scalar (sb, name, ext, n, block);
  }

  realloc_scalar (block*sizeof(double));
  trash (((scalar []){s, {-1}}));
  return s;
}

scalar new_scalar (const char * name)
{
  return new_block_scalar (name, "", 1);
}

scalar new_vertex_scalar (const char * name)
{
  return init_vertex_scalar (new_scalar (name), name);
}

static vector alloc_block_vector (const char * name, int block)
{
  vector v;
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 103

    v.x = new_block_scalar (name, ext.x, block);
#line 103

    v.y = new_block_scalar (name, ext.y, block);
#line 103

    v.z = new_block_scalar (name, ext.z, block);}
  return v;
}

vector new_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_vector (v, NULL);
  return v;
}

vector new_face_vector (const char * name)
{
  vector v = alloc_block_vector (name, 1);
  init_face_vector (v, NULL);
  return v;
}

vector new_block_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    {
#line 127

      vb.x.i = v.x.i + i;
#line 127

      vb.y.i = v.y.i + i;
#line 127

      vb.z.i = v.z.i + i;}
    init_vector (vb, NULL);
    {
#line 130

      _attribute[vb.x.i].block = - i;
#line 130

      _attribute[vb.y.i].block = - i;
#line 130

      _attribute[vb.z.i].block = - i;}
  }
  {
#line 133

    _attribute[v.x.i].block = block;
#line 133

    _attribute[v.y.i].block = block;
#line 133

    _attribute[v.z.i].block = block;}
  return v;
}

vector new_block_face_vector (const char * name, int block)
{
  vector v = alloc_block_vector (name, block);
  for (int i = 0; i < block; i++) {
    vector vb;
    {
#line 143

      vb.x.i = v.x.i + i;
#line 143

      vb.y.i = v.y.i + i;
#line 143

      vb.z.i = v.z.i + i;}
    init_face_vector (vb, NULL);
    {
#line 146

      _attribute[vb.x.i].block = - i;
#line 146

      _attribute[vb.y.i].block = - i;
#line 146

      _attribute[vb.z.i].block = - i;}
  }
  {
#line 149

    _attribute[v.x.i].block = block;
#line 149

    _attribute[v.y.i].block = block;
#line 149

    _attribute[v.z.i].block = block;}
  return v;
}

tensor new_tensor (const char * name)
{
  char cname[strlen(name) + 3];
  struct { char * x, * y, * z; } ext = {"%s.x", "%s.y", "%s.z"};
  tensor t;
  {
#line 159
 {
    sprintf (cname, ext.x, name);
    t.x = new_vector (cname);
  }
#line 159
 {
    sprintf (cname, ext.y, name);
    t.y = new_vector (cname);
  }
#line 159
 {
    sprintf (cname, ext.z, name);
    t.z = new_vector (cname);
  }}
  init_tensor (t, NULL);
  return t;
}

tensor new_symmetric_tensor (const char * name)
{
  char cname[strlen(name) + 5];
  struct { char * x, * y, * z; } ext = {"%s.x.x", "%s.y.y", "%s.z.z"};
  tensor t;
  {
#line 172
 {
    sprintf (cname, ext.x, name);
    t.x.x = new_scalar(cname);
  }
#line 172
 {
    sprintf (cname, ext.y, name);
    t.y.y = new_scalar(cname);
  }
#line 172
 {
    sprintf (cname, ext.z, name);
    t.z.z = new_scalar(cname);
  }}

    sprintf (cname, "%s.x.y", name);
    t.x.y = new_scalar(cname);
    t.y.x = t.x.y;


    sprintf (cname, "%s.x.z", name);
    t.x.z = new_scalar(cname);
    t.z.x = t.x.z;
    sprintf (cname, "%s.y.z", name);
    t.y.z = new_scalar(cname);
    t.z.y = t.y.z;




  init_tensor (t, NULL);
  return t;
}

static int nconst = 0;

void init_const_scalar (scalar s, const char * name, double val)
{
  if (s.i - _NVARMAX >= nconst) {
    nconst = s.i - _NVARMAX + 1;
    _constant = (double *) prealloc (_constant, (nconst)*sizeof(double),__func__,__FILE__,__LINE__);
  }
  _constant[s.i - _NVARMAX] = val;
}

scalar new_const_scalar (const char * name, int i, double val)
{
  scalar s = (scalar){i + _NVARMAX};
  init_const_scalar (s, name, val);
  return s;
}

void init_const_vector (vector v, const char * name, double * val)
{
  {
#line 216

    init_const_scalar (v.x, name, *val++);
#line 216

    init_const_scalar (v.y, name, *val++);
#line 216

    init_const_scalar (v.z, name, *val++);}
}

vector new_const_vector (const char * name, int i, double * val)
{
  vector v;
  {
#line 223

    v.x.i = _NVARMAX + i++;
#line 223

    v.y.i = _NVARMAX + i++;
#line 223

    v.z.i = _NVARMAX + i++;}
  init_const_vector (v, name, val);
  return v;
}

void scalar_clone (scalar a, scalar b)
{
  char * name = _attribute[a.i].name;
  double (** boundary) (Point, Point, scalar, void *) = _attribute[a.i].boundary;
  double (** boundary_homogeneous) (Point, Point, scalar, void *) =
    _attribute[a.i].boundary_homogeneous;
  assert (_attribute[b.i].block > 0 && _attribute[a.i].block == _attribute[b.i].block);
  _attribute[a.i] = _attribute[b.i];
  _attribute[a.i].name = name;
  _attribute[a.i].boundary = boundary;
  _attribute[a.i].boundary_homogeneous = boundary_homogeneous;
  for (int i = 0; i < nboundary; i++) {
    _attribute[a.i].boundary[i] = _attribute[b.i].boundary[i];
    _attribute[a.i].boundary_homogeneous[i] = _attribute[b.i].boundary_homogeneous[i];
  }
}

scalar * list_clone (scalar * l)
{
  scalar * list = NULL;
  int nvar = datasize/sizeof(double), map[nvar];
  for (int i = 0; i < nvar; i++)
    map[i] = -1;
  if (l) for (scalar s = *l, *_i12 = l; ((scalar *)&s)->i >= 0; s = *++_i12) {
    scalar c = _attribute[s.i].block > 1 ? new_block_scalar("c", "", _attribute[s.i].block) :
      new_scalar("c");
    scalar_clone (c, s);
    map[s.i] = c.i;
    list = list_append (list, c);
  }
  if (list) for (scalar s = *list, *_i13 = list; ((scalar *)&s)->i >= 0; s = *++_i13)
    {
#line 260

      if (_attribute[s.i].v.x.i >= 0 && map[_attribute[s.i].v.x.i] >= 0)
 _attribute[s.i].v.x.i = map[_attribute[s.i].v.x.i];
#line 260

      if (_attribute[s.i].v.y.i >= 0 && map[_attribute[s.i].v.y.i] >= 0)
 _attribute[s.i].v.y.i = map[_attribute[s.i].v.y.i];
#line 260

      if (_attribute[s.i].v.z.i >= 0 && map[_attribute[s.i].v.z.i] >= 0)
 _attribute[s.i].v.z.i = map[_attribute[s.i].v.z.i];}
  return list;
}

void delete (scalar * list)
{
  if (all == NULL)
    return;

  if (list) for (scalar f = *list, *_i14 = list; ((scalar *)&f)->i >= 0; f = *++_i14) {
    for (int i = 0; i < _attribute[f.i].block; i++) {
      scalar fb = {f.i + i};
      if (_attribute[f.i].delete)
 _attribute[f.i].delete (fb);
      pfree (_attribute[fb.i].name,__func__,__FILE__,__LINE__); _attribute[fb.i].name = NULL;
      pfree (_attribute[fb.i].boundary,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary = NULL;
      pfree (_attribute[fb.i].boundary_homogeneous,__func__,__FILE__,__LINE__); _attribute[fb.i].boundary_homogeneous = NULL;
      _attribute[fb.i].freed = true;
    }
  }

  if (list == all) {
    all[0].i = -1;
    return;
  }

  trash (list);
  if (list) for (scalar f = *list, *_i15 = list; ((scalar *)&f)->i >= 0; f = *++_i15) {
    if (_attribute[f.i].block > 0) {
      scalar * s = all;
      for (; s->i >= 0 && s->i != f.i; s++);
      if (s->i == f.i) {
 for (; s[_attribute[f.i].block].i >= 0; s++)
   s[0] = s[_attribute[f.i].block];
 s->i = -1;
      }
    }
  }
}

void free_solver()
{
  if (free_solver_funcs) {
    free_solver_func * a = (free_solver_func *) free_solver_funcs->p;
    for (int i = 0; i < free_solver_funcs->len/sizeof(free_solver_func); i++)
      a[i] ();
    array_free (free_solver_funcs);
  }

  delete (all);
  pfree (all,__func__,__FILE__,__LINE__); all = NULL;
  for (Event * ev = Events; !ev->last; ev++) {
    Event * e = ev->next;
    while (e) {
      Event * next = e->next;
      pfree (e,__func__,__FILE__,__LINE__);
      e = next;
    }
  }

  pfree (Events,__func__,__FILE__,__LINE__); Events = NULL;
  pfree (_attribute,__func__,__FILE__,__LINE__); _attribute = NULL;
  pfree (_constant,__func__,__FILE__,__LINE__); _constant = NULL;
  free_grid();
  qpclose_all();
#if TRACE
  trace_off();
#endif
#if MTRACE
  pmuntrace();
#endif
#if _CADNA
  cadna_end();
#endif
}



void (* boundary_level) (scalar *, int l);
void (* boundary_flux) (vector *);


void boundary (scalar * list)
{ trace ("boundary", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 345);
  if (list == NULL)
    { ; end_trace("boundary", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 347);  return; }
  vector * listf = NULL;
  if (list) for (scalar s = *list, *_i16 = list; ((scalar *)&s)->i >= 0; s = *++_i16)
    if (!is_constant(s) && _attribute[s.i].block > 0 && _attribute[s.i].face)
      listf = vectors_add (listf, _attribute[s.i].v);
  if (listf) {
    boundary_flux (listf);
    pfree (listf,__func__,__FILE__,__LINE__);
  }
  boundary_level (list, -1);
 end_trace("boundary", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 357); }

void cartesian_boundary_level (scalar * list, int l)
{
  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list, l); };
}

void cartesian_boundary_flux (vector * list)
{

}

static double symmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 370 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

  return val(s,0,0,0);
}

static double antisymmetry (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 375 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

  return -val(s,0,0,0);
}

double (* default_scalar_bc[]) (Point, Point, scalar, void *) = {
  symmetry, symmetry, symmetry, symmetry, symmetry, symmetry
};

scalar cartesian_init_scalar (scalar s, const char * name)
{

  char * pname;
  if (name) {
    pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
    pname = pstrdup (name,__func__,__FILE__,__LINE__);
  }
  else
    pname = _attribute[s.i].name;
  pfree (_attribute[s.i].boundary,__func__,__FILE__,__LINE__);
  pfree (_attribute[s.i].boundary_homogeneous,__func__,__FILE__,__LINE__);

  _attribute[s.i] = (const _Attributes){0};
  _attribute[s.i].block = 1;
  _attribute[s.i].name = pname;

  _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
    pmalloc (nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  for (int b = 0; b < nboundary; b++)
    _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] =
      b < 2*3 ? default_scalar_bc[b] : symmetry;
  _attribute[s.i].gradient = NULL;
  {
#line 408
 {
    _attribute[s.i].d.x = 0;
    _attribute[s.i].v.x.i = -1;
  }
#line 408
 {
    _attribute[s.i].d.y = 0;
    _attribute[s.i].v.y.i = -1;
  }
#line 408
 {
    _attribute[s.i].d.z = 0;
    _attribute[s.i].v.z.i = -1;
  }}
  _attribute[s.i].face = false;
  return s;
}

scalar cartesian_init_vertex_scalar (scalar s, const char * name)
{
  {
#line 418

    _attribute[s.i].d.x = -1;
#line 418

    _attribute[s.i].d.y = -1;
#line 418

    _attribute[s.i].d.z = -1;}
  for (int d = 0; d < nboundary; d++)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = NULL;
  return s;
}

double (* default_vector_bc[]) (Point, Point, scalar, void *) = {
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry,
  antisymmetry, antisymmetry
};

vector cartesian_init_vector (vector v, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 434
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_scalar (v.x, cname);
    }
    else
      init_scalar (v.x, NULL);
    _attribute[v.x.i].v = v;
  }
#line 434
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_scalar (v.y, cname);
    }
    else
      init_scalar (v.y, NULL);
    _attribute[v.y.i].v = v;
  }
#line 434
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_scalar (v.z, cname);
    }
    else
      init_scalar (v.z, NULL);
    _attribute[v.z.i].v = v;
  }}

  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] =
      d < 2*3 ? default_vector_bc[d] : antisymmetry;
  return v;
}

vector cartesian_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_vector (v, name);
  {
#line 454
 {
    _attribute[v.x.i].d.x = -1;
    _attribute[v.x.i].face = true;
  }
#line 454
 {
    _attribute[v.y.i].d.y = -1;
    _attribute[v.y.i].face = true;
  }
#line 454
 {
    _attribute[v.z.i].d.z = -1;
    _attribute[v.z.i].face = true;
  }}
  for (int d = 0; d < nboundary; d++)
    _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
  return v;
}

tensor cartesian_init_tensor (tensor t, const char * name)
{
  struct { char * x, * y, * z; } ext = {".x", ".y", ".z"};
  {
#line 466
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.x);
      init_vector (t.x, cname);
    }
    else
      init_vector (t.x, NULL);
  }
#line 466
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.y);
      init_vector (t.y, cname);
    }
    else
      init_vector (t.y, NULL);
  }
#line 466
 {
    if (name) {
      char cname[strlen(name) + 3];
      sprintf (cname, "%s%s", name, ext.z);
      init_vector (t.z, cname);
    }
    else
      init_vector (t.z, NULL);
  }}
#line 490 "/home/davidtang/basilisk/src/grid/cartesian-common.h"
    assert (false);

  return t;
}

struct OutputCells {
  FILE * fp;
  coord c;
  double size;
};

void output_cells (struct OutputCells p)
{
  if (!p.fp) p.fp = fout;
   { foreach(){

#line 504 "/home/davidtang/basilisk/src/grid/cartesian-common.h"
 {
    bool inside = true;
    coord o = {x,y,z};
    {
#line 507

      if (inside && p.size > 0. &&
   (o.x > p.c.x + p.size || o.x < p.c.x - p.size))
 inside = false;
#line 507

      if (inside && p.size > 0. &&
   (o.y > p.c.y + p.size || o.y < p.c.y - p.size))
 inside = false;
#line 507

      if (inside && p.size > 0. &&
   (o.z > p.c.z + p.size || o.z < p.c.z - p.size))
 inside = false;}
    if (inside) {
      Delta /= 2.;
#line 523 "/home/davidtang/basilisk/src/grid/cartesian-common.h"
      for (int i = -1; i <= 1; i += 2) {
 fprintf (p.fp, "%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n%g %g %g\n\n",
   x - Delta, y - Delta, z + i*Delta,
   x - Delta, y + Delta, z + i*Delta,
   x + Delta, y + Delta, z + i*Delta,
   x + Delta, y - Delta, z + i*Delta,
   x - Delta, y - Delta, z + i*Delta);
 for (int j = -1; j <= 1; j += 2)
   fprintf (p.fp, "%g %g %g\n%g %g %g\n\n",
     x + i*Delta, y + j*Delta, z - Delta,
     x + i*Delta, y + j*Delta, z + Delta);
      }

    }
  } } end_foreach(); }
  fflush (p.fp);
}
#line 548 "/home/davidtang/basilisk/src/grid/cartesian-common.h"
static char * replace_ (const char * vname)
{
  char * name = pstrdup (vname,__func__,__FILE__,__LINE__), * c = name;
  while (*c != '\0') {
    if (*c == '.')
      *c = '_';
    c++;
  }
  return name;
}

static void debug_plot (FILE * fp, const char * name, const char * cells,
   const char * stencil)
{
  char * vname = replace_ (name);
  fprintf (fp,
    "  load 'debug.plot'\n"
    "  v=%s\n"







    "  splot '%s' w l lc 0, "
    "'%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 1"
           " title columnhead(4+4*v)",

    vname, cells, stencil);
  pfree (vname,__func__,__FILE__,__LINE__);
}

void cartesian_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 582 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

  char name[80] = "cells";
  if (pid() > 0)
    sprintf (name, "cells-%d", pid());
  FILE * fp = fopen (name, "w");
  output_cells ((struct OutputCells){fp, (coord){x,y,z}, 4.*Delta});
  fclose (fp);

  char stencil[80] = "stencil";
  if (pid() > 0)
    sprintf (stencil, "stencil-%d", pid());
  fp = fopen (stencil, "w");
  if (all) for (scalar v = *all, *_i17 = all; ((scalar *)&v)->i >= 0; v = *++_i17)





    fprintf (fp, "x y z %s ", _attribute[v.i].name);

  fputc ('\n', fp);
#line 629 "/home/davidtang/basilisk/src/grid/cartesian-common.h"
    for (int k = -2; k <= 2; k++)
      for (int l = -2; l <= 2; l++)
 for (int m = -2; m <= 2; m++) {
   if (all) for (scalar v = *all, *_i18 = all; ((scalar *)&v)->i >= 0; v = *++_i18) {
     fprintf (fp, "%g %g %g ",
       x + k*Delta + _attribute[v.i].d.x*Delta/2.,
       y + l*Delta + _attribute[v.i].d.y*Delta/2.,
       z + m*Delta + _attribute[v.i].d.z*Delta/2.);
     if (allocated(k,l,m))
       fprintf (fp, "%g ", val(v,k,l,m));
     else
       fputs ("n/a ", fp);
   }
   fputc ('\n', fp);
 }

  fclose (fp);

  fp = fopen ("debug.plot", "w");
  fprintf (fp,
    "set term x11\n"
    "set size ratio -1\n"
    "set key outside\n");
  if (all) for (scalar s = *all, *_i19 = all; ((scalar *)&s)->i >= 0; s = *++_i19) {
    char * name = replace_ (_attribute[s.i].name);
    fprintf (fp, "%s = %d\n", name, s.i);
    pfree (name,__func__,__FILE__,__LINE__);
  }
  fclose (fp);

  fprintf (ferr, "Last point stencils can be displayed using (in gnuplot)\n");
  debug_plot (ferr, _attribute[0].name, name, stencil);
  fflush (ferr);

  fp = fopen ("plot", "w");
  debug_plot (fp, _attribute[0].name, name, stencil);
  fclose (fp);
}

void cartesian_methods()
{
  init_scalar = cartesian_init_scalar;
  init_vertex_scalar = cartesian_init_vertex_scalar;
  init_vector = cartesian_init_vector;
  init_tensor = cartesian_init_tensor;
  init_face_vector = cartesian_init_face_vector;
  boundary_level = cartesian_boundary_level;
  boundary_flux = cartesian_boundary_flux;
  debug = cartesian_debug;
}

struct _interpolate {
  scalar v;
  double x, y, z;
};

static double interpolate_linear (Point point, struct _interpolate p)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 686 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

  scalar v = p.v;
#line 703 "/home/davidtang/basilisk/src/grid/cartesian-common.h"
  x = (p.x - x)/Delta - _attribute[v.i].d.x/2.;
  y = (p.y - y)/Delta - _attribute[v.i].d.y/2.;
  z = (p.z - z)/Delta - _attribute[v.i].d.z/2.;
  int i = sign(x), j = sign(y), k = sign(z);
  x = fabs(x); y = fabs(y); z = fabs(z);

  return (((val(v,0,0,0)*(1. - x) + val(v,i,0,0)*x)*(1. - y) +
    (val(v,0,j,0)*(1. - x) + val(v,i,j,0)*x)*y)*(1. - z) +
   ((val(v,0,0,k)*(1. - x) + val(v,i,0,k)*x)*(1. - y) +
    (val(v,0,j,k)*(1. - x) + val(v,i,j,k)*x)*y)*z);

}


double interpolate (struct _interpolate p)
{ trace ("interpolate", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 718);
  Point point = locate ((struct _locate){p.x, p.y, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 719 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

  if (point.level < 0)
    { double _ret =  nodata; end_trace("interpolate", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 721);  return _ret; }
  { double _ret =  interpolate_linear (point, p); end_trace("interpolate", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 722);  return _ret; }
 end_trace("interpolate", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 723); }


void interpolate_array (scalar * list, coord * a, int n, double * v, bool linear)
{ trace ("interpolate_array", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 727);
  int j = 0;
  for (int i = 0; i < n; i++) {
    Point point = locate ((struct _locate){a[i].x, a[i].y, a[i].z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 730 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

    if (point.level >= 0) {
      if (list) for (scalar s = *list, *_i20 = list; ((scalar *)&s)->i >= 0; s = *++_i20)
 v[j++] = !linear ? val(s,0,0,0) :
   interpolate_linear (point,
         (struct _interpolate){s, a[i].x, a[i].y, a[i].z});
    }
    else
      if (list) for (scalar s = *list, *_i21 = list; ((scalar *)&s)->i >= 0; s = *++_i21)
 v[j++] = nodata;
  }
#if 1
  if (pid() == 0)
    MPI_Reduce (MPI_IN_PLACE, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
  else
    MPI_Reduce (v, v, n*list_len(list), MPI_DOUBLE,
  MPI_MIN, 0, MPI_COMM_WORLD);
#endif
 end_trace("interpolate_array", "/home/davidtang/basilisk/src/grid/cartesian-common.h", 749); }



typedef int bid;

bid new_bid()
{
  int b = nboundary++;
  if (all) for (scalar s = *all, *_i22 = all; ((scalar *)&s)->i >= 0; s = *++_i22) {
    _attribute[s.i].boundary = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
    _attribute[s.i].boundary_homogeneous = (double (**)(Point, Point, scalar, void *))
      prealloc (_attribute[s.i].boundary_homogeneous, nboundary*sizeof (void (*)()),__func__,__FILE__,__LINE__);
  }
  if (all) for (scalar s = *all, *_i23 = all; ((scalar *)&s)->i >= 0; s = *++_i23) {
    if (_attribute[s.i].v.x.i < 0)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b] = symmetry;
    else if (_attribute[s.i].v.x.i == s.i) {
      vector v = _attribute[s.i].v;
      {
#line 769

 _attribute[v.y.i].boundary[b] = _attribute[v.y.i].boundary_homogeneous[b] = symmetry;
#line 769

 _attribute[v.z.i].boundary[b] = _attribute[v.z.i].boundary_homogeneous[b] = symmetry;
#line 769

 _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] = symmetry;}
      _attribute[v.x.i].boundary[b] = _attribute[v.x.i].boundary_homogeneous[b] =
 _attribute[v.x.i].face ? NULL : antisymmetry;
    }
  }
  return b;
}



static double periodic_bc (Point point, Point neighbor, scalar s, void * data)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 781 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

  return nodata;
}

static void periodic_boundary (int d)
{

  if (all) for (scalar s = *all, *_i24 = all; ((scalar *)&s)->i >= 0; s = *++_i24)
    _attribute[s.i].boundary[d] = _attribute[s.i].boundary_homogeneous[d] = periodic_bc;

  if (all) for (scalar s = *all, *_i25 = all; ((scalar *)&s)->i >= 0; s = *++_i25)
    if (_attribute[s.i].face) {
      vector v = _attribute[s.i].v;
      _attribute[v.x.i].boundary[d] = _attribute[v.x.i].boundary_homogeneous[d] = NULL;
    }

  default_scalar_bc[d] = periodic_bc;
  default_vector_bc[d] = periodic_bc;
}

void periodic (int dir)
{





    assert (dir <= back);


  int c = dir/2;
  periodic_boundary (2*c);
  periodic_boundary (2*c + 1);
  (&Period.x)[c] = true;
}


double getvalue (Point point, scalar s, int i, int j, int k)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 819 "/home/davidtang/basilisk/src/grid/cartesian-common.h"

  return val(s,i,j,k);
}
#line 4 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

#ifndef foreach_level_or_leaf
# define foreach_level_or_leaf foreach_level
# define end_foreach_level_or_leaf end_foreach_level
#endif

#ifndef foreach_coarse_level
# define foreach_coarse_level foreach_level
# define end_foreach_coarse_level end_foreach_level
#endif










void (* restriction) (scalar *);

static inline void restriction_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 27 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

  double sum = 0.;
   { foreach_child()
    sum += val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3);
}

static inline void restriction_volume_average (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 35 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/(val_cm(cm,0,0,0) + 1e-30);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 35

  double sum = 0.;
   { foreach_child()
    sum += val_cm(cm,0,0,0)*val(s,0,0,0); end_foreach_child(); }
  val(s,0,0,0) = sum/(1 << 3)/(val_cm(cm,0,0,0) + 1e-30);
 }}

static inline void face_average (Point point, vector v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 43 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

  {
#line 44
 {







      val(v.x,0,0,0) = (fine(v.x,0,0,0) + fine(v.x,0,1,0) +
        fine(v.x,0,0,1) + fine(v.x,0,1,1))/4.;
      val(v.x,1,0,0) = (fine(v.x,2,0,0) + fine(v.x,2,1,0) +
  fine(v.x,2,0,1) + fine(v.x,2,1,1))/4.;

  }
#line 44
 {







      val(v.y,0,0,0) = (fine(v.y,0,0,0) + fine(v.y,0,0,1) +
        fine(v.y,1,0,0) + fine(v.y,1,0,1))/4.;
      val(v.y,0,1,0) = (fine(v.y,0,2,0) + fine(v.y,0,2,1) +
  fine(v.y,1,2,0) + fine(v.y,1,2,1))/4.;

  }
#line 44
 {







      val(v.z,0,0,0) = (fine(v.z,0,0,0) + fine(v.z,1,0,0) +
        fine(v.z,0,1,0) + fine(v.z,1,1,0))/4.;
      val(v.z,0,0,1) = (fine(v.z,0,0,2) + fine(v.z,1,0,2) +
  fine(v.z,0,1,2) + fine(v.z,1,1,2))/4.;

  }}
}

static inline void restriction_face (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 61 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

  face_average (point, _attribute[s.i].v);
}

static inline void restriction_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 66 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

  for (int i = 0; i <= 1; i++) {
    val(s,i,0,0) = fine(s,2*i,0,0);

    val(s,i,1,0) = fine(s,2*i,2,0);


    for (int j = 0; j <= 1; j++)
      val(s,i,j,1) = fine(s,2*i,2*j,2);

  }
}

static inline void no_restriction (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 78 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
}

static inline void no_data (Point point, scalar s) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 81 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = nodata; end_foreach_child(); }
}

void wavelet (scalar s, scalar w)
{
  restriction (((scalar []){s,{-1}}));
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level (l){

#line 90 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
 {
       { foreach_child()
        val(w,0,0,0) = val(s,0,0,0); end_foreach_child(); }
      _attribute[s.i].prolongation (point, s);
       { foreach_child() {
        double sp = val(s,0,0,0);
        val(s,0,0,0) = val(w,0,0,0);

        val(w,0,0,0) -= sp;
      } end_foreach_child(); }
    } } end_foreach_coarse_level(); }
    boundary_level (((scalar []){w,{-1}}), l + 1);
  }

   { foreach_level(0){

#line 104 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

    val(w,0,0,0) = val(s,0,0,0); } end_foreach_level(); }
  boundary_level (((scalar []){w,{-1}}), 0);
}

void inverse_wavelet (scalar s, scalar w)
{
   { foreach_level(0){

#line 111 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

    val(s,0,0,0) = val(w,0,0,0); } end_foreach_level(); }
  boundary_level (((scalar []){s,{-1}}), 0);
  for (int l = 0; l <= depth() - 1; l++) {
     { foreach_coarse_level (l){

#line 115 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
 {
      _attribute[s.i].prolongation (point, s);
       { foreach_child()
        val(s,0,0,0) += val(w,0,0,0); end_foreach_child(); }
    } } end_foreach_coarse_level(); }
    boundary_level (((scalar []){s,{-1}}), l + 1);
  }
}

static inline double bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 125 "/home/davidtang/basilisk/src/grid/multigrid-common.h"








    return (27.*coarse(s,0,0,0) +
     9.*(coarse(s,child.x,0,0) + coarse(s,0,child.y,0) +
  coarse(s,0,0,child.z)) +
     3.*(coarse(s,child.x,child.y,0) + coarse(s,child.x,0,child.z) +
  coarse(s,0,child.y,child.z)) +
     coarse(s,child.x,child.y,child.z))/64.;

}

static inline void refine_bilinear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 143 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = bilinear (point, s); end_foreach_child(); }
}

static inline double quadratic (double a, double b, double c)
{
  return (30.*a + 5.*b - 3.*c)/32.;
}

static inline double biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 154 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

#line 169 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
  assert (false);
  return 0.;

}

static inline double biquadratic_vertex (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 175 "/home/davidtang/basilisk/src/grid/multigrid-common.h"







  assert (false);
  return 0.;

}

static inline void refine_biquadratic (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 188 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(s,0,0,0) = biquadratic (point, s); end_foreach_child(); }
}

static inline void refine_linear (Point point, scalar s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 194 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 194

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 197

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 200

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 200

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 200

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 206

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 206

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 194

  coord g;
  if (_attribute[s.i].gradient)
    {
#line 197

      g.x = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0));
#line 197

      g.y = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0));
#line 197

      g.z = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1));}
  else
    {
#line 200

      g.x = (val(s,1,0,0) - val(s,-1,0,0))/2.;
#line 200

      g.y = (val(s,0,1,0) - val(s,0,-1,0))/2.;
#line 200

      g.z = (val(s,0,0,1) - val(s,0,0,-1))/2.;}

  double sc = val(s,0,0,0), cmc = 4.*val_cm(cm,0,0,0), sum = val_cm(cm,0,0,0)*(1 << 3);
   { foreach_child() {
    val(s,0,0,0) = sc;
    {
#line 206

      val(s,0,0,0) += child.x*g.x*val_cm(cm,-child.x,0,0)/cmc;
#line 206

      val(s,0,0,0) += child.y*g.y*val_cm(cm,0,-child.y,0)/cmc;
#line 206

      val(s,0,0,0) += child.z*g.z*val_cm(cm,0,0,-child.z)/cmc;}
    sum -= val_cm(cm,0,0,0);
  } end_foreach_child(); }
  assert (fabs(sum) < 1e-10);
 }}

static inline void refine_reset (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 214 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

   { foreach_child()
    val(v,0,0,0) = 0.; end_foreach_child(); }
}

static inline void refine_injection (Point point, scalar v)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 220 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

  double val = val(v,0,0,0);
   { foreach_child()
    val(v,0,0,0) = val; end_foreach_child(); }
}

static scalar multigrid_init_scalar (scalar s, const char * name)
{
  s = cartesian_init_scalar (s, name);
  _attribute[s.i].prolongation = refine_bilinear;
  _attribute[s.i].restriction = restriction_average;
  return s;
}

static scalar multigrid_init_vertex_scalar (scalar s, const char * name)
{
  s = cartesian_init_vertex_scalar (s, name);
  _attribute[s.i].restriction = restriction_vertex;
  return s;
}

static vector multigrid_init_face_vector (vector v, const char * name)
{
  v = cartesian_init_face_vector (v, name);
  {
#line 244

    _attribute[v.y.i].restriction = no_restriction;
#line 244

    _attribute[v.z.i].restriction = no_restriction;
#line 244

    _attribute[v.x.i].restriction = no_restriction;}
  _attribute[v.x.i].restriction = restriction_face;
  return v;
}

void multigrid_debug (Point point)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 251 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

  cartesian_debug (point);

  FILE * plot = fopen ("plot", "a");
  if (point.level > 0) {
    char name[80] = "coarse";
    if (pid() > 0)
      sprintf (name, "coarse-%d", pid());
    FILE * fp = fopen (name, "w");
#line 284 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
      double xc = x - child.x*Delta/2., yc = y - child.y*Delta/2.;
      double zc = z - child.z*Delta/2.;
      for (int k = 0; k <= 1; k++)
 for (int l = 0; l <= 1; l++)
   for (int m = 0; m <= 1; m++) {
     if (all) for (scalar v = *all, *_i26 = all; ((scalar *)&v)->i >= 0; v = *++_i26)
       fprintf (fp, "%g %g %g %g ",
         xc + k*child.x*Delta*2. + _attribute[v.i].d.x*Delta,
         yc + l*child.y*Delta*2. + _attribute[v.i].d.y*Delta,
         zc + m*child.z*Delta*2. + _attribute[v.i].d.z*Delta,
         coarse(v,k*child.x,l*child.y,m*child.z));
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 3 t ''",
        name);

    fclose (fp);
  }

  if (is_coarse()) {
    char name[80] = "fine";
    if (pid() > 0)
      sprintf (name, "fine-%d", pid());
    FILE * fp = fopen (name, "w");
#line 341 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
      double xf = x - Delta/4., yf = y - Delta/4., zf = z - Delta/4.;
      for (int k = -2; k <= 3; k++)
 for (int l = -2; l <= 3; l++)
   for (int m = -2; m <= 3; m++) {
     if (all) for (scalar v = *all, *_i27 = all; ((scalar *)&v)->i >= 0; v = *++_i27) {
       fprintf (fp, "%g %g %g ",
         xf + k*Delta/2. + _attribute[v.i].d.x*Delta/4.,
         yf + l*Delta/2. + _attribute[v.i].d.y*Delta/4.,
         zf + m*Delta/2. + _attribute[v.i].d.z*Delta/4.);
       if (allocated_child(k,l,m))
  fprintf (fp, "%g ", fine(v,k,l,m));
       else
  fputs ("n/a ", fp);
     }
     fputc ('\n', fp);
   }
      fprintf (ferr, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);
      fprintf (plot, ", '%s' u 1+4*v:2+4*v:3+4*v:4+4*v w labels tc lt 2 t ''",
        name);

    fclose (fp);
  }
  fflush (ferr);
  fclose (plot);
}

static void multigrid_restriction (scalar * list)
{
  scalar * listdef = NULL, * listc = NULL, * list2 = NULL;
  if (list) for (scalar s = *list, *_i28 = list; ((scalar *)&s)->i >= 0; s = *++_i28)
    if (!is_constant (s) && _attribute[s.i].block > 0) {
      if (_attribute[s.i].restriction == restriction_average) {
 listdef = list_add (listdef, s);
 list2 = list_add (list2, s);
      }
      else if (_attribute[s.i].restriction != no_restriction) {
 listc = list_add (listc, s);
 if (_attribute[s.i].face)
   {
#line 380

     list2 = list_add (list2, _attribute[s.i].v.x);
#line 380

     list2 = list_add (list2, _attribute[s.i].v.y);
#line 380

     list2 = list_add (list2, _attribute[s.i].v.z);}
 else
   list2 = list_add (list2, s);
      }
    }

  if (listdef || listc) {
    for (int l = depth() - 1; l >= 0; l--) {
       { foreach_coarse_level(l){

#line 389 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
 {
 if (listdef) for (scalar s = *listdef, *_i29 = listdef; ((scalar *)&s)->i >= 0; s = *++_i29)
  
     restriction_average (point, s);
 if (listc) for (scalar s = *listc, *_i30 = listc; ((scalar *)&s)->i >= 0; s = *++_i30) {
  
     _attribute[s.i].restriction (point, s);
 }
      } } end_foreach_coarse_level(); }
      { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->level) _b->level (_b, list2, l); };
    }
    pfree (listdef,__func__,__FILE__,__LINE__);
    pfree (listc,__func__,__FILE__,__LINE__);
    pfree (list2,__func__,__FILE__,__LINE__);
  }
}

void multigrid_methods()
{
  cartesian_methods();
  debug = multigrid_debug;
  init_scalar = multigrid_init_scalar;
  init_vertex_scalar = multigrid_init_vertex_scalar;
  init_face_vector = multigrid_init_face_vector;
  restriction = multigrid_restriction;
}







void subtree_size (scalar size, bool leaves)
{




   { foreach(){

#line 428 "/home/davidtang/basilisk/src/grid/multigrid-common.h"

    val(size,0,0,0) = 1; } end_foreach(); }





  { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), depth()); };
  for (int l = depth() - 1; l >= 0; l--) {
     { foreach_coarse_level(l){

#line 437 "/home/davidtang/basilisk/src/grid/multigrid-common.h"
 {
      double sum = !leaves;
       { foreach_child()
 sum += val(size,0,0,0); end_foreach_child(); }
      val(size,0,0,0) = sum;
    } } end_foreach_coarse_level(); }
    { Boundary ** _i = boundaries, * _b; while (_i && (_b = *_i++)) if (_b->restriction) _b->restriction (_b, ((scalar []){size,{-1}}), l); };
  }
}
#line 826 "/home/davidtang/basilisk/src/grid/multigrid.h"

struct Dimensions {
  int nx, ny, nz;
};

void dimensions (struct Dimensions p)
{

  for (int i = 0; i < 3; i++)
    mpi_dims[i] = (&p.nx)[i];

}



#if 3 == 1

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\

#line 847

#define end_foreach_slice_x() }

#elif 3 == 2

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\

#line 857

#define end_foreach_slice_x() }

#define foreach_slice_y(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = start; point.j < end; point.j++)\

#line 865

#define end_foreach_slice_y() }

#elif 3 == 3

#define foreach_slice_x(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = start; point.i < end; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\
      for (point.k = 0; point.k < point.n + 2*2; point.k++)\

#line 876

#define end_foreach_slice_x() }

#define foreach_slice_y(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = start; point.j < end; point.j++)\
      for (point.k = 0; point.k < point.n + 2*2; point.k++)\

#line 885

#define end_foreach_slice_y() }

#define foreach_slice_z(start, end, l) {\
  Point point = {0};\
  point.level = l; point.n = 1 << point.level;\
  for (point.i = 0; point.i < point.n + 2*2; point.i++)\
    for (point.j = 0; point.j < point.n + 2*2; point.j++)\
      for (point.k = start; point.k < end; point.k++)\

#line 894

#define end_foreach_slice_z() }

#endif

#line 1 "grid/multigrid-mpi.h"
#line 1 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"


typedef struct {
  Boundary b;
  MPI_Comm cartcomm;
} MpiBoundary;

#define BUF void *


#line 10

static BUF snd_x (int i, int dst, int tag, int level, scalar * list,
    MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  if (list) for (scalar s = *list, *_i31 = list; ((scalar *)&s)->i >= 0; s = *++_i31)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   { foreach_slice_x (i, i + 2, level){

#line 21 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"

    if (list) for (scalar s = *list, *_i32 = list; ((scalar *)&s)->i >= 0; s = *++_i32) {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_x(); }
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}
#line 10

static BUF snd_y (int i, int dst, int tag, int level, scalar * list,
    MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  if (list) for (scalar s = *list, *_i31 = list; ((scalar *)&s)->i >= 0; s = *++_i31)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   { foreach_slice_y (i, i + 2, level){

#line 21 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"

    if (list) for (scalar s = *list, *_i32 = list; ((scalar *)&s)->i >= 0; s = *++_i32) {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_y(); }
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}
#line 10

static BUF snd_z (int i, int dst, int tag, int level, scalar * list,
    MPI_Request * req)
{
  if (dst == MPI_PROC_NULL)
    return NULL;
  size_t size = 0;
  if (list) for (scalar s = *list, *_i31 = list; ((scalar *)&s)->i >= 0; s = *++_i31)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
   { foreach_slice_z (i, i + 2, level){

#line 21 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"

    if (list) for (scalar s = *list, *_i32 = list; ((scalar *)&s)->i >= 0; s = *++_i32) {
      memcpy (b, &val(s,0,0,0), sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_z(); }
  MPI_Isend (buf, size, MPI_BYTE, dst, tag, MPI_COMM_WORLD, req);
  return buf;
}


#line 30

static void rcv_x (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  if (list) for (scalar s = *list, *_i33 = list; ((scalar *)&s)->i >= 0; s = *++_i33)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   { foreach_slice_x (i, i + 2, level){

#line 42 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"

    if (list) for (scalar s = *list, *_i34 = list; ((scalar *)&s)->i >= 0; s = *++_i34) {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_x(); }
  pfree (buf,__func__,__FILE__,__LINE__);
}
#line 30

static void rcv_y (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  if (list) for (scalar s = *list, *_i33 = list; ((scalar *)&s)->i >= 0; s = *++_i33)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   { foreach_slice_y (i, i + 2, level){

#line 42 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"

    if (list) for (scalar s = *list, *_i34 = list; ((scalar *)&s)->i >= 0; s = *++_i34) {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_y(); }
  pfree (buf,__func__,__FILE__,__LINE__);
}
#line 30

static void rcv_z (int i, int src, int tag, int level, scalar * list)
{
  if (src == MPI_PROC_NULL)
    return;
  size_t size = 0;
  if (list) for (scalar s = *list, *_i33 = list; ((scalar *)&s)->i >= 0; s = *++_i33)
    size += _attribute[s.i].block;
  size *= pow((1 << level) + 2*2, 3 - 1)*2*sizeof(double);
  double * buf = (double *) pmalloc (size,__func__,__FILE__,__LINE__), * b = buf;
  MPI_Status s;
  MPI_Recv (buf, size, MPI_BYTE, src, tag, MPI_COMM_WORLD, &s);
   { foreach_slice_z (i, i + 2, level){

#line 42 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"

    if (list) for (scalar s = *list, *_i34 = list; ((scalar *)&s)->i >= 0; s = *++_i34) {
      memcpy (&val(s,0,0,0), b, sizeof(double)*_attribute[s.i].block);
      b += _attribute[s.i].block;
    } } end_foreach_slice_z(); }
  pfree (buf,__func__,__FILE__,__LINE__);
}


static void mpi_boundary_level (const Boundary * b, scalar * list, int level)
{ trace ("mpi_boundary_level", "/home/davidtang/basilisk/src/grid/multigrid-mpi.h", 52);
  scalar * list1 = NULL;
  if (list) for (scalar s = *list, *_i35 = list; ((scalar *)&s)->i >= 0; s = *++_i35)
    if (!is_constant(s) && _attribute[s.i].block > 0)
      list1 = list_add (list1, s);
  if (!list1)
    { ; end_trace("mpi_boundary_level", "/home/davidtang/basilisk/src/grid/multigrid-mpi.h", 58);  return; }

  prof_start ("mpi_boundary_level");

  if (level < 0) level = depth();
  MpiBoundary * mpi = (MpiBoundary *) b;
  struct { int x, y, z; } dir = {0,1,2};
  {
#line 65
 {
    int left, right;
    MPI_Cart_shift (mpi->cartcomm, dir.x, 1, &left, &right);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_x (npl - 2*2, right, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_x (2, left, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_x (0, left, 0, level, list1);
    rcv_x (npl - 2, right, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  }
#line 65
 {
    int bottom, top;
    MPI_Cart_shift (mpi->cartcomm, dir.y, 1, &bottom, &top);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_y (npl - 2*2, top, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_y (2, bottom, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_y (0, bottom, 0, level, list1);
    rcv_y (npl - 2, top, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  }
#line 65
 {
    int back, front;
    MPI_Cart_shift (mpi->cartcomm, dir.z, 1, &back, &front);
    MPI_Request reqs[2];
    void * buf[2];
    int npl = (1 << level) + 2*2, nr = 0;
    if ((buf[0] = snd_z (npl - 2*2, front, 0, level, list1, &reqs[nr])))
      nr++;
    if ((buf[1] = snd_z (2, back, 1, level, list1, &reqs[nr])))
      nr++;
    rcv_z (0, back, 0, level, list1);
    rcv_z (npl - 2, front, 1, level, list1);
    MPI_Status stats[nr];
    MPI_Waitall (nr, reqs, stats);
    pfree (buf[0],__func__,__FILE__,__LINE__); pfree (buf[1],__func__,__FILE__,__LINE__);
  }}

  pfree (list1,__func__,__FILE__,__LINE__);

  prof_stop();
 end_trace("mpi_boundary_level", "/home/davidtang/basilisk/src/grid/multigrid-mpi.h", 85); }

static void mpi_boundary_destroy (Boundary * b)
{
  MpiBoundary * m = (MpiBoundary *) b;
  MPI_Comm_free (&m->cartcomm);
  pfree (m,__func__,__FILE__,__LINE__);
}

Boundary * mpi_boundary_new()
{
  MpiBoundary * m = ((MpiBoundary *) pcalloc (1, sizeof(MpiBoundary),__func__,__FILE__,__LINE__));
  MPI_Dims_create (npe(), 3, mpi_dims);
  MPI_Cart_create (MPI_COMM_WORLD, 3,
     mpi_dims, &Period.x, 0, &m->cartcomm);
  MPI_Cart_coords (m->cartcomm, pid(), 3, mpi_coords);


  struct { int x, y, z; } dir = {0,1,2};
  {
#line 104
 {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.x, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (left);
    if (r != MPI_PROC_NULL)
      periodic_boundary (right);
  }
#line 104
 {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.y, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (bottom);
    if (r != MPI_PROC_NULL)
      periodic_boundary (top);
  }
#line 104
 {
    int l, r;
    MPI_Cart_shift (m->cartcomm, dir.z, 1, &l, &r);
    if (l != MPI_PROC_NULL)
      periodic_boundary (back);
    if (r != MPI_PROC_NULL)
      periodic_boundary (front);
  }}


  N /= mpi_dims[0];
  int r = 0;
  while (N > 1)
    N /= 2, r++;
  grid->depth = grid->maxdepth = r;
  N = mpi_dims[0]*(1 << r);
  grid->n = 1 << 3*depth();
  grid->tn = npe()*grid->n;


  Boundary * b = (Boundary *) m;
  b->level = mpi_boundary_level;
  b->destroy = mpi_boundary_destroy;
  add_boundary (b);

  return b;
}


double z_indexing (scalar index, bool leaves)
{ trace ("z_indexing", "/home/davidtang/basilisk/src/grid/multigrid-mpi.h", 134);
  long i;
  if (leaves)
    i = pid()*(1 << 3*depth());
  else
    i = pid()*((1 << 3*(depth() + 1)) - 1)/((1 << 3) - 1);
   { foreach_cell(){

#line 140 "/home/davidtang/basilisk/src/grid/multigrid-mpi.h"
 {
    if (!leaves || is_leaf(cell))
      val(index,0,0,0) = i++;
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary (((scalar []){index,{-1}}));
  { double _ret =  pid() == 0 ? i*npe() - 1 : -1; end_trace("z_indexing", "/home/davidtang/basilisk/src/grid/multigrid-mpi.h", 147);  return _ret; }
 end_trace("z_indexing", "/home/davidtang/basilisk/src/grid/multigrid-mpi.h", 148); }
#line 900 "/home/davidtang/basilisk/src/grid/multigrid.h"
#line 3 "/home/davidtang/basilisk/src/grid/multigrid3D.h"
#line 14 "isotropic.c"
#line 1 "navier-stokes/centered.h"
#line 1 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
#line 27 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
#line 1 "./run.h"
#line 1 "/home/davidtang/basilisk/src/run.h"
#line 9 "/home/davidtang/basilisk/src/run.h"
double dt = 1.;

#line 1 "./utils.h"
#line 1 "/home/davidtang/basilisk/src/utils.h"







double DT = 1e10, CFL = 0.5;




struct {

  long nc;

  long tnc;

  double t;

  double speed;

  timer gt;
} perf;





void update_perf() {
  perf.nc += grid->n;
  perf.tnc += grid->tn;
  perf.t = timer_elapsed (perf.gt);
  perf.speed = perf.tnc/perf.t;
}






typedef struct {
  double cpu;
  double real;
  double speed;
  double min;
  double avg;
  double max;
  size_t tnc;
  long mem;
} timing;






timing timer_timing (timer t, int i, size_t tnc, double * mpi)
{
  timing s;
#if 1
  s.avg = mpi_time - t.tm;
#endif
  clock_t end = clock();
  s.cpu = ((double) (end - t.c))/CLOCKS_PER_SEC;
  s.real = timer_elapsed (t);
  if (tnc == 0) {
    double n = 0;
     { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _n = n; 
#line 69
foreach(){

#line 69 "/home/davidtang/basilisk/src/utils.h"
 _n++; } end_foreach();OMP(omp critical) n += _n;
mpi_all_reduce_double (n, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 69
 }
    s.tnc = n;
    tnc = n*i;
  }
  else
    s.tnc = tnc;
#if _GNU_SOURCE
  struct rusage usage;
  getrusage (RUSAGE_SELF, &usage);
  s.mem = usage.ru_maxrss;
#else
  s.mem = 0;
#endif
#if 1
  if (mpi)
    MPI_Allgather (&s.avg, 1, MPI_DOUBLE, mpi, 1, MPI_DOUBLE, MPI_COMM_WORLD);
  s.max = s.min = s.avg;
  mpi_all_reduce (s.max, MPI_DOUBLE, MPI_MAX);
  mpi_all_reduce (s.min, MPI_DOUBLE, MPI_MIN);
  mpi_all_reduce (s.avg, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.real, MPI_DOUBLE, MPI_SUM);
  mpi_all_reduce (s.mem, MPI_LONG, MPI_SUM);
  s.real /= npe();
  s.avg /= npe();
  s.mem /= npe();
#else
  s.min = s.max = s.avg = 0.;
#endif
  s.speed = s.real > 0. ? tnc/s.real : -1.;
  return s;
}




void timer_print (timer t, int i, size_t tnc)
{
  timing s = timer_timing (t, i, tnc, NULL);
  fprintf (fout,
    "\n# " "Multigrid"
    ", %d steps, %g CPU, %.4g real, %.3g points.step/s, %d var\n",
    i, s.cpu, s.real, s.speed, (int) (datasize/sizeof(double)));
#if 1
  fprintf (fout,
    "# %d procs, MPI: min %.2g (%.2g%%) "
    "avg %.2g (%.2g%%) max %.2g (%.2g%%)\n",
    npe(),
    s.min, 100.*s.min/s.real,
    s.avg, 100.*s.avg/s.real,
    s.max, 100.*s.max/s.real);
#endif
}







typedef struct {
  double avg, rms, max, volume;
} norm;

norm normf (scalar f)
{
  double avg = 0., rms = 0., max = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; double _avg = avg; double _rms = rms; double _volume = volume; 
#line 135

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 135
foreach(){

#line 136 "/home/davidtang/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata && (cube(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 135
foreach(){

#line 136 "/home/davidtang/basilisk/src/utils.h"

    if (val(f,0,0,0) != nodata && (cube(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double v = fabs(val(f,0,0,0));
      if (v > _max) _max = v;
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _avg += (cube(Delta)*val_cm(cm,0,0,0))*v;
      _rms += (cube(Delta)*val_cm(cm,0,0,0))*sq(v);
    } } end_foreach(); }OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) avg += _avg;
mpi_all_reduce_double (avg, MPI_SUM);
OMP(omp critical) rms += _rms;
mpi_all_reduce_double (rms, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 143
 }
  norm n;
  n.avg = volume ? avg/volume : 0.;
  n.rms = volume ? sqrt(rms/volume) : 0.;
  n.max = max;
  n.volume = volume;
  return n;
}





typedef struct {
  double min, max, sum, stddev, volume;
} stats;

stats statsf (scalar f)
{
  double min = 1e100, max = -1e100, sum = 0., sum2 = 0., volume = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; double _sum2 = sum2; double _volume = volume; double _max = max; double _min = min; 
#line 163

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 163
foreach(){

#line 164 "/home/davidtang/basilisk/src/utils.h"

    if ((cube(Delta)*val_cm(cm,0,0,0)) > 0. && val(f,0,0,0) != nodata) {
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 163
foreach(){

#line 164 "/home/davidtang/basilisk/src/utils.h"

    if ((cube(Delta)*val_cm(cm,0,0,0)) > 0. && val(f,0,0,0) != nodata) {
      _volume += (cube(Delta)*val_cm(cm,0,0,0));
      _sum += (cube(Delta)*val_cm(cm,0,0,0))*val(f,0,0,0);
      _sum2 += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(f,0,0,0));
      if (val(f,0,0,0) > _max) _max = val(f,0,0,0);
      if (val(f,0,0,0) < _min) _min = val(f,0,0,0);
    } } end_foreach(); }OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);
OMP(omp critical) sum2 += _sum2;
mpi_all_reduce_double (sum2, MPI_SUM);
OMP(omp critical) volume += _volume;
mpi_all_reduce_double (volume, MPI_SUM);
OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);
OMP(omp critical) if (_min < min) min = _min;
mpi_all_reduce_double (min, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 171
 }
  stats s;
  s.min = min, s.max = max, s.sum = sum, s.volume = volume;
  if (volume > 0.)
    sum2 -= sum*sum/volume;
  s.stddev = sum2 > 0. ? sqrt(sum2/volume) : 0.;
  return s;
}
#line 187 "/home/davidtang/basilisk/src/utils.h"
static double generic_limiter (double r, double beta)
{
  double v1 = min (r, beta), v2 = min (beta*r, 1.);
  v1 = max (0., v1);
  return max (v1, v2);
}

double minmod (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.)*(s1 - s0);
}

double superbee (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 2.)*(s1 - s0);
}

double sweby (double s0, double s1, double s2) {
  return s1 == s0 ? 0. : generic_limiter ((s2 - s1)/(s1 - s0), 1.5)*(s1 - s0);
}
#line 213 "/home/davidtang/basilisk/src/utils.h"
double theta = 1.3;

double minmod2 (double s0, double s1, double s2)
{
  if (s0 < s1 && s1 < s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 < d1) d1 = d2;
    return min(d1, d3);
  }
  if (s0 > s1 && s1 > s2) {
    double d1 = theta*(s1 - s0), d2 = (s2 - s0)/2., d3 = theta*(s2 - s1);
    if (d2 > d1) d1 = d2;
    return max(d1, d3);
  }
  return 0.;
}
#line 237 "/home/davidtang/basilisk/src/utils.h"
void gradients (scalar * f, vector * g)
{
  assert (list_len(f) == vectors_len(g));
   { foreach(){

#line 240 "/home/davidtang/basilisk/src/utils.h"
 {
    scalar s; vector v;
    scalar * _i2 = f; vector * _i3 = g; if (f) for (s = *f, v = *g; ((scalar *)&s)->i >= 0; s = *++_i2, v = *++_i3) {
      if (_attribute[s.i].gradient)
 {
#line 244
 {





     val(v.x,0,0,0) = _attribute[s.i].gradient (val(s,-1,0,0), val(s,0,0,0), val(s,1,0,0))/Delta;
 }
#line 244
 {





     val(v.y,0,0,0) = _attribute[s.i].gradient (val(s,0,-1,0), val(s,0,0,0), val(s,0,1,0))/Delta;
 }
#line 244
 {





     val(v.z,0,0,0) = _attribute[s.i].gradient (val(s,0,0,-1), val(s,0,0,0), val(s,0,0,1))/Delta;
 }}
      else
 {
#line 253
 {





     val(v.x,0,0,0) = (val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta);
 }
#line 253
 {





     val(v.y,0,0,0) = (val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta);
 }
#line 253
 {





     val(v.z,0,0,0) = (val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta);
 }}
    }
  } } end_foreach(); }
  boundary ((scalar *) g);
}
#line 281 "/home/davidtang/basilisk/src/utils.h"
void vorticity (const vector u, scalar omega)
{
   { 
if (!is_constant(fm.x) && !is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 283
foreach(){

#line 283 "/home/davidtang/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
if (is_constant(fm.x) && !is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 283
foreach(){

#line 283 "/home/davidtang/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
if (!is_constant(fm.x) && is_constant(cm)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 283
foreach(){

#line 283 "/home/davidtang/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); }
if (is_constant(fm.x) && is_constant(cm)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 283
foreach(){

#line 283 "/home/davidtang/basilisk/src/utils.h"

    val(omega,0,0,0) = ((val_fm_x(fm.x,1,0,0) - val_fm_x(fm.x,0,0,0))*val(u.y,0,0,0) +
        val_fm_x(fm.x,1,0,0)*val(u.y,1,0,0) - val_fm_x(fm.x,0,0,0)*val(u.y,-1,0,0) -
        (val_fm_y(fm.y,0,1,0) - val_fm_y(fm.y,0,0,0))*val(u.x,0,0,0) +
        val_fm_y(fm.y,0,0,0)*val(u.x,0,-1,0) - val_fm_y(fm.y,0,1,0)*val(u.x,0,1,0))/(2.*val_cm(cm,0,0,0)*Delta + 0.); } end_foreach(); } }
  boundary (((scalar []){omega,{-1}}));
}





double change (scalar s, scalar sn)
{
  double max = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _max = max; 
#line 298

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 298
foreach(){

#line 298 "/home/davidtang/basilisk/src/utils.h"
 {
    if ((cube(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > _max)
 _max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 298
foreach(){

#line 298 "/home/davidtang/basilisk/src/utils.h"
 {
    if ((cube(Delta)*val_cm(cm,0,0,0)) > 0.) {
      double ds = fabs (val(s,0,0,0) - val(sn,0,0,0));
      if (ds > _max)
 _max = ds;
    }
    val(sn,0,0,0) = val(s,0,0,0);
  } } end_foreach(); }OMP(omp critical) if (_max > max) max = _max;
mpi_all_reduce_double (max, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 305
 }
  return max;
}





scalar lookup_field (const char * name)
{
  if (name)
    if (all) for (scalar s = *all, *_i36 = all; ((scalar *)&s)->i >= 0; s = *++_i36)
      if (!strcmp (_attribute[s.i].name, name))
 return s;
  return (scalar){-1};
}

vector lookup_vector (const char * name)
{
  if (name) {
    char component[strlen(name) + 3];
    strcpy (component, name);
    strcat (component, ".x");
    if (all) for (scalar s = *all, *_i37 = all; ((scalar *)&s)->i >= 0; s = *++_i37)
      if (!strcmp (_attribute[s.i].name, component))
 return _attribute[s.i].v;
  }
  return (vector){{-1}};
}







#define foreach_segment(_S,_p) {\
  coord t = {(_S)[1].x - (_S)[0].x, (_S)[1].y - (_S)[0].y};\
  double norm = sqrt(sq(t.x) + sq(t.y));\
  assert (norm > 0.);\
  t.x = t.x/norm + 1e-6, t.y = t.y/norm - 1.5e-6;\
  double alpha = ((_S)[0].x*((_S)[1].y - (_S)[0].y) -\
    (_S)[0].y*((_S)[1].x - (_S)[0].x))/norm;\
  foreach()\
    if (fabs(t.y*x - t.x*y - alpha) < 0.708*Delta) {\
      coord _o = {x,y}, _p[2];\
      int _n = 0;\
      {\
 if (t.x)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].x = _o.x + _i*Delta/2.;\
     double a = (_p[_n].x - (_S)[0].x)/t.x;\
     _p[_n].y = (_S)[0].y + a*t.y;\
     if (fabs(_p[_n].y - _o.y) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].x = (_S)[0].x + a*t.x, _p[_n].y = (_S)[0].y + a*t.y;\
       if (fabs(_p[_n].x - _o.x) <= Delta/2. &&\
    fabs(_p[_n].y - _o.y) <= Delta/2.)\
  _n++;\
     }\
   }\
 if (t.y)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].y = _o.y + _i*Delta/2.;\
     double a = (_p[_n].y - (_S)[0].y)/t.y;\
     _p[_n].z = (_S)[0].z + a*t.z;\
     if (fabs(_p[_n].z - _o.z) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].y = (_S)[0].y + a*t.y, _p[_n].z = (_S)[0].z + a*t.z;\
       if (fabs(_p[_n].y - _o.y) <= Delta/2. &&\
    fabs(_p[_n].z - _o.z) <= Delta/2.)\
  _n++;\
     }\
   }\
 if (t.z)\
   for (int _i = -1; _i <= 1 && _n < 2; _i += 2) {\
     _p[_n].z = _o.z + _i*Delta/2.;\
     double a = (_p[_n].z - (_S)[0].z)/t.z;\
     _p[_n].x = (_S)[0].x + a*t.x;\
     if (fabs(_p[_n].x - _o.x) <= Delta/2.) {\
       a = clamp (a, 0., norm);\
       _p[_n].z = (_S)[0].z + a*t.z, _p[_n].x = (_S)[0].x + a*t.x;\
       if (fabs(_p[_n].z - _o.z) <= Delta/2. &&\
    fabs(_p[_n].x - _o.x) <= Delta/2.)\
  _n++;\
     }\
   }}\
      if (_n == 2) {\

#line 367

#define end_foreach_segment() } } end_foreach(); }




void fields_stats()
{
  fprintf (ferr, "# t = %g, fields = {", t);
  if (all) for (scalar s = *all, *_i38 = all; ((scalar *)&s)->i >= 0; s = *++_i38)
    fprintf (ferr, " %s", _attribute[s.i].name);
  fputs (" }\n", ferr);
  fprintf (ferr, "# %12s: %12s %12s %12s %12s\n",
    "name", "min", "avg", "stddev", "max");
  if (all) for (scalar s = *all, *_i39 = all; ((scalar *)&s)->i >= 0; s = *++_i39) {
    stats ss = statsf (s);
    fprintf (ferr, "# %12s: %12g %12g %12g %12g\n",
      _attribute[s.i].name, ss.min, ss.sum/ss.volume, ss.stddev, ss.max);
  }
}

#line 1 "./output.h"
#line 1 "/home/davidtang/basilisk/src/output.h"
#line 37 "/home/davidtang/basilisk/src/output.h"
struct OutputField {
  scalar * list;
  FILE * fp;
  int n;
  bool linear;
  double box[2][2];
};


void output_field (struct OutputField p)
{ trace ("output_field", "/home/davidtang/basilisk/src/output.h", 47);
  if (!p.list) p.list = all;
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  p.n++;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }

  int len = list_len(p.list);
  double Delta = 0.999999*(p.box[1][0] - p.box[0][0])/(p.n - 1);
  int ny = (p.box[1][1] - p.box[0][1])/Delta + 1;
  double ** field = (double **) matrix_new (p.n, ny, len*sizeof(double));
  for (int i = 0; i < p.n; i++) {
    double x = Delta*i + p.box[0][0];
    for (int j = 0; j < ny; j++) {
      double y = Delta*j + p.box[0][1];
      if (p.linear) {
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i40 = p.list; ((scalar *)&s)->i >= 0; s = *++_i40)
   field[i][len*j + k++] = interpolate ((struct _interpolate){s, x, y});
      }
      else {
 Point point = locate ((struct _locate){x, y});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 72 "/home/davidtang/basilisk/src/output.h"

 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i41 = p.list; ((scalar *)&s)->i >= 0; s = *++_i41)
   field[i][len*j + k++] = point.level >= 0 ? val(s,0,0,0) : nodata;
      }
    }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, field[0], len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif
    fprintf (p.fp, "# 1:x 2:y");
    int i = 3;
    if (p.list) for (scalar s = *p.list, *_i42 = p.list; ((scalar *)&s)->i >= 0; s = *++_i42)
      fprintf (p.fp, " %d:%s", i++, _attribute[s.i].name);
    fputc('\n', p.fp);
    for (int i = 0; i < p.n; i++) {
      double x = Delta*i + p.box[0][0];
      for (int j = 0; j < ny; j++) {
 double y = Delta*j + p.box[0][1];

 fprintf (p.fp, "%g %g", x, y);
 int k = 0;
 if (p.list) for (scalar s = *p.list, *_i43 = p.list; ((scalar *)&s)->i >= 0; s = *++_i43)
   fprintf (p.fp, " %g", field[i][len*j + k++]);
 fputc ('\n', p.fp);
      }
      fputc ('\n', p.fp);
    }
    fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (field[0], NULL, len*p.n*ny, MPI_DOUBLE, MPI_MIN, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (field);
 end_trace("output_field", "/home/davidtang/basilisk/src/output.h", 112); }
#line 140 "/home/davidtang/basilisk/src/output.h"
struct OutputMatrix {
  scalar f;
  FILE * fp;
  int n;
  bool linear;
};


void output_matrix (struct OutputMatrix p)
{ trace ("output_matrix", "/home/davidtang/basilisk/src/output.h", 149);
  if (p.n == 0) p.n = N;
  if (!p.fp) p.fp = fout;
  float fn = p.n;
  float Delta = (float) L0/fn;
  fwrite (&fn, sizeof(float), 1, p.fp);
  for (int j = 0; j < p.n; j++) {
    float yp = (float) (Delta*j + X0 + Delta/2.);
    fwrite (&yp, sizeof(float), 1, p.fp);
  }
  for (int i = 0; i < p.n; i++) {
    float xp = (float) (Delta*i + X0 + Delta/2.);
    fwrite (&xp, sizeof(float), 1, p.fp);
    for (int j = 0; j < p.n; j++) {
      float yp = (float)(Delta*j + Y0 + Delta/2.), v;
      if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 167 "/home/davidtang/basilisk/src/output.h"

 assert (point.level >= 0);
 v = val(p.f,0,0,0);
      }
      fwrite (&v, sizeof(float), 1, p.fp);
    }
  }
  fflush (p.fp);
 end_trace("output_matrix", "/home/davidtang/basilisk/src/output.h", 175); }
#line 184 "/home/davidtang/basilisk/src/output.h"
typedef void (* colormap) (double cmap[127][3]);

void jet (double cmap[127][3])
{
  for (int i = 0; i < 127; i++) {
    cmap[i][0] =
      i <= 46 ? 0. :
      i >= 111 ? -0.03125*(i - 111) + 1. :
      i >= 78 ? 1. :
      0.03125*(i - 46);
    cmap[i][1] =
      i <= 14 || i >= 111 ? 0. :
      i >= 79 ? -0.03125*(i - 111) :
      i <= 46 ? 0.03125*(i - 14) :
      1.;
    cmap[i][2] =
      i >= 79 ? 0. :
      i >= 47 ? -0.03125*(i - 79) :
      i <= 14 ? 0.03125*(i - 14) + 1.:
      1.;
  }
}

void cool_warm (double cmap[127][3])
{






  static double basemap[33][3] = {
    {0.2298057, 0.298717966, 0.753683153},
    {0.26623388, 0.353094838, 0.801466763},
    {0.30386891, 0.406535296, 0.84495867},
    {0.342804478, 0.458757618, 0.883725899},
    {0.38301334, 0.50941904, 0.917387822},
    {0.424369608, 0.558148092, 0.945619588},
    {0.46666708, 0.604562568, 0.968154911},
    {0.509635204, 0.648280772, 0.98478814},
    {0.552953156, 0.688929332, 0.995375608},
    {0.596262162, 0.726149107, 0.999836203},
    {0.639176211, 0.759599947, 0.998151185},
    {0.681291281, 0.788964712, 0.990363227},
    {0.722193294, 0.813952739, 0.976574709},
    {0.761464949, 0.834302879, 0.956945269},
    {0.798691636, 0.849786142, 0.931688648},
    {0.833466556, 0.860207984, 0.901068838},
    {0.865395197, 0.86541021, 0.865395561},
    {0.897787179, 0.848937047, 0.820880546},
    {0.924127593, 0.827384882, 0.774508472},
    {0.944468518, 0.800927443, 0.726736146},
    {0.958852946, 0.769767752, 0.678007945},
    {0.96732803, 0.734132809, 0.628751763},
    {0.969954137, 0.694266682, 0.579375448},
    {0.966811177, 0.650421156, 0.530263762},
    {0.958003065, 0.602842431, 0.481775914},
    {0.943660866, 0.551750968, 0.434243684},
    {0.923944917, 0.49730856, 0.387970225},
    {0.89904617, 0.439559467, 0.343229596},
    {0.869186849, 0.378313092, 0.300267182},
    {0.834620542, 0.312874446, 0.259301199},
    {0.795631745, 0.24128379, 0.220525627},
    {0.752534934, 0.157246067, 0.184115123},
    {0.705673158, 0.01555616, 0.150232812}
  };

  for (int i = 0; i < 127; i++) {
    double x = i*(32 - 1e-10)/(127 - 1);
    int j = x; x -= j;
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (1. - x)*basemap[j][k] + x*basemap[j+1][k];
  }
}

void gray (double cmap[127][3])
{
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = i/(127 - 1.);
}

void randomap (double cmap[127][3])
{
  srand(0);
  for (int i = 0; i < 127; i++)
    for (int k = 0; k < 3; k++)
      cmap[i][k] = (noise() + 1.)/2.;
}

void blue_white_red (double cmap[127][3])
{
  for (int i = 0; i < (127 + 1)/2; i++) {
    cmap[i][0] = i/((127 - 1)/2.);
    cmap[i][1] = i/((127 - 1)/2.);
    cmap[i][2] = 1.;
  }
  for (int i = 0; i < (127 - 1)/2; i++) {
    cmap[i + (127 + 1)/2][0] = 1.;
    cmap[i + (127 + 1)/2][1] = cmap[(127 - 3)/2 - i][1];
    cmap[i + (127 + 1)/2][2] = cmap[(127 - 3)/2 - i][1];
  }
}





typedef struct {
  unsigned char r, g, b;
} color;

color colormap_color (double cmap[127][3],
        double val, double min, double max)
{
  color c;
  if (val == nodata) {
    c.r = c.g = c.b = 0;
    return c;
  }
  int i;
  double coef;
  if (max != min)
    val = (val - min)/(max - min);
  else
    val = 0.;
  if (val <= 0.) i = 0, coef = 0.;
  else if (val >= 1.) i = 127 - 2, coef = 1.;
  else {
    i = val*(127 - 1);
    coef = val*(127 - 1) - i;
  }
  assert (i < 127 - 1);
  unsigned char * c1 = (unsigned char *) &c;
  for (int j = 0; j < 3; j++)
    c1[j] = 255*(cmap[i][j]*(1. - coef) + cmap[i + 1][j]*coef);
  return c;
}
#line 335 "/home/davidtang/basilisk/src/output.h"
static const char * extension (const char * file, const char * ext) {
  int len = strlen(file);
  return len > 4 && !strcmp (file + len - 4, ext) ? file + len - 4 : NULL;
}

static const char * is_animation (const char * file) {
  const char * ext;
  if ((ext = extension (file, ".mp4")) ||
      (ext = extension (file, ".ogv")) ||
      (ext = extension (file, ".gif")))
    return ext;
  return NULL;
}

static struct {
  FILE ** fp;
  char ** names;
  int n;
} open_image_data = {NULL, NULL, 0};

static void open_image_cleanup()
{
  for (int i = 0; i < open_image_data.n; i++) {
    qpclose (open_image_data.fp[i]);
    pfree (open_image_data.names[i],__func__,__FILE__,__LINE__);
  }
  pfree (open_image_data.fp,__func__,__FILE__,__LINE__);
  pfree (open_image_data.names,__func__,__FILE__,__LINE__);
  open_image_data.fp = NULL;
  open_image_data.names = NULL;
  open_image_data.n = 0;
}

static FILE * open_image_lookup (const char * file)
{
  for (int i = 0; i < open_image_data.n; i++)
    if (!strcmp (file, open_image_data.names[i]))
      return open_image_data.fp[i];
  return NULL;
}

static bool which (const char * command)
{
  char * s = getenv ("PATH");
  if (!s)
    return false;
  char path[strlen(s) + 1];
  strcpy (path, s);
  s = strtok (path, ":");
  while (s) {
    char f[strlen(s) + strlen(command) + 2];
    strcpy (f, s);
    strcat (f, "/");
    strcat (f, command);
    FILE * fp = fopen (f, "r");
    if (fp) {
      fclose (fp);
      return true;
    }
    s = strtok (NULL, ":");
  }
  return false;
}

static FILE * ppm_fallback (const char * file, const char * mode)
{
  char filename[strlen(file) + 5];
  strcpy (filename, file);
  strcat (filename, ".ppm");
  FILE * fp = fopen (filename, mode);
  if (!fp) {
    perror (file);

    MPI_Abort (MPI_COMM_WORLD, 1);

    exit (1);
  }
  return fp;
}

FILE * open_image (const char * file, const char * options)
{
  assert (pid() == 0);
  const char * ext;
  if ((ext = is_animation (file))) {
    FILE * fp = open_image_lookup (file);
    if (fp)
      return fp;

    int len = strlen ("ppm2???    ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "ppm2"); strcat (command, ext + 1);

    static int has_ffmpeg = -1;
    if (has_ffmpeg < 0) {
      if (which (command) && (which ("ffmpeg") || which ("avconv")))
 has_ffmpeg = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find '%s' or 'ffmpeg'/'avconv'\n"
   "  falling back to raw PPM outputs\n", command);
 has_ffmpeg = false;
      }
    }
    if (!has_ffmpeg)
      return ppm_fallback (file, "a");

    static bool added = false;
    if (!added) {
      free_solver_func_add (open_image_cleanup);
      added = true;
    }
    open_image_data.n++;
    open_image_data.names = (char * *) prealloc (open_image_data.names, (open_image_data.n)*sizeof(char *),__func__,__FILE__,__LINE__);
    open_image_data.names[open_image_data.n - 1] = pstrdup (file,__func__,__FILE__,__LINE__);

    if (options) {
      strcat (command, " ");
      strcat (command, options);
    }
    strcat (command, !strcmp (ext, ".mp4") ? " " : " > ");
    strcat (command, file);
    open_image_data.fp = (FILE * *) prealloc (open_image_data.fp, (open_image_data.n)*sizeof(FILE *),__func__,__FILE__,__LINE__);
    return open_image_data.fp[open_image_data.n - 1] = qpopen (command, "w");
  }
  else {
    static int has_convert = -1;
    if (has_convert < 0) {
      if (which ("convert"))
 has_convert = true;
      else {
 fprintf (ferr,
   "open_image(): cannot find 'convert'\n"
   "  falling back to raw PPM outputs\n");
 has_convert = false;
      }
    }
    if (!has_convert)
      return ppm_fallback (file, "w");

    int len = strlen ("convert ppm:-   ") + strlen (file) +
      (options ? strlen (options) : 0);
    char command[len];
    strcpy (command, "convert ppm:- ");
    if (options) {
      strcat (command, options);
      strcat (command, " ");
    }
    strcat (command, file);
    return qpopen (command, "w");
  }
}

void close_image (const char * file, FILE * fp)
{
  assert (pid() == 0);
  if (is_animation (file)) {
    if (!open_image_lookup (file))
      fclose (fp);
  }
  else if (which ("convert"))
    qpclose (fp);
  else
    fclose (fp);
}
#line 566 "/home/davidtang/basilisk/src/output.h"
struct OutputPPM {
  scalar f;
  FILE * fp;
  int n;
  char * file;
  double min, max, spread, z;
  bool linear;
  double box[2][2];
  scalar mask;
  colormap map;
  char * opt;
};


void output_ppm (struct OutputPPM p)
{ trace ("output_ppm", "/home/davidtang/basilisk/src/output.h", 581);

  if (p.n == 0) p.n = N;
  if (p.min == 0 && p.max == 0) {
    stats s = statsf (p.f);
    if (p.spread < 0.)
      p.min = s.min, p.max = s.max;
    else {
      double avg = s.sum/s.volume, spread = (p.spread ? p.spread : 5.)*s.stddev;
      p.min = avg - spread; p.max = avg + spread;
    }
  }
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
  }
  if (!p.map)
    p.map = jet;

  double fn = p.n;
  double Delta = (p.box[1][0] - p.box[0][0])/fn;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;
  if (ny % 2) ny++;

  color ** ppm = (color **) matrix_new (ny, p.n, sizeof(color));
  double cmap[127][3];
  p.map (cmap);
  OMP_PARALLEL() {
    OMP(omp for schedule(static))
      for (int j = 0; j < ny; j++) {
 double yp = Delta*j + p.box[0][1] + Delta/2.;
 for (int i = 0; i < p.n; i++) {
   double xp = Delta*i + p.box[0][0] + Delta/2., v;
   if (p.mask.i) {
     if (p.linear) {
       double m = interpolate ((struct _interpolate){p.mask, xp, yp, p.z});
       if (m < 0.)
  v = nodata;
       else
  v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
     }
     else {
       Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 624 "/home/davidtang/basilisk/src/output.h"

       if (point.level < 0 || val(p.mask,0,0,0) < 0.)
  v = nodata;
       else
  v = val(p.f,0,0,0);
     }
   }
   else if (p.linear)
     v = interpolate ((struct _interpolate){p.f, xp, yp, p.z});
   else {
     Point point = locate ((struct _locate){xp, yp, p.z});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 634 "/home/davidtang/basilisk/src/output.h"

     v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
   }
   ppm[ny - 1 - j][i] = colormap_color (cmap, v, p.min, p.max);
 }
      }
  }

  if (pid() == 0) {
#if 1
    MPI_Reduce (MPI_IN_PLACE, ppm[0], 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif
    if (!p.fp) p.fp = fout;
    if (p.file)
      p.fp = open_image (p.file, p.opt);

    fprintf (p.fp, "P6\n%u %u 255\n", p.n, ny);
    fwrite (((void **) ppm)[0], sizeof(color), ny*p.n, p.fp);

    if (p.file)
      close_image (p.file, p.fp);
    else
      fflush (p.fp);
  }
#if 1
  else
    MPI_Reduce (ppm[0], NULL, 3*ny*p.n, MPI_UNSIGNED_CHAR, MPI_MAX, 0,
  MPI_COMM_WORLD);
#endif

  matrix_free (ppm);
 end_trace("output_ppm", "/home/davidtang/basilisk/src/output.h", 666); }
#line 698 "/home/davidtang/basilisk/src/output.h"
struct OutputGRD {
  scalar f;
  FILE * fp;
  double Delta;
  bool linear;
  double box[2][2];
  scalar mask;
};


void output_grd (struct OutputGRD p)
{ trace ("output_grd", "/home/davidtang/basilisk/src/output.h", 709);

  if (!p.fp) p.fp = fout;
  if (p.box[0][0] == 0. && p.box[0][1] == 0. &&
      p.box[1][0] == 0. && p.box[1][1] == 0.) {
    p.box[0][0] = X0; p.box[0][1] = Y0;
    p.box[1][0] = X0 + L0; p.box[1][1] = Y0 + L0;
    if (p.Delta == 0) p.Delta = L0/N;
  }

  double Delta = p.Delta;
  int nx = (p.box[1][0] - p.box[0][0])/Delta;
  int ny = (p.box[1][1] - p.box[0][1])/Delta;


  fprintf (p.fp, "ncols          %d\n", nx);
  fprintf (p.fp, "nrows          %d\n", ny);
  fprintf (p.fp, "xllcorner      %g\n", p.box[0][0]);
  fprintf (p.fp, "yllcorner      %g\n", p.box[0][1]);
  fprintf (p.fp, "cellsize       %g\n", Delta);
  fprintf (p.fp, "nodata_value   -9999\n");


  for (int j = ny-1; j >= 0; j--) {
    double yp = Delta*j + p.box[0][1] + Delta/2.;
    for (int i = 0; i < nx; i++) {
      double xp = Delta*i + p.box[0][0] + Delta/2., v;
      if (p.mask.i) {
 if (p.linear) {
   double m = interpolate ((struct _interpolate){p.mask, xp, yp});
   if (m < 0.)
     v = nodata;
   else
     v = interpolate ((struct _interpolate){p.f, xp, yp});
 }
 else {
   Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 745 "/home/davidtang/basilisk/src/output.h"

   if (point.level < 0 || val(p.mask,0,0,0) < 0.)
     v = nodata;
   else
     v = val(p.f,0,0,0);
 }
      }
      else if (p.linear)
 v = interpolate ((struct _interpolate){p.f, xp, yp});
      else {
 Point point = locate ((struct _locate){xp, yp});  int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 755 "/home/davidtang/basilisk/src/output.h"

 v = point.level >= 0 ? val(p.f,0,0,0) : nodata;
      }
      if (v == nodata)
 fprintf (p.fp, "-9999 ");
      else
 fprintf (p.fp, "%f ", v);
    }
    fprintf (p.fp, "\n");
  }

  fflush (p.fp);
 end_trace("output_grd", "/home/davidtang/basilisk/src/output.h", 767); }
#line 794 "/home/davidtang/basilisk/src/output.h"
struct OutputGfs {
  FILE * fp;
  scalar * list;
  double t;
  char * file;
  bool translate;
};

static char * replace (const char * input, int target, int with,
         bool translate)
{
  if (translate) {
    if (!strcmp (input, "u.x"))
      return pstrdup ("U",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.y"))
      return pstrdup ("V",__func__,__FILE__,__LINE__);
    if (!strcmp (input, "u.z"))
      return pstrdup ("W",__func__,__FILE__,__LINE__);
  }
  char * name = pstrdup (input,__func__,__FILE__,__LINE__), * i = name;
  while (*i != '\0') {
    if (*i == target)
      *i = with;
    i++;
  }
  return name;
}


void output_gfs (struct OutputGfs p)
{ trace ("output_gfs", "/home/davidtang/basilisk/src/output.h", 824);
  char * fname = p.file;

#if 1

  not_mpi_compatible();

  FILE * fp = p.fp;
  if (p.file == NULL) {
    long pid = getpid();
    MPI_Bcast (&pid, 1, MPI_LONG, 0, MPI_COMM_WORLD);
    fname = ((char *) pmalloc ((80)*sizeof(char),__func__,__FILE__,__LINE__));
    snprintf (fname, 80, ".output-%ld", pid);
    p.fp = NULL;
  }
#endif

  bool opened = false;
  if (p.fp == NULL) {
    if (fname == NULL)
      p.fp = fout;
    else if (!(p.fp = fopen (fname, "w"))) {
      perror (fname);
      exit (1);
    }
    else
      opened = true;
  }

  scalar * list = p.list ? p.list : list_copy (all);

  restriction (list);
  fprintf (p.fp,
    "1 0 GfsSimulation GfsBox GfsGEdge { binary = 1"
    " x = %g y = %g ",
    0.5 + X0/L0, 0.5 + Y0/L0);

  fprintf (p.fp, "z = %g ", 0.5 + Z0/L0);


  if (list != NULL && list[0].i != -1) {
    scalar s = list[0];
    char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
    fprintf (p.fp, "variables = %s", name);
    pfree (name,__func__,__FILE__,__LINE__);
    for (int i = 1; i < list_len(list); i++) {
      scalar s = list[i];
      if (_attribute[s.i].name) {
 char * name = replace (_attribute[s.i].name, '.', '_', p.translate);
 fprintf (p.fp, ",%s", name);
 pfree (name,__func__,__FILE__,__LINE__);
      }
    }
    fprintf (p.fp, " ");
  }
  fprintf (p.fp, "} {\n");
  fprintf (p.fp, "  Time { t = %g }\n", t);
  if (L0 != 1.)
    fprintf (p.fp, "  PhysicalParams { L = %g }\n", L0);
  fprintf (p.fp, "  VariableTracerVOF f\n");
  fprintf (p.fp, "}\nGfsBox { x = 0 y = 0 z = 0 } {\n");

#if 1
  long header;
  if ((header = ftell (p.fp)) < 0) {
    perror ("output_gfs(): error in header");
    exit (1);
  }
  int cell_size = sizeof(unsigned) + sizeof(double);
  if (list) for (scalar s = *list, *_i44 = list; ((scalar *)&s)->i >= 0; s = *++_i44)
    if (_attribute[s.i].name)
      cell_size += sizeof(double);
  scalar index = new_scalar("index");
  size_t total_size = header + (z_indexing (index, false) + 1)*cell_size;
#endif



   { foreach_cell(){

#line 902 "/home/davidtang/basilisk/src/output.h"
 {
#if 1
    if (is_local(cell))
#endif
    {
#if 1
      if (fseek (p.fp, header + val(index,0,0,0)*cell_size, SEEK_SET) < 0) {
 perror ("output_gfs(): error while seeking");
 exit (1);
      }
#endif
      unsigned flags =
 level == 0 ? 0 :
#line 923 "/home/davidtang/basilisk/src/output.h"
      child.x == -1 && child.y == -1 && child.z == -1 ? 0 :
 child.x == -1 && child.y == -1 && child.z == 1 ? 1 :
 child.x == -1 && child.y == 1 && child.z == -1 ? 2 :
 child.x == -1 && child.y == 1 && child.z == 1 ? 3 :
 child.x == 1 && child.y == -1 && child.z == -1 ? 4 :
 child.x == 1 && child.y == -1 && child.z == 1 ? 5 :
 child.x == 1 && child.y == 1 && child.z == -1 ? 6 :
 7;

      if (is_leaf(cell))
 flags |= (1 << 4);
      fwrite (&flags, sizeof (unsigned), 1, p.fp);
      double a = -1;
      fwrite (&a, sizeof (double), 1, p.fp);
      if (list) for (scalar s = *list, *_i45 = list; ((scalar *)&s)->i >= 0; s = *++_i45)
 if (_attribute[s.i].name) {
   if (_attribute[s.i].v.x.i >= 0) {




     if (_attribute[s.i].v.x.i == s.i) {
       s = _attribute[s.i].v.y;
       a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
     }
     else if (_attribute[s.i].v.y.i == s.i) {
       s = _attribute[s.i].v.x;
       a = is_local(cell) && val(s,0,0,0) != nodata ? - val(s,0,0,0) : (double) DBL_MAX;
     }


     else
       a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;

   }
   else
     a = is_local(cell) && val(s,0,0,0) != nodata ? val(s,0,0,0) : (double) DBL_MAX;
   fwrite (&a, sizeof (double), 1, p.fp);
 }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

#if 1
  delete (((scalar []){index,{-1}}));
  if (!pid() && fseek (p.fp, total_size, SEEK_SET) < 0) {
    perror ("output_gfs(): error while finishing");
    exit (1);
  }
  if (!pid())
#endif
    fputs ("}\n", p.fp);
  fflush (p.fp);

  if (!p.list)
    pfree (list,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);

#if 1
  if (p.file == NULL) {
    MPI_Barrier (MPI_COMM_WORLD);
    if (pid() == 0) {
      if (fp == NULL)
 fp = fout;
      p.fp = fopen (fname, "r");
      size_t l;
      unsigned char buffer[8192];
      while ((l = fread (buffer, 1, 8192, p.fp)) > 0)
 fwrite (buffer, 1, l, fp);
      fflush (fp);
      remove (fname);
    }
    pfree (fname,__func__,__FILE__,__LINE__);
  }
#endif
 end_trace("output_gfs", "/home/davidtang/basilisk/src/output.h", 1000); }
#line 1024 "/home/davidtang/basilisk/src/output.h"
struct Dump {
  char * file;
  scalar * list;
  FILE * fp;
  bool unbuffered;
};

struct DumpHeader {
  double t;
  long len;
  int i, depth, npe, version;
  coord n;
};

static const int dump_version =

  170901;

static scalar * dump_list (scalar * lista)
{
  scalar * list = is_constant(cm) ? NULL : list_concat (((scalar []){cm,{-1}}), NULL);
  if (lista) for (scalar s = *lista, *_i46 = lista; ((scalar *)&s)->i >= 0; s = *++_i46)
    if (!_attribute[s.i].face && !_attribute[s.i].nodump && s.i != cm.i)
      list = list_add (list, s);
  return list;
}

static void dump_header (FILE * fp, struct DumpHeader * header, scalar * list)
{
  if (fwrite (header, sizeof(struct DumpHeader), 1, fp) < 1) {
    perror ("dump(): error while writing header");
    exit (1);
  }
  if (list) for (scalar s = *list, *_i47 = list; ((scalar *)&s)->i >= 0; s = *++_i47) {
    unsigned len = strlen(_attribute[s.i].name);
    if (fwrite (&len, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing len");
      exit (1);
    }
    if (fwrite (_attribute[s.i].name, sizeof(char), len, fp) < len) {
      perror ("dump(): error while writing s.name");
      exit (1);
    }
  }
  double o[4] = {X0,Y0,Z0,L0};
  if (fwrite (o, sizeof(double), 4, fp) < 4) {
    perror ("dump(): error while writing coordinates");
    exit (1);
  }
}

#if !1

void dump (struct Dump p)
{ trace ("dump", "/home/davidtang/basilisk/src/output.h", 1078);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  char * name = NULL;
  if (file) {
    name = (char *) pmalloc (strlen(file) + 2,__func__,__FILE__,__LINE__);
    strcpy (name, file);
    if (!p.unbuffered)
      strcat (name, "~");
    if ((fp = fopen (name, "w")) == NULL) {
      perror (name);
      exit (1);
    }
  }
  assert (fp);

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };
  dump_header (fp, &header, list);

  subtree_size (size, false);

   { foreach_cell(){

#line 1104 "/home/davidtang/basilisk/src/output.h"
 {
    unsigned flags = is_leaf(cell) ? leaf : 0;
    if (fwrite (&flags, sizeof(unsigned), 1, fp) < 1) {
      perror ("dump(): error while writing flags");
      exit (1);
    }
    if (list) for (scalar s = *list, *_i48 = list; ((scalar *)&s)->i >= 0; s = *++_i48)
      if (fwrite (&val(s,0,0,0), sizeof(double), 1, fp) < 1) {
 perror ("dump(): error while writing scalars");
 exit (1);
      }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  pfree (list,__func__,__FILE__,__LINE__);
  if (file) {
    fclose (fp);
    if (!p.unbuffered)
      rename (name, file);
    pfree (name,__func__,__FILE__,__LINE__);
  }
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/davidtang/basilisk/src/output.h", 1126); }
#else

void dump (struct Dump p)
{ trace ("dump", "/home/davidtang/basilisk/src/output.h", 1130);
  FILE * fp = p.fp;
  char def[] = "dump", * file = p.file ? p.file : p.fp ? NULL : def;

  if (fp != NULL || file == NULL) {
    fprintf (ferr, "dump(): must specify a file name when using MPI\n");
    exit(1);
  }

  char name[strlen(file) + 2];
  strcpy (name, file);
  if (!p.unbuffered)
    strcat (name, "~");
  FILE * fh = fopen (name, "w");
  if (fh == NULL) {
    perror (name);
    exit (1);
  }

  scalar * dlist = dump_list (p.list ? p.list : all);
  scalar size= new_scalar("size");
  scalar * list = list_concat (((scalar []){size,{-1}}), dlist); pfree (dlist,__func__,__FILE__,__LINE__);
  struct DumpHeader header = { t, list_len(list), iter, depth(), npe(),
          dump_version };


  for (int i = 0; i < 3; i++)
    (&header.n.x)[i] = mpi_dims[i];
  MPI_Barrier (MPI_COMM_WORLD);


  if (pid() == 0)
    dump_header (fh, &header, list);

  scalar index = {-1};

  index = new_scalar("index");
  z_indexing (index, false);
  int cell_size = sizeof(unsigned) + header.len*sizeof(double);
  int sizeofheader = sizeof(header) + 4*sizeof(double);
  if (list) for (scalar s = *list, *_i49 = list; ((scalar *)&s)->i >= 0; s = *++_i49)
    sizeofheader += sizeof(unsigned) + sizeof(char)*strlen(_attribute[s.i].name);
  long pos = pid() ? 0 : sizeofheader;

  subtree_size (size, false);

   { foreach_cell(){

#line 1176 "/home/davidtang/basilisk/src/output.h"
 {

    if (is_local(cell)) {
      long offset = sizeofheader + val(index,0,0,0)*cell_size;
      if (pos != offset) {
 fseek (fh, offset, SEEK_SET);
 pos = offset;
      }
      unsigned flags = is_leaf(cell) ? leaf : 0;
      fwrite (&flags, 1, sizeof(unsigned), fh);
      if (list) for (scalar s = *list, *_i50 = list; ((scalar *)&s)->i >= 0; s = *++_i50)
 fwrite (&val(s,0,0,0), 1, sizeof(double), fh);
      pos += cell_size;
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }

  delete (((scalar []){index,{-1}}));

  pfree (list,__func__,__FILE__,__LINE__);
  fclose (fh);
  if (!p.unbuffered && pid() == 0)
    rename (name, file);
 delete (((scalar []){size,{-1}}));  end_trace("dump", "/home/davidtang/basilisk/src/output.h", 1200); }
#endif


bool restore (struct Dump p)
{ trace ("restore", "/home/davidtang/basilisk/src/output.h", 1205);
  FILE * fp = p.fp;
  char * file = p.file;
  if (file && (fp = fopen (file, "r")) == NULL)
    { bool _ret =  false; end_trace("restore", "/home/davidtang/basilisk/src/output.h", 1209);  return _ret; }
  assert (fp);

  struct DumpHeader header;
  if (fread (&header, sizeof(header), 1, fp) < 1) {
    fprintf (ferr, "restore(): error: expecting header\n");
    exit (1);
  }
#line 1227 "/home/davidtang/basilisk/src/output.h"
  if (header.npe != npe()) {
    fprintf (ferr,
      "restore(): error: the number of processes don't match:"
      " %d != %d\n",
      header.npe, npe());
    exit (1);
  }
  dimensions ((struct Dimensions){header.n.x, header.n.y, header.n.z});
  double n = header.n.x;
  int depth = header.depth;
  while (n > 1)
    depth++, n /= 2;
  init_grid (1 << depth);





  bool restore_all = (p.list == all);
  scalar * list = dump_list (p.list ? p.list : all);
  if (header.version == 161020) {
    if (header.len - 1 != list_len (list)) {
      fprintf (ferr,
        "restore(): error: the list lengths don't match: "
        "%ld (file) != %d (code)\n",
        header.len - 1, list_len (list));
      exit (1);
    }
  }
  else {
    if (header.version != dump_version) {
      fprintf (ferr,
        "restore(): error: file version mismatch: "
        "%d (file) != %d (code)\n",
        header.version, dump_version);
      exit (1);
    }

    scalar * input = NULL;
    for (int i = 0; i < header.len; i++) {
      unsigned len;
      if (fread (&len, sizeof(unsigned), 1, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting len\n");
 exit (1);
      }
      char name[len + 1];
      if (fread (name, sizeof(char), len, fp) < 1) {
 fprintf (ferr, "restore(): error: expecting s.name\n");
 exit (1);
      }
      name[len] = '\0';

      if (i > 0) {
 bool found = false;
 if (list) for (scalar s = *list, *_i51 = list; ((scalar *)&s)->i >= 0; s = *++_i51)
   if (!strcmp (_attribute[s.i].name, name)) {
     input = list_append (input, s);
     found = true; break;
   }
 if (!found) {
   if (restore_all) {
     scalar s = new_scalar("s");
     pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
     _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
     input = list_append (input, s);
   }
   else
     input = list_append (input, (scalar){INT_MAX});
 }
      }
    }
    pfree (list,__func__,__FILE__,__LINE__);
    list = input;

    double o[4];
    if (fread (o, sizeof(double), 4, fp) < 4) {
      fprintf (ferr, "restore(): error: expecting coordinates\n");
      exit (1);
    }
    origin ((struct _origin){o[0], o[1], o[2]});
    size (o[3]);
  }


  long cell_size = sizeof(unsigned) + header.len*sizeof(double);
  long offset = pid()*((1 << 3*(header.depth + 1)) - 1)/
    ((1 << 3) - 1)*cell_size;
  if (fseek (fp, offset, SEEK_CUR) < 0) {
    perror ("restore(): error while seeking");
    exit (1);
  }


  scalar * listm = is_constant(cm) ? NULL : (scalar *)((vector []){{fm.x,fm.y,fm.z},{{-1},{-1},{-1}}});



   { foreach_cell(){

#line 1324 "/home/davidtang/basilisk/src/output.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof(unsigned), 1, fp) != 1) {
      fprintf (ferr, "restore(): error: expecting 'flags'\n");
      exit (1);
    }

    fseek (fp, sizeof(double), SEEK_CUR);
    if (list) for (scalar s = *list, *_i52 = list; ((scalar *)&s)->i >= 0; s = *++_i52) {
      double val;
      if (fread (&val, sizeof(double), 1, fp) != 1) {
 fprintf (ferr, "restore(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX)
 val(s,0,0,0) = val;
    }
    if (!(flags & leaf) && is_leaf(cell))
      refine_cell (point, listm, 0, NULL);
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary (list);

  boundary (listm);

  scalar * other = NULL;
  if (all) for (scalar s = *all, *_i53 = all; ((scalar *)&s)->i >= 0; s = *++_i53)
    if (!list_lookup (list, s) && !list_lookup (listm, s))
      other = list_append (other, s);
  reset (other, 0.);
  pfree (other,__func__,__FILE__,__LINE__);

  pfree (list,__func__,__FILE__,__LINE__);
  if (file)
    fclose (fp);


  while (iter < header.i && events (false))
    iter = inext;
  events (false);
  while (t < header.t && events (false))
    t = tnext;
  t = header.t;
  events (false);

  { bool _ret =  true; end_trace("restore", "/home/davidtang/basilisk/src/output.h", 1370);  return _ret; }
 end_trace("restore", "/home/davidtang/basilisk/src/output.h", 1371); }
#line 389 "/home/davidtang/basilisk/src/utils.h"
#line 12 "/home/davidtang/basilisk/src/run.h"


void run (void)
{ trace ("run", "/home/davidtang/basilisk/src/run.h", 15);
  iter = 0, t = 0., dt = 1.;
  init_grid (N);

  perf.nc = perf.tnc = 0;
  perf.gt = timer_start();
  while (events (true)) {





    update_perf();
    iter = inext, t = tnext;
  }




  timer_print (perf.gt, iter, perf.tnc);

  free_grid();
 end_trace("run", "/home/davidtang/basilisk/src/run.h", 37); }




static int defaults_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults (const int i, const double t, Event * _ev) { trace ("defaults", "/home/davidtang/basilisk/src/run.h", 42);  {
  display ((struct _display){"box();"});
 end_trace("defaults", "/home/davidtang/basilisk/src/run.h", 44); } return 0; } 





static int cleanup_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t = 1234567890);   *ip = i; *tp = t;   return ret; } static int cleanup (const int i, const double t, Event * _ev) { trace ("cleanup", "/home/davidtang/basilisk/src/run.h", 50);  {
  display ((struct _display){"", true});
 end_trace("cleanup", "/home/davidtang/basilisk/src/run.h", 52); } return 0; } 
#line 28 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
#line 1 "./timestep.h"
#line 1 "/home/davidtang/basilisk/src/timestep.h"

double timestep (const vector u, double dtmax)
{
  static double previous = 0.;
  dtmax /= CFL;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _dtmax = dtmax; 
#line 6

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/davidtang/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/davidtang/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/davidtang/basilisk/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.z,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 16
 end_foreach_face(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 6
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 6
{

#line 6 "/home/davidtang/basilisk/src/timestep.h"

    if (val(u.x,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.x,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 6
{

#line 6 "/home/davidtang/basilisk/src/timestep.h"

    if (val(u.y,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.y,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < _dtmax) _dtmax = dt;
    } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 6
{

#line 6 "/home/davidtang/basilisk/src/timestep.h"

    if (val(u.z,0,0,0) != 0.) {
      double dt = Delta/fabs(val(u.z,0,0,0));




      dt *= val_cm(cm,0,0,0);

      if (dt < _dtmax) _dtmax = dt;
    } }  }}  end_foreach_face_generic()
#line 16
 end_foreach_face(); }OMP(omp critical) if (_dtmax < dtmax) dtmax = _dtmax;
mpi_all_reduce_double (dtmax, MPI_MIN);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 16
 }
  dtmax *= CFL;
  if (dtmax > previous)
    dtmax = (previous + 0.1*dtmax)/1.1;
  previous = dtmax;
  return dtmax;
}
#line 29 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
#line 1 "./bcg.h"
#line 1 "/home/davidtang/basilisk/src/bcg.h"
#line 11 "/home/davidtang/basilisk/src/bcg.h"
void tracer_fluxes (scalar f,
      vector uf,
      vector flux,
      double dt,
       scalar src)
{





  vector g= new_vector("g");
  gradients (((scalar []){f,{-1}}), ((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));




   { 
if (!is_constant(fm.x) && !is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_z(fm.z,i,0,0) && val_fm_z(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val_fm_z(fm.z,i,0,0) + val_fm_z(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_z(fm.z,0,i,0) && val_fm_z(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val_fm_z(fm.z,0,i,0) + val_fm_z(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val_fm_x(fm.x,0,0,i) && val_fm_x(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val_fm_x(fm.x,0,0,i) + val_fm_x(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_y(fm.y,0,0,i) && val_fm_y(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val_fm_y(fm.y,0,0,i) + val_fm_y(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_src
#define val_src(a,i,j,k) val(a,i,j,k)
#undef fine_src
#define fine_src(a,i,j,k) fine(a,i,j,k)
#undef coarse_src
#define coarse_src(a,i,j,k) coarse(a,i,j,k)
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_z(fm.z,i,0,0) && val_fm_z(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val_fm_z(fm.z,i,0,0) + val_fm_z(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_z(fm.z,0,i,0) && val_fm_z(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val_fm_z(fm.z,0,i,0) + val_fm_z(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val_fm_x(fm.x,0,0,i) && val_fm_x(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val_fm_x(fm.x,0,0,i) + val_fm_x(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_y(fm.y,0,0,i) && val_fm_y(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val_fm_y(fm.y,0,0,i) + val_fm_y(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(src)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_z(fm.z,i,0,0) && val_fm_z(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val_fm_z(fm.z,i,0,0) + val_fm_z(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_z(fm.z,0,i,0) && val_fm_z(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val_fm_z(fm.z,0,i,0) + val_fm_z(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val_fm_x(fm.x,0,0,i) && val_fm_x(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val_fm_x(fm.x,0,0,i) + val_fm_x(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_y(fm.y,0,0,i) && val_fm_y(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val_fm_y(fm.y,0,0,i) + val_fm_y(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(src)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const double _const_src = _constant[src.i -_NVARMAX];
NOT_UNUSED(_const_src);
#undef val_src
#define val_src(a,i,j,k) _const_src
#undef fine_src
#define fine_src(a,i,j,k) _const_src
#undef coarse_src
#define coarse_src(a,i,j,k) _const_src
#line 28
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.x,0,0,0)/(val_fm_x(fm.x,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,i,0,0) + (val_src(src,0,0,0) + val_src(src,-1,0,0))*dt/4. + s*(1. - s*un)*val(g.x,i,0,0)*Delta/2.;





    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double vn = (val(uf.y,i,0,0) + val(uf.y,i,1,0))/(val_fm_y(fm.y,i,0,0) + val_fm_y(fm.y,i,1,0));
      double fyy = vn < 0. ? val(f,i,1,0) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,-1,0);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_z(fm.z,i,0,0) && val_fm_z(fm.z,i,0,1)) {
      double wn = (val(uf.z,i,0,0) + val(uf.z,i,0,1))/(val_fm_z(fm.z,i,0,0) + val_fm_z(fm.z,i,0,1));
      double fzz = wn < 0. ? val(f,i,0,1) - val(f,i,0,0) : val(f,i,0,0) - val(f,i,0,-1);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.x,0,0,0) = f2*val(uf.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.y,0,0,0)/(val_fm_y(fm.y,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,i,0) + (val_src(src,0,0,0) + val_src(src,0,-1,0))*dt/4. + s*(1. - s*un)*val(g.y,0,i,0)*Delta/2.;





    if (val_fm_z(fm.z,0,i,0) && val_fm_z(fm.z,0,i,1)) {
      double vn = (val(uf.z,0,i,0) + val(uf.z,0,i,1))/(val_fm_z(fm.z,0,i,0) + val_fm_z(fm.z,0,i,1));
      double fyy = vn < 0. ? val(f,0,i,1) - val(f,0,i,0) : val(f,0,i,0) - val(f,0,i,-1);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double wn = (val(uf.x,0,i,0) + val(uf.x,1,i,0))/(val_fm_x(fm.x,0,i,0) + val_fm_x(fm.x,1,i,0));
      double fzz = wn < 0. ? val(f,1,i,0) - val(f,0,i,0) : val(f,0,i,0) - val(f,-1,i,0);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.y,0,0,0) = f2*val(uf.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 28
{

#line 28 "/home/davidtang/basilisk/src/bcg.h"
 {







    double un = dt*val(uf.z,0,0,0)/(val_fm_z(fm.z,0,0,0)*Delta + 0.), s = sign(un);
    int i = -(s + 1.)/2.;
    double f2 = val(f,0,0,i) + (val_src(src,0,0,0) + val_src(src,0,0,-1))*dt/4. + s*(1. - s*un)*val(g.z,0,0,i)*Delta/2.;





    if (val_fm_x(fm.x,0,0,i) && val_fm_x(fm.x,1,0,i)) {
      double vn = (val(uf.x,0,0,i) + val(uf.x,1,0,i))/(val_fm_x(fm.x,0,0,i) + val_fm_x(fm.x,1,0,i));
      double fyy = vn < 0. ? val(f,1,0,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,-1,0,i);
      f2 -= dt*vn*fyy/(2.*Delta);
    }


    if (val_fm_y(fm.y,0,0,i) && val_fm_y(fm.y,0,1,i)) {
      double wn = (val(uf.y,0,0,i) + val(uf.y,0,1,i))/(val_fm_y(fm.y,0,0,i) + val_fm_y(fm.y,0,1,i));
      double fzz = wn < 0. ? val(f,0,1,i) - val(f,0,0,i) : val(f,0,0,i) - val(f,0,-1,i);
      f2 -= dt*wn*fzz/(2.*Delta);
    }


    val(flux.z,0,0,0) = f2*val(uf.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 59
 end_foreach_face(); } }





  boundary_flux (((vector []){{flux.x,flux.y,flux.z},{{-1},{-1},{-1}}}));
 delete (((scalar []){g.x,g.y,g.z,{-1}})); }






struct Advection {
  scalar * tracers;
  vector u;
  double dt;
  scalar * src;
};

void advection (struct Advection p)
{




  scalar * lsrc = p.src;
  if (!lsrc) {
    scalar zero= new_const_scalar("zero", 8,  0.);
    if (p.tracers) for (scalar s = *p.tracers, *_i54 = p.tracers; ((scalar *)&s)->i >= 0; s = *++_i54)
      lsrc = list_append (lsrc, zero);
  }

  assert (list_len(p.tracers) == list_len(lsrc));
  scalar f, src;
  scalar * _i4 = p.tracers; scalar * _i5 = lsrc; if (p.tracers) for (f = *p.tracers, src = *lsrc; ((scalar *)&f)->i >= 0; f = *++_i4, src = *++_i5) {
    vector flux= new_face_vector("flux");
    tracer_fluxes (f, p.u, flux, p.dt, src);

     { 
if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 99
foreach(){

#line 99 "/home/davidtang/basilisk/src/bcg.h"

      {
#line 100

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 100

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 100

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 99
foreach(){

#line 99 "/home/davidtang/basilisk/src/bcg.h"

      {
#line 100

        val(f,0,0,0) += p.dt*(val(flux.x,0,0,0) - val(flux.x,1,0,0))/(Delta*val_cm(cm,0,0,0));
#line 100

        val(f,0,0,0) += p.dt*(val(flux.y,0,0,0) - val(flux.y,0,1,0))/(Delta*val_cm(cm,0,0,0));
#line 100

        val(f,0,0,0) += p.dt*(val(flux.z,0,0,0) - val(flux.z,0,0,1))/(Delta*val_cm(cm,0,0,0));}; } end_foreach(); } }



   delete (((scalar []){flux.x,flux.y,flux.z,{-1}})); }
  boundary (p.tracers);

  if (!p.src)
    pfree (lsrc,__func__,__FILE__,__LINE__);
}
#line 30 "/home/davidtang/basilisk/src/navier-stokes/centered.h"



#line 1 "./viscosity.h"
#line 1 "/home/davidtang/basilisk/src/viscosity.h"
#line 1 "./poisson.h"
#line 1 "/home/davidtang/basilisk/src/poisson.h"
#line 32 "/home/davidtang/basilisk/src/poisson.h"
void mg_cycle (scalar * a, scalar * res, scalar * da,
        void (* relax) (scalar * da, scalar * res,
          int depth, void * data),
        void * data,
        int nrelax, int minlevel, int maxlevel)
{




  restriction (res);





  minlevel = min (minlevel, maxlevel);
  for (int l = minlevel; l <= maxlevel; l++) {




    if (l == minlevel)
       { foreach_level_or_leaf (l){

#line 55 "/home/davidtang/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i55 = da; ((scalar *)&s)->i >= 0; s = *++_i55)
  
     val(s,0,0,0) = 0.; } end_foreach_level_or_leaf(); }





    else
       { foreach_level (l){

#line 65 "/home/davidtang/basilisk/src/poisson.h"

 if (da) for (scalar s = *da, *_i56 = da; ((scalar *)&s)->i >= 0; s = *++_i56)
  
     val(s,0,0,0) = bilinear (point, s); } end_foreach_level(); }





    boundary_level (da, l);
    for (int i = 0; i < nrelax; i++) {
      relax (da, res, l, data);
      boundary_level (da, l);
    }
  }




   { foreach(){

#line 84 "/home/davidtang/basilisk/src/poisson.h"
 {
    scalar s, ds;
    scalar * _i6 = a; scalar * _i7 = da; if (a) for (s = *a, ds = *da; ((scalar *)&s)->i >= 0; s = *++_i6, ds = *++_i7)
     
 val(s,0,0,0) += val(ds,0,0,0);
  } } end_foreach(); }
  boundary (a);
}
#line 103 "/home/davidtang/basilisk/src/poisson.h"
int NITERMAX = 100, NITERMIN = 1;
double TOLERANCE = 1e-3;




typedef struct {
  int i;
  double resb, resa;
  double sum;
  int nrelax;
  int minlevel;
} mgstats;
#line 126 "/home/davidtang/basilisk/src/poisson.h"
struct MGSolve {
  scalar * a, * b;
  double (* residual) (scalar * a, scalar * b, scalar * res,
         void * data);
  void (* relax) (scalar * da, scalar * res, int depth,
    void * data);
  void * data;

  int nrelax;
  scalar * res;
  int minlevel;
  double tolerance;
};

mgstats mg_solve (struct MGSolve p)
{





  scalar * da = list_clone (p.a), * res = p.res;
  if (!res)
    res = list_clone (p.b);






  for (int b = 0; b < nboundary; b++)
    if (da) for (scalar s = *da, *_i57 = da; ((scalar *)&s)->i >= 0; s = *++_i57)
      _attribute[s.i].boundary[b] = _attribute[s.i].boundary_homogeneous[b];




  mgstats s = {0};
  double sum = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _sum = sum; 
#line 165
foreach (){

#line 165 "/home/davidtang/basilisk/src/poisson.h"

    if (p.b) for (scalar s = *p.b, *_i58 = p.b; ((scalar *)&s)->i >= 0; s = *++_i58)
      _sum += val(s,0,0,0); } end_foreach();OMP(omp critical) sum += _sum;
mpi_all_reduce_double (sum, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 167
 }
  s.sum = sum;
  s.nrelax = p.nrelax > 0 ? p.nrelax : 4;




  double resb;
  resb = s.resb = s.resa = p.residual (p.a, p.b, res, p.data);






  if (p.tolerance == 0.)
    p.tolerance = TOLERANCE;
  for (s.i = 0;
       s.i < NITERMAX && (s.i < NITERMIN || s.resa > p.tolerance);
       s.i++) {
    mg_cycle (p.a, res, da, p.relax, p.data,
       s.nrelax,
       p.minlevel,
       grid->maxdepth);
    s.resa = p.residual (p.a, p.b, res, p.data);
#line 200 "/home/davidtang/basilisk/src/poisson.h"
    if (s.resa > p.tolerance) {
      if (resb/s.resa < 1.2 && s.nrelax < 100)
 s.nrelax++;
      else if (resb/s.resa > 10 && s.nrelax > 2)
 s.nrelax--;
    }







    resb = s.resa;
  }
  s.minlevel = p.minlevel;




  if (s.resa > p.tolerance) {
    scalar v = p.a[0];
    fprintf (ferr,
      "WARNING: convergence for %s not reached after %d iterations\n"
      "  res: %g sum: %g nrelax: %d\n", _attribute[v.i].name,
      s.i, s.resa, s.sum, s.nrelax), fflush (ferr);
  }




  if (!p.res)
    delete (res), pfree (res,__func__,__FILE__,__LINE__);
  delete (da), pfree (da,__func__,__FILE__,__LINE__);

  return s;
}
#line 263 "/home/davidtang/basilisk/src/poisson.h"
struct Poisson {
  scalar a, b;
   vector alpha;
   scalar lambda;
  double tolerance;
  int nrelax, minlevel;
  scalar * res;



};





static void relax (scalar * al, scalar * bl, int l, void * data)
{
  scalar a = al[0], b = bl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
#line 301 "/home/davidtang/basilisk/src/poisson.h"
  scalar c = a;






   { 
if (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 308
foreach_level_or_leaf (l){

#line 308 "/home/davidtang/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 310
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 310
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 310
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
#line 324 "/home/davidtang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 308
foreach_level_or_leaf (l){

#line 308 "/home/davidtang/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 310
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 310
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 310
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
#line 324 "/home/davidtang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 308
foreach_level_or_leaf (l){

#line 308 "/home/davidtang/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 310
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 310
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 310
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
#line 324 "/home/davidtang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); }
if (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 308
foreach_level_or_leaf (l){

#line 308 "/home/davidtang/basilisk/src/poisson.h"
 {
    double n = - sq(Delta)*val(b,0,0,0), d = - val_lambda(lambda,0,0,0)*sq(Delta);
    {
#line 310
 {
      n += val_alpha_x(alpha.x,1,0,0)*val(a,1,0,0) + val_alpha_x(alpha.x,0,0,0)*val(a,-1,0,0);
      d += val_alpha_x(alpha.x,1,0,0) + val_alpha_x(alpha.x,0,0,0);
    }
#line 310
 {
      n += val_alpha_y(alpha.y,0,1,0)*val(a,0,1,0) + val_alpha_y(alpha.y,0,0,0)*val(a,0,-1,0);
      d += val_alpha_y(alpha.y,0,1,0) + val_alpha_y(alpha.y,0,0,0);
    }
#line 310
 {
      n += val_alpha_z(alpha.z,0,0,1)*val(a,0,0,1) + val_alpha_z(alpha.z,0,0,0)*val(a,0,0,-1);
      d += val_alpha_z(alpha.z,0,0,1) + val_alpha_z(alpha.z,0,0,0);
    }}
#line 324 "/home/davidtang/basilisk/src/poisson.h"
      val(c,0,0,0) = n/d;
  } } end_foreach_level_or_leaf(); } }
#line 343 "/home/davidtang/basilisk/src/poisson.h"
}






static double residual (scalar * al, scalar * bl, scalar * resl, void * data)
{
  scalar a = al[0], b = bl[0], res = resl[0];
  struct Poisson * p = (struct Poisson *) data;
   vector alpha = p->alpha;
   scalar lambda = p->lambda;
  double maxres = 0.;
#line 378 "/home/davidtang/basilisk/src/poisson.h"
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 378

if (!is_constant(lambda) && !is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 378
foreach (){

#line 378 "/home/davidtang/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 380

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_z(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  val_alpha_z(alpha.z,0,0,1)*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
if (is_constant(lambda) && !is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 378
foreach (){

#line 378 "/home/davidtang/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 380

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_z(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  val_alpha_z(alpha.z,0,0,1)*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
if (!is_constant(lambda) && is_constant(alpha.x)) {
#undef val_lambda
#define val_lambda(a,i,j,k) val(a,i,j,k)
#undef fine_lambda
#define fine_lambda(a,i,j,k) fine(a,i,j,k)
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 378
foreach (){

#line 378 "/home/davidtang/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 380

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_z(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  val_alpha_z(alpha.z,0,0,1)*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }
if (is_constant(lambda) && is_constant(alpha.x)) {
const double _const_lambda = _constant[lambda.i -_NVARMAX];
NOT_UNUSED(_const_lambda);
#undef val_lambda
#define val_lambda(a,i,j,k) _const_lambda
#undef fine_lambda
#define fine_lambda(a,i,j,k) _const_lambda
#undef coarse_lambda
#define coarse_lambda(a,i,j,k) _const_lambda
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 378
foreach (){

#line 378 "/home/davidtang/basilisk/src/poisson.h"
 {
    val(res,0,0,0) = val(b,0,0,0) - val_lambda(lambda,0,0,0)*val(a,0,0,0);
    {
#line 380

      val(res,0,0,0) += (val_alpha_x(alpha.x,0,0,0)*((val(a,0,0,0) - val(a,0 -1,0,0))/Delta) -
  val_alpha_x(alpha.x,1,0,0)*((val(a,1,0,0) - val(a,1 -1,0,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_y(alpha.y,0,0,0)*((val(a,0,0,0) - val(a,0,0 -1,0))/Delta) -
  val_alpha_y(alpha.y,0,1,0)*((val(a,0,1,0) - val(a,0,1 -1,0))/Delta))/Delta;
#line 380

      val(res,0,0,0) += (val_alpha_z(alpha.z,0,0,0)*((val(a,0,0,0) - val(a,0,0,0 -1))/Delta) -
  val_alpha_z(alpha.z,0,0,1)*((val(a,0,0,1) - val(a,0,0,1 -1))/Delta))/Delta;}






    if (fabs (val(res,0,0,0)) > _maxres)
      _maxres = fabs (val(res,0,0,0));
  } } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 391
 }

  boundary (resl);
  return maxres;
}
#line 406 "/home/davidtang/basilisk/src/poisson.h"
mgstats poisson (struct Poisson p)
{






  if (!p.alpha.x.i)
    p.alpha = unityf;
  if (!p.lambda.i)
    p.lambda = zeroc;




  vector alpha = p.alpha;
  scalar lambda = p.lambda;
  restriction (((scalar []){alpha.x,alpha.y,alpha.z,lambda,{-1}}));





  double defaultol = TOLERANCE;
  if (p.tolerance)
    TOLERANCE = p.tolerance;

  scalar a = p.a, b = p.b;




  mgstats s = mg_solve ((struct MGSolve){((scalar []){a,{-1}}), ((scalar []){b,{-1}}), residual, relax,
   &p, p.nrelax, p.res, .minlevel = max(1, p.minlevel)});




  if (p.tolerance)
    TOLERANCE = defaultol;

  return s;
}
#line 468 "/home/davidtang/basilisk/src/poisson.h"
struct Project {
  vector uf;
  scalar p;
  vector alpha;
  double dt;
  int nrelax;
};


mgstats project (struct Project q)
{ trace ("project", "/home/davidtang/basilisk/src/poisson.h", 478);
  vector uf = q.uf;
  scalar p = q.p;
   vector alpha = q.alpha.x.i ? q.alpha : unityf;
  double dt = q.dt ? q.dt : 1.;
  int nrelax = q.nrelax ? q.nrelax : 4;






  scalar div= new_scalar("div");
   { foreach(){

#line 491 "/home/davidtang/basilisk/src/poisson.h"
 {
    val(div,0,0,0) = 0.;
    {
#line 493

      val(div,0,0,0) += val(uf.x,1,0,0) - val(uf.x,0,0,0);
#line 493

      val(div,0,0,0) += val(uf.y,0,1,0) - val(uf.y,0,0,0);
#line 493

      val(div,0,0,0) += val(uf.z,0,0,1) - val(uf.z,0,0,0);}
    val(div,0,0,0) /= dt*Delta;
  } } end_foreach(); }
#line 507 "/home/davidtang/basilisk/src/poisson.h"
  mgstats mgp = poisson ((struct Poisson){p, div, alpha,
    .tolerance = TOLERANCE/sq(dt), .nrelax = nrelax});




   { 
if (!is_constant(alpha.x)) {
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 513
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 513
{

#line 513 "/home/davidtang/basilisk/src/poisson.h"

    val(uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 513
{

#line 513 "/home/davidtang/basilisk/src/poisson.h"

    val(uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 513
{

#line 513 "/home/davidtang/basilisk/src/poisson.h"

    val(uf.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta); }  }}  end_foreach_face_generic()
#line 514
 end_foreach_face(); }
if (is_constant(alpha.x)) {
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 513
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 513
{

#line 513 "/home/davidtang/basilisk/src/poisson.h"

    val(uf.x,0,0,0) -= dt*val_alpha_x(alpha.x,0,0,0)*((val(p,0,0,0) - val(p,0 -1,0,0))/Delta); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 513
{

#line 513 "/home/davidtang/basilisk/src/poisson.h"

    val(uf.y,0,0,0) -= dt*val_alpha_y(alpha.y,0,0,0)*((val(p,0,0,0) - val(p,0,0 -1,0))/Delta); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 513
{

#line 513 "/home/davidtang/basilisk/src/poisson.h"

    val(uf.z,0,0,0) -= dt*val_alpha_z(alpha.z,0,0,0)*((val(p,0,0,0) - val(p,0,0,0 -1))/Delta); }  }}  end_foreach_face_generic()
#line 514
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));

  { mgstats _ret =  mgp; delete (((scalar []){div,{-1}}));  end_trace("project", "/home/davidtang/basilisk/src/poisson.h", 517);  return _ret; }
 delete (((scalar []){div,{-1}}));  end_trace("project", "/home/davidtang/basilisk/src/poisson.h", 518); }
#line 2 "/home/davidtang/basilisk/src/viscosity.h"

struct Viscosity {
  vector u;
  vector mu;
  scalar rho;
  double dt;
  int nrelax;
  scalar * res;
};
#line 25 "/home/davidtang/basilisk/src/viscosity.h"
static void relax_viscosity (scalar * a, scalar * b, int l, void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0])));




  vector w = u;


   { 
if (!is_constant(rho) && !is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/davidtang/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && !is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/davidtang/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (!is_constant(rho) && is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/davidtang/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); }
if (is_constant(rho) && is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 39
foreach_level_or_leaf (l){

#line 39 "/home/davidtang/basilisk/src/viscosity.h"
 {
    {
#line 40

      val(w.x,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*val(u.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)*val(u.x,-1,0,0)

      + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) +
     (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
     (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.x,0,-1,0) +
         (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
         (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


      + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) +
       (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
       (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.x,0,0,-1) +
         (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
         (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

      ) + val(r.x,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).x + dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0) + 2.*val_mu_x(mu.x,0,0,0)

          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)


          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)

        ));
#line 40

      val(w.y,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*val(u.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)*val(u.y,0,-1,0)

      + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) +
     (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
     (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
      - val_mu_z(mu.z,0,0,0)*(- val(u.y,0,0,-1) +
         (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
         (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


      + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) +
       (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
       (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.y,-1,0,0) +
         (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
         (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

      ) + val(r.y,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).y + dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0) + 2.*val_mu_y(mu.y,0,0,0)

          + val_mu_z(mu.z,0,0,1) + val_mu_z(mu.z,0,0,0)


          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)

        ));
#line 40

      val(w.z,0,0,0) = (dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*val(u.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)*val(u.z,0,0,-1)

      + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) +
     (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
     (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
      - val_mu_x(mu.x,0,0,0)*(- val(u.z,-1,0,0) +
         (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
         (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


      + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) +
       (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
       (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
      - val_mu_y(mu.y,0,0,0)*(- val(u.z,0,-1,0) +
         (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
         (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

      ) + val(r.z,0,0,0)*sq(Delta))/
    (sq(Delta)*((coord){1.,1.,1.}).z + dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1) + 2.*val_mu_z(mu.z,0,0,0)

          + val_mu_x(mu.x,1,0,0) + val_mu_x(mu.x,0,0,0)


          + val_mu_y(mu.y,0,1,0) + val_mu_y(mu.y,0,0,0)

        ));}
  } } end_foreach_level_or_leaf(); } }
#line 85 "/home/davidtang/basilisk/src/viscosity.h"
}

static double residual_viscosity (scalar * a, scalar * b, scalar * resl,
      void * data)
{
  struct Viscosity * p = (struct Viscosity *) data;
   vector mu = p->mu;
   scalar rho = p->rho;
  double dt = p->dt;
  vector u = (*((vector *)&(a[0]))), r = (*((vector *)&(b[0]))), res = (*((vector *)&(resl[0])));
  double maxres = 0.;
#line 127 "/home/davidtang/basilisk/src/viscosity.h"
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _maxres = maxres; 
#line 127

if (!is_constant(rho) && !is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 127
foreach (){

#line 127 "/home/davidtang/basilisk/src/viscosity.h"

    {
#line 128
 {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    }
#line 128
 {
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    }
#line 128
 {
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    }} } end_foreach(); }
if (is_constant(rho) && !is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
#undef val_mu_x
#define val_mu_x(a,i,j,k) val(a,i,j,k)
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_y
#define val_mu_y(a,i,j,k) val(a,i,j,k)
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) coarse(a,i,j,k)
#undef val_mu_z
#define val_mu_z(a,i,j,k) val(a,i,j,k)
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) coarse(a,i,j,k)
#line 127
foreach (){

#line 127 "/home/davidtang/basilisk/src/viscosity.h"

    {
#line 128
 {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    }
#line 128
 {
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    }
#line 128
 {
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    }} } end_foreach(); }
if (!is_constant(rho) && is_constant(mu.x)) {
#undef val_rho
#define val_rho(a,i,j,k) val(a,i,j,k)
#undef fine_rho
#define fine_rho(a,i,j,k) fine(a,i,j,k)
#undef coarse_rho
#define coarse_rho(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 127
foreach (){

#line 127 "/home/davidtang/basilisk/src/viscosity.h"

    {
#line 128
 {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    }
#line 128
 {
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    }
#line 128
 {
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    }} } end_foreach(); }
if (is_constant(rho) && is_constant(mu.x)) {
const double _const_rho = _constant[rho.i -_NVARMAX];
NOT_UNUSED(_const_rho);
#undef val_rho
#define val_rho(a,i,j,k) _const_rho
#undef fine_rho
#define fine_rho(a,i,j,k) _const_rho
#undef coarse_rho
#define coarse_rho(a,i,j,k) _const_rho
const struct { double x, y, z; } _const_mu = {_constant[mu.x.i -_NVARMAX], _constant[mu.y.i - _NVARMAX], _constant[mu.z.i - _NVARMAX]};
NOT_UNUSED(_const_mu);
#undef val_mu_x
#define val_mu_x(a,i,j,k) _const_mu.x
#undef fine_mu_x
#define fine_mu_x(a,i,j,k) _const_mu.x
#undef coarse_mu_x
#define coarse_mu_x(a,i,j,k) _const_mu.x
#undef val_mu_y
#define val_mu_y(a,i,j,k) _const_mu.y
#undef fine_mu_y
#define fine_mu_y(a,i,j,k) _const_mu.y
#undef coarse_mu_y
#define coarse_mu_y(a,i,j,k) _const_mu.y
#undef val_mu_z
#define val_mu_z(a,i,j,k) _const_mu.z
#undef fine_mu_z
#define fine_mu_z(a,i,j,k) _const_mu.z
#undef coarse_mu_z
#define coarse_mu_z(a,i,j,k) _const_mu.z
#line 127
foreach (){

#line 127 "/home/davidtang/basilisk/src/viscosity.h"

    {
#line 128
 {
      val(res.x,0,0,0) = val(r.x,0,0,0) - ((coord){1.,1.,1.}).x*val(u.x,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_x(mu.x,1,0,0)*(val(u.x,1,0,0) - val(u.x,0,0,0))
    - 2.*val_mu_x(mu.x,0,0,0)*(val(u.x,0,0,0) - val(u.x,-1,0,0))

    + val_mu_y(mu.y,0,1,0)*(val(u.x,0,1,0) - val(u.x,0,0,0) +
          (val(u.y,1,0,0) + val(u.y,1,1,0))/4. -
          (val(u.y,-1,0,0) + val(u.y,-1,1,0))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,-1,0) +
       (val(u.y,1,-1,0) + val(u.y,1,0,0))/4. -
       (val(u.y,-1,-1,0) + val(u.y,-1,0,0))/4.)


    + val_mu_z(mu.z,0,0,1)*(val(u.x,0,0,1) - val(u.x,0,0,0) +
     (val(u.z,1,0,0) + val(u.z,1,0,1))/4. -
     (val(u.z,-1,0,0) + val(u.z,-1,0,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.x,0,0,0) - val(u.x,0,0,-1) +
       (val(u.z,1,0,-1) + val(u.z,1,0,0))/4. -
       (val(u.z,-1,0,-1) + val(u.z,-1,0,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.x,0,0,0)) > _maxres)
 _maxres = fabs (val(res.x,0,0,0));
    }
#line 128
 {
      val(res.y,0,0,0) = val(r.y,0,0,0) - ((coord){1.,1.,1.}).y*val(u.y,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_y(mu.y,0,1,0)*(val(u.y,0,1,0) - val(u.y,0,0,0))
    - 2.*val_mu_y(mu.y,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,-1,0))

    + val_mu_z(mu.z,0,0,1)*(val(u.y,0,0,1) - val(u.y,0,0,0) +
          (val(u.z,0,1,0) + val(u.z,0,1,1))/4. -
          (val(u.z,0,-1,0) + val(u.z,0,-1,1))/4.)
    - val_mu_z(mu.z,0,0,0)*(val(u.y,0,0,0) - val(u.y,0,0,-1) +
       (val(u.z,0,1,-1) + val(u.z,0,1,0))/4. -
       (val(u.z,0,-1,-1) + val(u.z,0,-1,0))/4.)


    + val_mu_x(mu.x,1,0,0)*(val(u.y,1,0,0) - val(u.y,0,0,0) +
     (val(u.x,0,1,0) + val(u.x,1,1,0))/4. -
     (val(u.x,0,-1,0) + val(u.x,1,-1,0))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.y,0,0,0) - val(u.y,-1,0,0) +
       (val(u.x,-1,1,0) + val(u.x,0,1,0))/4. -
       (val(u.x,-1,-1,0) + val(u.x,0,-1,0))/4.)

    )/sq(Delta);
      if (fabs (val(res.y,0,0,0)) > _maxres)
 _maxres = fabs (val(res.y,0,0,0));
    }
#line 128
 {
      val(res.z,0,0,0) = val(r.z,0,0,0) - ((coord){1.,1.,1.}).z*val(u.z,0,0,0) +
        dt/val_rho(rho,0,0,0)*(2.*val_mu_z(mu.z,0,0,1)*(val(u.z,0,0,1) - val(u.z,0,0,0))
    - 2.*val_mu_z(mu.z,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,0,-1))

    + val_mu_x(mu.x,1,0,0)*(val(u.z,1,0,0) - val(u.z,0,0,0) +
          (val(u.x,0,0,1) + val(u.x,1,0,1))/4. -
          (val(u.x,0,0,-1) + val(u.x,1,0,-1))/4.)
    - val_mu_x(mu.x,0,0,0)*(val(u.z,0,0,0) - val(u.z,-1,0,0) +
       (val(u.x,-1,0,1) + val(u.x,0,0,1))/4. -
       (val(u.x,-1,0,-1) + val(u.x,0,0,-1))/4.)


    + val_mu_y(mu.y,0,1,0)*(val(u.z,0,1,0) - val(u.z,0,0,0) +
     (val(u.y,0,0,1) + val(u.y,0,1,1))/4. -
     (val(u.y,0,0,-1) + val(u.y,0,1,-1))/4.)
    - val_mu_y(mu.y,0,0,0)*(val(u.z,0,0,0) - val(u.z,0,-1,0) +
       (val(u.y,0,-1,1) + val(u.y,0,0,1))/4. -
       (val(u.y,0,-1,-1) + val(u.y,0,0,-1))/4.)

    )/sq(Delta);
      if (fabs (val(res.z,0,0,0)) > _maxres)
 _maxres = fabs (val(res.z,0,0,0));
    }} } end_foreach(); }OMP(omp critical) if (_maxres > maxres) maxres = _maxres;
mpi_all_reduce_double (maxres, MPI_MAX);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 151
 }

  return maxres;
}




mgstats viscosity (struct Viscosity p)
{ trace ("viscosity", "/home/davidtang/basilisk/src/viscosity.h", 160);
  vector u = p.u, r= new_vector("r");
   { foreach(){

#line 162 "/home/davidtang/basilisk/src/viscosity.h"

    {
#line 163

      val(r.x,0,0,0) = val(u.x,0,0,0);
#line 163

      val(r.y,0,0,0) = val(u.y,0,0,0);
#line 163

      val(r.z,0,0,0) = val(u.z,0,0,0);}; } end_foreach(); }

  vector mu = p.mu;
  scalar rho = p.rho;
  restriction (((scalar []){mu.x,mu.y,mu.z,rho,{-1}}));

  { mgstats _ret =  mg_solve ((struct MGSolve){(scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y,r.z},{{-1},{-1},{-1}}}),
     residual_viscosity, relax_viscosity, &p, p.nrelax, p.res}); delete (((scalar []){r.x,r.y,r.z,{-1}}));  end_trace("viscosity", "/home/davidtang/basilisk/src/viscosity.h", 171);  return _ret; }
 delete (((scalar []){r.x,r.y,r.z,{-1}}));  end_trace("viscosity", "/home/davidtang/basilisk/src/viscosity.h", 172); }


mgstats viscosity_explicit (struct Viscosity p)
{ trace ("viscosity_explicit", "/home/davidtang/basilisk/src/viscosity.h", 176);
  vector u = p.u, r= new_vector("r");
  mgstats mg = {0};
  mg.resb = residual_viscosity ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), (scalar *)((vector []){{r.x,r.y,r.z},{{-1},{-1},{-1}}}), &p);
   { foreach(){

#line 180 "/home/davidtang/basilisk/src/viscosity.h"

    {
#line 181

      val(u.x,0,0,0) += val(r.x,0,0,0);
#line 181

      val(u.y,0,0,0) += val(r.y,0,0,0);
#line 181

      val(u.z,0,0,0) += val(r.z,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
  { mgstats _ret =  mg; delete (((scalar []){r.x,r.y,r.z,{-1}}));  end_trace("viscosity_explicit", "/home/davidtang/basilisk/src/viscosity.h", 184);  return _ret; }
 delete (((scalar []){r.x,r.y,r.z,{-1}}));  end_trace("viscosity_explicit", "/home/davidtang/basilisk/src/viscosity.h", 185); }
#line 34 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
#line 44 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
scalar p= {0};
vector u= {{1},{2},{3}}, g= {{4},{5},{6}};
scalar pf= {7};
vector uf= {{8},{9},{10}};
#line 70 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 vector mu = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}}, a = {{_NVARMAX + 0},{_NVARMAX + 1},{_NVARMAX + 2}}, alpha = {{_NVARMAX + 3},{_NVARMAX + 4},{_NVARMAX + 5}};
 scalar rho = {(_NVARMAX + 6)};
mgstats mgp, mgpf, mgu;
bool stokes = false;
#line 91 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static void _set_boundary0 (void) { _attribute[p.i].boundary[right] = _boundary0; _attribute[p.i].boundary_homogeneous[right] = _boundary0_homogeneous; } 
static void _set_boundary1 (void) { _attribute[p.i].boundary[left] = _boundary1; _attribute[p.i].boundary_homogeneous[left] = _boundary1_homogeneous; } 
#line 101 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static void _set_boundary2 (void) { _attribute[p.i].boundary[top] = _boundary2; _attribute[p.i].boundary_homogeneous[top] = _boundary2_homogeneous; } 
static void _set_boundary3 (void) { _attribute[p.i].boundary[bottom] = _boundary3; _attribute[p.i].boundary_homogeneous[bottom] = _boundary3_homogeneous; } 


static void _set_boundary4 (void) { _attribute[p.i].boundary[front] = _boundary4; _attribute[p.i].boundary_homogeneous[front] = _boundary4_homogeneous; } 
static void _set_boundary5 (void) { _attribute[p.i].boundary[back] = _boundary5; _attribute[p.i].boundary_homogeneous[back] = _boundary5_homogeneous; } 
#line 126 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int defaults_0 (const int i, const double t, Event * _ev) { trace ("defaults_0", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 126); 
{

  CFL = 0.8;




  _attribute[p.i].nodump = _attribute[pf.i].nodump = true;




  if (alpha.x.i == unityf.x.i) {
    alpha = fm;
    rho = cm;
  }
  else if (!is_constant(alpha.x)) {
    vector alphav = alpha;
     { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 145
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 145
{

#line 145 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 145
{

#line 145 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 145
{

#line 145 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0); }  }}  end_foreach_face_generic()
#line 146
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 145
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 145
{

#line 145 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(alphav.x,0,0,0) = val_fm_x(fm.x,0,0,0); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 145
{

#line 145 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(alphav.y,0,0,0) = val_fm_y(fm.y,0,0,0); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 145
{

#line 145 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(alphav.z,0,0,0) = val_fm_z(fm.z,0,0,0); }  }}  end_foreach_face_generic()
#line 146
 end_foreach_face(); } }
    boundary ((scalar *)((vector []){{alpha.x,alpha.y,alpha.z},{{-1},{-1},{-1}}}));
  }
#line 173 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 end_trace("defaults_0", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 173); } return 0; } 





static int default_display_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int default_display (const int i, const double t, Event * _ev) { trace ("default_display", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 179); 
  display ((struct _display){"squares (color = 'u.x', spread = -1);"}); end_trace("default_display", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 180);  return 0; } 





double dtmax;

static int init_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init (const int i, const double t, Event * _ev) { trace ("init", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 188); 
{
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 192
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 192
{

#line 192 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 192
{

#line 192 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 192
{

#line 192 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.); }  }}  end_foreach_face_generic()
#line 193
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 192
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 192
{

#line 192 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 192
{

#line 192 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 192
{

#line 192 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.); }  }}  end_foreach_face_generic()
#line 193
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));




  event ("properties");





  dtmax = DT;
  event ("stability");
 end_trace("init", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 207); } return 0; } 
#line 216 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int set_dtmax (const int i, const double t, Event * _ev) { trace ("set_dtmax", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 216);  dtmax = DT; end_trace("set_dtmax", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 216);  return 0; } 

static int stability_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int stability (const int i, const double t, Event * _ev) { trace ("stability", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 218);  {
  dt = dtnext (stokes ? dtmax : timestep (uf, dtmax));
 end_trace("stability", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 220); } return 0; } 







static int vof_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int vof (const int i, const double t, Event * _ev) { trace ("vof", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 228); ; end_trace("vof", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 228);  return 0; } 
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_advection (const int i, const double t, Event * _ev) { trace ("tracer_advection", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 229); ; end_trace("tracer_advection", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 229);  return 0; } 
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int tracer_diffusion (const int i, const double t, Event * _ev) { trace ("tracer_diffusion", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 230); ; end_trace("tracer_diffusion", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 230);  return 0; } 






static int properties_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int properties (const int i, const double t, Event * _ev) { trace ("properties", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 237);  {
  boundary (((scalar []){alpha.x,alpha.y,alpha.z,mu.x,mu.y,mu.z,rho,{-1}}));
 end_trace("properties", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 239); } return 0; } 
#line 251 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
void prediction()
{
  vector du;
  {
#line 254
 {
    scalar s = new_scalar("s");
    du.x = s;
  }
#line 254
 {
    scalar s = new_scalar("s");
    du.y = s;
  }
#line 254
 {
    scalar s = new_scalar("s");
    du.z = s;
  }}

  if (_attribute[u.x.i].gradient)
     { foreach(){

#line 260 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      {
#line 261
 {





   val(du.x,0,0,0) = _attribute[u.x.i].gradient (val(u.x,-1,0,0), val(u.x,0,0,0), val(u.x,1,0,0))/Delta;
      }
#line 261
 {





   val(du.y,0,0,0) = _attribute[u.y.i].gradient (val(u.y,0,-1,0), val(u.y,0,0,0), val(u.y,0,1,0))/Delta;
      }
#line 261
 {





   val(du.z,0,0,0) = _attribute[u.z.i].gradient (val(u.z,0,0,-1), val(u.z,0,0,0), val(u.z,0,0,1))/Delta;
      }} } end_foreach(); }
  else
     { foreach(){

#line 270 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      {
#line 271
 {





   val(du.x,0,0,0) = (val(u.x,1,0,0) - val(u.x,-1,0,0))/(2.*Delta);
    }
#line 271
 {





   val(du.y,0,0,0) = (val(u.y,0,1,0) - val(u.y,0,-1,0))/(2.*Delta);
    }
#line 271
 {





   val(du.z,0,0,0) = (val(u.z,0,0,1) - val(u.z,0,0,-1))/(2.*Delta);
    }} } end_foreach(); }
  boundary ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));

  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 282
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 282
{

#line 282 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }


    if (val_fm_z(fm.z,i,0,0) && val_fm_z(fm.z,i,0,1)) {
      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);
    }

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 282
{

#line 282 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val_fm_z(fm.z,0,i,0) && val_fm_z(fm.z,0,i,1)) {
      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);
    }


    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);
    }

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 282
{

#line 282 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

    if (val_fm_x(fm.x,0,0,i) && val_fm_x(fm.x,1,0,i)) {
      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);
    }


    if (val_fm_y(fm.y,0,0,i) && val_fm_y(fm.y,0,1,i)) {
      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);
    }

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 299
 end_foreach_face(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 282
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 282
{

#line 282 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.x,0,0,0) + val(u.x,-1,0,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.x,0,0,0) = val(u.x,i,0,0) + (val(g.x,0,0,0) + val(g.x,-1,0,0))*dt/4. + s*(1. - s*un)*val(du.x,i,0,0)*Delta/2.;

    if (val_fm_y(fm.y,i,0,0) && val_fm_y(fm.y,i,1,0)) {
      double fyy = val(u.y,i,0,0) < 0. ? val(u.x,i,1,0) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,-1,0);
      val(uf.x,0,0,0) -= dt*val(u.y,i,0,0)*fyy/(2.*Delta);
    }


    if (val_fm_z(fm.z,i,0,0) && val_fm_z(fm.z,i,0,1)) {
      double fzz = val(u.z,i,0,0) < 0. ? val(u.x,i,0,1) - val(u.x,i,0,0) : val(u.x,i,0,0) - val(u.x,i,0,-1);
      val(uf.x,0,0,0) -= dt*val(u.z,i,0,0)*fzz/(2.*Delta);
    }

    val(uf.x,0,0,0) *= val_fm_x(fm.x,0,0,0);
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 282
{

#line 282 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.y,0,0,0) + val(u.y,0,-1,0))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.y,0,0,0) = val(u.y,0,i,0) + (val(g.y,0,0,0) + val(g.y,0,-1,0))*dt/4. + s*(1. - s*un)*val(du.y,0,i,0)*Delta/2.;

    if (val_fm_z(fm.z,0,i,0) && val_fm_z(fm.z,0,i,1)) {
      double fyy = val(u.z,0,i,0) < 0. ? val(u.y,0,i,1) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,0,i,-1);
      val(uf.y,0,0,0) -= dt*val(u.z,0,i,0)*fyy/(2.*Delta);
    }


    if (val_fm_x(fm.x,0,i,0) && val_fm_x(fm.x,1,i,0)) {
      double fzz = val(u.x,0,i,0) < 0. ? val(u.y,1,i,0) - val(u.y,0,i,0) : val(u.y,0,i,0) - val(u.y,-1,i,0);
      val(uf.y,0,0,0) -= dt*val(u.x,0,i,0)*fzz/(2.*Delta);
    }

    val(uf.y,0,0,0) *= val_fm_y(fm.y,0,0,0);
  } }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 282
{

#line 282 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
 {
    double un = dt*(val(u.z,0,0,0) + val(u.z,0,0,-1))/(2.*Delta), s = sign(un);
    int i = -(s + 1.)/2.;
    val(uf.z,0,0,0) = val(u.z,0,0,i) + (val(g.z,0,0,0) + val(g.z,0,0,-1))*dt/4. + s*(1. - s*un)*val(du.z,0,0,i)*Delta/2.;

    if (val_fm_x(fm.x,0,0,i) && val_fm_x(fm.x,1,0,i)) {
      double fyy = val(u.x,0,0,i) < 0. ? val(u.z,1,0,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,-1,0,i);
      val(uf.z,0,0,0) -= dt*val(u.x,0,0,i)*fyy/(2.*Delta);
    }


    if (val_fm_y(fm.y,0,0,i) && val_fm_y(fm.y,0,1,i)) {
      double fzz = val(u.y,0,0,i) < 0. ? val(u.z,0,1,i) - val(u.z,0,0,i) : val(u.z,0,0,i) - val(u.z,0,-1,i);
      val(uf.z,0,0,0) -= dt*val(u.y,0,0,i)*fzz/(2.*Delta);
    }

    val(uf.z,0,0,0) *= val_fm_z(fm.z,0,0,0);
  } }  }}  end_foreach_face_generic()
#line 299
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));

  delete ((scalar *)((vector []){{du.x,du.y,du.z},{{-1},{-1},{-1}}}));
}
#line 314 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static int advection_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int advection_term (const int i, const double t, Event * _ev) { trace ("advection_term", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 314); 
{
  if (!stokes) {
    prediction();
    mgpf = project ((struct Project){uf, pf, alpha, dt/2., mgpf.nrelax});
    advection ((struct Advection){(scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}), uf, dt, (scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}})});
  }
 end_trace("advection_term", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 321); } return 0; } 







static void correction (double dt)
{
   { foreach(){

#line 331 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    {
#line 332

      val(u.x,0,0,0) += dt*val(g.x,0,0,0);
#line 332

      val(u.y,0,0,0) += dt*val(g.y,0,0,0);
#line 332

      val(u.z,0,0,0) += dt*val(g.z,0,0,0);}; } end_foreach(); }
  boundary ((scalar *)((vector []){{u.x,u.y,u.z},{{-1},{-1},{-1}}}));
}
#line 344 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int viscous_term (const int i, const double t, Event * _ev) { trace ("viscous_term", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 344); 
{
  if (constant(mu.x) != 0.) {
    correction (dt);
    mgu = viscosity ((struct Viscosity){u, mu, rho, dt, mgu.nrelax});
    correction (-dt);
  }




  if (!is_constant(a.x)) {
    vector af = a;
    trash (((vector []){{af.x,af.y,af.z},{{-1},{-1},{-1}}}));
     { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 358
{

#line 358 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(af.x,0,0,0) = 0.; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 358
{

#line 358 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(af.y,0,0,0) = 0.; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 358
{

#line 358 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

      val(af.z,0,0,0) = 0.; }  }}  end_foreach_face_generic()
#line 359
 end_foreach_face(); }
  }
 end_trace("viscous_term", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 361); } return 0; } 
#line 380 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static int acceleration_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration (const int i, const double t, Event * _ev) { trace ("acceleration", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 380); 
{
  trash (((vector []){{uf.x,uf.y,uf.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x) && !is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 383
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val_a_z(a.z,0,0,0)); }  }}  end_foreach_face_generic()
#line 384
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#line 383
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val_a_z(a.z,0,0,0)); }  }}  end_foreach_face_generic()
#line 384
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(a.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 383
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val_a_z(a.z,0,0,0)); }  }}  end_foreach_face_generic()
#line 384
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(a.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#line 383
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*(((val(u.x,0,0,0) + val(u.x,0 -1,0,0))/2.) + dt*val_a_x(a.x,0,0,0)); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*(((val(u.y,0,0,0) + val(u.y,0,0 -1,0))/2.) + dt*val_a_y(a.y,0,0,0)); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 383
{

#line 383 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(uf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*(((val(u.z,0,0,0) + val(u.z,0,0,0 -1))/2.) + dt*val_a_z(a.z,0,0,0)); }  }}  end_foreach_face_generic()
#line 384
 end_foreach_face(); } }
  boundary ((scalar *)((vector []){{uf.x,uf.y,uf.z},{a.x,a.y,a.z},{{-1},{-1},{-1}}}));
 end_trace("acceleration", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 386); } return 0; } 
#line 395 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
void centered_gradient (scalar p, vector g)
{





  vector gf= new_face_vector("gf");
   { 
if (!is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(a.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(a.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); }
if (!is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); }
if (is_constant(fm.x) && !is_constant(a.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); }
if (!is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); }
if (is_constant(fm.x) && is_constant(a.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 403
foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.x,0,0,0) = val_fm_x(fm.x,0,0,0)*val_a_x(a.x,0,0,0) - val_alpha_x(alpha.x,0,0,0)*(val(p,0,0,0) - val(p,-1,0,0))/Delta; }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.y,0,0,0) = val_fm_y(fm.y,0,0,0)*val_a_y(a.y,0,0,0) - val_alpha_y(alpha.y,0,0,0)*(val(p,0,0,0) - val(p,0,-1,0))/Delta; }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 403
{

#line 403 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    val(gf.z,0,0,0) = val_fm_z(fm.z,0,0,0)*val_a_z(a.z,0,0,0) - val_alpha_z(alpha.z,0,0,0)*(val(p,0,0,0) - val(p,0,0,-1))/Delta; }  }}  end_foreach_face_generic()
#line 404
 end_foreach_face(); } }
  boundary_flux (((vector []){{gf.x,gf.y,gf.z},{{-1},{-1},{-1}}}));





  trash (((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
   { 
if (!is_constant(fm.x)) {
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#line 412
foreach(){

#line 412 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    {
#line 413

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0) + 0.);
#line 413

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0) + 0.);
#line 413

      val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(val_fm_z(fm.z,0,0,0) + val_fm_z(fm.z,0,0,1) + 0.);}; } end_foreach(); }
if (is_constant(fm.x)) {
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#line 412
foreach(){

#line 412 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

    {
#line 413

      val(g.x,0,0,0) = (val(gf.x,0,0,0) + val(gf.x,1,0,0))/(val_fm_x(fm.x,0,0,0) + val_fm_x(fm.x,1,0,0) + 0.);
#line 413

      val(g.y,0,0,0) = (val(gf.y,0,0,0) + val(gf.y,0,1,0))/(val_fm_y(fm.y,0,0,0) + val_fm_y(fm.y,0,1,0) + 0.);
#line 413

      val(g.z,0,0,0) = (val(gf.z,0,0,0) + val(gf.z,0,0,1))/(val_fm_z(fm.z,0,0,0) + val_fm_z(fm.z,0,0,1) + 0.);}; } end_foreach(); } }
  boundary ((scalar *)((vector []){{g.x,g.y,g.z},{{-1},{-1},{-1}}}));
 delete (((scalar []){gf.x,gf.y,gf.z,{-1}})); }






static int projection_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int projection (const int i, const double t, Event * _ev) { trace ("projection", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 423); 
{
  mgp = project ((struct Project){uf, p, alpha, dt, mgp.nrelax});
  centered_gradient (p, g);




  correction (dt);
 end_trace("projection", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 432); } return 0; } 





static int end_timestep_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int end_timestep (const int i, const double t, Event * _ev) { trace ("end_timestep", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 438); ; end_trace("end_timestep", "/home/davidtang/basilisk/src/navier-stokes/centered.h", 438);  return 0; } 
#line 15 "isotropic.c"





#line 1 "lambda2.h"
#line 1 "/home/davidtang/basilisk/src/lambda2.h"
static void eigsrt (double d[3],
      double v[3][3])
{
  int k, j, i;
  double p;

  for (i = 0; i < 3 - 1; i++) {
    p = d[k = i];

    for (j = i + 1; j < 3; j++)
      if (d[j] >= p)
 p = d[k = j];
    if (k != i) {
      d[k] = d[i];
      d[i] = p;
      for (j = 0; j < 3; j++) {
 p = v[j][i];
 v[j][i] = v[j][k];
 v[j][k] = p;
      }
    }
  }
}
#line 37 "/home/davidtang/basilisk/src/lambda2.h"
void eigenvalues (double a[3][3],
    double d[3],
    double v[3][3])
{
  int j, iq, ip, i;
  double tresh, theta, tau, t, sm, s, h, g, c, b[3], z[3];

  for (ip = 0; ip < 3; ip++) {
    for (iq = 0; iq < 3; iq++)
      v[ip][iq] = 0.0;
    v[ip][ip] = 1.0;
  }

  for (ip = 0; ip < 3; ip++) {
    b[ip] = d[ip] = a[ip][ip];
    z[ip] = 0.0;
  }

  for (i = 1; i <= 50; i++) {
    sm = 0.0;
    for (ip = 0; ip < 3 - 1; ip++) {
      for (iq = ip + 1; iq < 3; iq++)
 sm += fabs (a[ip][iq]);
    }
    if (sm == 0.0) {
      eigsrt (d, v);
      return;
    }
    if (i < 4)
      tresh = 0.2*sm/(3*3);
    else
      tresh = 0.0;
    for (ip = 0; ip < 3 - 1; ip++) {
      for (iq = ip + 1; iq < 3; iq++) {
 g = 100.0*fabs (a[ip][iq]);
 if (i > 4 && fabs(d[ip]) + g == fabs(d[ip]) &&
     fabs(d[iq]) + g == fabs(d[iq]))
   a[ip][iq] = 0.0;
 else if (fabs (a[ip][iq]) > tresh) {
   h = d[iq] - d[ip];
   if (fabs(h) + g == fabs(h))
     t = a[ip][iq]/h;
   else {
     theta = 0.5*h/a[ip][iq];
     t = 1.0/(fabs (theta) + sqrt (1.0 + theta*theta));
     if (theta < 0.0) t = -t;
   }
   c = 1.0/sqrt (1 + t*t);
   s = t*c;
   tau = s/(1.0 + c);
   h = t*a[ip][iq];
   z[ip] -= h;
   z[iq] += h;
   d[ip] -= h;
   d[iq] += h;
   a[ip][iq] = 0.0;
   for (j = 0; j <= ip - 1; j++)
     { g=a[j][ip];h=a[j][iq];a[j][ip]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};
   for (j = ip + 1; j <= iq - 1; j++)
     { g=a[ip][j];h=a[j][iq];a[ip][j]=g-s*(h+g*tau);a[j][iq]=h+s*(g-h*tau);};
   for (j = iq + 1; j < 3; j++)
     { g=a[ip][j];h=a[iq][j];a[ip][j]=g-s*(h+g*tau);a[iq][j]=h+s*(g-h*tau);};
   for (j = 0; j < 3; j++)
     { g=v[j][ip];h=v[j][iq];v[j][ip]=g-s*(h+g*tau);v[j][iq]=h+s*(g-h*tau);};
 }
      }
    }
    for (ip = 0; ip < 3; ip++) {
      b[ip] += z[ip];
      d[ip] = b[ip];
      z[ip] = 0.0;
    }
  }

  for (i = 0; i < 3; i++) {
    for (j = 0; j < 3; j++)
      fprintf (ferr, "%10.3g ", a[i][j]);
    fprintf (ferr, "\n");
  }
  assert (false);
}

void lambda2 (const vector u, scalar l2)
{
   { foreach(){

#line 121 "/home/davidtang/basilisk/src/lambda2.h"
 {
    double JJ[3][3];
    scalar s = u.x;
    int i = 0;
    {
#line 125

      JJ[0][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));
#line 125

      JJ[0][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));
#line 125

      JJ[0][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));}
    s = u.y; i = 0;
    {
#line 128

      JJ[1][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));
#line 128

      JJ[1][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));
#line 128

      JJ[1][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));}
    s = u.z; i = 0;
    {
#line 131

      JJ[2][i++] = ((val(s,1,0,0) - val(s,-1,0,0))/(2.*Delta));
#line 131

      JJ[2][i++] = ((val(s,0,1,0) - val(s,0,-1,0))/(2.*Delta));
#line 131

      JJ[2][i++] = ((val(s,0,0,1) - val(s,0,0,-1))/(2.*Delta));}
    double S2O2[3][3];
    for (int i = 0; i < 3; i++)
      for (int j = 0; j < 3; j++) {
 S2O2[i][j] = 0.;
 for (int k = 0; k < 3; k++)
   S2O2[i][j] += JJ[i][k]*JJ[k][j] + JJ[k][i]*JJ[j][k];
      }
    double lambda[3], ev[3][3];
    eigenvalues (S2O2, lambda, ev);
    val(l2,0,0,0) = lambda[1]/2.;
  } } end_foreach(); }
  boundary (((scalar []){l2,{-1}}));
}
#line 21 "isotropic.c"
#line 1 "view.h"
#line 1 "/home/davidtang/basilisk/src/view.h"
#line 67 "/home/davidtang/basilisk/src/view.h"
#include <gl/framebuffer.h>
#include <gl/trackball.h>
#include <gl/utils.h>


#line 1 "utils.h"
#line 73 "/home/davidtang/basilisk/src/view.h"
#line 1 "input.h"
#line 1 "/home/davidtang/basilisk/src/input.h"
#line 16 "/home/davidtang/basilisk/src/input.h"
struct InputPGM {

  scalar s;
  FILE * fp;

  double ox, oy, width;
};

void input_pgm (struct InputPGM p)
{
  scalar s = p.s;
  if (p.width == 0.) p.width = L0;

  char line[81];
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read magic number\n");
    exit (1);
  }
  if (strcmp (line, "P2\n") && strcmp (line, "P5\n")) {
    fprintf (ferr, "input_pgm: magic number '%s' does not match PGM\n",
      line);
    exit (1);
  }
  int binary = !strcmp (line, "P5\n");
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  int width, height;
  while (line[0] == '#' && fgets (line, 81, p.fp));
  if (line[0] == '#' || sscanf (line, "%d %d", &width, &height) != 2) {
    fprintf (ferr, "input_pgm: could not read width and height\n");
    exit (1);
  }
  if (!fgets (line, 81, p.fp)) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  int maxval;
  if (sscanf (line, "%d", &maxval) != 1) {
    fprintf (ferr, "input_pgm: could not read maxval\n");
    exit (1);
  }
  if (maxval < 256) {
    unsigned char * a = ((unsigned char *) pmalloc ((width*height)*sizeof(unsigned char),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 1, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { foreach(){

#line 73 "/home/davidtang/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
  else {
    unsigned short * a = ((unsigned short *) pmalloc ((width*height)*sizeof(unsigned short),__func__,__FILE__,__LINE__));
    size_t n = 0;
    if (binary)
      n = fread (a, 2, width*height, p.fp);
    else {
      int v;
      while (n < width*height && fscanf (p.fp, "%d ", &v) == 1)
 a[n++] = v;
    }
    if (n != width*height) {
      fprintf (ferr, "input_pgm: read only %ld values\n", n);
      exit (1);
    }
     { foreach(){

#line 96 "/home/davidtang/basilisk/src/input.h"
 {
      int i = (x - p.ox)*width/p.width, j = (y - p.oy)*width/p.width;
      if (i >= 0 && i < width && j >= 0 && j < height)
 val(s,0,0,0) = 1. - a[(height - 1 - j)*width + i]/(double)maxval;
      else
 val(s,0,0,0) = 0.;
    } } end_foreach(); }
    pfree (a,__func__,__FILE__,__LINE__);
  }
}

static void next_char (FILE * fp, int target)
{
  int c = fgetc(fp), para = 0;
  while (c != EOF && (c != target || para > 0)) {
    if (c == '{') para++;
    if (c == '}') para--;
    c = fgetc(fp);
  }
  if (c != target) {
    fprintf (ferr, "input_gfs(): error: expecting '%c'\n", target);
    exit (1);
  }
}

static int next_string (FILE * fp, const char * target)
{
  int slen = strlen (target), para = 0;
  char s[slen + 1];
  s[slen] = '\0';
  int len = 0, c = fgetc (fp);
  while (c != EOF && len < slen) {
    if (c == '{') para++;
    if (c == '}') para--;
    s[len++] = c;
    c = fgetc (fp);
  }
  while (c != EOF && para >= 0) {
    if (!strcmp (s, target) && para == 0)
      break;
    if (c == '{') para++;
    if (c == '}') para--;
    for (int i = 0; i < slen - 1; i++)
      s[i] = s[i+1];
    s[slen - 1] = c;
    c = fgetc (fp);
  }
  if (strcmp (s, target))
    c = -1;
  return c;
}
#line 166 "/home/davidtang/basilisk/src/input.h"

void input_gfs (struct OutputGfs p)
{ trace ("input_gfs", "/home/davidtang/basilisk/src/input.h", 168);
  not_mpi_compatible();

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }
  bool input_all = (p.list == all);
  if (p.list == NULL) p.list = all;





  next_char (p.fp, '{');

  char * s = ((char *) pmalloc ((1)*sizeof(char),__func__,__FILE__,__LINE__));
  int len = 0;
  int c = fgetc(p.fp);
  while (c != EOF && c != '}') {
    s[len++] = c;
    s = (char *) prealloc (s, (len + 1)*sizeof(char),__func__,__FILE__,__LINE__);
    s[len] = '\0';
    c = fgetc(p.fp);
  }
  if (c != '}') {
    fprintf (ferr, "input_gfs(): error: expecting '}'\n");
    exit (1);
  }

  char * s1 = strstr (s, "variables");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting 'variables'\n");
    exit (1);
  }

  s1 = strstr (s1, "=");
  if (!s1) {
    fprintf (ferr, "input_gfs(): error: expecting '='\n");
    exit (1);
  }
  s1++;

  while (strchr (" \t", *s1))
    s1++;

  scalar * input = NULL;
  s1 = strtok (s1, ", \t");
  while (s1) {
    char * name = replace (s1, '_', '.', false);
    bool found = false;
    if (p.list) for (scalar s = *p.list, *_i59 = p.list; ((scalar *)&s)->i >= 0; s = *++_i59)
      if (!is_constant(s) && _attribute[s.i].name && !strcmp (_attribute[s.i].name, name)) {
 input = list_append (input, s);
 found = true; break;
      }
    if (!found) {
      if (input_all) {
 scalar s = new_scalar("s");
 pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
 _attribute[s.i].name = pstrdup (name,__func__,__FILE__,__LINE__);
 input = list_append (input, s);
      }
      else
 input = list_append (input, (scalar){INT_MAX});
    }
    pfree (name,__func__,__FILE__,__LINE__);
    s1 = strtok (NULL, ", \t");
  }
  pfree (s,__func__,__FILE__,__LINE__);

  next_char (p.fp, '{');
  double t1 = 0.;
  if (next_string (p.fp, "Time") >= 0) {
    next_char (p.fp, '{');
    next_char (p.fp, 't');
    next_char (p.fp, '=');
    if (fscanf (p.fp, "%lf", &t1) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 't'\n");
      exit (1);
    }
    next_char (p.fp, '}');
    next_char (p.fp, '}');
  }

  if (next_string (p.fp, "Box") < 0) {
    fprintf (ferr, "input_gfs(): error: expecting 'GfsBox'\n");
    exit (1);
  }

  next_char (p.fp, '{');
  next_char (p.fp, '{');
  next_char (p.fp, '\n');

  scalar * listm = ((scalar []){cm,fm.x,fm.y,fm.z,{-1}});
  scalar * listr = !is_constant(cm) ? listm : NULL;
  NOT_UNUSED (listr);

   { foreach_cell(){

#line 273 "/home/davidtang/basilisk/src/input.h"
 {
    unsigned flags;
    if (fread (&flags, sizeof (unsigned), 1, p.fp) != 1) {
      fprintf (ferr, "input_gfs(): error: expecting 'flags'\n");
      exit (1);
    }
    if (!(flags & (1 << 4)) && is_leaf(cell))
      refine_cell (point, listr, 0, NULL);
    double a;
    if (fread (&a, sizeof (double), 1, p.fp) != 1 || a != -1) {
      fprintf (ferr, "input_gfs(): error: expecting '-1'\n");
      exit (1);
    }
    if (input) for (scalar s = *input, *_i60 = input; ((scalar *)&s)->i >= 0; s = *++_i60) {
      if (fread (&a, sizeof (double), 1, p.fp) != 1) {
 fprintf (ferr, "input_gfs(): error: expecting a scalar\n");
 exit (1);
      }
      if (s.i != INT_MAX) {
 if (_attribute[s.i].v.x.i >= 0) {



   if (_attribute[s.i].v.x.i == s.i) {
     s = _attribute[s.i].v.y;
     val(s,0,0,0) = a;
   }
   else if (_attribute[s.i].v.y.i == s.i) {
     s = _attribute[s.i].v.x;
     val(s,0,0,0) = - a;
   }


   else
     val(s,0,0,0) = a;

 }
 else
   val(s,0,0,0) = a;
      }
    }
    if (is_leaf(cell))
      continue;
  } } end_foreach_cell(); }
  boundary (listm);
  boundary (input);

  pfree (input,__func__,__FILE__,__LINE__);
  if (opened)
    fclose (p.fp);


  while (t < t1 && events (false))
    t = tnext;
  events (false);
 end_trace("input_gfs", "/home/davidtang/basilisk/src/input.h", 328); }
#line 367 "/home/davidtang/basilisk/src/input.h"
struct InputGRD {
  scalar s;
  FILE * fp;
  char * file;
  double nodatavalue;
  bool linear, periodic, zero;
  int smooth;
};

void input_grd (struct InputGRD p)
{
  scalar input = p.s;

  bool opened = false;
  if (p.fp == NULL) {
    if (p.file == NULL)
      p.fp = stdin;
    else if (!(p.fp = fopen (p.file, "r"))) {
      perror (p.file);
      exit (1);
    }
    else
      opened = true;
  }


  double DeltaGRD;
  int nx, ny;
  double XG0, YG0, ndv;


  char waste[100];
  fscanf (p.fp, "%s %d", waste, &nx);
  fscanf (p.fp, "%s %d", waste, &ny);
  fscanf (p.fp, "%s %lf", waste, &XG0);
  fscanf (p.fp, "%s %lf", waste, &YG0);
  fscanf (p.fp, "%s %lf", waste, &DeltaGRD);
  fscanf (p.fp, "%s %lf", waste, &ndv);


  if (!p.nodatavalue)
    p.nodatavalue = ndv;


  double * value = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
  for (int i = ny - 1; i >= 0; i--)
    for (int j = 0 ; j < nx; j++) {
      fscanf (p.fp, "%lf ", &value[j + i*nx]);
      if (p.zero && value[j + i*nx] == ndv)
 value[j + i*nx] = 0.;
    }


  if (p.smooth > 0) {
    double * smoothed = ((double *) pmalloc ((nx*ny)*sizeof(double),__func__,__FILE__,__LINE__));
    for (int s = 0; s < p.smooth; s++) {
      for (int i = 0; i < ny; i++)
 for (int j = 0 ; j < nx; j++) {
   int n = 0;
   smoothed[j + i*nx] = 0.;
   for (int k = -1; k <= 1; k++)
     for (int l = -1; l <= 1; l++)
       if ((l != 0 || k != 0) &&
    i + k >= 0 && i + k < ny &&
    j + l >= 0 && j + l < nx &&
    value[j + l + (i + k)*nx] != ndv)
  smoothed[j + i*nx] += value[j + l + (i + k)*nx], n++;
   if (n == 0)
     smoothed[j + i*nx] = p.zero ? 0. : ndv;
   else
     smoothed[j + i*nx] /= n;
 }
      swap (double *, value, smoothed);
    }
    pfree (smoothed,__func__,__FILE__,__LINE__);
  }

  bool warning = false;
   { 
#if _OPENMP
  #undef OMP
  #define OMP(x)
#endif
#line 445
foreach (){

#line 445 "/home/davidtang/basilisk/src/input.h"
 {
    if (p.periodic || _attribute[input.i].boundary[right] == periodic_bc) {
      if (x > XG0 + nx*DeltaGRD)
 x -= nx*DeltaGRD;
      else if (x < XG0)
 x += nx*DeltaGRD;
    }

    int j = (x - XG0 + DeltaGRD/2.)/DeltaGRD;
    int i = (y - YG0 + DeltaGRD/2.)/DeltaGRD;
    if (i >= 0 && i < ny && j >= 0 && j < nx) {
      double val;

      int j1 = (x - XG0)/DeltaGRD;
      int i1 = (y - YG0)/DeltaGRD;
      if (p.linear && i1 >= 0 && j1 >= 0 && i1 < ny - 1 && j1 < nx - 1 &&
   value[j1 + i1*nx] != ndv && value[j1 + 1 + i1*nx] != ndv &&
   value[j1 + (i1 + 1)*nx] != ndv && value[j1 + 1 + (i1 + 1)*nx] != ndv) {

 double dx = x - (j1*DeltaGRD + XG0);
 double dy = y - (i1*DeltaGRD + YG0);
 val = (value[j1 + i1*nx] +
        dx*(value[j1 + 1 + i1*nx] - value[j1 + i1*nx])/DeltaGRD +
        dy*(value[j1 + (i1 + 1)*nx] - value[j1 + i1*nx])/DeltaGRD +
        dx*dy*(value[j1 + i1*nx] + value[j1 + 1 + (i1 + 1)*nx] -
        value[j1 + (i1 + 1)*nx] - value[j1 + 1 + i1*nx])
        /sq(DeltaGRD));
      }
      else
 val = value[j + i*nx];
      if (val == ndv)
 val(input,0,0,0) = nodata;
      else
 val(input,0,0,0) = val;
    }
    else {
      val(input,0,0,0) = nodata;
      warning = true;
    }
  } } end_foreach();
#if _OPENMP
  #undef OMP
  #define OMP(x) _Pragma(#x)
#endif
#line 484
 }
  pfree (value,__func__,__FILE__,__LINE__);

  if (warning)
    fprintf (ferr,
      "input_grd(): Warning: Raster data is not covering all"
      " the simulation area\n");

  if (opened)
    fclose (p.fp);
}
#line 74 "/home/davidtang/basilisk/src/view.h"







typedef struct {
  char * expr;
  scalar s;
} cexpr;

static scalar get_cexpr (cexpr * cache, const char * expr)
{
  cexpr * c = cache;
  while (c->expr) {
    if (!strcmp (c->expr, expr)) {


      cexpr tmp = *c;
      while ((c + 1)->expr)
 *c = *(c + 1), c++;
      *c = tmp;
      return c->s;
    }
    c++;
  }
  return (scalar){-1};
}

static cexpr * add_cexpr (cexpr * cache, int maxlen,
     const char * expr, scalar s)
{
  cexpr * c = cache;
  while (c->expr) c++;
  int len = c - cache;
  if (len < maxlen) {
    cache = prealloc (cache, sizeof(cexpr)*(len + 2),__func__,__FILE__,__LINE__);
    c = &cache[len];
  }
  else {

    c = cache;
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar []){s,{-1}}));

    while ((c + 1)->expr)
      *c = *(c + 1), c++;
  }
  c->expr = pstrdup (expr,__func__,__FILE__,__LINE__);
  c->s = s;
  (c + 1)->expr = NULL;
  return cache;
}

static void free_cexpr (cexpr * cache)
{
  cexpr * c = cache;
  while (c->expr) {
    pfree (c->expr,__func__,__FILE__,__LINE__);
    scalar s = c->s;
    delete (((scalar []){s,{-1}}));
    c++;
  }
  pfree (cache,__func__,__FILE__,__LINE__);
}






struct _bview {
  float tx, ty, sx, sy, sz;
  float quat[4];
  float fov;
  float tz, near, far;

  bool gfsview;
  bool reversed;

  float bg[3];
  float lc;
  float res;

  unsigned width, height, samples;

  framebuffer * fb;
  Frustum frustum;

  void (* map) (coord *);

  int ni;

  bool active;

  cexpr * cache;
  int maxlen;
};

typedef struct _bview bview;




bview * bview_new()
{
  bview * p = ((bview *) pcalloc (1, sizeof(bview),__func__,__FILE__,__LINE__));

  p->tx = p->ty = 0;
  p->sx = p->sy = p->sz = 1.;
  p->quat[0] = p->quat[1] = p->quat[2] = 0; p->quat[3] = 1;
  p->fov = 24.;
  gl_trackball (p->quat, 0.0, 0.0, 0.0, 0.0);




  p->bg[0] = 0.3; p->bg[1] = 0.4; p->bg[2] = 0.6;

  p->res = 1.;
  p->lc = 0.001;

  p->samples = 4;
  p->width = 600*p->samples, p->height = 600*p->samples;


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  p->fb = framebuffer_new (p->width, p->height);

  init_gl();
  p->active = false;

  enable_fpe (FE_DIVBYZERO|FE_INVALID);

  return p;
}




void bview_destroy (bview * p)
{
  framebuffer_destroy (p->fb);
  if (p->cache)
    free_cexpr (p->cache);
  pfree (p,__func__,__FILE__,__LINE__);
}




static bview * _view = NULL;






static void destroy_view()
{
  assert (_view);
  bview_destroy (_view);
}

bview * get_view() {
  if (!_view) {
    _view = bview_new();
    free_solver_func_add (destroy_view);
  }
  return _view;
}




static void redraw() {
  bview * view = get_view();


  disable_fpe (FE_DIVBYZERO|FE_INVALID);

  glMatrixMode (GL_PROJECTION);
  glLoadIdentity ();

  if (view->far <= view->near) {
    double max = 2.;
    gluPerspective (view->fov, view->width/(float)view->height, 1., 1. + 2.*max);

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, - (1. + max));
  }
  else {
    gluPerspective (view->fov, view->width/(float)view->height,
      view->near, view->far);

    glMatrixMode (GL_MODELVIEW);
    glLoadIdentity ();
    glTranslatef (view->tx, view->ty, view->tz);
  }

  GLfloat m[4][4];
  gl_build_rotmatrix (m, view->quat);
  glMultMatrixf (&m[0][0]);

  if (view->gfsview) {
    m[0][0] = 0., m[0][1] = 0., m[0][2] = -1.;
    m[1][0] = 0., m[1][1] = -1., m[1][2] = 0.;
    m[2][0] = 1., m[2][1] = 0., m[2][2] = 0.;
    glMultMatrixf (&m[0][0]);
  }

  glScalef (view->sx/L0, view->sy/L0, view->sz/L0);

  glClearColor (view->bg[0], view->bg[1], view->bg[2], 0.);
  glClear (GL_COLOR_BUFFER_BIT|GL_DEPTH_BUFFER_BIT);

  gl_get_frustum (&view->frustum);

  view->active = true;
  view->ni = 0;
}




bview * draw() {
  bview * view = get_view();
  if (!view->active)
    redraw();
  else


    disable_fpe (FE_DIVBYZERO|FE_INVALID);
  return view;
}







typedef void * pointer;
#line 363 "/home/davidtang/basilisk/src/view.h"
typedef struct {
  GLubyte a[4];
  float depth;
} RGBA;

static void compose_image_op (void * pin, void * pout, int * len,
         MPI_Datatype * dptr)
{
  RGBA * in = pin, * out = pout;
  for (int i = 0; i < *len; i++,in++,out++)
    if (out->depth > in->depth)
      *out = *in;
}


static pointer compose_image (bview * view)
{ trace ("compose_image", "/home/davidtang/basilisk/src/view.h", 379);
  unsigned char * image = framebuffer_image (view->fb);
  if (npe() > 1) {
    MPI_Op op;
    MPI_Op_create (compose_image_op, true, &op);
    MPI_Datatype rgba;
    MPI_Type_create_struct (2,
       (int[]){4,1},
       (MPI_Aint[]){0,4},
       (MPI_Datatype[]){MPI_BYTE, MPI_FLOAT},
       &rgba);
    MPI_Type_commit (&rgba);
    fbdepth_t * depth = framebuffer_depth (view->fb);
    int size = view->width*view->height;
    RGBA * buf = pmalloc (size*sizeof(RGBA),__func__,__FILE__,__LINE__);
    unsigned char * ptr = image;
    fbdepth_t * dptr = depth;
    for (int i = 0; i < size; i++) {
      for (int j = 0; j < 4; j++)
 buf[i].a[j] = *ptr++;
      buf[i].depth = *dptr++;
    }
    if (pid() == 0) {
      MPI_Reduce (MPI_IN_PLACE, buf, size, rgba, op, 0, MPI_COMM_WORLD);
      unsigned char * ptr = image;
      for (int i = 0; i < size; i++)
 for (int j = 0; j < 4; j++)
   *ptr++ = buf[i].a[j];
    }
    else
      MPI_Reduce (buf, buf, size, rgba, op, 0, MPI_COMM_WORLD);
    pfree (buf,__func__,__FILE__,__LINE__);
    MPI_Op_free (&op);
    MPI_Type_free (&rgba);
  }
  { pointer _ret =  image; end_trace("compose_image", "/home/davidtang/basilisk/src/view.h", 414);  return _ret; }
 end_trace("compose_image", "/home/davidtang/basilisk/src/view.h", 415); }



#line 1 "vertexbuffer.h"
#line 1 "/home/davidtang/basilisk/src/vertexbuffer.h"
#line 14 "/home/davidtang/basilisk/src/vertexbuffer.h"
struct {

  Array * position, * normal, * color, * index;
  float modelview[16];
  int type;
  int dim;
  int vertex, nvertex;
  bool visible;


  int line_loop, lines, line_strip ;
  int quads, polygon, fan;
  int state;
} VertexBuffer = {
  .visible = false,
  .modelview = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1 }
};

static void vertex_buffer_push_index (unsigned int i)
{
  i -= VertexBuffer.vertex;
  array_append (VertexBuffer.index, &i, sizeof(unsigned int));
}

void vertex_buffer_setup()
{
  VertexBuffer.nvertex = 0;
  VertexBuffer.type = -1;
  VertexBuffer.dim = -1;
  VertexBuffer.position = array_new();
  VertexBuffer.normal = array_new();
  VertexBuffer.color = array_new();
  VertexBuffer.index = array_new();
}

void vertex_buffer_free()
{
  array_free (VertexBuffer.position);
  VertexBuffer.position = NULL;
  array_free (VertexBuffer.normal);
  VertexBuffer.normal = NULL;
  array_free (VertexBuffer.color);
  VertexBuffer.color = NULL;
  array_free (VertexBuffer.index);
  VertexBuffer.index = NULL;
}

static void vertex_buffer_glBegin (int state)
{
  if (VertexBuffer.index) {

    glGetFloatv (GL_MODELVIEW_MATRIX, VertexBuffer.modelview);

    bview * view = get_view();

    float q[16] = { 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0,
      - view->tx, - view->ty, 3, 1 };
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    gl_build_rotmatrix ((float (*)[4])q, view->quat);
    swap (float, q[1], q[4]);
    swap (float, q[2], q[8]);
    swap (float, q[6], q[9]);
    matrix_multiply (q, VertexBuffer.modelview);
    for (int i = 0; i < 16; i++)
      VertexBuffer.modelview[i] = q[i];

    VertexBuffer.state = state;
    switch (state) {
    case GL_LINE_LOOP:
      VertexBuffer.line_loop = VertexBuffer.nvertex;
      break;
    case GL_LINES:
      VertexBuffer.lines = VertexBuffer.nvertex;
      break;
    case GL_LINE_STRIP:
      VertexBuffer.line_strip = VertexBuffer.nvertex;
      break;
    case GL_QUADS:
      VertexBuffer.quads = VertexBuffer.nvertex;
      break;
    case GL_POLYGON:
      VertexBuffer.polygon = VertexBuffer.nvertex;
      break;
    case GL_TRIANGLE_FAN:
      VertexBuffer.fan = VertexBuffer.nvertex;
      break;
    default:
      fprintf (ferr, "glBegin (%d) not implemented yet\n", state);
      break;
    }
  }
  glBegin (state);
}

static void vertex_buffer_glEnd()
{
  glEnd();
  if (VertexBuffer.index) {
    int type = -1;
    switch (VertexBuffer.state) {

    case GL_LINE_LOOP:
      for (int i = VertexBuffer.line_loop; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      vertex_buffer_push_index (VertexBuffer.nvertex - 1);
      vertex_buffer_push_index (VertexBuffer.line_loop);
      type = 0;
      break;

    case GL_LINES:
      for (int i = VertexBuffer.lines; i < VertexBuffer.nvertex; i += 2) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case GL_LINE_STRIP:
      for (int i = VertexBuffer.line_strip; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 0;
      break;

    case GL_QUADS:
      for (int i = VertexBuffer.quads; i < VertexBuffer.nvertex; i += 4)
 for (int j = 1; j <= 2; j++) {
   vertex_buffer_push_index (i);
   vertex_buffer_push_index (i + j);
   vertex_buffer_push_index (i + j + 1);
 }
      type = 1;
      break;

    case GL_POLYGON:
      for (int j = 1; j <= VertexBuffer.nvertex - VertexBuffer.polygon - 2;
    j++) {
 vertex_buffer_push_index (VertexBuffer.polygon);
 vertex_buffer_push_index (VertexBuffer.polygon + j);
 vertex_buffer_push_index (VertexBuffer.polygon + j + 1);
      }
      type = 1;
      break;

    case GL_TRIANGLE_FAN:
      for (int i = VertexBuffer.fan + 1; i < VertexBuffer.nvertex - 1; i++) {
 vertex_buffer_push_index (VertexBuffer.fan);
 vertex_buffer_push_index (i);
 vertex_buffer_push_index (i + 1);
      }
      type = 1;
      break;

    default:
      break;
    }
    VertexBuffer.state = 0;
    if (VertexBuffer.type >= 0 && type >= 0)

      assert (VertexBuffer.type == type);
    else
      VertexBuffer.type = type;
  }
}

static void vertex_buffer_glColor3f (float r, float g, float b)
{
  glColor3f (r, g, b);
  if (VertexBuffer.color) {
    struct { float x, y, z; } color = {r, g, b};
    array_append (VertexBuffer.color, &color, 3*sizeof(float));
  }
}

static void vertex_buffer_glNormal3d (double nx, double ny, double nz)
{
  glNormal3d (nx, ny, nz);
  if (VertexBuffer.normal) {
    struct { float x, y, z; } normal = {nx, ny, nz};
    array_append (VertexBuffer.normal, &normal, 3*sizeof(float));
  }
}

static void vertex_buffer_glVertex3d (double x, double y, double z)
{
  glVertex3d (x, y, z);

  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 3)
      VertexBuffer.dim = 3;
    float v[4] = {x, y, z, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
}

static void vertex_buffer_glVertex2d (double x, double y)
{
  glVertex3d (x, y, 0.);

  if (VertexBuffer.position) {
    if (VertexBuffer.dim < 2)
      VertexBuffer.dim = 2;
    float v[4] = {x, y, 0, 1.};
    vector_multiply (v, VertexBuffer.modelview);
    array_append (VertexBuffer.position, v, 3*sizeof(float));
    VertexBuffer.nvertex++;
  }
}
#line 420 "/home/davidtang/basilisk/src/view.h"






#line 1 "draw.h"
#line 1 "/home/davidtang/basilisk/src/draw.h"




#line 1 "fractions.h"
#line 1 "/home/davidtang/basilisk/src/fractions.h"
#line 12 "/home/davidtang/basilisk/src/fractions.h"
#line 1 "geometry.h"
#line 1 "/home/davidtang/basilisk/src/geometry.h"
#line 28 "/home/davidtang/basilisk/src/geometry.h"
double line_alpha (double c, coord n)
{
  double alpha, n1, n2;

  n1 = fabs (n.x); n2 = fabs (n.y);
  if (n1 > n2)
    swap (double, n1, n2);

  c = clamp (c, 0., 1.);
  double v1 = n1/2.;
  if (c <= v1/n2)
    alpha = sqrt (2.*c*n1*n2);
  else if (c <= 1. - v1/n2)
    alpha = c*n2 + v1;
  else
    alpha = n1 + n2 - sqrt (2.*n1*n2*(1. - c));

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;

  return alpha - (n.x + n.y)/2.;
}



double plane_alpha (double c, coord n)
{
  double alpha;
  coord n1;

  n1.x = fabs (n.x); n1.y = fabs (n.y); n1.z = fabs (n.z);

  double m1, m2, m3;
  m1 = min(n1.x, n1.y);
  m3 = max(n1.x, n1.y);
  m2 = n1.z;
  if (m2 < m1) {
    double tmp = m1;
    m1 = m2;
    m2 = tmp;
  }
  else if (m2 > m3) {
    double tmp = m3;
    m3 = m2;
    m2 = tmp;
  }
  double m12 = m1 + m2;
  double pr = max(6.*m1*m2*m3, 1e-50);
  double V1 = m1*m1*m1/pr;
  double V2 = V1 + (m2 - m1)/(2.*m3), V3;
  double mm;
  if (m3 < m12) {
    mm = m3;
    V3 = (m3*m3*(3.*m12 - m3) + m1*m1*(m1 - 3.*m3) + m2*m2*(m2 - 3.*m3))/pr;
  }
  else {
    mm = m12;
    V3 = mm/(2.*m3);
  }

  c = clamp (c, 0., 1.);
  double ch = min(c, 1. - c);
  if (ch < V1)
    alpha = pow (pr*ch, 1./3.);
  else if (ch < V2)
    alpha = (m1 + sqrt(m1*m1 + 8.*m2*m3*(ch - V1)))/2.;
  else if (ch < V3) {
    double p12 = sqrt (2.*m1*m2);
    double q = 3.*(m12 - 2.*m3*ch)/(4.*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + m12;
  }
  else if (m12 < m3)
    alpha = m3*ch + mm/2.;
  else {
    double p = m1*(m2 + m3) + m2*m3 - 1./4., p12 = sqrt(p);
    double q = 3.*m1*m2*m3*(1./2. - ch)/(2.*p*p12);
    double teta = acos(clamp(q,-1.,1.))/3.;
    double cs = cos(teta);
    alpha = p12*(sqrt(3.*(1. - cs*cs)) - cs) + 1./2.;
  }
  if (c > 1./2.) alpha = 1. - alpha;

  if (n.x < 0.)
    alpha += n.x;
  if (n.y < 0.)
    alpha += n.y;
  if (n.z < 0.)
    alpha += n.z;

  return alpha - (n.x + n.y + n.z)/2.;;
}
#line 133 "/home/davidtang/basilisk/src/geometry.h"
double line_area (double nx, double ny, double alpha)
{
  double a, v, area;

  alpha += (nx + ny)/2.;
  if (nx < 0.) {
    alpha -= nx;
    nx = - nx;
  }
  if (ny < 0.) {
    alpha -= ny;
    ny = - ny;
  }

  if (alpha <= 0.)
    return 0.;

  if (alpha >= nx + ny)
    return 1.;

  if (nx < 1e-10)
    area = alpha/ny;
  else if (ny < 1e-10)
    area = alpha/nx;
  else {
    v = sq(alpha);

    a = alpha - nx;
    if (a > 0.)
      v -= a*a;

    a = alpha - ny;
    if (a > 0.)
      v -= a*a;

    area = v/(2.*nx*ny);
  }

  return clamp (area, 0., 1.);
}



double plane_volume (coord n, double alpha)
{
  double al = alpha + (n.x + n.y + n.z)/2. +
    max(0., -n.x) + max(0., -n.y) + max(0., -n.z);
  if (al <= 0.)
    return 0.;
  double tmp = fabs(n.x) + fabs(n.y) + fabs(n.z);
  if (al >= tmp)
    return 1.;
  if (tmp < 1e-10)
    return 0.;
  double n1 = fabs(n.x)/tmp;
  double n2 = fabs(n.y)/tmp;
  double n3 = fabs(n.z)/tmp;
  al = max(0., min(1., al/tmp));
  double al0 = min(al, 1. - al);
  double b1 = min(n1, n2);
  double b3 = max(n1, n2);
  double b2 = n3;
  if (b2 < b1) {
    tmp = b1;
    b1 = b2;
    b2 = tmp;
  }
  else if (b2 > b3) {
    tmp = b3;
    b3 = b2;
    b2 = tmp;
  }
  double b12 = b1 + b2;
  double bm = min(b12, b3);
  double pr = max(6.*b1*b2*b3, 1e-50);
  if (al0 < b1)
    tmp = al0*al0*al0/pr;
  else if (al0 < b2)
    tmp = 0.5*al0*(al0 - b1)/(b2*b3) + b1*b1*b1/pr;
  else if (al0 < bm)
    tmp = (al0*al0*(3.*b12 - al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0))/pr;
  else if (b12 < b3)
    tmp = (al0 - 0.5*bm)/b3;
  else
    tmp = (al0*al0*(3. - 2.*al0) + b1*b1*(b1 - 3.*al0) +
    b2*b2*(b2 - 3.*al0) + b3*b3*(b3 - 3.*al0))/pr;

  double volume = al <= 0.5 ? tmp : 1. - tmp;
  return clamp (volume, 0., 1.);
}
#line 237 "/home/davidtang/basilisk/src/geometry.h"
double rectangle_fraction (coord n, double alpha, coord a, coord b)
{
  coord n1;
  {
#line 240
 {
    alpha -= n.x*(b.x + a.x)/2.;
    n1.x = n.x*(b.x - a.x);
  }
#line 240
 {
    alpha -= n.y*(b.y + a.y)/2.;
    n1.y = n.y*(b.y - a.y);
  }
#line 240
 {
    alpha -= n.z*(b.z + a.z)/2.;
    n1.z = n.z*(b.z - a.z);
  }}
  return plane_volume (n1, alpha);
}
#line 277 "/home/davidtang/basilisk/src/geometry.h"
static coord cube_edge[12][2] = {
  {{0.,0.,0.},{1.,0.,0.}},{{0.,0.,1.},{1.,0.,1.}},
  {{0.,1.,1.},{1.,1.,1.}},{{0.,1.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,1.,0.}},{{0.,0.,1.},{0.,1.,1.}},
  {{1.,0.,1.},{1.,1.,1.}},{{1.,0.,0.},{1.,1.,0.}},
  {{0.,0.,0.},{0.,0.,1.}},{{1.,0.,0.},{1.,0.,1.}},
  {{1.,1.,0.},{1.,1.,1.}},{{0.,1.,0.},{0.,1.,1.}}
};




static int cube_connect[12][2][4] = {
  {{9, 1, 8}, {4, 3, 7}},
  {{6, 2, 5}, {8, 0, 9}},
  {{10, 3, 11}, {5, 1, 6}},
  {{7, 0, 4}, {11, 2, 10}},
  {{3, 7, 0}, {8, 5, 11}},
  {{11, 4, 8}, {1, 6, 2}},
  {{2, 5, 1}, {9, 7, 10}},
  {{10, 6, 9}, {0, 4, 3}},
  {{5, 11, 4}, {0, 9, 1}},
  {{1, 8, 0}, {7, 10, 6}},
  {{6, 9, 7}, {3, 11, 2}},
  {{2, 10, 3}, {4, 8, 5}}
};

int facets (coord n, double alpha, coord v[12], double h)
{
  coord a[12];
  int orient[12];

  for (int i = 0; i < 12; i++) {
    coord e, d;
    double den = 0., t = alpha;
    {
#line 312
 {
      d.x = h*(cube_edge[i][0].x - 0.5);
      e.x = h*(cube_edge[i][1].x - 0.5);
      den += n.x*(e.x - d.x);
      t -= n.x*d.x;
    }
#line 312
 {
      d.y = h*(cube_edge[i][0].y - 0.5);
      e.y = h*(cube_edge[i][1].y - 0.5);
      den += n.y*(e.y - d.y);
      t -= n.y*d.y;
    }
#line 312
 {
      d.z = h*(cube_edge[i][0].z - 0.5);
      e.z = h*(cube_edge[i][1].z - 0.5);
      den += n.z*(e.z - d.z);
      t -= n.z*d.z;
    }}
    orient[i] = -1;
    if (fabs (den) > 1e-10) {
      t /= den;
      if (t >= 0. && t < 1.) {
 double s = - alpha;
 {
#line 323
 {
   a[i].x = d.x + t*(e.x - d.x);
   s += n.x*e.x;
 }
#line 323
 {
   a[i].y = d.y + t*(e.y - d.y);
   s += n.y*e.y;
 }
#line 323
 {
   a[i].z = d.z + t*(e.z - d.z);
   s += n.z*e.z;
 }}
 orient[i] = (s > 0.);
      }
    }
  }

  for (int i = 0; i < 12; i++) {
    int nv = 0, e = i;
    while (orient[e] >= 0) {
      int m = 0, * ne = cube_connect[e][orient[e]];
      v[nv++] = a[e];
      orient[e] = -1;
      while (m < 3 && orient[e] < 0)
 e = ne[m++];
    }
    if (nv > 2)
      return nv;
  }
  return 0;
}






double line_length_center (coord m, double alpha, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  {
#line 357

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 357

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }}

  p->x = p->y = p->z = 0.;

  if (alpha <= 0. || alpha >= n.x + n.y)
    return 0.;

  {
#line 368

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = (m.y < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }
#line 368

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = (m.x < 0. ? 1. - alpha : alpha) - 0.5;
      return 1.;
    }}

  if (alpha >= n.x) {
    p->x += 1.;
    p->y += (alpha - n.x)/n.y;
  }
  else
    p->x += alpha/n.x;

  double ax = p->x, ay = p->y;
  if (alpha >= n.y) {
    p->y += 1.;
    ay -= 1.;
    p->x += (alpha - n.y)/n.x;
    ax -= (alpha - n.y)/n.x;
  }
  else {
    p->y += alpha/n.y;
    ay -= alpha/n.y;
  }

  {
#line 394
 {
    p->x /= 2.;
    p->x = clamp (p->x, 0., 1.);
    if (m.x < 0.)
      p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 394
 {
    p->y /= 2.;
    p->y = clamp (p->y, 0., 1.);
    if (m.y < 0.)
      p->y = 1. - p->y;
    p->y -= 0.5;
  }}

  return sqrt (ax*ax + ay*ay);
}




double plane_area_center (coord m, double alpha, coord * p)
{
  {
#line 410

    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      double length = line_length_center (n, alpha, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return length;
    }
#line 410

    if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      double length = line_length_center (n, alpha, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return length;
    }
#line 410

    if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      double length = line_length_center (n, alpha, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return length;
    }}

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  {
#line 424

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 424

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
#line 424

    if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }}

  double amax = n.x + n.y + n.z;
  if (alpha < 0. || alpha > amax) {
    p->x = p->y = p->z = 0.;
    return 0.;
  }

  double area = sq(alpha);
  p->x = p->y = p->z = area*alpha;

  {
#line 439
 {
    double b = alpha - n.x;
    if (b > 0.) {
      area -= b*b;
      p->x -= b*b*(2.*n.x + alpha);
      p->y -= b*b*b;
      p->z -= b*b*b;
    }
  }
#line 439
 {
    double b = alpha - n.y;
    if (b > 0.) {
      area -= b*b;
      p->y -= b*b*(2.*n.y + alpha);
      p->z -= b*b*b;
      p->x -= b*b*b;
    }
  }
#line 439
 {
    double b = alpha - n.z;
    if (b > 0.) {
      area -= b*b;
      p->z -= b*b*(2.*n.z + alpha);
      p->x -= b*b*b;
      p->y -= b*b*b;
    }
  }}

  amax = alpha - amax;
  {
#line 450
 {
    double b = amax + n.x;
    if (b > 0.) {
      area += b*b;
      p->y += b*b*(2.*n.y + alpha - n.z);
      p->z += b*b*(2.*n.z + alpha - n.y);
      p->x += b*b*b;
    }
  }
#line 450
 {
    double b = amax + n.y;
    if (b > 0.) {
      area += b*b;
      p->z += b*b*(2.*n.z + alpha - n.x);
      p->x += b*b*(2.*n.x + alpha - n.z);
      p->y += b*b*b;
    }
  }
#line 450
 {
    double b = amax + n.z;
    if (b > 0.) {
      area += b*b;
      p->x += b*b*(2.*n.x + alpha - n.y);
      p->y += b*b*(2.*n.y + alpha - n.x);
      p->z += b*b*b;
    }
  }}

  area *= 3.;
  {
#line 461
 {
    if (area) {
      p->x /= area*n.x;
      p->x = clamp (p->x, 0., 1.);
    }
    else
      p->x = 0.;
    if (m.x < 0.) p->x = 1. - p->x;
    p->x -= 0.5;
  }
#line 461
 {
    if (area) {
      p->y /= area*n.y;
      p->y = clamp (p->y, 0., 1.);
    }
    else
      p->y = 0.;
    if (m.y < 0.) p->y = 1. - p->y;
    p->y -= 0.5;
  }
#line 461
 {
    if (area) {
      p->z /= area*n.z;
      p->z = clamp (p->z, 0., 1.);
    }
    else
      p->z = 0.;
    if (m.z < 0.) p->z = 1. - p->z;
    p->z -= 0.5;
  }}

  return area*sqrt (1./(sq(n.x)*sq(n.y)) +
      1./(sq(n.x)*sq(n.z)) +
      1./(sq(n.z)*sq(n.y)))/6.;
}






void line_center (coord m, double alpha, double a, coord * p)
{
  alpha += (m.x + m.y)/2.;

  coord n = m;
  {
#line 487

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 487

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }}

  p->z = 0.;
  if (alpha <= 0.) {
    p->x = p->y = -0.5;
    return;
  }

  if (alpha >= n.x + n.y) {
    p->x = p->y = 0.;
    return;
  }

  {
#line 504

    if (n.x < 1e-4) {
      p->x = 0.;
      p->y = sign(m.y)*(a/2. - 0.5);
      return;
    }
#line 504

    if (n.y < 1e-4) {
      p->y = 0.;
      p->x = sign(m.x)*(a/2. - 0.5);
      return;
    }}

  p->x = p->y = cube(alpha);

  {
#line 513
 {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= sq(b)*(alpha + 2.*n.x);
      p->y -= cube(b);
    }
  }
#line 513
 {
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= sq(b)*(alpha + 2.*n.y);
      p->x -= cube(b);
    }
  }}

  {
#line 521
 {
    p->x /= 6.*sq(n.x)*n.y*a;
    p->x = sign(m.x)*(p->x - 0.5);
  }
#line 521
 {
    p->y /= 6.*sq(n.y)*n.x*a;
    p->y = sign(m.y)*(p->y - 0.5);
  }}
}
#line 534 "/home/davidtang/basilisk/src/geometry.h"
void plane_center (coord m, double alpha, double a, coord * p)
{
  {
#line 536

    if (fabs (m.x) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.y;
      ((double *)&n)[1] = m.z;
      line_center (n, alpha, a, &q);
      p->x = 0.;
      p->y = ((double *)&q)[0];
      p->z = ((double *)&q)[1];
      return;
    }
#line 536

    if (fabs (m.y) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.z;
      ((double *)&n)[1] = m.x;
      line_center (n, alpha, a, &q);
      p->y = 0.;
      p->z = ((double *)&q)[0];
      p->x = ((double *)&q)[1];
      return;
    }
#line 536

    if (fabs (m.z) < 1e-4) {
      coord n, q;
      ((double *)&n)[0] = m.x;
      ((double *)&n)[1] = m.y;
      line_center (n, alpha, a, &q);
      p->z = 0.;
      p->x = ((double *)&q)[0];
      p->y = ((double *)&q)[1];
      return;
    }}

  alpha += (m.x + m.y + m.z)/2.;
  coord n = m;
  {
#line 550

    if (n.x < 0.) {
      alpha -= n.x;
      n.x = - n.x;
    }
#line 550

    if (n.y < 0.) {
      alpha -= n.y;
      n.y = - n.y;
    }
#line 550

    if (n.z < 0.) {
      alpha -= n.z;
      n.z = - n.z;
    }}

  if (alpha <= 0. || a == 0.) {
    p->x = p->y = p->z = -0.5;
    return;
  }

  if (alpha >= n.x + n.y + n.z || a == 1.) {
    p->x = p->y = p->z = 0.;
    return;
  }

  p->x = p->y = p->z = sq(sq(alpha));
  {
#line 567
 {
    double b = alpha - n.x;
    if (b > 0.) {
      p->x -= cube(b)*(3.*n.x + alpha);
      p->y -= sq(sq(b));
      p->z -= sq(sq(b));
    }
  }
#line 567
 {
    double b = alpha - n.y;
    if (b > 0.) {
      p->y -= cube(b)*(3.*n.y + alpha);
      p->z -= sq(sq(b));
      p->x -= sq(sq(b));
    }
  }
#line 567
 {
    double b = alpha - n.z;
    if (b > 0.) {
      p->z -= cube(b)*(3.*n.z + alpha);
      p->x -= sq(sq(b));
      p->y -= sq(sq(b));
    }
  }}

  double amax = alpha - (n.x + n.y + n.z);
  {
#line 577
 {
    double b = amax + n.z;
    if (b > 0.) {
      p->x += cube(b)*(3.*n.x + alpha - n.y);
      p->y += cube(b)*(3.*n.y + alpha - n.x);
      p->z += sq(sq(b));
    }
  }
#line 577
 {
    double b = amax + n.x;
    if (b > 0.) {
      p->y += cube(b)*(3.*n.y + alpha - n.z);
      p->z += cube(b)*(3.*n.z + alpha - n.y);
      p->x += sq(sq(b));
    }
  }
#line 577
 {
    double b = amax + n.y;
    if (b > 0.) {
      p->z += cube(b)*(3.*n.z + alpha - n.x);
      p->x += cube(b)*(3.*n.x + alpha - n.z);
      p->y += sq(sq(b));
    }
  }}

  double b = 24.*n.x*n.y*n.z*a;
  {
#line 587
 {
    p->x /= b*n.x;
    p->x = sign(m.x)*(p->x - 0.5);
  }
#line 587
 {
    p->y /= b*n.y;
    p->y = sign(m.y)*(p->y - 0.5);
  }
#line 587
 {
    p->z /= b*n.z;
    p->z = sign(m.z)*(p->z - 0.5);
  }}
}
#line 13 "/home/davidtang/basilisk/src/fractions.h"
#line 21 "/home/davidtang/basilisk/src/fractions.h"
#line 1 "myc.h"
#line 1 "/home/davidtang/basilisk/src/myc.h"
#line 16 "/home/davidtang/basilisk/src/myc.h"
coord mycs (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 17 "/home/davidtang/basilisk/src/myc.h"

  double m1,m2,m[4][3],t0,t1,t2;
  int cn;



  m1 = val(c,-1,0,-1) + val(c,-1,0,1) + val(c,-1,-1,0) + val(c,-1,1,0) +
       val(c,-1,0,0);
  m2 = val(c,1,0,-1) + val(c,1,0,1) + val(c,1,-1,0) + val(c,1,1,0) +
       val(c,1,0,0);
  m[0][0] = m1 > m2 ? 1. : -1.;

  m1 = val(c,-1,-1,0)+ val(c,1,-1,0)+ val(c,0,-1,0);
  m2 = val(c,-1,1,0)+ val(c,1,1,0)+ val(c,0,1,0);
  m[0][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1)+ val(c,1,0,-1)+ val(c,0,0,-1);
  m2 = val(c,-1,0,1)+ val(c,1,0,1)+ val(c,0,0,1);
  m[0][2] = 0.5*(m1-m2);



  m1 = val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,0);
  m2 = val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,0);
  m[1][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1) + val(c,0,-1,1) + val(c,1,-1,0) + val(c,-1,-1,0) +
       val(c,0,-1,0);
  m2 = val(c,0,1,-1) + val(c,0,1,1) + val(c,1,1,0) + val(c,-1,1,0) +
       val(c,0,1,0);
  m[1][1] = m1 > m2 ? 1. : -1.;

  m1 = val(c,0,-1,-1)+ val(c,0,0,-1)+ val(c,0,1,-1);
  m2 = val(c,0,-1,1)+ val(c,0,0,1)+ val(c,0,1,1);
  m[1][2] = 0.5*(m1-m2);




  m1 = val(c,-1,0,-1)+ val(c,-1,0,1)+ val(c,-1,0,0);
  m2 = val(c,1,0,-1)+ val(c,1,0,1)+ val(c,1,0,0);
  m[2][0] = 0.5*(m1-m2);

  m1 = val(c,0,-1,-1)+ val(c,0,-1,1)+ val(c,0,-1,0);
  m2 = val(c,0,1,-1)+ val(c,0,1,1)+ val(c,0,1,0);
  m[2][1] = 0.5*(m1-m2);

  m1 = val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1) +
       val(c,0,0,-1);
  m2 = val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1) +
       val(c,0,0,1);
  m[2][2] = m1 > m2 ? 1. : -1.;


  t0 = fabs(m[0][0]) + fabs(m[0][1]) + fabs(m[0][2]);
  m[0][0] /= t0;
  m[0][1] /= t0;
  m[0][2] /= t0;

  t0 = fabs(m[1][0]) + fabs(m[1][1]) + fabs(m[1][2]);
  m[1][0] /= t0;
  m[1][1] /= t0;
  m[1][2] /= t0;

  t0 = fabs(m[2][0]) + fabs(m[2][1]) + fabs(m[2][2]);
  m[2][0] /= t0;
  m[2][1] /= t0;
  m[2][2] /= t0;


  t0 = fabs(m[0][0]);
  t1 = fabs(m[1][1]);
  t2 = fabs(m[2][2]);

  cn = 0;
  if (t1 > t0) {
    t0 = t1;
    cn = 1;
  }
  if (t2 > t0)
    cn = 2;


  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,-1,-1,1) + val(c,-1,1,1) +
       2.*(val(c,-1,-1,0) + val(c,-1,1,0) + val(c,-1,0,-1) + val(c,-1,0,1)) +
       4.*val(c,-1,0,0);
  m2 = val(c,1,-1,-1) + val(c,1,1,-1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,1,-1,0) + val(c,1,1,0) + val(c,1,0,-1) + val(c,1,0,1)) +
       4.*val(c,1,0,0);
  m[3][0] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,-1,1) + val(c,1,-1,-1) + val(c,1,-1,1) +
       2.*( val(c,-1,-1,0) + val(c,1,-1,0) + val(c,0,-1,-1) + val(c,0,-1,1)) +
       4.*val(c,0,-1,0);
  m2 = val(c,-1,1,-1) + val(c,-1,1,1) + val(c,1,1,-1) + val(c,1,1,1) +
       2.*(val(c,-1,1,0) + val(c,1,1,0) + val(c,0,1,-1) + val(c,0,1,1)) +
       4.*val(c,0,1,0);
  m[3][1] = m1 - m2;

  m1 = val(c,-1,-1,-1) + val(c,-1,1,-1) + val(c,1,-1,-1) + val(c,1,1,-1) +
       2.*(val(c,-1,0,-1) + val(c,1,0,-1) + val(c,0,-1,-1) + val(c,0,1,-1)) +
       4.*val(c,0,0,-1);
  m2 = val(c,-1,-1,1) + val(c,-1,1,1) + val(c,1,-1,1) + val(c,1,1,1) +
       2.*(val(c,-1,0,1) + val(c,1,0,1) + val(c,0,-1,1) + val(c,0,1,1)) +
       4.*val(c,0,0,1);
  m[3][2] = m1 - m2;


  t0 = fabs(m[3][0]) + fabs(m[3][1]) + fabs(m[3][2]);
  if (t0 < 1e-30) {
    coord mxyz = {1., 0., 0.};
    return mxyz;
  }

  m[3][0] /= t0;
  m[3][1] /= t0;
  m[3][2] /= t0;


  t0 = fabs (m[3][0]);
  t1 = fabs (m[3][1]);
  t2 = fabs (m[3][2]);
  if (t1 > t0)
    t0 = t1;
  if (t2 > t0)
    t0 = t2;

  if (fabs(m[cn][cn]) > t0)
    cn = 3;


  coord mxyz = {m[cn][0], m[cn][1], m[cn][2]};
  return mxyz;
}
#line 22 "/home/davidtang/basilisk/src/fractions.h"
#line 121 "/home/davidtang/basilisk/src/fractions.h"
struct Fractions {
  scalar Phi;
  scalar c;
  vector s;
  double val;
};


void fractions (struct Fractions a)
{ trace ("fractions", "/home/davidtang/basilisk/src/fractions.h", 130);
  scalar Phi = a.Phi;
  scalar c = a.c;
  vector s = (a.s).x.i ? (a.s) : new_face_vector("s");
  double val = a.val;
#line 143 "/home/davidtang/basilisk/src/fractions.h"
  vector p= new_vector("p");
#line 155 "/home/davidtang/basilisk/src/fractions.h"
   { foreach_vertex_aux(){

#line 155 "/home/davidtang/basilisk/src/fractions.h"
 {
#line 155
 if (_a.x < point.n + 2) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,1,0,0) - val) < 0.) {






      val(p.x,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,1,0,0));
      if (val(Phi,0,0,0) < val)
 val(p.x,0,0,0) = 1. - val(p.x,0,0,0);
    }
#line 180 "/home/davidtang/basilisk/src/fractions.h"
    else
      val(p.x,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,1,0,0) > val);
  }
#line 155
 if (_a.y < point.n + 2) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,1,0) - val) < 0.) {






      val(p.y,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,1,0));
      if (val(Phi,0,0,0) < val)
 val(p.y,0,0,0) = 1. - val(p.y,0,0,0);
    }
#line 180 "/home/davidtang/basilisk/src/fractions.h"
    else
      val(p.y,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,1,0) > val);
  }
#line 155
 if (_a.z < point.n + 2) {





    if ((val(Phi,0,0,0) - val)*(val(Phi,0,0,1) - val) < 0.) {






      val(p.z,0,0,0) = (val(Phi,0,0,0) - val)/(val(Phi,0,0,0) - val(Phi,0,0,1));
      if (val(Phi,0,0,0) < val)
 val(p.z,0,0,0) = 1. - val(p.z,0,0,0);
    }
#line 180 "/home/davidtang/basilisk/src/fractions.h"
    else
      val(p.z,0,0,0) = (val(Phi,0,0,0) > val || val(Phi,0,0,1) > val);
  }} } end_foreach_vertex_aux(); }
#line 194 "/home/davidtang/basilisk/src/fractions.h"
  scalar s_x = s.x, s_y = s.y, s_z = s.z;
   { foreach_face_generic() { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 195
{

#line 195 "/home/davidtang/basilisk/src/fractions.h"






  {
#line 233 "/home/davidtang/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 235
 {
      n.x = val(p.y,0,0,0) - val(p.y,1,0,0);
      nn += fabs(n.x);
    }
#line 235
 {
      n.y = val(p.x,0,0,0) - val(p.x,0,1,0);
      nn += fabs(n.y);
    }}





    if (nn == 0.)
      val(s_z,0,0,0) = val(p.x,0,0,0);
    else {





      {
#line 252

 n.x /= nn;
#line 252

 n.y /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 262

   if (val(p.x,0,i,0) > 0. && val(p.x,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.x,0,i,0) - 0.5);
     alpha += n.x*a + n.y*(i - 0.5);
     ni++;
   }
#line 262

   if (val(p.y,i,0,0) > 0. && val(p.y,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.y,i,0,0) - 0.5);
     alpha += n.y*a + n.x*(i - 0.5);
     ni++;
   }}
#line 276 "/home/davidtang/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_z,0,0,0) = max (val(p.x,0,0,0), val(p.y,0,0,0));
      else if (ni != 4)
 val(s_z,0,0,0) = line_area (n.x, n.y, alpha/ni);
      else {

 val(s_z,0,0,0) = (val(p.x,0,0,0) + val(p.x,0,1,0) + val(p.y,0,0,0) + val(p.y,1,0,0) > 2.);



      }
    }
  } }  }}  { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 195
{

#line 195 "/home/davidtang/basilisk/src/fractions.h"






  {
#line 233 "/home/davidtang/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 235
 {
      n.y = val(p.z,0,0,0) - val(p.z,0,1,0);
      nn += fabs(n.y);
    }
#line 235
 {
      n.z = val(p.y,0,0,0) - val(p.y,0,0,1);
      nn += fabs(n.z);
    }}





    if (nn == 0.)
      val(s_x,0,0,0) = val(p.y,0,0,0);
    else {





      {
#line 252

 n.y /= nn;
#line 252

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 262

   if (val(p.y,0,0,i) > 0. && val(p.y,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i) - val)*(val(p.y,0,0,i) - 0.5);
     alpha += n.y*a + n.z*(i - 0.5);
     ni++;
   }
#line 262

   if (val(p.z,0,i,0) > 0. && val(p.z,0,i,0) < 1.) {
     double a = sign(val(Phi,0,i,0) - val)*(val(p.z,0,i,0) - 0.5);
     alpha += n.z*a + n.y*(i - 0.5);
     ni++;
   }}
#line 276 "/home/davidtang/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_x,0,0,0) = max (val(p.y,0,0,0), val(p.z,0,0,0));
      else if (ni != 4)
 val(s_x,0,0,0) = line_area (n.y, n.z, alpha/ni);
      else {

 val(s_x,0,0,0) = (val(p.y,0,0,0) + val(p.y,0,0,1) + val(p.z,0,0,0) + val(p.z,0,1,0) > 2.);



      }
    }
  } }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 195
{

#line 195 "/home/davidtang/basilisk/src/fractions.h"






  {
#line 233 "/home/davidtang/basilisk/src/fractions.h"
    coord n;
    double nn = 0.;
    {
#line 235
 {
      n.z = val(p.x,0,0,0) - val(p.x,0,0,1);
      nn += fabs(n.z);
    }
#line 235
 {
      n.x = val(p.z,0,0,0) - val(p.z,1,0,0);
      nn += fabs(n.x);
    }}





    if (nn == 0.)
      val(s_y,0,0,0) = val(p.z,0,0,0);
    else {





      {
#line 252

 n.z /= nn;
#line 252

 n.x /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 {
#line 262

   if (val(p.z,i,0,0) > 0. && val(p.z,i,0,0) < 1.) {
     double a = sign(val(Phi,i,0,0) - val)*(val(p.z,i,0,0) - 0.5);
     alpha += n.z*a + n.x*(i - 0.5);
     ni++;
   }
#line 262

   if (val(p.x,0,0,i) > 0. && val(p.x,0,0,i) < 1.) {
     double a = sign(val(Phi,0,0,i) - val)*(val(p.x,0,0,i) - 0.5);
     alpha += n.x*a + n.z*(i - 0.5);
     ni++;
   }}
#line 276 "/home/davidtang/basilisk/src/fractions.h"
      if (ni == 0)
 val(s_y,0,0,0) = max (val(p.z,0,0,0), val(p.x,0,0,0));
      else if (ni != 4)
 val(s_y,0,0,0) = line_area (n.z, n.x, alpha/ni);
      else {

 val(s_y,0,0,0) = (val(p.z,0,0,0) + val(p.z,1,0,0) + val(p.x,0,0,0) + val(p.x,0,0,1) > 2.);



      }
    }
  } }  }}  end_foreach_face_generic()
#line 288
 end_foreach_face(); }







  boundary_flux (((vector []){{s.x,s.y,s.z},{{-1},{-1},{-1}}}));
   { foreach(){

#line 297 "/home/davidtang/basilisk/src/fractions.h"
 {




    coord n;
    double nn = 0.;
    {
#line 304
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 304
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
#line 304
 {
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }}
    if (nn == 0.)
      val(c,0,0,0) = val(s.x,0,0,0);
    else {
      {
#line 311

 n.x /= nn;
#line 311

 n.y /= nn;
#line 311

 n.z /= nn;}






      double alpha = 0., ni = 0.;
      for (int i = 0; i <= 1; i++)
 for (int j = 0; j <= 1; j++)
   {
#line 322

     if (val(p.x,0,i,j) > 0. && val(p.x,0,i,j) < 1.) {
       double a = sign(val(Phi,0,i,j) - val)*(val(p.x,0,i,j) - 0.5);
       alpha += n.x*a + n.y*(i - 0.5) + n.z*(j - 0.5);
       ni++;
     }
#line 322

     if (val(p.y,j,0,i) > 0. && val(p.y,j,0,i) < 1.) {
       double a = sign(val(Phi,j,0,i) - val)*(val(p.y,j,0,i) - 0.5);
       alpha += n.y*a + n.z*(i - 0.5) + n.x*(j - 0.5);
       ni++;
     }
#line 322

     if (val(p.z,i,j,0) > 0. && val(p.z,i,j,0) < 1.) {
       double a = sign(val(Phi,i,j,0) - val)*(val(p.z,i,j,0) - 0.5);
       alpha += n.z*a + n.x*(i - 0.5) + n.y*(j - 0.5);
       ni++;
     }}




      if (ni == 0)
 val(c,0,0,0) = val(s.x,0,0,0);
      else if (ni < 3 || ni > 6)
 val(c,0,0,0) = 0.;
      else
 val(c,0,0,0) = plane_volume (n, alpha/ni);
    }
  } } end_foreach(); }





  boundary (((scalar []){c,{-1}}));
 delete (((scalar []){p.x,p.y,p.z,{-1}}));  { if (!(a.s).x.i) delete (((scalar []){s.x,s.y,s.z,{-1}})); }  end_trace("fractions", "/home/davidtang/basilisk/src/fractions.h", 346); }
#line 384 "/home/davidtang/basilisk/src/fractions.h"
coord youngs_normal (Point point, scalar c)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 385 "/home/davidtang/basilisk/src/fractions.h"

  coord n;
  double nn = 0.;
  assert (3 == 2);
  {
#line 389
 {
    n.x = (val(c,-1,1,0) + 2.*val(c,-1,0,0) + val(c,-1,-1,0) -
    val(c,+1,1,0) - 2.*val(c,+1,0,0) - val(c,+1,-1,0));
    nn += fabs(n.x);
  }
#line 389
 {
    n.y = (val(c,0,-1,1) + 2.*val(c,0,-1,0) + val(c,0,-1,-1) -
    val(c,0,+1,1) - 2.*val(c,0,+1,0) - val(c,0,+1,-1));
    nn += fabs(n.y);
  }
#line 389
 {
    n.z = (val(c,1,0,-1) + 2.*val(c,0,0,-1) + val(c,-1,0,-1) -
    val(c,1,0,+1) - 2.*val(c,0,0,+1) - val(c,-1,0,+1));
    nn += fabs(n.z);
  }}

  if (nn > 0.)
    {
#line 396

      n.x /= nn;
#line 396

      n.y /= nn;
#line 396

      n.z /= nn;}
  else
    n.x = 1.;
  return n;
}





coord facet_normal (Point point, scalar c, vector s)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 408 "/home/davidtang/basilisk/src/fractions.h"

  if (s.x.i >= 0) {
    coord n;
    double nn = 0.;
    {
#line 412
 {
      n.x = val(s.x,0,0,0) - val(s.x,1,0,0);
      nn += fabs(n.x);
    }
#line 412
 {
      n.y = val(s.y,0,0,0) - val(s.y,0,1,0);
      nn += fabs(n.y);
    }
#line 412
 {
      n.z = val(s.z,0,0,0) - val(s.z,0,0,1);
      nn += fabs(n.z);
    }}
    if (nn > 0.)
      {
#line 417

 n.x /= nn;
#line 417

 n.y /= nn;
#line 417

 n.z /= nn;}
    else
      {
#line 420

 n.x = 1./3;
#line 420

 n.y = 1./3;
#line 420

 n.z = 1./3;}
    return n;
  }
  return mycs (point, c);
}
#line 434 "/home/davidtang/basilisk/src/fractions.h"

void reconstruction (const scalar c, vector n, scalar alpha)
{ trace ("reconstruction", "/home/davidtang/basilisk/src/fractions.h", 436);
   { foreach(){

#line 437 "/home/davidtang/basilisk/src/fractions.h"
 {





    if (val(c,0,0,0) <= 0. || val(c,0,0,0) >= 1.) {
      val(alpha,0,0,0) = 0.;
      {
#line 445

 val(n.x,0,0,0) = 0.;
#line 445

 val(n.y,0,0,0) = 0.;
#line 445

 val(n.z,0,0,0) = 0.;}
    }
    else {






      coord m = mycs (point, c);
      {
#line 456

 val(n.x,0,0,0) = m.x;
#line 456

 val(n.y,0,0,0) = m.y;
#line 456

 val(n.z,0,0,0) = m.z;}
      val(alpha,0,0,0) = plane_alpha (val(c,0,0,0), m);
    }
  } } end_foreach(); }
#line 484 "/home/davidtang/basilisk/src/fractions.h"
  boundary (((scalar []){n.x,n.y,n.z,alpha,{-1}}));
 end_trace("reconstruction", "/home/davidtang/basilisk/src/fractions.h", 485); }
#line 505 "/home/davidtang/basilisk/src/fractions.h"
struct OutputFacets {
  scalar c;
  FILE * fp;
  vector s;
};


void output_facets (struct OutputFacets p)
{ trace ("output_facets", "/home/davidtang/basilisk/src/fractions.h", 513);
  scalar c = p.c;
  vector s = p.s;
  if (!p.fp) p.fp = fout;
  if (!s.x.i) s.x.i = -1;

   { foreach(){

#line 519 "/home/davidtang/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = facet_normal (point, c, s);
      double alpha = plane_alpha (val(c,0,0,0), n);







      coord v[12];
      int m = facets (n, alpha, v, 1.);
      for (int i = 0; i < m; i++)
 fprintf (p.fp, "%g %g %g\n",
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
      if (m > 0)
 fputc ('\n', p.fp);

    } } end_foreach(); }

  fflush (p.fp);
 end_trace("output_facets", "/home/davidtang/basilisk/src/fractions.h", 541); }








double interface_area (scalar c)
{ trace ("interface_area", "/home/davidtang/basilisk/src/fractions.h", 551);
  double area = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _area = area; 
#line 553
foreach (){

#line 553 "/home/davidtang/basilisk/src/fractions.h"

    if (val(c,0,0,0) > 1e-6 && val(c,0,0,0) < 1. - 1e-6) {
      coord n = mycs (point, c), p;
      double alpha = plane_alpha (val(c,0,0,0), n);
      _area += pow(Delta, 3 - 1)*plane_area_center (n, alpha, &p);
    } } end_foreach();OMP(omp critical) area += _area;
mpi_all_reduce_double (area, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 558
 }
  { double _ret =  area; end_trace("interface_area", "/home/davidtang/basilisk/src/fractions.h", 559);  return _ret; }
 end_trace("interface_area", "/home/davidtang/basilisk/src/fractions.h", 560); }
#line 6 "/home/davidtang/basilisk/src/draw.h"
#line 1 "gl/font.h"
#line 1 "/home/davidtang/basilisk/src/gl/font.h"
#line 27 "/home/davidtang/basilisk/src/gl/font.h"
#include <stdio.h>
#line 1 "gl/og_font.h"
#line 1 "/home/davidtang/basilisk/src/gl/og_font.h"




typedef struct tagSOG_StrokeVertex SOG_StrokeVertex;
struct tagSOG_StrokeVertex
{
    GLfloat X, Y;
};

typedef struct tagSOG_StrokeStrip SOG_StrokeStrip;
struct tagSOG_StrokeStrip
{
    int Number;
    const SOG_StrokeVertex *Vertices;
};

typedef struct tagSOG_StrokeChar SOG_StrokeChar;
struct tagSOG_StrokeChar
{
    GLfloat Right;
    int Number;
    const SOG_StrokeStrip* Strips;
};

typedef struct tagSOG_StrokeFont SOG_StrokeFont;
struct tagSOG_StrokeFont
{
    char *Name;
    int Quantity;
    GLfloat Height;
    const SOG_StrokeChar **Characters;
};
#line 29 "/home/davidtang/basilisk/src/gl/font.h"
#line 39 "/home/davidtang/basilisk/src/gl/font.h"
extern SOG_StrokeFont ogStrokeMonoRoman;
#line 48 "/home/davidtang/basilisk/src/gl/font.h"
static SOG_StrokeFont *oghStrokeByID( void *font )
{


    if( font == ((void *)0x0001) )
        return &ogStrokeMonoRoman;

    fprintf (ferr, "stroke font %p not found", font );
    return 0;
}
#line 83 "/home/davidtang/basilisk/src/gl/font.h"
void gl_StrokeCharacter( int character )
{
    void *fontID = ((void *)0x0001);
    const SOG_StrokeChar *schar;
    const SOG_StrokeStrip *strip;
    int i, j;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( !font ||
        ( 1 > character ) ||
        ( font->Quantity < character ) )
        return;

    schar = font->Characters[ character ];
    if( schar )
    {
        strip = schar->Strips;

        for( i = 0; i < schar->Number; i++, strip++ )
        {
            vertex_buffer_glBegin( GL_LINE_STRIP );
            for( j = 0; j < strip->Number; j++ )
                vertex_buffer_glVertex2d( strip->Vertices[ j ].X, strip->Vertices[ j ].Y );
            vertex_buffer_glEnd( );
        }
        glTranslatef( schar->Right, 0.0, 0.0 );
    }
}
#line 147 "/home/davidtang/basilisk/src/gl/font.h"
void gl_StrokeString( const char *string )
{
    void *fontID = ((void *)0x0001);
    int i, j;
    float length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );
    unsigned char c;

    if( font && string )





        while(( c = *string++ ))
       if( c < font->Quantity ) {
                if( c == '\n' )
                {
                    glTranslatef ( -length, -( float )( font->Height ), 0.0 );
                    length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                    {
                        const SOG_StrokeStrip *strip = schar->Strips;

                        for( i = 0; i < schar->Number; i++, strip++ )
                        {
                            vertex_buffer_glBegin( GL_LINE_STRIP );

                            for( j = 0; j < strip->Number; j++ )
                                vertex_buffer_glVertex2d( strip->Vertices[ j ].X,
                                            strip->Vertices[ j ].Y);

                            vertex_buffer_glEnd( );
                        }

                        length += schar->Right;
                        glTranslatef( schar->Right, 0.0, 0.0 );
                    }
                }
     }
}
#line 226 "/home/davidtang/basilisk/src/gl/font.h"
float gl_StrokeWidth( int character )
{
    void *fontID = ((void *)0x0001);
    float ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font &&
        ( 0 < character ) &&
        ( font->Quantity > character ) )
    {
        const SOG_StrokeChar *schar = font->Characters[ character ];
        if( schar )
            ret = schar->Right;
    }

    return ret;
}
#line 269 "/home/davidtang/basilisk/src/gl/font.h"
float gl_StrokeLength( const char *string )
{
    void *fontID = ((void *)0x0001);
    unsigned char c;
    float length = 0.0;
    float this_line_length = 0.0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font && string )
        while(( c = *string++ ))
            if( c < font->Quantity )
            {
                if( c == '\n' )
                {
                    if( length < this_line_length )
                        length = this_line_length;
                    this_line_length = 0.0;
                }
                else
                {
                    const SOG_StrokeChar *schar =
                        font->Characters[ c ];
                    if( schar )
                        this_line_length += schar->Right;
                }
            }

    if( length < this_line_length )
        length = this_line_length;
    return length;
}
#line 321 "/home/davidtang/basilisk/src/gl/font.h"
GLfloat gl_StrokeHeight()
{
    void *fontID = ((void *)0x0001);
    GLfloat ret = 0;
    SOG_StrokeFont *font = oghStrokeByID( fontID );

    if( font )
        ret = font->Height;

    return ret;
}
#line 7 "/home/davidtang/basilisk/src/draw.h"




void clear()
{
  bview * view = get_view();
  if (view->active)
    view->active = false;
  draw();
}
#line 49 "/home/davidtang/basilisk/src/draw.h"
struct _view_set {
  float tx, ty;
  float fov;
  float quat[4];
  float sx, sy, sz;
  unsigned width, height, samples;
  float bg[3];
  float theta, phi, psi;
  bool relative;
  float tz, near, far;
  float res;
  char * camera;
  void (* map) (coord *);
  int cache;
  float p1x, p1y, p2x, p2y;
  bview * view;
};

void view (struct _view_set p)
{
  bview * v = p.view ? p.view : get_view();
  if (p.fov) {
    if (p.relative)
      v->fov += (0.1 + 3.*v->fov)*p.fov;
    else
      v->fov = p.fov;
    v->fov = clamp(v->fov,0.01,100.);
  }
  for (int i = 0; i < 4; i++)
    if (p.quat[i]) {
      for (int j = 0; j < 4; j++)
 v->quat[j] = p.quat[j];
      break;
    }
  if (p.tx) v->tx = p.relative ? v->tx + p.tx*0.02*(0.01 + 3.*v->fov) : p.tx;
  if (p.ty) v->ty = p.relative ? v->ty + p.ty*0.02*(0.01 + 3.*v->fov) : p.ty;
  if (p.sx) v->sx = p.sx;
  if (p.sy) v->sy = p.sy;
  if (p.sz) v->sz = p.sz;
  if (p.bg[0] || p.bg[1] || p.bg[2])
    for (int i = 0; i < 3; i++)
      v->bg[i] = p.bg[i];

  if (p.camera) {
    v->gfsview = false;
    if (strlen(p.camera) >= 4 &&
 !strcmp (&p.camera[strlen(p.camera) - 4], ".gfv")) {
      FILE * fp = fopen (p.camera, "r");
      if (!fp) {
 perror (p.camera);
 exit (1);
      }
      char s[81];
      float q[4], fov;
      int nq = 0, nf = 0;
      while (fgets (s, 81, fp) && (!nq || !nf)) {
 if (!nq)
   nq = sscanf (s, "  q0 = %f q1 = %f q2 = %f q3 = %f",
         &q[0], &q[1], &q[2], &q[3]);
 if (!nf)
   nf = sscanf (s, "  fov = %f", &fov);
      }
      if (nq != 4 || nf != 1) {
 fprintf (ferr, "%s: not a valid gfv file\n", p.camera);
 exit (1);
      }
      for (int j = 0; j < 4; j++)
 v->quat[j] = q[j];
      v->fov = fov;
      v->gfsview = true;
    }
    else if (!strcmp (p.camera, "left"))
      gl_axis_to_quat ((float[]){0,1,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "right"))
      gl_axis_to_quat ((float[]){0,1,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "top"))
      gl_axis_to_quat ((float[]){1,0,0}, - pi/2., v->quat);
    else if (!strcmp (p.camera, "bottom"))
      gl_axis_to_quat ((float[]){1,0,0}, pi/2., v->quat);
    else if (!strcmp (p.camera, "front"))
      gl_axis_to_quat ((float[]){0,0,1}, 0., v->quat);
    else if (!strcmp (p.camera, "back"))
      gl_axis_to_quat ((float[]){0,1,0}, pi, v->quat);
    else if (!strcmp (p.camera, "iso")) {
      gl_axis_to_quat ((float[]){0,1,0}, pi/4., v->quat);
      float q[4];
      gl_axis_to_quat ((float[]){1,0,0}, - pi/4., q);
      gl_add_quats(q, v->quat, v->quat);
    }
    else {
      fprintf (ferr, "view(): unknown camera '%s'\n", p.camera);
      exit (1);
    }
  }
  else if (p.theta || p.phi || p.psi) {
    v->gfsview = false;
    float q[4];
    gl_axis_to_quat ((float[]){1,0,0}, - p.phi, q);
    if (p.relative) {
      float q1[4];
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, q1);
      gl_add_quats(q, q1, q1);
      float q2[4];
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q2);
      gl_add_quats(q1, q2, q2);
      gl_add_quats(q2, v->quat, v->quat);
    }
    else {
      gl_axis_to_quat ((float[]){0,1,0}, p.theta, v->quat);
      gl_add_quats(q, v->quat, v->quat);
      gl_axis_to_quat ((float[]){0,0,1}, p.psi, q);
      gl_add_quats(q, v->quat, v->quat);
    }
  }

  if (p.map)
    v->map = p.map;

  if (p.p1x || p.p1y || p.p2x || p.p2y) {
    float q[4];
    gl_trackball(q, p.p1x, p.p1y, p.p2x, p.p2y);
    gl_add_quats (q, v->quat, v->quat);
  }

  if (p.far > p.near) {
    v->tz = p.tz;
    v->far = p.far;
    v->near = p.near;
  }

  if (p.res)
    v->res = p.res;

  if ((p.width && p.width != v->width) ||
      (p.height && p.height != v->height) ||
      (p.samples && p.samples != v->samples)) {
    v->width = v->width/v->samples;
    v->height = v->height/v->samples;
    if (p.width) v->width = p.width;
    if (p.height) v->height = p.height;
    if (p.samples) v->samples = p.samples;
    v->width *= v->samples;
    v->height *= v->samples;
    framebuffer_destroy (v->fb);
    v->fb = framebuffer_new (v->width, v->height);
    init_gl();
  }

  if (p.cache > 0) {
    v->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
    v->maxlen = p.cache;
  }

  clear();
}







struct _translate {
  float x, y, z;
};

void begin_translate (struct _translate p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glTranslatef (p.x, p.y, p.z);
  gl_get_frustum (&view->frustum);
}

void end_translate()
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  gl_get_frustum (&view->frustum);
}
#line 240 "/home/davidtang/basilisk/src/draw.h"
struct _mirror {
  coord n;
  double alpha;
};

void begin_mirror (struct _mirror p)
{
  bview * view = draw();
  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  normalize (&p.n);
  GLfloat s[16], t[16];
  s[0] = 1. - 2.*p.n.x*p.n.x;
  s[1] = - 2.*p.n.x*p.n.y; s[2] = - 2.*p.n.x*p.n.z;
  s[3] = 0.;
  s[4] = s[1];
  s[5] = 1. - 2.*p.n.y*p.n.y; s[6] = - 2.*p.n.y*p.n.z;
  s[7] = 0.;
  s[8] = s[2]; s[9] = s[6]; s[10] = 1. - 2.*p.n.z*p.n.z;
  s[11] = 0.;
  s[12] = 0.; s[13] = 0.; s[14] = 0.;
  s[15] = 1.;

  t[0] = 1.; t[1] = 0.; t[2] = 0.; t[3] = 0.;
  t[4] = 0.; t[5] = 1.; t[6] = 0.; t[7] = 0.;
  t[8] = 0.; t[9] = 0.; t[10] = 1.; t[11] = 0.;
  t[12] = - 2.*p.n.x*p.alpha;
  t[13] = - 2.*p.n.y*p.alpha;
  t[14] = - 2.*p.n.z*p.alpha;
  t[15] = 1.;
  matrix_multiply (s, t);
  glMultMatrixf (s);
  gl_get_frustum (&view->frustum);
  view->reversed = !view->reversed;
}

void end_mirror() {
  end_translate();
  bview * view = draw();
  view->reversed = !view->reversed;
}







static void mapped_position (bview * view, coord * p, double * r)
{
  double x = p->x, y = p->y, z = p->z, rm = 0.;
  view->map (p);
  for (int i = -1; i <= 1; i += 2)
    for (int j = -1; j <= 1; j += 2)
      for (int k = -1; k <= 1; k += 2) {
 coord q = {x + i**r, y + j**r, z + k**r};
 view->map (&q);
 double pq = sq(p->x - q.x) + sq(p->y - q.y) + sq(p->z - q.z);
 if (pq > rm)
   rm = pq;
      }
  *r = sqrt (rm);
}

#define foreach_visible(view)\
foreach_cell() {\
\
\
\
  double _r = Delta*0.87;\
\
  coord _p = {x, y, z};\
  if ((view)->map)\
    mapped_position (view, &_p, &_r);\
  if (VertexBuffer.visible &&\
      !sphere_in_frustum (_p.x, _p.y, _p.z, _r, &(view)->frustum))\
    continue;\
  if (is_leaf(cell) ||\
      (VertexBuffer.visible &&\
       sphere_diameter (_p.x, _p.y, _p.z, _r/L0, &(view)->frustum)\
       < (view)->res)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 322

#define end_foreach_visible()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 329

#line 340 "/home/davidtang/basilisk/src/draw.h"
static void glnormal3d (bview * view, double x, double y, double z) {

  if (view->gfsview || view->reversed)
    vertex_buffer_glNormal3d (- x, - y, - z);
  else
    vertex_buffer_glNormal3d (x, y, z);
}

#define foreach_visible_plane(view, n1, alpha1)\
coord n = {(n1).x, (n1).y, (n1).z};\
double _alpha = 0.9999999*(alpha1);\
{\
  double norm = sqrt(sq(n.x) + sq(n.y) + sq(n.z));\
  if (!norm)\
    n.z = 1.;\
  else\
    n.x /= norm, n.y /= norm, n.z /= norm, _alpha /= norm;\
}\
glnormal3d (view, n.x, n.y, n.z);\
foreach_cell() {\
\
  double _r = Delta*0.87, alpha = (_alpha - n.x*x - n.y*y - n.z*z)/Delta;\
  if (fabs(alpha) > 0.87 ||\
      (VertexBuffer.visible &&\
       !sphere_in_frustum (x, y, z, _r, &(view)->frustum)))\
    continue;\
  if (is_leaf(cell) ||\
      (VertexBuffer.visible &&\
       sphere_diameter (x, y, z, _r/L0, &(view)->frustum) < (view)->res)) {\
    if (is_active(cell) && is_local(cell)) {\

#line 370

#define end_foreach_visible_plane()\
    }\
    continue;\
  }\
}\
end_foreach_cell();\

#line 377



static bool _reversed = false;

static void begin_draw_lines (bview * view, float color[3], float lw)
{
  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glTranslatef (0., 0., view->lc*view->fov/24.);
  vertex_buffer_glColor3f (color[0], color[1], color[2]);
  glLineWidth (view->samples*(lw > 0. ? lw : 1.));
  _reversed = view->reversed;
  view->reversed = false;
}

static void end_draw_lines()
{
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();
  bview * view = draw();
  view->reversed = _reversed;
}

static inline double interp (Point point, coord p, scalar col) { int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 401 "/home/davidtang/basilisk/src/draw.h"

  struct _interpolate _r = { col, x + p.x*Delta, y + p.y*Delta, z + p.z*Delta };
  return interpolate_linear (point, _r);
}

static double evaluate_expression (Point point, Node * n)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 407 "/home/davidtang/basilisk/src/draw.h"

  assert (n);
  switch (n->type) {
  case '1': return n->d.value;
  case '+': return (evaluate_expression (point, n->e[0]) +
      evaluate_expression(point, n->e[1]));
  case '-': return (evaluate_expression (point, n->e[0]) -
      evaluate_expression(point, n->e[1]));
  case '*': return (evaluate_expression (point, n->e[0]) *
      evaluate_expression(point, n->e[1]));
  case '/': return (evaluate_expression (point, n->e[0]) /
      evaluate_expression(point, n->e[1]));
  case '^': return pow (evaluate_expression (point, n->e[0]),
   evaluate_expression(point, n->e[1]));
  case '>': return (evaluate_expression (point, n->e[0]) >
      evaluate_expression(point, n->e[1]));
  case '<': return (evaluate_expression (point, n->e[0]) <
      evaluate_expression(point, n->e[1]));
  case 'L': return (evaluate_expression (point, n->e[0]) <=
      evaluate_expression(point, n->e[1]));
  case 'G': return (evaluate_expression (point, n->e[0]) >=
      evaluate_expression(point, n->e[1]));
  case '=': return (evaluate_expression (point, n->e[0]) ==
      evaluate_expression(point, n->e[1]));
  case 'i': return (evaluate_expression (point, n->e[0]) !=
      evaluate_expression(point, n->e[1]));
  case 'O': return (evaluate_expression (point, n->e[0]) ||
      evaluate_expression(point, n->e[1]));
  case 'A': return (evaluate_expression (point, n->e[0]) &&
      evaluate_expression(point, n->e[1]));
  case '?': return (evaluate_expression (point, n->e[0]) ?
      evaluate_expression(point, n->e[1]) :
      evaluate_expression(point, n->e[2]));
  case 'm': return - evaluate_expression (point, n->e[0]);
  case 'f': return n->d.func (evaluate_expression (point, n->e[0]));
  case 'v': {
    scalar s = {n->s};
    int k[3] = {0,0,0};
    for (int i = 0; i < 3; i++)
      if (n->e[i])
 k[i] = evaluate_expression (point, n->e[i]);
    return val(s,k[0],k[1],k[2]);
  }
  case 'D': return Delta;
  case 'x': return x;
  case 'y': return y;
  case 'z': return z;
  default:
    fprintf (ferr, "unknown operation type '%c'\n", n->type);
    assert (false);
  }
  return undefined;
}

static bool assemble_node (Node * n)
{
  if (n->type == 'v') {
    char * id = n->d.id;
    scalar s = lookup_field (id);
    if (s.i >= 0)
      n->s = s.i;
    else {
      n->s = -1;
      if (!strcmp (id, "Delta"))
 reset_node_type (n, 'D');
      else if (!strcmp (id, "x"))
 reset_node_type (n, 'x');
      else if (!strcmp (id, "y"))
 reset_node_type (n, 'y');
      else if (!strcmp (id, "z"))
 reset_node_type (n, 'z');
      else {
 typedef struct { char * name; double val; } Constant;
 static Constant constants[] = {
   {"pi", pi },
   {"nodata", nodata },
   {"HUGE", HUGE },
   { NULL },
 };
 Constant * p = constants;
 while (p->name) {
   if (!strcmp (p->name, id)) {
     reset_node_type (n, '1');
     n->d.value = p->val;
     break;
   }
   p++;
 }
 if (n->type == 'v') {
   fprintf (ferr, "unknown identifier '%s'\n", id);
   return false;
 }
      }
    }
  }
  for (int i = 0; i < 3; i++)
    if (n->e[i] && !assemble_node (n->e[i]))
      return false;
  return true;
}

static scalar compile_expression (char * expr, bool * isexpr)
{
  *isexpr = false;
  if (!expr)
    return (scalar){-1};

  bview * view = get_view();
  scalar s;
  if (view->cache && (s = get_cexpr (view->cache, expr)).i >= 0)
    return s;

  Node * node = parse_node (expr);
  if (node == NULL) {
    fprintf (ferr, "'%s': syntax error\n", expr);
    return (scalar){-1};
  }
  if (!assemble_node (node)) {
    free_node (node);
    return (scalar){-1};
  }
  if (node->type == 'v' && node->e[0] == NULL) {
    scalar s = {node->s};
    free_node (node);
    return s;
  }
  s = new_scalar("s");
  pfree (_attribute[s.i].name,__func__,__FILE__,__LINE__);
  _attribute[s.i].name = pstrdup (expr,__func__,__FILE__,__LINE__);
   { foreach(){

#line 536 "/home/davidtang/basilisk/src/draw.h"

    val(s,0,0,0) = evaluate_expression (point, node); } end_foreach(); }
  boundary (((scalar []){s,{-1}}));
  restriction (((scalar []){s,{-1}}));
  free_node (node);

  if (view->cache)
    view->cache = add_cexpr (view->cache, view->maxlen, expr, s);
  else
    *isexpr = true;
  return s;
}
#line 603 "/home/davidtang/basilisk/src/draw.h"
static void begin_colorized (float fc[3], bool constant_color,
        double cmap[127][3], bool use_texture)
{

  if (use_texture) {
    GLfloat texture[3*256];
    for (int i = 0; i < 256; i++) {
      color j = colormap_color (cmap, i/255., 0, 1);
      texture[3*i] = j.r/255.;
      texture[3*i + 1] = j.g/255.;
      texture[3*i + 2] = j.b/255.;
    }
    glTexImage1D (GL_TEXTURE_1D, 0, GL_RGB, 256,0, GL_RGB, GL_FLOAT, texture);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MIN_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_MAG_FILTER, GL_LINEAR);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_S, GL_CLAMP_TO_EDGE);
    glTexParameteri (GL_TEXTURE_1D, GL_TEXTURE_WRAP_T, GL_CLAMP_TO_EDGE);
    glEnable (GL_TEXTURE_1D);
  }
  if (constant_color)
    vertex_buffer_glColor3f (fc[0], fc[1], fc[2]);
}

static void end_colorized() {
  glDisable (GL_TEXTURE_1D);
}
#line 660 "/home/davidtang/basilisk/src/draw.h"
struct _draw_vof {
  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};







static bool cfilter (Point point, scalar c, double cmin)
{ int ig = 0; NOT_UNUSED(ig); int jg = 0; NOT_UNUSED(jg); int kg = 0; NOT_UNUSED(kg); POINT_VARIABLES; 
#line 682 "/home/davidtang/basilisk/src/draw.h"

  double cmin1 = 4.*cmin;
  if (val(c,0,0,0) <= cmin) {
    {
#line 685

      if (val(c,1,0,0) >= 1. - cmin1 || val(c,-1,0,0) >= 1. - cmin1)
 return true;
#line 685

      if (val(c,0,1,0) >= 1. - cmin1 || val(c,0,-1,0) >= 1. - cmin1)
 return true;
#line 685

      if (val(c,0,0,1) >= 1. - cmin1 || val(c,0,0,-1) >= 1. - cmin1)
 return true;}
    return false;
  }
  if (val(c,0,0,0) >= 1. - cmin) {
    {
#line 691

      if (val(c,1,0,0) <= cmin1 || val(c,-1,0,0) <= cmin1)
 return true;
#line 691

      if (val(c,0,1,0) <= cmin1 || val(c,0,-1,0) <= cmin1)
 return true;
#line 691

      if (val(c,0,0,1) <= cmin1 || val(c,0,0,-1) <= cmin1)
 return true;}
    return false;
  }
  int n = 0;
  double min = HUGE, max = - HUGE;
   { foreach_neighbor(1) {
    if (val(c,0,0,0) > cmin && val(c,0,0,0) < 1. - cmin && ++n >= (1 << 3))
      return true;
    if (val(c,0,0,0) > max) max = val(c,0,0,0);
    if (val(c,0,0,0) < min) min = val(c,0,0,0);
  } end_foreach_neighbor(); }
  return max - min > 0.5;
}

static void glvertex3d (bview * view, double x, double y, double z) {
  if (view->map) {
    coord p = {x, y, z};
    view->map (&p);
    vertex_buffer_glVertex3d (p.x, p.y, p.z);
  }
  else
    vertex_buffer_glVertex3d (x, y, z);
}
#line 739 "/home/davidtang/basilisk/src/draw.h"

bool draw_vof (struct _draw_vof p)
{ trace ("draw_vof", "/home/davidtang/basilisk/src/draw.h", 741);
  scalar c = lookup_field (p.c);
  if (c.i < 0) {
    fprintf (ferr, "draw_vof(): no field named '%s'\n", p.c);
    { bool _ret =  false; end_trace("draw_vof", "/home/davidtang/basilisk/src/draw.h", 745);  return _ret; }
  }
  vector s = lookup_vector (p.s);

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) { bool _ret =  false; end_trace("draw_vof", "/home/davidtang/basilisk/src/draw.h", 749);  return _ret; } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((3 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;

  double cmin = 1e-3;
#line 762 "/home/davidtang/basilisk/src/draw.h"
  bview * view = draw();
#line 840 "/home/davidtang/basilisk/src/draw.h"
  double larger =
    p.larger ? p.larger : p.edges || (p.color && !p.linear) ? 1. : 1.1;
  if (p.edges)
    { begin_draw_lines (view, p.lc, p.lw); {
       { foreach_visible (view){

#line 844 "/home/davidtang/basilisk/src/draw.h"

 if (cfilter (point, c, cmin)) {
   coord n = facet_normal (point, c, s);
   double alpha = plane_alpha (val(c,0,0,0), n);
   coord v[12];
   int m = facets (n, alpha, v, larger);
   if (m > 2) {
     vertex_buffer_glBegin (GL_LINE_LOOP);
     for (int i = 0; i < m; i++)
       glvertex3d (view,
     x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
     vertex_buffer_glEnd ();
     view->ni++;
   }
 } } end_foreach_visible(); }
    } end_draw_lines(); }
  else
    { begin_colorized (p.fc, !VertexBuffer.color || !p.color, cmap, !VertexBuffer.color && p.color && p.linear && col.i >= 0); {
       { foreach_visible (view){

#line 862 "/home/davidtang/basilisk/src/draw.h"

 if (cfilter (point, c, cmin)) {
   coord n = facet_normal (point, c, s);
   double alpha = plane_alpha (val(c,0,0,0), n);
   coord v[12];
   int m = facets (n, alpha, v, larger);
   if (m > 2) {
     vertex_buffer_glBegin (GL_POLYGON);
     for (int i = 0; i < m; i++) {
       if (p.linear) {
  if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, interp (point, v[i], col), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[i], col); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
       }
       else {
  if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
       }
       glnormal3d (view, n.x, n.y, n.z);
       glvertex3d (view,
     x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
     }
     vertex_buffer_glEnd ();
     view->ni++;
   }
 } } end_foreach_visible(); }
    } end_colorized(); }
#line 896 "/home/davidtang/basilisk/src/draw.h"
  if (p.expr) delete(((scalar []){col,{-1}}));
  { bool _ret =  true; end_trace("draw_vof", "/home/davidtang/basilisk/src/draw.h", 897);  return _ret; }
 end_trace("draw_vof", "/home/davidtang/basilisk/src/draw.h", 898); }
#line 909 "/home/davidtang/basilisk/src/draw.h"
struct _isoline {
  char * phi;
  double val;
  int n;


  char * c;
  char * s;
  bool edges;
  double larger;
  int filled;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};


bool isoline (struct _isoline p)
{ trace ("isoline", "/home/davidtang/basilisk/src/draw.h", 931);
#line 952 "/home/davidtang/basilisk/src/draw.h"
  assert (false);

  { bool _ret =  true; end_trace("isoline", "/home/davidtang/basilisk/src/draw.h", 954);  return _ret; }
 end_trace("isoline", "/home/davidtang/basilisk/src/draw.h", 955); }
#line 968 "/home/davidtang/basilisk/src/draw.h"
struct _cells {
  coord n;
  double alpha;
  float lc[3], lw;
};


bool cells (struct _cells p)
{ trace ("cells", "/home/davidtang/basilisk/src/draw.h", 976);
  bview * view = draw();
  { begin_draw_lines (view, p.lc, p.lw); {
#line 990 "/home/davidtang/basilisk/src/draw.h"
     { foreach_visible_plane (view, p.n, p.alpha){

#line 990 "/home/davidtang/basilisk/src/draw.h"
 {
      coord v[12];
      int m = facets (n, alpha, v, 1.);
      if (m > 2) {
 vertex_buffer_glBegin (GL_LINE_LOOP);
 for (int i = 0; i < m; i++)
   glvertex3d (view, x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
 vertex_buffer_glEnd ();
 view->ni++;
      }
    } } end_foreach_visible_plane(); }

  } end_draw_lines(); }
  { bool _ret =  true; end_trace("cells", "/home/davidtang/basilisk/src/draw.h", 1003);  return _ret; }
 end_trace("cells", "/home/davidtang/basilisk/src/draw.h", 1004); }






struct _vectors {
  char * u;
  double scale;
  float lc[3], lw;
};


bool vectors (struct _vectors p)
{ trace ("vectors", "/home/davidtang/basilisk/src/draw.h", 1019);
#line 1052 "/home/davidtang/basilisk/src/draw.h"
  fprintf (ferr, "vectors() is not implemented in 3D yet\n");

  { bool _ret =  true; end_trace("vectors", "/home/davidtang/basilisk/src/draw.h", 1054);  return _ret; }
 end_trace("vectors", "/home/davidtang/basilisk/src/draw.h", 1055); }
#line 1075 "/home/davidtang/basilisk/src/draw.h"
struct _squares {
  char * color;
  char * z;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3];
  bool expr;

  coord n;
  double alpha;
};


bool squares (struct _squares p)
{ trace ("squares", "/home/davidtang/basilisk/src/draw.h", 1090);
#line 1106 "/home/davidtang/basilisk/src/draw.h"
  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) { bool _ret =  false; end_trace("squares", "/home/davidtang/basilisk/src/draw.h", 1106);  return _ret; } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((3 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;
  scalar f = col;

  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  if (p.linear) {
    { begin_colorized (p.fc, !VertexBuffer.color || !p.color, cmap, !VertexBuffer.color && p.color && p.linear && col.i >= 0); {
#line 1166 "/home/davidtang/basilisk/src/draw.h"
       { foreach_visible_plane (view, p.n, p.alpha){

#line 1166 "/home/davidtang/basilisk/src/draw.h"

 if (val(f,0,0,0) != nodata) {
   coord v[12];
   int m = facets (n, alpha, v, 1.);
   if (m > 2) {
     coord c = {0,0,0};
     for (int i = 0; i < m; i++)
       {
#line 1173

  c.x += v[i].x/m;
#line 1173

  c.y += v[i].y/m;
#line 1173

  c.z += v[i].z/m;}
     vertex_buffer_glBegin (GL_TRIANGLE_FAN);
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, interp (point, c, f), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, c, f); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex3d (view, x + c.x*Delta, y + c.y*Delta, z + c.z*Delta);
     for (int i = 0; i < m; i++) {
       if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, interp (point, v[i], f), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[i], f); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
       glvertex3d (view,
     x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
     }
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, interp (point, v[0], f), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v[0], f); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
     glvertex3d (view,
   x + v[0].x*Delta, y + v[0].y*Delta, z + v[0].z*Delta);
     vertex_buffer_glEnd ();
     view->ni++;
   }
 } } end_foreach_visible_plane(); }

    } end_colorized(); }
  }
  else {
#line 1211 "/home/davidtang/basilisk/src/draw.h"
     { foreach_visible_plane (view, p.n, p.alpha){

#line 1211 "/home/davidtang/basilisk/src/draw.h"

      if (val(f,0,0,0) != nodata) {
 coord v[12];
 int m = facets (n, alpha, v, 1.);
 if (m > 2) {
   vertex_buffer_glBegin (GL_POLYGON);
   for (int i = 0; i < m; i++) {
     if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
     glvertex3d (view,
   x + v[i].x*Delta, y + v[i].y*Delta, z + v[i].z*Delta);
   }
   vertex_buffer_glEnd ();
   view->ni++;
 }
      } } end_foreach_visible_plane(); }

  }
  if (p.expr) delete (((scalar []){col,{-1}}));




  { bool _ret =  true; end_trace("squares", "/home/davidtang/basilisk/src/draw.h", 1233);  return _ret; }
 end_trace("squares", "/home/davidtang/basilisk/src/draw.h", 1234); }
#line 1245 "/home/davidtang/basilisk/src/draw.h"
struct _box {
  bool notics;
  float lc[3], lw;
};


bool box (struct _box p)
{ trace ("box", "/home/davidtang/basilisk/src/draw.h", 1252);
  bview * view = draw();
  { begin_draw_lines (view, p.lc, p.lw); {

    float height = 0.5*gl_StrokeHeight();
    float width = gl_StrokeWidth ('1'), scale = L0/(60.*width), length;
    float Z1 = 3 == 2 ? 0. : Z0;
    char label[80];

    glMatrixMode (GL_MODELVIEW);

    if (!p.notics) {
      int nt = 8;
      for (int i = 0; i <= nt; i++) {
 glPushMatrix();
 glTranslatef (X0 + i*L0/nt - height/2.*scale, Y0 - width/3.*scale, Z1);
 glRotatef (-90, 0, 0, 1);
 glScalef (scale, scale, 1.);
 sprintf (label, "%g", X0 + i*L0/nt);
 gl_StrokeString (label);
 glPopMatrix();

 glPushMatrix();
 sprintf (label, "%g", Y0 + i*L0/nt);
 length = gl_StrokeLength (label);
 glTranslatef (X0 - (length + width/3.)*scale,
        Y0 + i*L0/nt - height/2.*scale, Z1);
 glScalef (scale, scale, 1.);
 gl_StrokeString (label);
 glPopMatrix();


 glPushMatrix();
 sprintf (label, "%g", Z0 + i*L0/nt);
 length = gl_StrokeLength (label);
 glTranslatef (X0 - (length + width/3.)*scale,
        Y0, Z0 + i*L0/nt + height/2.*scale);
 glRotatef (-90, 1, 0, 0);
 glScalef (scale, scale, 1.);
 gl_StrokeString (label);
 glPopMatrix();

      }

      glPushMatrix();
      sprintf (label, "%g", X0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 + L0/2 - height*scale, Y0 - (length + 4.*width)*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("X");
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Y0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + 4.*width)*scale,
      Y0 + L0/2. - height*scale, Z1);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("Y");
      glPopMatrix();


      glPushMatrix();
      sprintf (label, "%g", Z0 + L0/2.);
      length = gl_StrokeLength (label);
      glTranslatef (X0 - (length + 4.*width)*scale,
      Y0, Z0 + L0/2. + height*scale);
      glRotatef (-90, 1, 0, 0);
      glScalef (2.*scale, 2.*scale, 1.);
      gl_StrokeString ("Z");
      glPopMatrix();

    }
#line 1338 "/home/davidtang/basilisk/src/draw.h"
     { foreach_level (0){

#line 1338 "/home/davidtang/basilisk/src/draw.h"
 {
      for (int i = -1; i <= 1; i += 2) {
 vertex_buffer_glBegin (GL_LINE_LOOP);
 glvertex3d (view, x - Delta_x/2., y - Delta_y/2., z + i*Delta/2.);
 glvertex3d (view, x + Delta_x/2., y - Delta_y/2., z + i*Delta/2.);
 glvertex3d (view, x + Delta_x/2., y + Delta_y/2., z + i*Delta/2.);
 glvertex3d (view, x - Delta_x/2., y + Delta_y/2., z + i*Delta/2.);
 vertex_buffer_glEnd ();
 view->ni++;
 vertex_buffer_glBegin (GL_LINES);
 for (int j = -1; j <= 1; j += 2) {
   glvertex3d (view, x + i*Delta/2., y + j*Delta/2., z - Delta/2.);
   glvertex3d (view, x + i*Delta/2., y + j*Delta/2., z + Delta/2.);
 }
 vertex_buffer_glEnd ();
 view->ni++;
      }
    } } end_foreach_level(); }

  } end_draw_lines(); }
  { bool _ret =  true; end_trace("box", "/home/davidtang/basilisk/src/draw.h", 1358);  return _ret; }
 end_trace("box", "/home/davidtang/basilisk/src/draw.h", 1359); }
#line 1372 "/home/davidtang/basilisk/src/draw.h"
struct _isosurface {
  char * f;
  double v;

  char * color;
  double min, max, spread;
  bool linear;
  colormap map;
  float fc[3], lc[3], lw;
  bool expr;
};


bool isosurface (struct _isosurface p)
{ trace ("isosurface", "/home/davidtang/basilisk/src/draw.h", 1386);

  if (!p.f)
    { bool _ret =  false; end_trace("isosurface", "/home/davidtang/basilisk/src/draw.h", 1389);  return _ret; }

  scalar f = {-1};
  bool fexpr;
  if (strcmp (p.f, "level")) {
    f = compile_expression (p.f, &fexpr);
    if (f.i < 0)
      { bool _ret =  false; end_trace("isosurface", "/home/davidtang/basilisk/src/draw.h", 1396);  return _ret; }
  }

  scalar col = {-1}; if (p.color && strcmp (p.color, "level")) { col = compile_expression (p.color, &p.expr); if (col.i < 0) { bool _ret =  false; end_trace("isosurface", "/home/davidtang/basilisk/src/draw.h", 1399);  return _ret; } } double cmap[127][3]; if (p.color) { if (p.min == 0 && p.max == 0) { if (col.i < 0) p.min = 0, p.max = depth(); else { stats s = statsf (col); double avg = s.sum/s.volume; if (p.spread < 0.) p.min = s.min, p.max = s.max; else { if (!p.spread) p.spread = 5.; double spread = p.spread*s.stddev; p.min = avg - spread; p.max = avg + spread; } } } if (!p.map) p.map = jet; p.map (cmap); } if ((3 > 2 || p.linear) && !p.fc[0] && !p.fc[1] && !p.fc[2]) p.fc[0] = p.fc[1] = p.fc[2] = 1.;;

  scalar v= new_vertex_scalar("v");
   { foreach_vertex(){

#line 1402 "/home/davidtang/basilisk/src/draw.h"

    val(v,0,0,0) = (val(f,0,0,0) + val(f,-1,0,0) + val(f,0,-1,0) + val(f,-1,-1,0) +
    val(f,0,0,-1) + val(f,-1,0,-1) + val(f,0,-1,-1) + val(f,-1,-1,-1))/8.; } end_foreach_vertex(); }
  boundary (((scalar []){v,{-1}}));

  vector n= new_vector("n");
   { foreach(){

#line 1408 "/home/davidtang/basilisk/src/draw.h"

    {
#line 1409

      val(n.x,0,0,0) = ((val(f,1,0,0) - val(f,-1,0,0))/(2.*Delta));
#line 1409

      val(n.y,0,0,0) = ((val(f,0,1,0) - val(f,0,-1,0))/(2.*Delta));
#line 1409

      val(n.z,0,0,0) = ((val(f,0,0,1) - val(f,0,0,-1))/(2.*Delta));}; } end_foreach(); }
  boundary ((scalar *)((vector []){{n.x,n.y,n.z},{{-1},{-1},{-1}}}));

  bview * view = draw();
  glShadeModel (GL_SMOOTH);
  { begin_colorized (p.fc, !VertexBuffer.color || !p.color, cmap, !VertexBuffer.color && p.color && p.linear && col.i >= 0); {
     { foreach_visible (view){

#line 1416 "/home/davidtang/basilisk/src/draw.h"
 {
      double val[8] = {
 val(v,0,0,0), val(v,1,0,0), val(v,1,0,1), val(v,0,0,1),
 val(v,0,1,0), val(v,1,1,0), val(v,1,1,1), val(v,0,1,1)
      };
      double t[5][3][3];
      int nt = polygonize (val, p.v, t);
      for (int i = 0; i < nt; i++) {
 if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
 vertex_buffer_glBegin (GL_POLYGON);
 for (int j = 0; j < 3; j++) {
   coord v = {t[i][j][0], t[i][j][1], t[i][j][2]}, np;
   {
#line 1428

     np.x = interp (point, v, n.x);
#line 1428

     np.y = interp (point, v, n.y);
#line 1428

     np.z = interp (point, v, n.z);}
   glnormal3d (view, np.x, np.y, np.z);
   if (p.linear) {
     if (p.color && p.linear && col.i >= 0) { if (VertexBuffer.color) { color b = colormap_color (cmap, interp (point, v, col), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); } else { double _v = interp (point, v, col); glTexCoord1d (clamp(((_v) - p.min)/(p.max - p.min), 0., 1.)); } };
   }
   else {
     if (p.color && (!p.linear || col.i < 0)) { color b = colormap_color (cmap, col.i < 0 ? (double) level : val(col,0,0,0), p.min, p.max); vertex_buffer_glColor3f (b.r/255., b.g/255., b.b/255.); };
   }
   glvertex3d (view, x + v.x*Delta_x, y + v.y*Delta_y, z + v.z*Delta_z);
 }
 vertex_buffer_glEnd ();
 view->ni++;
      }
    } } end_foreach_visible(); }
  } end_colorized(); }
  if (p.expr) delete (((scalar []){col,{-1}}));
  if (fexpr) delete (((scalar []){f,{-1}}));

  { bool _ret =  true; delete (((scalar []){n.x,n.y,n.z,v,{-1}}));  end_trace("isosurface", "/home/davidtang/basilisk/src/draw.h", 1447);  return _ret; }
 delete (((scalar []){n.x,n.y,n.z,v,{-1}}));  end_trace("isosurface", "/home/davidtang/basilisk/src/draw.h", 1448); }
#line 1458 "/home/davidtang/basilisk/src/draw.h"
struct _travelling {
  double start, end;
  float tx, ty, quat[4], fov;
};




void travelling (struct _travelling p)
{
  static float tx, ty, quat[4], fov;
  static double told = -1.;
  if (told < p.start && t >= p.start) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  if (t >= p.start && t <= p.end)
    view ((struct _view_set){.tx = (!p.tx ? tx : ((t - p.start)*(p.tx) + (p.end - t)*(tx))/(p.end - p.start)), .ty = (!p.ty ? ty : ((t - p.start)*(p.ty) + (p.end - t)*(ty))/(p.end - p.start)),
   .fov = (!p.fov ? fov : ((t - p.start)*(p.fov) + (p.end - t)*(fov))/(p.end - p.start)),
   .quat = {(!p.quat[0] ? quat[0] : ((t - p.start)*(p.quat[0]) + (p.end - t)*(quat[0]))/(p.end - p.start)), (!p.quat[1] ? quat[1] : ((t - p.start)*(p.quat[1]) + (p.end - t)*(quat[1]))/(p.end - p.start)),
           (!p.quat[2] ? quat[2] : ((t - p.start)*(p.quat[2]) + (p.end - t)*(quat[2]))/(p.end - p.start)), (!p.quat[3] ? quat[3] : ((t - p.start)*(p.quat[3]) + (p.end - t)*(quat[3]))/(p.end - p.start))}});
  if (told < p.end && t >= p.end) {
    bview * view = get_view();
    tx = view->tx, ty = view->ty, fov = view->fov;
    for (int i = 0; i < 4; i++)
      quat[i] = view->quat[i];
  }
  told = t;
}
#line 1505 "/home/davidtang/basilisk/src/draw.h"
struct _draw_string {
  char * str;
  int pos;
  float size;
  float lc[3], lw;
};


bool draw_string (struct _draw_string p)
{ trace ("draw_string", "/home/davidtang/basilisk/src/draw.h", 1514);
  bview * view = draw();

  glMatrixMode (GL_PROJECTION);
  glPushMatrix();
  glLoadIdentity();

  glMatrixMode (GL_MODELVIEW);
  glPushMatrix();
  glLoadIdentity();

  vertex_buffer_glColor3f (p.lc[0], p.lc[1], p.lc[2]);
  glLineWidth (view->samples*(p.lw > 0. ? p.lw : 1.));

  float width = gl_StrokeWidth ('1'), height = gl_StrokeHeight();
  if (!p.size)
    p.size = 40;
  float hscale = 2./(p.size*width), vscale = hscale*view->width/view->height;
  float vmargin = width/2.*vscale;
  if (p.pos == 0)
    glTranslatef (-1., -1. + vmargin, 0.);
  else if (p.pos == 1)
    glTranslatef (-1., 1. - height*vscale, 0.);
  else if (p.pos == 2)
    glTranslatef (1. - strlen(p.str)*width*hscale, 1. - height*vscale, 0.);
  else
    glTranslatef (1. - strlen(p.str)*width*hscale, -1. + vmargin, 0.);
  glScalef (hscale, vscale, 1.);
  gl_StrokeString (p.str);

  glMatrixMode (GL_MODELVIEW);
  glPopMatrix();
  glMatrixMode (GL_PROJECTION);
  glPopMatrix();

  { bool _ret =  true; end_trace("draw_string", "/home/davidtang/basilisk/src/draw.h", 1549);  return _ret; }
 end_trace("draw_string", "/home/davidtang/basilisk/src/draw.h", 1550); }




struct _labels {
  char * f;
  float lc[3], lw;
};


bool labels (struct _labels p)
{ trace ("labels", "/home/davidtang/basilisk/src/draw.h", 1562);
#line 1591 "/home/davidtang/basilisk/src/draw.h"
  fprintf (ferr, "labels() is not implemented in 3D yet\n");
  { bool _ret =  false; end_trace("labels", "/home/davidtang/basilisk/src/draw.h", 1592);  return _ret; }

 end_trace("labels", "/home/davidtang/basilisk/src/draw.h", 1594); }







#line 1 "draw_json.h"
#line 1 "/home/davidtang/basilisk/src/draw_json.h"



int _view_set_json (void * q, char * s, int len) {
  struct _view_set * p = (struct _view_set *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"view_set\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->ty);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->fov);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", p->quat[0], p->quat[1], p->quat[2], p->quat[3]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sy\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sy);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"sz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->sz);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"width\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->width);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"height\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->height);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"samples\": { \"type\": \"punsigned\", \"cardinality\": 1, \"value\": \"%d\" }", p->samples);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"bg\": { \"type\": \"pfloat\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->bg[0], p->bg[1], p->bg[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"theta\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->theta);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"phi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->phi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"psi\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->psi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"relative\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->relative);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tz\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tz);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"near\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->near);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"far\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->far);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"res\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->res);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"camera\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->camera);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"cache\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->cache);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p1x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p1y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p1y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p2x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"p2y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->p2y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _translate_json (void * q, char * s, int len) {
  struct _translate * p = (struct _translate *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"translate\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"x\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->x);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"y\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->y);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _mirror_json (void * q, char * s, int len) {
  struct _mirror * p = (struct _mirror *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"mirror\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _draw_vof_json (void * q, char * s, int len) {
  struct _draw_vof * p = (struct _draw_vof *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_vof\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->c);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->s);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->edges);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->larger);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->filled);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _isoline_json (void * q, char * s, int len) {
  struct _isoline * p = (struct _isoline *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isoline\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"phi\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->phi);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"val\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->val);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->n);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"c\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->c);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"s\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->s);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"edges\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->edges);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"larger\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->larger);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"filled\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->filled);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _cells_json (void * q, char * s, int len) {
  struct _cells * p = (struct _cells *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"cells\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _vectors_json (void * q, char * s, int len) {
  struct _vectors * p = (struct _vectors *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"vectors\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"u\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->u);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"scale\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->scale);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _squares_json (void * q, char * s, int len) {
  struct _squares * p = (struct _squares *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"squares\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"z\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"n\": { \"type\": \"coord\", \"cardinality\": 3, \"value\": [%lf,%lf,%lf] }", p->n.x, p->n.y, p->n.z);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"alpha\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->alpha);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _box_json (void * q, char * s, int len) {
  struct _box * p = (struct _box *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"box\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"notics\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->notics);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _isosurface_json (void * q, char * s, int len) {
  struct _isosurface * p = (struct _isosurface *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"isosurface\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->f);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"v\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->v);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"color\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->color);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"min\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->min);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"max\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->max);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"spread\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->spread);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"linear\": { \"type\": \"pbool\", \"cardinality\": 1, \"value\": \"%d\" }", p->linear);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->fc[0], p->fc[1], p->fc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _travelling_json (void * q, char * s, int len) {
  struct _travelling * p = (struct _travelling *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"travelling\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"start\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->start);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"end\": { \"type\": \"pdouble\", \"cardinality\": 1, \"value\": \"%lf\" }", p->end);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"tx\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->tx);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"ty\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->ty);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"quat\": { \"type\": \"pfloat\", \"cardinality\": 4, \"value\": [%f,%f,%f,%f] }", p->quat[0], p->quat[1], p->quat[2], p->quat[3]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"fov\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->fov);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _draw_string_json (void * q, char * s, int len) {
  struct _draw_string * p = (struct _draw_string *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"draw_string\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"str\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->str);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"pos\": { \"type\": \"pint\", \"cardinality\": 1, \"value\": \"%d\" }", p->pos);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"size\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->size);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}

int _labels_json (void * q, char * s, int len) {
  struct _labels * p = (struct _labels *) q;
  int i, len1 = 0;
  i = snprintf (s, len, "  \"labels\" : {");
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n    \"f\": { \"type\": \"pstring\", \"cardinality\": 1, \"value\": \"%s\" }", p->f);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lc\": { \"type\": \"color\", \"cardinality\": 3, \"value\": [%f,%f,%f] }", p->lc[0], p->lc[1], p->lc[2]);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, ",\n    \"lw\": { \"type\": \"pfloat\", \"cardinality\": 1, \"value\": \"%f\" }", p->lw);
  s += i, len -= i, len1 += i;
  i = snprintf (s, len, "\n  }");
  s += i, len -= i, len1 += i;
  return len1;
}
#line 1603 "/home/davidtang/basilisk/src/draw.h"

struct {
  int (* json) (void * q, char * s, int len);
} bview_interface[] = {
  { _draw_vof_json },
  { _squares_json },
  { _cells_json },
  { _box_json },





  { _isosurface_json },

  { NULL }
};
#line 427 "/home/davidtang/basilisk/src/view.h"
#line 442 "/home/davidtang/basilisk/src/view.h"
struct _load {
  FILE * fp;
  char * file;
  Array * buf;
};

bool load (struct _load p);
#line 496 "/home/davidtang/basilisk/src/view.h"
struct _save {
  char * file, * format, * opt;
  FILE * fp;
  float lw;
  int sort, options;

  bview * view;
};

static void bview_draw (bview * view)
{
  if (!view->active)
    return;
  view->active = false;
  glFinish ();
  enable_fpe (FE_DIVBYZERO|FE_INVALID);
}


bool save (struct _save p)
{ trace ("save", "/home/davidtang/basilisk/src/view.h", 516);
  char ppm[] = "ppm";
  if (!p.format) {
    p.format = ppm;
    if (p.file) {
      char * s = strchr (p.file, '.'), * dot = s;
      while (s) {
 dot = s;
 s = strchr (s + 1, '.');
      }
      if (dot)
 p.format = dot + 1;
    }
  }

  bview * view = p.view ? p.view : get_view();

  if (!strcmp (p.format, "png") ||
      !strcmp (p.format, "jpg") ||
      (p.file && is_animation (p.file))) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0) {
      FILE * fp = open_image (p.file, p.opt);
      if (!fp) {
 perror (p.file);
 { bool _ret =  false; end_trace("save", "/home/davidtang/basilisk/src/view.h", 542);  return _ret; }
      }
      gl_write_image (fp, image, view->width, view->height, view->samples);
      close_image (p.file, fp);
    }
    { bool _ret =  true; end_trace("save", "/home/davidtang/basilisk/src/view.h", 547);  return _ret; }
  }

  if (p.file && (p.fp = fopen (p.file, "w")) == NULL) {
    perror (p.file);
    { bool _ret =  false; end_trace("save", "/home/davidtang/basilisk/src/view.h", 552);  return _ret; }
  }
  if (!p.fp)
    p.fp = fout;

  if (!strcmp (p.format, "ppm")) {
    bview_draw (view);
    unsigned char * image = (unsigned char *) compose_image (view);
    if (pid() == 0)
      gl_write_image (p.fp, image, view->width, view->height, view->samples);
  }

  else if (!strcmp (p.format, "bv")) {

    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      p.format);
#line 584 "/home/davidtang/basilisk/src/view.h"
  }

  else if (!strcmp (p.format, "gnu") ||
    !strcmp (p.format, "obj") ||
    !strcmp (p.format, "kml") ||
    !strcmp (p.format, "ps") ||
    !strcmp (p.format, "eps") ||
    !strcmp (p.format, "tex") ||
    !strcmp (p.format, "pdf") ||
    !strcmp (p.format, "svg") ||
    !strcmp (p.format, "pgf"))
    fprintf (ferr, "save(): error: the '%s' format is no longer supported\n",
      p.format);

  else {
    fprintf (ferr, "save(): unknown format '%s'\n", p.format);
    if (p.file) {
      fclose (p.fp);
      remove (p.file);
    }
    { bool _ret =  false; end_trace("save", "/home/davidtang/basilisk/src/view.h", 604);  return _ret; }
  }

  fflush (p.fp);
  if (p.file)
    fclose (p.fp);

  { bool _ret =  true; end_trace("save", "/home/davidtang/basilisk/src/view.h", 611);  return _ret; }
 end_trace("save", "/home/davidtang/basilisk/src/view.h", 612); }







static char * remove_blanks (char * line)
{
  while (strchr (" \t", *line)) line++;
  char * s = line, * cur = line;
  bool instring = false;
  while (*s != '\0' && *s != '#') {
    if (*s == '"')
      instring = !instring;
    if (instring || !strchr (" \t", *s))
      *cur++ = *s;
    s++;
  }
  *cur = '\0';
  return line;
}






#line 1 "draw_get.h"
#line 1 "/home/davidtang/basilisk/src/draw_get.h"

#line 1 "parse.h"
#line 1 "/home/davidtang/basilisk/src/parse.h"




enum ParamsType { pstring, pint, punsigned, pbool, pfloat, pdouble, pcolormap };

typedef struct {
  char * key;
  enum ParamsType type;
  void * val;
  int n;
} Params;

static bool atobool (char * s)
{
  if (!strcmp (s, "true"))
    return true;
  if (!strcmp (s, "false"))
    return false;
  return atoi (s) != 0;
}

static bool args (Params * p, char * val)
{
  static char * name[] = { "string", "int", "unsigned",
      "bool", "float", "double", "colormap" };
  switch (p->type) {

  case pstring:
    if (val[0] != '"') {
      fprintf (ferr, "expecting a string for '%s' got '%s'\n", p->key, val);
      return false;
    }
    if (val[strlen(val) - 1] != '"') {
      fprintf (ferr, "unterminated quoted string '%s'\n", val);
      return false;
    }
    val[strlen(val) - 1] = '\0';
    char * s = &val[1];
    int nc = 0;
    while (*s != '\0') {
      if (!strchr (" \t\n\r", *s))
 nc++;
      s++;
    }
    *((char **)p->val) = nc > 0 ? &val[1] : NULL;
    break;

  case pcolormap:
    if (!strcmp (val, "jet"))
      *((colormap *)p->val) = jet;
    else if (!strcmp (val, "cool_warm"))
      *((colormap *)p->val) = cool_warm;
    else if (!strcmp (val, "gray"))
      *((colormap *)p->val) = gray;
    else if (!strcmp (val, "randomap"))
      *((colormap *)p->val) = randomap;
    else {
      fprintf (ferr, "unknown colormap '%s'\n", val);
      return false;
    }
    break;

  case pint: case punsigned: case pbool: case pdouble: case pfloat:
    if (val[0] == '"') {
      fprintf (ferr, "expecting a %s for '%s' got %s\n",
        name[p->type], p->key, val);
      return false;
    }
    if (!p->n) {
      switch (p->type) {
      case pint: *((int *)p->val) = atoi(val); break;
      case punsigned: *((unsigned *)p->val) = atoi(val); break;
      case pbool: *((bool *)p->val) = atobool(val); break;
      case pfloat: *((float *)p->val) = atof(val); break;
      case pdouble: *((double *)p->val) = atof(val); break;
      default: assert (false);
      }
    }
    else {
      if (val[0] != '{') {
 fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
 return false;
      }
      val++;
      int i = 0;
      char c = ',';
      while (i < p->n && c != '}') {
 char * s = strchr (val, ',');
 if (!s)
   s = strchr (val, '}');
 if (!s) {
   fprintf (ferr, "expecting an array for '%s' got %s\n", p->key, val);
   return false;
 }
 c = *s;
 *s++ = '\0';
 switch (p->type) {
 case pint: ((int *)p->val)[i++] = atoi (val); break;
 case punsigned: ((unsigned *)p->val)[i++] = atoi (val); break;
 case pbool: ((bool *)p->val)[i++] = atobool (val); break;
 case pfloat: ((float *)p->val)[i++] = atof (val); break;
 case pdouble: ((double *)p->val)[i++] = atof (val); break;
 default: assert (false);
 }
 val = s;
      }
      if (c != '}') {
 fprintf (ferr, "expecting '}' for '%s' got %s\n", p->key, val);
 return false;
      }
    }
    break;

  default:
    assert (false);
  }
  return true;
}

static char * find_comma (char * s)
{
  int par = 0;
  while (*s != '\0') {
    if (*s == ',' && par == 0) {
      *s = '\0';
      return s + 1;
    }
    if (*s == '{')
      par++;
    else if (*s == '}')
      par--;
    s++;
  }
  return NULL;
}

static char * mystrtok (char * str, const char * delim)
{
  static char * s = NULL;
  char * start = str ? str : s;
  bool string = false;
  s = start;
  while (*s != '\0') {
    if (*s == '"')
      string = !string;
    if (!string && strchr(delim, *s))
      break;
    s++;
  }
  if (*s != '\0')
    *s++ = '\0';
  return start;
}

int parse_params (Params * params)
{
  char * s;
  int i = 0, n = 0;
  Params * p = params;
  while (p->key) p++, n++;
  if (!(s = mystrtok (NULL, ");")) || s[0] == '\n')
    return false;
  while (s && *s != '\0') {
    char * next = find_comma (s), * key = s;
    if ((s = strchr (key, '='))) {
      s[0] = '\0', s++;
      i = -1;
      Params * p = params;
      while (p->key && strcmp(p->key, key)) p++;
      if (!p->key) {
 fprintf (ferr, "unknown key '%s'\n", key);
 return false;
      }
      if (!args (p, s))
 return false;
    }
    else {
      if (i < 0) {
 fprintf (ferr, "anonymous value '%s' after keys\n", key);
 return false;
      }
      if (i >= n) {
 fprintf (ferr, "too many parameters: '%s' %d %d\n", key, i, n);
 return false;
      }
      if (!args (&params[i], key))
 return false;
      i++;
    }
    s = next;
  }
  return true;
}
#line 3 "/home/davidtang/basilisk/src/draw_get.h"

bool _view_set_get (struct _view_set * p) {
  Params params[] = {
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"fov", pfloat, &p->fov},
    {"quat", pfloat, p->quat, 4},
    {"sx", pfloat, &p->sx},
    {"sy", pfloat, &p->sy},
    {"sz", pfloat, &p->sz},
    {"width", punsigned, &p->width},
    {"height", punsigned, &p->height},
    {"samples", punsigned, &p->samples},
    {"bg", pfloat, p->bg, 3},
    {"theta", pfloat, &p->theta},
    {"phi", pfloat, &p->phi},
    {"psi", pfloat, &p->psi},
    {"relative", pbool, &p->relative},
    {"tz", pfloat, &p->tz},
    {"near", pfloat, &p->near},
    {"far", pfloat, &p->far},
    {"res", pfloat, &p->res},
    {"camera", pstring, &p->camera},
    {"cache", pint, &p->cache},
    {"p1x", pfloat, &p->p1x},
    {"p1y", pfloat, &p->p1y},
    {"p2x", pfloat, &p->p2x},
    {"p2y", pfloat, &p->p2y},
    {NULL}
  };
  return parse_params (params);
}

bool _translate_get (struct _translate * p) {
  Params params[] = {
    {"x", pfloat, &p->x},
    {"y", pfloat, &p->y},
    {"z", pfloat, &p->z},
    {NULL}
  };
  return parse_params (params);
}

bool _mirror_get (struct _mirror * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_vof_get (struct _draw_vof * p) {
  Params params[] = {
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _isoline_get (struct _isoline * p) {
  Params params[] = {
    {"phi", pstring, &p->phi},
    {"val", pdouble, &p->val},
    {"n", pint, &p->n},
    {"c", pstring, &p->c},
    {"s", pstring, &p->s},
    {"edges", pbool, &p->edges},
    {"larger", pdouble, &p->larger},
    {"filled", pint, &p->filled},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _cells_get (struct _cells * p) {
  Params params[] = {
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _vectors_get (struct _vectors * p) {
  Params params[] = {
    {"u", pstring, &p->u},
    {"scale", pdouble, &p->scale},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _squares_get (struct _squares * p) {
  Params params[] = {
    {"color", pstring, &p->color},
    {"z", pstring, &p->z},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"expr", pbool, &p->expr},
    {"n", pdouble, &p->n, 3},
    {"alpha", pdouble, &p->alpha},
    {NULL}
  };
  return parse_params (params);
}

bool _box_get (struct _box * p) {
  Params params[] = {
    {"notics", pbool, &p->notics},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _isosurface_get (struct _isosurface * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"v", pdouble, &p->v},
    {"color", pstring, &p->color},
    {"min", pdouble, &p->min},
    {"max", pdouble, &p->max},
    {"spread", pdouble, &p->spread},
    {"linear", pbool, &p->linear},
    {"map", pcolormap, &p->map},
    {"fc", pfloat, p->fc, 3},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {"expr", pbool, &p->expr},
    {NULL}
  };
  return parse_params (params);
}

bool _travelling_get (struct _travelling * p) {
  Params params[] = {
    {"start", pdouble, &p->start},
    {"end", pdouble, &p->end},
    {"tx", pfloat, &p->tx},
    {"ty", pfloat, &p->ty},
    {"quat", pfloat, p->quat, 4},
    {"fov", pfloat, &p->fov},
    {NULL}
  };
  return parse_params (params);
}

bool _draw_string_get (struct _draw_string * p) {
  Params params[] = {
    {"str", pstring, &p->str},
    {"pos", pint, &p->pos},
    {"size", pfloat, &p->size},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}

bool _labels_get (struct _labels * p) {
  Params params[] = {
    {"f", pstring, &p->f},
    {"lc", pfloat, p->lc, 3},
    {"lw", pfloat, &p->lw},
    {NULL}
  };
  return parse_params (params);
}
#line 642 "/home/davidtang/basilisk/src/view.h"

bool process_line (char * line)
{
  if (line[0] == '\0')
    return true;
  char * buf = pstrdup (line,__func__,__FILE__,__LINE__);
  char * s = mystrtok (remove_blanks (line), "(");
  if (!s) {
    pfree (buf,__func__,__FILE__,__LINE__);
    return true;
  }

  if (!strcmp (s, "restore")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      bview * view = get_view();
      if (view->cache) {
 free_cexpr (view->cache);
 view->cache = pcalloc (1, sizeof (cexpr),__func__,__FILE__,__LINE__);
      }
      if (!restore ((struct Dump){.file = file, .list = all}))
 fprintf (ferr, "could not restore from '%s'\n", file);
      else {
 restriction (all);
 fields_stats();
 clear();
      }
    }
  }

  else if (!strcmp (s, "dump")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    dump ((struct Dump){.file = file});
  }

  else if (!strcmp (s, "input_gfs")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file) {
      input_gfs ((struct OutputGfs){.file = file, .list = all});
      restriction (all);
      fields_stats();
      clear();
    }
  }

  else if (!strcmp (s, "save")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      save ((struct _save){.file = file});
  }

  else if (!strcmp (s, "load")) {
    char * file = NULL;
    parse_params ((Params[]){{"file", pstring, &file}, {NULL}});
    if (file)
      load ((struct _load){.file = file});
  }

  else if (!strcmp (s, "cells")) {
    struct _cells p = {{0}};
    if (!_cells_get (&p) || !cells (p))
      return false;
  }

  else if (!strcmp (s, "vectors")) {
    struct _vectors p = {0};
    if (!_vectors_get (&p) || !vectors (p))
      return false;
  }

  else if (!strcmp (s, "draw_vof")) {
    struct _draw_vof p = {0};
    if (!_draw_vof_get (&p) || !draw_vof (p))
      return false;
  }

  else if (!strcmp (s, "isoline")) {
    struct _isoline p = {0};
    if (!_isoline_get (&p) || !isoline (p))
      return false;
  }

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    if (!_squares_get (&p) || !squares (p))
      return false;
  }

  else if (!strcmp (s, "begin_translate")) {
    struct _translate p = {0};
    _translate_get (&p);
    begin_translate (p);
  }

  else if (!strcmp (s, "end_translate"))
    end_translate();

  else if (!strcmp (s, "begin_mirror")) {
    struct _mirror p = {{0}};
    _mirror_get (&p);
    begin_mirror (p);
  }

  else if (!strcmp (s, "end_mirror"))
    end_mirror();

  else if (!strcmp (s, "squares")) {
    struct _squares p = {0};
    if (!_squares_get (&p) || !squares (p))
      return false;
  }

  else if (!strcmp (s, "isosurface")) {
    struct _isosurface p = {0};
    if (!_isosurface_get (&p) || !isosurface (p))
      return false;
  }

  else if (!strcmp (s, "draw_string")) {
    struct _draw_string p = {0};
    if (!_draw_string_get (&p) || !draw_string (p))
      return false;
  }

  else if (!strcmp (s, "labels")) {
    struct _labels p = {0};
    if (!_labels_get (&p) || !labels (p))
      return false;
  }

  else if (!strcmp (s, "clear"))
    clear();

  else if (!strcmp (s, "box")) {
    struct _box p = {0};
    if (!_box_get (&p) || !box (p))
      return false;
  }

  else if (!strcmp (s, "view")) {
    struct _view_set p = {0};
    _view_set_get (&p);
    view (p);
  }

  else if (s[0] != '\n' && s[0] != '\0')
    fprintf (ferr, "load(): syntax error: '%s'\n", s);

  pfree (buf,__func__,__FILE__,__LINE__);
  return true;
}

bool load (struct _load p) {
  if (p.file) {
    p.fp = fopen (p.file, "r");
    if (!p.fp) {
      perror (p.file);
      return false;
    }
  }

  if (p.fp) {
    char line[256];
    while (fgets (line, 256, p.fp) && process_line (line));
  }
  else if (p.buf) {
    int i = 0;
    char * s = (char *) p.buf->p;
    while (i < p.buf->len) {
      char * start = s;
      while (i < p.buf->len && *s != '\n')
 s++, i++;
      if (*s == '\n' && ++s > start) {
 char line[s - start + 1];
 strncpy (line, start, s - start);
 line[s - start] = '\0';
 process_line (line);
      }
    }
  }
  return true;
}
#line 22 "isotropic.c"




#line 1 "navier-stokes/perfs.h"
#line 1 "/home/davidtang/basilisk/src/navier-stokes/perfs.h"






static int perfs_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i += 1);   *ip = i; *tp = t;   return ret; } static int perfs (const int i, const double t, Event * _ev) { trace ("perfs", "/home/davidtang/basilisk/src/navier-stokes/perfs.h", 7);  {
  static FILE * fp =NULL; if (!fp || i == 0) fp = pid() > 0 ? fopen("/dev/null", "w") :  fopen ("perfs", "w");
  if (i == 0)
    fprintf (fp,
      "t dt mgp.i mgp.nrelax mgpf.i mgpf.nrelax mgu.i mgu.nrelax "
      "grid->tn perf.t perf.speed npe\n");
  fprintf (fp, "%g %g %d %d %d %d %d %d %ld %g %g %d\n",
    t, dt, mgp.i, mgp.nrelax, mgpf.i, mgpf.nrelax, mgu.i, mgu.nrelax,
    grid->tn, perf.t, perf.speed, npe());
  fflush (fp);
 end_trace("perfs", "/home/davidtang/basilisk/src/navier-stokes/perfs.h", 17); } return 0; } 






static int perf_plot_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 10);   *ip = i; *tp = t;   return ret; } static int perf_plot (const int i, const double t, Event * _ev) { trace ("perf_plot", "/home/davidtang/basilisk/src/navier-stokes/perfs.h", 24);  {
  if (getenv ("DISPLAY"))
    qpopen ("gnuplot -e 'set term x11 noraise title perfs' "
    "$BASILISK/navier-stokes/perfs.plot "
    "& read dummy; kill $!", "w");
 end_trace("perf_plot", "/home/davidtang/basilisk/src/navier-stokes/perfs.h", 29); } return 0; } 
#line 27 "isotropic.c"
#line 1 "maxruntime.h"
#line 1 "/home/davidtang/basilisk/src/maxruntime.h"
#line 18 "/home/davidtang/basilisk/src/maxruntime.h"
static double _maxruntime = HUGE;

static int runtime_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i += 10);   *ip = i; *tp = t;   return ret; } static int runtime (const int i, const double t, Event * _ev) { trace ("runtime", "/home/davidtang/basilisk/src/maxruntime.h", 20);  {
  mpi_all_reduce (perf.t, MPI_DOUBLE, MPI_MAX);
  if (perf.t >= _maxruntime - 300) {
    dump ((struct Dump){.file = "restart"});
    { int _ret =  1; end_trace("runtime", "/home/davidtang/basilisk/src/maxruntime.h", 24);  return _ret; }
  }
 end_trace("runtime", "/home/davidtang/basilisk/src/maxruntime.h", 26); } return 0; } 

void maxruntime (int * argc, char * argv[])
{
  for (int i = 0; i < *argc; i++)
    if (!strcmp (argv[i], "--maxruntime") || !strcmp (argv[i], "-m")) {
      if (i + 1 < *argc) {
 char * s = strtok (argv[i + 1], ":");
 int n = 0;
 _maxruntime = 0;
 do {
   _maxruntime = 60*_maxruntime + atoi(s);
   n++;
 } while ((s = strtok (NULL, ":")));
 if (n > 3) {
   fprintf (ferr, "maxruntime: TIME format must be H:M:S\n");
   exit (1);
 }
      }
      else {
 fprintf (ferr, "usage: %s TIME\n", argv[i]);
 exit (1);
      }
      *argc -= 2;
      for (int j = i; j < *argc; j++)
 argv[j] = argv[j + 2];
    }
}
#line 28 "isotropic.c"






vector av= {{11},{12},{13}};





int maxlevel = 7;

int main (int argc, char * argv[])
{ _init_solver();
  maxruntime (&argc, argv);
  if (argc > 1)
    maxlevel = atoi(argv[1]);




  L0 = 2.*pi;
  {
#line 52

    periodic (right);
#line 52

    periodic (top);
#line 52

    periodic (front);}





  vector muc= new_const_vector("muc", 9, (double []) {0.01,0.01,0.01});
  mu = muc;
  a = av;
  N = 1 << maxlevel;
  run();
 free_solver(); }







static int init_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i = 0);   *ip = i; *tp = t;   return ret; } static int init_0 (const int i, const double t, Event * _ev) { trace ("init_0", "isotropic.c", 72);  {
  if (!restore ((struct Dump){.file = "restart"}))
     { foreach(){

#line 74 "isotropic.c"
 {
      val(u.x,0,0,0) = cos(y) + sin(z);
      val(u.y,0,0,0) = sin(x) + cos(z);
      val(u.z,0,0,0) = cos(x) + sin(y);
    } } end_foreach(); }
 end_trace("init_0", "isotropic.c", 79); } return 0; } 







static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int acceleration_0 (const int i, const double t, Event * _ev) { trace ("acceleration_0", "isotropic.c", 87);  {
  coord ubar;
  {
#line 89
 {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
#line 89
 {
    stats s = statsf(u.y);
    ubar.y = s.sum/s.volume;
  }
#line 89
 {
    stats s = statsf(u.z);
    ubar.z = s.sum/s.volume;
  }}
   { foreach_face_generic() { int ig = -1; VARIABLES;  if (is_face_x()) {
#line 93
{

#line 93 "isotropic.c"

    val(av.x,0,0,0) += 0.1*((val(u.x,0,0,0) + val(u.x,-1,0,0))/2. - ubar.x); }  }}  { int jg = -1; VARIABLES;  if (is_face_y()) {
#line 93
{

#line 93 "isotropic.c"

    val(av.y,0,0,0) += 0.1*((val(u.y,0,0,0) + val(u.y,0,-1,0))/2. - ubar.y); }  }}  { int kg = -1; VARIABLES;  if (is_face_z()) {
#line 93
{

#line 93 "isotropic.c"

    val(av.z,0,0,0) += 0.1*((val(u.z,0,0,0) + val(u.z,0,0,-1))/2. - ubar.z); }  }}  end_foreach_face_generic()
#line 94
 end_foreach_face(); }
 end_trace("acceleration_0", "isotropic.c", 95); } return 0; } 







static int logfile_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (i++);   *ip = i; *tp = t;   return ret; } static int logfile_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t <= 300);   *ip = i; *tp = t;   return ret; } static int logfile (const int i, const double t, Event * _ev) { trace ("logfile", "isotropic.c", 103);  {
  coord ubar;
  {
#line 105
 {
    stats s = statsf(u.x);
    ubar.x = s.sum/s.volume;
  }
#line 105
 {
    stats s = statsf(u.y);
    ubar.y = s.sum/s.volume;
  }
#line 105
 {
    stats s = statsf(u.z);
    ubar.z = s.sum/s.volume;
  }}

  double ke = 0., vd = 0., vol = 0.;
   { 
#undef OMP_PARALLEL
#define OMP_PARALLEL()
OMP(omp parallel) {
double _ke = ke; double _vd = vd; double _vol = vol; 
#line 111

if (!is_constant(cm)) {
#undef val_cm
#define val_cm(a,i,j,k) val(a,i,j,k)
#undef fine_cm
#define fine_cm(a,i,j,k) fine(a,i,j,k)
#undef coarse_cm
#define coarse_cm(a,i,j,k) coarse(a,i,j,k)
#line 111
foreach(){

#line 111 "isotropic.c"
 {
    _vol += (cube(Delta)*val_cm(cm,0,0,0));
    {
#line 113
 {

      _ke += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(u.x,0,0,0) - ubar.x);

      _vd += (cube(Delta)*val_cm(cm,0,0,0))*(sq(val(u.x,1,0,0) - val(u.x,-1,0,0)) +
    sq(val(u.x,0,1,0) - val(u.x,0,-1,0)) +
    sq(val(u.x,0,0,1) - val(u.x,0,0,-1)))/sq(2.*Delta);
    }
#line 113
 {

      _ke += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(u.y,0,0,0) - ubar.y);

      _vd += (cube(Delta)*val_cm(cm,0,0,0))*(sq(val(u.y,0,1,0) - val(u.y,0,-1,0)) +
    sq(val(u.y,0,0,1) - val(u.y,0,0,-1)) +
    sq(val(u.y,1,0,0) - val(u.y,-1,0,0)))/sq(2.*Delta);
    }
#line 113
 {

      _ke += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(u.z,0,0,0) - ubar.z);

      _vd += (cube(Delta)*val_cm(cm,0,0,0))*(sq(val(u.z,0,0,1) - val(u.z,0,0,-1)) +
    sq(val(u.z,1,0,0) - val(u.z,-1,0,0)) +
    sq(val(u.z,0,1,0) - val(u.z,0,-1,0)))/sq(2.*Delta);
    }}
  } } end_foreach(); }
if (is_constant(cm)) {
const double _const_cm = _constant[cm.i -_NVARMAX];
NOT_UNUSED(_const_cm);
#undef val_cm
#define val_cm(a,i,j,k) _const_cm
#undef fine_cm
#define fine_cm(a,i,j,k) _const_cm
#undef coarse_cm
#define coarse_cm(a,i,j,k) _const_cm
#line 111
foreach(){

#line 111 "isotropic.c"
 {
    _vol += (cube(Delta)*val_cm(cm,0,0,0));
    {
#line 113
 {

      _ke += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(u.x,0,0,0) - ubar.x);

      _vd += (cube(Delta)*val_cm(cm,0,0,0))*(sq(val(u.x,1,0,0) - val(u.x,-1,0,0)) +
    sq(val(u.x,0,1,0) - val(u.x,0,-1,0)) +
    sq(val(u.x,0,0,1) - val(u.x,0,0,-1)))/sq(2.*Delta);
    }
#line 113
 {

      _ke += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(u.y,0,0,0) - ubar.y);

      _vd += (cube(Delta)*val_cm(cm,0,0,0))*(sq(val(u.y,0,1,0) - val(u.y,0,-1,0)) +
    sq(val(u.y,0,0,1) - val(u.y,0,0,-1)) +
    sq(val(u.y,1,0,0) - val(u.y,-1,0,0)))/sq(2.*Delta);
    }
#line 113
 {

      _ke += (cube(Delta)*val_cm(cm,0,0,0))*sq(val(u.z,0,0,0) - ubar.z);

      _vd += (cube(Delta)*val_cm(cm,0,0,0))*(sq(val(u.z,0,0,1) - val(u.z,0,0,-1)) +
    sq(val(u.z,1,0,0) - val(u.z,-1,0,0)) +
    sq(val(u.z,0,1,0) - val(u.z,0,-1,0)))/sq(2.*Delta);
    }}
  } } end_foreach(); }OMP(omp critical) ke += _ke;
mpi_all_reduce_double (ke, MPI_SUM);
OMP(omp critical) vd += _vd;
mpi_all_reduce_double (vd, MPI_SUM);
OMP(omp critical) vol += _vol;
mpi_all_reduce_double (vol, MPI_SUM);

#undef OMP_PARALLEL
#define OMP_PARALLEL() OMP(omp parallel)
}
#line 121
 }
  ke /= 2.*vol;
  vd *= 0.01/vol;

  if (i == 0)
    fprintf (ferr, "t dissipation energy Reynolds\n");
  fprintf (ferr, "%g %g %g %g\n",
    t, vd, ke, 2./3.*ke/0.01*sqrt(15.*0.01/vd));
 end_trace("logfile", "isotropic.c", 129); } return 0; } 
#line 138 "isotropic.c"
static int movie_expr0 (int * ip, double * tp, Event * _ev) {  int i = *ip; double t = *tp;  int ret = (t += 0.25);   *ip = i; *tp = t;   return ret; } static int movie_expr1 (int * ip, double * tp, Event * _ev) {   int i = *ip; double t = *tp;   int ret = ( t <= 150);   *ip = i; *tp = t;   return ret; } static int movie (const int i, const double t, Event * _ev) { trace ("movie", "isotropic.c", 138); 
{
  view ((struct _view_set){.fov = 44, .camera = "iso", .ty = .2,
 .width = 600, .height = 600, .bg = {1,1,1}, .samples = 4});
  clear();
  squares ((struct _squares){"u.y", .linear = true});
  squares ((struct _squares){"u.x", .linear = true, .n = {1,0,0}});
  scalar omega= new_scalar("omega");
  vorticity (u, omega);
  squares ((struct _squares){"omega", .linear = true, .n = {0,1,0}});
  scalar l2= new_scalar("l2");
  lambda2 (u, l2);
  isosurface ((struct _isosurface){"l2", -1});
  save ((struct _save){"movie.mp4"});
 delete (((scalar []){l2,omega,{-1}}));  end_trace("movie", "isotropic.c", 152); } return 0; } 
#line 91 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static double _boundary0 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 90 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann ((val_a_x(a.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_x(fm.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_x(alpha.x,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); } return 0.; } static double _boundary0_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 90 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 91
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 91
return  neumann_homogeneous(); } return 0.; }
#line 92 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static double _boundary1 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 91 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann (- (val_a_x(a.x,0,0,0)*val_fm_x(fm.x,0,0,0)/val_alpha_x(alpha.x,0,0,0))); } return 0.; } static double _boundary1_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 91 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 92
return  neumann_homogeneous(); } return 0.; }
#line 101 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static double _boundary2 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 100 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann ((val_a_y(a.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_y(fm.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_y(alpha.y,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); } return 0.; } static double _boundary2_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 100 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 101
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 101
return  neumann_homogeneous(); } return 0.; }
#line 102 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static double _boundary3 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 101 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann (- (val_a_y(a.y,0,0,0)*val_fm_y(fm.y,0,0,0)/val_alpha_y(alpha.y,0,0,0))); } return 0.; } static double _boundary3_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 101 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 102
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 102
return  neumann_homogeneous(); } return 0.; }
#line 105 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static double _boundary4 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 104 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann ((val_a_z(a.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))*val_fm_z(fm.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0))/val_alpha_z(alpha.z,(ig > 0 ? 1 : ig < 0 ? -1 : 0),(jg > 0 ? 1 : jg < 0 ? -1 : 0),(kg > 0 ? 1 : kg < 0 ? -1 : 0)))); } return 0.; } static double _boundary4_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 104 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 105
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 105
return  neumann_homogeneous(); } return 0.; }
#line 106 "/home/davidtang/basilisk/src/navier-stokes/centered.h"
static double _boundary5 (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 105 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann (- (val_a_z(a.z,0,0,0)*val_fm_z(fm.z,0,0,0)/val_alpha_z(alpha.z,0,0,0))); } return 0.; } static double _boundary5_homogeneous (Point point, Point neighbor, scalar _s, void * data) { int ig = neighbor.i - point.i;  if (ig == 0) ig = _attribute[_s.i].d.x;  NOT_UNUSED(ig); int jg = neighbor.j - point.j;  if (jg == 0) jg = _attribute[_s.i].d.y;  NOT_UNUSED(jg); int kg = neighbor.k - point.k;  if (kg == 0) kg = _attribute[_s.i].d.z;  NOT_UNUSED(kg); POINT_VARIABLES; 
#line 105 "/home/davidtang/basilisk/src/navier-stokes/centered.h"

if (!is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && !is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) coarse(a,i,j,k)
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) val(a,i,j,k)
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) coarse(a,i,j,k)
#line 106
return  neumann_homogeneous(); }
if (!is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann_homogeneous(); }
if (is_constant(a.x) && !is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
#undef val_fm_x
#define val_fm_x(a,i,j,k) val(a,i,j,k)
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_y
#define val_fm_y(a,i,j,k) val(a,i,j,k)
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) coarse(a,i,j,k)
#undef val_fm_z
#define val_fm_z(a,i,j,k) val(a,i,j,k)
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann_homogeneous(); }
if (!is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
#undef val_a_x
#define val_a_x(a,i,j,k) val(a,i,j,k)
#undef fine_a_x
#define fine_a_x(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) coarse(a,i,j,k)
#undef val_a_y
#define val_a_y(a,i,j,k) val(a,i,j,k)
#undef fine_a_y
#define fine_a_y(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) coarse(a,i,j,k)
#undef val_a_z
#define val_a_z(a,i,j,k) val(a,i,j,k)
#undef fine_a_z
#define fine_a_z(a,i,j,k) fine(a,i,j,k)
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) coarse(a,i,j,k)
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann_homogeneous(); }
if (is_constant(a.x) && is_constant(fm.x) && is_constant(alpha.x)) {
const struct { double x, y, z; } _const_a = {_constant[a.x.i -_NVARMAX], _constant[a.y.i - _NVARMAX], _constant[a.z.i - _NVARMAX]};
NOT_UNUSED(_const_a);
#undef val_a_x
#define val_a_x(a,i,j,k) _const_a.x
#undef fine_a_x
#define fine_a_x(a,i,j,k) _const_a.x
#undef coarse_a_x
#define coarse_a_x(a,i,j,k) _const_a.x
#undef val_a_y
#define val_a_y(a,i,j,k) _const_a.y
#undef fine_a_y
#define fine_a_y(a,i,j,k) _const_a.y
#undef coarse_a_y
#define coarse_a_y(a,i,j,k) _const_a.y
#undef val_a_z
#define val_a_z(a,i,j,k) _const_a.z
#undef fine_a_z
#define fine_a_z(a,i,j,k) _const_a.z
#undef coarse_a_z
#define coarse_a_z(a,i,j,k) _const_a.z
const struct { double x, y, z; } _const_fm = {_constant[fm.x.i -_NVARMAX], _constant[fm.y.i - _NVARMAX], _constant[fm.z.i - _NVARMAX]};
NOT_UNUSED(_const_fm);
#undef val_fm_x
#define val_fm_x(a,i,j,k) _const_fm.x
#undef fine_fm_x
#define fine_fm_x(a,i,j,k) _const_fm.x
#undef coarse_fm_x
#define coarse_fm_x(a,i,j,k) _const_fm.x
#undef val_fm_y
#define val_fm_y(a,i,j,k) _const_fm.y
#undef fine_fm_y
#define fine_fm_y(a,i,j,k) _const_fm.y
#undef coarse_fm_y
#define coarse_fm_y(a,i,j,k) _const_fm.y
#undef val_fm_z
#define val_fm_z(a,i,j,k) _const_fm.z
#undef fine_fm_z
#define fine_fm_z(a,i,j,k) _const_fm.z
#undef coarse_fm_z
#define coarse_fm_z(a,i,j,k) _const_fm.z
const struct { double x, y, z; } _const_alpha = {_constant[alpha.x.i -_NVARMAX], _constant[alpha.y.i - _NVARMAX], _constant[alpha.z.i - _NVARMAX]};
NOT_UNUSED(_const_alpha);
#undef val_alpha_x
#define val_alpha_x(a,i,j,k) _const_alpha.x
#undef fine_alpha_x
#define fine_alpha_x(a,i,j,k) _const_alpha.x
#undef coarse_alpha_x
#define coarse_alpha_x(a,i,j,k) _const_alpha.x
#undef val_alpha_y
#define val_alpha_y(a,i,j,k) _const_alpha.y
#undef fine_alpha_y
#define fine_alpha_y(a,i,j,k) _const_alpha.y
#undef coarse_alpha_y
#define coarse_alpha_y(a,i,j,k) _const_alpha.y
#undef val_alpha_z
#define val_alpha_z(a,i,j,k) _const_alpha.z
#undef fine_alpha_z
#define fine_alpha_z(a,i,j,k) _const_alpha.z
#undef coarse_alpha_z
#define coarse_alpha_z(a,i,j,k) _const_alpha.z
#line 106
return  neumann_homogeneous(); } return 0.; }
size_t datasize = 14*sizeof (double);
static int defaults (const int i, const double t, Event * _ev);
static int defaults_expr0 (int * ip, double * tp, Event * _ev);
static int cleanup (const int i, const double t, Event * _ev);
static int cleanup_expr0 (int * ip, double * tp, Event * _ev);
static int defaults_0 (const int i, const double t, Event * _ev);
static int defaults_0_expr0 (int * ip, double * tp, Event * _ev);
static int default_display (const int i, const double t, Event * _ev);
static int default_display_expr0 (int * ip, double * tp, Event * _ev);
static int init (const int i, const double t, Event * _ev);
static int init_expr0 (int * ip, double * tp, Event * _ev);
static int set_dtmax (const int i, const double t, Event * _ev);
static int set_dtmax_expr0 (int * ip, double * tp, Event * _ev);
static int stability (const int i, const double t, Event * _ev);
static int stability_expr0 (int * ip, double * tp, Event * _ev);
static int vof (const int i, const double t, Event * _ev);
static int vof_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_advection (const int i, const double t, Event * _ev);
static int tracer_advection_expr0 (int * ip, double * tp, Event * _ev);
static int tracer_diffusion (const int i, const double t, Event * _ev);
static int tracer_diffusion_expr0 (int * ip, double * tp, Event * _ev);
static int properties (const int i, const double t, Event * _ev);
static int properties_expr0 (int * ip, double * tp, Event * _ev);
static int advection_term (const int i, const double t, Event * _ev);
static int advection_term_expr0 (int * ip, double * tp, Event * _ev);
static int viscous_term (const int i, const double t, Event * _ev);
static int viscous_term_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration (const int i, const double t, Event * _ev);
static int acceleration_expr0 (int * ip, double * tp, Event * _ev);
static int projection (const int i, const double t, Event * _ev);
static int projection_expr0 (int * ip, double * tp, Event * _ev);
static int end_timestep (const int i, const double t, Event * _ev);
static int end_timestep_expr0 (int * ip, double * tp, Event * _ev);
static int perfs (const int i, const double t, Event * _ev);
static int perfs_expr0 (int * ip, double * tp, Event * _ev);
static int perf_plot (const int i, const double t, Event * _ev);
static int perf_plot_expr0 (int * ip, double * tp, Event * _ev);
static int runtime (const int i, const double t, Event * _ev);
static int runtime_expr0 (int * ip, double * tp, Event * _ev);
static int init_0 (const int i, const double t, Event * _ev);
static int init_0_expr0 (int * ip, double * tp, Event * _ev);
static int acceleration_0 (const int i, const double t, Event * _ev);
static int acceleration_0_expr0 (int * ip, double * tp, Event * _ev);
static int logfile (const int i, const double t, Event * _ev);
static int logfile_expr0 (int * ip, double * tp, Event * _ev);
static int logfile_expr1 (int * ip, double * tp, Event * _ev);
static int movie (const int i, const double t, Event * _ev);
static int movie_expr0 (int * ip, double * tp, Event * _ev);
static int movie_expr1 (int * ip, double * tp, Event * _ev);
static void _set_boundary0 (void);
static void _set_boundary1 (void);
static void _set_boundary2 (void);
static void _set_boundary3 (void);
static void _set_boundary4 (void);
static void _set_boundary5 (void);
void _init_solver (void) {
  void init_solver();
  init_solver();
  Events = (Event *) pmalloc (sizeof (Event), __func__, __FILE__, __LINE__);
  Events[0].last = 1;
  event_register ((Event){ 0, 1, defaults, {defaults_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/run.h", 42, "defaults"});
  event_register ((Event){ 0, 1, defaults_0, {defaults_0_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 126, "defaults"});
  event_register ((Event){ 0, 1, default_display, {default_display_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 179, "default_display"});
  event_register ((Event){ 0, 1, init, {init_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 188, "init"});
  event_register ((Event){ 0, 1, init_0, {init_0_expr0}, ((int *)0), ((double *)0),
    "isotropic.c", 72, "init"});
  event_register ((Event){ 0, 1, perfs, {perfs_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/perfs.h", 7, "perfs"});
  event_register ((Event){ 0, 1, perf_plot, {perf_plot_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/perfs.h", 24, "perf_plot"});
  event_register ((Event){ 0, 1, runtime, {runtime_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/maxruntime.h", 20, "runtime"});
  event_register ((Event){ 0, 2, logfile, {logfile_expr0, logfile_expr1}, ((int *)0), ((double *)0),
    "isotropic.c", 103, "logfile"});
  event_register ((Event){ 0, 2, movie, {movie_expr0, movie_expr1}, ((int *)0), ((double *)0),
    "isotropic.c", 138, "movie"});
  event_register ((Event){ 0, 1, cleanup, {cleanup_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/run.h", 50, "cleanup"});
  event_register ((Event){ 0, 1, set_dtmax, {set_dtmax_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 216, "set_dtmax"});
  event_register ((Event){ 0, 1, stability, {stability_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 218, "stability"});
  event_register ((Event){ 0, 1, vof, {vof_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 228, "vof"});
  event_register ((Event){ 0, 1, tracer_advection, {tracer_advection_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 229, "tracer_advection"});
  event_register ((Event){ 0, 1, tracer_diffusion, {tracer_diffusion_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 230, "tracer_diffusion"});
  event_register ((Event){ 0, 1, properties, {properties_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 237, "properties"});
  event_register ((Event){ 0, 1, advection_term, {advection_term_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 314, "advection_term"});
  event_register ((Event){ 0, 1, viscous_term, {viscous_term_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 344, "viscous_term"});
  event_register ((Event){ 0, 1, acceleration, {acceleration_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 380, "acceleration"});
  event_register ((Event){ 0, 1, acceleration_0, {acceleration_0_expr0}, ((int *)0), ((double *)0),
    "isotropic.c", 87, "acceleration"});
  event_register ((Event){ 0, 1, projection, {projection_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 423, "projection"});
  event_register ((Event){ 0, 1, end_timestep, {end_timestep_expr0}, ((int *)0), ((double *)0),
    "/home/davidtang/basilisk/src/navier-stokes/centered.h", 438, "end_timestep"});
  _attribute = (_Attributes *) pcalloc (datasize/sizeof(double), sizeof (_Attributes), __func__, __FILE__, __LINE__);
  all = (scalar *) pmalloc (sizeof (scalar)*15,__func__, __FILE__, __LINE__);
  for (int i = 0; i < 14; i++)
    all[i].i = i;
  all[14].i = -1;
  set_fpe();
  multigrid_methods();
  init_face_vector ((vector){{11},{12},{13}}, "av");
  init_face_vector ((vector){{8},{9},{10}}, "uf");
  init_scalar ((scalar){7}, "pf");
  init_vector ((vector){{4},{5},{6}}, "g");
  init_vector ((vector){{1},{2},{3}}, "u");
  init_scalar ((scalar){0}, "p");
  init_const_scalar ((scalar){_NVARMAX+7}, "zeroc",  0.);
  init_const_scalar ((scalar){_NVARMAX+6}, "unity",  1.);
  init_const_vector ((vector){{_NVARMAX+3},{_NVARMAX+4},{_NVARMAX+5}}, "unityf", (double []) {1.,1.,1.});
  init_const_vector ((vector){{_NVARMAX+0},{_NVARMAX+1},{_NVARMAX+2}}, "zerof", (double []) {0.,0.,0.});
  _set_boundary0();
  _set_boundary1();
  _set_boundary2();
  _set_boundary3();
  _set_boundary4();
  _set_boundary5();
}
