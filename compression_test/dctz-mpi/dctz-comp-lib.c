/**
 * @file dctz-comp-lib.c
 * @author Seung Woo Son
 * @date July 2019
 * @brief DCTZ compression library routine
 * (C) 2019 University of Massachuetts Lowell.
       See LICENSE in top-level directory.
*/

#include <stdlib.h>
#include <memory.h>
#include <string.h>
#ifdef TIME_DEBUG
#include <sys/time.h>
#endif /* TIME_DEBUG */
#include <unistd.h>
#ifdef _OPENMP
#include <omp.h>
#endif
#include <pthread.h>
#include "zlib.h"
#include "dctz.h"
#include "dct.h"

#define DEF_MEM_LEVEL 8

union
{
  float *f;
  double *d;
} a, a_x;

void *compress_thread(void *arg)
{
  z_stream *defstream = (z_stream *)arg;
#ifdef DEBUG
  printf("compress started ...\n");
#endif
  deflate(defstream, Z_FINISH);
#ifdef DEBUG
  printf("done! compression...\n");
#endif
  uLong ret = defstream->total_out;
  deflateEnd(defstream);
  pthread_exit((void *)ret);
}

int dctz_compress(t_var *var, int N, size_t *outSize, t_var *var_z, double error_bound)
{
  int i, j, nblk, rem;
#ifdef TIME_DEBUG
  struct timeval start_t, end_t, gstart_t;
  double sf_t, dct_t, DC_AC_t, zlib_t, comp_t, malloc_t, genbin_t;
#endif
  double *bin_maxes, *bin_center, bin_width, range_min, range_max;
  t_bin_id *bin_index, *bin_indexz, *bin_indexz2;
#ifdef USE_TRUNCATE
  float *DC, *DCz, *DCz2, *AC_exact, *AC_exactz, *AC_exactz2;
#else
  double *DC, *DCz, *DCz2, *AC_exact, *AC_exactz, *AC_exactz2;
#endif
  struct header h;
  t_bstat bs;
  size_t type_size = 0;
#ifdef USE_QTABLE
  double *qtable; /* Quantizer Table */
#endif

  if (var->datatype == DOUBLE)
    type_size = sizeof(double);
  else /* FLOAT */
    type_size = sizeof(float);

  if (var->datatype == DOUBLE)
  {
    if (NULL == (a_x.d = (double *)malloc(N * type_size)))
    {
      fprintf(stderr, "Out of memory: a_x\n");
      exit(1);
    }
  }
  else
  { /* FLOAT */
    if (NULL == (a_x.f = (float *)malloc(N * type_size)))
    {
      fprintf(stderr, "Out of memory: a_x\n");
      exit(1);
    }
  }

  if (error_bound < 1E-6)
  {
    printf("ERROR BOUND is not acceptable");
    exit(1);
  }

  /* TODO: bin_maxes and bin_center should be double or float? */
  if (NULL == (bin_maxes = (double *)malloc(NBINS * sizeof(double))))
  {
    fprintf(stderr, "Out of memory: bin_maxes\n");
    exit(1);
  }

  if (NULL == (bin_center = (double *)malloc(NBINS * sizeof(double))))
  {
    fprintf(stderr, "Out of memory: bin_center\n");
    exit(1);
  }

#ifdef DEBUG
  for (i = 0; i < BLK_SZ; i++)
  { /* show the first block */
    printf("a[%d] = %e\n", i, var->buf.d[i]);
    if (i % BLK_SZ == 0 && i != 0)
      printf("\n");
  }
#endif

#ifdef USE_QTABLE
  /* Start of Initialize  Quantizer Table */
  if (NULL == (qtable = (double *)malloc(BLK_SZ * sizeof(double))))
  {
    fprintf(stderr, "Out of memory: qtable\n");
    exit(1);
  }

  for (i = 0; i < BLK_SZ; i++)
  {
    qtable[i] = 0.0;
  }
  if (NULL == (bin_index = (t_bin_id *)malloc(2 * N * sizeof(t_bin_id))))
  {
    fprintf(stderr, "Out of memory: bin_index[]\n");
    exit(1);
  }
  memset(bin_index, 0, sizeof(t_bin_id) * 2 * N);
#ifdef DEBUG
  for (i = 0; i < BLK_SZ; i++)
  {
    printf("qtable[%d] = %e\n", i, qtable[i]);
  }
#endif
  /* End of Initialize  Quantizer Table */
#else
  if (NULL == (bin_index = (t_bin_id *)malloc(N * sizeof(t_bin_id))))
  {
    fprintf(stderr, "Out of memory: bin_index[]\n");
    exit(1);
  }
  memset(bin_index, 0, sizeof(t_bin_id) * N);

#endif /* USE_QTABLE */

#ifdef TIME_DEBUG
  gettimeofday(&start_t, NULL);
  gstart_t = start_t;
#endif

  /* determine scaling factor */
  calc_data_stat(var, &bs, N);

  if (var->datatype == DOUBLE)
  {
#ifdef DEBUG
    printf("scaling factor = %f\n", bs.sf.d);
#endif
    double xscale = pow(10, (bs.sf.d) - 1);
    /* apply scaling factor */
    if (bs.sf.d != 1.0)
    {
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(bs.sf.d)
#endif
      for (i = 0; i < N; i++)
        var->buf.d[i] /= xscale;
    }
  }
  else
  { /* FLOAT */
#ifdef DEBUG
    printf("scaling factor = %f\n", bs.sf.f);
#endif
    float xscale = pow(10, (bs.sf.f) - 1);
    /* apply scaling factor */
    if (bs.sf.f != 1.0)
    {
#ifdef _OPENMP
#pragma omp parallel for private(i) shared(bs.sf.f)
#endif
      for (i = 0; i < N; i++)
        var->buf.f[i] /= xscale;
    }
  }

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  sf_t = (double)((end_t.tv_sec * 1000000 + end_t.tv_usec) - (start_t.tv_sec * 1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

  /* DCT over decomposed blocks */
  nblk = CEIL(N, BLK_SZ);
  rem = N % BLK_SZ;
#ifdef DEBUG
  printf("\nnumber of blocks = %d, remainder = %d\n", nblk, rem);
#endif

#ifdef USE_TRUNCATE
  if (NULL == (DC = (float *)malloc(nblk * sizeof(float))))
  {
    fprintf(stderr, "Out of memory: DC[]\n");
    exit(1);
  }
#else
  if (NULL == (DC = (double *)malloc(nblk * sizeof(double))))
  {
    fprintf(stderr, "Out of memory: DC[]\n");
    exit(1);
  }
#endif

#ifdef USE_TRUNCATE
  if (NULL == (DCz = (float *)malloc(nblk * sizeof(float))))
  {
    fprintf(stderr, "Out of memory: DCz[]\n");
    exit(1);
  }
  memset(DCz, 0, sizeof(float) * nblk); /* TODO: is it necessary? */
#else
  if (NULL == (DCz = (double *)malloc(nblk * sizeof(double))))
  {
    fprintf(stderr, "Out of memory: DCz[]\n");
    exit(1);
  }
#endif

  if (NULL == (bin_indexz = (t_bin_id *)malloc(N * sizeof(t_bin_id))))
  {
    fprintf(stderr, "Out of memory: bin_indexz[]\n");
    exit(1);
  }
  memset(bin_indexz, 0, sizeof(t_bin_id) * N);

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  malloc_t = (double)((end_t.tv_sec * 1000000 + end_t.tv_usec) - (start_t.tv_sec * 1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

  if (var->datatype == DOUBLE)
    gen_bins(bs.min.d, bs.max.d, bin_maxes, bin_center, NBINS, error_bound);
  else /* FLOAT */
    gen_bins(bs.min.f, bs.max.f, bin_maxes, bin_center, NBINS, error_bound);

  int half = NBINS / 2;
  bin_width = error_bound * 2 * BRSF;
  range_min = -(half * 2 + 1) * (error_bound * BRSF);
  range_max = (half * 2 + 1) * (error_bound * BRSF);

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  genbin_t = (double)((end_t.tv_sec * 1000000 + end_t.tv_usec) - (start_t.tv_sec * 1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

#ifdef USE_TRUNCATE
  if (NULL == (AC_exact = (float *)malloc(N * sizeof(float))))
  {
    fprintf(stderr, "Out of memory: AC_exact\n");
    exit(1);
  }
  memset(AC_exact, 0, sizeof(float) * N); /* TODO: is it necessary? */
#else
  if (NULL == (AC_exact = (double *)malloc(N * sizeof(double))))
  {
    fprintf(stderr, "Out of memory: AC_exact\n");
    exit(1);
  }
  memset(AC_exact, 0, sizeof(double) * N); /* TODO: is it necessary? */
#endif

#ifdef USE_TRUNCATE
  if (NULL == (AC_exactz = (float *)malloc(N * sizeof(float))))
  {
    fprintf(stderr, "Out of memory: AC_exactz[]\n");
    exit(1);
  }
  memset(AC_exactz, 0, sizeof(float) * N); /* TODO: is it necessary? */
#else
  if (NULL == (AC_exactz = (double *)malloc(N * sizeof(double))))
  {
    fprintf(stderr, "Out of memory: AC_exactz[]\n");
    exit(1);
  }
  memset(AC_exactz, 0, sizeof(double) * N); /* TODO: is it necessary? */
#endif

  if (var->datatype == DOUBLE)
    dct_init(BLK_SZ);
  else /* FLOAT */
    dct_init_f(BLK_SZ);

  int tot_AC_exact_count = 0;
  /* DCT block decomposed */
  for (i = 0; i < nblk; i++)
  { /* for each decomposed blk */
    int l_blk_sz = ((i == nblk - 1) && (rem != 0)) ? rem : BLK_SZ;
    if ((i == nblk - 1) && (rem != 0))
    {
      if (var->datatype == DOUBLE)
      {
        dct_finish();
        dct_init(rem);
      }
      else
      { /* FLOAT */
        dct_finish_f();
        dct_init_f(rem);
      }
    }
    if (var->datatype == DOUBLE)
      dct_fftw(var->buf.d + i * BLK_SZ, a_x.d + i * BLK_SZ, l_blk_sz, nblk);
    else /* FLOAT */
      dct_fftw_f(var->buf.f + i * BLK_SZ, a_x.f + i * BLK_SZ, l_blk_sz, nblk);

#ifdef DEBUG
    printf("block %d: after DCT:\n", i);
    for (j = 0; j < BLK_SZ && (i < 3); j++)
    { /* show the first block only */
      printf("a_x[%d] = %e \n", i * BLK_SZ + j, a_x.d[i * BLK_SZ + j]);
    }
    printf("\n");
#endif

#ifdef USE_TRUNCATE
    DC[i] = (float)(var->datatype == DOUBLE ? a_x.d[i * BLK_SZ] : a_x.f[i * BLK_SZ]); /* save DC component in truncated*/
#else
    DC[i] = (double)a_x.d[i * BLK_SZ];      /* save DC component */
#endif
    bin_index[i * BLK_SZ] = NBINS; /* store as it is */

    for (j = 1; j < l_blk_sz; j++)
    {
      if (var->datatype == DOUBLE)
      {
        double item = a_x.d[i * BLK_SZ + j];
        t_bin_id bin_id;
        if (item < range_min || item > range_max)
        {
          bin_id = NBINS;
#ifdef USE_QTABLE
          /* The Start  of Making Quantizer Table -QT applied to block coefficients  */
          if (fabs(item) >= qtable[j])
            qtable[j] = fabs(item);
#endif /* USE_QTABLE */
        }
        else
          bin_id = (t_bin_id)((item - range_min) / bin_width);
#ifdef DEBUG
        printf("bin_id = %d\n", bin_id);
#endif
        bin_index[i * BLK_SZ + j] = bin_id;
#ifdef DEBUG
        printf("a_x[%d]=%e => %d\n", i * BLK_SZ + j, (double)item, bin_id);
#endif
      } /* DOUBLE */
      else
      { /* FLOAT */
        float item = a_x.f[i * BLK_SZ + j];
        t_bin_id bin_id;
        if (item < range_min || item > range_max)
        {
          bin_id = NBINS;
#ifdef USE_QTABLE
          /* The Start  of Making Quantizer Table -QT applied to block coefficients  */
          if (fabs(item) >= qtable[j])
            qtable[j] = fabs(item);
#endif /* USE_QTABLE */
        }
        else
          bin_id = (t_bin_id)((item - range_min) / bin_width);
#ifdef DEBUG
        printf("bin_id = %d\n", bin_id);
#endif
        bin_index[i * BLK_SZ + j] = bin_id;
#ifdef DEBUG
        printf("a_x[%d]=%e => %d\n", i * BLK_SZ + j, (float)item, bin_id);
#endif
      } /* else: FLOAT */
    }
    /* The End of of Making Quantizer Table  */
  }
  if (var->datatype == DOUBLE)
    dct_finish();
  else /* FLOAT */
    dct_finish_f();

#ifdef DCT_FILE_DEBUG
  FILE *fp = fopen("dct_result.bin", "w+");
  if (var->datatype == DOUBLE)
    fwrite(a_x.d, sizeof(double), N, fp);
  else /* FLOAT */
    fwrite(a_x.f, sizeof(float), N, fp);
  fclose(fp);
#endif

#ifdef USE_QTABLE
#ifdef DEBUG
  printf("Quantizer Table:\n");
  for (j = 0; j < BLK_SZ; j++)
  { /* Show Quantizer Table */
    printf("before qtable[%d] = %e \n", j, qtable[j]);
  }
#endif

  for (j = 1; j < BLK_SZ; j++)
  { /* Show Quantizer Table */
    // if (qtable[j] < bin_maxes[NBINS-1]) {
    if (qtable[j] < 1.0)
    {
      qtable[j] = 1.0;
    }
  }

#ifdef DEBUG
  printf("Quantizer Table:\n");
  for (j = 0; j < BLK_SZ; j++)
  { /* Show Quantizer Table */
    printf("after qtable[%d] = %e \n", j, qtable[j]);
  }
#endif
#endif

  unsigned int k = N;
  double qt_factor = (NBINS == 255 ? 10.0 : 2000.0);

  for (i = 0; i < nblk; i++)
  {
    int l_blk_sz = ((i == nblk - 1) && (rem != 0)) ? rem : BLK_SZ;
    for (j = 1; j < l_blk_sz; j++)
    {
      t_bin_id bin_id;
      bin_id = bin_index[i * BLK_SZ + j];
      if (bin_id == NBINS)
      {
#ifdef USE_QTABLE
        double item = (var->datatype == DOUBLE) ? a_x.d[i * BLK_SZ + j] : a_x.f[i * BLK_SZ + j]; /* in case of FLOAT, it will be cast to double */
                                                                                                 /* if out of bin area, normalize it to the area from range_max/range_min to range_max/range_min +/- error_bound */
        if (item < range_min)
        {
          item = (item / qtable[j]) * error_bound * qt_factor + range_min;
        }
        else if (item > range_max)
        {
          item = (item / qtable[j]) * error_bound * qt_factor + range_max;
        }
        if (var->datatype == DOUBLE)
          a_x.d[i * BLK_SZ + j] = item; /* update a_x with updated value */
        else                            /* FLOAT */
          a_x.f[i * BLK_SZ + j] = item; /* update a_x with updated value */
        if (item < range_min || item > range_max)
        {
          bin_id = NBINS;
#ifdef USE_TRUNCATE
          AC_exact[tot_AC_exact_count++] = (float)(var->datatype == DOUBLE ? a_x.d[i * BLK_SZ + j] : a_x.f[i * BLK_SZ + j]);
#else
    AC_exact[tot_AC_exact_count++] = (double)(var->datatype == DOUBLE ? a_x.d[i*BLK_SZ+j] | a_x.f[i*BLK_SZ+j]);
    ;
#endif /* USE_TRUNCATE */
        }
        else
          bin_id = (t_bin_id)((item - range_min) / bin_width);
        bin_index[k++] = bin_id;
#ifdef DEBUG
        printf("a_x[%d]=%e => %d\n", i * BLK_SZ + j, (double)item, bin_id);
#endif /* DEBUG */
#else
#ifdef USE_TRUNCATE
        AC_exact[tot_AC_exact_count++] = (float)(var->datatype == DOUBLE ? a_x.d[i * BLK_SZ + j] : a_x.f[i * BLK_SZ + j]);
#else
  AC_exact[tot_AC_exact_count++] = (double)(var->datatype == DOUBLE ? a_x.d[i*BLK_SZ+j] | a_x.f[i*BLK_SZ+j]); /* in case of FLOAT, casting to double makes sense? */
#endif
#endif /* USE_QTABLE */
      }
    }
  }

#ifdef DEBUG
  printf("total AC_exact_count = %d\n", tot_AC_exact_count);
#endif

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  dct_t = (double)((end_t.tv_sec * 1000000 + end_t.tv_usec) - (start_t.tv_sec * 1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

  free(bin_maxes);

#ifdef DEBUG
  int bin_freq[NBINS + 1] = {0};

  i = 0;
  while (i < N)
  {
    bin_freq[(int)bin_index[i++]]++;
  }
  printf("i=%d\n", i);

  int sum = 0;
  printf("bin_freq: ");
  for (i = 0; i < NBINS + 1; i++)
  {
    printf("%d, ", bin_freq[i]);
    sum += bin_freq[i];
  }
  printf("sum=%d\n", sum);
#endif

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);

  DC_AC_t = (double)((end_t.tv_sec * 1000000 + end_t.tv_usec) - (start_t.tv_sec * 1000000 + start_t.tv_usec));

  gettimeofday(&start_t, NULL);
#endif

  // char bin_index_file[640];
  // FILE *fp_index;
  // sprintf(bin_index_file, "bin_index.bin");
  // fp_index = fopen(bin_index_file, "wb");
  // fwrite(bin_index, N, 1, fp_index);
  // fclose(fp_index);

#ifdef DEBUG
  printf("tot_AC_exact_count=%d\n", tot_AC_exact_count);
#ifdef USE_QTABLE
  printf("bin_index before compression = %lu\n", k * sizeof(t_bin_id));
#else
  printf("bin_index before compression = %lu\n", N * sizeof(t_bin_id));
#endif
#ifdef USE_TRUNCATE
  printf("DC before compression = %lu\n", nblk * sizeof(float));
  printf("AC_exact before compression = %lu\n", tot_AC_exact_count * sizeof(float));
#else
  printf("DC before compression = %lu\n", nblk * sizeof(double));
  printf("AC_exact before compression = %lu\n", tot_AC_exact_count * sizeof(double));
#endif
#endif

  pthread_t thread[3];
  pthread_attr_t attr; /* thread attributes (left at defaults) */

  /* set defaults (not all pthread implementations default to joinable) */
  pthread_attr_init(&attr);
  pthread_attr_setdetachstate(&attr, PTHREAD_CREATE_JOINABLE);

  /* setup for compress */
  z_stream defstream[3];

  defstream[0].zalloc = Z_NULL;
  defstream[0].zfree = Z_NULL;
  defstream[0].opaque = Z_NULL;

  /* compress bin_index */
#ifdef USE_QTABLE
  uLong ucompSize_binindex = k * sizeof(t_bin_id);
#else
  uLong ucompSize_binindex = N * sizeof(t_bin_id);
#endif
  uLong compSize_binindex = compressBound(ucompSize_binindex);

  int windowBits = 14;
  deflateInit2(&defstream[0], 1, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);

  defstream[0].avail_in = ucompSize_binindex;
  defstream[0].next_in = (Bytef *)bin_index;
  defstream[0].avail_out = compSize_binindex;
  defstream[0].next_out = (Bytef *)bin_indexz;
  defstream[0].data_type = Z_UNKNOWN; /* Z_ASCII, Z_BINARY, Z_UNKNOWN */

  if (pthread_create(&thread[0], &attr, compress_thread, (void *)&defstream[0]))
  {
    fprintf(stderr, "Error creating thread\n");
    exit(0);
  }

  /* compress DC */
  defstream[1].zalloc = Z_NULL;
  defstream[1].zfree = Z_NULL;
  defstream[1].opaque = Z_NULL;

#ifdef USE_TRUNCATE
  uLong ucompSize_DC = nblk * sizeof(float);
  uLong compSize_DC = compressBound(ucompSize_DC);
#else
  uLong ucompSize_DC = nblk * sizeof(double);
  uLong compSize_DC = compressBound(ucompSize_DC);
#endif

  deflateInit2(&defstream[1], 1, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);

  defstream[1].avail_in = ucompSize_DC;
  defstream[1].next_in = (Bytef *)DC;
  defstream[1].avail_out = compSize_DC;
  defstream[1].next_out = (Bytef *)DCz;
  defstream[1].data_type = Z_UNKNOWN;

  if (pthread_create(&thread[1], &attr, compress_thread, (void *)&defstream[1]))
  {
    fprintf(stderr, "Error creating thread\n");
    exit(0);
  }

  /* compress AC_exact */
  defstream[2].zalloc = Z_NULL;
  defstream[2].zfree = Z_NULL;
  defstream[2].opaque = Z_NULL;

#ifdef USE_TRUNCATE
  uLong ucompSize_AC_exact = N * sizeof(float);
  uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#else
  uLong ucompSize_AC_exact = N * sizeof(double);
  uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#endif

  deflateInit2(&defstream[2], 1, Z_DEFLATED, windowBits, DEF_MEM_LEVEL, Z_DEFAULT_STRATEGY);

  defstream[2].avail_in = ucompSize_AC_exact;
  defstream[2].next_in = (Bytef *)AC_exact;
  defstream[2].avail_out = compSize_AC_exact;
  defstream[2].next_out = (Bytef *)AC_exactz;
  defstream[2].data_type = Z_UNKNOWN;

  if (pthread_create(&thread[2], &attr, compress_thread, (void *)&defstream[2]))
  {
    fprintf(stderr, "Error creating thread\n");
    exit(0);
  }

#ifdef USE_TRUNCATE
  // uLong ucompSize_AC_exact = tot_AC_exact_count*sizeof(float);
  //   uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#else
  // uLong ucompSize_AC_exact = tot_AC_exact_count*sizeof(double);
  //  uLong compSize_AC_exact = compressBound(ucompSize_AC_exact);
#endif
  void *ret;
  for (i = 0; i < 3; i++)
  {
    pthread_join(thread[i], &ret);
#ifdef DEBUG
    printf("thread %d joined\n", i);
#endif
    switch (i)
    {
    case 0:
      compSize_binindex = (uLong)ret;
      break;
    case 1:
      compSize_DC = (uLong)ret;
      break;
    case 2:
      compSize_AC_exact = (uLong)ret;
      break;
    }
  }

  pthread_attr_destroy(&attr);
#if 0
  compSize_binindex = defstream[0].total_out; /* update with actual size */
  deflateEnd(&defstream[0]);

  compSize_DC = defstream[1].total_out; /* update with actual size */
  deflateEnd(&defstream[1]);

  compSize_AC_exact_count = defstream[2].total_out; /* update with actual size */
  deflateEnd(&defstream[2]);
#endif

  bin_indexz2 = (t_bin_id *)realloc(bin_indexz, compSize_binindex); /* TODO: check error */

#ifdef SIZE_DEBUG
  printf("Compressed bin_index size is: %lu\n", compSize_binindex);
#endif

  DCz2 = realloc(DCz, compSize_DC); /* TODO: check error */
#ifdef SIZE_DEBUG
  printf("Compressed DC size is: %lu\n", compSize_DC);
#endif

  AC_exactz2 = realloc(AC_exactz, compSize_AC_exact); /* TODO: check error */
#ifdef SIZE_DEBUG
  printf("Compressed AC_exact size is: %lu\n", compSize_AC_exact);
#endif

#ifdef TIME_DEBUG
  gettimeofday(&end_t, NULL);
  double comp_rate;

  zlib_t = (double)((end_t.tv_sec * 1000000 + end_t.tv_usec) - (start_t.tv_sec * 1000000 + start_t.tv_usec));
  comp_t = (double)((end_t.tv_sec * 1000000 + end_t.tv_usec) - (gstart_t.tv_sec * 1000000 + gstart_t.tv_usec));
  comp_rate = (N * sizeof(double) / (double)(1024 * 1024)) / (comp_t / 1000000);

  printf("sf_t=%f(s), dct_t=%f(s), zlib_t(compress)=%f(s)\n", sf_t / 1000000, dct_t / 1000000, zlib_t / 1000000);
  printf("malloc_t=%f(s), genbin=%f(s), DC_AC_t=%f(s)\n", malloc_t / 1000000, genbin_t / 1000000, DC_AC_t / 1000000);
  printf("comp_time = %f (s), compression rate = %f (MB/s)\n", comp_t / 1000000, comp_rate);
#endif

  *outSize = sizeof(struct header) + compSize_binindex + compSize_DC + compSize_AC_exact;

  h.datatype = var->datatype;
  h.num_elements = N;
  h.error_bound = error_bound;
  h.tot_AC_exact_count = tot_AC_exact_count;
  if (var->datatype == DOUBLE)
    h.scaling_factor.d = bs.sf.d;
  else /* FLOAT */
    h.scaling_factor.f = bs.sf.f;
  h.bindex_sz_compressed = compSize_binindex;
  h.DC_sz_compressed = compSize_DC;
  h.AC_exact_sz_compressed = compSize_AC_exact;
#ifdef USE_QTABLE
  h.bindex_count = k;
#endif
  // h.AC_exact_count_sz_compressed = compSize_AC_exact_count;

  unsigned char *cur_p;
  if (var->datatype == DOUBLE)
    cur_p = (unsigned char *)(var_z->buf.d);
  else
    cur_p = (unsigned char *)(var_z->buf.f);
  memcpy(cur_p, &h, sizeof(struct header));
  cur_p += sizeof(struct header);
  memcpy(cur_p, bin_indexz2, compSize_binindex);
  cur_p += compSize_binindex;
  // memcpy (cur_p, AC_exact_countz2, compSize_AC_exact_count);
  // cur_p += compSize_AC_exact_count;
  memcpy(cur_p, DCz2, compSize_DC);
  cur_p += compSize_DC;
  memcpy(cur_p, AC_exactz2, compSize_AC_exact);
#ifdef USE_QTABLE
  cur_p += compSize_AC_exact;
  memcpy(cur_p, qtable, BLK_SZ * sizeof(double));
#endif /* USE_QTABLE */

  if (var->datatype == DOUBLE)
    free(a_x.d);
  else /* FLOAT */
    free(a_x.f);
  // free(a_x);
  free(DC);
  // free(DCz);
  free(DCz2);
  free(bin_center);
  // free(AC_exact_count);
  // free(AC_exact_countz2);
  free(AC_exact);
  // free(AC_exactz);
  free(AC_exactz2);
  free(bin_index);
  // free(bin_indexz);
  free(bin_indexz2);
#ifdef USE_QTABLE
  free(qtable);
#endif

#ifndef SIZE_DEBUG
  printf("outSize = %zu\n", *outSize);
#endif

  return (1);
}
