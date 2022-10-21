#include <cstdio>
#include <cstdlib>
#include <cmath>

#include "netcdf.h"
#include "cJSON.h"
#include <mpi.h>

#include <iostream>
#include <string>
#include <iomanip>
#include <stdlib.h>
#include <typeinfo>
#include "tthresh.hpp"
#include "compress_revise.hpp"
#include "decompress_revise.hpp"
#include "Slice.hpp"

using namespace std;

int main(int argc, char *argv[])
{
  // bool binaryOutput = true;
  // int printCmpResults = 0;
  // bool compression = true;
  // bool decompression = false;
  // // int dataType = SZ_FLOAT;
  // char *inPath = nullptr;
  // char *cmpPath = nullptr;
  // char *conPath = "sz.config";
  // char *decPath = nullptr;
  // bool delCmpPath = false;

  // char *errBoundMode = nullptr;
  // char *errBound = nullptr;
  // char *absErrorBound = nullptr;
  // char *relErrorBound = nullptr;
  // char *pwrErrorBound = nullptr;
  // char *psnrErrorBound = nullptr;
  // char *normErrorBound = nullptr;

  // bool sz2mode = false;

  size_t r4 = 0;
  size_t r3 = 0;
  size_t r2 = 0;
  size_t r1 = 0;

  // size_t i = 0;
  int status;

  int nprocx, nprocy;

  int myid, mpi_size;

  char init_root[1024];
  char out_root[1024];

  char fnm_init[1024];

  double tolerance;

  MPI_Init(&argc, &argv);
  MPI_Comm_rank(MPI_COMM_WORLD, &myid);

  char wave_json[1024];
  int velocity, stress, strain;
  int reconstructed;
  sprintf(wave_json, "test.json");

  FILE *fp;
  fp = fopen(wave_json, "r");
  fseek(fp, 0, SEEK_END);
  long len = ftell(fp);
  fseek(fp, 0, SEEK_SET);
  char *wave_tomo = (char *)malloc(len + 1);
  fread(wave_tomo, 1, len, fp);
  fclose(fp);

  cJSON *root_wave = cJSON_Parse(wave_tomo);
  if (NULL == root_wave)
  {
    printf("Error at parsing json");
    exit(-1);
  }
  cJSON *item;
  if (item = cJSON_GetObjectItem(root_wave, "number_of_mpiprocs_x"))
  {
    nprocx = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root_wave, "number_of_mpiprocs_y"))
  {
    nprocy = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root_wave, "sz_timestep"))
  {

    r4 = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root_wave, "input_root"))
  {
    sprintf(init_root, "%s", item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root_wave, "output_root"))
  {
    sprintf(out_root, "%s", item->valuestring);
  }
  if (item = cJSON_GetObjectItem(root_wave, "velocity"))
  {
    velocity = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root_wave, "stress"))
  {
    stress = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root_wave, "strain"))
  {
    strain = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root_wave, "reconstructed"))
  {
    reconstructed = item->valueint;
  }
  if (item = cJSON_GetObjectItem(root_wave, "tolerance"))
  {
    tolerance = item->valuedouble;
  }
  //////////////////////////////////////////////////////////////////////////////
  int pdims[2] = {nprocx, nprocy};
  int periods[2] = {0, 0};
  int myid2[2];

  MPI_Comm comm;

  MPI_Cart_create(MPI_COMM_WORLD, 2, pdims, periods, 0, &comm);

  MPI_Cart_coords(comm, myid, 2, myid2);

  int proc_x = myid2[0];
  int proc_y = myid2[1];

  printf("myid = %d, i = %d, j = %d\n", myid, proc_x, proc_y);
  //////////////////////////////////////////////////////////////////////////////
  int ncid;
  int Vxid, Vyid, Vzid, Txxid, Tyyid, Tzzid, Txyid, Txzid, Tyzid;
  int Exxid, Eyyid, Ezzid, Exyid, Exzid, Eyzid;
  int xi_id, eta_id, zeta_id, time_id;
  size_t xi_length, eta_length, zeta_length, time_length, sum_length, compress_length;
  int compress_it_num;

  int ncid_r;
  int dim_id_r[4];
  int varid_r[3];

  // double compression_time, write_time, reconstruct_time;

  double compress_time[1] = {0.0};
  double decompress_time[1] = {0.0};

  if (velocity == 1)
  {
    sprintf(fnm_init, "%s/volume_velocity_px%d_py%d.nc", init_root, proc_x, proc_y);
    printf("%s\n", fnm_init);
    nc_open(fnm_init, NC_NOWRITE, &ncid);

    nc_inq_varid(ncid, "Vx", &Vxid);
    nc_inq_varid(ncid, "Vy", &Vyid);
    nc_inq_varid(ncid, "Vz", &Vzid);

    nc_inq_dimid(ncid, "xi", &xi_id);
    nc_inq_dimid(ncid, "eta", &eta_id);
    nc_inq_dimid(ncid, "zeta", &zeta_id);
    nc_inq_dimid(ncid, "time", &time_id);

    nc_inq_dimlen(ncid, xi_id, &xi_length);
    nc_inq_dimlen(ncid, eta_id, &eta_length);
    nc_inq_dimlen(ncid, zeta_id, &zeta_length);
    nc_inq_dimlen(ncid, time_id, &time_length);

    char rec_file[1024];

    r3 = zeta_length;
    r2 = eta_length;
    r1 = xi_length;

    compress_length = r1 * r2 * r3 * r4;

    if (reconstructed == 1)
    {
      sprintf(rec_file, "reconstructed/reconstructed_velocity_px%d_py%d.nc", proc_x, proc_y);

      nc_create(rec_file, NC_CLOBBER | NC_64BIT_OFFSET, &ncid_r);
      nc_def_dim(ncid_r, "time", NC_UNLIMITED, &dim_id_r[0]);
      nc_def_dim(ncid_r, "zeta", zeta_length, &dim_id_r[1]);
      nc_def_dim(ncid_r, "eta", eta_length, &dim_id_r[2]);
      nc_def_dim(ncid_r, "xi", xi_length, &dim_id_r[3]);

      nc_def_var(ncid_r, "Vx", NC_FLOAT, 4, dim_id_r, &varid_r[0]);
      nc_def_var(ncid_r, "Vy", NC_FLOAT, 4, dim_id_r, &varid_r[1]);
      nc_def_var(ncid_r, "Vz", NC_FLOAT, 4, dim_id_r, &varid_r[2]);
      nc_enddef(ncid_r);
    }

    // size_t outSize_Vx, outSize_Vy, outSize_Vz;
    char outputFilePath[1024];

    compress_it_num = time_length / r4;

    for (int ii = 0; ii < compress_it_num; ii++)
    {
      size_t startp[] = {ii * r4, 0, 0, 0};
      size_t countp[] = {r4, r3, r2, r1};

      MPI_Barrier(MPI_COMM_WORLD);

      printf("startp = %d * %d\n", ii, r4);

      float *compress_Vx = (float *)malloc(sizeof(float) * compress_length);
      float *compress_Vy = (float *)malloc(sizeof(float) * compress_length);
      float *compress_Vz = (float *)malloc(sizeof(float) * compress_length);
      float *compress = (float *)malloc(sizeof(float) * compress_length * 3);

      nc_get_vara_float(ncid, Vxid, startp, countp, compress_Vx);
      nc_get_vara_float(ncid, Vyid, startp, countp, compress_Vy);
      nc_get_vara_float(ncid, Vzid, startp, countp, compress_Vz);

      memcpy(compress, compress_Vx, sizeof(float) * compress_length);
      memcpy(compress + compress_length, compress_Vy, sizeof(float) * compress_length);
      memcpy(compress + compress_length * 2, compress_Vz, sizeof(float) * compress_length);

      sprintf(outputFilePath, "%s/compressed_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);

      string outfile = outputFilePath;

      string io_type = "float";
      Target target = eps;
      double target_value = tolerance;
      size_t skip_bytes = 0;
      bool verbose_flag = false;
      bool debug_flag = false;
      int dims = 5;
      int base[5] = {r1, r2, r3, r4, 3};

      // if (myid == 10)
      // {
      //   verbose_flag = true;
      // }

      MPI_Barrier(MPI_COMM_WORLD);

      compress_revised(compress, dims, base, outfile, io_type, target, target_value, skip_bytes, verbose_flag, debug_flag, compress_time);

      // MPI_Barrier(MPI_COMM_WORLD);
      // double compress_end = MPI_Wtime();
      // compression_time = compression_time + compress_end - compress_start;

      MPI_Barrier(MPI_COMM_WORLD);
      // double write_end = MPI_Wtime();
      // write_time = write_time + write_end - compress_end;

      if (reconstructed == 1)
      {
        bool autocrop = false;

        float *reconstructed_data = (float *)calloc(compress_length * 3, sizeof(float));

        // sprintf(rec_file, "reconstructed/reconstructed_velocity_px%d_py%d_index%d.nc", proc_x, proc_y, ii);

        // string reconstruct_file = rec_file;

        MPI_Barrier(MPI_COMM_WORLD);
        double decompress_start = MPI_Wtime();

        decompress(outfile, reconstructed_data, autocrop, verbose_flag, debug_flag);

        MPI_Barrier(MPI_COMM_WORLD);
        double decompress_end = MPI_Wtime();
        decompress_time[0] = decompress_time[0] + decompress_end - decompress_start;

        nc_put_vara_float(ncid_r, varid_r[0], startp, countp, reconstructed_data);
        nc_put_vara_float(ncid_r, varid_r[1], startp, countp, reconstructed_data + compress_length);
        nc_put_vara_float(ncid_r, varid_r[2], startp, countp, reconstructed_data + compress_length * 2);

        MPI_Barrier(MPI_COMM_WORLD);
        // printf("reconstructed successful");
        free(reconstructed_data);
        // free(reconstructed_Vy);
        // free(reconstructed_Vz);
      }

      free(compress);
      free(compress_Vx);
      free(compress_Vy);
      free(compress_Vz);
    }

    nc_close(ncid);
    nc_close(ncid_r);
  }
  /*
      if (stress == 1)
      {
          sprintf(fnm_init, "%s/volume_stress_px%d_py%d.nc", init_root, proc_x, proc_y);
          printf("%s\n", fnm_init);
          nc_open(fnm_init, NC_NOWRITE, &ncid);

          nc_inq_varid(ncid, "Txx", &Txxid);
          nc_inq_varid(ncid, "Tyy", &Tyyid);
          nc_inq_varid(ncid, "Tzz", &Tzzid);
          nc_inq_varid(ncid, "Txy", &Txyid);
          nc_inq_varid(ncid, "Txz", &Txzid);
          nc_inq_varid(ncid, "Tyz", &Tyzid);

          nc_inq_dimid(ncid, "xi", &xi_id);
          nc_inq_dimid(ncid, "eta", &eta_id);
          nc_inq_dimid(ncid, "zeta", &zeta_id);
          nc_inq_dimid(ncid, "time", &time_id);

          nc_inq_dimlen(ncid, xi_id, &xi_length);
          nc_inq_dimlen(ncid, eta_id, &eta_length);
          nc_inq_dimlen(ncid, zeta_id, &zeta_length);
          nc_inq_dimlen(ncid, time_id, &time_length);

          compress_length = xi_length * eta_length * zeta_length * r4;
          sum_length = xi_length * eta_length * zeta_length * time_length;

          float *Txx = (float *)malloc(sizeof(float) * sum_length);
          float *Tyy = (float *)malloc(sizeof(float) * sum_length);
          float *Tzz = (float *)malloc(sizeof(float) * sum_length);
          float *Txy = (float *)malloc(sizeof(float) * sum_length);
          float *Txz = (float *)malloc(sizeof(float) * sum_length);
          float *Tyz = (float *)malloc(sizeof(float) * sum_length);
          float *compress_Txx = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Tyy = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Tzz = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Txy = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Tyz = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Txz = (float *)malloc(sizeof(float) * compress_length);

          nc_get_var_float(ncid, Txxid, Txx);
          nc_get_var_float(ncid, Tyyid, Tyy);
          nc_get_var_float(ncid, Tzzid, Tzz);
          nc_get_var_float(ncid, Txyid, Txy);
          nc_get_var_float(ncid, Txzid, Txz);
          nc_get_var_float(ncid, Tyzid, Tyz);

          nc_close(ncid);

          r3 = zeta_length;
          r2 = eta_length;
          r1 = xi_length;

          SZ::Config conf;
          conf = SZ::Config(r4, r3, r2, r1);

          conf.loadcfg(conPath);
          size_t outSize_Txx, outSize_Tyy, outSize_Tzz;
          size_t outSize_Txy, outSize_Txz, outSize_Tyz;

          char outputFilePath_Txx[1024], outputFilePath_Tyy[1024], outputFilePath_Tzz[1024];
          char outputFilePath_Txy[1024], outputFilePath_Txz[1024], outputFilePath_Tyz[1024];

          compress_it_num = time_length / r4;

          for (int ii = 0; ii < compress_it_num; ii++)
          {
              memcpy(compress_Txx, Txx + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Tyy, Tyy + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Tzz, Tzz + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Txy, Txy + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Txz, Txz + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Tyz, Tyz + ii * compress_length, compress_length * sizeof(float));

              MPI_Barrier(MPI_COMM_WORLD);
              double compress_start = MPI_Wtime();

              char *bytes_Txx = SZ_compress<float>(conf, compress_Txx, outSize_Txx);
              char *bytes_Tyy = SZ_compress<float>(conf, compress_Tyy, outSize_Tyy);
              char *bytes_Tzz = SZ_compress<float>(conf, compress_Tzz, outSize_Tzz);
              char *bytes_Txy = SZ_compress<float>(conf, compress_Txy, outSize_Txy);
              char *bytes_Txz = SZ_compress<float>(conf, compress_Txz, outSize_Txz);
              char *bytes_Tyz = SZ_compress<float>(conf, compress_Tyz, outSize_Tyz);

              MPI_Barrier(MPI_COMM_WORLD);
              double compress_end = MPI_Wtime();
              compression_time = compression_time + compress_end - compress_start;

              sprintf(outputFilePath_Txx, "%s/compresse_Txx_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Txx, bytes_Txx, outSize_Txx);

              sprintf(outputFilePath_Tyy, "%s/compresse_Tyy_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Tyy, bytes_Tyy, outSize_Tyy);

              sprintf(outputFilePath_Tzz, "%s/compresse_Tzz_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Tzz, bytes_Tzz, outSize_Tzz);

              sprintf(outputFilePath_Txy, "%s/compresse_Txy_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Txy, bytes_Txy, outSize_Txy);

              sprintf(outputFilePath_Txz, "%s/compresse_Txz_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Txz, bytes_Txz, outSize_Txz);

              sprintf(outputFilePath_Tyz, "%s/compresse_Tyz_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Tyz, bytes_Tyz, outSize_Tyz);
              delete[] bytes_Txx;
              delete[] bytes_Tyy;
              delete[] bytes_Tzz;
              delete[] bytes_Txy;
              delete[] bytes_Txz;
              delete[] bytes_Tyz;

              MPI_Barrier(MPI_COMM_WORLD);
              double write_end = MPI_Wtime();
              write_time = write_time + write_end - compress_end;
          }
          free(Txx);
          free(Tyy);
          free(Tzz);
          free(Txy);
          free(Txz);
          free(Tyz);
          free(compress_Txx);
          free(compress_Tyy);
          free(compress_Tzz);
          free(compress_Txy);
          free(compress_Txz);
          free(compress_Tyz);
      }

      if (strain == 1)
      {
          sprintf(fnm_init, "%s/volume_strain_px%d_py%d.nc", init_root, proc_x, proc_y);
          printf("%s\n", fnm_init);
          nc_open(fnm_init, NC_NOWRITE, &ncid);

          nc_inq_varid(ncid, "Txx", &Exxid);
          nc_inq_varid(ncid, "Tyy", &Eyyid);
          nc_inq_varid(ncid, "Tzz", &Ezzid);
          nc_inq_varid(ncid, "Txy", &Exyid);
          nc_inq_varid(ncid, "Txz", &Exzid);
          nc_inq_varid(ncid, "Tyz", &Eyzid);

          nc_inq_dimid(ncid, "xi", &xi_id);
          nc_inq_dimid(ncid, "eta", &eta_id);
          nc_inq_dimid(ncid, "zeta", &zeta_id);
          nc_inq_dimid(ncid, "time", &time_id);

          nc_inq_dimlen(ncid, xi_id, &xi_length);
          nc_inq_dimlen(ncid, eta_id, &eta_length);
          nc_inq_dimlen(ncid, zeta_id, &zeta_length);
          nc_inq_dimlen(ncid, time_id, &time_length);

          compress_length = xi_length * eta_length * zeta_length * r4;
          sum_length = xi_length * eta_length * zeta_length * time_length;

          float *Exx = (float *)malloc(sizeof(float) * sum_length);
          float *Eyy = (float *)malloc(sizeof(float) * sum_length);
          float *Ezz = (float *)malloc(sizeof(float) * sum_length);
          float *Exy = (float *)malloc(sizeof(float) * sum_length);
          float *Exz = (float *)malloc(sizeof(float) * sum_length);
          float *Eyz = (float *)malloc(sizeof(float) * sum_length);
          float *compress_Exx = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Eyy = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Ezz = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Exy = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Eyz = (float *)malloc(sizeof(float) * compress_length);
          float *compress_Exz = (float *)malloc(sizeof(float) * compress_length);

          nc_get_var_float(ncid, Exxid, Exx);
          nc_get_var_float(ncid, Eyyid, Eyy);
          nc_get_var_float(ncid, Ezzid, Ezz);
          nc_get_var_float(ncid, Exyid, Exy);
          nc_get_var_float(ncid, Exzid, Exz);
          nc_get_var_float(ncid, Eyzid, Eyz);

          nc_close(ncid);

          r3 = zeta_length;
          r2 = eta_length;
          r1 = xi_length;

          SZ::Config conf;
          conf = SZ::Config(r4, r3, r2, r1);

          conf.loadcfg(conPath);
          size_t outSize_Exx, outSize_Eyy, outSize_Ezz;
          size_t outSize_Exy, outSize_Exz, outSize_Eyz;

          char outputFilePath_Exx[1024], outputFilePath_Eyy[1024], outputFilePath_Ezz[1024];
          char outputFilePath_Exy[1024], outputFilePath_Exz[1024], outputFilePath_Eyz[1024];

          compress_it_num = time_length / r4;

          for (int ii = 0; ii < compress_it_num; ii++)
          {
              memcpy(compress_Exx, Exx + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Eyy, Eyy + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Ezz, Ezz + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Exy, Exy + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Exz, Exz + ii * compress_length, compress_length * sizeof(float));
              memcpy(compress_Eyz, Eyz + ii * compress_length, compress_length * sizeof(float));

              MPI_Barrier(MPI_COMM_WORLD);
              double compress_start = MPI_Wtime();

              char *bytes_Exx = SZ_compress<float>(conf, compress_Exx, outSize_Exx);
              char *bytes_Eyy = SZ_compress<float>(conf, compress_Eyy, outSize_Eyy);
              char *bytes_Ezz = SZ_compress<float>(conf, compress_Ezz, outSize_Ezz);
              char *bytes_Exy = SZ_compress<float>(conf, compress_Exy, outSize_Exy);
              char *bytes_Exz = SZ_compress<float>(conf, compress_Exz, outSize_Exz);
              char *bytes_Eyz = SZ_compress<float>(conf, compress_Eyz, outSize_Eyz);

              MPI_Barrier(MPI_COMM_WORLD);
              double compress_end = MPI_Wtime();
              compression_time = compression_time + compress_end - compress_start;

              sprintf(outputFilePath_Exx, "%s/compresse_Exx_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Exx, bytes_Exx, outSize_Exx);

              sprintf(outputFilePath_Eyy, "%s/compresse_Eyy_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Eyy, bytes_Eyy, outSize_Eyy);

              sprintf(outputFilePath_Ezz, "%s/compresse_Ezz_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Ezz, bytes_Ezz, outSize_Ezz);

              sprintf(outputFilePath_Exy, "%s/compresse_Exy_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Exy, bytes_Exy, outSize_Exy);

              sprintf(outputFilePath_Exz, "%s/compresse_Exz_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Exz, bytes_Exz, outSize_Exz);

              sprintf(outputFilePath_Eyz, "%s/compresse_Eyz_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
              SZ::writefile(outputFilePath_Eyz, bytes_Eyz, outSize_Eyz);
              delete[] bytes_Exx;
              delete[] bytes_Eyy;
              delete[] bytes_Ezz;
              delete[] bytes_Exy;
              delete[] bytes_Exz;
              delete[] bytes_Eyz;

              MPI_Barrier(MPI_COMM_WORLD);
              double write_end = MPI_Wtime();
              write_time = write_time + write_end - compress_end;
          }
          free(Exx);
          free(Eyy);
          free(Ezz);
          free(Exy);
          free(Exz);
          free(Eyz);
          free(compress_Exx);
          free(compress_Eyy);
          free(compress_Ezz);
          free(compress_Exy);
          free(compress_Exz);
          free(compress_Eyz);
      }
  */

  if (myid == 0)
  {
    printf("compression time is %f\n", compress_time[0]);
    printf("decompression time is %f\n", decompress_time[0]);
  }
  // printf("written time is %f\n", write_time);
  // printf("reconstruct time is %f\n", reconstruct_time);
  MPI_Finalize();
  return 0;
}