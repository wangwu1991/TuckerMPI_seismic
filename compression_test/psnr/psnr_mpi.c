#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "netcdf.h"
#include "cJSON.h"
#include <mpi.h>

#include <unistd.h>

int main(int argc, char *argv[])
{
    int nprocx, nprocy;

    int myid, mpi_size;

    char init_root[1024];
    char out_root[1024];

    char fnm_init[1024];

    int r4;

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
    // if (item = cJSON_GetObjectItem(root_wave, "tolerance"))
    // {
    //     tolerance = item->valuedouble;
    // }
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
    int Vxid, Vyid, Vzid;
    int Vxid_r, Vyid_r, Vzid_r;
    int Txxid, Tyyid, Tzzid, Txyid, Txzid, Tyzid;
    int Exxid, Eyyid, Ezzid, Exyid, Exzid, Eyzid;
    int xi_id, eta_id, zeta_id, time_id;
    size_t xi_length, eta_length, zeta_length, time_length, sum_length, compress_length;
    int compress_it_num;

    int ncid_r;
    int dim_id_r[4];
    int varid_r[3];

    // double compression_time, write_time, reconstruct_time;

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

        compress_length = xi_length * eta_length * zeta_length * r4;

        float *init_Vx = (float *)malloc(sizeof(float) * compress_length);
        float *init_Vy = (float *)malloc(sizeof(float) * compress_length);
        float *init_Vz = (float *)malloc(sizeof(float) * compress_length);

        float *rect_Vx = (float *)malloc(sizeof(float) * compress_length);
        float *rect_Vy = (float *)malloc(sizeof(float) * compress_length);
        float *rect_Vz = (float *)malloc(sizeof(float) * compress_length);

        char rec_file[1024];
        sprintf(rec_file, "reconstructed/reconstructed_velocity_px%d_py%d.nc", proc_x, proc_y);

        nc_open(rec_file, NC_NOWRITE, &ncid_r);
        nc_inq_varid(ncid_r, "Vx", &Vxid_r);
        nc_inq_varid(ncid_r, "Vy", &Vyid_r);
        nc_inq_varid(ncid_r, "Vz", &Vzid_r);

        compress_it_num = time_length / r4;

        float min_Vx = 0, max_Vx = 0, global_max_Vx, global_min_Vx;
        float min_Vy = 0, max_Vy = 0, global_max_Vy, global_min_Vy;
        float min_Vz = 0, max_Vz = 0, global_max_Vz, global_min_Vz;
        long double sum = 0, sum_Vz = 0, sum_Vy = 0, sum_Vx = 0;

        long double global_sum_Vx = 0;
        long double global_sum_Vy = 0;
        long double global_sum_Vz = 0;

        float err_x, err_y, err_z;
        unsigned long int nbEle, global_nbEle;

        for (int ii = 0; ii < compress_it_num; ii++)
        {

            size_t startp[] = {ii * r4, 0, 0, 0};
            size_t countp[] = {r4, zeta_length, eta_length, xi_length};

            printf("startp = %d * %d\n", ii, r4);

            nc_get_vara_float(ncid, Vxid, startp, countp, init_Vx);
            nc_get_vara_float(ncid, Vyid, startp, countp, init_Vy);
            nc_get_vara_float(ncid, Vzid, startp, countp, init_Vz);

            nc_get_vara_float(ncid_r, Vxid_r, startp, countp, rect_Vx);
            nc_get_vara_float(ncid_r, Vyid_r, startp, countp, rect_Vy);
            nc_get_vara_float(ncid_r, Vzid_r, startp, countp, rect_Vz);

            for (int jj = 0; jj < compress_length; jj++)
            {
                if (min_Vx > init_Vx[jj])
                {
                    min_Vx = init_Vx[jj];
                }
                if (max_Vx < init_Vx[jj])
                {
                    max_Vx = init_Vx[jj];
                }

                if (min_Vy > init_Vy[jj])
                {
                    min_Vy = init_Vy[jj];
                }
                if (max_Vy < init_Vy[jj])
                {
                    max_Vy = init_Vy[jj];
                }

                if (min_Vz > init_Vz[jj])
                {
                    min_Vz = init_Vz[jj];
                }
                if (max_Vz < init_Vz[jj])
                {
                    max_Vz = init_Vz[jj];
                }

                if (isnormal(rect_Vx[jj]) != 1)
                {
                    rect_Vx[jj] = 0;
                }
                if (isnormal(rect_Vy[jj]) != 1)
                {
                    rect_Vy[jj] = 0;
                }
                if (isnormal(rect_Vz[jj]) != 1)
                {
                    rect_Vz[jj] = 0;
                }

                float err_x = rect_Vx[jj] - init_Vx[jj];
                float err_y = rect_Vy[jj] - init_Vy[jj];
                float err_z = rect_Vz[jj] - init_Vz[jj];

                sum_Vx += (long double)(err_x * err_x);
                sum_Vy += (long double)(err_y * err_y);
                sum_Vz += (long double)(err_z * err_z);
            }
        }

        nbEle = compress_it_num * compress_length;

        MPI_Barrier(MPI_COMM_WORLD);

        printf("compress_it_num = %d\n", compress_it_num);
        printf("compress_length = %d\n", compress_length);
        printf("nbEle = %ld\n", nbEle);

        printf("sum_Vx = %Lf, min_Vx = %f, max_Vx = %f\n", sum_Vx, min_Vx, max_Vx);
        printf("sum_Vy = %Lf, min_Vy = %f, max_Vy = %f\n", sum_Vy, min_Vy, max_Vy);
        printf("sum_Vz = %Lf, min_Vz = %f, max_Vz = %f\n", sum_Vz, min_Vz, max_Vz);
        printf("\n\n");

        MPI_Barrier(MPI_COMM_WORLD);

        MPI_Reduce(&max_Vx, &global_max_Vx, 1, MPI_FLOAT, MPI_MAX, 0, comm);
        MPI_Reduce(&max_Vy, &global_max_Vy, 1, MPI_FLOAT, MPI_MAX, 0, comm);
        MPI_Reduce(&max_Vz, &global_max_Vz, 1, MPI_FLOAT, MPI_MAX, 0, comm);

        MPI_Reduce(&min_Vx, &global_min_Vx, 1, MPI_FLOAT, MPI_MIN, 0, comm);
        MPI_Reduce(&min_Vy, &global_min_Vy, 1, MPI_FLOAT, MPI_MIN, 0, comm);
        MPI_Reduce(&min_Vz, &global_min_Vz, 1, MPI_FLOAT, MPI_MIN, 0, comm);

        MPI_Reduce(&sum_Vx, &global_sum_Vx, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&sum_Vy, &global_sum_Vy, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, comm);
        MPI_Reduce(&sum_Vz, &global_sum_Vz, 1, MPI_LONG_DOUBLE, MPI_SUM, 0, comm);

        MPI_Reduce(&nbEle, &global_nbEle, 1, MPI_UNSIGNED_LONG, MPI_SUM, 0, comm);

        if (myid == 0)
        {
            printf("##################\n");
            // printf("global_sum = %lf, global_min = %lf, global_max = %lf\n", global_sum, global_min, global_max);
            printf("global_nbEle = %ld\n", global_nbEle);
            printf("global_sum_Vx = %Lf, global_min_Vx = %f, global_max_Vx = %f\n", global_sum_Vx, global_min_Vx, global_max_Vx);
            printf("global_sum_Vy = %Lf, global_min_Vy = %f, global_max_Vy = %f\n", global_sum_Vy, global_min_Vy, global_max_Vy);
            printf("global_sum_Vz = %Lf, global_min_Vz = %f, global_max_Vz = %f\n", global_sum_Vz, global_min_Vz, global_max_Vz);

            long double range_Vx = (long double)(global_max_Vx - global_min_Vx);
            long double range_Vy = (long double)(global_max_Vy - global_min_Vy);
            long double range_Vz = (long double)(global_max_Vz - global_min_Vz);

            long double mse_Vx = (global_sum_Vx / global_nbEle);
            long double mse_Vy = (global_sum_Vy / global_nbEle);
            long double mse_Vz = (global_sum_Vz / global_nbEle);

            long double global_sum = global_sum_Vx + global_sum_Vy + global_sum_Vz;

            float global_max = global_max_Vx, global_min = global_min_Vx;
            if (global_max <= global_max_Vy)
            {
                global_max = global_max_Vy;
            }
            if (global_max <= global_max_Vz)
            {
                global_max = global_max_Vz;
            }

            if (global_min >= global_min_Vy)
            {
                global_min = global_min_Vy;
            }
            if (global_min >= global_min_Vz)
            {
                global_min = global_min_Vz;
            }

            long double range = (long double)(global_max - global_min);
            long double mse = global_sum / global_nbEle / 3;

            long double psnr_Vx = 20 * log10l(range_Vx) - 10 * log10l(mse_Vx);
            long double psnr_Vy = 20 * log10l(range_Vy) - 10 * log10l(mse_Vy);
            long double psnr_Vz = 20 * log10l(range_Vz) - 10 * log10l(mse_Vz);

            long double psnr = 20 * log10l(range) - 10 * log10l(mse);

            long double nrmse_Vx = sqrtl(mse_Vx) / range_Vx;
            long double nrmse_Vy = sqrtl(mse_Vy) / range_Vy;
            long double nrmse_Vz = sqrtl(mse_Vz) / range_Vz;

            long double nrmse = sqrt(mse) / range;

            printf("PSNR_Vx = %Lf, NRMSE_Vx = %Lf\n", psnr_Vx, nrmse_Vx);
            printf("PSNR_Vy = %Lf, NRMSE_Vy = %Lf\n", psnr_Vy, nrmse_Vy);
            printf("PSNR_Vz = %Lf, NRMSE_Vz = %Lf\n", psnr_Vz, nrmse_Vz);
            printf("PSNR = %Lf, NRMSE = %Lf\n", psnr, nrmse);
        }

        free(init_Vx);
        free(init_Vy);
        free(init_Vz);

        free(rect_Vx);
        free(rect_Vy);
        free(rect_Vz);

        nc_close(ncid);
        nc_close(ncid_r);
    }

    MPI_Finalize();
    return 0;
}