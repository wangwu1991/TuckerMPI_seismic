#include <sys/time.h>
#include <stdio.h>
#include <stdlib.h>
#include "sz.h"
#include "rw.h"

#include "netcdf.h"
#include "cJSON.h"
#include <mpi.h>

int main(int argc, char *argv[])
{
    size_t r5 = 0;
    size_t r4 = 0;
    size_t r3 = 0;
    size_t r2 = 0;
    size_t r1 = 0;

    int status;

    int nprocx, nprocy;

    int myid, mpi_size;

    char init_root[1024];
    char out_root[1024];

    char fnm_init[1024];

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

    double compression_time, write_time, reconstruct_time;

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

        if (reconstructed == 1)
        {
            char rec_file[1024];
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

        r5 = 0;
        r3 = zeta_length;
        r2 = eta_length;
        r1 = xi_length;

        char cfgFile[1024] = "sz.config";

        int status = SZ_Init(cfgFile);

        size_t outSize_Vx, outSize_Vy, outSize_Vz;
        char outputFilePath_Vx[1024], outputFilePath_Vy[1024], outputFilePath_Vz[1024];

        compress_it_num = time_length / r4;

        compress_length = r1 * r2 * r3 * r4;
        float *compress_Vx = (float *)malloc(sizeof(float) * compress_length);
        float *compress_Vy = (float *)malloc(sizeof(float) * compress_length);
        float *compress_Vz = (float *)malloc(sizeof(float) * compress_length);

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
            size_t countp[] = {r4, r3, r2, r1};

            printf("startp = %d * %d\n", ii, r4);

            nc_get_vara_float(ncid, Vxid, startp, countp, compress_Vx);
            nc_get_vara_float(ncid, Vyid, startp, countp, compress_Vy);
            nc_get_vara_float(ncid, Vzid, startp, countp, compress_Vz);

            MPI_Barrier(MPI_COMM_WORLD);
            double compress_start = MPI_Wtime();

            unsigned char *bytes_Vx = SZ_compress(SZ_FLOAT, compress_Vx, &outSize_Vx, r5, r4, r3, r2, r1);
            unsigned char *bytes_Vy = SZ_compress(SZ_FLOAT, compress_Vy, &outSize_Vy, r5, r4, r3, r2, r1);
            unsigned char *bytes_Vz = SZ_compress(SZ_FLOAT, compress_Vz, &outSize_Vz, r5, r4, r3, r2, r1);

            MPI_Barrier(MPI_COMM_WORLD);
            double compress_end = MPI_Wtime();
            compression_time = compression_time + compress_end - compress_start;

            sprintf(outputFilePath_Vx, "%s/compressed_Vx_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
            writeByteData(bytes_Vx, outSize_Vx, outputFilePath_Vx, &status);

            sprintf(outputFilePath_Vy, "%s/compressed_Vy_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
            writeByteData(bytes_Vy, outSize_Vy, outputFilePath_Vy, &status);

            sprintf(outputFilePath_Vz, "%s/compressed_Vz_px%d_py%d_index%d.sz", out_root, proc_x, proc_y, ii);
            writeByteData(bytes_Vz, outSize_Vz, outputFilePath_Vz, &status);

            MPI_Barrier(MPI_COMM_WORLD);
            double write_end = MPI_Wtime();
            write_time = write_time + write_end - compress_end;

            if (reconstructed == 1)
            {
                float *reconstructed_Vx = SZ_decompress(SZ_FLOAT, bytes_Vx, outSize_Vx, r5, r4, r3, r2, r1);
                float *reconstructed_Vy = SZ_decompress(SZ_FLOAT, bytes_Vy, outSize_Vy, r5, r4, r3, r2, r1);
                float *reconstructed_Vz = SZ_decompress(SZ_FLOAT, bytes_Vz, outSize_Vz, r5, r4, r3, r2, r1);

                MPI_Barrier(MPI_COMM_WORLD);
                double reconstruct_end = MPI_Wtime();
                reconstruct_time = reconstruct_time + reconstruct_end - write_end;

                nc_put_vara_float(ncid_r, varid_r[0], startp, countp, reconstructed_Vx);
                nc_put_vara_float(ncid_r, varid_r[1], startp, countp, reconstructed_Vy);
                nc_put_vara_float(ncid_r, varid_r[2], startp, countp, reconstructed_Vz);

                for (int jj = 0; jj < compress_length; jj++)
                {
                    if (min_Vx > compress_Vx[jj])
                    {
                        min_Vx = compress_Vx[jj];
                    }
                    if (max_Vx < compress_Vx[jj])
                    {
                        max_Vx = compress_Vx[jj];
                    }

                    if (min_Vy > compress_Vy[jj])
                    {
                        min_Vy = compress_Vy[jj];
                    }
                    if (max_Vy < compress_Vy[jj])
                    {
                        max_Vy = compress_Vy[jj];
                    }

                    if (min_Vz > compress_Vz[jj])
                    {
                        min_Vz = compress_Vz[jj];
                    }
                    if (max_Vz < compress_Vz[jj])
                    {
                        max_Vz = compress_Vz[jj];
                    }

                    float err_x = reconstructed_Vx[jj] - compress_Vx[jj];
                    float err_y = reconstructed_Vy[jj] - compress_Vy[jj];
                    float err_z = reconstructed_Vz[jj] - compress_Vz[jj];

                    sum_Vx += (long double)(err_x * err_x);
                    sum_Vy += (long double)(err_y * err_y);
                    sum_Vz += (long double)(err_z * err_z);
                }

                free(reconstructed_Vx);
                free(reconstructed_Vy);
                free(reconstructed_Vz);
            }
            free(bytes_Vx);
            free(bytes_Vy);
            free(bytes_Vz);
        }
        free(compress_Vx);
        free(compress_Vy);
        free(compress_Vz);

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

            long double mse_Vx = global_sum_Vx / global_nbEle;
            long double mse_Vy = global_sum_Vy / global_nbEle;
            long double mse_Vz = global_sum_Vz / global_nbEle;

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

            long double names = sqrtl(mse) / range;

            printf("PSNR_Vx = %Lf, NRMSE_Vx = %Lf\n", psnr_Vx, nrmse_Vx);
            printf("PSNR_Vy = %Lf, NRMSE_Vy = %Lf\n", psnr_Vy, nrmse_Vy);
            printf("PSNR_Vz = %Lf, NRMSE_Vz = %Lf\n", psnr_Vz, nrmse_Vz);
            printf("PSNR = %Lf, NRMSE = %Lf\n", psnr, names);
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

    printf("compression time is %f\n", compression_time);
    printf("written time is %f\n", write_time);
    printf("reconstruct time is %f\n", reconstruct_time);
    MPI_Finalize();
    return 0;
}
