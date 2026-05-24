#include "ldem.h"

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <gdal.h>     /* sudo apt install libgdal-dev */
#include <cpl_conv.h>

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD (M_PI / 180.0)

/* --------------------------------------------------------------------------
 * Internal helpers
 * -------------------------------------------------------------------------- */

static float ldem_get_pixel(const LDEMGrid *dem, int col, int row)
{
    return dem->z[row * dem->width + col];
}

/* --------------------------------------------------------------------------
 * Load / free
 * -------------------------------------------------------------------------- */

int ldem_load(const char *path, LDEMGrid *dem)
{
    GDALDatasetH  ds   = NULL;
    GDALRasterBandH band = NULL;
    double        gt[6];
    double        px_w, px_h;
    CPLErr        err;

    if (dem == NULL) {
        return 0;
    }

    dem->width      = 0;
    dem->height     = 0;
    dem->origin_x   = 0.0;
    dem->origin_y   = 0.0;
    dem->pixel_size = 0.0;
    dem->z          = NULL;

    /*
     * Initialise the GDAL driver registry.
     * Calling this multiple times is safe (no-op after the first call).
     */
    GDALAllRegister();

    ds = GDALOpen(path, GA_ReadOnly);
    if (ds == NULL) {
        fprintf(stderr, "Error: GDAL failed to open file: %s\n", path);
        return 0;
    }

    /* ------------------------------------------------------------------
     * Validate band count and data type.
     * ------------------------------------------------------------------ */
    if (GDALGetRasterCount(ds) != 1) {
        fprintf(stderr, "Error: expected a single-band raster, got %d bands.\n",
                GDALGetRasterCount(ds));
        GDALClose(ds);
        return 0;
    }

    band = GDALGetRasterBand(ds, 1);

    if (GDALGetRasterDataType(band) != GDT_Float32) {
        fprintf(stderr, "Error: expected Float32 raster, got %s.\n",
                GDALGetDataTypeName(GDALGetRasterDataType(band)));
        GDALClose(ds);
        return 0;
    }

    /* ------------------------------------------------------------------
     * Read and validate the geotransform.
     *
     * GDAL geotransform convention:
     *   gt[0]  X of the upper-left corner of the upper-left pixel
     *   gt[1]  pixel width  (positive eastward)
     *   gt[2]  row rotation (0 for north-up)
     *   gt[3]  Y of the upper-left corner of the upper-left pixel
     *   gt[4]  column rotation (0 for north-up)
     *   gt[5]  pixel height (negative for north-up rasters)
     * ------------------------------------------------------------------ */
    if (GDALGetGeoTransform(ds, gt) != CE_None) {
        fprintf(stderr, "Error: failed to read geotransform from: %s\n", path);
        GDALClose(ds);
        return 0;
    }

    if (gt[2] != 0.0 || gt[4] != 0.0) {
        fprintf(stderr, "Error: rotated or sheared rasters are not supported "
                "(gt[2]=%g, gt[4]=%g).\n", gt[2], gt[4]);
        GDALClose(ds);
        return 0;
    }

    if (gt[1] <= 0.0) {
        fprintf(stderr, "Error: non-positive pixel width: %g\n", gt[1]);
        GDALClose(ds);
        return 0;
    }

    /*
     * Require square pixels: pixel width and the absolute pixel height
     * must agree to within a small relative tolerance.
     */
    px_w = gt[1];
    px_h = fabs(gt[5]);
    if (fabs(px_w - px_h) / px_w > 1e-6) {
        fprintf(stderr, "Error: non-square pixels are not supported "
                "(pixel_width=%.6f, pixel_height=%.6f).\n", px_w, px_h);
        GDALClose(ds);
        return 0;
    }

    dem->width      = GDALGetRasterXSize(ds);
    dem->height     = GDALGetRasterYSize(ds);
    dem->origin_x   = gt[0];
    dem->origin_y   = gt[3];
    dem->pixel_size = px_w;

    /* ------------------------------------------------------------------
     * Allocate height array and read all pixels in one call.
     * ------------------------------------------------------------------ */
    size_t n = (size_t)dem->width * (size_t)dem->height;

    dem->z = (float *)malloc(n * sizeof(float));
    if (dem->z == NULL) {
        fprintf(stderr, "Error: failed to allocate DEM memory: %.1f MB\n",
                (double)(n * sizeof(float)) / (1024.0 * 1024.0));
        GDALClose(ds);
        return 0;
    }

    /*
     * GDALRasterIO reads the entire band into the pre-allocated buffer.
     * nPixelSpace = 0 and nLineSpace = 0 let GDAL choose the tightest
     * packing, which matches the row-major float layout of dem->z.
     */
    err = GDALRasterIO(band,
                       GF_Read,
                       0, 0,                /* xoff, yoff   */
                       dem->width, dem->height, /* xsize, ysize */
                       dem->z,
                       dem->width, dem->height, /* buf_xsize, buf_ysize */
                       GDT_Float32,
                       0, 0);              /* nPixelSpace, nLineSpace */

    if (err != CE_None) {
        fprintf(stderr, "Error: GDALRasterIO failed to read raster data.\n");
        free(dem->z);
        dem->z = NULL;
        GDALClose(ds);
        return 0;
    }

    GDALClose(ds);

    fprintf(stderr,
            "Loaded DEM: %d x %d pixels, pixel_size=%.3f m, "
            "origin=(%.3f, %.3f), %.1f MB\n",
            dem->width, dem->height,
            dem->pixel_size,
            dem->origin_x, dem->origin_y,
            (double)(n * sizeof(float)) / (1024.0 * 1024.0));

    return 1;
}

void ldem_free(LDEMGrid *dem)
{
    if (dem == NULL) {
        return;
    }

    free(dem->z);
    dem->z      = NULL;
    dem->width  = 0;
    dem->height = 0;
    dem->origin_x   = 0.0;
    dem->origin_y   = 0.0;
    dem->pixel_size = 0.0;
}

/* --------------------------------------------------------------------------
 * Coordinate conversion
 * -------------------------------------------------------------------------- */

int ldem_latlon_to_xy(double lat_deg, double lon_deg,
                      double *x, double *y)
{
    if (x == NULL || y == NULL) {
        return 0;
    }

    /*
     * LDEM_80S_80MPP_ADJ.TIF covers the south polar region.
     * The map corners correspond to approximately 75.89 deg S.
     * A small margin is allowed here.
     */
    if (lat_deg < -90.0 || lat_deg > -75.0) {
        return 0;
    }

    double lat = lat_deg * DEG2RAD;
    double lon = lon_deg * DEG2RAD;

    /*
     * Moon south polar stereographic projection:
     *
     *   latitude of natural origin  = -90 deg
     *   longitude of natural origin =   0 deg
     *   scale factor                =   1
     *   false easting/northing      =   0 m
     *
     * X positive direction corresponds to 90 deg E.
     * Y positive direction corresponds to 0 deg longitude.
     */
    double rho = 2.0 * LDEM_MOON_RADIUS * tan((M_PI / 2.0 + lat) / 2.0);

    *x = rho * sin(lon);
    *y = rho * cos(lon);

    return 1;
}

int ldem_xy_to_pixel(const LDEMGrid *dem,
                     double x, double y,
                     double *pixel, double *line)
{
    if (dem == NULL || pixel == NULL || line == NULL) {
        return 0;
    }

    if (dem->width <= 0 || dem->height <= 0 || dem->pixel_size <= 0.0) {
        return 0;
    }

    double x_min = dem->origin_x;
    double x_max = dem->origin_x + dem->pixel_size * dem->width;
    double y_max = dem->origin_y;
    double y_min = dem->origin_y - dem->pixel_size * dem->height;

    if (x < x_min || x > x_max || y < y_min || y > y_max) {
        return 0;
    }

    /*
     * GeoTIFF origin is the upper-left corner of the upper-left pixel.
     * DEM values are treated as pixel-centered samples.
     */
    *pixel = (x - dem->origin_x) / dem->pixel_size - 0.5;
    *line  = (dem->origin_y - y) / dem->pixel_size - 0.5;

    return 1;
}

/* --------------------------------------------------------------------------
 * Height query
 * -------------------------------------------------------------------------- */

int ldem_height_xy(const LDEMGrid *dem,
                   double x, double y,
                   double *height)
{
    if (dem == NULL || dem->z == NULL || height == NULL) {
        return 0;
    }

    double pixel, line;

    if (!ldem_xy_to_pixel(dem, x, y, &pixel, &line)) {
        return 0;
    }

    int col0 = (int)floor(pixel);
    int row0 = (int)floor(line);
    int col1 = col0 + 1;
    int row1 = row0 + 1;

    if (col0 < 0 || row0 < 0 ||
        col1 >= dem->width || row1 >= dem->height) {
        return 0;
    }

    float z00_f = ldem_get_pixel(dem, col0, row0);
    float z10_f = ldem_get_pixel(dem, col1, row0);
    float z01_f = ldem_get_pixel(dem, col0, row1);
    float z11_f = ldem_get_pixel(dem, col1, row1);

    if (isnan(z00_f) || isnan(z10_f) ||
        isnan(z01_f) || isnan(z11_f)) {
        return 0;
    }

    double dx = pixel - (double)col0;
    double dy = line  - (double)row0;

    double z0 = (double)z00_f * (1.0 - dx) + (double)z10_f * dx;
    double z1 = (double)z01_f * (1.0 - dx) + (double)z11_f * dx;

    *height = z0 * (1.0 - dy) + z1 * dy;

    return 1;
}

int ldem_height_latlon(const LDEMGrid *dem,
                       double lat_deg, double lon_deg,
                       double *height)
{
    double x, y;

    if (!ldem_latlon_to_xy(lat_deg, lon_deg, &x, &y)) {
        return 0;
    }

    return ldem_height_xy(dem, x, y, height);
}

/* --------------------------------------------------------------------------
 * Surface normal
 * -------------------------------------------------------------------------- */

int ldem_surface_normal_xy(const LDEMGrid *dem,
                           double x, double y,
                           double step,
                           LDEMVec3 *n_out)
{
    if (dem == NULL || n_out == NULL || step <= 0.0) {
        return 0;
    }

    double hx_p, hx_m, hy_p, hy_m;

    /*
     * Central-difference terrain gradients:
     *
     *   dh/dx = [h(x+step, y) - h(x-step, y)] / (2 * step)
     *   dh/dy = [h(x, y+step) - h(x, y-step)] / (2 * step)
     */
    if (!ldem_height_xy(dem, x + step, y, &hx_p)) return 0;
    if (!ldem_height_xy(dem, x - step, y, &hx_m)) return 0;
    if (!ldem_height_xy(dem, x, y + step, &hy_p)) return 0;
    if (!ldem_height_xy(dem, x, y - step, &hy_m)) return 0;

    double dhdx = (hx_p - hx_m) / (2.0 * step);
    double dhdy = (hy_p - hy_m) / (2.0 * step);

    /*
     * Surface: z = h(x, y)
     * Upward normal: n = normalize(-dh/dx, -dh/dy, 1)
     */
    double nx = -dhdx;
    double ny = -dhdy;
    double nz = 1.0;

    double norm = sqrt(nx * nx + ny * ny + nz * nz);
    if (norm <= 0.0) {
        return 0;
    }

    n_out->x = nx / norm;
    n_out->y = ny / norm;
    n_out->z = nz / norm;

    return 1;
}

int ldem_surface_normal_latlon(const LDEMGrid *dem,
                               double lat_deg, double lon_deg,
                               double step,
                               LDEMVec3 *n_out,
                               double *x_out,
                               double *y_out)
{
    double x, y;

    if (!ldem_latlon_to_xy(lat_deg, lon_deg, &x, &y)) {
        return 0;
    }

    if (x_out) *x_out = x;
    if (y_out) *y_out = y;

    return ldem_surface_normal_xy(dem, x, y, step, n_out);
}

/* --------------------------------------------------------------------------
 * Utility
 * -------------------------------------------------------------------------- */

double ldem_slope_from_normal(const LDEMVec3 *n)
{
    if (n == NULL) {
        return 0.0;
    }

    double nz = n->z;

    if      (nz >  1.0) nz =  1.0;
    else if (nz < -1.0) nz = -1.0;

    return acos(nz);
}
