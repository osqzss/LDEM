#ifndef LDEM_H
#define LDEM_H

#ifdef __cplusplus
extern "C" {
#endif

/*
 * Lunar south polar stereographic projection parameters.
 *
 * The Moon radius matches the value used by the LOLA LDEM products.
 * Georeferencing (origin, pixel size) is read from the GeoTIFF metadata
 * by ldem_load() and stored in LDEMGrid; it is not hardcoded here.
 */
#define LDEM_MOON_RADIUS 1737400.0

typedef struct {
    double x;
    double y;
    double z;
} LDEMVec3;

typedef struct {
    int width;
    int height;

    /*
     * Georeferencing parameters read from the GeoTIFF metadata.
     *
     * origin_x   : X coordinate of the upper-left corner of the upper-left pixel [m]
     * origin_y   : Y coordinate of the upper-left corner of the upper-left pixel [m]
     * pixel_size : ground sampling distance (square pixels assumed) [m]
     */
    double origin_x;
    double origin_y;
    double pixel_size;

    /*
     * DEM height array in meters.
     * Row-major order:
     *
     *   z[row * width + col]
     */
    float *z;
} LDEMGrid;

/*
 * Load a single-band Float32 GeoTIFF LDEM file into memory using GDAL.
 *
 * Georeferencing parameters (origin_x, origin_y, pixel_size) are read
 * from the file's embedded geotransform.  The geotransform must describe
 * an axis-aligned, north-up raster with square pixels (no rotation or shear).
 *
 * Returns 1 on success, 0 on failure.
 */
int ldem_load(const char *path, LDEMGrid *dem);

/*
 * Release memory allocated by ldem_load().
 */
void ldem_free(LDEMGrid *dem);

/*
 * Convert lunar latitude/longitude to LOLA south polar stereographic X/Y.
 *
 * Input:
 *   lat_deg : planetocentric latitude [deg]
 *   lon_deg : east-positive longitude [deg]
 *
 * Output:
 *   x, y    : projected coordinates [m]
 *
 * Returns 1 on success, 0 on failure.
 */
int ldem_latlon_to_xy(double lat_deg, double lon_deg,
                      double *x, double *y);

/*
 * Convert projected X/Y coordinates to fractional pixel/line coordinates.
 *
 * pixel and line are expressed with respect to pixel-centered samples.
 *
 * Returns 1 on success, 0 if the point is outside the DEM.
 */
int ldem_xy_to_pixel(const LDEMGrid *dem,
                     double x, double y,
                     double *pixel, double *line);

/*
 * Get DEM height at projected X/Y coordinates using bilinear interpolation.
 *
 * Input:
 *   x, y   : LOLA south polar stereographic coordinates [m]
 *
 * Output:
 *   height : terrain height [m]
 *
 * Returns 1 on success, 0 on failure.
 */
int ldem_height_xy(const LDEMGrid *dem,
                   double x, double y,
                   double *height);

/*
 * Get DEM height at latitude/longitude using bilinear interpolation.
 *
 * Returns 1 on success, 0 on failure.
 */
int ldem_height_latlon(const LDEMGrid *dem,
                       double lat_deg, double lon_deg,
                       double *height);

/*
 * Compute the local surface normal at projected X/Y coordinates.
 *
 * The normal vector is defined in the local LDEM coordinate system:
 *
 *   n = normalize(-dh/dx, -dh/dy, 1)
 *
 * where x and y are the LDEM polar stereographic coordinates and h is
 * the DEM height.
 *
 * step is the central-difference interval [m]. For the 80 m/pixel LDEM,
 * a value of 160 to 400 m is usually reasonable.
 *
 * Returns 1 on success, 0 on failure.
 */
int ldem_surface_normal_xy(const LDEMGrid *dem,
                           double x, double y,
                           double step,
                           LDEMVec3 *n_out);

/*
 * Compute the local surface normal at latitude/longitude.
 *
 * If x_out and y_out are not NULL, the projected coordinates are also returned.
 *
 * Returns 1 on success, 0 on failure.
 */
int ldem_surface_normal_latlon(const LDEMGrid *dem,
                               double lat_deg, double lon_deg,
                               double step,
                               LDEMVec3 *n_out,
                               double *x_out,
                               double *y_out);

/*
 * Compute the terrain slope angle from a surface normal.
 *
 * Returns the slope angle in radians.
 */
double ldem_slope_from_normal(const LDEMVec3 *n);

#ifdef __cplusplus
}
#endif

#endif /* LDEM_H */
