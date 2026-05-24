/*
 * shackleton_edge_mask.c
 *
 * Generate an azimuth-dependent elevation mask CSV for a receiver located
 * near the center of Shackleton crater, using LDEM_80S_80MPP_ADJ.TIF and
 * the provided ldem.c / ldem.h functions.
 *
 * Output CSV format for the MATLAB skyplot function:
 *   azimuth_deg,elevation_mask_deg
 *
 * Optional verbose output:
 *   azimuth_deg,elevation_mask_deg,edge_range_m,edge_lat_deg,edge_lon_deg,
 *   edge_dem_height_m,rx_dem_height_m
 *
 * Default receiver:
 *   rx_lat = -89.67 deg
 *   rx_lon = 129.78 deg
 *   h_ant  = 2.0 m above the DEM surface
 *
 * Method:
 *   For each azimuth, sample terrain points radially outward from the receiver.
 *   The terrain point with the maximum apparent elevation angle is used as the
 *   effective crater-edge horizon/elevation mask in that azimuth direction.
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <float.h>

#include "ldem.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

#define DEG2RAD (M_PI / 180.0)
#define RAD2DEG (180.0 / M_PI)

static const double RX_LAT_DEG_DEFAULT  = -89.67;
static const double RX_LON_DEG_DEFAULT  = 129.78;
static const double H_ANT_M_DEFAULT     = 2.0;
static const double AZ_STEP_DEG_DEFAULT = 5.0;

/* Shackleton radius is roughly 10 km. Search slightly beyond it. */
static const double R_MIN_M_DEFAULT  = 100.0;
static const double R_MAX_M_DEFAULT  = 20000.0;
static const double R_STEP_M_DEFAULT = 20.0;

static double clamp(double x, double lo, double hi)
{
    if (x < lo) return lo;
    if (x > hi) return hi;
    return x;
}

static double wrap_lon_0_360(double lon_deg)
{
    while (lon_deg < 0.0) lon_deg += 360.0;
    while (lon_deg >= 360.0) lon_deg -= 360.0;
    return lon_deg;
}

static void llh_to_moon_fixed(double lat_deg, double lon_deg, double h_m,
                              double r[3])
{
    const double lat = lat_deg * DEG2RAD;
    const double lon = lon_deg * DEG2RAD;
    const double radius = LDEM_MOON_RADIUS + h_m;
    const double clat = cos(lat);

    r[0] = radius * clat * cos(lon);
    r[1] = radius * clat * sin(lon);
    r[2] = radius * sin(lat);
}

static void local_enu_basis(double lat_deg, double lon_deg,
                            double east[3], double north[3], double up[3])
{
    const double lat = lat_deg * DEG2RAD;
    const double lon = lon_deg * DEG2RAD;
    const double slat = sin(lat);
    const double clat = cos(lat);
    const double slon = sin(lon);
    const double clon = cos(lon);

    east[0] = -slon;
    east[1] =  clon;
    east[2] =  0.0;

    north[0] = -slat * clon;
    north[1] = -slat * slon;
    north[2] =  clat;

    up[0] = clat * clon;
    up[1] = clat * slon;
    up[2] = slat;
}

static double dot3(const double a[3], const double b[3])
{
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

/* Destination point on a sphere from lat/lon, azimuth clockwise from north,
 * and surface distance.  This is sufficiently accurate for 10--20 km scales.
 */
static void destination_point_sphere(double lat1_deg, double lon1_deg,
                                     double az_deg, double dist_m,
                                     double *lat2_deg, double *lon2_deg)
{
    const double lat1 = lat1_deg * DEG2RAD;
    const double lon1 = lon1_deg * DEG2RAD;
    const double az = az_deg * DEG2RAD;
    const double d = dist_m / LDEM_MOON_RADIUS;

    const double sin_lat1 = sin(lat1);
    const double cos_lat1 = cos(lat1);
    const double sin_d = sin(d);
    const double cos_d = cos(d);

    const double sin_lat2 = sin_lat1 * cos_d + cos_lat1 * sin_d * cos(az);
    const double lat2 = asin(clamp(sin_lat2, -1.0, 1.0));

    const double y = sin(az) * sin_d * cos_lat1;
    const double x = cos_d - sin_lat1 * sin(lat2);
    const double lon2 = lon1 + atan2(y, x);

    *lat2_deg = lat2 * RAD2DEG;
    *lon2_deg = wrap_lon_0_360(lon2 * RAD2DEG);
}

static double apparent_elevation_deg(const double rx_xyz[3],
                                     const double east[3],
                                     const double north[3],
                                     const double up[3],
                                     const double pt_xyz[3])
{
    double v[3];
    v[0] = pt_xyz[0] - rx_xyz[0];
    v[1] = pt_xyz[1] - rx_xyz[1];
    v[2] = pt_xyz[2] - rx_xyz[2];

    const double ve = dot3(v, east);
    const double vn = dot3(v, north);
    const double vu = dot3(v, up);
    const double vh = hypot(ve, vn);

    return atan2(vu, vh) * RAD2DEG;
}

static void usage(const char *prog)
{
    fprintf(stderr,
        "Usage:\n"
        "  %s <LDEM_80S_80MPP_ADJ.TIF> <output.csv> [options]\n\n"
        "Options:\n"
        "  --rx-lat deg       Receiver latitude [default %.8f]\n"
        "  --rx-lon deg       Receiver longitude, east-positive [default %.8f]\n"
        "  --h-ant m          Antenna height above DEM surface [default %.3f]\n"
        "  --az-step deg      Azimuth interval [default %.3f]\n"
        "  --r-min m          Minimum search range [default %.1f]\n"
        "  --r-max m          Maximum search range [default %.1f]\n"
        "  --r-step m         Radial search step [default %.1f]\n"
        "  --allow-negative   Do not clamp negative terrain horizons to 0 deg\n"
        "  --verbose          Add diagnostic columns to output CSV\n\n"
        "Example:\n"
        "  %s LDEM_80S_80MPP_ADJ.TIF shackleton_edge_mask.csv --verbose\n",
        prog,
        RX_LAT_DEG_DEFAULT, RX_LON_DEG_DEFAULT, H_ANT_M_DEFAULT,
        AZ_STEP_DEG_DEFAULT, R_MIN_M_DEFAULT, R_MAX_M_DEFAULT,
        R_STEP_M_DEFAULT,
        prog);
}

int main(int argc, char **argv)
{
    if (argc < 3) {
        usage(argv[0]);
        return 1;
    }

    const char *dem_path = argv[1];
    const char *out_path = argv[2];

    double rx_lat_deg = RX_LAT_DEG_DEFAULT;
    double rx_lon_deg = RX_LON_DEG_DEFAULT;
    double h_ant_m = H_ANT_M_DEFAULT;
    double az_step_deg = AZ_STEP_DEG_DEFAULT;
    double r_min_m = R_MIN_M_DEFAULT;
    double r_max_m = R_MAX_M_DEFAULT;
    double r_step_m = R_STEP_M_DEFAULT;
    int verbose = 0;
    int allow_negative = 0;

    for (int i = 3; i < argc; i++) {
        if (strcmp(argv[i], "--rx-lat") == 0 && i + 1 < argc) {
            rx_lat_deg = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--rx-lon") == 0 && i + 1 < argc) {
            rx_lon_deg = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--h-ant") == 0 && i + 1 < argc) {
            h_ant_m = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--az-step") == 0 && i + 1 < argc) {
            az_step_deg = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--r-min") == 0 && i + 1 < argc) {
            r_min_m = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--r-max") == 0 && i + 1 < argc) {
            r_max_m = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--r-step") == 0 && i + 1 < argc) {
            r_step_m = atof(argv[++i]);
        }
        else if (strcmp(argv[i], "--allow-negative") == 0) {
            allow_negative = 1;
        }
        else if (strcmp(argv[i], "--verbose") == 0) {
            verbose = 1;
        }
        else {
            fprintf(stderr, "Unknown or incomplete option: %s\n", argv[i]);
            usage(argv[0]);
            return 1;
        }
    }

    if (az_step_deg <= 0.0 || az_step_deg > 360.0 ||
        r_min_m < 0.0 || r_max_m <= r_min_m || r_step_m <= 0.0) {
        fprintf(stderr, "Error: invalid search parameters.\n");
        return 1;
    }

    LDEMGrid dem;
    if (!ldem_load(dem_path, &dem)) {
        fprintf(stderr, "Error: failed to load DEM: %s\n", dem_path);
        return 1;
    }

    double rx_dem_h_m = 0.0;
    if (!ldem_height_latlon(&dem, rx_lat_deg, rx_lon_deg, &rx_dem_h_m)) {
        fprintf(stderr,
                "Error: failed to read receiver DEM height at lat=%.8f lon=%.8f\n",
                rx_lat_deg, rx_lon_deg);
        ldem_free(&dem);
        return 1;
    }

    const double rx_h_m = rx_dem_h_m + h_ant_m;

    double rx_xyz[3], east[3], north[3], up[3];
    llh_to_moon_fixed(rx_lat_deg, rx_lon_deg, rx_h_m, rx_xyz);
    local_enu_basis(rx_lat_deg, rx_lon_deg, east, north, up);

    FILE *fp = fopen(out_path, "w");
    if (fp == NULL) {
        fprintf(stderr, "Error: failed to open output CSV: %s\n", out_path);
        ldem_free(&dem);
        return 1;
    }

    if (verbose) {
        fprintf(fp,
                "azimuth_deg,elevation_mask_deg,edge_range_m,edge_lat_deg,"
                "edge_lon_deg,edge_dem_height_m,rx_dem_height_m\n");
    }
    else {
        fprintf(fp, "azimuth_deg,elevation_mask_deg\n");
    }

    fprintf(stderr,
            "Receiver: lat=%.8f deg, lon=%.8f deg, DEM height=%.3f m, h_ant=%.3f m\n",
            rx_lat_deg, rx_lon_deg, rx_dem_h_m, h_ant_m);
    fprintf(stderr,
            "Search: az_step=%.3f deg, range=%.1f..%.1f m, dr=%.1f m\n",
            az_step_deg, r_min_m, r_max_m, r_step_m);

    for (double az_deg = 0.0; az_deg < 360.0 - 1e-12; az_deg += az_step_deg) {
        double best_el_deg = -DBL_MAX;
        double best_r_m = NAN;
        double best_lat_deg = NAN;
        double best_lon_deg = NAN;
        double best_dem_h_m = NAN;

        for (double r_m = r_min_m; r_m <= r_max_m + 1e-12; r_m += r_step_m) {
            double lat_deg, lon_deg;
            double dem_h_m;
            double pt_xyz[3];

            destination_point_sphere(rx_lat_deg, rx_lon_deg, az_deg, r_m,
                                     &lat_deg, &lon_deg);

            if (!ldem_height_latlon(&dem, lat_deg, lon_deg, &dem_h_m)) {
                continue;
            }

            llh_to_moon_fixed(lat_deg, lon_deg, dem_h_m, pt_xyz);

            const double el_deg = apparent_elevation_deg(rx_xyz,
                                                         east, north, up,
                                                         pt_xyz);

            if (el_deg > best_el_deg) {
                best_el_deg = el_deg;
                best_r_m = r_m;
                best_lat_deg = lat_deg;
                best_lon_deg = lon_deg;
                best_dem_h_m = dem_h_m;
            }
        }

        if (!isfinite(best_r_m)) {
            fprintf(stderr, "Warning: no valid DEM sample for az=%.3f deg\n", az_deg);
            best_el_deg = 0.0;
        }

        if (!allow_negative && best_el_deg < 0.0) {
            best_el_deg = 0.0;
        }

        if (verbose) {
            fprintf(fp, "%.6f,%.6f,%.3f,%.9f,%.9f,%.3f,%.3f\n",
                    az_deg, best_el_deg, best_r_m, best_lat_deg,
                    best_lon_deg, best_dem_h_m, rx_dem_h_m);
        }
        else {
            fprintf(fp, "%.6f,%.6f\n", az_deg, best_el_deg);
        }
    }

    fclose(fp);
    ldem_free(&dem);

    fprintf(stderr, "Wrote: %s\n", out_path);
    return 0;
}
