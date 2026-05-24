#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

#include "ldem.h"

#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

int main(int argc, char **argv)
{
    LDEMGrid dem;
    double   h;
    double   x, y;
    LDEMVec3 n;

    /*
     * normal_step should be larger than one DEM pixel to avoid
     * excessive sensitivity to local pixel noise.
     * For the 80 m/pixel LDEM, 160 to 400 m is a reasonable initial value.
     */
    const double normal_step = 240.0;

    if (argc != 5) {
        fprintf(stderr, "Usage:\n");
        fprintf(stderr, "  %s LDEM.TIF latlon latitude_deg longitude_deg\n", argv[0]);
        fprintf(stderr, "  %s LDEM.TIF xy x_m y_m\n", argv[0]);
        fprintf(stderr, "\nExamples:\n");
        fprintf(stderr, "  %s LDEM_80S_80MPP_ADJ.TIF latlon -89.67 129.78\n", argv[0]);
        fprintf(stderr, "  %s LDEM_80S_80MPP_ADJ.TIF xy 7690.244 -6402.723\n", argv[0]);
        return 1;
    }

    const char *dem_path = argv[1];
    const char *mode     = argv[2];

    if (!ldem_load(dem_path, &dem)) {
        fprintf(stderr, "Error: failed to load LDEM file.\n");
        return 1;
    }

    if (strcmp(mode, "latlon") == 0) {
        double lat_deg = atof(argv[3]);
        double lon_deg = atof(argv[4]);

        if (!ldem_latlon_to_xy(lat_deg, lon_deg, &x, &y)) {
            fprintf(stderr, "Error: failed to convert latitude/longitude to X/Y.\n");
            ldem_free(&dem);
            return 1;
        }

        if (!ldem_height_latlon(&dem, lat_deg, lon_deg, &h)) {
            fprintf(stderr, "Error: failed to get DEM height.\n");
            ldem_free(&dem);
            return 1;
        }

        if (!ldem_surface_normal_latlon(&dem, lat_deg, lon_deg,
                                        normal_step, &n, NULL, NULL)) {
            fprintf(stderr, "Error: failed to compute surface normal.\n");
            ldem_free(&dem);
            return 1;
        }

        printf("Input mode      : latlon\n");
        printf("Latitude        : %.10f deg\n", lat_deg);
        printf("Longitude       : %.10f deg\n", lon_deg);
    }
    else if (strcmp(mode, "xy") == 0) {
        x = atof(argv[3]);
        y = atof(argv[4]);

        if (!ldem_height_xy(&dem, x, y, &h)) {
            fprintf(stderr, "Error: failed to get DEM height.\n");
            ldem_free(&dem);
            return 1;
        }

        if (!ldem_surface_normal_xy(&dem, x, y, normal_step, &n)) {
            fprintf(stderr, "Error: failed to compute surface normal.\n");
            ldem_free(&dem);
            return 1;
        }

        printf("Input mode      : xy\n");
    }
    else {
        fprintf(stderr, "Error: unknown mode: %s\n", mode);
        fprintf(stderr, "Use either 'latlon' or 'xy'.\n");
        ldem_free(&dem);
        return 1;
    }

    double slope_deg = ldem_slope_from_normal(&n) * 180.0 / M_PI;

    printf("Projected X     : %.3f m\n", x);
    printf("Projected Y     : %.3f m\n", y);
    printf("Surface height  : %.3f m\n", h);
    printf("Normal step     : %.3f m\n", normal_step);
    printf("Surface normal  : nx = %.10f, ny = %.10f, nz = %.10f\n",
           n.x, n.y, n.z);
    printf("Slope angle     : %.6f deg\n", slope_deg);

    ldem_free(&dem);

    return 0;
}
