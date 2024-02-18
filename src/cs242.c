/**
 * @file cs242.c
 * @brief Main file for cs242 model
 * @version 1.0
 *
 * @section DESCRIPTION
 *
 * Delivers the cs242 model which consists base on the mesh created with
 * ucvm2mesh using model tiling : sfcvm,cs242,cvmsi(background)
 *
 */

#include "limits.h"
#include "ucvm_model_dtypes.h"
#include "cs242.h"
#include <assert.h>

int cs242_debug=0;

/** The config of the model */
char *cs242_config_string=NULL;
int cs242_config_sz=0;

// Constants
/** The version of the model. */
const char *cs242_version_string = "cs242";

// Variables
/** Set to 1 when the model is ready for query. */
int cs242_is_initialized = 0;

char cs242_data_directory[128];

/** Configuration parameters. */
cs242_configuration_t *cs242_configuration;
/** Holds pointers to the velocity model data OR indicates it can be read from file. */
cs242_model_t *cs242_velocity_model;

/** Proj coordinate transformation objects. can go from geo <-> utm */
PJ *cs242_geo2utm = NULL;

/** The cosine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double cs242_cos_rotation_angle = 0;
/** The sine of the rotation angle used to rotate the box and point around the bottom-left corner. */
double cs242_sin_rotation_angle = 0;

/** The height of this model's region, in meters. */
double cs242_total_height_m = 0;
/** The width of this model's region, in meters. */
double cs242_total_width_m = 0;
/**
 * Initializes the cs242 plugin model within the UCVM framework. In order to initialize
 * the model, we must provide the UCVM install path and optionally a place in memory
 * where the model already exists.
 *
 * @param dir The directory in which UCVM has been installed.
 * @param label A unique identifier for the velocity model.
 * @return Success or failure, if initialization was successful.
 */
int cs242_init(const char *dir, const char *label) {
    int tempVal = 0;
    char configbuf[512];
    double north_height_m = 0, east_width_m = 0, rotation_angle = 0;

    // Initialize variables.
    cs242_configuration = calloc(1, sizeof(cs242_configuration_t));
    cs242_velocity_model = calloc(1, sizeof(cs242_model_t));

    cs242_config_string = calloc(CS242_CONFIG_MAX, sizeof(char));
    cs242_config_string[0]='\0';
    cs242_config_sz=0;

    // Configuration file location.
    sprintf(configbuf, "%s/model/%s/data/config", dir, label);

    // Read the cs242_configuration file.
    if (cs242_read_configuration(configbuf, cs242_configuration) != SUCCESS)
        return FAIL;

    // Set up the data directory.
    sprintf(cs242_data_directory, "%s/model/%s/data/%s/", dir, label, cs242_configuration->model_dir);

    // Can we allocate the model, or parts of it, to memory. If so, we do.
    tempVal = cs242_try_reading_model(cs242_velocity_model);

    if (tempVal == SUCCESS) {
        fprintf(stderr, "WARNING: Could not load model into memory. Reading the model from the\n");
        fprintf(stderr, "hard disk may result in slow performance.\n");
    } else if (tempVal == FAIL) {
        cs242_print_error("No model file was found to read from.");
        return FAIL;
    }

    // We need to convert the point from lat, lon to UTM, let's set it up.
    char cs242_projstr[64];
    snprintf(cs242_projstr, 64, "+proj=utm +zone=%d +datum=NAD27 +units=m +no_defs", cs242_configuration->utm_zone);
    if (!(cs242_geo2utm = proj_create_crs_to_crs(PJ_DEFAULT_CTX, "EPSG:4326", cs242_projstr, NULL))) {
        cs242_print_error("Could not set up Proj transformation from EPSG:4326 to UTM.");
        cs242_print_error((char  *)proj_context_errno_string(PJ_DEFAULT_CTX, proj_context_errno(PJ_DEFAULT_CTX)));
        return (UCVM_CODE_ERROR);
    }


    // In order to simplify our calculations in the query, we want to rotate the box so that the bottom-left
    // corner is at (0m,0m). Our box's height is total_height_m and total_width_m. We then rotate the
    // point so that is is somewhere between (0,0) and (total_width_m, total_height_m). How far along
    // the X and Y axis determines which grid points we use for the interpolation routine.

    // Calculate the rotation angle of the box.
    north_height_m = cs242_configuration->top_left_corner_n - cs242_configuration->bottom_left_corner_n;
    east_width_m = cs242_configuration->top_left_corner_e - cs242_configuration->bottom_left_corner_e;

    // Rotation angle. Cos, sin, and tan are expensive computationally, so calculate once.
    rotation_angle = atan(east_width_m / north_height_m);

    cs242_cos_rotation_angle = cos(rotation_angle);
    cs242_sin_rotation_angle = sin(rotation_angle);

    cs242_total_height_m = sqrt(pow(cs242_configuration->top_left_corner_n - cs242_configuration->bottom_left_corner_n, 2.0f) +
          pow(cs242_configuration->top_left_corner_e - cs242_configuration->bottom_left_corner_e, 2.0f));
    cs242_total_width_m  = sqrt(pow(cs242_configuration->top_right_corner_n - cs242_configuration->top_left_corner_n, 2.0f) +
          pow(cs242_configuration->top_right_corner_e - cs242_configuration->top_left_corner_e, 2.0f));

    if(cs242_debug) {
      fprintf(stderr,"height %lf width %lf\n", north_height_m, east_width_m);
      fprintf(stderr,"totol height %lf total width %lf\n", cs242_total_height_m, cs242_total_width_m);
      fprintf(stderr,"cos angle %lf sin angle %lf\n", cs242_cos_rotation_angle, cs242_sin_rotation_angle);
    }

    // setup config_string 
    sprintf(cs242_config_string,"config = %s\n",configbuf);
    cs242_config_sz=1;


    // Let everyone know that we are initialized and ready for business.
    cs242_is_initialized = 1;

    return SUCCESS;
}

/*
#define PROJ_GEO_IPJ "+proj=latlong +datum=WGS84"
#define PROJ_GEO_OPJ "+proj=utm +zone=11 +ellps=WGS84"
static int to_utm(double *lon, double *lat) {
    projPJ ipj= pj_init_plus(PROJ_GEO_IPJ);
    projPJ opj = pj_init_plus(PROJ_GEO_OPJ);
    int p = pj_transform(ipj, opj, 1, 1, lon, lat, NULL );
    return p;
}
*/

/*** transform to UTM zone 32, then back to geographical 
         https://proj.org/en/9.3/development/quickstart.html
14    P = proj_create_crs_to_crs(
15        C, "EPSG:4326", "+proj=utm +zone=32 +datum=WGS84", 
16        NULL);

41    b = proj_trans(P, PJ_FWD, a);
42    printf("easting: %.3f, northing: %.3f\n", b.enu.e, b.enu.n);
43
44    b = proj_trans(P, PJ_INV, b);
45    printf("longitude: %g, latitude: %g\n", b.lp.lam, b.lp.phi);

and
    PJ_COORD c_in;
33    c_in.lpzt.z = 0.0;
34    c_in.lpzt.t = HUGE_VAL; // important only for time-dependent projections
35    c_in.lp.lam = lon;
36    c_in.lp.phi = lat;
37

and

typedef union {
    double v[4];
    PJ_XYZT xyzt;
    PJ_UVWT uvwt;
    PJ_LPZT lpzt;
    PJ_GEOD geod;
    PJ_OPK opk;
    PJ_ENU enu;
    PJ_XYZ  xyz;
    PJ_UVW  uvw;
    PJ_LPZ  lpz;
    PJ_XY   xy;
    PJ_UV   uv;
    PJ_LP   lp;
} PJ_COORD ;

examples:
https://proj.org/en/5.0/development/migration.html

*/

static int to_utm(double lon, double lat, double *point_u, double *point_v) {
    PJ_COORD xyzSrc = proj_coord(lat, lon, 0.0, HUGE_VAL);
    PJ_COORD xyzDest = proj_trans(cs242_geo2utm, PJ_FWD, xyzSrc);
    int err = proj_context_errno(PJ_DEFAULT_CTX);
    if (err) {
       fprintf(stderr, "Error occurred while transforming latitude=%.4f, longitude=%.4f to UTM.\n",
              lat, lon);
        fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
        return UCVM_CODE_ERROR;
    }
    *point_u = xyzDest.xyzt.x;
    *point_v = xyzDest.xyzt.y;
    return err;
}

static int to_geo(double point_u, double point_v, double *lon, double *lat) {
    PJ_COORD xyzSrc;
    xyzSrc.xyzt.x=point_u;
    xyzSrc.xyzt.y=point_v;
    PJ_COORD xyzDest = proj_trans(cs242_geo2utm, PJ_INV, xyzSrc);
    
    int err = proj_context_errno(PJ_DEFAULT_CTX);
    if (err) {
       fprintf(stderr, "Error occurred while transforming u=%.4f, v=%.4f to Geo.\n",
              point_u, point_v);
        fprintf(stderr, "Proj error: %s\n", proj_context_errno_string(PJ_DEFAULT_CTX, err));
        return UCVM_CODE_ERROR;
    }
    *lon=xyzDest.lp.lam;
    *lat=xyzDest.lp.phi;
    return err;
}


/**
 * Queries cs242 at the given points and returns the data that it finds.
 *
 * @param points The points at which the queries will be made.
 * @param data The data that will be returned (Vp, Vs, rho, Qs, and/or Qp).
 * @param numpoints The total number of points to query.
 * @return SUCCESS or FAIL.
 */
int cs242_query(cs242_point_t *points, cs242_properties_t *data, int numpoints) {
    int i = 0;

    double point_u = 0, point_v = 0;
    double point_x = 0, point_y = 0; 
				   
    int load_x_coord = 0, load_y_coord = 0, load_z_coord = 0;
    double x_percent = 0, y_percent = 0, z_percent = 0;

    cs242_properties_t surrounding_points[8];
    int zone = cs242_configuration->utm_zone;

    for (i = 0; i < numpoints; i++) {

        // We need to be below the surface to service this query.
        if (points[i].depth < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

	// lon,lat,u,v			     
	to_utm(points[i].longitude, points[i].latitude, &point_u, &point_v);

//if(cs242_debug) { fprintf(stderr,"lon %lf lat %lf\n", points[i].longitude, points[i].latitude); }
//if(cs242_debug) { fprintf(stderr,"point_u %lf point_v %lf\n", point_u, point_v); }

        // Point within rectangle.
        point_u -= cs242_configuration->bottom_left_corner_e;
        point_v -= cs242_configuration->bottom_left_corner_n;

        // We need to rotate that point, the number of degrees we calculated above.
        point_x = cs242_cos_rotation_angle * point_u - cs242_sin_rotation_angle * point_v;
        point_y = cs242_sin_rotation_angle * point_u + cs242_cos_rotation_angle * point_v;

        // Which point base point does that correspond to?
        load_x_coord = floor(point_x / cs242_total_width_m * (cs242_configuration->nx - 1));
        load_y_coord = floor(point_y / cs242_total_height_m * (cs242_configuration->ny - 1));

        // And on the Z-axis?
        load_z_coord = (cs242_configuration->depth / cs242_configuration->depth_interval - 1) -
                       floor(points[i].depth / cs242_configuration->depth_interval);

//if(cs242_debug) { fprintf(stderr,"load_x_coord %d load_y_coord %d load_z_coord %d\n", load_x_coord,load_y_coord,load_z_coord); }

        // Are we outside the model's X and Y boundaries?
        if (load_x_coord > cs242_configuration->nx - 2 || load_y_coord > cs242_configuration->ny - 2 || load_x_coord < 0 || load_y_coord < 0) {
            data[i].vp = -1;
            data[i].vs = -1;
            data[i].rho = -1;
            data[i].qp = -1;
            data[i].qs = -1;
            continue;
        }

        if(cs242_configuration->interpolation) {

          // Get the X, Y, and Z percentages for the bilinear or trilinear interpolation below.
          double x_interval=(cs242_configuration->nx > 1) ?
                     cs242_total_width_m / (cs242_configuration->nx-1):cs242_total_width_m;
          double y_interval=(cs242_configuration->ny > 1) ?
                     cs242_total_height_m / (cs242_configuration->ny-1):cs242_total_height_m;

          x_percent = fmod(point_u, x_interval) / x_interval;
          y_percent = fmod(point_v, y_interval) / y_interval;
          z_percent = fmod(points[i].depth, cs242_configuration->depth_interval) / cs242_configuration->depth_interval;

          if (load_z_coord < 1) {
              // We're below the model boundaries. Bilinearly interpolate the bottom plane and use that value.
              data[i].vp = -1;
              data[i].vs = -1;
              data[i].rho = -1;
              data[i].qp = -1;
              data[i].qs = -1;
              continue;
          } else {
              // Read all the surrounding point properties.
              cs242_read_properties(load_x_coord, load_y_coord, load_z_coord, &(surrounding_points[0]));    // Orgin.
              cs242_read_properties(load_x_coord + 1, load_y_coord, load_z_coord, &(surrounding_points[1]));    // Orgin + 1x
              cs242_read_properties(load_x_coord, load_y_coord + 1, load_z_coord, &(surrounding_points[2]));    // Orgin + 1y
              cs242_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord, &(surrounding_points[3]));    // Orgin + x + y, forms top plane.
              cs242_read_properties(load_x_coord, load_y_coord, load_z_coord - 1, &(surrounding_points[4]));    // Bottom plane origin
              cs242_read_properties(load_x_coord + 1, load_y_coord, load_z_coord - 1, &(surrounding_points[5]));    // +1x
              cs242_read_properties(load_x_coord, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[6]));    // +1y
              cs242_read_properties(load_x_coord + 1, load_y_coord + 1, load_z_coord - 1, &(surrounding_points[7]));    // +x +y, forms bottom plane.
  
              cs242_trilinear_interpolation(x_percent, y_percent, z_percent, surrounding_points, &(data[i]));
          }
          } else {
              cs242_read_properties(load_x_coord, load_y_coord, load_z_coord, &(data[i]));    // Orgin.
        }

        // Calculate Qp and Qs.
        if (data[i].vs < 1500) 
            data[i].qs = data[i].vs * 0.02;
        else
            data[i].qs = data[i].vs * 0.10;

        data[i].qp = data[i].qs * 1.5;
    }

    return SUCCESS;
}

/**
 * Retrieves the material properties (whatever is available) for the given data point, expressed
 * in x, y, and z co-ordinates.
 *
 * @param x The x coordinate of the data point.
 * @param y The y coordinate of the data point.
 * @param z The z coordinate of the data point.
 * @param data The properties struct to which the material properties will be written.
 */
void cs242_read_properties(int x, int y, int z, cs242_properties_t *data) {
    // Set everything to -1 to indicate not found.
    data->vp = -1;
    data->vs = -1;
    data->rho = -1;
    data->qp = -1;
    data->qs = -1;

    float *ptr = NULL;
    FILE *fp = NULL;
        long location = 0;

        // the z is inverted at line #145
        if ( strcmp(cs242_configuration->seek_axis, "fast-y") == 0 ||
                 strcmp(cs242_configuration->seek_axis, "fast-Y") == 0 ) { // fast-y,  cs242 
            if(strcmp(cs242_configuration->seek_direction, "bottom-up") == 0) { 
                location = ((long) z * cs242_configuration->nx * cs242_configuration->ny) + (x * cs242_configuration->ny) + y;
//if(cs242_debug) {fprintf(stderr,"LOCATION==%d(fast-y, bottom-up)\n", location); }
                } else { // nz starts from 0 up to nz-1
                    location = ((long)((cs242_configuration->nz -1) - z) * cs242_configuration->nx * cs242_configuration->ny) + (x * cs242_configuration->ny) + y;
//if(cs242_debug) {fprintf(stderr,"LOCATION==%d(fast-y, not bottom-up)\n", location); }
            }
        } else {  // fast-X, cca data
            if ( strcmp(cs242_configuration->seek_axis, "fast-x") == 0 ||
                     strcmp(cs242_configuration->seek_axis, "fast-X") == 0 ) { // fast-y,  cs242 
                if(strcmp(cs242_configuration->seek_direction, "bottom-up") == 0) { 
                    location = ((long)z * cs242_configuration->nx * cs242_configuration->ny) + (y * cs242_configuration->nx) + x;
//if(cs242_debug) {fprintf(stderr,"LOCATION==%d(fast-x, bottom-up)\n", location); }
                    } else { // bottom-up
                        location = ((long)(cs242_configuration->nz - z) * cs242_configuration->nx * cs242_configuration->ny) + (y * cs242_configuration->nx) + x;
//if(cs242_debug) {fprintf(stderr,"LOCATION==%d(fast-x, not bottom-up)\n", location); }
                }
            }
        }

    // Check our loaded components of the model.
    if (cs242_velocity_model->vs_status == 2) {
        // Read from memory.
        ptr = (float *)cs242_velocity_model->vs;
        data->vs = ptr[location];
    } else if (cs242_velocity_model->vs_status == 1) {
        // Read from file.
        fp = (FILE *)cs242_velocity_model->vs;
        fseek(fp, location * sizeof(float), SEEK_SET);
                float temp;
        fread(&(temp), sizeof(float), 1, fp);
                data->vs = temp;
    }

    // Check our loaded components of the model.
    if (cs242_velocity_model->vp_status == 2) {
        // Read from memory.
        ptr = (float *)cs242_velocity_model->vp;
        data->vp = ptr[location];
    } else if (cs242_velocity_model->vp_status == 1) {
        // Read from file.
        fp = (FILE *)cs242_velocity_model->vp;
        fseek(fp, location * sizeof(float), SEEK_SET);
                float temp;
        fread(&(temp), sizeof(float), 1, fp);
                data->vp=temp;
    }

    // Check our loaded components of the model.
    if (cs242_velocity_model->rho_status == 2) {
        // Read from memory.
        ptr = (float *)cs242_velocity_model->rho;
        data->rho = ptr[location];
    } else if (cs242_velocity_model->rho_status == 1) {
        // Read from file.
        fp = (FILE *)cs242_velocity_model->rho;
        fseek(fp, location * sizeof(float), SEEK_SET);
                float temp;
        fread(&(temp), sizeof(float), 1, fp);
                data->rho=temp;
    }
}

/**
 * Trilinearly interpolates given a x percentage, y percentage, z percentage and a cube of
 * data properties in top origin format (top plane first, bottom plane second).
 *
 * @param x_percent X percentage
 * @param y_percent Y percentage
 * @param z_percent Z percentage
 * @param eight_points Eight surrounding data properties
 * @param ret_properties Returned data properties
 */
void cs242_trilinear_interpolation(double x_percent, double y_percent, double z_percent,
                             cs242_properties_t *eight_points, cs242_properties_t *ret_properties) {
    cs242_properties_t *temp_array = calloc(2, sizeof(cs242_properties_t));
    cs242_properties_t *four_points = eight_points;

    cs242_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[0]);

    // Now advance the pointer four "cvms5_properties_t" spaces.
    four_points += 4;

    // Another interpolation.
    cs242_bilinear_interpolation(x_percent, y_percent, four_points, &temp_array[1]);

    // Now linearly interpolate between the two.
    cs242_linear_interpolation(z_percent, &temp_array[0], &temp_array[1], ret_properties);

    free(temp_array);
}

/**
 * Bilinearly interpolates given a x percentage, y percentage, and a plane of data properties in
 * origin, bottom-right, top-left, top-right format.
 *
 * @param x_percent X percentage.
 * @param y_percent Y percentage.
 * @param four_points Data property plane.
 * @param ret_properties Returned data properties.
 */
void cs242_bilinear_interpolation(double x_percent, double y_percent, cs242_properties_t *four_points, cs242_properties_t *ret_properties) {
    cs242_properties_t *temp_array = calloc(2, sizeof(cs242_properties_t));
    cs242_linear_interpolation(x_percent, &four_points[0], &four_points[1], &temp_array[0]);
    cs242_linear_interpolation(x_percent, &four_points[2], &four_points[3], &temp_array[1]);
    cs242_linear_interpolation(y_percent, &temp_array[0], &temp_array[1], ret_properties);
    free(temp_array);
}

/**
 * Linearly interpolates given a percentage from x0 to x1, a data point at x0, and a data point at x1.
 *
 * @param percent Percent of the way from x0 to x1 (from 0 to 1 interval).
 * @param x0 Data point at x0.
 * @param x1 Data point at x1.
 * @param ret_properties Resulting data properties.
 */
void cs242_linear_interpolation(double percent, cs242_properties_t *x0, cs242_properties_t *x1, cs242_properties_t *ret_properties) {
    ret_properties->vp  = (1 - percent) * x0->vp  + percent * x1->vp;
    ret_properties->vs  = (1 - percent) * x0->vs  + percent * x1->vs;
    ret_properties->rho = (1 - percent) * x0->rho + percent * x1->rho;
    ret_properties->qp  = (1 - percent) * x0->qp  + percent * x1->qp;
    ret_properties->qs  = (1 - percent) * x0->qs  + percent * x1->qs;
}

/**
 * Called when the model is being discarded. Free all variables.
 *
 * @return SUCCESS
 */
int cs242_finalize() {

    proj_destroy(cs242_geo2utm);
    cs242_geo2utm = NULL;

    if (cs242_velocity_model) free(cs242_velocity_model);
    if (cs242_configuration) free(cs242_configuration);

    return SUCCESS;
}

/**
 * Returns the version information.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int cs242_version(char *ver, int len)
{
  int verlen;
  verlen = strlen(cs242_version_string);
  if (verlen > len - 1) {
    verlen = len - 1;
  }
  memset(ver, 0, len);
  strncpy(ver, cs242_version_string, verlen);
  return 0;
}

/**
 * Returns the model config information.
 *
 * @param key Config key string to return.
 * @param sz number of config terms.
 * @return Zero
 */
int cs242_config(char **config, int *sz)
{
  int len=strlen(cs242_config_string);
  if(len > 0) {
    *config=cs242_config_string;
    *sz=cs242_config_sz;
    return SUCCESS;
  }
  return FAIL;
}


/**
 * Reads the cs242_configuration file describing the various properties of CVM-S5 and populates
 * the cs242_configuration struct. This assumes cs242_configuration has been "calloc'ed" and validates
 * that each value is not zero at the end.
 *
 * @param file The cs242_configuration file location on disk to read.
 * @param config The cs242_configuration struct to which the data should be written.
 * @return Success or failure, depending on if file was read successfully.
 */
int cs242_read_configuration(char *file, cs242_configuration_t *config) {
    FILE *fp = fopen(file, "r");
    char key[40];
    char value[80];
    char line_holder[128];

    // If our file pointer is null, an error has occurred. Return fail.
    if (fp == NULL) {
        cs242_print_error("Could not open the cs242_configuration file.");
        return FAIL;
    }

    // Read the lines in the cs242_configuration file.
    while (fgets(line_holder, sizeof(line_holder), fp) != NULL) {
        if (line_holder[0] != '#' && line_holder[0] != ' ' && line_holder[0] != '\n') {
            sscanf(line_holder, "%s = %s", key, value);

            // Which variable are we editing?
            if (strcmp(key, "utm_zone") == 0) config->utm_zone = atoi(value);
            if (strcmp(key, "model_dir") == 0) sprintf(config->model_dir, "%s", value);
            if (strcmp(key, "nx") == 0) config->nx = atoi(value);
            if (strcmp(key, "ny") == 0) config->ny = atoi(value);
            if (strcmp(key, "nz") == 0) config->nz = atoi(value);
            if (strcmp(key, "depth") == 0) config->depth = atof(value);
            if (strcmp(key, "top_left_corner_e") == 0) config->top_left_corner_e = atof(value);
            if (strcmp(key, "top_left_corner_n") == 0) config->top_left_corner_n = atof(value);
            if (strcmp(key, "top_right_corner_e") == 0) config->top_right_corner_e = atof(value);
            if (strcmp(key, "top_right_corner_n") == 0) config->top_right_corner_n = atof(value);
            if (strcmp(key, "bottom_left_corner_e") == 0) config->bottom_left_corner_e = atof(value);
            if (strcmp(key, "bottom_left_corner_n") == 0) config->bottom_left_corner_n = atof(value);
            if (strcmp(key, "bottom_right_corner_e") == 0) config->bottom_right_corner_e = atof(value);
            if (strcmp(key, "bottom_right_corner_n") == 0) config->bottom_right_corner_n = atof(value);
            if (strcmp(key, "depth_interval") == 0) config->depth_interval = atof(value);
            if (strcmp(key, "seek_axis") == 0) sprintf(config->seek_axis, "%s", value);
            if (strcmp(key, "seek_direction") == 0) sprintf(config->seek_direction, "%s", value);
            if (strcmp(key, "interpolation") == 0) { 
                config->interpolation=0;
                if (strcmp(value,"on") == 0) config->interpolation=1;
            }
        }
    }

    // Have we set up all cs242_configuration parameters?
    if (config->utm_zone == 0 || config->nx == 0 || config->ny == 0 || config->nz == 0 || config->model_dir[0] == '\0' || 
        config->seek_direction[0] == '\0' || config->seek_axis[0] == '\0' ||
        config->top_left_corner_e == 0 || config->top_left_corner_n == 0 || config->top_right_corner_e == 0 ||
        config->top_right_corner_n == 0 || config->bottom_left_corner_e == 0 || config->bottom_left_corner_n == 0 ||
        config->bottom_right_corner_e == 0 || config->bottom_right_corner_n == 0 || config->depth == 0 ||
        config->depth_interval == 0) {
        cs242_print_error("One cs242_configuration parameter not specified. Please check your cs242_configuration file.");
        return FAIL;
    }

    fclose(fp);

    return SUCCESS;
}

/**
 * Prints the error string provided.
 *
 * @param err The error string to print out to stderr.
 */
void cs242_print_error(char *err) {
    fprintf(stderr, "An error has occurred while executing cs242. The error was: %s\n",err);
    fprintf(stderr, "\n\nPlease contact software@scec.org and describe both the error and a bit\n");
    fprintf(stderr, "about the computer you are running cs242 on (Linux, Mac, etc.).\n");
}

/**
 * Check if the data is too big to be loaded internally (exceed maximum
 * allowable by a INT variable)
 *
 */
static int too_big() {
        long max_size= (long) (cs242_configuration->nx) * cs242_configuration->ny * cs242_configuration->nz;
        long delta= max_size - INT_MAX;

    if( delta > 0) {
        return 1;
        } else {
        return 0;
        }
}

/**
 * Tries to read the model into memory.
 *
 * @param model The model parameter struct which will hold the pointers to the data either on disk or in memory.
 * @return 2 if all files are read to memory, SUCCESS if file is found but at least 1
 * is not in memory, FAIL if no file found.
 */
int cs242_try_reading_model(cs242_model_t *model) {
    double base_malloc = cs242_configuration->nx * cs242_configuration->ny * cs242_configuration->nz * sizeof(float);
    int file_count = 0;
    int all_read_to_memory =0;
    char current_file[128];
    FILE *fp;

    // Let's see what data we actually have.
    sprintf(current_file, "%s/vp.dat", cs242_data_directory);
    if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
            model->vp = malloc(base_malloc);
            if (model->vp != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->vp, 1, base_malloc, fp);
                        all_read_to_memory++;
            fclose(fp);
            model->vp_status = 2;
            } else {
            model->vp = fopen(current_file, "rb");
            model->vp_status = 1;
                    }
        } else {
            model->vp = fopen(current_file, "rb");
            model->vp_status = 1;
                }
        file_count++;
    }

    sprintf(current_file, "%s/vs.dat", cs242_data_directory);
    if (access(current_file, R_OK) == 0) {
                if( !too_big( )) { // only if fit
            model->vs = malloc(base_malloc);
            if (model->vs != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->vs, 1, base_malloc, fp);
                        all_read_to_memory++;
            fclose(fp);
            model->vs_status = 2;
            } else {
            model->vs = fopen(current_file, "rb");
            model->vs_status = 1;
            }
        } else {
            model->vs = fopen(current_file, "rb");
            model->vs_status = 1;
                }
        file_count++;
    }

    sprintf(current_file, "%s/rho.dat", cs242_data_directory);
    if (access(current_file, R_OK) == 0) {
                if(!too_big() ) { // only if fit
            model->rho = malloc(base_malloc);
            if (model->rho != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->rho, 1, base_malloc, fp);
                        all_read_to_memory++;
            fclose(fp);
            model->rho_status = 2;
            } else {
            model->rho = fopen(current_file, "rb");
            model->rho_status = 1;
            }
        } else {
            model->rho = fopen(current_file, "rb");
            model->rho_status = 1;
        }
        file_count++;
    }

    sprintf(current_file, "%s/qp.dat", cs242_data_directory);
    if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
            model->qp = malloc(base_malloc);
            if (model->qp != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->qp, 1, base_malloc, fp);
                        all_read_to_memory++;
            fclose(fp);
            model->qp_status = 2;
            } else {
            model->qp = fopen(current_file, "rb");
            model->qp_status = 1;
            }
        } else {
            model->qp = fopen(current_file, "rb");
            model->qp_status = 1;
        }
        file_count++;
    }

    sprintf(current_file, "%s/qs.dat", cs242_data_directory);
    if (access(current_file, R_OK) == 0) {
                if( !too_big() ) { // only if fit
            model->qs = malloc(base_malloc);
            if (model->qs != NULL) {
            // Read the model in.
            fp = fopen(current_file, "rb");
            fread(model->qs, 1, base_malloc, fp);
                        all_read_to_memory++;
            model->qs_status = 2;
            } else {
            model->qs = fopen(current_file, "rb");
            model->qs_status = 1;
            }
        } else {
            model->qs = fopen(current_file, "rb");
            model->qs_status = 1;
        }
        file_count++;
    }

    if (file_count == 0)
        return FAIL;
    else if (file_count > 0 && all_read_to_memory != file_count)
        return SUCCESS;
    else
        return 2;
}

// The following functions are for dynamic library mode. If we are compiling
// a static library, these functions must be disabled to avoid conflicts.
#ifdef DYNAMIC_LIBRARY

/**
 * Init function loaded and called by the UCVM library. Calls cs242_init.
 *
 * @param dir The directory in which UCVM is installed.
 * @return Success or failure.
 */
int model_init(const char *dir, const char *label) {
    return cs242_init(dir, label);
}

/**
 * Query function loaded and called by the UCVM library. Calls cs242_query.
 *
 * @param points The basic_point_t array containing the points.
 * @param data The basic_properties_t array containing the material properties returned.
 * @param numpoints The number of points in the array.
 * @return Success or fail.
 */
int model_query(cs242_point_t *points, cs242_properties_t *data, int numpoints) {
    return cs242_query(points, data, numpoints);
}

/**
 * Finalize function loaded and called by the UCVM library. Calls cs242_finalize.
 *
 * @return Success
 */
int model_finalize() {
    return cs242_finalize();
}

/**
 * Version function loaded and called by the UCVM library. Calls cs242_version.
 *
 * @param ver Version string to return.
 * @param len Maximum length of buffer.
 * @return Zero
 */
int model_version(char *ver, int len) {
    return cs242_version(ver, len);
}

/**
 * Version function loaded and called by the UCVM library. Calls cs242_config.
 *
 * @param config Config string to return.
 * @param sz number of config terms
 * @return Zero
 */
int model_config(char **config, int *sz) {
    return cs242_config(config, sz);
}


int (*get_model_init())(const char *, const char *) {
        return &cs242_init;
}
int (*get_model_query())(cs242_point_t *, cs242_properties_t *, int) {
         return &cs242_query;
}
int (*get_model_finalize())() {
         return &cs242_finalize;
}
int (*get_model_version())(char *, int) {
         return &cs242_version;
}
int (*get_model_config())(char **, int*) {
    return &cs242_config;
}



#endif
