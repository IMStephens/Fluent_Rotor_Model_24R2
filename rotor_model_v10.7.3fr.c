#include "udf.h"
#include "seem.h"
#include "global.h"
#include "mem.h"
#include "thread_mem.h"
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define NUM_UDM              15    /* Number of UDM locations to pre-allocate. */
#define max_rotor_zones      10    /* Maximum number of rotor zones. */
#define max_sections         20    /* Maximum number of blade sections. */
#define cell_volume_samples   3    /* Number of cells sampled for volume calc. */
#define pb_vert_angle_lim    (M_PI / 4.0) /* Vertical angle limit (rad). */
#define dbg                   0    /* Set to 1 to enable verbose debug output. */
#define my_mem_loop(m)        for ((m) = my_mem; (m) != NULL; (m) = (m)->next)
#define dimport               1    /* Import rotor data from file: 1=YES 0=NO. */
#define gauss_routine         0    /* Enable Gaussian optimisation routine. */


/*---------------------------------------------------------------------------*/
/* GUI helper: returns the n-th element of a Scheme list.                    */
/*---------------------------------------------------------------------------*/
#if !RP_NODE
static Pointer
List_Ref(Pointer l, int n)
{
    int i;
    for (i = 0; i < n; ++i)
        l = CDR(l);
    return CAR(l);
}
#endif


/*---------------------------------------------------------------------------*/
/* v1.0 does not allow for rotation in math. negative azimuthal direction    */
/* v2.0 allows for math. negative azimuthal rotation                         */
/* v3.0 URF for source terms, ramp up of rotor rpm                           */
/* v4.0 Thrust trimming in analog to lz                                      */
/* v5.0 Compute momentum coefficients, own trimming                          */
/* v5.1 Initialize xsource_old, ysource_old, and zsource_old with current    */
/*      source strength instead of 0                                         */
/* v6.0 GUI implementation. twist and type linear interpolation              */
/* v7.0 NON ideal gas for density allowed                                    */
/* v8.0 Torque calculation                                                   */
/* v9.0 Parallelized and Fluent 6.2 (boolean --> cxboolean)                  */
/*---------------------------------------------------------------------------*/
/*                 BEGIN GLOBAL INPUT VARIABLES                              */
/*---------------------------------------------------------------------------*/
/*Max. number of rotor zones 10*/
/*Max. number of blade sections 20*/
/*Max. name length for profile types 30 characters*/
/*Max. name length for airfoil tables 30 characters*/
/*Max. number of different airfoil tables 100*/
/*Max. number of AOA and coeff per cl or cd 80*/

/*GUI*/
int ac_mass = 7000;             /* Aircraft mass for force balancing, kg. */
int nrtz;                       /* Number of Rotor zones. */
int nbld[max_rotor_zones];      /* Number of blades. */
real rrl[max_rotor_zones];      /* Rotor radius, m. */
real rspe[max_rotor_zones];     /* Rotor speed, rad/s. */
real teff[max_rotor_zones];     /* Tip effect. */
real dskco[max_rotor_zones][3]; /* Rotor disk origin 0x 1y 2z, m. */
real dpit[max_rotor_zones];     /* Rotor disk pitch angle, rad. */
real tdpit[max_rotor_zones];    /* Rotor disk trim pitch angle, rad. */
real dban[max_rotor_zones];     /* Rotor disk bank angle, rad. */
real tdban[max_rotor_zones];    /* Rotor disk trim bank angle, rad. */
int fzon[max_rotor_zones];      /* Rotor face zone ID. */

real bcop[max_rotor_zones];     /* Blade pitch collective, rad. */
real bcys[max_rotor_zones];     /* Blade pitch cyclic sin*/
real bcyc[max_rotor_zones];     /* Blade pitch cyclic cos*/
real bflco[max_rotor_zones];    /* Blade flapping cone, rad. */
real bfls[max_rotor_zones];     /* Blade flapping cyclic sin*/
real bflc[max_rotor_zones];     /* Blade flapping cyclic cos*/

int nsec[max_rotor_zones];      /*Number of blade sections along span. */
real rin[max_rotor_zones][max_sections];   /* Normalized inner radius. */
real rout[max_rotor_zones][max_sections];  /* Normalized outer radius. */
real rsec[max_rotor_zones][max_sections];  /* Normalizes radius from GUI. */
real cin[max_rotor_zones][max_sections];   /* Chord at nd. inner radius. */
real cout[max_rotor_zones][max_sections];  /* Chord at nd. outer radius. */
real csec[max_rotor_zones][max_sections];  /* Chord at nd. sectional radius. */
real twst[max_rotor_zones][max_sections];  /* Twist*/
char type[max_rotor_zones][max_sections][30];  /*Profile type name */

int trmco[max_rotor_zones]={0,0,0,0,0,0,0,0,0,0};  /*Trimming collective pitch 1(ON) 0(OFF)*/
int trmcy[max_rotor_zones]={0,0,0,0,0,0,0,0,0,0};  /*Trimming cyclic pitch 1(ON) 0(OFF)*/
int forba[max_rotor_zones]={0,0,0,0,0,0,0,0,0,0};  /*Force balancing 1(ON) 0(OFF)*/
int momba[max_rotor_zones]={0,0,0,0,0,0,0,0,0,0};  /*Moment balancing 1(ON) 0(OFF)*/
int fuseid[max_rotor_zones];  /* Fuselage face thread ID. */
real cmgx[max_rotor_zones];   /* x-location of fuselage CG. */
real cmgy[max_rotor_zones];   /* y-location of fuselage CG. */
real cmgz[max_rotor_zones];   /* z-location of fuselage CG. */
int errchk;                   /* Error check for NaN values. */    
   
real trdf[max_rotor_zones];   /* Trimming damping factor. */
int trufq[max_rotor_zones];   /* Trimming update frequency. */
real ctd[max_rotor_zones];    /* Trimming desired ct. */
real cmxd[max_rotor_zones];   /* Trimming desired cmx. */
real cmyd[max_rotor_zones];   /* Trimming desired cmy. */
int trnw = 0;                 /* Currently trimming. */

real CT[max_rotor_zones];   /* Rotor thrust coeff. */
real CMX[max_rotor_zones];  /* Rotor moment coeff x. */
real CMY[max_rotor_zones];  /* Rotor moment coeff y. */

real sum_thrust[max_rotor_zones];  /* Summation of Rotor Thrust. */
real sum_torque[max_rotor_zones];  /* Summation of Rotor Torque. */
real power[max_rotor_zones];       /* Summation of Rotor Power. */

real force[ND_ND];                 /* Force on thread object. */
real moment[ND_ND];                /* Moment on thread object. */
real cg[ND_ND];                    /* CG coordinates of thread object. */
real main_rot_thrust_target_init;  /* Initial main rotor thrust coeff target. */

int ktot;                          /* Number of different airfoil tables read in. */
char file_name[100][30];           /* Airfoil table filenames, e.g. "airfoil.dat". */
char check_name[100][30];          /* Airfoil names, e.g. "naca0012". */
char clorcd[100][50][max_rotor_zones]; /* CL or CD identifier per table entry. */
float RE[100][50];                 /* Airfoil Reynolds numbers. */
float MA[100][50];                 /* Airfoil Mach numbers. */
int itot[100];                     /* Number of entries per .dat file. */
int jtot[100][50];                 /* Number of AOA entries per set. */
float aoa[100][50][80];            /* Angle of attack table values. */
float coeff[100][50][80];          /* CL or CD table values. */
cxboolean rho_const;               /* TRUE if density is constant, FALSE otherwise. */

/*TUI*/
real urf_source = 0.5;              /* Under-relaxation factor for source calculation. */
int  i_start    = 10;               /* Iterations for linear rotor speed ramp-up. */
real dalpha     = 10.0 * M_PI / 180.0; /* Trimming: Jacobian perturbation angle (rad). */
real limiter    = 0.1;              /* Max angle change relative to dalpha. */

char newfilename[256];              /* Filename for reading/writing rotor data. */
char time_str[100];                 /* Time string for filename tagging. */
int  hdr = 1;                       /* Output file header flag: 1=write header. */

int  n_r = 15;                       /* Number of ALM nodes. */
real sigma_az_old = 0.0;    /* Previous value of sigma_az. */
real sigma_vert_old = 0.0;  /* Previous value of sigma_vert. */
real gauss_denom_old = 0.0; /* Previous value of Gaussian denominator. */
int global_cells_ahead_behind_rotor_old[max_rotor_zones] = {0}; /* Previous value of cells ahead/behind rotor. */
int total_cells_above_rotor_old[max_rotor_zones] = {0}; /* Previous value of cells above/below rotor. */

real sigma_az_mplr; /* Following Gauss opt routine, calculated sigma_az multiplier. */
real sigma_vert_mplr; /* Following Gauss opt routine, calculated sigma_vert multiplier. */
real gauss_denom_mplr; /* Following Gauss opt routine, calculated Gauss denom multiplier. */

const char *data_fname[] = {
"coll_data.txt",
"cycc_data.txt",
"cycs_data.txt",
"thrust_data.txt",
"torque_data.txt",
"power_data.txt",
"thrust_coeff_data.txt",
"roll_mom_coeff_data.txt",
"pitch_mom_coeff_data.txt",
};

const char *data_vname[] = {
"Collective Pitch",
"Cyclic Pitch Cos",
"Cyclic Pitch Sine",
"Rotor Thrust",
"Rotor Torque",
"Rotor Power",
"Rotor Thrust Coefficient",
"Rotor Roll Moment Coefficient",
"Rotor Pitch Moment Coefficient",
};

/*---------------------------------------------------------------------------*/
/*                 END GLOBAL INPUT VARIABLES                                */
/*---------------------------------------------------------------------------*/
int istflag = 1; /* 1 at start; set to 0 after first ADJUST run. */
/*---------------------------------------------------------------------------*/
/*                 BEGIN GLOBAL TRIMMING VARIABLES                           */
/*---------------------------------------------------------------------------*/
int trufq_co[max_rotor_zones], trufq_cy[max_rotor_zones];
real trdf_co[max_rotor_zones], trdf_cy[max_rotor_zones];
int up_co[max_rotor_zones], up_cy[max_rotor_zones];
real sldity[max_rotor_zones];
/*---------------------------------------------------------------------------*/
/*                   END GLOBAL TRIMMING VARIABLES                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*                 BEGIN GLOBAL MEMORY                                       */
/*---------------------------------------------------------------------------*/
#if !RP_HOST
  int czon[max_rotor_zones];                      /*Rotor cell zone ID*/
  char factor[max_rotor_zones][30];               /*Memory: Geometry weighting factor*/
  char xsource[max_rotor_zones][30];        /*Memory: Source x-component*/
  char ysource[max_rotor_zones][30];        /*Memory: Source y-component*/
  char zsource[max_rotor_zones][30];        /*Memory: Source z-component*/

  char xsource_old[max_rotor_zones][30];  /*Memory: Source x-component previous timestep UR*/
  char ysource_old[max_rotor_zones][30];  /*Memory: Source y-component previous timestep UR*/
  char zsource_old[max_rotor_zones][30];  /*Memory: Source z-component previous timestep UR*/
#endif
/*---------------------------------------------------------------------------*/
/*                 END GLOBAL MEMORY                                         */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*                 PUT ALL THE FUNCTIONS HERE                                */
/*---------------------------------------------------------------------------*/

typedef struct thread_mem_struct 
{
  int thread_id;
  char mem_name[64];
  void *mem;
  struct thread_mem_struct *next;
} 
Thread_Mem_Struct;

Thread_Mem_Struct *my_mem = NULL;

void * Get_Thread_Memory(Thread *t, char *mem_name)
{
  Thread_Mem_Struct *m;
  int id = THREAD_ID(t);
  my_mem_loop(m)
  {
    if ((id == m->thread_id) && STREQ(mem_name, m->mem_name))
      return m->mem;
  }
  return NULL;
}

static void add_thread_memory(Thread *t, char *mem_name, void *mem)
{
  Thread_Mem_Struct *new_mem = malloc(sizeof(Thread_Mem_Struct));
  if (NULL == new_mem)
    Error("Unable to add new thread memory");

  strncpy(new_mem->mem_name,mem_name,64);
  new_mem->mem = mem;
  new_mem->thread_id = THREAD_ID(t);
  new_mem->next = my_mem;
  my_mem = new_mem;
}

void Alloc_Thread_Memory(Thread *t, char *mem_name, size_t size)
{
  void *new_mem;
  if (NNULLP(new_mem = Get_Thread_Memory(t,mem_name))) {
    Message("Warning: Memory Already Allocated");
    return;
  }
  new_mem = malloc(size);
  if (NULL == new_mem)
    Error("Unable to add new thread memory");

  add_thread_memory(t, mem_name, new_mem);
  return;
}

void Free_Thread_Memory(Thread *t, char *mem_name)
{
  Thread_Mem_Struct *m, *pre;
  pre = NULL;
  my_mem_loop(m) {
    if ((THREAD_ID(t) == m->thread_id) && STREQ(mem_name, m->mem_name)) {
      if (NNULLP(pre))
        pre->next = m->next;

      if (m == my_mem)
        my_mem = m->next;

      free(m->mem);
      free(m);
      break;
    }
    pre = m;
  }
}

#if !RP_HOST
/*---------------------------------------------------------------------------*/
  void my_get_solidity()
/*---------------------------------------------------------------------------*/
/* Computes the solidity for each rotor zone.                                */
/*---------------------------------------------------------------------------*/
  {
    int i, j;
    real area, area_tot;

    for (i = 0; i < nrtz; i++) {
      area_tot = 0.0;
      for (j = 0; j < nsec[i]; j++) {
        area = 0.5 * (cout[i][j] + cin[i][j]) * (rout[i][j] - rin[i][j]) * rrl[i];
        area_tot += area;
      }

      if ((M_PI * rrl[i] * rrl[i]) == 0.0) {
        Message("Error: Rotor %d radius is undefined.\n", i + 1);
        sldity[i] = 0.0;
      }
      else {
        sldity[i] = nbld[i] * area_tot / (M_PI * rrl[i] * rrl[i]);
      }
    }
  }

/*---------------------------------------------------------------------------*/
  void my_obtain_cell_id()
/*---------------------------------------------------------------------------*/
/* Identifies the fluid cell zone ID for each rotor zone.                    */
/*---------------------------------------------------------------------------*/
  {
    int i;
    int fl_zone;
    static int fluid_zone_id = -1;
    Thread *t;
    Domain *domain = Get_Domain(1);

    if (domain == NULL) {
      Message0("Error: Domain is NULL.\n");
      return;
    }

    if (dbg == 1)
      Message0("%d obtain_cell_id\n", myid);

    fl_zone = 0;
    for (i = 0; i < nrtz; i++) {
      if (fluid_zone_id < 0) {
        thread_loop_c(t, domain) {
          if (FLUID_THREAD_P(t)) {
            if (fl_zone == 1) {
              Message0("Warning: More than one fluid zone detected. "
                "Using first detected zone: %d\n", czon[i]);
              continue;
            }
            fl_zone = 1;
            czon[i] = THREAD_ID(t);
            Message0("Fluid zone detected: %d\n", czon[i]);
          }
        }
      }
    }
  }


/*---------------------------------------------------------------------------*/
  void my_allocate_memory(char factor[][30], char xsource[][30],
    char ysource[][30], char zsource[][30], char xsource_old[][30],
    char ysource_old[][30], char zsource_old[][30])
/*---------------------------------------------------------------------------*/
/* Allocates per-cell thread memory for force source arrays and initialises  */
/* the previous-timestep (old) source terms to zero.                         */
/*---------------------------------------------------------------------------*/
  {
    int i, count, sum_count, local_idx;
    Thread *tc = NULL;
    Domain *domain;
    cell_t cc;
    real *xsource_old_ptr, *ysource_old_ptr, *zsource_old_ptr;

    domain = Get_Domain(1);
    if (domain == NULL) {
      Message0("Error: Domain is NULL.\n");
      return;
    }

    if (dbg == 1)
      Message0("%d allocate_memory\n", myid);

    /* Build memory name strings for each rotor zone. */
    for (i = 0; i < max_rotor_zones; i++) {
      sprintf(factor[i],      "factor%d",      i);
      sprintf(xsource[i],     "xsource%d",     i);
      sprintf(ysource[i],     "ysource%d",     i);
      sprintf(zsource[i],     "zsource%d",     i);
      sprintf(xsource_old[i], "xsource_old%d", i);
      sprintf(ysource_old[i], "ysource_old%d", i);
      sprintf(zsource_old[i], "zsource_old%d", i);
    }

    /* Allocate per-cell arrays for each rotor zone. */
    for (i = 0; i < nrtz; i++) {
      tc = Lookup_Thread(domain, czon[i]);
      if (tc == NULL) {
        Message0("Error: tc is NULL in my_allocate_memory\n");
        return;
      }

      count = 0;
      begin_c_loop_int(cc, tc)
        count++;
      end_c_loop_int(cc, tc)

      sum_count = PRF_GISUM1(count);
      if (dbg == 1)
        Message0("  %d Rotor %d: cells in volume = %d (on this node = %d)\n",
          myid, i + 1, sum_count, count);

      Alloc_Thread_Memory(tc, factor[i],      sizeof(real) * count);
      Alloc_Thread_Memory(tc, xsource[i],     sizeof(real) * count);
      Alloc_Thread_Memory(tc, ysource[i],     sizeof(real) * count);
      Alloc_Thread_Memory(tc, zsource[i],     sizeof(real) * count);
      Alloc_Thread_Memory(tc, xsource_old[i], sizeof(real) * count);
      Alloc_Thread_Memory(tc, ysource_old[i], sizeof(real) * count);
      Alloc_Thread_Memory(tc, zsource_old[i], sizeof(real) * count);
    }

    /* Initialise old-timestep source arrays to zero. */
    for (i = 0; i < nrtz; i++) {
      tc = Lookup_Thread(domain, czon[i]);
      if (tc == NULL) {
        Message0("Error: tc is NULL in my_allocate_memory (init pass)\n");
        return;
      }

      xsource_old_ptr = Get_Thread_Memory(tc, xsource_old[i]);
      ysource_old_ptr = Get_Thread_Memory(tc, ysource_old[i]);
      zsource_old_ptr = Get_Thread_Memory(tc, zsource_old[i]);

      if (!xsource_old_ptr || !ysource_old_ptr || !zsource_old_ptr) {
        Message0("Error: NULL old-source pointer in my_allocate_memory\n");
        continue;
      }

      local_idx = 0;
      begin_c_loop_int(cc, tc)
        xsource_old_ptr[local_idx] = 0.0;
        ysource_old_ptr[local_idx] = 0.0;
        zsource_old_ptr[local_idx] = 0.0;
        local_idx++;
      end_c_loop_int(cc, tc)
    }
  }

/*---------------------------------------------------------------------------*/
  void my_geom_factor()
/*---------------------------------------------------------------------------*/
/* Computes the geometric weighting factor for each cell in each rotor zone. */
/*---------------------------------------------------------------------------*/
  {
    int i, local_indx;
    Thread *tc = NULL;
    cell_t cc;
    real xrc[ND_ND], rc;
    real *area_ptr;

    Domain *domain = Get_Domain(1);
    if (domain == NULL) {
      Message0("Error: Domain is NULL.\n");
      return;
    }

    if (dbg == 1)
      Message0("%d geom_factor\n", myid);

    for (i = 0; i < nrtz; i++) {
      tc = Lookup_Thread(domain, czon[i]);
      if (tc == NULL) {
        Message0("Error: tc is NULL in my_geom_factor (Rotor %d)\n", i + 1);
        continue;
      }

      area_ptr = Get_Thread_Memory(tc, factor[i]);
      if (area_ptr == NULL) {
        Message0("Error: area_ptr is NULL in my_geom_factor (Rotor %d)\n", i + 1);
        continue;
      }

      local_indx = 0;
      begin_c_loop_int(cc, tc)
        C_CENTROID(xrc, cc, tc);

        /* Distance from rotor disk origin. */
        rc = sqrt((xrc[0] - dskco[i][0]) * (xrc[0] - dskco[i][0]) +
                  (xrc[1] - dskco[i][1]) * (xrc[1] - dskco[i][1]) +
                  (xrc[2] - dskco[i][2]) * (xrc[2] - dskco[i][2]));

        if (rc <= 1e-12) {
          Message0("Error: Rotor %d origin coincides with cell centroid.\n", i + 1);
          area_ptr[local_indx] = 0.0;
        }
        else {
          area_ptr[local_indx] = ((real)nbld[i]) * ((rrl[i] * rrl[i]) /
                                  (n_r * n_r)) / (2.0 * M_PI * rc);
        }
        local_indx++;
      end_c_loop_int(cc, tc)
    }
  }


/*---------------------------------------------------------------------------*/
  void my_start_trimming(int trufq_co[], int trufq_cy[],
    real trdf_co[], real trdf_cy[],
    int up_co[], int up_cy[])
/*---------------------------------------------------------------------------*/
/* Initialises trimming update frequencies and damping factors.              */
/*---------------------------------------------------------------------------*/
  {
    int i;

    if (dbg == 1)
      Message0("%d start_trimming\n", myid);

    for (i = 0; i < nrtz; i++) {
      trufq_co[i] = trufq[i];
      trufq_cy[i] = trufq[i];
      trdf_co[i]  = trdf[i];
      trdf_cy[i]  = trdf[i];
      up_co[i]    = trufq_co[i];
      up_cy[i]    = trufq_cy[i];
    }
  }

#endif

#if !RP_NODE
/*---------------------------------------------------------------------------*/
  void my_read_in_airfoil_tables()
/*---------------------------------------------------------------------------*/
/* Reads CL/CD polars from .dat files for each airfoil type.                 */
/* Expects AOA entries ordered from -180 to +180 degrees.                    */
/*---------------------------------------------------------------------------*/
  {
    FILE *f_ptr;
    int i, j, k;

    if (dbg == 1) {
      Message("\n read_in_airfoil_tables\n");
      Message("     -180 < AOA < 180\n");
    }

    for (k = 0; k < ktot; k++) {
      f_ptr = fopen(file_name[k], "r");
      if (f_ptr == NULL) {
        Message("\nNull file pointer for %s\n", file_name[k]);
        continue;
      }

      fscanf(f_ptr, "%s\n", check_name[k]);
      fscanf(f_ptr, "%d\n", &itot[k]);

      for (i = 0; i < itot[k]; i++) {
        fscanf(f_ptr, "%s\n",  clorcd[k][i]);
        fscanf(f_ptr, "%f\n",  &RE[k][i]);
        fscanf(f_ptr, "%f\n",  &MA[k][i]);
        fscanf(f_ptr, "%d\n",  &jtot[k][i]);

        for (j = 0; j < jtot[k][i]; j++)
          fscanf(f_ptr, "%f %f\n", &aoa[k][i][j], &coeff[k][i][j]);

        /* Verify AOA entries are in ascending order. */
        for (j = 1; j < jtot[k][i]; j++) {
          if (aoa[k][i][j - 1] > aoa[k][i][j])
            Message("ERROR: Airfoil table entries must be ordered with "
              "increasing AOA\n");
        }

        /* Sanity-check CL/CD signs. */
        for (j = 0; j < jtot[k][i]; j++) {
          if (strcmp(clorcd[k][i], "cd") == 0) {
            if (coeff[k][i][j] < 0.0)
              Message("WARNING: CD is negative in %s at AOA %f?\n",
                check_name[k], aoa[k][i][j]);
          }
          if (strcmp(clorcd[k][i], "cl") == 0) {
            if (aoa[k][i][j] < -90.0 && coeff[k][i][j] < 0.0)
              Message("WARNING: CL is negative in %s at AOA %f?\n",
                check_name[k], aoa[k][i][j]);
            if (aoa[k][i][j] > 90.0 && coeff[k][i][j] > 0.0)
              Message("WARNING: CL is positive in %s at AOA %f?\n",
                check_name[k], aoa[k][i][j]);
            if (aoa[k][i][j] >= 0.0 && aoa[k][i][j] <= 90.0 && coeff[k][i][j] < 0.0)
              Message("WARNING: CL is negative in %s at AOA %f?\n",
                check_name[k], aoa[k][i][j]);
            if (aoa[k][i][j] >= -90.0 && aoa[k][i][j] < 0.0 && coeff[k][i][j] > 0.0)
              Message("WARNING: CL is positive in %s at AOA %f?\n",
                check_name[k], aoa[k][i][j]);
          }
        }
      }
      fclose(f_ptr);
    }
  }
  
#endif

#if !RP_HOST
/*---------------------------------------------------------------------------*/
  void my_get_radial_position(int i, cell_t cc, Thread *tc, int count,
    real *r_pb_ptr, real *psi_pb_ptr, real *r_pb2_ptr,
    real *x_pb_ptr, real *y_pb_ptr, real *z_pb_ptr, int *rotor_cell,
    real *r_total, real *xx, real *yy, real *zz, real *beta_pb_ptr,
    real *dpitc, real *dpits, real *dbanc, real *dbans)
/*---------------------------------------------------------------------------*/
/* Transforms the global cell centroid position into rotor-disc coordinates  */
/* and returns normalised radius, azimuth, vertical angle, and trig values.  */
/* Sets *rotor_cell = 0 if the cell is outside the rotor influence region.   */
/*---------------------------------------------------------------------------*/
  {
    real xrc[ND_ND];
    real x_pb, y_pb, z_pb;
    real tmp_radius, tmp_angle;
  
    if (count == 0 && dbg == 1)
        Message0("  %d Rotor %d - get_radial_position \n", myid, i+1);
  
    /* Distance from rotor centre in global x, y, and z. */
    C_CENTROID(xrc,cc,tc);
    *xx = xrc[0] - dskco[i][0];
    *yy = xrc[1] - dskco[i][1];
    *zz = xrc[2] - dskco[i][2];
    *r_total = sqrt((*xx) * (*xx) + (*yy) * (*yy) + (*zz) * (*zz));
      
    /* If cell is beyond rotor radius, move on to next cell. */
    if (*r_total > rrl[i]) {
      *rotor_cell = 0;
      return;
    }

    /* Pitch angle components. */
    *dpits = sin(dpit[i]);
    *dpitc = cos(dpit[i]);
  
    /* Bank angle components. */
    *dbans = sin(dban[i]);
    *dbanc = cos(dban[i]);

    /* Distance from rotor centre in rotor x, y, and z. */
    x_pb = (*xx) * (*dpitc) - (*zz) * (*dpits);
    y_pb = (*xx) * (*dpits) * (*dbans) + (*yy) * (*dbanc) + (*zz) * (*dpitc) * (*dbans);
    z_pb = (*xx) * (*dpits) * (*dbanc) - (*yy) * (*dbans) + (*zz) * (*dpitc) * (*dbanc);

    *x_pb_ptr = x_pb;
    *y_pb_ptr = y_pb;
    *z_pb_ptr = z_pb;
  
    /* If cell is 45 deg above rotor, move on to next cell. */
    if (fabs(z_pb) > pb_vert_angle_lim * rrl[i]) {
      *rotor_cell = 0;
      return;
    }

    /* Radius in rotor x and y co-ordinates. */
    *r_pb2_ptr = sqrt(x_pb * x_pb + y_pb * y_pb);

    /* Total Distance from rotor centre. */
    tmp_radius  = sqrt(x_pb * x_pb + y_pb * y_pb + z_pb * z_pb);
    if (rrl[i] <= 1e-12) {
      Message("Error: Rotor %d radius is undefined.\n", i+1);
      tmp_radius = 0.0;
    }
    else {
      /* Normalised radius. */
      tmp_radius /= rrl[i];
    }
    *r_pb_ptr = tmp_radius;
  
    /*Azimuthal angle psi in pitch and bank plane. */
    tmp_angle  = atan2(y_pb, x_pb);
    if(tmp_angle  < 0.0)
      tmp_angle  += 2.0 * M_PI;
    *psi_pb_ptr = tmp_angle;

    /* Vertical angle beta from rotor plane */
    if (*r_pb2_ptr <= 1e-12) {
        *beta_pb_ptr = 0.0;
    } else {
        *beta_pb_ptr = atan2(z_pb, *r_pb2_ptr);
    }

    *rotor_cell = 1; /* Mark as valid rotor cell */
  }


/*---------------------------------------------------------------------------*/
  void my_vel_xyz_to_blade(int i, cell_t cc, Thread *tc, int count,
    real psi_pb, real r_pb2,
    real *Usc_ptr, real *Utc_ptr, real *Unc_ptr,
    real *aeps_ptr, real *Utotal_ptr, real *beta_ptr,
    real dpitc, real dpits, real dbanc, real dbans,
    real *psi_pbs, real *psi_pbc, real *betac, real *betas)
/*---------------------------------------------------------------------------*/
/* Transforms global-frame velocities into the local blade-frame and         */
/* computes the total 2D velocity and inflow angle at the blade section.     */
/*                                                                           */
/*   Usc   - radial velocity in blade frame                                  */
/*   Utc   - tangential velocity in blade frame                              */
/*   Unc   - normal (out-of-plane) velocity in blade frame                   */
/*   aeps  - inflow angle epsilon (rad)                                      */
/*   Utotal - total local 2D velocity (tangential + normal)                  */
/*   beta  - total flapping angle (rad)                                      */
/*---------------------------------------------------------------------------*/
  {
    real Ux, Uy, Uz;
    real Ux_pb, Uy_pb, Uz_pb;
    real Ur_pb, Ut_pb;
    real beta;
    real tmp, tmp1, Utc, Unc, aeps;
    real Ut_omega, Ut_total, Utotal;

    if (count == 0 && dbg == 1)
      Message0("  %d Rotor %d - vel_xyz_to_blade\n", myid, i + 1);

    /* Flow-field velocities in global x, y, z. */
    Ux = C_U(cc, tc);
    Uy = C_V(cc, tc);
    Uz = C_W(cc, tc);

    /* Velocities transformed into the pitch/bank plane. */
    Ux_pb = Ux * dpitc - Uz * dpits;
    Uy_pb = Ux * dpits * dbans + Uy * dbanc + Uz * dpitc * dbans;
    Uz_pb = Ux * dpits * dbanc - Uy * dbans + Uz * dpitc * dbanc;

    /* Azimuth angle components. */
    *psi_pbs = sin(psi_pb);
    *psi_pbc = cos(psi_pb);

    /* Radial and tangential velocities in rotor plane. */
    Ur_pb =  Ux_pb * (*psi_pbc) + Uy_pb * (*psi_pbs);
    Ut_pb = -Ux_pb * (*psi_pbs) + Uy_pb * (*psi_pbc);

    /* Total flapping angle. */
    beta = bflco[i] - bflc[i] * (*psi_pbc) - bfls[i] * (*psi_pbs);
    if (isnan(beta))
      Message0("Error: Rotor %d beta returns NaN.\n", i + 1);
    *beta_ptr = beta;

    /* Flapping angle components. */
    *betas = sin(beta);
    *betac = cos(beta);

    /* Sign convention: +1 for negative rotation direction, -1 for positive. */
    tmp1 = 0.0;
    if      (rspe[i] > 0.0) tmp1 = -1.0;
    else if (rspe[i] < 0.0) tmp1 =  1.0;

    /* Radial, tangential, and normal velocities in flapping plane. */
    tmp       = Ur_pb * (*betac) + Uz_pb * (*betas);
    *Usc_ptr  = tmp;
    Utc       = tmp1 * Ut_pb;
    *Utc_ptr  = Utc;
    Unc       = -Ur_pb * (*betas) + Uz_pb * (*betac);
    *Unc_ptr  = Unc;

    /* Rotor tangential speed (m/s), ramped up linearly over i_start iters. */
    if (N_ITER >= i_start)
      Ut_omega = fabs(rspe[i]) * r_pb2;
    else
      Ut_omega = fabs(rspe[i]) * r_pb2 * ((real)N_ITER / (real)i_start);

    /* Total tangential velocity and resultant 2D velocity. */
    Ut_total  = Ut_omega + Utc;
    Utotal    = sqrt(Ut_total * Ut_total + Unc * Unc);
    *Utotal_ptr = Utotal;

    /* Inflow angle epsilon, wrapped to [-pi, pi]. */
    aeps = atan2(Unc, fabs(Ut_total));
    if      (aeps < -M_PI) aeps += 2.0 * M_PI;
    else if (aeps >  M_PI) aeps -= 2.0 * M_PI;
    *aeps_ptr = aeps;
  }

/*---------------------------------------------------------------------------*/
  void my_get_pitch(int i, cell_t cc, Thread *tc, int count, real r_pb,
    real psi_pb, real *theta_ptr, real psi_pbs, real psi_pbc)
/*---------------------------------------------------------------------------*/
/* Computes the total blade pitch angle at the given normalised radius,      */
/* including collective, cyclic, and linearly-interpolated twist.            */
/*---------------------------------------------------------------------------*/
  {
    int j;
    real theta, tmp;

    /* Find the blade section containing r_pb (linear interpolation of twist). */
    for (j = 0; j < nsec[i]; j++) {
      if      (r_pb > rout[i][j]) continue;
      else if (r_pb > rin[i][j])  break;
      else    Message("ERROR: Blade section radius not found.\n");
    }

    if ((rout[i][j] - rin[i][j]) == 0.0)
      tmp = 0.0;
    else
      tmp = twst[i][j] + (twst[i][j + 1] - twst[i][j]) /
            (rout[i][j] - rin[i][j]) * (r_pb - rin[i][j]);

    theta = bcop[i] + tmp - bcyc[i] * psi_pbc - bcys[i] * psi_pbs;

    /* Wrap to (-pi, pi]. */
    if      (theta < -M_PI) theta += 2.0 * M_PI;
    else if (theta >  M_PI) theta -= 2.0 * M_PI;
    *theta_ptr = theta;
  }


/*---------------------------------------------------------------------------*/
  void my_get_chord_alm_area(int i, cell_t cc, Thread *tc, int count, real r_pb,
    real *chord_ptr)
/*---------------------------------------------------------------------------*/
/* Returns the linearly interpolated chord length at normalised radius r_pb. */
/*---------------------------------------------------------------------------*/
  {
    int j;

    if (count == 0 && dbg == 1)
      Message0("  %d Rotor %d - get_chord\n", myid, i + 1);

    /* Find section and interpolate chord linearly. */
    for (j = 0; j < nsec[i]; j++) {
      if      (r_pb > rout[i][j]) continue;
      else if (r_pb > rin[i][j])  break;
      else    Message("ERROR: Blade section radius not found.\n");
    }

    if ((rout[i][j] - rin[i][j]) == 0.0)
      *chord_ptr = 0.0;
    else
      *chord_ptr = cin[i][j] + (cout[i][j] - cin[i][j]) /
                  (rout[i][j] - rin[i][j]) * (r_pb - rin[i][j]);
  }

/*---------------------------------------------------------------------------*/
  void my_get_re(int i, cell_t cc, Thread *tc, int count, real chord,
    real Utotal, real *Re_t_ptr)
/*---------------------------------------------------------------------------*/
/* Computes the local Reynolds number at the blade section.                  */
/*---------------------------------------------------------------------------*/
  {
    if (count == 0 && dbg == 1)
      Message0("  %d Rotor %d - get_re\n", myid, i + 1);

    if (C_MU_L(cc, tc) == 0.0) {
      Message("Error: Laminar viscosity is zero.\n");
      *Re_t_ptr = 0.0;
    }
    else {
      *Re_t_ptr = Utotal * chord * C_R(cc, tc) / C_MU_L(cc, tc);
    }
  }


/*---------------------------------------------------------------------------*/
  void my_get_ma(int i, cell_t cc, Thread *tc, int count, real Utotal,
    real *Ma_t_ptr)
/*---------------------------------------------------------------------------*/
/* Computes the local Mach number. Returns 0.2 for constant-density cases.  */
/*---------------------------------------------------------------------------*/
  {
    real kappa;

    if (rho_const) {
      *Ma_t_ptr = 0.2;
    }
    else {
      if ((C_CP(cc, tc) - C_RGAS(cc, tc)) == 0.0) {
        Message("Error: gamma for rotor %d is undefined (zero denominator).\n", i + 1);
        kappa = 0.0;
      }
      else {
        kappa = C_CP(cc, tc) / (C_CP(cc, tc) - C_RGAS(cc, tc));
      }

      if (sqrt(kappa * C_RGAS(cc, tc) * C_T(cc, tc)) == 0.0) {
        Message("Error: Speed of sound for rotor %d is zero.\n", i + 1);
        *Ma_t_ptr = 0.0;
      }
      else {
        *Ma_t_ptr = Utotal / sqrt(kappa * C_RGAS(cc, tc) * C_T(cc, tc));
      }
    }
  }


/*---------------------------------------------------------------------------*/
  void my_get_aoa(int i, cell_t cc, Thread *tc, int count, real aeps,
    real theta, real *alpha_t_ptr)
/*---------------------------------------------------------------------------*/
/* Computes the local angle of attack as epsilon + theta.                    */
/*---------------------------------------------------------------------------*/
  {
    *alpha_t_ptr = aeps + theta;

    if (isnan(*alpha_t_ptr) && i == 1 && errchk != 1) {
      Message0("Warning: alpha_t returns NaN at iter %d\n", count);
      errchk = 1;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_get_cl(int i, cell_t cc, Thread *tc, int count, real r_pb,
    real Ma_t, real Re_t, real alpha_t, real *CL_ptr)
/*---------------------------------------------------------------------------*/
/* Looks up and interpolates the lift coefficient from the airfoil tables    */
/* for the local Reynolds number, Mach number, and angle of attack.          */
/*---------------------------------------------------------------------------*/
  {
    int j,k,it,jt;
    real Re_min, Re_max, Ma_min, Ma_max;
    int ire_min = 0, ire_max = 0, ima_min = 0, ima_max = 0;
    real clre_min = 0, clre_max = 0, clre = 0, clma_min = 0, clma_max = 0, 
      clma = 0, cl = 0, rfac = 0;
    int iflag;
    real cl_m = 0, cl_p = 0;
  
    rfac=0.5;   /*0 --> cl=clma, 1 --> cl=clre*/
  
    /*Re_t=90000;*/
    /*Ma_t=0.96;*/
    /*alpha_t=3.0;*/
  
    alpha_t=alpha_t*180.0/M_PI;
  
    /*Re_min=min(1.2,Re_t);
    Re_max=max(1.3,Re_t);
    Message("Re_min %f Re_max %f \n",Re_min,Re_max);*/
  
    Re_min=1000000000.0;
    Re_max=-1000000000.0;
    Ma_min=1000000000.0;
    Ma_max=-1000000000.0;
  
    /*Assumed to be linear per blade section*/
    /*Loop through all blade sections*/
    j=0;
    iflag=0;
    while (j < nsec[i]) {
      if (r_pb > rout[i][j])
        j += 1;
      else if (r_pb > rin[i][j]) {
        /*Message("%f %f %f \n",rin[i][j],r_pb,rout[i][j]);*/
        break;
      }
      else
        Message("ERROR: Blade section radius not found. \n");
    }
  
    /*Message("%f %f %f %s \n",rin[i][j],r_pb,rout[i][j],type[i][j]);*/
  
    k=0;
    while (k < ktot) {
      if (strcmp(check_name[k], type[i][j]) == 0 ||
        strcmp(check_name[k], type[i][j+1]) == 0) {
        /*Message("  check %s type %s \n",check_name[k], type[i][j]);*/
        /*Message("  %d \n",itot[k]);*/
  
        it=0;
        while (it < itot[k]) {
          if (strcmp(clorcd[k][it],"cl") == 0) {
            /*Message("  %s \n",clorcd[k][it]);*/
            if (RE[k][it] < Re_min) {
              Re_min=RE[k][it];
              ire_min=it;
            }
            if (RE[k][it] > Re_max) {
              Re_max=RE[k][it];
              ire_max=it;
            }
            if (MA[k][it] < Ma_min) {
              Ma_min=MA[k][it];
              ima_min=it;
            }
            if (MA[k][it] > Ma_max) {
              Ma_max=MA[k][it];
              ima_max=it;
            }
          }
          it+=1;
        }
  
        jt=0;
  
        if (aoa[k][ire_min][jt] > alpha_t)
          Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
            k, ire_min);
  
        if (aoa[k][ire_max][jt] > alpha_t)
          Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
            k, ire_max);
  
        if (aoa[k][ima_min][jt] > alpha_t)
          Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
            k, ima_min);
  
        if (aoa[k][ima_max][jt] > alpha_t)
          Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
            k, ima_max);
  
        jt=jtot[k][ire_min]-1;
        if (aoa[k][ire_min][jt] < alpha_t)
          Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
            k, ire_min);
  
        jt=jtot[k][ire_max]-1;
        if (aoa[k][ire_max][jt] < alpha_t)
          Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
            k, ire_max);
  
        jt=jtot[k][ima_min]-1;
        if (aoa[k][ima_min][jt] < alpha_t)
          Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
            k, ima_min);
  
        jt=jtot[k][ima_max]-1;
        if (aoa[k][ima_max][jt] < alpha_t)
          Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
            k, ima_max);
  
      /*------------------------------------------------------------------*/
        if ((Re_t > Re_min) && (Re_t < Re_max)) {
          /*Message("RE interpolation \n");*/
  
          jt=0;
          while (jt < jtot[k][ire_min]) {
            if (aoa[k][ire_min][jt] > alpha_t) {
              clre_min=coeff[k][ire_min][jt]+
                (aoa[k][ire_min][jt] - alpha_t)/
                (aoa[k][ire_min][jt] - aoa[k][ire_min][jt-1]) *
                (coeff[k][ire_min][jt-1] - coeff[k][ire_min][jt]);
              break;
            }
            jt+=1;
          }
  
          jt=0;
          while (jt < jtot[k][ire_max]) {
            if (aoa[k][ire_max][jt] > alpha_t) {
              clre_max=coeff[k][ire_max][jt]+
                (aoa[k][ire_max][jt] - alpha_t)/
                (aoa[k][ire_max][jt] - aoa[k][ire_max][jt-1]) *
                (coeff[k][ire_max][jt-1] - coeff[k][ire_max][jt]);
              break;
            }
            jt+=1;
          }
          /*Message("clre_min %f clre_max %f \n",clre_min, clre_max);*/
  
          clre=clre_max+(Re_max-Re_t)/(Re_max-Re_min)*(clre_min-clre_max);
  
          /*Message("clre %f \n",clre);*/
        }
        /*------------------------------------------------------------------*/
        if (Re_t <= Re_min) {
        /*Message("Re_t smaller/equal \n");*/
  
          jt=0;
          while (jt < jtot[k][ire_min]) {
            if (aoa[k][ire_min][jt] > alpha_t) {
              clre=coeff[k][ire_min][jt]+
                (aoa[k][ire_min][jt] - alpha_t)/
                (aoa[k][ire_min][jt] - aoa[k][ire_min][jt-1]) *
                (coeff[k][ire_min][jt-1] - coeff[k][ire_min][jt]);
              break;
            }
            jt+=1;
          }
        }
      /*------------------------------------------------------------------*/
        if (Re_t >= Re_max) {
          /*Message("Re_t greater/equal \n");*/
          jt=0;
          while (jt < jtot[k][ire_max]) {
            if (aoa[k][ire_max][jt] > alpha_t) {
              clre=coeff[k][ire_max][jt]+
                (aoa[k][ire_max][jt] - alpha_t)/
                (aoa[k][ire_max][jt] - aoa[k][ire_max][jt-1]) *
                (coeff[k][ire_max][jt-1] - coeff[k][ire_max][jt]);
              break;
            }
            jt+=1;
          }
          /*Message("clre = clre_max = %f \n",clre);*/
        }
      
        /*------------------------------------------------------------------*/
        if ((Ma_t > Ma_min) && (Ma_t < Ma_max)) {
          /*Message("MA interpolation \n");*/
          jt=0;
          while (jt < jtot[k][ima_min]) {
            if (aoa[k][ima_min][jt] > alpha_t) {
              clma_min=coeff[k][ima_min][jt]+
                (aoa[k][ima_min][jt] - alpha_t)/
                (aoa[k][ima_min][jt] - aoa[k][ima_min][jt-1]) *
                (coeff[k][ima_min][jt-1] - coeff[k][ima_min][jt]);
              break;
            }
            jt+=1;
          }
  
          jt=0;
          while (jt < jtot[k][ima_max]) {
            if (aoa[k][ima_max][jt] > alpha_t) {
              clma_max=coeff[k][ima_max][jt]+
                (aoa[k][ima_max][jt] - alpha_t)/
                (aoa[k][ima_max][jt] - aoa[k][ima_max][jt-1]) *
                (coeff[k][ima_max][jt-1] - coeff[k][ima_max][jt]);
              break;
            }
            jt+=1;
          }
          /*Message("clma_min %f clma_max %f \n",clma_min, clma_max);*/
  
          clma=clma_max+(Ma_max-Ma_t)/(Ma_max-Ma_min)*(clma_min-clma_max);
  
          /*Message("clma %f \n",clma);*/
        }
        /*------------------------------------------------------------------*/
        if (Ma_t <= Ma_min) {
          /*Message("Ma_t smaller/equal \n");*/
  
          jt=0;
          while (jt < jtot[k][ima_min]) {
            if (aoa[k][ima_min][jt] > alpha_t) {
              clma=coeff[k][ima_min][jt]+
                (aoa[k][ima_min][jt] - alpha_t)/
                (aoa[k][ima_min][jt] - aoa[k][ima_min][jt-1]) *
                (coeff[k][ima_min][jt-1] - coeff[k][ima_min][jt]);
              break;
            }
            jt+=1;
          }
          /*Message("clma = clma_min = %f \n",clma);*/
  
        }
        /*------------------------------------------------------------------*/
        if (Ma_t >= Ma_max) {
          /*Message("Ma_t greater/equal \n");*/
          jt=0;
          while (jt < jtot[k][ima_max]) {
            if (aoa[k][ima_max][jt] > alpha_t) {
              clma=coeff[k][ima_max][jt]+
                (aoa[k][ima_max][jt] - alpha_t)/
                (aoa[k][ima_max][jt] - aoa[k][ima_max][jt-1]) *
                (coeff[k][ima_max][jt-1] - coeff[k][ima_max][jt]);
              break;
            }
            jt+=1;
          }
        }
        /*------------------------------------------------------------------*/
  
        if (strcmp(check_name[k], type[i][j]) == 0) {
          cl_m=clre*rfac+(1.0-rfac)*clma;
          iflag=iflag+1;
          /*Message("  check %s type %s \n",check_name[k], type[i][j]);*/
        }
  
        if (strcmp(check_name[k], type[i][j+1]) == 0) {
          cl_p=clre*rfac+(1.0-rfac)*clma;
          iflag=iflag+1;
          /*Message("  check %s type %s \n",check_name[k], type[i][j+1]);*/
        }
        if (iflag == 2) 
          break;
      }
      k+=1;
    }
  
    if(isnan(cl_m)) {
      if(i==1) {
        if(errchk !=1) {
          Message0("Warning, cl_m returns NaN at iter %d \n",count);
          errchk=1;
        }
      }
    }
  
    if(isnan(cl_p)) {
      if(i==1) {
        if(errchk !=1) {
          Message0("Warning, cl_p returns NaN at iter %d \n",count);
          errchk=1;
        }
      }
    }
  
    cl = cl_m + (cl_p-cl_m)/(rout[i][j]-rin[i][j])*(r_pb-rin[i][j]);
  
    if(isnan(cl)) {
      if(i==1) {
        if(errchk !=1) {
          Message0("Warning, cl returns NaN at iter %d \n",count);
          errchk=1;
        }
      }
    }
  
    *CL_ptr=cl;
  
    if(isnan(*CL_ptr)) {
      if(i==1) {
        if(errchk !=1) {
          Message0("Warning, *CL_ptr returns NaN at iter %d \n",count);
          errchk=1;
        }
      }
    } 
    
    /*Message("  rout %g rin %g alpha_t %g iflag %d \n",
              rout[i][j], rin[i][j], alpha_t,iflag);
    Message("  cl_m %g cl %g cl_p %g \n",cl_m,cl,cl_p);*/
  }

/*---------------------------------------------------------------------------*/
  void my_get_cd(int i, cell_t cc, Thread *tc, int count, real r_pb,
    real Ma_t, real Re_t, real alpha_t, real *CD_ptr)
/*---------------------------------------------------------------------------*/
/* Looks up and interpolates the drag coefficient from the airfoil tables    */
/* for the local Reynolds number, Mach number, and angle of attack.          */
/*---------------------------------------------------------------------------*/
{
  int j,k,it,jt;
  real Re_min, Re_max, Ma_min, Ma_max;
  int ire_min = 0, ire_max = 0, ima_min = 0, ima_max = 0;
  real cdre_min = 0, cdre_max = 0, cdre = 0, cdma_min = 0, cdma_max = 0, 
    cdma = 0, cd = 0, rfac = 0;
  int iflag;
  real cd_m = 0, cd_p = 0;


  rfac=0.5;   /*0 --> cd=cdma, 1 --> cd=cdre*/


  /*Re_t=600000;*/
  /*Ma_t=0.95;*/
  /*alpha_t=3.0;*/

  alpha_t=alpha_t*180.0/M_PI;


  /*Re_min=min(1.2,Re_t);
  Re_max=max(1.3,Re_t);
  Message("Re_min %f Re_max %f \n",Re_min,Re_max);*/

  Re_min=1000000000.0;
  Re_max=-1000000000.0;
  Ma_min=1000000000.0;
  Ma_max=-1000000000.0;


  /*Assumed to be linear per blade section*/
  /*Loop through all blade sections*/
  j=0;
  iflag=0;
  while (j < nsec[i]) {
    if (r_pb > rout[i][j])
      j += 1;
    else if (r_pb > rin[i][j])
      break;
    else
      Message("ERROR: Blade section radius not found. \n");
  }

  /*Message("%f %f %f %s \n",rin[i][j],r_pb,rout[i][j],type[i][j]);*/


  k=0;
  while (k < ktot) {
    if (strcmp(check_name[k], type[i][j]) == 0 ||
      strcmp(check_name[k], type[i][j+1]) == 0) {

      it=0;
      while (it < itot[k]) {
        if (strcmp(clorcd[k][it],"cd") == 0) {
          if (RE[k][it] < Re_min) {
            Re_min=RE[k][it];
            ire_min=it;
          }
          if (RE[k][it] > Re_max) {
            Re_max=RE[k][it];
            ire_max=it;
          }
          if (MA[k][it] < Ma_min) {
            Ma_min=MA[k][it];
            ima_min=it;
          }
          if (MA[k][it] > Ma_max) {
            Ma_max=MA[k][it];
            ima_max=it;
          }
        }
        it+=1;
      }

      jt=0;
      if (aoa[k][ire_min][jt] > alpha_t)
        Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
          k, ire_min);

      if (aoa[k][ire_max][jt] > alpha_t)
        Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
          k, ire_max);

      if (aoa[k][ima_min][jt] > alpha_t)
        Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
          k, ima_min);

      if (aoa[k][ima_max][jt] > alpha_t)
        Message("ERROR: Start AOA in airfoil table too big k %d i %d \n",
          k, ima_max);

      jt=jtot[k][ire_min]-1;
      if (aoa[k][ire_min][jt] < alpha_t)
        Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
          k, ire_min);

      jt=jtot[k][ire_max]-1;
      if (aoa[k][ire_max][jt] < alpha_t)
        Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
          k, ire_max);

      jt=jtot[k][ima_min]-1;
      if (aoa[k][ima_min][jt] < alpha_t)
        Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
          k, ima_min);

      jt=jtot[k][ima_max]-1;
      if (aoa[k][ima_max][jt] < alpha_t)
        Message("ERROR: End AOA in airfoil table too small k %d i %d \n",
          k, ima_max);

      /*------------------------------------------------------------------*/
      if ((Re_t > Re_min) && (Re_t < Re_max)) {
        /*Message("RE interpolation \n");*/

        jt=0;
        while (jt < jtot[k][ire_min]) {
          if (aoa[k][ire_min][jt] > alpha_t) {
            cdre_min=coeff[k][ire_min][jt]+
              (aoa[k][ire_min][jt] - alpha_t)/
              (aoa[k][ire_min][jt] - aoa[k][ire_min][jt-1]) *
              (coeff[k][ire_min][jt-1] - coeff[k][ire_min][jt]);
            break;
          }
          jt+=1;
        }

        jt=0;
        while (jt < jtot[k][ire_max]) {
          if (aoa[k][ire_max][jt] > alpha_t) {
            cdre_max=coeff[k][ire_max][jt]+
              (aoa[k][ire_max][jt] - alpha_t)/
              (aoa[k][ire_max][jt] - aoa[k][ire_max][jt-1]) *
              (coeff[k][ire_max][jt-1] - coeff[k][ire_max][jt]);
            break;
          }
          jt+=1;
        }
        /*Message("cdre_min %f cdre_max %f \n",cdre_min, cdre_max);*/

        cdre=cdre_max+(Re_max-Re_t)/(Re_max-Re_min)*(cdre_min-cdre_max);

      }
      /*------------------------------------------------------------------*/
      if (Re_t <= Re_min) {
        /*Message("Re_t smaller/equal \n");*/

        jt=0;
        while (jt < jtot[k][ire_min]) {
          if (aoa[k][ire_min][jt] > alpha_t) {
            cdre=coeff[k][ire_min][jt]+
              (aoa[k][ire_min][jt] - alpha_t)/
              (aoa[k][ire_min][jt] - aoa[k][ire_min][jt-1]) *
              (coeff[k][ire_min][jt-1] - coeff[k][ire_min][jt]);
            break;
          }
          jt+=1;
        }
        /*Message("cdre = cdre_min = %f \n",cdre);*/
      }
      /*------------------------------------------------------------------*/
      if (Re_t >= Re_max) {
        /*Message("Re_t greater/equal \n");*/

        jt=0;
        while (jt < jtot[k][ire_max]) {
          if (aoa[k][ire_max][jt] > alpha_t) {
            cdre=coeff[k][ire_max][jt]+
              (aoa[k][ire_max][jt] - alpha_t)/
              (aoa[k][ire_max][jt] - aoa[k][ire_max][jt-1]) *
              (coeff[k][ire_max][jt-1] - coeff[k][ire_max][jt]);
            break;
          }
          jt+=1;
        }
        /*Message("cdre = cdre_max = %f \n",cdre);*/
      }
      /*------------------------------------------------------------------*/
      if ((Ma_t > Ma_min) && (Ma_t < Ma_max)) {
        /*Message("MA interpolation \n");*/

        jt=0;
        while (jt < jtot[k][ima_min]) {
          if (aoa[k][ima_min][jt] > alpha_t) {
            cdma_min=coeff[k][ima_min][jt]+
              (aoa[k][ima_min][jt] - alpha_t)/
              (aoa[k][ima_min][jt] - aoa[k][ima_min][jt-1]) *
              (coeff[k][ima_min][jt-1] - coeff[k][ima_min][jt]);
            break;
          }
          jt+=1;
        }
        jt=0;
        while (jt < jtot[k][ima_max]) {
          if (aoa[k][ima_max][jt] > alpha_t) {
            cdma_max=coeff[k][ima_max][jt]+
              (aoa[k][ima_max][jt] - alpha_t)/
              (aoa[k][ima_max][jt] - aoa[k][ima_max][jt-1]) *
              (coeff[k][ima_max][jt-1] - coeff[k][ima_max][jt]);
            break;
          }
          jt+=1;
        }
        cdma=cdma_max+(Ma_max-Ma_t)/(Ma_max-Ma_min)*(cdma_min-cdma_max);
      }
      /*------------------------------------------------------------------*/
      if (Ma_t <= Ma_min) {
        /*Message("Ma_t smaller/equal \n");*/
        jt=0;
        while (jt < jtot[k][ima_min]) {
          if (aoa[k][ima_min][jt] > alpha_t) {
            cdma=coeff[k][ima_min][jt]+
              (aoa[k][ima_min][jt] - alpha_t)/
              (aoa[k][ima_min][jt] - aoa[k][ima_min][jt-1]) *
              (coeff[k][ima_min][jt-1] - coeff[k][ima_min][jt]);
            break;
          }
          jt+=1;
        }
        /*Message("cdma = cdma_min = %f \n",cdma);*/
      }
      /*------------------------------------------------------------------*/
      if (Ma_t >= Ma_max) {
        /*Message("Ma_t greater/equal \n");*/
        jt=0;
        while (jt < jtot[k][ima_max]) {
          if (aoa[k][ima_max][jt] > alpha_t) {
            cdma=coeff[k][ima_max][jt]+
              (aoa[k][ima_max][jt] - alpha_t)/
              (aoa[k][ima_max][jt] - aoa[k][ima_max][jt-1]) *
              (coeff[k][ima_max][jt-1] - coeff[k][ima_max][jt]);
            break;
          }
          jt+=1;
        }
      }
      /*------------------------------------------------------------------*/

      if (strcmp(check_name[k], type[i][j]) == 0) {
        cd_m=cdre*rfac+(1.0-rfac)*cdma;
        iflag=iflag+1;
        /*Message("  check %s type %s \n",check_name[k], type[i][j]);*/
      }

      if (strcmp(check_name[k], type[i][j+1]) == 0) {
        cd_p=cdre*rfac+(1.0-rfac)*cdma;
        iflag=iflag+1;
        /*Message("  check %s type %s \n",check_name[k], type[i][j+1]);*/
      }

      if (iflag == 2) 
        break;

    }
    k+=1;
  }

  cd = cd_m + (cd_p-cd_m)/(rout[i][j]-rin[i][j])*(r_pb-rin[i][j]);
  *CD_ptr=cd;
  if(isnan(*CD_ptr)) {
    if(i==1) {
      if(errchk !=1) {
        Message0("Warning, *CD_ptr returns NaN at iter %d \n",count);
        errchk=1;
      }
    }
  }
}

/*---------------------------------------------------------------------------*/
  void my_get_force(int i, cell_t cc, Thread *tc, real *factor_ptr,
    real *sum_area, int count, real CL, real CD,
    real chord, real aeps, real Utotal, real r_pb,
    real *Ftc_ptr, real *Fnc_ptr, real *Fsc_ptr, int isunst, int local_indx,
    real sample_vol)
/*---------------------------------------------------------------------------*/
/* Computes the blade tangential, normal, and radial force components at the */
/* current cell, applying a Gaussian (ALM) or geometric (BEM) distribution. */
/* Applies tip loss cutoff based on teff.                                    */
/*---------------------------------------------------------------------------*/
  {
    real FL, FD;
    real aepss, aepsc;
    real tmp;

    if (isunst == 1)
      tmp = 0.5 * Utotal * Utotal * C_R(cc, tc) * pow(cbrt(C_VOLUME(cc, tc)), 2);
    else
      tmp = chord * 0.5 * Utotal * Utotal * C_R(cc, tc) * factor_ptr[local_indx];

    if (isnan(tmp))
      Message0("Warning: force scaling factor is NaN at iter %d\n", count);

    FL = CL * tmp;
    FD = CD * tmp;

    /* Apply tip loss: zero lift beyond the tip effect radius. */
    if (r_pb > teff[i])
      FL = 0.0;

    aepss = sin(aeps);
    aepsc = cos(aeps);

    *Ftc_ptr = FD * aepsc - FL * aepss;
    *Fnc_ptr = FD * aepss + FL * aepsc;
    *Fsc_ptr = 0.0;

    if (isnan(*Fnc_ptr)) {
      Message0("Warning: Fnc is NaN at iter %d\n", count);
      errchk = 1;
    }
  }


/*---------------------------------------------------------------------------*/
  void my_force_blade_to_xyz(int i, cell_t cc, Thread *tc, int count,
    real Ftc, real Fnc, real Fsc, real beta, real psi_pb,
    real *Fx_ptr, real *Fy_ptr, real *Fz_ptr,
    real *Fz_pb_ptr, real *Ft_pb_ptr, real chord, real r_pb,
    real *gwf, int isunst, int *rotor_cell,
    real beta_pb, int *gauss_active, real *sum_gwf,
    real dpitc, real dpits, real dbanc, real dbans,
    real psi_pbs, real psi_pbc, real betas, real betac, real *deltabeta,
    int *cellsaboverotor, int *cells_ahead_behind_rotor,
    real *sigma_az, real *sigma_vert, real *gauss_denom,
    real *saz_fac, real *sve_fac, real *gde_fac, real *deltapsi,
    real *az_dist, real *vert_dist)
/*---------------------------------------------------------------------------*/
/* Rotates blade-frame forces back into global XYZ. In unsteady (ALM) mode, */
/* applies a Gaussian distribution based on the nearest blade position.      */
/*---------------------------------------------------------------------------*/
  {
    real Fr_pb, Ft_pb, Fz_pb;
    real Fx_pb, Fy_pb;
    real tmp1;
    real blade_psi_curr; /* Current blade azimuth angle (rad). */
    real blade_psi_init; /* Initial blade azimuth angle (rad). */
    real bsa;            /* Blade spacing angle (rad). */

    if (count == 0 && dbg == 1)
      Message0("  %d Rotor %d - force_blade_to_xyz\n", myid, i + 1);

    tmp1 = (rspe[i] > 0.0) ? -1.0 : 1.0;
  
    Fr_pb = Fsc * betac -          Fnc*betas;
    Ft_pb =             tmp1*Ftc;
    Fz_pb = Fsc  *betas +          Fnc*betac;
  
    Fx_pb = Fr_pb*psi_pbc - Ft_pb*psi_pbs;
    Fy_pb = Fr_pb*psi_pbs + Ft_pb*psi_pbc;

    /*If sim is unsteady, perform rotor position check, then remove source 
      terms where applicable. */
    if(isunst == 1){

      (*sigma_az) = chord / 2.0;
      (*sigma_vert) = (*sigma_az) / 2.0;
      (*gauss_denom) = 0.8;

      /* Calculate current position of blades, blade_psi_curr. */
      bsa = ((2.0 * M_PI)/nbld[i]);
      
      /* Loop through each blade position. */
      for(blade_psi_init = 0; blade_psi_init < (2.0 * M_PI); blade_psi_init = 
        blade_psi_init + bsa) {

        /* First, find where the rotor blades are at the current time-step, 
          blade_psi_curr. */
        blade_psi_curr = fmod(blade_psi_init + (rspe[i] * CURRENT_TIME), 2 * M_PI);

        /* Find vertical difference between current cell and current blade. */
        *deltabeta = fabs(beta_pb - beta); 

        /* Find angular difference between current cell and current blade. 
          Adding 2 pi to each psi angle prevents issues with wrapping 
          around 0.*/
        (*deltapsi) = fabs(psi_pb - blade_psi_curr);
        if ((*deltapsi) > M_PI) {
          (*deltapsi) = 2.0 * M_PI - (*deltapsi);
        }

        /* Continue only if the current blade is closest to this cell. There 
          are 3 logical conditions here. The first is to see if the difference 
          between the psi angles is less than half the blade spacing (bsa). 
          The second and third are for wrapping the angular differences around 2pi*/
        if((*deltapsi) < bsa/2.0){
          
          /* Multiply the radial position of the cell by the angular difference 
          between the cell and the blade's current location. */
          *az_dist = fabs(r_pb * rrl[i] * (*deltapsi));
          /* Similar to above for vertical distance. */
          *vert_dist = fabs(r_pb * rrl[i] * (*deltabeta));

          /* Calculate the Gaussian distribution multiplier. */
          *gwf = (1.0 / (*gauss_denom)) * exp( ( -(*az_dist) * (*az_dist) /
            ( (*sigma_az) * (*sigma_az) ) ) - ( (*vert_dist) * (*vert_dist) /
            ( (*sigma_vert) * (*sigma_vert) ) ) );
          *sum_gwf += *gwf;

          /* Multiply forces by Gaussian distrbution factor. */
          *Fx_ptr  = ( Fx_pb* dpitc + Fy_pb * dpits * dbans + Fz_pb * dpits * dbanc ) * (*gwf);
          *Fy_ptr  = (                Fy_pb * dbanc         - Fz_pb * dbans   ) * (*gwf);
          *Fz_ptr  = (-Fx_pb* dpits + Fy_pb * dpitc * dbans + Fz_pb * dpitc * dbanc ) * (*gwf);
          
          *Fz_pb_ptr = Fz_pb * ( *gwf );
          *Ft_pb_ptr = Ft_pb * ( *gwf );

          if( (*gwf) >= 1e-1 ) {
            (*gauss_active)++;
            if(abs( *vert_dist > (chord/2.0) ) ) {
              (*cellsaboverotor)++;
            }
            if(abs( *az_dist > (chord/2.0) ) ) {
              (*cells_ahead_behind_rotor)++;
            }
          }
        }
      }
    }
    else {
      *Fx_ptr  =  Fx_pb * dpitc + Fy_pb * dpits * dbans + Fz_pb * dpits * dbanc;
      *Fy_ptr  =                  Fy_pb * dbanc         - Fz_pb * dbans ; 
      *Fz_ptr  = -Fx_pb * dpits + Fy_pb * dpitc * dbans + Fz_pb * dpitc * dbanc;
    
      *Fz_pb_ptr = Fz_pb;
      *Ft_pb_ptr = Ft_pb;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_gauss_opt(int i, real sum_sample_vol, real sum_gwf, int gauss_active,
    int vol_sample_entries, int cellsaboverotor, int cells_ahead_behind_rotor,
    real *gauss_denom, real *sigma_vert, real *sigma_az,
    real *saz_fac, real *sve_fac, real *gde_fac)
/*---------------------------------------------------------------------------*/
/* Optional Gaussian parameter optimisation routine. Adjusts sigma_az to    */
/* minimise cells ahead/behind the rotor blade. Only active when             */
/* gauss_routine == 1.                                                       */
/*---------------------------------------------------------------------------*/
  {
    real gauss_sldty_old; /* Old Gaussian solidity. */
    real gauss_sldty; /* Gaussian solidity. */
    real rot_zone_cell_count; /* Number of cells in rotor zone. */
    real avg_vol; /* Average cell volume in rotor zone. */
    real avg_cell_size; /* Average cell size in rotor zone. */
    real global_sum_gwf; /* Global sum of Gaussian Weighting Factor. */
    real global_sum_sample_vol; /* Global sum of all sample's volumes. */
    real global_sum_gauss_active; /* Global sum of active Gaussian cells. */
    int global_vol_sample_entries; /* Total cell vol entries across nodes. */
    int total_cells_above_rotor; /* Total cells above/below rotor across 
      nodes. */
    int global_cells_ahead_behind_rotor; /* Total cells ahead/behind rotor across 
      nodes. */
    int delta_c_ah_bh_r; /* Difference in cell count ahead/behind rotor. */
    int delta_c_ab_r; /* Difference in cell count above/below rotor. */
    real delta_sigma_az; /* Difference in sigma_az value. */
    real delta_sigma_vert; /* Difference in sigma_vert value. */
    real delta_gauss_denom; /* Difference in Gaussian denominator value. */
    real delta_cabr_sigvert; /* Gradient of delta cells above/below rotor to 
      delta sigma_vert. */
    real delta_cahr_sigaz; /* Gradient of delta cells ahead/behind rotor to 
      delta sigma_az. */
    real sigma_az_new;
    real sigma_vert_new;
    real gauss_denom_new; 


    /* Unsure what sigma and width to use? Try activating this optimisation 
      routine. */
    if(gauss_routine){


      /* Cells located above/below and ahead/behind the rotor. 
      total_cells_above_rotor = PRF_GRSUM1(cellsaboverotor);*/
      global_cells_ahead_behind_rotor = PRF_GRSUM1(cells_ahead_behind_rotor);
      Message0("  global_cells_ahead_behind_rotor = %d \n", 
        global_cells_ahead_behind_rotor);

      /* Total cell volume sample across all nodes. 
      global_vol_sample_entries = PRF_GRSUM1(vol_sample_entries);*/

      /* Global sum of all sample's volumes. 
      global_sum_sample_vol = PRF_GRSUM1(sum_sample_vol);*/

      /* Calculate the rotor-local cell sizing using the samples taken earlier. 
      avg_vol = global_sum_sample_vol / global_vol_sample_entries;
      avg_cell_size = cbrt(avg_vol);
      Message0("  Average cell size in rotor zone = %f cm\n", 100.0 * avg_cell_size);*/

      /* Calculate the global sum of the Gaussian Weighting Factor. 
      global_sum_gwf = PRF_GRSUM1(sum_gwf);*/

      /* Calculate the global sum of the active Gaussian cells. 
      global_sum_gauss_active = PRF_GRSUM1(gauss_active);*/

      /* Calculate the number of cells in current rotor zone. 
      rot_zone_cell_count = M_PI * rrl[i] * rrl[i] / (avg_cell_size * avg_cell_size);*/

      /* Calculate the Gaussian solidity.
      gauss_sldty = global_sum_gauss_active / rot_zone_cell_count;
      Message0("  Gaussian solidity = %f, Actual solidity = %f \n", 
        gauss_sldty, sldity[i]); */
        
      /* Adjust the sigma values to lower the number of active Gaussian cells
        above/below and ahead/behind of the rotor. */

      /* First, find differece between new and old cell count above/below, 
        and ahead/behind rotor. 
      delta_c_ab_r = total_cells_above_rotor - total_cells_above_rotor_old[i];*/
      delta_c_ah_bh_r = global_cells_ahead_behind_rotor - global_cells_ahead_behind_rotor_old[i];
      Message0("  delta_c_ah_bh_r = %d \n", delta_c_ah_bh_r);

      /* Find difference between new and old sigma and denom values. 
      delta_sigma_vert = (*sigma_vert) - sigma_vert_old;
      delta_gauss_denom = (*gauss_denom) - gauss_denom_old;*/
      delta_sigma_az = (*sigma_az) - sigma_az_old;

      /* Save sigma and denom values. 
      sigma_vert_old = (*sigma_vert);
      gauss_denom_old = (*gauss_denom);*/
      sigma_az_old *= (*saz_fac);

      /* Save cell count above/below and ahead/below rotor. 
      total_cells_above_rotor_old[i] = total_cells_above_rotor;*/
      global_cells_ahead_behind_rotor_old[i] = global_cells_ahead_behind_rotor;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_get_source_terms(int i, cell_t cc, Thread *tc, int count,
    real *xsource_ptr, real *ysource_ptr, real *zsource_ptr,
    real *xsource_old_ptr, real *ysource_old_ptr, real *zsource_old_ptr,
    real Fx, real Fy, real Fz, real urf_source, int local_indx)
/*---------------------------------------------------------------------------*/
/* Writes under-relaxed momentum source terms to per-cell memory arrays.     */
/* On the first iteration (istflag == 1) no URF blending is applied.         */
/*---------------------------------------------------------------------------*/
  {
    real vol;

    if (count == 0 && dbg == 1)
      Message0("  %d Rotor %d - get_source_terms\n", myid, i + 1);

    vol = C_VOLUME(cc, tc);
    if (vol == 0.0) {
      Message0("Warning: Rotor %d cell volume is zero.\n", i + 1);
      xsource_ptr[local_indx] = 0.0;
      ysource_ptr[local_indx] = 0.0;
      zsource_ptr[local_indx] = 0.0;
    }
    else if (istflag != 1 || N_ITER == 0) {
      /* Blended update (all iterations after the first ADJUST run). */
      xsource_ptr[local_indx] = (-Fx / vol) * urf_source +
        (1.0 - urf_source) * xsource_old_ptr[local_indx];
      ysource_ptr[local_indx] = (-Fy / vol) * urf_source +
        (1.0 - urf_source) * ysource_old_ptr[local_indx];
      zsource_ptr[local_indx] = (-Fz / vol) * urf_source +
        (1.0 - urf_source) * zsource_old_ptr[local_indx];
    }
    else {
      /* First iteration: no blending. */
      xsource_ptr[local_indx] = -Fx / vol;
      ysource_ptr[local_indx] = -Fy / vol;
      zsource_ptr[local_indx] = -Fz / vol;
    }

    /* Store for next iteration's under-relaxation. */
    xsource_old_ptr[local_indx] = xsource_ptr[local_indx];
    ysource_old_ptr[local_indx] = ysource_ptr[local_indx];
    zsource_old_ptr[local_indx] = zsource_ptr[local_indx];
  }

/*---------------------------------------------------------------------------*/
  void my_get_cone(int i, real thrust)
/*---------------------------------------------------------------------------*/
/* Placeholder for a blade coning angle calculation. Not yet implemented.    */
/*---------------------------------------------------------------------------*/
  {
  /*real lb, lbx;                                                                                               
    real wb, wbx;                                                                                               
    real cf, cfx;                                                                                               
    real conum, coden;                                                                                          
  																											  
    lb = thrust / nbld[i]; Lift force per blade (thrust assumed as lift, N)
    lbx = (3/4) * rrl[i]; Location where the total lift force acts (assumed parabolic, m)
    
    wb = 116.5; Weight of each blade (kg)
    wbx = (1/2) * rrl[i]; Centre of mass per blade (m)
    cf = wb * ((rspe[i]*rrl[i])*(rspe[i]*rrl[i])) / rrl[i]; Centrifugal force (N)
    cfx = rrl[i]/sqrt(2); Location where centrifugal force acts (m)
    conum = lb * lbx; Numerator
    coden = sqrt((cf*cfx)*(cf*cfx)+(wb*wbx)*(wb*wbx));
    bflco[i] = (conum/coden) - (wb*wbx)/(cf*cfx);
    Message0("  %d Rotor %d: new cone angle = %d \n",myid,i+1,bflco[i]);*/
  }
 

/*---------------------------------------------------------------------------*/
  void my_get_ct(int i, real thrust, real coeff_denom)
/*---------------------------------------------------------------------------*/
/* Computes the rotor thrust coefficient CT[i].                              */
/*---------------------------------------------------------------------------*/
  {
    if (coeff_denom == 0.0) {
      Message("Error: Thrust coeff denominator for rotor %d is zero.\n", i + 1);
      CT[i] = 0.0;
    }
    else {
      CT[i] = thrust / coeff_denom;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_get_cmx(int i, real Mx_pb, real coeff_denom)
/*---------------------------------------------------------------------------*/
/* Computes the rotor roll moment coefficient CMX[i].                        */
/*---------------------------------------------------------------------------*/
  {
    if (coeff_denom == 0.0) {
      Message("Error: X-moment coeff denominator for rotor %d is zero.\n", i + 1);
      CMX[i] = 0.0;
    }
    else {
      CMX[i] = Mx_pb / coeff_denom;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_get_cmy(int i, real My_pb, real coeff_denom)
/*---------------------------------------------------------------------------*/
/* Computes the rotor pitch moment coefficient CMY[i].                       */
/*---------------------------------------------------------------------------*/
  {
    if (coeff_denom == 0.0) {
      Message("Error: Y-moment coeff denominator for rotor %d is zero.\n", i + 1);
      CMY[i] = 0.0;
    }
    else {
      CMY[i] = My_pb / coeff_denom;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_force_moment_balancing(int i, real thr)
/*---------------------------------------------------------------------------*/
/* Adjusts the main and tail rotor thrust targets to balance fuselage forces */
/* and moments. Only active when forba/momba flags are set.                  */
/*---------------------------------------------------------------------------*/
  {
    Domain *d = Get_Domain(1);
    Thread *t_object = Lookup_Thread(d,fuseid[0]);
    real cttr_new, deltacttr, cttr_old;
    real ctmr_new, deltactmr, ctmr_old;
    real force_x, force_z;
    real moment_x;
    real coeff_denom;

    if (d == NULL) {
      Message0("Error: Fuselage thread not found for rotor %d.\n", i+1);
      return;
    }
  
    if (t_object == NULL) {
      Message0("Error: Fuselage thread not found for rotor %d.\n", i+1);
      return;
    }
  
    coeff_denom = 0.5 * (RP_Get_Real("reference-density")) * rspe[i] * 
      rspe[i] * rrl[i] * rrl[i] * rrl[i] * rrl[i] * M_PI;
    
  
    if(forba[i] == 1) {
      if (i == 0) {
        
        if (dbg ==1){
          Message0("Compensating main rotor for fuselage downwash \n");
          Message0("Old main rotor %d target thrust coeff = %f \n", i+1, ctd[i]);
        }
  
        ctmr_old = ctd[i];
  
        Compute_Force_And_Moment(d, t_object, cg, force, moment, FALSE);
  
        force_x = force[0];
        force_z = force[2];
        moment_x = moment[0];
  
        if (dbg ==1){
          Message0("X Force on Fuselage = %f N\n",force_x);
          Message0("Z Force on Fuselage = %f N\n",force_z);
          Message0("Current X Force Produced by Rotor = %f N \n", 
            thr*sin(dpit[i]));
          Message0("Required disk pitch angle = %f deg\n",
            (180/M_PI) * atan((-1.0 * force_x)/((ac_mass * 9.81) - force_z)));
        }
        
        ctmr_new = ((ac_mass * 9.81 - force_z)/cos(dpit[i])) /
          coeff_denom;
  
        deltactmr = ctmr_new - ctmr_old;
  
        ctd[i] = ctmr_old + trdf[i]*deltactmr;
        
        if (dbg ==1)
          Message0("New main rotor thrust coeff target = %f\n", ctd[i]);
      }
      if (i == 1 && momba[i] == 1) {
        cttr_old = ctd[i];
        if (dbg ==1){
          Message0("Correcting tail rotor thrust for main rotor "
            "torque... \n");
          Message0("Old tail rotor thrust coeff target = %f\n", 
            cttr_old);
        }
        
        cttr_new =  (fabs(sum_torque[i-1]) / ((dskco[i][0] - dskco[i-1][0]) * 
          coeff_denom)) / sin(-1.0*dban[i]);
        /*The line above corrects tail rotor thrust target with main rotor torque. Disable when using single rotor. */
  
        deltacttr = cttr_new - cttr_old;
        ctd[i] = cttr_old + trdf[i]*deltacttr;
  
        if (dbg ==1)
          Message0("New tail rotor thrust coeff target = %f\n", ctd[i]);
      }
    }
  }

/*---------------------------------------------------------------------------*/
  void my_update_trimming(int i, int trufq_co[], int trufq_cy[],
    real trdf_co[], real trdf_cy[],
    int up_co[], int up_cy[],
    int *flag_jac_ptr,
    real *ct_l_ptr,  real *ct_s_ptr,  real *ct_h_ptr,
    real *bcop_l_ptr, real *bcop_s_ptr, real *bcop_h_ptr,
    real *cmx_lc_ptr, real *cmx_ls_ptr, real *cmx_s_ptr,
    real *cmx_hc_ptr, real *cmx_hs_ptr,
    real *cmy_lc_ptr, real *cmy_ls_ptr, real *cmy_s_ptr,
    real *cmy_hc_ptr, real *cmy_hs_ptr,
    real *bcyc_l_ptr, real *bcyc_s_ptr, real *bcyc_h_ptr,
    real *bcys_l_ptr, real *bcys_s_ptr, real *bcys_h_ptr,
    real *ct_lc_ptr,  real *ct_ls_ptr,
    real *ct_hc_ptr,  real *ct_hs_ptr,
    real *cmx_l_ptr,  real *cmx_h_ptr,
    real *cmy_l_ptr,  real *cmy_h_ptr,
    real thr)
/*---------------------------------------------------------------------------*/
/* Runs the Jacobian-based trimming algorithm for collective and/or cyclic   */
/* pitch, targeting desired CT, CMX, and CMY values.                         */
/*---------------------------------------------------------------------------*/
  {
    real del_co;
    real Ctmp, Btmp;
    real ct_l, ct_s, ct_h, bcop_l, bcop_s, bcop_h;
    real cmx_lc, cmx_ls, cmx_s, cmx_hc, cmx_hs;
    real cmy_lc, cmy_ls, cmy_s, cmy_hc, cmy_hs;
    real bcyc_l, bcyc_s, bcyc_h, bcys_l, bcys_s, bcys_h;
    real del_cyc, del_cys;
    real dcmx_bcys, dcmx_bcyc, dcmy_bcys, dcmy_bcyc, det_jac;
    real ct_lc, ct_ls, ct_hc, ct_hs, cmx_l, cmx_h, cmy_l, cmy_h;
    real dct_bcop, dct_bcyc, dct_bcys, dcmx_bcop, dcmy_bcop;
  
    /*Message0("  %d update_trimming \n",myid); */
  
    /*THRUST ONLY*/
    if (trmco[i] == 1 && trmcy[i] == 0) {
      if (up_co[i] < trufq_co[i]) {
        up_co[i] += 1;
      }
      else {
        trnw = 1;
        if (*flag_jac_ptr == 0) {
          if (dbg == 1)
            Message0("Trimming collective pitch: ON \n");
  
          /*modify_bcop();*/
          /*store ct_s*/
          *ct_s_ptr=CT[i];
          *bcop_s_ptr=bcop[i];
  
          if (dbg == 1){
            Message0("CT(actual)= %f CT(desired)= %f \n", CT[i], ctd[i]);
            Message0("%d CO ALPHA PLUS \n",myid);
          }

          bcop[i]=*bcop_s_ptr+dalpha;
          *bcop_h_ptr=bcop[i];
          *flag_jac_ptr = 1;
          goto COEND;
        }
  
        if (*flag_jac_ptr == 1) {
          if (dbg == 1)
            Message0("%d CO ALPHA MINUS \n",myid);
  
          /*store ct_h*/
          *ct_h_ptr=CT[i];
          bcop[i]=*bcop_s_ptr-dalpha;
          *bcop_l_ptr=bcop[i];
  
          *flag_jac_ptr = -1;
  
          goto COEND;
        }
  
        if (*flag_jac_ptr == -1) {
          /*store ct_l*/
          *ct_l_ptr=CT[i];
  
          /*determine new angle*/
          ct_l   = *ct_l_ptr;
          ct_s   = *ct_s_ptr;
          ct_h   = *ct_h_ptr;
          bcop_l = *bcop_l_ptr;
          bcop_s = *bcop_s_ptr;
          bcop_h = *bcop_h_ptr;
  
  
          Ctmp=(ct_s-ct_h-(bcop_s-bcop_h)/(bcop_s-bcop_l)*(ct_s-ct_l))/
            ((bcop_s-bcop_h)*(bcop_h-bcop_l));
  
          Btmp=(ct_s-ct_l-Ctmp*(bcop_s*bcop_s - bcop_l*bcop_l))/
            (bcop_s-bcop_l);
  
          del_co=trdf_co[i]*(ctd[i]-ct_s)/(Btmp+2.*Ctmp*bcop_s);

          /*limiter*/
          if(del_co > dalpha*limiter) 
            del_co=dalpha*limiter;
          if(del_co < -dalpha*limiter) 
            del_co=-dalpha*limiter;
          bcop[i]=bcop_s+del_co;
  
          if (dbg == 1){
            Message0("bcop_l %f bcop_s %f bcop_h %f \n", bcop_l*180/M_PI, 
              bcop_s*180/M_PI, bcop_h*180/M_PI);
            Message0("ct_l %f ct_s %f ct_h %f \n", ct_l, ct_s, ct_h);
  
            Message0("old collective pitch %f \n", bcop_s*180/M_PI);
            Message0("%f = new collective for rotor %d, \n ",
              bcop[i]*180/M_PI, i+1);
            Message0("difference %f \n", del_co*180/M_PI);
          }

          *flag_jac_ptr = 10;
  
          up_co[i]=1;
          goto COEND;
        }
      }
    }
    COEND:
  
    /*MOMENT ONLY*/
    if (trmco[i] == 0 && trmcy[i] == 1) {
      if (up_cy[i] < trufq_cy[i]) {
        up_cy[i] += 1;
      }
      else {
        trnw = 1;
        if (*flag_jac_ptr == 0) {
          Message0("    %d Trimming cyclic pitch: ON \n",myid);
          /*modify_bcys_and_bcyc();*/
  
  
          /*store cmx_s cmy_s*/
          *cmx_s_ptr=CMX[i];
          *cmy_s_ptr=CMY[i];
  
          *bcys_s_ptr = bcys[i];
          *bcyc_s_ptr = bcyc[i];
  
  
          Message0("    %d CMX(actual)= %f CMX(desired)= %f \n",myid, CMX[i], cmxd[i]);
          Message0("    %d CMY(actual)= %f CMY(desired)= %f \n",myid, CMY[i], cmyd[i]);
  
  
          Message0("%d CYC ALPHA PLUS \n",myid);
  
          bcys[i]=*bcys_s_ptr;
          bcyc[i]=*bcyc_s_ptr+dalpha;
  
          /**bcys_h_ptr = bcys[i];*/
          *bcyc_h_ptr = bcyc[i];
  
          *flag_jac_ptr = 2;
  
          goto CYEND;
        }
  
        if (*flag_jac_ptr == 2) {
  
          /*store cmx_h cmy_h*/
          *cmx_hc_ptr=CMX[i];
          *cmy_hc_ptr=CMY[i];
  
          Message0("%d CYC ALPHA MINUS \n",myid);
  
          bcys[i]=*bcys_s_ptr;
          bcyc[i]=*bcyc_s_ptr-dalpha;
  
          /**bcys_l_ptr = bcys[i];*/
          *bcyc_l_ptr = bcyc[i];
  
  
          *flag_jac_ptr = -2;
  
          goto CYEND;
        }
  
  
        if (*flag_jac_ptr == -2) {
          /*store cmx_h cmy_h*/
          *cmx_lc_ptr=CMX[i];
          *cmy_lc_ptr=CMY[i];
  
          Message0("%d CYS ALPHA PLUS \n",myid);
  
          bcys[i]=*bcys_s_ptr+dalpha;
          bcyc[i]=*bcyc_s_ptr;
  
          *bcys_h_ptr = bcys[i];
          /**bcyc_l_ptr = bcyc[i];*/
  
  
          *flag_jac_ptr = 3;
  
          goto CYEND;
        }
  
  
        if (*flag_jac_ptr == 3) {
          /*store cmx_h cmy_h*/
          *cmx_hs_ptr=CMX[i];
          *cmy_hs_ptr=CMY[i];
  
          Message0("%d CYS ALPHA MINUS \n",myid);
  
          bcys[i]=*bcys_s_ptr-dalpha;
          bcyc[i]=*bcyc_s_ptr;
  
          *bcys_l_ptr = bcys[i];
          /**bcyc_l_ptr = bcyc[i];*/
  
  
          *flag_jac_ptr = -3;
  
          goto CYEND;
        }
  
  
        if (*flag_jac_ptr == -3) {
  
  
          /*store cmx_l cmy_l*/
          *cmx_ls_ptr=CMX[i];
          *cmy_ls_ptr=CMY[i];
  
  
          /*determine new angle*/
  
  
          cmx_lc   = *cmx_lc_ptr;
          cmx_ls   = *cmx_ls_ptr;
          cmy_lc   = *cmy_lc_ptr;
          cmy_ls   = *cmy_ls_ptr;
          cmx_s    = *cmx_s_ptr;
          cmy_s    = *cmy_s_ptr;
          cmx_hc   = *cmx_hc_ptr;
          cmx_hs   = *cmx_hs_ptr;
          cmy_hc   = *cmy_hc_ptr;
          cmy_hs   = *cmy_hs_ptr;
  
          bcyc_l = *bcyc_l_ptr;
          bcyc_s = *bcyc_s_ptr;
          bcyc_h = *bcyc_h_ptr;
          bcys_l = *bcys_l_ptr;
          bcys_s = *bcys_s_ptr;
          bcys_h = *bcys_h_ptr;
  
  
          /*d(cmx)/d(bcys)*/
          Ctmp=(cmx_s-cmx_hs-(bcys_s-bcys_h)/(bcys_s-bcys_l)*(cmx_s-cmx_ls))/
            ((bcys_s-bcys_h)*(bcys_h-bcys_l));
          Btmp=(cmx_s-cmx_ls-Ctmp*(bcys_s*bcys_s - bcys_l*bcys_l))/
            (bcys_s-bcys_l);
          dcmx_bcys = Btmp+2.*Ctmp*bcys_s;
  
  
          /*d(cmx)/d(bcyc)*/
          Ctmp=(cmx_s-cmx_hc-(bcyc_s-bcyc_h)/(bcyc_s-bcyc_l)*(cmx_s-cmx_lc))/
            ((bcyc_s-bcyc_h)*(bcyc_h-bcyc_l));
          Btmp=(cmx_s-cmx_lc-Ctmp*(bcyc_s*bcyc_s - bcyc_l*bcyc_l))/
            (bcyc_s-bcyc_l);
          dcmx_bcyc = Btmp+2.*Ctmp*bcyc_s;
  
  
          /*d(cmy)/d(bcys)*/
          Ctmp=(cmy_s-cmy_hs-(bcys_s-bcys_h)/(bcys_s-bcys_l)*(cmy_s-cmy_ls))/
            ((bcys_s-bcys_h)*(bcys_h-bcys_l));
          Btmp=(cmy_s-cmy_ls-Ctmp*(bcys_s*bcys_s - bcys_l*bcys_l))/
            (bcys_s-bcys_l);
          dcmy_bcys = Btmp+2.*Ctmp*bcys_s;
  
  
          /*d(cmy)/d(bcyc)*/
          Ctmp=(cmy_s-cmy_hc-(bcyc_s-bcyc_h)/(bcyc_s-bcyc_l)*(cmy_s-cmy_lc))/
            ((bcyc_s-bcyc_h)*(bcyc_h-bcyc_l));
          Btmp=(cmy_s-cmy_lc-Ctmp*(bcyc_s*bcyc_s - bcyc_l*bcyc_l))/
            (bcyc_s-bcyc_l);
          dcmy_bcyc = Btmp+2.*Ctmp*bcyc_s;
  
  
          /*det_jac*/
          det_jac = dcmy_bcyc*dcmx_bcys - dcmy_bcys*dcmx_bcyc;
          /*if(det_jac <=0.0) det_jac=-1.0;
          if(det_jac > 0.0) det_jac=1.0;*/
          det_jac=1./det_jac;
  
  
          /*Message("%f %f %f %f \n",dcmy_bcyc,dcmx_bcys,dcmy_bcys,dcmx_bcyc);*/
  
  
          del_cyc = trdf_cy[i]*det_jac * (dcmx_bcys*(cmyd[i]-cmy_s) -
            dcmy_bcys*(cmxd[i]-cmx_s));
  
          /*Message("%f %f %f %f %f %f \n",trdf_cy[i],det_jac,dcmx_bcys,
            (cmyd[i]-cmy_s),
            dcmy_bcys,(cmxd[i]-cmx_s));*/
          /*Message("%f \n", dcmx_bcys*(cmyd[i]-cmy_s) -
            dcmy_bcys*(cmxd[i]-cmx_s));*/
  
  
          del_cys = trdf_cy[i]*det_jac * (-dcmx_bcyc*(cmyd[i]-cmy_s) +
            dcmy_bcyc*(cmxd[i]-cmx_s));
  
          /*limiter*/
          if(del_cyc > dalpha*limiter) 
            del_cyc=dalpha*limiter;
          if(del_cys > dalpha*limiter) 
            del_cys=dalpha*limiter;
          if(del_cyc < -dalpha*limiter) 
            del_cyc=-dalpha*limiter;
          if(del_cys < -dalpha*limiter) 
            del_cys=-dalpha*limiter;
  
  
          bcyc[i]=bcyc_s+del_cyc;
          bcys[i]=bcys_s+del_cys;
  
  
          *flag_jac_ptr = 10;
  
          up_cy[i]=1;
          goto CYEND;
        }
      }
    }
    CYEND:
  
  
    /*THRUST AND MOMENT*/
    if ((trmco[i] == 1) && (trmcy[i] == 1)) {
      if (up_co[i] < trufq_co[i]) /*will be updated at the collective pitch update frequency*/ {
        up_co[i] += 1;
      }
      else {
        trnw = 1;
        /*Message("flag_jac %d",*flag_jac_ptr);*/
        if (*flag_jac_ptr == 0) {
          /*modify_bcop();*/
  
          /*store ct_s cmx_s cmy_s*/
          *ct_s_ptr  = CT[i];
          *cmx_s_ptr = CMX[i];
          *cmy_s_ptr = CMY[i];
  
          *bcop_s_ptr = bcop[i];
          *bcys_s_ptr = bcys[i];
          *bcyc_s_ptr = bcyc[i];
  
          if (dbg == 1) {
            Message0("Trimming collective and cyclic pitch: ON \n");
            Message0("CT(actual)= %f CT(desired)= %f \n", CT[i], ctd[i]);
            Message0("CMX(actual)= %f CMX(desired)= %f \n", CMX[i], cmxd[i]);
            Message0("CMY(actual)= %f CMY(desired)= %f \n", CMY[i], cmyd[i]);
          }
          
          /*Message0("%d CO ALPHA PLUS \n",myid);*/
          bcop[i]=*bcop_s_ptr+dalpha;
          bcys[i]=*bcys_s_ptr;
          bcyc[i]=*bcyc_s_ptr;
  
          *bcop_h_ptr=bcop[i];
  
          *flag_jac_ptr = 1;
  
          goto COYEND;
        }
  
        if (*flag_jac_ptr == 1) {
  
          /*store ct_h cmx_h cmy_h*/
          *ct_h_ptr  = CT[i];
          *cmx_h_ptr = CMX[i];
          *cmy_h_ptr = CMY[i];
  
  
          /*Message0("%d CO ALPHA MINUS \n",myid);*/
          bcop[i]=*bcop_s_ptr-dalpha;
          bcys[i]=*bcys_s_ptr;
          bcyc[i]=*bcyc_s_ptr;
  
          *bcop_l_ptr=bcop[i];
  
          *flag_jac_ptr = -1;
  
          goto COYEND;
        }
  
  
        if (*flag_jac_ptr == -1) {
          /*store ct_l cmx_l cmy_l*/
          *ct_l_ptr  = CT[i];
          *cmx_l_ptr = CMX[i];
          *cmy_l_ptr = CMY[i];
  
  
          /*Message0("%d CYC ALPHA PLUS \n",myid); */
          bcop[i]=*bcop_s_ptr;
          bcys[i]=*bcys_s_ptr;
          bcyc[i]=*bcyc_s_ptr+dalpha;
  
          *bcyc_h_ptr = bcyc[i];
  
          *flag_jac_ptr = 2;
  
          goto COYEND;
        }
  
  
        if (*flag_jac_ptr == 2) {
          /*store ct_hc cmx_hc cmy_hc*/
          *ct_hc_ptr  = CT[i];
          *cmx_hc_ptr = CMX[i];
          *cmy_hc_ptr = CMY[i];
  
  
          /*Message0("%d CYC ALPHA MINUS \n",myid); */
          bcop[i]=*bcop_s_ptr;
          bcys[i]=*bcys_s_ptr;
          bcyc[i]=*bcyc_s_ptr-dalpha;
  
          *bcyc_l_ptr = bcyc[i];
  
          *flag_jac_ptr = -2;
  
          goto COYEND;
        }
  
  
        if (*flag_jac_ptr == -2) {
          /*store ct_lc cmx_lc cmy_lc*/
          *ct_lc_ptr  = CT[i];
          *cmx_lc_ptr = CMX[i];
          *cmy_lc_ptr = CMY[i];
  
  
          /*Message0("%d CYS ALPHA PLUS \n",myid);*/
          bcop[i]=*bcop_s_ptr;
          bcys[i]=*bcys_s_ptr+dalpha;
          bcyc[i]=*bcyc_s_ptr;
  
          *bcys_h_ptr = bcys[i];
  
          *flag_jac_ptr = 3;
  
          goto COYEND;
        }
  
  
        if (*flag_jac_ptr == 3) {
          /*store ct_hs cmx_hs cmy_hs*/
          *ct_hs_ptr  = CT[i];
          *cmx_hs_ptr = CMX[i];
          *cmy_hs_ptr = CMY[i];
  
  
          /*Message0("%d CYS ALPHA MINUS \n", myid);*/
          bcop[i]=*bcop_s_ptr;
          bcys[i]=*bcys_s_ptr-dalpha;
          bcyc[i]=*bcyc_s_ptr;
  
          *bcys_l_ptr = bcys[i];
  
          *flag_jac_ptr = -3;
  
          goto COYEND;
        }
  
        if (*flag_jac_ptr == -3) {
          /*store ct_ls cmx_ls cmy_ls*/
          *ct_ls_ptr  = CT[i];
          *cmx_ls_ptr = CMX[i];
          *cmy_ls_ptr = CMY[i];
  
  
          /*determine new angle*/
  
          ct_s     = *ct_s_ptr;
          cmx_s    = *cmx_s_ptr;
          cmy_s    = *cmy_s_ptr;
          ct_l     = *ct_l_ptr;
          cmx_l    = *cmx_l_ptr;
          cmy_l    = *cmy_l_ptr;
          ct_h     = *ct_h_ptr;
          cmx_h    = *cmx_h_ptr;
          cmy_h    = *cmy_h_ptr;
          ct_lc    = *ct_lc_ptr;
          cmx_lc   = *cmx_lc_ptr;
          cmy_lc   = *cmy_lc_ptr;
          ct_hc    = *ct_hc_ptr;
          cmx_hc   = *cmx_hc_ptr;
          cmy_hc   = *cmy_hc_ptr;
          ct_ls    = *ct_ls_ptr;
          cmx_ls   = *cmx_ls_ptr;
          cmy_ls   = *cmy_ls_ptr;
          ct_hs    = *ct_hs_ptr;
          cmx_hs   = *cmx_hs_ptr;
          cmy_hs   = *cmy_hs_ptr;
  
          bcop_l = *bcop_l_ptr;
          bcop_s = *bcop_s_ptr;
          bcop_h = *bcop_h_ptr;
          bcyc_l = *bcyc_l_ptr;
          bcyc_s = *bcyc_s_ptr;
          bcyc_h = *bcyc_h_ptr;
          bcys_l = *bcys_l_ptr;
          bcys_s = *bcys_s_ptr;
          bcys_h = *bcys_h_ptr;
  
  
          /*d(ct)/d(bcop)*/
          Ctmp=(ct_s-ct_h-(bcop_s-bcop_h)/(bcop_s-bcop_l)*(ct_s-ct_l))/
            ((bcop_s-bcop_h)*(bcop_h-bcop_l));
          Btmp=(ct_s-ct_l-Ctmp*(bcop_s*bcop_s - bcop_l*bcop_l))/
            (bcop_s-bcop_l);
          dct_bcop = Btmp+2.*Ctmp*bcop_s;
  
  
          /*d(cmx)/d(bcop)*/
          Ctmp=(cmx_s-cmx_h-(bcop_s-bcop_h)/(bcop_s-bcop_l)*(cmx_s-cmx_l))/
            ((bcop_s-bcop_h)*(bcop_h-bcop_l));
          Btmp=(cmx_s-cmx_l-Ctmp*(bcop_s*bcop_s - bcop_l*bcop_l))/
            (bcop_s-bcop_l);
          dcmx_bcop = Btmp+2.*Ctmp*bcop_s;
  
  
          /*d(cmy)/d(bcop)*/
          Ctmp=(cmy_s-cmy_h-(bcop_s-bcop_h)/(bcop_s-bcop_l)*(cmy_s-cmy_l))/
            ((bcop_s-bcop_h)*(bcop_h-bcop_l));
          Btmp=(cmy_s-cmy_l-Ctmp*(bcop_s*bcop_s - bcop_l*bcop_l))/
            (bcop_s-bcop_l);
          dcmy_bcop = Btmp+2.*Ctmp*bcop_s;
  
  
          /*d(ct)/d(bcys)*/
          Ctmp=(ct_s-ct_hs-(bcys_s-bcys_h)/(bcys_s-bcys_l)*(ct_s-ct_ls))/
            ((bcys_s-bcys_h)*(bcys_h-bcys_l));
          Btmp=(ct_s-ct_ls-Ctmp*(bcys_s*bcys_s - bcys_l*bcys_l))/
            (bcys_s-bcys_l);
          dct_bcys = Btmp+2.*Ctmp*bcys_s;
  
  
          /*d(cmx)/d(bcys)*/
          Ctmp=(cmx_s-cmx_hs-(bcys_s-bcys_h)/(bcys_s-bcys_l)*(cmx_s-cmx_ls))/
            ((bcys_s-bcys_h)*(bcys_h-bcys_l));
          Btmp=(cmx_s-cmx_ls-Ctmp*(bcys_s*bcys_s - bcys_l*bcys_l))/
            (bcys_s-bcys_l);
          dcmx_bcys = Btmp+2.*Ctmp*bcys_s;
  
  
          /*d(cmy)/d(bcys)*/
          Ctmp=(cmy_s-cmy_hs-(bcys_s-bcys_h)/(bcys_s-bcys_l)*(cmy_s-cmy_ls))/
            ((bcys_s-bcys_h)*(bcys_h-bcys_l));
          Btmp=(cmy_s-cmy_ls-Ctmp*(bcys_s*bcys_s - bcys_l*bcys_l))/
            (bcys_s-bcys_l);
          dcmy_bcys = Btmp+2.*Ctmp*bcys_s;
  
  
          /*d(ct)/d(bcyc)*/
          Ctmp=(ct_s-ct_hc-(bcyc_s-bcyc_h)/(bcyc_s-bcyc_l)*(ct_s-ct_lc))/
            ((bcyc_s-bcyc_h)*(bcyc_h-bcyc_l));
          Btmp=(ct_s-ct_lc-Ctmp*(bcyc_s*bcyc_s - bcyc_l*bcyc_l))/
            (bcyc_s-bcyc_l);
          dct_bcyc = Btmp+2.*Ctmp*bcyc_s;
  
  
          /*d(cmx)/d(bcyc)*/
          Ctmp=(cmx_s-cmx_hc-(bcyc_s-bcyc_h)/(bcyc_s-bcyc_l)*(cmx_s-cmx_lc))/
            ((bcyc_s-bcyc_h)*(bcyc_h-bcyc_l));
          Btmp=(cmx_s-cmx_lc-Ctmp*(bcyc_s*bcyc_s - bcyc_l*bcyc_l))/
            (bcyc_s-bcyc_l);
          dcmx_bcyc = Btmp+2.*Ctmp*bcyc_s;
  
  
          /*d(cmy)/d(bcyc)*/
          Ctmp=(cmy_s-cmy_hc-(bcyc_s-bcyc_h)/(bcyc_s-bcyc_l)*(cmy_s-cmy_lc))/
            ((bcyc_s-bcyc_h)*(bcyc_h-bcyc_l));
          Btmp=(cmy_s-cmy_lc-Ctmp*(bcyc_s*bcyc_s - bcyc_l*bcyc_l))/
            (bcyc_s-bcyc_l);
          dcmy_bcyc = Btmp+2.*Ctmp*bcyc_s;
  
  
  
          /*det_jac*/
          det_jac = dct_bcop*dcmy_bcyc*dcmx_bcys +
          dct_bcyc*dcmy_bcys*dcmx_bcop +
          dct_bcys*dcmy_bcop*dcmx_bcyc -
          dct_bcys*dcmy_bcyc*dcmx_bcop -
          dct_bcyc*dcmy_bcop*dcmx_bcys -
          dct_bcop*dcmy_bcys*dcmx_bcyc;
  
          det_jac=1./det_jac;
  
  
          del_co=trdf_co[i]*det_jac * (
            (dcmy_bcyc*dcmx_bcys - dcmy_bcys*dcmx_bcyc)*(ctd[i]-ct_s) +
            (dct_bcys*dcmx_bcyc - dct_bcyc*dcmx_bcys)*(cmyd[i]-cmy_s) +
            (dct_bcyc*dcmy_bcys - dct_bcys*dcmy_bcyc)*(cmxd[i]-cmx_s));
  
          del_cyc=trdf_cy[i]*det_jac * (
            (dcmy_bcys*dcmx_bcop - dcmy_bcop*dcmx_bcys)*(ctd[i]-ct_s) +
            (dct_bcop*dcmx_bcys - dct_bcys*dcmx_bcop)*(cmyd[i]-cmy_s) +
            (dct_bcys*dcmy_bcop - dct_bcop*dcmy_bcys)*(cmxd[i]-cmx_s));
  
          del_cys=trdf_cy[i]*det_jac * (
            (dcmy_bcop*dcmx_bcyc - dcmy_bcyc*dcmx_bcop)*(ctd[i]-ct_s) +
            (dct_bcyc*dcmx_bcop - dct_bcop*dcmx_bcyc)*(cmyd[i]-cmy_s) +
            (dct_bcop*dcmy_bcyc - dct_bcyc*dcmy_bcop)*(cmxd[i]-cmx_s));
  
  
          /*limiter*/
          if(del_co  > dalpha*limiter) 
            del_co = dalpha*limiter;
          if(del_cyc > dalpha*limiter) 
            del_cyc= dalpha*limiter;
          if(del_cys > dalpha*limiter) 
            del_cys= dalpha*limiter;
          if(del_co  < -dalpha*limiter) 
            del_co = -dalpha*limiter;
          if(del_cyc < -dalpha*limiter) 
            del_cyc= -dalpha*limiter;
          if(del_cys < -dalpha*limiter) 
            del_cys= -dalpha*limiter;

          bcop[i] = bcop_s + del_co;
          bcyc[i] = bcyc_s + del_cyc;
          bcys[i] = bcys_s + del_cys;
  
          if (dbg == 1){
            Message0("%f = new collective for rotor %d, \n "
              "%f = new cyclic cos for rotor %d, \n "
              "%f = new cyclic sin for rotor %d\n",
              bcop[i]*180/M_PI, i+1, bcyc[i]*180/M_PI, i+1, bcys[i]*180/M_PI,
              i+1);
          }
  
          *flag_jac_ptr = 10;
      
          up_co[i]=1;

          /*my_force_moment_balancing(i, sum_thrust[i]);*/

          goto COYEND;
        }
      }
    }
    COYEND:;
  }

#endif

/*---------------------------------------------------------------------------*/
DEFINE_EXECUTE_ON_LOADING(my_on_loading, libname)
/*---------------------------------------------------------------------------*/
/* Registers UDM names when the UDF library is loaded.                       */
/*---------------------------------------------------------------------------*/
{
  Message0("!! Info: Ensure you have pre-allocated at least %d UDM locations "
    "in the GUI/TUI !!\n", NUM_UDM + 2);
  Set_User_Memory_Name(0,  "Distance from Rotor Centre");
  Set_User_Memory_Name(1,  "Azimuth Angle");
  Set_User_Memory_Name(2,  "Pitch Angle at Section");
  Set_User_Memory_Name(3,  "Total Velocity at Section");
  Set_User_Memory_Name(4,  "Chord Length at Section");
  Set_User_Memory_Name(5,  "Cell in Rotor Zone");
  Set_User_Memory_Name(6,  "Azimuthal Sigma Value");
  Set_User_Memory_Name(7,  "Angle of Attack at Section");
  Set_User_Memory_Name(8,  "Angle Between Cell and Rotor Blade");
  Set_User_Memory_Name(9,  "Azimuthal Distance to Rotor Blade");
  Set_User_Memory_Name(10, "Gaussian Weighting Factor");
  Set_User_Memory_Name(11, "Fz term at Section");
  Set_User_Memory_Name(12, "Local Index");
}

/*---------------------------------------------------------------------------*/
DEFINE_ON_DEMAND(my_rotor_inputs)
/*---------------------------------------------------------------------------*/
/* Reads rotor geometry, blade, and trimming parameters from the GUI (scheme */
/* panel) and broadcasts them from host to solver nodes.                     */
/*---------------------------------------------------------------------------*/
{

  int number_of_sections, ic;
  char exst_type[100][30];
  const char *type_new[max_rotor_zones][max_sections];
  int col_pitch, cyc_pitch, fuse_forces, fuse_moments;
  real bcop_udf[max_rotor_zones], bcyc_udf[max_rotor_zones], bcys_udf[max_rotor_zones];
  int trmco_dat[max_rotor_zones], trmcy_dat[max_rotor_zones];
  register int j,k,m;
  char exfname[100];
  FILE *inp_fp;
  int ch, itrm, inpt_it, pos, ffound;
  real inpt_data1[max_rotor_zones][max_rotor_zones], inpt_data2[max_rotor_zones][max_rotor_zones];
  real dskco_0[max_rotor_zones], dskco_1[max_rotor_zones], dskco_2[max_rotor_zones];
  real rsec_flatten[200], csec_flatten[200], twst_flatten[200];
  char type_flatten[6000];
  real RE_flatten[5000], MA_flatten[5000], jtot_flatten[5000];
  char check_name_flatten[3000], clorcd_flatten[50000];
  real aoa_flatten[400000], coeff_flatten[400000];



  #if !RP_NODE
    Pointer l = rpgetvar ("my-rotor");
  
    Domain *d;
    Thread *t;
    d = Get_Domain(1);

    if (d == NULL) {
      Message0("Error: Domain is NULL.\n");
      return; // Exit the function to prevent further execution
    }
  
    if (dbg == 1)
      Message("\n\n\nIn Define on demand - rotor_inputs \n");
  
    nrtz = RP_Get_List_Length("my-rotor");
  
    if(RP_Get_Integer("my-list-length-check") == 1) {
      if(dbg == 1)
        Message("\n\nPassing variables to the solver");
    }
    else
      Message("\n\nPlease click on change/create... button to pass the variables to solver");
  
  
    for (j=0; j<nrtz; j++) {
      /*save old pitch angles and trimming flags from .dat file*/
      bcop_udf[j] = bcop[j];
      bcyc_udf[j] = bcyc[j];
      bcys_udf[j] = bcys[j];
      trmco_dat[j] = trmco[j];
      trmcy_dat[j] = trmcy[j];
  
      nbld[j] = fixnum_arg(List_Ref (List_Ref (List_Ref (l,j),1),0), "");
      rrl[j]  = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),1), "");
      rspe[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),2), "");
      teff[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),3), "");
      teff[j] = teff[j]*0.01;
  
      dskco[j][0] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),4), "");
      dskco_0[j] = dskco[j][0];
      dskco[j][1] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),5), "");
      dskco_1[j] = dskco[j][1];
      dskco[j][2] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),6), "");
      dskco_2[j] = dskco[j][2];
      dpit[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),7), "");
      dban[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),8), "");

      bcop[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),9), "");
      bcys[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),10), "");
      bcyc[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),11), "");
      bflco[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),12), "");
      bfls[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),13), "");
      bflc[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),1),14), "");
      fzon[j] = 0;
      
      nsec[j] = fixnum_arg(List_Ref (List_Ref (List_Ref (l,j),2),0), "");
  
      number_of_sections = nsec[j] ;
  
      for (k=0; k<number_of_sections; k++) {
        rsec[j][k] = flonum_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),2),2),k), "");
        csec[j][k] = flonum_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),2),3),k), "");
        twst[j][k] = flonum_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),2),4),k), "");

        type_new[j][k] = string_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),2),5),k), "");
        (void)strcpy(type[j][k],type_new[j][k]);

        rsec_flatten[j*20 + k] = rsec[j][k];
        csec_flatten[j*20 + k] = csec[j][k];
        twst_flatten[j*20 + k] = twst[j][k];
        strncpy(&type_flatten[(j*20 + k)*30], type[j][k], 30);
      }
  
      if (fixnum_arg(List_Ref (List_Ref (List_Ref (l,j), 3),0), "") == 1) {
        itrm = 1;
        col_pitch = fixnum_arg(List_Ref (List_Ref (List_Ref (l,j),3),1), "");
        if (col_pitch == 1) {
          trmco[j] = 1;
          ctd[j]   = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),3),7), "");
          if (dbg == 1)
            Message("\n Desired Thrust Coeff in Model = %f",ctd[j]);
        }
        else {
          trmco[j] = 0;
        }
        cyc_pitch = fixnum_arg(List_Ref (List_Ref (List_Ref (l,j),3),2), "");
  
        if (cyc_pitch == 1) {
          trmcy[j] = 1;
          if (col_pitch == 1) {
            cmxd[j]  = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),3),8), "");
            cmyd[j]  = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),3),9), "");
            fuse_forces = fixnum_arg(List_Ref (List_Ref (List_Ref (l,j),3),3), "");
            if (fuse_forces == 1) {
              forba[j] = 1;
              for (m = 1; m < 10; m++) {
                forba[m] = 1;
              }
              fuse_moments = fixnum_arg(List_Ref (List_Ref (List_Ref (l,j),3),4), "");
              if (fuse_moments == 1) {
                momba[j] = 1;
                for(m = 1; m < 10; m++) {
                  momba[m] = 1;
                }
                cmgx[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),3),11), "");
                cmgy[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),3),12), "");
                cmgz[j] = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),3),13), "");
              }
              if (fuse_moments == 0) {
                momba[j] = 0;
              }
              if (fixnum_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),3),10),0),"") == 0) {
                Message("\n Please select a surface for fuselage face thread");
                fuseid[j] = 0;
              }
              else {
                fuseid[j] = fixnum_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),3),10),0), "");
              }
            }
            else {
              forba[j] = 0;
            }
          }
        }
        else {
          trmcy[j] = 0;
        }
        trufq[j] = fixnum_arg(List_Ref (List_Ref (List_Ref (l,j),3),5), "");
        trdf[j]  = flonum_arg(List_Ref (List_Ref (List_Ref (l,j),3),6), "");
      }
    }
    /*Message("DONE");*/
    
    int num_data_files = sizeof(data_fname) / sizeof(data_fname[0]);
    for (j=0; j<nrtz; j++) {
      for (int h = 0; h < num_data_files; h++) {
        Message("Searching for rotor %d %s...\n", j+1, data_fname[h]);
        sprintf(exfname, "rotor_%d_%s", j+1, data_fname[h]);
        inp_fp = fopen(exfname, "r");

        if (inp_fp == NULL) {
          Message("File %s for rotor %d does not exist or cannot be opened!\n",
            data_fname[h], j+1);
          ffound = 0;
        }
        else {
          if (dimport == 0){
            break;
          }
          Message("Data file %s for rotor %d found, importing...\n", 
            data_fname[h], j+1);
          ffound = 1;
          fseek(inp_fp, -1, SEEK_END);
          pos = ftell(inp_fp);
          ch = fgetc(inp_fp);
          if (isdigit(ch)){
            while (ch != '\n') {
              if (dbg == 1)
                Message("The current character in the file for rotor %d "
                  "is: '%c'\n", j+1, ch);
              
              pos--;
              fseek(inp_fp, pos, SEEK_SET);
              ch = fgetc(inp_fp);
            }
  
            if (dbg == 1)
              Message("Last number found at position %d in the file "
                "for rotor %d.\n", pos, j+1);
       
            fscanf(inp_fp, "%d %lf %lf", &inpt_it, &inpt_data1[h][j], 
              &inpt_data2[h][j]);
          }
          if (isalpha(ch)){
            ffound = 0;
            Message("File was previously initialized without writing data.\n");
          }
          fclose(inp_fp);
        }
      }
      
      if (ffound == 1) {
        bcop[j] = inpt_data1[0][j]* M_PI / 180.0;
        bcyc[j] = inpt_data1[1][j]* M_PI / 180.0;
        bcys[j] = inpt_data1[2][j]* M_PI / 180.0;
        ctd[j]  = inpt_data2[6][j];
        cmxd[j] = inpt_data2[7][j];
        cmyd[j] = inpt_data2[8][j];
        Message("Imported values for rotor %d: bcop = %f deg, bcyc = %f deg, "
          "bcys = %f deg, ctd = %f, cmxd = %f, cmyd = %f\n",
          j+1, bcop[j]/ (M_PI / 180.0), bcyc[j]/ (M_PI / 180.0), 
          bcys[j]/ (M_PI / 180.0), ctd[j], cmxd[j], cmyd[j]);
      }
    }
  #endif

  /* ALL the information resides on the host. Now distribute to nodes */
  Message0("Sending data to nodes...\n");

  /*Only host sends*/
  host_to_node_int_2(nrtz,itrm);
  host_to_node_int(nbld,10);
  host_to_node_real(rrl,10);
  host_to_node_real(rspe,10);
  host_to_node_real(teff,10);
  host_to_node_real(dskco_0,10);
  host_to_node_real(dskco_1,10);
  host_to_node_real(dskco_2,10);
  /*host_to_node_real(dskco,30);*/
  host_to_node_real(dpit,10);
  host_to_node_real(dban,10);

  host_to_node_int(trmco,10);
  host_to_node_int(trmcy,10);
  host_to_node_int(trmco_dat,10);
  host_to_node_int(trmcy_dat,10);
  host_to_node_int(forba,10);
  host_to_node_int(momba,10);
  host_to_node_int(fuseid,10);
  host_to_node_real(trdf,10);
  host_to_node_int(trufq,10);
  host_to_node_real(ctd,10);
  host_to_node_real(cmxd,10);
  host_to_node_real(cmyd,10);
  host_to_node_real(cmgx,10);
  host_to_node_real(cmgy,10);
  host_to_node_real(cmgz,10);

  host_to_node_real(bcop,10);
  host_to_node_real(bcop_udf,10);

  host_to_node_real(bcyc,10);
  host_to_node_real(bcyc_udf,10);
  host_to_node_real(bcys,10);
  host_to_node_real(bcys_udf,10);

  host_to_node_real(bflco,10);
  host_to_node_real(bfls,10);
  host_to_node_real(bflc,10);

  host_to_node_int(fzon,10);
  host_to_node_int(nsec,10);

  /*host_to_node_string(type,6000);*/

  host_to_node_real(rsec_flatten, 200);
  host_to_node_real(csec_flatten, 200);
  host_to_node_real(twst_flatten, 200);
  host_to_node_string(type_flatten, 6000);

  /*Only nodes receive*/

  #if !RP_HOST
    if (dbg == 1) {
      Message0 ("\n\n %d Rotor Model Inputs\n",myid);
      Message0 (" %d ------------------\n",myid);
      Message0 (" %d General Inputs \n",myid);
      Message0 (" %d ------------------\n",myid);
      Message0 (" %d Number of Rotor zones = %d\n", myid, nrtz);
    }
    for (j=0; j<nrtz; j++) {
      dskco[j][0] = dskco_0[j];
      dskco[j][1] = dskco_1[j];
      dskco[j][2] = dskco_2[j];
      if (dbg == 1) {
        Message0 ("\n\n %d Rotor zone %d\n", myid, j+1);
        Message0 (" %d Number of Blades = %d\n",myid, nbld[j]) ;
        Message0 (" %d Rotor Radius = %g\n",myid, rrl[j]) ;
        Message0 (" %d Rotor Speed = %g\n",myid, rspe[j]) ;
        Message0 (" %d Tip Effect = %g\n",myid, teff[j]) ;
        
        Message0 (" %d Rotor Disk Origin - X = %g\n",myid, dskco[j][0]) ;
        Message0 (" %d Rotor Disk Origin - Y = %g\n",myid, dskco[j][1]) ;
        Message0 (" %d Rotor Disk Origin - Z = %g\n",myid, dskco[j][2]) ;

        Message0 (" %d Rotor Disk Pitch Angle = %g\n",myid, dpit[j]*180/M_PI) ;
        Message0 (" %d Rotor Disk Bank Angle = %g\n",myid, dban[j]*180/M_PI) ;
        /*if(trmco[j] == 1 && istflag != 1 || trmcy[j] == 1 && istflag != 1 ||
           trmco_dat[j] == 1 || trmcy_dat[j] == 1)
        {
          bcop[j]=bcop_udf[j];
          Message0 (" %d   GUI COLLECTIVE BLADE PITCH IGNORED\n",myid);
        }*/
        Message0 (" %d Blade Pitch - Collective = %g\n",myid, bcop[j]*180/M_PI) ;
        /*if(trmco[j] == 1 && istflag != 1 || trmcy[j] == 1 && istflag != 1 ||
           trmco_dat[j] == 1 || trmcy_dat[j] == 1)
        {
          bcyc[j]=bcyc_udf[j];
          bcys[j]=bcys_udf[j];
          Message0 (" %d   GUI CYCLIC BLADE PITCH IGNORED\n",myid);
        }*/
        Message0 (" %d Blade Pitch - Cyclic Sin = %g\n",myid, bcys[j]*180/M_PI) ;
        Message0 (" %d Blade Pitch - Cyclic Cos = %g\n",myid, bcyc[j]*180/M_PI) ;
        Message0 (" %d Blade Flapping - Collective = %g\n",myid, bflco[j]*180/M_PI) ;
        Message0 (" %d Blade Flapping - Cyclic Sin = %g\n",myid, bfls[j]*180/M_PI) ;
        Message0 (" %d Blade Flapping - Cyclic Cos = %g\n",myid, bflc[j]*180/M_PI) ;
        Message0 (" %d Rotor Face Zone ID = %d\n",myid, fzon[j]) ;
        Message0 (" %d ------------------\n",myid);
        Message0 (" %d Geometry Inputs\n",myid);
        Message0 (" %d ------------------\n",myid);
        Message0 (" %d Number of Blade sections = %d\n",myid, nsec[j]);
      }
      number_of_sections = nsec[j] ;
      for (k=0; k<number_of_sections; k++) {
        rsec[j][k] = rsec_flatten[j*20 + k];
        csec[j][k] = csec_flatten[j*20 + k];
        twst[j][k] = twst_flatten[j*20 + k];
        strncpy(type[j][k], &type_flatten[(j*20 + k)*30], 30);
        if (dbg == 1) {
          Message0 ("\n %d Radius of section = %g\n",myid, rsec[j][k]) ;
          Message0 (" %d Chord of section = %g\n",myid, csec[j][k]) ;
          Message0 (" %d Twist of section = %g\n",myid, twst[j][k]*180/M_PI) ;
          /*(void)strcpy(type[j][k],type_new[j][k]);*/
            /*Message ("Name of the file = %s\n", type_new[j][k]);*/
          Message0 (" %d Name of the file = %s\n",myid, type[j][k]);
        }
      }

      if(itrm == 1) {
        if (dbg == 1) {
          Message0 (" %d ------------------\n",myid);
          Message0 (" %d Trimming Inputs\n",myid);
          Message0 (" %d ------------------\n",myid);
          Message0 (" %d Collective Pitch = %d\n",myid, trmco[j]);
          Message0 (" %d Cyclic Pitch = %d\n",myid, trmcy[j]);
          Message0 (" %d Update Frequency = %d\n",myid, trufq[j]);
          Message0 (" %d Damping Factor = %g\n",myid, trdf[j]);
        }
        if (trmco[j] == 1) {
          main_rot_thrust_target_init = ctd[0];
          if (dbg == 1) {
            Message0 (" %d Desired thrust coefficient = %g\n",myid, ctd[j]);
            Message0 (" %d ------------------\n",myid);
          }
          Message0 (" %d Rotor %d thrust target = %g kN\n",myid, j+1, 
            0.5*RP_Get_Real("reference-density")*
            rspe[j] * rspe[j] * rrl[j] * rrl[j] * M_PI * rrl[j] * rrl[j] * ctd[j]/1000);
        }
        else {
          if (dbg == 1) {
            Message0 (" %d ------------------\n",myid);
            Message0 (" %d Collective Pitch not selected\n",myid);
            Message0 (" %d ------------------\n",myid);
          }
        }
        if (dbg == 1) {
          if (trmcy[j] == 1) {
            Message0 (" %d Desired x-momentum coefficient = %g\n",myid, cmxd[j]);
            Message0 (" %d Desired y-momentum coefficient = %g\n",myid, cmyd[j]);
            if (j == 0) {
              for (m=0; m<10; m++) {                         
                if (forba[m] == 1) {
                  if (momba[m] == 1) {
                    Message0 (" %d Fuselage forces and moments to be "
                      "balanced\n", myid);
                    Message0 (" %d Co-ordinates of fuselage CG = [%f,%f,%f]\n",
                      myid, cmgx[m], cmgy[m], cmgz[m]);
                  }
                  else {
                    Message0 (" %d Fuselage forces to be balanced\n",myid);
                  }
                  Message0 (" %d Fuselage Face Thread ID = %d\n",myid, 
                    fuseid[m]);
                  break;        
                }                                                                        
              }   
            }
            Message0 (" %d ------------------\n",myid);
          }
          else {
            Message0 (" %d ------------------\n",myid);
            Message0 (" %d Cyclic pitch not selected\n",myid);
            Message0 (" %d -------------------------\n",myid);
          }
        }
      }
      else
        Message0 ("\n %d Trimming is not selected\n",myid);
    }

    /*CALCULATE rout, rin, cout, cin from rsec, csec*/
    for (j=0; j<nrtz; j++) {
      number_of_sections = nsec[j]-1 ;
      for (k=0; k<number_of_sections; k++) {
        rin[j][k]=rsec[j][k];
        cin[j][k]=csec[j][k];

        rout[j][k]=rsec[j][k+1];
        cout[j][k]=csec[j][k+1];
      }
    }
    for (j=0; j<nrtz; j++) {
      number_of_sections = nsec[j]-1 ;
      for (k=0; k<number_of_sections; k++) {
        if (dbg == 1) {
          Message0("%d Rotor %d Section %d | rin %g rout %g | cin %g cout %g \n",myid,
            j+1,k,rin[j][k], rout[j][k], cin[j][k], cout[j][k]);
        }
      }
    }  
  
    /*input_airfoil_tables(&ktot,file_name);*/
    /*CALCULATE ktot*/
    ktot=1;
    for (j=0; j<nrtz; j++) {
      number_of_sections = nsec[j] ;
      for (k=0; k<number_of_sections; k++) {
        /*Message ("TEST\n") ;*/
        for(ic = 0; ic < ktot; ic++) {
          if (strcmp(type[j][k],exst_type[ic]) == 0) {
            /*Message("yes\n");*/
            goto LABEL2;
          }
        }
        (void)strcpy(exst_type[ktot-1],type[j][k]);
        ktot=ktot+1;
        LABEL2:;
      }
    }
  
    ktot=ktot-1;
    if (dbg == 1) {
      Message0("%d ktot %d \n",myid,ktot);
      Message0("%d nrtz %d \n",myid,nrtz);
    }
    
    for(ic=0; ic<ktot; ic++){
      if (dbg == 1)
        Message0("%d Airfoils found %s\n",myid, exst_type[ic]);

      (void)strcpy(file_name[ic],exst_type[ic]);
      (void)strcat(file_name[ic],".dat");
      (void)strcpy(exst_type[ic],"");
      if (dbg == 1)
        Message0("%d Airfoil files existing in wkdir: %s\n",myid, file_name[ic]);
    }

    my_get_solidity();
  
    /*deg_to_rad(dpit, dban, bcop, bcys, bcyc, bflco, bfls, bflc, twst, &dalpha);*/
    /*dalpha=dalpha*M_PI/180.0;*/
  
    my_obtain_cell_id();

    if (istflag == 1) 
      my_allocate_memory(factor, xsource, ysource, zsource,
        xsource_old, ysource_old, zsource_old);
  
    my_geom_factor();
  
    my_start_trimming(trufq_co, trufq_cy, trdf_co, trdf_cy, up_co, up_cy);

  #endif

  node_to_host_int_1(ktot);
  for (j=0; j<ktot; j++) {
    node_to_host_string(file_name[j],30);
  }

  #if !RP_NODE
    my_read_in_airfoil_tables(); /*-180deg to +180deg*/

    /*istflag=1;*/

    /*determine if rho=const*/
    thread_loop_c(t,d) {
      if(FLUID_THREAD_P(t)) {
        if (DENSITY_METHOD(t) == RHO_CONSTANT) {
          if (dbg == 1) {
            Message("DENSITY CONST \n");
            Message("Mach number interpolation ignored \n");
          }
          rho_const = TRUE;
        }
        else {
          rho_const = FALSE;
        }
      }
    }
  #endif


  if (dbg == 0) 
    Message0("\nRotor passed to solver.\n");

  #if !RP_NODE
    for (int k = 0; k < ktot; k++) {
      for (int i = 0; i < itot[k]; i++) {
        RE_flatten[k*50 + i] = RE[k][i];
        MA_flatten[k*50 + i] = MA[k][i];
        jtot_flatten[k*50 + i] = jtot[k][i];
        for (int j=0; j<jtot[k][i]; j++) {
          aoa_flatten[k*50 + i*80 + j] = aoa[k][i][j];
          coeff_flatten[k*50 + i*80 + j] = coeff[k][i][j];
        }
      }
    }


    for(int k = 0; k < ktot; k++){
      strncpy(&check_name_flatten[k*30], check_name[k], 30);
      for(int i = 0; i < itot[k]; i++) 
        strncpy(&clorcd_flatten[k*500 + i*10], clorcd[k][i], 10);
    }

  #endif

  host_to_node_real(RE_flatten, 5000);
  host_to_node_real(MA_flatten, 5000);
  host_to_node_real(jtot_flatten, 5000);
  host_to_node_real(aoa_flatten, 400000);
  host_to_node_real(coeff_flatten, 400000);
  host_to_node_int(itot,100);
  host_to_node_boolean_1(rho_const);
  host_to_node_string(check_name_flatten,3000);
  host_to_node_string(clorcd_flatten,50000);

  #if !RP_HOST

    for (int k = 0; k < ktot; k++) {
      for (int i=0; i<itot[k]; i++) {
        RE[k][i] = RE_flatten[k*50 + i];
        MA[k][i] = MA_flatten[k*50 + i] ;
        jtot[k][i] = jtot_flatten[k*50 + i];
        for (int j=0; j<jtot[k][i]; j++) {
          aoa[k][i][j] = aoa_flatten[k*50 + i*80 + j];
          coeff[k][i][j] = coeff_flatten[k*50 + i*80 + j];
        }
      }
    }    

    for(int k = 0; k < ktot; k++){
      strncpy(check_name[k], &check_name_flatten[k*30], 30);

      for(int i = 0; i < itot[k]; i++){
        strncpy(clorcd[k][i], &clorcd_flatten[k*500 + i*10], 10);
      }
    }

  #endif

}

/*---------------------------------------------------------------------------*/
DEFINE_ADJUST(my_SrcComp, domain)
/*---------------------------------------------------------------------------*/
/* Main computation routine. Called each iteration. Loops over all rotor     */
/* zones, computes aerodynamic forces per cell, applies source terms to the  */
/* momentum equations, and runs the trimming algorithm if enabled.           */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST /* This portion of the code is only run on nodes*/
    int i, a, count;
    Thread *tc;
    cell_t cc;
    real *factor_ptr, sum_area, sum_sum_area;
    real dpitc, dpits; /* Sine and cosine rotor disk pitch angle components in 
      degrees. */
    real dbanc, dbans; /* Sine and cosine rotor disk bank angle components in 
      degrees. */
    real r_pb, r_pb2; /* Radial distance angle in the pitch/bank plane. */
    real psi_pb, psi_pbs, psi_pbc; /* Psi angle in the pitch/bank plane and 
      their sine and cosine components. */
    real Usc, Utc, Unc;  
    real aeps;
    real Utotal;
    real beta, beta_pb, betac, betas; 
    real theta, chord;
    real CL, CD;
    real Ftc, Fnc, Fsc;
    real Fx, Fy, Fz;
    real *xsource_ptr, *ysource_ptr, *zsource_ptr;
    real *xsource_old_ptr, *ysource_old_ptr, *zsource_old_ptr;
    real Re_t, Ma_t, alpha_t;
    real thrust, Fz_pb;
    real x_pb, y_pb, z_pb, Mx_pb, My_pb;
    int flag_jac;
    real ct_s, ct_h, ct_l, bcop_l, bcop_s, bcop_h;
    real cmx_lc, cmx_ls, cmx_s, cmx_hc, cmx_hs;
    real cmy_lc, cmy_ls, cmy_s, cmy_hc, cmy_hs;
    real bcyc_l, bcyc_s, bcyc_h, bcys_l, bcys_s, bcys_h; 
    real ct_lc, ct_ls, ct_hc, ct_hs, cmx_l, cmx_h, cmy_l, cmy_h; /*Low and high
      thrust and moment coeffs. */
    real alpha_tmax, alpha_tmin; /* Max and min angle of attack */
    real torque, Ft_pb; /* Torque force component */
    real sum_Mx_pb, sum_My_pb; /* Sum of moments about pb */
    real rho_ref; /* Reference density */
    real coeff_denom; /* Coefficient denominator */
    FILE *fp; /* File pointer for output */
    int sum_count = 0; /* Counter for active cells */
    int istr_target = 0; /* Flag to indicate if output data has trim target. */
    int isunst = 0; /* Flag for unsteady simulation*/
    real gwf; /* Gaussian weighting factor. */
    int local_indx; /* Local index for cell in rotor zone. */
    int rotor_cell; /* Flag to indicate if the cell is within rotor 
      influence. */
    real r_total; /* Radial distance from rotor center to cell center. */
    real xx, yy, zz; /* Coordinates of the cell center. */
    int gauss_active; /* Flag to indicate if Gaussian weighting is active in 
      the cell. */
    real sample_vol; /* Volume of the sample cell. */
    real sum_sample_vol; /* Summation of sample volumes across all active 
      cells. */
    real sum_gwf; /* Summation of the Gaussian weighting fn across all active 
      cells. */
    int vol_sample_entries; /* Number of volume samples. */
    real deltabeta; /* Vertical difference between blade and cell (rad). */
    int cellsaboverotor; /* Flag to indicate if cell is above rotor. */
    int cells_ahead_behind_rotor; /* Flag to indicate if cell is ahead of rotor. */
    real gauss_denom; /* Current value of Gaussian denominator. */
    real sigma_az; /* Azimuthal sigma value for Gaussian fn. */
    real sigma_vert; /* Vertical sigma value for Gaussian fn. */
    real saz_fac; /* Sigma azimuthal multiplication factor. */
    real sve_fac; /* Sigma vertical multiplication factor. */
    real gde_fac; /* Gaussian denom multiplication factor. */
    real deltapsi; /* Angular difference between blade and cell (rad). */
    real az_dist; /* Azimuthal distance from rotor blade (m). */
    real vert_dist; /* Vertical distance from blade (m). */

    /* If in debug mode, print to console where the code currently is at. */
    if (dbg == 1){
      for (i = 0; i < nrtz; i++) {
        Message("In DEFINE_ADJUST, rotor %d \n", i+1);
      }
    }

    /* If Fluent is set to transient, the unsteady rotor code (ALM) activates. 
       This will be printed to console if in debug mode. */
    if (CURRENT_TIME > 0.0){
      isunst = 1;
      /* hdr is only true if run for the first time. */
      if (hdr == 1)
        Message0("VBM in unsteady mode. \n");
    }
    else{
      isunst = 0;
      if (hdr == 1)
        Message0("VBM in steady mode. \n");
    }

    /* Upon first run after creating/loading, create header lines to the output
       files. */
    if ((hdr == 1 || N_ITER == 0) && myid == 0) {
      for (i = 0; i < nrtz; i++) {
        for (int h = 0; h < (sizeof(data_fname) / sizeof(data_fname[0])); h++){

          /* Only some output variables have trim targets. The ones that do will
             have their targets printed to output file. */
          if (h > 5) {
            istr_target = 1;
          }
          if (h < 5) {
            istr_target = 0;
          }

          sprintf(newfilename, "rotor_%d_%s", i + 1, data_fname[h]);

          fp = fopen(newfilename, "a");
          if (fp == NULL) {
            Message0("Error: Unable to open file for writing data.\n");
          }
          else {
            fprintf(fp, "\n BEGIN FILE");
            if (istr_target == 1)
              fprintf(fp, "\niteration\tdata\ttarget");
            if (istr_target == 0)
              fprintf(fp, "\niteration\tdata");
            
            fclose(fp);
          }
        }
      }
    }
    
    hdr = 0;

    if (domain == NULL) {
      Message0("Error: Domain is NULL.\n");
      return; /* Exit the function to prevent further execution */
    }

    /*Loop through all zones*/
    i = 0;
    while (i < nrtz) {
      if (dbg == 1){
        Message0("%d Rotor %d of %d \n",myid, i+1, nrtz);
      }
      
      tc = Lookup_Thread(domain, czon[i]);
      if (tc == NULL) {
        Message0("Error: Cell thread is NULL in define_adjust.\n");
        return; 
      }

      factor_ptr = Get_Thread_Memory(tc, factor[i]);
      if (factor_ptr == NULL) {
        Message0("Error: factor_ptr is NULL in define_adjust.\n");
        return; 
      }
      xsource_ptr = Get_Thread_Memory(tc, xsource[i]);
      if (xsource_ptr == NULL) {
        Message0("Error: xsource_ptr %d is NULL in define_adjust.\n", i);
      }
      ysource_ptr = Get_Thread_Memory(tc, ysource[i]);
      if (ysource_ptr == NULL) {
        Message0("Error: ysource_ptr %d is NULL in define_adjust.\n", i);
      }
      zsource_ptr = Get_Thread_Memory(tc, zsource[i]);
      if (zsource_ptr == NULL) {
        Message0("Error: zsource_ptr %d is NULL in define_adjust.\n", i);
      }

      zsource_old_ptr = Get_Thread_Memory(tc, zsource_old[i]);
      if (zsource_old_ptr == NULL) {
        Message0("Error: zsource_old_ptr %d is NULL in define_adjust.\n", i);
      }
      ysource_old_ptr = Get_Thread_Memory(tc, ysource_old[i]);
      if (ysource_old_ptr == NULL) {
        Message0("Error: ysource_old_ptr %d is NULL in define_adjust.\n", i);
      }
      xsource_old_ptr = Get_Thread_Memory(tc, xsource_old[i]);
      if (xsource_old_ptr == NULL) {
        Message0("Error: xsource_old_ptr %d is NULL in define_adjust.\n", i);
      }

      rho_ref = RP_Get_Real("reference-density");

      coeff_denom = 0.5 * rho_ref * rrl[i] * rrl[i] * M_PI * rspe[i] * 
        rspe[i] * rrl[i] * rrl[i];

      flag_jac=0;
      LABEL1:
      thrust = 0.0;
      torque = 0.0;
      Mx_pb = 0.0;
      My_pb = 0.0;
      sum_area = 0.0;
      sum_sum_area = 0.0;
      alpha_tmax = -1000.0;
      alpha_tmin = 1000.0;
      local_indx = 0;
      count = 0;
      gauss_active = 0;
      sum_sample_vol = 0.0;
      sum_gwf = 0.0;
      vol_sample_entries = 0;
      sample_vol = 0.0;
      deltabeta = 0.0;
      cells_ahead_behind_rotor = 0;
      cellsaboverotor = 0;

      begin_c_loop_int(cc,tc)

        rotor_cell = 1;
        gwf = 0.0;

        my_get_radial_position(i, cc, tc, count, 
          &r_pb, &psi_pb, &r_pb2,
          &x_pb, &y_pb, &z_pb, &rotor_cell, 
          &r_total, &xx, &yy, &zz, &beta_pb,
          &dpitc, &dpits, &dbanc, &dbans);

        /* If cell outside rotor, move to another cell. */
        if(rotor_cell == 0){
          continue;
        }

        /* Calculate near-rotor cell volume, and summate. */
        if ((count < cell_volume_samples) && (sample_vol <= 1e-12)){
          sample_vol = C_VOLUME(cc,tc);
          sum_sample_vol += sample_vol;
          vol_sample_entries++;
        }

        my_vel_xyz_to_blade(i, cc, tc, count, 
          psi_pb, r_pb2, 
          &Usc, &Utc, &Unc,
          &aeps, &Utotal, &beta, 
          dpitc, dpits, dbanc, dbans,
          &psi_pbs, &psi_pbc, &betac, &betas);

        my_get_pitch(i, cc, tc, count, r_pb, 
          psi_pb, &theta, psi_pbs, psi_pbc);

        my_get_chord_alm_area(i, cc, tc, count, r_pb, 
          &chord);

        my_get_re(i, cc, tc, count, chord, 
          Utotal, &Re_t);
        /*Message0("   Re_t %f \n",Re_t);*/

        my_get_ma(i, cc, tc, count, Utotal, 
          &Ma_t);
        /*Message0("   Ma_t %f \n",Ma_t);*/

        my_get_aoa(i, cc, tc, count, aeps, 
          theta, &alpha_t);
        /*Message0("   alpha_t %f \n",alpha_t*180/M_PI);*/        

        my_get_cl(i, cc, tc, count, r_pb, 
          Ma_t, Re_t, alpha_t, &CL);
        /*Message0("   CL %f \n",CL);*/   

        my_get_cd(i, cc, tc, count, r_pb, 
          Ma_t, Re_t, alpha_t, &CD);
        /*Message0("   CD %f \n",CD);*/

        my_get_force(i, cc, tc, factor_ptr, 
          &sum_area, count, CL, CD, 
          chord, aeps, Utotal, r_pb, 
          &Ftc, &Fnc, &Fsc, isunst, local_indx, 
          sample_vol);
        /*Message0("Ftc=%f Fnc=%f Fsc=%f  \n",Ftc, Fnc, Fsc);*/

        my_force_blade_to_xyz(i, cc, tc, count, 
          Ftc, Fnc, Fsc, beta, psi_pb,
          &Fx, &Fy, &Fz,
          &Fz_pb, &Ft_pb, chord, r_pb, 
          &gwf, isunst, &rotor_cell,
          beta_pb, &gauss_active, &sum_gwf,
          dpitc, dpits, dbanc, dbans,
          psi_pbs, psi_pbc, betas, betac, &deltabeta,
          &cellsaboverotor, &cells_ahead_behind_rotor,
          &sigma_az, &sigma_vert, &gauss_denom,
          &saz_fac, &sve_fac, &gde_fac, &deltapsi,
          &az_dist, &vert_dist); 

        /* If cell outside rotor, move to another cell. 
        if(rotor_cell == 0){
          continue;
        }*/
        
        if (flag_jac == 0) 
          my_get_source_terms(i, cc, tc, count, 
            xsource_ptr, ysource_ptr, zsource_ptr,
            xsource_old_ptr, ysource_old_ptr, zsource_old_ptr,
            Fx, Fy, Fz, urf_source, local_indx);

        thrust=thrust+Fz_pb;

        Mx_pb = Mx_pb + Fz_pb*y_pb;
        My_pb = My_pb + Fz_pb*x_pb;
        torque = torque + Ft_pb*sqrt(x_pb*x_pb + y_pb*y_pb);
        
        C_UDMI(cc, tc, 0) = r_pb;
        C_UDMI(cc, tc, 1) = psi_pb;
        C_UDMI(cc, tc, 2) = theta * (180.0 / M_PI);
        C_UDMI(cc, tc, 3) = Utotal;
        C_UDMI(cc, tc, 4) = chord;
        C_UDMI(cc, tc, 5) = rotor_cell;
        C_UDMI(cc, tc, 6) = sigma_az;
        C_UDMI(cc, tc, 7) = alpha_t* (180.0 / M_PI);
        C_UDMI(cc, tc, 8) = deltapsi;
        C_UDMI(cc, tc, 9) = az_dist;
        C_UDMI(cc, tc, 10) = gwf;
        C_UDMI(cc, tc, 11) = Fz;
        C_UDMI(cc, tc, 12) = local_indx; /* DO NOT ALTER THIS */

        alpha_tmax = fmax(alpha_tmax, alpha_t);
        alpha_tmin = fmin(alpha_tmin, alpha_t);
      
        local_indx++;
        count++;
      end_c_loop_int(cc,tc)
      
      if(isunst)
        my_gauss_opt(i, sum_sample_vol, sum_gwf, gauss_active, 
          vol_sample_entries, cellsaboverotor, cells_ahead_behind_rotor,
          &gauss_denom, &sigma_vert, &sigma_az,
          &saz_fac, &sve_fac, &gde_fac);

      /*Message("    sum_area %f \n",sum_area);*/
      sum_thrust[i] = PRF_GRSUM1(thrust);

      sum_Mx_pb = PRF_GRSUM1(Mx_pb);
      sum_My_pb = PRF_GRSUM1(My_pb);
      sum_torque[i] = PRF_GRSUM1(torque);

      /*get_cone(i, sum_thrust);*/

      my_get_ct(i, sum_thrust[i], coeff_denom);
      my_get_cmx(i, sum_Mx_pb, coeff_denom);
      my_get_cmy(i, sum_My_pb, coeff_denom);

      if (flag_jac == 0){
        power[i]=sum_torque[i]*rspe[i];

        if (dbg == 1) {
          Message0("Rotor %d Thrust %f ct= %f \n",i+1,sum_thrust[i],CT[i]);
          /*Message("    %d Rotor %d Thrust on node %f \n",myid,i+1,thrust);*/
          Message0("Rotor %d Torque %f \n",i+1,sum_torque[i]);
          /*Message("    %d Rotor %d Torque on node %f \n",myid,i+1,torque);*/
          Message0("Rotor %d Power %f \n",i+1,power[i]);
          /*Message0("    %d Rotor %d Mx_pb %f cmx= %f \n",myid,i+1,sum_Mx_pb,
          CMX); */
          /*Message("    %d Rotor %d Mx_pb on node %f \n",myid,i+1,Mx_pb);*/
          /*Message0("    %d Rotor %d My_pb %f cmy= %f \n",myid,i+1,sum_My_pb,
          CMY); */
          /*Message("    %d Rotor %d My_pb on node %f \n",myid,i+1,My_pb);*/
        }
      }

      my_update_trimming(i, trufq_co, trufq_cy, 
        trdf_co, trdf_cy,
        up_co, up_cy, 
        &flag_jac,
        &ct_l, &ct_s, &ct_h,
        &bcop_l, &bcop_s, &bcop_h,
        &cmx_lc, &cmx_ls, &cmx_s, 
        &cmx_hc, &cmx_hs,
        &cmy_lc, &cmy_ls, &cmy_s, 
        &cmy_hc, &cmy_hs,
        &bcyc_l, &bcyc_s, &bcyc_h,
        &bcys_l, &bcys_s, &bcys_h,
        &ct_lc, &ct_ls, 
        &ct_hc, &ct_hs,
        &cmx_l, &cmx_h, 
        &cmy_l, &cmy_h,
        sum_thrust[i]);

      if (flag_jac == 1 || flag_jac == -1 ||
        flag_jac == 2 || flag_jac == -2 ||
        flag_jac == 3 || flag_jac == -3) 
        goto LABEL1;

      i += 1;
    }
    istflag=0;
    
  #endif

}



/*---------------------------------------------------------------------------*/
/* MOMENTUM SOURCE TERM UDFs                                                 */
/*                                                                           */
/* Each rotor zone requires three momentum source UDFs (x, y, z). All       */
/* twelve functions share the same structure; only the rotor index and the   */
/* source array differ. The helper macro ROTOR_SRC_BODY reduces repetition. */
/*---------------------------------------------------------------------------*/
#if !RP_HOST
#define ROTOR_SRC_BODY(src_array, rotor_idx, fn_name)                         \
  do {                                                                        \
    real *src_ptr;                                                            \
    real source   = 0.0;                                                      \
    int local_indx = (int)C_UDMI(cc, tc, 12);                                \
    src_ptr = Get_Thread_Memory(tc, (src_array)[rotor_idx]);                  \
    if (src_ptr == NULL)                                                      \
      Message0("Error: source ptr NULL in DEFINE_SOURCE(" #fn_name ")\n");    \
    else                                                                      \
      source = src_ptr[local_indx];                                           \
    return source;                                                            \
  } while (0)
#endif

/* Rotor 1 (index 0) */
DEFINE_SOURCE(my_xmom_src_1, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(xsource, 0, my_xmom_src_1); 
        #endif 
    }
        
DEFINE_SOURCE(my_ymom_src_1, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(ysource, 0, my_ymom_src_1); 
        #endif 
    }

DEFINE_SOURCE(my_zmom_src_1, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(zsource, 0, my_zmom_src_1); 
        #endif 
    }

/* Rotor 2 (index 1) */
DEFINE_SOURCE(my_xmom_src_2, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(xsource, 1, my_xmom_src_2); 
        #endif 
    }
DEFINE_SOURCE(my_ymom_src_2, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(ysource, 1, my_ymom_src_2); 
        #endif 
    }
DEFINE_SOURCE(my_zmom_src_2, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(zsource, 1, my_zmom_src_2); 
        #endif 
    }

/* Rotor 3 (index 2) */
DEFINE_SOURCE(my_xmom_src_3, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(xsource, 2, my_xmom_src_3); 
        #endif 
    }
DEFINE_SOURCE(my_ymom_src_3, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(ysource, 2, my_ymom_src_3); 
        #endif 
    }
DEFINE_SOURCE(my_zmom_src_3, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(zsource, 2, my_zmom_src_3); 
        #endif 
    }

/* Rotor 4 (index 3) */
DEFINE_SOURCE(my_xmom_src_4, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(xsource, 3, my_xmom_src_4); 
        #endif 
    }
DEFINE_SOURCE(my_ymom_src_4, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(ysource, 3, my_ymom_src_4); 
        #endif 
    }
DEFINE_SOURCE(my_zmom_src_4, cc, tc, dS, eqn) 
    { 
        #if !RP_HOST 
            ROTOR_SRC_BODY(zsource, 3, my_zmom_src_4); 
        #endif 
    }

/*---------------------------------------------------------------------------*/
DEFINE_EXECUTE_AT_END(my_write_rotor_data_to_file)
/*---------------------------------------------------------------------------*/
/* Appends per-rotor performance data (pitch angles, thrust, torque, power,  */
/* and force/moment coefficients) to output files at the end of each         */
/* iteration. Angle quantities are converted from radians to degrees.        */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    if (myid == 0) {
      FILE *fp = NULL;
      int i, h;
      int num_outputs = (int)(sizeof(data_fname) / sizeof(data_fname[0]));
      real data_save[9], tr_target[3];

      for (i = 0; i < nrtz; i++) {
        data_save[0] = bcop[i];
        data_save[1] = bcyc[i];
        data_save[2] = bcys[i];
        data_save[3] = sum_thrust[i];
        data_save[4] = sum_torque[i];
        data_save[5] = power[i];
        data_save[6] = CT[i];
        data_save[7] = CMX[i];
        data_save[8] = CMY[i];

        tr_target[0] = ctd[i];
        tr_target[1] = cmxd[i];
        tr_target[2] = cmyd[i];

        for (h = 0; h < num_outputs; h++) {
          int isangle     = (h < 3) ? 1 : 0;
          int istr_target = (h > 5) ? 1 : 0;

          sprintf(newfilename, "rotor_%d_%s", i + 1, data_fname[h]);
          fp = fopen(newfilename, "a");
          if (fp == NULL) {
            Message0("Error: Unable to open %s for writing.\n", newfilename);
            continue;
          }

          if (isangle) {
            fprintf(fp, "\n%d\t%f", N_ITER, data_save[h] * 180.0 / M_PI);
            if (dbg == 1)
              Message0("Writing Rotor %d %s = %f deg to %s\n",
                i + 1, data_vname[h], data_save[h] * 180.0 / M_PI, newfilename);
          }
          else if (istr_target) {
            fprintf(fp, "\n%d\t%f\t%f", N_ITER, data_save[h], tr_target[h - 6]);
            if (dbg == 1)
              Message0("Writing Rotor %d %s = %f (target %f) to %s\n",
                i + 1, data_vname[h], data_save[h], tr_target[h - 6], newfilename);
          }
          else {
            fprintf(fp, "\n%d\t%f", N_ITER, data_save[h]);
            if (dbg == 1)
              Message0("Writing Rotor %d %s = %f to %s\n",
                i + 1, data_vname[h], data_save[h], newfilename);
          }

          fclose(fp);
          if (dbg == 1)
            Message0("%s written to %s\n", data_vname[h], newfilename);
        }
      }
    }
    trnw = 0;
  #endif
}
