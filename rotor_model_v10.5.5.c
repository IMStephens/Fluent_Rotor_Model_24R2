#include "udf.h"
#include "seem.h"
#include "global.h"
#include "mem.h"
#include "thread_mem.h"
#include <stdio.h>
#include <time.h>
#include <string.h>
#include <math.h>

#define NUM_UDM 11
#define dbg 0 /*Set to 1 for debugging. Will print more messages.*/
#define my_mem_loop(m)for(m = my_mem; m != NULL; m = m->next)


/*GUI start*/
#if !RP_NODE
static Pointer
List_Ref(Pointer l, int n){
  register int i;
  for (i=0; i<n; ++i){
    l = CDR(l);
  }
  return CAR(l);
}
#endif
/*GUI end*/


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
int ac_mass = 7000;                   /*Aircraft mass for force balancing, kg*/
int nrtz;                             /*Number of Rotor zones*/
int nbld[10];                         /*Number of blades*/
real rrl[10];                         /*Rotor radius*/
real rspe[10];                        /*Rotor speed*/
real teff[10];                        /*Tip effect*/
real dskco[10][3];                    /*Rotor disk origin 0x 1y 2z*/
real dpit[10];                        /*Rotor disk pitch angle*/
real tdpit[10];                       /*Rotor disk trim pitch angle*/
real dban[10];                        /*Rotor disk bank angle*/
real tdban[10];                       /*Rotor disk trim bank angle*/
int fzon[10];                         /*Rotor face zone ID*/

real bcop[10];                        /*Blade pitch collective*/
real bcys[10];                        /*Blade pitch cyclic sin*/
real bcyc[10];                        /*Blade pitch cyclic cos*/
real bflco[10];                       /*Blade flapping cone*/
real bfls[10];                        /*Blade flapping cyclic sin*/
real bflc[10];                        /*Blade flapping cyclic cos*/

int nsec[10];                         /*Number of blade sections along span*/
real rin[10][20];                     /*Normalized inner radius*/
real rout[10][20];                    /*Normalized outer radius*/
real rsec[10][20];				            /*Normalizes radius from GUI*/
real cin[10][20];                     /*Chord at normalized inner radius*/
real cout[10][20];                    /*Chord at normalized outer radius*/
real csec[10][20];				            /*Chord at normalized sectional radius (GUI)*/
real twst[10][20];                    /*Twist*/
char type[10][20][30];                /*Profile type name */

int trmco[10]={0,0,0,0,0,0,0,0,0,0};  /*Trimming collective pitch 1(ON) 0(OFF)*/
int trmcy[10]={0,0,0,0,0,0,0,0,0,0};  /*Trimming cyclic pitch 1(ON) 0(OFF)*/
int forba[10]={0,0,0,0,0,0,0,0,0,0};  /*Force balancing 1(ON) 0(OFF)*/
int momba[10]={0,0,0,0,0,0,0,0,0,0};  /*Moment balancing 1(ON) 0(OFF)*/
int fuseid[10];					              /*Fuselage face thread ID*/
real cmgx[10];                        /*x-location of fuselage CG*/
real cmgy[10];					              /*y-location of fuselage CG*/
real cmgz[10];                        /*z-location of fuselage CG*/
int errchk;                           /*Error check for NaN values*/    
  
real trdf[10];                        /*Trimming damping factor*/
int trufq[10];                        /*Trimming update frequency*/
real ctd[10];                         /*Trimming desired ct*/
real cmxd[10];                        /*Trimming desired cmx*/
real cmyd[10];                        /*Trimming desired cmy*/
int trnw = 0;                         /*Currently trimming*/

real CT[10];                          /*Rotor thrust coeff*/
real CMX[10];                         /*Rotor moment coeff x*/
real CMY[10];                         /*Rotor moment coeff y*/

real sum_thrust[10];                  /*Summation of Rotor Thrust*/
real sum_torque[10];                  /*Summation of Rotor Torque*/
real power[10];                       /*Summation of Rotor Power*/

real force[ND_ND];				            /*Force on thread obj*/
real moment[ND_ND];				            /*Moment on thread obj*/
real cg[ND_ND];					              /*cg on thread obj*/
real main_rot_thrust_target_init;     /*Initial main rotor thrust coeff target*/

int ktot;							                /*Number of different airfoiltables read in*/
char file_name[100][30];			        /*Name of airfoil tables, eg. airfoil.dat*/
char check_name[100][30];			        /*Name of airfoil, eg. naca0012*/
char clorcd[100][50][10];			        /*cl or cd entry*/
float RE[100][50];				            /*Airfoil Reynolds number*/
float MA[100][50];				            /*Airfoil Mach number*/
int itot[100];					              /*Number of entries per .dat file*/
int jtot[100][50];				            /*Number of aoa entries per set*/
float aoa[100][50][80];			          /*Angle of attack*/
float coeff[100][50][80];			        /*cl or cd*/
cxboolean rho_const;                  /*true if density const, false otherwise*/

/*TUI*/
real urf_source=0.5;				          /*URF for source calculation*/
int i_start=10;					              /*No iterations for linear rspe ramp up*/
real coef_fac=0.5;				            /*Factor to calculate the thrust, momentum,
									                      and torque coefficients.
									                      coef_fac=1.0   US
									                      coef_fac=0.5   EU*/
real dalpha=10.0*M_PI/180.0;		      /*Trimming: perturbation angle for Jacobian
                                        in deg*/
real limiter=0.1;				              /*Max. angle change relative to dalpha*/

char newfilename[256];                /*Filename for reading and writing rotor data*/
char time_str[100];                   /*String for time in filename*/  
/*time_t now;                           Current time for filename*/ 
int hdr = 1;                          /*Header for writing data file, 1=on, 0=off*/ 

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
int istflag=1;					/*flag 1 at start
								       0 after first ADJUST run              */
/*---------------------------------------------------------------------------*/
/*                 BEGIN GLOBAL TRIMMING VARIABLES                           */
/*---------------------------------------------------------------------------*/
int trufq_co[10], trufq_cy[10];
real trdf_co[10], trdf_cy[10];
int up_co[10], up_cy[10];
real sldity[10];
/*---------------------------------------------------------------------------*/
/*                   END GLOBAL TRIMMING VARIABLES                           */
/*---------------------------------------------------------------------------*/

/*---------------------------------------------------------------------------*/
/*                 BEGIN GLOBAL MEMORY                                       */
/*---------------------------------------------------------------------------*/
#if !RP_HOST
  int czon[10];                      /*Rotor cell zone ID*/
  char factor[10][30];               /*Memory: Geometry weighting factor*/
  char xsource[10][30];				/*Memory: Source x-component*/
  char ysource[10][30];				/*Memory: Source y-component*/
  char zsource[10][30];				/*Memory: Source z-component*/

  char xsource_old[10][30];	/*Memory: Source x-component previous timestep UR*/
  char ysource_old[10][30];	/*Memory: Source y-component previous timestep UR*/
  char zsource_old[10][30];	/*Memory: Source z-component previous timestep UR*/
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
  void my_get_solidity(real sldity_ptr[])
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    int i,j;
    real area, area_tot;

    /*Loop through all zones*/
    i=0;
    while (i < nrtz) {
      area_tot=0.0;
	    j=0;
	    while (j < nsec[i]) {
        area=0.5*(cout[i][j]+cin[i][j])*(rout[i][j]-rin[i][j])*rrl[i];
        area_tot=area_tot + area;
	      j += 1;
	    }

      /*Below is an error check to ensure no division by zero*/
      if ((M_PI*rrl[i]*rrl[i]) == 0.0) {
        Message("Error: Rotor %d radius is undefined.\n",i+1);
        sldity_ptr[i] = 0.0;
      } 
      else {
        sldity_ptr[i]=nbld[i]*area_tot/(M_PI*rrl[i]*rrl[i]);
      }
	    
	    i += 1;
    }
  }


/*---------------------------------------------------------------------------*/
  void my_obtain_cell_id()
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    int i, icell_f,icell_c0, icell_c1;
    Thread *tff = NULL;
    Thread *tc0 = NULL, *tc1 = NULL;
    face_t ff;
    cell_t cc;
    int sum_icell_f,sum_icell_c0, sum_icell_c1;

    Domain *domain;
    domain = Get_Domain(1);

    if (domain == NULL) {
      Message0("Error: Domain is NULL.\n");
      return; // Exit the function to prevent further execution
    }

    if (dbg == 1)
      Message0("%d obtain_cell_id \n",myid);
    
    i=0;
    while (i < nrtz){
      tff = Lookup_Thread(domain,fzon[i]);
      if (tff == NULL) {
        Message("Error: tff is NULL in my_obtain_cell_id\n");
        return; /* Exit function to prevent further execution */
      }

      icell_f=0;
      begin_f_loop(ff,tff)
        if PRINCIPAL_FACE_P(ff,tff) icell_f += 1;
      end_f_loop(ff,tff)

      tc0=F_C0_THREAD(fzon[i],tff);
      tc1=F_C1_THREAD(fzon[i],tff);

      icell_c0=0;
      begin_c_loop_int(cc,tc0)
        icell_c0 += 1;
      end_c_loop_int(cc,tc0)
  
      icell_c1=0;
      begin_c_loop_int(cc,tc1)
        icell_c1 += 1;
      end_c_loop_int(cc,tc1)
  
      sum_icell_f = PRF_GISUM1(icell_f);
      sum_icell_c0 = PRF_GISUM1(icell_c0);
      sum_icell_c1 = PRF_GISUM1(icell_c1);

      if (dbg == 1){
        Message0("  %d sum_icell_f = %d, sum_icell_c0 = %d, "
          "sum_icell_c1 = %d \n",
          myid, sum_icell_f, sum_icell_c0, sum_icell_c1);
      }

      if (sum_icell_c0 == sum_icell_f) {
        czon[i]=THREAD_ID(tc0);
        if (dbg == 1)
          Message0("Rotor %d: czon %d \n", i+1, czon[i]);
      }
      else if (sum_icell_c1 == sum_icell_f) {
        czon[i]=THREAD_ID(tc1);
        if (dbg == 1)
          Message0("Rotor %d: czon %d \n", i+1, czon[i]);
      }
      else{
        if (dbg == 1)
          Message0("ERROR: Cell ID not determined \n");
      }
      i += 1;
    }
  }


/*---------------------------------------------------------------------------*/
  void my_allocate_memory(char factor[][30], char xsource[][30],
	  char ysource[][30], char zsource[][30], char xsource_old[][30],
    char ysource_old[][30], char zsource_old[][30])
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
  
    int i, count, sum_count;
    Thread *tc = NULL;
    Domain *domain;
    cell_t cc;
    real *xsource_old_ptr,*ysource_old_ptr,*zsource_old_ptr;
    
    domain = Get_Domain(1);
    
    if (domain == NULL) {
      Message0("Error: Domain is NULL.\n");
      return;
    }
  
    if (dbg == 1)
      Message0("%d allocate_memory \n",myid);
  
    (void)strcpy(factor[0],"factor0");
    (void)strcpy(factor[1],"factor1");
    (void)strcpy(factor[2],"factor2");
    (void)strcpy(factor[3],"factor3");
    (void)strcpy(factor[4],"factor4");
    (void)strcpy(factor[5],"factor5");
    (void)strcpy(factor[6],"factor6");
    (void)strcpy(factor[7],"factor7");
    (void)strcpy(factor[8],"factor8");
    (void)strcpy(factor[9],"factor9");
  
    (void)strcpy(xsource[0],"xsource0");
    (void)strcpy(xsource[1],"xsource1");
    (void)strcpy(xsource[2],"xsource2");
    (void)strcpy(xsource[3],"xsource3");
    (void)strcpy(xsource[4],"xsource4");
    (void)strcpy(xsource[5],"xsource5");
    (void)strcpy(xsource[6],"xsource6");
    (void)strcpy(xsource[7],"xsource7");
    (void)strcpy(xsource[8],"xsource8");
    (void)strcpy(xsource[9],"xsource9");
  
    (void)strcpy(ysource[0],"ysource0");
    (void)strcpy(ysource[1],"ysource1");
    (void)strcpy(ysource[2],"ysource2");
    (void)strcpy(ysource[3],"ysource3");
    (void)strcpy(ysource[4],"ysource4");
    (void)strcpy(ysource[5],"ysource5");
    (void)strcpy(ysource[6],"ysource6");
    (void)strcpy(ysource[7],"ysource7");
    (void)strcpy(ysource[8],"ysource8");
    (void)strcpy(ysource[9],"ysource9");
  
    (void)strcpy(zsource[0],"zsource0");
    (void)strcpy(zsource[1],"zsource1");
    (void)strcpy(zsource[2],"zsource2");
    (void)strcpy(zsource[3],"zsource3");
    (void)strcpy(zsource[4],"zsource4");
    (void)strcpy(zsource[5],"zsource5");
    (void)strcpy(zsource[6],"zsource6");
    (void)strcpy(zsource[7],"zsource7");
    (void)strcpy(zsource[8],"zsource8");
    (void)strcpy(zsource[9],"zsource9");
  
  
    (void)strcpy(xsource_old[0],"xsource_old0");
    (void)strcpy(xsource_old[1],"xsource_old1");
    (void)strcpy(xsource_old[2],"xsource_old2");
    (void)strcpy(xsource_old[3],"xsource_old3");
    (void)strcpy(xsource_old[4],"xsource_old4");
    (void)strcpy(xsource_old[5],"xsource_old5");
    (void)strcpy(xsource_old[6],"xsource_old6");
    (void)strcpy(xsource_old[7],"xsource_old7");
    (void)strcpy(xsource_old[8],"xsource_old8");
    (void)strcpy(xsource_old[9],"xsource_old9");
  
    (void)strcpy(ysource_old[0],"ysource_old0");
    (void)strcpy(ysource_old[1],"ysource_old1");
    (void)strcpy(ysource_old[2],"ysource_old2");
    (void)strcpy(ysource_old[3],"ysource_old3");
    (void)strcpy(ysource_old[4],"ysource_old4");
    (void)strcpy(ysource_old[5],"ysource_old5");
    (void)strcpy(ysource_old[6],"ysource_old6");
    (void)strcpy(ysource_old[7],"ysource_old7");
    (void)strcpy(ysource_old[8],"ysource_old8");
    (void)strcpy(ysource_old[9],"ysource_old9");
  
    (void)strcpy(zsource_old[0],"zsource_old0");
    (void)strcpy(zsource_old[1],"zsource_old1");
    (void)strcpy(zsource_old[2],"zsource_old2");
    (void)strcpy(zsource_old[3],"zsource_old3");
    (void)strcpy(zsource_old[4],"zsource_old4");
    (void)strcpy(zsource_old[5],"zsource_old5");
    (void)strcpy(zsource_old[6],"zsource_old6");
    (void)strcpy(zsource_old[7],"zsource_old7");
    (void)strcpy(zsource_old[8],"zsource_old8");
    (void)strcpy(zsource_old[9],"zsource_old9");
  
    i=0;
    while (i < nrtz) {
      tc = Lookup_Thread(domain, czon[i]);
  
      if (tc == NULL) {
        Message0("Error: tc is NULL in my_allocate_memory\n");
        return;
      }
  
      count = 0;
      begin_c_loop_int(cc,tc)
        count++;
      end_c_loop_int(cc,tc)

  	  sum_count = PRF_GISUM1(count);
  
      if (dbg == 1){
        Message0("  %d Rotor %d: Number of cells in volume=%d "
          "on this node =%d\n", myid, i+1, sum_count, count);
      }

      Alloc_Thread_Memory(tc,factor[i],sizeof(real)*count);
      Alloc_Thread_Memory(tc,xsource[i],sizeof(real)*count);
      Alloc_Thread_Memory(tc,ysource[i],sizeof(real)*count);
      Alloc_Thread_Memory(tc,zsource[i],sizeof(real)*count);
  
      Alloc_Thread_Memory(tc,xsource_old[i],sizeof(real)*count);
      Alloc_Thread_Memory(tc,ysource_old[i],sizeof(real)*count);
      Alloc_Thread_Memory(tc,zsource_old[i],sizeof(real)*count);
  
      i += 1;
    }
  
  
    /*Initialize old timestep for source terms*/
    /*if (istflag == 1)
      {*/
    i=0;
    while (i < nrtz) {
      tc = Lookup_Thread(domain, czon[i]);
  
      if (tc == NULL) {
        Message0("Error: tc is NULL in my_allocate_memory\n");
        return;
      }
  
      xsource_old_ptr= Get_Thread_Memory(tc,xsource_old[i]);
      if (xsource_old_ptr == NULL) {
        Message0("Error: xsource_old_ptr is NULL in my_allocate_memory\n");
        xsource_old_ptr = 0;
      }
  
      ysource_old_ptr= Get_Thread_Memory(tc,ysource_old[i]);
      if (ysource_old_ptr == NULL) {
        Message0("Error: ysource_old_ptr is NULL in my_allocate_memory\n");
        ysource_old_ptr = 0;
      }
      
      zsource_old_ptr= Get_Thread_Memory(tc,zsource_old[i]);
      if (zsource_old_ptr == NULL) {
        Message0("Error: zsource_old_ptr is NULL in my_allocate_memory\n");
        zsource_old_ptr = 0;
      }
  
      begin_c_loop_int(cc,tc)
        xsource_old_ptr[cc]=0.0;
        ysource_old_ptr[cc]=0.0;
        zsource_old_ptr[cc]=0.0;
      end_c_loop_int(cc,tc)
  
      i += 1;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_geom_factor()
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    int i, count, n;
    Thread *tc = NULL, *tf = NULL, *t = NULL;
    cell_t cc;
    face_t ff;
    real sum_area = 0, area = 0, A[ND_ND], xrc[ND_ND], rc;
    real *area_ptr;
    real sum_sum_area;
  
    Domain *domain;
    domain = Get_Domain(1);

    if (domain == NULL) {
      Message0("Error: Domain is NULL.\n");
      return;
    }
  
    if (dbg == 1)
      Message0("%d geom_factor \n",myid);
  
    i=0;
    while (i < nrtz){
      tc = Lookup_Thread(domain,czon[i]);
  
      if (tc == NULL) {
        Message0("Error: tc is NULL in my_geom_factor\n");
        return;
      }
  
      tf = Lookup_Thread(domain,fzon[i]);
      if (tf == NULL) {
        Message0("Error: tf is NULL in my_geom_factor\n");
        return;
      }
  
      if (dbg == 1)
        Message0("Rotor %d: czon %d fzon %d \n", i+1, czon[i], fzon[i]);
      
      area_ptr = Get_Thread_Memory(tc,factor[i]);
      if (area_ptr == NULL) {
        Message0("Error: area_ptr is NULL in my_geom_factor\n");
        return;
      }
  
      sum_area=0.0;
      begin_c_loop_int(cc,tc) {
  	    C_CENTROID(xrc,cc,tc);
        rc=sqrt( (xrc[0]-dskco[i][0])*(xrc[0]-dskco[i][0]) +
                   (xrc[1]-dskco[i][1])*(xrc[1]-dskco[i][1]) +
                   (xrc[2]-dskco[i][2])*(xrc[2]-dskco[i][2])
  	      	    );
  	    /*Message("rc %f",rc);*/
  
        count=0;
        c_face_loop(cc,tc,n) {
          t=C_FACE_THREAD(cc,tc,n);
          if(t == tf) {
            ff = C_FACE(cc,tc,n);
            F_AREA(A,ff,tf);
            area=NV_MAG(A);
  
            count=count+1;
          }
        }
  	    if(count == 0)
  	      Message0("ERROR: No area found");
  	    if((count =! 1))
  	      Message0("ERROR: Multiple areas found");
  
  	    sum_area=sum_area+area;
  
        /*Below is an error check to avoid division by zero*/
        if (2.0*M_PI*rc == 0.0) {
          Message("Error: Rotor %d origin is undefined.\n",i+1);
          area_ptr[cc] = 0.0;
        }
        else {
          area_ptr[cc]=((real) nbld[i])*area/(2.0*M_PI*rc);
        }
  	    
  	    /*area_ptr[cc]=area;*/
  	  }
      end_c_loop_int(cc,tc)
      sum_sum_area = PRF_GRSUM1(sum_area);
      if (dbg == 1)
        Message0("  %d Rotor %d: Face area =%f on this node=%f \n",
        myid, i+1, sum_sum_area, sum_area);
      
      i += 1;
    }
  }


/*---------------------------------------------------------------------------*/
  void my_start_trimming(int trufq_co[], int trufq_cy[],
                           real trdf_co[], real trdf_cy[],
                           int up_co[], int up_cy[])
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    int i;
  
    if (dbg == 1)
      Message0("%d start_trimming \n",myid);
  
   /*Set update frequencies and damping factors for cyclic and colective as
      same*/
    i=0;
    while (i < nrtz) {
      trufq_co[i]=trufq[i];
      trufq_cy[i]=trufq[i];
  
      trdf_co[i]=trdf[i];
      trdf_cy[i]=trdf[i];
  
      up_co[i]=trufq_co[i];
      up_cy[i]=trufq_cy[i];
  
      i += 1;
    }
  }

#endif

#if !RP_NODE
/*---------------------------------------------------------------------------*/
  void my_read_in_airfoil_tables()
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    FILE *f_ptr;
    int i, j, k;
  
    if (dbg == 1) {
      Message("\n read_in_airfoil_tables \n");
      Message("     -180 < AOA < 180 \n");
    }
  
    k=0;
    while (k < ktot) {
  	  f_ptr = fopen(file_name[k],"r");
      fscanf(f_ptr,"%s \n",check_name[k]);
      fscanf(f_ptr,"%d \n",&itot[k]);
  
      if(f_ptr == NULL) {
      	Message("\nNull file pointer\n");
      	continue;
  	  }
      i=0;
      while (i < itot[k])	{
        fscanf(f_ptr,"%s \n",clorcd[k][i]);
        fscanf(f_ptr,"%f \n",&RE[k][i]);
        fscanf(f_ptr,"%f \n",&MA[k][i]);
        fscanf(f_ptr,"%d \n",&jtot[k][i]);
  
        j=0;
        while (j < jtot[k][i]) {
          fscanf(f_ptr,"%f %f \n",&aoa[k][i][j], &coeff[k][i][j]);
          j+=1;
  	    }
  
  	    j=1;
        while (j < jtot[k][i]) {
          if(aoa[k][i][j-1] > aoa[k][i][j])
            Message("ERROR: Airfoil table entries must be ordered with "
              "increasing AOA \n");
  	      
          j+=1;
  	    }
  
  	    j=0;
        while (j < jtot[k][i]) {
          if (strcmp(clorcd[k][i],"cd") == 0)	{
  	        /*Message("  %s \n",clorcd[k][i]);*/
            if(coeff[k][i][j] < 0.0)
              Message("WARNING: Are you sure that cd is negative in %s, "
                "AOA %f ? \n", check_name[k], aoa[k][i][j]);
          }
  
          if (strcmp(clorcd[k][i],"cl") == 0) {
  	        /*Message("  %s \n",clorcd[k][i]);*/
  		      if(aoa[k][i][j] < -90.0 &&  coeff[k][i][j] < 0.0)
              Message("WARNING: Are you sure that cl is negative in %s, "
                "AOA %f ? \n", check_name[k], aoa[k][i][j]);
  
  		      if(aoa[k][i][j] > 90.0 &&  coeff[k][i][j] > 0.0)
              Message("WARNING: Are you sure that cl is positive in %s, "
                "AOA %f ? \n", check_name[k], aoa[k][i][j]);
  
            if(aoa[k][i][j] <= 90.0 && aoa[k][i][j] >= 0.0 && coeff[k][i][j] < 0.0)
              Message("WARNING: Are you sure that cl is negative in %s, AOA %f ? \n",
  	            check_name[k],aoa[k][i][j]);
  
            if(aoa[k][i][j] >= -90.0 && aoa[k][i][j] < 0.0 && coeff[k][i][j] > 0.0)
              Message("WARNING: Are you sure that cl is positive in %s, AOA %f ? \n",
  	            check_name[k],aoa[k][i][j]);
  		    }
  	      j+=1;
  	    }
  	    i+=1;
  	  }
      k+=1;
    }
    fclose(f_ptr);
  }
  
#endif

#if !RP_HOST
/*---------------------------------------------------------------------------*/
  void my_get_radial_position(int i, cell_t cc, Thread *tc, int count,
    real *r_pb_ptr, real *psi_pb_ptr, real *r_pb2_ptr, 
    real *x_pb_ptr, real *y_pb_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    real xrc[ND_ND], xx, yy, zz;
    real dpits, dpitc;
    real dbans, dbanc;
    real x_pb, y_pb, z_pb;
    real tmp1, tmp2;
  
    if(count == 0) {
      if (dbg == 1)
        Message0("  %d Rotor %d - get_radial_position \n",myid,i+1);
    }
  
    /*x y z coordinates*/
    C_CENTROID(xrc,cc,tc);
    xx = xrc[0]-dskco[i][0];
    yy = xrc[1]-dskco[i][1];
    zz = xrc[2]-dskco[i][2];
  
    /*Pitch "dpit[i]"*/
    dpits=sin(dpit[i]);
    dpitc=cos(dpit[i]);
  
    /*Bank "dban[i]"*/
    dbans=sin(dban[i]);
    dbanc=cos(dban[i]);
  
  
    x_pb= xx*dpitc                  - zz*dpits;
    *x_pb_ptr=x_pb;
    y_pb= xx*dpits*dbans + yy*dbanc + zz*dpitc*dbans;
    *y_pb_ptr=y_pb;
    z_pb= xx*dpits*dbanc - yy*dbans + zz*dpitc*dbanc;
  
    tmp1=sqrt (x_pb*x_pb + y_pb*y_pb + z_pb*z_pb);
  
    if (rrl[i] == 0.0) {
      Message("Error: Rotor %d radius is undefined.\n",i+1);
      tmp1 = 0.0;
    } 
    else {
      tmp1=tmp1/rrl[i];;
    }
    
    *r_pb_ptr=tmp1;
  
    *r_pb2_ptr=sqrt (x_pb*x_pb + y_pb*y_pb);
  
    tmp2=atan2(y_pb,x_pb);
    if(tmp2 < 0.0)
      tmp2=2.0*M_PI+tmp2;
    *psi_pb_ptr=tmp2;
    /*Now azimuthal angle [0,2pi]*/
    /*Message("  r_pb %f     psi_pb %f \n",r_pb_ptr, psi_pb_ptr*180/M_PI);*/
  
  }


/*---------------------------------------------------------------------------*/
  void my_vel_xyz_to_blade(int i, cell_t cc, Thread *tc, int count,
    real psi_pb, real r_pb2,
	  real *Usc_ptr, real *Utc_ptr, real *Unc_ptr,
	  real *aeps_ptr, real *Utotal_ptr, real *beta_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    real Ux, Uy, Uz;
    real dpits, dpitc;
    real dbans, dbanc;
    real Ux_pb, Uy_pb, Uz_pb;
    real psi_pbs, psi_pbc;
    real Ur_pb, Ut_pb;
    real beta, betas, betac;
    real tmp,tmp1,Utc,Unc,aeps;
    real Ut_omega, Ut_total, Utotal;
    int i_count;
    real tmp2,tmp3,tmp4;
  
    if(count == 0) {
      if (dbg == 1)
        Message0("  %d Rotor %d - vel_xyz_to_blade \n",myid,i+1);
    }
  
    /*Velocities in x y z coordinates*/
    Ux=C_U(cc,tc);
    Uy=C_V(cc,tc);
    Uz=C_W(cc,tc);
  
    /*Pitch "dpit[i]"*/
    dpits=sin(dpit[i]);
    dpitc=cos(dpit[i]);
  
    /*Bank "dban[i]"*/
    dbans=sin(dban[i]);
    dbanc=cos(dban[i]);
  
    Ux_pb= Ux*dpitc                  - Uz*dpits;
    Uy_pb= Ux*dpits*dbans + Uy*dbanc + Uz*dpitc*dbans;
    Uz_pb= Ux*dpits*dbanc - Uy*dbans + Uz*dpitc*dbanc;
  
    /*Polar coordinates in pitch/bank plane*/
    /*Message("   psi_pb %f \n", psi_pb*180/M_PI);*/
    psi_pbs=sin(psi_pb);
    psi_pbc=cos(psi_pb);
  
    Ur_pb= Ux_pb*psi_pbc  + Uy_pb*psi_pbs;
    Ut_pb= -Ux_pb*psi_pbs + Uy_pb*psi_pbc;
  
    /*Flapping and coning plane*/
    /*Message("   bflco %f bfls %f bflc %f \n", bflco[i], bfls[i], bflc[i]);*/
    beta= bflco[i] - bflc[i]*psi_pbc - bfls[i]*psi_pbs;
  
    if(isnan(beta)) {
      Message0("Error: Rotor %d beta returns NaN.\n",i+1);
    }
    
    *beta_ptr=beta;
    betas=sin(beta);
    betac=cos(beta);
    if(rspe[i] > 0.0)
  	  tmp1=-1.0;
    if(rspe[i] < 0.0)
  	  tmp1=1.0;
  
    tmp= Ur_pb*betac  +       + Uz_pb*betas;
    *Usc_ptr=tmp;
    Utc=              tmp1*Ut_pb;
    *Utc_ptr=Utc;
    Unc= -Ur_pb*betas         + Uz_pb*betac;
    *Unc_ptr=Unc;
  
    i_count=N_ITER;
    /*Message("i_start %d i_count %d \n",i_start,i_count);*/
    if(i_count >= i_start) {
  	  Ut_omega= fabs(rspe[i])*r_pb2;
  	  /*Message("i_count > i_start");*/
    }
    else {
      tmp2=i_count;
  	  tmp3=i_start;
  	  tmp4=tmp2/tmp3;
  	  Ut_omega= fabs(rspe[i])*r_pb2*tmp4;
      /*Message("i_count < i_start %f",tmp4);*/
    }
  
    Ut_total= Ut_omega+Utc;
  
    Utotal= sqrt(Ut_total*Ut_total + Unc*Unc);
    *Utotal_ptr=Utotal;
  
    /*Angle epsilon*/
    /*-pi < aeps < pi*/
    aeps=atan2(Unc,fabs(Ut_total));
    *aeps_ptr=aeps;
    if(aeps < -M_PI)
      *aeps_ptr=2.0*M_PI+aeps;
    if(aeps > M_PI)
      *aeps_ptr=aeps - 2.0*M_PI;
  
  }

/*---------------------------------------------------------------------------*/
  void my_get_pitch(int i, cell_t cc, Thread *tc, int count, real r_pb,
	  real psi_pb, real *theta_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    int j;
    real psi_pbs, psi_pbc;
    real theta, tmp;
  
    if(count == 0) {
      if (dbg == 1)
        Message0("  %d Rotor %d - get_pitch \n",myid,i+1);
    }
    psi_pbs=sin(psi_pb);
    psi_pbc=cos(psi_pb);
  
    /*Twist*/
    /*Assumed to be linear between blade section*/
    /*Loop through all blade sections*/
    j=0;
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
    /*Message("%f %f %f %f \n",rin[i][j],r_pb,rout[i][j],twst[i][j]*180/M_PI);*/
    /*Message("%f %f %f %f %f \n",rsec[i][j],r_pb,rsec[i][j+1],
  	                          twst[i][j]*180/M_PI,twst[i][j+1]*180/M_PI);*/
  
    if ((rout[i][j]-rin[i][j]) == 0) {
      /*Message("Error: Rotor %d radial section %d is undefined.\n",i+1,j);*/
      tmp = 0.0;
    }
    else {
      tmp=twst[i][j]+(twst[i][j+1]-twst[i][j])/(rout[i][j]-rin[i][j])*
      (r_pb-rin[i][j]);
    }
    
    /*Message("%f %f %f %f %f %f \n",rsec[i][j],r_pb,rsec[i][j+1],
  	                 twst[i][j]*180/M_PI,tmp*180/M_PI,twst[i][j+1]*180/M_PI);*/
  
    theta= bcop[i] + tmp - bcyc[i]*psi_pbc - bcys[i]*psi_pbs;
    *theta_ptr=theta;
    
    /*-pi < theta < pi*/
    if(theta < -M_PI)
      *theta_ptr=2.0*M_PI+theta;
    if(theta > M_PI)
      *theta_ptr=theta - 2.0*M_PI;
  }


/*---------------------------------------------------------------------------*/
  void my_get_chord(int i, cell_t cc, Thread *tc, int count, real r_pb,
	  real *chord_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    int j;
  
    if(count == 0) {
      if (dbg == 1)
        Message0("  %d Rotor %d - get_chord \n",myid,i+1);
    }
  
    /*Chord*/
    /*Assumed to be linear per blade section*/
    /*Loop through all blade sections*/
    j=0;
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
  
    if ((rout[i][j]-rin[i][j]) == 0.0) {
      /*Message("Error: Rotor %d chord at section %d is zero.\n", i+1,j);*/
      *chord_ptr = 0.0;
    } 
    else {
      *chord_ptr=cin[i][j]+(cout[i][j]-cin[i][j])/(rout[i][j]-rin[i][j])*
        (r_pb-rin[i][j]);
    }
  
    /*Message("%f %f %f | %f %f %f \n",rin[i][j],r_pb,rout[i][j],
  	       cin[i][j],chord,cout[i][j]);*/
  }

/*---------------------------------------------------------------------------*/
  void my_get_re(int i, cell_t cc, Thread *tc, int count, real chord, 
    real Utotal, real *Re_t_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
  
    if(count == 0) {
      if (dbg == 1)
        Message0("  %d Rotor %d - get_re \n",myid,i+1);
    }
  
  /*C_MU_L(cc,tc) laminar viscosity (in properies panel)
  C_MU_T(cc,tc) turbulent viscosity (computed using turbulence model)
  C_MU_EFF(cc,tc) effective viscosity (sum of laminar and turbulent)
  C_R(cc,tc) density*/
  /*Message("mul %f mut %f mueff %f rho %f \n",
  		C_MU_L(cc,tc), C_MU_T(cc,tc), C_MU_EFF(cc,tc), C_R(cc,tc));*/
  
  
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
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    real kappa;
  
  
  /*Message("R %f cp %f \n",C_RGAS(cc,tc), C_CP(cc,tc));*/
  
    if (rho_const) {
  	  *Ma_t_ptr= 0.2;
      /*Message("rho is const");*/
    }
    else {
      if ((C_CP(cc,tc) - C_RGAS(cc,tc)) == 0.0) {
        Message("Error: kappa for rotor %d is zero.\n",i+1);
        kappa = 0.0;
      } 
      else {
        kappa=C_CP(cc,tc) / (C_CP(cc,tc) - C_RGAS(cc,tc));
      }
      if ((sqrt(kappa*C_RGAS(cc,tc)*C_T(cc,tc))) == 0.0) {
        Message("Error: Mach no for rotor %d is zero.\n",i+1);
        *Ma_t_ptr = 0.0;
      } 
      else {
        *Ma_t_ptr=Utotal/(sqrt(kappa*C_RGAS(cc,tc)*C_T(cc,tc)));
      }
      
      /*Message("rho is ideal gas");*/
    }
  }


/*---------------------------------------------------------------------------*/
  void my_get_aoa(int i, cell_t cc, Thread *tc, int count, real aeps,
	  real theta, real *alpha_t_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
  
    *alpha_t_ptr=aeps+theta;
  
    if(isnan(*alpha_t_ptr)) {
      if(i==1) {
        if(errchk !=1) {
          Message0("Warning, *alpha_t_ptr returns NaN at iter %d \n",count);
          errchk=1;
        }
      }
    }
  }

/*---------------------------------------------------------------------------*/
  void my_get_cl(int i, cell_t cc, Thread *tc, int count, real r_pb,
	  real Ma_t, real Re_t, real alpha_t, real *CL_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    int j,k,it,jt;
    real Re_min, Re_max, Ma_min, Ma_max;
    int ire_min = 0, ire_max = 0, ima_min = 0, ima_max = 0;
    real clre_min = 0, clre_max = 0, clre = 0, clma_min = 0, clma_max = 0, 
      clma = 0, cl = 0, rfac = 0;
    int iflag;
    real cl_m = 0, cl_p = 0;
  
    rfac=0.5;		/*0 --> cl=clma, 1 --> cl=clre*/
  
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
            if (MA[k][it] < Ma_min)	{
  			      Ma_min=MA[k][it];
  			      ima_min=it;
  			    }
            if (MA[k][it] > Ma_max)	{
  			      Ma_max=MA[k][it];
  			      ima_max=it;
  			    }
  		    }
          it+=1;
        }
  
        /*Message("Re_min %f Re_t %f Re_max %f \n",Re_min,Re_t,Re_max);*/
        /*Message("Ma_min %f Ma_t %f Ma_max %f \n",Ma_min,Ma_t,Ma_max);*/
        /*Message("Re_min %f Re_t %f Re_max %f \n",*/
        /*  RE[k][ire_min],Re_t,RE[k][ire_max]);*/
        /*Message("Ma_min %f Ma_t %f Ma_max %f \n",*/
  			/*  MA[k][ima_min],Ma_t,MA[k][ima_max]);*/
  
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
        if (Re_t <= Re_min)	{
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
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
{
  int j,k,it,jt;
  real Re_min, Re_max, Ma_min, Ma_max;
  int ire_min = 0, ire_max = 0, ima_min = 0, ima_max = 0;
  real cdre_min = 0, cdre_max = 0, cdre = 0, cdma_min = 0, cdma_max = 0, 
    cdma = 0, cd = 0, rfac = 0;
  int iflag;
  real cd_m = 0, cd_p = 0;


  rfac=0.5;		/*0 --> cd=cdma, 1 --> cd=cdre*/


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
    real *sum_area_ptr, real sum_area, int count, real CL, real CD, 
    real chord, real aeps, real Utotal, real r_pb,
	  real *Ftc_ptr, real *Fnc_ptr, real *Fsc_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    real FL, FD;
    real aepss, aepsc;
    real tmp;
  
  
    if(count == 0) {
      if (dbg == 1)
        Message0("  %d Rotor %d - get_force \n",myid,i+1);
    }
  
    tmp=chord*0.5*Utotal*Utotal*C_R(cc,tc)*factor_ptr[cc];
  
    if(isnan(tmp)) {
      if(i==1) {
        if(errchk !=1) {
          Message0("Warning, tmp returns NaN at iter %d \n",count);
          errchk=1;
        }
      }
    }
  
    FL=CL*tmp;
    FD=CD*tmp;
  
    /*Tip effect*/
    if (r_pb > teff[i]) {
	    FL=0.0;
      /*Message("Tip effect %f \n",r_pb);*/
    }
  
    aepss=sin(aeps);
    aepsc=cos(aeps);
  
    *Ftc_ptr= FD*aepsc - FL*aepss;
    *Fnc_ptr= FD*aepss + FL*aepsc;
    *Fsc_ptr= 0.0;
  
    if(isnan(*Fnc_ptr)) {
      if(i==1) {
        if(errchk !=1) {
          Message0("Warning, Fnc_ptr returns NaN at iter %d \n",count);
          errchk=1;
        }
      }
    }
  }

/*---------------------------------------------------------------------------*/
  void my_force_blade_to_xyz(int i, cell_t cc, Thread *tc, int count,
	  real Ftc, real Fnc, real Fsc, real beta, real psi_pb,
	  real *Fx_ptr, real *Fy_ptr, real *Fz_ptr, real *Fz_pb_ptr, real *Ft_pb_ptr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    real betas, betac, Fr_pb, Ft_pb, Fz_pb;
    real psi_pbs, psi_pbc, Fx_pb, Fy_pb;
    real dpits, dpitc, dbans, dbanc;
    real tmp1;
  
    if(count == 0) {
      if (dbg == 1) 
        Message0("  %d Rotor %d - force_blade_to_xyz \n",myid,i+1);
    }
  
    betas=sin(beta);
    betac=cos(beta);
    if(rspe[i] > 0.0)
	    tmp1=-1.0;
    if(rspe[i] < 0.0)
	    tmp1=1.0;
  
    Fr_pb= Fsc*betac -          Fnc*betas;
    Ft_pb=             tmp1*Ftc;
    Fz_pb= Fsc*betas +          Fnc*betac;
  
    psi_pbs=sin(psi_pb);
    psi_pbc=cos(psi_pb);
  
    Fx_pb= Fr_pb*psi_pbc - Ft_pb*psi_pbs;
    Fy_pb= Fr_pb*psi_pbs + Ft_pb*psi_pbc;
  
  
    /*Pitch "dpit[i]"*/
    dpits=sin(dpit[i]);
    dpitc=cos(dpit[i]);
    /*Bank "dban[i]"*/
    dbans=sin(dban[i]);
    dbanc=cos(dban[i]);
  
    *Fx_ptr= Fx_pb*dpitc + Fy_pb*dpits*dbans + Fz_pb*dpits*dbanc;
    *Fy_ptr=               Fy_pb*dbanc       - Fz_pb*dbans;
    *Fz_ptr=-Fx_pb*dpits + Fy_pb*dpitc*dbans + Fz_pb*dpitc*dbanc;
  
    *Fz_pb_ptr=Fz_pb;
    *Ft_pb_ptr=Ft_pb;
  }

/*---------------------------------------------------------------------------*/
  void my_get_source_terms(int i, cell_t cc, Thread *tc, int count,
	  real *xsource_ptr, real *ysource_ptr, real *zsource_ptr,
	  real *xsource_old_ptr, real *ysource_old_ptr, real *zsource_old_ptr,
	  real Fx, real Fy, real Fz, real urf_source)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    real vol;
    int i_count;
  
    if(count == 0) {
      if (dbg == 1)
        Message0("  %d Rotor %d - get_source_terms \n",myid,i+1);
    }
  
    vol= C_VOLUME(cc,tc);
  
    if(vol == 0)
      Message0("Warning: Rotor %d volume is zero.\n",i+1);
  
    i_count=N_ITER;
    /*Message("i_count %d istflag %d \n",i_count,istflag);*/
    if(istflag != 1 || i_count == 0) {
      if (vol == 0.0) {
        xsource_ptr[cc] = 0.0;
        ysource_ptr[cc] = 0.0;
        zsource_ptr[cc] = 0.0;
      } 
      else {
        xsource_ptr[cc]=(-Fx/vol)*urf_source + (1.0-urf_source)*xsource_old_ptr[cc];
        ysource_ptr[cc]=(-Fy/vol)*urf_source + (1.0-urf_source)*ysource_old_ptr[cc];
        zsource_ptr[cc]=(-Fz/vol)*urf_source + (1.0-urf_source)*zsource_old_ptr[cc];
      }
  
      /*Message("start or typical");*/
    }
    else {
      if (vol == 0.0) {
        xsource_ptr[cc] = 0.0;
        ysource_ptr[cc] = 0.0;
        zsource_ptr[cc] = 0.0;
      } 
      else {
       xsource_ptr[cc]=(-Fx/vol);
       ysource_ptr[cc]=(-Fy/vol);
       zsource_ptr[cc]=(-Fz/vol);
  	  /*Message("restart");*/
      }
    }
  
    xsource_old_ptr[cc]=xsource_ptr[cc];
    ysource_old_ptr[cc]=ysource_ptr[cc];
    zsource_old_ptr[cc]=zsource_ptr[cc];
  }

/*---------------------------------------------------------------------------*/
  void my_get_cone(int i, real thrust)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/

  {
  /*  real lb, lbx;                                                                                               
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
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    /*Message("rho_ref = %f %f %f %f %f \n",rho_ref,
    coef_fac, rrl[i], rspe[i], M_PI);*/
  
    if (coeff_denom == 0.0) {
      Message("Error: Thrust coeff denominator for rotor %d is zero.\n",i+1);
      CT[i] = 0.0;
    } 
    else {
      CT[i] = thrust/coeff_denom;
    }
  }

/*---------------------------------------------------------------------------*/
  void my_get_cmx(int i, real Mx_pb, real coeff_denom)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    /*Message("rho_ref = %f %f %f %f %f \n",rho_ref,
    coef_fac, rrl[i], rspe[i], M_PI);*/
  
    if (coeff_denom == 0.0) {
  
      Message("Error: X-moment coeff denominator for rotor %d is zero.\n",i+1);
      CMX[i] = 0.0;
    } 
    else {
      CMX[i]= Mx_pb/coeff_denom;
    }  
  }

/*---------------------------------------------------------------------------*/
  void my_get_cmy(int i, real My_pb, real coeff_denom)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
/*---------------------------------------------------------------------------*/
  {
    /*Message("rho_ref = %f %f %f %f %f \n",rho_ref,
    coef_fac, rrl[i], rspe[i], M_PI);*/
  
    if (coeff_denom == 0.0) {
  
      Message("Error: Y-moment coeff denominator for rotor %d is zero.\n",i+1);
      CMY[i] = 0.0;
    } 
    else {
      CMY[i] = My_pb/coeff_denom;
    }  
  }

/*---------------------------------------------------------------------------*/
  void my_force_moment_balancing(int i, real thr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
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
    real *ct_l_ptr, real *ct_s_ptr, real *ct_h_ptr,
    real *bcop_l_ptr, real *bcop_s_ptr, real *bcop_h_ptr,
    real *cmx_lc_ptr, real *cmx_ls_ptr, real *cmx_s_ptr,
    real *cmx_hc_ptr, real *cmx_hs_ptr,
    real *cmy_lc_ptr, real *cmy_ls_ptr, real *cmy_s_ptr,
    real *cmy_hc_ptr, real *cmy_hs_ptr,
    real *bcyc_l_ptr, real *bcyc_s_ptr, real *bcyc_h_ptr,
    real *bcys_l_ptr, real *bcys_s_ptr, real *bcys_h_ptr,
    real *ct_lc_ptr, real *ct_ls_ptr,
    real *ct_hc_ptr, real *ct_hs_ptr,
    real *cmx_l_ptr, real *cmx_h_ptr,
    real *cmy_l_ptr, real *cmy_h_ptr,
    real thr)
/*---------------------------------------------------------------------------*/
/*     								             */
/*     								             */
/*Version	Date	Name			Remarks		             */
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
      if (up_co[i] < trufq_co[i])	/*will be updated at the collective pitch update frequency*/ {
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

          my_force_moment_balancing(i, sum_thrust[i]);

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

{
  Message0("!! Info: Ensure you have pre-allocated at least %d UDM locations in the GUI/TUI !!", NUM_UDM+2);
	Set_User_Memory_Name(0, "Incident Velocity");
	Set_User_Memory_Name(1, "Pitch Angle");
	Set_User_Memory_Name(2, "Angle of Attack");
	Set_User_Memory_Name(3, "Coefficient of Lift");
	Set_User_Memory_Name(4, "Coefficient of Drag");
	Set_User_Memory_Name(5, "Normal Force");
	Set_User_Memory_Name(6, "Tangential Force");
	Set_User_Memory_Name(7, "Time-Averaged X-Force");
	Set_User_Memory_Name(8, "Time-Averaged Z-Force");
  Set_User_Memory_Name(9, "Rotor 1 CT");
}

/*---------------------------------------------------------------------------*/
DEFINE_ON_DEMAND(my_rotor_inputs)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{

  int number_of_sections, ic;
  char exst_type[100][30];
  const char *type_new[10][20];
  int col_pitch, cyc_pitch, fuse_forces, fuse_moments;
  real bcop_udf[10], bcyc_udf[10], bcys_udf[10];
  int trmco_dat[10], trmcy_dat[10];
  register int j,k,m;
  char exfname[100];
  FILE *inp_fp;
  int ch, itrm, inpt_it, pos, ffound;
  real inpt_data1[10][10], inpt_data2[10][10];
  real dskco_0[10], dskco_1[10], dskco_2[10];
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
  
  
      if(fixnum_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),1),15),0), "") == 0) {
        Message("\n Please select a surface for rotor face zone");
        fzon[j] = 0;
      }
      else {
        fzon[j] = fixnum_arg(List_Ref (List_Ref (List_Ref (List_Ref (l,j),1),15),0), "");
      }
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
        itrm=1;
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
          Message("Data file %s for rotor %d found, importing...\n", 
            data_fname[h], j+1);
          ffound = 1;
          fseek(inp_fp, -1, SEEK_END);
          pos = ftell(inp_fp);
          ch = fgetc(inp_fp);
          if (isdigit(ch)){
            while (ch != '\n') {
              
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
	      for(ic=0; ic<ktot; ic++) {
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

    my_get_solidity(sldity);
  
    /*deg_to_rad(dpit, dban, bcop, bcys, bcyc, bflco, bfls, bflc, twst, &dalpha);*/
    /*dalpha=dalpha*M_PI/180.0;*/
  
    my_obtain_cell_id();
  
    if (dbg == 1) {
      for (int i = 0; i < nrtz; i++) {
        Message0("Rotor inputs chkpt1, czon[%d] = %d \n", i, czon[i]);
      }
    }

    if (istflag == 1) 
      my_allocate_memory(factor,xsource,ysource,zsource,
        xsource_old,ysource_old,zsource_old);
  
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
    if(dbg == 1) {
      for (int i = 0; i < nrtz; i++) {
        Message0("Rotor inputs chkpt4, czon[%d] = %d \n", i, czon[i]);
        Message0("Rotor inputs chkpt4, rspe[%d] = %f \n", i, rspe[i]);
      }
    }

 
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
    
    if (dbg == 1){
      for (int i = 0; i < 10; i++) {
        Message0("Rotor inputs chkpt5, czon[%d] = %d \n", i, czon[i]);
        Message0("Rotor inputs chkpt5, rspe[%d] = %f \n", i, rspe[i]);
      }
    }
    

    for(int k = 0; k < ktot; k++){
      strncpy(check_name[k], &check_name_flatten[k*30], 30);

      for(int i = 0; i < itot[k]; i++){
        strncpy(clorcd[k][i], &clorcd_flatten[k*500 + i*10], 10);
      }
    }

    if (dbg == 1){
      for (int i = 0; i < 10; i++) {
        Message0("Rotor inputs chkpt6, czon[%d] = %d \n", i, czon[i]);
        Message0("Rotor inputs chkpt6, rspe[%d] = %f \n", i, rspe[i]);
      }
    } 
  #endif

}

/*---------------------------------------------------------------------------*/
DEFINE_ADJUST(my_SrcComp, domain)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    int i, count;
    Thread *tc;
    cell_t cc;
    real *factor_ptr, sum_area;
    real r_pb,psi_pb,r_pb2;
    real Usc, Utc, Unc, aeps, Utotal, beta;
    real theta, chord;
    real CL, CD;
    real Ftc, Fnc, Fsc;
    real Fx, Fy, Fz;
    real *xsource_ptr, *ysource_ptr, *zsource_ptr;
    real *xsource_old_ptr, *ysource_old_ptr, *zsource_old_ptr;
    real Re_t, Ma_t, alpha_t;
    real thrust, Fz_pb;
    real x_pb, y_pb, Mx_pb, My_pb;
    int flag_jac;
    real ct_s, ct_h, ct_l, bcop_l, bcop_s, bcop_h;
    real cmx_lc, cmx_ls, cmx_s, cmx_hc, cmx_hs;
    real cmy_lc, cmy_ls, cmy_s, cmy_hc, cmy_hs;
    real bcyc_l, bcyc_s, bcyc_h, bcys_l, bcys_s, bcys_h;
    real ct_lc, ct_ls, ct_hc, ct_hs, cmx_l, cmx_h, cmy_l, cmy_h;
    real alpha_tmax, alpha_tmin;
    real torque, Ft_pb;
    real sum_Mx_pb, sum_My_pb;
    real rho_ref;
    real coeff_denom;
    FILE *fp;
    int istr_target = 0;
    /*real x_sum, y_sum, z_sum;  Sums for coordinates 
    real x_avg, y_avg, z_avg;       Averages for coordinates 
    real x[ND_ND];                  Array to store cell centroid */
  
    if (dbg == 1){
      for (i = 0; i < nrtz; i++) {
        Message("In DEFINE_ADJUST, czon[%d] = %d \n", i, czon[i]);
      }
    }


    if ((hdr == 1 || N_ITER == 0) && myid == 0) {
      for (i = 0; i < nrtz; i++) {
        for (int h = 0; h < (sizeof(data_fname) / sizeof(data_fname[0])); h++){

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

    if (dbg == 1)
      Message0("%d In DEFINE_ADJUST - SrcComp \n", myid);

    /*Loop through all zones*/
    i=0;
    while (i < nrtz) {
      if (dbg == 1){
	      Message0("%d Rotor %d of %d \n",myid, i+1, nrtz);
        Message0("Rotor %d czon = %d \n",i+1, czon[i]);
      }
      
      tc = Lookup_Thread(domain, czon[i]);
      if (tc == NULL) {
        Message0("Error: Cell thread is NULL in define_adjust.\n");
        return; 
      }
      factor_ptr = Get_Thread_Memory(tc,factor[i]);
      if (factor_ptr == NULL) {
        Message0("Error: factor_ptr is NULL in define_adjust.\n");
        return; 
      }

      xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
      if (xsource_ptr == NULL) {
        Message0("Error: xsource_ptr %d is NULL in define_adjust.\n", i);
      }
      ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
      if (ysource_ptr == NULL) {
        Message0("Error: ysource_ptr %d is NULL in define_adjust.\n", i);
      }
      zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
      if (zsource_ptr == NULL) {
        Message0("Error: zsource_ptr %d is NULL in define_adjust.\n", i);
      }

      zsource_old_ptr= Get_Thread_Memory(tc,zsource_old[i]);
      if (zsource_old_ptr == NULL) {
        Message0("Error: zsource_old_ptr %d is NULL in define_adjust.\n", i);
      }
      ysource_old_ptr= Get_Thread_Memory(tc,ysource_old[i]);
      if (ysource_old_ptr == NULL) {
        Message0("Error: ysource_old_ptr %d is NULL in define_adjust.\n", i);
      }
      xsource_old_ptr= Get_Thread_Memory(tc,xsource_old[i]);
      if (xsource_old_ptr == NULL) {
        Message0("Error: xsource_old_ptr %d is NULL in define_adjust.\n", i);
      }

	    /*count=0;
      begin_c_loop_int(cc,tc)
        count++;
      end_c_loop_int(cc,tc)*/

      rho_ref = RP_Get_Real("reference-density");

      coeff_denom = coef_fac * rho_ref * rrl[i] * rrl[i] * M_PI * rspe[i] * 
        rspe[i] * rrl[i] * rrl[i];

	    flag_jac=0;
      LABEL1:
	    thrust=0.0;
	    torque=0.0;
	    Mx_pb=0.0;
	    My_pb=0.0;
      sum_area=0.0;
	    alpha_tmax=-1000.0;
	    alpha_tmin=1000.0;
      count=0;

      begin_c_loop_int(cc,tc)

        my_get_radial_position(i, cc, tc, count, 
          &r_pb, &psi_pb, &r_pb2,
	  	    &x_pb,&y_pb);
	      /*Message0("  r_pb %f     psi_pb %f  x_pb %f   y_pb %f \n",
	  	    r_pb, psi_pb*180/M_PI, x_pb, y_pb); */

	      my_vel_xyz_to_blade(i, cc, tc, count, 
          psi_pb, r_pb2, 
          &Usc, &Utc, &Unc,
          &aeps, &Utotal, &beta);
        /*Message0("  Usc %f Utc %f Unc %f \n",Usc, Utc, Unc);*/
	  	  /*Message0("aeps %f \n", aeps*180.0/M_PI);*/

	  	  my_get_pitch(i, cc, tc, count, r_pb, 
          psi_pb, &theta);
	  	  /*Message0("  psi_pb %f theta %f \n",psi_pb*180/M_PI, theta*180/M_PI);*/

	  	  my_get_chord(i, cc, tc, count, r_pb, 
          &chord);
	  	  /*Message0("%f \n",chord);*/

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
          &sum_area, sum_area, count, CL, CD, 
          chord, aeps, Utotal, r_pb, 
          &Ftc, &Fnc, &Fsc);
	  	  /*Message0("Ftc=%f Fnc=%f Fsc=%f  \n",Ftc, Fnc, Fsc);*/

        my_force_blade_to_xyz(i, cc, tc, count, 
          Ftc, Fnc, Fsc, beta, psi_pb,
	  	    &Fx, &Fy, &Fz, &Fz_pb, &Ft_pb);

        if (flag_jac == 0) 
          my_get_source_terms(i, cc, tc, count, 
            xsource_ptr,ysource_ptr,zsource_ptr,
	  	      xsource_old_ptr,ysource_old_ptr,zsource_old_ptr,
	  	      Fx, Fy, Fz, urf_source);

	  	  thrust=thrust+Fz_pb;

	  	  Mx_pb=Mx_pb+Fz_pb*y_pb;
	  	  My_pb=My_pb+Fz_pb*x_pb;
        torque=torque+Ft_pb*sqrt(x_pb*x_pb + y_pb*y_pb);

        C_UDMI(cc, tc, 0) = Utotal;
        C_UDMI(cc, tc, 1) = theta * (180.0 / M_PI);
        C_UDMI(cc, tc, 2) = (aeps + theta) * (180.0 / M_PI);
        C_UDMI(cc, tc, 3) = CL;
        C_UDMI(cc, tc, 4) = CD;
	  	  C_UDMI(cc, tc, 5) = Fz_pb;
	  	  C_UDMI(cc, tc, 6) = Ft_pb;
        C_UDMI(cc, tc, 7) = Fx;
	  	  C_UDMI(cc, tc, 8) = Fz;
        C_UDMI(cc, tc, 9) = Fz_pb/coeff_denom;

	  	  alpha_tmax=fmax(alpha_tmax, alpha_t);
	  	  alpha_tmin=fmin(alpha_tmin, alpha_t);

        /*C_CENTROID(x, cc, tc); Get the centroid of the cell*/
        /*x_sum += x[0];*/
        /*y_sum += x[1];*/
        /*z_sum += x[2];*/
      
	      count++;
      end_c_loop_int(cc,tc)

      /*x_avg = x_sum / count;*/
      /*y_avg = y_sum / count;*/
      /*z_avg = z_sum / count;*/
      /*Message0("Average coordinates for Rotor %d: X = %f, Y = %f, Z = %f\n", i+1, x_avg, y_avg, z_avg);*/

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
	  	    /*Message0("    %d Rotor %d Mx_pb %f cmx= %f \n",myid,i+1,sum_Mx_pb,CMX); */
          /*Message("    %d Rotor %d Mx_pb on node %f \n",myid,i+1,Mx_pb);*/
	  	    /*Message0("    %d Rotor %d My_pb %f cmy= %f \n",myid,i+1,sum_My_pb,CMY); */
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
DEFINE_SOURCE(my_xmom_src_1,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=0;				/*Rotor 1*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
  
    if (xsource_ptr == NULL){
      Message0("Error: xsource_ptr is NULL in DEFINE_SOURCE(my_xmom_src_1)\n");
      source = 0.0;
    }
    else{
      source=xsource_ptr[cc];
    }
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_1,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=0;				/*Rotor 1*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
  
    if (ysource_ptr == NULL) {
      Message0("Error: ysource_ptr is NULL in DEFINE_SOURCE(my_ymom_src_1)\n");
      source = 0.0;
    }
    else{
      source=ysource_ptr[cc];
    }
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_1,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=0;				/*Rotor 1*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
  
    if (zsource_ptr == NULL) {
      Message0("Error: zsource_ptr is NULL in DEFINE_SOURCE(my_zmom_src_1)\n");
      source = 0.0;
    }
    else{
      source=zsource_ptr[cc];
    }
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_2,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=1;				/*Rotor 2*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
  
    if (xsource_ptr == NULL) {
      Message0("Error: xsource_ptr is NULL in DEFINE_SOURCE(my_xmom_src_1)\n");
      source = 0.0;
    }
    else{
      source=xsource_ptr[cc];
    }
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_2,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=1;				/*Rotor 2*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
  
    if (ysource_ptr == NULL) {
      Message0("Error: ysource_ptr is NULL in DEFINE_SOURCE(my_ymom_src_2)\n");
      source = 0.0;
    }
    else{
      source=ysource_ptr[cc];
    }
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_2,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=1;				/*Rotor 2*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
  
    if (zsource_ptr == NULL) {
      Message0("Error: zsource_ptr is NULL in DEFINE_SOURCE(my_zmom_src_2)\n");
      source = 0.0;
    }
    else{
      source=zsource_ptr[cc];
    }
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_3,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=2;				/*Rotor 3*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
  
    if (xsource_ptr == NULL) {
      Message("Error: xsource_ptr is NULL in DEFINE_SOURCE(my_xmom_src_3)\n");
      source = 0.0;
    }
    else{
      source=xsource_ptr[cc];
    }
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_3,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=2;				/*Rotor 3*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_3,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=2;				/*Rotor 3*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_4,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=3;				/*Rotor 4*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
    source=xsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_4,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=3;				/*Rotor 4*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_4,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=3;				/*Rotor 4*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_5,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=4;				/*Rotor 5*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
    source=xsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_5,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=4;				/*Rotor 5*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_5,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=4;				/*Rotor 5*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_6,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=5;				/*Rotor 6*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
    source=xsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_6,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=5;				/*Rotor 6*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_6,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=5;				/*Rotor 6*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_7,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=6;				/*Rotor 7*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
    source=xsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_7,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=6;				/*Rotor 7*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_7,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=6;				/*Rotor 7*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_8,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=7;				/*Rotor 8*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
    source=xsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_8,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=7;				/*Rotor 8*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_8,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=7;				/*Rotor 8*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_9,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=8;				/*Rotor 9*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
    source=xsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_9,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=8;				/*Rotor 9*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_9,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=8;				/*Rotor 9*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_xmom_src_10,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *xsource_ptr;
    real source = 0.0;
    int i=9;				/*Rotor 10*/
  
    xsource_ptr= Get_Thread_Memory(tc,xsource[i]);
    source=xsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_ymom_src_10,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *ysource_ptr;
    real source = 0.0;
    int i=9;				/*Rotor 10*/
  
    ysource_ptr= Get_Thread_Memory(tc,ysource[i]);
    source=ysource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_SOURCE(my_zmom_src_10,cc,tc,dS,eqn)
/*---------------------------------------------------------------------------*/
/*				       				       */
/*				       				       */
/*Version	Date	Name			Remarks	       	       */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    
    real *zsource_ptr;
    real source = 0.0;
    int i=9;				/*Rotor 10*/
  
    zsource_ptr= Get_Thread_Memory(tc,zsource[i]);
    source=zsource_ptr[cc];
  
    return source;
  #endif
}

/*---------------------------------------------------------------------------*/
DEFINE_EXECUTE_AT_END(my_write_rotor_data_to_file)
/*---------------------------------------------------------------------------*/
/* When calling the file saving function, inputs are in the following order: */
/* char File_name, char Data_name, real data_to_write, char time_str,        */ 
/* int is_angle? The purpose of the final input, is_angle, is to convert */
/* the data to degrees.                                                      */
/*---------------------------------------------------------------------------*/
{
  #if !RP_HOST
    if(trnw == 1 && N_ITER > 1 && myid == 0) {
      FILE *fp = NULL;
      int i = 0, h = 0;
      int isangle = 0, istr_target = 0;
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

        for (h = 0; h < (sizeof(data_fname) / sizeof(data_fname[0])); h++) {
          isangle = 0;
          istr_target = 0;

          if (h < 3) {
            isangle = 1;
          }
          if (h > 5) {
            istr_target = 1;
          }

          sprintf(newfilename, "rotor_%d_%s", i + 1, data_fname[h]);

          fp = fopen(newfilename, "a");
          if (fp == NULL) {
            Message0("Error: Unable to open file for writing data.\n");
          }
          else {
            if (isangle == 1){
              if (dbg == 1)
                Message0("Writing Rotor %d %s = %f to %s\n", i+1, 
                  data_vname[h], data_save[h]*180/M_PI, newfilename);

              fprintf(fp, "\n%d\t%f", N_ITER, data_save[h]*180/M_PI);
            }
            if (isangle == 0){
              if (istr_target == 1){
                if (dbg == 1)
                  Message0("Writing Rotor %d %s = %f to %s\n", i+1, 
                    data_vname[h], data_save[h], newfilename);

                fprintf(fp, "\n%d\t%f\t%f", N_ITER, data_save[h], 
                  tr_target[h-6]);
              }
              if (istr_target == 0){
                if (dbg == 1)
                  Message0("Writing Rotor %d %s = %f to %s\n", i+1, 
                    data_vname[h], data_save[h], newfilename);

                fprintf(fp, "\n%d\t%f", N_ITER, data_save[h]);
              }
            }
            fclose(fp);
            if (dbg == 1)
              Message0("%s values written to %s\n", data_vname[h], 
                newfilename);
          }
        }   
      }
    }
    trnw = 0;
  #endif
}
