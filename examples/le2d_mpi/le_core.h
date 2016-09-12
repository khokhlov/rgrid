#ifndef LE_CORE_H
#define LE_CORE_H

#include <sys/time.h>

#include "rgrid/darrayscatter.h"

static inline double timer()
{
	struct timeval theStartTime;
	gettimeofday(&theStartTime, NULL);
	return theStartTime.tv_sec + 1e-6 * theStartTime.tv_usec;
}

/* Real type. */
#ifdef USE_DOUBLE
typedef double real;
#else
typedef float real;
#endif

/* Integer type. */
typedef int int_t;

typedef rgrid::DArrayScatter<real, int_t> DAS;
typedef rgrid::DArrayContainer<real, int_t> DAC;
typedef rgrid::DArray<real, int_t> DA;

typedef rgrid::DArrayScatter<float, int_t> DAS_F;
typedef rgrid::DArrayContainer<float, int_t> DAC_F;
typedef rgrid::DArray<float, int_t> DA_F;

/* 2d vector struct. */
typedef struct {
	real x, y;
} le_vec2;

/*
 * Elastic node.
 * Store velocity and stress tensor (http://en.wikipedia.org/wiki/Stress_%28mechanics%29).
 */
enum le_node {
	/* 2d vector. */
	VX,
	VY,
	/* 
	 * Symmetrical tensor 2d.
	 * sxx sxy
	 * sxy sxx
	 */
	SXX,
	SXY,
	SYY,
	NODE_SIZE
};

/* Integer 2d point. */
typedef struct {
	int_t x, y;
} le_point2;

/*
 * Riemann invariant structure (http://en.wikipedia.org/wiki/Riemann_invariant).
 */
typedef struct {
	real w1, w2, w3, w4/*, w5*/;
} le_w;

/*
 * Elastic material structure.
 * We use simple approach when all parameters are constant on whole region.
 * c1 - speed of P-wave http://en.wikipedia.org/wiki/P-wave.
 * c2 - speed of S-wave http://en.wikipedia.org/wiki/S-wave.
 * rho - density.
 */
typedef struct {
	real c1, c2, rho;
	
	/*
	 * Some cached values to speedup calculations.
	 */
	real irhoc1; // 1.0 / (c1 * rho)
	real irhoc2; // 1.0 / (c2 * rho)
	real rhoc1; // c1 * rho
	real rhoc2; // c2 * rho
	real rhoc3; // c3 * rho
} le_material;

/*
 * Structure for storing all parameters of task.
 */
typedef struct {
	/* Time step. */
	real dt;
	
	/* Grid spacing. */
	le_vec2 h;
	
	/* Number of nodes in grid on each axis. */
	le_point2 n;
	
	/* Material. */
	le_material mat;
	
	/* Grid data (nodes). */
	DA da; // initialization
	DAS das;
	DAS_F dasSave;
} le_task;

/*
 * Create material and init all fields of structure.
 */
void le_init_material(const real c1, const real c2, const real rho, le_material &m);

/* Create task with given parameters. Allocate memory for nodes. */
void le_init_task(le_task &task, const real dt, const le_vec2 h, const le_material mat, const le_point2 n);

/* Free memory. */
void le_free_task(le_task &task);

/*
 * Set initial disturbance on the grid.
 */
void le_set_ball(le_task &t, const le_vec2 c, const real r, const real s);

/*
 * Save grid to legasy VTK format (http://www.vtk.org/VTK/img/file-formats.pdf).
 * You can use ParaView (http://www.paraview.org/),
 * MayaVi (http://mayavi.sourceforge.net/) or VisIt (https://wci.llnl.gov/codes/visit/)
 * to visualize results.
 * Return: 0 - all ok, 1 - error.
 */
int le_save_task(le_task &task, const char *file);

/*
 * One time step of difference scheme.
 */
void le_step(le_task &task);

#endif // LE_CORE_H
