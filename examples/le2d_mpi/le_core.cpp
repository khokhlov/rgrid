#include <cstdlib>
#include <cstdio>
#include <cmath>
#include <cassert>
#include <cstring>

#include "le_core.h"

#include "rgrid/darrayscatter.h"

inline real le_min(real a, real b) { return a > b ? b : a; }
inline real le_max(real a, real b) { return a > b ? a : b; }
inline real le_max3(real a, real b, real c) { return le_max(a, le_max(b, c)); }

#define limiter_superbee(r) (le_max3(0.0, le_min(1.0, 2.0 * r), le_min(2.0, r)))

/*
 * Set TVD limiter (http://en.wikipedia.org/wiki/Flux_limiters).
 */
#define limiter limiter_superbee

#define TVD2_EPS 1e-6

inline real& gind(real* grid, int_t x, int_t y, int_t c, int_t nx, int_t ny, int_t cn) {
#ifdef USE_AOS
	return grid[c + x * cn + y * cn * nx];
#else // USE_SOA
	return grid[x + y * nx + c * nx * ny];
#endif
}
#define ind(x,y,nx) ((x) + (y) * (nx))

/*
 * Second order TVD scheme (http://en.wikipedia.org/wiki/Total_variation_diminishing).
 */
inline real tvd2(const real c, const real u_2, const real u_1, const real u, const real u1)
{
	real r1 = (u  - u_1);
	real r2 = (u1 - u);
	//if (r2 == 0.0) {
		r1 += TVD2_EPS;
		r2 += TVD2_EPS;
	//}
	const real r = r1 / r2;
	r1 = (u_1 - u_2);
	r2 = (u   - u_1);
	//if (r2 == 0.0) {
		r1 += TVD2_EPS;
		r2 += TVD2_EPS;
	//}
	const real r_1 = r1 / r2;

	const real f12  = u   + limiter(r)   / 2.0 * (1.0 - c) * (u1 - u);
	const real f_12 = u_1 + limiter(r_1) / 2.0 * (1.0 - c) * (u  - u_1);

	return u - c * (f12 - f_12);
}

void le_set_ball(le_task &t, const le_vec2 c, const real r, const real s)
{
	int i, j;
	for (i = 0; i < t.n.x; i++) {
		for (j = 0; j < t.n.y; j++) {
			le_vec2 x = {t.h.x * i, t.h.y * j};
			le_vec2 d = {x.x - c.x, x.y - c.y};
			if (sqrt(d.x*d.x + d.y*d.y) < r) {
				/*
				 * Set pressure disturbance
				 */
				t.da.val(i, j, 0, SXX) = s;
				t.da.val(i, j, 0, SYY) = s;
			}
		}
	}
}

void le_init_task(le_task &task, const real dt, const le_vec2 h, const le_material mat, const le_point2 n)
{
	task.dt  = dt;
	task.h   = h;
	task.mat = mat;
	task.n   = n;
}

void le_free_task(le_task& task)
{
}

int le_save_task(le_task &t, const char *file)
{	
	std::stringstream ss;

	ss.precision(6); // diff with C version of this code
	ss << std::fixed;

	ss << "# vtk DataFile Version 3.0" << std::endl;
	ss << "Created by le_save_task" << std::endl;
	ss << "BINARY" << std::endl;
	ss << "DATASET STRUCTURED_POINTS" << std::endl;
	ss << "DIMENSIONS " << t.n.x << " " << t.n.y << " 1" << std::endl;
	ss << "SPACING " << t.h.x << " " << t.h.y << " 0.0" << std::endl;
	ss << "ORIGIN 0.0 0.0 0.0" << std::endl;
	ss << "POINT_DATA " << t.n.x * t.n.y << std::endl;
	
	/* velocity */
	ss << "SCALARS v float 1" << std::endl;
	ss << "LOOKUP_TABLE v_table" << std::endl;

	DA_F& ds = t.dasSave.getDArrayPart(0, 0, 0);
	DAC& dac = t.das.getLocalContainer();
	for (int_t k = 0; k != dac.numParts(rgrid::Z); ++k)
	for (int_t j = 0; j != dac.numParts(rgrid::Y); ++j)
	for (int_t i = 0; i != dac.numParts(rgrid::X); ++i) {
		DA& d = dac.getDArrayPart(i, j, k);
		for (int_t k2 = 0; k2 != d.localSize(rgrid::Z); ++k2)
		for (int_t j2 = 0; j2 != d.localSize(rgrid::Y); ++j2)
		for (int_t i2 = 0; i2 != d.localSize(rgrid::X); ++i2) {
			ds(i2 + dac.partOrigin(rgrid::X, i), j2 + dac.partOrigin(rgrid::Y, j), k2 + dac.partOrigin(rgrid::Z, k), 0) =
				sqrt(d(i2, j2, k2, VX)*d(i2, j2, k2, VX)+d(i2, j2, k2, VY)*d(i2, j2, k2, VY));
		}
	}
	ds.inverseBytes();
	t.dasSave.saveDataBegin(std::string(file), rgrid::rgio::CUSTOM_HEADER, &ss);
	t.dasSave.saveDataEnd();

	return 0;
}

void le_init_material(const real c1, const real c2, const real rho, le_material& m)
{
	m.c1 = c1;
	m.c2 = c2;
	m.rho = rho;

	/*
	 * Cached values.
	 */
	m.irhoc1 = 1.0 / (c1 * rho);
	m.irhoc2 = 1.0 / (c2 * rho);
	m.rhoc1 = c1 * rho;
	m.rhoc2 = c2 * rho;
	real mu = rho * c2 * c2;
	real la = rho * c1 * c1 - 2.0 * mu;
	m.rhoc3 = rho * c1 * la / (la + 2.0 * mu);
}

/*
 * w = OmegaR * u.
 */
inline void omega_x(const le_material &m, 
                    le_w& w,
                    const real vx,
                    const real vy,
                    const real sxx,
                    const real sxy,
                    const real syy
)
{
	const real nv = vx;
	const real N00T = sxx * m.irhoc1;

	const real n1v = vy;
	const real N01T = sxy * m.irhoc2;

	w.w1 = nv  - N00T;
	w.w2 = nv  + N00T;
	w.w3 = n1v - N01T;
	w.w4 = n1v + N01T;
}

inline void omega_y(const le_material &m,
                    le_w& w,
                    const real vx,
                    const real vy,
                    const real sxx,
                    const real sxy,
                    const real syy
                   )
{
	const real nv = vy;
	const real N00T = syy * m.irhoc1;

	const real n1v = vx;
	const real N01T = sxy * m.irhoc2;

	w.w1 = nv  - N00T;
	w.w2 = nv  + N00T;
	w.w3 = n1v - N01T;
	w.w4 = n1v + N01T;
}

inline void inc_x(const le_material &m,
                  const le_w &d,
                  real& vx,
                  real& vy,
                  real& sxx,
                  real& sxy,
                  real& syy
                 )
{
	const real d1 = 0.5 * d.w1;
	const real d2 = 0.5 * d.w2;
	const real d3 = 0.5 * d.w3;
	const real d4 = 0.5 * d.w4;

	vx += d1 + d2;
	vy += d3 + d4;

	sxx += (d2 - d1) * m.rhoc1;
	syy += (d2 - d1) * m.rhoc3;
	sxy += m.rhoc2 * (d4 - d3);
}

inline void inc_y(const le_material &m,
                  const le_w &d,
                  real& vx,
                  real& vy,
                  real& sxx,
                  real& sxy,
                  real& syy
                 )
{
	const real d1 = 0.5 * d.w1;
	const real d2 = 0.5 * d.w2;
	const real d3 = 0.5 * d.w3;
	const real d4 = 0.5 * d.w4;

	vy += d1 + d2;
	vx += d3 + d4;

	syy += (d2 - d1) * m.rhoc1;
	sxx += (d2 - d1) * m.rhoc3;
	sxy += m.rhoc2 * (d4 - d3);
}

inline void reconstruct(const le_w& ppu, const le_w& pu, const le_w& u, const le_w& nu, const le_w& nnu, const real k1, const real k2, le_w& d)
{
	d.w1 = tvd2(k1, ppu.w1, pu.w1, u.w1, nu.w1) - u.w1; // c1
	d.w2 = tvd2(k1, nnu.w2, nu.w2, u.w2, pu.w2) - u.w2; // -c1
	d.w3 = tvd2(k2, ppu.w3, pu.w3, u.w3, nu.w3) - u.w3; // c2
	d.w4 = tvd2(k2, nnu.w4, nu.w4, u.w4, pu.w4) - u.w4; // -c2
}

void le_step_x(le_task& t)
{
	/*
	 * Due to our system of pde is linear, we can use some simple way to solve it.
	 * du/dt + A * du/dx = 0.
	 * Vector u = {vx, vy, sxx, sxy, syy}.
	 * Matrix A could be represent in form OmegaL * Lambda * OmegaR,
	 * where Lambda - diagonal matrix of eigen values of matrix A.
	 * In our case Lambda = diag{c1, -c1, c2, -c2, 0}.
	 * OmegaR and OmegaL - metrices from eigen vectors of matrix A,
	 * OmegaR * OmegaL = E, where E = diag{1, 1, 1, 1, 1}.
	 * 
	 * We can rewrite out system in form:
	 * du/dt + OmegaL * Lambda * OmegaR du/dx = 0, multiply on matrix OmegaR:
	 * 
	 * OmegaR * du/dt + OmegaR * OmegaL * Lambda * OmegaR du/dx = 0.
	 * 
	 * Introduce new variables (http://en.wikipedia.org/wiki/Riemann_invariant):
	 * w = {w1, w2, w3, w4, w5},
	 * w = OmegaR * u, then we got:
	 * 
	 * dw/dt + Lambda * dw/dx = 0.
	 * And we get system of independent advection equations, that we can solve separatly.
	 * 
	 * So we get next algorithm:
	 * 1. Introduce new variables w = OmegaR * u;
	 * 2. Solve 5 equations of linear advection (in real we solve only 4, because of in fifth equation speed is 0);
	 * 3. Make inverse transformation u = OmegaL * w. 
	 */
	int i, j;
	
	// Courant number (http://en.wikipedia.org/wiki/Courant%E2%80%93Friedrichs%E2%80%93Lewy_condition).
	const real k1 = t.dt * t.mat.c1 / t.h.x;
	const real k2 = t.dt * t.mat.c2 / t.h.x;

	for (int_t partZ = 0; partZ != t.das.getLocalContainer().numParts(rgrid::Z); ++partZ)
	for (int_t partY = 0; partY != t.das.getLocalContainer().numParts(rgrid::Y); ++partY)
	for (int_t partX = 0; partX != t.das.getLocalContainer().numParts(rgrid::X); ++partX) {
		DA& lda = t.das.getDArrayPart(partX, partY, partZ);
		for (j = 0; j < lda.localSize(rgrid::Y); j++) {
			/*
			 * Riemann invariants for 5-point sctencil difference scheme.
			 */
			le_w w_2, w_1, w, w1, w2;

			omega_x(t.mat, w_2,
			        lda.val(-2, j, VX),
			        lda.val(-2, j, VY),
			        lda.val(-2, j, SXX),
			        lda.val(-2, j, SXY),
			        lda.val(-2, j, SYY)
			       );
			omega_x(t.mat, w_1,
			        lda.val(-1, j, VX),
			        lda.val(-1, j, VY),
			        lda.val(-1, j, SXX),
			        lda.val(-1, j, SXY),
			        lda.val(-1, j, SYY)
			       );
			omega_x(t.mat, w,
			        lda.val(0, j, VX),
			        lda.val(0, j, VY),
			        lda.val(0, j, SXX),
			        lda.val(0, j, SXY),
			        lda.val(0, j, SYY)
			       );
			omega_x(t.mat, w1,
			        lda.val(1, j, VX),
			        lda.val(1, j, VY),
			        lda.val(1, j, SXX),
			        lda.val(1, j, SXY),
			        lda.val(1, j, SYY)
			       );
			
			for (i = 0; i < lda.localSize(rgrid::X); i++) {
				omega_x(t.mat, w2,
				        lda.val(i + 2, j, VX),
				        lda.val(i + 2, j, VY),
				        lda.val(i + 2, j, SXX),
				        lda.val(i + 2, j, SXY),
				        lda.val(i + 2, j, SYY)
				       );
				le_w d;
				reconstruct(w_2, w_1, w, w1, w2, k1, k2, d);
				inc_x(t.mat, d,
				      lda.val(i, j, VX),
				      lda.val(i, j, VY),
				      lda.val(i, j, SXX),
				      lda.val(i, j, SXY),
				      lda.val(i, j, SYY)
				     );
				w_2 = w_1;
				w_1 = w;
				w   = w1;
				w1  = w2;
			}
		}
	}
}

void le_step_y(le_task &t)
{
	int i, j;
	const real k1 = t.dt * t.mat.c1 / t.h.y;
	const real k2 = t.dt * t.mat.c2 / t.h.y;
	
	for (int_t partZ = 0; partZ != t.das.getLocalContainer().numParts(rgrid::Z); ++partZ)
	for (int_t partY = 0; partY != t.das.getLocalContainer().numParts(rgrid::Y); ++partY)
	for (int_t partX = 0; partX != t.das.getLocalContainer().numParts(rgrid::X); ++partX) {
		DA& lda = t.das.getDArrayPart(partX, partY, partZ);
		for (i = 0; i < lda.localSize(rgrid::X); i++) {
			le_w w_2, w_1, w, w1, w2;
			
			omega_y(t.mat, w_2,
			        lda.val(i, -2, VX),
			        lda.val(i, -2, VY),
			        lda.val(i, -2, SXX),
			        lda.val(i, -2, SXY),
			        lda.val(i, -2, SYY)
			       );
			omega_y(t.mat, w_1,
			        lda.val(i, -1, VX),
			        lda.val(i, -1, VY),
			        lda.val(i, -1, SXX),
			        lda.val(i, -1, SXY),
			        lda.val(i, -1, SYY)
			       );
			omega_y(t.mat, w,
			        lda.val(i, 0, VX),
			        lda.val(i, 0, VY),
			        lda.val(i, 0, SXX),
			        lda.val(i, 0, SXY),
			        lda.val(i, 0, SYY)
			       );
			omega_y(t.mat, w1,
			        lda.val(i, 1, VX),
			        lda.val(i, 1, VY),
			        lda.val(i, 1, SXX),
			        lda.val(i, 1, SXY),
			        lda.val(i, 1, SYY)
			       );
			
			for (j = 0; j < lda.localSize(rgrid::Y); j++) {
				omega_y(t.mat, w2,
				        lda.val(i, j + 2, VX),
				        lda.val(i, j + 2, VY),
				        lda.val(i, j + 2, SXX),
				        lda.val(i, j + 2, SXY),
				        lda.val(i, j + 2, SYY)
				       );
				le_w d;
				reconstruct(w_2, w_1, w, w1, w2, k1, k2, d);
				inc_y(t.mat, d,
				      lda.val(i, j, VX),
				      lda.val(i, j, VY),
				      lda.val(i, j, SXX),
				      lda.val(i, j, SXY),
				      lda.val(i, j, SYY)
				     );
				w_2 = w_1;
				w_1 = w;
				w   = w1;
				w1  = w2;
			}
		}
	}
}


void le_step(le_task &task)
{
	/*
	 * We solve regular hyperbolic system of PDE (http://en.wikipedia.org/wiki/Hyperbolic_partial_differential_equation) in form:
	 * du/dt + Ax * du/dx + Ay * du/dy = 0.
	 * 
	 * During time integration we use dimension split method:
	 * 1. Step:
	 * Integrate system dv/dt + Ax * dv/dx = 0, get u = v^(n + 1).
	 * 2. Step:
	 * Integrate system du/dt + Ay * du/dy = 0, get on next time step u^(n + 1).
	 */
	
	
	task.das.externalSyncStart();
	task.das.internalSync();
	task.das.fillGhost();
	task.das.externalSyncEnd();
	le_step_x(task);

	task.das.externalSyncStart();
	task.das.internalSync();
	task.das.fillGhost();
	task.das.externalSyncEnd();
	le_step_y(task);
}
