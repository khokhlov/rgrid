#include <cstdio>
#include <cstdlib>

#include "le_core.h"

#include "rgrid/darrayscatter.h"

int main(int argc, char **argv) {
	rgmpi::init(&argc, &argv);
	{ // this brace because we want all destructors to be executed before MPI_Finalize
		if (argc != 8) {
			if (rgmpi::worldRank() == 0)
				printf("usage: %s nx ny steps gx gy lx ly\n", argv[0]);
			rgmpi::forceFinalize();
			return 0;
		}

		int i;
		le_point2 n = {atoi(argv[1]), atoi(argv[2])};
		int steps = atoi(argv[3]);
		le_task task;
		le_material mat;
		le_vec2 h = {1.0, 1.0};
		real dt = 0.3;
		le_vec2 center;
		center.x = n.x / 2;
		center.y = n.y / 2;
		double t;
		
		int_t globalPt[3] = {atoi(argv[4]), atoi(argv[5]), 1};
		int_t localPt[3] = {atoi(argv[6]), atoi(argv[7]), 1};
		if (task.das.getInternalRank() == 1) {
			localPt[rgrid::X] = 3;
		}
		if (task.das.getInternalRank() == 4) {
			localPt[rgrid::Y] = 3;
		}
		int_t ghost[3] = {2, 2, 0};
		task.das.setParts(globalPt, localPt, ghost);
		
		{
			int_t size[3] = {n.x, n.y, 1};
			int_t globalPt[3] = {task.das.numParts(rgrid::X), task.das.numParts(rgrid::Y), task.das.numParts(rgrid::Z)};
			int_t localPt[3] = {1, 1, 1};
			int_t ghost[3] = {0, 0, 0};
			int_t nc = 1;
			task.dasSave.setSizes(size, globalPt, localPt, ghost, nc);
		}
		
		if (task.das.getInternalRank() == 0) {
			task.da.resize(n.x, n.y, 1, ghost[rgrid::X], ghost[rgrid::Y], ghost[rgrid::Z], (int_t)NODE_SIZE);
		}
		/*
		 * Init material.
		*/
		le_init_material(2.0, 1.0, 1.5, mat);
		
		/*
		 * Init task.
		 */
		le_init_task(task, dt, h, mat, n);
		
		/*
		 * Initial condition.
		 */
		if (task.das.getInternalRank() == 0)
			le_set_ball(task, center, 10.0, 1.0);

		task.das.setAndScatter(0, task.da);
		
		t = timer();
		for (i = 0; i < steps; i++) {
			#ifdef SAVE_EVERY_STEP
				char name[1000];
				sprintf(name, "out-%06d.vtk", i);
				le_save_task(task, name);
			#endif
			le_step(task);
		}
		t = timer() - t;
		
		if (rgmpi::worldRank() == 0) {
			printf("%d %d %d %f\n", n.x, n.y, steps, t);
		}

		/*
		 * Save last step.
		 */
		le_save_task(task, "result.vtk");
		
		/*
		 * Free memory.
		 */
		le_free_task(task);
	}
	rgmpi::forceFinalize();
	return 0;
}

