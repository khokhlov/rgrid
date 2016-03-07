#ifndef DARRAY_SCATTER_H
#define DARRAY_SCATTER_H

#include "rgrid/darraycontainer.h"
#include "rgrid/pdim.h"
#include "rgrid/rgmpi.h"
#include "rgrid/darray.h"

#ifdef USE_MPI
namespace rgrid {

/*
 * DArrayScatter splits DArray into parts, send them to different processes and synchronize them
 * 1. call setSizes() or create by constructor with the same params, params must be equal on all processes
 * 2. call setAndScatter() on one process
 * 3. use getLocalContainer() and fillGhost() on all processes to work with local part
 */
template <typename T, typename I>
class DArrayScatter : public RGCut<I> {
public:
	/* 
	 * Use setSizes() if you create object by this constructor
	 */
	DArrayScatter() {};
	DArrayScatter(I const size[ALL_DIRS], 
	              I const globalPt[ALL_DIRS], 
	              I const localPt[ALL_DIRS], 
	              I const ghost[ALL_DIRS],
	              I const nc)
	{
		setSizes(size, globalPt, localPt, ghost, nc); 
	};
	~DArrayScatter();
	/*
	 * size is number of nodes in each direction
	 * globalPt is number of process parts, must be the same on all processes and 
	 * globalPt[X] * globalPt[Y] * globalPt[Z] == size of processes in MPI_COMM_WORLD
	 * localPt - number of local parts in each direction, can be different on defferent processes
	 * ghost - number of ghost nodes
	 * nc - number of components
	 */
	void setSizes(I const size[ALL_DIRS], 
	              I const globalPt[ALL_DIRS], 
	              I const localPt[ALL_DIRS], 
	              I const ghost[ALL_DIRS],
	              I const nc);
	/*
	 * get local container, use this when DArray scattered
	 */
	DArrayContainer<T, I>& getLocalContainer() { return dac; }
	/*
	 * Scatter DArray from process that calls this functon
	 * Only one process call this function
	 * da - array to scatter
	 */
	void setAndScatter(DArray<T, I> const& da);
	/*
	 * Scatter DArray from process with rank == processRank in cart comm
	 * All processes except one call this function
	 * processRank - rank of process to receive
	 */
	void setAndScatter(int const processRank = MPI_ANY_SOURCE);
	/*
	 * Gather DArray from process that calls this functon
	 * Only one process call this function 
	 */
	void gatherAndGet(DArray<T, I>& da);
	/*
	 * Gather DArray to process with rank == processRank in cart comm
	 * All processes except one call this function
	 * processRank - rank to send parts
	 */
	void gatherAndGet(int processRank);
	/*
	 * fill ghost nodes of local DArrays in local DArrayContainer by 
	 * 1) copying them from other processes
	 * 2) copying them from other DArrays inside local DArrayContainer
	 * 3) filling them with values from boundary nodes
	 */
	void fillGhost();
	/* get rank in internal cart comm */
	int getInternalRank() { return cartRank; }
private:
	
	enum mpiTag {
		SCATTER = 0,
		GATHER
	};
	
	DArrayScatter(DArrayScatter const&);
	DArrayScatter& operator=(DArrayScatter const&);
	
	/* generate types for scatter, gather */
	void generateSubarrayTypes();
	/* release types created by generateSubarrayTypes() */
	void releaseSubarrayTypes(std::vector<MPI_Datatype>& dt);
	
	/* container for local darrays */
	DArrayContainer<T, I> dac;
	/* position of local DArray */
	I cartPos[ALL_DIRS];
	/* cart comm for current decomposition */
	MPI_Comm cartComm;
	/* rank of current process in cart comm */
	int cartRank;
	/* neighbours ranks in cart comm */
	int neigh[SIDE_ALL][ALL_DIRS];
	/* value in neighbours array, if there is no neighbour */
	static const int NO_NEIGHBOUR = -1;
	/* parts in DArrayContainer */
	I localPt[ALL_DIRS];
	/* size of ghost */
	I ghost[ALL_DIRS];
	/* number of componentes */
	I nc;
	/* types for every subarray in big array */
	std::vector<MPI_Datatype> subArrayDt;
	/* types for every subarray itself */
	std::vector<MPI_Datatype> arrayDt;
};

template <typename T, typename I>
void DArrayScatter<T, I>::setSizes(I const size[ALL_DIRS], 
                                   I const globalPt[ALL_DIRS], 
                                   I const localPt[ALL_DIRS], 
                                   I const ghost[ALL_DIRS],
                                   I const nc) 
{
	RGCut<I>::setCutParams(size, globalPt);
	rgmpi::cartCreate(cartComm, globalPt);
	cartRank = rgmpi::commRank(cartComm);
	rgmpi::cartCoords(cartComm, cartRank, cartPos);
	// set neighbours ranks
	for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1)) {
		if (cartPos[d] == 0) {
			neigh[SIDE_LEFT][d] = NO_NEIGHBOUR;
		} else {
			--cartPos[d];
			neigh[SIDE_LEFT][d] = rgmpi::cartRank(cartComm, cartPos);
			++cartPos[d];
		}
		if (cartPos[d] == RGCut<I>::numParts(d) - 1) {
			neigh[SIDE_RIGHT][d] = NO_NEIGHBOUR;
		} else {
			++cartPos[d];
			neigh[SIDE_RIGHT][d] = rgmpi::cartRank(cartComm, cartPos);
			--cartPos[d];
		}
	}
	for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1)) {
		this->localPt[d] = localPt[d];
		this->ghost[d] = ghost[d];
	}
	this->nc = nc;
	generateSubarrayTypes();
}

template <typename T, typename I>
void DArrayScatter<T, I>::setAndScatter(const DArray<T, I>& da) {
	RG_ASSERT(RGCut<I>::numNodes() == da.localSize(), "Wrong number of nodes");
	std::vector<MPI_Request> req;
	req.resize(RGCut<I>::numParts());
	// send parts to all processes
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
	for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
	for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
		I destCoord[ALL_DIRS] = { i, j, k };
		int destRank = rgmpi::cartRank(cartComm, destCoord);
		MPI_CHECK(MPI_Isend(&da[0], 1, subArrayDt[destRank], destRank, SCATTER, cartComm, &req[destRank]));
	}
	// receive own part
	setAndScatter(cartRank);
	// wait while all parts will be received
	MPI_CHECK(MPI_Waitall(RGCut<I>::numParts(), &req[0], MPI_STATUSES_IGNORE));
}

template <typename T, typename I>
void DArrayScatter<T, I>::setAndScatter(int const processRank) {
	DArray<T, I> da; // temporary darray
	da.resize(RGCut<I>::numNodes(X),
	          RGCut<I>::numNodes(Y),
	          RGCut<I>::numNodes(Z),
	          RGCut<I>::partNodes(X, cartPos[X]), 
	          RGCut<I>::partNodes(Y, cartPos[Y]), 
	          RGCut<I>::partNodes(Z, cartPos[Z]),
	          RGCut<I>::partOrigin(X, cartPos[X]),
	          RGCut<I>::partOrigin(Y, cartPos[Y]),
	          RGCut<I>::partOrigin(Z, cartPos[Z]),
	          ghost[X], ghost[Y], ghost[Z]);
	da.alloc(nc);
	MPI_CHECK(MPI_Recv(&da[0], 1, arrayDt[cartRank], processRank, SCATTER, cartComm, MPI_STATUS_IGNORE));
	dac.setDArray(da, localPt[X], localPt[Y], localPt[Z]);
}

template <typename T, typename I>
void DArrayScatter<T, I>::generateSubarrayTypes() {
	releaseSubarrayTypes(subArrayDt);
	releaseSubarrayTypes(arrayDt);
	subArrayDt.resize(RGCut<I>::numParts());
	arrayDt.resize(RGCut<I>::numParts());
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
	for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
	for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
		int cartRank3[ALL_DIRS] = { i, j, k };
		int rank = rgmpi::cartRank(cartComm, cartRank3);
		{
			// generate subArrayDt
			int array_of_sizes[4] = { RGCut<I>::numNodes(X) + 2 * ghost[X], 
			                          RGCut<I>::numNodes(Y) + 2 * ghost[Y], 
			                          RGCut<I>::numNodes(Z) + 2 * ghost[Z], 
			                          nc };
			int array_of_subsizes[4] = { RGCut<I>::partNodes(X, i),
			                             RGCut<I>::partNodes(Y, j),
			                             RGCut<I>::partNodes(Z, k),
			                             nc };
			int array_of_starts[4] = { ghost[X] + RGCut<I>::partOrigin(X, i), 
			                           ghost[Y] + RGCut<I>::partOrigin(Y, j),
			                           ghost[Z] + RGCut<I>::partOrigin(Z, k),
			                           0 };
			subArrayDt[rank] = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
		}
		{
			// generate arrayDt
			int array_of_sizes[4] = { RGCut<I>::partNodes(X, i) + 2 * ghost[X],
			                          RGCut<I>::partNodes(Y, j) + 2 * ghost[Y], 
			                          RGCut<I>::partNodes(Z, k) + 2 * ghost[Z], 
			                          nc };
			int array_of_subsizes[4] = { RGCut<I>::partNodes(X, i),
			                             RGCut<I>::partNodes(Y, j),
			                             RGCut<I>::partNodes(Z, k),
			                             nc };
			int array_of_starts[4] = { ghost[X], ghost[Y], ghost[Z], 0 };
			arrayDt[rank] = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
		}
	}
}

template <typename T, typename I>
void DArrayScatter<T, I>::releaseSubarrayTypes(std::vector<MPI_Datatype>& dt) {
	for(std::vector<MPI_Datatype>::iterator it = dt.begin(); it != dt.end(); ++it) {
		rgmpi::freeSubarrayType(*it);
	}
	dt.resize(0);
}

template<typename T, typename I>
DArrayScatter<T, I>::~DArrayScatter() {
	releaseSubarrayTypes(subArrayDt);
	releaseSubarrayTypes(arrayDt);
}

template<typename T, typename I>
void DArrayScatter<T, I>::gatherAndGet(DArray<T, I>& da ) {
	da.resize(RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
	          RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
	          0, 0, 0,
	          ghost[X], ghost[Y], ghost[Z]);
	da.alloc(nc);
	std::vector<MPI_Request> req;
	req.resize(RGCut<I>::numParts());
	// recv parts from all processes
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
	for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
	for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
		I srcCoord[ALL_DIRS] = { i, j, k };
		int srcRank = rgmpi::cartRank(cartComm, srcCoord);
		MPI_CHECK(MPI_Irecv(&da[0], 1, subArrayDt[srcRank], srcRank, GATHER, cartComm, &req[srcRank]));
	}
	// gather own part
	gatherAndGet(cartRank);
	// wait while all parts will be received
	MPI_CHECK(MPI_Waitall(RGCut<I>::numParts(), &req[0], MPI_STATUSES_IGNORE));
}

template<typename T, typename I>
void DArrayScatter<T, I>::gatherAndGet(int processRank) {
	DArray<T, I> da;
	dac.getDArray(da);
	MPI_CHECK(MPI_Send(&da[0], 1, arrayDt[cartRank], processRank, GATHER, cartComm));
}

} /* namespace rgrid */

#endif // USE_MPI

#endif /* DARRAY_SCATTER_H */
