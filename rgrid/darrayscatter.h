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
 */
template <typename T, typename I>
class DArrayScatter : public RGCut<I> {
public:
	/* 
	 * Use setSizes() if you create object by this constructor
	 */
	DArrayScatter() {
		this->init();
	};
	DArrayScatter(I const size[ALL_DIRS], 
	              I const globalPt[ALL_DIRS], 
	              I const localPt[ALL_DIRS], 
	              I const ghost[ALL_DIRS],
	              I const nc)
	{
		this->init();
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
	 * Scatter DArray from process with rank == processRank in cart comm
	 * All processes in cart comm call this function
	 * da - array to scatter
	 * processes with rank != processRank set dummy DArray in da that remains unchanged
	 */
	void setAndScatter(int const processRank, DArray<T, I> const& da);
	/*
	 * Gather DArray to process with rank == processRank in cart comm
	 * All processes in cart comm call this function
	 * process with rank == processRank receive full DArray in da,
	 * da on other processes remains unchanged
	 */
	void gatherAndGet(int const processRank, DArray<T, I>& da);
	/*
	 * fill ghost nodes of local DArrays in local DArrayContainer by 
	 * copying them from other DArrays inside local DArrayContainer
	 */
	void internalSync();
	/*
	 * start to receive ghost ghost for local DArrays in local
	 * DArrayContainer from other processes
	 */
	void externalSyncStart();
	/*
	 * wait until all data will be received and written in ghost nodes of local DArrays
	 */
	void externalSyncEnd();
	/*
	 * fill ghost nodes of local DArrays in local DArrayContainer with 
	 * values from boundary nodes
	 */
	void fillGhost();
	/* get rank in internal cart comm */
	int getInternalRank() const { return cartRank; }
	/* get number of components */
	I getNC() const { return nc; }
	/* get number of ghost nodes */
	I getGhost(CartDir d) const { return ghost[d]; }
	/* get type for file view in IO operations */
	MPI_Datatype fileViewType() const { return viewType; }
	/* return local DArrayContainer */
	DArrayContainer<T, I>& getDAC() { return dac; }
	/* get type for darray with indexes x, y, z */
	MPI_Datatype getDArrayDt(I x, I y, I z) const { return arrayDt.at(RGCut<I>::linInd(x, y, z)); }
	/* get DArray with indexes x, y, z in local DArrayContainer */
	DArray<T, I>& getDArrayPart(I x, I y, I z) { return dac.getDArrayPart(x, y, z); }
	/* get buffer of DArray with indexes x, y, z in local DArrayContainer */
	T* getDArrayBuffer(I x, I y, I z) { return &dac.getDArrayPart(x, y, z)[0]; }
	/* start saving data to file */
	void saveDataBegin(const std::string filename, const rgio::format fmt) {
		std::stringstream ss;
		I size[ALL_DIRS] = { RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z) };
		rgio::writeHeader(ss, size, getNC(), fmt);
		std::string header = ss.str();
		// write header
		if (rgmpi::commRank(cartComm) == 0) {
			MPI_CHECK(MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
			MPI_CHECK(MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL));
			MPI_CHECK(MPI_File_write(fh, header.c_str(), header.size(), MPI_CHAR, MPI_STATUS_IGNORE));
			MPI_CHECK(MPI_File_close(&fh));
		}
		// write data
		MPI_CHECK(MPI_File_open(cartComm, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
		MPI_CHECK(MPI_File_set_view(fh, header.size() * sizeof(char), rgmpi::getMPItype<T>(), fileViewType(), "native", MPI_INFO_NULL));
		MPI_Request req;
		saveDataRequest.clear();
		saveDataRequest.reserve(RGCut<I>::numParts());
		for (I k = 0; k != dac.numParts(Z); ++k)
		for (I j = 0; j != dac.numParts(Y); ++j)
		for (I i = 0; i != dac.numParts(X); ++i) {
			MPI_CHECK(MPI_File_iwrite(fh, &dac.getDArrayPart(i, j, k)[0], 1, locArrayDt.at(dac.linInd(i, j, k)), &req));
			saveDataRequest.push_back(req);
		}
	}
	/* Wait while all data will be saved */
	void saveDataEnd() {
		MPI_CHECK(MPI_Waitall(saveDataRequest.size(), &saveDataRequest.at(0), MPI_STATUSES_IGNORE));
		saveDataRequest.clear();
		MPI_CHECK(MPI_File_close(&fh));
	}
private:
	MPI_File fh;
	std::vector<MPI_Request> saveDataRequest;
	std::vector<MPI_Datatype> locArrayDt;
	
private:
	
	enum mpiTag {
		SCATTER = 0,
		GATHER,
		SYNC
	};
	
	DArrayScatter(DArrayScatter const&);
	DArrayScatter& operator=(DArrayScatter const&);
	
	/* initialize common things for all constructors */
	void init() {
		viewType = MPI_DATATYPE_NULL;
	}
	
	/* generate types for scatter, gather */
	void generateSubarrayTypes();
	/* release MPI subarray types */
	void releaseSubarrayTypes(std::vector<MPI_Datatype>& dt);
	/* send boundary nodes to negihbours */
	void sendNodesStart(int rank, CartSide side, CartDir dir);
	void sendNodesEnd();
	
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
	/* external sync buffers */
	std::vector<T> sendBuf[SIDE_ALL][ALL_DIRS];
	std::vector<T> recvBuf[SIDE_ALL][ALL_DIRS];
	/* request for external sync start and end */
	MPI_Request syncSendReq[SIDE_ALL * ALL_DIRS];
	MPI_Request syncRecvReq[SIDE_ALL * ALL_DIRS];
	/* value in neighbours array, if there is no neighbour */
	static const I NO_NEIGHBOUR = -1;
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
	/* view type */
	MPI_Datatype viewType;
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
void DArrayScatter<T, I>::setAndScatter(int const processRank, DArray<T, I> const& da) {
	std::vector<MPI_Request> req;
	if (processRank == cartRank) {
		RG_ASSERT(RGCut<I>::numNodes() == da.localSize(), "Wrong number of nodes");
		req.resize(RGCut<I>::numParts());
		// send parts to all processes
		for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
		for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
			I destCoord[ALL_DIRS] = { i, j, k };
			int destRank = rgmpi::cartRank(cartComm, destCoord);
			MPI_CHECK(MPI_Isend(&da[0], 1, subArrayDt[destRank], destRank, SCATTER, cartComm, &req[destRank]));
		}
	}
	
	DArray<T, I> tda; // temporary darray
	tda.resize(
		RGCut<I>::numNodes(X),
		RGCut<I>::numNodes(Y),
		RGCut<I>::numNodes(Z),
		RGCut<I>::partNodes(X, cartPos[X]), 
		RGCut<I>::partNodes(Y, cartPos[Y]), 
		RGCut<I>::partNodes(Z, cartPos[Z]),
		RGCut<I>::partOrigin(X, cartPos[X]),
		RGCut<I>::partOrigin(Y, cartPos[Y]),
		RGCut<I>::partOrigin(Z, cartPos[Z]),
		ghost[X], ghost[Y], ghost[Z]);
	tda.alloc(nc);
	MPI_CHECK(MPI_Recv(&tda[0], 1, arrayDt[cartRank], processRank, SCATTER, cartComm, MPI_STATUS_IGNORE));
	dac.setDArray(tda, localPt[X], localPt[Y], localPt[Z]);
	
	if (processRank == cartRank) {
		// wait while all parts will be received
		MPI_CHECK(MPI_Waitall(RGCut<I>::numParts(), &req[0], MPI_STATUSES_IGNORE));
	}
}

template <typename T, typename I>
void DArrayScatter<T, I>::generateSubarrayTypes() {
	releaseSubarrayTypes(subArrayDt);
	releaseSubarrayTypes(arrayDt);
	releaseSubarrayTypes(locArrayDt);
	subArrayDt.resize(RGCut<I>::numParts());
	arrayDt.resize(RGCut<I>::numParts());
	for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
	for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
	for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
		int cartRank3[ALL_DIRS] = { i, j, k };
		int rank = rgmpi::cartRank(cartComm, cartRank3);
		{
			// generate subArrayDt
			int array_of_sizes[4] = { 
				RGCut<I>::numNodes(X) + 2 * ghost[X], 
				RGCut<I>::numNodes(Y) + 2 * ghost[Y], 
				RGCut<I>::numNodes(Z) + 2 * ghost[Z], 
				nc };
			int array_of_subsizes[4] = { 
				RGCut<I>::partNodes(X, i),
				RGCut<I>::partNodes(Y, j),
				RGCut<I>::partNodes(Z, k),
				nc };
			int array_of_starts[4] = { 
				ghost[X] + RGCut<I>::partOrigin(X, i), 
				ghost[Y] + RGCut<I>::partOrigin(Y, j),
				ghost[Z] + RGCut<I>::partOrigin(Z, k),
				0 };
			subArrayDt[rank] = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
		}
		{
			// generate arrayDt
			int array_of_sizes[4] = { 
				RGCut<I>::partNodes(X, i) + 2 * ghost[X],
				RGCut<I>::partNodes(Y, j) + 2 * ghost[Y], 
				RGCut<I>::partNodes(Z, k) + 2 * ghost[Z], 
				nc };
			int array_of_subsizes[4] = { 
				RGCut<I>::partNodes(X, i),
				RGCut<I>::partNodes(Y, j),
				RGCut<I>::partNodes(Z, k),
				nc };
			int array_of_starts[4] = { ghost[X], ghost[Y], ghost[Z], 0 };
			arrayDt[rank] = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
		}
	}
	// make view type
	int array_of_sizes[4] = { 
		RGCut<I>::numNodes(X), 
		RGCut<I>::numNodes(Y), 
		RGCut<I>::numNodes(Z), 
		nc
	};
	int array_of_subsizes[4] = { 
		RGCut<I>::partNodes(X, cartPos[X]),
		RGCut<I>::partNodes(Y, cartPos[Y]),
		RGCut<I>::partNodes(Z, cartPos[Z]),
		nc
	};
	int array_of_starts[4] = {
		RGCut<I>::partOrigin(X, cartPos[X]), 
		RGCut<I>::partOrigin(Y, cartPos[Y]),
		RGCut<I>::partOrigin(Z, cartPos[Z]),
		0
	};
	viewType = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
	// type for local DArrays
	locArrayDt.resize(dac.numParts());
	for (I k = 0; k != dac.numParts(Z); ++k)
	for (I j = 0; j != dac.numParts(Y); ++j)
	for (I i = 0; i != dac.numParts(X); ++i) {
		// generate subArrayDt
		int array_of_sizes[4] = { 
			dac.numNodes(X), 
			dac.numNodes(Y),
			dac.numNodes(Z), 
			nc };
		int array_of_subsizes[4] = { 
			dac.partNodes(X, i),
			dac.partNodes(Y, j),
			dac.partNodes(Z, k),
			nc };
		int array_of_starts[4] = { 
			dac.partOrigin(X, i), 
			dac.partOrigin(Y, j),
			dac.partOrigin(Z, k),
			0 };
		locArrayDt.at(dac.linInd(i, j, k)) = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
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
	releaseSubarrayTypes(locArrayDt);
	if (viewType != MPI_DATATYPE_NULL) {
		rgmpi::freeSubarrayType(viewType);
	}
}

template<typename T, typename I>
void DArrayScatter<T, I>::gatherAndGet(int processRank, DArray<T, I>& da) {
	std::vector<MPI_Request> req;
	if (processRank == cartRank) {
		da.resize(
			RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
			RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
			0, 0, 0,
			ghost[X], ghost[Y], ghost[Z]);
		da.alloc(nc);
		req.resize(RGCut<I>::numParts());
		// recv parts from all processes
		for (I k = 0; k != RGCut<I>::numParts(Z); ++k)
		for (I j = 0; j != RGCut<I>::numParts(Y); ++j)
		for (I i = 0; i != RGCut<I>::numParts(X); ++i) {
			I srcCoord[ALL_DIRS] = { i, j, k };
			int srcRank = rgmpi::cartRank(cartComm, srcCoord);
			MPI_CHECK(MPI_Irecv(&da[0], 1, subArrayDt[srcRank], srcRank, GATHER, cartComm, &req[srcRank]));
		}
	}
	
	DArray<T, I> tda;
	dac.getDArray(tda);
	MPI_CHECK(MPI_Send(&tda[0], 1, arrayDt[cartRank], processRank, GATHER, cartComm));
	
	if (processRank == cartRank) {
		// wait while all parts will be received
		MPI_CHECK(MPI_Waitall(RGCut<I>::numParts(), &req[0], MPI_STATUSES_IGNORE));
	}
}

template<typename T, typename I>
void DArrayScatter<T, I>::internalSync() {
	dac.sync();
}

template<typename T, typename I>
void DArrayScatter<T, I>::fillGhost() {
	for (CartSide s = SIDE_LEFT; s != SIDE_ALL; s = static_cast<CartSide>(s + 1))
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1)) 
			if (neigh[s][d] == NO_NEIGHBOUR) {
				CartDir ort1, ort2;
				ortDirs(d, ort1, ort2);
				I k = (s == SIDE_LEFT) ? 0 : RGCut<I>::numParts(d) - 1;
				for (I i = 0; i < RGCut<I>::numParts(ort1); i++) 
				for (I j = 0; j < RGCut<I>::numParts(ort2); j++) {
					I coord[ALL_DIRS];
					coord[d] = k;
					coord[ort1] = i;
					coord[ort2] = j;
					I x = coord[X], y = coord[Y], z = coord[Z];
					DArray<T, I> da = dac.getDArrayPart(x, y, z);
					da.fillGhost(d, s);
				}
			}
}

template<typename T, typename I>
void DArrayScatter<T, I>::externalSyncStart() {
	for (CartSide s = SIDE_LEFT; s != SIDE_ALL; s = static_cast<CartSide>(s + 1))
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1)) 
			if (neigh[s][d] != NO_NEIGHBOUR) {
				CartDir ort1, ort2;
				ortDirs(d, ort1, ort2);
				I bufSize = RGCut<I>::partNodes(ort1) * RGCut<I>::partNodes(ort2) * ghost[d];
				// start recv
				recvBuf[s][d].resize(bufSize);
				MPI_CHECK(MPI_Irecv(&recvBuf[s][d].front(), bufSize, rgmpi::getMPItype<T>(), neigh[s][d], SYNC, cartComm, &syncRecvReq[s + d * SIDE_ALL]));
				// start send
				// sub array size and position
				I orig[ALL_DIRS] = { 0, 0, 0 };
				I size[ALL_DIRS] = {
					RGCut<I>::partNodes(X, cartPos[X]),
					RGCut<I>::partNodes(Y, cartPos[Y]), 
					RGCut<I>::partNodes(Z, cartPos[Z]),
				};
				size[d] = ghost[d];
				orig[d] = (s == SIDE_LEFT) ? 0 : RGCut<I>::numNodes(d) - ghost(d);
				dac.getSubArray(orig[X], orig[Y], orig[Z], size[X], size[Y], size[Z], sendBuf[s][d]);
				MPI_CHECK(MPI_Isend(&sendBuf[s][d].front(), bufSize, rgmpi::getMPItype<T>(), neigh[s][d], SYNC, cartComm, &syncSendReq[s + d * SIDE_ALL]));
			} else {
				syncRecvReq[s + d * SIDE_ALL] = MPI_REQUEST_NULL;
			}
}

template<typename T, typename I>
void DArrayScatter<T, I>::externalSyncEnd() {
	while (1) {
		int index;
		MPI_CHECK(MPI_Waitany(SIDE_ALL * ALL_DIRS, syncRecvReq, &index, MPI_STATUS_IGNORE));
		if (index == MPI_UNDEFINED) break;
		CartSide s = index % SIDE_ALL;
		CartDir d = index / ALL_DIRS;
		CartDir ort1, ort2;
		ortDirs(d, ort1, ort2);
		// sub array size and position
		I orig[ALL_DIRS] = { 0, 0, 0 };
		I size[ALL_DIRS] = {
			RGCut<I>::partNodes(X, cartPos[X]),
			RGCut<I>::partNodes(Y, cartPos[Y]),
			RGCut<I>::partNodes(Z, cartPos[Z]),
		};
		size[d] = ghost[d];
		orig[d] = (s == SIDE_LEFT) ? -ghost[d] : RGCut<I>::numNodes(d);
		dac.setSubArray(orig[X], orig[Y], orig[Z], size[X], size[Y], size[Z], recvBuf[s][d]);
	}
	MPI_CHECK(MPI_Waitall(SIDE_ALL * ALL_DIRS, syncSendReq, MPI_STATUSES_IGNORE));
}

} /* namespace rgrid */

#endif // USE_MPI

#endif /* DARRAY_SCATTER_H */
