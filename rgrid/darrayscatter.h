/**
 * \file
 * \brief Two partitioning steps (local and global) of DArray
 */

#ifndef DARRAY_SCATTER_H
#define DARRAY_SCATTER_H

#include <cstddef>

#include "rgrid/darraycontainer.h"
#include "rgrid/pdim.h"
#include "rgrid/rgmpi.h"
#include "rgrid/darray.h"

#ifdef USE_MPI
namespace rgrid {

/**
 * \brief DArrayScatter represents grid with two partitioning steps
 * 
 * First step is global partitioning - MPI.
 * DArrayScatter contains one DArrayContainer on each process.
 * 
 * Second step is local partitioning - OpenMP, OpenCL, etc.
 * Each local DArrayContainer contains DArrays.
 * 
 * \tparam T type of every grid node (i.e. double, float)
 * \tparam I type of grid indexes (i.e. int, long)
 */
template <typename T, typename I>
class DArrayScatter : public RGCut<I> {
public:
	/**
	 * \brief Default constructor, one global and one local part
	 * 
	 * To make more parts call setParts() before loading grid into DArrayScatter
	 */
	DArrayScatter() : cartComm(MPI_COMM_NULL), viewType(MPI_DATATYPE_NULL) 
	{
		I const size[ALL_DIRS] = {1, 1, 1};
		I const globalPt[ALL_DIRS] = {1, 1, 1};
		I const localPt[ALL_DIRS] = {1, 1, 1};
		I const ghost[ALL_DIRS] = {0, 0, 0};
		I const nc = 1;
		setSizes(size, globalPt, localPt, ghost, nc);
	};
	
	/**
	 * \brief The most verbose constructor: grid size and partitioning
	 * 
	 * \note It's better to use DArrayScatter(I const globalPt[ALL_DIRS], I const localPt[ALL_DIRS], I const ghost[ALL_DIRS]) instead
	 * \param[in] size number of nodes of entire grid
	 * \param[in] globalPt number of global parts (MPI processes in each direction)
	 * \param[in] localPt number of local parts (DArray's in local DArrayContainer), can be different on different MPI processes
	 * \param[in] ghost number of ghost nodes on each side of grid
	 * \param[in] nc number of type T variables in each grid node
	 */
	DArrayScatter(I const size[ALL_DIRS],
	              I const globalPt[ALL_DIRS],
	              I const localPt[ALL_DIRS],
	              I const ghost[ALL_DIRS],
	              I const nc) 
		: cartComm(MPI_COMM_NULL), viewType(MPI_DATATYPE_NULL)
	{
		setSizes(size, globalPt, localPt, ghost, nc);
	};
	
	/**
	 * \brief Constructor with grid partitioning
	 * 
	 * \param[in] globalPt number of global parts (MPI processes in each direction)
	 * \param[in] localPt number of local parts (DArray's in local DArrayContainer), can be different on different MPI processes
	 * \param[in] ghost number of ghost nodes on each side of grid
	 */
	DArrayScatter(I const globalPt[ALL_DIRS],
	              I const localPt[ALL_DIRS],
	              I const ghost[ALL_DIRS])
		: cartComm(MPI_COMM_NULL), viewType(MPI_DATATYPE_NULL)
	{
		setParts(globalPt, localPt, ghost);
	}
	
	~DArrayScatter();
	
	/**
	 * \brief Set grid partitioning
	 * 
	 * Call to this function causes creation of MPI communicator with cartesian topology
	 * with dimensions specified in globalPt. All processes in MPI_COMM_WORLD call this function.
	 * 
	 * \param[in] globalPt number of global parts (MPI processes in each direction)
	 * \param[in] localPt number of local parts (DArray's in local DArrayContainer), can be different on different MPI processes
	 * \param[in] ghost number of ghost nodes on each side of grid
	 */
	void setParts(I const globalPt[ALL_DIRS],
	              I const localPt[ALL_DIRS],
	              I const ghost[ALL_DIRS]) {
		I size[3]; // stub
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1)) {
			size[d] = globalPt[d] * localPt[d];
		}
		setSizes(size, globalPt, localPt, ghost, nc);
	}
	
	/**
	 * \brief Set grid size and partitioning
	 * 
	 * \note It's better to use setParts(I const globalPt[ALL_DIRS], I const localPt[ALL_DIRS], I const ghost[ALL_DIRS]) instead
	 * \param[in] size number of nodes of entire grid
	 * \param[in] globalPt number of global parts (MPI processes in each direction)
	 * \param[in] localPt number of local parts (DArray's in local DArrayContainer), can be different on different MPI processes
	 * \param[in] ghost number of ghost nodes on each side of grid
	 * \param[in] nc number of type T variables in each grid node
	 */
	void setSizes(I const size[ALL_DIRS],
	              I const globalPt[ALL_DIRS],
	              I const localPt[ALL_DIRS],
	              I const ghost[ALL_DIRS],
	              I const nc);
	/**
	 * \brief Get local part of grid
	 * 
	 * Get local container when grid loaded
	 * 
	 * \return DArrayContainer which represents local part of entire grid
	 */
	DArrayContainer<T, I> &getLocalContainer() {
		return dac;
	}
	/**
	 * \brief Scatter DArray from one process to all
	 * 
	 * Process which rank == processRank divides own DArray da into parts and
	 * sends them to other processes. 
	 * All processes call this function, but processes with rank != processRank 
	 * set dummy DArray in da which is not used.
	 * 
	 * \note Load grid by loadDataBegin when there are 
	 * not enough memory to save entire grid on one process
	 * 
	 * \param[in] processRank rank of process to scatter from
	 * \param[in] da DArray to scatter
	 */
	void setAndScatter(int const processRank, DArray<T, I> const &da);
	/**
	 * \brief Gather DArray from all processes to one
	 * 
	 * Gather DArray to process with rank == processRank 
	 * All processes call this function. Process with rank == processRank 
	 * receive full DArray in da, da on other processes remains unchanged.
	 * 
	 * \param[in] processRank rank of process to gather to
	 * \param[out] da DArray to gather
	 */
	void gatherAndGet(int const processRank, DArray<T, I> &da);
	/**
	 * \brief Synchronization inside local DArrayContainer
	 * 
	 * Fill ghost nodes of local DArrays in local DArrayContainer by
	 * copying them from other DArrays inside local DArrayContainer
	 */
	void internalSync();
	/**
	 * \brief Synchronization between DArrayContainers on different processes
	 * 
	 * Start to receive ghost nodes for local DArrays in local
	 * DArrayContainer from other processes
	 * 
	 * \warning Don't use ghost nodes before call to rgrid::DArrayScatter< T, I >::externalSyncEnd()
	 * 
	 * \sa externalSyncEnd
	 */
	void externalSyncStart();
	/**
	 * \brief Synchronization between DArrayContainers on different processes
	 * 
	 * Wait until all data will be received and written in ghost nodes of local DArrays
	 * 
	 * \sa externalSyncStart
	 */
	void externalSyncEnd();
	/**
	 * \brief Fill local ghost nodes with values from boundary nodes
	 * 
	 * Local ghost nodes on the edge of entire grid will be filled.
	 */
	void fillGhost();
	/**
	 * \brief Get rank in internal communicator
	 * \return rank of current process in internal cartesian communicator
	 */
	int getInternalRank() const {
		return cartRank;
	}
	/**
	 * \brief Get position in internal communicator
	 * \param[in] d direction
	 * \return position of current process in cartesian communicator in direction d
	 */
	int getInternalPos(CartDir d) {
		return cartPos[d];
	}
	/**
	 * \brief Get number of components
	 * \return number of components in each node
	 */
	I getNC() const {
		return nc;
	}
	/**
	 * \brief Get number of ghost nodes
	 * \param[in] d direction
	 * \return number of ghost nodes in direction d
	 */
	I getGhost(CartDir d) const {
		return ghost[d];
	}
	/**
	 * \brief Get type for file view in MPI IO operations
	 * \return file view type
	 */
	MPI_Datatype fileViewType() const {
		return viewType;
	}
	/**
	 * \brief Get MPI type for DArray.
	 * \details This type is subarray type of DArray with coordinates (x, y, z) in entire grid.
	 * \return type for darray with indexes (x, y, z)
	 */
	MPI_Datatype getDArrayDt(I x, I y, I z) const {
		return arrayDt.at(RGCut<I>::linInd(x, y, z));
	}
	/**
	 * \brief Get DArray from local DarrayContainer
	 * \return DArray with indexes (x, y, z) in local DArrayContainer 
	 */
	DArray<T, I> &getDArrayPart(I x, I y, I z) {
		return dac.getDArrayPart(x, y, z);
	}
	/**
	 * \brief Get buffer of DArray with indexes (x, y, z) in local DArrayContainer
	 * \return pointer to raw data
	 */
	T *getDArrayBuffer(I x, I y, I z) {
		return &dac.getDArrayPart(x, y, z)[0];
	}
	/**
	 * \brief Start saving entire DArray to file
	 * \param[in] filename name of output file
	 * \param[in] fmt format of file
	 * \param[in] ss is used when fmt == CUSTOM_HEADER
	 * \todo support for TEXT format
	 * 
	 * Use this function to save data from computational DArrayScatter when: 
	 * <br> 1. you don't have to reverse bytes
	 * <br> 2. all DArray components have to be saved
	 * 
	 * In all other cases you should create new DArrayScatter with no local partitioning,
	 * create data to be saved from computational DArrayScatter and place to new. 
	 * During this process chose components to sace convert bytes.
	 * Next, save data from new DArrayScatter by calling this function.
	 * 
	 * \warning Don't modify DArray data before call to rgrid::DArrayScatter< T, I >::saveDataEnd
	 * \sa saveDataEnd
	 */
	void saveDataBegin(std::string filename, const rgio::format fmt, std::stringstream* ss = NULL) {
		if (cartComm == MPI_COMM_NULL) return;
		Dim3D<I> size(RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z));
		if (fmt != rgio::CUSTOM_HEADER) {
			ss = new std::stringstream();
			rgio::writeHeader(*ss, size, getNC(), fmt);
		}
		std::string header = ss->str();
		if (fmt != rgio::CUSTOM_HEADER) {
			delete ss;
		}
		// write header
		if (rgmpi::commRank(cartComm) == 0) {
			MPI_CHECK(MPI_File_open(MPI_COMM_SELF, filename.c_str(), MPI_MODE_CREATE | MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
			MPI_CHECK(MPI_File_set_size(fh, 0));
			MPI_CHECK(MPI_File_set_view(fh, 0, MPI_CHAR, MPI_CHAR, "native", MPI_INFO_NULL));
			MPI_CHECK(MPI_File_write(fh, header.c_str(), header.size(), MPI_CHAR, MPI_STATUS_IGNORE));
			MPI_CHECK(MPI_File_close(&fh));
		}
		// write data
		MPI_CHECK(MPI_File_open(cartComm, filename.c_str(), MPI_MODE_WRONLY, MPI_INFO_NULL, &fh));
		MPI_CHECK(MPI_File_set_view(fh, header.size() * sizeof(char), rgmpi::getMPItype<T>(), fileViewType(), "native", MPI_INFO_NULL));
		if (dac.numParts() == 1) {
			DArray<T, I>& singleDA = dac.getDArrayPart(0);
			MPI_CHECK(MPI_File_write_all_begin(fh, singleDA.getDataRaw(), 1, arrayDt.at(cartRank)));
		} else {
			// use slow implementation, make a copy of local dac in single DArray
			dac.getDArray(tda);
			MPI_CHECK(MPI_File_write_all_begin(fh, tda.getDataRaw(), 1, arrayDt.at(cartRank)));
		}
	}
	/** 
	 * \brief Wait while all data will be saved
	 * \sa saveDataBegin
	 */
	void saveDataEnd() {
		if (cartComm == MPI_COMM_NULL) return;
		if (dac.numParts() == 1) {
			MPI_CHECK(MPI_File_write_all_end(fh, dac.getDArrayPart(0).getDataRaw(), MPI_STATUS_IGNORE));
		} else {
			MPI_CHECK(MPI_File_write_all_end(fh, tda.getDataRaw(), MPI_STATUS_IGNORE));
		}
		MPI_CHECK(MPI_File_close(&fh));
	}
	/** 
	 * \brief Start loading data from file 
	 * \todo support for formats other than BINARY
	 * \param[in] filename name of input file
	 * \warning Don't modify grid data before call to rgrid::DArrayScatter< T, I >::loadDataEnd
	 * \sa loadDataEnd
	 */
	void loadDataBegin(std::string filename) {
		if (cartComm == MPI_COMM_NULL) return;
		std::fstream fs;
		struct newParams {
			I size[ALL_DIRS];
			I components;
			long offset;
			int format;
		} np;
		// read header by first process
		if (rgmpi::commRank(cartComm) == 0) {
			rgio::format fmt;
			Dim3D<I> size;
			np.offset = rgio::loadHeaderMPI(filename, size, np.components, fmt);
			np.size[X] = size[X];
			np.size[Y] = size[Y];
			np.size[Z] = size[Z];
			np.format = fmt;
		}
		// send header info to other processes
		const int nitems = 4;
		const int blocklengths[nitems] = {ALL_DIRS, 1, 1, 1};
		const MPI_Datatype types[nitems] = {rgmpi::getMPItype<I>(), rgmpi::getMPItype<I>(), MPI_LONG, MPI_INT};
		MPI_Datatype mpi_np_type;
		MPI_Aint offsets[nitems] = { offsetof(newParams, size),
		                             offsetof(newParams, components),
		                             offsetof(newParams, offset),
		                             offsetof(newParams, format) };
		
		MPI_CHECK(MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_np_type));
		MPI_CHECK(MPI_Type_commit(&mpi_np_type));
		MPI_CHECK(MPI_Bcast(&np, 1, mpi_np_type, 0, cartComm));
		MPI_CHECK(MPI_Type_free(&mpi_np_type));
		// resize
		const I glPt[ALL_DIRS] = { RGCut<I>::numParts(X), RGCut<I>::numParts(Y), RGCut<I>::numParts(Z) };
		setSizes(np.size, glPt, localPt, ghost, np.components);
		// read data
		MPI_CHECK(MPI_File_open(cartComm, filename.c_str(), MPI_MODE_RDONLY, MPI_INFO_NULL, &fh));
		MPI_CHECK(MPI_File_set_view(fh, np.offset * sizeof(char), rgmpi::getMPItype<T>(), fileViewType(), "native", MPI_INFO_NULL));
		tda.resize(
			RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
			RGCut<I>::partNodes(X, cartPos[X]), RGCut<I>::partNodes(Y, cartPos[Y]), RGCut<I>::partNodes(Z, cartPos[Z]),
			RGCut<I>::partOrigin(X, cartPos[X]), RGCut<I>::partOrigin(Y, cartPos[Y]), RGCut<I>::partOrigin(Z, cartPos[Z]),
			ghost[X], ghost[Y], ghost[Z],
			nc
		);
		MPI_CHECK(MPI_File_read_all_begin(fh, tda.getDataRaw(), 1, arrayDt.at(cartRank)));
		
	}
	/**
	 * \brief Wait while all data will be loaded
	 * \sa loadDataBegin
	 */
	void loadDataEnd() {
		if (cartComm == MPI_COMM_NULL) return;
		MPI_CHECK(MPI_File_read_all_end(fh, tda.getDataRaw(), MPI_STATUS_IGNORE));
		dac.setDArray(tda, localPt[X], localPt[Y], localPt[Z]);
		MPI_CHECK(MPI_File_close(&fh));
		tda.resize(0,0,0, 0,0,0, 0,0,0, 0,0,0, 0); // free memory
	}
private:
	MPI_File fh;
	// datatypes for local array
	std::vector<MPI_Datatype> locArrayDt;

private:

	enum mpiTag {
		SCATTER = 0,
		GATHER,
		SYNC
	};

	DArrayScatter(DArrayScatter const &);
	DArrayScatter &operator=(DArrayScatter const &);

	/* generate types for scatter, gather */
	void generateSubarrayTypes();
	/* release MPI subarray types */
	void releaseSubarrayTypes(std::vector<MPI_Datatype> &dt);
	/* send boundary nodes to negihbours */
	void sendNodesStart(int rank, CartSide side, CartDir dir);
	void sendNodesEnd();

	/* container for local darrays */
	DArrayContainer<T, I> dac;
	/* position of local DArrayContainer */
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
	/* RGCut for local DArrayContainer */
	RGCut<I> dacCut;
	/* temp DArray for save/load */
	DArray<T, I> tda;
};

template <typename T, typename I>
void DArrayScatter<T, I>::setSizes(I const size[ALL_DIRS],
                                   I const globalPt[ALL_DIRS],
                                   I const localPt[ALL_DIRS],
                                   I const ghost[ALL_DIRS],
                                   I const nc) {
	for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1)) {
		RG_ASSERT(size[d] >= globalPt[d] * localPt[d], "Wrong grid partitioning");
	}
	this->nc = nc;
	RGCut<I>::setCutParams(size, globalPt);
	MPI_CHECK(MPI_Barrier(MPI_COMM_WORLD)); // wait before destroy communicator
	if (cartComm != MPI_COMM_NULL) {
		rgmpi::commFree(cartComm);
	}
	rgmpi::cartCreate(cartComm, globalPt);
	if (cartComm == MPI_COMM_NULL) return; // make stubs if we choose to take more processes later		
	cartRank = rgmpi::commRank(cartComm);
	rgmpi::cartCoords(cartComm, cartRank, cartPos);

	// set neighbours ranks
	for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1)) {
		if (cartPos[d] == 0) {
			neigh[SIDE_LEFT][d] = NO_NEIGHBOUR;
		} else {
			--cartPos[d];			// use slow implementation, make a copy of local dac in single DArray
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
	dacCut = RGCut<I>(
	             RGCut<I>::partNodes(X, cartPos[X]),
	             RGCut<I>::partNodes(Y, cartPos[Y]),
	             RGCut<I>::partNodes(Z, cartPos[Z]),
	             localPt[X],
	             localPt[Y],
	             localPt[Z]
	         );
	generateSubarrayTypes();
	// allocate memory for all DArray's
	Dim3D<I> globalSize(size[X], size[Y], size[Z]);
	Dim3D<I> localSize(RGCut<I>::partNodes(X, cartPos[X]), RGCut<I>::partNodes(Y, cartPos[Y]), RGCut<I>::partNodes(Z, cartPos[Z]));
	Dim3D<I> localParts(localPt[X], localPt[Y], localPt[Z]);
	Dim3D<I> dimGhost(ghost[X], ghost[Y], ghost[Z]);
	Dim3D<I> origin(RGCut<I>::partOrigin(X, cartPos[X]), RGCut<I>::partOrigin(Y, cartPos[Y]), RGCut<I>::partOrigin(Z, cartPos[Z]));
	
	dac.setParts(
		globalSize,
		localSize,
		localParts,
		dimGhost,
		origin,
		nc
	); 
}

template <typename T, typename I>
void DArrayScatter<T, I>::setAndScatter(int const processRank, DArray<T, I> const &da) {
	/// \todo for every da create new subarray types
	if (cartRank == processRank) {
		RG_ASSERT(da.ghost(X) == ghost[X] && da.ghost(Z) == ghost[Z] && da.ghost(Z) == ghost[Z], "ghost nodes have to be equal");
	}
	if (cartComm == MPI_COMM_NULL) return;
	std::vector<MPI_Request> req;
	// send info about sizes to other processes
	struct newSizes {
		I size[ALL_DIRS];
		I nc;
	} ns;
	if (processRank == cartRank) {
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d+1)) {
			ns.size[d] = da.localSize(d);
		}
		ns.nc = da.getNC();
	}
	// send header info to other processes
	const int nitems = 2;
	const int blocklengths[nitems] = {ALL_DIRS, 1};
	const MPI_Datatype types[nitems] = {rgmpi::getMPItype<I>(), rgmpi::getMPItype<I>()};
	MPI_Datatype mpi_ns_type;
	MPI_Aint offsets[nitems] = { offsetof(newSizes, size),
	                             offsetof(newSizes, nc)};
	MPI_CHECK(MPI_Type_create_struct(nitems, blocklengths, offsets, types, &mpi_ns_type));
	MPI_CHECK(MPI_Type_commit(&mpi_ns_type));
	MPI_CHECK(MPI_Bcast(&ns, 1, mpi_ns_type, processRank, cartComm));
	MPI_CHECK(MPI_Type_free(&mpi_ns_type));
	// resize local DArrayContainers and create types
	I globalPt[ALL_DIRS] = { RGCut<I>::numParts(X),
	                         RGCut<I>::numParts(Y),
	                         RGCut<I>::numParts(Z) };
	I localPt[ALL_DIRS] = { dacCut.numParts(X),
	                        dacCut.numParts(Y),
	                        dacCut.numParts(Z) };
	setSizes(ns.size, globalPt, localPt, ghost, ns.nc);
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
	    ghost[X], ghost[Y], ghost[Z],
	    nc);
	MPI_CHECK(MPI_Recv(&tda[0], 1, arrayDt[cartRank], processRank, SCATTER, cartComm, MPI_STATUS_IGNORE));
	dac.setDArray(tda, localPt[X], localPt[Y], localPt[Z]);
	if (processRank == cartRank) {
		// wait while all parts will be received
		MPI_CHECK(MPI_Waitall(RGCut<I>::numParts(), &req[0], MPI_STATUSES_IGNORE));
	}
}

template <typename T, typename I>
void DArrayScatter<T, I>::generateSubarrayTypes() {
	if (cartComm == MPI_COMM_NULL) return;
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
						nc
					};
					int array_of_subsizes[4] = {
						RGCut<I>::partNodes(X, i),
						RGCut<I>::partNodes(Y, j),
						RGCut<I>::partNodes(Z, k),
						nc
					};
					int array_of_starts[4] = {
						ghost[X] + RGCut<I>::partOrigin(X, i),
						ghost[Y] + RGCut<I>::partOrigin(Y, j),
						ghost[Z] + RGCut<I>::partOrigin(Z, k),
						0
					};
					subArrayDt[rank] = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
				}
				{
					// generate arrayDt
					int array_of_sizes[4] = {
						RGCut<I>::partNodes(X, i) + 2 * ghost[X],
						RGCut<I>::partNodes(Y, j) + 2 * ghost[Y],
						RGCut<I>::partNodes(Z, k) + 2 * ghost[Z],
						nc
					};
					int array_of_subsizes[4] = {
						RGCut<I>::partNodes(X, i),
						RGCut<I>::partNodes(Y, j),
						RGCut<I>::partNodes(Z, k),
						nc
					};
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
	locArrayDt.resize(dacCut.numParts());
	for (I k = 0; k != dacCut.numParts(Z); ++k)
		for (I j = 0; j != dacCut.numParts(Y); ++j)
			for (I i = 0; i != dacCut.numParts(X); ++i) {
				// type for local DArrays
				int array_of_sizes[4] = {
					dacCut.partNodes(X, i) + 2 * ghost[X],
					dacCut.partNodes(Y, j) + 2 * ghost[Y],
					dacCut.partNodes(Z, k) + 2 * ghost[Z],
					nc
				};
				int array_of_subsizes[4] = {
					dacCut.partNodes(X, i),
					dacCut.partNodes(Y, j),
					dacCut.partNodes(Z, k),
					nc
				};
				int array_of_starts[4] = {
					ghost[X],
					ghost[Y],
					ghost[Z],
					0
				};
				locArrayDt.at(dacCut.linInd(i, j, k)) = rgmpi::createSubarrayType<T>(array_of_sizes, array_of_subsizes, array_of_starts);
			}
}

template <typename T, typename I>
void DArrayScatter<T, I>::releaseSubarrayTypes(std::vector<MPI_Datatype> &dt) {
	for (std::vector<MPI_Datatype>::iterator it = dt.begin(); it != dt.end(); ++it) {
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
	rgmpi::commFree(cartComm);
}

template<typename T, typename I>
void DArrayScatter<T, I>::gatherAndGet(int processRank, DArray<T, I> &da) {
	if (cartComm == MPI_COMM_NULL) return;
	std::vector<MPI_Request> req;
	if (processRank == cartRank) {
		da.resize(
		    RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
		    RGCut<I>::numNodes(X), RGCut<I>::numNodes(Y), RGCut<I>::numNodes(Z),
		    0, 0, 0,
		    ghost[X], ghost[Y], ghost[Z],
		    nc);
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
	if (cartComm == MPI_COMM_NULL) return;
	dac.sync();
}

template<typename T, typename I>
void DArrayScatter<T, I>::fillGhost() {
	if (cartComm == MPI_COMM_NULL) return;
	for (CartSide s = SIDE_LEFT; s != SIDE_ALL; s = static_cast<CartSide>(s + 1))
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1))
			if (neigh[s][d] == NO_NEIGHBOUR) {
				CartDir ort1, ort2;
				ortDirs(d, ort1, ort2);
				I k = (s == SIDE_LEFT) ? 0 : dac.numParts(d) - 1;
				for (I i = 0; i < dac.numParts(ort1); i++)
					for (I j = 0; j < dac.numParts(ort2); j++) {
						I coord[ALL_DIRS];
						coord[d] = k;
						coord[ort1] = i;
						coord[ort2] = j;
						I x = coord[X], y = coord[Y], z = coord[Z];
						DArray<T, I>& da = dac.getDArrayPart(x, y, z);
						da.fillGhost(d, s);
					}
			}
}

template<typename T, typename I>
void DArrayScatter<T, I>::externalSyncStart() {
	if (cartComm == MPI_COMM_NULL) return;
	for (CartSide s = SIDE_LEFT; s != SIDE_ALL; s = static_cast<CartSide>(s + 1))
		for (CartDir d = X; d != ALL_DIRS; d = static_cast<CartDir>(d + 1))
			if (neigh[s][d] != NO_NEIGHBOUR) {
				CartDir ort1, ort2;
				ortDirs(d, ort1, ort2);
				if (dac.numParts(ort1) == 1 && dac.numParts(ort2) == 1) {
					// optimized synchronization without temporary buffer when there are no local partitioning
					
					// current DArray
					Dim3D<I> locIdx(0, 0, 0);
					locIdx[d] = (s == SIDE_LEFT) ? 0 : dac.numParts(d) - 1;
					DArray<T,I>& cda = dac.getDArrayPart(locIdx.x, locIdx.y, locIdx.z);
					// TODO don't generate types in every sync
					int sizes[4] = {
						cda.localGhostSize(X),
						cda.localGhostSize(Y),
						cda.localGhostSize(Z),
						nc
					};
					int subsizes[4] = { 0, 0, 0, nc };
					subsizes[d] = ghost[d];
					subsizes[ort1] = cda.localSize(ort1);
					subsizes[ort2] = cda.localSize(ort2);
					int starts[4] = { ghost[X], ghost[Y], ghost[Z], 0 };
					// types for ghost nodes in current DArray (recv to this nodes)
					starts[d] = (s == SIDE_LEFT) ? 0 : cda.localGhostSize(d) - ghost[d];
					MPI_Datatype ghostDT = rgmpi::createSubarrayType<T>(sizes, subsizes, starts);
					// types for boundary (not ghost) nodes in current DArray (send from this nodes)
					starts[d] = (s == SIDE_LEFT) ? ghost[d] : cda.localGhostSize(d) - 2 * ghost[d];
					MPI_Datatype boundDT = rgmpi::createSubarrayType<T>(sizes, subsizes, starts);
					
					MPI_CHECK(MPI_Irecv(cda.getDataRaw(), 1, ghostDT, neigh[s][d], SYNC, cartComm, &syncRecvReq[s + d * SIDE_ALL]));
					MPI_CHECK(MPI_Isend(cda.getDataRaw(), 1, boundDT, neigh[s][d], SYNC, cartComm, &syncSendReq[s + d * SIDE_ALL]));
				} else {
					I bufSize = dac.numNodes(ort1) * dac.numNodes(ort2) * ghost[d] * nc;
					// start recv
					recvBuf[s][d].resize(bufSize);
					MPI_CHECK(MPI_Irecv(&recvBuf[s][d].front(), bufSize, rgmpi::getMPItype<T>(), neigh[s][d], SYNC, cartComm, &syncRecvReq[s + d * SIDE_ALL]));
					// start send
					
					// now we have many DA inside DAC and we want to get ghost plane from DAC,
					// so we have to combine ghost planes from every DA placed on the border of local DAC
					
					// sub array size and position
					I orig[ALL_DIRS] = { 0, 0, 0 };
					I size[ALL_DIRS] = {
						dac.numNodes(X),
						dac.numNodes(Y),
						dac.numNodes(Z),
					};
					size[d] = ghost[d];
					orig[d] = (s == SIDE_LEFT) ? 0 : dac.numNodes(d) - ghost[d];
					dac.getSubArray(orig[X], orig[Y], orig[Z], size[X], size[Y], size[Z], sendBuf[s][d]);
					MPI_CHECK(MPI_Isend(&sendBuf[s][d].front(), bufSize, rgmpi::getMPItype<T>(), neigh[s][d], SYNC, cartComm, &syncSendReq[s + d * SIDE_ALL]));
				}
			} else {
				syncRecvReq[s + d * SIDE_ALL] = MPI_REQUEST_NULL;
				syncSendReq[s + d * SIDE_ALL] = MPI_REQUEST_NULL;
			}
}

template<typename T, typename I>
void DArrayScatter<T, I>::externalSyncEnd() {
	if (cartComm == MPI_COMM_NULL) return;
	while (1) {
		int index;
		MPI_CHECK(MPI_Waitany(SIDE_ALL * ALL_DIRS, syncRecvReq, &index, MPI_STATUS_IGNORE));
		if (index == MPI_UNDEFINED) {
			break;
		}
		CartSide s = static_cast<CartSide>(index % SIDE_ALL);
		CartDir d = static_cast<CartDir>(index / SIDE_ALL);
		CartDir ort1, ort2;
		ortDirs(d, ort1, ort2);
		
		if (dac.numParts(ort1) == 1 && dac.numParts(ort2) == 1) {
			// nothing to do, optimized synchronization
		} else {
			// sub array size and position
			I orig[ALL_DIRS] = { 0, 0, 0 };
			I size[ALL_DIRS] = {
				dac.numNodes(X),
				dac.numNodes(Y),
				dac.numNodes(Z),
			};
			size[d] = ghost[d];
			orig[d] = (s == SIDE_LEFT) ? -ghost[d] : dac.numNodes(d);
			dac.setSubArrayWithGhost(orig[X], orig[Y], orig[Z], size[X], size[Y], size[Z], recvBuf[s][d]);
		}
	}
	MPI_CHECK(MPI_Waitall(SIDE_ALL * ALL_DIRS, syncSendReq, MPI_STATUSES_IGNORE));
}

} /* namespace rgrid */

#endif // USE_MPI

#endif /* DARRAY_SCATTER_H */
