#ifndef GRID_OPERATOR_H
#define GRID_OPERATOR_H

#include <cmath>
#include "rgrid/darrayscatter.h"
#include "rgrid/nodeoperator.h"

namespace rgrid {

namespace operators {

template <typename T, typename I>
struct GridOp {
private:
	typedef DArrayScatter<T, I> DAS;
	typedef DArray<T, I> DA;
	
public:
	typedef VarOp<T, I> VarType;
	typedef ConstOp<T, I> ConstType;
	
	template <typename Expr>
	void apply(const Range<I>& r, Expr expr, DAS& out)
	{
		for (I gk = 0; gk != out.numParts(Z); ++gk)
		for (I gj = 0; gj != out.numParts(Y); ++gj)
		for (I gi = 0; gi != out.numParts(X); ++gi) {
			DA& nu = out.getDArrayPart(gi, gj, gk);
			
			for (VarsIt it = m_var_op.begin(); it != m_var_op.end(); ++it) {
				expr.attachInput(it->first, it->second->getDArrayPart(gi,gj,gk));
			}
			
			Range<I> glob = intersect(r, nu.getRange());
			Range<I> loc = glob.shift(-nu.origin());
			
			for (I k = loc.get(Z, SIDE_LEFT); k != loc.get(Z, SIDE_RIGHT); ++k)
			for (I j = loc.get(Y, SIDE_LEFT); j != loc.get(Y, SIDE_RIGHT); ++j)
			for (I i = loc.get(X, SIDE_LEFT); i != loc.get(X, SIDE_RIGHT); ++i) {
				nu(i, j, k, 0) = expr.eval(i, j, k);
			}
		}
	}
	
	VarType setVar(const std::string& name, DAS& in) {
		m_var_op.erase(name);
		m_var_op.insert(std::make_pair(name, &in));
		return VarType(name);
	}
	
private:
	typedef std::map<std::string, DAS*> Vars;
	typedef std::map<std::string, ConstOp<T, I> > Consts;
	
	typedef typename Vars::iterator VarsIt;
	typedef typename Consts::iterator ConstsIt;
	
	Vars m_var_op;
};

/*
 * if Side == SIDE_LEFT, rightmost value has pml parameter 0 and leftmost pml_max[D]
 * if SIde == SIDE_RIGHT, rightmost value has pml parameter pml_max[D] and leftmost 0
 * With PMLOffset you can specify shift of nodes relative unshifted grid
 * where PML is specified
 */
// template <
// 	typename T, typename I,
// 	template <typename> class Func,
// 	CartSide Side, CartDir Dir, PMLType PT, StaggeredOffset NodeOffset,
// 	template <typename,typename,CartDir> class NodeOperator>
// struct GridPMLOperator {
// 	typedef DArrayScatter<T, I> DAS;
// 	typedef DArray<T, I> DA;
// 	typedef NodeOperator<T, I, Dir> NO;
// 	
// 	GridPMLOperator(
// 		const TaskParams<T, I>& tp,
// 		const Range<I>& range,
// 		const DAS& in, DAS& out)
// 	: m_tp(tp), m_in(in), m_out(out), m_range(range)
// 	{
// 		m_pml_len = range.get(Dir, SIDE_RIGHT) - range.get(Dir, SIDE_LEFT);
// 		
// 		// create DArrays with pml values for every local DArray
// 		m_pml.resize(in.numParts(X) * in.numParts(Y) * in.numParts(Z));
// 
// 		// create arrays with pml coefficients
// 		m_pml_b.resize(m_pml_len);
// 		m_pml_a.resize(m_pml_len);
// 		
// 		// calculate offset
// 		T node_offset = 0;
// 		if (NodeOffset == OFFSET_NONE) {
// 			node_offset = 0;
// 		} else if (NodeOffset == OFFSET_PLUS_HALF) {
// 			node_offset = 0.5;
// 		} else if (NodeOffset == OFFSET_MINUS_HALF) {
// 			node_offset = -0.5;
// 		}
// 		
// 		T pml_offset = Side == SIDE_RIGHT ? node_offset : -node_offset; 
// 		
// 		// calculate all coefficients
// 		for (I i = 0; i != m_pml_len; ++i) {
// 			T pml_depth = i + pml_offset;
// 			if (pml_depth < 0) pml_depth = 0;
// 
// 			T pml_d = tp.pml_d_max[Dir] * pow(pml_depth * 1.0 / m_pml_len, 2);
// 			
// 			m_pml_b.at(i) = exp(-pml_d * tp.time_step);
// 			m_pml_a.at(i) = m_pml_b.at(i) - 1;
// 		}
// 		
// 		// for every local DA in DAC
// 		for (I gk = 0; gk != in.numParts(Z); ++gk)
// 		for (I gj = 0; gj != in.numParts(Y); ++gj)
// 		for (I gi = 0; gi != in.numParts(X); ++gi) {
// 			// get its local DA
// 			const DA& u = in.getDArrayPart(gi, gj, gk);
// 			//DA& nu = out.getDArrayPart(gi, gj, gk);
// 			
// 			// find intersection in global indexes
// 			Range<I> glob = intersect(range, u.getRange());
// 			// shift intersection origin to origin of current DA
// 			// so it will be local DA intersection
// 			Range<I> loc = glob.shift(-u.origin());
// 			// get size of intersection
// 			Dim3D<I> loc_size = loc.shift(-loc.get(SIDE_LEFT)).get(SIDE_RIGHT);
// 			
// 			// index of pml values
// 			int ind = gi + gj * in.numParts(X) + gk * in.numParts(X) * in.numParts(Y);
// 			// set size of pml value as size of current DA intersection
// 			// now pml indices lay from 0 to loc_size
// 			m_pml.at(ind).resize(
// 				loc_size, loc_size, Dim3D<I>(0,0,0), Dim3D<I>(0,0,0), 1);
// 		}
// 	}
// 	void apply() {
// 		STATIC_ASSERT(Dir == X || Dir == Y  || Dir == Z, GridPMLOperator_Wrong_CartSide_param);
// 		STATIC_ASSERT(Side == SIDE_LEFT || Side == SIDE_RIGHT, GridPMLOperator_Wrong_CartSide_param);
// 		
// 		for (I gk = 0; gk != m_in.numParts(Z); ++gk)
// 		for (I gj = 0; gj != m_in.numParts(Y); ++gj)
// 		for (I gi = 0; gi != m_in.numParts(X); ++gi) {
// 			const DA& u = m_in.getDArrayPart(gi, gj, gk);
// 			DA& nu = m_out.getDArrayPart(gi, gj, gk);
// 			
// 			Range<I> glob = intersect(m_range, u.getRange());
// 			Dim3D<I> pml_offset = glob.get(SIDE_LEFT) - m_range.get(SIDE_LEFT); // index of current pml coefficient
// 			Range<I> loc = glob.shift(-u.origin());
// 			//Dim3D<I> loc_size = loc.shift(-loc.get(SIDE_LEFT)).get(SIDE_RIGHT);
// 			
// 			int ind = gi + gj * m_in.numParts(X) + gk * m_in.numParts(X) * m_in.numParts(Y);
// 			
// 			for (I k = loc.get(Z, SIDE_LEFT); k != loc.get(Z, SIDE_RIGHT); ++k)
// 			for (I j = loc.get(Y, SIDE_LEFT); j != loc.get(Y, SIDE_RIGHT); ++j)
// 			for (I i = loc.get(X, SIDE_LEFT); i != loc.get(X, SIDE_RIGHT); ++i) {
// 				
// 				T temp_result = NO::apply(m_tp, i, j, k, u);
// 				
// 				I k_pml = k - loc.get(Z, SIDE_LEFT);
// 				I j_pml = j - loc.get(Y, SIDE_LEFT);
// 				I i_pml = i - loc.get(X, SIDE_LEFT);
// 				
// 				I ind_pml = pml_offset[Dir];
// 				if (Dir == X) {
// 					ind_pml += i_pml;
// 				} else if (Dir == Y) {
// 					ind_pml += j_pml;
// 				} else if (Dir == Z) {
// 					ind_pml += k_pml;
// 				}
// 				
// 				if (Side == SIDE_LEFT) {
// 					ind_pml = m_pml_len - 1 - ind_pml;
// 					// lower index -> less pml_d
// 				} else if (Side == SIDE_RIGHT) {
// 					// lower index -> more pml_d
// 					// do nothing
// 				}
// 				
// 				T pml_b = m_pml_b.at(ind_pml);
// 				T pml_a = m_pml_a.at(ind_pml);
// 				
// 				m_pml.at(ind)(i_pml, j_pml, k_pml, 0) =
// 					pml_b * m_pml.at(ind)(i_pml, j_pml, k_pml, 0) +
// 					pml_a * temp_result;
// 				temp_result = temp_result + m_pml.at(ind)(i_pml, j_pml, k_pml, 0);
// 				
// 				Func<T>::f(nu(i, j, k, 0), temp_result);
// 			}
// 			
// 		}
// 	}
// private:
// 	std::vector<T> m_pml_a;
// 	std::vector<T> m_pml_b;
// 	
// 	std::vector<DA> m_pml;
// 	
// 	I m_pml_len;
// 	const TaskParams<T, I>& m_tp;
// 	const DAS& m_in;
// 	DAS& m_out;
// 	Range<I> m_range;
// };

} // namespace operators

} // namespace rgrid

#endif


