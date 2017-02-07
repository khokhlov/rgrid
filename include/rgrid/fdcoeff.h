#ifndef FD_COEFF_H
#define FD_COEFF_H

#include <cstdlib>
#include <cmath>

#include <vector>

#include "rgrid/debug.h"

namespace rgrid { 

/**
 * \brief compute n!
 */
inline double fact(const unsigned int n) {
	double res = 1;
	for (unsigned int i = 2; i <= n; ++i) res *= i;
	return res;
}

/**
 * \brief Compute end!/start!
 */
inline double partFact(const unsigned int start, const unsigned int end) {
	double res = 1;
	DEBUG_ASSERT(start+1 <= end, "start+1 > end !!!");
	for (unsigned int i = start+1; i <= end; ++i) res *= i;
	return res;
}

class FDCoeff {
public:
	/**
	 * \brief Calculate FD coefficients
	 * \param N half of order
	 */
	FDCoeff(int N) {
		// first derivative
		m_c1.resize(N+1, 0);
		m_c1.at(0) = 0;
		for (int n = 1; n <= N; ++n) {
			double mul = 1.0;
			for (int n_mul = 1; n_mul <= n; ++n_mul) {
				mul *= static_cast<double>(N - n + n_mul) / static_cast<double>(N + n_mul);
			}
			m_c1.at(n) = (n % 2 ? 1.0 : -1.0) * mul / n;
		}
		// second derivative
		m_c2.resize(N+1, 0);
		for (int n = 1; n <= N; ++n) {
			m_c2.at(n) = m_c1.at(n) * 2.0 / n;
			m_c2.at(0) += m_c2.at(n);
		}
		m_c2.at(0) *= -2;
		// first derivative staggered
		m_sc1.resize(N, 0);
		if (N == 1) m_sc1.at(0) = 1;
		else {
			for (int n = 0; n < N; ++n) {
				m_sc1.at(n) = pow(2,5) * pow(-1, n) * (2*N-1);
				m_sc1.at(n) /= pow(2*n+1, 2) * pow(2, 4*N);
				m_sc1.at(n) *= partFact(N-1, 2*N-1) * partFact(N-2, 2*N-3);
				m_sc1.at(n) /= fact(N-n-1) * fact(N+n);
			}
		}
	}
	/**
	 * Get specified FD coefficient
	 * \param n number of coefficient
	 */
	inline double c1(int n) {
		DEBUG_ASSERT(!m_c1.empty(), "Call calc() first");
		DEBUG_ASSERT(0 <= n && n < static_cast<int>(m_c1.size()), "Wrong number of coefficient");
		if (n >= 0) {
			return m_c1.at(n);
		} else {
			return - m_c1.at(-n);
		}
	}
	inline double c2(int n) {
		n = abs(n);
		DEBUG_ASSERT(!m_c2.empty(), "Call calc() first");
		DEBUG_ASSERT(0 <= n && n < static_cast<int>(m_c2.size()), "Call calc() first");
		return m_c2[abs(n)];
	}
	inline double sc1(const int n) const {
		DEBUG_ASSERT(!m_sc1.empty(), "Call calc() first");
		DEBUG_ASSERT(0 < n && n <= static_cast<int>(m_sc1.size()), "Call calc() first");
		return m_sc1[n-1];
	}
	
private:
	std::vector<double> m_c1; // first derivatice coefficients
	std::vector<double> m_c2; // second derivative coefficients
	std::vector<double> m_sc1; // first derivatice staggered coefficients
};

} // namespace rgrid

#endif

