#ifndef FDC_H
#define FDC_H

namespace fdc {

static const int MAX_HALF_ORDER = 20;

extern const double C1S[MAX_HALF_ORDER][MAX_HALF_ORDER];
extern const double C1[MAX_HALF_ORDER][MAX_HALF_ORDER];
extern const double C2[MAX_HALF_ORDER][MAX_HALF_ORDER+1];

} // namespace fdc

#endif
