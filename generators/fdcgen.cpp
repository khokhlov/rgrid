/**
 * Generator of finite difference coefficients
 */

#include <iostream>

#include "rgrid/fdcoeff.h"
#include "rgrid/fdc.h"

using namespace std;
using namespace fdcoeff;
using namespace fdc;

int main() {
	cout << "#include \"rgrid/fdc.h\"" << endl;
	cout << endl;
	cout << "namespace fdc {" << endl;
	cout << endl;

	cout << "const double SC1[MAX_HALF_ORDER][MAX_HALF_ORDER] = {" << endl;
	for (int i = 1; i != MAX_HALF_ORDER; ++i) {
		FDCoeff fdc(i);
		cout << "\t{";
		for (int j = 1; j < i; ++j) {
			cout << fdc.sc1(j) << ", ";
		}
		cout << fdc.sc1(i) << "}," << endl;
	}
	cout << "};" << endl;
	cout << endl;

	cout << "const double C1[MAX_HALF_ORDER][MAX_HALF_ORDER] = {" << endl;
	for (int i = 1; i != MAX_HALF_ORDER; ++i) {
		FDCoeff fdc(i);
		cout << "\t{";
		for (int j = 1; j < i; ++j) {
			cout << fdc.c1(j) << ", ";
		}
		cout << fdc.c1(i) << "}," << endl;
	}
	cout << "};" << endl;
	cout << endl;

	cout << "const double C2[MAX_HALF_ORDER][MAX_HALF_ORDER+1] = {" << endl;
	for (int i = 1; i != MAX_HALF_ORDER; ++i) {
		FDCoeff fdc(i);
		cout << "\t{";
		for (int j = 0; j < i; ++j) {
			cout << fdc.c2(j) << ", ";
		}
		cout << fdc.c2(i) << "}," << endl;
	}
	cout << "};" << endl;
	cout << endl;

	cout << "} // namespace fdc" << endl;
	return 0;
}
