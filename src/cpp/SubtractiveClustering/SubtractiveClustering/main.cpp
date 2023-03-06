#include "SC.h"
#define console(o) std::cout << o << std::endl

int main(int argc, char** argv[])
{
	double R1 = 5.0; double R2 = R1 * 1.25; double eps = 0.001;
	SC sc(R1, R2, eps);
	std::string file_name = "data.txt";
	Matrix<double>* data = new Matrix<double>;
	*data = Matrix<double>::read_from_file(file_name);
	sc.fit(data);
	sc.centers->print();
	std::cin.get();

}