#pragma once
/*
*	Implemetation Subtractive Clustering Method for computing the cluster centroids,
 *	which belongs to unsupervised learning and can quickly determine the number of clusters and
*	cluster centroids based on the raw data.
*	
*	Author: Alex Povod
*	Version: 0.1
*/

#include <iostream>
#include "matrix.h"
#include <cassert>

class SC
{
public:
	SC();
	SC(double R1, double R2, double eps);
	~SC();

	//readly-only function
	double get_R1(void) const; 
	double get_R2(void) const;
	double get_eps(void) const;
	void read_data() const;
	void test();


	// A write function; can't be const
	void set_R1(double R1);
	void set_R2(double R2);
	void set_eps(double eps);

	// Base function algorithm's
	void fit(Matrix<double>* data);
	Matrix<double>* centers = new Matrix<double>;
	

private:
	double m_R1;
	double m_R2;
	double m_eps;
	Matrix<double>* m_data = new Matrix<double>;

	Matrix<double>* ONES(int h, int w);
	Matrix<double>* DISTANCE_MATRIX(Matrix<double>* A, Matrix<double>* B, bool squared);
	Matrix<double>* m_dist_matrix = new Matrix<double>;

	Matrix<double>* compute_mountain_function(Matrix<double>* data);
	double maximum(Matrix<double>* data);
	int argmax(double maximum, Matrix<double>* data);
	Matrix<double>* find_center(int index, Matrix<double>* data);
	Matrix<double>* vstack(Matrix<double>* total, Matrix<double>* center);
	Matrix<double>* update_mountain_function(Matrix<double>* potential, double max, Matrix<double>* center, Matrix<double>* data);
	Matrix<double>* subtr2Dfrom1D(Matrix<double>* Matrix2D, Matrix<double>* Matrix1D);
	Matrix<double>* euclidean_distance_axis_1(Matrix<double>* A);

};

