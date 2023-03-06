#include "SC.h"

SC::SC(): m_R1(5), m_R2(1.25 * 5), m_eps(0.003) { }

SC::SC(double R1, double R2, double eps): 
	m_R1(R1), 
	m_R2(R2), 
	m_eps(eps)
{

}

SC::~SC(){

	delete m_data;
}

double SC::get_R1() const
{
	/* Get equidistant lines R1 */
	return m_R1;
}

double SC::get_R2() const
{
	/* Get equidistant lines R2 */
	return m_R2;
}

double SC::get_eps() const
{
	/* Get stopping constant 'epsilon' */
	return m_eps;
}

void SC::set_R1(double R1)
{
	/* Set equidistant lines R1 */
	m_R1 = R1;
}

void SC::set_R2(double R2)
{
	/* Set equidistant lines R2 */
	m_R2 = R2;
}

void SC::set_eps(double eps)
{
	/* Set stopping constant 'epsilon' */
	m_eps = eps;
}

void SC::read_data() const
{
	m_data->print();
}

void SC::fit(Matrix<double>* data)
{
    Matrix<double>* center = new Matrix<double>;
    
	m_data = data;

    m_dist_matrix = DISTANCE_MATRIX(data, data, false);
    Matrix<double>* MF = new Matrix<double>;
    MF = compute_mountain_function(data);
    
    double max = maximum(MF);
    double max_star = 0;
    int it = 1;
    do 
    {
        int indmax = argmax(max, MF);
        center = find_center(indmax, data);
        if (it == 1)
            centers = center;
        else
            centers = vstack(centers, center);
        MF = update_mountain_function(MF, max, center, data);
        max_star = max;
        max = maximum(MF);
        it += 1;
    } while (max >= m_eps * max_star);

}

void SC::test()
{
    m_dist_matrix->print();
}


Matrix<double>* SC::ONES(int h, int w)
{
    /* Return a new array of given shape and type, filled with ones. */

    Matrix<double>* output = new Matrix<double>;
    double** array = new double* [h];
    for (auto i = 0; i < h; i++)
    {
        array[i] = new double[w];
        for (auto j = 0; j < w; j++)
            array[i][j] = 1;
    }
    output->set(h, w, array);
    return output;
}

Matrix<double>* SC::DISTANCE_MATRIX(Matrix<double>* A, Matrix<double>* B, bool squared)
{
    /* Compute the distance matrix. Returns the matrix of all pair-wise distances.*/

    Matrix<double>* AA = new Matrix<double>;
    Matrix<double>* O = new Matrix<double>;
    Matrix<double>* B_T = new Matrix<double>;
    *B_T = B->transpose();
    // double *sum = new double;
    int M = A->get_height();
    int N = B->get_height();
    double sum = 0.0;
    double** array = new double* [M];
    for (auto i = 0; i < M; i++)
    {
        array[i] = new double[0];
        for (auto j = 0; j < A->get_width(); j++)
        {
            sum += (double)A->get_data()[i][j] * (double)A->get_data()[i][j];
        }
        array[i][0] = sum;
        sum = 0.0;
    }
    AA->set(A->get_height(), 1, array);
    O = ONES(1, N);
    Matrix<double>* A_dots = new Matrix<double>;
    *A_dots = *AA * *O; // 3 x 2


    Matrix<double>* BB = new Matrix<double>;
    sum = 0.0;
    double** arrayB = new double* [N];
    for (auto i = 0; i < N; i++)
    {
        arrayB[i] = new double[B->get_width()];
        for (auto j = 0; j < B->get_width(); j++)
        {
            sum += (double)B->get_data()[i][j] * (double)B->get_data()[i][j];
        }
        arrayB[i][0] = sum;
        sum = 0.0;
    }
    BB->set(N, 1, arrayB);
    O = ONES(1, M);
    Matrix<double>* B_dots = new Matrix<double>;
    *B_dots = *BB * *O; 
    Matrix<double>* B_dots_T = new Matrix<double>;
    *B_dots_T = B_dots->transpose(); 
    Matrix<double>* D_squred = new Matrix<double>;
    *D_squred = *A_dots + *B_dots_T - 2 * *A * *B_T;

    if (!squared)
    {
        Matrix<double>* D = new Matrix<double>;
        auto h = D_squred->get_height(); auto w = D_squred->get_width();
        double** arrayD = new double* [h];
        for (auto i = 0; i < h; i++)
        {
            arrayD[i] = new double[w];
            for (auto j = 0; j < w; j++)
                arrayD[i][j] = std::sqrt((double)D_squred->get_data()[i][j]);
        }
        D->set(h, w, arrayD);
        return D;
    }
    return D_squred;
}


Matrix<double>* SC::compute_mountain_function(Matrix<double>* data)
{
    /* Step 1 of algoritm compute the potential (mountain function) for each point or sample.
        returns: the potential for that given row given no previous potential known  */

    Matrix<double>* P = new Matrix<double>;
    Matrix<double>* dist = new Matrix<double>;
    double denomerator = std::pow(m_R1 / 2, 2);
    dist = DISTANCE_MATRIX(data, data, false);
    auto h = data->get_height();
    auto w = data->get_width();
    double** array = new double* [h];
    double sum;
    for (auto i = 0; i < h; i++)
    {
        sum = 0;
        array[i] = new double[h];
        for (auto j = 0; j < h; j++)
        {
            sum += exp(-dist->get_data()[i][j] / denomerator);
        }
        array[i][0] = sum;
    }
    P->set(h, 1, array);
    return P;
}


double SC::maximum(Matrix<double>* data)
{
    /* Find maximum in array */
    double max = data->get_data()[0][0];
    for (auto i = 0; i < data->get_height(); i++)
    {
        for (auto j = 0; j < data->get_width(); j++)
        {
            if (data->get_data()[i][j] > max)
                max = data->get_data()[i][j];
        }
    }
    return max;
}


int SC::argmax(double maximum, Matrix<double>* data)
{
    /* Return index maximum's from data */
    auto h = data->get_height();
    auto w = data->get_width();
    int m = -1;

    for (auto i = 0; i < h; i++)
    {
        for (auto j = 0; j < w; j++)
        {
            if ((double)data->get_data()[i][j] == maximum)
            {
                m = i;
                return m;
            }
        }
    }
    if (m == -1)
        std::cout << "Such an element " << maximum << " does not exist in this array" << "\n";
    return m;
}


Matrix<double>*SC::find_center(int index, Matrix<double>* data)
{

    Matrix<double>* output = new Matrix<double>;
    double** array = new double* [0];
    array[0] = new double[data->get_width()];
    for (auto j = 0; j < data->get_width(); j++)
        array[0][j] = data->get_data()[index][j];
    output->set(1, data->get_width(), array);
    return output;
}


Matrix<double>* SC::vstack(Matrix<double>* total, Matrix<double>* center)
{
    /* Adding in stack new centroid of clusters */
    Matrix<double>* output = new Matrix<double>;
    int h = total->get_height();
    int w = total->get_width();
    double** array = new double* [h + 1];
    for (auto i = 0; i < h + 1; i++)
    {
        array[i] = new double[w];
        for (auto j = 0; j < w; j++)
        {
            if (i == h)
                array[i][j] = center->get_data()[0][j];
            else
                array[i][j] = total->get_data()[i][j];

        }
    }
    output->set(h + 1, w, array);
    return output;
}


Matrix<double>* SC::update_mountain_function(Matrix<double>* potential, double max, Matrix<double>* center, Matrix<double>* data)
{
    /* Step 3 of algorithm update the mountain function of each data. */
    
    Matrix<double>* output = new Matrix<double>;
    Matrix<double>* diff = new Matrix<double>;
    diff = subtr2Dfrom1D(data, center);
    Matrix<double>* numerator = new Matrix<double>;
    numerator = euclidean_distance_axis_1(diff);
    double denomerator = std::pow(m_R2 / 2, 2);
    Matrix<double>* ratio = new Matrix<double>;
    double** array = new double* [numerator->get_height()];
    for (auto i = 0; i < numerator->get_height(); i++)
    {
        array[i] = new double[numerator->get_width()];
        for (auto j = 0; j < numerator->get_width(); j++)
        {
            array[i][j] = std::exp(-numerator->get_data()[i][j] / denomerator);
        }
    }
    ratio->set(numerator->get_height(), numerator->get_width(), array);
    Matrix<double>* dot = new Matrix<double>;
    *dot = max * *ratio;
    *output = *potential - *dot;
    return output;
}


Matrix<double>* SC::subtr2Dfrom1D(Matrix<double>* Matrix2D, Matrix<double>* Matrix1D)
{
    /* Subtract 2D vector from 1D vector */

    Matrix<double>* output = new Matrix<double>;
    assert(Matrix2D->get_width() == Matrix1D->get_width());
    double** array = new double* [Matrix2D->get_height()];
    for (auto i = 0; i < Matrix2D->get_height(); i++)
    {
        array[i] = new double[Matrix2D->get_width()];
        for (auto j = 0; j < Matrix2D->get_width(); j++)
            array[i][j] = Matrix2D->get_data()[i][j] - Matrix1D->get_data()[0][j];
    }

    output->set(Matrix2D->get_height(), Matrix2D->get_width(), array);
    return output;
}


Matrix<double>* SC::euclidean_distance_axis_1(Matrix<double>* A)
{
    /* Compute Euclidean distance of array by axis 1 */

    Matrix<double>* output = new Matrix<double>;
    double** array = new double* [A->get_height()];
    double sum = 0.0;
    double distance;
    for (auto i = 0; i < A->get_height(); i++)
    {
        array[i] = new double[A->get_width()];
        for (auto j = 0; j < A->get_width(); j++)
        {
            sum = sum + std::pow(A->get_data()[i][j], 2.0);
            distance = std::sqrt(sum);
        }
        array[i][0] = distance; sum = 0.0; distance = 0.0;
    }
    output->set(A->get_height(), 1, array);
    return output;
}
