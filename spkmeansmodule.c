#define PY_SSIZE_T_CLEAN
#include <Python.h>
#include "helper.h"
#include "nsclustering.h"
/*
we get PyObject that represents 2 dimensional array and we convert it
to 2 dimensional double array (C object)
 */

void convertPython2DArray(PyObject *python_arr,double **arr,int rows,int columns){
    int i,j;
    PyObject *current_list,*item;
    for (i = 0; i < rows; i++) {
        current_list = PyList_GetItem(python_arr, i);
        for(j=0;j<columns;j++){
            item = PyList_GetItem(current_list, j);
            if (!PyFloat_Check(item)){
                arr[i][j] = 0.0;
            }
            arr[i][j] = PyFloat_AsDouble(item);
        }
    }
}

/*
we get 2 dimensional array and num of rows and columns,
we return 2 dimensional PyObject
*/

PyObject* createPyObjectFrom2DArray(double **arr,int rows,int columns, int k){
    int i,j;
    PyObject *current,*result,*val;

    result = PyList_New(k);
    for(i=0;i<rows;i++){
        current = PyList_New(columns);
        for(j=0;j<columns;j++){
            val = PyFloat_FromDouble(arr[i][j]);
            PyList_SetItem(current, j, val);
        }
        PyList_SetItem(result, i, current);
    }
    return result;
}

/*
 * This actually defines the geo function using a wrapper C API function
 * The wrapping function needs a PyObject* self argument.
 * This is a requirement for all functions and methods in the C API.
 * It has input PyObject *args from Python.
 */

/*
We perceive the variables from the Python program, initialize the
*/

static PyObject* kmeans(PyObject *self, PyObject *args){

    PyObject *py_centroids,*py_data_points,*result;
    int k, max_iter, num_rows, d;
    double **data_points, **centroids;
    double epsilon;

    /* This parses the Python arguments into a double (d)  variable named z and int (i) variable named n*/
    if(!PyArg_ParseTuple(args, "iiiidOO", &k, &max_iter,&d,&num_rows,&epsilon,&py_centroids,&py_data_points)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    centroids = allocate_array_2d(k, d);
    data_points = allocate_array_2d(num_rows, d);

    convertPython2DArray(py_centroids,centroids,k,d);
    convertPython2DArray(py_data_points,data_points,num_rows,d);

    algorithm(k,d,num_rows,max_iter,data_points, centroids);

    result = createPyObjectFrom2DArray(centroids,k,d,k);

    free_array_2d(centroids, k);
    free_array_2d(data_points,num_rows);

    return result;
}

static PyObject* spk(PyObject *self, PyObject *args){

    PyObject *py_data_points,*result;
    double **data_points;
    int num_rows,k,d;
    Eigen * eigen;

    if(!PyArg_ParseTuple(args, "iiiO", &k,&d,&num_rows,&py_data_points)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }
    data_points = allocate_array_2d(num_rows, d);
    convertPython2DArray(py_data_points,data_points,num_rows,d);

    double** weighted_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_weighted_matrix(weighted_matrix,data_points,num_rows,d);
    double** diagonal_degree_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_diagonal_degree_matrix(diagonal_degree_matrix,weighted_matrix,num_rows);

    double** lnorm_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_lnorm_matrix(lnorm_matrix,weighted_matrix,diagonal_degree_matrix,num_rows);

    eigen = calloc(num_rows, sizeof (Eigen));
    error_occurred(eigen == NULL);

    init_Eigen_struct(eigen,num_rows);
    Jacobi_algorithm(lnorm_matrix,num_rows,eigen);

    if (k == 0){
        k = calculate_eigengap_heuristic(eigen,num_rows);
    }
    error_occurred(k == 1);

    double** U_matrix = allocate_array_2d(num_rows, k);
    set_U_matrix(U_matrix,eigen,num_rows,k);

    double** T_matrix = allocate_array_2d(num_rows, k);
    calculate_T_matrix(T_matrix,U_matrix,num_rows,k);

    result = createPyObjectFrom2DArray(T_matrix,num_rows,k,num_rows);

    free_array_2d(lnorm_matrix, num_rows);
    free_array_2d(diagonal_degree_matrix,num_rows);
    free_array_2d(weighted_matrix, num_rows);
    free_Eigen_struct(eigen,num_rows);
    free_array_2d(U_matrix,num_rows);
    free_array_2d(T_matrix,num_rows);

    return result;
}

static PyObject* wam(PyObject *self, PyObject *args) {

    PyObject *py_data_points,*result;
    double **data_points;
    int num_rows,d;

    if(!PyArg_ParseTuple(args, "iiO", &d,&num_rows,&py_data_points)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    data_points = allocate_array_2d(num_rows, d);
    convertPython2DArray(py_data_points,data_points,num_rows,d);

    double** weighted_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_weighted_matrix(weighted_matrix,data_points,num_rows,d);

    result = createPyObjectFrom2DArray(weighted_matrix,num_rows,num_rows,num_rows);

    free_array_2d(weighted_matrix,num_rows);
    free_array_2d(data_points,num_rows);
    return result;
}

static PyObject* ddg(PyObject *self, PyObject *args) {

    PyObject *py_data_points,*result;
    double **data_points;
    int num_rows,d;

    if(!PyArg_ParseTuple(args, "iiO" ,&d,&num_rows,&py_data_points)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    data_points = allocate_array_2d(num_rows, d);
    convertPython2DArray(py_data_points,data_points,num_rows,d);

    double** weighted_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_weighted_matrix(weighted_matrix,data_points,num_rows,d);
    double** diagonal_degree_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_diagonal_degree_matrix(diagonal_degree_matrix,weighted_matrix,num_rows);

    result = createPyObjectFrom2DArray(diagonal_degree_matrix,num_rows,num_rows,num_rows);

    free_array_2d(weighted_matrix,num_rows);
    free_array_2d(diagonal_degree_matrix,num_rows);
    free_array_2d(data_points,num_rows);
    return result;

}

static PyObject* lnorm(PyObject *self, PyObject *args) {

    PyObject *py_data_points,*result;
    double **data_points;
    int num_rows,d;

    if(!PyArg_ParseTuple(args, "iiO" ,&d,&num_rows,&py_data_points)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    data_points = allocate_array_2d(num_rows, d);
    convertPython2DArray(py_data_points,data_points,num_rows,d);

    double** weighted_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_weighted_matrix(weighted_matrix,data_points,num_rows,d);
    double** diagonal_degree_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_diagonal_degree_matrix(diagonal_degree_matrix,weighted_matrix,num_rows);
    double** lnorm_matrix = allocate_array_2d(num_rows, num_rows);
    calculate_lnorm_matrix(lnorm_matrix,weighted_matrix,diagonal_degree_matrix,num_rows);

    result = createPyObjectFrom2DArray(lnorm_matrix,num_rows,num_rows,num_rows);

    free_array_2d(weighted_matrix,num_rows);
    free_array_2d(diagonal_degree_matrix,num_rows);
    free_array_2d(lnorm_matrix,num_rows);
    free_array_2d(data_points,num_rows);
    return result;

}

static PyObject* jacobi(PyObject *self, PyObject *args) {

    PyObject *py_symmetric_matrix,*result;
    double **symmetric_matrix, **eigen_matrix;
    int n,j;
    Eigen *eigen;

    if(!PyArg_ParseTuple(args, "iO" , &n, &py_symmetric_matrix)) {
        return NULL; /* In the CPython API, a NULL value is never valid for a
                        PyObject* so it is used to signal that an error has occurred. */
    }

    symmetric_matrix = allocate_array_2d(n, n);
    convertPython2DArray(py_symmetric_matrix, symmetric_matrix, n, n);

    eigen = calloc(n, sizeof (Eigen));
    error_occurred(eigen == NULL);

    init_Eigen_struct(eigen,n);
    Jacobi_algorithm(symmetric_matrix,n,eigen);

    eigen_matrix = allocate_array_2d(n+1, n); /*the first nXn is the vectors and the last row is the values*/
    set_U_matrix(eigen_matrix,eigen,n,n);

    for (j=0;j<n;j++){
        eigen_matrix[n][j] = eigen[j].value;
    }

    result = createPyObjectFrom2DArray(eigen_matrix,n+1,n,n+1);


    free_Eigen_struct(eigen,n);
    free_array_2d(symmetric_matrix,n);
    free_array_2d(eigen_matrix,n);
    return result;
}

static PyMethodDef capiMethods[] = {
        {"kmeans",(PyCFunction) kmeans,METH_VARARGS,PyDoc_STR("Kmeans algorithm to find the clusters")},
        {"spk",(PyCFunction) spk,METH_VARARGS,PyDoc_STR("Perform full spectral kmeans")},
        {"wam",(PyCFunction) wam,METH_VARARGS,PyDoc_STR("Calculate and output the Weighted Adjacency Matrix")},
        {"ddg",(PyCFunction) ddg,METH_VARARGS,PyDoc_STR("Calculate and output the Diagonal Degree Matrix ")},
        {"lnorm",(PyCFunction) lnorm,METH_VARARGS,PyDoc_STR("Calculate and output the Normalized Graph Laplacian")},
        {"jacobi",(PyCFunction) jacobi,METH_VARARGS,PyDoc_STR("Calculate and output the eigenvalues and eigenvectors")},
        {NULL, NULL, 0, NULL}
};


/* This initiates the module using the above definitions. */
static struct PyModuleDef moduledef = {
        PyModuleDef_HEAD_INIT,
        "spkmeans", /* name of module */
        NULL, /* module documentation, may be NULL */
        -1,  /* size of per-interpreter state of the module, or -1 if the module keeps state in global variables. */
        capiMethods /* the PyMethodDef array from before containing the methods of the extension */
};


/*
 * The PyModuleDef structure, in turn, must be passed to the interpreter in the module’s initialization function.
 * The initialization function must be named PyInit_name(), where name is the name of the module and should match
 * what we wrote in struct PyModuleDef.
 * This should be the only non-static item defined in the module file
 */
PyMODINIT_FUNC
PyInit_spkmeans(void)
{
    PyObject *m;
    m = PyModule_Create(&moduledef);
    if (!m) {
        return NULL;
    }
    return m;
}