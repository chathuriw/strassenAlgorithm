#include "Product.h"
#include <stdlib.h>
#include <stdio.h>
#include <math.h>
/*
  Product of two matrices
  + Product_basic(...) implements Algo 1.9
  + Product(...) makes use of Algo 1.9 and 5.3 to maximize speed

  @input
  A: pointer of int, the starting address of an Am * An matrix
  B: pointer of int, the starting address of an An * Bn matrix
  C: pointer of int, the starting address of an Am * Bn matrix

  @output
  void. the matrix C keeps the result of the product AB

  @other info
  - Each entry of A or B is an integer with value between -16 and +16

  @hints
  - For matrix A, (A + i * An  + j) will be the address of A[i][j]
*/

// cutoff value to call Product_basic
#define CUTOFF 64

void strassen_static_padding(int n, int *A, int as, int *B, int bs, int *C, int cs);
void strassen_dynamic_padding(int Am, int An, int Bn, int *A, int *B, int *C);
int find_closest_even(int Am, int An, int c);
void generate_square(int size, int *A, int Am, int An, int *S);
void copy_matrix(int *S, int *R, int rr, int rc, int stride);
void strassen_dynamic_peeling(int Am, int An, int Bn, int *A, int *B, int *C, int as, int bs, int cs);
void static_padding_init(int *A, int Am, int An, int *B, int Bn, int *C);

/**
 * Basic matrix multiplication algorithm
 */
void Product_basic(int *A, int Am, int An, int *B, int Bn,  int *C){
	int i, j, k;
	for (i=0; i < Am; i++){
		for (j=0 ; j < Bn ; j++){
			*(C + i * Bn +j) = 0;
			for (k=0; k < An; k++){
				*(C + i * Bn + j) += (*(A + i * An +k))*(*(B + k * Bn +j));
			}
		}
	}
}

/**
 * This method calls different versions of the Strassen's algorithm
 */
void Product(int *A, int Am, int An, int *B, int Bn, int *C) {
    // calling strassen_static_padding approach
    //static_padding_init(A, Am, An, B, Bn, C);
    // calling strassen_dynamic_padding approach
    //strassen_dynamic_padding(Am, An, Bn, A, B, C);
    // calling strassen's dynamic peeling approach
    strassen_dynamic_peeling(Am, An, Bn, A, B, C, An, Bn, Bn);
    // basic matrix multiplication
    //Product_basic(A, Am, An, B, Bn,C);
}

/**
 * This method calls Strassen's static padding method. Before calling static padding
 * it will find the closest even number according to the provided cut off and pad the
 * input matrices according to the closest even number.
 * At the end of calculation, it will copy the correct matrix to output matrix and
 * free the memory allocated in the process.
 */
void static_padding_init(int *A, int Am, int An, int *B, int Bn, int *C) {
	int sizeA = find_closest_even(Am, An, CUTOFF);
	int size = find_closest_even(sizeA, Bn, CUTOFF);

	int *SA = (int *) malloc((size * size) * sizeof(int));
	int *SB = (int *) malloc((size * size) * sizeof(int));
	int *SC = (int *) malloc((size * size) * sizeof(int));
	generate_square(size, A, Am, An, SA);
	generate_square(size, B, An, Bn, SB);
	generate_square(size, C, Am, Bn, SC);

	strassen_static_padding(size, SA, size, SB, size, SC, size);
	copy_matrix(SC, C, Am, Bn, size);
	free(SA);
	free(SB);
	free(SC);
}

/**
 * This is to copy source square matrix in to another matrix.
 *
 * */
void copy_matrix(int *S, int *R, int rr, int rc, int stride) {
    int i, k = 0;
    for (i = 0; i < rr; i++) {
        for (k = 0; k < rc; k++) {
            *(R + i * rc + k) = *(S + i * stride + k);
        }
    }
}

/**
 * This method is to find the closest even number for given size of the matrix. Closest
 * even number is calculated according to the equation c/2  < q ≤ c where c is the
 * cut off
 *
 * */
int find_closest_even(int Am, int An, int c) {
    int k;
    int max = Am > An ? Am : An;
    int p = 1;
    int q = c / 2 + 1;
    for (k = 0; ; k++) {
        for (q = c / 2 + 1;q <= c; q++) {
    	    if (q * p >= max) {
                return q*p;
    	    }
        }
        p = p * 2;
    }
    return 0; // will not reach
}

/**
 * This method is to padding 0s according to the size of the resultant 
 * matrix. A is the original matrix, S is the resultant matrix
 *
 *
 * */
void generate_square(int size, int *A, int Am, int An, int *S) {
    int i = 0;
    int k = 0;
    for (i = 0; i < size; i++) {
        for (k = 0; k < size; k++) {
            if (i >= Am) {
                *(S + i * size + k) = 0;
            } else if (k >= An) {
                *(S + i * size + k) = 0;
            } else {
                int temp = *(A + i * An + k);
                *(S + i * size + k) = temp;
            }

        }
    }
}

/**
 * This method is to add two matrices. Added results store in Matrix C. Correct
 * locations are calculated using the stride.
 * */
void add_matrix(int m, int n, int *A, int as, int *B, int bs, int *C, int cs) {
    int i = 0, j = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            *(C + i * cs + j) = *(A + i * as + j) + *(B + i * bs + j);
        }
    }
}

/**
 * This method is to subtract two matrices. Subtracted results store in Matrix C. Correct
 * locations are calculated using the stride.
 * */
void subtract_matrix(int m, int n, int *A, int as, int *B, int bs, int *C, int cs) {
    int i = 0, j = 0;
    for (i = 0; i < m; i++) {
        for (j = 0; j < n; j++) {
            *(C + i * cs + j) = *(A + i * as + j) - *(B + i * bs + j);
        }
    }
}


/**
 * This method is to return the index of the input matrix according to the
 * sub-matrix (11, 12, 21, 22).
 * */
int idx(int i, int j, int stride, int rows, int cols) {
    if (i == 1 && j == 1) {
        return 0;
    }
    if (i == 2 && j == 1) {
        return stride * rows;
    }
    if (i == 1 && j == 2) {
        return cols;
    }
    if (i == 2 && j == 2) {
        return stride * (rows) + cols;
    }
    return 0 ;
}

/**
 * Modified version of product basic where you can give the stride of the matrices
 */
void product_basic_with_stride (int Am, int An, int *A, int as, int *B, int Bn, int bs, int *C, int cs){
	int i, j, k, sum = 0;
	        for (i= 0; i < Am; i++){
	            for (j= 0 ; j < Bn ; j++){
	               sum = 0;
	               for (k= 0; k < An; k++){
	                   sum += (*(A + i * as +k))*(*(B + k * bs +j));
	               }
	               *(C + i * cs + j) = sum;
	            }
	        }
}

/**
 * This method is also a modified version of the product basic algorithm. Difference is
 * the product is added to existing values of the resultant matrix.
 */
void product_basic_with_stride_adding (int Am, int An, int *A, int as, int *B, int Bn, int bs, int *C, int cs){
	int i, j, k, sum = 0;
	        for (i= 0; i < Am; i++){
	            for (j= 0 ; j < Bn ; j++){
	               sum = 0;
	               for (k= 0; k < An; k++){
	                   sum += (*(A + i * as +k))*(*(B + k * bs +j));
	               }
	               *(C + i * cs + j) += sum;
	            }
	        }
}

/**
 * This is the strassen's algorithm with static padding
 *
 * */
void strassen_static_padding(int n, int *A, int as, int *B, int bs, int *C, int cs) {
	if (n <= CUTOFF) {
	        int i, j, k, sum = 0;
	        for (i= 0; i < n; i++){
	            for (j= 0 ; j < n ; j++){
	               sum = 0;
	               for (k= 0; k < n; k++){
	                   sum += (*(A + i * as +k))*(*(B + k * bs +j));
	               }
	               *(C + i * cs + j) = sum;
	            }
	        }
	        return;
	}
    int div = n / 2;
    int *S1 = (int *)malloc(sizeof(int) * div * div);
    add_matrix(div, div, A + idx(2, 1, as, div, div), as, A + idx(2, 2, as, div, div), as, S1, div);
    int *S2 = (int *)malloc(sizeof(int) * div * div);
    subtract_matrix(div, div, S1, div, A + idx(1, 1, as, div, div), as, S2, div);
    int *S3 = (int *)malloc(sizeof(int) * div * div);
    subtract_matrix(div, div, A + idx(1, 1, as, div, div), as, A + idx(2, 1, as, div, div), as, S3, div);
    int *S4 = (int *)malloc(sizeof(int) * div * div);
    subtract_matrix(div, div, A + idx(1, 2, as, div, div), as, S2, div, S4, div);
    int *S5 = (int *)malloc(sizeof(int) * div * div);
    subtract_matrix(div, div, B + idx(1, 2, bs, div, div), bs, B + idx(1, 1, bs, div, div), bs, S5, div);
    int *S6 = (int *)malloc(sizeof(int) * div * div);
    subtract_matrix(div, div, B + idx(2, 2, bs, div, div), bs, S5, div, S6, div);
    int *S7 = (int *)malloc(sizeof(int) * div * div);
    subtract_matrix(div, div, B + idx(2, 2, bs, div, div), bs, B + idx(1, 2, bs, div, div), bs, S7, div);
    int *S8 = (int *)malloc(sizeof(int) * div * div);
    subtract_matrix(div, div, S6, div, B + idx(2, 1, bs, div, div), bs, S8, div);

    int *M1 = (int *)malloc(sizeof(int) * div * div);
    strassen_static_padding(div, S2, div, S6, div, M1, div);
    int *M2 = (int *)malloc(sizeof(int) * div * div);
    strassen_static_padding(div, A + idx(1, 1, as, div, div), as, B + idx(1, 1, bs, div, div), bs, M2, div);
    int *M3 = (int *)malloc(sizeof(int) * div * div);
    strassen_static_padding(div, A + idx(1, 2, as, div, div), as, B + idx(2, 1, bs, div, div), bs, M3, div);
    int *M4 = (int *)malloc(sizeof(int) * div * div);
    strassen_static_padding(div, S3, div, S7, div, M4, div);
    int *M5 = (int *)malloc(sizeof(int) * div * div);
    strassen_static_padding(div, S1, div, S5, div, M5, div);
    int *M6 = (int *)malloc(sizeof(int) * div * div);
    strassen_static_padding(div, S4, div, B + idx(2, 2, bs, div, div), bs, M6, div);
    int *M7 = (int *)malloc(sizeof(int) * div * div);
    strassen_static_padding(div, A + idx(2, 2, as, div, div), as, S8, div, M7, div);
      
    int *T1 = (int *)malloc(sizeof(int) * div * div);
    add_matrix(div, div, M1, div, M2, div, T1, div);
    int *T2 = (int *)malloc(sizeof(int) * div * div);
    add_matrix(div, div, T1, div, M4, div, T2, div);
    
    add_matrix(div, div, M2, div, M3, div, C + idx(1, 1, cs, div, div), cs);
    add_matrix(div, div, T1, div, M5, div, C + idx(1, 2, cs, div, div), cs);
    add_matrix(div, div, C + idx(1, 2, cs, div, div), cs, M6, div, C + idx(1, 2, cs, div, div), cs);
    subtract_matrix(div, div, T2, div, M7, div, C + idx(2, 1, cs, div, div), cs);
    add_matrix(div, div, T2, div, M5, div, C + idx(2, 2, cs, div, div), cs);
    free(S1); free(S2); free(S3);free(S4); free(S5);free(S6);free(S7);free(S8);
    free(M1);free(M2);free(M3);free(M4);free(M5);free(M6);free(M7);
    free(T1);free(T2);
}

/**
 * This method is to get closest even value for the size of the
 * sub matrix
 *
 * */
int get_even(int n) {
    if (n % 2 == 0) {
        return n / 2;
    } else {
        return (n + 1) / 2;
    }
}

/**
 * This method is to do the dynamic padding of sub matrices
 *
 * */
void fill_matrices(int rows, int cols, int Am, int An, int *A, int *A11, int *A12, int *A21, int *A22) {
    int tempCols = 0;
    int tempRows  = 0;
    int i, j;
    for (i = 0; i < 2 * rows; i++) {
        for (j = 0; j < 2 * cols; j++) {
            if (i < rows) {
                 if (j < cols) {
                     *(A11 + cols * i + j) = *(A + i * An + j);
                 } else {
                     tempCols = j - cols;    
                     if (j < An) {
                         *(A12 + cols * i + tempCols) = *(A + i * An + j);
                     } else {
                         *(A12 + cols * i + tempCols) = 0;
                     }
                 }
            } else {
                 tempRows = i - rows; 
                 if (i < Am) {
                     if (j < cols) {
                         *(A21 + cols * tempRows + j) = *(A + i * An + j);         
                     } else {
                         if (j < An) {
                             *(A22 + cols * tempRows + j - cols) = *(A + i * An + j);  
                         } else {
                             *(A22 + cols * tempRows + j - cols) = 0;  
                         }
                     }
                 } else {
                     if (j < cols) {
                         *(A21 + cols * tempRows + j) = 0;         
                     } else {
                         *(A22 + cols * tempRows + j - cols) = 0;  
                     }
                 }
            }
        }
    }
}

/**
 * This method is to combine sub matrices to get the resultant matrix C
 *
 * */
void combine_matrices(int rows, int cols, int m, int n, int *A11, int *A12, int *A21, int *A22, int *C) {
    int i, j;
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            *(C + i * n + j) = *(A11 + i * cols + j);              
        }
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (j + cols < n) {
                *(C + i * n + j + cols) = *(A12 + i * cols + j);              
            }
        }
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (i + rows < m) {
                *(C + (i + rows)* n + j) = *(A21 + i * cols + j);              
            }
        }
    }
    for (i = 0; i < rows; i++) {
        for (j = 0; j < cols; j++) {
            if (i + rows < m && j + cols < n) {
                *(C + (i + rows)* n + j + cols) = *(A22 + i * cols + j);              
            }
        }
    }
}

/**
 * This is the Strassen's algorithm with dynamic padding
 *
 * */
void strassen_dynamic_padding(int Am, int An, int Bn, int *A, int *B, int *C) {
    if (Am <= CUTOFF || An < CUTOFF || Bn < CUTOFF) {
        Product_basic(A, Am, An, B, Bn, C);
        return;
    }

    int row = get_even(Am);
    int col = get_even(An);
    int *A11 = (int *)malloc(sizeof(int) * row * col);
    int *A12 = (int *)malloc(sizeof(int) * row * col);
    int *A21 = (int *)malloc(sizeof(int) * row * col);
    int *A22 = (int *)malloc(sizeof(int) * row * col);

    int bRow = get_even(An);
    int bCol = get_even(Bn);
    int *B11 = (int *)malloc(sizeof(int) * bRow * bCol);
    int *B12 = (int *)malloc(sizeof(int) * bRow * bCol);
    int *B21 = (int *)malloc(sizeof(int) * bRow * bCol);
    int *B22 = (int *)malloc(sizeof(int) * bRow * bCol);

    fill_matrices(row, col, Am, An, A, A11, A12, A21, A22);
    fill_matrices(bRow, bCol, An, Bn, B, B11, B12, B21, B22);
    int *S1 = (int *)malloc(sizeof(int) * row * col);
    add_matrix(row, col, A21, col, A22, col, S1, col);
    int *S2 = (int *)malloc(sizeof(int) * row * col);
    subtract_matrix(row, col, S1, col, A11, col, S2, col);
    int *S3 = (int *)malloc(sizeof(int) * row * col);
    subtract_matrix(row, col, A11, col, A21, col, S3, col);
    int *S4 = (int *)malloc(sizeof(int) * row * col);
    subtract_matrix(row, col, A12, col, S2, col, S4, col);

    int *S5 = (int *)malloc(sizeof(int) * bRow * bCol);
    subtract_matrix(bRow, bCol, B12, bCol, B11, bCol, S5, bCol);
    int *S6 = (int *)malloc(sizeof(int) * bRow * bCol);
    subtract_matrix(bRow, bCol, B22, bCol, S5, bCol, S6, bCol);
    int *S7 = (int *)malloc(sizeof(int) * bRow * bCol);
    subtract_matrix(bRow, bCol, B22, bCol, B12, bCol, S7, bCol);
    int *S8 = (int *)malloc(sizeof(int) * bRow * bCol);
    subtract_matrix(bRow, bCol, S6, bCol, B21, bCol, S8, bCol);

    int *M1 = (int *)malloc(sizeof(int) * row * bCol);
    strassen_dynamic_padding(row, col, bCol, S2, S6, M1);
    int *M2 = (int *)malloc(sizeof(int) * row * bCol);
    strassen_dynamic_padding(row, col, bCol, A11, B11, M2);
    int *M3 = (int *)malloc(sizeof(int) * row * bCol);
    strassen_dynamic_padding(row, col, bCol, A12, B21, M3);
    int *M4 = (int *)malloc(sizeof(int) * row * bCol);
    strassen_dynamic_padding(row, col, bCol, S3, S7, M4);
    int *M5 = (int *)malloc(sizeof(int) * row * bCol);
    strassen_dynamic_padding(row, col, bCol, S1, S5, M5);
    int *M6 = (int *)malloc(sizeof(int) * row * bCol);
    strassen_dynamic_padding(row, col, bCol, S4, B22, M6);
    int *M7 = (int *)malloc(sizeof(int) * row * bCol);
    strassen_dynamic_padding(row, col, bCol, A22, S8, M7);
      
    int *T1 = (int *)malloc(sizeof(int) * row * bCol);
    add_matrix(row, bCol, M1, bCol, M2, bCol, T1, bCol);
    int *T2 = (int *)malloc(sizeof(int) * row * bCol);
    add_matrix(row, bCol, T1, bCol, M4, bCol, T2, bCol);

    int *C11 = (int *)malloc(sizeof(int) * row * bCol);
    int *C12 = (int *)malloc(sizeof(int) * row * bCol);
    int *C21 = (int *)malloc(sizeof(int) * row * bCol);
    int *C22 = (int *)malloc(sizeof(int) * row * bCol);
    
    add_matrix(row, bCol, M2, bCol, M3, bCol, C11, bCol);
    add_matrix(row, bCol, T1, bCol, M5, bCol, C12, bCol);
    add_matrix(row, bCol, C12, bCol, M6, bCol, C12, bCol);
    subtract_matrix(row, bCol, T2, bCol, M7, bCol, C21, bCol);
    add_matrix(row, bCol, T2, bCol, M5, bCol, C22, bCol);
    combine_matrices(row, bCol, Am, Bn, C11, C12, C21, C22, C);
    free(A11); free(A12); free(A21); free(A22);
    free(B11); free(B12); free(B21); free(B22);
    free(S1); free(S2); free(S3);free(S4);free(S5);free(S6);free(S7);free(S8);
    free(M1);free(M2);free(M3);free(M4);free(M5);free(M6);free(M7);
    free(T1);free(T2);
    free(C11);free(C12);free(C21);free(C22);
}

/**
 * Before peeling, we need to split the original matrix in to 4 parts. We use the same input matrix
 * for the matrix with the size of m-1 * n-1 which is the biggest matrix. Two of them are
 * thin matrices and the other one is a single element matrix
 */
void fill_matrices_peeling_pre(int rows, int cols, int Am, int An, int *A, int as, int *A12, int *A21, int *A22) {
    int tempCols = 0;
    int tempRows  = 0;
    int i, j;
    for (i = 0; i < Am ; i++) {
		for (j = 0; j < An; j++) {
			if (i < rows) {
				if (j >= cols) {
					tempCols = j - cols;
					*(A12 + i) = *(A + i * as + j);
				}
			} else {
				tempRows = i - rows;
				if (i < Am) {
					if (j < cols) {
						*(A21 + cols * tempRows + j) = *(A + i * as + j);
					} else {
						*(A22 + cols * tempRows + j - cols) = *(A + i * as + j);

					}
				} else {
					if (j < cols) {
						*(A21 + cols * tempRows + j) = 0;
					} else {
						*(A22 + cols * tempRows + j - cols) = 0;
					}
				}
			}
		}
    }
}

/**
 * This method does strassen's algorithm according to dynamic peeling approach.
 */
void strassen_dynamic_peeling(int Am, int An, int Bn, int *A, int *B, int *C, int as, int bs, int cs) {
	if (Am <= CUTOFF || An < CUTOFF || Bn < CUTOFF) {
		product_basic_with_stride(Am, An, A, as, B, Bn,bs, C,cs);
		return;
	}
	int a11row, a11col, a12row, a12col, a21row, a21col, a22row, a22col;
	int b11row, b11col, b12row, b12col, b21row, b21col, b22row, b22col;
	if (Am % 2 == 0) {
		a11row = Am; a12row = Am;
		a21row = 0; a22row = 0;
	} else {
		a11row = Am - 1; a12row = Am - 1;
		a21row = 1; a22row = 1;
	}

	if (An % 2 == 0) {
		a11col = An; a12col = 0;
		a21col = An; a22col = 0;
		b11row = An; b12row = An;
		b21row = 0; b22row = 0;
	} else {
		a11col = An - 1; a12col = 1;
		a21col = An - 1; a22col = 1;
		b11row = An - 1; b12row = An - 1;
		b21row = 1; b22row = 1;
	}

	if (Bn % 2 == 0) {
		b11col = Bn; b12col = 0;
		b21col = Bn; b22col = 0;
	} else {
		b11col = Bn - 1; b12col = 1;
		b21col = Bn - 1; b22col = 1;
	}

	int *A12 = (int *) malloc(sizeof(int) * a12row * a12col);
	int *A21 = (int *) malloc(sizeof(int) * a21row * a21col);
	int *A22 = (int *) malloc(sizeof(int) * a22row * a22col);

	int *B12 = (int *) malloc(sizeof(int) * b12row * b12col);
	int *B21 = (int *) malloc(sizeof(int) * b21row * b21col);
	int *B22 = (int *) malloc(sizeof(int) * b22row * b22col);

	int *ABC3 = (int *) malloc(sizeof(int) * a11row * b12col);
	int *ABC4 = (int *) malloc(sizeof(int) * a11row * b22col);
	int *ABC5 = (int *) malloc(sizeof(int) * a21row * b11col);
	int *ABC6 = (int *) malloc(sizeof(int) * a22row * b21col);
	int *ABC7 = (int *) malloc(sizeof(int) * a21row * b12col);
	int *ABC8 = (int *) malloc(sizeof(int) * a22row * b22col);

	fill_matrices_peeling_pre(a11row, a11col, Am, An, A, as, A12, A21, A22);
	fill_matrices_peeling_pre(b11row, b11col, An, Bn, B, bs, B12, B21, B22);

	int row = a11row / 2;
	int col = a11col / 2;
	int bRow = b11row / 2;
	int bCol = b11col / 2;

	int *S1 = (int *) malloc(sizeof(int) * row * col);
	add_matrix(row, col, A + idx(2, 1, as, row, col), as, A + idx(2, 2, as, row, col), as, S1, col);
	int *S2 = (int *) malloc(sizeof(int) * row * col);
	subtract_matrix(row, col, S1, col, A + idx(1, 1, as, row, col), as, S2, col);
	int *S3 = (int *) malloc(sizeof(int) * row * col);
	subtract_matrix(row, col, A + idx(1, 1, as, row, col), as,  A + idx(2, 1, as, row, col), as, S3, col);
	int *S4 = (int *) malloc(sizeof(int) * row * col);
	subtract_matrix(row, col, A + idx(1, 2, as, row, col), as, S2, col, S4, col);

	int *S5 = (int *) malloc(sizeof(int) * bRow * bCol);
	subtract_matrix(bRow, bCol, B + idx(1, 2, bs, bRow, bCol), bs, B + idx(1, 1, bs, bRow, bCol), bs, S5, bCol);
	int *S6 = (int *) malloc(sizeof(int) * bRow * bCol);
	subtract_matrix(bRow, bCol,  B + idx(2, 2, bs,bRow, bCol), bs, S5, bCol, S6, bCol);
	int *S7 = (int *) malloc(sizeof(int) * bRow * bCol);
	subtract_matrix(bRow, bCol, B + idx(2, 2, bs, bRow, bCol), bs, B + idx(1, 2, bs, bRow, bCol), bs, S7, bCol);
	int *S8 = (int *) malloc(sizeof(int) * bRow * bCol);
	subtract_matrix(bRow, bCol, S6, bCol, B + idx(2, 1, bs, bRow, bCol), bs, S8, bCol);

	int *M1 = (int *) malloc(sizeof(int) * row * bCol);
	strassen_dynamic_peeling(row, col, bCol, S2, S6, M1, col, bCol, bCol);
	int *M2 = (int *) malloc(sizeof(int) * row * bCol);
	strassen_dynamic_peeling(row, col, bCol, A + idx(1, 1, as, row, col), B + idx(1, 1, bs, bRow, bCol), M2, as, bs, bCol);
	int *M3 = (int *) malloc(sizeof(int) * row * bCol);
	strassen_dynamic_peeling(row, col, bCol, A + idx(1, 2, as, row, col), B + idx(2, 1, bs, bRow, bCol), M3, as, bs, bCol);
	int *M4 = (int *) malloc(sizeof(int) * row * bCol);
	strassen_dynamic_peeling(row, col, bCol, S3, S7, M4, col, bCol, bCol);
	int *M5 = (int *) malloc(sizeof(int) * row * bCol);
	strassen_dynamic_peeling(row, col, bCol, S1, S5, M5, col, bCol, bCol);
	int *M6 = (int *) malloc(sizeof(int) * row * bCol);
	strassen_dynamic_peeling(row, col, bCol, S4, B + idx(2, 2, bs, bRow, bCol), M6, col, bs, bCol);
	int *M7 = (int *) malloc(sizeof(int) * row * bCol);
	strassen_dynamic_peeling(row, col, bCol, A + idx(2, 2, as, row, col), S8, M7, as, bCol, bCol);

	int *T1 = (int *) malloc(sizeof(int) * row * bCol);
	add_matrix(row, bCol, M1, bCol, M2, bCol, T1, bCol);
	int *T2 = (int *) malloc(sizeof(int) * row * bCol);
	add_matrix(row, bCol, T1, bCol, M4, bCol, T2, bCol);

	add_matrix(row, bCol, M2, bCol, M3, bCol,  C + idx(1, 1, cs, row, bCol), cs);
	add_matrix(row, bCol, T1, bCol, M5, bCol, C + idx(1, 2, cs, row, bCol), cs);
	add_matrix(row, bCol, C + idx(1, 2, cs, row, bCol), cs, M6, bCol,  C + idx(1, 2, cs, row, bCol), cs);
	subtract_matrix(row, bCol, T2, bCol, M7, bCol, C + idx(2, 1, cs, row, bCol), cs);
	add_matrix(row, bCol, T2, bCol, M5, bCol, C + idx(2, 2, cs, row, bCol), cs);

	product_basic_with_stride_adding(a12row, a12col, A12, a12col, B21, b21col,b21col, C, Bn);
	product_basic_with_stride(a11row, a11col, A, as, B12, b12col,b12col, ABC3, b12col);
	Product_basic(A12, a12row, a12col, B22, b22col, ABC4);
	product_basic_with_stride(a21row,a21col, A21, a21col, B, b11col,bs, ABC5, b11col);
	Product_basic(A22, a22row, a22col, B21, b21col, ABC6);
	Product_basic(A21, a21row, a21col, B12, b12col, ABC7);
	Product_basic(A22, a22row, a22col, B22, b22col, ABC8);

	add_matrix(a11row, b12col, ABC3, b12col, ABC4, b22col, C + Bn - 1,cs);
	add_matrix(a21row, b11col, ABC5, b11col, ABC6, b21col, C + a11row * Bn ,cs);
	add_matrix(a21row, b22col, ABC7, b12col, ABC8, b22col, C + Am * Bn - 1 ,cs);

	free(A12); free(A21); free(A22);
	free(B12); free(B21); free(B22);
	free(ABC3); free(ABC4);free(ABC5); free(ABC6); free(ABC7); free(ABC8);
	free(S1); free(S2); free(S3); free(S4); free(S5); free(S6); free(S7); free(S8);
	free(M1); free(M2); free(M3); free(M4); free(M5); free(M6);	free(M7);
	free(T1); free(T2);
}



