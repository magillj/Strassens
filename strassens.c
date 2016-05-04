/*
 * Developed by: Derek Coley & Jake Magill
 * Provides an implementation of Strassen's algorithm
 * and uses the best threshold value for a hybrid
 * algorithm of matrix multiplication over integers
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>
#define THRESHOLD 2

/* Asserts that the passed in matrix is the same
 * as a matrix computed using the naive algorithm */
void checkCorrect(int**, int**, int**, int, int);

/* Uses Strassen's algorithm to multiply together the
 * two passed in matrices and return a result matrix */
int** strassens(int**, int, int**, int);

/* Uses the naive multiplication algorithm to multiply together 
 * the two passed in matrices and return a result matrix */
int** naive(int**, int, int**, int);

/* Prints the passed in matrix to stdout */
void print(int**, int, int);

/* Creates a matrix with values 1, ..., n and returns it */
int** createMatrix(int, int);

/* Helper function to compute the value for a single
 * cell in the matrix for the naive algorithm */
int computeCell(int**, int**, int, int, int);

/* Adds or subtracts the two passed in matrices based
 * on the flag passed in and returns a result matrix */
int** matrixOp(int**, int**, int, int);

int main(int argc, char** argv) {
  int size = 32;
  int** a = createMatrix(size, size);
  int** b = createMatrix(size, size);
  int** c = strassens(a, size, b, size);
  checkCorrect(a, b, c, size, size);
  return EXIT_SUCCESS;
}

int** matrixOp(int** m1, int** m2, int size, int add) {
  int** m3 = createMatrix(size, size);
  for (int i = 0; i < size; i++) {
    for (int j = 0; j < size; j++) {
      m3[i][j] = (add) ? m1[i][j] + m2[i][j] :
          m1[i][j] - m2[i][j];
    }
  }
  return m3;
}

int** createMatrix(int row, int col) {
  int** a = (int**) malloc(row * sizeof(int*));
  for (int i = 0; i < col; i++) {
    a[i] = (int *) malloc(col * sizeof(int));
  }
  int count = 0;
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      a[i][j] = count++;
    }
  }
  return a;
}

int** naive(int** m1, int m1Size, int** m2, int m2Size) {
  int** m3 = createMatrix(m1Size, m1Size);
  for (int i = 0; i < m1Size; i++) {
    for (int j = 0; j < m1Size; j++) {
      m3[i][j] = computeCell(m1, m2, m1Size, i, j);
    }
  }
  return m3;
}

int computeCell(int** m1, int** m2, int len,
                int row, int col) {
  int result = 0;
  for (int i = 0; i < len; i++) {
    result += m1[row][i] * m2[i][col];
  }
  return result;
}

void print(int** grid, int row, int col) {
  for (int i = 0; i < row; i++) {
    for (int j = 0; j < col; j++) {
      printf("%d ", grid[i][j]);
    }
    printf("\n");
  }
  printf("\n");
}

void checkCorrect(int** m1, int** m2, int** test, int m1Size, int m2Size) {
  int** m3 = naive(m1, m1Size, m2, m2Size);
  for (int i = 0; i < m1Size; i++) {
    for (int j = 0; j < m1Size; j++) {
      if (test[i][j] != m3[i][j]) {
        printf("expected: \n");
        print(m3, m1Size, m1Size);
        printf("actual: \n");
        print(test, m1Size, m1Size);
        printf("failed :( \n");
      }
      assert(test[i][j] == m3[i][j]);
    }
  }
  printf("passed :) \n");
}

int** strassens(int** m1, int n2, int** m2, int n) {
  if (n == THRESHOLD) {
    return naive(m1, n2, m2, n);
  } else {
    /* Split m1 and m2 into four quadrants */
    int** a11 = createMatrix(n/2, n/2);
    int** a12 = createMatrix(n/2, n/2);
    int** a21 = createMatrix(n/2, n/2);
    int** a22 = createMatrix(n/2, n/2);
    int** b11 = createMatrix(n/2, n/2);
    int** b12 = createMatrix(n/2, n/2);
    int** b21 = createMatrix(n/2, n/2);
    int** b22 = createMatrix(n/2, n/2);
    int row;
    int col;
    row = col = 0;

    /* Fill quadrant matrices with appropriate values */
    for (int i = 0; i < n/2; i++) {
      for (int j = 0; j < n/2; j++) {
        a11[row][col] = m1[i][j];
        b11[row][col] = m2[i][j];
        col++;
      }
      col = 0;
      row++;
    }
    row = col = 0;
    
    for (int i = 0; i < n/2; i++) {
      for (int j = n/2; j < n; j++) {
        a12[row][col] = m1[i][j];
        b12[row][col++] = m2[i][j];
      }
      col = 0;
      row++;
    }
    row = col = 0;
    
    for (int i = n/2; i < n; i++) {
      for (int j = 0; j < n/2; j++) {
        a21[row][col] = m1[i][j];
        b21[row][col++] = m2[i][j];
      }
      col = 0;
      row++;
    }
    row = col = 0;
    
    for (int i = n/2; i < n; i++) {
      for (int j = n/2; j < n; j++) {
        a22[row][col] = m1[i][j];
        b22[row][col++] = m2[i][j];
      }
      col = 0;
      row++;
    }

    /* Use four quadrants and seven multiplications to compute p1, ..., p7 */
    int** p1 = strassens(a12, n/2, matrixOp(b11, b21, n/2, 1), n/2);
    int** p2 = strassens(a21, n/2, matrixOp(b12, b22, n/2, 1), n/2);
    int** p3 = strassens(matrixOp(a11, a12, n/2, 0), n/2, b11, n/2);
    int** p4 = strassens(matrixOp(a22, a21, n/2, 0), n/2, b22, n/2);
    int** p5 = strassens(matrixOp(a22, a12, n/2, 0), n/2, matrixOp(b21, b22, n/2, 0), n/2);
    int** p6 = strassens(matrixOp(a11, a21, n/2, 0), n/2, matrixOp(b12, b11, n/2, 0), n/2);
    int** p7 = strassens(matrixOp(a21, a12, n/2, 0), n/2, matrixOp(b11, b22, n/2, 1), n/2);

    /* Build four quadrants of result matrix */
    int** c = createMatrix(n, n);
    int** c11 = matrixOp(p1, p3, n/2, 1);
    int** c12 = matrixOp(matrixOp(p2, p3, n/2, 1),
                         matrixOp(p6, p7, n/2, 0),
                         n/2, 1);
    int** c21 = matrixOp(matrixOp(p1, p4, n/2, 1),
                         matrixOp(p5, p7, n/2, 1),
                         n/2, 1);
    int** c22 = matrixOp(p2, p4, n/2, 1);
    row = col = 0;

    /* Fill result matrix to return */
    for (int i = 0; i < n/2; i++) {
      for (int j = 0; j < n/2; j++) {
        c[i][j] = c11[row][col++];
      }
      col = 0;
      row++;
    }
    row = col = 0;
    
    for (int i = 0; i < n/2; i++) {
      for (int j = n/2; j < n; j++) {
        c[i][j] = c12[row][col++];
      }
      col = 0;
      row++;
    }
    row = col = 0;
    
    for (int i = n/2; i < n; i++) {
      for (int j = 0; j < n/2; j++) {
        c[i][j] = c21[row][col++];
      }
      col = 0;
      row++;
    }
    row = col = 0;
    
    for (int i = n/2; i < n; i++) {
      for (int j = n/2; j < n; j++) {
        c[i][j] = c22[row][col++];
      }
      col = 0;
      row++;
    }

    /* Return result :) */
    return c;
  }
}











