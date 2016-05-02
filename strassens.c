/*
 * Developed by: Derek Coley & Jake Magill
 * Provides an implementation of Strassen's algorithm
 * and uses the best threshold value for a hybrid
 * algorithm of matrix multiplication over integers
 */
#include <stdlib.h>
#include <stdio.h>
#include <assert.h>

void checkCorrect(int**, int**, int**, int, int);
int** strassens(int**, int**, int**, int, int, int);
int** naive(int**, int, int**, int);
void print(int**, int, int);
int** createMatrix(int, int);
int computeCell(int**, int**, int, int, int);
int** matrixOp(int**, int**, int, int);

int main(int argc, char** argv) {
  int size = 4;
  int** a = createMatrix(size, size);
  int** b = createMatrix(size, size);
  int** c = createMatrix(size, size);
  strassens(a, b, c, 0, 0, size);
  print(c, size, size);
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
}

void checkCorrect(int** m1, int** m2, int** test, int m1Size, int m2Size) {
  int** m3 = naive(m1, m1Size, m2, m2Size);
  print(m3, m1Size, m1Size);
  for (int i = 0; i < m1Size; i++) {
    for (int j = 0; j < m1Size; j++) {
      assert(test[i][j] == m3[i][j]);
    }
  }
}

int** strassens(int** m1, int** m2, int** result, int i, int j, int n) {
  if (n == 2 || 1) {
    int** result = createMatrix(2, 2);
    int p1 = m1[i][j+1] * (m2[i][j] + m2[i+1][j]);
    int p2 = m1[i+1][j] * (m2[i][j+1] + m2[i+1][j+1]);
    int p3 = (m1[i][j] - m1[i][j+1]) * m2[i][j];
    int p4 = (m1[i+1][j+1] - m1[i+1][j]) * m2[i+1][j+1];
    int p5 = (m1[i+1][j+1] - m1[i][j+1]) * (m2[i+1][j] - m2[i+1][j+1]);
    int p6 = (m1[i][j] - m1[i+1][j]) * (m2[i][j+1] - m2[i][j]);
    int p7 = (m1[i+1][j] - m1[i][j+1]) * (m2[i][j] + m2[i+1][j+1]);

    result[i][j] = p1 + p3;
    result[i][j+1] = p2 + p3 + p6 - p7;
    result[i+1][j] = p1 + p4 + p5 + p7;
    result[i+1][j+1] = p2 + p4;
    return result;
  } else {
    // Top left
    strassens(m1, m2, result, i, j, n/2);
    // Top right
    strassens(m1, m2, result, i + n/2, j, n/2);
    // Bottom left
    strassens(m1, m2, result, i, j + n/2, n/2);
    // Bottom right
    strassens(m1, m2, result, i + n/2, j + n/2, n/2);

    // Merge
  }
}











