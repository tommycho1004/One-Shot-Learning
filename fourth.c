#include <stdio.h>
#include <stdlib.h>

void transposeMat(int m, int n, double mat[m][n], double transpose[n][m]) 
{
	//will transpose any mxn matrix into its nxm transpose
        if (m < 1 || n < 1) 
	{
                printf("cannot transpose an empty matrix\n");
                return;
        }
        int i, j;
        for (i = 0; i < m; i++) 
	{
                for (j = 0; j < n; j++) 
		{
                        transpose[j][i] = mat[i][j];
                }
        }
}

void matrixMult(int m1, int n1, double mat1[m1][n1], int m2, int n2, double mat2[m2][n2], double product[m1][n2])
{
	//multiplies an mxn matrix by an nxk matrix to create an mxk matrix in an inputted product matrix
        if (n1 != m2 || m1 < 1 || m2 < 1 || n1 < 1 || n2 < 1) 
	{
                printf("incompatible matrices\n");
                return;
        }
        int i, j, k;
        for (i = 0; i < m1; i++)
	{
        	for (i = 0; i < m1; i++)
		{
                	for (j = 0; j < n2; j++) 
			{
                        	product[i][j] = 0;
                        	for (k = 0; k < n1; k++)
				{
                                	product[i][j] += mat1[i][k] * mat2[k][j];
                        	}
                	}
        	}
	}
}



void getIdentity(int size, double mat[size][size])
{
	//will find identity matrix of any given size n (nxn identity)
        if (size < 1) 
	{
                printf("no identity for an empty matrix\n");
        }
        int i, j;
        for (i = 0; i < size; i++) 
	{
                for (j = 0; j < size; j++) 
		{
			if (i==j)
			{
				mat[i][j] = 1;
			}
			else
			{
				mat[i][j] = 0;
			}	                         
                }
        }
}


void getInverse(int size, double mat[size][size], double inverse[size][size]) 
{
        if (size < 1)
	{
                printf("no inverse for an empty matrix\n");
                return;
        }
	//we need to get an identity matrix first to do gaussian elimination to find inverse
        getIdentity(size, inverse);
	//create local vars
        int pRow, cRow, cCol;
	int pCol = 0;
        for (pRow = 0; pRow < size; pRow++, pCol++) 
	{
                double oPivot = (mat[pRow][pCol] == 0) ? 1 : mat[pRow][pCol];
                for (cCol = 0; cCol < size; cCol++)
		{
                        mat[pRow][cCol] /= oPivot;
                        inverse[pRow][cCol] /= oPivot;
                }
                for (cRow = pRow+1; cRow < size; cRow++) 
		{
                        double c = - (mat[cRow][pCol] / mat[pRow][pCol]);
                        for (cCol = 0; cCol < size; cCol++)
			{
                                mat[cRow][cCol] += c*mat[pRow][cCol];
                                inverse[cRow][cCol] += c*inverse[pRow][cCol];
                        }
                }
        }
        for (pCol= size-1, pRow = size-1; pRow >= 0; pCol--,pRow--)
	{
                for (cRow = pRow-1; cRow >= 0; cRow--)
		{
                        double c = - (mat[cRow][pCol] / mat[pRow][pCol]);
                        for (cCol = size-1; cCol >= 0; cCol--)
			{
                                mat[cRow][cCol] += c*mat[pRow][cCol];
                                inverse[cRow][cCol] += c*inverse[pRow][cCol];
                        }
                }
        }
}

void printPrices(int size, double vector[size][1]) 
{
        //hopefully will print prices found in the final prediction vector
	int i;
        for (i = 0; i < size; i++)
	{
                printf("%0.0lf\n", vector[i][0]);
        }
}


int main(int argc, char *argv[])
{
	//check for all complete files
        if (argc != 3) 
	{
                return -1;
        }
	//open training file
        FILE *fptr = fopen(argv[1], "r");
        if (fptr == NULL) 
	{
                return -1;
        }
	//w=# of weights, a = # of attributes, 
        int w, a, num;
        fscanf(fptr, "%d\n", &a);
        w = a+1;
        fscanf(fptr, "%d\n", &num);
	//create matrices to hold attributes, prices, and eventually weights
        double temp[num][w];
	double prices[num][1];
        double weight[w][1];
        int i, j;
	//create temp matrix with given values
        for (i = 0; i < num; i++)
	{
                for (j = 0; j < w; j++)
		{
                        double n;
                        if (j == 0)
			{
                                temp[i][j] = 1;
                                fscanf(fptr, "%lf ,", &n);
                                prices[i][0] = n;
                        }
                        else if (j == a)
			{
                                fscanf(fptr, "%lf\n", &n);
                                temp[i][j] = n;
                        }
                        else 
			{
                                fscanf(fptr, "%lf ,", &n);
                                temp[i][j] = n;
                        }
                }
        }
	//now we're done with this file
        fclose(fptr);

	//open test file
        FILE *tptr = fopen(argv[2], "r");
        if (tptr == NULL) 
	{
                return -1;
        }
        int num2;
        fscanf(tptr, "%d\n", &num2);
	//create test matrix to store test values
        double test[num2][w];
	//similar process to above
        for (i = 0; i < num2; i++)
	{
                for (j = 0; j < w; j++) 
		{
                        double n;
                        if (j == 0) 
			{
                                test[i][j] = 1;
                        }
                        else if (j == a)
			{
                                fscanf(tptr, "%lf\n", &n);
                                test[i][j] = n;
                        }
                        else 
			{
                                fscanf(tptr, "%lf ,", &n);
                                test[i][j] = n;
                        }
                }
        }
	//close test file
        fclose(tptr);
	
	// W=(X^T X)^-1 X^T  Y
	//this is the equation we want to work towards
	//create a transpose matrix with opposite configuration of temp matrix
        double transpose[w][num];
	//now we transpose our temp matrix onto transpose
        transposeMat(num, w, temp, transpose);
	//to find inner product
        double inner[w][w];
	//multiply transpose matrix with temp matrix (we're doing X^T X here)
        matrixMult(w, num, transpose, num, w, temp, inner);
	//now get inverse of this whole thing ((X^T X)^-1)
        double inverse[w][w];
        getInverse(w, inner, inverse);
	//now we want to get (X^T X)^-1 X^T
        double product1[w][num];
        matrixMult(w, w, inverse, w, num, transpose, product1);
	//to finish we reintroduce the weight matrix we made earlier and put our new values in there
        matrixMult(w, num, product1, num, 1, prices, weight);
	//create matrix for our predictions
        double prediction[num2][1];
	//multiply our test points with our found weights to get predicted prices
        matrixMult(num2, w, test, w, 1, weight, prediction);
	//print prices
        printPrices(num2, prediction);

        return 0;
}
