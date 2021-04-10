//==========================//
//Student Name: Niyazi Ulke
//Student Number: 2017400114
//Compile Status: Compiling
//Program Status: Working
//==========================//
#include <iostream>
#include <fstream>
#include "mpi.h"
using namespace std;

/*
 * Merges two sorted arrays: arr1 and arr2, then places the result to merged array.
 * size1 and size2 are number of elements in arr1 and arr2.
 * Returns number of elements (length) of merged array.
 */
int merge(int arr1[], int arr2[], int merged[], int size1, int size2) {
	int i1 = 0; // index of arr1
	int i2 = 0; // index of arr2
	int j = 0;  // index of merged

	while (i1 < size1 && i2 < size2) {
		if (arr1[i1] < arr2[i2]) {
			merged[j++] = arr1[i1++];

		} else if (arr1[i1] == arr2[i2]) {
			merged[j++] = arr1[i1++];
			i2++;

		} else {
			merged[j++] = arr2[i2++];

		}
	}
	while (i1 < size1) { // If some elements are left in arr1, append them to merged array.
		merged[j++] = arr1[i1++];

	}
	while (i2 < size2) { // If some elements are left in arr2, append them to merged array.
		merged[j++] = arr2[i2++];

	}

	return j;
}
/*
 * Sorts weights and chooses the maximum T weighted features.
 * Then sorts feature ids in ascending order.
 * Places the highest weighted T feature ids to choose array.
 * A is number of features.
 */
void sortChoose(pair<double, int> weight[], int A, int T, int choose[]) {
	/*
	 * Sort by weights in descending order, insertion sort
	 */
	for (int i = 1; i < A; i++) {
		int j = i - 1;
		pair<double, int> compare(weight[i]);
		while (j >= 0) {
			if (weight[j].first < compare.first) {
				weight[j + 1] = weight[j];
				j--;
			} else {
				break;
			}
		}
		weight[j + 1] = compare;

	}
	/*
	 * Fill choose array with feature ids of the highest weighted T features.
	 * Insert feature ids in a sorted manner.
	 */
	for (int i = 0; i < T; i++) {
		choose[i] = weight[i].second;
		int compare = choose[i];
		int j = i - 1;
		while (j >= 0) {
			if (choose[j] > compare) {
				choose[j + 1] = choose[j];
				j--;
			} else {
				break;
			}
		}
		choose[j + 1] = compare;
	}

}

/*
 * Returns absolute value of a double.
 */
double absval(double num) {
	return (num < -num) ? -num : num;
}
/*
 * Calculates manhattanDistance between arr1 and arr2
 * A is length of arr1 and arr2, equivalently number of features.
 */
double manhattanDistance(double arr1[], double arr2[], int A) {

	double sum = 0;
	for (int i = 0; i < A; i++) {
		sum += absval(arr1[i] - arr2[i]);
	}
	return sum;
}
int main(int argc, char **argv) {
	int rank; // rank of the current processor
	int size; // total number of processors
	MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank); // gets the rank of the current processor
	MPI_Comm_size(MPI_COMM_WORLD, &size); // gets the total number of processors

	if (rank == 0) {
		string inputfile = argv[1];
		ifstream infile(inputfile);
		int P, N, A, M, T; // read parameters from the file.
		/*
		 * P is number of processors
		 * N is number of lines
		 * A is number of features
		 * M is number of iterations in relief algorithm
		 * T is number of features that each slave processor will choose
		 */
		infile >> P >> N >> A >> M >> T;
		int numLines = N / (P - 1); // number of lines that will be sended to each processor.
		/*
		 * Broadcast important parameters to all processors
		 */
		MPI_Bcast(&numLines, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&A, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD);

		double val;
		for (int i = 0; i < N; i++) {
			for (int j = 0; j < A + 1; j++) {
				// read values and send to the corresponding processors.
				infile >> val;
				MPI_Send(&val, 1, MPI_DOUBLE, (i / numLines + 1), 0,
						MPI_COMM_WORLD);
			}

		}
		int store[P * T];
		int size1 = T;
		int merged[P * T];
		for (int i = 0; i < T; i++)
			MPI_Recv(&store[i], 1, MPI_INT, 1, 0, MPI_COMM_WORLD,
					MPI_STATUS_IGNORE); // receive T integers as feature ids from slave 1
		int features[T];
		for (int i = 2; i < P; i++) {
			for (int j = 0; j < T; j++)
				MPI_Recv(&features[j], 1, MPI_INT, i, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE); // receive T integers as feature ids from slaves 2 to P-1

			size1 = merge(store, features, merged, size1, T); // size1 = number of elements in merged array
			for (int i = 0; i < size1; i++)
				store[i] = merged[i]; // copy from merged to store

		}
		string acc = "Master P" + to_string(rank) + " :"; // output string for the master

		for (int i = 0; i < size1; i++)
			acc += " " + to_string(store[i]); // prepare the output string

		cout << acc << endl;
	}

	else {
		int numLines, A, M, T;
		/* receive broadcasted information */
		MPI_Bcast(&numLines, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&A, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&M, 1, MPI_INT, 0, MPI_COMM_WORLD);
		MPI_Bcast(&T, 1, MPI_INT, 0, MPI_COMM_WORLD);
		double **myLines = new double*[numLines]; // stores instances sent to the processor

		for (int i = 0; i < numLines; i++)
			myLines[i] = new double[A + 1];

		for (int i = 0; i < numLines; i++) {
			for (int j = 0; j < A + 1; j++) {
				double val;
				MPI_Recv(&val, 1, MPI_DOUBLE, 0, 0, MPI_COMM_WORLD,
						MPI_STATUS_IGNORE); // receive numbers from the master

				myLines[i][j] = val;
			}

		}

		pair<double, int> weight[A]; // stores weight-featureid pairs
		double featureRanges[A]; // stores maximum - minimum of each feature
		for (int k = 0; k < A; k++) {
			double minval = myLines[0][k];
			double maxval = myLines[0][k];
			for (int i = 1; i < numLines; i++) { // find maximum and minimum values for each feature

				if (myLines[i][k] > maxval)
					maxval = myLines[i][k];

				if (myLines[i][k] < minval)
					minval = myLines[i][k];

			}
			featureRanges[k] = maxval - minval; // calculate maximum - minimum for each feature
			weight[k] = pair<double, int>(0, k); // initialize weights as 0 and initialize feature ids
		}

		for (int i = 0; i < M; i++) { // M is number of iterations of relief algorithm.
			int nearHit = -1; // index of nearest hit instance
			int nearMiss = -1; // index of nearest miss instance
			double hitDist = 0; // distance to nearest hit
			double missDist = 0; // distance to nearest miss

			for (int j = 0; j < numLines; j++) { // traverse all instances
				if (i == j) // skip the target instance
					continue;

				double dist = manhattanDistance(myLines[i], myLines[j], A);

				if (myLines[i][A] == myLines[j][A]) { // true if the instances are of the same class.

					if (hitDist == 0 || dist < hitDist) { // if it's the first hit distance calculation or the new minimum, update the variables
						hitDist = dist;
						nearHit = j;
					}

				} else {
					if (missDist == 0 || dist < missDist) { // if it's the first miss distance calculation or the new minimum, update the variables
						missDist = dist;
						nearMiss = j;
					}
				}

			}

			for (int k = 0; k < A; k++) { // update weights according to the relief algorithm

				weight[k].first -= absval(myLines[i][k] - myLines[nearHit][k])
						/ (featureRanges[k] * M);
				weight[k].first += absval(myLines[i][k] - myLines[nearMiss][k])
						/ (featureRanges[k] * M);
			}

		}
		int choose[T];
		string acc = "Slave P" + to_string(rank) + " :"; // output string for the slave
		sortChoose(weight, A, T, choose);
		for (int i = 0; i < T; i++) // prepare the output string for the slave
			acc += " " + to_string(choose[i]);

		cout << acc << endl; // print the output
		for (int i = 0; i < T; i++)
			MPI_Send(&choose[i], 1, MPI_INT, 0, 0, MPI_COMM_WORLD); // send T maximum weight features to the master

		for (int i = 0; i < numLines; i++)
			delete[] myLines[i]; // free dynamically allocated memory

		delete[] myLines; // free dynamically allocated memory
	}
	MPI_Barrier (MPI_COMM_WORLD);
	MPI_Finalize();
	return 0;

}

