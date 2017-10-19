#include "CSR_STRUCTURE.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>
#include <vector>
#include <algorithm>
#include <time.h>

//returns pointer to new transpose matrix
CSR_STRUCTURE* CreateTransposeMatrix(CSR_STRUCTURE &csr_matrix)
{

	int transpose_rows = csr_matrix.columns_;
	int transpose_cols = csr_matrix.rows_;
	int transpose_nnzs = csr_matrix.nonzeroes_;
	int *row_counts = new int[transpose_rows+1];

	/*overwrites any trash values that are in the array during allocation*/
	CSR_STRUCTURE transpose_matrix(transpose_rows, transpose_cols, transpose_nnzs);

	for (int i = 0; i < transpose_matrix.rows_ + 1; i++)
	{
		transpose_matrix.row_ptr_[i] = 0;
		row_counts[i] = 0;
	}//end for

	 /*counts number of non zero values in a column by incrementing the value in the row_ptr everytime it sees the same column*/
	 /*it also updates the number of non zeros in that column which turns into the # of nonzeros in that row of the transpose*/
	for (int i = 0; i < csr_matrix.rows_ + 1; i++)
	{
		for (int j = csr_matrix.row_ptr_[i]; j < csr_matrix.row_ptr_[i + 1]; j++)
		{
			//std::cout << "at column " << csr_matrix.col_ptr_[j] << " when j = " << j << " add one" << std::endl;
			int col_nnz = csr_matrix.col_ptr_[j];
			transpose_matrix.row_ptr_[col_nnz + 1]++;
		}
	}

	//std::cout << "before constructing row_ptr2" << std::endl;
	/*constructs the transpose matrix row_ptr and updates its values*/
	for (int i = 0; i < transpose_matrix.rows_ + 1; i++)
	{
		transpose_matrix.row_ptr_[i+1] += transpose_matrix.row_ptr_[i];
	}

	/*populates the new csr matrix with the transpose information*/
	for (int i = 0; i < csr_matrix.rows_; i++)
	{
		for (int j = csr_matrix.row_ptr_[i]; j < csr_matrix.row_ptr_[i + 1]; j++)
		{
			int i2 = csr_matrix.col_ptr_[j];
			transpose_matrix.col_ptr_[transpose_matrix.row_ptr_[i2] + row_counts[i2]] = i;
			transpose_matrix.nnz_val_[transpose_matrix.row_ptr_[i2] + row_counts[i2]] = csr_matrix.nnz_val_[j];
			row_counts[i2]++;
		}

	}

	return new CSR_STRUCTURE(transpose_matrix);
}

/*IO function that outputs any matrix to an output file*/
void OutputMatrix(CSR_STRUCTURE &csr_matrix, std::ofstream &matrix_of)
{
	matrix_of << csr_matrix.rows_ << " " << csr_matrix.columns_ << " " << csr_matrix.nonzeroes_ << "\n";

	for (int i = 0; i < csr_matrix.rows_; i++) {
		for (int j = csr_matrix.row_ptr_[i]; j < csr_matrix.row_ptr_[i + 1]; j++)
		{
			matrix_of << csr_matrix.col_ptr_[j] + 1 << " " << csr_matrix.nnz_val_[j] << " ";
		}
		matrix_of << "\n";
	}

	matrix_of.close();
}

//returns pointer to new csr_matrix
CSR_STRUCTURE* InputMatrixToCSR(std::ifstream &in_file)
{
	int counter = 0;
	int matrix_info[3];
	int info_val, int_val;
	std::string str_entire_line, str_info, str_val;
	getline(in_file, str_entire_line);
	std::istringstream first_line_ss(str_entire_line);

	/*grab first line which contains matrix information*/
	do {
		first_line_ss >> str_info;
		info_val = atoi(str_info.c_str());
		matrix_info[counter] = info_val;
		counter++;
	} while (first_line_ss);

	CSR_STRUCTURE newCSR(matrix_info[0], matrix_info[1], matrix_info[2]);

	int col_nnz_index_ptr = 0;
	int col_ptr_counter = 0;
	int row_ptr_counter = 0;
	int nnz_values_in_row = 0;

	/*read the entire file*/
	while (!in_file.eof())
	{
		getline(in_file, str_entire_line);

		std::istringstream entire_line(str_entire_line);

		while (entire_line)
		{

			getline(entire_line, str_val, ' ');

			/*description says odd numbers in each line are columns, start from 0, evens are columns odds are nnz values*/
			if (col_ptr_counter % 2 == 0) {
				int_val = atoi(str_val.c_str());
				newCSR.col_ptr_[col_nnz_index_ptr] = int_val - 1;
			}
			else {
				int_val = atof(str_val.c_str());
				newCSR.nnz_val_[col_nnz_index_ptr] = int_val;
				col_nnz_index_ptr++;
				nnz_values_in_row++;
			}

			col_ptr_counter++;

		}

		newCSR.row_ptr_[row_ptr_counter + 1] = nnz_values_in_row + newCSR.row_ptr_[row_ptr_counter];
		col_ptr_counter -= 1;
		nnz_values_in_row = 0;
		row_ptr_counter++;

	}//end while in file

	return new CSR_STRUCTURE(newCSR);
}

/*this function overallocates for the csr nnz_val[] and col_ptr, but the rest of the code doesnt go past it, very wasteful on memory*/
CSR_STRUCTURE* SimilarityCSR(CSR_STRUCTURE &csr_matrix)
{

	
	//double *sim_nnz = new double[csr_matrix.rows_*csr_matrix.columns_];
	//int *sim_col = new int[csr_matrix.rows_*csr_matrix.columns_];
	std::vector<double> sim_nnz;
	std::vector<int> sim_col;
	int *row_counts = new int[csr_matrix.rows_ + 1];
	CSR_STRUCTURE sim_CSR(csr_matrix.rows_, csr_matrix.rows_, csr_matrix.rows_*csr_matrix.rows_);

	for (int i = 0; i < csr_matrix.rows_ + 1; i++) {
		row_counts[i] = 0;
	}

	int col_ctr_index = 0;
	int col_ctr = 0;

	for (int i = 0; i < csr_matrix.rows_; i++)
	{
		for (int j = 0; j < csr_matrix.rows_; j++)
		{
			col_ctr = j;
			int row_i = csr_matrix.row_ptr_[i + 1] - csr_matrix.row_ptr_[i];
			int row_j = csr_matrix.row_ptr_[j + 1] - csr_matrix.row_ptr_[j];
			int ni = 0, nj = 0;
			double cosine = 0, lengthi = 0, lengthj = 0;

			/*go from 0 to however many nonzero values are in each row*/
			while (ni < row_i && nj < row_j)
			{
				int col_i = csr_matrix.row_ptr_[i] + ni;
				int col_j = csr_matrix.row_ptr_[j] + nj;

				/*check to see if the nonzero values are in the same column*/
				/*if the columns are the same, compute the dot product*/
				if (csr_matrix.col_ptr_[col_i] == csr_matrix.col_ptr_[col_j])
				{
					cosine += csr_matrix.nnz_val_[col_i] * csr_matrix.nnz_val_[col_j];
					lengthi += csr_matrix.nnz_val_[col_i] * csr_matrix.nnz_val_[col_i];
					lengthj += csr_matrix.nnz_val_[col_j] * csr_matrix.nnz_val_[col_j];
					ni++;
					nj++;
				}//end if
				 /*if one column [i] is higher than the other [j], then the other [j] needs to catch up*/
				else if (csr_matrix.col_ptr_[col_i] > csr_matrix.col_ptr_[col_j])
				{
					lengthj += csr_matrix.nnz_val_[col_j] * csr_matrix.nnz_val_[col_j];
					nj++;
				}
				/*if neither of the two above, then opposite case of the else if*/
				else
				{
					lengthi += csr_matrix.nnz_val_[col_i] * csr_matrix.nnz_val_[col_i];
					ni++;
				}

			}//end while ni && nj

			/*check to see which row finished first because the while loop exits before calculating the last length at the last index of the nonzero value*/
			while (ni < row_i)
			{
				lengthi += csr_matrix.nnz_val_[ni + csr_matrix.row_ptr_[i]] * csr_matrix.nnz_val_[ni + csr_matrix.row_ptr_[i]];
				ni++;
			}//end while

			while (nj < row_j)
			{
				lengthj += csr_matrix.nnz_val_[nj + csr_matrix.row_ptr_[j]] * csr_matrix.nnz_val_[nj + csr_matrix.row_ptr_[j]];
				nj++;
			}//end while

			/*calculate the cos(i,j)*/
			//only stores the non-zero calculations in the CSR structure and counts non-zero values
			if (lengthi * lengthj)
			{
				cosine /= sqrt(lengthi * lengthj);
				if (cosine > 0) 
				{
					//sim_CSR.increase_nnz_val_();
					//sim_nnz[col_ctr_index] = cosine;
					//sim_col[col_ctr_index] = col_ctr;
					//sim_nnz.push_back(cosine);
					//sim_col.push_back(col_ctr);
					sim_CSR.nnz_val_[col_ctr_index] = cosine;
					sim_CSR.col_ptr_[col_ctr_index] = col_ctr;
					sim_CSR.row_ptr_[i+1]++;
					col_ctr_index++;
				}
				else {
					col_ctr++;
				}
			}
			else {
				cosine = 0;
				col_ctr++;
			}
			//std::cout << "for row i=" << i << "nnz_val_[" << j << "] = " << cosine << ", but actual value at nnz_val_[" << j << "] is " << sim_CSR.nnz_val_[j] << std::endl;
			//std::cout << col_ctr_index << std::endl;
		}//end for j
		sim_CSR.row_ptr_[i+1] += sim_CSR.row_ptr_[i];
		//for (int k = 0; k < sim_CSR.nonzeroes_; k++) {
			//std::cout << "actual value at nnz_val_[" << k << "] is " << sim_CSR.nnz_val_[k] << std::endl;
		//}
	}//end for i
	//std::cout << "out of cosine for loops" << std::endl;
	
	sim_CSR.nnz_val_vec = sim_nnz;
	sim_CSR.col_ptr_vec = sim_col;

	//sim_CSR.setNonZeros(col_ctr_index);
	sim_CSR.nonzeroes_ = col_ctr_index;
	
	//could initialize sim_CSR nonzeros to row*row, which is the max # of entries
	//then resize it, reduce the linear complexity of having to copy the values over
	//for (int i = 0; i < sim_CSR.nonzeroes_; i++) {
	//	sim_CSR.nnz_val_[i] = sim_nnz[i];
	//	sim_CSR.col_ptr_[i] = sim_col[i];
	//}

	/*
	std::cout << "values in sim_CSR row_ptr:" << std::endl;
	for (int i = 0; i < sim_CSR.rows_ + 1; i++) {
		std::cout << sim_CSR.row_ptr_[i] << ", ";
	}
	std::cout << std::endl;
	*/

	/*
	std::cout << "values in sim_CSR nnz_val:";
	for (int i = 0; i < sim_CSR.nonzeroes_; i++) {
		std::cout << sim_CSR.col_ptr_[i] << ", ";
	}
	std::cout << std::endl;
	*/

	return new CSR_STRUCTURE(sim_CSR);
}

double element_in_S(CSR_STRUCTURE csr_matrix, int i, int j) 
{

	return csr_matrix.nnz_val_[csr_matrix.row_ptr_[i] + j];

}

//, double threshold_in
void CalculateSimilarity(CSR_STRUCTURE &csr_matrix, std::ofstream &sim_of, double threshold_in)
{
	/*checks all pairs including the similarity with itself and (j,i), which sim(i,j) == sim(j,i)*/
	for (int i = 0; i < csr_matrix.rows_; i++)
	{
		for (int j = 0; j < csr_matrix.rows_; j++)
		{
			
			int row_i = csr_matrix.row_ptr_[i + 1] - csr_matrix.row_ptr_[i];
			int row_j = csr_matrix.row_ptr_[j + 1] - csr_matrix.row_ptr_[j];
			int ni = 0, nj = 0;
			double cosine = 0, lengthi = 0, lengthj = 0;

			/*go from 0 to however many nonzero values are in each row*/
			while (ni < row_i && nj < row_j)
			{
				int col_i = csr_matrix.row_ptr_[i] + ni;
				int col_j = csr_matrix.row_ptr_[j] + nj;

				/*check to see if the nonzero values are in the same column*/
				/*if the columns are the same, compute the dot product*/
				if (csr_matrix.col_ptr_[col_i] == csr_matrix.col_ptr_[col_j])
				{
					cosine += csr_matrix.nnz_val_[col_i] * csr_matrix.nnz_val_[col_j];
					lengthi += csr_matrix.nnz_val_[col_i] * csr_matrix.nnz_val_[col_i];
					lengthj += csr_matrix.nnz_val_[col_j] * csr_matrix.nnz_val_[col_j];
					ni++;
					nj++;
				}//end if
				/*if one column [i] is higher than the other [j], then the other [j] needs to catch up*/
				else if (csr_matrix.col_ptr_[col_i] > csr_matrix.col_ptr_[col_j])
				{
					lengthj += csr_matrix.nnz_val_[col_j] * csr_matrix.nnz_val_[col_j];
					nj++;
				}
				/*if neither of the two above, then opposite case of the else if*/
				else
				{
					lengthi += csr_matrix.nnz_val_[col_i] * csr_matrix.nnz_val_[col_i];
					ni++;
				}

			}//end while ni && nj

			 
			 /*check to see which row finished first because the while loop exits before calculating the last length at the last index of the nonzero value*/
			
			 while (ni < row_i)
			 {
				 lengthi += csr_matrix.nnz_val_[ni + csr_matrix.row_ptr_[i]] * csr_matrix.nnz_val_[ni + csr_matrix.row_ptr_[i]];
				 ni++;
			 }//end while

			 while (nj < row_j) 
			 {
				 lengthj += csr_matrix.nnz_val_[nj + csr_matrix.row_ptr_[j]] * csr_matrix.nnz_val_[nj + csr_matrix.row_ptr_[j]];
				 nj++;
			 }//end while


			/*calculate the cos(i,j)*/
			if (lengthi * lengthj)
			{
				cosine /= sqrt(lengthi * lengthj);
			}
			else {
				cosine = 0;
			}

			/*output it to the file in i j sim(i,j) format*/
			if (cosine > threshold_in) {
				sim_of << i + 1 << " " << j + 1 << " " << cosine << "\n";
			}

		}//end for j

	}//end for i

	sim_of.close();

}//end cosine sim

/*calculates item similarites and returns the vector*/
std::vector<double>* CalculateVecMatrix(CSR_STRUCTURE &csr_R_in, int u, CSR_STRUCTURE &csr_S_in, int N_rec)
{
	//int j = 0;
	double *uvec = new double[csr_R_in.columns_];
	//double *x = new double[csr_R_in.columns_];
	std::vector<double> x(csr_R_in.columns_);
	int nnz_vals_in_row = csr_R_in.row_ptr_[u] - csr_R_in.row_ptr_[u + 1];

	//initialize both arrays to zeros, biproduct of doing so will only leave us with putting the nonzero values in the array
	for (int i = 0; i < csr_R_in.columns_; i++) 
	{
		//x[i] = 0;
		uvec[i] = 0;
	}

	//throw non-zero values in there respective place in the row
	for (int i = csr_R_in.row_ptr_[u]; i < csr_R_in.row_ptr_[u + 1]; i++) {
		uvec[csr_R_in.col_ptr_[i]] = csr_R_in.nnz_val_[csr_R_in.col_ptr_[i]];
	}

	//for (int i = 0; i < csr_S_in.nonzeroes_; i++) {
	//	std::cout << "i: " << i << " nonzero: " << csr_S_in.nnz_val_[i] << ", " << std::endl;
	//}

	for (int i = 0; i < csr_S_in.rows_; i++) {
		//std::cout << "for row " << i << " in S: ";
		for (int k = csr_S_in.row_ptr_[i]; k < csr_S_in.row_ptr_[i + 1]; k++) {
			int j = csr_S_in.col_ptr_[k];
			double val = csr_S_in.nnz_val_[k];
			//double num = 0;
			//num += val * uvec[i];
			//std::cout << val << ", nonzero val at k: " << csr_S_in.nnz_val_[k] << std::endl;

			//x[j] += val * uvec[i];

			x.at(j) += val * uvec[i];
			
		}
		//std::cout << std::endl;

	}

	//get rid of purchased items for the user
	int nnz_ctr = 0;
	while (nnz_ctr < csr_R_in.row_ptr_[u + 1] - csr_R_in.row_ptr_[u]) 
	{
		if (csr_R_in.nnz_val_[nnz_ctr] > 0) {
			x[csr_R_in.col_ptr_[nnz_ctr]] = 0;
			nnz_ctr++;
		}
		else {
			nnz_ctr++;
		}
	}
	
	//for (std::vector<double>::iterator it = x.begin(); it != x.end(); it++)
	//{
		//if (x[i] > 0) {

	//	std::cout << *it << ", ";

		//}
	//}

	//std::vector<double> copy_x(x);

	//std::sort(x.begin(), x.end());

	//std::cout << "for user " << u << ": ";
	//grab the values for the top N and find their index, which is what item they are

	//std::vector<double>::iterator it = x.end();
	//std::advance(it, csr_R_in.columns_ - N_rec);


	//get top N items
	/*
	int n_ctr = 0;
	double *item_sim_val = new double[N_rec];
	//std::cout << "x.at(csr_R_columns)" << x.at(1681) << std::endl;
	while (n_ctr < N_rec) {
		item_sim_val[n_ctr] = x.at(csr_R_in.columns_ - n_ctr - 1);
		n_ctr++;
	}
	*/
	//for (int i = 0; i < N_rec; i++) {
	//	std::cout << item_sim_val[i] << ", ";
	//}

	//find their item number
	/*
	int *item = new int[N_rec];
	//go through the each val in item_sim array and find its column number ("the item")
	for (int j = 0; j < N_rec; j++) {
		for (int i = 0; i < csr_R_in.columns_; i++) {
			if (item_sim_val[j] == copy_x[i]) {
				item[j] = i+1;
			}
		}
	}
	*/
	//std::cout << std::endl;

	//for (int i = 0; i < N_rec; i++) {
	//	std::cout << item[i] << ", ";
	//}

	return new std::vector<double>(x);

}

//takes in copy of vec_matrix
double* TopNItemScores(std::vector<double> vec_matrix, CSR_STRUCTURE &csr_R_in, int N_rec)
{
	//std::vector<double> copy_vec_matrix(vecmatrix);
	std::sort(vec_matrix.begin(), vec_matrix.end());
	int n_ctr = 0;
	double *item_sim_val = new double[N_rec];

	while (n_ctr < N_rec)
	{
		item_sim_val[n_ctr] = vec_matrix.at(csr_R_in.columns_ - n_ctr - 1);
		n_ctr++;
	}

	return item_sim_val;
}

int* TopNItems(std::vector<double> vec_matrix, CSR_STRUCTURE &csr_R_in, double *item_scores, int N_rec)
{
	int *item = new int[N_rec];
	for (int i = 0; i < N_rec; i++) {
		for (int j = 0; j < csr_R_in.columns_; j++) {
			if (item_scores[i] == vec_matrix[j]) {
				item[i] = j + 1;
			}
		}
	}

	return item;
}

void OutputRecommendations(double *score_matrix, int *item_matrix, std::ofstream &rec_of, int N_rec)
{
	for (int i = 0; i < N_rec; i++) 
	{
		rec_of << item_matrix[i] << " " << score_matrix[i] << " ";
	}
	rec_of << "\n";
}


int main(int argc, char * argv[])
{
	int *top_N;
	double *top_N_scores;
	std::vector<double> vec_matrix_calc;
	int hit_ctr = 0;
	double hit_rate = 0;
	double ar_hit_rate = 0;
	double rank_i = 0;
	clock_t t;
	/*makes sure there are exactly 5 arguments and the threshold*/
	if (argc == 5 && atoi(argv[3]) > 0)
	{
		std::ifstream train_file(argv[1]);
		std::ifstream test_file(argv[2]);

		int N = atoi(argv[3]);
		std::ofstream rec_outfile(argv[4], std::ofstream::out);

		//std::ofstream test_outfile(argv[4], std::ofstream::out);
		//std::ofstream cos_sim_outfile(argv[3], std::ofstream::out);

		//gets reference of new csr_structure created by the functions
		CSR_STRUCTURE csr_M = *InputMatrixToCSR(train_file);
		CSR_STRUCTURE csr_T = *InputMatrixToCSR(test_file);

		t = clock();

		CSR_STRUCTURE M_transpose = *CreateTransposeMatrix(csr_M);
		CSR_STRUCTURE csr_S = *SimilarityCSR(M_transpose);


		//for (int i = 0; i < csr_S.nonzeroes_; i++) {
		//	std::cout << "actual value at nnz_val_[" << i << "] is " << csr_S.nnz_val_[i] << " actual value at col_ptr_[" << i << "] is " << csr_S.col_ptr_[i] << std::endl;
		//}

		//OutputMatrix(csr_T, S_outfile);
		//
		for (int i = 0; i < csr_M.rows_; i++) {
			vec_matrix_calc = *CalculateVecMatrix(csr_M, i, csr_S, N);
			top_N_scores = TopNItemScores(vec_matrix_calc, csr_M, N);
			top_N = TopNItems(vec_matrix_calc, csr_M, top_N_scores, N);
			OutputRecommendations(top_N_scores, top_N, rec_outfile, N);
			
			//std::cout << "user " << i+1 << ": ";
			//std::cout << csr_T.row_ptr_[i] << ", " << csr_T.row_ptr_[i + 1];
			for (int j = 0; j < N; j++) {
				for (int k = csr_T.row_ptr_[i]; k <= csr_T.row_ptr_[i+1]; k++) {
				
					//std::cout << top_N[j] << ", ";
				
					if (top_N[j] == csr_T.col_ptr_[k] + 1) {
						hit_ctr++;
						//std::cout << "hit for user " << i << " item " << top_N[j] << std::endl;
						rank_i += double(j);
					}
				}
			}
			//std::cout << "rank" << rank_i << std::endl;

			//for (int k = csr_T.row_ptr_[i]; k < csr_T.row_ptr_[i + 1]; k++) {
			//	std::cout << csr_T.col_ptr_[k] << ", ";
			//}
			//std::cout << std::endl;
			
		}
		//output_matrix(M_transpose, test_outfile);
		
		hit_rate = hit_ctr / double(csr_T.rows_);
		ar_hit_rate = (1 / double(csr_T.rows_)) * rank_i;

		t = clock() - t;

		std::cout << "HR = " << hit_rate << std::endl;
		std::cout << "ARHR = " << ar_hit_rate << std::endl;
		std::cout << (((float)t) / CLOCKS_PER_SEC) * 1000 << " ms" << std::endl;
		//output_matrix(csr_T, test_outfile);
		//calculate_similarity(M_transpose, test_outfile, 0);

	}//end if
	else
	{
		std::cout << "This program takes 4 arguments" << std::endl;
	}

}//end main
