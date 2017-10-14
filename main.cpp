#include "CSR_STRUCTURE.h"
#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstdlib>
#include <cmath>


CSR_STRUCTURE create_transpose_matrix(CSR_STRUCTURE &csr_matrix)
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

	return transpose_matrix;
}

/*IO function that outputs any matrix to an output file*/
void output_matrix(CSR_STRUCTURE &csr_matrix, std::ofstream &transpose_of)
{
	transpose_of << csr_matrix.rows_ << " " << csr_matrix.columns_ << " " << csr_matrix.nonzeroes_ << "\n";

	for (int i = 0; i < csr_matrix.rows_; i++) {
		for (int j = csr_matrix.row_ptr_[i]; j < csr_matrix.row_ptr_[i + 1]; j++)
		{
			transpose_of << csr_matrix.col_ptr_[j] + 1 << " " << csr_matrix.nnz_val_[j] << " ";
		}
		transpose_of << "\n";
	}

	transpose_of.close();
}

//returns 
CSR_STRUCTURE* input_matrix_to_CSR(std::ifstream &in_file)
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



void cosine_simularity(CSR_STRUCTURE &csr_matrix, std::ofstream &sim_of, double threshold_in)
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
			if (nj == row_j && ni < row_i) {
				lengthi += csr_matrix.nnz_val_[ni + csr_matrix.row_ptr_[i]] * csr_matrix.nnz_val_[ni + csr_matrix.row_ptr_[i]];
			}
			else if (ni == row_i && nj < row_j) {
				lengthj += csr_matrix.nnz_val_[nj + csr_matrix.row_ptr_[j]] * csr_matrix.nnz_val_[nj + csr_matrix.row_ptr_[j]];
			}

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
				sim_of << i+1 << " " << j+1 << " " << cosine << "\n";
			}

		}//end for j

	}//end for i

	sim_of.close();

}//end cosine sim




int main(int argc, char * argv[])
{

	int info_val, int_val;
	double threshold = atof(argv[4]);
	
	/*makes sure there are exactly 5 aguments and the threshold is a value between [0,1]*/
	if (argc == 5 && threshold <= 1 && threshold >= 0)
	{
		std::ifstream input_file(argv[1]);
		std::ofstream rec_outfile(argv[2], std::ofstream::out);
		//std::ofstream cos_sim_outfile(argv[3], std::ofstream::out);

		CSR_STRUCTURE *current_matrix = input_matrix_to_CSR(input_file);
		CSR_STRUCTURE trans_matrix = create_transpose_matrix(*current_matrix);

		//output_matrix(trans_matrix, transpose_outfile);
		//cosine_simularity(trans_matrix, cos_sim_outfile, threshold);

	}//end if
	else
	{
		std::cout << "This program takes 4 arguments" << std::endl;
	}

}//end main
