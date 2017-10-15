#include "CSR_STRUCTURE.h"

CSR_STRUCTURE::CSR_STRUCTURE(void)
	: nnz_val_(new double [0]), col_ptr_(new int [0]), row_ptr_(new int [0]), rows_(0), columns_(0), nonzeroes_(0)
{
	for (int i = 0; i < rows_; i++) {
		row_ptr_[i] = 0;
	}
	for (int i = 0; i < nonzeroes_; i++)
	{
		col_ptr_[i] = 0;
		nnz_val_[i] = 0;
	}
}

//csr representation pairs nonzero values with their respective column index, which means the size of the array will be how many nonzero values are in the sparse matrix
//csr representation also iterates on how many nonzero values are in each ROW, so the size of row_ptr will be how many rows+1 (the 0th value in the row_ptr[] will be 1)  
CSR_STRUCTURE::CSR_STRUCTURE(int row_in, int col_in, int nnz_in)
	: nnz_val_(new double [nnz_in]), col_ptr_(new int [nnz_in]), row_ptr_(new int [row_in+1]), rows_(row_in), columns_(col_in), nonzeroes_(nnz_in)
{
	for (int i = 0; i < rows_; i++) {
		row_ptr_[i] = 0;
	}
	for (int i = 0; i <	nonzeroes_; i++)
	{
		col_ptr_[i] = 0;
		nnz_val_[i] = 0;
	}
	//std::cout << "rows: " << rows << " columns: " << columns << " nonzeros: " << nonzeroes << std::endl;
}

CSR_STRUCTURE::CSR_STRUCTURE(const CSR_STRUCTURE &csr) 
	: nnz_val_(new double [csr.nonzeroes_]), col_ptr_(new int[csr.nonzeroes_]), row_ptr_(new int [csr.rows_+1]), rows_(csr.rows_), columns_(csr.columns_), nonzeroes_(csr.nonzeroes_)
{
	for (int i = 0; i < csr.nonzeroes_; i++) {
		nnz_val_[i] = csr.nnz_val_[i];
		col_ptr_[i] = csr.col_ptr_[i];
		//std::cout << " in copy constructor i = " << i << std::endl;
	}

	for (int i = 0; i < csr.rows_+1; i++) {
		row_ptr_[i] = csr.row_ptr_[i];
		//std::cout << " in copy constructor i = " << i << std::endl;
	}
}

CSR_STRUCTURE::~CSR_STRUCTURE(void) 
{
	
	delete[] nnz_val_;
	delete[] col_ptr_;
	delete[] row_ptr_;
	
}

void CSR_STRUCTURE::increase_nnz_val_(void)
{
	nonzeroes_ += 1;

	int *new_col_ptr = new int[nonzeroes_];
	double *new_nnz_val = new double[nonzeroes_];

	//std::cout << "col_ptr and nnz_val resize by one: " << nonzeroes_ << std::endl;

	for (int i = 0; i < nonzeroes_; i++) {
		new_nnz_val[i] = 0;
	}

	//std::cout << "resize at " << nonzeroes_ - 1 << ": ";
	for (int i = 0; i < nonzeroes_; i++) {
		new_col_ptr[i] = col_ptr_[i];
		new_nnz_val[i] = nnz_val_[i];

		//std::cout << "i=" << i << ", new_nnz_val[i]=" << new_nnz_val[i] << ", ";
	}
	//std::cout << std::endl;
	delete[] col_ptr_;
	delete[] nnz_val_;

	//col_ptr_ = new_col_ptr;
	//nnz_val_ = new_nnz_val;
	
	nnz_val_ = new double[nonzeroes_];
	col_ptr_ = new int[nonzeroes_];

	for (int i = 0; i < nonzeroes_; i++) {
		nnz_val_[i] = 0;
	}

	for (int i = 0; i < nonzeroes_; i++) {
		col_ptr_[i] = new_col_ptr[i];
		nnz_val_[i] = new_nnz_val[i];
	}
	

	delete[] new_col_ptr;
	delete[] new_nnz_val;

}


