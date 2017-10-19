a)
Compile and run commands(assumes everything is in the same directory):

make -f makefile.mak
./pa2 "train input file" "test input file" N "recommendation output file"

*anything between two quotations ("") is a file name and N is an integer > 0. for example:

make -f makefile.mak
./pa2 train1.txt test1.txt 5 recommendations.out


b)
SimilarityCSR calculations the item-item similarities AND stores it in CSR representation
InputMatrixToCSR takes input files with sparse matrix format and stores it in CSR representation
CalculateVecMatrix gives the list of all non purchased items
TopNItemScores sorts the items based off of item scores
TopNItems matches the scores from TopNItemScores to its item "id" or "number"
Recommendations are done in main(), because of the way I designed my program to return all the information needed to recommend
	without generating r(u,*) for all users, so it does it iteratively, as specified in the PA2 description
	
**I've included my own train1.txt test1.txt (if these are the files you test on) with the correct formatting
	(i cant remember if these data sets were formatted weird)**