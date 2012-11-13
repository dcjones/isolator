
#ifndef ISOLATOR_SPARSE_MAT
#define ISOLATOR_SPARSE_MAT

/* An interface for initilazing a sparse matrixs with a stream of sorted
 * entries. */
class SparseMatEntryStream
{
    public:
        virtual bool finished() = 0;
        virtual void next() = 0;
        virtual unsigned int row() = 0;
        virtual unsigned int col() = 0;
        virtual float value() = 0;
};



/* An extremely simple sparse matrix implementation. */
class SparseMat
{
    public:
        /* Create an empty sparse matrix.
         *
         * Args:
         *   nrow: Number of rows.
         *   ncol: Number of columns.
         *   nnz: Number of non-zero elements to expect.
         */
        SparseMat(unsigned int nrow, unsigned int ncol,
                  const unsigned int* rowlens,
                  SparseMatEntryStream& entries);

        ~SparseMat();

        /* A want a stream of sorted entries with which to build this thing. */

        /* Return a pointer directly to the data at row i. */
        float* rowdata(unsigned int i);
        const float* rowdata(unsigned int i) const;

        unsigned int* rowindexes(unsigned int i);
        const unsigned int* rowindexes(unsigned int i) const;


    private:
        unsigned int nrow, ncol;

        /* rowlens[i] is the number of non-zero entries in row i */
        unsigned int* rowlens;

        float* data;
        unsigned int* indexes;

        /* Row pointers into data and indexes, resp. */
        float** data_rows;
        unsigned int** index_rows;
};

#endif

