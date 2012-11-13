
#include <algorithm>

#include "sparse_mat.hpp"
#include "linalg.hpp"


SparseMat::SparseMat(unsigned int nrow, unsigned int ncol,
                     const unsigned int* rowlens,
                     SparseMatEntryStream& entries)
    : nrow(nrow)
    , ncol(ncol)
{
    /* We want the start of each row to be 16-byte aligned. */
    size_t data_size_needed = 0, indexes_size_needed = 0;
    size_t s, r;
    for (unsigned int i = 0; i < nrow; ++i) {
        s = rowlens[i] * sizeof(float);
        r = s % 16;
        data_size_needed += s + (r == 0 ? 0 : (16 - r));

        s = rowlens[i] * sizeof(unsigned int);
        r = s % 16;
        indexes_size_needed += s + (r == 0 ? 0 : (16 - r));
    }

    data = reinterpret_cast<float*>(aalloc(data_size_needed));
    indexes = reinterpret_cast<unsigned int*>(aalloc(indexes_size_needed));

    /* Set row pointers. */
    unsigned int dataoff = 0, indexoff = 0;
    for (unsigned int i = 0; i < nrow; ++i) {
        if (dataoff % 16 != 0) {
            dataoff += 16 - (dataoff % 16);
        }

        data_rows[i] = data + dataoff / sizeof(float);
        dataoff += rowlens[i] * sizeof(float);

        if (indexoff % 16 != 0) {
            indexoff += 16 - (indexoff % 16);
        }

        index_rows[i] = indexes + indexoff / sizeof(unsigned int);
        indexoff += rowlens[i] * sizeof(unsigned int);
    }

    this->rowlens = new unsigned int [nrow];
    std::copy(rowlens, rowlens + nrow, this->rowlens);

    unsigned int prev_row = UINT_MAX;
    unsigned int coloff = 0;
    while (!entries.finished()) {
        if (entries.row() != prev_row) {
            coloff = 0;
            prev_row = entries.row();
        }

        index_rows[entries.row()][coloff] = entries.col();
        data_rows[entries.row()][coloff] = entries.value();
        ++coloff;

        entries.next();
    }
}


SparseMat::~SparseMat()
{
    afree(data);
    afree(indexes);
    delete rowlens;
}


float* SparseMat::rowdata(unsigned int i) { return data_rows[i]; }
const float* SparseMat::rowdata(unsigned int i) const { return data_rows[i]; }


unsigned int* SparseMat::rowindexes(unsigned int i) { return index_rows[i]; }
const unsigned int* SparseMat::rowindexes(unsigned int i) const { return index_rows[i]; }


