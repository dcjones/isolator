
#include "constants.hpp"

unsigned int         constants::num_threads             = 1;
size_t               constants::max_estimate_queue_size = 500;
size_t               constants::seqbias_num_reads       = 50000;
size_t               constants::seqbias_left_pos        = 20;
size_t               constants::seqbias_right_pos       = 20;
constants::libtype_t constants::libtype                 = constants::LIBTYPE_FR;
pos_t                constants::max_frag_len            = 1000;

