
#include "constants.hpp"

unsigned int         constants::num_threads              = 1;
size_t               constants::max_estimate_queue_size  = 500;
pos_t                constants::min_estimate_exon_length = 200;
size_t               constants::seqbias_num_reads        = 25000;
size_t               constants::seqbias_left_pos         = 15;
size_t               constants::seqbias_right_pos        = 15;
constants::libtype_t constants::libtype                  = constants::LIBTYPE_FR;
unsigned int         constants::max_alignments           = 10;
pos_t                constants::max_frag_len             = 1000;
float                constants::min_frag_len_pr          = 1e-8;
float                constants::frag_len_dist_smoothing  = 0.1;
size_t               constants::frag_len_min_pe_reads    = 10000;
double               constants::frag_len_mu              = 200.0;
double               constants::frag_len_sd              = 20.0;
float                constants::min_frag_weight          = 1e-8;
float                constants::zero_eps                 = 1e-8;
float                constants::simplex_constraint_eps   = 1e-3;
float                constants::max_post_objf_tolerance  = 1e-2;
float                constants::max_post_x_tolerance     = 1e-4;
unsigned int         constants::max_post_max_eval        = 10; 
float                constants::min_transcript_weight    = 1e-1;

