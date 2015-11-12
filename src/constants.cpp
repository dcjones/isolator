
#include "constants.hpp"

unsigned int         constants::num_threads                  = 1;
unsigned int         constants::num_opt_rounds               = 30;
size_t               constants::max_estimate_queue_size      = 500;
size_t               constants::seqbias_num_collected_reads  = 500000;
size_t               constants::seqbias_num_reads            = 50000;
size_t               constants::seqbias_left_pos             = 15;
size_t               constants::seqbias_right_pos            = 15;
pos_t                constants::gcbias_min_seq_len           = 500;
size_t               constants::gcbias_num_bins              = 15;
float                constants::gcbias_bins[]                = {0.29, 0.32, 0.35, 0.37, 0.39,
                                                                0.41, 0.43, 0.45, 0.48, 0.51,
                                                                0.53, 0.57, 0.60, 0.64, 1.00};
double               constants::gcbias_max_bias              = 100.0;
pos_t                constants::tpbias_min_tlen              = 5000;
pos_t                constants::tpbias_max_tlen              = 10000;
size_t               constants::tpbias_max_transcripts       = 200000;
constants::libtype_t constants::libtype                      = constants::LIBTYPE_FR;
pos_t                constants::max_frag_len                 = 1000;
float                constants::min_frag_len_pr              = 1e-5;
float                constants::transcript_len_min_frag_pr   = 1e-3;
size_t               constants::frag_len_min_pe_reads        = 10000;
double               constants::frag_len_mu                  = 200.0;
double               constants::frag_len_sd                  = 100.0;
float                constants::min_frag_weight              = 1e-6;
float                constants::zero_eps                     = 1e-9;
pos_t                constants::alt_exon_flank_length        = 200;
pos_t                constants::transcript_5p_extension      = 0;
pos_t                constants::transcript_3p_extension      = 0;
float                constants::min_transcript_weight        = 25.0;
unsigned int         constants::sampler_component_block_size = 1;
double               constants::sample_scaling_quantile      = 0.9;
size_t               constants::sample_scaling_truncation    = 100000;
float                constants::frag_prob_epsilon            = 1e-14;

// These numbers are scaled by the number of transcripts
double               constants::analyze_experiment_mean0     = 2e-4;
float                constants::min_expr                     = 2e-7;

double               constants::analyze_experiment_shape0    = 0.005;
double               constants::analyze_experiment_shape     = 1.0;
double               constants::analyze_experiment_splice_mu0    = 0.5;
double               constants::analyze_experiment_splice_sigma0 = 1.0;
double               constants::analyze_experiment_splice_nu     = 100.0;
double               constants::analyze_min_splice_sigma         = 0.03;
pos_t                constants::max_intron_overlap           = 2;
double               constants::aligned_mismatch_pr          = 0.01;
double               constants::misaligned_mismatch_pr       = 0.05;
double               constants::misalign_prior               = 0.05;
pos_t                constants::fragbias_endlen              = 100;
float                constants::fragbias_scale               = 3e5;
double               constants::min_align_pr                 = 0.1;
size_t               constants::max_newton_steps             = 10;
