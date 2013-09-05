

// Let's try to sketch out how this thing might look.

// The data needs to be something from which I can efficiently compute a smooth
// estimate of abundance. Look at binned kde maybe. The bandwidth will need to
// vary dramatically.

data {
    // number of sequenced samples
    int<lower=1> K;

    // number of conditions
    int<lower=1> C;

    // number of transcripts
    int<lower=1> N;

    // number of transcription groups (typically, transcription start sites)
    int<lower=1> T;

    // samples generated during quantification
    /*matrix[N, K] quantification[M];*/

    // the condition of each sample
    int<lower=1, upper=C> condition[K];

    // the transcription group of each transcripts
    int<lower=1, upper=T> tgroup[N];

    // log tgroup transcription rate
    real ts[K, T];

    // relative abundance of transcript i, in sample j
    matrix<lower=0>[N, K] xs;

    // degrees of freedom for condition tgroup splice rate
    real tgroup_nu;

    // hyper-hyper-parameters (or whatever) for the prior on gene-level
    // variance.
    real tgroup_alpha_alpha;
    real tgroup_beta_alpha;
    real tgroup_alpha_beta;
    real tgroup_beta_beta;
}


parameters  {
    // log-normal parameters for condition-spcefic expression
    matrix<upper=0>[C, T] tgroup_mu;
    matrix<lower=0>[C, T] tgroup_sigma;
    vector<lower=0>[T] tgroup_alpha;
    vector<lower=0>[T] tgroup_beta;
}


model {
    // Transcription
    // -------------

    // Transcript group abundance
    for (i in 1:K) {
        ts[i] ~ student_t(tgroup_nu,
                          tgroup_mu[condition[i]],
                          tgroup_sigma[condition[i]]);
    }

    // Pooled estimate of variances
    for (i in 1:C) {
        tgroup_sigma[i] ~ inv_gamma(tgroup_alpha, tgroup_beta);
    }

    tgroup_alpha ~ inv_gamma(tgroup_alpha_alpha, tgroup_beta_alpha);
    tgroup_beta ~ inv_gamma(tgroup_alpha_beta, tgroup_beta_beta);


    // Splicing
    // --------
    // TODO
}




