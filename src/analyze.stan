

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

    // number of samples generator in the quantification phase
    int<lower=1> M;

    // number of transcription groups (typically, transcription start sites)
    int<lower=1> T;

    // samples generated during quantification
    real quantification[M, N, K];

    // kde bandwidth for using quantification samples as a density estimate
    real bandwidth[N, K];

    // the condition of each sample
    int<lower=1, upper=C> cond[K];

    // the transcription group of each transcripts
    int<lower=1, upper=T> tss[N];

    // normalization factor across samples
    real depth[K];
}


parameters  {
    // relative abundance of transcript i, in sample j
    matrix[N, K] xs;

    // log-normal parameters for condition-spcefic expression
    matrix[T, C] mu;
    matrix[T, C] sigma;
}


model {

    // transcription group usage
    real ts[T, K];

    // temporary for kde computation
    real ps[M];

    // Transcription
    // -------------

    // Kernel density estimate of transcript abundance marginal likelihood with
    // a gaussian kernel, which just becomes a uniform mixture of standard normals.
    for (i in 1:N) {
        for (j in 1:K) {
            for (l in 1:M) {
                ps[l] <- normal_log(
                    (xs[i, j] - quantification[l, i, j]) / bandwidth[i, j], 0, 1);
            }
            increment_log_prob(log_sum_exp(ps) - log(M * bandwidth[i, j]));
        }
    }

    // Transcript group abundance
    for (i in 1:N) {
        for (j in 1:K) {
            ts[tss[i], j] <- ts[tss[i], j] + xs[i, j];
        }
    }

    for (i in 1:T) {
        for (j in 1:K) {
            ts[i, j] / depth[j] ~ lognormal(mu[i, cond[j]], sigma[i, cond[j]]);
        }
    }

    // pooled estimate of variances
    // TODO:
    // Let's concentrate on this stuff



    // Splicing
    // --------
    // TODO
}




