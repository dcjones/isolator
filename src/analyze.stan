

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
    matrix[N, K] quantification[M];

    // kde bandwidth for using quantification samples as a density estimate
    matrix[N, K]  bandwidth;

    // the condition of each sample
    int<lower=1, upper=C> condition[K];

    // the transcription group of each transcripts
    int<lower=1, upper=T> tss[N];

    // normalization factor across samples
    real depth[K];
}


parameters  {
    // relative abundance of transcript i, in sample j
    matrix<lower=0>[N, K] xs;

    // log-normal parameters for condition-spcefic expression
    matrix<upper=0>[T, C] mu;
    matrix<lower=0>[T, C] sigma;
}


model {

    // transcription group usage
    real ts[T, K];

    // temporary for kde computation
    real p;

    // Transcription
    // -------------

    // Kernel density estimate of transcript abundance marginal likelihood with
    // a gaussian kernel, which just becomes a uniform mixture of standard normals.
    for (l in 1:M) {
        p <- normal_log(to_vector((xs - quantification[l]) ./ bandwidth), 0, 1);
        increment_log_prob(p);
    }
    increment_log_prob(sum(log(M * bandwidth)));

    /*
     *for (i in 1:N) {
     *    for (j in 1:K) {
     *        for (l in 1:M) {
     *            ps[l] <- normal_log(
     *                (xs[i, j] - quantification[l, i, j]) / bandwidth[i, j], 0, 1);
     *        }
     *        increment_log_prob(log_sum_exp(ps) - log(M * bandwidth[i, j]));
     *    }
     *}
     */

    // Transcript group abundance
    /*for (i in 1:N) {*/
        /*for (j in 1:K) {*/
            /*ts[tss[i], j] <- ts[tss[i], j] + xs[i, j];*/
        /*}*/
    /*}*/

    /*for (i in 1:T) {*/
        /*for (j in 1:K) {*/
            /*ts[i, j] / depth[j] ~ lognormal(mu[i, condition[j]], sigma[i, condition[j]]);*/

            /*// log determinant of the jacobian for the depth transformation*/
            /*increment_log_prob(-K * log(depth[j]));*/
        /*}*/
    /*}*/

    // pooled estimate of variances
    // TODO:
    // Let's concentrate on this stuff



    // Splicing
    // --------
    // TODO
}




