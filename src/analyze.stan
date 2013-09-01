

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
    /*
        A number of problems with this model.

        1. The variable xs does not work. We need to generate tss usage samples.

        2. This thing uses way too much memory and takes way to long to run.

           I think having a gazillion normal_log terms just breaks the bank when
           computing the gradient. Possible solutions.

           a. Try to manually sample over xs
           b. Fit a normal distribution to the samples.

            Honostly, b seems like the best bet even though it disgusts me.
            The whole point of this exercise was to take into account complex
            variance inherent in these models. I guess that is intractable
            though. Fuck.

            How will we do the same for splicing?



        Another possibility is to actually do the quantification across all
        samples at the same time.

        That way we can actually add in a term for the mean condition
        expression.

        That's kind of a nice idea, but I still really want to try to use stan
        for the complex heirarchical stuff. Is it at all feasable to switch back
        and forth between the two?

        I guess so, since parameters are stored in the sample structure.


    */

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




