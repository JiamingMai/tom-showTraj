package GmmModel;

/**
 * Created by Tom.fu on 22/7/2015.
 */
public class VlGMM_double {

    static double VL_GMM_MIN_VARIANCE = 1e-6;
    static double VL_GMM_MIN_POSTERIOR = 1e-2;
    static double VL_GMM_MIN_PRIOR = 1e-6;

    static enum VlGMMInitialization {VlGMMKMeans, VlGMMRand, VlGMMCustom};

    int dimension;
    int numClusters;
    int numData;
    int maxNumIterations;
    int numRepetitions;
    int verbosity;
    double[] means;
    double[] covariances;
    double[] priors;
    double[] posteriors;
    double[] sigmaLowBound;
    VlGMMInitialization initialization;
    double LL ;
    //VlKMeans * kmeansInit;
    //vl_bool kmeansInitIsOwner;

    public VlGMM_double(int dimension, int numComponents){
        this.numClusters = numComponents;
        this.numData = 0;
        this.dimension = dimension;
        this.initialization = VlGMMInitialization.VlGMMRand;
        this.verbosity = 0;
        this.maxNumIterations = 50;
        this.numRepetitions = 1;
        this.posteriors = null;

        this.priors = new double[numComponents];
        this.means = new double[numComponents * dimension];
        this.covariances = new double[numComponents * dimension];
        this.sigmaLowBound = new double[dimension];

        for (int i = 0; i < dimension; i ++){
            this.sigmaLowBound[i] = 1e-4;
        }
    }

}
