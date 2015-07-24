package GmmModel;

import sun.nio.cs.ext.MacThai;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.util.ArrayList;
import java.util.List;
import java.util.Objects;

/**
 * Created by Tom.fu on 16/7/2015.
 */
public class General {

    static double VL_GMM_MIN_VARIANCE = 1e-6;
    static double VL_GMM_MIN_POSTERIOR = 1e-2;
    static double VL_GMM_MIN_PRIOR = 1e-6;


    /** @def VL_FISHER_FLAG_SQUARE_ROOT
     ** @brief Use signed squared-root (@ref fisher-normalization).
     **/
    /** @def VL_FISHER_FLAG_NORMALIZED
     ** @brief Gobally normalize the Fisher vector in L2 norm (@ref fisher-normalization).
     **/
    /** @def VL_FISHER_FLAG_IMPROVED
     ** @brief Improved Fisher vector.
     ** This is the same as @c VL_FISHER_FLAG_SQUARE_ROOT|VL_FISHER_FLAG_NORMALIZED.
     **/
    /** @def VL_FISHER_FLAG_FAST
     ** @brief Fast but more approximate calculations (@ref fisher-fast).
     ** Keep only the larges data to cluster assignment (posterior).
     **/
    /**
     * @}
     */
    static int VL_FISHER_FLAG_SQUARE_ROOT = (0x1 << 0);
    static int VL_FISHER_FLAG_NORMALIZED = (0x1 << 1);
    static int VL_FISHER_FLAG_IMPROVED = (VL_FISHER_FLAG_NORMALIZED | VL_FISHER_FLAG_SQUARE_ROOT);
    static int VL_FISHER_FLAG_FAST = (0x1 << 2);

    public static void main(String[] args) {
        String inputFilePath = "C:\\Users\\Tom.fu\\Desktop\\fromPeiYong\\";
        String fileMeans = inputFilePath + "encode_codebook\\codebook256_f\\validation50.score.gmm.means";
        String fileCovs = inputFilePath + "encode_codebook\\codebook256_f\\validation50.score.gmm.covs";
        String filePriors = inputFilePath + "encode_codebook\\codebook256_f\\validation50.score.gmm.priors";
        String modelFile = inputFilePath + "kth_model_row.model";

        String dataFile = inputFilePath + "trajsTest\\person21_boxing_d1_uncomp.txt";

        int numDimension = 288;
        int numCluster = 256;

        int fvLength = 2 * numCluster * numDimension;
        int numClass = 6;

        double[] means = getData(fileMeans, numCluster * numDimension);
        double[] covs = getData(fileCovs, numCluster * numDimension);
        double[] priors = getData(filePriors, numCluster);
        List<double[]> trainingResult = getTrainingResult(modelFile, fvLength);
/*
    ~isempty(strfind(datalist(i).name,'boxing'))		classids(i)=1;
	~isempty(strfind(datalist(i).name,'handclapping'))  classids(i)=2;
	elseif ~isempty(strfind(datalist(i).name,'handwaving'))classids(i)=3;
	elseif ~isempty(strfind(datalist(i).name,'jogging'))  classids(i)=4;
	elseif ~isempty(strfind(datalist(i).name,'running')) classids(i)=5;
	elseif ~isempty(strfind(datalist(i).name,'walking')) classids(i)=6;
 */
        for (int i = 0; i < 6; i++) {
            double dis = innerProduct(fv_data, trainingResult.get(i));
            System.out.println("i: " + i + ", dis: " + dis);
        }
    }

    public static double check(String testFileFolder, List<double[]> trainingResult,
                               double[] means, int numDimension, int numCluster, double[] covs, double[] priors, int flags) {
        double acc = 0.0;
        int totalCnt = 0;
        int accCnt = 0;

        File folder = new File(testFileFolder);
        File[] listOfFiles = folder.listFiles();
        for (int i = 0; i < listOfFiles.length; i++) {
            File f = listOfFiles[i];
            if (f.isFile()) {
                totalCnt++;
                int orgClassID = getClassID(f.getName());
                double[] fv_data = getFvData(f.getName(), means, numDimension, numCluster, covs, priors, flags);
                int testClassID =
            }
        }

        return acc;
    }

    public static double[] getFvData(
            String dataFile, double[] means, int numDimension, int numCluster, double[] covs, double[] priors, int flags){
        double[] data = getTrajDataWithNormalization(dataFile, 327);
        int numData = data.length / numDimension;

        Object[] enc = vl_fisher_encode_double(means, numDimension, numCluster, covs, priors, data, numData, flags);
        return (double[]) enc[0];
    }

    public static int getClassID(String fileName){
        if (fileName.contains("boxing")) {
            return 0;
        } else if (fileName.contains("handclapping")) {
            return 1;
        } else if (fileName.contains("handwaving")) {
            return 2;
        } else if (fileName.contains("jogging")) {
            return 3;
        }else if (fileName.contains("running")) {
            return 4;
        }else if (fileName.contains("walking")) {
            return 5;
        }else return -1;
    }

    public static double innerProduct(double[] a, double[] b) {
        if (a.length != b.length) {
            System.out.println("warning, a.length != b.length, " + a.length + ", " + b.length);
            return -1;
        }
        double retVal = 0.0;
        for (int i = 0; i < a.length; i++) {
            retVal += a[i] * b[i];
        }
        return retVal;
    }

    public static double[] getData(String fileName, int checkLength) {
        List<Double> dataList = new ArrayList<>();
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String[] rdLineSplit = myfile.readLine().split(" ");
            for (int i = 0; i < rdLineSplit.length; i++) {
                if (!rdLineSplit[i].trim().isEmpty()) {
                    dataList.add(Double.parseDouble(rdLineSplit[i].trim()));
                }
            }
            myfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (checkLength > 0 && dataList.size() != checkLength) {
            System.out.println("warning, dataList.size() != checkLength, " + dataList.size() + ", " + checkLength);
        }

        double[] retVal = new double[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).doubleValue();
        }

        return retVal;
    }

    public static double[] getTrajData(String fileName, int checkLength) {
        List<Double> dataList = new ArrayList<>();
        int lineCount = 0;
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
                lineCount++;
                List<Double> dList = new ArrayList<>();
                String[] rdLineSplit = rdLine.split(" ");
                for (int i = 0; i < rdLineSplit.length; i++) {
                    if (!rdLineSplit[i].trim().isEmpty()) {
                        dList.add(Double.parseDouble(rdLineSplit[i].trim()));
                    }
                }

                if (checkLength > 0 &&  dList.size() != checkLength) {
                    System.out.println("warning, dList.size() != checkLength, " + dList.size() + ", " + checkLength);
                }

                ///only remain the
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        dataList.add(dList.get(i));
                    }
                }
            }
            myfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("LineCount: " + lineCount + ", total data size: " + dataList.size());
        double[] retVal = new double[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).doubleValue();
        }
        return retVal;
    }

    public static double[] getTrajDataWithNormalization(String fileName, int checkLength) {
        List<Double> dataList = new ArrayList<>();
        int lineCount = 0;
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
                lineCount++;
                List<Double> dList = new ArrayList<>();
                String[] rdLineSplit = rdLine.split(" ");
                for (int i = 0; i < rdLineSplit.length; i++) {
                    if (!rdLineSplit[i].trim().isEmpty()) {
                        dList.add(Double.parseDouble(rdLineSplit[i].trim()));
                    }
                }

                if (checkLength > 0 && dList.size() != checkLength) {
                    System.out.println("warning, dList.size() != checkLength, " + dList.size() + ", " + checkLength);
                }

                double sum = 0.0;
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        double p = Math.sqrt(Math.abs(dList.get(i)));
                        p = dList.get(i) > 0.0 ? p : -p;
                        sum += p * p;
                        dList.set(i, p);
                    }
                }

                double base = Math.sqrt(sum);
                ///only remain the
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        if (base > 0.0) {
                            dataList.add(dList.get(i) / base);
                        } else {
                            dataList.add(dList.get(i));
                        }
                    }
                }
            }
            myfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        System.out.println("LineCount: " + lineCount + ", total data size: " + dataList.size());
        double[] retVal = new double[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).doubleValue();
        }
        return retVal;
    }

    public static List<double[]> getTrainingResult(String fileName, int checkLength) {
        List<double[]> retVal = new ArrayList<>();
        int lineCount = 0;
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
                lineCount++;
                List<Double> dList = new ArrayList<>();
                String[] rdLineSplit = rdLine.split(" ");
                for (int i = 0; i < rdLineSplit.length; i++) {
                    if (!rdLineSplit[i].trim().isEmpty()) {
                        dList.add(Double.parseDouble(rdLineSplit[i].trim()));
                    }
                }

                if (checkLength > 0 && dList.size() != checkLength) {
                    System.out.println("warning, dList.size() != checkLength, " + dList.size() + ", " + checkLength);
                }

                double[] v = new double[dList.size()];
                for (int i = 0; i < dList.size(); i++) {
                    v[i] = dList.get(i);
                }
                retVal.add(v);
            }
            myfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        return retVal;
    }

    public static double VlDistanceMahalanobis(
            List<Double> x, List<Double> mu, List<Double> s) {
        double acc = 0.0;
        for (int i = 0; i < x.size(); i++) {
            double diff = x.get(i) - mu.get(i);
            acc += (diff * diff * s.get(i));
        }
        return acc;
    }

    public static double VlDistanceMahalanobis(
            double[] x, double[] mu, double[] s) {
        double acc = 0.0;
        for (int i = 0; i < x.length; i++) {
            double diff = x[i] - mu[i];
            acc += (diff * diff * s[i]);
        }
        return acc;
    }

    ///Object[0] is the output of posteriors, List<Double>
    ///Object[1] is the log-likelihood
    ///

    /**
     * @fn vl_get_gmm_data_posterior_f(float*, vl_size, vl_size, float const*, float const*, vl_size, float const*, float const*)
     * * @brief Get Gaussian modes posterior probabilities
     * *  posteriors posterior probabilities (output)/
     * *  numClusters number of modes in the GMM model.
     * *  numData number of data elements.
     * *  priors prior mode probabilities of the GMM model.
     * *  means means of the GMM model.
     * *  dimension data dimension.
     * *  covariances diagonal covariances of the GMM model.
     * *  data data.
     * *  data log-likelihood.
     * *
     * * This is a helper function that does not require a ::VlGMM_double object
     * * instance to operate.
     **/
    public static Object[] vl_get_gmm_data_posteriors_double(
            int numClusters, int numData, double[] priors,
            double[] means, int dimension, double[] covariances, double[] data) {

        double LL = 0;
        double halfDimLog2Pi = (dimension / 2.0) * Math.log(2.0 * Math.PI);

        //TODO: check the size of this matrix
        double[] posteriors = new double[numClusters * numData];

        double[] logCovariances = new double[numClusters];
        double[] invCovariances = new double[numClusters * dimension];
        double[] logWeights = new double[numClusters];

        for (int i_cl = 0; i_cl < numClusters; i_cl++) {
            double logSigma = 0.0;
            if (priors[i_cl] < VL_GMM_MIN_PRIOR) {
                logWeights[i_cl] = Double.NEGATIVE_INFINITY;
            } else {
                logWeights[i_cl] = Math.log(priors[i_cl]);
            }

            for (int dim = 0; dim < dimension; dim++) {
                logSigma += Math.log(covariances[i_cl * dimension + dim]);
                invCovariances[i_cl * dimension + dim] = 1.0 / covariances[i_cl * dimension + dim];
            }

            logCovariances[i_cl] = logSigma;
        }

        for (int i_d = 0; i_d < numData; i_d++) {
            double clusterPosteriorsSum = 0;
            double maxPosterior = Double.NEGATIVE_INFINITY;

            for (int i_cl = 0; i_cl < numClusters; i_cl++) {
                double[] x = new double[dimension];
                double[] mu = new double[dimension];
                double[] s = new double[dimension];
                for (int k = 0; k < dimension; k++) {
                    x[k] = data[i_d * dimension + k];
                    mu[k] = means[i_cl * dimension + k];
                    s[k] = invCovariances[i_cl * dimension + k];
                }

                double disMeasure = VlDistanceMahalanobis(x, mu, s);
                double p = logWeights[i_cl] - halfDimLog2Pi - 0.5 * logCovariances[i_cl] - 0.5 * disMeasure;
                posteriors[i_cl + i_d * numClusters] = p;
                if (p > maxPosterior) {
                    maxPosterior = p;
                }
            }

            for (int i_cl = 0; i_cl < numClusters; i_cl++) {
                double p = posteriors[i_cl + i_d * numClusters];
                p = Math.exp(p - maxPosterior);
                posteriors[i_cl + i_d * numClusters] = p;
                clusterPosteriorsSum += p;
            }

            LL += (Math.log(clusterPosteriorsSum) + maxPosterior);

            for (int i_cl = 0; i_cl < numClusters; i_cl++) {
                posteriors[i_cl + i_d * numClusters] /= clusterPosteriorsSum;
            }
        }

        return new Object[]{posteriors, LL};
    }


    public static Object[] vl_fisher_encode_double(
            double[] means, int dimension, int numClusters,
            double[] covariances, double[] priors, double[] data, int numData, int flags) {

        int numTerms = 0;

        if (numClusters < 1) {
            System.out.println("vl_fisher_encode_double, warning:  numClusters < 1 !");
        }
        if (dimension < 1) {
            System.out.println("vl_fisher_encode_double, warning:  dimension < 1 !");
        }

        //double[] posteriors = new double[numClusters * numData];
        double[] sqrtInvSigma = new double[dimension * numClusters];

        double[] enc = new double[2 * dimension * numClusters];
        for (int i = 0; i < enc.length; i++) {
            enc[i] = 0.0;
        }

        for (int i_cl = 0; i_cl < numClusters; i_cl++) {
            for (int dim = 0; dim < dimension; dim++) {
                sqrtInvSigma[i_cl * dimension + dim] = Math.sqrt(1.0 / covariances[i_cl * dimension + dim]);
            }
        }

        Object[] gmmPosteriors = vl_get_gmm_data_posteriors_double(
                numClusters, numData, priors, means, dimension, covariances, data);

        double[] posteriors = (double[]) gmmPosteriors[0];
        //double LL = (double)gmmPosteriors[1];

        // sparsify posterior assignments with the FAST option
        if ((flags & VL_FISHER_FLAG_FAST) > 0) {
            System.out.println("vl_fisher_encode_double, enters: flags & VL_FISHER_FLAG_FAST) > 0");
            for (int i_d = 0; i_d < numData; i_d++) {
                // find largest posterior assignment for datum i_d
                int best = 0;
                double bestValue = posteriors[i_d * numClusters];
                for (int i_cl = 0; i_cl < numClusters; i_cl++) {
                    double p = posteriors[i_cl + i_d * numClusters];
                    if (p > bestValue) {
                        bestValue = p;
                        best = i_cl;
                    }
                }
                // make all posterior assignments zero but the best one
                for (int i_cl = 0; i_cl < numClusters; i_cl++) {
                    posteriors[i_cl + i_d * numClusters] = (i_cl == best ? 1.0 : 0.0);
                }
            }
        }

        for (int i_cl = 0; i_cl < numClusters; i_cl++) {
            double uprefix = 0.0;
            double vprefix = 0.0;

             /* If the GMM component is degenerate and has a null prior, then it
            must have null posterior as well. Hence it is safe to skip it.  In
            practice, we skip over it even if the prior is very small; if by
            any chance a feature is assigned to such a mode, then its weight
            would be very high due to the division by priors[i_cl] below.
            */
            if (priors[i_cl] < 1e-6) {
                continue;
            }

            for (int i_d = 0; i_d < numData; i_d++) {
                double p = posteriors[i_cl + i_d * numClusters];
                if (p < 1e-6) {
                    continue;
                }

                numTerms++;
                for (int dim = 0; dim < dimension; dim++) {
                    double diff = data[i_d * dimension + dim] - means[i_cl * dimension + dim];
                    diff *= sqrtInvSigma[i_cl * dimension + dim];
                    enc[i_cl * dimension + dim] += (p * diff);
                    enc[i_cl * dimension + numClusters * dimension + dim] += (p * (diff * diff - 1));
                }
            }

            if (numData > 0) {
                uprefix = 1.0 / ((double) numData * Math.sqrt(priors[i_cl]));
                vprefix = 1.0 / ((double) numData * Math.sqrt(2.0 * priors[i_cl]));
                for (int dim = 0; dim < dimension; dim++) {
                    enc[i_cl * dimension + dim] *= uprefix;
                    enc[i_cl * dimension + numClusters * dimension + dim] *= vprefix;
                }
            }
        }

        if ((flags & VL_FISHER_FLAG_SQUARE_ROOT) > 0) {
            System.out.println("vl_fisher_encode_double, enters: flags & VL_FISHER_FLAG_SQUARE_ROOT) > 0");
            for (int dim = 0; dim < 2 * dimension * numClusters; dim++) {
                double z = enc[dim];
                if (z >= 0) {
                    enc[dim] = Math.sqrt(z);
                } else {
                    enc[dim] = -Math.sqrt(-z);
                }
            }
        }

        if ((flags & VL_FISHER_FLAG_NORMALIZED) > 0) {
            System.out.println("vl_fisher_encode_double, enters: flags & VL_FISHER_FLAG_NORMALIZED) > 0");
            double n = 0.0;
            for (int dim = 0; dim < 2 * dimension * numClusters; dim++) {
                double z = enc[dim];
                n += (z * z);
            }
            n = Math.sqrt(n);
            n = Math.max(n, 1e-12);

            for (int dim = 0; dim < 2 * dimension * numClusters; dim++) {
                enc[dim] /= n;
            }
        }

        return new Object[]{enc, numTerms};
    }


}























