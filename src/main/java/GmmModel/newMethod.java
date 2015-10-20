package GmmModel;

import com.google.common.collect.HashMultiset;
import redis.clients.jedis.Jedis;

import javax.swing.text.html.parser.Entity;
import java.io.*;
import java.util.*;

/**
 * Created by Tom.fu on 19/8/2015.
 */
public class newMethod {

    static int VL_FISHER_FLAG_SQUARE_ROOT = (0x1 << 0);
    static int VL_FISHER_FLAG_NORMALIZED = (0x1 << 1);
    static int VL_FISHER_FLAG_IMPROVED = (VL_FISHER_FLAG_NORMALIZED | VL_FISHER_FLAG_SQUARE_ROOT);
    static int VL_FISHER_FLAG_FAST = (0x1 << 2);

    static double VL_GMM_MIN_PRIOR_F = 1e-6;

//    public static void mainTest(String[] args) {
//        Map<String, Double> results = new HashMap();
//        results.put("v_Skiing_g21_c05", 0.04);
//        results.put("v_Skiing_g23_c02", 0.03);
//        results.put("v_Biking_g24_c04", 0.01);
//        results.put("v_BoxingPunchingBag_g21_c01", 0.09);
//        results.put("v_Drumming_g24_c04", 0.08);
//        results.put("v_JumpingJack_g22_c02", 0.07);
//
//        results.entrySet().stream().sorted((e1, e2) -> e1.getKey().compareTo(e2.getKey())).filter(e -> e.getValue() > 0.02).forEach(e -> System.out.println(e));
//    }

    public static void main(String[] args) {

        String inputFilePath = "C:\\Users\\Tom.fu\\Desktop\\fromPeiYong\\";

        String demoFolder = "ucf101-Demo";

        String hogPcaFile = inputFilePath + demoFolder + "\\codebook\\hog_pca.mat";
        String mbhxPcaFile = inputFilePath + demoFolder + "\\codebook\\mbhx_pca.mat";
        String mbhyPcaFile = inputFilePath + demoFolder + "\\codebook\\mbhy_pca.mat";
        String hogGmmFile = inputFilePath + demoFolder + "\\codebook\\hog_gmm.codebook";
        String mbhxGmmFile = inputFilePath + demoFolder + "\\codebook\\mbhx_gmm.codebook";
        String mbhyGmmFile = inputFilePath + demoFolder + "\\codebook\\mbhy_gmm.codebook";

//        String modelFile = inputFilePath + demoFolder + "\\KTHDemo2_row.model";
        String modelFile = inputFilePath + demoFolder + "\\ucf101_demo_row.model";
        String classNameFile = inputFilePath + demoFolder + "\\ucf101_subset_labels.txt";

        List<String> classNameDic = getClassNames(classNameFile);
        for (int i = 0; i < classNameDic.size(); i ++){
            System.out.println("i: " + i + ", " + classNameDic.get(i));
        }

//        String dataFileFolder = inputFilePath + "trajsTest-new";
        String dataFileFolder = inputFilePath + "testForUcf101";

        int numDimension = 288;
        int numCluster = 256;

        int fvLength = 2 * numCluster * numDimension / 2;
        int numClass = 6;

        PcaData hogPca = new PcaData(hogPcaFile);
        PcaData mbhxPca = new PcaData(mbhxPcaFile);
        PcaData mbhyPca = new PcaData(mbhyPcaFile);

        GmmData hogGmm = new GmmData(hogGmmFile);
        GmmData mbhxGmm = new GmmData(mbhxGmmFile);
        GmmData mbhyGmm = new GmmData(mbhyGmmFile);

        List<float[]> trainingResult = getTrainingResult_float(modelFile, fvLength);

//        String testFile = inputFilePath + "KTHDemo2\\person25_boxing_d3_uncomp.txt";
//        float[] rawData = getTrajDataFromFileNew_float(testFile, 327);

//        List<float[]> rawInputData = getTrajDataNew_float("192.168.0.30", 6379, "tomQ", 0, -1);

//        float[] rawData = getTrajDataNew_float("192.168.0.30", 6379, "tomQ", 0, -1, 288);
//        int numData = rawData.length / numDimension;
//        System.out.println("rawData.length: " + rawData.length + ", numDimension: " + numDimension + ", numData: " + numData);

//        int classResult = checkNew_float(rawInputData, trainingResult, numDimension, hogPca, mbhxPca, mbhyPca, hogGmm, mbhxGmm, mbhyGmm, true);
//        float[] hogPart = extractPartWithNormalization(rawData, numData, numDimension, 0, 96);
//        float[] mbhxPart = extractPartWithNormalization(rawData, numData, numDimension, 96, 96);
//        float[] mbhyPart = extractPartWithNormalization(rawData, numData, numDimension, 192, 96);
//
//        float[] hogPcaPart = applyPca(hogPart, numData, hogPca.aveDim, hogPca);
//        float[] mbhxPcaPart = applyPca(mbhxPart, numData, mbhxPca.aveDim, mbhxPca);
//        float[] mbhyPcaPart = applyPca(mbhyPart, numData, mbhyPca.aveDim, mbhyPca);
//
//        float[] hogFV = getFVDataNew_float(hogPcaPart, numData, hogPca.eigDim, hogGmm);
//        float[] mbhxFV = getFVDataNew_float(mbhxPcaPart, numData, mbhxPca.eigDim, mbhxGmm);
//        float[] mbhyFV = getFVDataNew_float(mbhyPcaPart, numData, mbhyPca.eigDim, mbhyGmm);
//
//        Object[] classifyRestult = getClassificationResult_float(trainingResult, hogFV, mbhxFV, mbhyFV, true);
//        int maxIndex = (int) classifyRestult[0];
//        float maxSimilarity = (float) classifyRestult[1];
//        System.out.println("numData: " + numData + ", maxIndex: " + maxIndex + ", maxSim: " + maxSimilarity);

        Map<String, Float> results = checkNew_float(dataFileFolder, trainingResult, numDimension, hogPca, mbhxPca, mbhyPca, hogGmm, mbhxGmm, mbhyGmm, classNameDic);

        results.entrySet().stream().sorted((e1, e2) -> e1.getKey().compareTo(e2.getKey())).filter(e -> e.getValue() > 0.6f).forEach(e -> System.out.println(e));

//        System.out.println("results > 0.6f");
//        results.forEach((k, v)-> {if ((float)v > 0.6f) {
//            System.out.println(k + "," +  v);
//        }});

}

    public static int checkNew_float(List<float[]> rawInputData, List<float[]> trainingResult, int numDimension,
                                     PcaData hogPca, PcaData mbhxPca, PcaData mbhyPca, GmmData hogGmm, GmmData mbhxGmm, GmmData mbhyGmm, boolean toPrint) {


//        float[] rawData = new float[numDimension * rawInputData.size()];
//        for (int i = 0; i < rawInputData.size(); i ++){
//            for (int j = 0; j < numDimension; j ++){
//                rawData[i * numDimension + j] = rawInputData.get(i)[j];
//            }
//        }
        int numData = rawInputData.size();
        float[] hogPart = extractPartWithNormalization(rawInputData, numData, numDimension, 0, 96);
        float[] mbhxPart = extractPartWithNormalization(rawInputData, numData, numDimension, 96, 96);
        float[] mbhyPart = extractPartWithNormalization(rawInputData, numData, numDimension, 192, 96);

        float[] hogPcaPart = applyPca(hogPart, numData, hogPca.aveDim, hogPca);
        float[] mbhxPcaPart = applyPca(mbhxPart, numData, mbhxPca.aveDim, mbhxPca);
        float[] mbhyPcaPart = applyPca(mbhyPart, numData, mbhyPca.aveDim, mbhyPca);

        float[] hogFV = getFVDataNew_float(hogPcaPart, numData, hogPca.eigDim, hogGmm);
        float[] mbhxFV = getFVDataNew_float(mbhxPcaPart, numData, mbhxPca.eigDim, mbhxGmm);
        float[] mbhyFV = getFVDataNew_float(mbhyPcaPart, numData, mbhyPca.eigDim, mbhyGmm);

        Object[] classifyRestult = getClassificationResult_float(trainingResult, hogFV, mbhxFV, mbhyFV, toPrint);

        int testClassID = (int) classifyRestult[0];
        float similarity = (float) classifyRestult[1];

        return testClassID;
    }

    public static Map checkNew_float(String testFileFolder, List<float[]> trainingResult, int numDimension,
                                       PcaData hogPca, PcaData mbhxPca, PcaData mbhyPca, GmmData hogGmm, GmmData mbhxGmm, GmmData mbhyGmm, List<String> classNameDictionary) {
        int totalCnt = 0;
        int accCnt = 0;

        Map<String, Float> retVal = new HashMap();

        int[] classTotalCnt = new int[classNameDictionary.size()];
        int[] classAccCnt = new int[classNameDictionary.size()];

        File folder = new File(testFileFolder);
        File[] listOfFiles = folder.listFiles();
        for (int i = 0; i < listOfFiles.length; i++) {
            File f = listOfFiles[i];
            if (f.isFile()) {
                int orgClassID = getClassID(f.getName(), classNameDictionary);
                if (orgClassID < 0) {
                    System.out.println("Warning, in check_double, orgClassID < 0, fileName: " + f.getName());
                    continue;
                }

                totalCnt++;
                classTotalCnt[orgClassID] ++;

                float[] rawData = getTrajDataFromFileNew_float(f.getAbsolutePath(), 327);
                int numData = rawData.length / numDimension;
                float[] hogPart = extractPartWithNormalization(rawData, numData, numDimension, 0, 96);
                float[] mbhxPart = extractPartWithNormalization(rawData, numData, numDimension, 96, 96);
                float[] mbhyPart = extractPartWithNormalization(rawData, numData, numDimension, 192, 96);

                float[] hogPcaPart = applyPca(hogPart, numData, hogPca.aveDim, hogPca);
                float[] mbhxPcaPart = applyPca(mbhxPart, numData, mbhxPca.aveDim, mbhxPca);
                float[] mbhyPcaPart = applyPca(mbhyPart, numData, mbhyPca.aveDim, mbhyPca);

                float[] hogFV = getFVDataNew_float(hogPcaPart, numData, hogPca.eigDim, hogGmm);
                float[] mbhxFV = getFVDataNew_float(mbhxPcaPart, numData, mbhxPca.eigDim, mbhxGmm);
                float[] mbhyFV = getFVDataNew_float(mbhyPcaPart, numData, mbhyPca.eigDim, mbhyGmm);

                Object[] classifyRestult = getClassificationResult_float(trainingResult, hogFV, mbhxFV, mbhyFV, true);

                int testClassID = (int) classifyRestult[0];
                float similarity = (float) classifyRestult[1];
                if (orgClassID == testClassID) {
                    accCnt++;
                    classAccCnt[orgClassID]++;
                    System.out.println("File: " + f.getName() + ", orgClassID: " + orgClassID + ", testID: " + testClassID + ", sim: " + similarity);

                    retVal.put(f.getName(), similarity);

                } else {
                    System.out.println("Wrong classification, File: " + f.getName() + ", orgClassID: " + orgClassID + ", testID: " + testClassID + ", sim: " + similarity);
                }
            }
        }
        double acc = totalCnt == 0 ? 0.0 : (double) accCnt / (double) totalCnt;
        System.out.println("check_float, totalCnt: " + totalCnt + ", accCnt: " + accCnt + ", accuracy: " + acc);
        for (int i = 0; i < classNameDictionary.size(); i ++){
            double classAcc = classTotalCnt[i] == 0 ? 0.0 : (double) classAccCnt[i] / (double) classTotalCnt[i];
            System.out.println("class: " + classNameDictionary.get(i) + ", classCnt: " + classTotalCnt[i] + ", classAccCnt: " + classAccCnt[i] + ", classAcc: " + classAcc);
        }
        return retVal;
    }

    public static void dataCheck(float[] input, int numData, int showStart, int showCnt) {
        int dim = input.length / numData;
        float max = Float.MIN_VALUE;
        float min = Float.MAX_VALUE;

        for (int i = 0; i < input.length; i++) {
            if (input[i] > max) {
                max = input[i];
            }

            if (input[i] < min) {
                min = input[i];
            }
        }

        System.out.println("total length: " + input.length + ", numData: " + numData + ", dim: " + dim + ", max: " + max + ", min: " + min);
        for (int i = 0; i < showCnt; i++) {
            System.out.printf(input[showStart + i] + " ");
        }
        System.out.println();
    }

    public static float[] applyPca(float[] input, int numData, int numDim, PcaData pcaData) {

        if (input.length != numData * numDim) {
            System.out.println("in applyPca, warning input.length != numData * numDim, input.Length: " + input.length
                    + ",numDim: " + numDim + ", numData: " + numData + ", " + numData * numDim);
            return null;
        }

        if (numDim != pcaData.aveDim) {
            System.out.println("in applyPca, warning numDim != pcaData.aveDim, numDim: " + numDim + ", pcaData.aveDim: " + pcaData.aveDim);
            return null;
        }

        float[] retVal = new float[numData * pcaData.eigDim];
        for (int i = 0; i < numData; i++) {
            for (int j = 0; j < numDim; j++) {
                int inputIndex = i * numDim + j;
                input[inputIndex] -= pcaData.aveVec[j];
            }
            for (int k = 0; k < pcaData.eigDim; k++) {
                float tempS = 0;
                for (int j = 0; j < numDim; j++) {
                    int inputIndex = i * numDim + j;
                    int pcaDataIndex = k * numDim + j;
                    tempS += (input[inputIndex] * pcaData.projMat[pcaDataIndex]);
                }
                int outputIndex = i * pcaData.eigDim + k;
                retVal[outputIndex] = tempS;
            }
        }
        return retVal;
    }

    ///output: [288][288]...traceNum...[288]
    public static float[] getTrajDataFromFileNew_float(String fileName, int checkLength) {
        List<Float> dataList = new ArrayList<>();
        int lineCount = 0;
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
                lineCount++;
                List<Float> dList = new ArrayList<>();
                String[] rdLineSplit = rdLine.split(" ");
                for (int i = 0; i < rdLineSplit.length; i++) {
                    if (!rdLineSplit[i].trim().isEmpty()) {
                        dList.add(Float.parseFloat(rdLineSplit[i].trim()));
                    }
                }

                if (checkLength > 0 && dList.size() != checkLength) {
                    System.out.println("warning, dList.size() != checkLength, " + dList.size() + ", " + checkLength);
                }

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
        float[] retVal = new float[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).floatValue();
        }
        return retVal;
    }

    ///input: [288][288]...traceNum...[288]
    ///output: [normalized(96)][normalized(96)]...traceNum...[normalized(96)] for each of hog, mbhx, mhby
    public static float[] extractPartWithNormalization(float[] source, int numData, int numDim, int start, int count) {
        if (source.length != numData * numDim) {
            System.out.println("in extractPart, warning input.length != numData * numDim, input.Length: " + source.length
                    + ",numDim: " + numDim + ", numData: " + numData + ", " + numData * numDim);
            return null;
        }

        float[] retVal = new float[numData * count];

        for (int i = 0; i < numData; i++) {
            float sum = 0;
            for (int j = 0; j < count; j++) {
                int sourceIndex = i * numDim + start + j;
                sum += (source[sourceIndex] * source[sourceIndex]);
            }
            float base = (float) Math.sqrt(sum);
            for (int j = 0; j < count; j++) {
                int sourceIndex = i * numDim + start + j;
                int newIndex = i * count + j;
                retVal[newIndex] = source[sourceIndex] / base;
            }
        }
        return retVal;
    }

    ///input: [288][288]...traceNum...[288]
    ///output: [normalized(96)][normalized(96)]...traceNum...[normalized(96)] for each of hog, mbhx, mhby
    public static float[] extractPartWithNormalization(List<float[]> source, int numData, int numDim, int start, int count) {
        if (source.size() != numData) {
            System.out.println("in extractPart, warning source.size() != numData, source.size: " + source.size()
                    + ",numDim: " + numDim + ", numData: " + numData);
            return null;
        }

        float[] retVal = new float[numData * count];

        for (int i = 0; i < numData; i++) {
            float sum = 0;
            for (int j = 0; j < count; j++) {
                sum += (source.get(i)[start + j] * source.get(i)[start + j]);
            }
            float base = (float) Math.sqrt(sum);
            for (int j = 0; j < count; j++) {
                int newIndex = i * count + j;
                retVal[newIndex] = source.get(i)[start + j] / base;
            }
        }
        return retVal;
    }

    public static float[] getFVDataNew_float(float[] rawData, int numData, int numDim, GmmData gmmData) {
//        Object[] enc = vl_fisher_encode_float(means, numDimension, numCluster, covs, priors, rawData, numData, 0);
//        float[] fvData = (float[]) enc[0];
        if (rawData.length != numData * numDim) {
            System.out.println("in getFVDataNew_float, warning input.length != numData * numDim, rawData.Length: " + rawData.length
                    + ",numDim: " + numDim + ", numData: " + numData + ", " + numData * numDim);
            return null;
        }

        if (numDim != gmmData.numDimension) {
            System.out.println("in getFVDataNew_float, warning numDim != gmmData.numDimension, numDim: "
                    + numDim + ", gmmData.numDimension: " + gmmData.numDimension);
            return null;
        }

        Object[] enc = vl_fisher_encode_float(gmmData.means, gmmData.numDimension, gmmData.numCluster,
                gmmData.covs, gmmData.priors, rawData, numData, 0);

        float[] rawFVData = (float[]) enc[0];
        if (rawFVData.length != gmmData.fvLength) {
            System.out.println("in getFVDataNew_float, warning rawFVData.length != gmmData.fvLength, rawFVData.Length: "
                    + rawFVData.length + ",numDim: " + numDim + ", gmmData.numCluster: " + gmmData.numCluster + ", " + gmmData.fvLength);
            return null;
        }

        int half = rawFVData.length / 2;
        float sum = 0;
        for (int i = 0; i < half; i++) {
            int sqrtPriIndex = i / numDim;
            float q = rawFVData[i] / gmmData.sqrtPriors[sqrtPriIndex];
            float p = (float) Math.sqrt(Math.abs(q));
            p = q > 0 ? p : -p;
            sum += p * p;
            rawFVData[i] = p;
        }
        float base = (float) Math.sqrt(sum);
        for (int i = 0; i < half; i++) {
            rawFVData[i] /= base;
        }

        sum = 0;
        base = 0;
        for (int i = half; i < rawFVData.length; i++) {
            int sqrtPriIndex = (i - half) / numDim;
            float q = rawFVData[i] / gmmData.sqrtPriors[sqrtPriIndex];
            float p = (float) Math.sqrt(Math.abs(q));
            p = q > 0 ? p : -p;
            sum += p * p;
            rawFVData[i] = p;
        }
        base = (float) Math.sqrt(sum);
        for (int i = half; i < rawFVData.length; i++) {
            rawFVData[i] /= base;
        }

        return rawFVData;
    }

    public static float innerProduct(float[] x1, float[] y1, float[] y2, float[] y3) {
        if (x1.length != (y1.length + y2.length + y3.length)) {
            System.out.println("warning, x1.length != (y1.length + y2.length + y3.length), "
                    + x1.length + ", " + y1.length + ", " + y2.length + ", " + y3.length);
            return -1;
        }
        float retVal = 0;
        for (int i = 0; i < y1.length; i++) {
            retVal += x1[i] * y1[i];
        }
        for (int i = 0; i < y2.length; i++) {
            retVal += x1[y1.length + i] * y2[i];
        }
        for (int i = 0; i < y3.length; i++) {
            retVal += x1[y1.length + y2.length + i] * y3[i];
        }
        return retVal;
    }

    public static Object[] getClassificationResult_float(List<float[]> trainingResult, float[] fv1, float[] fv2, float[] fv3, boolean toPrint) {
        int maxIndex = 0;
        float maxSimilarity = Float.NEGATIVE_INFINITY;

        for (int i = 0; i < trainingResult.size(); i++) {
            float sim = innerProduct(trainingResult.get(i), fv1, fv2, fv3);
            if (toPrint) {
                System.out.println("i: " + i + ", sim: " + sim);
            }
            if (sim > maxSimilarity) {
                maxSimilarity = sim;
                maxIndex = i;
            }
        }
        return new Object[]{maxIndex, maxSimilarity};
    }

    public static Object[] vl_fisher_encode_float(
            float[] means, int dimension, int numClusters,
            float[] covariances, float[] priors, float[] data, int numData, int flags) {

        int numTerms = 0;

        if (numClusters < 1) {
            System.out.println("vl_fisher_encode_double, warning:  numClusters < 1 !");
        }
        if (dimension < 1) {
            System.out.println("vl_fisher_encode_double, warning:  dimension < 1 !");
        }

        //double[] posteriors = new double[numClusters * numData];
        float[] sqrtInvSigma = new float[dimension * numClusters];

        float[] enc = new float[2 * dimension * numClusters];
        for (int i = 0; i < enc.length; i++) {
            enc[i] = 0;
        }

        for (int i_cl = 0; i_cl < numClusters; i_cl++) {
            for (int dim = 0; dim < dimension; dim++) {
                sqrtInvSigma[i_cl * dimension + dim] = (float) Math.sqrt(1.0f / covariances[i_cl * dimension + dim]);
            }
        }

        Object[] gmmPosteriors = vl_get_gmm_data_posteriors_float(
                numClusters, numData, priors, means, dimension, covariances, data);

        float[] posteriors = (float[]) gmmPosteriors[0];
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
                    posteriors[i_cl + i_d * numClusters] = (i_cl == best ? 1.0f : 0.0f);
                }
            }
        }

        for (int i_cl = 0; i_cl < numClusters; i_cl++) {
            float uprefix = 0;
            float vprefix = 0;

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
                float p = posteriors[i_cl + i_d * numClusters];
                if (p < 1e-6) {
                    continue;
                }

                numTerms++;
                for (int dim = 0; dim < dimension; dim++) {
                    float diff = data[i_d * dimension + dim] - means[i_cl * dimension + dim];
                    diff *= sqrtInvSigma[i_cl * dimension + dim];
                    enc[i_cl * dimension + dim] += (p * diff);
                    enc[i_cl * dimension + numClusters * dimension + dim] += (p * (diff * diff - 1));
                }
            }

            if (numData > 0) {
                uprefix = 1.0f / ((float) numData * (float) Math.sqrt(priors[i_cl]));
                vprefix = 1.0f / ((float) numData * (float) Math.sqrt(2.0 * priors[i_cl]));
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
                    enc[dim] = (float) Math.sqrt(z);
                } else {
                    enc[dim] = (float) -Math.sqrt(-z);
                }
            }
        }

        if ((flags & VL_FISHER_FLAG_NORMALIZED) > 0) {
            System.out.println("vl_fisher_encode_double, enters: flags & VL_FISHER_FLAG_NORMALIZED) > 0");
            float n = 0;
            for (int dim = 0; dim < 2 * dimension * numClusters; dim++) {
                float z = enc[dim];
                n += (z * z);
            }
            n = (float) Math.sqrt(n);
            n = (float) Math.max(n, 1e-12);

            for (int dim = 0; dim < 2 * dimension * numClusters; dim++) {
                enc[dim] /= n;
            }
        }

        return new Object[]{enc, numTerms};
    }

    public static Object[] vl_get_gmm_data_posteriors_float(
            int numClusters, int numData, float[] priors,
            float[] means, int dimension, float[] covariances, float[] data) {

        float LL = 0;
        float halfDimLog2Pi = (dimension / 2.0f) * (float) Math.log(2.0 * Math.PI);

        //TODO: check_double the size of this matrix
        float[] posteriors = new float[numClusters * numData];

        float[] logCovariances = new float[numClusters];
        float[] invCovariances = new float[numClusters * dimension];
        float[] logWeights = new float[numClusters];

        for (int i_cl = 0; i_cl < numClusters; i_cl++) {
            float logSigma = 0;
            if (priors[i_cl] < VL_GMM_MIN_PRIOR_F) {
                logWeights[i_cl] = Float.NEGATIVE_INFINITY;
            } else {
                logWeights[i_cl] = (float) Math.log(priors[i_cl]);
            }

            for (int dim = 0; dim < dimension; dim++) {
                logSigma += Math.log(covariances[i_cl * dimension + dim]);
                invCovariances[i_cl * dimension + dim] = 1.0f / covariances[i_cl * dimension + dim];
            }
            logCovariances[i_cl] = logSigma;
        }

        for (int i_d = 0; i_d < numData; i_d++) {
            float clusterPosteriorsSum = 0;
            float maxPosterior = Float.NEGATIVE_INFINITY;

            for (int i_cl = 0; i_cl < numClusters; i_cl++) {
                float[] x = new float[dimension];
                float[] mu = new float[dimension];
                float[] s = new float[dimension];
                for (int k = 0; k < dimension; k++) {
                    x[k] = data[i_d * dimension + k];
                    mu[k] = means[i_cl * dimension + k];
                    s[k] = invCovariances[i_cl * dimension + k];
                }

                float disMeasure = VlDistanceMahalanobis(x, mu, s);
                float p = (float) logWeights[i_cl] - halfDimLog2Pi - 0.5f * logCovariances[i_cl] - 0.5f * disMeasure;
                posteriors[i_cl + i_d * numClusters] = p;
                if (p > maxPosterior) {
                    maxPosterior = p;
                }
            }

            for (int i_cl = 0; i_cl < numClusters; i_cl++) {
                float p = posteriors[i_cl + i_d * numClusters];
                p = (float) Math.exp(p - maxPosterior);
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

    public static List<String> getClassNames(String classNameFile){
        List<String> retVal = new ArrayList<>();
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(classNameFile));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
                if (rdLine.length() > 0) {
                    retVal.add(rdLine);
                }
            }
            myfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
        return retVal;
    }

    public static int getClassID(String fileName) {
        if (fileName.contains("boxing")) {
            return 0;
        } else if (fileName.contains("handclapping")) {
            return 1;
        } else if (fileName.contains("handwaving")) {
            return 2;
        } else if (fileName.contains("jogging")) {
            return 3;
        } else if (fileName.contains("running")) {
            return 4;
        } else if (fileName.contains("walking")) {
            return 5;
        } else return -1;
    }

    public static int getClassID(String fileName, List<String> classNameDictionary) {
        for (int i = 0; i < classNameDictionary.size(); i ++){
            if (fileName.contains(classNameDictionary.get(i))){
                return i;
            }
        }
        return -1;
    }

    public static List<float[]> getTrainingResult_float(String fileName, int checkLength) {
        List<float[]> retVal = new ArrayList<>();
        int lineCount = 0;
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
                lineCount++;
                List<Float> dList = new ArrayList<>();
                String[] rdLineSplit = rdLine.split(" ");
                for (int i = 0; i < rdLineSplit.length; i++) {
                    if (!rdLineSplit[i].trim().isEmpty()) {
                        dList.add(Float.parseFloat(rdLineSplit[i].trim()));
                    }
                }

                if (checkLength > 0 && dList.size() != checkLength) {
                    System.out.println("warning, dList.size() != checkLength, " + dList.size() + ", " + checkLength);
                }

                float[] v = new float[dList.size()];
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

    public static float VlDistanceMahalanobis(
            float[] x, float[] mu, float[] s) {
        double acc = 0;
        for (int i = 0; i < x.length; i++) {
            double diff = x[i] - mu[i];
            acc += (diff * diff * s[i]);
        }
        return (float) acc;
    }

    /**
     * For get data from redis
     *
     * @param host
     * @param port
     * @param queue
     * @param startIndex
     * @param endIndex
     * @param checkLength
     * @return
     */
    public static float[] getTrajDataNew_float(String host, int port, String queue, int startIndex, int endIndex, int checkLength) {

        ///String host = "192.168.0.30";
        ///int port = 6379;
        ///String queue = "tomQ";
        byte[] queueName = queue.getBytes();
        Jedis jedis = new Jedis(host, port);
        List<byte[]> getData = jedis.lrange(queueName, startIndex, endIndex);
        if (getData == null) {
            System.out.println("Warning! in getTrajDataWithNormalization_float, getData == null, queueName: " + queue);
        }

        List<Float> dataList = new ArrayList<>();
        int lineCount = 0;
        for (int j = 0; j < getData.size(); j++) {
            List<float[]> traces = readArrays(getData.get(j));

            for (int k = 0; k < traces.size(); k++) {
                float[] trace = traces.get(k);
                if (checkLength > 0 && trace.length != checkLength) {
                    System.out.println("warning, trace.length != checkLength, " + trace.length + ", " + checkLength);
                }
                lineCount++;
                for (int i = 0; i < trace.length; i++) {
                    dataList.add(trace[i]);
                }
            }
        }

        System.out.println("DataSize: " + getData.size() + ", lineCount: " + lineCount + ", total data size: " + dataList.size());
        float[] retVal = new float[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).floatValue();
        }

        return retVal;
    }

    /**
     * For get data from redis
     *
     * @param host
     * @param port
     * @param queue
     * @param startIndex
     * @param endIndex
     * @return
     */
    public static List<float[]> getTrajDataNew_float(String host, int port, String queue, int startIndex, int endIndex) {

        ///String host = "192.168.0.30";
        ///int port = 6379;
        ///String queue = "tomQ";
        byte[] queueName = queue.getBytes();
        Jedis jedis = new Jedis(host, port);
        List<byte[]> getData = jedis.lrange(queueName, startIndex, endIndex);
        if (getData == null) {
            System.out.println("Warning! in getTrajDataWithNormalization_float, getData == null, queueName: " + queue);
        }

        List<float[]> retVal = new ArrayList<>();
        int lineCount = 0;
        for (int j = 0; j < getData.size(); j++) {
            List<float[]> traces = readArrays(getData.get(j));
            retVal.addAll(traces);
        }
        return retVal;
    }

    public static byte[] toBytes(List<float[]> data) {
        ByteArrayOutputStream bout = new ByteArrayOutputStream();
        DataOutputStream dout = new DataOutputStream(bout);
        try {
            dout.writeInt(data.size());
        } catch (IOException e) {
            // never arrive here
        }
        data.stream().forEach(arr -> {
            try {
                dout.writeInt(arr.length);
                for (int i = 0; i < arr.length; i++) {
                    dout.writeFloat(arr[i]);
                }
            } catch (IOException e) {
                e.printStackTrace();
            }
        });
        return bout.toByteArray();
    }

    public static List<float[]> readArrays(byte[] in) {
        DataInputStream input = new DataInputStream(new ByteArrayInputStream(in));
        List<float[]> data = null;
        try {
            int size = input.readInt();
            data = new ArrayList<>(size);
            for (int i = 0; i < size; i++) {
                float[] arr = new float[input.readInt()];
                for (int j = 0; j < arr.length; j++) {
                    arr[j] = input.readFloat();
                }
                data.add(arr);
            }
        } catch (IOException e) {
            // never arrive here
        }
        return data;
    }


    public static String getClassificationString(int result) {

        switch (result){
            case 0: return "boxing";
            case 1: return "handclapping";
            case 2: return "handwaving";
            case 3: return "jogging";
            case 4: return "running";
            case 5: return "walking";
            default: return "unknown";
        }
    }
}


