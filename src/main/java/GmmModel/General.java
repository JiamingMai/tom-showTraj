package GmmModel;

import redis.clients.jedis.Jedis;

import java.io.*;
import java.util.ArrayList;
import java.util.List;

/**
 * Created by Tom.fu on 16/7/2015.
 */
public class General {

    static double VL_GMM_MIN_VARIANCE_D = 1e-6;
    static double VL_GMM_MIN_POSTERIOR_D = 1e-2;
    static double VL_GMM_MIN_PRIOR_D = 1e-6;

    static double VL_GMM_MIN_VARIANCE_F = 1e-6;
    static double VL_GMM_MIN_POSTERIOR_F = 1e-2;
    static double VL_GMM_MIN_PRIOR_F = 1e-6;

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
//        String modelFile = inputFilePath + "kth_model_row.model";
        String modelFile = inputFilePath + "traj-new-aveDim.model";

        String dataFileFolder = inputFilePath + "trajsTest-new";

        int numDimension = 288;
        int numCluster = 256;

        int fvLength = 2 * numCluster * numDimension;
        int numClass = 6;

//        double[] means = getData_double(fileMeans, numCluster * numDimension);
//        double[] covs = getData_double(fileCovs, numCluster * numDimension);
//        double[] priors = getData_double(filePriors, numCluster);
//        List<double[]> trainingResult = getTrainingResult_double(modelFile, fvLength);
//        check_double(dataFileFolder, trainingResult, means, numDimension, numCluster, covs, priors, 0);

        float[] means = getData_float(fileMeans, numCluster * numDimension);
        float[] covs = getData_float(fileCovs, numCluster * numDimension);
        float[] priors = getData_float(filePriors, numCluster);
        List<float[]> trainingResult = getTrainingResult_float(modelFile, fvLength);

        check_float(dataFileFolder, trainingResult, means, numDimension, numCluster, covs, priors, 0);
//
//////        String testFile = inputFilePath + "trajsTest\\person21_boxing_d3_uncomp.txt";
////        String testFile = inputFilePath + "person25_boxing_d3_uncomp.txt";
////        float[] rawData = checkGetTrajDataWithNormalization_float("192.168.0.30", 6379, "tomQ", 0, 0);
////        float[] rawData = checkGetTrajDataWithNormalization_float(testFile);
//        float[] rawData = getTrajDataWithNormalization_float(testFile, 327);
//
//        String oFile = inputFilePath + "person25_boxing_d3_redis.txt";
////        getTrajDataWithNormalizationToFile_float("192.168.0.30", 6379, "tomQ", 0, -1, 288, oFile);
        float[] rawData = getTrajDataWithNormalization_float("192.168.0.30", 6379, "tomQ", 0, -1, 288);
//        float[] rawData = getTrajDataWithNormalizationFromFile_float(oFile, 288);
//        int numData = rawData.length / numDimension;
//
//        Object[] enc = vl_fisher_encode_float(means, numDimension, numCluster, covs, priors, rawData, numData, 0);
//        float[] fvData = (float[]) enc[0];
//
//        Object[] classifyResult = getClassificationResult_float(trainingResult, fvData);
//        int maxIndex = (int) classifyResult[0];
//        float maxSimilarity = (float) classifyResult[1];
//        System.out.println("numData: " + numData + ", fvData_length: " + fvData.length + ", maxIndex: " + maxIndex + ", maxSim: " + maxSimilarity);
    }

    /*
    ~isempty(strfind(datalist(i).name,'boxing'))		classids(i)=1;
	~isempty(strfind(datalist(i).name,'handclapping'))  classids(i)=2;
	elseif ~isempty(strfind(datalist(i).name,'handwaving'))classids(i)=3;
	elseif ~isempty(strfind(datalist(i).name,'jogging'))  classids(i)=4;
	elseif ~isempty(strfind(datalist(i).name,'running')) classids(i)=5;
	elseif ~isempty(strfind(datalist(i).name,'walking')) classids(i)=6;
 */
    public static double check_double(String testFileFolder, List<double[]> trainingResult,
                                      double[] means, int numDimension, int numCluster, double[] covs, double[] priors, int flags) {
        double acc = 0.0;
        int totalCnt = 0;
        int accCnt = 0;

        File folder = new File(testFileFolder);
        File[] listOfFiles = folder.listFiles();
        for (int i = 0; i < listOfFiles.length; i++) {
            File f = listOfFiles[i];
            if (f.isFile()) {
                int orgClassID = getClassID(f.getName());
                if (orgClassID < 0) {
                    System.out.println("Warning, in check_double, orgClassID < 0, fileName: " + f.getName());
                    continue;
                }

                totalCnt++;
                double[] fv_data = getFvData_double(f.getAbsolutePath(), means, numDimension, numCluster, covs, priors, flags, 167, 243);
                Object[] classifyRestult = getClassificationResult_double(trainingResult, fv_data);
                int testClassID = (int) classifyRestult[0];
                double similarity = (double) classifyRestult[1];
                if (orgClassID == testClassID) {
                    accCnt++;
                    System.out.println("File: " + f.getName() + ", orgClassID: " + orgClassID + ", testID: " + testClassID + ", sim: " + similarity);
                } else {
                    System.out.println("Wrong classification, File: " + f.getName() + ", orgClassID: " + orgClassID + ", testID: " + testClassID + ", sim: " + similarity);
                }
            }
        }
        acc = totalCnt == 0 ? 0.0 : (double) accCnt / (double) totalCnt;
        System.out.println("totalCnt: " + totalCnt + ", accCnt: " + accCnt + ", accuracy: " + acc);
        return acc;
    }

    public static float check_float(String testFileFolder, List<float[]> trainingResult,
                                    float[] means, int numDimension, int numCluster, float[] covs, float[] priors, int flags) {
        double acc = 0.0;
        int totalCnt = 0;
        int accCnt = 0;

        File folder = new File(testFileFolder);
        File[] listOfFiles = folder.listFiles();
        for (int i = 0; i < listOfFiles.length; i++) {
            File f = listOfFiles[i];
            if (f.isFile()) {
                int orgClassID = getClassID(f.getName());
                if (orgClassID < 0) {
                    System.out.println("Warning, in check_double, orgClassID < 0, fileName: " + f.getName());
                    continue;
                }

                totalCnt++;
                //float[] fv_data = getFvData_float(f.getAbsolutePath(), means, numDimension, numCluster, covs, priors, flags, 167, 243);
                float[] fv_data = getFvData_float(f.getAbsolutePath(), means, numDimension, numCluster, covs, priors, flags);
                Object[] classifyRestult = getClassificationResult_float(trainingResult, fv_data);
                int testClassID = (int) classifyRestult[0];
                float similarity = (float) classifyRestult[1];
                if (orgClassID == testClassID) {
                    accCnt++;
                    System.out.println("File: " + f.getName() + ", orgClassID: " + orgClassID + ", testID: " + testClassID + ", sim: " + similarity);
                } else {
                    System.out.println("Wrong classification, File: " + f.getName() + ", orgClassID: " + orgClassID + ", testID: " + testClassID + ", sim: " + similarity);
                }
            }
        }
        acc = totalCnt == 0 ? 0.0 : (double) accCnt / (double) totalCnt;
        System.out.println("check_float, totalCnt: " + totalCnt + ", accCnt: " + accCnt + ", accuracy: " + acc);
        return (float) acc;
    }

    public static Object[] getClassificationResult_double(List<double[]> trainingResult, double[] fv_data) {
        int maxIndex = 0;
        double maxSimilarity = Double.NEGATIVE_INFINITY;

        for (int i = 0; i < trainingResult.size(); i++) {
            double sim = innerProduct(fv_data, trainingResult.get(i));
            if (sim > maxSimilarity) {
                maxSimilarity = sim;
                maxIndex = i;
            }
        }

        return new Object[]{maxIndex, maxSimilarity};
    }

    public static Object[] getClassificationResult_float(List<float[]> trainingResult, float[] fv_data) {
        int maxIndex = 0;
        float maxSimilarity = Float.NEGATIVE_INFINITY;

        for (int i = 0; i < trainingResult.size(); i++) {
            float sim = innerProduct(fv_data, trainingResult.get(i));
            System.out.println("i: " + i + ", sim: " + sim);
            if (sim > maxSimilarity) {
                maxSimilarity = sim;
                maxIndex = i;
            }
        }
        return new Object[]{maxIndex, maxSimilarity};
    }

    public static double[] getFvData_double(
            String dataFile, double[] means, int numDimension, int numCluster, double[] covs, double[] priors, int flags) {
        double[] data = getTrajDataWithNormalization_double(dataFile, 327);
        int numData = data.length / numDimension;

        Object[] enc = vl_fisher_encode_double(means, numDimension, numCluster, covs, priors, data, numData, flags);
        return (double[]) enc[0];
    }

    public static float[] getFvData_float(
            String dataFile, float[] means, int numDimension, int numCluster, float[] covs, float[] priors, int flags) {
        float[] data = getTrajDataWithNormalization_float(dataFile, 327);
        int numData = data.length / numDimension;

        Object[] enc = vl_fisher_encode_float(means, numDimension, numCluster, covs, priors, data, numData, flags);
        return (float[]) enc[0];
    }

    public static double[] getFvData_double(
            String dataFile, double[] means, int numDimension, int numCluster, double[] covs, double[] priors, int flags, int startFrameID, int endFrameID) {
        double[] data = getTrajDataWithNormalization_double(dataFile, 327, startFrameID, endFrameID);
        int numData = data.length / numDimension;

        Object[] enc = vl_fisher_encode_double(means, numDimension, numCluster, covs, priors, data, numData, flags);
        return (double[]) enc[0];
    }

    public static float[] getFvData_float(
            String dataFile, float[] means, int numDimension, int numCluster, float[] covs, float[] priors, int flags, int startFrameID, int endFrameID) {
        float[] data = getTrajDataWithNormalization_float(dataFile, 327, startFrameID, endFrameID);
        int numData = data.length / numDimension;

        Object[] enc = vl_fisher_encode_float(means, numDimension, numCluster, covs, priors, data, numData, flags);
        return (float[]) enc[0];
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

    public static float innerProduct(float[] a, float[] b) {
        if (a.length != b.length) {
            System.out.println("warning, a.length != b.length, " + a.length + ", " + b.length);
            return -1;
        }
        float retVal = 0;
        for (int i = 0; i < a.length; i++) {
            retVal += a[i] * b[i];
        }
        return retVal;
    }

    public static double[] getData_double(String fileName, int checkLength) {
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

    public static float[] getData_float(String fileName, int checkLength) {
        List<Float> dataList = new ArrayList<>();
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String[] rdLineSplit = myfile.readLine().split(" ");
            for (int i = 0; i < rdLineSplit.length; i++) {
                if (!rdLineSplit[i].trim().isEmpty()) {
                    dataList.add(Float.parseFloat(rdLineSplit[i].trim()));
                }
            }
            myfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }

        if (checkLength > 0 && dataList.size() != checkLength) {
            System.out.println("warning, dataList.size() != checkLength, " + dataList.size() + ", " + checkLength);
        }

        float[] retVal = new float[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).floatValue();
        }
        return retVal;
    }

    public static double[] getTrajData_double(String fileName, int checkLength) {
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

    public static float[] getTrajData_float(String fileName, int checkLength) {
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
        float[] retVal = new float[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).floatValue();
        }
        return retVal;
    }

    public static double[] getTrajDataWithNormalization_double(String fileName, int checkLength) {
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

    public static float[] getTrajDataWithNormalization_float(String fileName, int checkLength) {
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

                float sum = 0;
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        float p = (float) Math.sqrt(Math.abs(dList.get(i)));
                        p = dList.get(i) > 0 ? p : -p;
                        sum += p * p;
                        dList.set(i, p);
                    }
                }

                float base = (float) Math.sqrt(sum);
                ///only remain the
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        if (base > 0) {
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
        float[] retVal = new float[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).floatValue();
        }
        return retVal;
    }

    public static float[] getTrajDataWithNormalizationFromFile_float(String fileName, int checkLength) {
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

                float sum = 0;
                for (int i = 0; i < dList.size(); i++) {
                    float p = (float) Math.sqrt(Math.abs(dList.get(i)));
                    p = dList.get(i) > 0 ? p : -p;
                    sum += p * p;
                    dList.set(i, p);
                }

                float base = (float) Math.sqrt(sum);
                ///only remain the
                for (int i = 0; i < dList.size(); i++) {
                    if (base > 0) {
                        dataList.add(dList.get(i) / base);
                    } else {
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

    public static float[] checkGetTrajDataWithNormalization_float(String fileName) {
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

                float sum = 0;
                System.out.println(lineCount + "-" + dList.size());
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        float tmp = dList.get(i);
                        float p = (float) Math.sqrt(Math.abs(dList.get(i)));
                        p = dList.get(i) > 0 ? p : -p;
                        sum += p * p;
                        dList.set(i, p);
                        System.out.print(tmp + " ");
                    }
                }
                System.out.println();

//                if (lineCount == 2) return null;
                float base = (float) Math.sqrt(sum);
                ///only remain the
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        if (base > 0) {
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
        float[] retVal = new float[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).floatValue();
        }
        return retVal;
    }

    public static float[] checkGetTrajDataWithNormalization_float(String host, int port, String queue, int startIndex, int endIndex) {

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
                lineCount++;
                float sum = 0;
                System.out.println(lineCount + "-" + trace.length);
                for (int i = 0; i < trace.length; i++) {
                    float tmp = trace[i];
                    float p = (float) Math.sqrt(Math.abs(trace[i]));
                    p = trace[i] > 0 ? p : -p;
                    sum += p * p;
                    trace[i] = p;

                    System.out.print(tmp + " ");
                }
                System.out.println();

                float base = (float) Math.sqrt(sum);
                ///only remain the
                for (int i = 0; i < trace.length; i++) {
                    if (base > 0) {
                        dataList.add(trace[i] / base);
                    } else {
                        dataList.add(trace[i]);
                    }

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
     * @param checkLength
     * @return
     */
    public static float[] getTrajDataWithNormalization_float(String host, int port, String queue, int startIndex, int endIndex, int checkLength) {

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
                float sum = 0;
                for (int i = 0; i < trace.length; i++) {
                    float p = (float) Math.sqrt(Math.abs(trace[i]));
                    p = trace[i] > 0 ? p : -p;
                    sum += p * p;
                    trace[i] = p;
                }

                float base = (float) Math.sqrt(sum);
                ///only remain the
                for (int i = 0; i < trace.length; i++) {
                    if (base > 0) {
                        dataList.add(trace[i] / base);
                    } else {
                        dataList.add(trace[i]);
                    }
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

    public static void getTrajDataWithNormalizationToFile_float(String host, int port, String queue, int startIndex, int endIndex, int checkLength, String fileName) {

        ///String host = "192.168.0.30";
        ///int port = 6379;
        ///String queue = "tomQ";
        byte[] queueName = queue.getBytes();
        Jedis jedis = new Jedis(host, port);
        List<byte[]> getData = jedis.lrange(queueName, startIndex, endIndex);
        if (getData == null) {
            System.out.println("Warning! in getTrajDataWithNormalization_float, getData == null, queueName: " + queue);
        }

        try {
            BufferedWriter myfile = new BufferedWriter(new FileWriter(fileName));


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
                    float sum = 0;
                    for (int i = 0; i < trace.length; i++) {
                        myfile.write(trace[i] + " ");
                    }
                    myfile.newLine();
                }
            }
            System.out.println("DataSize: " + getData.size() + ", lineCount: " + lineCount + ", total data size: " + dataList.size());
            myfile.close();
        } catch (IOException e) {
            e.printStackTrace();
        }
    }

    public static double[] getTrajDataWithNormalization_double(String fileName, int checkLength, int startFrameID, int endFrameID) {
        List<Double> dataList = new ArrayList<>();
        int lineCount = 0;
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
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

                int frameID = (int) dList.get(0).intValue();
                if (frameID < startFrameID || frameID > endFrameID) {
                    continue;
                }

                lineCount++;

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

    public static float[] getTrajDataWithNormalization_float(String fileName, int checkLength, int startFrameID, int endFrameID) {
        List<Float> dataList = new ArrayList<>();
        int lineCount = 0;
        try {
            BufferedReader myfile = new BufferedReader(new FileReader(fileName));
            String rdLine = null;
            while ((rdLine = myfile.readLine()) != null) {
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

                int frameID = (int) dList.get(0).intValue();
                if (frameID < startFrameID || frameID > endFrameID) {
                    continue;
                }

                lineCount++;

                float sum = 0;
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        float p = (float) Math.sqrt(Math.abs(dList.get(i)));
                        p = dList.get(i) > 0 ? p : -p;
                        sum += p * p;
                        dList.set(i, p);
                    }
                }

                float base = (float) Math.sqrt(sum);
                ///only remain the
                for (int i = 0; i < dList.size(); i++) {
                    if (i >= 39) {
                        if (base > 0) {
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
        float[] retVal = new float[dataList.size()];
        for (int i = 0; i < dataList.size(); i++) {
            retVal[i] = dataList.get(i).floatValue();
        }
        return retVal;
    }

    public static List<double[]> getTrainingResult_double(String fileName, int checkLength) {
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

    public static float VlDistanceMahalanobis(
            float[] x, float[] mu, float[] s) {
        double acc = 0;
        for (int i = 0; i < x.length; i++) {
            double diff = x[i] - mu[i];
            acc += (diff * diff * s[i]);
        }
        return (float) acc;
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

        //TODO: check_double the size of this matrix
        double[] posteriors = new double[numClusters * numData];

        double[] logCovariances = new double[numClusters];
        double[] invCovariances = new double[numClusters * dimension];
        double[] logWeights = new double[numClusters];

        for (int i_cl = 0; i_cl < numClusters; i_cl++) {
            double logSigma = 0.0;
            if (priors[i_cl] < VL_GMM_MIN_PRIOR_D) {
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

    public static float[] transpose(float[] org, int orgRow, int orgCol) {
        if (org.length != orgRow * orgCol) {
            System.out.println("warning!! input matrix dimension is incorrect!!, org.length: "
                    + org.length + ", orgRow: " + orgRow + ", orgCol: " + orgCol);
            return null;
        }

        float[] retVal = new float[orgCol * orgRow];
        for (int i = 0; i < orgRow; i++) {
            for (int j = 0; j < orgCol; j++) {
                int orgIndex = i * orgCol + j;
                int transIndex = j * orgRow + i;
                retVal[transIndex] = org[orgIndex];
            }
        }
        return retVal;
    }

    public static float[] l2Normalization(float[] input) {
        float[] retVal = new float[input.length];
        float sum = 0;
        for (int i = 0; i < input.length; i++) {
            sum += (input[i] * input[i]);
        }
        float base = (float) Math.sqrt(sum);
        for (int i = 0; i < input.length; i++) {
            input[i] = input[i] / base;
        }
        return retVal;
    }

}























