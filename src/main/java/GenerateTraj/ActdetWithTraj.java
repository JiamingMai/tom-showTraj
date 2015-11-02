package GenerateTraj;

import GmmModel.GmmData;
import GmmModel.PcaData;
import org.bytedeco.javacpp.opencv_core;
import org.bytedeco.javacpp.opencv_imgproc;
import org.bytedeco.javacpp.opencv_video;
import org.bytedeco.javacv.FrameRecorder;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.BufferedWriter;
import java.io.File;
import java.io.FileWriter;
import java.io.IOException;
import java.nio.FloatBuffer;
import java.util.ArrayList;
import java.util.Arrays;
import java.util.LinkedList;
import java.util.List;

import static GenerateTraj.FeatureExtra.*;
import static GmmModel.newMethod.*;
import static org.bytedeco.javacpp.opencv_core.*;
import static org.bytedeco.javacpp.opencv_highgui.*;


/**
 * Created by Tom.fu on 28/11/2014.
 */
public class ActdetWithTraj {
    static int patch_size = 32;
    static int nxy_cell = 2;
    static int nt_cell = 3;
    static int dimension = 32;
    static boolean fullOrientation = true;
    static float epsilon = 0.05f;
    static float min_flow = 0.4f * 0.4f;

    static int start_frame = 0;
    static int end_frame = 1000000;
    static double quality = 0.001;
    static double min_distance = 5; //13;
    static int init_gap = 1; //3;

    static int track_length = 15; //15; update (bingbing longer traj)

    static int scale_num = 1;
    static float scale_stride = (float) Math.sqrt(2.0);


    public static void main(String[] args) throws FrameRecorder.Exception {

        int drawIndicator = 1;


        long startTime = System.currentTimeMillis();
        //System.out.println("new test, start at: " + startTime);

        String inputFolder = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\UCF_combineAction\\";
        String outputFileName = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\empty.txt";
        String filePrefix = "frame";

        int winPosX = Integer.parseInt(args[0]);
        int winPosY = Integer.parseInt(args[1]);

        ///for comparison demo, winPosX = 1200, winPosY = 113
        ActdetWithTrajDrawing(drawIndicator, inputFolder, outputFileName, filePrefix,
                320, 240, 160, 120, 480, 360, winPosX, winPosY, null, false, null);

        long endTime = System.currentTimeMillis();
        System.out.println("finished with duration: " + (endTime - startTime));
    }

    ///show_track == 1, draw in window and in file
    ///show_track > 1, only generate frames, not drawing in a window
    public static void ActdetWithTrajDrawing(
            int show_track, String sourceVideoFolder, String outputFile, String filePrefix,
            int inW, int inH, int procW, int procH, int outW, int outH, int winPosX, int winPosY,
            String recordFile, boolean toOutputTrajInfo, String outputFrameFolder) throws FrameRecorder.Exception {

        String inputFilePath = "C:\\Users\\Tom.fu\\Desktop\\fromPeiYong\\";

        String demoFolder = "ucf101-Demo";

        String hogPcaFile = inputFilePath + demoFolder + "\\ucf_codebook\\hog_pca.mat";
        String mbhxPcaFile = inputFilePath + demoFolder + "\\ucf_codebook\\mbhx_pca.mat";
        String mbhyPcaFile = inputFilePath + demoFolder + "\\ucf_codebook\\mbhy_pca.mat";
        String hogGmmFile = inputFilePath + demoFolder + "\\ucf_codebook\\hog_gmm.codebook";
        String mbhxGmmFile = inputFilePath + demoFolder + "\\ucf_codebook\\mbhx_gmm.codebook";
        String mbhyGmmFile = inputFilePath + demoFolder + "\\ucf_codebook\\mbhy_gmm.codebook";

//        String modelFile = inputFilePath + demoFolder + "\\KTHDemo2_row.model";
        String modelFile = inputFilePath + demoFolder + "\\ucf101_demo_row.model";
        String classNameFile = inputFilePath + demoFolder + "\\ucf101_actNames.txt";

        List<String> actionNameList = getClassNames(classNameFile);
//        for (int i = 0; i < actionNameList.size(); i ++){
//            System.out.println("i: " + i + ", " + actionNameList.get(i));
//        }

//        String dataFileFolder = inputFilePath + "trajsTest-new";
        String dataFileFolder = inputFilePath + "testForUcf101";

        int numDimension = 288;
        int numCluster = 256;

        int fvLength = 2 * numCluster * numDimension / 2;

        int frameRate = 25;
        int windowInSeconds = 4;
        int windowInFrames = windowInSeconds * frameRate;

        int resultLastSeconds = 2;
        int countDownSeconds = windowInSeconds - resultLastSeconds;

        PcaData hogPca = new PcaData(hogPcaFile);
        PcaData mbhxPca = new PcaData(mbhxPcaFile);
        PcaData mbhyPca = new PcaData(mbhyPcaFile);

        GmmData hogGmm = new GmmData(hogGmmFile);
        GmmData mbhxGmm = new GmmData(mbhxGmmFile);
        GmmData mbhyGmm = new GmmData(mbhyGmmFile);

        List<float[]> trainingResult = getTrainingResult_float(modelFile, fvLength);

        int CurrentClass = 0;

        IplImage frame = null;
        IplImage image = null;
        IplImage trajImage = null;
        IplImage prev_image = null;
        IplImage actImage = null;

        IplImage grey = null;
        IplImage prev_grey = null;
        IplImagePyramid grey_pyramid = null;
        IplImagePyramid prev_grey_pyramid = null;
        IplImagePyramid eig_pyramid = null;

        float[] fscales = new float[scale_num];
        List<LinkedList<Track>> xyScaleTracks = new ArrayList<>();

        DescInfo hogInfo = new DescInfo(8, 0, 1, patch_size, nxy_cell, nt_cell, min_flow);
        DescInfo mbhInfo = new DescInfo(8, 0, 1, patch_size, nxy_cell, nt_cell, min_flow);
        TrackerInfo tracker = new TrackerInfo(track_length, init_gap);
        int nChannel = 3;
        int nDepth = 8;

        float adjX = (float) inW / (float) procW;
        float adjY = (float) inH / (float) procH;

        File folder = new File(sourceVideoFolder);
        File[] listOfFiles = folder.listFiles();
        int start_frame_no = 1;
        int end_frame_no = listOfFiles.length;
        //System.out.println("the input folder: " + sourceVideoFolder + " has totally files: " + end_frame_no + ", outputFile: " + outputFile);

        FrameRecorder recorder = null;
        if (recordFile != null) {
            recorder = FrameRecorder.createDefault(recordFile, outW, outH);
            recorder.setFrameRate(15);
            recorder.setVideoQuality(1.0);
            recorder.start();
        }

        List<float[]> traceFeatures = new ArrayList<>();

        int frameFileIndex = start_frame_no;
        int init_counter = 0;
        int frameNum = 0;

        if (show_track == 1) {
            cvNamedWindow("DenseTraj-ActionDet");
            cvMoveWindow("DenseTraj-ActionDet", winPosX, winPosY);
        }

        long startTime = System.currentTimeMillis();
        System.out.println("Start stand alone version of action detection with drawing dense trajectories, at: " + startTime);
        try {
            BufferedWriter myfile = new BufferedWriter(new FileWriter(outputFile));
            while (frameFileIndex <= end_frame_no) {
                int i, j, c;

                String fileName = sourceVideoFolder + String.format("%s%06d.jpg", filePrefix, frameFileIndex);
                IplImage ppimage = cvLoadImage(fileName);
                IplImage inputImage = cvLoadImage(fileName);
                frame = cvCreateImage(cvSize(procW, procH), nDepth, nChannel);
                opencv_imgproc.cvResize(ppimage, frame, opencv_imgproc.CV_INTER_AREA);

                if (frameNum >= start_frame && frameNum <= end_frame) {
                    if (image == null) {
                        image = cvCreateImage(cvGetSize(frame), 8, 3);
                        image.origin(frame.origin());
                        prev_image = cvCreateImage(cvGetSize(frame), 8, 3);
                        prev_image.origin(frame.origin());

                        grey = cvCreateImage(cvGetSize(frame), 8, 1);
                        grey_pyramid = new IplImagePyramid(scale_stride, scale_num, cvGetSize(frame), 8, 1);
                        prev_grey = cvCreateImage(cvGetSize(frame), 8, 1);
                        prev_grey_pyramid = new IplImagePyramid(scale_stride, scale_num, cvGetSize(frame), 8, 1);
                        eig_pyramid = new IplImagePyramid(scale_stride, scale_num, cvGetSize(frame), 32, 1);

                        cvCopy(frame, image, null);
                        opencv_imgproc.cvCvtColor(image, grey, opencv_imgproc.CV_BGR2GRAY);
                        grey_pyramid.rebuild(grey);

                        for (int ixyScale = 0; ixyScale < scale_num; ++ixyScale) {
                            LinkedList<Track> tracks = new LinkedList<>();
                            fscales[ixyScale] = (float) Math.pow(scale_stride, ixyScale);

                            IplImage grey_temp = cvCloneImage(grey_pyramid.getImage(ixyScale));
                            IplImage eig_temp = cvCloneImage(eig_pyramid.getImage(ixyScale));
                            LinkedList<CvPoint2D32f> points = cvDenseSample(grey_temp, eig_temp, quality, min_distance);

                            for (i = 0; i < points.size(); i++) {
                                Track track = new Track(tracker.trackLength);
                                PointDesc point = new PointDesc(mbhInfo, hogInfo, points.get(i));
                                track.addPointDesc(point);
                                tracks.addLast(track);
                            }
                            xyScaleTracks.add(tracks);
                            cvReleaseImage(grey_temp);
                            cvReleaseImage(eig_temp);
                        }
                    }

                    cvCopy(frame, image, null);
                    opencv_imgproc.cvCvtColor(image, grey, opencv_imgproc.CV_BGR2GRAY);
                    grey_pyramid.rebuild(grey);

                    if (frameNum > 0) {
                        init_counter++;
                        for (int ixyScale = 0; ixyScale < scale_num; ixyScale++) {
                            LinkedList<CvPoint2D32f> points_in = new LinkedList<>();
                            LinkedList<Track> tracks = xyScaleTracks.get(ixyScale);

                            for (int tempI = 0; tempI < tracks.size(); tempI++) {
                                Track iTrack = tracks.get(tempI);
                                CvPoint2D32f point = new CvPoint2D32f(iTrack.pointDescs.getLast().point);
                                points_in.addLast(point);
                            }

                            IplImage prev_grey_temp = cvCloneImage(prev_grey_pyramid.getImage(ixyScale));
                            IplImage grey_temp = cvCloneImage(grey_pyramid.getImage(ixyScale));

                            IplImage flow = cvCreateImage(cvGetSize(grey_temp), IPL_DEPTH_32F, 2);
                            opencv_video.cvCalcOpticalFlowFarneback(prev_grey_temp, grey_temp, flow,
                                    Math.sqrt(2.0) / 2.0, 5, 10, 2, 7, 1.5, opencv_video.OPTFLOW_FARNEBACK_GAUSSIAN);

                            Object[] pOutWithStatus = OpticalFlowTrackerSimple(flow, points_in);
                            LinkedList<CvPoint2D32f> points_out = (LinkedList<CvPoint2D32f>) pOutWithStatus[0];
                            int[] status = (int[]) pOutWithStatus[1];

                            int width = grey_temp.width();
                            int height = grey_temp.height();

                            DescMat[] mbhMatXY = MbhComp(flow, mbhInfo, width, height);
                            DescMat mbhMatX = mbhMatXY[0];
                            DescMat mbhMatY = mbhMatXY[1];

                            DescMat hogMat = HogComp(prev_grey_temp, hogInfo, width, height);

                            int tSize = tracks.size();
                            int tSizeCnt = 0;
                            for (i = 0; i < tSize; i++) {
                                Track iTrack = tracks.get(tSizeCnt);
                                if (status[i] == 1) {
                                    CvPoint2D32f prev_point = points_in.get(i);
                                    CvScalar rect = getRect(prev_point, cvSize(width, height), hogInfo);

                                    PointDesc pointDesc = iTrack.pointDescs.getLast();
                                    pointDesc.mbhX = getDesc(mbhMatX, rect, mbhInfo);
                                    pointDesc.mbhY = getDesc(mbhMatY, rect, mbhInfo);
                                    pointDesc.hog = getDesc(hogMat, rect, hogInfo);

                                    PointDesc point = new PointDesc(mbhInfo, hogInfo, points_out.get(i));
                                    iTrack.addPointDesc(point);

                                    if (show_track > 0) {
                                        float length = iTrack.pointDescs.size();

                                        float point0_x = fscales[ixyScale] * iTrack.pointDescs.get(0).point.x();
                                        float point0_y = fscales[ixyScale] * iTrack.pointDescs.get(0).point.y();
                                        CvPoint2D32f point0 = new CvPoint2D32f();
                                        point0.x(point0_x * adjX);
                                        point0.y(point0_y * adjY);

                                        Integer t = (int) point0_x * 1023 + (int) point0_y * 59;
                                        if (Math.abs(t.hashCode()) % 5 == 0) {
                                            continue;
                                        }

                                        float jIndex = 0;
                                        for (int jj = 1; jj < length; jj++, jIndex++) {
                                            float point1_x = fscales[ixyScale] * iTrack.pointDescs.get(jj).point.x();
                                            float point1_y = fscales[ixyScale] * iTrack.pointDescs.get(jj).point.y();
                                            CvPoint2D32f point1 = new CvPoint2D32f();
                                            point1.x(point1_x * adjX);
                                            point1.y(point1_y * adjY);


//                                            cvLine(image, cvPointFrom32f(point0), cvPointFrom32f(point1),
//                                                    CV_RGB(0, cvFloor(255.0 * (jIndex + 1.0) / length), 0), 1, 8, 0);
                                            cvLine(ppimage, cvPointFrom32f(point0), cvPointFrom32f(point1),
                                                    CV_RGB(0, cvFloor(255.0 * (jIndex + 1.0) / length), 0), 1, 8, 0);
                                            point0 = point1;
                                        }
                                    }
                                    tSizeCnt++;
                                } else {
                                    tracks.remove(iTrack);
                                }
                            }
                            cvReleaseImage(prev_grey_temp);
                            cvReleaseImage(grey_temp);
                            cvReleaseImage(flow);
                        }

                        //output
                        for (int ixyScale = 0; ixyScale < scale_num; ixyScale++) {
                            LinkedList<Track> tracks = xyScaleTracks.get(ixyScale);
                            int tSizeBefore = tracks.size();
                            int tSizeCnt = 0;
                            int indCnt = 0;
                            int removeSize = 0;
                            for (i = 0; i < tSizeBefore; i++) {
                                Track iTrack = tracks.get(tSizeCnt);
                                // if the trajectory achieves the length we want
                                if (iTrack.pointDescs.size() >= tracker.trackLength + 1) {
                                    CvPoint2D32f[] trajectory = new CvPoint2D32f[tracker.trackLength + 1];
                                    for (int count = 0; count <= tracker.trackLength; count++) {
                                        PointDesc iDesc = iTrack.pointDescs.get(count);
                                        trajectory[count] = new CvPoint2D32f();
                                        trajectory[count].x(iDesc.point.x() * fscales[ixyScale]);
                                        trajectory[count].y(iDesc.point.y() * fscales[ixyScale]);
                                    }

                                    ///**
                                    //return new Object[]{bk2, new float[]{mean_x, mean_y, var_x, var_y, length}, 1};
                                    Object[] ret = isValid(trajectory);
                                    CvPoint2D32f[] trajectory_backup = (CvPoint2D32f[]) ret[0];
                                    float[] vals = (float[]) ret[1];
                                    int indicator = (int) ret[2];

                                    if (indicator == 1) {
                                        indCnt++;

                                        float[] Track_Info = new float[]{
                                                (float) frameNum, vals[0], vals[1], vals[2], vals[3], vals[4], fscales[ixyScale]};

                                        float[] XYs = new float[(tracker.trackLength + 1) * 2];
                                        for (int count = 0; count < tracker.trackLength + 1; count++) {
                                            XYs[count * 2] = trajectory[count].x();
                                            XYs[count * 2 + 1] = trajectory[count].y();
                                        }

                                        List<Float> hogFeature = new ArrayList<>();
                                        int t_stride = cvFloor(tracker.trackLength / hogInfo.ntCells);
                                        int iDescIndex = 0;
                                        for (int n = 0; n < hogInfo.ntCells; n++) {
                                            float[] vec = new float[hogInfo.dim];
                                            for (int m = 0; m < hogInfo.dim; m++) {
                                                vec[m] = 0;
                                            }
                                            for (int t = 0; t < t_stride; t++, iDescIndex++) {
                                                PointDesc iDesc = iTrack.pointDescs.get(iDescIndex);
                                                for (int m = 0; m < hogInfo.dim; m++) {
                                                    vec[m] += iDesc.hog[m];
                                                }
                                            }
                                            for (int m = 0; m < mbhInfo.dim; m++) {
                                                hogFeature.add(vec[m] / (float) t_stride);
                                            }
                                        }

                                        List<Float> mbhFeature = new ArrayList<>();
                                        t_stride = cvFloor(tracker.trackLength / mbhInfo.ntCells);
                                        iDescIndex = 0;
                                        for (int n = 0; n < mbhInfo.ntCells; n++) {
                                            float[] vec = new float[mbhInfo.dim];
                                            for (int m = 0; m < mbhInfo.dim; m++) {
                                                vec[m] = 0;
                                            }
                                            for (int t = 0; t < t_stride; t++, iDescIndex++) {
                                                PointDesc iDesc = iTrack.pointDescs.get(iDescIndex);
                                                for (int m = 0; m < mbhInfo.dim; m++) {
                                                    vec[m] += iDesc.mbhX[m];
                                                }
                                            }
                                            for (int m = 0; m < mbhInfo.dim; m++) {
                                                mbhFeature.add(vec[m] / (float) t_stride);
                                            }
                                        }

                                        t_stride = cvFloor(tracker.trackLength / mbhInfo.ntCells);
                                        iDescIndex = 0;
                                        for (int n = 0; n < mbhInfo.ntCells; n++) {
                                            float[] vec = new float[mbhInfo.dim];
                                            for (int m = 0; m < mbhInfo.dim; m++) {
                                                vec[m] = 0;
                                            }
                                            for (int t = 0; t < t_stride; t++, iDescIndex++) {
                                                PointDesc iDesc = iTrack.pointDescs.get(iDescIndex);
                                                for (int m = 0; m < mbhInfo.dim; m++) {
                                                    vec[m] += iDesc.mbhY[m];
                                                }
                                            }
                                            for (int m = 0; m < mbhInfo.dim; m++) {
                                                mbhFeature.add(vec[m] / (float) t_stride);
                                            }
                                        }

                                        float[] allFeatures = new float[nt_cell * dimension * 3];
                                        for (int iii = 0; iii < hogFeature.size(); iii++) {
                                            allFeatures[iii] = hogFeature.get(iii);
                                        }
                                        for (int iii = 0; iii < mbhFeature.size(); iii++) {
                                            allFeatures[hogFeature.size() + iii] = mbhFeature.get(iii);
                                        }

                                        traceFeatures.add(allFeatures);

                                        if (toOutputTrajInfo) {
                                            WriteTrajFeature2Txt(myfile, Track_Info, XYs, hogFeature, mbhFeature);
                                        }
                                    }
                                    tracks.remove(iTrack);
                                    removeSize++;
                                } else {
                                    tSizeCnt++;
                                }
                            }
                        }

                        if (init_counter == tracker.initGap) {
                            init_counter = 0;
                            for (int ixyScale = 0; ixyScale < scale_num; ixyScale++) {
                                LinkedList<Track> tracks = xyScaleTracks.get(ixyScale);
                                LinkedList<CvPoint2D32f> points_in = new LinkedList<>();
                                for (i = 0; i < tracks.size(); i++) {
                                    Track iTrack = tracks.get(i);
                                    CvPoint2D32f point = new CvPoint2D32f(iTrack.pointDescs.getLast().point);
                                    points_in.addLast(point);
                                }

                                IplImage grey_temp = cvCloneImage(grey_pyramid.getImage(ixyScale));
                                IplImage eig_temp = cvCloneImage(eig_pyramid.getImage(ixyScale));

                                LinkedList<CvPoint2D32f> points_out =
                                        cvDenseSample(grey_temp, eig_temp, points_in, quality, min_distance);

                                for (i = 0; i < points_out.size(); i++) {
                                    Track track = new Track(tracker.trackLength);
                                    PointDesc point = new PointDesc(mbhInfo, hogInfo, points_out.get(i));
                                    track.addPointDesc(point);
                                    tracks.addLast(track);
                                }
                                cvReleaseImage(grey_temp);
                                cvReleaseImage(eig_temp);
                            }
                        }

                        if (frameNum % windowInFrames == 0) {
                            CurrentClass = checkNew_float(traceFeatures, trainingResult,
                                    numDimension, hogPca, mbhxPca, mbhyPca, hogGmm, mbhxGmm, mbhyGmm, false);
                            traceFeatures.clear();
                        }
                    }

                    cvCopy(frame, prev_image, null);
                    opencv_imgproc.cvCvtColor(prev_image, prev_grey, opencv_imgproc.CV_BGR2GRAY);
                    prev_grey_pyramid.rebuild(prev_grey);
                    frameNum++;
                }

                trajImage = cvCreateImage(cvSize(outW, outH), nDepth, nChannel);
                opencv_imgproc.cvResize(ppimage, trajImage, opencv_imgproc.CV_INTER_AREA);

                actImage = cvCreateImage(cvSize(outW, outH), nDepth, nChannel);
                opencv_imgproc.cvResize(inputImage, actImage, opencv_imgproc.CV_INTER_AREA);

                CvFont font = new CvFont();
                cvInitFont(font, CV_FONT_VECTOR0, 1.2f, 1.2f, 0, 2, 8);
                CvPoint showPos = cvPoint(5, 40);
                CvScalar showColor = CvScalar.RED;
                CvPoint showPos2 = cvPoint(5, outH - 15);

                if (frameNum < track_length + resultLastSeconds * frameRate) {
                    cvPutText(actImage, "Action Detection", showPos, font, showColor);
                } else {
                    int adjFrameID = frameNum - track_length - resultLastSeconds * frameRate; ///window is 75, 0-14, 15-29, 30-44, 45-59, 60-74
                    //int winIndex = adjFrameID / windowInFrames;
                    int secPos = (adjFrameID % windowInFrames) / frameRate;

                    //3, 2, 1, x, x, 3, 2, 1, x, x,
                    if (secPos < countDownSeconds) {//
                        int showSecondInfo = countDownSeconds - secPos;
                        int t = windowInSeconds - showSecondInfo;
                        int percent = t * 100 / windowInSeconds;
                        cvPutText(actImage, "Detecting action... " + percent + "%", showPos2, font, showColor);
                    } else {
                        cvPutText(actImage, "Action: " + getClassificationString(CurrentClass, actionNameList), showPos, font, showColor);
                    }
                }

                opencv_core.Mat trajMat = new opencv_core.Mat(trajImage);
                opencv_core.Mat actMat = new opencv_core.Mat(actImage);

                opencv_core.Mat combineMat = new opencv_core.Mat();
                opencv_core.Size size = new opencv_core.Size(outW, 2 * outH);
                opencv_imgproc.resize(actMat, combineMat, size);

                opencv_core.Mat dst_roi1 = new opencv_core.Mat(combineMat, new opencv_core.Rect(0, 0, outW, outH));
                trajMat.copyTo(dst_roi1);

                opencv_core.Mat dst_roi2 = new opencv_core.Mat(combineMat, new opencv_core.Rect(0, outH, outW, outH));
                actMat.copyTo(dst_roi2);

                if (show_track == 1) {
                    cvShowImage("DenseTraj-ActionDet", combineMat.asIplImage());
                    c = cvWaitKey(1);
                    if ((char) c == 27) {
                        break;
                    }
                }
                if (recordFile != null) {
                    recorder.record(trajImage);
                }
                if (show_track > 0 && outputFrameFolder != null) {
                    String outputFrameFile = outputFrameFolder + String.format("%s%06d.jpg", filePrefix, frameFileIndex);
                    File initialImage = new File(outputFrameFile);
                    try {
                        Mat mat = new Mat(trajImage);
                        BufferedImage bufferedImage = mat.getBufferedImage();
                        ImageIO.write(bufferedImage, "JPEG", initialImage);
                    } catch (IOException e) {
                        e.printStackTrace();
                    }
                }
                frameFileIndex++;
                //long current = System.currentTimeMillis();

                System.out.println("frameID: " + frameNum + ", timeStamp: " + System.currentTimeMillis());
            }
            myfile.close();
            if (recordFile != null) {
                recorder.stop();
                recorder.release();
            }

            if (show_track == 1) {
                cvDestroyWindow("DenseTraj-ActionDet");
            }
        } catch (IOException e) {
            e.printStackTrace();
        } catch (Exception e) {
            e.printStackTrace();
        }
    }

}
