package tools;

import org.bytedeco.javacpp.opencv_core;
import org.bytedeco.javacpp.opencv_imgproc;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;
import java.util.ArrayList;
import java.util.HashMap;
import java.util.List;
import java.util.Map;
import java.util.stream.Collectors;

import static org.bytedeco.javacpp.opencv_highgui.*;
import static org.bytedeco.javacpp.opencv_highgui.cvWaitKey;

/**
 * Created by Tom.fu on 24/8/2015.
 */
public class functionTest {
    public static void main(String[] args) {

        Map<String, Integer> a = new HashMap<>();

        a.put("a", 1);
        a.put("ab", 2);
        a.put("b", 3);
        a.put("c", 4);

        System.out.println(a.toString());

        List<Integer>x = new ArrayList<>(a.entrySet().stream().filter(e -> e.getKey().contains("a"))
                .collect(Collectors.toMap(Map.Entry::getKey, Map.Entry::getValue)).values());

        System.out.println(x.toString());
        System.out.println(a.toString());

        String fileName1 = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\UCF_combineAction_2\\frame000001.jpg";
        opencv_core.IplImage frame1 = cvLoadImage(fileName1);
        opencv_core.Mat mat1 = new opencv_core.Mat(frame1);

        String fileName2 = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\UCF_combineAction_2\\frame000002.jpg";
        opencv_core.IplImage frame2 = cvLoadImage(fileName2);
        opencv_core.Mat mat2 = new opencv_core.Mat(frame2);


        int w = mat1.cols();
        int h = mat1.rows();

        opencv_core.Mat mat3 = new opencv_core.Mat();
        opencv_core.Size size = new opencv_core.Size(w, 2*h);
        opencv_imgproc.resize(mat1, mat3, size);

        opencv_core.Mat dst_roi1 = new opencv_core.Mat(mat3, new opencv_core.Rect(0, 0, w, h));
        mat1.copyTo(dst_roi1);

        opencv_core.Mat dst_roi2 = new opencv_core.Mat(mat3, new opencv_core.Rect(0, h, w, h));
        mat2.copyTo(dst_roi2);


        cvNamedWindow("DenseTrack");
        cvMoveWindow("DenseTrack", 520, 30);

        cvShowImage("DenseTrack", mat3.asIplImage());
        //cvShowImage("DenseTrack", newFrame);
        cvWaitKey();

//        Mat imgPanel(100, 250, CV_8UC1, Scalar(0));
//        Mat imgPanelRoi(imgPanel, Rect(0, 0, imgSrc.cols, imgSrc.rows));
//        imgSrc.copyTo(imgPanelRoi);
//
//        imshow("imgPanel", imgPanel);
//        waitKey();

        //Mat dst_roi = dst(Rect(left, top, src.cols, src.rows));
        //src.copyTo(dst_roi);
    }
}
