package tools;

import org.bytedeco.javacpp.opencv_core;
import org.bytedeco.javacpp.opencv_imgproc;

import javax.imageio.ImageIO;
import java.awt.image.BufferedImage;
import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;

import static org.bytedeco.javacpp.opencv_highgui.cvLoadImage;

/**
 * Created by Tom.fu on 24/8/2015.
 */
public class FileRenameAndArrangement {
    public static void main(String[] args) {
        String inputFileName = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\combine-for-test.txt";

        String outputFolder = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\UCF_combineAction\\";

        String sourceFolder = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\subset\\";

        String filePrefix = "frame";
        int totalLineCnt = 0;
        int overallIndex = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(inputFileName));
            String rdLine = null;

            while ((rdLine = br.readLine())!=null){
                totalLineCnt ++;
                String[] rdLineSplit = rdLine.split("_");
                String folderName = rdLineSplit[1];

                String src = sourceFolder + folderName + "\\" + rdLine + "\\";
                //String src = sourceFolder + rdLine + "\\";
                File folder = new File(src);
                File[] listOfFiles = folder.listFiles();
                int start_frame_no = 1;
                int end_frame_no = listOfFiles.length;

                for (int i = start_frame_no; i <= end_frame_no; i ++){
                    String orgFileName = src + String.format("%s%06d.jpg", filePrefix, i);
                    int targetFileNameIndex = i + overallIndex;
                    String targetFileName = outputFolder + String.format("%s%06d.jpg", filePrefix, targetFileNameIndex);

//                    opencv_core.IplImage image = cvLoadImage(orgFileName);
//                    opencv_core.Mat matOrg = new opencv_core.Mat(image);
//                    opencv_core.Mat matNew = new opencv_core.Mat();
//                    opencv_core.Size size = new opencv_core.Size(160, 120);
//                    opencv_imgproc.resize(matOrg, matNew, size);
//
//                    File initialImage = new File(targetFileName);
//                    BufferedImage bufferedImage = matNew.getBufferedImage();
//                    ImageIO.write(bufferedImage, "JPEG", initialImage);

                    File from = new File(orgFileName);
                    File to = new File(targetFileName);
                    Files.copy(from.toPath(), to.toPath());
                }
                overallIndex += listOfFiles.length;
                System.out.println("Finished folder: " + rdLine + ", with total files: " + listOfFiles.length + ", overallIndex: "  + overallIndex);
            }
        } catch (IOException ioE){
            ioE.printStackTrace();
        }
    }

}
