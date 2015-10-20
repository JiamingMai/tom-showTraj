package tools;

import java.io.BufferedReader;
import java.io.File;
import java.io.FileReader;
import java.io.IOException;
import java.nio.file.Files;

/**
 * Created by Tom.fu on 24/8/2015.
 */
public class FileRenameAndArrangement {
    public static void main(String[] args) {
        String inputFileName = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\combine-for-test.txt";

        String outputFolder = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\combineAction\\";

        String sourceFolder = "C:\\Users\\Tom.fu\\Desktop\\ucf101subset\\subset\\";

        String filePrefix = "frame";
        int totalLineCnt = 0;
        int overallIndex = 0;
        try {
            BufferedReader br = new BufferedReader(new FileReader(inputFileName));
            String rdLine = null;

            while ((rdLine = br.readLine())!=null){
                totalLineCnt ++;
//                String[] rdLineSplit = rdLine.split("_");
//                String folderName = rdLineSplit[1];

//                String src = sourceFolder + folderName + "\\" + rdLine + "\\";
                String src = sourceFolder + rdLine + "\\";
                File folder = new File(src);
                File[] listOfFiles = folder.listFiles();
                int start_frame_no = 1;
                int end_frame_no = listOfFiles.length;

                for (int i = start_frame_no; i <= end_frame_no; i ++){
                    String orgFileName = src + String.format("%s%06d.jpg", filePrefix, i);
                    int targetFileNameIndex = i + overallIndex;
                    String targetFileName = outputFolder + String.format("%s%06d.jpg", filePrefix, targetFileNameIndex);

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
