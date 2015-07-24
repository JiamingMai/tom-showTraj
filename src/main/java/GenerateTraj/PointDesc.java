package GenerateTraj;

import org.bytedeco.javacpp.opencv_core;

/**
 * Created by Tom.fu on 2/12/2014.
 */
public class PointDesc {
    public float[] mbhX;
    public float[] mbhY;
    public float[] hog;
    public opencv_core.CvPoint2D32f point;

    PointDesc(DescInfo mbhInfo, DescInfo hogInfo, opencv_core.CvPoint2D32f point){
        mbhX = new float[mbhInfo.nxCells * mbhInfo.nyCells * mbhInfo.nBins];
        mbhY = new float[mbhInfo.nxCells * mbhInfo.nyCells * mbhInfo.nBins];
        hog = new float[hogInfo.nxCells * hogInfo.nyCells * hogInfo.nBins];
        this.point = new opencv_core.CvPoint2D32f(point);
        this.point.put(point);
    }
}
