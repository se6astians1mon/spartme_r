package com.example.magnetometerrawdata;

import android.util.Log;

import org.apache.commons.math3.linear.Array2DRowRealMatrix;
import org.apache.commons.math3.linear.EigenDecomposition;
import org.apache.commons.math3.linear.LUDecomposition;
import org.apache.commons.math3.linear.RealMatrix;
import org.apache.commons.math3.linear.RealVector;

import java.util.Queue;

public class PCAEllipse {
    public static RealVector findMajorAxis(double[][] data) {
        int n = data.length;
        double[] mean = calculateMean(data);
        RealMatrix covarianceMatrix = getRealMatrix(data, mean, n);

        // Check for degenerate data
        if (isDegenerateData(covarianceMatrix)) {
            Log.i("matrix","Degenerate data - no ellipse found.");
            return null;  // Handle the degenerate case
        }

        // Perform eigenvalue decomposition
        EigenDecomposition eigenDecomposition = new EigenDecomposition(covarianceMatrix);
        RealVector majorAxis = eigenDecomposition.getEigenvector(0);  // Largest eigenvalue's eigenvector
        return majorAxis;
    }
    public static double calculateMajorAxisLength(double[][] data) {
        // Eigen decomposition
        int n = data.length;
        double[] mean = calculateMean(data);
        RealMatrix covarianceMatrix = getRealMatrix(data, mean, n);
        EigenDecomposition eigenDecomposition = new EigenDecomposition(covarianceMatrix);

        // Get the eigenvalues (these are already sorted in decreasing order)
        double[] eigenvalues = eigenDecomposition.getRealEigenvalues();

        // The largest eigenvalue corresponds to the major axis
        double majorEigenvalue = eigenvalues[0]; // Assuming eigenvalues[0] is the largest

        // Major axis length is 2 times the square root of the major eigenvalue
        double majorAxisLength = 2 * Math.sqrt(majorEigenvalue);

        return majorAxisLength;
    }
    private static RealMatrix getRealMatrix(double[][] data, double[] mean, int n) {
        RealMatrix covarianceMatrix = new Array2DRowRealMatrix(2, 2);  // 2x2 matrix for covariance

        // Calculate covariance matrix
        for (double[] point : data) {
            double[] centered = {point[0] - mean[0], point[1] - mean[1]};

            // Create a 2x1 matrix (column vector) for the centered point
            RealMatrix centeredMatrix = new Array2DRowRealMatrix(new double[][]{{centered[0]}, {centered[1]}});

            // Multiply the column vector by its transpose to get a 2x2 matrix
            RealMatrix product = centeredMatrix.multiply(centeredMatrix.transpose());

            // Add the product to the covariance matrix
            covarianceMatrix = covarianceMatrix.add(product);
        }

        covarianceMatrix = covarianceMatrix.scalarMultiply(1.0 / n);
        return covarianceMatrix;
    }

    public static double[][] queueTo2DArray(Queue<float[]> dataQueue, int[] start_end) {
        int size = dataQueue.size();
        double[][] points = new double[start_end[1]-start_end[0]][2];  // 2 columns for x and y

        int index = 0;
        for (float[] point : dataQueue) {
            if(index>=start_end[0] && index<start_end[1]) {
                points[index-start_end[0]][0] = point[0];  // x value
                points[index-start_end[0]][1] = point[1];  // y value
            }
            index++;
        }

        return points;  // Returns 2D array of (x, y) pairs
    }

    public static double[] calculateMean(double[][] data) {
        double[] mean = new double[2];
        for (double[] point : data) {
            mean[0] += point[0];
            mean[1] += point[1];
        }
        mean[0] /= data.length;
        mean[1] /= data.length;
        return mean;
    }

    private static boolean isDegenerateData(RealMatrix covarianceMatrix) {
        // Use LUDecomposition to compute the determinant
        LUDecomposition luDecomposition = new LUDecomposition(covarianceMatrix);
        double determinant = luDecomposition.getDeterminant();

        if (Math.abs(determinant) < 1e-10) {
            // The determinant is very close to 0, indicating singularity or near-singularity
            return true;
        }

        // Optional: You can also check the variances of x and y
        double varianceX = covarianceMatrix.getEntry(0, 0);  // Variance in x direction
        double varianceY = covarianceMatrix.getEntry(1, 1);  // Variance in y direction

        if (varianceX < 1e-10 || varianceY < 1e-10) {
            return true;
        }

        return false;
    }


    public static double[] findMajor(double[][] points) {

        RealVector majorAxis = findMajorAxis(points);
        double[] vector = {0.0, 0.0};
        if(majorAxis !=null){
            vector = majorAxis.toArray();
            //Log.i("ellipse", "found vector! "+majorAxis);
        }

        return vector;
    }

    public static double[] getCenter(double[][] points) {
        double sumX = 0, sumY = 0;
        for (double[]point: points
             ) {
            sumX+=point[0];
            sumY+=point[1];
        }
        return (new double[]{sumX / points.length, sumY / points.length});

    }
}
