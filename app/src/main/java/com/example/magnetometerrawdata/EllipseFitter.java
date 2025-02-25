package com.example.magnetometerrawdata;

import org.ejml.data.Complex_F64;
import org.ejml.data.DMatrixRMaj;
import org.ejml.dense.row.CommonOps_DDRM;
import org.ejml.dense.row.factory.DecompositionFactory_DDRM;
import org.ejml.interfaces.decomposition.EigenDecomposition;
import org.ejml.dense.row.decomposition.eig.SwitchingEigenDecomposition_DDRM;

import java.util.Arrays;

public class EllipseFitter {

    public static double[] fitEllipse(double[] x, double[] y) {
        /**
         * Fit the coefficients a,b,c,d,e,f, representing an ellipse described by
         * the formula F(x,y) = ax^2 + bxy + cy^2 + dx + ey + f = 0 to the provided
         * arrays of data points x=[x1, x2, ..., xn] and y=[y1, y2, ..., yn].
         * Based on the algorithm of Halir and Flusser, "Numerically stable direct
         * least squares fitting of ellipses'.
         */

        int n = x.length;
        if (n != y.length || n < 5) {
            throw new IllegalArgumentException("Input arrays x and y must have the same length and at least 5 points.");
        }

        DMatrixRMaj D1 = new DMatrixRMaj(n, 3);
        DMatrixRMaj D2 = new DMatrixRMaj(n, 3);

        for (int i = 0; i < n; ++i) {
            D1.set(i, 0, x[i] * x[i]);
            D1.set(i, 1, x[i] * y[i]);
            D1.set(i, 2, y[i] * y[i]);
            D2.set(i, 0, x[i]);
            D2.set(i, 1, y[i]);
            D2.set(i, 2, 1.0);
        }

        DMatrixRMaj D1_T = CommonOps_DDRM.transpose(D1, null);
        DMatrixRMaj D2_T = CommonOps_DDRM.transpose(D2, null);

        DMatrixRMaj S1 = new DMatrixRMaj(3, 3);
        CommonOps_DDRM.mult(D1_T, D1, S1);

        DMatrixRMaj S2 = new DMatrixRMaj(3, 3); //Actually 3x3 but will only use 3x3 part from mult(D1_T, D2, S2); which returns 3x3 part
        CommonOps_DDRM.mult(D1_T, D2, S2);

        DMatrixRMaj S3 = new DMatrixRMaj(3, 3);
        CommonOps_DDRM.mult(D2_T, D2, S3);


        DMatrixRMaj S3_inv = new DMatrixRMaj(3, 3);
        if (!CommonOps_DDRM.invert(S3, S3_inv)) {
            throw new ArithmeticException("S3 is singular or non-invertible."); // Handle non-invertible case
        }


        DMatrixRMaj S2_T = CommonOps_DDRM.transpose(S2, null);
        DMatrixRMaj T = new DMatrixRMaj(3, 3); // Actually 3x3, but will only use 3x3 part from mult
        CommonOps_DDRM.mult(S3_inv, S2_T, T);
        CommonOps_DDRM.scale(-1, T); // T = -S3_inv * S2_T


        DMatrixRMaj M_term = new DMatrixRMaj(3, 3);
        CommonOps_DDRM.mult(S2, T, M_term);

        DMatrixRMaj M = new DMatrixRMaj(3, 3);
        CommonOps_DDRM.add(S1, M_term, M); // M = S1 + S2 * T


        DMatrixRMaj C_inv = new DMatrixRMaj(3, 3);
        C_inv.set(0, 2, 0.5);
        C_inv.set(1, 1, -1.0); // Already -1, no need to scale
        C_inv.set(2, 0, 0.5);


        DMatrixRMaj temp = new DMatrixRMaj(M.getNumRows(), M.getNumCols());
        CommonOps_DDRM.mult(C_inv, M, temp);
        System.arraycopy(temp.data, 0, M.data, 0, M.getNumElements());


        EigenDecomposition<DMatrixRMaj> eigenDecomp = DecompositionFactory_DDRM.eig(3, true);
        eigenDecomp.decompose(M);
        SwitchingEigenDecomposition_DDRM eigenDecompConcrete = (SwitchingEigenDecomposition_DDRM) eigenDecomp;
        double[] eigval = new double[3];
        DMatrixRMaj[] eigvecs = new DMatrixRMaj[3];
        for (int i = 0; i < eigenDecompConcrete.getNumberOfEigenvalues(); i++) {
            Complex_F64 eigenvalue = eigenDecompConcrete.getEigenvalue(i);
            eigval[i] = eigenvalue.getReal();
            eigvecs[i] = eigenDecomp.getEigenVector(i);
        }

        double[] con = new double[3];
        for(int i=0; i<3; ++i) {
            con[i] = 4 * eigvecs[i].get(0) * eigvecs[i].get(2) - Math.pow(eigvecs[i].get(1), 2);
        }


        int bestIndex = -1;
        for(int i=0; i<3; ++i){
            if(con[i] > 0){
                bestIndex = i;
                break; // Assuming only one eigenvector will satisfy the condition.
            }
        }

        if (bestIndex == -1) {
            throw new ArithmeticException("No eigenvector found corresponding to a valid ellipse.");
        }


        DMatrixRMaj ak_matrix = eigvecs[bestIndex];

        DMatrixRMaj ak = new DMatrixRMaj(3,1); // Column vector
        ak.set(0,0, ak_matrix.get(0,0));
        ak.set(1,0, ak_matrix.get(1,0));
        ak.set(2,0, ak_matrix.get(2,0));


        DMatrixRMaj Tak = new DMatrixRMaj(3, 1);
        CommonOps_DDRM.mult(T, ak, Tak);


        double[] coeffs = new double[6];
        coeffs[0] = ak.get(0); // a
        coeffs[1] = ak.get(1); // b
        coeffs[2] = ak.get(2); // c
        coeffs[3] = Tak.get(0); // d
        coeffs[4] = Tak.get(1); // e
        coeffs[5] = Tak.get(2); // f


        return coeffs;
    }


    public static double[] cartToPol(double[] coeffs) {
        /**
         * Convert the cartesian conic coefficients, (a, b, c, d, e, f), to the
         * ellipse parameters, where F(x, y) = ax^2 + bxy + cy^2 + dx + ey + f = 0.
         * The returned parameters are x0, y0, ap, bp, e, phi, where (x0, y0) is the
         * ellipse centre; (ap, bp) are the semi-major and semi-minor axes,
         * respectively; e is the eccentricity; and phi is the rotation of the semi-
         * major axis from the x-axis.
         */
        double a = coeffs[0];
        double b = coeffs[1] / 2.0;
        double c = coeffs[2];
        double d = coeffs[3] / 2.0;
        double f = coeffs[4] / 2.0;
        double g = coeffs[5];

        double den = b * b - a * c;
        if (den > 0) {
            throw new IllegalArgumentException("coeffs do not represent an ellipse: b^2 - ac must be negative!"); // Original was b^2 - 4ac, but here it's b^2 - ac as b is already b/2
        }

        // The location of the ellipse centre.
        double x0 = (c * d - b * f) / den;
        double y0 = (a * f - b * d) / den;

        double num = 2 * (a * f * f + c * d * d + g * b * b - 2 * b * d * f - a * c * g);
        double fac = Math.sqrt(Math.pow(a - c, 2) + 4 * b * b);
        // The semi-major and semi-minor axis lengths (these are not sorted).
        double ap = Math.sqrt(num / den / (fac - a - c));
        double bp = Math.sqrt(num / den / (-fac - a - c));

        // Sort the semi-major and semi-minor axis lengths but keep track of
        // the original relative magnitudes of width and height.
        boolean width_gt_height = true;
        if (ap < bp) {
            width_gt_height = false;
            double temp = ap;
            ap = bp;
            bp = temp;
        }

        // The eccentricity.
        double r = (bp / ap) * (bp / ap);
        if (r > 1) {
            r = 1 / r;
        }
        double e = Math.sqrt(1 - r);

        // The angle of anticlockwise rotation of the major-axis from x-axis.
        double phi;
        if (b == 0) {
            phi = (a < c) ? 0 : Math.PI / 2.0;
        } else {
            phi = Math.atan((2.0 * b) / (a - c)) / 2.0;
            if (a > c) {
                phi += Math.PI / 2.0;
            }
        }
        if (!width_gt_height) {
            // Ensure that phi is the angle to rotate to the semi-major axis.
            phi += Math.PI / 2.0;
        }
        phi = phi % Math.PI;

        return new double[]{x0, y0, ap, bp, e, phi};
    }

    public static double[][] getEllipsePts(double[] params, int npts, double tmin, double tmax) {
        /**
         * Return npts points on the ellipse described by the params = x0, y0, ap,
         * bp, e, phi for values of the parametric variable t between tmin and tmax.
         */
        double x0 = params[0];
        double y0 = params[1];
        double ap = params[2];
        double bp = params[3];
        double phi = params[5]; //e is params[4], but not used here

        double[] t = new double[npts];
        for(int i=0; i<npts; ++i) {
            t[i] = tmin + (tmax - tmin) * i / (npts - 1.0); //Equivalent of linspace
        }


        double[] x = new double[npts];
        double[] y = new double[npts];

        for (int i = 0; i < npts; ++i) {
            x[i] = x0 + ap * Math.cos(t[i]) * Math.cos(phi) - bp * Math.sin(t[i]) * Math.sin(phi);
            y[i] = y0 + ap * Math.cos(t[i]) * Math.sin(phi) + bp * Math.sin(t[i]) * Math.cos(phi);
        }

        return new double[][]{x, y};
    }

    public static void main(String[] args) {
        int npts = 250;
        double tmin = Math.PI / 6.0;
        double tmax = 4.0 * Math.PI / 3.0;
        double x0_exact = 4;
        double y0_exact = -3.5;
        double ap_exact = 7;
        double bp_exact = 3;
        double phi_exact = Math.PI / 4.0;

        double[][] ellipse_pts = getEllipsePts(new double[]{x0_exact, y0_exact, ap_exact, bp_exact, 0, phi_exact}, npts, tmin, tmax);
        double[] x = ellipse_pts[0];
        double[] y = ellipse_pts[1];

        double noise = 0.1;
        java.util.Random random = new java.util.Random();
        for (int i = 0; i < npts; ++i) {
            x[i] += noise * random.nextGaussian();
            y[i] += noise * random.nextGaussian();
        }

        double[] coeffs = fitEllipse(x, y);
        System.out.println("Exact parameters:");
        System.out.println("x0, y0, ap, bp, phi = " + x0_exact + ", " + y0_exact + ", " + ap_exact + ", " + bp_exact + ", " + phi_exact);
        System.out.println("Fitted parameters (coeffs a, b, c, d, e, f):");
        System.out.println("a, b, c, d, e, f = " + Arrays.toString(coeffs));


        double[] fittedParams = cartToPol(coeffs);
        System.out.println("Fitted parameters (x0, y0, ap, bp, e, phi):");
        System.out.println("x0, y0, ap, bp, e, phi = " + fittedParams[0] + ", " + fittedParams[1] + ", " + fittedParams[2] + ", " + fittedParams[3] + ", " + fittedParams[4] + ", " + fittedParams[5]);
    }
}