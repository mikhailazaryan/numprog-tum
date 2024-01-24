import java.util.Arrays;

public class Gauss {
    public static double[] backSubst(double[][] matrix, double[] b) {
        int n = matrix.length - 1;
        double[] solution = new double[n+1];

        solution[n] = b[n]/matrix[n][n];
        double dividend;
        for (int i = n-1; i >= 0; i--) {
            dividend = b[i];
            for (int j = n; j >= i+1; j--) {
                dividend -= matrix[i][j]*solution[j];
            }
            solution[i] = dividend/matrix[i][i];
        }
        return solution;
    }

    public static double[] solve(double[][] matrix, double[] b) {
        int n = matrix.length;
        double[][] nm = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                nm[i][j] = matrix[i][j];
            }
        }
        double[] bcopy = new double[n];
        for (int j = 0; j < n; j++) {
            bcopy[j]= b[j];
        }

        int maxLine;
        double max;

        for (int i = 0; i < n-1; i++){
            max = Math.abs(nm[i][i]);
            maxLine = i;
            for (int j = i+1; j < n; j++) {
                double cur = Math.abs(nm[j][i]);
                if (cur > max) {
                    max = cur;
                    maxLine = j;
                }
            }

            double[] temp = nm[maxLine];
            nm[maxLine] = nm[i];
            nm[i] = temp;

            double t = bcopy[maxLine];
            bcopy[maxLine] = bcopy[i];
            bcopy[i] = t;



            for (int j = i + 1; j < n; j++) {
                double fac = nm[j][i] / nm[i][i];
                for (int k = i; k < n; k++) {
                    nm[j][k] -= fac * nm[i][k];
                }
                bcopy[j] -= fac*bcopy[i];
            }
        }
        return backSubst(nm, bcopy);
    }
    public static double[] solveSing(double[][] matrix) {
        int n = matrix.length;
        double[][] nm = new double[n][n];
        for (int i = 0; i < n; i++) {
            for (int j = 0; j < n; j++) {
                nm[i][j] = matrix[i][j];
            }
        }

        int maxLine;
        double max;
        int lastColumnOfT = -1;
        for (int i = 0; i < n-1; i++){

            max = Math.abs(nm[i][i]);
            maxLine = i;
            for (int j = i+1; j < n; j++) {
                double cur = Math.abs(nm[j][i]);
                if (cur > max) {
                    max = cur;
                    maxLine = j;
                }
            }
            if (max == 0.0) lastColumnOfT = i;
            if (lastColumnOfT != -1) break;

            double[] temp = nm[maxLine];
            nm[maxLine] = nm[i];
            nm[i] = temp;


            for (int j = i + 1; j < n; j++) {
                double fac = nm[j][i] / nm[i][i];
                for (int k = i; k < n; k++) {
                    nm[j][k] -= fac * nm[i][k];
                }
            }
        }
        if (lastColumnOfT == -1) {
            double[] ans = new double[n];
            Arrays.fill(ans, 0);
            return ans;
        } else {

            double[][] m1 = new double[lastColumnOfT][lastColumnOfT];
            double[] v1 = new double[lastColumnOfT];

            for (int i = 0; i < lastColumnOfT; i++) {
                v1[i] = nm[i][lastColumnOfT] == 0? 0 : -nm[i][lastColumnOfT];
                for (int j = 0; j < lastColumnOfT; j++) {
                    m1[i][j] = nm[i][j];
                }
            }

            double[] p = new double[n];
            double[] res = backSubst(m1, v1);

            for (int i = 0; i < res.length; i++) {
                p[i] = res[i];
            }
            p[res.length] = 1.0;
            for (int i = res.length+1; i < n; i++){
                p[i] = 0.0;
            }
            return p;
        }
    }

    public static double[][] buildProbabilityMatrix(int[][] linkMatrix, double rho) {
        int n = linkMatrix.length;
        double k = 1.0-rho, s = rho/(1.0*n);
        double[][] ans = new double[n][n];
        for (int i = 0; i < n; i++){
            int sitesCount = 0;
            for (int j = 0; j < n; j++) {
                sitesCount += linkMatrix[j][i];
            }
            for (int j = 0; j < linkMatrix[i].length; j++) {
                if (linkMatrix[j][i] == 1) ans[j][i] = (1.0 / sitesCount)*k+s;
                else ans[j][i] = s;
            }
        }
        return ans;
    }

}
