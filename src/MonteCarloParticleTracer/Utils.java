package MonteCarloParticleTracer;

import org.apache.commons.math3.analysis.interpolation.LinearInterpolator;
import org.apache.commons.math3.analysis.polynomials.PolynomialSplineFunction;
import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.io.File;
import java.util.ArrayList;
import java.util.Scanner;

/**
 * Created by lahmann on 2016-06-15.
 */
public class Utils {

    private static final NormalDistribution NORMAL_DISTRIBUTION = new NormalDistribution();

    static double[] removeElement(double[] array, int index){

        double[] newArray = new double[array.length - 1];
        for (int i = 0; i < array.length; i++){

            if (i < index){
                newArray[i] = array[i];
            }

            else if (i > index){
                newArray[i-1] = array[i];
            }
        }
        return newArray;
    }





    public static double[][] parseCSV(File file) throws Exception {

        ArrayList<ArrayList<Double>> list = new ArrayList<>();

        Scanner s = new Scanner(file);
        while (s.hasNextLine()) {

            String[] entries = s.nextLine().split(",");
            while (list.size() < entries.length) {
                list.add(new ArrayList<>());
            }

            for (int i = 0; i < entries.length; i++) {
                list.get(i).add(Double.parseDouble(entries[i]));
            }

        }
        s.close();

        double[][] data = new double[list.size()][list.get(0).size()];
        for (int i = 0; i < data.length; i++) {
            for (int j = 0; j < data[i].length; j++) {
                data[i][j] = list.get(i).get(j);
            }
        }

        return data;

    }

    static double getNormalizedRadius(Particle particle, Plasma plasma){
        double[] coordinates = Utils.getSphericalFromVector(particle.getPosition());

        double R     = coordinates[0];
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double rMax = plasma.getOuterRadiusBound(theta, phi);
        double rMin = plasma.getInnerRadiusBound(theta, phi);
        return (R - rMin) / (rMax - rMin);
    }

    static double[] getSphericalFromVector(Vector3D vector){
        SphericalCoordinates coordinates = new SphericalCoordinates(vector);

        // Apache is in Math notation whilst we're in Physics notation
        double r     = coordinates.getR();
        double theta = coordinates.getPhi();
        double phi   = coordinates.getTheta();

        if (Double.isNaN(theta)){
            theta = 0.0;
        }

        if (Double.isNaN(phi)){
            phi = 0.0;
        }

        return new double[] {r, theta, phi};
    }

    public static Vector3D getVectorFromSpherical(double r, double theta, double phi){
        // Apache is in Math notation whilst we're in Physics notation
        return new SphericalCoordinates(r, phi, theta).getCartesian();
    }

    static Vector3D sampleRandomNormalizedVector(){
        double[] xR = NORMAL_DISTRIBUTION.sample(3);
        return new Vector3D(xR[0], xR[1], xR[2]).normalize();
    }

}
