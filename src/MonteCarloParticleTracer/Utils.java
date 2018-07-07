package MonteCarloParticleTracer;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Created by lahmann on 2016-06-15.
 */
public class Utils {

    private static final NormalDistribution NORMAL_DISTRIBUTION = new NormalDistribution();







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
