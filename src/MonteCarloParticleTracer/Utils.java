package MonteCarloParticleTracer;

import org.apache.commons.math3.distribution.NormalDistribution;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Created by lahmann on 2016-06-15.
 */
public class Utils {



    public static final double MINIMUM_LI_PETRASSO_ENERGY_MeV = 0.01;


    public static final String D3He_ENDF_XS_FILE = "./data/d3He_BoschHale.dat";
    public static final String DT_ENDF_XS_FILE   = "./data/dT_BoschHale.dat";

    public static final String DDp_REACTIVITY_FILE = "./data/DDp_Reactivities.dat";
    public static final String DDn_REACTIVITY_FILE = "./data/DDn_Reactivities.dat";
    public static final String D3Hep_REACTIVITY_FILE = "./data/D3Hep_Reactivities.dat";


    public static double[] linspace(double a, double b, int N){
        double[] nodes = new double[N];

        for (int i = 0; i < nodes.length; i++){
            nodes[i] = a + (b-a)*i/(nodes.length - 1);
        }

        return nodes;
    }

    public static double[] logspace(double a, double b, int N){
        double[] nodes = new double[N];

        for (int i = 0; i < nodes.length; i++){
            double power =  a + (b-a)*i/(nodes.length - 1);
            nodes[i] = Math.pow(10, power);
        }

        return nodes;
    }

    public static double getNormalizedRadius(Particle particle, PlasmaLayer plasmaLayer){
        double[] coordinates = Utils.getSphericalFromVector(particle.getPosition());

        double R     = coordinates[0];
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double rMax = plasmaLayer.getOuterRadiusBound(theta, phi);
        double rMin = plasmaLayer.getInnerRadiusBound(theta, phi);
        return (R - rMin) / (rMax - rMin);
    }

    public static double[] getSphericalFromVector(Vector3D vector){
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

    public static Vector3D sampleRandomNormalizedVector(){
        double[] xR = new NormalDistribution().sample(3);
        return new Vector3D(xR[0], xR[1], xR[2]).normalize();
    }

}
