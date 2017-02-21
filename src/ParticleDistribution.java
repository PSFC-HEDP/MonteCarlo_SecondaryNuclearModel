import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Created by lahmann on 2016-06-15.
 */
public class ParticleDistribution {

    private int A;      // Mass number of Particles in this Distribution
    private int Z;      // Atomic number of Particles in the Distribution

    private Distribution radialDistribution;        // Normalized radius distribution [0,1]
    private Distribution energyDistribution;        // Energy Distribution in MeV



    /**
     * Built in pre-defined particle distributions
     */

    public static ParticleDistribution ThermalDDtDistribution(Distribution radialDistribution, double temperature){

        temperature *= 1e-3;                                                            // keV -> MeV

        double Q = (Utils.DD_T_birthEnergy_MeV + Utils.DD_p_birthEnergy_MeV);           // Q values in keV
        double mu = Utils.DD_T_birthEnergy_MeV;                                         // Mean birth energy (MeV)
        double factor = 2.0 * (3.0*1.0) / (3.0+1.0) / (2.0+2.0);                        // Factor = 2 (mC*mD) / (mC+mD) / (mA+mB)
        double sigma = Math.sqrt(factor * Q * temperature);                             // Spectral width (keV)

        Distribution energyDistribution = Distribution.normDistribution(mu, sigma);
        return new ParticleDistribution(1, 3, radialDistribution, energyDistribution);
    }

    public static ParticleDistribution ThermalDD3HeDistribution(Distribution radialDistribution, double temperature){

        temperature *= 1e-3;                                                            // keV -> MeV

        double Q = (Utils.DD_3He_birthEnergy_MeV + Utils.DD_n_birthEnergy_MeV);         // Q values in keV
        double mu = Utils.DD_3He_birthEnergy_MeV;                                       // Mean birth energy (MeV)
        double factor = 2.0 * (3.0*1.0) / (3.0+1.0) / (2.0+2.0);                        // Factor = 2 (mC*mD) / (mC+mD) / (mA+mB)
        double sigma = Math.sqrt(factor * Q * temperature);                             // Spectral width (keV)

        Distribution energyDistribution = Distribution.normDistribution(mu, sigma);
        return new ParticleDistribution(2, 3, radialDistribution, energyDistribution);
    }


    /**
     * Constructors
     */

    public ParticleDistribution(int Z, int A,  Distribution radialDistribution, Distribution energyDistribution) {
        this.Z = Z;
        this.A = A;

        this.radialDistribution = radialDistribution;
        this.energyDistribution = energyDistribution;
    }


    /**
     * Sampling function
     * Depends on the boundaries of the plasma since this Distribution is defined in terms of normalized radii [0,1]
     */

    public Particle sample(Plasma plasma){

        Particle particle = new Particle(energyDistribution.sample(), this.Z, this.A);

        // Sample a random position vector (we'll ignore the magnitude)
        double dx = Math.random() - 0.5;
        double dy = Math.random() - 0.5;
        double dz = Math.random() - 0.5;

        Vector3D position = new Vector3D(dx, dy, dz).normalize();

        // Convert it to spherical
        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        // Use our radial distribution and the plasma bounds to determine the magnitude of
        // our position vector
        double r = radialDistribution.sample() * plasma.getRadiusBound(theta, phi);

        // Rescale our position vector to this magnitude
        position = position.scalarMultiply(r);
        particle.setPosition(position);

        // Sample a random direction vector
        dx = Math.random() - 0.5;
        dy = Math.random() - 0.5;
        dz = Math.random() - 0.5;

        Vector3D direction = new Vector3D(dx, dy, dz).normalize();
        particle.setDirection(direction);

        return particle;
    }

    public int getA() {
        return A;
    }

    public int getZ() {
        return Z;
    }

    public double getMaxEnergy(){
        double[] energyNodes = energyDistribution.getNodes();
        return energyNodes[energyNodes.length - 1];
    }

}
