package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

/**
 * Created by lahmann on 2016-06-15.
 */
public class ParticleDistribution {

    private double particleWeight = 1.0;
    private ParticleType type;                      // MonteCarloParticleTracer.Particle type in this MonteCarloParticleTracer.Distribution
    private Distribution radialDistribution;        // Normalized radius distribution [0,1]
    private Distribution energyDistribution;        // Energy MonteCarloParticleTracer.Distribution in MeV



    /**
     * Built in pre-defined particle distributions
     */

    // D + D -> T + p
    public static ParticleDistribution ThermalDDtDistribution(Distribution radialDistribution, double temperature){

        temperature *= 1e-3;                                                            // keV -> MeV

        double factor = 2.0;                                                            // Factor = 2 (mC*mD) / (mC+mD) / (mA+mB)
        factor *= ParticleType.proton.getMass();
        factor *= ParticleType.triton.getMass();
        factor /= (ParticleType.proton.getMass() + ParticleType.triton.getMass());
        factor /= (ParticleType.deuteron.getMass() + ParticleType.deuteron.getMass());

        double Q = Constants.DD_P_ENERGY_RELEASE;                                         // Q values in MeV
        double mu = Constants.DD_T_BIRTH_ENERGY_MEV;                                   // Mean birth energy (MeV)
        double sigma = Math.sqrt(factor * Q * temperature);                             // Spectral width (MeV)

        Distribution energyDistribution = Distribution.normDistribution(mu, sigma);
        return new ParticleDistribution(ParticleType.triton, radialDistribution, energyDistribution);
    }

    // D + D -> 3He + n
    public static ParticleDistribution ThermalDD3HeDistribution(Distribution radialDistribution, double temperature){

        temperature *= 1e-3;                                                            // keV -> MeV

        double factor = 2.0;                                                            // Factor = 2 (mC*mD) / (mC+mD) / (mA+mB)
        factor *= ParticleType.neutron.getMass();
        factor *= ParticleType.helium3.getMass();
        factor /= (ParticleType.neutron.getMass() + ParticleType.helium3.getMass());
        factor /= (ParticleType.deuteron.getMass() + ParticleType.deuteron.getMass());

        double Q = Constants.DD_N_ENERGY_RELEASE;                                         // Q values in MeV
        double mu = Constants.DD_3HE_BIRTH_ENERGY_MEV;                                   // Mean birth energy (MeV)
        double sigma = Math.sqrt(factor * Q * temperature);                             // Spectral width (MeV)

        Distribution energyDistribution = Distribution.normDistribution(mu, sigma);
        return new ParticleDistribution(ParticleType.helium3, radialDistribution, energyDistribution);
    }




    /**
     * Constructors
     */

    public ParticleDistribution(int Z, double mass,  Distribution radialDistribution, Distribution energyDistribution) {
        this(new ParticleType(Z, mass), radialDistribution, energyDistribution);
    }

    public ParticleDistribution(ParticleType type, Distribution radialDistribution, Distribution energyDistribution) {
        this.type = type;
        this.radialDistribution = radialDistribution;
        this.energyDistribution = energyDistribution;
    }


    /**
     * Sampling function
     * Depends on the boundaries of the plasma since this Distribution is defined in terms of normalized radii [0,1]
     */

    public Particle sample(PlasmaLayer plasmaLayer){

        // Sample the position
        Vector3D position = Utils.sampleRandomNormalizedVector();

        // Convert it to spherical
        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        // Use our radial distribution and the plasma bounds to determine the magnitude of
        // our position vector
        double rNorm = radialDistribution.sample();
        double rMax  = plasmaLayer.getOuterRadiusBound(theta, phi);
        double rMin  = plasmaLayer.getInnerRadiusBound(theta, phi);

        double r = (rMax - rMin) * rNorm + rMin;

        // Rescale our position vector to this magnitude
        position = position.scalarMultiply(r);


        // Sample the direction
        Vector3D direction = Utils.sampleRandomNormalizedVector();


        // Sample the energy
        double energy = energyDistribution.sample();

        // Return the MonteCarloParticleTracer.Particle Object

        Particle particle = new Particle(getZ(), getMass(), position, direction, energy, 0.0);
        particle.setWeight(particleWeight);
        return particle;
    }

    public ParticleType getType() {
        return type;
    }

    public double getMass() {
        return type.getMass();
    }

    public int getZ() {
        return type.getZ();
    }

    public double getMaxEnergy(){
        double[] energyNodes = energyDistribution.getNodes();
        return energyNodes[energyNodes.length - 1];
    }

    public void setParticleWeight(double particleWeight) {
        this.particleWeight = particleWeight;
    }

    public String toString(){
        String string = "";
        string += "Particle type: " + type.toString() + "\n";
        string += "         Radial Distribution        \n";
        string += "------------------------------------\n";
        string += radialDistribution.toString();
        string += "\n";
        string += "         Energy Distribution        \n";
        string += "------------------------------------\n";
        string += energyDistribution.toString();
        return string;
    }



}
