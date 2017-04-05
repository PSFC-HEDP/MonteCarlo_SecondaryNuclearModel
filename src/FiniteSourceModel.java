import MonteCarloParticleTracer.Particle;
import MonteCarloParticleTracer.ParticleDistribution;
import MonteCarloParticleTracer.PlasmaLayer;
import MonteCarloParticleTracer.Utils;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolatingFunction;
import org.apache.commons.math3.analysis.interpolation.BicubicInterpolator;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;


/**
 * Created by lahmann on 2017-02-27.
 */
public class FiniteSourceModel {

    // MonteCarloParticleTracer.Constants used when generating the stopping power function
    final double MINIMUM_STOPPING_POWER_ENERGY = Utils.MINIMUM_LI_PETRASSO_ENERGY_MeV;       // MeV
    final int NUM_ENERGY_NODES = 200;


    private final int NUM_STEPS_PER_SYSTEM_SIZE = 200;
    static final double ENERGY_CUTOFF = 0.01;                   // MeV
    private final double ACCEPTABLE_ENERGY_ERROR = 0.01;        // Fraction


    private ParticleDistribution testParticleDistribution;
    private BicubicInterpolatingFunction stoppingPower;     // Stopping power profile defined in f(E,r) space (MeV / cm)
    private PlasmaLayer plasma;


    private Vector3D diagnosticLocation;


    public FiniteSourceModel(ParticleDistribution testParticleDistribution, PlasmaLayer plasma, Vector3D diagnosticLocation) {
        this.testParticleDistribution = testParticleDistribution;
        this.plasma = plasma;
        this.diagnosticLocation = diagnosticLocation;
        this.stoppingPower = calculateStoppingPowerFunction();
    }


    public Particle[] rangeParticles(int N){

        Particle[] particles = new Particle[N];
        for (int i = 0; i < particles.length; i++){

            Particle particle = testParticleDistribution.sample(plasma);

            //Vector3D direction = particle.getPosition().normalize();
            //Vector3D direction = diagnosticLocation.subtract(particle.getPosition()).normalize();
            //particle.setDirection(direction);

            // Determining a step size is non-trivial since our plasma could be any shape with any profiles
            // A cheap calculation we can do is to take the plasmaBound R(theta = 0, phi = 0) to be somewhat
            // representative of the system size and take our step size to be some fraction of that.
            // This approximation SHOULD hold for plasmas with relatively small low mode asymmetries
            // High amplitude and/or high mode asymmetries have the potential to break this
            double plasmaBound = plasma.getRadiusBound(0.0, 0.0);
            double dx = plasmaBound / NUM_STEPS_PER_SYSTEM_SIZE;


            double energyToLose = particle.getEnergy() - ENERGY_CUTOFF;


            // While the particle is inside this plasma
            while (plasma.getIsInside(particle.getPosition()) && particle.getEnergy() >  ENERGY_CUTOFF) {

                // Parameters at the starting point
                double r1 = getNormalizedRadius(particle);
                double E1 = particle.getEnergy();
                double dEdx1 = stoppingPower.value(E1, r1);

                // If the particle is losing a significant amount of energy, we need to take smaller steps
                dx = Math.min(dx, -energyToLose / dEdx1 / NUM_STEPS_PER_SYSTEM_SIZE);

                // Attempt to step and verify we're still in the plasma
                Particle steppedParticle = particle.step(dx, dEdx1);
                if (!plasma.getIsInside(steppedParticle.getPosition()) || steppedParticle.getEnergy() < ENERGY_CUTOFF){
                    break;
                }

                // Parameters at the destination point
                double r2 = getNormalizedRadius(steppedParticle);
                double E2 = steppedParticle.getEnergy();
                double dEdx2 = stoppingPower.value(E2, r2);

                // Use trapezoidal method to iterate onto the particle's final energy
                double energyError = Double.MAX_VALUE;
                while (energyError > ACCEPTABLE_ENERGY_ERROR){
                    double averageStopping = 0.5*(dEdx1 + dEdx2);
                    steppedParticle = particle.step(dx, averageStopping);

                    double newEnergy = steppedParticle.getEnergy();
                    if (newEnergy < ENERGY_CUTOFF)  break;

                    energyError = FastMath.abs(newEnergy - E2)/E2;
                    E2 = newEnergy;
                    dEdx2 = stoppingPower.value(E2, r2);
                }

                particle = steppedParticle;
            }

            particles[i] = particle;
        }

        return particles;
    }

    private BicubicInterpolatingFunction calculateStoppingPowerFunction(){

        final double maxEnergy = testParticleDistribution.getMaxEnergy();
        double[] energyNodes = Utils.linspace(MINIMUM_STOPPING_POWER_ENERGY, maxEnergy, NUM_ENERGY_NODES);
        double[] radiusNodes = plasma.getRadiusNodes();

        double[][] dEdx = new double[energyNodes.length][radiusNodes.length];
        for (int i = 0; i < energyNodes.length; i++){
            //MonteCarloParticleTracer.Particle p = new MonteCarloParticleTracer.Particle(energyNodes[i], testParticleDistribution.getZ(), testParticleDistribution.getA());
            //dEdx[i] = plasma.getStoppingPower(p);
        }

        return new BicubicInterpolator().interpolate(energyNodes, radiusNodes, dEdx);
    }

    private double getNormalizedRadius(Particle particle){
        Vector3D position = particle.getPosition();

        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        // Position and energy at our starting point
        return position.getNorm() / plasma.getRadiusBound(theta, phi);
    }
}
