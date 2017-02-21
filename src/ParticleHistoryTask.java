import org.apache.commons.math3.analysis.interpolation.BicubicInterpolatingFunction;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

import java.util.concurrent.ExecutionException;
import java.util.concurrent.RunnableFuture;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Created by lahmann on 2016-09-10.
 */
public class ParticleHistoryTask implements RunnableFuture {

    private CrossSection crossSection;
    private BicubicInterpolatingFunction stoppingPower;     // Stopping power profile defined in f(E,r) space (MeV / cm)
    private ParticleDistribution particleDistribution;
    private Plasma plasma;
    private int numParticles;

    private boolean debugMode = false;

    private Double reactionProbability = 0.0;

    /**
     * Step parameters
     */
    private final int NUM_STEPS_PER_SYSTEM_SIZE = 200;

    static final double ENERGY_CUTOFF = 0.01;                   // MeV
    private final double ACCEPTABLE_ENERGY_ERROR = 0.01;        // Fraction



    public ParticleHistoryTask(CrossSection crossSection, ParticleDistribution particleDistribution, Plasma plasma,
                               BicubicInterpolatingFunction stoppingPower, int numParticles) {
        this.crossSection = crossSection;
        this.particleDistribution = particleDistribution;
        this.plasma = plasma;
        this.stoppingPower = stoppingPower;
        this.numParticles = numParticles;
    }

    public void setDebugMode(boolean debugMode) {
        this.debugMode = debugMode;
    }

    @Override
    public void run() {

        // Determining a step size is non-trivial since our plasma could be any shape with any profiles
        // A cheap calculation we can do is to take the plasmaBound R(theta = 0, phi = 0) to be somewhat
        // representative of the system size and take our step size to be some fraction of that.
        // This approximation SHOULD hold for plasmas with relatively small low mode asymmetries
        // High amplitude and/or high mode asymmetries have the potential to break this

        double plasmaBound = plasma.getRadiusBound(0.0, 0.0);
        double dx = plasmaBound / NUM_STEPS_PER_SYSTEM_SIZE;


        // In reality we should be sampling these from some thermal distribution
        // Doing this assumes that those effects would be small
        // (i.e E ~ 0 relative to the energy of the particles we're tracking)

        Particle backgroundDeuteron = Particle.deuteron(0.0);


        for (int i = 0; i < numParticles; i++) {

            int totalSteps = 0;
            double historyReactionProbability = 0.0;
            Particle particle = particleDistribution.sample(plasma);
            double energyToLose = particle.getE() - ENERGY_CUTOFF;


            // DEBUG MODE LINE
            if (debugMode){
                Vector3D position = particle.getPosition();
                Vector3D direction = particle.getDirection();
                System.out.printf("Starting history %d with:\n" +
                        "rx = %+.4e cm, ry = %+.4e cm, rz = %+.4e cm\n" +
                        "dx = %+.4e cm, dy = %+.4e cm, dz = %+.4e cm\n" +
                        "E  = %+.4e MeV\n",

                        (i+1), position.getX(), position.getY(), position.getZ(),
                        direction.getX(), direction.getY(), direction.getZ(), particle.getE());
            }


            // While the particle is inside this plasma
            while (plasma.getIsInside(particle.getPosition()) && particle.getE() >  ENERGY_CUTOFF) {

                // Parameters at the starting point
                double r1 = getNormalizedRadius(particle);
                double E1 = particle.getE();
                double dEdx1 = stoppingPower.value(E1, r1);

                // If the particle is losing a significant amount of energy, we need to take smaller steps
                dx = Math.min(dx, -energyToLose / dEdx1 / NUM_STEPS_PER_SYSTEM_SIZE);

                // Attempt to step and verify we're still in the plasma
                Particle steppedParticle = particle.step(dx, dEdx1);
                if (!plasma.getIsInside(steppedParticle.getPosition()) || steppedParticle.getE() < ENERGY_CUTOFF){
                    break;
                }

                // Parameters at the destination point
                double r2 = getNormalizedRadius(steppedParticle);
                double E2 = steppedParticle.getE();
                double dEdx2 = stoppingPower.value(E2, r2);

                // Use trapezoidal method to iterate onto the particle's final energy
                double energyError = Double.MAX_VALUE;
                while (energyError > ACCEPTABLE_ENERGY_ERROR){
                    double averageStopping = 0.5*(dEdx1 + dEdx2);
                    steppedParticle = particle.step(dx, averageStopping);

                    double newEnergy = steppedParticle.getE();
                    if (newEnergy < ENERGY_CUTOFF)  break;

                    energyError = FastMath.abs(newEnergy - E2)/E2;
                    E2 = newEnergy;
                    dEdx2 = stoppingPower.value(E2, r2);
                }

                // Grab the average deuteron number density between these two positions
                double nD1 = plasma.getDeuteronNumberDensity(particle.getPosition());
                double nD2 = plasma.getDeuteronNumberDensity(steppedParticle.getPosition());
                double nD = 0.5*(nD1 + nD2);

                // Grab the average cross section between these two positions
                double sigma1 = 1e-24 * crossSection.evaluate(particle, backgroundDeuteron);
                double sigma2 = 1e-24 * crossSection.evaluate(steppedParticle, backgroundDeuteron);
                double sigma = 0.5*(sigma1 + sigma2);

                // Calculate the reaction probability
                historyReactionProbability += nD * sigma * dx * (1 - historyReactionProbability);
                particle = steppedParticle;

                totalSteps++;
            }
            reactionProbability += historyReactionProbability;
        }
        reactionProbability /= numParticles;
    }

    /**
     * Various getters
     */

    public Double getReactionProbability() {
        return reactionProbability;
    }


    /**
     * Private convenience methods
     */

    private double getNormalizedRadius(Particle particle){
        Vector3D position = particle.getPosition();

        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        // Position and energy at our starting point
        return position.getNorm() / plasma.getRadiusBound(theta, phi);
    }


    /**
     * Inherited methods that we don't currently use
     */

    @Override
    public Object get() throws InterruptedException, ExecutionException {
        return null;
    }

    @Override
    public boolean cancel(boolean mayInterruptIfRunning) {
        return false;
    }

    @Override
    public boolean isCancelled() {
        return false;
    }

    @Override
    public boolean isDone() {
        return this.reactionProbability == null;
    }

    @Override
    public Object get(long timeout, TimeUnit unit) throws InterruptedException, ExecutionException, TimeoutException {
        return null;
    }
}
