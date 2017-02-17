import org.apache.commons.math3.analysis.interpolation.BicubicInterpolatingFunction;
import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

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
    private final int NUM_STEPS_PER_SYSTEM_SIZE    = 100;
    private final double ENERGY_CUTOFF = 0.01;



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

            if (debugMode)      System.out.printf("Starting history %d ... \n", (i+1));

            int totalSteps = 0;
            double historyReactionProbability = 0.0;
            Particle particle = particleDistribution.sample(plasma);

            // While the particle is inside this plasma
            while (plasma.getIsInside(particle.getPosition()) && particle.getE() >  ENERGY_CUTOFF) {

                // We need to know where we are relative to the plasma boundary
                // In order to use our stopping power function which is defined in normalized r
                Vector3D position = particle.getPosition();

                SphericalCoordinates coordinates = new SphericalCoordinates(position);

                // Apache is in Math notation whilst we're in Physics notation
                double theta = coordinates.getPhi();
                double phi = coordinates.getTheta();

                double r = position.getNorm() / plasma.getRadiusBound(theta, phi);      // Normalized r


                // Grab all of the plasma quantities at this position
                double nD = plasma.getDeuteronNumberDensity(particle.getPosition());
                double sigma = 1e-24 * crossSection.evaluate(particle, backgroundDeuteron);
                double dEdx = stoppingPower.value(particle.getE(), r);


                // Calculate the reaction probability
                historyReactionProbability += nD * sigma * dx * (1 - historyReactionProbability);
                particle.step(dx, dEdx);
            }
            reactionProbability += historyReactionProbability;
        }
        reactionProbability /= numParticles;
    }

    @Override
    public Object get() throws InterruptedException, ExecutionException {
        return this.reactionProbability;
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
