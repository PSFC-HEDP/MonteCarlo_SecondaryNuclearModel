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

    private NuclearReaction nuclearReaction;
    private StoppingPowerModel stoppingPower;     // Stopping power profile defined in f(E,r) space (MeV / cm)
    private ParticleDistribution particleDistribution;
    private ParticleType reactingPlasmaSpecies;
    private Plasma plasma;
    private int numParticles;

    private boolean debugMode = false;

    private Double reactionProbability = 0.0;

    /**
     * Step parameters
     */
    private final int MIN_NUM_STEPS = 200;

    private final double ACCEPTABLE_ENERGY_ERROR   = 1e-3;      // Fraction
    private final double ACCEPTABLE_DISTANCE_ERROR = 1e-3;      // Fraction



    public ParticleHistoryTask(NuclearReaction nuclearReaction, ParticleDistribution particleDistribution, Plasma plasma,
                               StoppingPowerModel stoppingPower, int numParticles) {
        this.nuclearReaction = nuclearReaction;
        this.particleDistribution = particleDistribution;
        this.plasma = plasma;
        this.stoppingPower = stoppingPower;
        this.numParticles = numParticles;

        if (particleDistribution.getType().equals(nuclearReaction.getReactantParticleTypeA())){
            if (plasma.containsSpecies(nuclearReaction.getReactantParticleTypeB())){
                reactingPlasmaSpecies = nuclearReaction.getReactantParticleTypeB();
            }else{
                System.exit(-1);
            }
        }
        else if (particleDistribution.getType().equals(nuclearReaction.getReactantParticleTypeB())){
            if (plasma.containsSpecies(nuclearReaction.getReactantParticleTypeA())){
                reactingPlasmaSpecies = nuclearReaction.getReactantParticleTypeA();
            }else{
                System.exit(-1);
            }
        }
        else {
            System.exit(-1);
        }

    }

    public void setDebugMode(boolean debugMode) {
        this.debugMode = debugMode;
    }

    @Override
    public void run() {

        for (int i = 0; i < numParticles; i++) {

            // Sample a particle
            Particle particle = particleDistribution.sample(plasma);

            // Determine how far this particle needs to travel
            double distance = getDistanceToPlasmaBoundary(particle);

            // Calculate a step size
            double dx = distance / MIN_NUM_STEPS;

            // Initialize some variables we'll need
            int totalSteps = 0;                             // This is for debugging only
            double historyReactionProbability = 0.0;
            double initialEnergy = particle.getEnergy();


            // DEBUG MODE LINE
            if (debugMode){
                Vector3D position = particle.getPosition();
                Vector3D direction = particle.getDirection();
                System.out.printf("Starting history %d with:\n" +
                        "rx = %+.4e cm, ry = %+.4e cm, rz = %+.4e cm\n" +
                        "dx = %+.4e cm, dy = %+.4e cm, dz = %+.4e cm\n" +
                        "E  = %+.4e MeV\n",

                        (i+1), position.getX(), position.getY(), position.getZ(),
                        direction.getX(), direction.getY(), direction.getZ(), particle.getEnergy());
            }


            // While the particle is inside this plasma
            while (distance > dx && particle.getEnergy() >  0.0) {

                // Parameters at the starting point
                double r1 = getNormalizedRadius(particle);
                double E1 = particle.getEnergy();
                double dEdx1 = stoppingPower.evaluate(E1, r1);

                // If the particle is losing a significant amount of energy, we need to take smaller steps
                dx = Math.min(dx, -initialEnergy / dEdx1 / MIN_NUM_STEPS);

                // Calculate r2 before the energy loop to save time
                double r2 = getNormalizedRadius(particle.step(dx, 0.0));


                // Use trapezoidal method to iterate onto the particle's final energy
                double E2 = E1;
                double dEdx2 = dEdx1;
                Particle steppedParticle = particle;
                double energyError = Double.MAX_VALUE;

                while (energyError > ACCEPTABLE_ENERGY_ERROR){
                    double averageStopping = 0.5*(dEdx1 + dEdx2);
                    steppedParticle = particle.step(dx, averageStopping);

                    double newEnergy = steppedParticle.getEnergy();

                    // If we've hit a negative energy, we need to take a smaller step size
                    if (newEnergy < 0.0){
                        newEnergy = 0.0;
                        dx = -E1 / averageStopping;
                        r2 = getNormalizedRadius(particle.step(dx, 0.0));
                    }

                    energyError = FastMath.abs(newEnergy - E2)/E2;
                    E2 = newEnergy;
                    dEdx2 = stoppingPower.evaluate(E2, r2);
                }

                // Grab the average number density of the species we're interacting with between these two positions
                double n1 = plasma.getSpeciesNumberDensity(particle.getPosition(), reactingPlasmaSpecies);
                double n2 = plasma.getSpeciesNumberDensity(steppedParticle.getPosition(), reactingPlasmaSpecies);
                double n = 0.5*(n1 + n2);


                // Grab the average ion temperature between these two positions
                double Ti1 = plasma.getIonTemperature(particle.getPosition());
                double Ti2 = plasma.getIonTemperature(steppedParticle.getPosition());
                double Ti = 0.5*(Ti1 + Ti2);


                // Create a background particle who's energy in T_ion
                Particle backgroundParticle = new Particle(reactingPlasmaSpecies, 1e-3*Ti);


                // Generate the resultant particle
                Particle reactantParticle = nuclearReaction.getProductParticleC(particle, backgroundParticle, Utils.sampleRandomNormalizedVector());


                // Grab the average cross section between these two positions
                double sigma1 = 1e-24 * nuclearReaction.getCrossSection(particle, backgroundParticle);
                double sigma2 = 1e-24 * nuclearReaction.getCrossSection(steppedParticle, backgroundParticle);
                double sigma = 0.5*(sigma1 + sigma2);


                // Calculate the reaction probability
                historyReactionProbability += n * sigma * dx * (1 - historyReactionProbability);
                particle = steppedParticle;

                totalSteps++;
                distance -= dx;
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

    private double getDistanceToPlasmaBoundary(Particle particle){

        // Create a clone that we'll use for the calculation
        Particle tempParticle = particle.clone();

        // We need some kind of guess of what the step size should be to get started
        double plasmaBound = plasma.getRadiusBound(0.0, 0.0);
        double dx = plasmaBound / MIN_NUM_STEPS;

        // Step our particle until we're no longer in the plasma
        double distance = 0.0;
        while (plasma.getIsInside(tempParticle.getPosition())){
            tempParticle = tempParticle.step(dx, 0.0);
            distance += dx;
        }

        // Step back into the plasma
        tempParticle = tempParticle.step(-dx, 0.0);
        distance -= dx;


        // Bisection method the refine our distance calculation
        double distanceError = Double.MAX_VALUE;
        while (distanceError > ACCEPTABLE_DISTANCE_ERROR){
            dx /= 2;
            tempParticle = tempParticle.step(dx, 0.0);


            // To avoid "out of bound" errors, it's very important that we always stay inside the plasma
            // As a result, we'll always underestimate the TRUE distance. Never overestimate
            if (plasma.getIsInside(tempParticle.getPosition())){
                double newDistance = distance + dx;
                distanceError = FastMath.abs(newDistance - distance)/distance;
                distance = newDistance;
            }else{
                tempParticle = tempParticle.step(-dx, 0.0);    // Undo the previous step
            }
        }

        return distance;
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
