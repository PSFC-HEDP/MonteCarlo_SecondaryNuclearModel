package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.RunnableFuture;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Created by lahmann on 2016-09-10.
 */
public class ParticleHistoryTask implements RunnableFuture {

    private NuclearReaction nuclearReaction;
    private ArrayList<StoppingPowerModel> stoppingPowerModels;     // Stopping power profile defined in f(E,r) space (MeV / cm)
    private ParticleDistribution particleDistribution;
    private ParticleType reactingPlasmaSpecies;
    private ArrayList<PlasmaLayer> plasmaLayers;
    private Vector3D detectorLineOfSight;
    private EnergyTally energyTally;
    private int numParticles;
    private boolean debugMode = false;

    private Double taskReactionProbability = 0.0;

    /**
     * Step parameters
     */
    private final int MIN_NUM_STEPS = 200;

    private final double ACCEPTABLE_ENERGY_ERROR   = 1e-3;      // Fraction
    private final double ACCEPTABLE_DISTANCE_ERROR = 1e-3;      // Fraction



    public ParticleHistoryTask(ParticleDistribution particleDistribution, NuclearReaction nuclearReaction,
                               ArrayList<PlasmaLayer> plasmaLayers, ArrayList<StoppingPowerModel> stoppingPowerModels,
                               int numParticles, Vector3D detectorLineOfSight, double[] energyNodes) {

        this.particleDistribution = particleDistribution;
        this.nuclearReaction = nuclearReaction;

        this.plasmaLayers = plasmaLayers;
        this.stoppingPowerModels = stoppingPowerModels;

        this.numParticles = numParticles;
        this.detectorLineOfSight = detectorLineOfSight;
        this.energyTally = new EnergyTally(energyNodes);

        /**
        // Verify that we have a combination of distribution / plasma / nuclear reaction that makes sense
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
         */

    }

    public void setDebugMode(boolean debugMode) {
        this.debugMode = debugMode;
    }

    @Override
    public void run() {

        for (int i = 0; i < numParticles; i++) {

            // TODO: We are currently assuming that particles always get born in the center most layer
            int currentLayer = 0;

            // Sample a particle
            Particle particle = particleDistribution.sample(plasmaLayers.get(currentLayer));

            // DEBUG MODE LINE
            if (debugMode)  System.out.printf("Starting history %d:\n" + particle.toString(), i+1);


            // Determine how far this particle needs to travel
            double distance = getDistanceToPlasmaBoundary(particle, plasmaLayers.get(currentLayer));


            // Calculate a step size
            double dx = distance / MIN_NUM_STEPS;


            // Initialize some variables we'll need
            int totalSteps = 0;                             // This is for debugging only
            double particleReactionProb = 0.0;
            double initialEnergy = particle.getEnergy();


            // While the particle is inside this plasma
            while (distance > dx && particle.getEnergy() >  0.0) {

                // Parameters at the starting point
                double r1 = getNormalizedRadius(particle, plasmaLayers.get(currentLayer));
                double E1 = particle.getEnergy();
                double dEdx1 = stoppingPowerModels.get(currentLayer).evaluate(E1, r1);

                // If the particle is losing a significant amount of energy, we need to take smaller steps
                dx = Math.min(dx, -initialEnergy / dEdx1 / MIN_NUM_STEPS);

                // Calculate r2 before the energy loop to save time
                double r2 = getNormalizedRadius(particle.step(dx, 0.0), plasmaLayers.get(currentLayer));


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
                        r2 = getNormalizedRadius(particle.step(dx, 0.0), plasmaLayers.get(currentLayer));
                    }

                    energyError = FastMath.abs(newEnergy - E2)/E2;
                    E2 = newEnergy;
                    dEdx2 = stoppingPowerModels.get(currentLayer).evaluate(E2, r2);
                }

                // DEBUG MODE LINE
                //if (debugMode)  System.out.printf("-->  After step %d:\n" + particle.toString(), totalSteps+1);

                /**
                 * TODO: I'm commenting this all out as a quick way to get Maria's data.
                 * TODO: Future implementations need to have this be an option

                // Grab the average number density of the species we're interacting with between these two positions
                double n1 = plasma.getSpeciesNumberDensity(particle.getPosition(), reactingPlasmaSpecies);
                double n2 = plasma.getSpeciesNumberDensity(steppedParticle.getPosition(), reactingPlasmaSpecies);
                double n = 0.5*(n1 + n2);


                // Grab the average ion temperature between these two positions
                double Ti1 = plasma.getIonTemperature(particle.getPosition());
                double Ti2 = plasma.getIonTemperature(steppedParticle.getPosition());
                double Ti = 0.5*(Ti1 + Ti2);


                // Create a background particle whose energy in T_ion
                Particle backgroundParticle = new Particle(reactingPlasmaSpecies, 1e-3*Ti);


                // Generate the resultant particle
                Particle productParticle = nuclearReaction.getProductParticleC(particle, backgroundParticle, detectorLineOfSight);
                productParticle.multiplyWeight(particle.getWeight());       // Normalize it to the weigh of it's parent


                // Grab the average cross section between these two positions
                double sigma1 = nuclearReaction.getCrossSection(particle, backgroundParticle);
                double sigma2 = nuclearReaction.getCrossSection(steppedParticle, backgroundParticle);
                double sigma = 0.5*(sigma1 + sigma2);


                // Calculate the reaction probability at this step
                double stepReactionProbability = n * sigma * dx * (1 - particleReactionProb);
                productParticle.multiplyWeight(stepReactionProbability);


                // TODO: The multiplying by energy is an adhoc way to account for the non-isotropic lab frame distribution
                // TODO: I don't currently have a justification for why it works beyond the fact that others have done this in the past
                // Add this to the tally
                energyTally.addTally(productParticle.getEnergy(), productParticle.getWeight() * productParticle.getEnergy());

                 particleReactionProb += stepReactionProbability;
                 */

                // Initialize the next step
                particle = steppedParticle;

                totalSteps++;
                distance -= dx;

                // If we're at the edge of a layer check and see if we're moving into a new one
                if (distance < dx){

                    // If we're closer to the outer boundary of this layer
                    if (plasmaLayers.get(currentLayer).getIsCloserToOuterBoundary(particle.getPosition())){
                        // AND there's another layer beyond this one
                        if (currentLayer < plasmaLayers.size() - 1){
                            currentLayer++;

                            // Place the particle on the inner boundary of the next layer
                            double[] coordinates = Utils.getSphericalFromVector(particle.getPosition());
                            double theta  = coordinates[1];
                            double phi    = coordinates[2];
                            double radius = plasmaLayers.get(currentLayer).getInnerRadiusBound(theta, phi);

                            particle.setPosition(Utils.getVectorFromSpherical(radius, theta, phi));

                            // Calculate the new distance we need to travel
                            distance = getDistanceToPlasmaBoundary(particle, plasmaLayers.get(currentLayer));
                        }
                    }

                    // Otherwise, we're close to the inner boundary
                    else{
                        currentLayer--;

                        // Place the particle on the outer boundary of the previous layer
                        double[] coordinates = Utils.getSphericalFromVector(particle.getPosition());
                        double theta  = coordinates[1];
                        double phi    = coordinates[2];
                        double radius = plasmaLayers.get(currentLayer).getOuterRadiusBound(theta, phi);

                        particle.setPosition(Utils.getVectorFromSpherical(radius, theta, phi));

                        // Calculate the new distance we need to travel
                        distance = getDistanceToPlasmaBoundary(particle, plasmaLayers.get(currentLayer));
                    }
                }
            }

            // TODO: Temp line
            energyTally.addTally(particle.getEnergy(), 1.0);

            /**
            particle.multiplyWeight(particleReactionProb);
            taskReactionProbability += particle.getWeight();
             */
        }

        /**
        // We need to re-normalize the energy tally due to the ad hoc weighting we do
        energyTally.setSum(taskReactionProbability);
         */
    }

    /**
     * Various getters
     */

    public Double getReactionProbability() {
        return taskReactionProbability;
    }

    public EnergyTally getEnergyTally() {
        return energyTally;
    }

    /**
     * Private convenience methods
     */

    private double getNormalizedRadius(Particle particle, PlasmaLayer plasmaLayer){
        double[] coordinates = Utils.getSphericalFromVector(particle.getPosition());

        double R     = coordinates[0];
        double theta = coordinates[1];
        double phi   = coordinates[2];

        double rMax = plasmaLayer.getOuterRadiusBound(theta, phi);
        double rMin = plasmaLayer.getInnerRadiusBound(theta, phi);
        return (R - rMin) / (rMax - rMin);
    }

    private double getDistanceToPlasmaBoundary(Particle particle, PlasmaLayer plasmaLayer){

        // Create a clone that we'll use for the calculation
        Particle tempParticle = particle.clone();

        // We need some kind of guess of what the step size should be to get started
        double plasmaRmax = plasmaLayer.getOuterRadiusBound(0.0, 0.0);
        double plasmaRmin = plasmaLayer.getInnerRadiusBound(0.0, 0.0);
        double dx = (plasmaRmax - plasmaRmin) / MIN_NUM_STEPS;

        // Step our particle until we're no longer in the plasma
        double distance = 0.0;
        while (plasmaLayer.getIsInside(tempParticle.getPosition())){
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
            // We do this by never updating the error calculation if we're outside the plasma
            // As a result, we'll always underestimate the TRUE distance. Never overestimate
            if (plasmaLayer.getIsInside(tempParticle.getPosition())){
                double newDistance = distance + dx;
                distanceError = FastMath.abs(newDistance - distance)/distance;
                distance = newDistance;
            }else{
                tempParticle = tempParticle.step(-2*dx, 0.0);    // Undo last move and step backwards
                distance -= dx;
            }
        }

        return distance;
    }

    public class EnergyTally {

        double[] nodes;
        double[] tallies;

        public EnergyTally(double[] nodes) {
            this.nodes = nodes;
            this.tallies = new double[nodes.length];
        }

        public void addTally(double energy, double weight){
            for (int i = 0; i < nodes.length; i++){
                if (energy < nodes[i]){
                    tallies[i-1] += weight;
                    return;
                }
            }
        }

        public void setSum(double sum){
            double currentSum = 0.0;

            for (double tally : tallies){
                currentSum += tally;
            }

            for (int i = 0; i < tallies.length; i++){
                tallies[i] *= (sum / currentSum);
            }
        }

        public String toString(){
            String string = "   Nodes   |   Value\n";

            for (int i = 0; i < nodes.length; i++){
                string += String.format("%.4e | %.4e\n", nodes[i], tallies[i]);
            }

            return string;
        }
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
        return this.taskReactionProbability == null;
    }

    @Override
    public Object get(long timeout, TimeUnit unit) throws InterruptedException, ExecutionException, TimeoutException {
        return null;
    }
}
