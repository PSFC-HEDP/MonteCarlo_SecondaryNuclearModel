package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

import java.lang.reflect.Array;
import java.util.ArrayList;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.RunnableFuture;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Created by lahmann on 2016-09-10.
 */
public class ParticleHistoryTask implements RunnableFuture {

    // TODO: Consider making this a user option
    private final double ENERGY_NODE_WIDTH = 0.1;      // MeV

    private ArrayList<Model.LayerData> layerData = new ArrayList<>();
    private ParticleDistribution particleDistribution;

    private Vector3D detectorLineOfSight;
    private int numParticles;

    private Tally escapingParticleEnergyTally;

    private boolean debugMode = false;



    private Double taskReactionProbability = 0.0;

    /**
     * Step parameters
     */
    private final int MIN_NUM_STEPS = 200;

    private final double ACCEPTABLE_ENERGY_ERROR   = 1e-3;      // Fraction
    private final double ACCEPTABLE_DISTANCE_ERROR = 1e-3;      // Fraction



    public ParticleHistoryTask(ParticleDistribution particleDistribution, int numParticles) {

        this.particleDistribution = particleDistribution;
        this.numParticles = numParticles;

        // Create the energy tally for the escaping particles
        double maxEnergy = 1.05 * particleDistribution.getMaxEnergy();
        int numNodes = (int) Math.ceil(maxEnergy / ENERGY_NODE_WIDTH);
        escapingParticleEnergyTally = new Tally(Utils.linspace(0.0, maxEnergy, numNodes));
    }

    public void addPlasmaLayer(Model.LayerData data){
        this.layerData.add(data);
    }

    public void setDetectorLineOfSight(Vector3D detectorLineOfSight) {
        this.detectorLineOfSight = detectorLineOfSight;
    }

    public void setDebugMode(boolean debugMode) {
        this.debugMode = debugMode;
    }

    @Override
    public void run() {

        for (int i = 0; i < numParticles; i++) {

            // TODO: We are currently assuming that particles always get born in the center most layer
            int currentLayerIndex = 0;

            // Grab the current layer data
            Model.LayerData currentLayerData = layerData.get(currentLayerIndex);

            // Sample a particle
            Particle particle = particleDistribution.sample(currentLayerData.layer);

            boolean finished = false;
            while (!finished) {

                // Trace the particle through the current layer
                particle = traceThroughPlasmaLayer(particle, currentLayerData);

                // If the particle is dead, we're done
                if (particle.getEnergy() <= 0) {
                    finished = true;
                }

                // Check to see if we're near the outer boundary
                else if (currentLayerData.layer.getIsCloserToOuterBoundary(particle.getPosition())) {

                    // If so, and this is the last layer we're done
                    if (currentLayerIndex == layerData.size() - 1) {
                        finished = true;

                        // Add the escaped particle to our escaping tally
                        escapingParticleEnergyTally.addValue(particle.getEnergy(), particle.getWeight());
                    }

                    // Otherwise, move on to the next layer
                    else {
                        currentLayerIndex++;
                        currentLayerData = layerData.get(currentLayerIndex);
                    }
                }

                // Otherwise, we must be at the inner boundary so we'll move to the previous layer
                else {
                    currentLayerIndex--;
                    currentLayerData = layerData.get(currentLayerIndex);
                }
            }

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

    public Tally getEscapingParticleEnergyTally() {
        return escapingParticleEnergyTally;
    }

    /**
     * Private convenience methods
     */

    private Particle traceThroughPlasmaLayer(Particle particle, Model.LayerData layerData){

        // Rename for convenience
        PlasmaLayer plasmaLayer = layerData.layer;
        StoppingPowerModel stoppingPowerModel = layerData.stoppingPowerModel;

        // Determine how far this particle needs to travel
        double distance = getDistanceToPlasmaBoundary(particle, plasmaLayer);

        // Calculate a step size
        double dx = distance / MIN_NUM_STEPS;

        int totalSteps = 0;                             // This is for debugging only
        double particleReactionProb = 0.0;
        double initialEnergy = particle.getEnergy();


        // While the particle is inside this plasma
        while (distance > dx && particle.getEnergy() >  0.0) {

            // Parameters at the starting point
            double r1 = Utils.getNormalizedRadius(particle, plasmaLayer);
            double E1 = particle.getEnergy();
            double dEdx1 = stoppingPowerModel.evaluate(E1, r1);

            // If the particle is losing a significant amount of energy, we need to take smaller steps
            dx = Math.min(dx, -initialEnergy / dEdx1 / MIN_NUM_STEPS);

            // Calculate r2 before the energy loop to save time
            double r2 = Utils.getNormalizedRadius(particle.step(dx, 0.0), plasmaLayer);


            // Use trapezoidal method to iterate onto the particle's final energy
            double E2 = E1;
            double dEdx2 = dEdx1;
            Particle steppedParticle = particle;

            double energyError = Double.MAX_VALUE;
            while (energyError > ACCEPTABLE_ENERGY_ERROR) {
                double averageStopping = 0.5 * (dEdx1 + dEdx2);
                steppedParticle = particle.step(dx, averageStopping);

                double newEnergy = steppedParticle.getEnergy();

                // If we've hit a negative energy, we need to take a smaller step size
                if (newEnergy < 0.0) {
                    newEnergy = 0.0;
                    dx = -E1 / averageStopping;
                    r2 = Utils.getNormalizedRadius(particle.step(dx, 0.0), plasmaLayer);
                }

                energyError = FastMath.abs(newEnergy - E2) / E2;
                E2 = newEnergy;
                dEdx2 = stoppingPowerModel.evaluate(E2, r2);
            }

            // DEBUG MODE LINE
            //if (debugMode)  System.out.printf("-->  After step %d:\n" + particle.toString(), totalSteps+1);

            // Loop through all of the reactions relevant to this layer
            for (Model.ReactionData data : layerData.reactionData) {

                NuclearReaction nuclearReaction = data.nuclearReaction;
                ParticleType reactingPlasmaSpecies = data.reactingSpecies;

                // Grab the average number density of the species we're interacting with between these two positions
                double n1 = plasmaLayer.getSpeciesNumberDensity(particle.getPosition(), reactingPlasmaSpecies);
                double n2 = plasmaLayer.getSpeciesNumberDensity(steppedParticle.getPosition(), reactingPlasmaSpecies);
                double n = 0.5 * (n1 + n2);


                // Grab the average ion temperature between these two positions
                double Ti1 = plasmaLayer.getIonTemperature(particle.getPosition());
                double Ti2 = plasmaLayer.getIonTemperature(steppedParticle.getPosition());
                double Ti = 0.5 * (Ti1 + Ti2);


                // Create a background particle whose energy in T_ion
                Particle backgroundParticle = new Particle(reactingPlasmaSpecies, 1e-3 * Ti);


                // Generate the resultant particle
                Particle productParticle = nuclearReaction.getProductParticle(particle, backgroundParticle, detectorLineOfSight);
                productParticle.multiplyWeight(particle.getWeight());       // Normalize it to the weigh of it's parent


                // Grab the average cross section between these two positions
                double sigma1 = nuclearReaction.getCrossSection(particle, backgroundParticle);
                double sigma2 = nuclearReaction.getCrossSection(steppedParticle, backgroundParticle);
                double sigma = 0.5 * (sigma1 + sigma2);


                //TODO: Fix the tallying of nuclear products
                /**
                // Calculate the reaction probability at this step
                double stepReactionProbability = n * sigma * dx * (1 - particleReactionProb);
                productParticle.multiplyWeight(stepReactionProbability);


                // TODO: The multiplying by energy is an adhoc way to account for the non-isotropic lab frame distribution
                // TODO: I don't currently have a justification for why it works beyond the fact that others have done this in the past
                // Add this to the tally
                energyTally.addTally(productParticle.getEnergy(), productParticle.getWeight() * productParticle.getEnergy());
                particleReactionProb += stepReactionProbability;
                 */
            }

            // Initialize the next step
            particle = steppedParticle;

            totalSteps++;
            distance -= dx;
        }

        return particle;
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
