package MonteCarloParticleTracer;

import org.apache.commons.math3.analysis.function.Constant;
import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.HashMap;
import java.util.concurrent.ExecutionException;
import java.util.concurrent.RunnableFuture;
import java.util.concurrent.TimeUnit;
import java.util.concurrent.TimeoutException;

/**
 * Created by lahmann on 2016-09-10.
 */
public class ParticleHistoryTask implements RunnableFuture {

    // TODO: Consider making this a user option
    private final double ENERGY_NODE_WIDTH = 0.01;      // MeV

    private ArrayList<Model.LayerData> layerDataList = new ArrayList<>();
    private ParticleDistribution sourceParticleDistribution;

    private Tally[] sourceParticlePositionTallies;
    private Tally[] sourceParticleEnergyTallies;
    private Tally[] sourceParticleTimeTallies;

    private HashMap<NuclearReaction, Tally[]> productParticlePositionTallyMap = new HashMap<>();
    private HashMap<NuclearReaction, Tally[]> productParticleEnergyTallyMap   = new HashMap<>();
    private HashMap<NuclearReaction, Tally[]> productParticleTimeTallyMap     = new HashMap<>();


    private Vector3D detectorLineOfSight;
    private int numParticles;

    private boolean debugMode = false;

    private Double taskReactionProbability = 0.0;

    /**
     * Step parameters
     */
    private final int MIN_NUM_STEPS = 500;

    private final double ACCEPTABLE_ENERGY_ERROR   = 1e-3;      // Fraction
    private final double ACCEPTABLE_DISTANCE_ERROR = 1e-3;      // Fraction



    public ParticleHistoryTask(ParticleDistribution sourceParticleDistribution, int numParticles) {

        this.sourceParticleDistribution = sourceParticleDistribution;
        this.numParticles = numParticles;

        // Create the energy tally for the escaping particles
        double maxEnergy = 1.05 * sourceParticleDistribution.getMaxEnergy();
        int numNodes = (int) Math.ceil(maxEnergy / ENERGY_NODE_WIDTH);
    }

    public void addPlasmaLayer(Model.LayerData data){
        this.layerDataList.add(data);
    }

    public void setDetectorLineOfSight(Vector3D detectorLineOfSight) {
        this.detectorLineOfSight = detectorLineOfSight;
    }

    public void setDebugMode(boolean debugMode) {
        this.debugMode = debugMode;
    }

    @Override
    public void run() {

        setUpTallies();

        for (int i = 0; i < numParticles; i++) {

            //System.out.println(i);

            // TODO: We are currently assuming that particles always get born in the center most layer
            int startingLayerIndex = 0;

            // Sample a particle
            Particle particle = sourceParticleDistribution.sample(layerDataList.get(startingLayerIndex).layer);

            // Trace this particle through the plasma
            traceThroughPlasma(particle, startingLayerIndex, sourceParticlePositionTallies,
                    sourceParticleEnergyTallies, sourceParticleTimeTallies);

        }


    }


    /**
     * Particle tracing methods
     */

    private Particle traceThroughPlasma(Particle particle, int startingLayerIndex, Tally[] positionTallies,
                                        Tally[] energyTallies, Tally[] timeTallies){

        int currentLayerIndex = startingLayerIndex;

        // Add the starting particle to our 0th index array. This will serve as a source tally
        positionTallies[0].addValue(particle.getPosition().getNorm(), particle.getWeight());
        energyTallies  [0].addValue(particle.getEnergy()            , particle.getWeight());
        timeTallies    [0].addValue(particle.getTime()              , particle.getWeight());

        // Grab the current layer data
        Model.LayerData currentLayerData = layerDataList.get(currentLayerIndex);

        boolean finished = false;
        while (!finished) {

            // TODO: Bandage Fix for neutrons having dEdx == 0
            if (particle.getZ() == 0) {
                return particle;
            }

            // Trace the particle through the current layer
            particle = traceThroughPlasmaLayer(particle, currentLayerIndex);
            double weight = particle.getWeight();
            double bias   = weight;

            // TODO: Adhoc weighting
            if (detectorLineOfSight != null) {
                bias *= particle.getEnergy();
            }

            // If the particle is dead, we're done
            if (particle.getEnergy() <= 0) {
                finished = true;
            }

            // Else, check to see if we're near the outer boundary
            else if (currentLayerData.layer.getIsCloserToOuterBoundary(particle.getPosition())) {

                // Tally this particle as crossing the outer surface (+1 is because we're using the 0th index)
                positionTallies[currentLayerIndex+1].addBiasedValue(particle.getPosition().getNorm(), bias, weight);
                energyTallies  [currentLayerIndex+1].addBiasedValue(particle.getEnergy()            , bias, weight);
                timeTallies    [currentLayerIndex+1].addBiasedValue(particle.getTime()              , bias, weight);

                // If so, and this is the last layer we're done
                if (currentLayerIndex == layerDataList.size() - 1) {
                    finished = true;
                }

                // Otherwise, move on to the next layer
                else {
                    currentLayerIndex++;
                    currentLayerData = layerDataList.get(currentLayerIndex);

                    // Move the particle inside the new layer
                    particle = currentLayerData.layer.moveInside(particle);
                }
            }

            // Otherwise, we must be at the inner boundary so we'll move to the previous layer
            else {

                // Tally this particle as crossing the outer surface of the previous layer (there's a +1 and a -1 that cancel out in the indexing)
                positionTallies[currentLayerIndex].addBiasedValue(particle.getPosition().getNorm(), bias, weight);
                energyTallies  [currentLayerIndex].addBiasedValue(particle.getEnergy()            , bias, weight);
                timeTallies    [currentLayerIndex].addBiasedValue(particle.getTime()              , bias, weight);

                currentLayerIndex--;
                currentLayerData = layerDataList.get(currentLayerIndex);

                // Move the particle inside the new layer
                particle = currentLayerData.layer.moveInside(particle);
            }
        }

        return particle;
    }

    private Particle traceThroughPlasmaLayer(Particle particle, int layerIndex){

        // Debug message
        if (debugMode){
            System.out.printf("\nStarting particle with %.4e MeV in layer %d: \n",
                    particle.getEnergy(), layerIndex);
            System.out.println(" Step :    x (cm)   |    y (cm)   |    z (cm)   |  r (norm)  |   E (MeV)  ");
        }

        // Rename for convenience
        Model.LayerData layerData = this.layerDataList.get(layerIndex);

        PlasmaLayer plasmaLayer = layerData.layer;
        StoppingPowerModel stoppingPowerModel = layerData.stoppingPowerModels.get(particle.getType());

        // Determine how far this particle needs to travel
        double distance = getDistanceToPlasmaBoundary(particle, plasmaLayer);

        // TODO: Bandage Fix for neutrons having dEdx == 0
        if (particle.getZ() == 0) {
            return particle.step(distance, 0.0);
        }

        // Calculate a step size
        double maxStepSize = distance / MIN_NUM_STEPS;
        double dx = maxStepSize;

        int totalSteps = 0;                             // This is for debugging only
        double particleReactionProb = 0.0;
        double initialEnergy = particle.getEnergy();


        // While the particle is inside this plasma
        while (distance > dx && particle.getEnergy() >  0.0) {

            // Parameters at the starting point
            double r1 = Utils.getNormalizedRadius(particle, plasmaLayer);
            double E1 = particle.getEnergy();

            // Debug message
            if (debugMode){
                System.out.printf("Stepping %.4e Velocity = %.4e Time = %.4e dT = %.4e\n", dx,
                        particle.getVelocity().getNorm()*Constants.SPEED_OF_LIGHT_CM_PER_SEC,
                        particle.getTime(),
                        dx / (particle.getVelocity().getNorm()*Constants.SPEED_OF_LIGHT_CM_PER_SEC));
                /*System.out.printf("%5d : %+.4e | %+.4e | %+.4e | %.4e | %.4e | %.4e\n",
                        totalSteps,
                        particle.getPosition().getX(),
                        particle.getPosition().getY(),
                        particle.getPosition().getZ(),
                        r1,
                        particle.getEnergy(),
                        particle.getTime());*/
            }


            // If the particle is losing a significant amount of energy, we need to take smaller steps
            double dEdx1 = stoppingPowerModel.evaluate(E1, r1);
            dx = Math.min(maxStepSize, -initialEnergy / dEdx1 / MIN_NUM_STEPS);

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

            // Loop through all of the reactions relevant to this layer
            ArrayList<Model.ReactionData> reactions = layerData.reactionDataMap.get(particle.getType());
            for (Model.ReactionData data : reactions) {

                NuclearReaction nuclearReaction = data.nuclearReaction;
                ParticleType reactingPlasmaSpecies = data.backgroundParticle;

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


                // Grab the average cross section between these two positions
                double sigma1 = nuclearReaction.getCrossSection(particle, backgroundParticle);
                double sigma2 = nuclearReaction.getCrossSection(steppedParticle, backgroundParticle);
                double sigma = 0.5 * (sigma1 + sigma2);


                // Generate the resultant particle
                Particle productParticle = nuclearReaction.getProductParticle(particle, backgroundParticle, detectorLineOfSight);
                productParticle.multiplyWeight(particle.getWeight());       // Normalize it to the weigh of it's parent


                 // Calculate the reaction probability at this step
                 double stepReactionProbability = n * sigma * dx * (1 - particleReactionProb);
                 productParticle.multiplyWeight(stepReactionProbability);

                 // TODO: The multiplying by energy is an adhoc way to account for the non-isotropic lab frame distribution
                 // TODO: I don't currently have a justification for why it works beyond the fact that others have done this in the past
                if (detectorLineOfSight != null) {
                    //productParticle.multiplyWeight(productParticle.getEnergy());
                }

                // Grab the tallies for this product
                Tally[] positionTallies = productParticlePositionTallyMap.get(nuclearReaction);
                Tally[] energyTallies   = productParticleEnergyTallyMap  .get(nuclearReaction);
                Tally[] timeTallies     = productParticleTimeTallyMap    .get(nuclearReaction);

                // Trace this product through the plasma
                traceThroughPlasma(productParticle, layerIndex, positionTallies, energyTallies, timeTallies);

                // Add the updated tallies to the tally map
                productParticlePositionTallyMap.put(nuclearReaction, positionTallies);
                productParticleEnergyTallyMap  .put(nuclearReaction, energyTallies);
                productParticleTimeTallyMap    .put(nuclearReaction, timeTallies);

                // Update the reaction probability
                particleReactionProb += stepReactionProbability;
            }

            // Initialize the next step
            particle = steppedParticle;

            totalSteps++;
            distance -= dx;
        }

        return particle;
    }



    /**
     * Private convenience methods
     */

    private double getDistanceToPlasmaBoundary(Particle particle, PlasmaLayer plasmaLayer){

        // Create a clone that we'll use for the calculation
        Particle tempParticle = particle.clone();

        // We need some kind of guess of what the step size should be to get started
        double plasmaRmax = plasmaLayer.getOuterRadiusBound(0.0, 0.0);
        double plasmaRmin = plasmaLayer.getInnerRadiusBound(0.0, 0.0);
        double dx = (plasmaRmax - plasmaRmin) / MIN_NUM_STEPS;

        // Step our particle until we're no longer in the plasma
        double distance = 0.0;
        do {
            tempParticle = tempParticle.step(dx, 0.0);
            distance += dx;
        }while (plasmaLayer.getIsInside(tempParticle.getPosition()));

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

    private void setUpTallies(){

        // Build the source position tallies
        double[] positionNodes = Utils.linspace(0.0, 75.0*1e-4, 200);      // TODO HARD CODED!
        sourceParticlePositionTallies = new Tally[layerDataList.size() + 1];
        for (int i =0; i < sourceParticlePositionTallies.length; i++){
            sourceParticlePositionTallies[i] = new Tally(positionNodes);
        }


        // Build the source time tallies
        double[] timeNodes = Utils.linspace(0.0, 100.0*1e-12, 200);         // TODO HARD CODED!
        sourceParticleTimeTallies = new Tally[layerDataList.size() + 1];
        for (int i =0; i < sourceParticleTimeTallies.length; i++){
            sourceParticleTimeTallies[i] = new Tally(timeNodes);
        }


        // Build the source energy tallies
        double maxSourceEnergy = 1.1*sourceParticleDistribution.getMaxEnergy();
        int numSourceNodes = (int) Math.ceil(maxSourceEnergy / ENERGY_NODE_WIDTH);
        double[] sourceNodes = Utils.linspace(0.0, maxSourceEnergy, numSourceNodes);
        sourceParticleEnergyTallies   = new Tally[layerDataList.size() + 1];
        for (int i =0; i < sourceParticlePositionTallies.length; i++){
            sourceParticleEnergyTallies[i] = new Tally(sourceNodes);
        }


        Particle maxEnergySourceParticle = new Particle(sourceParticleDistribution.getType(), maxSourceEnergy);

        // For every layer in this model
        for (Model.LayerData layerData : this.layerDataList){

            // For each list of reactions
            for (ParticleType key : layerData.reactionDataMap.keySet()){

                // For each reaction
                for (Model.ReactionData reactionData : layerData.reactionDataMap.get(key)){

                    // Get the nuclear reaction
                    NuclearReaction nuclearReaction = reactionData.nuclearReaction;


                    // Build the product particle position tallies
                    Tally[] productParticlePositionTallies = new Tally[layerDataList.size() + 1];
                    for (int i =0; i < productParticlePositionTallies.length; i++){
                        productParticlePositionTallies[i] = new Tally(positionNodes);
                    }
                    productParticlePositionTallyMap.put(nuclearReaction, productParticlePositionTallies);


                    // Build the product particle time tallies
                    Tally[] productParticleTimeTallies = new Tally[layerDataList.size() + 1];
                    for (int i =0; i < productParticleTimeTallies.length; i++){
                        productParticleTimeTallies[i] = new Tally(timeNodes);
                    }
                    productParticleTimeTallyMap.put(nuclearReaction, productParticleTimeTallies);


                    // Build the product particle energy tallies
                    Particle backgroundParticle = new Particle(reactionData.backgroundParticle, 0.0);
                    Particle maxEnergyProductParticle = nuclearReaction.getMaxEnergyProductParticle(maxEnergySourceParticle, backgroundParticle);

                    double maxProductEnergy = 1.1*maxEnergyProductParticle.getEnergy();
                    int numProductNodes = (int) Math.ceil(maxProductEnergy / ENERGY_NODE_WIDTH);
                    double[] productNodes = Utils.linspace(0.0, maxProductEnergy, numProductNodes);

                    Tally[] productParticleEnergyTallies = new Tally[layerDataList.size() + 1];
                    for (int i = 0; i < productParticleEnergyTallies.length; i++){
                        productParticleEnergyTallies[i] = new Tally(productNodes);
                    }
                    productParticleEnergyTallyMap.put(nuclearReaction, productParticleEnergyTallies);

                }
            }
        }
    }



    /**
     * Various getters
     */

    public Tally[] getSourceParticlePositionTallies() {
        return sourceParticlePositionTallies;
    }

    public Tally[] getSourceParticleEnergyTallies() {
        return sourceParticleEnergyTallies;
    }

    public Tally[] getSourceParticleTimeTallies() {
        return sourceParticleTimeTallies;
    }

    public HashMap<NuclearReaction, Tally[]> getProductParticlePositionTallyMap() {
        return productParticlePositionTallyMap;
    }

    public HashMap<NuclearReaction, Tally[]> getProductParticleEnergyTallyMap() {
        return productParticleEnergyTallyMap;
    }

    public HashMap<NuclearReaction, Tally[]> getProductParticleTimeTallyMap() {
        return productParticleTimeTallyMap;
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
