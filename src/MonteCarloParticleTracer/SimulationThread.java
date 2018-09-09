package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;
import org.apache.commons.math3.util.FastMath;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by lahmann on 2016-09-10.
 */
public class SimulationThread extends Thread {

    /**
     * Magic constants
     */

    // Width between energy nodes used in tallies                       TODO: Make a user input
    private final double ENERGY_NODE_WIDTH = 0.01;                      // MeV

    // Number of steps each source particle will take through any plasma
    private final int MIN_NUM_STEPS = 200;

    // Tolerance for calculating the particle energy of the next step
    private final double ACCEPTABLE_ENERGY_ERROR   = 1e-3;              // Fraction

    // Tolerance for calculating the distance that a particle must travel
    private final double ACCEPTABLE_DISTANCE_ERROR = 1e-3;              // Fraction



    /**
     * Everything else
     */

    // ORDERED Array list with all of the plasma information
    private ArrayList<Model.PlasmaData> plasmaData = new ArrayList<>();

    // Plasma from which the source particles originate
    private Plasma sourcePlasma;

    // Radial source distribution
    private Distribution radialSourceDistribution;

    // Polar source distribution
    private Distribution polarSourceDistribution;

    // Nuclear reaction that we sample the source particles from
    private NuclearReaction sourceReaction;

    // DoubleArray of tallies for the source particles (1 for birth and 1 for each surface)
    private Tally[] sourceParticlePositionTallies;
    private Tally[] sourceParticleEnergyTallies;
    private Tally[] sourceParticleTimeTallies;

    // Hashmap of tally arrays for every nuclear reaction we're modelling
    private HashMap<NuclearReaction, Tally[]> productParticlePositionTallyMap = new HashMap<>();
    private HashMap<NuclearReaction, Tally[]> productParticleEnergyTallyMap   = new HashMap<>();
    private HashMap<NuclearReaction, Tally[]> productParticleTimeTallyMap     = new HashMap<>();

    // The vector that all source particles will be forced to travel
    // Secondary reactions cannot be trusted when this is used
    private Vector3D sourceDirection = null;

    // The vector that all secondary reactions will be forced to travel
    private Vector3D productDirection = null;

    // Number of particles simulated by this thread
    private int numParticles;

    // Weight that all particles are born with by default
    private double birthWeight;

    // Logger object for keeping track of this threads process
    private Logger logger = new Logger(false);

    // Debug mode flag
    private boolean debugMode = false;


    /**
     * Default constructor
     * @param numThreadParticles number of particles to be simulated by this thread
     * @param numTotalParticles total number of particles in the entire simulation (for proper weighting)
     */
    SimulationThread(int numThreadParticles, int numTotalParticles) {
        this.numParticles = numThreadParticles;
        this.birthWeight  = 1.0 / numTotalParticles;
    }


    /**
     * Method for setting the source information
     * @param sourcePlasma Plasma where the source particles originate
     * @param radialSourceDistribution Radial distribution to same source particles from
     * @param sourceReaction Reaction that we use to sample the source particles
     */
    void setSourceInformation(Plasma sourcePlasma, NuclearReaction sourceReaction, Distribution radialSourceDistribution){

        // If the plasma bounds vary in phi, we need to proper 2D distribution for the angular dependence
        if (sourcePlasma.boundsVaryInPhi()){
            System.err.println("Code does not currently support sampling a 3D plasma");
            System.exit(-1);
        }


        // If the plasma only varies in theta, we only need a 1D polar distribution
        else if (sourcePlasma.boundsVaryInTheta()){
            polarSourceDistribution = sourcePlasma.getPolarDistribution();

        }

        // Set the radial distribution
        this.radialSourceDistribution = radialSourceDistribution;



        // Set the source plasma and reaction (make unique copies)
        this.sourcePlasma   = sourcePlasma.copy();
        this.sourceReaction = sourceReaction.copy();

    }


    void setRadialSourceDistribution(Distribution radialSourceDistribution){
        this.radialSourceDistribution = radialSourceDistribution;
    }


    /**
     * Method for adding a plasma to the thread. Plasmas must be added in order (inside -> outside)
     * @param data Data object that stores everything relevant to the simulation
     */
    void addPlasma(Model.PlasmaData data){
        this.plasmaData.add(data);
    }



    public void setSourceDirection(Vector3D sourceDirection) {
        this.sourceDirection = sourceDirection;
    }



    public void setProductDirection(Vector3D productDirection) {
        this.productDirection = productDirection;
    }



    @Override
    /**
     * Main method of this thread
     */
    public void run() {

        for (int i = 0; i < numParticles; i++) {

            //System.out.println(i + " - " + Thread.currentThread().getId());

            if (debugMode)  logger.startTimer("Source Particle Lifespan");

            // We're assuming particles are always born in the center plasma TODO: Fix this
            int startingLayerIndex = 0;

            // Sample a particle
            if (debugMode)  logger.startTimer("Sample source particle");
            Particle particle = sample(birthWeight, sourceDirection);
            if (debugMode)  logger.stopTimer("Sample source particle");

            // Simulate this particle
            traceParticle(particle, startingLayerIndex, sourceParticlePositionTallies,
                    sourceParticleEnergyTallies, sourceParticleTimeTallies);

            if (debugMode)  logger.stopTimer("Source Particle Lifespan");
        }

        if (debugMode)  logger.dumpLogToConsole();
    }


    /**
     * Method that handles tracing and tallying a particle through the entire geometry
     * Moves the particle between plasmas as surfaces are crossed
     * @param particle Particle that we want to simulate
     * @param startingLayerIndex Index of the plasma that this particle started in
     * @param positionTallies Position tallies corresponding to this particle
     * @param energyTallies Energy tallies corresponding to this particle
     * @param timeTallies Time tallies corresponding to this particle
     * @return This particle after it escapes or dies
     */
    private Particle traceParticle(Particle particle, int startingLayerIndex, Tally[] positionTallies,
                                   Tally[] energyTallies, Tally[] timeTallies){

        int currentLayerIndex = startingLayerIndex;

        // Add the starting particle to our 0th index array. This will serve as a source tally
        positionTallies[0].addValue(particle.getPosition().getNorm(), particle.getWeight());
        energyTallies  [0].addValue(particle.getEnergy()            , particle.getWeight());
        timeTallies    [0].addValue(particle.getTime()              , particle.getWeight());

        // Grab the current plasma data
        Model.PlasmaData currentPlasmaData = plasmaData.get(currentLayerIndex);


        // Loop until the particle dies or escapes
        while ( true ) {

            // We don't want to bother simulating neutrons for now TODO: Fix this
            if (particle.getZ() == 0) {
                return particle;
            }


            // Trace the particle through the current plasma
            particle = traceThroughPlasma(particle, currentLayerIndex);


            // If the particle is dead, we're done
            if (particle.getEnergy() <= 0) {
                return particle;
            }


            // Figure out the weight and bias of this particle for tallies
            double weight = particle.getWeight();
            double bias   = weight;


            // If we're at the inner boundary, we need to move down to the next inner plasma
            if (currentPlasmaData.plasma.getIsCloserToInnerBoundary(particle.getPosition())){

                // Tally this particle as crossing the outer surface of the previous plasma (there's a +1 and a -1 that cancel out in the indexing)
                positionTallies[currentLayerIndex].addBiasedValue(particle.getPosition().getNorm(), bias, weight);
                energyTallies  [currentLayerIndex].addBiasedValue(particle.getEnergy()            , bias, weight);
                timeTallies    [currentLayerIndex].addBiasedValue(particle.getTime()              , bias, weight);

                // Decrease the index by 1
                currentLayerIndex--;
                currentPlasmaData = plasmaData.get(currentLayerIndex);

                // Move the particle inside the new plasma
                particle = currentPlasmaData.plasma.moveInside(particle);

            }

            // Otherwise we're at the outer boundary
            else {

                // Tally this particle as crossing the outer surface (+1 is because we're using the 0th index)
                positionTallies[currentLayerIndex+1].addBiasedValue(particle.getPosition().getNorm(), bias, weight);
                energyTallies  [currentLayerIndex+1].addBiasedValue(particle.getEnergy()            , bias, weight);
                timeTallies    [currentLayerIndex+1].addBiasedValue(particle.getTime()              , bias, weight);


                // If this is the outer most plasma, we've escaped
                if (currentLayerIndex == plasmaData.size() - 1) {
                    return particle;
                }


                // If not, we need to move to the next plasma
                currentLayerIndex++;
                currentPlasmaData = plasmaData.get(currentLayerIndex);


                // Move the particle inside the new plasma
                particle = currentPlasmaData.plasma.moveInside(particle);

            }
        }
    }


    /**
     * Method that traces the particle through a single plasma
     * @param particle Particle that we want to simulate
     * @param layerIndex Index of the plasma we're tracing through
     * @return This particle after it escapes or dies
     */
    private Particle traceThroughPlasma(Particle particle, int layerIndex){

        // Get the plasma data
        Model.PlasmaData plasmaData = this.plasmaData.get(layerIndex);


        // Rename the fields for convenience
        Plasma plasma = plasmaData.plasma;
        StoppingPowerModel stoppingPowerModel = plasmaData.stoppingPowerModels.get(particle.getType());


        // Determine how far this particle needs to travel
        if (debugMode)  logger.startTimer("Get Distance to Boundary");
        double distance = getDistanceToPlasmaBoundary(particle, plasma);
        if (debugMode)  logger.stopTimer("Get Distance to Boundary");


        // Calculate a step size
        double dx = distance / MIN_NUM_STEPS;

        int totalSteps = 0;                             // This is for debugging only
        double particleReactionProb = 0.0;
        double initialEnergy = particle.getEnergy();


        // While the particle is inside this plasma
        while (distance > dx && particle.getEnergy() >  0.0) {

            if (debugMode)  logger.startTimer("Particle Step");

            // Parameters at the starting point
            double r1 = Utils.getNormalizedRadius(particle, plasma);
            double E1 = particle.getEnergy();


            // If the particle is losing a significant amount of energy, we need to take smaller steps
            if (debugMode)  logger.startTimer("Determine next position");
            double dEdx1 = stoppingPowerModel.evaluate(E1, r1);
            dx = Math.min(dx, -initialEnergy / dEdx1 / MIN_NUM_STEPS);


            // We need to ensure that the particle doesn't step outside the plasma
            dx = Math.min(dx, distance);


            // Calculate r2 before the energy loop to save time
            double r2 = Utils.getNormalizedRadius(particle.step(dx, 0.0), plasma);
            if (debugMode)  logger.stopTimer("Determine next position");


            // Use trapezoidal method to iterate onto the particle's final energy
            if (debugMode)  logger.startTimer("Determine next energy");
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
                    r2 = Utils.getNormalizedRadius(particle.step(dx, 0.0), plasma);
                }

                energyError = FastMath.abs(newEnergy - E2) / E2;
                E2 = newEnergy;
                dEdx2 = stoppingPowerModel.evaluate(E2, r2);
            }
            if (debugMode)  logger.stopTimer("Determine next energy");

            // Loop through all of the reactions relevant to this plasma
            ArrayList<Model.ReactionData> reactions = plasmaData.reactionDataMap.get(particle.getType());
            for (Model.ReactionData data : reactions) {

                NuclearReaction nuclearReaction = data.nuclearReaction;
                ParticleType reactingPlasmaSpecies = data.backgroundParticle;

                // Grab the average number density of the species we're interacting with between these two positions
                if (debugMode)  logger.startTimer("Get average number density");
                double n1 = plasma.getSpeciesNumberDensity(particle.getPosition(), reactingPlasmaSpecies);
                double n2 = plasma.getSpeciesNumberDensity(steppedParticle.getPosition(), reactingPlasmaSpecies);
                double n = 0.5 * (n1 + n2);
                if (debugMode)  logger.stopTimer("Get average number density");


                // Grab the average ion temperature between these two positions
                if (debugMode)  logger.startTimer("Get average ion temperature");
                double Ti1 = plasma.getIonTemperature(particle.getPosition());
                double Ti2 = plasma.getIonTemperature(steppedParticle.getPosition());
                double Ti = 0.5 * (Ti1 + Ti2);
                if (debugMode)  logger.stopTimer("Get average ion temperature");


                // Create a background particle whose energy in T_ion
                Particle backgroundParticle = new Particle(reactingPlasmaSpecies, 1e-3 * Ti);


                // Grab the average cross section between these two positions
                if (debugMode)  logger.startTimer("Get average cross section");
                double sigma1 = nuclearReaction.getCrossSection(particle, backgroundParticle);
                double sigma2 = nuclearReaction.getCrossSection(steppedParticle, backgroundParticle);
                double sigma = 0.5 * (sigma1 + sigma2);
                if (debugMode)  logger.stopTimer("Get average cross section");


                // Generate the resultant particle
                if (debugMode)  logger.startTimer("Generate secondary particle");
                Particle productParticle = nuclearReaction.getProductParticle(particle, backgroundParticle, productDirection);
                productParticle.multiplyWeight(particle.getWeight());       // Normalize it to the weigh of it's parent
                if (debugMode)  logger.stopTimer("Generate secondary particle");


                // Calculate the reaction probability at this step
                double stepReactionProbability = n * sigma * dx * (1 - particleReactionProb);
                productParticle.multiplyWeight(stepReactionProbability);

                // Grab the tallies for this product
                if (debugMode)  logger.startTimer("Tally secondary particle");
                Tally[] positionTallies = productParticlePositionTallyMap.get(nuclearReaction);
                Tally[] energyTallies   = productParticleEnergyTallyMap  .get(nuclearReaction);
                Tally[] timeTallies     = productParticleTimeTallyMap    .get(nuclearReaction);

                // Trace this product through the plasma
                traceParticle(productParticle, layerIndex, positionTallies, energyTallies, timeTallies);

                // Add the updated tallies to the tally map
                productParticlePositionTallyMap.put(nuclearReaction, positionTallies);
                productParticleEnergyTallyMap  .put(nuclearReaction, energyTallies);
                productParticleTimeTallyMap    .put(nuclearReaction, timeTallies);
                if (debugMode)  logger.stopTimer("Tally secondary particle");

                // Update the reaction probability
                particleReactionProb += stepReactionProbability;

            }

            // Initialize the next step
            particle = steppedParticle;

            totalSteps++;
            distance -= dx;

            if (debugMode)  logger.stopTimer("Particle Step");
        }

        return particle;
    }



    /**
     * Private convenience methods
     */

    private double getDistanceToPlasmaBoundary(Particle particle, Plasma plasma){

        // Create a clone that we'll use for the calculation
        Particle tempParticle = particle.clone();

        // We need some kind of initial guess of the distance TODO: Sample random positions?
        double plasmaRmax = plasma.getOuterRadiusBound(0.0, 0.0);
        double plasmaRmin = plasma.getInnerRadiusBound(0.0, 0.0);
        double dx = (plasmaRmax - plasmaRmin);


        // Bisection method the refine our distance calculation
        double distance = 0;
        double distanceError = Double.MAX_VALUE;
        while (distanceError > ACCEPTABLE_DISTANCE_ERROR && (dx/distance) > ACCEPTABLE_DISTANCE_ERROR){

            /*
             *      To avoid "out of bound" errors, it's very important that we always stay inside the plasma
             *      We do this by never updating the error calculation if we're outside the plasma
             *      As a result, we'll always underestimate the TRUE distance. Never overestimate
             */

            // Try stepping dx
            tempParticle = tempParticle.step(dx, 0.0);

            // If we're still in the plasma update the distance
            if (plasma.getIsInside(tempParticle.getPosition())){

                double newDistance = distance + dx;
                if (distance != 0){
                    distanceError = FastMath.abs(newDistance - distance)/distance;
                }
                distance = newDistance;

            }

            // If we left the plasma, undo the step and reduce the step size
            else {
                tempParticle = tempParticle.step(-dx, 0.0);
                dx *= 0.5;
            }
        }

        return distance;
    }

    void setUpTallies(){

        // Build the source position tallies
        double[] positionNodes = DoubleArray.linspace(0.0, sourcePlasma.getOuterRadiusBound(0.0, 0.0), 200).getValues();      // TODO HARD CODED!
        sourceParticlePositionTallies = new Tally[plasmaData.size() + 1];
        for (int i =0; i < sourceParticlePositionTallies.length; i++){
            sourceParticlePositionTallies[i] = new Tally(positionNodes);
        }


        // Build the source time tallies
        double[] timeNodes = DoubleArray.linspace(0.0, 100.0*1e-12, 200).getValues();         // TODO HARD CODED!
        sourceParticleTimeTallies = new Tally[plasmaData.size() + 1];
        for (int i =0; i < sourceParticleTimeTallies.length; i++){
            sourceParticleTimeTallies[i] = new Tally(timeNodes);
        }


        // Build the source energy tallies
        double maxSourceEnergy = 2.0*sourceReaction.getZeroTemperatureMeanEnergy();     // TODO: HARD CODED!
        int numSourceNodes = (int) Math.ceil(maxSourceEnergy / ENERGY_NODE_WIDTH);
        double[] sourceNodes = DoubleArray.linspace(0.0, maxSourceEnergy, numSourceNodes).getValues();
        sourceParticleEnergyTallies   = new Tally[plasmaData.size() + 1];
        for (int i =0; i < sourceParticlePositionTallies.length; i++){
            sourceParticleEnergyTallies[i] = new Tally(sourceNodes);
        }


        Particle maxEnergySourceParticle = new Particle(sourceReaction.getProducts()[0], maxSourceEnergy);

        // For every plasma in this model
        for (Model.PlasmaData layerData : this.plasmaData){

            // For each list of reactions
            for (ParticleType key : layerData.reactionDataMap.keySet()){

                // For each reaction
                for (Model.ReactionData reactionData : layerData.reactionDataMap.get(key)){

                    // Get the nuclear reaction
                    NuclearReaction nuclearReaction = reactionData.nuclearReaction;


                    // Build the product particle position tallies
                    Tally[] productParticlePositionTallies = new Tally[plasmaData.size() + 1];
                    for (int i =0; i < productParticlePositionTallies.length; i++){
                        productParticlePositionTallies[i] = new Tally(positionNodes);
                    }
                    productParticlePositionTallyMap.put(nuclearReaction, productParticlePositionTallies);


                    // Build the product particle time tallies
                    Tally[] productParticleTimeTallies = new Tally[plasmaData.size() + 1];
                    for (int i =0; i < productParticleTimeTallies.length; i++){
                        productParticleTimeTallies[i] = new Tally(timeNodes);
                    }
                    productParticleTimeTallyMap.put(nuclearReaction, productParticleTimeTallies);


                    // Build the product particle energy tallies
                    Particle backgroundParticle = new Particle(reactionData.backgroundParticle, 0.0);
                    Particle maxEnergyProductParticle = nuclearReaction.getMaxEnergyProductParticle(maxEnergySourceParticle, backgroundParticle);

                    double maxProductEnergy = 1.1*maxEnergyProductParticle.getEnergy();
                    int numProductNodes = (int) Math.ceil(maxProductEnergy / ENERGY_NODE_WIDTH);
                    double[] productNodes = DoubleArray.linspace(0.0, maxProductEnergy, numProductNodes).getValues();

                    Tally[] productParticleEnergyTallies = new Tally[plasmaData.size() + 1];
                    for (int i = 0; i < productParticleEnergyTallies.length; i++){
                        productParticleEnergyTallies[i] = new Tally(productNodes);
                    }
                    productParticleEnergyTallyMap.put(nuclearReaction, productParticleEnergyTallies);

                }
            }
        }
    }

    private Particle sample(double birthWeight, Vector3D direction){

        // Start by sampling a spherically random normal vector
        Vector3D position = Utils.sampleRandomNormalizedVector();

        // Convert it to spherical
        double[] coordinates = Utils.getSphericalFromVector(position);
        double theta = coordinates[1];
        double phi   = coordinates[2];


        // If we have a polar distribution, use it instead
        if (polarSourceDistribution != null){
            theta = polarSourceDistribution.sample();
            position = Utils.getVectorFromSpherical(1.0, theta, phi);
        }


        // Use our radial distribution and the plasma bounds to determine the magnitude of our position vector
        double rNorm = Math.abs(radialSourceDistribution.sample());         // TODO: Temp bug fix
        double rMax  = sourcePlasma.getOuterRadiusBound(theta, phi);
        double rMin  = sourcePlasma.getInnerRadiusBound(theta, phi);

        double r = (rMax - rMin) * rNorm + rMin;


        // Rescale our position vector to this magnitude
        position = position.scalarMultiply(r);


        // Sample the direction
        if (direction == null) {
            direction = Utils.sampleRandomNormalizedVector();
        }

        // Build the energy distribution based on our position in the plasma
        double mu = sourceReaction.getZeroTemperatureMeanEnergy();              // TODO: Include temperature up-shift
        double sigma = sourcePlasma.getThermalSigma(position, sourceReaction);


        // Sample the energy
        double energy = Distribution.normDistribution(mu, sigma).sample();


        // Get the particle type
        ParticleType type = sourceReaction.getProducts()[0];


        // Return the sample particle
        Particle particle = new Particle(type, position, direction, energy, 0.0);
        particle.setWeight(birthWeight);                                        // TODO: Use actual XS?
        return particle;
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


}
