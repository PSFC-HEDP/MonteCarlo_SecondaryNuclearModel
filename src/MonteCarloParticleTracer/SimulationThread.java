package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.SphericalCoordinates;
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
    private final int MIN_NUM_STEPS = 500;

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

    // Nuclear reaction that we sample the source particles from
    private NuclearReaction sourceReaction;

    // Array of tallies for the source particles (1 for birth and 1 for each surface)
    private Tally[] sourceParticlePositionTallies;
    private Tally[] sourceParticleEnergyTallies;
    private Tally[] sourceParticleTimeTallies;

    // Hashmap of tally arrays for every nuclear reaction we're modelling
    private HashMap<NuclearReaction, Tally[]> productParticlePositionTallyMap = new HashMap<>();
    private HashMap<NuclearReaction, Tally[]> productParticleEnergyTallyMap   = new HashMap<>();
    private HashMap<NuclearReaction, Tally[]> productParticleTimeTallyMap     = new HashMap<>();

    // Position where the data will be tallied (averages over 4pi if null)
    private Vector3D detectorLineOfSight;

    // Number of particles simulated by this thread
    private int numParticles;

    // Logger object for keeping track of this threads process
    private Logger logger = new Logger(false);



    /**
     * Default constructor
     * @param numParticles Number of particles simulated by this thread
     */
    SimulationThread(int numParticles) {
        this.numParticles = numParticles;
    }


    /**
     * Method for setting the source information
     * @param sourcePlasma Plasma where the source particles originate
     * @param reactivity Reactivity of the reaction the generates the source particles
     * @param sourceReaction Reaction that we use to sample the source particles
     */
    void setSourceInformation(Plasma sourcePlasma, Reactivity reactivity, NuclearReaction sourceReaction){

        // Generate the radial distribution using the reactivity
        this.radialSourceDistribution = sourcePlasma.getSpatialBurnDistribution(reactivity);

        // Set the source plasma and reaction
        this.sourcePlasma = sourcePlasma;
        this.sourceReaction = sourceReaction;

    }


    /**
     * Method for adding a plasma to the thread. Plasmas must be added in order (inside -> outside)
     * @param data Data object that stores everything relevant to the simulation
     */
    void addPlasma(Model.PlasmaData data){
        this.plasmaData.add(data);
    }


    /**
     * Method for setting the detector line of sight (not required)
     * @param detectorLineOfSight Physical location of the tallies (can be null)
     */
    void setDetectorLineOfSight(Vector3D detectorLineOfSight) {
        this.detectorLineOfSight = detectorLineOfSight;
    }


    @Override
    /**
     * Main method of this thread
     */
    public void run() {

        logger.addLog("Setting up tallies ... ");
        setUpTallies();
        logger.logTaskCompletion();

        for (int i = 0; i < numParticles; i++) {

            logger.addLog("Starting particle num " + i);

            // We're assuming particles are always born in the center plasma TODO: Fix this
            int startingLayerIndex = 0;

            // Sample a particle
            Particle particle = sample();

            // Simulate this particle
            traceParticle(particle, startingLayerIndex, sourceParticlePositionTallies,
                    sourceParticleEnergyTallies, sourceParticleTimeTallies);

            logger.logTaskCompletion();
        }
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

            if (detectorLineOfSight != null) {
                bias *= particle.getEnergy();   // TODO: Adhoc correction
            }


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
        double distance = getDistanceToPlasmaBoundary(particle, plasma);\


        // Calculate a step size
        double maxStepSize = distance / MIN_NUM_STEPS;
        double dx = maxStepSize;

        int totalSteps = 0;                             // This is for debugging only
        double particleReactionProb = 0.0;
        double initialEnergy = particle.getEnergy();


        // While the particle is inside this plasma
        while (distance > dx && particle.getEnergy() >  0.0) {

            // Parameters at the starting point
            double r1 = Utils.getNormalizedRadius(particle, plasma);
            double E1 = particle.getEnergy();


            // If the particle is losing a significant amount of energy, we need to take smaller steps
            double dEdx1 = stoppingPowerModel.evaluate(E1, r1);
            dx = Math.min(maxStepSize, -initialEnergy / dEdx1 / MIN_NUM_STEPS);

            // Calculate r2 before the energy loop to save time
            double r2 = Utils.getNormalizedRadius(particle.step(dx, 0.0), plasma);


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
                    r2 = Utils.getNormalizedRadius(particle.step(dx, 0.0), plasma);
                }

                energyError = FastMath.abs(newEnergy - E2) / E2;
                E2 = newEnergy;
                dEdx2 = stoppingPowerModel.evaluate(E2, r2);
            }

            // Loop through all of the reactions relevant to this plasma
            ArrayList<Model.ReactionData> reactions = plasmaData.reactionDataMap.get(particle.getType());
            for (Model.ReactionData data : reactions) {

                NuclearReaction nuclearReaction = data.nuclearReaction;
                ParticleType reactingPlasmaSpecies = data.backgroundParticle;

                // Grab the average number density of the species we're interacting with between these two positions
                double n1 = plasma.getSpeciesNumberDensity(particle.getPosition(), reactingPlasmaSpecies);
                double n2 = plasma.getSpeciesNumberDensity(steppedParticle.getPosition(), reactingPlasmaSpecies);
                double n = 0.5 * (n1 + n2);


                // Grab the average ion temperature between these two positions
                double Ti1 = plasma.getIonTemperature(particle.getPosition());
                double Ti2 = plasma.getIonTemperature(steppedParticle.getPosition());
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
                traceParticle(productParticle, layerIndex, positionTallies, energyTallies, timeTallies);

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

    private double getDistanceToPlasmaBoundary(Particle particle, Plasma plasma){

        // Create a clone that we'll use for the calculation
        Particle tempParticle = particle.clone();

        // We need some kind of guess of what the step size should be to get started
        double plasmaRmax = plasma.getOuterRadiusBound(0.0, 0.0);
        double plasmaRmin = plasma.getInnerRadiusBound(0.0, 0.0);
        double dx = (plasmaRmax - plasmaRmin) / MIN_NUM_STEPS;

        // Step our particle until we're no longer in the plasma
        double distance = 0.0;
        do {
            tempParticle = tempParticle.step(dx, 0.0);
            distance += dx;
        }while (plasma.getIsInside(tempParticle.getPosition()));

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
            if (plasma.getIsInside(tempParticle.getPosition())){
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
        sourceParticlePositionTallies = new Tally[plasmaData.size() + 1];
        for (int i =0; i < sourceParticlePositionTallies.length; i++){
            sourceParticlePositionTallies[i] = new Tally(positionNodes);
        }


        // Build the source time tallies
        double[] timeNodes = Utils.linspace(0.0, 100.0*1e-12, 200);         // TODO HARD CODED!
        sourceParticleTimeTallies = new Tally[plasmaData.size() + 1];
        for (int i =0; i < sourceParticleTimeTallies.length; i++){
            sourceParticleTimeTallies[i] = new Tally(timeNodes);
        }


        // Build the source energy tallies
        double maxSourceEnergy = 2.0*sourceReaction.getZeroTemperatureMeanEnergy();     // TODO: HARD CODED!
        int numSourceNodes = (int) Math.ceil(maxSourceEnergy / ENERGY_NODE_WIDTH);
        double[] sourceNodes = Utils.linspace(0.0, maxSourceEnergy, numSourceNodes);
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
                    double[] productNodes = Utils.linspace(0.0, maxProductEnergy, numProductNodes);

                    Tally[] productParticleEnergyTallies = new Tally[plasmaData.size() + 1];
                    for (int i = 0; i < productParticleEnergyTallies.length; i++){
                        productParticleEnergyTallies[i] = new Tally(productNodes);
                    }
                    productParticleEnergyTallyMap.put(nuclearReaction, productParticleEnergyTallies);

                }
            }
        }
    }

    private Particle sample(){

        // Sample the direction of our position vector
        // TODO: I don't think this is valid for a non-symmetric plasma...
        Vector3D position = Utils.sampleRandomNormalizedVector();

        // Convert it to spherical
        SphericalCoordinates coordinates = new SphericalCoordinates(position);

        // Apache is in Math notation whilst we're in Physics notation
        double theta = coordinates.getPhi();
        double phi = coordinates.getTheta();

        // Use our radial distribution and the plasma bounds to determine the magnitude of our position vector
        double rNorm = radialSourceDistribution.sample();
        double rMax  = sourcePlasma.getOuterRadiusBound(theta, phi);
        double rMin  = sourcePlasma.getInnerRadiusBound(theta, phi);

        double r = (rMax - rMin) * rNorm + rMin;

        // Rescale our position vector to this magnitude
        position = position.scalarMultiply(r);

        // Sample the direction
        Vector3D direction = Utils.sampleRandomNormalizedVector();

        // Build the energy distribution based on our position in the plasma
        double mu = sourceReaction.getZeroTemperatureMeanEnergy();      // TODO: Include temperature upshift...
        double sigma = sourcePlasma.getThermalSigma(position, sourceReaction);

        // TODO: Very temp...
        /*
        double Tion = sourcePlasma.getIonTemperature(position);
        double deltaE = -9.5302E-04*Math.pow(Tion, 4) +
                4.2759E-02*Math.pow(Tion, 3) -
                6.8887E-01*Math.pow(Tion, 2) +
                8.6333E+00*Math.pow(Tion, 1);
         */


        // Sample the energy
        double energy = Distribution.normDistribution(mu, sigma).sample();

        // Get the particle type
        ParticleType type = sourceReaction.getProducts()[0];

        // Return the sample particle
        Particle particle = new Particle(type, position, direction, energy, 0.0);
        particle.setWeight(1.0 / numParticles);     // TODO: Use actual XS?
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
