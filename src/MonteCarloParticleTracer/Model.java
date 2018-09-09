package MonteCarloParticleTracer;

import org.apache.commons.math3.geometry.euclidean.threed.Vector3D;

import java.util.ArrayList;
import java.util.HashMap;

/**
 * Created by lahmann on 2016-06-19.
 *
 * Class the holds all of the information required to sample particles from some ParticleDistribution
 * and track them through some Plasma. The actual step by step logic of the simulation is handled within
 * individual ParticleHistoryTasks, this class just initiates them
 *
 */


public class Model {

    // File that all of the output data will be dumped too
    private OutputFile outputFile;

    // Nuclear reaction from which the source particles will originate
    private NuclearReaction sourceNuclearReaction;

    // Reactivity of the source particle (ALTERNATIVE WAY TO SET THE RADIAL DISTRIBUTION)
    private Reactivity sourceReactivity;

    // Radial distribution from which to sample the source particle (ALTERNATIVE TO REACTIVITY METHOD)
    private Distribution sourceRadialDistribution;

    // An ORDERED array list of all the plasmas in this model
    private ArrayList<Plasma> plasmas = new ArrayList<>();

    // An array list of all nuclear reactions we'll simulate in this model
    private ArrayList<NuclearReaction> nuclearReactions = new ArrayList<>();

    // The vector that all source particles will be forced to travel
    // Secondary reactions cannot be trusted when this is used
    private Vector3D sourceDirection = null;

    // The vector that all secondary reactions will be forced to travel
    private Vector3D productDirection = null;

    // Threads that this model interacts with
    private SimulationThread[] threads;

    // Logger object that will keep track of progress / errors
    private Logger logger = new Logger(false);

    // DoubleArray of tallies for the source particles (1 tally for birth, and 1 for every surface)
    Tally[] sourceParticlePositionTallies = null;
    Tally[] sourceParticleEnergyTallies = null;
    Tally[] sourceParticleTimeTallies = null;

    // HashMap of tally arrays for every nuclear reaction we model
    HashMap<NuclearReaction, Tally[]> productParticlePositionTallyMap = null;
    HashMap<NuclearReaction, Tally[]> productParticleEnergyTallyMap = null;
    HashMap<NuclearReaction, Tally[]> productParticleTimeTallyMap = null;



    /**
     * Default constructor
     * @param name Filename for the file where the output data will be saved
     */
    public Model(String name) {
        this.outputFile = new OutputFile(name);
    }



    /**
     * Method to set all of the required source information. This method much be called to have a valid model
     * @param sourceNuclearReaction Nuclear reaction of the source (primary) particle
     * @param sourceReactivity Reactivity of the source particle reaction
     */
    public void setSourceInformation(NuclearReaction sourceNuclearReaction, Reactivity sourceReactivity) {
        this.sourceNuclearReaction = sourceNuclearReaction;
        this.sourceReactivity = sourceReactivity;
        this.sourceRadialDistribution = null;
    }

    public void setSourceInformation(NuclearReaction sourceNuclearReaction, Distribution radialDistribution) {
        this.sourceNuclearReaction = sourceNuclearReaction;
        this.sourceRadialDistribution = radialDistribution;
    }


    /**
     * Method to add a plasma plasma to the model. Plasma layers should be added in order from inside to outside.
     * At least one plasma plasma needs to be added for a valid model.
     * @param plasma Plasma object to add to this model
     */
    public void addPlasmaLayer(Plasma plasma){

        // Force the new plasma to share a boundary with the out most plasma
        if (plasmas.size() > 0){
            Plasma outerMostLayer = plasmas.get(plasmas.size() - 1);
            plasma.setInnerBoundaryLegendreModes(outerMostLayer.getOuterBoundaryLegendreModes());
        }

        plasmas.add(plasma);
    }


    /**
     * Method to add a nuclear reaction to tally. No nuclear reactions are required for a working model.
     * @param reaction Nuclear reaction to add to this model
     */
    public void addNuclearReaction(NuclearReaction reaction){
        nuclearReactions.add(reaction);
    }


    /**
     * Method to force all source particles to go in a particular direction
     * @param sourceDirection direction that ALL source particles will travel
     */
    public void setSourceDirection(Vector3D sourceDirection) {
        this.sourceDirection = sourceDirection;
    }

    /**
     * Method to force all product particles to go in a particular direction
     * @param productDirection direction that ALL product particles will travel
     */
    public void setProductDirection(Vector3D productDirection) {
        this.productDirection = productDirection;
    }


    /**
     * Main method that sets up all of the individual Tasks and collects the tallys upon completion
     * @param totalParticles The total number of source particles to be sampled
     * @param numCPUs The number of nodes to spread the calculation across
     * @throws Exception
     */
    public void runSimulation(int totalParticles, int numCPUs) {

        try {

            /**
             * Verification and data organization
             */

            // Verify that we have a plasma
            if (plasmas.size() == 0)
                throw new Exceptions.NoPlasmaSpecifiedException();


            // Verify that the source reaction has been set
            if (sourceNuclearReaction == null)
                throw new Exceptions.NoSourceInformationException();

            // Verify that the radial distribution has been set
            if (sourceRadialDistribution == null){

                if (sourceReactivity != null){
                    sourceRadialDistribution = plasmas.get(0).getRadialBurnDistribution(sourceReactivity);
                }
                else {
                    throw new Exceptions.NoSourceInformationException();
                }
            }

            // Verify that no nuclear reactions are being used if we're forcing the source direction
            if (sourceDirection != null && nuclearReactions.size() > 0)
                throw new Exceptions.InvalidNuclearReactionModelException();


            // Verify this machine has enough processors
            if (numCPUs > Runtime.getRuntime().availableProcessors() - 1 || numCPUs < 1)
                throw new Exceptions.InvalidNumberProcessorsRequestedException(numCPUs, Runtime.getRuntime().availableProcessors() - 1);


            // Load the stopping power libraries
            String osName = System.getProperty("os.name");
            if (osName.contains("Windows")) {
                System.load(DataFiles.windows_StopPow_Lib.getAbsolutePath());
            }
            else if (osName.contains("Linux")) {
                System.load(DataFiles.linux_StopPow_Lib.getAbsolutePath());
            }
            else {
                throw new Exceptions.UnsupportedOperatingSystemException(osName);
            }


            // Calculate the number of particles we'll simulate in each thread
            int particlesPerThread = Math.floorDiv(totalParticles, numCPUs) + 1;


            // Organize the all of the plasma data to save runtime
            logger.addLog("Starting model setup ... ");
            PlasmaData[] plasmaData = new PlasmaData[plasmas.size()];
            for (int i = 0; i < plasmaData.length; i++) {
                plasmaData[i] = new PlasmaData(plasmas.get(i));
            }
            logger.logTaskCompletion();


            // Clear out the tallies (in case this simulation is being re-ran)
            sourceParticlePositionTallies = null;
            sourceParticleEnergyTallies = null;
            sourceParticleTimeTallies = null;
            productParticlePositionTallyMap = null;
            productParticleEnergyTallyMap = null;
            productParticleTimeTallyMap = null;



            /**
             * Running the actual simulation
             */

            // Set up the simulation threads
            logger.addLog(String.format("Building %d threads with %d particles each ... ", numCPUs, particlesPerThread));
            threads = new SimulationThread[numCPUs];
            for (int i = 0; i < threads.length; i++) {

                // Make the thread
                threads[i] = new SimulationThread(particlesPerThread, numCPUs*particlesPerThread);
                logger.addLog("Started thread " + threads[i].getId());

                // Add the source and tally info
                threads[i].setSourceInformation(plasmas.get(0), sourceNuclearReaction, sourceRadialDistribution);
                threads[i].setSourceDirection(sourceDirection);
                threads[i].setProductDirection(productDirection);

                // Add the plasma data
                for (PlasmaData data : plasmaData) {
                    threads[i].addPlasma(data);
                }

                // Set up the tallies internally
                threads[i].setUpTallies();

                // Start the thread
                threads[i].start();
            }
            logger.logTaskCompletion();



            /**
             * Post-simulation clean-up
             */

            // Rejoin with the threads and sum the tallies
            logger.addLog("Waiting on threads .... ");
            for (SimulationThread thread : threads) {

                // Join
                thread.join();
                logger.addLog("Rejoined with thread " + thread.getId());

                // This means this is our first thread
                if (sourceParticlePositionTallies == null) {

                    sourceParticlePositionTallies = thread.getSourceParticlePositionTallies();
                    sourceParticleEnergyTallies = thread.getSourceParticleEnergyTallies();
                    sourceParticleTimeTallies = thread.getSourceParticleTimeTallies();

                    productParticlePositionTallyMap = thread.getProductParticlePositionTallyMap();
                    productParticleEnergyTallyMap = thread.getProductParticleEnergyTallyMap();
                    productParticleTimeTallyMap = thread.getProductParticleTimeTallyMap();
                }

                // Otherwise we need to start adding them
                else {

                    // Source position tally
                    Tally[] threadTallies = thread.getSourceParticlePositionTallies();
                    for (int i = 0; i < threadTallies.length; i++) {
                        sourceParticlePositionTallies[i].addTally(threadTallies[i]);
                    }

                    // Source energy tally
                    threadTallies = thread.getSourceParticleEnergyTallies();
                    for (int i = 0; i < threadTallies.length; i++) {
                        sourceParticleEnergyTallies[i].addTally(threadTallies[i]);
                    }

                    // Source time tally
                    threadTallies = thread.getSourceParticleTimeTallies();
                    for (int i = 0; i < threadTallies.length; i++) {
                        sourceParticleTimeTallies[i].addTally(threadTallies[i]);
                    }

                    for (NuclearReaction key : productParticlePositionTallyMap.keySet()) {

                        // Product position tally
                        threadTallies = thread.getProductParticlePositionTallyMap().get(key);
                        for (int i = 0; i < threadTallies.length; i++) {
                            productParticlePositionTallyMap.get(key)[i].addTally(threadTallies[i]);
                        }

                        // Product energy tally
                        threadTallies = thread.getProductParticleEnergyTallyMap().get(key);
                        for (int i = 0; i < threadTallies.length; i++) {
                            productParticleEnergyTallyMap.get(key)[i].addTally(threadTallies[i]);
                        }

                        // Product energy tally
                        threadTallies = thread.getProductParticleTimeTallyMap().get(key);
                        for (int i = 0; i < threadTallies.length; i++) {
                            productParticleTimeTallyMap.get(key)[i].addTally(threadTallies[i]);
                        }
                    }
                }
            }
            logger.logTaskCompletion();


            // Finally try to write the output file
            writeOutputFile();

        }

        catch (Exception e){
            e.printStackTrace();
        }

        finally {

        }

    }

    private void writeOutputFile() throws Exception {

        // Clear any previous contents
        outputFile.clear();

        // Add the log info
        outputFile.addHeader("Logs");
        outputFile.addString(logger.toString());


        // Add all of the plasma data
        outputFile.addHeader("Simulated Plasma Data");
        for (int i = 0; i < plasmas.size(); i++){
            outputFile.addString("Plasma Layer " + i);
            outputFile.addString(plasmas.get(i).toString());
        }


        // Add source information to the output file
        outputFile.addHeader("Source Information");
        outputFile.addString(sourceNuclearReaction.toString());


        // Add all of the modeled nuclear reactions to the output file
        outputFile.addHeader("Modelled Nuclear Reactions");
        for (int i = 0; i < nuclearReactions.size(); i++){
            outputFile.addString("Reaction " + i);
            outputFile.addString(nuclearReactions.get(i).toString());
        }


        // Add the source particle tallies to the output
        for (int i = 0; i < sourceParticlePositionTallies.length; i++) {

            // Birth tallies
            if (i == 0) {

                outputFile.addHeader("Source Particle Birth Tallies");

                outputFile.addSubheader("Source Birth Position Tally");
                outputFile.addString(sourceParticlePositionTallies[0].toString());

                outputFile.addSubheader("Source Birth Spectrum Tally");
                outputFile.addString(sourceParticleEnergyTallies[0].toString());

                outputFile.addSubheader("Source Birth Time Tally");
                outputFile.addString(sourceParticleTimeTallies[0].toString());

            }

            // Surface tallies
            else {

                outputFile.addHeader("Source Particle Surface Tallies (Plasma Layer " + (i - 1) + ")");

                outputFile.addSubheader("Source Position Surface Tally");
                outputFile.addString(sourceParticlePositionTallies[i].toString());

                outputFile.addSubheader("Source Spectrum Surface Tally");
                outputFile.addString(sourceParticleEnergyTallies[i].toString());

                outputFile.addSubheader("Source Time Surface Tally");
                outputFile.addString(sourceParticleTimeTallies[i].toString());

            }
        }


        // Add the product particle tallies to the output
        for (NuclearReaction key : productParticlePositionTallyMap.keySet()) {

            Tally[] productParticlePositionTallies = productParticlePositionTallyMap.get(key);
            Tally[] productParticleEnergyTallies   = productParticleEnergyTallyMap  .get(key);
            Tally[] productParticleTimeTallies     = productParticleTimeTallyMap    .get(key);

            for (int i = 0; i < productParticlePositionTallies.length; i++) {

                // Birth tallies
                if (i == 0) {

                    outputFile.addHeader("Particle Birth Tallies");
                    outputFile.addString(key.toString());

                    outputFile.addSubheader("Birth Position Tally");
                    outputFile.addString(productParticlePositionTallies[0].toString());

                    outputFile.addSubheader("Birth Spectrum Tally");
                    outputFile.addString(productParticleEnergyTallies[0].toString());

                    outputFile.addSubheader("Birth Time Tally");
                    outputFile.addString(productParticleTimeTallies[0].toString());

                }

                // Surface tallies
                else {

                    outputFile.addHeader("Particle Surface Tallies (Plasma Layer " + (i - 1) + ")");
                    outputFile.addString(key.toString());

                    outputFile.addSubheader("Position Surface Tally");
                    outputFile.addString(productParticlePositionTallies[i].toString());

                    outputFile.addSubheader("Spectrum Surface Tally");
                    outputFile.addString(productParticleEnergyTallies[i].toString());

                    outputFile.addSubheader("Time Surface Tally");
                    outputFile.addString(productParticleTimeTallies[i].toString());

                }
            }
        }



    }

    public OutputFile getOutputFile() {
        return outputFile;
    }

    class PlasmaData {

        Plasma plasma;
        HashMap<ParticleType, StoppingPowerModel> stoppingPowerModels = new HashMap<>();
        HashMap<ParticleType, ArrayList<ReactionData>> reactionDataMap = new HashMap<>();


        PlasmaData(Plasma plasma) {

            // Store the plasma
            this.plasma = plasma;

            // Determine the source particle
            ParticleType sourceParticle = sourceNuclearReaction.getProducts()[0];

            // Build the stopping power model for the source particle
            stoppingPowerModels.put(sourceParticle, new StoppingPowerModel(sourceParticle, plasma));

            // Init the reaction data map with an empty array
            reactionDataMap.put(sourceParticle, new ArrayList<>());


            // Loop through all of the nuclear reactions
            for (NuclearReaction reaction : nuclearReactions){

                // Copy the reaction
                reaction = reaction.copy();

                // Get the product particle associated with this reaction
                ParticleType productParticle = reaction.getProducts()[0];

                // Build the stopping power model for this particle
                stoppingPowerModels.put(productParticle, new StoppingPowerModel(productParticle, plasma));

                // Determine if our source particle can cause this reaction in this plasma
                ParticleType A = reaction.getReactants()[0];
                ParticleType B = reaction.getReactants()[1];

                if (sourceParticle.equals(A) && plasma.containsSpecies(B)){
                    reactionDataMap.get(sourceParticle).add(new ReactionData(reaction, B));
                }

                if (sourceParticle.equals(B) && plasma.containsSpecies(A)){
                    reactionDataMap.get(sourceParticle).add(new ReactionData(reaction, A));
                }
            }
        }


        PlasmaData copy(){
            return new PlasmaData(plasma.copy());
        }
    }

    class ReactionData{

        NuclearReaction nuclearReaction;
        ParticleType backgroundParticle;

        public ReactionData(NuclearReaction nuclearReaction, ParticleType backgroundParticle) {
            this.nuclearReaction = nuclearReaction;
            this.backgroundParticle = backgroundParticle;
        }

        public ReactionData copy(){
            NuclearReaction nuclearReaction = this.nuclearReaction.copy();
            ParticleType backgroundParticle = this.backgroundParticle.copy();
            return new ReactionData(nuclearReaction, backgroundParticle);
        }
    }

}
